#include <mc_rtc/gui/Form.h>

#include <BaselineWalkingController/BaselineWalkingController.h>
#include <BaselineWalkingController/FootManager.h>
#include <BaselineWalkingController/CentroidalManager.h>
#include <BaselineWalkingController/SwingTraj.h>
#include <BaselineWalkingController/states/SoftFootState.h>
#include <BaselineWalkingController/swing/SwingTrajLandingSearch.h>

#include <variable_stiffness/connectionFile.h>

#include <numeric>
#include <functional>
#include <cmath>

#include <mc_tasks/CoMTask.h>
#include <mc_mujoco/devices/RangeSensor.h>

#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/QhullPoint.h"

using namespace BWC;

void SoftFootState::start(mc_control::fsm::Controller & _ctl)
{
  State::start(_ctl);

  // Create client
  client_ = mc_rtc::ROSBridge::get_node_handle()->serviceClient<variable_stiffness::connectionFile>("connectionFile");

  // Initialize map
  foot_data_[Foot::Left] = FootData{};
  foot_data_[Foot::Right] = FootData{};

  _ctl.gui()->addElement({"SoftFoot"},
    mc_rtc::gui::Label("cost", [this]() { return this->cost_; })
  );
  
  _ctl.logger().addLogEntry("cost", [this]() { return cost_; });

  _ctl.gui()->addElement({"SoftFoot"}, mc_rtc::gui::Label("PhalangesStiffness", [this]() { return this->PhalangesStiffness_; }));
  _ctl.logger().addLogEntry("PhalangesStiffness", [this]() { return PhalangesStiffness_; });

  output("OK");
}

bool SoftFootState::run(mc_control::fsm::Controller & ctl)
{
  // Cast ctl to BaselineWalkingController
  auto & ctrl = static_cast<BWC::BaselineWalkingController&>(ctl);

  // Compute cost
  calculateCost(ctl);

  // Check if we are in single support or not
  if(ctrl.footManager_->supportPhase() == BWC::SupportPhase::DoubleSupport)
  {
    foot_data_[Foot::Right].needReset = true;
    foot_data_[Foot::Left].needReset = true; 
    return false;
  }
  // Get the current moving foot
  Foot current_moving_foot;
  if(ctrl.footManager_->supportPhase() == BWC::SupportPhase::LeftSupport) // This is the name "LeftFootCenter" or "RightFootCenter"
  {
    current_moving_foot = Foot::Right;
  }
  else
  {
    current_moving_foot = Foot::Left;
  }
  
  if(foot_data_[current_moving_foot].needReset)
  {
    reset(ctl, current_moving_foot);
  }

  // Estimate ground from sensors
  estimateGround(ctl, current_moving_foot);

  // Check foot position with respect to desired landing pose
  const auto & ground = foot_data_[current_moving_foot].ground;
  // Get landing: a bit dirty for the moment
  const sva::PTransformd & X_0_landing = ctrl.footManager_->swingTraj_->endPose_;
  // If we saw more than the landing pose + half of the foot, we have enough data to perform all the computations
  if(!foot_data_[current_moving_foot].areComputationDone && ground.back().x() >= X_0_landing.translation().x() + foot_length_ * 0.5)
  {
    mc_rtc::log::success("Accumulated enough data");
    // We need to do these steps only one time
    foot_data_[current_moving_foot].areComputationDone = true;
    // Extract ground segment from ground data
    extractGroundSegment(ctl, current_moving_foot, X_0_landing.translation());
    // Continue only if the sensor measured more than 10 values
    if(ground_segment_[current_moving_foot].raw.size() > 5)
    {
      // Compute the altitude profile
      extractAltitudeProfileFromGroundSegment(current_moving_foot);
      // call the server and update the variable stiffness
      updateVariableStiffness(ctl, current_moving_foot);
      // Compute convex hull of the segment -> right now it does not work 
      computeSegmentConvexHull(ctl, current_moving_foot);
      // Now with the convex hull we can compute the angle
      computeFootLandingAngle(current_moving_foot, X_0_landing.translation());
      // Update targeted pose
      updateFootSwingPose(ctl, current_moving_foot, X_0_landing);
    }
  }

  return false;
}

void SoftFootState::teardown(mc_control::fsm::Controller &)
{
  // Clean up GUI
  ctl().gui()->removeCategory({ctl().name(), "SoftFoot"});
  
  ctl().logger().removeLogEntry("cost");
  ctl().logger().removeLogEntry("PhalangesStiffness");
}

void SoftFootState::calculateCost(mc_control::fsm::Controller & ctl)
{
  // Cast ctl to BaselineWalkingController
  auto & ctrl = static_cast<BWC::BaselineWalkingController&>(ctl);

  double zmp = 0.0;
  // ZMP Error
  {
    // Get ZMP from Centroidal Manager
    const Eigen::Vector3d & zmp_ref = ctrl.centroidalManager_->refZmp_;
    // Compute Measured ZMP
    std::unordered_map<Foot, sva::ForceVecd> sensorWrenchList;
    for(const auto & foot : ctrl.footManager_->getCurrentContactFeet())
    {
      const auto & surfaceName = ctrl.footManager_->surfaceName(foot);
      const auto & sensorName = ctrl.robot().indirectSurfaceForceSensor(surfaceName).name();
      const auto & sensor = ctrl.robot().forceSensor(sensorName);
      const auto & sensorWrench = sensor.worldWrenchWithoutGravity(ctrl.robot());
      sensorWrenchList.emplace(foot, sensorWrench);
    }
    const Eigen::Vector3d & zmp_mes = ctrl.centroidalManager_->calcZmp(sensorWrenchList, zmp_ref.z());
    const double zmp_error = (zmp_ref - zmp_mes).norm();

    zmp_.push_back(zmp_error);
    double sum = std::accumulate(zmp_.begin(), zmp_.end(), 0.0);
    double mean = sum / zmp_.size();
    std::vector<double> diff(zmp_.size());
    std::transform(zmp_.begin(), zmp_.end(), diff.begin(), [mean](double x) { return x - mean; });
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / zmp_.size());
    zmp = lambda_zmp_ * mean + lambda_zmp_ * lambda_zmp_ * stdev;
  }

  double CoM = 0.0;
  // CoM error
  {
    const Eigen::Vector3d & CoM_ref = ctrl.comTask_->com();
    const Eigen::Vector3d & CoM_mes = ctrl.realRobot().com();
    const double CoM_error = (CoM_ref - CoM_mes).norm();

    CoM_.push_back(CoM_error);
    double sum = std::accumulate(CoM_.begin(), CoM_.end(), 0.0);
    double mean = sum / CoM_.size();
    std::vector<double> diff(CoM_.size());
    std::transform(CoM_.begin(), CoM_.end(), diff.begin(), [mean](double x) { return x - mean; });
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / CoM_.size());
    CoM = lambda_CoM_ * mean + lambda_CoM_ * lambda_CoM_ * stdev;
  }

  double footstep = 0.0;
  // Footstep
  {
    // Current number of steps to do
    footstep_error_ =  ctrl.footManager_->footstepQueue().size() ;
    footstep = lambda_footstep_ * footstep_error_;
  } 
  
  cost_ = - footstep - zmp - CoM;
  //mc_rtc::log::success("Foot_Target_Pose : {}, Foot_Robot_Pose : {}, error_distance : {}", Foot_Target_Pose.translation().head<2>(), Foot_Robot_Pose.translation().head<2>(), error_distance);
}

void SoftFootState::estimateGround(mc_control::fsm::Controller & ctl, const Foot & current_moving_foot)
{
  // Cast ctl to BaselineWalkingController
  auto & ctrl = static_cast<BWC::BaselineWalkingController&>(ctl);

  // Select data/string based on current_moving_foot
  std::string sensor_name = range_sensor_name_[current_moving_foot];

  // From here do not need to worry about which foot it is
  FootData & data = foot_data_[current_moving_foot];

  // Return the parent body of the sensor (phalanx)
  const std::string& BodyOfSensor = ctl.robot().device<mc_mujoco::RangeSensor>(sensor_name).parent(); 
  // Access the position of body name in world coordinates (phalanx position)
  sva::PTransformd X_0_ph = ctl.realRobot().bodyPosW(BodyOfSensor); 
  // Returns the transformation from the parent body to the sensor
  const sva::PTransformd& X_ph_s = ctl.robot().device<mc_mujoco::RangeSensor>(sensor_name).X_p_s();
  // Sensor position in global frame Z coordinate
  data.range = ctl.robot().device<mc_mujoco::RangeSensor>(sensor_name).data(); 
  const sva::PTransformd X_s_m = sva::PTransformd(Eigen::Vector3d(0, 0, data.range));
  sva::PTransformd X_0_m = X_s_m*X_ph_s*X_0_ph;
  // Keep the estimated 3d point for the ground
  data.ground.push_back(X_0_m.translation());
}

void SoftFootState::extractGroundSegment(mc_control::fsm::Controller & ctl, const Foot & current_moving_foot, const Eigen::Vector3d & landing)
{ 
  FootData & data = foot_data_[current_moving_foot];
  // Do a copy to sort the vector
  std::vector<Eigen::Vector3d> ground = data.ground;
  // Sort alongside x
  std::sort(ground.begin(), ground.end(), 
    [](const Eigen::Vector3d & a, const Eigen::Vector3d & b)
    {
      return a.x() < b.x();
    }
  );

  // Find beginning and ending of segment to extract
  const auto begin_iterator = std::find_if(ground.begin(), ground.end(),
    [&](const Eigen::Vector3d & v) { return v.x() >= landing.x() - foot_length_ * 0.5; });
  const auto end_iterator = std::find_if(ground.begin(), ground.end(),
    [&](const Eigen::Vector3d & v) { return v.x() >= landing.x() + foot_length_ * 0.5; });
  
  // Save the selected segment in raw data of ground segment structure
  // ground_segment_[current_moving_foot].raw.clear();
  std::transform(begin_iterator, end_iterator, std::back_inserter(ground_segment_[current_moving_foot].raw),
    [](const Eigen::Vector3d & v) { return v; }); 

  mc_rtc::log::error("ground_segment_[{}].raw {}", current_moving_foot == Foot::Left ? "Left" : "Right", ground_segment_[current_moving_foot].raw.size());

  // Add trajectory
  if(ground_segment_[current_moving_foot].raw.size() > 5)
  {
    const std::string name = current_moving_foot == Foot::Left ? "left" : "right";  
    ctl.gui()->addElement({"SoftFoot"},
      mc_rtc::gui::Trajectory(name + "_segment", {mc_rtc::gui::Color::Red}, [this, current_moving_foot]() { return ground_segment_[current_moving_foot].raw; })
    );
  }
}

void SoftFootState::extractAltitudeProfileFromGroundSegment(const Foot & current_moving_foot)
{
  // Get the segment
  const auto & raw_segment = ground_segment_[current_moving_foot].raw;
  // Save the selected segment in raw data of ground segment structure
  std::transform(raw_segment.begin(), raw_segment.end(), std::back_inserter(foot_data_[current_moving_foot].altitude),
    [](const Eigen::Vector3d & v) { return v.z(); }); 

  mc_rtc::log::info("extractAltitudeProfileFromGroundSegment");
}

void SoftFootState::updateVariableStiffness(mc_control::fsm::Controller & ctl, const Foot & current_moving_foot)
{
  variable_stiffness::connectionFile srv;
  srv.request.profile = foot_data_[current_moving_foot].altitude;

  if (client_.call(srv))
  {
    mc_rtc::log::success("Stiffness: {}", srv.response.stiffness);
    foot_data_[current_moving_foot].k = srv.response.stiffness;
    // Solution to modify the variable stiffness
    auto stiffnessToAngle = [this](double VarStiff) 
      {
        double angle_low = 0;
        double angle_high = 1;
        double stiffness_low = 0;
        double stiffness_high = 100;
        return angle_low+(VarStiff-stiffness_low)*(angle_high-angle_low)/(stiffness_high-stiffness_low);
      };
    auto postureTask = ctl.getPostureTask(ctl.robot().name());
    // Reset stifness for both feet gains
    postureTask->jointGains(ctl.solver(), {tasks::qp::JointGains("R_VARSTIFF", 350), tasks::qp::JointGains("L_VARSTIFF", 350)});
    // Set computed stiffness for current moving foot
    postureTask->target({{variable_stiffness_jointname_[current_moving_foot], std::vector<double>{stiffnessToAngle(foot_data_[current_moving_foot].k)}}});
  }
  else
  {
    mc_rtc::log::error("Failed to call service connectionFile to compute the variable stiffness");
  }

  mc_rtc::log::info("updateVariableStiffness");
}

void SoftFootState::computeSegmentConvexHull(mc_control::fsm::Controller & ctl, const Foot & current_moving_foot)
{ 
  const auto & raw_segment = ground_segment_[current_moving_foot].raw;

  // Build the input for qhull
  std::vector<double> points_in;
  double min = std::numeric_limits<double>::max();
  double max = std::numeric_limits<double>::min();
  for(size_t i = 0; i < raw_segment.size(); i+=5)
  {
    const auto& p = raw_segment[i];
    points_in.push_back(p.x());
    points_in.push_back(p.z());

    min = std::min(min, p.z());
    max = std::max(max, p.z());
  }

  std::vector<Eigen::Vector2d> convex_hull;
  if(max - min > 0.001)
  {
    try
    {
      // Run qhull
      orgQhull::Qhull qhull;
      qhull.runQhull("", 2, points_in.size() / 2, points_in.data(), "QJ");

      orgQhull::QhullFacet face = qhull.firstFacet();
      orgQhull::QhullVertex v;
      size_t prev_id = std::numeric_limits<size_t>::max();

      for(size_t i = 0; i < qhull.facetCount(); ++i)
      {
        face = face.nextFacet2d(&v);
        if(v.point().id() < prev_id)
        {
          convex_hull.emplace_back(v.point().coordinates()[0], v.point().coordinates()[1]);
        }
        prev_id = v.point().id();
      }

      // Sort alongside x
      std::sort(convex_hull.begin(), convex_hull.end(), [](const Eigen::Vector2d& a, const Eigen::Vector2d& b) { return a.x() < b.x(); });

    }
    catch (std::exception& e)
    {
      mc_rtc::log::error("Error during qhull ! {}", e.what());
      convex_hull.push_back(Eigen::Vector2d(raw_segment.front().x(), raw_segment.front().z()));
      convex_hull.push_back(Eigen::Vector2d(raw_segment.back().x(), raw_segment.back().z()));
    }
  }
  else
  {
    mc_rtc::log::error("Too flat to run qhull !");
    convex_hull.push_back(Eigen::Vector2d(raw_segment.front().x(), raw_segment.front().z()));
    convex_hull.push_back(Eigen::Vector2d(raw_segment.back().x(), raw_segment.back().z()));
  }

  // For display
  std::vector<Eigen::Vector3d> convex_display;
  for(size_t i = 0; i < convex_hull.size(); ++i)
  {
    convex_display.push_back(Eigen::Vector3d(convex_hull[i].x(), 0., convex_hull[i].y()));
  }

  // Add trajectory
  const std::string name = current_moving_foot == Foot::Left ? "left" : "right";  
  ctl.gui()->addElement({"SoftFoot"},
    mc_rtc::gui::Trajectory(name + "_convex_display", {mc_rtc::gui::Color::Blue}, [this, convex_display]() { return convex_display; })
  );
 
  // Interpolate points in-between point from convex hull
  auto lerp = [](const Eigen::Vector2d & A, const Eigen::Vector2d& B, double y, double t)
  {
    const Eigen::Vector2d vec = B * t + A * (1. - t);
    return Eigen::Vector3d(vec.x(), y, vec.y());
  };

  ground_segment_[current_moving_foot].convex.clear();
  for(size_t i = 0; i < convex_hull.size() - 1; ++i)
  {
    for(double t = 0.; t <= 1.0; t += 0.1)
    {
      ground_segment_[current_moving_foot].convex.push_back(lerp(convex_hull[i], convex_hull[i+1], raw_segment.front().y(), t));
    }
  }

  // Add trajectory
  ctl.gui()->addElement({"SoftFoot"},
    mc_rtc::gui::Trajectory(name + "_convex_segment", {mc_rtc::gui::Color::Magenta}, [this, current_moving_foot]() { return ground_segment_[current_moving_foot].convex; })
  );
}

void SoftFootState::computeFootLandingAngle(const Foot & current_moving_foot, const Eigen::Vector3d & landing)
{
  // Get closest point from landing in convex
  const auto & convex = ground_segment_[current_moving_foot].convex;  
  auto convex_iterator = std::find_if(convex.begin(), convex.end(), [&](const Eigen::Vector3d & v) { return v.x() >= landing.x(); });

  const Eigen::Vector3d & p_1 = *convex_iterator;
  convex_iterator = --convex_iterator;
  const Eigen::Vector3d & p_0 = *(convex_iterator);

  double dz = p_1.z() - p_0.z();
  double dx = p_1.x() - p_0.x();

  foot_data_[current_moving_foot].angle = -std::atan(dz / dx);

  mc_rtc::log::info("angle {} [rad] {} [deg]", foot_data_[current_moving_foot].angle, foot_data_[current_moving_foot].angle * 180. / M_PI);
}

void SoftFootState::updateFootSwingPose(mc_control::fsm::Controller & ctl, const Foot & current_moving_foot, const sva::PTransformd & X_0_landing)
{
  // Cast ctl to BaselineWalkingController
  auto & ctrl = static_cast<BWC::BaselineWalkingController&>(ctl);
  // Get the angle
  double desired_angle = foot_data_[current_moving_foot].angle;
  // Dirty access to swing traj
  
  // Follow Murooka-san code to update, sorry for the dirty access
  dynamic_cast<BWC::SwingTrajLandingSearch*>(ctrl.footManager_->swingTraj_.get())->updatePitch(desired_angle);
}


void SoftFootState::reset(mc_control::fsm::Controller & ctl, const Foot & foot)
{
  foot_data_[foot].needReset = false;
  // Delete all the data that are behind the foot
  const auto X_0_p = ctl.robot().surfacePose(surface_name_[foot]);
  auto & ground = foot_data_[foot].ground;
  const auto ground_iterator = std::find_if(ground.begin(), ground.end(), [&](const Eigen::Vector3d & v) { return v.x() >= X_0_p.translation().x() - 0.5 * foot_length_; });
  // Delete from the beginning of the vector up to the back part of the foot
  ground.erase(ground.begin(), ground_iterator);
  
  // Reset altitude
  foot_data_[foot].altitude.clear();

  // Reset boolean
  foot_data_[foot].areComputationDone = false;

  // Reset ground segment data
  ground_segment_[foot].raw.clear();
  ground_segment_[foot].convex.clear();

  // Get name for the feet
  const Foot other_foot = foot == Foot::Left ? Foot::Right : Foot::Left;  
  const std::string other_name = other_foot == Foot::Left ? "left" : "right";
  const std::string name = foot == Foot::Left ? "left" : "right";  
  // Reset logger
  ctl.logger().removeLogEntry("MyMeasures_" + other_name + "_range");
  ctl.logger().removeLogEntry("MyMeasures_" + other_name + "_ground");

  // Reset GUI
  ctl.gui()->removeElement({"SoftFoot"}, other_name + "_point_ground");
  ctl.gui()->removeElement({"SoftFoot"}, other_name + "_ground");

  ctl.gui()->removeElement({"SoftFoot"}, other_name + "_segment");
  ctl.gui()->removeElement({"SoftFoot"}, name + "_segment");

  ctl.gui()->removeElement({"SoftFoot"}, other_name + "_convex_segment");
  ctl.gui()->removeElement({"SoftFoot"}, name + "_convex_segment");

  ctl.gui()->removeElement({"SoftFoot"}, other_name + "_convex_display");
  ctl.gui()->removeElement({"SoftFoot"}, name + "_convex_display");

  // Handle logger and gui
  ctl.logger().addLogEntry("MyMeasures_" + name + "_ground", [this, foot]() { return foot_data_[foot].ground.back();} );
  ctl.logger().addLogEntry("MyMeasures_" + name + "_range", [this, foot]() { return foot_data_[foot].range;} );

  ctl.gui()->addElement({"SoftFoot"},
    mc_rtc::gui::Point3D(name + "_point_ground", {mc_rtc::gui::Color::Green}, [this, foot](){ return foot_data_[foot].ground.back(); }),
    mc_rtc::gui::Trajectory(name + "_ground", {mc_rtc::gui::Color::Green}, [this, foot]() { return foot_data_[foot].ground; })
  );
}



EXPORT_SINGLE_STATE("BWC::SoftFoot", SoftFootState)
