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

#include <boost/numeric/odeint.hpp>

namespace
{

// https://stackoverflow.com/questions/27028226/python-linspace-in-c
template<typename T>
std::vector<T> linspace(T start_in, T end_in, int num_in)
{
  std::vector<T> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
  {
    linspaced.push_back(start);
    return linspaced;
  }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
  {
    linspaced.push_back(start + delta * i);
  }
  linspaced.push_back(end); // I want to ensure that start and end are exactly the same as the input
  return linspaced;
}

}

using namespace BWC;

void SoftFootState::start(mc_control::fsm::Controller & _ctl)
{
  State::start(_ctl);

  with_variable_stiffness_ = _ctl.config()("with_variable_stiffness", false);
  with_ankle_rotation_ = config_("with_ankle_rotation", false);
  with_foot_adjustment_ = config_("with_foot_adjustment", false);

  nr_phalanxes_ = config_("nr_phalanx", 10);
  phalanx_length_ = 0.27742 / static_cast<double>(nr_phalanxes_);

  // Display configuration
  mc_rtc::log::info("nr_phalanxes is set to {}", nr_phalanxes_);
  mc_rtc::log::info("phalanx_length is set to {}", phalanx_length_);
  mc_rtc::log::info("with_variable_stiffness is set to {}", with_variable_stiffness_);
  mc_rtc::log::info("with_ankle_rotation is set to {}", with_ankle_rotation_);
  mc_rtc::log::info("with_foot_adjustment is set to {}", with_foot_adjustment_);

  // Create client
  client_ = mc_rtc::ROSBridge::get_node_handle()->serviceClient<variable_stiffness::connectionFile>("connectionFile");

  if(nr_remaining_footstep_ == 0.)
  {
    nr_footstep_ = 0.;
  }
  
  // Initialize map
  foot_data_[Foot::Left] = FootData{};
  foot_data_[Foot::Right] = FootData{};

  _ctl.gui()->addElement({"SoftFoot"}, mc_rtc::gui::Label("cost", [this]() { return this->cost_; }));
  _ctl.logger().addLogEntry("cost", [this]() { return cost_; });

  _ctl.gui()->addElement({"SoftFoot"}, mc_rtc::gui::Label("nr_remaining_footstep", [this]() { return this->nr_remaining_footstep_; }));
  _ctl.logger().addLogEntry("nr_remaining_footstep", [this]() { return nr_remaining_footstep_; });

  _ctl.gui()->addElement({"SoftFoot"}, mc_rtc::gui::Label("nr_footstep", [this]() { return this->nr_footstep_; }));
  _ctl.logger().addLogEntry("nr_footstep", [this]() { return nr_footstep_; });

  _ctl.gui()->addElement({"SoftFoot"}, mc_rtc::gui::Label("PhalangesStiffness", [this]() { return this->PhalangesStiffness_; }));
  _ctl.logger().addLogEntry("PhalangesStiffness", [this]() { return PhalangesStiffness_; });

  _ctl.gui()->addXYPlot(
    "SoftFoot",
    // mc_rtc::gui::plot::Polygon("Ground",
    //   [&_ctl, this]()
    //   {
    //     // Construct Ground Polygon
    //     std::vector<std::array<double, 2>> points;
    //     // Get current moving foot
    //     Foot current_moving_foot = getCurrentMovingFoot(_ctl);
    //     // Get the segment
    //     const auto & raw_segment = ground_segment_[current_moving_foot].raw;
    //     // Save the selected segment in raw data of ground segment structure
    //     std::transform(raw_segment.begin(), raw_segment.end(), std::back_inserter(points),
    //       [](const Eigen::Vector3d & v) { return std::array<double, 2>{v.x(), v.z()}; }); 
    //     // Create polygon
    //     auto polygon = mc_rtc::gui::plot::PolygonDescription(points, mc_rtc::gui::Color::Magenta);
    //     polygon.closed(false);
    //     return polygon;
    //   }),
    mc_rtc::gui::plot::Polygon("Filtered_Ground",
      [&_ctl, this]()
      {
        // Construct Ground Polygon
        std::vector<std::array<double, 2>> points;
        // Get current moving foot
        Foot current_moving_foot = getCurrentMovingFoot(_ctl);
        // Get the segment
        const auto & filtered_segment = ground_segment_[current_moving_foot].filtered;
        // Save the selected segment in raw data of ground segment structure
        std::transform(filtered_segment.begin(), filtered_segment.end(), std::back_inserter(points),
          [](const Eigen::Vector3d & v) { return std::array<double, 2>{v.x(), v.z()}; }); 
        // Create polygon
        auto polygon = mc_rtc::gui::plot::PolygonDescription(points, mc_rtc::gui::Color::Blue);
        polygon.closed(false);
        return polygon;
      }),
    mc_rtc::gui::plot::Polygon("Convex",
      [&_ctl, this]()
      {
        // Construct Ground Polygon
        std::vector<std::array<double, 2>> points;
        // Get current moving foot
        Foot current_moving_foot = getCurrentMovingFoot(_ctl);
        // Get the segment
        const auto & convex_segment = ground_segment_[current_moving_foot].convex;
        // Save the selected segment in convex data of ground segment structure
        std::transform(convex_segment.begin(), convex_segment.end(), std::back_inserter(points),
          [](const Eigen::Vector3d & v) { return std::array<double, 2>{v.x(), v.z()}; }); 
        // Create polygon
        auto polygon = mc_rtc::gui::plot::PolygonDescription(points, mc_rtc::gui::Color::Green);
        polygon.closed(false);
        return polygon;
      }),
    mc_rtc::gui::plot::Polygon("Phalanxes",
      [&_ctl, this]()
      {
        // Get current moving foot
        Foot current_moving_foot = getCurrentMovingFoot(_ctl);
        // Create polygon
        auto polygon = mc_rtc::gui::plot::PolygonDescription(foot_data_[current_moving_foot].phalanxes, mc_rtc::gui::Color::Red);
        polygon.closed(false);
        return polygon;
      })
  );

  // Subscriber for the range sensors
  right_foot_range_sensor_sub_ = mc_rtc::ROSBridge::get_node_handle()->subscribe("range_sensor/distance", 1, &SoftFootState::rightFRSCallback, this);
  left_foot_range_sensor_sub_ = mc_rtc::ROSBridge::get_node_handle()->subscribe("range_sensor/distance", 1, &SoftFootState::leftFRSCallback, this);
}

bool SoftFootState::run(mc_control::fsm::Controller & ctl)
{
  ros::spinOnce();

  // Cast ctl to BaselineWalkingController
  auto & ctrl = static_cast<BWC::BaselineWalkingController&>(ctl);

  // Compute cost
  calculateCost(ctl);

  // Check if we are in single support or not
  if(ctrl.footManager_->supportPhase() == BWC::SupportPhase::DoubleSupport)
  {
    foot_data_[Foot::Right].need_reset = true;
    foot_data_[Foot::Left].need_reset = true;
    return false;
  }
  
  // Get the current moving foot
  Foot current_moving_foot = getCurrentMovingFoot(ctl);
  
  // Scope to handle reset/GUI
  {
    const Foot foot = current_moving_foot;
    const std::string name = foot == Foot::Left ? "left" : "right";
    if(foot_data_[current_moving_foot].need_reset)
    {
      ++ nr_footstep_;
      reset(ctl, current_moving_foot);
    }
    else if (!ctl.gui()->hasElement({"SoftFoot"}, name + "_point_ground"))
    {
      // Handle logger and gui
      ctl.logger().addLogEntry("MyMeasures_" + name + "_ground", [this, foot]() { return foot_data_[foot].ground.back();} );
      ctl.logger().addLogEntry("MyMeasures_" + name + "_filtered_x", [this, foot]()
        {
          std::vector<double> d;
          std::transform(ground_segment_[foot].filtered.begin(), ground_segment_[foot].filtered.end(), std::back_inserter(d),
            [](const Eigen::Vector3d & v) { return v.x(); }); 
          return d;
        });
      ctl.logger().addLogEntry("MyMeasures_" + name + "_filtered_y", [this, foot]()
        {
          std::vector<double> d;
          std::transform(ground_segment_[foot].filtered.begin(), ground_segment_[foot].filtered.end(), std::back_inserter(d),
            [](const Eigen::Vector3d & v) { return v.y(); }); 
          return d;
        });
      ctl.logger().addLogEntry("MyMeasures_" + name + "_convex_x", [this, foot]()
        {
          std::vector<double> d;
          std::transform(ground_segment_[foot].convex.begin(), ground_segment_[foot].convex.end(), std::back_inserter(d),
            [](const Eigen::Vector3d & v) { return v.x(); }); 
          return d;
        });
      ctl.logger().addLogEntry("MyMeasures_" + name + "_convex_y", [this, foot]()
        {
          std::vector<double> d;
          std::transform(ground_segment_[foot].convex.begin(), ground_segment_[foot].convex.end(), std::back_inserter(d),
            [](const Eigen::Vector3d & v) { return v.y(); }); 
          return d;
        });
      ctl.logger().addLogEntry("MyMeasures_" + name + "_range", [this, foot]() { return foot_data_[foot].range;} );
      ctl.logger().addLogEntry("MyMeasures_" + name + "_k", [this, foot]() { return foot_data_[foot].k;} );
      ctl.logger().addLogEntry("MyMeasures_" + name + "_angle", [this, foot]() { return foot_data_[foot].angle;} );
      ctl.logger().addLogEntry("MyMeasures_" + name + "_min_max_phalanxes_angle", [this, foot]() { return foot_data_[foot].min_max_phalanxes_angle;} );
      ctl.logger().addLogEntry("MyMeasures_" + name + "_position_offset", [this, foot]() { return foot_data_[foot].position_offset_x;} );
      ctl.logger().addLogEntry("MyMeasures_" + name + "_phalanxes_x", [this, foot]()
        {
          std::vector<double> d;
          std::transform(foot_data_[foot].phalanxes.begin(), foot_data_[foot].phalanxes.end(), std::back_inserter(d),
            [](const std::array<double, 2> & v) { return v[0]; }); 
          return d;
        });
      ctl.logger().addLogEntry("MyMeasures_" + name + "_phalanxes_y", [this, foot]()
        {
          std::vector<double> d;
          std::transform(foot_data_[foot].phalanxes.begin(), foot_data_[foot].phalanxes.end(), std::back_inserter(d),
            [](const std::array<double, 2> & v) { return v[1]; }); 
          return d;
        });

      ctl.gui()->addElement({"SoftFoot"},
        mc_rtc::gui::Point3D(name + "_point_ground", {mc_rtc::gui::Color::Green}, [this, foot](){ return foot_data_[foot].ground.back(); }),
        mc_rtc::gui::Trajectory(name + "_ground", {mc_rtc::gui::Color::Green}, [this, foot]() { return foot_data_[foot].ground; })
      );
    }
  }

  // Estimate ground from sensors
  estimateGround(ctl, current_moving_foot);

  // Check foot position with respect to desired landing pose
  const auto & ground = foot_data_[current_moving_foot].ground;
  // TODO: Get landing: a bit dirty for the moment
  sva::PTransformd X_0_landing = ctrl.footManager_->swingTraj_->endPose_;
  // If we saw more than the landing pose + half of the foot, we have enough data to perform all the computations
  if(!foot_data_[current_moving_foot].computation_done && ground.back().x() >= X_0_landing.translation().x() + foot_length_ * 0.5 + extra_to_compute_best_position_)
  {
    mc_rtc::log::success("Accumulated enough data");
    // We need to do these steps only one time
    foot_data_[current_moving_foot].computation_done = true;
    // Extract ground segment from ground data
    extractGroundSegment(ctl, current_moving_foot, X_0_landing.translation());
    // Continue only if we have enough data such as the foot is smaller than the estimated ground
    if(ground_segment_[current_moving_foot].raw.back().x() - ground_segment_[current_moving_foot].raw.front().x() > 0.95 * foot_length_)
    {
      // Compute the altitude profile
      extractAltitudeProfileFromGroundSegment(current_moving_foot);
      // call the server and update the variable stiffness
      if(with_variable_stiffness_)
      {
        updateVariableStiffness(ctl, current_moving_foot);
      }
      // Compute convex hull of the segment -> right now it does not work 
      computeSegmentConvexHull(ctl, current_moving_foot);
      // Find best landing position
      if(with_foot_adjustment_)
      {
        // Compute min/max phalanxes angle depending of stiffness 
        computeMinMaxAngle(current_moving_foot);
        // Compute the landing position
        computeFootLandingPosition(current_moving_foot, X_0_landing.translation());
        // Update landing pose
        X_0_landing = sva::PTransformd(Eigen::Vector3d(foot_data_[current_moving_foot].position_offset_x, 0., 0.)) * X_0_landing;
      }
      if(with_ankle_rotation_)
      {
        // Now with the convex hull we can compute the angle
        computeFootLandingAngle(current_moving_foot, X_0_landing.translation());
      }
      // Update targeted pose
      updateFootSwingPose(ctl, current_moving_foot, X_0_landing);
    }
    else
    {
      mc_rtc::log::warning("Too small {} > {}", ground_segment_[current_moving_foot].raw.back().x() - ground_segment_[current_moving_foot].raw.front().x(), foot_length_);
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
    nr_remaining_footstep_ =  ctrl.footManager_->footstepQueue().size() ;
    footstep = lambda_footstep_ * nr_remaining_footstep_;
  } 
  
  cost_ = - footstep - zmp - CoM;
  //mc_rtc::log::success("Foot_Target_Pose : {}, Foot_Robot_Pose : {}, error_distance : {}", Foot_Target_Pose.translation().head<2>(), Foot_Robot_Pose.translation().head<2>(), error_distance);
  // mc_rtc::log::info("footstep {}", footstep);
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
  {
    const std::lock_guard<std::mutex> lock(range_sensor_mutex_);
    data.range = ctl.robot().device<mc_mujoco::RangeSensor>(sensor_name).data();
  }
  // mc_rtc::log::error("data.range {}", data.range);
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
    [&](const Eigen::Vector3d & v) { return v.x() >= landing.x() - foot_length_ * 0.5 - extra_to_compute_best_position_; });
  const auto end_iterator = std::find_if(ground.begin(), ground.end(),
    [&](const Eigen::Vector3d & v) { return v.x() >= landing.x() + foot_length_ * 0.5 + extra_to_compute_best_position_; });
  
  // Save the selected segment in raw data of ground segment structure
  // ground_segment_[current_moving_foot].raw.clear();
  std::vector<Eigen::Vector3d> raw;
  std::transform(begin_iterator, end_iterator, std::back_inserter(raw), [](const Eigen::Vector3d & v) { return v; }); 
  mc_rtc::log::error("raw.size() {}", raw.size());
  // Check if there are missing data in collected data due to obstruction
  // Interpolate points in-between point from convex hull
  auto lerp = [](const Eigen::Vector3d & A, const Eigen::Vector3d& B, double t)
  {
    const Eigen::Vector3d vec = B * t + A * (1. - t);
    return vec;
  };
  // If yes, interpolate data
  for(size_t i = 0; i < raw.size() - 1; ++i)
  {
    ground_segment_[current_moving_foot].raw.push_back(raw[i]);
    // Check if next element is far
    if(raw[i+1].x() - raw[i].x() > 0.01)
    {
      // We need to interpolate data
      for(double t = 0.; t <= 1.0; t += 0.05)
      {
        ground_segment_[current_moving_foot].raw.push_back(lerp(raw[i], raw[i+1], t));
      }
    }
  }
  
  // Filter the trajectory
  {
    const auto& raw = ground_segment_[current_moving_foot].raw;
    auto& filtered = ground_segment_[current_moving_foot].filtered;
    const double ePow = 1. - std::exp(-0.005 * 2 * M_PI * 5.);
    filtered.push_back(raw.front());
    for(size_t i = 1; i < raw.size(); ++i)
    {
      filtered.emplace_back(raw[i].x(), raw[i].y(), filtered.back().z() + (raw[i].z() - filtered.back().z()) * ePow);
    }

    // Add small variation noise
    const double noise = 0.00075;
    double direction = 1.0;
    for(auto& f: filtered)
    {
      f.z() += direction * noise;
      direction *= -1.;
    }
  }

  mc_rtc::log::info("ground_segment_[{}].raw {}", current_moving_foot == Foot::Left ? "Left" : "Right", ground_segment_[current_moving_foot].raw.size());
  mc_rtc::log::info("ground_segment_[{}].filtered {}", current_moving_foot == Foot::Left ? "Left" : "Right", ground_segment_[current_moving_foot].filtered.size());
  
  // Add trajectory
  if(ground_segment_[current_moving_foot].raw.size() > 5)
  {
    const std::string name = current_moving_foot == Foot::Left ? "left" : "right";  
    ctl.gui()->addElement({"SoftFoot"},
      mc_rtc::gui::Trajectory(name + "_segment", {mc_rtc::gui::Color::Red}, [this, current_moving_foot]() { return ground_segment_[current_moving_foot].raw; })
    );
    ctl.gui()->addElement({"SoftFoot"},
      mc_rtc::gui::Trajectory(name + "_segment_filtered", {mc_rtc::gui::Color::Magenta}, [this, current_moving_foot]() { return ground_segment_[current_moving_foot].filtered; })
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
    PhalangesStiffness_ = foot_data_[current_moving_foot].k;
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
  const auto & raw_segment = ground_segment_[current_moving_foot].filtered;

  const std::vector<Eigen::Vector3d> convex_hull = computeConvexHull(raw_segment);

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
 
  ground_segment_[current_moving_foot].convex = convex_hull;

  // Add trajectory
  ctl.gui()->addElement({"SoftFoot"},
    mc_rtc::gui::Trajectory(name + "_convex_segment", {mc_rtc::gui::Color::Magenta}, [this, current_moving_foot]() { return ground_segment_[current_moving_foot].convex; })
  );
}

std::vector<Eigen::Vector3d> SoftFootState::computeConvexHull(const std::vector<Eigen::Vector3d>& data)
{
  std::vector<Eigen::Vector3d> ret;
  
  // Build the input for qhull
  std::vector<double> points_in;
  double min = std::numeric_limits<double>::max();
  double max = std::numeric_limits<double>::lowest();
  // Skip some points
  for(size_t i = 0; i < data.size(); i += 2)
  {
    const auto& p = data[i];
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
      for(size_t i = 0; i < qhull.facetCount(); ++i)
      {
        face = face.nextFacet2d(&v);
        convex_hull.emplace_back(v.point().coordinates()[0], v.point().coordinates()[1]);
      }

      // Select only the upper part of the convex hull
      size_t idx_min_x = std::numeric_limits<size_t>::max();
      double min_x = std::numeric_limits<double>::max();
      size_t idx_max_x = std::numeric_limits<size_t>::lowest();
      double max_x = std::numeric_limits<double>::lowest();
      for(size_t i = 0; i < convex_hull.size(); ++i)
      {
        if(convex_hull[i].x() < min_x)
        {
          min_x = convex_hull[i].x();
          idx_min_x = i;
        } 

        if(convex_hull[i].x() > max_x)
        {
          max_x = convex_hull[i].x();
          idx_max_x = i;
        } 
      }
      //
      std::vector<Eigen::Vector2d> selected_part_0;
      double max_part_0 = std::numeric_limits<double>::lowest();
      for(size_t i = idx_min_x; i <= idx_max_x; ++i)
      {
        selected_part_0.push_back(convex_hull[i]);
        max_part_0 = selected_part_0.back().y() > max_part_0 ? selected_part_0.back().y() : max_part_0;
      }
      //
      std::vector<Eigen::Vector2d> selected_part_1;
      double max_part_1 = std::numeric_limits<double>::lowest();
      int idx = idx_max_x;
      while (idx != idx_min_x)
      {
        selected_part_1.push_back(convex_hull[idx]);
        idx = (idx + 1) % convex_hull.size();
        max_part_1 = selected_part_1.back().y() > max_part_0 ? selected_part_1.back().y() : max_part_0;
      }
      selected_part_1.push_back(convex_hull[idx_min_x]);
      //
      if(max_part_0 > max_part_1)
      {
        convex_hull = selected_part_0;
      }
      else
      {
        convex_hull = selected_part_1;
      }
      // Sort alongside x
      std::sort(convex_hull.begin(), convex_hull.end(), [](const Eigen::Vector2d& a, const Eigen::Vector2d& b) { return a.x() < b.x()  ; });

    }
    catch (std::exception& e)
    {
      mc_rtc::log::error("Error during qhull ! {}", e.what());
      convex_hull.push_back(Eigen::Vector2d(data.front().x(), data.front().z()));
      convex_hull.push_back(Eigen::Vector2d(data.back().x(), data.back().z()));
    }
  }
  else
  {
    mc_rtc::log::error("Too flat to run qhull !");
    convex_hull.push_back(Eigen::Vector2d(data.front().x(), data.front().z()));
    convex_hull.push_back(Eigen::Vector2d(data.back().x(), data.back().z()));
  }

  // Interpolate points in-between point from convex hull
  auto lerp = [](const Eigen::Vector2d & A, const Eigen::Vector2d& B, double y, double t)
  {
    const Eigen::Vector2d vec = B * t + A * (1. - t);
    return Eigen::Vector3d(vec.x(), y, vec.y());
  };

  for(size_t i = 0; i < convex_hull.size() - 1; ++i)
  {
    for(double t = 0.; t <= 1.0; t += 0.1)
    {
      ret.push_back(lerp(convex_hull[i], convex_hull[i+1], data.front().y(), t));
    }
  }

  return ret;
}

void SoftFootState::computeMinMaxAngle(const Foot & current_moving_foot)
{
  if(nr_phalanxes_ > 1.)
  {
    const double MASS_PHALANX = 0.2;
    const double MASS_ROBOT = 43.;
    const double nr_phalanxes = static_cast<double>(nr_phalanxes_);
    const double k = foot_data_[current_moving_foot].k;
    const double d = 2.;
    const double f0 = MASS_ROBOT*9.81*(foot_length_/2.); 
    const double I = (10.*MASS_PHALANX*phalanx_length_*phalanx_length_)/3. ;
    
    auto odefcn = [=](const std::vector<double>& x, std::vector<double>& dx, const double t)
    {
        dx[0] = x[1];
        dx[1] = (- (nr_phalanxes - 1.) * k * x[0] - (nr_phalanxes - 1.) * d * x[1] + f0) / I;
    };

    double angle = 0.;
    auto observer = [&angle](const std::vector<double>& x, const double t)
    {
        angle = x[0] * 180. / M_PI;
    };

    typedef boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>> stepper_type;
    const double dt = 0.1;
    std::vector<double> x(2);
    x[0] = 0.0;
    x[1] = 0.0;
    boost::numeric::odeint::integrate_const(boost::numeric::odeint::make_dense_output<stepper_type>( 1E-6, 1E-3 ),
          odefcn, x, 0.0, 10.0, dt , observer);

    angle = angle / (nr_phalanxes - 1.);

    mc_rtc::log::info("Min/max angle is {}", angle);

    foot_data_[current_moving_foot].min_max_phalanxes_angle = angle;
  }
  else
  {
    foot_data_[current_moving_foot].min_max_phalanxes_angle = M_PI;
  }
}

void SoftFootState::computeFootLandingPosition(const Foot & current_moving_foot, const Eigen::Vector3d & default_landing_pos)
{
  //
  auto computeThetaBasedOnQhull = [&](const Eigen::Vector3d & start, const Eigen::Vector3d & end)
  {
    const auto & ground = ground_segment_[current_moving_foot].filtered;

    std::vector<Eigen::Vector3d> ground_under_phalanx;

    const auto begin_iterator = std::find_if(ground.begin(), ground.end(),
      [&](const Eigen::Vector3d & v) { return v.x() >= start.x() - 0.01; });
    const auto end_iterator = std::find_if(ground.begin(), ground.end(),
      [&](const Eigen::Vector3d & v) { return v.x() >= end.x() + 0.01; });
    std::transform(begin_iterator, end_iterator, std::back_inserter(ground_under_phalanx), [](const Eigen::Vector3d & v) { return v; });
    
    const std::vector<Eigen::Vector3d> convex = computeConvexHull(ground_under_phalanx);

    const auto start_convex_iterator = std::find_if(convex.begin(), convex.end(), [&](const Eigen::Vector3d & v) { return v.x() >= start.x(); });
    const auto end_convex_iterator = std::find_if(convex.begin(), convex.end(), [&](const Eigen::Vector3d & v) { return v.x() >= end.x(); });

    const Eigen::Vector3d & p_0 = *start_convex_iterator;
    const Eigen::Vector3d & p_1 = *end_convex_iterator;

    const double dz = p_1.z() - p_0.z();
    const double dx = p_1.x() - p_0.x();

    const double theta = std::atan(dz / dx);

    return theta;
  };

  //
  auto whichPhalanxItIs = [&](const Eigen::Vector3d & landing_pos, const Eigen::Vector3d & highest_point_on_convex) -> size_t
  {
    const auto & convex = ground_segment_[current_moving_foot].convex;

    // Compute phalanxes pos
    const double start_x = (sva::PTransformd(Eigen::Vector3d(-foot_length_ * 0.5, 0., 0.)) * landing_pos).translation().x();
    const double end_x = (sva::PTransformd(Eigen::Vector3d(+foot_length_ * 0.5, 0., 0.)) * landing_pos).translation().x();

    double ratio = (highest_point_on_convex.x() - start_x) / phalanx_length_;
    return static_cast<size_t>(ratio);
  };

  //
  auto addLeftPhalanxToPhalanxes = [&](const Eigen::Vector3d& phalanx_pos, double theta)
  {
    const std::array<double, 2> start = {phalanx_pos.x(), phalanx_pos.z()};
    foot_data_[current_moving_foot].phalanxes.push_back(start);

    const std::array<double, 2> end = {
      phalanx_pos.x() + phalanx_length_ * std::cos(theta),
      phalanx_pos.z() + phalanx_length_ * std::sin(theta)
    };
    foot_data_[current_moving_foot].phalanxes.push_back(end);
  };
  //
  auto addRightPhalanxToPhalanxes = [&](const Eigen::Vector3d& phalanx_pos, double theta)
  {
    const std::array<double, 2> start = {
      phalanx_pos.x() - phalanx_length_ * std::cos(theta),
      phalanx_pos.z() - phalanx_length_ * std::sin(theta)
    };
    foot_data_[current_moving_foot].phalanxes.push_back(start);

    const std::array<double, 2> end = {phalanx_pos.x(), phalanx_pos.z()};
    foot_data_[current_moving_foot].phalanxes.push_back(end);
  };
  
  //
  auto computeSide = [&](int start, int end, double direction, const Eigen::Vector3d & phalanx_pos,
    double theta, std::vector<double>& ret_phalanxes_theta, std::vector<Eigen::Vector2d>& ret_contact_points) -> std::vector<Eigen::Vector3d>
  {
    std::vector<Eigen::Vector3d> ret;
    if(start == end)
    {
      return ret;
    }
    const std::vector<Eigen::Vector3d> & ground = ground_segment_[current_moving_foot].filtered;

    // std::cout << "phalanx_pos " << phalanx_pos.transpose() << std::endl;

    Eigen::Vector3d next_phalanx_side_pos(
      phalanx_pos.x() + direction * (0.5 * phalanx_length_) * std::cos(theta),
      phalanx_pos.y(),
      phalanx_pos.z() + direction * (0.5 * phalanx_length_) * std::sin(theta)
    );

    // std::cout << "phalanx_pos_side " << next_phalanx_side_pos.transpose() << std::endl;

    // Create linspace of theta to evaluate using the foot_data_[current_moving_foot].min_max_phalanxes_angle
    const double min_max_phalanxes_angle = std::min(M_PI * 0.5, foot_data_[current_moving_foot].min_max_phalanxes_angle);
    const double start_theta = direction * min_max_phalanxes_angle - direction * 0.01;
    const double end_theta = -direction * min_max_phalanxes_angle + direction * 0.01;
    // const double start_theta = direction * M_PI * 0.5 - direction * 0.01;
    // const double end_theta = -direction * M_PI * 0.5  + direction * 0.01;
    // std::cout << "direction " << direction << std::endl; 
    std::cout << "start_theta " << start_theta << " end_theta " << end_theta << " - " << static_cast<int>(min_max_phalanxes_angle * 180. / M_PI) << std::endl;
    const std::vector<double> evaluated_thetas = linspace(start_theta, end_theta, static_cast<int>(min_max_phalanxes_angle * 180. / M_PI)); // [rad]
    //
    const std::vector<int> phalanx_idxs = linspace(start, end, std::abs(end - start) + 1);
    for(const auto& phalanx_i: phalanx_idxs)
    {
      double theta_intersection = 0.;
      bool is_intersection = false;
      Eigen::Vector2d contact_found;

      if(phalanx_i < nr_phalanxes_)
      {        
        Eigen::Vector3d phalanx_pos_left;
        Eigen::Vector3d phalanx_pos_right;

        std::vector<Eigen::Vector3d> ground_under_phalanx;
        // Find beginning and ending of segment to extract to select ground under phalanx
        if(direction == -1)
        {
          phalanx_pos_right = next_phalanx_side_pos;
          phalanx_pos_left = Eigen::Vector3d(
            phalanx_pos_right.x() - phalanx_length_,
            phalanx_pos_right.y(),
            phalanx_pos_right.z()
          );
        }
        else
        {
          phalanx_pos_left = next_phalanx_side_pos;
          phalanx_pos_right = Eigen::Vector3d(
            phalanx_pos_left.x() + phalanx_length_,
            phalanx_pos_left.y(),
            phalanx_pos_left.z()
          );
        }
        const auto begin_iterator = std::find_if(ground.begin(), ground.end(),
          [&](const Eigen::Vector3d & v) { return v.x() >= phalanx_pos_left.x(); });
        const auto end_iterator = std::find_if(ground.begin(), ground.end(),
          [&](const Eigen::Vector3d & v) { return v.x() >= phalanx_pos_right.x(); });
        std::transform(begin_iterator, end_iterator, std::back_inserter(ground_under_phalanx), [](const Eigen::Vector3d & v) { return v; }); 
        // std::cout << ground_under_phalanx.size() << std::endl;
        if(!ground_under_phalanx.empty())
        {
          for(const auto& ev_theta: evaluated_thetas)
          {
            // std::cout << "phalanx_pos_left " << phalanx_pos_left.transpose() << " phalanx_pos_right " << phalanx_pos_right.transpose() << std::endl;
            is_intersection = evaluateTheta(ev_theta, direction, next_phalanx_side_pos, ground_under_phalanx, contact_found);
            if(is_intersection)
            {
              const double step_theta = 2. * M_PI / 180.;
              const double start_theta = ev_theta + direction * step_theta;
              const double end_theta = ev_theta - direction * step_theta;
              const std::vector<double> evaluated_precise_thetas = linspace(start_theta, end_theta, 10); // [rad]
              for(const auto& ev_pre_theta: evaluated_precise_thetas)
              {
                is_intersection = evaluateTheta(ev_pre_theta, direction, next_phalanx_side_pos, ground_under_phalanx, contact_found);
                if(is_intersection)
                {
                  // std::cout << "ev_pre_theta: " << ev_pre_theta << std::endl;
                  theta_intersection = ev_pre_theta;
                  break;
                }
              }
              break;             
            }
          }
        }

        if(is_intersection)
        {
          // std::cout << "intersection found" << std::endl;
          ret_contact_points.push_back(contact_found);
        }
        else
        {
          std::cout << "no intersection found" << std::endl;
          theta_intersection = evaluated_thetas.back();
        }

        // Update pos
        next_phalanx_side_pos = Eigen::Vector3d(
          next_phalanx_side_pos.x() + direction * phalanx_length_ * std::cos(theta_intersection),
          next_phalanx_side_pos.y(),
          next_phalanx_side_pos.z() + direction * phalanx_length_ * std::sin(theta_intersection)
        );

        // std::cout << next_phalanx_side_pos.transpose() << std::endl;

        ret_phalanxes_theta.push_back(theta_intersection);
        ret.push_back(next_phalanx_side_pos);
        // break;
      }
    }
    
    return ret;
  };

  const std::vector<double> position_offsets_x = {0., -0.05, -0.03, -0.01, 0.01, 0.03, 0.05};
  // const std::vector<double> position_offsets_x = {0.};
  const auto & convex = ground_segment_[current_moving_foot].convex;
  double maximized_distance_between_phalanxes = std::numeric_limits<double>::lowest();
  for(const auto & position_offset_x: position_offsets_x)
  {
    // Compute landing pos to evaluate
    const sva::PTransformd landing_pos = sva::PTransformd(Eigen::Vector3d(position_offset_x, 0., 0.)) * default_landing_pos;
    // Get min and max foot x
    const double min_x_foot = (sva::PTransformd(Eigen::Vector3d(-foot_length_ * 0.5, 0., 0.)) * landing_pos).translation().x();
    const double max_x_foot = (sva::PTransformd(Eigen::Vector3d(+foot_length_ * 0.5, 0., 0.)) * landing_pos).translation().x();   
    // Find highest point on hull
    Eigen::Vector3d highest_point_on_convex(0., 0., std::numeric_limits<double>::lowest());
    for(const auto & pos: convex)
    {
      if(pos.x() >= min_x_foot && pos.x() <= max_x_foot && pos.z() > highest_point_on_convex.z())
      {
        highest_point_on_convex = pos;
      }
    }

    // std::cout << "Highest point on convex hull is " << highest_point_on_convex.transpose() << std::endl;
    // Compute which n-th phalanx will land on the highest point
    size_t n_phalanx = whichPhalanxItIs(landing_pos.translation(), highest_point_on_convex);
    // Keep highest altitude
    // foot_data_[current_moving_foot].position_offset_z = std::max(foot_data_[current_moving_foot].position_offset_z, highest_point_on_convex.z() + 0.015 - 0.04); // 0.04 is altitude of foot from configuration file
    // std::cout << "The n-th phalanx is " <<  n_phalanx << std::endl;
    // phalanx_pos is the middle of the phalanx
    const Eigen::Vector3d phalanx_pos_start = (sva::PTransformd(Eigen::Vector3d(-foot_length_ * 0.5 + n_phalanx * phalanx_length_, 0., 0.)) * landing_pos).translation();
    const Eigen::Vector3d phalanx_pos_end = (sva::PTransformd(Eigen::Vector3d(-foot_length_ * 0.5 + (n_phalanx + 1) * phalanx_length_, 0., 0.)) * landing_pos).translation();
    
    // Compute landing theta at this point for the n-th phalanx
    const double theta = computeThetaBasedOnQhull(phalanx_pos_start, phalanx_pos_end);
    // std::cout << "theta " << theta << std::endl;
    // std::cout << highest_point_on_convex.transpose() << std::endl;
  
    // std::cout << "min_x_foot " << min_x_foot << " max_x_foot " << max_x_foot << std::endl;
    Eigen::Vector3d phalanx_pos = (sva::PTransformd(Eigen::Vector3d(-foot_length_ * 0.5 + (n_phalanx + 0.5) * phalanx_length_, 0., 0.)) * landing_pos).translation();
    phalanx_pos.z() = highest_point_on_convex.z();
    // std::cout << "phalanx_pos " << phalanx_pos.transpose() << std::endl;
    // std::cout << "start phalanx " << min_x_foot + phalanx_length_ * n_phalanx << std::endl;

    const Eigen::Vector3d phalanx_pos_left(
      phalanx_pos.x() - (0.5 * phalanx_length_) * std::cos(theta),
      phalanx_pos.y(),
      phalanx_pos.z() - (0.5 * phalanx_length_) * std::sin(theta)
    );
    const Eigen::Vector3d phalanx_pos_right(
      phalanx_pos.x() + (0.5 * phalanx_length_) * std::cos(theta),
      phalanx_pos.y(),
      phalanx_pos.z() + (0.5 * phalanx_length_) * std::sin(theta)
    );

    // std::cout << "phalanx_pos_left " << phalanx_pos_left.transpose() << std::endl;
    // std::cout << "phalanx_pos_right " << phalanx_pos_right.transpose() << std::endl;

    // std::cout << "n_phalanx " << n_phalanx << " nr_phalanxes_ " << nr_phalanxes_ << std::endl;
    std::vector<double> right_side_phalanxes_theta;
    std::vector<Eigen::Vector2d> right_contact_points;
    std::vector<Eigen::Vector3d> right_side_phalanxes_pos;

    std::vector<double> left_side_phalanxes_theta;
    std::vector<Eigen::Vector2d> left_contact_points;
    std::vector<Eigen::Vector3d> left_side_phalanxes_pos;

    double distance_between_phalanxes = 0.;

    if(nr_phalanxes_ > 2)
    {
      // Compute right part
      // std::cout << "Right side" << std::endl;
      right_side_phalanxes_pos = computeSide(n_phalanx, nr_phalanxes_, +1., phalanx_pos, theta, right_side_phalanxes_theta, right_contact_points);

      // Compute left part
      // std::cout << "Left side" << std::endl;
      left_side_phalanxes_pos = computeSide(n_phalanx, 0, -1., phalanx_pos, theta, left_side_phalanxes_theta, left_contact_points);

      if(left_contact_points.size() > 0 && right_contact_points.size() > 0)
      {
        distance_between_phalanxes = right_contact_points.back().x() - left_contact_points.back().x();
      }
      else if(left_contact_points.size() > 0)
      {
        distance_between_phalanxes = phalanx_pos.x() - left_contact_points.back().x();
      }
      else if(right_contact_points.size() > 0)
      {
        distance_between_phalanxes = right_contact_points.back().x() - phalanx_pos.x();
      }
    }
    else
    {
      // To check if the intersection belongs to the segment we are evaluating
      auto isPointOnSegment = [](const Eigen::Vector2d& A, const Eigen::Vector2d& B, const Eigen::Vector2d& C) -> bool
      {
        const Eigen::Vector2d AB = B - A;
        const Eigen::Vector2d AC = C - A;

        // Check the they are collinear
        if(std::abs(AB.x() * AC.y() - AC.x() * AB.y()) < 1e-3)
        {
          // Check if C is in between A and B
          const double dot_AB = AB.dot(AB);
          const double dot_AC = AB.dot(AC);

          // // std::cout << "dot_AB " << dot_AB << " dot_AC " << dot_AC << std::endl;
          if (dot_AC == 0. || dot_AC == dot_AB || (dot_AC > 0. && dot_AC < dot_AB))
          {
            // // std::cout << "Point on segment" << A.transpose() << " " << B.transpose() << " " << C.transpose() << std::endl;
            return true;
          }
        }
        return false;
      };

      const std::vector<Eigen::Vector3d> & ground = ground_segment_[current_moving_foot].filtered;
      std::vector<Eigen::Vector3d> ground_under_phalanx;
      const auto begin_iterator = std::find_if(ground.begin(), ground.end(),
        [&](const Eigen::Vector3d & v) { return v.x() >= phalanx_pos_left.x(); });
      const auto end_iterator = std::find_if(ground.begin(), ground.end(),
        [&](const Eigen::Vector3d & v) { return v.x() >= phalanx_pos_right.x(); });
      std::transform(begin_iterator, end_iterator, std::back_inserter(ground_under_phalanx), [](const Eigen::Vector3d & v) { return v; }); 
      
      // const std::array<double, 2> start = {phalanx_pos_left.x(), phalanx_pos_left.z()};
      // const std::array<double, 2> end = {phalanx_pos_right.x(), phalanx_pos_right.z()};
      
      const Eigen::Vector2d start(phalanx_pos_left.x(), phalanx_pos_left.z());
      const Eigen::Vector2d end(phalanx_pos_right.x(), phalanx_pos_right.z());
      
      // Compute line equation
      const double slope = (end.y() - start.y()) / (end.x() - start.x());
      const double intercept = end.y() - slope * end.x();

      auto f = [slope, intercept](double x) -> double
      {
        return slope * x + intercept;
      };

      std::vector<Eigen::Vector2d> contact_points;
      // Compute points for line that matches
      if(!ground_under_phalanx.empty())
      {
        double previous_diff_y = f(ground_under_phalanx[0].x()) - ground_under_phalanx[0].z();
        for(size_t i = 1; i < ground_under_phalanx.size(); ++i)
        {
          const double line_y = f(ground_under_phalanx[i].x());
          const double diff_y = line_y - ground_under_phalanx[i].z();
          // We check if it is a possible 0
          if(previous_diff_y * diff_y <= 0.)
          {
            // Check if the point is on the segment
            Eigen::Vector2d output(ground_under_phalanx[i].x(), ground_under_phalanx[i].z());
            if(isPointOnSegment(start, end, output))
            {
              contact_points.push_back(output);
            }
          }
          previous_diff_y = diff_y;
        }
      }
      std::cout << "nr contact points found " << contact_points.size() << std::endl;
      // distance_between_phalanxes = TODO do something here to compute something
    }

    if(foot_data_[current_moving_foot].phalanxes.empty() || distance_between_phalanxes > maximized_distance_between_phalanxes)
    {
      maximized_distance_between_phalanxes = distance_between_phalanxes;
      foot_data_[current_moving_foot].position_offset_x = position_offset_x;
      mc_rtc::log::info("New maximized found for {}", position_offset_x);
      mc_rtc::log::info("distance_between_phalanxes {}", distance_between_phalanxes);
      foot_data_[current_moving_foot].phalanxes.clear();
      // Fill phalanxes
      for(int i = left_side_phalanxes_pos.size() - 1; i >= 0 ; --i)
      {
        addLeftPhalanxToPhalanxes(left_side_phalanxes_pos[i], left_side_phalanxes_theta[i]);
      }

      const std::array<double, 2> start = {phalanx_pos_left.x(), phalanx_pos_left.z()};
      foot_data_[current_moving_foot].phalanxes.push_back(start);

      const std::array<double, 2> end = {phalanx_pos_right.x(), phalanx_pos_right.z()};
      foot_data_[current_moving_foot].phalanxes.push_back(end);

      for(int i = 0; i < right_side_phalanxes_pos.size() ; ++i)
      {
        addRightPhalanxToPhalanxes(right_side_phalanxes_pos[i], right_side_phalanxes_theta[i]);
      }
    }
  }
}

bool SoftFootState::evaluateTheta(double theta, double direction, const Eigen::Vector3d& pos, const std::vector<Eigen::Vector3d>& ground_under_phalanx, Eigen::Vector2d& output)
{
  // To check if the intersection belongs to the segment we are evaluating
  auto isPointOnSegment = [](const Eigen::Vector2d& A, const Eigen::Vector2d& B, const Eigen::Vector2d& C) -> bool
  {
    const Eigen::Vector2d AB = B - A;
    const Eigen::Vector2d AC = C - A;

    // Check the they are collinear
    if(std::abs(AB.x() * AC.y() - AC.x() * AB.y()) < 1e-3)
    {
      // Check if C is in between A and B
      const double dot_AB = AB.dot(AB);
      const double dot_AC = AB.dot(AC);

      // // std::cout << "dot_AB " << dot_AB << " dot_AC " << dot_AC << std::endl;
      if (dot_AC == 0. || dot_AC == dot_AB || (dot_AC > 0. && dot_AC < dot_AB))
      {
        // // std::cout << "Point on segment" << A.transpose() << " " << B.transpose() << " " << C.transpose() << std::endl;
        return true;
      }
    }
    return false;
  };

  //
  const Eigen::Vector2d start_offset(
    pos.x() + direction * 0.005 * std::cos(theta),
    pos.z() + direction * 0.005 * std::sin(theta)
  );

  
  const Eigen::Vector2d end(
    pos.x() + direction * phalanx_length_ * std::cos(theta),
    pos.z() + direction * phalanx_length_ * std::sin(theta)
  );

  // Compute line equation
  const double slope = (end.y() - pos.z()) / (end.x() - pos.x());
  const double intercept = end.y() - slope * end.x();

  auto f = [slope, intercept](double x) -> double
  {
    return slope * x + intercept;
  };

  // std::cout << "-----------------------" << std::endl;
  // std::cout << "For theta: " << theta << std::endl;
  // // std::cout << pos.z() << " = " << f(pos.x()) << " -- " << end.y() << " = " << f(end.x()) << std::endl;
  // std::cout << "start_offset " << start_offset.transpose() << " end " << end.transpose() << std::endl;
  // std::cout << "gxs " << ground_under_phalanx.front().x() << " gxe " << ground_under_phalanx.back().x() << std::endl;
  // Compute points for line that matches
  double previous_diff_y = f(ground_under_phalanx[0].x()) - ground_under_phalanx[0].z();
  for(size_t i = 1; i < ground_under_phalanx.size(); ++i)
  {
    const double line_y = f(ground_under_phalanx[i].x());
    const double diff_y = line_y - ground_under_phalanx[i].z();
    // std::cout.precision(4);
    // std::cout << ground_under_phalanx[i].x() << " (" << line_y << " - " << ground_under_phalanx[i].z() << " = " << diff_y << ") ";
    // We check if it is a possible 0
    if(previous_diff_y * diff_y <= 0.)
    {
      // Check if the point is on the segment
      // // std::cout << "Found a zero at " << i << std::endl;
      // // std::cout << "Ground " << Eigen::Vector2d(ground_under_phalanx[i].x(), ground_under_phalanx[i].z()).transpose() << std::endl;
      output = Eigen::Vector2d(ground_under_phalanx[i].x(), ground_under_phalanx[i].z());
      if(isPointOnSegment(start_offset, end, output))
      {
        return true;
      }
    }
    previous_diff_y = diff_y;
  }
  // std::cout << std::endl;

  return false;
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
  double desired_offset_position_x = foot_data_[current_moving_foot].position_offset_x;
  double desired_offset_position_z = foot_data_[current_moving_foot].position_offset_z;
  // Dirty access to swing traj
  
  // Follow Murooka-san code to update, sorry for the dirty access
  // The order is important
  if(with_foot_adjustment_)
  {
    dynamic_cast<BWC::SwingTrajLandingSearch*>(ctrl.footManager_->swingTraj_.get())->updatePosXZ(desired_offset_position_x, desired_offset_position_z);
  }
  if(with_ankle_rotation_)
  {
    dynamic_cast<BWC::SwingTrajLandingSearch*>(ctrl.footManager_->swingTraj_.get())->updatePitch(desired_angle);
  }
}


void SoftFootState::reset(mc_control::fsm::Controller & ctl, const Foot & foot)
{
  foot_data_[foot].need_reset = false;
  foot_data_[foot].computation_done = false;
  foot_data_[foot].min_max_phalanxes_angle = 0.;
  foot_data_[foot].position_offset_x = 0.;
  foot_data_[foot].position_offset_z = 0.;
  foot_data_[foot].k = 0.;
  // Delete all the data that are behind the foot
  const auto X_0_p = ctl.robot().surfacePose(surface_name_[foot]);
  auto & ground = foot_data_[foot].ground;
  const auto ground_iterator = std::find_if(ground.begin(), ground.end(), [&](const Eigen::Vector3d & v) { return v.x() >= X_0_p.translation().x() - 0.5 * foot_length_ - extra_to_compute_best_position_; });
  // Delete from the beginning of the vector up to the back part of the foot
  ground.erase(ground.begin(), ground_iterator);
  
  // Reset altitude
  foot_data_[foot].altitude.clear();
  foot_data_[foot].phalanxes.clear();

  // Reset ground segment data
  ground_segment_[foot].raw.clear();
  ground_segment_[foot].filtered.clear();
  ground_segment_[foot].convex.clear();

  // Get name for the feet
  const Foot other_foot = foot == Foot::Left ? Foot::Right : Foot::Left;  
  const std::string other_name = other_foot == Foot::Left ? "left" : "right";
  const std::string name = foot == Foot::Left ? "left" : "right";  
  // Reset logger
  ctl.logger().removeLogEntry("MyMeasures_" + other_name + "_range");
  ctl.logger().removeLogEntry("MyMeasures_" + other_name + "_ground");
  ctl.logger().removeLogEntry("MyMeasures_" + other_name + "_filtered_x");
  ctl.logger().removeLogEntry("MyMeasures_" + other_name + "_convex_x");
  ctl.logger().removeLogEntry("MyMeasures_" + other_name + "_filtered_y");
  ctl.logger().removeLogEntry("MyMeasures_" + other_name + "_convex_y");
  ctl.logger().removeLogEntry("MyMeasures_" + other_name + "_k");
  ctl.logger().removeLogEntry("MyMeasures_" + other_name + "_angle");
  ctl.logger().removeLogEntry("MyMeasures_" + other_name + "_min_max_phalanxes_angle");
  ctl.logger().removeLogEntry("MyMeasures_" + other_name + "_position_offset");
  ctl.logger().removeLogEntry("MyMeasures_" + other_name + "_phalanxes_x");
  ctl.logger().removeLogEntry("MyMeasures_" + other_name + "_phalanxes_y");



  // Reset GUI
  ctl.gui()->removeElement({"SoftFoot"}, other_name + "_point_ground");
  ctl.gui()->removeElement({"SoftFoot"}, other_name + "_ground");

  ctl.gui()->removeElement({"SoftFoot"}, other_name + "_segment");
  ctl.gui()->removeElement({"SoftFoot"}, name + "_segment");

  ctl.gui()->removeElement({"SoftFoot"}, other_name + "_segment_filtered");
  ctl.gui()->removeElement({"SoftFoot"}, name + "_segment_filtered");

  ctl.gui()->removeElement({"SoftFoot"}, other_name + "_convex_segment");
  ctl.gui()->removeElement({"SoftFoot"}, name + "_convex_segment");

  ctl.gui()->removeElement({"SoftFoot"}, other_name + "_convex_display");
  ctl.gui()->removeElement({"SoftFoot"}, name + "_convex_display");

  // Reset log
}

void SoftFootState::rightFRSCallback(const std_msgs::Float64::ConstPtr& data)
{
  // Select data/string based on current_moving_foot
  const std::string sensor_name = range_sensor_name_[Foot::Right];
  const std::lock_guard<std::mutex> lock(range_sensor_mutex_);
  ctl().robot().device<mc_mujoco::RangeSensor>(sensor_name).update(data->data * 0.001);
}

void SoftFootState::leftFRSCallback(const std_msgs::Float64::ConstPtr& data)
{
  // Select data/string based on current_moving_foot
  const std::string sensor_name = range_sensor_name_[Foot::Left];
  const std::lock_guard<std::mutex> lock(range_sensor_mutex_);
  ctl().robot().device<mc_mujoco::RangeSensor>(sensor_name).update(data->data * 0.001);
}

EXPORT_SINGLE_STATE("BWC::SoftFoot", SoftFootState)
