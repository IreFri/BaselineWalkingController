#pragma once

#include <BaselineWalkingController/State.h>

#include <mc_rtc/ros.h>
#include <ros/ros.h>
#include <std_msgs/Float64.h>
#include <mutex>

namespace BWC
{
/** \brief FSM state to calculate ground profile and calculate sole stiffness online. */
struct SoftFootState : State
{
public:

  /** \brief Start. */
  void start(mc_control::fsm::Controller & ctl) override;

  /** \brief Run. */
  bool run(mc_control::fsm::Controller & ctl) override;

  /** \brief Teardown. */
  void teardown(mc_control::fsm::Controller & ctl) override;

protected:

  void calculateCost(mc_control::fsm::Controller & ctl);

  void estimateGround(mc_control::fsm::Controller & ctl, const Foot & current_moving_foot);

  void extractGroundSegment(mc_control::fsm::Controller & ctl, const Foot & current_moving_foot, const Eigen::Vector3d & landing);

  void extractAltitudeProfileFromGroundSegment(const Foot & current_moving_foot);

  void updateVariableStiffness(mc_control::fsm::Controller & ctl, const Foot & current_moving_foot);

  void computeSegmentConvexHull(mc_control::fsm::Controller & ctl, const Foot & current_moving_foot);

  std::vector<Eigen::Vector3d> computeConvexHull(const std::vector<Eigen::Vector3d>& data);

  void computeMinMaxAngle(const Foot & current_moving_foot);

  void computeFootLandingPosition(const Foot & current_moving_foot, const Eigen::Vector3d & landing);

  void computeFootLandingAngle(const Foot & current_moving_foot, const Eigen::Vector3d & landing);

  void updateFootSwingPose(mc_control::fsm::Controller & ctl, const Foot & current_moving_foot, const sva::PTransformd & X_0_landing);

  Foot getCurrentMovingFoot(mc_control::fsm::Controller & ctl) const
  {
    // Cast ctl to BaselineWalkingController
    auto & ctrl = static_cast<BWC::BaselineWalkingController&>(ctl);

    Foot current_moving_foot; 
    if(ctrl.footManager_->supportPhase() == BWC::SupportPhase::LeftSupport) // This is the name "LeftFootCenter" or "RightFootCenter"
    {
      current_moving_foot = Foot::Right;
    }
    else
    {
      current_moving_foot = Foot::Left;
    }
    return current_moving_foot;
  }

  bool evaluateTheta(double theta, double direction, const Eigen::Vector3d& pos, const std::vector<Eigen::Vector3d>& ground, Eigen::Vector2d& output);

protected:
  // Used to compute the cost to understand which is the best stiffness
  double lambda_zmp_ = 3.0;
  std::vector<double> zmp_;

  double lambda_CoM_ = 3.0;
  std::vector<double> CoM_; 

  double nr_remaining_footstep_ = 0.0;
  double nr_footstep_ = 0.0;
  double lambda_footstep_ = 10.0;

  double cost_ = 0.0;
  double PhalangesStiffness_ = 0.0;
  double extra_to_compute_best_position_ = 0.07;

  // FootData contains data used to estimate the ground profile
  struct FootData
  {
    double range;
    std::vector<Eigen::Vector3d> ground;
    std::vector<std::array<double, 2>> phalanxes;
    std::vector<double> altitude;
    double k;
    double angle;
    double min_max_phalanxes_angle;
    double position_offset_x;
    double position_offset_z;
    bool need_reset;
    bool computation_done;
  };
  std::unordered_map<Foot, FootData> foot_data_;
  
  // Data for GroundSegment to compute best position/orientation
  struct GroundSegment
  {
    std::vector<Eigen::Vector3d> raw; // means all the data
    std::vector<Eigen::Vector3d> filtered;
    std::vector<Eigen::Vector3d> convex;
  };
  std::unordered_map<Foot, GroundSegment> ground_segment_; 

  std::unordered_map<Foot, std::string> variable_stiffness_jointname_ = {
    {Foot::Left, "L_VARSTIFF"},
    {Foot::Right, "R_VARSTIFF"}
  };

  std::unordered_map<Foot, std::string> range_sensor_name_ = {
    {Foot::Left, "LeftFootRangeSensor"},
    {Foot::Right, "RightFootRangeSensor"}
  };

  std::unordered_map<Foot, std::string> surface_name_ = {
    {Foot::Left, "LeftFootCenter"},
    {Foot::Right, "RightFootCenter"}
  };

  // TODO: Ugly hardcoded value
  double foot_length_ = 0.27742;
  size_t nr_phalanxes_;
  double phalanx_length_;

  // Reset data
  void reset(mc_control::fsm::Controller & ctl, const Foot & foot);
 
  // Client is here to call the service to compute the stiffness based on the ground profile
  ros::ServiceClient client_; 

  // options
  bool with_variable_stiffness_ = false;
  bool with_ankle_rotation_ = false;
  bool with_foot_adjustment_ = false;

  // Subscriber to update the range sensor
  ros::Subscriber right_foot_range_sensor_sub_;
  void rightFRSCallback(const std_msgs::Float64::ConstPtr& data);

  ros::Subscriber left_foot_range_sensor_sub_;
  void leftFRSCallback(const std_msgs::Float64::ConstPtr& data);

  std::mutex range_sensor_mutex_;
};
} // namespace BWC
