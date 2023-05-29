#pragma once

#include <BaselineWalkingController/State.h>

#include <mc_rtc/ros.h>
#include <ros/ros.h>

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

  void computeFootLandingAngle(const Foot & current_moving_foot, const Eigen::Vector3d & landing);

  void updateFootSwingPose(mc_control::fsm::Controller & ctl, const Foot & current_moving_foot, const sva::PTransformd & X_0_landing);

protected:
  // Used to compute the cost to understand which is the best stiffness
  double lambda_zmp_ = 3.0;
  std::vector<double> zmp_;

  double lambda_CoM_ = 3.0;
  std::vector<double> CoM_; 

  double footstep_error_ = 0.0;
  double lambda_footstep_ = 10.0;

  double cost_ = 0.0;
  double PhalangesStiffness_ = 0.0;

  // FootData contains data used to estimate the ground profile
  struct FootData
  {
    double range;
    std::vector<Eigen::Vector3d> ground;
    std::vector<double> altitude;
    double profile;
    double profileFiltered;
    double k;
    double angle;
    bool areComputationDone;
    bool needReset;
  };
  std::unordered_map<Foot, FootData> foot_data_;
  
  // Data for GroundSegment to compute best position/orientation
  struct GroundSegment
  {
    std::vector<Eigen::Vector3d> raw; // means all the data
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

  double foot_length_ = 0.27742;

  // Reset data
  void reset(mc_control::fsm::Controller & ctl, const Foot & foot);
 
  // Client is here to call the service to compute the stiffness based on the ground profile
  ros::ServiceClient client_; 
};
} // namespace BWC
