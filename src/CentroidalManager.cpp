#include <RBDyn/Momentum.h>

#include <mc_rtc/gui/Checkbox.h>
#include <mc_rtc/gui/Label.h>
#include <mc_rtc/gui/NumberInput.h>
#include <mc_tasks/CoMTask.h>
#include <mc_tasks/FirstOrderImpedanceTask.h>

#include <CCC/Constants.h>
#include <ForceColl/WrenchDistribution.h>

#include <BaselineWalkingController/BaselineWalkingController.h>
#include <BaselineWalkingController/CentroidalManager.h>
#include <BaselineWalkingController/FootManager.h>

using namespace BWC;

void CentroidalManager::Configuration::load(const mc_rtc::Configuration & mcRtcConfig)
{
  mcRtcConfig("name", name);
  mcRtcConfig("method", method);
  mcRtcConfig("useActualStateForMpc", useActualStateForMpc);
  mcRtcConfig("enableZmpFeedback", enableZmpFeedback);
  mcRtcConfig("enableComZFeedback", enableComZFeedback);
  mcRtcConfig("dcmGainP", dcmGainP);
  mcRtcConfig("zmpVelGain", zmpVelGain);
  mcRtcConfig("comZGainP", comZGainP);
  mcRtcConfig("comZGainD", comZGainD);
  mcRtcConfig("refComZ", refComZ);
  mcRtcConfig("useTargetPoseForControlRobotAnchorFrame", useTargetPoseForControlRobotAnchorFrame);
  mcRtcConfig("useActualComForWrenchDist", useActualComForWrenchDist);
  mcRtcConfig("wrenchDistConfig", wrenchDistConfig);
}

CentroidalManager::CentroidalManager(BaselineWalkingController * ctlPtr, const mc_rtc::Configuration & // mcRtcConfig
                                     )
: ctlPtr_(ctlPtr)
{
}

void CentroidalManager::reset()
{
  robotMass_ = ctl().robot().mass();
}

void CentroidalManager::update()
{
  // Check if enable
  if(!enabled_)
  {
    return;
  }

  // Check if the robot is in the air
  bool isInTheAir = true;
  for(const auto & foot : Feet::Both)
  {
    isInTheAir = isInTheAir && ctl().footTasks_.at(foot)->measuredWrench().force().z() < 10;
  }
  
  if(isInTheAir)
  {
    mc_rtc::log::warning("Robot is in the air");
    ctl().footManager_->disable();
    wasInTheAir_ = true;
    return;
  }
  else if(wasInTheAir_)
  {
    ctl().footManager_->reset();
    ctl().footManager_->enable();
    wasInTheAir_ = false;
  }
  
  // Set MPC state
  if(config().useActualStateForMpc)
  {
    mpcCom_ = ctl().realRobot().com();
    mpcComVel_ = ctl().realRobot().comVelocity();
  }
  else
  {
    // Task targets are the planned state in the previous step
    mpcCom_ = ctl().comTask_->com();
    mpcComVel_ = ctl().comTask_->refVel();
  }
  refZmp_ = ctl().footManager_->calcRefZmp(ctl().t());

  // Run MPC
  runMpc();

  // Calculate target wrench
  {
    controlZmp_ = plannedZmp_;
    controlForceZ_ = plannedForceZ_;

    // Compensate ZMP delay
    // See equation (10) of https://ieeexplore.ieee.org/abstract/document/6094838
    Eigen::Vector3d refZmpVel = ctl().footManager_->calcRefZmp(ctl().t(), 1);
    controlZmp_.head<2>() += config().zmpVelGain * refZmpVel.head<2>();

    // Apply DCM feedback
    if(config().enableZmpFeedback)
    {
      double omega = std::sqrt(plannedForceZ_ / (robotMass_ * (mpcCom_.z() - refZmp_.z())));
      Eigen::Vector3d plannedDcm = ctl().comTask_->com() + ctl().comTask_->refVel() / omega;
      Eigen::Vector3d actualDcm = ctl().realRobot().com() + ctl().realRobot().comVelocity() / omega;
      controlZmp_.head<2>() += config().dcmGainP * (actualDcm - plannedDcm).head<2>();
    }

    // Apply ForceZ feedback
    if(config().enableComZFeedback)
    {
      double plannedComZ = ctl().comTask_->com().z();
      double actualComZ = ctl().realRobot().com().z();
      double plannedComVelZ = ctl().comTask_->refVel().z();
      double actualComVelZ = ctl().realRobot().comVelocity().z();
      controlForceZ_ -=
          config().comZGainP * (actualComZ - plannedComZ) + config().comZGainD * (actualComVelZ - plannedComVelZ);
    }

    // Convert ZMP to wrench and distribute
    contactList_ = ctl().footManager_->calcCurrentContactList();
    wrenchDist_ = std::make_shared<ForceColl::WrenchDistribution>(ForceColl::getContactVecFromMap(contactList_),
                                                                  config().wrenchDistConfig);
    Eigen::Vector3d comForWrenchDist =
        (config().useActualComForWrenchDist ? ctl().realRobot().com() : ctl().comTask_->com());
    sva::ForceVecd controlWrench;
    controlWrench.force() << controlForceZ_ / (comForWrenchDist.z() - refZmp_.z())
                                 * (comForWrenchDist.head<2>() - controlZmp_.head<2>()),
        controlForceZ_;
    controlWrench.moment().setZero(); // Moment is represented around CoM
    wrenchDist_->run(controlWrench, comForWrenchDist);
  }

  // Set target of tasks
  {
    // Set target of CoM task
    Eigen::Vector3d plannedComAccel = calcPlannedComAccel();
    Eigen::Vector3d nextPlannedCom =
        mpcCom_ + ctl().dt() * mpcComVel_ + 0.5 * std::pow(ctl().dt(), 2) * plannedComAccel;
    Eigen::Vector3d nextPlannedComVel = mpcComVel_ + ctl().dt() * plannedComAccel;
    if(isConstantComZ())
    {
      nextPlannedCom.z() = config().refComZ + ctl().footManager_->calcRefGroundPosZ(ctl().t());
      nextPlannedComVel.z() = ctl().footManager_->calcRefGroundPosZ(ctl().t(), 1);
      plannedComAccel.z() = ctl().footManager_->calcRefGroundPosZ(ctl().t(), 2);
    }
    ctl().comTask_->com(nextPlannedCom);
    ctl().comTask_->refVel(nextPlannedComVel);
    ctl().comTask_->refAccel(plannedComAccel);

    // Set target wrench of foot tasks
    const auto & targetWrenchList = ForceColl::calcWrenchList(contactList_, wrenchDist_->resultWrenchRatio_);
    for(const auto & foot : Feet::Both)
    {
      sva::ForceVecd targetWrench = sva::ForceVecd::Zero();
      if(targetWrenchList.count(foot))
      {
        targetWrench = targetWrenchList.at(foot);
      }
      ctl().footTasks_.at(foot)->targetWrenchW(targetWrench);
    }
  }

  // Update force visualization
  {
    ctl().gui()->removeCategory({ctl().name(), config().name, "ForceMarker"});
    wrenchDist_->addToGUI(*ctl().gui(), {ctl().name(), config().name, "ForceMarker"});
  }
}

void CentroidalManager::stop()
{
  removeFromGUI(*ctl().gui());
  removeFromLogger(ctl().logger());
}

void CentroidalManager::addToGUI(mc_rtc::gui::StateBuilder & gui)
{
  gui.addElement(
      {ctl().name(), config().name, "Config"}, mc_rtc::gui::Label("method", [this]() { return config().method; }),
      mc_rtc::gui::Checkbox(
          "useActualStateForMpc", [this]() { return config().useActualStateForMpc; },
          [this]() { config().useActualStateForMpc = !config().useActualStateForMpc; }),
      mc_rtc::gui::Checkbox(
          "enableZmpFeedback", [this]() { return config().enableZmpFeedback; },
          [this]() { config().enableZmpFeedback = !config().enableZmpFeedback; }),
      mc_rtc::gui::Checkbox(
          "enableComZFeedback", [this]() { return config().enableComZFeedback; },
          [this]() { config().enableComZFeedback = !config().enableComZFeedback; }),
      mc_rtc::gui::NumberInput(
          "dcmGainP", [this]() { return config().dcmGainP; }, [this](double v) { config().dcmGainP = v; }),
      mc_rtc::gui::NumberInput(
          "zmpVelGain", [this]() { return config().zmpVelGain; }, [this](double v) { config().zmpVelGain = v; }),
      mc_rtc::gui::NumberInput(
          "comZGainP", [this]() { return config().comZGainP; }, [this](double v) { config().comZGainP = v; }),
      mc_rtc::gui::NumberInput(
          "comZGainD", [this]() { return config().comZGainD; }, [this](double v) { config().comZGainD = v; }),
      mc_rtc::gui::NumberInput(
          "refComZ", [this]() { return config().refComZ; }, [this](double v) { config().refComZ = v; }),
      mc_rtc::gui::Checkbox(
          "useTargetPoseForControlRobotAnchorFrame",
          [this]() { return config().useTargetPoseForControlRobotAnchorFrame; },
          [this]() {
            config().useTargetPoseForControlRobotAnchorFrame = !config().useTargetPoseForControlRobotAnchorFrame;
          }),
      mc_rtc::gui::Checkbox(
          "useActualComForWrenchDist", [this]() { return config().useActualComForWrenchDist; },
          [this]() { config().useActualComForWrenchDist = !config().useActualComForWrenchDist; }));

    using Style = mc_rtc::gui::plot::Style;
    using Side = mc_rtc::gui::plot::Side;
    gui.addElement({ctl().name(), config().name, "Debug"}, mc_rtc::gui::ElementsStacking::Horizontal,
      mc_rtc::gui::Button("Plot DCM-ZMP Tracking (x)",
            [this, &gui]()
            {
              gui.addPlot(
                  "DCM-ZMP Tracking (x)", mc_rtc::gui::plot::X("t", [this]() { static double t = 0.; return t += ctl().solver().dt(); }),
                  mc_rtc::gui::plot::Y(
                      "support_min",
                      [this]()
                      {
                        Eigen::Vector2d minPos = Eigen::Vector2d::Constant(std::numeric_limits<double>::max());
                        for(const auto & contactKV : contactList_)
                        {
                          for(const auto & vertexWithRidge : contactKV.second->vertexWithRidgeList_)
                          {
                            minPos = minPos.cwiseMin(vertexWithRidge.vertex.head<2>());
                          }
                        }
                        return minPos.x();
                      }, mc_rtc::gui::Color::Red),
                  mc_rtc::gui::plot::Y(
                      "support_max",
                      [this]()
                      {
                        Eigen::Vector2d maxPos = Eigen::Vector2d::Constant(std::numeric_limits<double>::lowest());
                        for(const auto & contactKV : contactList_)
                        {
                          for(const auto & vertexWithRidge : contactKV.second->vertexWithRidgeList_)
                          {
                            maxPos = maxPos.cwiseMax(vertexWithRidge.vertex.head<2>());
                          }
                        }
                        return maxPos.x();
                      }, mc_rtc::gui::Color::Red),
                  mc_rtc::gui::plot::Y(
                      "refZmp", [this]() { return refZmp_.x(); }, mc_rtc::gui::Color::Cyan),
                  mc_rtc::gui::plot::Y(
                      "plannedZmp", [this]() { return plannedZmp_.x(); }, mc_rtc::gui::Color::Cyan, Style::Dashed),
                  mc_rtc::gui::plot::Y(
                      "measuredZmp",
                      [this]()
                      {
                        std::unordered_map<Foot, sva::ForceVecd> sensorWrenchList;
                        for(const auto & foot : ctl().footManager_->getCurrentContactFeet())
                        {
                          const auto & surfaceName = ctl().footManager_->surfaceName(foot);
                          const auto & sensorName = ctl().robot().indirectSurfaceForceSensor(surfaceName).name();
                          const auto & sensor = ctl().robot().forceSensor(sensorName);
                          const auto & sensorWrench = sensor.worldWrenchWithoutGravity(ctl().robot());
                          sensorWrenchList.emplace(foot, sensorWrench);
                        }
                        return calcZmp(sensorWrenchList, refZmp_.z()).x();
                      }, mc_rtc::gui::Color::Blue, Style::Dashed),
                  mc_rtc::gui::plot::Y(
                      "controlZmp", [this]() { return controlZmp_.x(); }, mc_rtc::gui::Color::Blue),
                  mc_rtc::gui::plot::Y(
                      "plannedDcm",
                      [this]()
                      {
                        double omega = std::sqrt(plannedForceZ_ / (robotMass_ * (mpcCom_.z() - refZmp_.z())));
                        return (ctl().comTask_->com() + ctl().comTask_->refVel() / omega).x();
                      }, mc_rtc::gui::Color::Magenta, Style::Solid, Side::Left),
                  mc_rtc::gui::plot::Y(
                      "actualDcm",
                      [this]()
                      {
                        double omega = std::sqrt(plannedForceZ_ / (robotMass_ * (mpcCom_.z() - refZmp_.z())));
                        return (ctl().realRobot().com() + ctl().realRobot().comVelocity() / omega).x();
                      }, mc_rtc::gui::Color::Magenta, Style::Dashed, Side::Left)
              );
            }),
      mc_rtc::gui::Button("Stop DCM-ZMP (x)", [&gui]() { gui.removePlot("DCM-ZMP Tracking (x)"); }));

  gui.addElement({ctl().name(), config().name, "Debug"}, mc_rtc::gui::ElementsStacking::Horizontal,
      mc_rtc::gui::Button("Plot DCM-ZMP Tracking (y)",
            [this, &gui]()
            {
              gui.addPlot(
                  "DCM-ZMP Tracking (y)", mc_rtc::gui::plot::X("t", [this]() { static double t = 0.; return t += ctl().solver().dt(); }),
                  mc_rtc::gui::plot::Y(
                      "support_min",
                      [this]()
                      {
                        Eigen::Vector2d minPos = Eigen::Vector2d::Constant(std::numeric_limits<double>::max());
                        for(const auto & contactKV : contactList_)
                        {
                          for(const auto & vertexWithRidge : contactKV.second->vertexWithRidgeList_)
                          {
                            minPos = minPos.cwiseMin(vertexWithRidge.vertex.head<2>());
                          }
                        }
                        return minPos.y();
                      }, mc_rtc::gui::Color::Red),
                  mc_rtc::gui::plot::Y(
                      "support_max",
                      [this]()
                      {
                        Eigen::Vector2d maxPos = Eigen::Vector2d::Constant(std::numeric_limits<double>::lowest());
                        for(const auto & contactKV : contactList_)
                        {
                          for(const auto & vertexWithRidge : contactKV.second->vertexWithRidgeList_)
                          {
                            maxPos = maxPos.cwiseMax(vertexWithRidge.vertex.head<2>());
                          }
                        }
                        return maxPos.y();
                      }, mc_rtc::gui::Color::Red),
                  mc_rtc::gui::plot::Y(
                      "refZmp", [this]() { return refZmp_.y(); }, mc_rtc::gui::Color::Cyan),
                  mc_rtc::gui::plot::Y(
                      "plannedZmp", [this]() { return plannedZmp_.y(); }, mc_rtc::gui::Color::Cyan, Style::Dashed),
                  mc_rtc::gui::plot::Y(
                      "measuredZmp",
                      [this]()
                      {
                        std::unordered_map<Foot, sva::ForceVecd> sensorWrenchList;
                        for(const auto & foot : ctl().footManager_->getCurrentContactFeet())
                        {
                          const auto & surfaceName = ctl().footManager_->surfaceName(foot);
                          const auto & sensorName = ctl().robot().indirectSurfaceForceSensor(surfaceName).name();
                          const auto & sensor = ctl().robot().forceSensor(sensorName);
                          const auto & sensorWrench = sensor.worldWrenchWithoutGravity(ctl().robot());
                          sensorWrenchList.emplace(foot, sensorWrench);
                        }
                        return calcZmp(sensorWrenchList, refZmp_.z()).y();
                      }, mc_rtc::gui::Color::Blue, Style::Dashed),
                  mc_rtc::gui::plot::Y(
                      "controlZmp", [this]() { return controlZmp_.y(); }, mc_rtc::gui::Color::Blue),
                  mc_rtc::gui::plot::Y(
                      "plannedDcm",
                      [this]()
                      {
                        double omega = std::sqrt(plannedForceZ_ / (robotMass_ * (mpcCom_.z() - refZmp_.z())));
                        return (ctl().comTask_->com() + ctl().comTask_->refVel() / omega).y();
                      }, mc_rtc::gui::Color::Magenta, Style::Solid, Side::Left),
                  mc_rtc::gui::plot::Y(
                      "actualDcm",
                      [this]()
                      {
                        double omega = std::sqrt(plannedForceZ_ / (robotMass_ * (mpcCom_.z() - refZmp_.z())));
                        return (ctl().realRobot().com() + ctl().realRobot().comVelocity() / omega).y();
                      }, mc_rtc::gui::Color::Magenta, Style::Dashed, Side::Left)
              );
            }),
      mc_rtc::gui::Button("Stop DCM-ZMP (y)", [&gui]() { gui.removePlot("DCM-ZMP Tracking (y)"); }));

  gui.addElement({ctl().name(), config().name, "Debug"}, mc_rtc::gui::ElementsStacking::Horizontal,
      mc_rtc::gui::Button("Plot CoM Tracking (x)",
            [this, &gui]()
            {
              gui.addPlot("CoM Tracking (x)", mc_rtc::gui::plot::X("t", [this]() { static double t = 0.; return t += ctl().solver().dt(); }),
                mc_rtc::gui::plot::Y(
                    "com_ref", [this]() { return ctl().robot().com().x(); }, mc_rtc::gui::Color::Red),
                mc_rtc::gui::plot::Y(
                    "com_mes", [this]() { return ctl().realRobot().com().x(); }, mc_rtc::gui::Color::Magenta));
            }),
      mc_rtc::gui::Button("Stop CoM (x)", [&gui]() { gui.removePlot("CoM Tracking (x)"); }));

  gui.addElement({ctl().name(), config().name, "Debug"}, mc_rtc::gui::ElementsStacking::Horizontal,
      mc_rtc::gui::Button("Plot CoM Tracking (y)",
            [this, &gui]()
            {
              gui.addPlot("CoM Tracking (y)", mc_rtc::gui::plot::X("t", [this]() { static double t = 0.; return t += ctl().solver().dt(); }),
              mc_rtc::gui::plot::Y(
                    "com_ref", [this]() { return ctl().robot().com().y(); }, mc_rtc::gui::Color::Red),
                mc_rtc::gui::plot::Y(
                    "com_mes", [this]() { return ctl().realRobot().com().y(); }, mc_rtc::gui::Color::Magenta));
            }),
      mc_rtc::gui::Button("Stop CoM (y)", [&gui]() { gui.removePlot("CoM Tracking (y)"); }));
}

void CentroidalManager::removeFromGUI(mc_rtc::gui::StateBuilder & gui)
{
  gui.removeCategory({ctl().name(), config().name});
}

void CentroidalManager::addToLogger(mc_rtc::Logger & logger)
{
  logger.addLogEntry(config().name + "_Config_method", this, [this]() { return config().method; });
  logger.addLogEntry(config().name + "_Config_useActualStateForMpc", this,
                     [this]() { return config().useActualStateForMpc; });
  logger.addLogEntry(config().name + "_Config_enableZmpFeedback", this,
                     [this]() { return config().enableZmpFeedback; });
  logger.addLogEntry(config().name + "_Config_enableComZFeedback", this,
                     [this]() { return config().enableComZFeedback; });
  logger.addLogEntry(config().name + "_Config_dcmGainP", this, [this]() { return config().dcmGainP; });
  logger.addLogEntry(config().name + "_Config_zmpVelGain", this, [this]() { return config().zmpVelGain; });
  logger.addLogEntry(config().name + "_Config_comZGainP", this, [this]() { return config().comZGainP; });
  logger.addLogEntry(config().name + "_Config_comZGainD", this, [this]() { return config().comZGainD; });
  logger.addLogEntry(config().name + "_Config_refComZ", this, [this]() { return config().refComZ; });
  logger.addLogEntry(config().name + "_Config_useTargetPoseForControlRobotAnchorFrame", this,
                     [this]() { return config().useTargetPoseForControlRobotAnchorFrame; });
  logger.addLogEntry(config().name + "_Config_useActualComForWrenchDist", this,
                     [this]() { return config().useActualComForWrenchDist; });

  MC_RTC_LOG_HELPER(config().name + "_CoM_MPC", mpcCom_);
  logger.addLogEntry(config().name + "_CoM_planned", this, [this]() { return ctl().comTask_->com(); });
  logger.addLogEntry(config().name + "_CoM_controlRobot", this, [this]() { return ctl().robot().com(); });
  logger.addLogEntry(config().name + "_CoM_realRobot", this, [this]() { return ctl().realRobot().com(); });

  MC_RTC_LOG_HELPER(config().name + "_forceZ_planned", plannedForceZ_);
  MC_RTC_LOG_HELPER(config().name + "_forceZ_control", controlForceZ_);

  logger.addLogEntry(config().name + "_ZMP_ref", this, [this]() { return refZmp_; });
  MC_RTC_LOG_HELPER(config().name + "_ZMP_planned", plannedZmp_);
  MC_RTC_LOG_HELPER(config().name + "_ZMP_control", controlZmp_);
  logger.addLogEntry(config().name + "_ZMP_controlWrenchDist", this, [this]() {
    return wrenchDist_ ? calcZmp(ForceColl::calcWrenchList(contactList_, wrenchDist_->resultWrenchRatio_), refZmp_.z())
                       : Eigen::Vector3d::Zero();
  });
  logger.addLogEntry(config().name + "_ZMP_measured", this, [this]() {
    std::unordered_map<Foot, sva::ForceVecd> sensorWrenchList;
    for(const auto & foot : ctl().footManager_->getCurrentContactFeet())
    {
      const auto & surfaceName = ctl().footManager_->surfaceName(foot);
      const auto & sensorName = ctl().robot().indirectSurfaceForceSensor(surfaceName).name();
      const auto & sensor = ctl().robot().forceSensor(sensorName);
      const auto & sensorWrench = sensor.worldWrenchWithoutGravity(ctl().robot());
      sensorWrenchList.emplace(foot, sensorWrench);
    }
    return calcZmp(sensorWrenchList, refZmp_.z());
  });
  logger.addLogEntry(config().name + "_ZMP_SupportRegion_min", this, [this]() {
    Eigen::Vector2d minPos = Eigen::Vector2d::Constant(std::numeric_limits<double>::max());
    for(const auto & contactKV : contactList_)
    {
      for(const auto & vertexWithRidge : contactKV.second->vertexWithRidgeList_)
      {
        minPos = minPos.cwiseMin(vertexWithRidge.vertex.head<2>());
      }
    }
    return minPos;
  });
  logger.addLogEntry(config().name + "_ZMP_SupportRegion_max", this, [this]() {
    Eigen::Vector2d maxPos = Eigen::Vector2d::Constant(std::numeric_limits<double>::lowest());
    for(const auto & contactKV : contactList_)
    {
      for(const auto & vertexWithRidge : contactKV.second->vertexWithRidgeList_)
      {
        maxPos = maxPos.cwiseMax(vertexWithRidge.vertex.head<2>());
      }
    }
    return maxPos;
  });

  logger.addLogEntry(config().name + "_CentroidalMomentum_controlRobot", this, [this]() {
    return rbd::computeCentroidalMomentum(ctl().robot().mb(), ctl().robot().mbc(), ctl().robot().com());
  });
}

void CentroidalManager::removeFromLogger(mc_rtc::Logger & logger)
{
  logger.removeLogEntries(this);
}

void CentroidalManager::setAnchorFrame()
{
  std::string anchorName = "KinematicAnchorFrame::" + ctl().robot().name();
  if(ctl().datastore().has(anchorName))
  {
    ctl().datastore().remove(anchorName);
  }
  ctl().datastore().make_call(anchorName, [this](const mc_rbdyn::Robot & robot) { return calcAnchorFrame(robot); });
}

sva::PTransformd CentroidalManager::calcAnchorFrame(const mc_rbdyn::Robot & robot) const
{
  double leftFootSupportRatio = ctl().footManager_->leftFootSupportRatio();
  bool isControlRobot = (&(ctl().robot()) == &robot);

  if(isControlRobot && config().useTargetPoseForControlRobotAnchorFrame)
  {
    return sva::interpolate(ctl().footManager_->targetFootPose(Foot::Right),
                            ctl().footManager_->targetFootPose(Foot::Left), leftFootSupportRatio);
  }
  else
  {
    return sva::interpolate(robot.surfacePose(ctl().footManager_->surfaceName(Foot::Right)),
                            robot.surfacePose(ctl().footManager_->surfaceName(Foot::Left)), leftFootSupportRatio);
  }
}

Eigen::Vector3d CentroidalManager::calcZmp(const std::unordered_map<Foot, sva::ForceVecd> & wrenchList,
                                           double zmpPlaneHeight,
                                           const Eigen::Vector3d & zmpPlaneNormal) const
{
  sva::ForceVecd totalWrench = sva::ForceVecd::Zero();
  for(const auto & wrenchKV : wrenchList)
  {
    totalWrench += wrenchKV.second;
  }

  Eigen::Vector3d zmpPlaneOrigin = Eigen::Vector3d(0, 0, zmpPlaneHeight);
  Eigen::Vector3d zmp = zmpPlaneOrigin;

  if(totalWrench.force().z() > 0)
  {
    Eigen::Vector3d momentInZmpPlane = totalWrench.moment() - zmpPlaneOrigin.cross(totalWrench.force());
    zmp += zmpPlaneNormal.cross(momentInZmpPlane) / totalWrench.force().z();
  }

  return zmp;
}

Eigen::Vector3d CentroidalManager::calcPlannedComAccel() const
{
  Eigen::Vector3d plannedComAccel;
  plannedComAccel << plannedForceZ_ / (robotMass_ * (mpcCom_.z() - refZmp_.z()))
                         * (mpcCom_.head<2>() - plannedZmp_.head<2>()),
      plannedForceZ_ / robotMass_;
  plannedComAccel.z() -= CCC::constants::g;
  return plannedComAccel;
}
