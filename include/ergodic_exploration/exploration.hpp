/**
 * @file exploration.hpp
 * @author Boston Cleek
 * @date 4 Dec 2020
 * @brief Ergodic exploration with a user defined target and dynamic window
 * approach for collision avoidance
 */

#include <iostream>

#include <ros/ros.h>
#include <ros/console.h>
#include <geometry_msgs/Twist.h>
#include <geometry_msgs/TransformStamped.h>
#include <nav_msgs/OccupancyGrid.h>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/LaserScan.h>
#include <tf2_ros/transform_listener.h>
#include <tf2/LinearMath/Quaternion.h>
#include <tf2/utils.h>
#include <gazebo_msgs/ModelStates.h>

#include <ergodic_exploration/ergodic_control.hpp>
#include <ergodic_exploration/dynamic_window.hpp>

// TODO: copy constructors for dwa, collision

namespace ergodic_exploration
{
constexpr char LOGNAME[] = "ergodic exploration";

using arma::eye;
using arma::mat;
using arma::span;
using arma::vec;

/** @brief Exploration template */
template <class ModelT>
class Exploration
{
public:
  /**
   * @brief Constructor
   * @param nh - NodeHandle
   * @param ergodic_control - ergodic controller
   * @param collision - collision detector
   * @param dwa - dynamic window local planner
   */
  Exploration(ros::NodeHandle& nh, ErgodicControl<ModelT>& ergodic_control,
              const Collision& collision, const DynamicWindow& dwa);

  /** @brief Initialize subscribers and publishers */
  void init();

  /**
   * @brief Odometry callback
   * @param msg - odometry message
   * @details updates the robot's twist in the body frame [vx, vy, w]
   */
  void odomCallback(const nav_msgs::Odometry& msg);

  /**
   * @brief Occupancy grid callback
   * @param msg - occupancy grid  message
   */
  void mapCallback(const nav_msgs::OccupancyGrid::ConstPtr& msg);

  /**
   * @brief Executes the exploration stack
   * @param target - target distribution
   * @param map_frame_id - map frame
   * @param base_frame_id - robot's base link frame
   * @param frequency - control loop frequency
   * @param dt - time step in integration
   * @param horizon - length of integration
   * @details publishes a body twist at a fixed frequency
   */
  void control(const Target& target, const std::string& map_frame_id,
               const std::string& base_frame_id, double frequency, double val_dt,
               double val_horizon);

private:
  ros::NodeHandle nh_;  // NodeHandle

  ErgodicControl<ModelT> ergodic_control_;  // ergodic controller
  Collision collision_;                     // collision detector
  DynamicWindow dwa_;                       // dynamic window local planner

  bool map_received_;  // occupancy grid received

  vec vb_;        // current body twist [vx, vy, w]
  GridMap grid_;  // occupancy grid

  tf2_ros::Buffer tfBuffer_;               // tf2 buffer
  tf2_ros::TransformListener tfListener_;  // tf2 listener

  ros::Subscriber map_sub_;   // occupancy grid subscriber
  ros::Subscriber odom_sub_;  // odometry grid subscriber

  ros::Publisher cmd_pub_;       // twist publisher
  ros::Publisher opt_traj_pub_;  // ergodic controller trajectory publisher
  ros::Publisher dwa_path_pub_;  // dynamic window trajectory publisher
  ros::Publisher target_pub_;    // target distribution publisher
};

template <class ModelT>
Exploration<ModelT>::Exploration(ros::NodeHandle& nh,
                                 ErgodicControl<ModelT>& ergodic_control,
                                 const Collision& collision, const DynamicWindow& dwa)
  : nh_(nh)
  , ergodic_control_(ergodic_control)
  , collision_(collision)
  , dwa_(dwa)
  , map_received_(false)
  , vb_(3, arma::fill::zeros)
  , tfListener_(tfBuffer_)
{
  init();
}

template <class ModelT>
void Exploration<ModelT>::init()
{
  map_sub_ = nh_.subscribe("map", 1, &Exploration<ModelT>::mapCallback, this);
  odom_sub_ = nh_.subscribe("odom", 1, &Exploration<ModelT>::odomCallback, this);

  cmd_pub_ = nh_.advertise<geometry_msgs::Twist>("cmd_vel", 1);
  opt_traj_pub_ = nh_.advertise<nav_msgs::Path>("trajectory", 1, true);
  dwa_path_pub_ = nh_.advertise<nav_msgs::Path>("dwa_trajectory", 1, true);
  target_pub_ = nh_.advertise<visualization_msgs::MarkerArray>("target", 1, true);
}

template <class ModelT>
void Exploration<ModelT>::odomCallback(const nav_msgs::Odometry& msg)
{
  vb_(0) = msg.twist.twist.linear.x;
  vb_(1) = msg.twist.twist.linear.y;
  vb_(2) = msg.twist.twist.angular.z;
}

template <class ModelT>
void Exploration<ModelT>::mapCallback(const nav_msgs::OccupancyGrid::ConstPtr& msg)
{
  grid_.update(msg);
  // grid_.print();
  map_received_ = true;
}

template <class ModelT>
void Exploration<ModelT>::control(const Target& target, const std::string& map_frame_id,
                                  const std::string& base_frame_id, double frequency,
                                  double val_dt, double val_horizon)
{
  // Initialize the target distribution
  ergodic_control_.setTarget(target);

  visualization_msgs::MarkerArray marker_array;
  target.markers(marker_array, map_frame_id);
  target_pub_.publish(marker_array);

  const auto dwa_steps = dwa_.steps();
  mat opt_traj;
  vec u(3, arma::fill::zeros);
  vec pose(3, arma::fill::zeros);

  bool pose_known = false;
  bool follow_dwa = false;

  ros::Rate rate(frequency);
  unsigned int i = 0;

  while (nh_.ok())
  {
    try
    {
      const geometry_msgs::TransformStamped t_map_base =
          tfBuffer_.lookupTransform(map_frame_id, base_frame_id, ros::Time(0));
      pose(0) = t_map_base.transform.translation.x;
      pose(1) = t_map_base.transform.translation.y;
      pose(2) = normalize_angle_PI(tf2::getYaw(t_map_base.transform.rotation));

      // Add state to memory
      ergodic_control_.addStateMemory(pose);

      pose_known = true;
    }

    catch (tf2::TransformException& ex)
    {
      ROS_WARN_ONCE_NAMED(LOGNAME, "%s", ex.what());
      // continue;
    }

    // Contol loop
    if (map_received_ && pose_known)
    {
      if (follow_dwa)
      {
        i++;

        if (i == dwa_steps)
        {
          follow_dwa = false;
        }

        // else
        // {
        //   ROS_INFO_NAMED(LOGNAME, "Following DWA! %u ", i);
        //   // u.print("u");
        // }
      }

      else
      {
        u = ergodic_control_.control(grid_, pose);

        nav_msgs::Path trajectory;
        trajectory.header.frame_id = map_frame_id;
        ergodic_control_.path(trajectory);

        // ROS_INFO_STREAM_NAMED(LOGNAME, "Publish ergodic trajectory");
        opt_traj_pub_.publish(trajectory);
      }

      if (!validate_control(collision_, grid_, pose, u, val_dt, val_horizon))
      {
        ROS_INFO_STREAM_NAMED(LOGNAME, "Collision detected! Enabling DWA!");

        // Collision is a result from the previous dwa twist
        if (follow_dwa)
        {
          // ROS_WARN_STREAM_NAMED(LOGNAME, "Collision from DWA!");
          // If dwa fails twist is set to zeros
          dwa_.control(u, grid_, pose, vb_, u);

          // Whether or not dwa is successful the
          // ergodic controller will replan next iteration
          follow_dwa = false;
        }

        // Collision is a result from the ergodic controller
        else
        {
          // ROS_WARN_STREAM_NAMED(LOGNAME, "Collision from Ergodic Control!");
          ergodic_control_.optTraj(opt_traj);

          // If dwa fails twist is set to zeros
          if (!dwa_.control(u, grid_, pose, vb_, opt_traj, ergodic_control_.timeStep()))
          {
            // ergodic controller will replan next iteration
            // ROS_WARN_STREAM_NAMED(LOGNAME, "Trigger ergodic replan");
            follow_dwa = false;
          }

          else
          {
            // dwa found at a solution and the robot will follow it
            follow_dwa = true;
            i = 0;
            // ROS_INFO_NAMED(LOGNAME, "Following DWA! %u ", i);
          }
        }

        nav_msgs::Path dwa_traj;
        dwa_traj.header.frame_id = map_frame_id;
        constTwistPath(dwa_traj, pose, u, dwa_.timeStep(), dwa_.horizon());

        dwa_path_pub_.publish(dwa_traj);
      }  // end validate control

      geometry_msgs::Twist twist_msg;
      twist_msg.linear.x = u(0);
      twist_msg.linear.y = u(1);
      twist_msg.angular.z = u(2);

      cmd_pub_.publish(twist_msg);
    }  // end control loop

    rate.sleep();
  }  // end while loop
}

}  // namespace ergodic_exploration
