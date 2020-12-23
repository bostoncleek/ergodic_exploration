/*********************************************************************
 * BSD 3-Clause License
 *
 * Copyright (c) 2020 Northwestern University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 *  * Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/
/**
 * @file exploration.hpp
 * @author Boston Cleek
 * @date 4 Dec 2020
 * @brief Ergodic exploration with a user defined target and dynamic window
 * approach for collision avoidance
 */

#include <ros/ros.h>
#include <ros/console.h>
#include <geometry_msgs/Twist.h>
#include <geometry_msgs/TransformStamped.h>
#include <nav_msgs/OccupancyGrid.h>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/LaserScan.h>
#include <tf2_ros/transform_listener.h>
#include <tf2/LinearMath/Quaternion.h>

#include <ergodic_exploration/ergodic_control.hpp>
#include <ergodic_exploration/dynamic_window.hpp>

#include <ergodic_exploration/models/omni.hpp>

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
  Exploration(ros::NodeHandle& nh, const ErgodicControl<ModelT>& ergodic_control,
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
                                 const ErgodicControl<ModelT>& ergodic_control,
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
  opt_traj_pub_ = nh_.advertise<nav_msgs::Path>("trajectory", 1);
  dwa_path_pub_ = nh_.advertise<nav_msgs::Path>("dwa_trajectory", 1);
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

  const visualization_msgs::MarkerArray marker_array = target.markers(map_frame_id);
  target_pub_.publish(marker_array);

  const auto dwa_steps = dwa_.steps();
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
      pose(2) = getYaw(t_map_base.transform.rotation.x, t_map_base.transform.rotation.y,
                       t_map_base.transform.rotation.z, t_map_base.transform.rotation.w);

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
        // May need to replan
        follow_dwa = (i != dwa_steps) ? true : false;
      }

      if (!follow_dwa)
      {
        u = ergodic_control_.control(grid_, pose);

        const nav_msgs::Path trajectory = ergodic_control_.path(map_frame_id);
        opt_traj_pub_.publish(trajectory);
      }

      if (!validate_control(collision_, grid_, pose, u, val_dt, val_horizon))
      {
        // ROS_INFO_STREAM_NAMED(LOGNAME, "Collision detected! Enabling DWA!");

        // Collision is a result from the previous dwa twist
        if (follow_dwa)
        {
          // If dwa fails twist is set to zeros
          u = std::get<vec>(dwa_.control(grid_, pose, vb_, u));

          // Whether or not dwa is successful the
          // ergodic controller will replan next iteration
          follow_dwa = false;
        }

        // Collision is a result from the ergodic controller
        else
        {
          // Reference trajectory for dwa
          const mat opt_traj = ergodic_control_.optTraj();

          // If dwa fails twist is set to zeros
          const tuple<bool, vec> dwa_state =
              dwa_.control(grid_, pose, vb_, opt_traj, ergodic_control_.timeStep());

          u = std::get<vec>(dwa_state);

          // if dwa found at a solution and the robot will follow it
          // else ergodic controller will replan next iteration
          follow_dwa = std::get<bool>(dwa_state);

          if (follow_dwa)
          {
            i = 0;

            // Visualize dwa predicted trajectory
            const nav_msgs::Path dwa_traj =
                constTwistPath(map_frame_id, pose, u, dwa_.timeStep(), dwa_.horizon());
            dwa_path_pub_.publish(dwa_traj);
          }
        }
      }  // end validate control

      geometry_msgs::Twist twist_msg;
      twist_msg.linear.x = u(0);
      twist_msg.linear.y = u(1);
      twist_msg.angular.z = u(2);

      cmd_pub_.publish(twist_msg);
      // u.print();

    }  // end control loop

    rate.sleep();
  }  // end while loop
}
}  // namespace ergodic_exploration
