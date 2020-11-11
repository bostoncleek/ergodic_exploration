/**
 * @file exploration_node.cpp
 * @author Boston Cleek
 * @date 23 Oct 2020
 * @brief Ergodic exploration
 */

#include <iostream>
#include <ctime>
#include <chrono>

#include <ros/ros.h>
#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/OccupancyGrid.h>
#include <tf2/LinearMath/Quaternion.h>

#include <ergodic_exploration/cart.hpp>
#include <ergodic_exploration/omni.hpp>
#include <ergodic_exploration/ergodic_control.hpp>
#include <ergodic_exploration/dynamic_window.hpp>

using arma::eye;
using arma::mat;
using arma::vec;
using arma::span;
using ergodic_exploration::PI;

static ergodic_exploration::GridMap grid;
static bool received;

void mapCallback(const nav_msgs::OccupancyGrid::ConstPtr& msg)
{
  grid.update(msg);
  received = true;
}

int main(int argc, char** argv)
{
  ROS_INFO_STREAM_NAMED("main", "Starting exploration");
  ros::init(argc, argv, "exploration");
  ros::NodeHandle nh;

  ros::Subscriber map_sub = nh.subscribe("map", 1, mapCallback);
  ros::Publisher pose_pub =
      nh.advertise<geometry_msgs::PoseStamped>("robot", 1, true);
  ros::Publisher goal_pub =
      nh.advertise<geometry_msgs::PoseStamped>("goal", 1, true);


  received = false;



  // TODO: ADD noise to moition model when forward propagating

  // const auto wheel_radius = 0.2;
  // const auto wheel_base_x = 0.5;
  // const auto wheel_base_y = 0.5;
  // ergodic_exploration::Mecanum mecanum(wheel_radius, wheel_base_x, wheel_base_y);

  ergodic_exploration::Omni omni;
  // ergodic_exploration::SimpleCart cart;

  // const auto wheel_radius = 0.033;
  // const auto wheel_base = 0.08;
  // ergodic_exploration::Cart cart(wheel_radius, wheel_base);

  //////////////////////////////////////////////////////////////////////////////
  const auto dt = 0.1;
  const auto horizon = 0.5;
  const auto num_samples = 1000;
  const auto buffer_size = 1e6;
  const auto batch_size = 100;
  // const mat R = eye<mat>(4, 4) * 0.001;
  const mat R = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 0.01 } };
  const mat Sigma = eye<mat>(2, 2) * 0.0025;

  ergodic_exploration::ErgodicControl ergodic_control(omni, dt, horizon, num_samples,
                                                      buffer_size, batch_size, R, Sigma);
 //////////////////////////////////////////////////////////////////////////////

 //////////////////////////////////////////////////////////////////////////////
 // double boundary_radius = 0.05;
 // double search_radius = 0.1;
 // double obstacle_threshold = 0.01;
 // double occupied_threshold = 0.9;
 // ergodic_exploration::Collision collision(boundary_radius, search_radius,
 //                                          obstacle_threshold, occupied_threshold);
 //
 // double dt = 0.1;
 // double horizon = 0.5;
 //
 // double acc_lim_x = 1.0;
 // double acc_lim_y = 1.0;
 // double acc_lim_th = 2.0;
 //
 // double max_vel_trans = 1.0;
 // double min_vel_trans = -1.0;
 //
 // double max_vel_x = 1.0;
 // double min_vel_x = -1.0;
 //
 // double max_vel_y = 1.0;
 // double min_vel_y = -1.0;
 //
 // double max_rot_vel = 1.0;
 // double min_rot_vel = -1.0;
 //
 // unsigned int vx_samples = 5;
 // unsigned int vy_samples = 5;
 // unsigned int vth_samples = 10;
 //
 // ergodic_exploration::DynamicWindow dwa(omni, dt, horizon, acc_lim_x, acc_lim_y,
 //                                        acc_lim_th, max_vel_trans, min_vel_trans,
 //                                        max_vel_x, min_vel_x, max_vel_y, min_vel_y,
 //                                        max_rot_vel, min_rot_vel, vx_samples, vy_samples,
 //                                        vth_samples);
 //////////////////////////////////////////////////////////////////////////////

  ergodic_exploration::RungeKutta rk4(dt);

  // vec x = { 3.0, 3.0, 3.0/2.0 * PI};
  // vec x = { 3.0, 1.5, 0.0};
  vec x = { 2.5, 1.5, PI};

  vec vb = { 0.0, 0.0, 0.0 };

  vec goal = { 2.0, 1.5, PI };

  auto i = 0;
  while (nh.ok())
  {
    ros::spinOnce();

    if (received)
    {
      // auto c_start = std::clock();
      // auto t_start = std::chrono::high_resolution_clock::now();

      vec u = ergodic_control.control(grid, x);
      // vec u = dwa.control(collision, grid, x, goal, vb);
      // u.print("u(t):");

      x = rk4.step(omni, x, u);
      vb = u;

      // auto c_end = std::clock();
      // auto t_end = std::chrono::high_resolution_clock::now();
      //
      // std::cout
      //     << "Clock frequency (Hz): "
      //     << 1.0 / (static_cast<double>((c_end - c_start)) /
      //               static_cast<double>(CLOCKS_PER_SEC))
      //     << " Wall frequency (Hz): "
      //     << 1.0 / (std::chrono::duration<double, std::milli>(t_end - t_start).count() /
      //               1000.0)
      //     << std::endl;

      // x.print("pose");
      // std::cout << "Heading: " << x(2) << std::endl;

      //////////////////////////////////////////////////////////////////////////////
      tf2::Quaternion quat_bot;
      quat_bot.setRPY(0.0, 0.0, x(2));

      geometry_msgs::PoseStamped pose_bot;
      pose_bot.header.frame_id = "map";
      pose_bot.pose.position.x = x(0);
      pose_bot.pose.position.y = x(1);
      pose_bot.pose.position.z = 0.0;
      pose_bot.pose.orientation.x = quat_bot.x();
      pose_bot.pose.orientation.y = quat_bot.y();
      pose_bot.pose.orientation.z = quat_bot.z();
      pose_bot.pose.orientation.w = quat_bot.w();
      //////////////////////////////////////////////////////////////////////////////
      tf2::Quaternion quat_goal;
      quat_goal.setRPY(0.0, 0.0, goal(2));

      geometry_msgs::PoseStamped pose_goal;
      pose_goal.header.frame_id = "map";
      pose_goal.pose.position.x = goal(0);
      pose_goal.pose.position.y = goal(1);
      pose_goal.pose.position.z = 0.0;
      pose_goal.pose.orientation.x = quat_goal.x();
      pose_goal.pose.orientation.y = quat_goal.y();
      pose_goal.pose.orientation.z = quat_goal.z();
      pose_goal.pose.orientation.w = quat_goal.w();
      //////////////////////////////////////////////////////////////////////////////


      pose_pub.publish(pose_bot);
      goal_pub.publish(pose_goal);

      ros::Duration(0.1).sleep();

      if (i == 5000)
      {
        break;
      }

      i++;
    }
  }

  ros::waitForShutdown();
  ROS_INFO_STREAM_NAMED("main", "Shutting down.");
  return 0;
}
