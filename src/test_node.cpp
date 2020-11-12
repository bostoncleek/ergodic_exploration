/**
 * @file test_node.cpp
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
using arma::span;
using arma::vec;
using ergodic_exploration::Collision;
using ergodic_exploration::GridMap;
using ergodic_exploration::PI;
using ergodic_exploration::RungeKutta;

static ergodic_exploration::GridMap grid;
static bool received;

template <typename modelT>
bool validateControl(const modelT& model, const Collision& collision, const GridMap& grid,
                     const vec x0, const vec& u, double dt, double horizon)
{
  vec x = x0;
  RungeKutta rk4(dt);
  const auto steps = static_cast<unsigned int>(horizon / std::abs(dt));

  for (unsigned int i = 0; i < steps; i++)
  {
    x = rk4.step(model, x, u);
    if (collision.collisionCheck(grid, x))
    {
      return false;
    }
  }

  return true;
}

void mapCallback(const nav_msgs::OccupancyGrid::ConstPtr& msg)
{
  grid.update(msg);
  grid.print();
  received = true;
}

int main(int argc, char** argv)
{
  ROS_INFO_STREAM_NAMED("main", "Starting exploration");
  ros::init(argc, argv, "exploration");
  ros::NodeHandle nh;

  ros::Subscriber map_sub = nh.subscribe("map", 1, mapCallback);
  ros::Publisher pose_pub = nh.advertise<geometry_msgs::PoseStamped>("robot", 1, true);
  ros::Publisher goal_pub = nh.advertise<geometry_msgs::PoseStamped>("goal", 1, true);

  received = false;

  // TODO: ADD noise to moition model when forward propagating
  ergodic_exploration::Omni omni;

  const auto boundary_radius = 0.1;
  const auto search_radius = 0.2;
  const auto obstacle_threshold = 0.05;
  const auto occupied_threshold = 0.9;
  ergodic_exploration::Collision collision(boundary_radius, search_radius,
                                           obstacle_threshold, occupied_threshold);
  //////////////////////////////////////////////////////////////////////////////
  const auto dt = 0.05;
  const auto horizon = 3.0;
  const auto num_samples = 1e3;
  const auto buffer_size = 1e6;
  const auto batch_size = 100;
  // const mat R = eye<mat>(3, 3);
  mat R(3, 3, arma::fill::zeros);
  R(0, 0) = 1.0;
  R(1, 1) = 1.0;
  R(2, 2) = 0.1;

  const mat Sigma = eye<mat>(2, 2) * 0.0025;

  ergodic_exploration::ErgodicControl ergodic_control(omni, dt, horizon, num_samples,
                                                      buffer_size, batch_size, R, Sigma);
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  double dwa_dt = 0.1;
  double dwa_horizon = 0.5;
  double dwa_frequency = 20.0;

  double acc_lim_x = 1.0;
  double acc_lim_y = 1.0;
  double acc_lim_th = 2.0;

  double max_vel_trans = 1.0;
  double min_vel_trans = -1.0;

  double max_vel_x = 1.0;
  double min_vel_x = -1.0;

  double max_vel_y = 1.0;
  double min_vel_y = -1.0;

  double max_rot_vel = 1.0;
  double min_rot_vel = -1.0;

  unsigned int vx_samples = 5;
  unsigned int vy_samples = 10;
  unsigned int vth_samples = 5;

  ergodic_exploration::DynamicWindow dwa(omni, dwa_dt, dwa_horizon, dwa_frequency,
                                         acc_lim_x, acc_lim_y, acc_lim_th, max_vel_trans,
                                         min_vel_trans, max_vel_x, min_vel_x, max_vel_y,
                                         min_vel_y, max_rot_vel, min_rot_vel, vx_samples,
                                         vy_samples, vth_samples);
  //////////////////////////////////////////////////////////////////////////////

  // vec x = { 3.0, 3.0, 3.0/2.0 * PI};
  // vec x = { 3.0, 1.5, 0.0};
  vec x = { 2.5, 1.5, PI };

  vec vb = { 0.0, 0.0, 0.0 };

  vec u = { 0.0, 0.0, 0.0 };

  // vec vref = { 0.5, 0.0, 0.2 };

  auto loop_freq = 50.0;
  auto control_freq = 20.0;
  ros::Rate rate(loop_freq);

  ergodic_exploration::RungeKutta rk4(dt);

  auto i = 0;
  auto ctr = 0;
  while (nh.ok())
  {
    ros::spinOnce();

    if (received)
    {
      if (ctr == static_cast<int>(loop_freq / control_freq))
      {
        u = ergodic_control.control(collision, grid, x);
        if (!validateControl(omni, collision, grid, x, u, 0.1, 0.5))
        {
          std::cout << "DWA" << std::endl;
          u = dwa.control(collision, grid, x, vb, u);
        }
        ctr = 0;
      }

      // u = ergodic_control.control(collision, grid, x);
      //
      // if (!validateControl(omni, collision, grid, x, u, 0.1, 0.5))
      // {
      //   std::cout << "DWA" << std::endl;
      //   u = dwa.control(collision, grid, x, vb, u);
      // }

      // auto t_start = std::chrono::high_resolution_clock::now();
      //
      // // u = dwa.control(collision, grid, x, vb, vref);
      // // u = ergodic_control.control(collision, grid, x);
      //
      //
      // auto t_end = std::chrono::high_resolution_clock::now();
      // std::cout << "Hz: " << 1.0 / (std::chrono::duration<double, std::milli>(t_end -
      // t_start).count() / 1000.0) << std::endl;

      x = rk4.step(omni, x, u);
      vb = u;

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

      pose_pub.publish(pose_bot);

      if (i == 5000)
      {
        break;
      }

      i++;
      ctr++;
    }

    // rate.sleep();
  }

  std::cout << " END " << std::endl;

  ros::waitForShutdown();
  ROS_INFO_STREAM_NAMED("main", "Shutting down.");
  return 0;
}
