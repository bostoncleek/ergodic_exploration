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
#include <geometry_msgs/Twist.h>
#include <nav_msgs/OccupancyGrid.h>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/LaserScan.h>
#include <tf2/LinearMath/Quaternion.h>
#include <tf2/LinearMath/Matrix3x3.h>
#include <tf2/utils.h>

#include <ergodic_exploration/models/cart.hpp>
#include <ergodic_exploration/models/omni.hpp>
#include <ergodic_exploration/ergodic_control.hpp>
#include <ergodic_exploration/dynamic_window.hpp>
#include <ergodic_exploration/mapping.hpp>

using arma::eye;
using arma::mat;
using arma::span;
using arma::vec;

using namespace ergodic_exploration;


int main(int argc, char** argv)
{
  ROS_INFO_STREAM_NAMED("main", "Starting exploration");
  ros::init(argc, argv, "exploration");
  ros::NodeHandle nh;

  ros::Publisher map_pub = nh.advertise<nav_msgs::OccupancyGrid>("map", 1, true);
  ros::Publisher pose_pub = nh.advertise<geometry_msgs::PoseStamped>("robot", 1, true);
  ros::Publisher path_pub = nh.advertise<nav_msgs::Path>("trajectory", 1, true);
  ros::Publisher target_pub =
      nh.advertise<visualization_msgs::MarkerArray>("target", 1, true);

  arma::arma_rng::set_seed_random();

  // received = false;
  //////////////////////////////////////////////////////////////////////////////
  // TODO: Add noise to motion model
  // motion model
  Omni omni;
  SimpleCart cart;
  //////////////////////////////////////////////////////////////////////////////
  const auto loop_freq = 20.0;
  const auto frequency = 10.0;

  const auto max_vel_x = 1.0;
  const auto min_vel_x = -1.0;

  const auto max_vel_y = 1.0;
  const auto min_vel_y = -1.0;

  const auto max_rot_vel = 2.0;
  const auto min_rot_vel = -2.0;

  const vec umin = { min_vel_x, min_vel_y, min_rot_vel };
  const vec umax = { max_vel_x, max_vel_y, max_rot_vel };
  //////////////////////////////////////////////////////////////////////////////
  // grid
  const auto xmin = 0.0;
  const auto xmax = 1.0;
  const auto ymin = 0.0;
  const auto ymax = 1.0;
  const auto resolution = 0.05;
  const auto xsize = ergodic_exploration::axis_length(xmin, xmax, resolution);
  const auto ysize = ergodic_exploration::axis_length(ymin, ymax, resolution);
  std::vector<int8_t> map_data(xsize * ysize, 50);

  GridMap grid(xmin, xmax, ymin, ymax, resolution, map_data);
  //////////////////////////////////////////////////////////////////////////////
  // ** Testing **
  nav_msgs::OccupancyGrid grid_msg;
  grid_msg.header.frame_id = "map";
  grid_msg.info.resolution = resolution;
  grid_msg.info.width = xsize;
  grid_msg.info.height = ysize;
  grid_msg.info.origin.position.x = xmin;
  grid_msg.info.origin.position.y = ymin;
  grid_msg.data = grid.gridData();
  //////////////////////////////////////////////////////////////////////////////
  // collision
  const auto boundary_radius = 0.1;
  const auto search_radius = 0.2;
  const auto obstacle_threshold = 0.05;
  const auto occupied_threshold = 0.9;
  Collision collision(boundary_radius, search_radius, obstacle_threshold,
                      occupied_threshold);
  //////////////////////////////////////////////////////////////////////////////
  const auto ec_dt = 0.1;
  const auto ec_horizon = 2.0;
  const auto target_resolution = 0.02;
  const auto expl_weight = 1.0;
  const auto num_basis = 5;
  const auto buffer_size = 1e6;
  const auto batch_size = 100;
  // const mat R = eye<mat>(3, 3);
  mat R(3, 3, arma::fill::zeros);
  R(0, 0) = 10.;
  R(1, 1) = 10.;
  R(2, 2) = 1.;

  ErgodicControl ergodic_control(omni, collision, ec_dt, ec_horizon, target_resolution,
                                 expl_weight, num_basis, buffer_size, batch_size, R, umin,
                                 umax);

  // double dwa_dt = 0.1;
  // double dwa_horizon = 1.0;
  // double acc_dt = 0.2;
  //
  // double acc_lim_x = 1.0;
  // double acc_lim_y = 0.0;
  // double acc_lim_th = 2.0;
  //
  // unsigned int vx_samples = 5;
  // unsigned int vy_samples = 0;
  // unsigned int vth_samples = 5;
  //
  // ergodic_exploration::DynamicWindow dwa(collision, dwa_dt, dwa_horizon, acc_dt, acc_lim_x,
  //                                       acc_lim_y, acc_lim_th, max_vel_x, min_vel_x,
  //                                       max_vel_y, min_vel_y, max_rot_vel, min_rot_vel,
  //                                       vx_samples, vy_samples, vth_samples);



  Gaussian g1({ 0.7, 0.7 }, { 0.1, 0.1 });
  Gaussian g2({ 0.3, 0.3 }, { 0.1, 0.1 });
  Target target({ g1, g2 });

  ergodic_control.setTarget(target);

  visualization_msgs::MarkerArray marker_array;
  target.markers(marker_array, "map");

  map_pub.publish(grid_msg);
  target_pub.publish(marker_array);

  vec x = { 0.2, 0.3, 0.0 };
  vec u = { 0.0, 0.0, 0.0 };

  // vec uref = { 0.7, 0.0, 0.0 };
  // vec vb = { 0.1, 0.0, 0.0 };
  //
  // u = dwa.control(grid, x, vb, uref);
  // u.print("u");

  RungeKutta rk4(0.1);
  ros::Rate rate(loop_freq);

  while (nh.ok())
  {
    ros::spinOnce();

    u = ergodic_control.control(grid, x);

    x = rk4.step(omni, x, u);

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

    pose_pub.publish(pose_bot);

    nav_msgs::Path trajectory;
    ergodic_control.path(trajectory, "map");

    path_pub.publish(trajectory);

    rate.sleep();
  }

  ROS_INFO_STREAM_NAMED("main", "Shutting down.");
  return 0;
}
