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

#include <ergodic_exploration/cart.hpp>
#include <ergodic_exploration/omni.hpp>
#include <ergodic_exploration/ergodic_control.hpp>
#include <ergodic_exploration/dynamic_window.hpp>
#include <ergodic_exploration/mapping.hpp>

using arma::eye;
using arma::mat;
using arma::span;
using arma::vec;

using namespace ergodic_exploration;


// static ergodic_exploration::GridMap grid;
// static bool received;

// void pdf(std::vector<int8_t>& map_data, unsigned int xsize, unsigned int ysize,
//          double resolution, double xmin, double ymin)
// {
//   double xm1 = 0.7;
//   double ym1 = 0.7;
//   double x1var = 0.1 * 0.1;
//   double y1var = 0.1 * 0.1;
//
//   double xm2 = 0.3;
//   double ym2 = 0.3;
//   double x2var = 0.1 * 0.1;
//   double y2var = 0.1 * 0.1;
//
//   unsigned int size = map_data.size();
//   for (unsigned int k = 0; k < size; k++)
//   {
//     const auto i = static_cast<unsigned int>(k / xsize);
//     const auto j = k - i * xsize;
//
//     const auto x = static_cast<double>(j * resolution) + resolution / 2.0 + xmin;
//     const auto y = static_cast<double>(i * resolution) + resolution / 2.0 + ymin;
//
//     double pdf1 = std::pow((x - xm1), 2) / x1var + std::pow((y - ym1), 2) / y1var;
//     double pdf2 = std::pow((x - xm2), 2) / x2var + std::pow((y - ym2), 2) / y2var;
//
//     map_data.at(k) = (std::exp(-0.5 * pdf1) + std::exp(-0.5 * pdf2)) * 500;
//     // std::cout << map_data.at(k) << std::endl;
//   }
// }

// template <typename modelT>
// bool validateControl(const modelT& model, const Collision& collision, const GridMap& grid,
//                      const vec x0, const vec& u, double dt, double horizon)
// {
//   vec x = x0;
//   RungeKutta rk4(dt);
//   const auto steps = static_cast<unsigned int>(horizon / std::abs(dt));
//
//   for (unsigned int i = 0; i < steps; i++)
//   {
//     x = rk4.step(model, x, u);
//     if (collision.collisionCheck(grid, x))
//     {
//       return false;
//     }
//   }
//
//   return true;
// }

// void mapCallback(const nav_msgs::OccupancyGrid::ConstPtr& msg)
// {
//   // grid.update(msg);
//   // grid.print();
//   // received = true;
// }

int main(int argc, char** argv)
{
  ROS_INFO_STREAM_NAMED("main", "Starting exploration");
  ros::init(argc, argv, "exploration");
  ros::NodeHandle nh;

  // ros::Subscriber map_sub = nh.subscribe("map", 1, mapCallback);

  ros::Publisher map_pub = nh.advertise<nav_msgs::OccupancyGrid>("map", 1, true);
  ros::Publisher pose_pub = nh.advertise<geometry_msgs::PoseStamped>("robot", 1, true);
  ros::Publisher path_pub = nh.advertise<nav_msgs::Path>("trajectory", 1, true);
  ros::Publisher target_pub = nh.advertise<visualization_msgs::MarkerArray>("target", 1, true);
  ros::Publisher sample_pub = nh.advertise<visualization_msgs::Marker>("samples", 1, true);


  arma::arma_rng::set_seed_random();

  // received = false;
  //////////////////////////////////////////////////////////////////////////////
  // TODO: Add noise to motion model
  // motion model
  Omni omni;
  //////////////////////////////////////////////////////////////////////////////
  const auto loop_freq = 20.0;
  const auto frequency = 10.0;

  // EC validation
  // const auto dt = 0.1;
  // const auto horizon = 0.5;
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
  Collision collision(boundary_radius, search_radius,
                        obstacle_threshold, occupied_threshold);
  //////////////////////////////////////////////////////////////////////////////
  const auto ec_dt = 0.1;
  const auto ec_horizon = 0.5;
  const auto num_samples = 1e3;
  const auto buffer_size = 1e6;
  const auto batch_size = 100;
  // const mat R = eye<mat>(3, 3);
  mat R(3, 3, arma::fill::zeros);
  R(0, 0) = 100.0;
  R(1, 1) = 100.0;
  R(2, 2) = 50.0;

  const mat Sigma = eye<mat>(2, 2) * 0.001;

  ErgodicControl ergodic_control(omni, ec_dt, ec_horizon, num_samples,
                                                      buffer_size, batch_size, R, Sigma);
  //////////////////////////////////////////////////////////////////////////////
  // dwa
  // double dwa_dt = 0.1;
  // double dwa_horizon = 1.0;
  // double dwa_frequency = 20.0;
  //
  // double acc_lim_x = 1.0;
  // double acc_lim_y = 1.0;
  // double acc_lim_th = 2.0;
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
  // unsigned int vy_samples = 10;
  // unsigned int vth_samples = 5;
  //
  // ergodic_exploration::DynamicWindow dwa(dwa_dt, dwa_horizon, dwa_frequency, acc_lim_x,
  //                                        acc_lim_y, acc_lim_th, max_vel_x, min_vel_x,
  //                                        max_vel_y, min_vel_y, max_rot_vel, min_rot_vel,
  //                                        vx_samples, vy_samples, vth_samples);
  //////////////////////////////////////////////////////////////////////////////

  RungeKutta rk4(1.0/loop_freq);

  Gaussian g1({ 0.7, 0.7 }, { 0.1, 0.1 });
  Gaussian g2({ 0.3, 0.3 }, { 0.1, 0.1 });
  Target target({ g1, g2 });

  visualization_msgs::MarkerArray marker_array;
  target.markers(marker_array, "map");

  map_pub.publish(grid_msg);
  target_pub.publish(marker_array);


  vec x = { 0.3, 0.4, 0.0 };
  vec u = { 0.0, 0.0, 0.0 };

  auto ctr = 0;
  ros::Rate rate(loop_freq);

  while (nh.ok())
  {
    // ros::spinOnce();

    if (ctr == static_cast<int>(loop_freq / frequency))
    {
      u = ergodic_control.control(collision, grid, target, x);
      // if (!validateControl(omni, collision, grid, x, u, 0.1, 0.5))
      // {
      //   std::cout << "DWA" << std::endl;
      //   u = dwa.control(collision, grid, x, vb, u);
      // }
      ctr = 0;
    }

    // u = ergodic_control.control(collision, grid, target, x);
    // u.print("u_t:");

    // if (!validateControl(omni, collision, grid, x, u, 0.1, 0.5))
    // {
    //   std::cout << "DWA" << std::endl;
    //   u = dwa.control(collision, grid, x, vb, u);
    // }

    // auto t_start = std::chrono::high_resolution_clock::now();
    // auto t_end = std::chrono::high_resolution_clock::now();
    // std::cout << "Hz: " << 1.0 / (std::chrono::duration<double, std::milli>(t_end -
    // t_start).count() / 1000.0) << std::endl;

    x = rk4.step(omni, x, u);

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

    nav_msgs::Path trajectory;
    ergodic_control.path(trajectory, "map");

    path_pub.publish(trajectory);

    visualization_msgs::Marker marker;
    ergodic_control.samples(marker, "map");

    sample_pub.publish(marker);

    ctr++;

    rate.sleep();
  }

  std::cout << " END " << std::endl;

  ros::waitForShutdown();
  ROS_INFO_STREAM_NAMED("main", "Shutting down.");
  return 0;
}
