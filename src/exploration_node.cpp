/**
 * @file exploration_node.cpp
 * @author Boston Cleek
 * @date 23 Oct 2020
 * @brief Ergodic exploration
 */

#include <iostream>
#include <typeinfo>

#include <ros/ros.h>
#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/OccupancyGrid.h>
#include <tf2/LinearMath/Quaternion.h>

#include <ergodic_exploration/collision.hpp>
#include <ergodic_exploration/grid.hpp>
#include <ergodic_exploration/integrator.hpp>
#include <ergodic_exploration/cart.hpp>
#include <ergodic_exploration/omni.hpp>
#include <ergodic_exploration/types.hpp>
#include <ergodic_exploration/numerics.hpp>
#include <ergodic_exploration/ergodic_control.hpp>

using arma::eye;
using arma::mat;
using arma::vec;

void pdf(std::vector<int8_t>& map_data, unsigned int xsize, unsigned int ysize,
         double resolution, double xmin, double ymin)
{
  double xm1 = 0.7;
  double ym1 = 0.7;
  double x1var = 0.2 * 0.2;
  double y1var = 0.2 * 0.2;

  double xm2 = 0.3;
  double ym2 = 0.3;
  double x2var = 0.3 * 0.3;
  double y2var = 0.3 * 0.3;

  unsigned int size = map_data.size();
  for (unsigned int k = 0; k < size; k++)
  {
    const auto i = static_cast<unsigned int>(k / xsize);
    const auto j = k - i * xsize;

    const auto x = static_cast<double>(j * resolution) + resolution / 2.0 + xmin;
    const auto y = static_cast<double>(i * resolution) + resolution / 2.0 + ymin;

    double pdf1 = std::pow((x - xm1), 2) / x1var + std::pow((y - ym1), 2) / y1var;
    double pdf2 = std::pow((x - xm2), 2) / x2var + std::pow((y - ym2), 2) / y2var;

    map_data.at(k) = (std::exp(-0.5 * pdf1) + std::exp(-0.5 * pdf2)) * 100;
    // std::cout << map_data.at(k) << std::endl;
  }
}

int main(int argc, char** argv)
{
  ROS_INFO_STREAM_NAMED("main", "Starting exploration");
  ros::init(argc, argv, "exploration");
  ros::NodeHandle nh;

  ros::Publisher map_pub = nh.advertise<nav_msgs::OccupancyGrid>("map", 1, true);
  ros::Publisher pose_pub =
      nh.advertise<geometry_msgs::PoseStamped>("robot_pose", 1, true);

  const auto xmin = 0.0;
  const auto xmax = 1.0;
  const auto ymin = 0.0;
  const auto ymax = 1.0;
  const auto resolution = 0.02;
  const auto xsize = ergodic_exploration::axis_length(xmin, xmax, resolution);
  const auto ysize = ergodic_exploration::axis_length(ymin, ymax, resolution);
  std::vector<int8_t> map_data(xsize * ysize, 50);

  pdf(map_data, xsize, ysize, resolution, xmin, ymin);

  nav_msgs::OccupancyGrid map;
  map.header.frame_id = "map";
  map.info.resolution = resolution;
  map.info.origin.position.x = xmin;
  map.info.origin.position.y = ymin;
  map.info.width = xsize;
  map.info.height = ysize;
  map.data = map_data;

  ergodic_exploration::GridMap grid_map(xmin, xmax, ymin, ymax, resolution, map_data);

  // TODO: ADD noise to moition model when forward propagating

  const auto wheel_radius = 0.2;
  const auto wheel_base_x = 0.5;
  const auto wheel_base_y = 0.5;
  ergodic_exploration::Mecanum mecanum(wheel_radius, wheel_base_x, wheel_base_y);

  // ergodic_exploration::Omni omni;

  // const auto wheel_radius = 0.033;
  // const auto wheel_base = 0.08;
  // ergodic_exploration::Cart cart(wheel_radius, wheel_base);

  const auto dt = 0.1;
  const auto horizon = 0.5;
  const auto num_samples = 1000;
  const auto buffer_size = 1e6;
  const auto batch_size = 100;
  const mat R = eye<mat>(4, 4) * 0.001;

  // const mat R = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };

  const mat Sigma = eye<mat>(2, 2) * 0.0025;

  ergodic_exploration::ErgodicControl ergodic_control(mecanum, dt, horizon, num_samples,
                                                      buffer_size, batch_size, R, Sigma);

  ergodic_exploration::RungeKutta rk4(dt);

  vec x = { 0.5, 0.5, 0.0 };
  auto i = 0;

  while (nh.ok())
  {
    vec u = ergodic_control.control(grid_map, x);
    x = rk4.step(mecanum, x, u);

    tf2::Quaternion quat;
    quat.setRPY(0.0, 0.0, x(2));

    // std::cout << "Heading: " << x(2) << std::endl;

    geometry_msgs::PoseStamped pose;
    pose.header.frame_id = "map";
    pose.pose.position.x = x(0);
    pose.pose.position.y = x(1);
    pose.pose.position.z = 0.0;
    pose.pose.orientation.x = quat.x();
    pose.pose.orientation.y = quat.y();
    pose.pose.orientation.z = quat.z();
    pose.pose.orientation.w = quat.w();

    map_pub.publish(map);
    pose_pub.publish(pose);

    ros::Duration(0.01).sleep();

    if (i == 5000)
    {
      break;
    }

    i++;
  }

  ros::waitForShutdown();
  ROS_INFO_STREAM_NAMED("main", "Shutting down.");
  return 0;
}
