/**
 * @file mutual_information_node.cpp
 * @author Boston Cleek
 * @date 23 Oct 2020
 * @brief Publish mututal information given 2D occupancy grid
 */

#include <iostream>
#include <typeinfo>

#include <ros/ros.h>
#include <nav_msgs/OccupancyGrid.h>

#include <ergodic_exploration/collision.hpp>
#include <ergodic_exploration/grid.hpp>
#include <ergodic_exploration/integrator.hpp>
#include <ergodic_exploration/cart.hpp>
#include <ergodic_exploration/omni.hpp>

// #include <armadillo>

int main(int argc, char** argv)
{
  ROS_INFO_STREAM_NAMED("main", "Starting exploration");
  ros::init(argc, argv, "exploration");
  ros::NodeHandle nh;

  ros::AsyncSpinner spinner(1);
  // spinner.start();

  const auto xmin = 0.0;
  const auto xmax = 10.0;
  const auto ymin = 0.0;
  const auto ymax = 10.0;
  const double resolution = 0.05;
  const unsigned int xsize = ergodic_exploration::axis_length(xmin, xmax, resolution);
  const unsigned int ysize = ergodic_exploration::axis_length(ymin, ymax, resolution);
  std::vector<int8_t> map_data(xsize * ysize, 100);

  ergodic_exploration::GridMap grid_map(xmin, xmax, ymin, ymax, resolution, map_data);

  double boundary_radius = 0.5;
  double search_radius = 1.0;
  double obstacle_threshold = 0.5;
  int8_t occupied_threshold = 90;
  ergodic_exploration::Collision collision_detector(
      boundary_radius, search_radius, obstacle_threshold, occupied_threshold);

  // std::cout << "xsize: " << xsize << std::endl;
  // std::cout << "ysize: " << ysize << std::endl;
  // std::cout << "size of map: " << map_data.size() << std::endl;

  // ergodic_exploration::Cart cart(0.1, 2.0);

  const auto wheel_radius = 0.1;
  const auto wheel_base_x = 0.5;
  const auto wheel_base_y = 0.5;
  ergodic_exploration::Mecanum mecanum(wheel_radius, wheel_base_x, wheel_base_y);

  const auto horizon = 0.4;
  const auto dt = 0.1;
  const auto N = static_cast<int>(horizon / dt);

  arma::mat ut(4, N);
  ut.row(0).fill(0.0);
  ut.row(1).fill(0.0);
  ut.row(2).fill(0.0);
  ut.row(3).fill(0.0);

  arma::vec x0 = { 0.0, 0.0, 0.0 };

  // arma::vec xt = mecanum(x0, ut.col(0));

  ergodic_exploration::RungeKutta rk4(0.1);

  // arma::vec x1 = rk4.step(mecanum, x0, ut.col(0));
  // x1.print("x1:");

  arma::mat xt = rk4.solve(mecanum, x0, ut, horizon);
  xt.print("xt:");

  collision_detector.obstacleCells(grid_map, 100, 100);
  ergodic_exploration::CollisionMap collision_map = collision_detector.collisionMap();

  for (const auto item : collision_map)
  {
    map_data.at(item.first) = 0;
  }

  nav_msgs::OccupancyGrid map;
  map.header.frame_id = "map";
  map.info.resolution = resolution;
  map.info.origin.position.x = xmin;
  map.info.origin.position.y = ymin;
  map.info.width = xsize;
  map.info.height = ysize;
  map.data = map_data;

  ros::Publisher map_pub = nh.advertise<nav_msgs::OccupancyGrid>("map", 1);

  while (nh.ok())
  {
    ros::spinOnce();
    map_pub.publish(map);
  }

  // ros::waitForShutdown();
  ROS_INFO_STREAM_NAMED("main", "Shutting down.");
  return 0;
}
