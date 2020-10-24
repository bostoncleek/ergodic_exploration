/**
 * @file mutual_information_node.cpp
 * @author Boston Cleek
 * @date 23 Oct 2020
 * @brief Publish mututal information given 2D occupancy grid
 */

#include <iostream>

#include <ros/ros.h>
#include <nav_msgs/OccupancyGrid.h>

#include <ergodic_exploration/grid.hpp>

int main(int argc, char** argv)
{
  ROS_INFO_STREAM_NAMED("main", "Starting mututal_information");
  ros::init(argc, argv, "mututal_information");
  ros::NodeHandle nh;

  ros::AsyncSpinner spinner(1);
  // spinner.start();

  double resolution = 0.5;
  unsigned int xsize = 4;
  unsigned int ysize = 3;
  std::vector<int8_t> map_data(xsize * ysize, 0.0);

  map_data.at(0) = 100;
  map_data.at(1) = 100;
  map_data.at(2) = 100;
  map_data.at(3) = 100;

  map_data.at(8) = 100;
  map_data.at(9) = 100;
  map_data.at(10) = 100;
  map_data.at(11) = 100;

  nav_msgs::OccupancyGrid map;
  map.header.frame_id = "map";
  map.info.resolution = resolution;
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
