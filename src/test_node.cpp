

#include <iostream>
#include <typeinfo>
#include <ctime>
#include <chrono>
#include <vector>
#include <memory>

#include <ros/ros.h>
#include <nav_msgs/OccupancyGrid.h>

#include <ergodic_exploration/collision.hpp>
#include <ergodic_exploration/grid.hpp>
#include <ergodic_exploration/integrator.hpp>
#include <ergodic_exploration/cart.hpp>
#include <ergodic_exploration/omni.hpp>
#include <ergodic_exploration/types.hpp>
#include <ergodic_exploration/numerics.hpp>
#include <ergodic_exploration/ergodic_control.hpp>
#include <ergodic_exploration/mapping.hpp>
#include <ergodic_exploration/dynamic_window.hpp>

// #include <armadillo>

using arma::eye;
using arma::mat;
using arma::span;
using arma::vec;
using ergodic_exploration::PI;

int main(int argc, char** argv)
{
  ROS_INFO_STREAM_NAMED("main", "Starting test");
  ros::init(argc, argv, "test");
  ros::NodeHandle nh;

  // ros::AsyncSpinner spinner(1);
  // spinner.start();

  const auto xmin = 0.0;
  const auto xmax = 2.0;
  const auto ymin = 0.0;
  const auto ymax = 3.0;
  const double resolution = 0.5;
  const unsigned int xsize = ergodic_exploration::axis_length(xmin, xmax, resolution);
  const unsigned int ysize = ergodic_exploration::axis_length(ymin, ymax, resolution);
  std::vector<int8_t> map_data(xsize * ysize, 50);

  map_data.at(13) = 100;

  ergodic_exploration::GridMap grid(xmin, xmax, ymin, ymax, resolution, map_data);

  double boundary_radius = 1.0;
  double search_radius = 3.0;
  double obstacle_threshold = 0.5;
  double occupied_threshold = 0.9;
  ergodic_exploration::Collision collision(boundary_radius, search_radius,
                                           obstacle_threshold, occupied_threshold);

  ergodic_exploration::Omni omni;

  double dt = 0.1;
  double horizon = 0.5;
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
  unsigned int vy_samples = 5;
  unsigned int vth_samples = 10;

  ergodic_exploration::DynamicWindow dwa(omni, dt, horizon, acc_lim_x, acc_lim_y,
                                         acc_lim_th, max_vel_trans, min_vel_trans,
                                         max_vel_x, min_vel_x, max_vel_y, min_vel_y,
                                         max_rot_vel, min_rot_vel, vx_samples, vy_samples,
                                         vth_samples);

  vec x = { 0.0, 0.0, PI / 2.0 };
  vec vb = { 0.0, 0.0, 0.0 };
  vec goal = { 7.0, 5.0, PI / 2.0 };

  //////////////////////////////////////////////////////////////////////////////
  clock_t c_start = clock();
  auto t_start = std::chrono::high_resolution_clock::now();

  // dwa.control(collision, grid, x, goal, vb);

  double dmin = 0.0;
  if (collision.minDistance(dmin, grid, x) == ergodic_exploration::CollisionMsg::crash)
  {
    std::cout << "collision" << std::endl;
  }
  std::cout << "dis: " << dmin << std::endl;

  vec disp;
  if (collision.minDirection(disp, grid, x) == ergodic_exploration::CollisionMsg::obstacle)
  {
    disp.print();
  }

  clock_t c_end = clock();
  auto t_end = std::chrono::high_resolution_clock::now();

  std::cout << "Clock frequency (Hz): "
            << 1.0 / (static_cast<double>((c_end - c_start)) /
                      static_cast<double>(CLOCKS_PER_SEC))
            << " Wall frequency (Hz): "
            << 1.0 / (std::chrono::duration<double, std::milli>(t_end - t_start).count() /
                      1000.0)
            << std::endl;
  //////////////////////////////////////////////////////////////////////////////

  // // std::cout <<  << std::endl;

  // std::cout << "xsize: " << xsize << std::endl;
  // std::cout << "ysize: " << ysize << std::endl;
  // std::cout << "size of map: " << map_data.size() << std::endl;

  // ergodic_exploration::Cart cart(0.1, 2.0);

  // const auto wheel_radius = 0.2;
  // const auto wheel_base_x = 0.5;
  // const auto wheel_base_y = 0.5;
  // ergodic_exploration::Mecanum mecanum(wheel_radius, wheel_base_x, wheel_base_y);

  // ergodic_exploration::Omni omni;
  //
  // // const auto wheel_radius = 0.033;
  // // const auto wheel_base = 0.08;
  // // ergodic_exploration::Cart cart(wheel_radius, wheel_base);
  //
  // const auto dt = 0.1;
  // const auto horizon = 0.5;
  // const auto num_samples = 100;
  // const auto buffer_size = 1e6;
  // const auto batch_size = 1000;
  // // const mat R = eye<mat>(4, 4) * 0.001;
  // // const mat Sigma = eye<mat>(2, 2) * 0.0025;
  //
  // const mat R = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
  // const mat Sigma = eye<mat>(2, 2) * 0.0025;
  //
  // ergodic_exploration::ErgodicControl ergodic_control(omni, dt, horizon, num_samples,
  //                                                     buffer_size, batch_size, R, Sigma);
  //
  // // ergodic_control.sample(grid_map);
  // const vec x0 = { 0.7, 0.5, 0.0 };
  // ergodic_control.control(grid_map, x0);

  // ergodic_exploration::RungeKutta rk4(dt);
  //
  // vec x = x0;
  // for (unsigned int i = 0; i < 100; i++)
  // {
  //   // clock_t c_start = clock();
  //   // auto t_start = std::chrono::high_resolution_clock::now();
  //   //
  //   // vec u = ergodic_control.control(grid_map, x);
  //   //
  //   // clock_t c_end = clock();
  //   // auto t_end = std::chrono::high_resolution_clock::now();
  //   //
  //   // std::cout
  //   //     << "Clock frequency (Hz): "
  //   //     << 1.0 / (static_cast<double>((c_end - c_start)) /
  //   //               static_cast<double>(CLOCKS_PER_SEC))
  //   //     << " Wall frequency (Hz): "
  //   //     << 1.0 / (std::chrono::duration<double, std::milli>(t_end - t_start).count() /
  //   //               1000.0)
  //   //     << std::endl;
  //
  //   vec u = ergodic_control.control(grid_map, x);
  //   x = rk4.step(omni, x, u);
  //   x = rk4.step(mecanum, x, u);
  //   // x.print();
  // }
  // x.print();

  // vec rho = {1.0, 0.2, 0.5};
  // vec kldx = {0.1, 0.2, 0.0};
  // vec u = {1.0, 0.2, 0.5, 1.0};
  // mat fdx = mecanum.fdx(x0, u);
  //
  // vec drho = ergodic_exploration::rhodot(rho, kldx, fdx);

  // std::cout << "entropy: " << ergodic_exploration::entropy(0.0) << std::endl;
  // std::cout << "entropy: " << ergodic_exploration::entropy(-1.0/100) << std::endl;
  // std::cout << "entropy: " << ergodic_exploration::entropy(0.5) << std::endl;

  // unsigned int num_pts = 1080;
  // float angle_min = 0.0;
  // float angle_max = 2.0 * PI;
  // float angle_increment = (angle_max - angle_min) / static_cast<float>(num_pts);
  // std::vector<float> ranges(num_pts, 2.0);
  //
  // // float angle = angle_min;
  // // for (unsigned int i = 0; i < num_pts; i++)
  // // {
  // //   ranges.at(i) += angle_min;
  // //   angle_min += 0.1;
  // // }
  //
  // sensor_msgs::LaserScan scan_msg;
  // scan_msg.angle_min = angle_min;
  // scan_msg.angle_max = angle_max;
  // scan_msg.angle_increment = angle_increment;
  // scan_msg.ranges = ranges;
  //
  // auto scan_ptr = boost::make_shared<const sensor_msgs::LaserScan>(scan_msg);
  //
  // // vec x = { 5.0, 5.0, PI / 2.0 };
  // vec x = { 5.0, 5.0, 0.0 };
  // mat Tbs = ergodic_exploration::transform2d(-PI / 2.0);
  // // Tbs.print();
  //
  // ergodic_exploration::OccupancyMapper mapper(Tbs);
  //
  // auto c_start = std::clock();
  // auto t_start = std::chrono::high_resolution_clock::now();
  //
  // mapper.updateMap(grid_map, scan_ptr, x);
  // // mapper.updateMap(grid_map, scan_ptr, x);
  // // mapper.updateMap(grid_map, scan_ptr, x);
  // // mapper.updateMap(grid_map, scan_ptr, x);
  //
  // auto c_end = std::clock();
  // auto t_end = std::chrono::high_resolution_clock::now();
  //
  // std::cout << "Clock frequency (Hz): "
  //           << 1.0 / (static_cast<double>((c_end - c_start)) /
  //                     static_cast<double>(CLOCKS_PER_SEC))
  //           << " Wall frequency (Hz): "
  //           << 1.0 / (std::chrono::duration<double, std::milli>(t_end - t_start).count() /
  //                     1000.0)
  //           << std::endl;
  //
  nav_msgs::OccupancyGrid map;
  map.header.frame_id = "map";
  map.info.resolution = resolution;
  map.info.origin.position.x = xmin;
  map.info.origin.position.y = ymin;
  map.info.width = xsize;
  map.info.height = ysize;
  // map.data = map_data;
  map.data = grid.gridData();

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
