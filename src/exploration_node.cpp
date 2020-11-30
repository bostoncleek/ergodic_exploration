/**
 * @file exploration_node.cpp
 * @author Boston Cleek
 * @date 23 Oct 2020
 * @brief Ergodic exploration
 */

#include <iostream>
#include <variant>
#include <utility>
#include <ctime>
#include <chrono>
#include <memory>

#include <ros/ros.h>
#include <ros/console.h>
#include <geometry_msgs/Twist.h>
#include <geometry_msgs/TransformStamped.h>
#include <nav_msgs/OccupancyGrid.h>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/LaserScan.h>
#include <tf2_ros/transform_listener.h>
#include <tf2/LinearMath/Quaternion.h>
#include <tf2/LinearMath/Matrix3x3.h>
#include <tf2/utils.h>
#include <gazebo_msgs/ModelStates.h>

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

constexpr char LOGNAME[] = "ergodic exploration";

static GridMap grid;
static vec pose = { 0.0, 0.0, 0.0 };
static vec vb = { 0.0, 0.0, 0.0 };

static bool odom_update = false;
static bool map_received = false;

void odomCallback(const nav_msgs::Odometry& msg)
{
  vb(0) = msg.twist.twist.linear.x;
  vb(1) = msg.twist.twist.linear.y;
  vb(2) = msg.twist.twist.angular.z;
}

// void modelCallBack(const gazebo_msgs::ModelStates& msg)
// {
//   // store names of all items in gazebo
//   std::vector<std::string> names = msg.name;
//
//   // index of robot
//   int robot_index = 0;
//
//   // find diff_drive robot
//   int ctr = 0;
//   for (const auto& item : names)
//   {
//     // check for robot
//     if (item == "nuridgeback")
//     {
//       robot_index = ctr;
//     }
//
//     ctr++;
//   }
//
//   // pose(0) = msg.pose[robot_index].position.x;
//   // pose(1) = msg.pose[robot_index].position.y;
//   // pose(2) = tf2::getYaw(msg.pose[robot_index].orientation);
//
//   // std::cout << "Pose gazebo: " << msg.pose[robot_index].position.x <<
//   // " " << msg.pose[robot_index].position.y << " " <<
//   // tf2::getYaw(msg.pose[robot_index].orientation) << std::endl;
// }

void mapCallback(const nav_msgs::OccupancyGrid::ConstPtr& msg)
{
  grid.update(msg);
  // grid.print();
  map_received = true;
}

int main(int argc, char** argv)
{
  ROS_INFO_STREAM_NAMED(LOGNAME, "Starting exploration");
  ros::init(argc, argv, "exploration");
  ros::NodeHandle nh;
  ros::NodeHandle pnh("~");

  ros::Subscriber map_sub = nh.subscribe("map", 1, mapCallback);
  ros::Subscriber odom_sub = nh.subscribe("odom", 1, odomCallback);
  // ros::Subscriber model_sub = nh.subscribe("/gazebo/model_states", 1, modelCallBack);

  // ros::Publisher map_pub = nh.advertise<nav_msgs::OccupancyGrid>("map_update", 1, true);
  ros::Publisher cmd_pub = nh.advertise<geometry_msgs::Twist>("cmd_vel", 1);
  ros::Publisher path_pub = nh.advertise<nav_msgs::Path>("trajectory", 1, true);
  ros::Publisher target_pub =
      nh.advertise<visualization_msgs::MarkerArray>("target", 1, true);

  tf2_ros::Buffer tfBuffer;
  tf2_ros::TransformListener tfListener(tfBuffer);
  geometry_msgs::TransformStamped t_map_base;

  arma::arma_rng::set_seed_random();

  //////////////////////////////////////////////////////////////////////////////
  // TODO: Add noise to motion model
  // motion model
  Omni omni;
  SimpleCart cart;

  // std::variant<Omni, SimpleCart> motion_model;
  // bool holonomic = true;
  // // if(!pnh.getParam("holonomic", holonomic))
  // // {
  // //   ROS_ERROR_STREAM_NAMED(LOGNAME, "Vechilce type not specified");
  // //   ros::shutdown();
  // // }
  //
  // if (!holonomic)
  // {
  //   Omni omni;
  //   motion_model = omni;
  // }
  //
  // else
  // {
  //   SimpleCart cart;
  //   motion_model = cart;
  // }

  const auto map_frame_id = pnh.param<std::string>("map_frame_id", "map");
  const auto base_frame_id = pnh.param<std::string>("base_frame_id", "base_link");

  // publish on cmd_vel at a constant frequency
  const double frequency = pnh.param("frequency", 10.0);
  // EC validation
  const double val_dt = pnh.param("val_dt", 0.1);
  const double val_horizon = pnh.param("val_horizon", 0.5);

  const double max_vel_x = pnh.param("max_vel_x", 1.0);
  const double max_vel_y = pnh.param("max_vel_y", 1.0);
  const double max_rot_vel = pnh.param("max_rot_vel", 1.0);

  const double min_vel_x = pnh.param("min_vel_x", -1.0);
  const double min_vel_y = pnh.param("min_vel_y", -1.0);
  const double min_rot_vel = pnh.param("min_rot_vel", -1.0);

  const double acc_lim_x = pnh.param("acc_lim_x", 1.0);
  const double acc_lim_y = pnh.param("acc_lim_y", 1.0);
  const double acc_lim_th = pnh.param("acc_lim_th", 1.0);

  const vec umin = { min_vel_x, min_vel_y, min_rot_vel };
  const vec umax = { max_vel_x, max_vel_y, max_rot_vel };

  // collision
  const double boundary_radius = pnh.param("boundary_radius", 0.7);
  const double search_radius = pnh.param("search_radius", 1.0);
  const double obstacle_threshold = pnh.param("obstacle_threshold", 0.2);
  const double occupied_threshold = pnh.param("occupied_threshold", 0.8);

  // ergodic control
  const double ec_dt = pnh.param("ec_dt", 0.1);
  const double ec_horizon = pnh.param("ec_horizon", 2.0);
  const double target_resolution = pnh.param("target_resolution", 0.1);
  const double expl_weight = pnh.param("expl_weight", 1.0);
  const unsigned int num_basis = pnh.param("num_basis", 10);
  const unsigned int buffer_size = pnh.param("buffer_size", 1e6);
  const unsigned int batch_size = pnh.param("batch_size", 100);

  std::vector<double> control_weights = { 1.0, 1.0, 1.0 };
  pnh.getParam("control_weights", control_weights);

  mat R(3, 3, arma::fill::zeros);
  R(0, 0) = control_weights.at(0);
  R(1, 1) = control_weights.at(1);
  R(2, 2) = control_weights.at(2);

  // dwa
  const double dwa_dt = pnh.param("dwa_dt", 0.1);
  const double dwa_horizon = pnh.param("dwa_horizon", 1.0);
  const double acc_dt = pnh.param("acc_dt", 0.2);
  const unsigned int vx_samples = pnh.param("vx_samples", 3);
  const unsigned int vy_samples = pnh.param("vy_samples", 8);
  const unsigned int vth_samples = pnh.param("vth_samples", 5);

  // target
  XmlRpc::XmlRpcValue means;
  XmlRpc::XmlRpcValue sigmas;
  pnh.getParam("means", means);
  pnh.getParam("sigmas", sigmas);
  const auto num_targets = static_cast<unsigned int>(means.size());

  GaussianList gaussians(num_targets);
  for (unsigned int i = 0; i < num_targets; i++)
  {
    gaussians.at(i) = { { means[i][0], means[i][1] }, { sigmas[i][0], sigmas[i][1] } };
  }

  Target target(gaussians);

  // std::cout << buffer_size << std::endl;
  // std::cout << num_basis << std::endl;
  // //////////////////////////////////////////////////////////////////////////////
  // // grid
  // const auto xmin = -25.0;
  // const auto xmax = 25.0;
  // const auto ymin = -25.0;
  // const auto ymax = 25.0;
  // const auto resolution = 0.05;
  // const auto xsize = ergodic_exploration::axis_length(xmin, xmax, resolution);
  // const auto ysize = ergodic_exploration::axis_length(ymin, ymax, resolution);
  // std::vector<int8_t> map_data(xsize * ysize, 50);
  //
  // GridMap my_grid(xmin, xmax, ymin, ymax, resolution, map_data);
  // // //////////////////////////////////////////////////////////////////////////////
  // nav_msgs::OccupancyGrid grid_msg;
  // grid_msg.header.frame_id = "map";
  // grid_msg.info.resolution = resolution;
  // grid_msg.info.width = xsize;
  // grid_msg.info.height = ysize;
  // grid_msg.info.origin.position.x = xmin;
  // grid_msg.info.origin.position.y = ymin;
  //
  // // TODO: use TF listener to get transform or param server
  // const mat Tbs = ergodic_exploration::transform2d(0.393, 0.0);
  // OccupancyMapper mapper(Tbs);
  // //////////////////////////////////////////////////////////////////////////////
  Collision collision(boundary_radius, search_radius, obstacle_threshold,
                      occupied_threshold);

  ErgodicControl ergodic_control(omni, collision, ec_dt, ec_horizon, target_resolution,
                                 expl_weight, num_basis, buffer_size, batch_size, R, umin,
                                 umax);

  DynamicWindow dwa(collision, dwa_dt, dwa_horizon, acc_dt, acc_lim_x, acc_lim_y,
                    acc_lim_th, max_vel_x, min_vel_x, max_vel_y, min_vel_y, max_rot_vel,
                    min_rot_vel, vx_samples, vy_samples, vth_samples);
  // //////////////////////////////////////////////////////////////////////////////

  ergodic_control.setTarget(target);

  visualization_msgs::MarkerArray marker_array;
  target.markers(marker_array, map_frame_id);
  target_pub.publish(marker_array);

  // // vec u = ergodic_control.control(collision, grid, target, pose);
  vec u = { 0.0, 0.0, 0.0 };
  // vec uref = { 0.7, 0.0, 0.0 };

  bool pose_known = false;
  // bool first_map = false;
  ros::Rate rate(frequency);
  unsigned int i = 0;
  while (nh.ok())
  {
    ros::spinOnce();

    // Update pose
    try
    {
      t_map_base = tfBuffer.lookupTransform(map_frame_id, base_frame_id, ros::Time(0));
      pose(0) = t_map_base.transform.translation.x;
      pose(1) = t_map_base.transform.translation.y;
      pose(2) = tf2::getYaw(t_map_base.transform.rotation);  // wrapped -PI to PI ?
      pose_known = true;
    }

    catch (tf2::TransformException& ex)
    {
      ROS_WARN_NAMED(LOGNAME, "%s", ex.what());
      // continue;
    }

    // pose_known = true;

    // Contol loop
    if (map_received && pose_known)
    {
      // auto t_start = std::chrono::high_resolution_clock::now();

      // if (!first_map)
      // {
      //   ergodic_control.configTarget(grid);
      //   first_map = true;
      // }
      //
      // if (i == 5)
      // {
      //   ergodic_control.configTarget(grid);
      //   i = 0;
      // }

      u = ergodic_control.control(grid, pose);
      if (!validate_control(collision, grid, pose, u, val_dt, val_horizon))
      {
        ROS_INFO_STREAM_NAMED(LOGNAME, "Collision detected! Enabling DWA!");

        // If dwa fails twist is set to zeros
        u = dwa.control(grid, pose, vb, u);
      }

      // ROS_INFO_STREAM_NAMED(LOGNAME, "DWA!");

      // vec u = dwa.control(grid, pose, vb, uref);

      // u.print("u_t:");

      // auto t_end = std::chrono::high_resolution_clock::now();
      // std::cout
      //     << "Hz: "
      //     << 1.0 / (std::chrono::duration<double, std::milli>(t_end - t_start).count() /
      //               1000.0)
      //     << std::endl;

      geometry_msgs::Twist twist_msg;
      twist_msg.linear.x = u(0);
      twist_msg.linear.y = u(1);
      twist_msg.angular.z = u(2);

      cmd_pub.publish(twist_msg);

      nav_msgs::Path trajectory;
      ergodic_control.path(trajectory, map_frame_id);

      path_pub.publish(trajectory);

      i++;

      // break;
    }

    // if (scan_update)
    // {
    //   if (mapper.updateMap(my_grid, scan, pose))
    //   {
    //     grid_msg.data = my_grid.gridData();
    //     map_pub.publish(grid_msg);
    //   }
    //   else
    //   {
    //     ROS_ERROR_NAMED(LOGNAME, "Failed to update map");
    //   }
    //
    //   scan_update = false;
    // }

    rate.sleep();
  }

  ros::waitForShutdown();
  ROS_INFO_STREAM_NAMED(LOGNAME, "Shutting down.");
  return 0;
}
