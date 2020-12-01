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
#include <std_srvs/Empty.h>
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

using arma::eye;
using arma::mat;
using arma::span;
using arma::vec;
using namespace ergodic_exploration;

constexpr char LOGNAME[] = "ergodic exploration";

static GridMap grid, mi_grid;
static vec pose = { 0.0, 0.0, 0.0 };
static vec vb = { 0.0, 0.0, 0.0 };
static double distance_traveled = 0.0;

static bool map_received = false;
static bool mi_received = false;

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
//   distance_traveled += distance(pose(0), pose(1), msg.pose[robot_index].position.x,
//                                 msg.pose[robot_index].position.y);
//
//   pose(0) = msg.pose[robot_index].position.x;
//   pose(1) = msg.pose[robot_index].position.y;
//   pose(2) = normalize_angle_PI(tf2::getYaw(msg.pose[robot_index].orientation));
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

void miCallback(const nav_msgs::OccupancyGrid::ConstPtr& msg)
{
  mi_grid.update(msg);
  // mi_grid.print();
  mi_received = true;
}

int main(int argc, char** argv)
{
  ROS_INFO_STREAM_NAMED(LOGNAME, "Starting exploration");
  ros::init(argc, argv, "exploration");
  ros::NodeHandle nh;
  ros::NodeHandle pnh("~");

  ros::Subscriber map_sub = nh.subscribe("map", 1, mapCallback);
  ros::Subscriber mi_map_sub = nh.subscribe("mi_map", 1, miCallback);
  ros::Subscriber odom_sub = nh.subscribe("odom", 1, odomCallback);
  // ros::Subscriber model_sub = nh.subscribe("/gazebo/model_states", 1, modelCallBack);

  ros::Publisher cmd_pub = nh.advertise<geometry_msgs::Twist>("cmd_vel", 1);
  ros::Publisher path_pub = nh.advertise<nav_msgs::Path>("trajectory", 1);
  ros::Publisher dwa_path_pub = nh.advertise<nav_msgs::Path>("dwa_trajectory", 1);

  ros::Publisher target_pub =
      nh.advertise<visualization_msgs::MarkerArray>("target", 1, true);

  ros::ServiceClient mi_client = nh.serviceClient<std_srvs::Empty>("start_mi");


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

  //////////////////////////////////////////////////////////////////////////////
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

  vec u = { 0.0, 0.0, 0.0 };

  bool pose_known = false;
  bool target_set = false;
  bool update_target = false;
  bool update_mi = true;

  ros::Rate rate(frequency);
  while (nh.ok())
  {
    ros::spinOnce();

    // pose_known = true;

    // Update pose
    try
    {
      t_map_base = tfBuffer.lookupTransform(map_frame_id, base_frame_id, ros::Time(0));

      distance_traveled += distance(pose(0), pose(1), t_map_base.transform.translation.x,
                                    t_map_base.transform.translation.y);

      pose(0) = t_map_base.transform.translation.x;
      pose(1) = t_map_base.transform.translation.y;
      pose(2) = normalize_angle_PI(tf2::getYaw(t_map_base.transform.rotation));  // wrapped -PI to PI ?
      pose_known = true;
    }

    catch (tf2::TransformException& ex)
    {
      ROS_WARN_NAMED(LOGNAME, "%s", ex.what());
      // continue;
    }

    // TODO: Check if map has grown and if so update mi
    if (distance_traveled > 5.0)
    {
      ROS_INFO_NAMED(LOGNAME, "Robot traveled: %f", distance_traveled);
      update_mi = true;
    }

    // TODO: publish 0 twist when updating mi
    if (map_received && update_mi)
    {
      ROS_INFO_STREAM_NAMED(LOGNAME, "Checking if mi service exists");
      if (mi_client.exists())
      {
        ROS_INFO_STREAM_NAMED(LOGNAME, "Calling mi service");

        std_srvs::Empty srv;
        mi_client.call(srv);

        update_target = true;
        update_mi = false;
        distance_traveled = 0.0;
      }
    }

    // Contol loop
    if (map_received && mi_received && pose_known)
    {
      if (!target_set || update_target)
      {
        ergodic_control.configTarget(mi_grid);
        target_set = true;
        update_target = false;
      }

      if (target_set)
      {
        u = ergodic_control.control(grid, pose);
        if (!validate_control(collision, grid, pose, u, val_dt, val_horizon))
        {
          ROS_INFO_STREAM_NAMED(LOGNAME, "Collision detected! Enabling DWA!");

          // If dwa fails twist is set to zeros
          u = dwa.control(grid, pose, vb, u);

          nav_msgs::Path dwa_traj;
          dwa.path(dwa_traj, pose, u, map_frame_id);

          dwa_path_pub.publish(dwa_traj);

        }
      }

      geometry_msgs::Twist twist_msg;
      twist_msg.linear.x = u(0);
      twist_msg.linear.y = u(1);
      twist_msg.angular.z = u(2);

      cmd_pub.publish(twist_msg);

      nav_msgs::Path trajectory;
      ergodic_control.path(trajectory, map_frame_id);

      path_pub.publish(trajectory);
    }

    rate.sleep();
  }

  ros::waitForShutdown();
  ROS_INFO_STREAM_NAMED(LOGNAME, "Shutting down.");
  return 0;
}
