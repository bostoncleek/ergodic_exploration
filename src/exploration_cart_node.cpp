/*********************************************************************
 * BSD 3-Clause License
 *
 * Copyright (c) 2020 Northwestern University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 *  * Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/
/**
 * @file exploration_node.cpp
 * @author Boston Cleek
 * @date 4 Dec 2020
 * @brief Ergodic exploration using a differential drive cart
 *
 * @PARAMETERS:
 *  * map_frame_id - map frame id
 *  * base_frame_id - base link frame id
 *  * frequency - control loop frequency (Hz)
 *  * val_dt - control validation time step for collision detection (s)
 *  * val_horizon - control validation horizon for collision detection (s)
 *  * max_vel_x - max x velcocity (m/s)
 *  * max_rot_vel - max ratotional velcocity (rad/s)
 *  * min_vel_x - min x velcocity (m/s)
 *  * min_rot_vel - min ratotional velcocity (rad/s)
 *  * acc_lim_x - x acceleration limit (m/s^2)
 *  * acc_lim_th - rotational acceleration limit (rad/s^2)
 *  * boundary_radius - bounding radius around robot (m)
 *  * search_radius - max search radius for collision detection (m)
 *  * obstacle_threshold - obstacles within radius from boundary are cosidered collisions (m)
 *  * occupied_threshold - occupancy grid cell probability to be considered an obstacle [0 1]
 *  * ec_dt - time step used in integration (s)
 *  * ec_horizon - control horizon (s)
 *  * target_resolution - target grid resolution (m)
 *  * expl_weight - ergodic exploration weight
 *  * num_basis - number of basis functions
 *  * buffer_size - total number of past states stored
 *  * batch_size - number of past states randomly sampled in each control loop
 *  * control_weights - weights on twist [vx w]
 *  * dwa_dt - time step used in integration (s)
 *  * dwa_horizon - control horizon (s)
 *  * acc_dt - time step the acceleration limits are applied (s)
 *  * vx_samples - number of x velcocity samples
 *  * vth_samples - number of rotational velcocity samples
 *  * means - target x and y means (m)
 *  * sigmas - target x and y standard deviations (m)
 *
 * @PUBLISHES:
 *  * cmd_vel (geometry_msgs/Twist) - body twist
 *  * trajectory (nav_msgs/Path) - ergodic controller optimzed trajectory
 *  * dwa_trajectory (nav_msgs/Path) - dynamic window trajectory
 *  * target (visualization_msgs/MarkerArray) - target distribution
 *
 * @SUBSCRIBES:
 *  * map (nav_msgs/OccupancyGrid) - occupancy grid
 *  * odom (nav_msgs/Odometry) - robot's odometry
 */

#include <iostream>

#include <ergodic_exploration/models/cart.hpp>
#include <ergodic_exploration/exploration.hpp>

using arma::mat;
using arma::vec;

using ergodic_exploration::Collision;
using ergodic_exploration::DynamicWindow;
using ergodic_exploration::ErgodicControl;
using ergodic_exploration::Exploration;
using ergodic_exploration::GaussianList;
using ergodic_exploration::Target;
using ergodic_exploration::models::SimpleCart;

constexpr char LOGNAME[] = "cart exploration";

int main(int argc, char** argv)
{
  ROS_INFO_STREAM_NAMED(LOGNAME, "Starting exploration");
  ros::init(argc, argv, "exploration");
  ros::NodeHandle nh;
  ros::NodeHandle pnh("~");

  ros::AsyncSpinner spinner(1);
  spinner.start();

  arma::arma_rng::set_seed_random();

  // TODO: Add noise to motion model
  // motion model
  const SimpleCart cart;

  //////////////////////////////////////////////////////////////////////////////
  const auto map_frame_id = pnh.param<std::string>("map_frame_id", "map");
  const auto base_frame_id = pnh.param<std::string>("base_frame_id", "base_link");

  // publish on cmd_vel at a constant frequency
  const double frequency = pnh.param("frequency", 10.0);
  // EC validation
  const double val_dt = pnh.param("val_dt", 0.1);
  const double val_horizon = pnh.param("val_horizon", 0.5);

  const double max_vel_x = pnh.param("max_vel_x", 1.0);
  const double max_rot_vel = pnh.param("max_rot_vel", 1.0);

  const double min_vel_x = pnh.param("min_vel_x", -1.0);
  const double min_rot_vel = pnh.param("min_rot_vel", -1.0);

  const double acc_lim_x = pnh.param("acc_lim_x", 1.0);
  const double acc_lim_th = pnh.param("acc_lim_th", 1.0);

  const vec umin = { min_vel_x, 0.0, min_rot_vel };
  const vec umax = { max_vel_x, 0.0, max_rot_vel };

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

  std::vector<double> control_weights = { 1.0, 1.0 };
  pnh.getParam("control_weights", control_weights);

  // this is technically not the inverse but the penalty for vy is not used
  mat Rinv(3, 3, arma::fill::zeros);
  Rinv(0, 0) = 1.0 / control_weights.at(0);
  Rinv(1, 1) = 0.0;
  Rinv(2, 2) = 1.0 / control_weights.at(1);

  // dwa
  const double dwa_dt = pnh.param("dwa_dt", 0.1);
  const double dwa_horizon = pnh.param("dwa_horizon", 1.0);
  const double acc_dt = pnh.param("acc_dt", 0.2);
  const unsigned int vx_samples = pnh.param("vx_samples", 3);
  const unsigned int vth_samples = pnh.param("vth_samples", 5);

  if (dwa_horizon > ec_horizon)
  {
    ROS_ERROR_STREAM_NAMED(
        LOGNAME, "Dynamic window horizon is greater than the ergodic control horizon");
    ros::shutdown();
  }

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

  const Target target(gaussians);

  //////////////////////////////////////////////////////////////////////////////
  const Collision collision(boundary_radius, search_radius, obstacle_threshold,
                            occupied_threshold);

  const ErgodicControl ergodic_control(cart, collision, ec_dt, ec_horizon,
                                       target_resolution, expl_weight, num_basis,
                                       buffer_size, batch_size, Rinv, umin, umax);

  const DynamicWindow dwa(collision, dwa_dt, dwa_horizon, acc_dt, acc_lim_x, 0.0,
                          acc_lim_th, max_vel_x, min_vel_x, 0.0, 0.0, max_rot_vel,
                          min_rot_vel, vx_samples, 0.0, vth_samples);

  ergodic_exploration::Exploration<SimpleCart> exploration(nh, ergodic_control, collision,
                                                           dwa);

  // start exploring
  exploration.control(target, map_frame_id, base_frame_id, frequency, val_dt, val_horizon);

  //////////////////////////////////////////////////////////////////////////////

  // ros::waitForShutdown();
  ROS_INFO_STREAM_NAMED(LOGNAME, "Shutting down.");
  return 0;
}
