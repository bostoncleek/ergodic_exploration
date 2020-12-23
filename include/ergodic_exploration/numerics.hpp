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
 * @file numerics.hpp
 * @author Boston Cleek
 * @date 30 Oct 2020
 * @brief Useful numerical utilities
 */
#ifndef NUMERICS_HPP
#define NUMERICS_HPP

#include <cmath>
#include <algorithm>
// #include <numbers>

#include <armadillo>

#include <nav_msgs/Path.h>
#include <tf2/LinearMath/Quaternion.h>

#include <ergodic_exploration/collision.hpp>

namespace ergodic_exploration
{
using arma::mat;
using arma::vec;

// TODO: why gcc cant find numbers
constexpr double PI = 3.14159265358979323846;

/**
 * @brief approximately compare two floating-point numbers
 * @param d1 - a number to compare
 * @param d2 - a second number to compare
 * @param epsilon - absolute threshold required for equality
 * @return true if abs(d1 - d2) < epsilon
 */
inline bool almost_equal(double d1, double d2, double epsilon = 1.0e-12)
{
  return std::fabs(d1 - d2) < epsilon ? true : false;
}

/**
 * @brief Wraps angle between -pi and pi
 * @param rad - angle in radians
 * @return wrapped angle in radians
 */
inline double normalize_angle_PI(double rad)
{
  // floating point remainder essentially this is fmod
  const auto q = std::floor((rad + PI) / (2.0 * PI));
  rad = (rad + PI) - q * 2.0 * PI;

  if (rad < 0.0)
  {
    rad += 2.0 * PI;
  }

  return (rad - PI);
}

/**
 * @brief Wraps angle between 0 and 2pi or 0 to -2pi
 * @param rad - angle in radians
 * @return wrapped angle in radians
 */
inline double normalize_angle_2PI(double rad)
{
  // floating point remainder essentially this is fmod
  const auto q = std::floor(rad / (2.0 * PI));
  rad = (rad)-q * 2.0 * PI;

  if (rad < 0.0)
  {
    rad += 2.0 * PI;
  }

  return rad;
}

/**
 * @brief Get yaw from quaternion
 * @param qx - x-axis rotation component
 * @param qy - y-axis rotation component
 * @param qz - z-axis rotation component
 * @param qw - rotation magnitude
 * @return yaw (rad)
 */
inline double getYaw(double qx, double qy, double qz, double qw)
{
  // TODO: add unit test
  // Source: https://github.com/ros/geometry2/blob/noetic-devel/tf2/include/tf2/impl/utils.h#L122
  const auto sqx = qx * qx;
  const auto sqy = qy * qy;
  const auto sqz = qz * qz;
  const auto sqw = qw * qw;

  // Normalization added from urdfom_headers
  const auto sarg = -2.0 * (qx * qz - qw * qy) / (sqx + sqy + sqz + sqw);

  // Cases derived from https://orbitalstation.wordpress.com/tag/quaternion/
  if (sarg < -0.99999 || almost_equal(sarg, -0.99999))
  {
    return -2.0 * std::atan2(qy, qx);
  }

  else if (sarg > 0.99999 || almost_equal(sarg, 0.99999))
  {
    return 2.0 * std::atan2(qy, qx);
  }

  return std::atan2(2.0 * (qx * qy + qw * qz), sqw + sqx - sqy - sqz);
}

/**
 * @brief Euclidean distance between two points
 * @param x0 - x-position point 0
 * @param y0 - y-position point 0
 * @param x1 - x-position point 1
 * @param y1 - y-position point 1
 * @return euclidean distance
 */
inline double distance(double x0, double y0, double x1, double y1)
{
  const auto dx = x1 - x0;
  const auto dy = y1 - y0;
  return std::sqrt(dx * dx + dy * dy);
}

/**
 * @brief Entropy of a single grid cell
 * @param p - probability grid cell is occupied represented as a decimal
 * @return entropy
 */
inline double entropy(double p)
{
  // Assign zero information gain
  if (almost_equal(0.0, p) || almost_equal(1.0, p) /*|| p < 0.0*/)
  {
    return 1e-3;
  }

  // unknowm: p = -1 => entropy(0.5) = 0.7
  else if (p < 0.0)
  {
    return 0.7;
  }

  return -p * std::log(p) - (1.0 - p) * std::log(1.0 - p);
}

/**
 * @brief Convert polar to cartesian coordinates
 * @param angle - angle in radians
 * @param range - range measurement
 */
inline vec polar2Cartesian(double angle, double range)
{
  const auto x = range * std::cos(angle);
  const auto y = range * std::sin(angle);
  return { x, y };
}

/**
 * @brief Convert polar to cartesian homogenous coordinates
 * @param angle - angle in radians
 * @param range - range measurement
 */
inline vec polar2CartesianHomo(double angle, double range)
{
  const auto x = range * std::cos(angle);
  const auto y = range * std::sin(angle);
  return { x, y, 1.0 };
}

/**
 * @brief Construct 2D transformation matrix
 * @param x - x position
 * @param y - y position
 * @param angle - yaw in radians
 * @details 2D transformation
 */
inline mat transform2d(double x, double y, double angle)
{
  const mat trans2d = { { std::cos(angle), -std::sin(angle), x },
                        { std::sin(angle), std::cos(angle), y },
                        { 0.0, 0.0, 1.0 } };

  return trans2d;
}

/**
 * @brief Construct 2D transformation matrix
 * @param x - x position
 * @param y - y position
 * @details 2D transformation
 */
inline mat transform2d(double x, double y)
{
  const mat trans2d = { { 1.0, 0.0, x }, { 0.0, 1.0, y }, { 0.0, 0.0, 1.0 } };

  return trans2d;
}

/**
 * @brief Construct 2D transformation
 * @param angle - yaw in radians
 * @details 2D transformation
 */
inline mat transform2d(double angle)
{
  const mat trans2d = { { std::cos(angle), -std::sin(angle), 0.0 },
                        { std::sin(angle), std::cos(angle), 0.0 },
                        { 0.0, 0.0, 1.0 } };
  return trans2d;
}

/**
 * @brief Construct 2D transformation inverse
 * @param trans2d - 2D transformation
 * @details 2D transformation inverse
 */
inline mat transform2dInv(const mat& trans2d)
{
  // R^T flip sign in sin
  const auto stheta = -trans2d(1, 0);
  const auto ctheta = trans2d(0, 0);
  const auto theta = std::atan2(stheta, ctheta);

  // p' = -R^T * p
  const auto x = -(ctheta * trans2d(0, 2) - stheta * trans2d(1, 2));
  const auto y = -(stheta * trans2d(0, 2) + ctheta * trans2d(1, 2));

  return transform2d(theta, x, y);
}

/**
 * @brief Integrate a constant twist
 * @param x - current state [x, y, theta]
 * @param vb - current twist [vx, vy, w]
 * @param dt - time step
 * @return new pose
 */
inline vec integrate_twist(const vec& x, const vec& u, double dt)
{
  // Eqn. 13.35 and 13.36 pg 471 Modern Robotics
  // displacement b to b' (dx, dy, dth)
  vec dqb(3);

  // no rotation
  if (almost_equal(u(2), 0.0))
  {
    dqb(0) = u(0) * dt;
    dqb(1) = u(1) * dt;
    dqb(2) = 0.0;
  }

  else
  {
    const vec vb = u * dt;
    dqb(0) = (vb(0) * std::sin(vb(2)) + vb(1) * (std::cos(vb(2)) - 1.0)) / vb(2);

    dqb(1) = (vb(1) * std::sin(vb(2)) + vb(0) * (1.0 - std::cos(vb(2)))) / vb(2);

    dqb(2) = vb(2);
  }

  return x + transform2d(x(2)) * dqb;
}

/**
 * @brief Determine if control will cause a collision
 * @param collision - collision detector
 * @param grid - grid map
 * @param x0 - initial state
 * @param u - twist [vx, vy, w]
 * @param dt - time step in integration
 * @param horizon - length of integration
 * @return true if the control is collision free
 * @details The control is assumed to be constant and a twist is
 * integrated for a fixed amout of time
 */
inline bool validate_control(const Collision& collision, const GridMap& grid,
                             const vec& x0, const vec& u, double dt, double horizon)
{
  vec x = x0;
  const auto steps = static_cast<unsigned int>(std::abs(horizon / dt));

  for (unsigned int i = 0; i < steps; i++)
  {
    x = integrate_twist(x, u, dt);
    x(2) = normalize_angle_PI(x(2));

    if (collision.collisionCheck(grid, x))
    {
      return false;
    }
  }

  return true;
}

/**
 * @brief Visualize path from following a constant twist
 * @param x0 - current state
 * @param u - twist [vx, vy, w]
 * @param dt - time step
 * @param horizon - control horizon
 * @return trajectory
 */
inline nav_msgs::Path constTwistPath(const std::string& map_frame_id, const vec& x0,
                                     const vec& u, double dt, double horizon)
{
  nav_msgs::Path path;
  path.header.frame_id = map_frame_id;

  const auto steps = static_cast<unsigned int>(std::abs(horizon / dt));
  path.poses.resize(steps);

  vec x = x0;
  for (unsigned int i = 0; i < steps; i++)
  {
    x = integrate_twist(x, u, dt);

    path.poses.at(i).pose.position.x = x(0);
    path.poses.at(i).pose.position.y = x(1);

    tf2::Quaternion quat;
    quat.setRPY(0.0, 0.0, normalize_angle_PI(x(2)));

    path.poses.at(i).pose.orientation.x = quat.x();
    path.poses.at(i).pose.orientation.y = quat.y();
    path.poses.at(i).pose.orientation.z = quat.z();
    path.poses.at(i).pose.orientation.w = quat.w();
  }

  return path;
}

}  // namespace ergodic_exploration
#endif
