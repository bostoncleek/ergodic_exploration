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
 * @file omni.hpp
 * @author Boston Cleek
 * @date 30 Oct 2020
 * @brief Kinematic omni directional models control wheel velocities or body twist
 */
#ifndef OMNI_HPP
#define OMNI_HPP

#include <cmath>
#include <armadillo>

#include <ergodic_exploration/numerics.hpp>

namespace ergodic_exploration
{
namespace models
{
using arma::mat;
using arma::vec;

/**
 * @brief Kinematic model of 4 mecanum wheel robot
 * @details The state is [x, y, theta] and controls are the angular velocities of each
 * wheel [u0, u1, u2, u3] corresponding to (front left, front right, rear right, rear
 * left). Assumes the mecanum wheel rollers are at +/- 45 degrees
 */
struct Mecanum
{
  /**
   * @brief Constructor
   * @param wheel_radius - radius of wheel
   * @param wheel_base_x - distance from chassis center to wheel center along x-axis
   * @param wheel_base_y - distance from chassis center to wheel center along y-axis
   */
  Mecanum(double wheel_radius, double wheel_base_x, double wheel_base_y)
    : wheel_radius(wheel_radius)
    , wheel_base_x(wheel_base_x)
    , wheel_base_y(wheel_base_y)
    , state_space(3)
  {
  }

  /**
   * @brief Convert wheel velocities to a body frame twist
   * @param u - control [u0, u1, u2, u3]
   * @return twist in body frame Vb = [vx, vy, w]
   */
  vec wheels2Twist(const vec u) const
  {
    const auto l = 1.0 / (wheel_base_x + wheel_base_y);

    // pseudo inverse of jacobian matrix
    const mat Hp = { { 1.0, 1.0, 1.0, 1.0 }, { -1.0, 1.0, -1.0, 1.0 }, { -l, l, l, -l } };

    const vec vb = (wheel_radius / 4.0) * Hp * u;

    return { vb(0), vb(1), vb(2) };
  }

  /**
   * @brief Kinematic model of 4 mecanum wheel robot
   * @param x - state [x, y, theta]
   * @param u - control [u0, u1, u2, u3]
   * @return [xdot, ydot, thetadot] = f(x,u)
   */
  vec operator()(const vec x, const vec u) const
  {
    vec xdot(3);
    const auto s = (wheel_radius / 4.0) * std::sin(x(2));
    const auto c = (wheel_radius / 4.0) * std::cos(x(2));
    const auto l = wheel_radius / (4.0 * (wheel_base_x + wheel_base_y));

    xdot(0) = u(0) * (s + c) + u(1) * (-s + c) + u(2) * (s + c) + u(3) * (-s + c);
    xdot(1) = u(0) * (s - c) + u(1) * (s + c) + u(2) * (s - c) + u(3) * (s + c);
    xdot(2) = -u(0) * l + u(1) * l + u(2) * l - u(3) * l;

    return xdot;
  }

  /**
   * @brief Jacobian of the model with respect to the state
   * @param x - state [x, y, theta]
   * @param u - control [u0, u1, u2, u3]
   * @return A = D1(f(x,u)) of shape (3x3)
   */
  mat fdx(const vec x, const vec u) const
  {
    mat A(3, 3, arma::fill::zeros);

    const auto s = (wheel_radius / 4.0) * std::sin(x(2));
    const auto c = (wheel_radius / 4.0) * std::cos(x(2));

    const auto df0dth =
        u(0) * (-s + c) + u(1) * (-s - c) + u(2) * (-s + c) + u(3) * (-s - c);
    const auto df1dth =
        u(0) * (s + c) + u(1) * (-s + c) + u(2) * (s + c) + u(3) * (-s + c);

    A(0, 2) = df0dth;
    A(1, 2) = df1dth;

    return A;
  }

  /**
   * @brief Jacobian of the model with respect to the control
   * @param x - state [x, y, theta]
   * @return B = D2(f(x,u)) of shape (3x4)
   */
  mat fdu(const vec x) const
  {
    const auto s = (wheel_radius / 4.0) * std::sin(x(2));
    const auto c = (wheel_radius / 4.0) * std::cos(x(2));
    const auto l = wheel_radius / (4.0 * (wheel_base_x + wheel_base_y));

    const mat B = { { s + c, -s + c, s + c, -s + c },
                    { s - c, s + c, s - c, s + c },
                    { -l, l, l, -l } };
    return B;
  }

  double wheel_radius;       // radius of wheel
  double wheel_base_x;       // distance from chassis center to wheel center along x-axis
  double wheel_base_y;       // distance from chassis center to wheel center along y-axis
  unsigned int state_space;  // states space dimension
};

/**
 * @brief Kinematic model of omni directonal robot
 * @details The state is [x, y, theta] and controls are the linear and
 * angular velocities [vx, vy, w] (body twist)
 */
struct Omni
{
  /** @brief Constructor */
  Omni() : state_space(3)
  {
  }

  /**
   * @brief Kinematic model of 4 mecanum wheel robot
   * @param x - state [x, y, theta]
   * @param u - body twist control [vx, vy, w]
   * @return [xdot, ydot, thetadot] = f(x,u)
   */
  vec operator()(const vec x, const vec u) const
  {
    const auto xdot = u(0) * std::cos(x(2)) - u(1) * std::sin(x(2));
    const auto ydot = u(0) * std::sin(x(2)) + u(1) * std::cos(x(2));
    return { xdot, ydot, u(2) };
    // return { xdot, ydot, 0.0 };
    // return { u(0), u(1), u(2) };
  }

  /**
   * @brief Jacobian of the model with respect to the state
   * @param x - state [x, y, theta]
   * @param u - body twist control [vx, vy, w]
   * @return A = D1(f(x,u)) of shape (3x3)
   */
  mat fdx(const vec x, const vec u) const
  {
    mat A(3, 3, arma::fill::zeros);
    A(0, 2) = -u(0) * std::sin(x(2)) - u(1) * std::cos(x(2));
    A(1, 2) = u(0) * std::cos(x(2)) - u(1) * std::sin(x(2));
    return A;
  }

  /**
   * @brief Jacobian of the model with respect to the control
   * @param x - state [x, y, theta]
   * @return B = D2(f(x,u)) of shape (3x3)
   */
  mat fdu(const vec x) const
  {
    // mat B(3,3, arma::fill::eye);
    const mat B = { { std::cos(x(2)), -std::sin(x(2)), 0.0 },
                    { std::sin(x(2)), std::cos(x(2)), 0.0 },
                    { 0.0, 0.0, 1.0 } };
    return B;
  }

  unsigned int state_space;  // states space dimension
};
}  // namespace models
}  // namespace ergodic_exploration
#endif
