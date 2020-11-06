/**
 * @file omni.hpp
 * @author Boston Cleek
 * @date 30 Oct 2020
 * @brief Kinematic 4 wheel omni directional robot
 */

#pragma once

#include <cmath>
#include <armadillo>

#include <ergodic_exploration/types.hpp>

namespace ergodic_exploration
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
    , action_space_(4)
    , state_space_(3)
  {
  }

  /**
   * @brief Convert wheel velocities to a body frame twist
   * @param u - control [u0, u1, u2, u3] (column vector)
   * @return twist in body frame Vb = [vx, vy, w]
   */
  Twist2D wheels2Twist(const vec u) const
  {
    const auto l = 1.0 / (wheel_base_x + wheel_base_y);

    // pseudo inverse of jacobian matrix
    const mat Hp = { { 1.0, 1.0, 1.0, 1.0 }, { -1.0, 1.0, -1.0, 1.0 }, { -l, l, l, -l } };

    const vec v = (wheel_radius / 4.0) * Hp * u;

    return { v(0), v(1), v(2) };
  }

  /**
   * @brief Kinematic model of 4 mecanum wheel robot
   * @param x - state [x, y, theta] (column vector)
   * @param u - control [u0, u1, u2, u3] (column vector)
   * @return [xdot, ydot, thetadot] = f(x,u) (column vector)
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
   * @param x - state [x, y, theta] (column vector)
   * @param u - control [u0, u1, u2, u3] (column vector)
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
   * @param x - state [x, y, theta] (column vector)
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

  double wheel_radius;  // radius of wheel
  double wheel_base_x;  // distance from chassis center to wheel center along x-axis
  double wheel_base_y;  // distance from chassis center to wheel center along y-axis
  unsigned int action_space_;  // control space dimension
  unsigned int state_space_;   // states space dimension
};

}  // namespace ergodic_exploration
