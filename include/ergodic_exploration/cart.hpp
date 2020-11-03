/**
 * @file cart.hpp
 * @author Boston Cleek
 * @date 28 Oct 2020
 * @brief Kinematic 2 wheel cart
 */

#pragma once

#include <cmath>
#include <armadillo>

namespace ergodic_exploration
{
using arma::mat;
using arma::vec;

/**
 * @brief Kinematic model of 2 wheel differential drive robot
 * @details The state is [x, y, theta] and controls are the velocities of each wheel [uL, uR]
 */
struct Cart
{
  /**
   * @brief Constructor
   * @param wheel_radius - radius of wheel
   * @param wheel_base - distance from point at the center in between both wheels to the
   * center of a wheel
   */
  Cart(double wheel_radius, double wheel_base)
    : wheel_radius(wheel_radius), wheel_base(wheel_base)
  {
  }

  /**
   * @brief Kinematic model of a 2 wheel differential drive robot
   * @param x - state [x, y, theta] (column vector)
   * @param u - control [uL, uR] (column vector)
   * @return [xdot, ydot, thetadot] = f(x,u) (column vector)
   */
  vec operator()(const vec x, const vec u) const
  {
    vec xdot(3);
    xdot(0) = (u(0) + u(1)) * std::cos(x(2));
    xdot(1) = (u(0) + u(1)) * std::sin(x(2));
    xdot(2) = (u(1) - u(0)) / wheel_base;

    return (wheel_radius / 2.0) * xdot;
  }

  /**
   * @brief Jacobian of the model with respect to the state
   * @param x - state [x, y, theta] (column vector)
   * @param u - control [uL, uR] (column vector)
   * @return A = D1(f(x,u)) of shape (3x3)
   */
  mat fdx(const vec x, const vec u) const
  {
    mat A(3, 3, arma::fill::zeros);

    const auto df0dth = -(wheel_radius / 2.0) * (u(0) + u(1)) * std::sin(x(2));
    const auto df1dth = (wheel_radius / 2.0) * (u(0) + u(1)) * std::cos(x(2));

    A(0, 2) = df0dth;
    A(1, 2) = df1dth;

    return A;
  }

  /**
   * @brief Jacobian of the model with respect to the control
   * @param x - state [x, y, theta] (column vector)
   * @return B = D2(f(x,u)) of shape (3x2)
   */
  mat fdu(const vec x) const
  {
    mat B(3, 2);

    B(0, 0) = std::cos(x(2));
    B(0, 1) = std::cos(x(2));

    B(1, 0) = std::sin(x(2));
    B(1, 1) = std::sin(x(2));

    B(2, 0) = -1.0 / wheel_base;
    B(2, 1) = 1.0 / wheel_base;

    return (wheel_radius / 2.0) * B;
  }

  double wheel_radius;  // radius of wheel
  double wheel_base;    // distance from robot center to wheel center
};

/**
 * @brief Kinematic model of 2 wheel differential drive robot
 * @details The state is [x, y, theta] and controls are the linear and angular velocities [v, w]
 */
struct SimpleCart
{
  /**
   * @brief Kinematic model of a 2 wheel differential drive robot
   * @param x - state [x, y, theta] (column vector)
   * @param u - control [v, w] (column vector)
   * @return xdot = f(x,u)
   */
  vec operator()(const vec x, const vec u) const
  {
    return { u(0) * std::cos(x(2)), u(0) * std::sin(x(2)), u(1) };
  }

  /**
   * @brief Jacobian of the model with respect to the state
   * @param x - state [x, y, theta] (column vector)
   * @param u - control [v, w](column vector)
   * @return A = D1(f(x,u)) of shape (3x3)
   */
  mat fdx(const vec x, const vec u) const
  {
    mat A(3, 3, arma::fill::zeros);
    A(0, 2) = -u(0) * std::sin(x(2));
    A(1, 2) = u(0) * std::cos(x(2));
    return A;
  }

  /**
   * @brief Jacobian of the model with respect to the control
   * @param x - state [x, y, theta] (column vector)
   * @return B = D2(f(x,u)) of shape (3x2)
   */
  mat fdu(const vec x) const
  {
    mat B(3, 2, arma::fill::zeros);

    B(0, 0) = std::cos(x(2));
    B(1, 0) = std::sin(x(2));
    B(2, 1) = 1.0;

    return B;
  }
};
}  // namespace ergodic_exploration