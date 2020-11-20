/**
 * @file numerics.hpp
 * @author Boston Cleek
 * @date 30 Oct 2020
 * @brief Useful numerical utilities
 */
#ifndef NUMERICS_HPP
#define NUMERICS_HPP

#include <cmath>
// #include <numbers>

#include <armadillo>

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
 * @param prob_occu - probability grid cell is occupied represented as a decimal
 * @return entropy
 */
inline double entropy(double p)
{
  // Assign zero information gain
  if (almost_equal(0.0, p) || almost_equal(1.0, p) /*|| p < 0.0*/)
  {
    return 1e-3;
  }

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
 * @return new pose [x, y, theta]
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

}  // namespace ergodic_exploration
#endif
