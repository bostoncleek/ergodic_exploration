/**
 * @file numerics.hpp
 * @author Boston Cleek
 * @date 30 Oct 2020
 * @brief Useful numerical utilities
 */

#pragma once

#include <cmath>
// #include <numbers>

namespace ergodic_exploration
{
// TODO: why gcc cant fine numbers
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

  if (rad < 0)
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

  if (rad < 0)
  {
    rad += 2.0 * PI;
  }

  return rad;
}

}  // namespace ergodic_exploration
