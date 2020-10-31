/**
 * @file numerics.hpp
 * @author Boston Cleek
 * @date 30 Oct 2020
 * @brief Useful numerical utilities and types
 */

#pragma once

#include <cmath>

namespace ergodic_exploration
{
/** @ brief A 2-Dimensional twist */
struct Twist2D
{
  Twist2D(double vx, double vy, double w) : vx(vx), vy(vy), w(w)
  {
  }

  Twist2D() : vx(0.0), vy(0.0), w(0.0)
  {
  }

  double vx;  // linear x velocity
  double vy;  // linear y velocity
  double w;   // rotation about z-axis
};

/**
 * @brief approximately compare two floating-point numbers
 * @param d1 - a number to compare
 * @param d2 - a second number to compare
 * @param epsilon - absolute threshold required for equality
 * @return true if abs(d1 - d2) < epsilon
 */
constexpr bool almost_equal(double d1, double d2, double epsilon = 1.0e-12)
{
  return std::fabs(d1 - d2) < epsilon ? true : false;
}

}  // namespace ergodic_exploration
