/**
 * @file numerics.hpp
 * @author Boston Cleek
 * @date 30 Oct 2020
 * @brief Useful numerical utilities
 */

#pragma once

#include <cmath>
#include <random>

namespace ergodic_exploration
{
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

/** @brief Returns a random number engine */
std::mt19937_64& get_twister();

/**
 * @brief Samples a real uniform distribution
 * @param min - lower bound
 * @param max - upper bound
 * @return a random sample
 */
double sampleUniformDistribution(const double min, const double max);

}  // namespace ergodic_exploration
