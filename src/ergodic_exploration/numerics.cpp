/**
 * @file numerics.cpp
 * @author Boston Cleek
 * @date 30 Oct 2020
 * @brief Useful numerical utilities and types
 */

#include <ergodic_exploration/numerics.hpp>

namespace ergodic_exploration
{
std::mt19937_64& get_twister()
{
  static std::random_device rd;
  static std::mt19937_64 gen(rd());
  return gen;
}

double sampleUniformDistribution(double min, double max)
{
  std::uniform_real_distribution<double> dis(min, max);
  return dis(get_twister());
}

}  // namespace ergodic_exploration
