/**
 * @file buffer.hpp
 * @author Boston Cleek
 * @date 8 Nov 2020
 * @brief Stores past states
 */

#pragma once

#include <unordered_map>
#include <armadillo>

namespace ergodic_exploration
{
using arma::mat;
using arma::vec;

/** @brief Store and smaple past states */
class ReplayBuffer
{
public:
  /**
   * @brief Constructor
   * @param buffer_size - max number of states stored
   * @param batch_size - number of states randomly sampled from memory
   */
  ReplayBuffer(unsigned int buffer_size, unsigned int batch_size);

  /**
   * @brief Add current state to memory
   * @param x - current state
   */
  void append(const vec& x);

  /**
   * @brief Sample states from memory
   * @param xt_total[out] - predicted trajectory + sampled states
   * @paramxt - predicted trajectory
   */
  void sampleMemory(mat& xt_total, const mat& xt);

private:
  unsigned int buffer_size_;                      // total number of past states in memory
  unsigned int batch_size_;                       // number of states sampled from memory
  std::unordered_map<unsigned int, vec> memory_;  // past states
};

}  // namespace ergodic_exploration
