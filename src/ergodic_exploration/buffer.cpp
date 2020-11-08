/**
 * @file buffer.cpp
 * @author Boston Cleek
 * @date 8 Nov 2020
 * @brief Stores past states
 */

#include <iostream>

#include <ergodic_exploration/buffer.hpp>

namespace ergodic_exploration
{
using arma::distr_param;
using arma::ivec;
using arma::randi;

ReplayBuffer::ReplayBuffer(unsigned int buffer_size, unsigned int batch_size)
  : buffer_size_(buffer_size), batch_size_(batch_size)
{
}

void ReplayBuffer::append(const vec& x)
{
  if (memory_.size() < buffer_size_)
  {
    memory_.emplace(memory_.size(), x);
    return;
  }
  std::cout << "WARNING: Buffer is full" << std::endl;
}

void ReplayBuffer::sampleMemory(mat& xt_total, const mat& xt)
{
  if (memory_.empty())
  {
    xt_total = xt;
  }

  // Concatenate the current store states with predicted trajectory
  else if (memory_.size() <= batch_size_)
  {
    const auto num_stored = memory_.size();
    const auto num_states = xt.n_cols + num_stored;
    xt_total.resize(xt.n_rows, num_states);

    for (unsigned int i = 0; i < num_stored; i++)
    {
      // Index is the key
      xt_total.col(i) = memory_.at(i);
    }

    // Copy predicted trajectory to end
    xt_total.cols(num_stored, num_states - 1) = xt;
  }

  // Randomly sample memory and concatenate with predicted trajectory
  else
  {
    const auto num_states = xt.n_cols + batch_size_;
    xt_total.resize(xt.n_rows, num_states);

    // random ints on interval [a b]
    const ivec rand_ints = randi<ivec>(batch_size_, distr_param(0, memory_.size() - 1));

    for (unsigned int i = 0; i < batch_size_; i++)
    {
      // Index is the key
      xt_total.col(i) = memory_.at(rand_ints(i));
    }

    // Copy predicted trajectory to end
    xt_total.cols(batch_size_, num_states - 1) = xt;
  }
}
}  // namespace ergodic_exploration
