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

mat ReplayBuffer::sampleMemory(const mat& xt) const
{
  if (memory_.empty())
  {
    return xt;
  }

  // total states
  mat xt_total;

  // Concatenate the current store states with predicted trajectory
  if (memory_.size() <= batch_size_)
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

  return xt_total;
}
}  // namespace ergodic_exploration
