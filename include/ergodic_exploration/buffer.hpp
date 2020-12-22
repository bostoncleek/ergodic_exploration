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
 * @file buffer.hpp
 * @author Boston Cleek
 * @date 8 Nov 2020
 * @brief Stores past states
 */
#ifndef BUFFER_HPP
#define BUFFER_HPP

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
   * @param xt - predicted trajectory
   * @return predicted trajectory + sampled states
   */
  mat sampleMemory(const mat& xt) const;

private:
  unsigned int buffer_size_;                      // total number of past states in memory
  unsigned int batch_size_;                       // number of states sampled from memory
  std::unordered_map<unsigned int, vec> memory_;  // past states
};
}  // namespace ergodic_exploration
#endif
