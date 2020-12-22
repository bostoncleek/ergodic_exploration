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
 * @file dynamic_window.hpp
 * @author Boston Cleek
 * @date 5 Nov 2020
 * @brief Dynamic window control
 */
#ifndef DYNAMIC_WINDOW_HPP
#define DYNAMIC_WINDOW_HPP

#include <cmath>

#include <ergodic_exploration/grid.hpp>
#include <ergodic_exploration/collision.hpp>
#include <ergodic_exploration/numerics.hpp>

namespace ergodic_exploration
{
using arma::mat;

/** @brief Dynamic window approach */
class DynamicWindow
{
public:
  /**
   * @brief Constructor
   * @param collision - collision detector
   * @param dt - time step in integration
   * @param horizon - control horizon
   * @param acc_dt - time acceleration limit is applied for determining window limits
   * @param acc_lim_x - x acceleration limit
   * @param acc_lim_y - y acceleration limit
   * @param acc_lim_th - rotational acceleration limit
   * @param max_vel_x - max x velocity
   * @param min_vel_x - min x velocity
   * @param max_vel_y - max y velocity
   * @param min_vel_y - min y velocity
   * @param max_rot_vel - max rotational velocity
   * @param min_rot_vel - min rotational velocity
   * @param vx_samples - number of x velocity samples
   * @param vy_samples - number of y velocity samples
   * @param vth_samples - number of rotational velocity samples
   * @details Given a reference control finds the closest collision free control
   */
  DynamicWindow(const Collision& collision, double dt, double horizon, double acc_dt,
                double acc_lim_x, double acc_lim_y, double acc_lim_th, double max_vel_x,
                double min_vel_x, double max_vel_y, double min_vel_y, double max_rot_vel,
                double min_rot_vel, unsigned int vx_samples, unsigned int vy_samples,
                unsigned int vth_samples);

  /**
   * @brief Compose best control
   * @param grid - grid map
   * @param x0 - current state [x, y, theta]
   * @param vb - current twist [vx, vy, w]
   * @param vref - desired twist to follow [vx, vy, w]
   * @return true if at least 1 solution found and optimal twist [vx, vy, w]
   * @details u_opt is set to zeros if no solution found
   */
  tuple<bool, vec> control(const GridMap& grid, const vec& x0, const vec& vb,
                           const vec& vref) const;

  /**
   * @brief Compose best control
   * @param grid - grid map
   * @param x0 - current state [x, y, theta]
   * @param vb - current twist [vx, vy, w]
   * @param xt_ref - reference trajectory
   * @param dt_ref - reference trajectory time discretization
   * @return true if at least one solution is found and optimal twist [vx, vy, w]
   * @details u_opt is set to zeros if no solution found
   */
  tuple<bool, vec> control(const GridMap& grid, const vec& x0, const vec& vb,
                           const mat& xt_ref, double dt_ref) const;

  /** @brief return time step */
  double timeStep() const
  {
    return dt_;
  }

  /** @brief return control horizon */
  double horizon() const
  {
    return horizon_;
  }

  /** @brief return number of steps in horizon */
  unsigned int steps() const
  {
    return steps_;
  }

private:
  /**
   * @brief Compose window size and discretization
   * @param vb - current twist [vx, vy, w]
   * @return twist lower limits and twist discretization
   */
  tuple<vec, vec> window(const vec& vb) const;

  /**
   * @brief Objective function for reference twist
   * @param grid - grid map
   * @param x0 - current state [x, y, theta]
   * @param vref - desired twist to follow [vx, vy, w]
   * @param u - twist from dynamic window
   * @return cost of the current sampled twist
   * @details - Finds best twist to reference twist that is collision free. If there
   * is a collision the cost is the max numerical value for a double.
   */
  double objective(const GridMap& grid, const vec& x0, const vec& vref,
                   const vec& u) const;

  /**
   * @brief Objective function for reference trajectory
   * @param grid - grid map
   * @param x0 - current state [x, y, theta]
   * @param u - twist from dynamic window
   * @param xt_ref - reference trajectory
   * @param tf - reference trajectory length in time
   * @return cost of the current sampled twist
   * @details - Finds best twist that generates a path closest to the reference
   * trajectory that is collision free. If there is a collision the cost is the
   * max numerical value for a double.
   */
  double objective(const GridMap& grid, const vec& x0, const vec& u, const mat& xt_ref,
                   double tf) const;

private:
  Collision collision_;                        // collision detection
  double dt_;                                  // time step in integration
  double horizon_;                             // control horizon
  double acc_dt_;                              // time acceleration limit is applied
  double acc_lim_x_, acc_lim_y_, acc_lim_th_;  // acceleration limits
  double max_vel_x_, min_vel_x_;               // velocity limits in x-direction
  double max_vel_y_, min_vel_y_;               // velocity limits in y-direction
  double max_rot_vel_, min_rot_vel_;           // rotational velocity limits
  unsigned int vx_samples_;   // number of velocity samples in x-direction
  unsigned int vy_samples_;   // number of velocity samples in y-direction
  unsigned int vth_samples_;  // number of angular velocity samples
  unsigned int steps_;        // number of steps in each rollout
};
}  // namespace ergodic_exploration
#endif
