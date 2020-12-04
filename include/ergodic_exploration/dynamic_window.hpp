/**
 * @file dynamic_window.hpp
 * @author Boston Cleek
 * @date 5 Nov 2020
 * @brief Dynamic window control
 */
#ifndef DYNAMIC_WINDOW_HPP
#define DYNAMIC_WINDOW_HPP

#include <cmath>
#include <limits>

#include <ergodic_exploration/grid.hpp>
#include <ergodic_exploration/collision.hpp>

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
   * @param u_opt[out] - optimal twist [vx, vy, w]
   * @param grid - grid map
   * @param x0 - current state [x, y, theta]
   * @param vb - current twist [vx, vy, w]
   * @param vref - desired twist to follow [vx, vy, w]
   * @return true if at least one solution is found
   * @details u_opt is set to zeros if no solution found
   */
  bool control(vec& u_opt, const GridMap& grid, const vec& x0, const vec& vb,
               const vec& vref) const;

  /**
   * @brief Compose best control
   * @param u_opt[out] - optimal twist [vx, vy, w]
   * @param grid - grid map
   * @param x0 - current state [x, y, theta]
   * @param vb - current twist [vx, vy, w]
   * @param xt_ref - reference trajectory
   * @param dt_ref - reference trajectory time discretization
   * @return true if at least one solution is found
   * @details u_opt is set to zeros if no solution found
   */
  bool control(vec& u_opt, const GridMap& grid, const vec& x0, const vec& vb,
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
   * @param vel_lower[out] - twist lower limits
   * @param delta_vb[out] - twist discretization
   * @param vb - current twist [vx, vy, w]
   */
  void window(vec& vel_lower, vec& delta_vb, const vec& vb) const;

  /**
   * @brief Objective function
   * @param cost[out] - cost of objective function
   * @param grid - grid map
   * @param x0 - current state [x, y, theta]
   * @param vref - desired twist to follow [vx, vy, w]
   * @param u - twist from dynamic window
   * @return true if the twist is collision free
   * @details - finds best twist to reference twist that is collision free
   */
  bool objective(double& cost, const GridMap& grid, const vec& x0, const vec& vref,
                 const vec& u) const;

  /**
   * @brief Objective function
   * @param cost[out] - cost of objective function
   * @param grid - grid map
   * @param x0 - current state [x, y, theta]
   * @param u - twist from dynamic window
   * @param xt_ref - reference trajectory
   * @param tf - reference trajectory length in time
   * @return true if the twist is collision free
   * @details - finds best twist that generates a path closest to the reference
   * trajectory that is collision free
   */
  bool objective(double& cost, const GridMap& grid, const vec& x0, const vec& u,
                 const mat& xt_ref, double tf) const;

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
