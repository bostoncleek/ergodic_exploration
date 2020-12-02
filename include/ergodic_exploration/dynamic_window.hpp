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

#include <nav_msgs/Path.h>
#include <tf2/LinearMath/Quaternion.h>

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

  // TODO: return bool indicating failure so recovery behviors can be used
  /**
   * @brief Compose best control
   * @param grid - grid map
   * @param x0 - current state [x, y, theta]
   * @param vb - current twist [vx, vy, w]
   * @param vref - desired twist to follow [vx, vy, w]
   * @return optimal twist [vx, vy, w]
   * @details u_opt is set to zeros if no solution found
   */
  vec control(const GridMap& grid, const vec& x0, const vec& vb, const vec& vref) const;

  vec control(const GridMap& grid, const vec& x0, const vec& vb, const mat& xt_ref,
              double delta) const;

  /**
   * @brief Objective function to minimize
   * @param cost[out] - cost of objective function
   * @param grid - grid map
   * @param x0 - current state [x, y, theta]
   * @param vref - desired twist to follow [vx, vy, w]
   * @param u - twist from dynamic window
   * @return true if the twist is collision free
   */
  bool objective(double& cost, const GridMap& grid, const vec& x0, const vec& vref,
                 const vec& u) const;

  bool objective(double& cost, const GridMap& grid, const vec& x0, const vec& u,
                 const mat& xt_ref, double tf) const;

  /**
   * @brief Visulaize path from following twist
   * @param path - trajectory
   * @param x - current state
   * @param u - twist [vx, vy, w]
   * @param frame - trajectory frame
   */
  void path(nav_msgs::Path& path, const vec& x, const vec& u, std::string frame) const;

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
