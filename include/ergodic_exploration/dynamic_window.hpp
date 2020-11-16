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
class DynamicWindow
{
public:
  DynamicWindow(double dt, double horizon, double frequency, double acc_lim_x,
                double acc_lim_y, double acc_lim_th, double max_vel_x, double min_vel_x,
                double max_vel_y, double min_vel_y, double max_rot_vel,
                double min_rot_vel, unsigned int vx_samples, unsigned int vy_samples,
                unsigned int vth_samples);

  /**
   * @brief Update control
   * @param collision - collision detector
   * @param grid - grid map
   * @param x0 - current state [x, y, theta]
   * @param vb - current twist [vx, vy, w]
   * @param vref - desired twist to follow [vx, vy, w]
   * @return twist
   */
  vec control(const Collision& collision, const GridMap& grid, const vec& x0,
              const vec& vb, const vec& vref);

  bool objective(double& loss, const Collision& collision, const GridMap& grid, const vec& x,
                   const vec& vref, const vec& u);
                   
private:
  double dt_;                                  // time step in integration
  double horizon_;                             // control horizon
  double frequency_;                           // control loop frequency
  double acc_lim_x_, acc_lim_y_, acc_lim_th_;  // acceleration limits
  double max_vel_x_, min_vel_x_;               // velocity limits in x-direction
  double max_vel_y_, min_vel_y_;               // velocity limits in y-direction
  double max_rot_vel_, min_rot_vel_;           // rotational velocity limits
  unsigned int vx_samples_;   // number of velocity samples in x-direction
  unsigned int vy_samples_;   // number of velocity samples in y-direction
  unsigned int vth_samples_;  // number of angular velocity samples
};
}  // namespace ergodic_exploration
#endif
