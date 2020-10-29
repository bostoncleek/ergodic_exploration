/**
 * @file integrator.hpp
 * @author Boston Cleek
 * @date 28 Oct 2020
 * @brief Numerical integration methods
 */

#include <functional>
#include <armadillo>

namespace ergodic_exploration
{
using arma::mat;
using arma::vec;

/**
 * @brief Function representing a dynamic system
 * @details xdot = f(x,u)
 */
typedef std::function<vec(const vec&, const vec&)> DynaFunc;

/** @brief 4th order Runge-Kutta integration */
class RungeKutta
{
public:
  /**
   * @brief Constructor
   * @param dt - time step
   */
  RungeKutta(double dt);

  /**
   * @brief Simulate the dynamics
   * @param func - function to integrate
   * @param x0 - initial state (column vector)
   * @param u - control signal (column vector)
   * @param horizon - length of trajectory in time
   * @return new state (column vector)
   */
  mat solve(const DynaFunc& func, const vec& x0, const mat& ut, double horizon);

  /**
   * @brief Performs one step of RK4
   * @param func - function to integrate
   * @param x - state (column vector)
   * @param u - control (column vector)
   * @return new state (column vector)
   */
  vec step(const DynaFunc& func, const vec& x, const vec& u);

private:
  double dt_;
};

}  // namespace ergodic_exploration
