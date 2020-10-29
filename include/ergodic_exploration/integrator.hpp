/**
 * @file integrator.hpp
 * @author Boston Cleek
 * @date 28 Oct 2020
 * @brief Numerical integration methods
 */

#include <functional>
#include <armadillo>

#include <ergodic_exploration/collision.hpp>

namespace ergodic_exploration
{
using arma::mat;
using arma::vec;

/**
 * @brief Function representing a dynamic system
 * @details xdot = f(x,u)
 */
typedef std::function<vec(const vec&, const vec&)> DynaFunc;

// template<class T>
// using CoState = std::function<void(T)>;

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

  /**
   * @brief Performs one step of RK4
   * @param model - robot model
   * @param collision - collision detector
   * @param xt - forward propagated trajectory
   * @param ut - control signal
   * @param edx - derivatve of the ergodic measure approximated by KL divergence
   * @return co-state variable
   * @details The robot model is used to compose A = D1[f(x,u)]. The collision
   * detector composes the boundary derivatives (db/dx), collision derivatives (dc/dx),
   * and range finder collisions (dr/dx).
   */
  template <class T>
  vec step(const T& model, const Collision& collision, const mat& xt, const mat& ut, const vec& edx);

private:
  double dt_;
};

template <class T>
vec RungeKutta::step(const T& model, const Collision& collision, const mat& xt,
                     const mat& ut, const vec& edx)
{
}

}  // namespace ergodic_exploration
