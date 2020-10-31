/**
 * @file integrator.hpp
 * @author Boston Cleek
 * @date 28 Oct 2020
 * @brief Numerical integration methods
 */

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
// typedef std::function<vec(const vec&, const vec&)> DynaFunc;

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
  template <class T>
  mat solve(const T& model, const vec& x0, const mat& ut, double horizon);

  /**
   * @brief Performs one step of RK4
   * @param func - function to integrate
   * @param x - state (column vector)
   * @param u - control (column vector)
   * @return new state (column vector)
   */
  template <class T>
  vec step(const T& model, const vec& x, const vec& u);

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
  // template <class T>
  // vec step(const T& model, const Collision& collision, const mat& xt, const mat& ut, const vec& edx);

private:
  double dt_;
};

RungeKutta::RungeKutta(double dt) : dt_(dt)
{
}

template <class T>
mat RungeKutta::solve(const T& model, const vec& x0, const mat& ut, double horizon)
{
  vec x = x0;
  const auto steps = static_cast<unsigned int>(horizon / dt_);
  mat xt(x0.n_rows, steps);

  for (unsigned int i = 0; i < steps; i++)
  {
    x = step(model, x, ut.col(i));
    xt.col(i) = x;
  }
  return xt;
}

template <class T>
vec RungeKutta::step(const T& model, const vec& x, const vec& u)
{
  const vec k1 = model(x, u);
  const vec k2 = model(x + dt_ * (0.5 * k1), u);
  const vec k3 = model(x + dt_ * (0.5 * k2), u);
  const vec k4 = model(x + dt_ * k3, u);
  return x + dt_ / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

// template <class T>
// vec RungeKutta::step(const T& model, const Collision& collision, const mat& xt,
//                      const mat& ut, const vec& edx)
// {
// }

}  // namespace ergodic_exploration
