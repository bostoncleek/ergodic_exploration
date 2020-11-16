/**
 * @file integrator.hpp
 * @author Boston Cleek
 * @date 28 Oct 2020
 * @brief Numerical integration methods
 */
#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include <cmath>
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
// typedef std::function<vec(const vec&, const vec&)> DynaFunc;

/**
 * @brief Function representing the time derivatve of the co-state variable
 * @details inputs are co-state, ergodic measure derivatve, barrier derivatve,
 * and the jacobian of the dynamics w.r.t state
 */
typedef std::function<vec(const vec&, const vec&, const vec&, const mat&)> CoStateFunc;

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
   * @brief Simulate the dynamics forward in time
   * @param xt[out] - trajectory
   * @param model - dynamic model
   * @param x0 - initial state (column vector)
   * @param ut - control signal (each column is applied at a single time step)
   * @param horizon - length of trajectory in time
   */
  template <class ModelT>
  void solve(mat& xt, const ModelT& model, const vec& x0, const mat& ut, double horizon);

  /**
   * @brief Performs one step of RK4 forward in time
   * @param model - dynamic model
   * @param x - state (column vector)
   * @param u - control (column vector)
   * @return new state (column vector)
   */
  template <class ModelT>
  vec step(const ModelT& model, const vec& x, const vec& u);

  /**
   * @brief Solve the co-state variable backwards in time
   * @param rhot[out] - co-state variable solution
   * @param func - time derivatve of co-state variable
   * @param model - dynamic model
   * @param rhoT - co-state variable terminal condition (zero vector)
   * @param xt - forward porpagated dynamic model trajectory
   * @param ut - control signal
   * @param edx - derivatve of the ergodic measure for each state in xt
   * @param bdx - derivatve of barrier function for each state in xt
   * @param horizon - length of trajectory in time
   * @details co-state is sorted from [t0 tf] no need to index backwards
   */
  template <class ModelT>
  void solve(mat& rhot, const CoStateFunc& func, const ModelT& model, const vec& rhoT,
             const mat& xt, const mat& ut, const mat& edx, const mat& bdx,
             double horizon);

  /**
   * @brief Performs one step of RK4 backwards in time
   * @param func - time derivatve of the co-state variable
   * @param rho - co-state variable
   * @param kldx - derivatve of the ergodic measure approximated by KL divergence
   * @param dbar - derivatve of barrier function for a state
   * @param fdx - jacobian of the model with respect to the control
   * @return co-state variable
   * @details The robot model is used to compose A = D1[f(x,u)]. The collision
   * detector composes the boundary derivatives (db/dx), collision derivatives (dc/dx),
   * and range finder collisions (dr/dx). The columns of xt, ut, and edx correspond to
   * the state, control, or derivative at a given time.
   */
  vec step(const CoStateFunc& func, const vec& rho, const vec& kldx, const vec& dbar,
           const mat& fdx);

private:
  double dt_;
};

RungeKutta::RungeKutta(double dt) : dt_(dt)
{
}

template <class ModelT>
void RungeKutta::solve(mat& xt, const ModelT& model, const vec& x0, const mat& ut,
                       double horizon)
{
  vec x = x0;
  const auto steps = static_cast<unsigned int>(horizon / std::abs(dt_));
  xt.resize(x.n_rows, steps);

  for (unsigned int i = 0; i < steps; i++)
  {
    x = step(model, x, ut.col(i));
    xt.col(i) = x;
  }
}

template <class ModelT>
vec RungeKutta::step(const ModelT& model, const vec& x, const vec& u)
{
  const vec k1 = model(x, u);
  const vec k2 = model(x + dt_ * (0.5 * k1), u);
  const vec k3 = model(x + dt_ * (0.5 * k2), u);
  const vec k4 = model(x + dt_ * k3, u);
  return x + (dt_ / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

template <class ModelT>
void RungeKutta::solve(mat& rhot, const CoStateFunc& func, const ModelT& model,
                       const vec& rhoT, const mat& xt, const mat& ut, const mat& edx,
                       const mat& bdx, double horizon)
{
  vec rho = rhoT;
  const auto steps = static_cast<unsigned int>(horizon / std::abs(dt_));
  rhot.resize(rho.n_rows, steps);

  // Iterate backwards
  // this way rhot from t0 to tf in the returned matrix
  for (unsigned int i = steps; i-- > 0;)
  {
    // Index states, controls, and ergodic measures at end of array
    rho = step(func, rho, edx.col(i), bdx.col(i), model.fdx(xt.col(i), ut.col(i)));
    rhot.col(i) = rho;
  }
}

vec RungeKutta::step(const CoStateFunc& func, const vec& rho, const vec& kldx,
                     const vec& dbar, const mat& fdx)
{
  const vec k1 = func(rho, kldx, dbar, fdx);
  const vec k2 = func(rho - dt_ * (0.5 * k1), kldx, dbar, fdx);
  const vec k3 = func(rho - dt_ * (0.5 * k2), kldx, dbar, fdx);
  const vec k4 = func(rho - dt_ * k3, kldx, dbar, fdx);
  return rho - dt_ / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}
}  // namespace ergodic_exploration
#endif
