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
#include <ergodic_exploration/numerics.hpp>

namespace ergodic_exploration
{
using arma::linspace;
using arma::mat;
using arma::span;
using arma::vec;

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
   * @param model - dynamic model
   * @param x0 - initial state
   * @param ut - control signal (each column is applied at a single time step)
   * @param horizon - length of trajectory in time
   * @return trajectory
   * @details the boundary condition is not added to the trajectory
   */
  template <class ModelT>
  mat solve(const ModelT& model, const vec& x0, const mat& ut, double horizon) const;

  /**
   * @brief Solve the co-state variable backwards in time
   * @param func - time derivatve of co-state variable
   * @param model - dynamic model
   * @param rhoT - co-state variable terminal condition (zero vector)
   * @param xt - forward porpagated dynamic model trajectory
   * @param ut - control signal
   * @param edx - gradient of the ergodic metric for each state in xt
   * @param bdx - derivatve of barrier function for each state in xt
   * @param horizon - length of trajectory in time
   * @return co-state variable solution
   * @details co-state is sorted from [t0 tf] no need to index backwards and the
   * boundary condition is not added to the trajectory
   */
  template <class ModelT>
  mat solve(const CoStateFunc& func, const ModelT& model, const vec& rhoT, const mat& xt,
            const mat& ut, const mat& edx, const mat& bdx, double horizon) const;

  /**
   * @brief Performs one step of RK4 forward in time
   * @param model - dynamic model
   * @param x - state
   * @param u - control
   * @return new state
   */
  template <class ModelT>
  vec step(const ModelT& model, const vec& x, const vec& u) const;

  /**
   * @brief Performs one step of RK4 backwards in time
   * @param func - time derivatve of the co-state variable
   * @param rho - co-state variable
   * @param gdx - gradient of the ergodic metric
   * @param dbar - derivatve of barrier function for a state
   * @param fdx - jacobian of the model with respect to the control
   * @return co-state variable
   * @details The robot model is used to compose A = D1[f(x,u)].
   * The columns of xt, ut, and edx correspond to
   * the state, control, or derivative at a given time.
   */
  vec step(const CoStateFunc& func, const vec& rho, const vec& gdx, const vec& dbar,
           const mat& fdx) const;

private:
  double dt_;  // time step
};

RungeKutta::RungeKutta(double dt) : dt_(dt)
{
}

template <class ModelT>
mat RungeKutta::solve(const ModelT& model, const vec& x0, const mat& ut,
                      double horizon) const
{
  // TODO: Add terminal x0?
  vec x = x0;
  const auto steps = static_cast<unsigned int>(std::abs(horizon / dt_));
  mat xt(x.n_rows, steps);

  for (unsigned int i = 0; i < steps; i++)
  {
    x = step(model, x, ut.col(i));
    x(2) = normalize_angle_PI(x(2));
    xt.col(i) = x;
  }

  return xt;
}

template <class ModelT>
mat RungeKutta::solve(const CoStateFunc& func, const ModelT& model, const vec& rhoT,
                      const mat& xt, const mat& ut, const mat& edx, const mat& bdx,
                      double horizon) const
{
  // TODO: Add terminal p(T)?
  vec rho = rhoT;
  const auto steps = static_cast<unsigned int>(std::abs(horizon / dt_));
  mat rhot(rho.n_rows, steps);

  // Iterate backwards
  // this way rhot from t0 to tf in the returned matrix
  for (unsigned int i = steps; i-- > 0;)
  {
    // Index states, controls, and ergodic measures at end of array
    rho = step(func, rho, edx.col(i), bdx.col(i), model.fdx(xt.col(i), ut.col(i)));
    rhot.col(i) = rho;
  }

  return rhot;
}

template <class ModelT>
vec RungeKutta::step(const ModelT& model, const vec& x, const vec& u) const
{
  const vec k1 = model(x, u);
  const vec k2 = model(x + dt_ * (0.5 * k1), u);
  const vec k3 = model(x + dt_ * (0.5 * k2), u);
  const vec k4 = model(x + dt_ * k3, u);
  return x + (dt_ / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

vec RungeKutta::step(const CoStateFunc& func, const vec& rho, const vec& gdx,
                     const vec& dbar, const mat& fdx) const
{
  const vec k1 = func(rho, gdx, dbar, fdx);
  const vec k2 = func(rho - dt_ * (0.5 * k1), gdx, dbar, fdx);
  const vec k3 = func(rho - dt_ * (0.5 * k2), gdx, dbar, fdx);
  const vec k4 = func(rho - dt_ * k3, gdx, dbar, fdx);
  return rho - dt_ / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

class RungeKutta45
{
public:
  RungeKutta45(double hmax, double hmin, double epsilon, unsigned int max_iter);

  template <class ModelT>
  bool solve(mat& xt, const ModelT& model, const vec& x0, const mat& ut, double dt,
             double horizon) const;

  template <class ModelT>
  bool solve(mat& rhot, const CoStateFunc& func, const ModelT& model, const vec& rhoT,
             const mat& xt, const mat& ut, const mat& edx, const mat& bdx, double dt,
             double horizon) const;

  template <class ModelT>
  double step(vec& x_new, const ModelT& model, const vec& x, const vec& u, double h) const;

  double step(vec& rho_new, const CoStateFunc& func, const vec& rho, const vec& gdx,
              const vec& dbar, const mat& fdx, double h) const;

private:
  double hmin_, hmax_;
  double epsilon_;
  unsigned int max_iter_;
};

RungeKutta45::RungeKutta45(double hmin, double hmax, double epsilon, unsigned int max_iter)
  : hmin_(hmin), hmax_(hmax), epsilon_(epsilon), max_iter_(max_iter)
{
}

template <class ModelT>
bool RungeKutta45::solve(mat& xt, const ModelT& model, const vec& x0, const mat& ut,
                         double dt, double horizon) const
{
  // TODO: how to initialize h?
  auto h = (hmax_ + hmin_) / 2.0;
  // auto h = hmin_;
  vec x = x0;

  // desired length of trajectory
  const auto steps = static_cast<unsigned int>(std::abs(horizon / dt));
  const vec tvec = linspace(0.0, horizon, steps);

  // allocate memory for trajectory based on max possible steps
  const auto max_steps = static_cast<unsigned int>(horizon / hmin_) + 1;
  mat xt_max(x0.n_rows, max_steps);
  vec tvec_max(max_steps);

  // add x0
  xt_max.col(0) = x0;
  tvec_max(0) = 0.0;

  vec x_new;
  auto t = 0.0;
  unsigned int stored = 1;
  unsigned int iter = 0;
  bool flag = 0;
  while (!flag)
  {
    const auto i = static_cast<unsigned int>(
        std::floor(std::round(static_cast<double>(steps - 1) * std::abs(t / horizon))));
    // std::cout << "i: " << i << std::endl;

    const auto r = step(x_new, model, x, ut.col(i), h);
    // std::cout << "r: " << r << std::endl;

    if (r < epsilon_ || almost_equal(r, epsilon_))
    {
      x = x_new;
      xt_max.col(stored) = x;

      t += h;
      tvec_max(stored) = t;

      stored++;
    }

    h *= 0.84 * std::pow(epsilon_ / r, 0.25);
    h = std::clamp(h, hmin_, hmax_);

    if (t > horizon || almost_equal(t, horizon))
    {
      flag = 1;
    }

    // detect final time
    else if (t + h > horizon)
    {
      h = horizon - t;
    }

    // std::cout << "t: " << t << std::endl;
    // std::cout << "h: " << h << std::endl;

    if (iter == max_iter_)
    {
      std::cout << "WARNING: max iterations reached " << std::endl;
      return false;
    }

    iter++;
  }

  // xt_max.cols(0, stored-1).t().print("traj");
  // tvec_max.rows(0, stored - 1).print("t rk45");
  // tvec.print("t desired");
  //
  // std::cout << "tvec_max: " << tvec_max.n_rows << std::endl;
  // std::cout << "tvec: " << tvec.n_rows << std::endl;
  //
  // std::cout << "steps: " << steps << std::endl;
  // std::cout << "stored: " << stored << std::endl;

  // fit quartic polynomials
  xt.resize(x0.n_rows, steps);
  for (unsigned int i = 0; i < x0.n_rows; i++)
  {
    const vec p =
        polyfit(tvec_max.rows(0, stored - 1), xt_max(span(i, i), span(0, stored - 1)), 4);
    // p.print("p");
    xt.row(i) = polyval(p, tvec).t();
  }

  return true;
}

template <class ModelT>
double RungeKutta45::step(vec& x_new, const ModelT& model, const vec& x, const vec& u,
                          double h) const
{
  const vec k1 = h * model(x, u);

  const vec k2 = h * model(x + ((1.0 / 4.0) * k1), u);

  const vec k3 = h * model(x + ((3.0 / 32.0) * k1) + ((9.0 / 32.0) * k2), u);

  const vec k4 = h * model(x + ((1932.0 / 2197.0) * k1) - ((7200.0 / 2197.0) * k2) +
                               ((7296.0 / 2197.0) * k3),
                           u);
  const vec k5 = h * model(x + ((439.0 / 216.0) * k1) - (8.0 * k2) +
                               ((3680.0 / 513.0) * k3) - ((845.0 / 4104.0) * k4),
                           u);
  const vec k6 =
      h * model(x - ((8.0 / 27.0) * k1) + (2.0 * k2) - ((3544.0 / 2565.0) * k3) +
                    ((1859.0 / 4104.0) * k4) - ((11.0 / 40.0) * k5),
                u);

  x_new = x + ((25.0 / 216.0) * k1) + ((1408.0 / 2565.0) * k3) +
          ((2197.0 / 4101.0) * k4) - ((1.0 / 5.0) * k5);

  const vec z = x + ((16.0 / 135.0) * k1) + ((6656.0 / 12825.0) * k3) +
                ((28561.0 / 56430.0) * k4) - ((9.0 / 50.0) * k5) + ((2.0 / 55.0) * k6);

  // std::cout << "max diff: " << max(abs(z - x_new)) << std::endl;
  // std::cout << "error norm: " << norm(z - x_new, 2) << std::endl;

  return max(abs(z - x_new)) / h;
  // return min(abs(z - x_new)) / h;
}

template <class ModelT>
bool RungeKutta45::solve(mat& rhot, const CoStateFunc& func, const ModelT& model,
                         const vec& rhoT, const mat& xt, const mat& ut, const mat& edx,
                         const mat& bdx, double dt, double horizon) const
{
  // TODO: how to initialize h?
  // auto h = (hmax_ + hmin_) / 2.0;
  auto h = hmin_;
  vec rho = rhoT;

  // desired length of trajectory
  const auto steps = static_cast<unsigned int>(std::abs(horizon / dt));
  const vec tvec = linspace(0.0, horizon, steps);

  // allocate memory for trajectory based on max possible steps
  const auto max_steps = static_cast<unsigned int>(horizon / hmin_) + 1;
  mat rhot_max(rhoT.n_rows, max_steps);
  vec tvec_max(max_steps);

  // add x0
  rhot_max.col(0) = rhoT;
  tvec_max(0) = horizon;

  vec rho_new;
  auto t = horizon;
  unsigned int stored = 1;
  unsigned int iter = 0;
  bool flag = 0;
  while (!flag)
  {
    const auto i = static_cast<unsigned int>(
        std::floor(std::round(static_cast<double>(steps - 1) * std::abs(t / horizon))));
    // std::cout << "i: " << i << std::endl;

    const auto r = step(rho_new, func, rho, edx.col(i), bdx.col(i),
                        model.fdx(xt.col(i), ut.col(i)), h);
    // std::cout << "r: " << r << std::endl;

    if (r < epsilon_ || almost_equal(r, epsilon_))
    {
      rho = rho_new;
      rhot_max.col(stored) = rho;

      t -= h;
      tvec_max(stored) = t;

      stored++;
    }

    h *= 0.84 * std::pow(epsilon_ / r, 0.25);
    h = std::clamp(h, hmin_, hmax_);

    if (t < 0.0 || almost_equal(t, 0.0))
    {
      flag = 1;
    }

    // detect final time
    else if (t - h < 0.0)
    {
      h = t;
    }

    // std::cout << "t: " << t << std::endl;
    // std::cout << "h: " << h << std::endl;

    if (iter == max_iter_)
    {
      std::cout << "WARNING: max iterations reached " << std::endl;
      return false;
    }

    iter++;
  }

  // rhot_max.cols(0, stored-1).t().print("rho traj");
  // tvec_max.rows(0, stored - 1).print("t rk45");
  // tvec.print("t desired");
  //
  // std::cout << "tvec_max: " << tvec_max.n_rows << std::endl;
  // std::cout << "tvec: " << tvec.n_rows << std::endl;
  //
  // std::cout << "steps: " << steps << std::endl;
  // std::cout << "stored: " << stored << std::endl;

  // fit quartic polynomials
  rhot.resize(rho.n_rows, steps);
  for (unsigned int i = 0; i < rhoT.n_rows; i++)
  {
    const vec p = polyfit(tvec_max.rows(0, stored - 1),
                          rhot_max(span(i, i), span(0, stored - 1)), 4);
    // p.print("p");
    rhot.row(i) = polyval(p, tvec).t();
  }

  return true;
}

double RungeKutta45::step(vec& rho_new, const CoStateFunc& func, const vec& rho,
                          const vec& gdx, const vec& dbar, const mat& fdx, double h) const
{
  const vec k1 = -h * func(rho, gdx, dbar, fdx);

  const vec k2 = -h * func(rho + ((1.0 / 4.0) * k1), gdx, dbar, fdx);

  const vec k3 =
      -h * func(rho + ((3.0 / 32.0) * k1) + ((9.0 / 32.0) * k2), gdx, dbar, fdx);

  const vec k4 = -h * func(rho + ((1932.0 / 2197.0) * k1) - ((7200.0 / 2197.0) * k2) +
                               ((7296.0 / 2197.0) * k3),
                           gdx, dbar, fdx);

  const vec k5 = -h * func(rho + ((439.0 / 216.0) * k1) - (8.0 * k2) +
                               ((3680.0 / 513.0) * k3) - ((845.0 / 4104.0) * k4),
                           gdx, dbar, fdx);

  const vec k6 =
      -h * func(rho - ((8.0 / 27.0) * k1) + (2.0 * k2) - ((3544.0 / 2565.0) * k3) +
                    ((1859.0 / 4104.0) * k4) - ((11.0 / 40.0) * k5),
                gdx, dbar, fdx);

  rho_new = rho + ((25.0 / 216.0) * k1) + ((1408.0 / 2565.0) * k3) +
            ((2197.0 / 4101.0) * k4) - ((1.0 / 5.0) * k5);

  const vec z = rho + ((16.0 / 135.0) * k1) + ((6656.0 / 12825.0) * k3) +
                ((28561.0 / 56430.0) * k4) - ((9.0 / 50.0) * k5) + ((2.0 / 55.0) * k6);

  // std::cout << "max diff: " << max(abs(z - rho_new)) << std::endl;
  // std::cout << "error norm: " << norm(z - rho_new, 2) << std::endl;

  return max(abs(z - rho_new)) / h;
}

}  // namespace ergodic_exploration
#endif
