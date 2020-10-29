/**
 * @file integrator.cpp
 * @author Boston Cleek
 * @date 28 Oct 2020
 * @brief Numerical integration methods
 */

#include <ergodic_exploration/integrator.hpp>

namespace ergodic_exploration
{
RungeKutta::RungeKutta(double dt) : dt_(dt)
{
}

mat RungeKutta::solve(const DynaFunc& func, const vec& x0, const mat& ut, double horizon)
{
  vec x = x0;
  const auto steps = static_cast<unsigned int>(horizon / dt_);
  mat xt(x0.n_rows, steps);

  for (unsigned int i = 0; i < steps; i++)
  {
    x = step(func, x, ut.col(i));
    xt.col(i) = x;
  }
  return xt;
}

vec RungeKutta::step(const DynaFunc& func, const vec& x, const vec& u)
{
  const vec k1 = func(x, u);
  const vec k2 = func(x + dt_ * (0.5 * k1), u);
  const vec k3 = func(x + dt_ * (0.5 * k2), u);
  const vec k4 = func(x + dt_ * k3, u);
  return x + dt_ / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

}  // namespace ergodic_exploration
