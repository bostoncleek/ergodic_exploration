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
 * @file basis.cpp
 * @author Boston Cleek
 * @date 20 Nov 2020
 * @brief Fourier cosine basis
 */

#include <cmath>

#include <ergodic_exploration/basis.hpp>
#include <ergodic_exploration/numerics.hpp>

namespace ergodic_exploration
{
using arma::span;

Basis::Basis(double lx, double ly, unsigned int num_basis)
  : lx_(lx)
  , ly_(ly)
  , total_basis_(num_basis * num_basis)
  , lamdak_(total_basis_)
  , k_(2, total_basis_)
{
  // TODO: add hk

  // basis numbers
  auto col = 0;
  for (unsigned int i = 0; i < num_basis; i++)
  {
    for (unsigned int j = 0; j < num_basis; j++)
    {
      k_(0, col) = j;
      k_(1, col) = i;
      col++;
    }
  }
  // k_.t().print("k");

  // TODO: verify power on demoninator
  // frequency coefficients weights
  for (unsigned int i = 0; i < total_basis_; i++)
  {
    lamdak_(i) = 1.0 / std::pow((1.0 + std::sqrt(sum(square(k_.col(i))))), 1.5);
  }
  // lamdak_.print("lamdak");
}

vec Basis::fourierBasis(const vec& x) const
{
  vec fk(total_basis_);
  for (unsigned int i = 0; i < total_basis_; i++)
  {
    fk(i) =
        std::cos(k_(0, i) * (PI / lx_) * x(0)) * std::cos(k_(1, i) * (PI / ly_) * x(1));
  }

  return fk;
}

mat Basis::gradFourierBasis(const vec& x) const
{
  mat dfk(2, total_basis_);
  auto k1 = 0.0;
  auto k2 = 0.0;

  for (unsigned int i = 0; i < total_basis_; i++)
  {
    k1 = k_(0, i) * (PI / lx_);
    k2 = k_(1, i) * (PI / ly_);

    dfk(0, i) = -k1 * std::sin(k1 * x(0)) * std::cos(k2 * x(1));
    dfk(1, i) = -k2 * std::cos(k1 * x(0)) * std::sin(k2 * x(1));
  }

  return dfk;
}

vec Basis::trajCoeff(const mat& xt) const
{
  mat fk_mat(total_basis_, xt.n_cols);

  for (unsigned int i = 0; i < xt.n_cols; i++)
  {
    fk_mat.col(i) = fourierBasis(xt(span(0, 1), span(i, i)));
  }

  // sum of each row
  return (1.0 / xt.n_cols) * sum(fk_mat, 1);
}

vec Basis::spatialCoeff(const vec& phi_vals, const mat& phi_grid) const
{
  mat fk_mat(total_basis_, phi_grid.n_cols);

  for (unsigned int i = 0; i < phi_grid.n_cols; i++)
  {
    fk_mat.col(i) = fourierBasis(phi_grid.col(i)) * phi_vals(i);
  }

  // sum of each row
  return sum(fk_mat, 1);
}
}  // namespace ergodic_exploration
