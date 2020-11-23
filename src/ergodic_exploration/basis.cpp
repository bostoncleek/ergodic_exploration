/**
 * @file basis.cpp
 * @author Boston Cleek
 * @date 20 Nov 2020
 * @brief Fourier cosine basis
 */

#include <iostream>
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

  // TODO: verify norm on k
  // TODO: verify power on demoninator
  // frequency coefficients weights}
  for (unsigned int i = 0; i < total_basis_; i++)
  {
    lamdak_(i) = 1.0 / std::pow((1.0 + std::sqrt(sum(square(k_.col(i))))), 1.5);
  }
  // lamdak_.print("lamdak");
}

void Basis::fourierBasis(vec& fk, const vec& x)
{
  fk.resize(total_basis_);
  for (unsigned int i = 0; i < total_basis_; i++)
  {
    fk(i) =
        std::cos(k_(0, i) * (PI / lx_) * x(0)) * std::cos(k_(1, i) * (PI / ly_) * x(1));
  }
}

void Basis::gradFourierBasis(mat& dfk, const vec& x)
{
  dfk.resize(2, total_basis_);
  auto k1 = 0.0;
  auto k2 = 0.0;

  for (unsigned int i = 0; i < total_basis_; i++)
  {
    k1 = k_(0, i) * (PI / lx_);
    k2 = k_(1, i) * (PI / ly_);

    dfk(0, i) = -k1 * std::sin(k1 * x(0)) * std::cos(k2 * x(1));
    dfk(1, i) = -k2 * std::cos(k1 * x(0)) * std::sin(k2 * x(1));
  }
}

void Basis::trajCoeff(vec& ck, const mat& xt)
{
  ck.resize(total_basis_);
  mat fk_mat(total_basis_, xt.n_cols);

  vec fk;
  for (unsigned int i = 0; i < xt.n_cols; i++)
  {
    fourierBasis(fk, xt(span(0, 1), span(i, i)));
    fk_mat.col(i) = fk;
  }

  // sum accross each column
  ck = (1.0 / xt.n_cols) * sum(fk_mat, 1);
}

void Basis::spatialCoeff(vec& phik, const vec& phi_vals, const mat& phi_grid)
{
  phik.resize(total_basis_);
  mat fk_mat(total_basis_, phi_grid.n_cols);

  vec fk;
  for (unsigned int i = 0; i < phi_grid.n_cols; i++)
  {
    fourierBasis(fk, phi_grid.col(i));
    fk_mat.col(i) = fk * phi_vals(i);
  }

  // sum accross each column
  phik = sum(fk_mat, 1);
}

}  // namespace ergodic_exploration
