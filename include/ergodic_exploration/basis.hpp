/**
 * @file basis.hpp
 * @author Boston Cleek
 * @date 20 Nov 2020
 * @brief Fourier cosine basis
 */
#ifndef BASIS_HPP
#define BASIS_HPP

#include <armadillo>

namespace ergodic_exploration
{
using arma::imat;
using arma::mat;
using arma::vec;

class Basis
{
private:
  template <class ModelT>
  friend class ErgodicControl;

public:
  Basis(double lx, double ly, unsigned int num_basis);

  void fourierBasis(vec& fk, const vec& x);

  void gradFourierBasis(mat& dfk, const vec& x);

  void trajCoeff(vec& ck, const mat& xt);

  void spatialCoeff(vec& phik, const vec& phi_vals, const mat& phi_grid);

private:
  double lx_, ly_;            // length of domain
  unsigned int total_basis_;  // number of basis functions
  vec lamdak_;                // frequency coefficients weights
  imat k_;                    // basis number
};

}  // namespace ergodic_exploration
#endif
