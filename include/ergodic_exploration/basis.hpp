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

/**
 * @brief Fourier cosine basis */
class Basis
{
private:
  template <class ModelT>
  friend class ErgodicControl;

public:
  /**
   * @brief Constructor
   * @param lx - x-axis lenght of domain
   * @param ly - y-axis lenght of domain
   * @param num_basis - number of each basis functions per dimension
   */
  Basis(double lx, double ly, unsigned int num_basis);

  /**
   * @brief Compose cosine basis functions given positon
   * @param x - robot position [x y]
   * @return cosine basis functions (num_basis^2 x 1)
   * @details x must be on domain [0 lx] x [0 ly]
   */
  vec fourierBasis(const vec& x) const;

  /**
   * @brief Compose gradient cosine basis functions
   * @param x - robot position [x y]
   * @return gradient of each basis function (2 x num_basis^2)
   * @details x must be on domain [0 lx] x [0 ly]
   */
  mat gradFourierBasis(const vec& x) const;

  /**
   * @brief Compose trajectory fourier coefficients
   * @param xt - trajectory
   * @return trajectory fourier coefficients (num_basis^2 x 1)
   * @details xt must be on domain [0 lx] x [0 ly] and robot heading is not required
   */
  vec trajCoeff(const mat& xt) const;

  /**
   * @brief Compose spatial fourier coefficients
   * @param phi_vals - target evaluated at each grid cell in phi_grid
   * @param phi_grid - discretization of fourier domain
   * @return spatial fourier coefficients (num_basis^2 x 1)
   * @details phi_grid is not the occupancy grid, the first row contains the
   * x-cordinates and the second row contains the corresponding y-coordinates.
   * The index of the elements in phi_vals must correspond to the columns in phi_grid.
   */
  vec spatialCoeff(const vec& phi_vals, const mat& phi_grid) const;

private:
  double lx_, ly_;            // length of domain
  unsigned int total_basis_;  // number of basis functions
  vec lamdak_;                // frequency coefficients weights
  imat k_;                    // basis number
};
}  // namespace ergodic_exploration
#endif
