/**
 * @file types.hpp
 * @author Boston Cleek
 * @date 4 Nov 2020
 * @brief Types
 */
#ifndef TYPES_HPP
#define TYPES_HPP

#include <iosfwd>

namespace ergodic_exploration
{
/** @ brief A 2-Dimensional twist */
struct Twist2D
{
  /**
   * @brief Constructor
   * @param vx - velocity in x-direction
   * @param vy - velocity in y-direction
   * @param w - angular velocity about z-axis
   */
  Twist2D(double vx, double vy, double w) : vx(vx), vy(vy), w(w)
  {
  }

  /** @brief Constructor */
  Twist2D() : vx(0.0), vy(0.0), w(0.0)
  {
  }

  double vx;  // linear x velocity
  double vy;  // linear y velocity
  double w;   // rotation about z-axis
};

/**
 * @brief Output a 2 dimensional twisy as [vx vy w]
 * @param os[out] - stream to output to
 * @param twist - the twist to display
 * @return output steam
 */
std::ostream& operator<<(std::ostream& os, const Twist2D& twist);
}  // namespace ergodic_exploration
#endif
