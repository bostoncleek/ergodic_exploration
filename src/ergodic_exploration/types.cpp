/**
 * @file types.cpp
 * @author Boston Cleek
 * @date 4 Nov 2020
 * @brief Types
 */

#include <iostream>

#include <ergodic_exploration/types.hpp>

namespace ergodic_exploration
{
std::ostream& operator<<(std::ostream& os, const Twist2D& twist)
{
  os << "[" << twist.vx << " " << twist.vy << " " << twist.w << "]\n";
  return os;
}
}  // namespace ergodic_exploration
