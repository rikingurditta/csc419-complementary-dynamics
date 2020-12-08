#include <Eigen/Core>

#ifndef COMPLEMENTARY_DYNAMICS_RIG_JACOBIAN_H
#define COMPLEMENTARY_DYNAMICS_RIG_JACOBIAN_H

#endif //COMPLEMENTARY_DYNAMICS_RIG_JACOBIAN_H

// Inputs:
//   V  #V by 3 list of rest (initial pose) mesh vertex positions
//   W  #V by k matrix of weights, where W(i, j) is how much effect bone j has on vertex i
// Outputs:
//   J  #V*3 by k*12 rig Jacobian of linear blend skinning based rig
void rig_jacobian(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXd & W,
  Eigen::MatrixXd & J);