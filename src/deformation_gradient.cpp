#include "deformation_gradient.h"

void deformation_gradient(
  const Eigen::MatrixXd & V,
  const Eigen::RowVector4i & tet,
  const Eigen::MatrixXd & U,
  Eigen::Matrix3d & F) {
  int v0 = tet(0);
  int v1 = tet(1);
  int v2 = tet(2);
  int v3 = tet(3);
  // get deformed shape matrix
  Eigen::MatrixXd D_def(3, 4);
  D_def.col(0) = V.row(v0) + U.row(v0);
  D_def.col(1) = V.row(v1) + U.row(v1);
  D_def.col(2) = V.row(v2) + U.row(v2);
  D_def.col(3) = V.row(v3) + U.row(v3);
  // get undeformed shape matrix stuff
  Eigen::Matrix3d Te = Eigen::Matrix3d::Zero();
  Te.col(0) = V.row(v1) - V.row(v0);
  Te.col(1) = V.row(v2) - V.row(v0);
  Te.col(2) = V.row(v3) - V.row(v0);
  Eigen::Matrix3d Te_inv = Te.inverse();
  Eigen::Matrix43d dphi;
  dphi.row(0) = -Eigen::RowVector3d::Ones() * Te_inv;
  dphi.block<3, 3>(1, 0) = Te_inv;
  // get deformation gradient
  F = D_def * dphi;
}