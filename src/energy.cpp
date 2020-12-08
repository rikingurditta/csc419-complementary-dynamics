#include "energy.h"

double energy(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXd &T,
  const Eigen::VectorXd &Ur,
  const Eigen::VectorXd &Uc,
  const double dt,
  const double neohookean_C,
  const double neohookean_D) {
  // sum potential energy for each individual tetrahedron to get potential of entire mesh
  double e = 0;
  for (int t = 0; t < T.rows(); t++) {
    int v0 = T(t, 0);
    int v1 = T(t, 1);
    int v2 = T(t, 2);
    int v3 = T(t, 3);
    // get deformed shape matrix
    Eigen::Matrix34d D_def;
    D_def.col(0) = Ur.row(v0) + Uc.row(v0);
    D_def.col(1) = Ur.row(v1) + Uc.row(v1);
    D_def.col(2) = Ur.row(v2) + Uc.row(v2);
    D_def.col(3) = Ur.row(v3) + Uc.row(v3);
    // get undeformed shape matrix
    Eigen::Matrix3d Te = Eigen::Matrix3d::Zero();
    Te.col(0) = V.row(v1) - V.row(v0);
    Te.col(1) = V.row(v2) - V.row(v0);
    Te.col(2) = V.row(v3) - V.row(v0);
    Eigen::Matrix3d Te_inv = T.inverse();
    Eigen::Matrix4d D_und;
    D_und.row(0) = -Eigen::RowVector3d::Ones() * Te_inv;
    D_und.block<3, 3>(1, 0) = Te_inv;
    // get deformation gradient
    Eigen::Matrix3d = D_def * D_und;
    // calculate energy for current tetrahedron using neohookean elasticity potential formula
    e += neohookean_C * (pow(F.determinant(), -2. / 3.) * (F.transpose() * F).trace() - 3.)
         + neohookean_D * pow(F.determinant() - 1, 2.);
  }
  return e;
}