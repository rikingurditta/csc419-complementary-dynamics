#include "energy.h"
#include "deformation_gradient.h"

double energy(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &T,
  const Eigen::VectorXd &Ur,
  const Eigen::VectorXd &Uc,
  const double dt,
  const double neohookean_C,
  const double neohookean_D) {
  // sum potential energy for each individual tetrahedron to get potential of entire mesh
  double e = 0;
  for (int t = 0; t < T.rows(); t++) {
    Eigen::Matrix3d F;
    deformation_gradient(V, T.row(t), Ur, Uc, F);
    // calculate energy for current tetrahedron using neohookean elasticity potential formula
    e += neohookean_C * (pow(F.determinant(), -2. / 3.) * (F.transpose() * F).trace() - 3.)
         + neohookean_D * pow(F.determinant() - 1, 2.);
  }
  return e;
}