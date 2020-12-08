#include "gradient.h"
#include "dF_dU_flattened.h"
#include "dphi_neo_hookean_dF.h"

// Inputs:
//   V  #V by 3 list of rest (initial pose) mesh vertex positions
//   T  #T by 4 list tetrahedra indices into rows of V
//   Ur #V by 3 rig displacement at the current frame
//   Uc #V by 3 complementary displacement at the current frame
//   dt time step
// Outputs:
//   g  #V*3 by 1 gradient of total energy of the system
double gradient(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &T,
  const Eigen::VectorXd &Ur,
  const Eigen::VectorXd &Uc,
  const double dt,
  const double neohookean_C,
  const double neohookean_D,
  Eigen::VectorXd g) {
  g = Eigen::VectorXd::Zero(V.rows() * 3);
  for (int t = 0; t < T.rows(); t++) {
    // get deformation gradient
    Eigen::Matrix3d F;
    deformation_gradient(V, T.row(t), Ur, Uc, F);
    // calculate dF/dU, i.e. derivative of deformation gradient with respect to displacements
    Eigen::Matrix912d B;
    dF_dU_flattened(V, T.row(t), B);
    // calculate derivative of neohookean potential energy with respect to deformation gradient
    Eigen::Vector9d dphi;
    dphi_neo_hookean_dF(dphi, F, neohookean_C, neohookean_D)
    // calculate derivative of neohookean potential energy with respect to displacements for current tet (chain rule)
    Eigen::Vector12d g_tet = B.transpose() * dphi;
    // distribute to global gradient vector
    g.segment(T(t, 0) * 3, 3) += g_tet.segment(0, 3);
    g.segment(T(t, 1) * 3, 3) += g_tet.segment(3, 3);
    g.segment(T(t, 2) * 3, 3) += g_tet.segment(6, 3);
    g.segment(T(t, 3) * 3, 3) += g_tet.segment(9, 3);
  }
}