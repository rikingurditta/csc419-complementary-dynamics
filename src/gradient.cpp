#include "gradient.h"
#include "dF_dU_flattened.h"
#include "dpsi_neo_hookean_dF.h"

// Inputs:
//   V  #V by 3 list of rest (initial pose) mesh vertex positions
//   T  #T by 4 list tetrahedra indices into rows of V
//   Ur #V by 3 rig displacement at the current frame
//   Uc #V by 3 complementary displacement at the current frame
//   dt time step
// Outputs:
//   g  #V by 3 gradient of total energy of the system
double energy(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & T,
  const Eigen::VectorXd & Ur,
  const Eigen::VectorXd & Uc,
  const double dt,
  const double neohookean_C,
  const double neohookean_D,
  Eigen::MatrixXd g) {
  g = Eigen::MatrixXd::Zero(V.rows(), 3);
  for (int t = 0; t < T.rows(); t++) {
    // get deformation gradient
    Eigen::Matrix3d F;
    deformation_gradient(V, T.row(t), Ur, Uc, F);
    // calculate dF/dU, i.e. derivative of deformation gradient with respect to displacements (flattened)
    Eigen::Matrix912d B;
    dF_dU_flattened(V, T.row(t), B);
    // calculate derivative of neohookean potential energy with respect to deformation gradient (flattened)
    Eigen::Vector9d dpsi;
    dpsi_neo_hookean_dF(dpsi, F, neohookean_C, neohookean_D)
    // calculate derivative of neohookean potential energy with respect to displacements (flattened)
    Eigen::VectorXd g_flattened = B.transpose() * dpsi;
    // unflatten
    g.row(T(t, 0)) += g_flattened.segment(0, 3);
    g.row(T(t, 1)) += g_flattened.segment(3, 3);
    g.row(T(t, 2)) += g_flattened.segment(6, 3);
    g.row(T(t, 3)) += g_flattened.segment(9, 3);
  }
}