#include "hessian.h"
#include "dF_dU_flattened.h"
#include "d2phi_neo_hookean_dF2.h"

double hessian(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &T,
  const Eigen::VectorXd &Ur,
  const Eigen::VectorXd &Uc,
  const double dt,
  const double neohookean_C,
  const double neohookean_D,
  Eigen::SparseMatrixd H) {
  H.setZero();
  H.resize(V.rows() * 3, V.rows() * 3);
  std::vector<Eigen::Triplet<double>> tl;
  tl.reserve(144 * T.size());
  H.reserve(144 * T.size());
  for (int t = 0; t < T.rows(); t++) {
    // get deformation gradient
    Eigen::Matrix3d F;
    deformation_gradient(V, T.row(t), Ur, Uc, F);
    // calculate dF/dU, i.e. derivative of deformation gradient with respect to displacements
    Eigen::Matrix912d B;
    dF_dU_flattened(V, T.row(t), B);
    // calculate hessian of neohookean potential energy with respect to deformation gradient
    Eigen::Matrix99d d2phi;
    d2phi_neo_hookean_dF2(d2phi, F, neohookean_C, neohookean_D);
    // calculate hessian of neohookean potential energy with respect to displacements for current tet (chain rule)
    Eigen::Matrix1212d H_tet = B.transpose() * d2phi * B;
    // distribute to global hessian matrix
    for (int el_row = 0; el_row < 4; el_row++) {
      for (int el_col = 0; el_col < 4; el_col++) {
        // i, j loops: iterate over diagonal of 3x3 block for el_row, el_col in tetrahedron hessian
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            tl.emplace_back(T(t, el_row) * 3 + i, T(t, el_col) * 3 + j,
                            H_tet(el_row * 3 + i, el_col * 3 + j));
          }
        }
      }
    }
  }
  H.setFromTriplets(tl.begin(), tl.end())
}