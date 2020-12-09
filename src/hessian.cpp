#include "hessian.h"
#include "deformation_gradient.h"
#include "dF_dU_flattened.h"
#include <igl/volume.h>

void d2phi_neo_hookean_dF2(Eigen::Matrix99d &ddw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D);

void hessian(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &T,
  const Eigen::MatrixXd &U,
  const double dt,
  const double neohookean_C,
  const double neohookean_D,
  Eigen::SparseMatrixd H) {
  H.setZero();
  H.resize(V.rows() * 3, V.rows() * 3);
  std::vector<Eigen::Triplet<double>> tl;
  tl.reserve(144 * T.size());
  H.reserve(144 * T.size());
  Eigen::VectorXd tet_volumes(T.rows());
  igl::volume(V, T, tet_volumes);
  for (int t = 0; t < T.rows(); t++) {
    // get deformation gradient
    Eigen::Matrix3d F;
    deformation_gradient(V, T.row(t), U, F);
    // calculate dF/dU, i.e. derivative of deformation gradient with respect to displacements
    Eigen::Matrix912d B;
    dF_dU_flattened(V, T.row(t), B);
    // calculate hessian of neohookean potential energy with respect to deformation gradient
    Eigen::Matrix99d d2phi;
    d2phi_neo_hookean_dF2(d2phi, F, neohookean_C, neohookean_D);
    // calculate hessian of neohookean potential energy with respect to displacements for current tet (chain rule)
    Eigen::Matrix1212d H_tet = tet_volumes(t) * B.transpose() * d2phi * B;
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
  H.setFromTriplets(tl.begin(), tl.end());
}

void d2phi_neo_hookean_dF2(Eigen::Matrix99d &ddw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
  // painstakingly used matlab
  double F1_1 = F(0, 0);
  double F1_2 = F(0, 1);
  double F1_3 = F(0, 2);
  double F2_1 = F(1, 0);
  double F2_2 = F(1, 1);
  double F2_3 = F(1, 2);
  double F3_1 = F(2, 0);
  double F3_2 = F(2, 1);
  double F3_3 = F(2, 2);

  // factor out common terms
  double J = F.determinant();
  double J_neg_2_3 = pow(J, -2. / 3.);
  double J_neg_5_3 = pow(J, -5. / 3.);
  double J_neg_8_3 = pow(J, -8. / 3.);
  double F_sqnorm = F.squaredNorm();

  double F_11_22_m = F1_1*F2_2-F1_2*F2_1;
  double F_11_32_m = F1_1*F3_2-F1_2*F3_1;
  double F_21_32_m = F2_1*F3_2-F2_2*F3_1;
  double F_11_23_m = F1_1*F2_3-F1_3*F2_1;
  double F_12_23_m = F1_2*F2_3-F1_3*F2_2;
  double F_11_33_m = F1_1*F3_3-F1_3*F3_1;
  double F_21_33_m = F2_1*F3_3-F2_3*F3_1;
  double F_12_33_m = F1_2*F3_3-F1_3*F3_2;
  double F_22_33_m = F2_2*F3_3-F2_3*F3_2;

  // F was flattened row-wise, as in lecture
  ddw(0, 0) = C*(J_neg_2_3*2.0+pow(F_22_33_m,2.0)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0)-F1_1*(F_22_33_m)*J_neg_5_3*(8.0/3.0))+D*pow(F_22_33_m,2.0)*2.0;
  ddw(0, 1) = -C*(F1_1*(F_21_33_m)*J_neg_5_3*(-4.0/3.0)+F1_2*(F_22_33_m)*J_neg_5_3*(4.0/3.0)+(F_21_33_m)*(F_22_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_21_33_m)*(F_22_33_m)*2.0;
  ddw(0, 2) = -C*(F1_1*(F_21_32_m)*J_neg_5_3*(4.0/3.0)+F1_3*(F_22_33_m)*J_neg_5_3*(4.0/3.0)-(F_21_32_m)*(F_22_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_21_32_m)*(F_22_33_m)*2.0;
  ddw(0, 3) = -C*(F1_1*(F_12_33_m)*J_neg_5_3*(-4.0/3.0)+F2_1*(F_22_33_m)*J_neg_5_3*(4.0/3.0)+(F_12_33_m)*(F_22_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_12_33_m)*(F_22_33_m)*2.0;
  ddw(0, 4) = -C*(F3_3*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F1_1*(F_11_33_m)*J_neg_5_3*(4.0/3.0)+F2_2*(F_22_33_m)*J_neg_5_3*(4.0/3.0)-(F_11_33_m)*(F_22_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_33_m)*(F_22_33_m)*2.0-D*F3_3*(-J+1.0)*2.0;
  ddw(0, 5) = C*(F3_2*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F1_1*(F_11_32_m)*J_neg_5_3*(4.0/3.0)-F2_3*(F_22_33_m)*J_neg_5_3*(4.0/3.0)-(F_11_32_m)*(F_22_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_32_m)*(F_22_33_m)*2.0+D*F3_2*(-J+1.0)*2.0;
  ddw(0, 6) = -C*(F1_1*(F_12_23_m)*J_neg_5_3*(4.0/3.0)+F3_1*(F_22_33_m)*J_neg_5_3*(4.0/3.0)-(F_12_23_m)*(F_22_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_12_23_m)*(F_22_33_m)*2.0;
  ddw(0, 7) = C*(F2_3*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F1_1*(F_11_23_m)*J_neg_5_3*(4.0/3.0)-F3_2*(F_22_33_m)*J_neg_5_3*(4.0/3.0)-(F_11_23_m)*(F_22_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_23_m)*(F_22_33_m)*2.0+D*F2_3*(-J+1.0)*2.0;
  ddw(0, 8) = -C*(F2_2*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F1_1*(F_11_22_m)*J_neg_5_3*(4.0/3.0)+F3_3*(F_22_33_m)*J_neg_5_3*(4.0/3.0)-(F_11_22_m)*(F_22_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_22_m)*(F_22_33_m)*2.0-D*F2_2*(-J+1.0)*2.0;
  ddw(1, 0) = -C*(F1_1*(F_21_33_m)*J_neg_5_3*(-4.0/3.0)+F1_2*(F_22_33_m)*J_neg_5_3*(4.0/3.0)+(F_21_33_m)*(F_22_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_21_33_m)*(F_22_33_m)*2.0;
  ddw(1, 1) = C*(J_neg_2_3*2.0+pow(F_21_33_m,2.0)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0)+F1_2*(F_21_33_m)*J_neg_5_3*(8.0/3.0))+D*pow(F_21_33_m,2.0)*2.0;
  ddw(1, 2) = -C*(F1_2*(F_21_32_m)*J_neg_5_3*(4.0/3.0)-F1_3*(F_21_33_m)*J_neg_5_3*(4.0/3.0)+(F_21_32_m)*(F_21_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_21_32_m)*(F_21_33_m)*2.0;
  ddw(1, 3) = C*(F3_3*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F1_2*(F_12_33_m)*J_neg_5_3*(4.0/3.0)+F2_1*(F_21_33_m)*J_neg_5_3*(4.0/3.0)+(F_12_33_m)*(F_21_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_12_33_m)*(F_21_33_m)*2.0+D*F3_3*(-J+1.0)*2.0;
  ddw(1, 4) = -C*(F1_2*(F_11_33_m)*J_neg_5_3*(4.0/3.0)-F2_2*(F_21_33_m)*J_neg_5_3*(4.0/3.0)+(F_11_33_m)*(F_21_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_33_m)*(F_21_33_m)*2.0;
  ddw(1, 5) = C*(F3_1*J_neg_5_3*(F_sqnorm)*(-2.0/3.0)+F1_2*(F_11_32_m)*J_neg_5_3*(4.0/3.0)+F2_3*(F_21_33_m)*J_neg_5_3*(4.0/3.0)+(F_11_32_m)*(F_21_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_32_m)*(F_21_33_m)*2.0-D*F3_1*(-J+1.0)*2.0;
  ddw(1, 6) = -C*(F2_3*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F1_2*(F_12_23_m)*J_neg_5_3*(4.0/3.0)-F3_1*(F_21_33_m)*J_neg_5_3*(4.0/3.0)+(F_12_23_m)*(F_21_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_12_23_m)*(F_21_33_m)*2.0-D*F2_3*(-J+1.0)*2.0;
  ddw(1, 7) = C*(F1_2*(F_11_23_m)*J_neg_5_3*(4.0/3.0)+F3_2*(F_21_33_m)*J_neg_5_3*(4.0/3.0)+(F_11_23_m)*(F_21_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_23_m)*(F_21_33_m)*2.0;
  ddw(1, 8) = C*(F2_1*J_neg_5_3*(F_sqnorm)*(2.0/3.0)-F1_2*(F_11_22_m)*J_neg_5_3*(4.0/3.0)+F3_3*(F_21_33_m)*J_neg_5_3*(4.0/3.0)-(F_11_22_m)*(F_21_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_22_m)*(F_21_33_m)*2.0+D*F2_1*(-J+1.0)*2.0;
  ddw(2, 0) = -C*(F1_1*(F_21_32_m)*J_neg_5_3*(4.0/3.0)+F1_3*(F_22_33_m)*J_neg_5_3*(4.0/3.0)-(F_21_32_m)*(F_22_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_21_32_m)*(F_22_33_m)*2.0;
  ddw(2, 1) = -C*(F1_2*(F_21_32_m)*J_neg_5_3*(4.0/3.0)-F1_3*(F_21_33_m)*J_neg_5_3*(4.0/3.0)+(F_21_32_m)*(F_21_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_21_32_m)*(F_21_33_m)*2.0;
  ddw(2, 2) = C*(J_neg_2_3*2.0+pow(F_21_32_m,2.0)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0)-F1_3*(F_21_32_m)*J_neg_5_3*(8.0/3.0))+D*pow(F_21_32_m,2.0)*2.0;
  ddw(2, 3) = -C*(F3_2*J_neg_5_3*(F_sqnorm)*(2.0/3.0)-F1_3*(F_12_33_m)*J_neg_5_3*(4.0/3.0)+F2_1*(F_21_32_m)*J_neg_5_3*(4.0/3.0)+(F_12_33_m)*(F_21_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_12_33_m)*(F_21_32_m)*2.0-D*F3_2*(-J+1.0)*2.0;
  ddw(2, 4) = C*(F3_1*J_neg_5_3*(F_sqnorm)*(2.0/3.0)-F1_3*(F_11_33_m)*J_neg_5_3*(4.0/3.0)-F2_2*(F_21_32_m)*J_neg_5_3*(4.0/3.0)+(F_11_33_m)*(F_21_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_33_m)*(F_21_32_m)*2.0+D*F3_1*(-J+1.0)*2.0;
  ddw(2, 5) = -C*(F1_3*(F_11_32_m)*J_neg_5_3*(-4.0/3.0)+F2_3*(F_21_32_m)*J_neg_5_3*(4.0/3.0)+(F_11_32_m)*(F_21_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_32_m)*(F_21_32_m)*2.0;
  ddw(2, 6) = C*(F2_2*J_neg_5_3*(F_sqnorm)*(2.0/3.0)-F1_3*(F_12_23_m)*J_neg_5_3*(4.0/3.0)-F3_1*(F_21_32_m)*J_neg_5_3*(4.0/3.0)+(F_12_23_m)*(F_21_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_12_23_m)*(F_21_32_m)*2.0+D*F2_2*(-J+1.0)*2.0;
  ddw(2, 7) = -C*(F2_1*J_neg_5_3*(F_sqnorm)*(2.0/3.0)-F1_3*(F_11_23_m)*J_neg_5_3*(4.0/3.0)+F3_2*(F_21_32_m)*J_neg_5_3*(4.0/3.0)+(F_11_23_m)*(F_21_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_23_m)*(F_21_32_m)*2.0-D*F2_1*(-J+1.0)*2.0;
  ddw(2, 8) = -C*(F1_3*(F_11_22_m)*J_neg_5_3*(4.0/3.0)+F3_3*(F_21_32_m)*J_neg_5_3*(4.0/3.0)-(F_11_22_m)*(F_21_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_22_m)*(F_21_32_m)*2.0;
  ddw(3, 0) = -C*(F1_1*(F_12_33_m)*J_neg_5_3*(-4.0/3.0)+F2_1*(F_22_33_m)*J_neg_5_3*(4.0/3.0)+(F_12_33_m)*(F_22_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_12_33_m)*(F_22_33_m)*2.0;
  ddw(3, 1) = C*(F3_3*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F1_2*(F_12_33_m)*J_neg_5_3*(4.0/3.0)+F2_1*(F_21_33_m)*J_neg_5_3*(4.0/3.0)+(F_12_33_m)*(F_21_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_12_33_m)*(F_21_33_m)*2.0+D*F3_3*(-J+1.0)*2.0;
  ddw(3, 2) = -C*(F3_2*J_neg_5_3*(F_sqnorm)*(2.0/3.0)-F1_3*(F_12_33_m)*J_neg_5_3*(4.0/3.0)+F2_1*(F_21_32_m)*J_neg_5_3*(4.0/3.0)+(F_12_33_m)*(F_21_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_12_33_m)*(F_21_32_m)*2.0-D*F3_2*(-J+1.0)*2.0;
  ddw(3, 3) = C*(J_neg_2_3*2.0+pow(F_12_33_m,2.0)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0)+F2_1*(F_12_33_m)*J_neg_5_3*(8.0/3.0))+D*pow(F_12_33_m,2.0)*2.0;
  ddw(3, 4) = -C*(F2_1*(F_11_33_m)*J_neg_5_3*(4.0/3.0)-F2_2*(F_12_33_m)*J_neg_5_3*(4.0/3.0)+(F_11_33_m)*(F_12_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_33_m)*(F_12_33_m)*2.0;
  ddw(3, 5) = C*(F2_1*(F_11_32_m)*J_neg_5_3*(4.0/3.0)+F2_3*(F_12_33_m)*J_neg_5_3*(4.0/3.0)+(F_11_32_m)*(F_12_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_32_m)*(F_12_33_m)*2.0;
  ddw(3, 6) = -C*(F2_1*(F_12_23_m)*J_neg_5_3*(4.0/3.0)-F3_1*(F_12_33_m)*J_neg_5_3*(4.0/3.0)+(F_12_23_m)*(F_12_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_12_23_m)*(F_12_33_m)*2.0;
  ddw(3, 7) = C*(F1_3*J_neg_5_3*(F_sqnorm)*(-2.0/3.0)+F2_1*(F_11_23_m)*J_neg_5_3*(4.0/3.0)+F3_2*(F_12_33_m)*J_neg_5_3*(4.0/3.0)+(F_11_23_m)*(F_12_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_23_m)*(F_12_33_m)*2.0-D*F1_3*(-J+1.0)*2.0;
  ddw(3, 8) = C*(F1_2*J_neg_5_3*(F_sqnorm)*(2.0/3.0)-F2_1*(F_11_22_m)*J_neg_5_3*(4.0/3.0)+F3_3*(F_12_33_m)*J_neg_5_3*(4.0/3.0)-(F_11_22_m)*(F_12_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_22_m)*(F_12_33_m)*2.0+D*F1_2*(-J+1.0)*2.0;
  ddw(4, 0) = -C*(F3_3*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F1_1*(F_11_33_m)*J_neg_5_3*(4.0/3.0)+F2_2*(F_22_33_m)*J_neg_5_3*(4.0/3.0)-(F_11_33_m)*(F_22_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_33_m)*(F_22_33_m)*2.0-D*F3_3*(-J+1.0)*2.0;
  ddw(4, 1) = -C*(F1_2*(F_11_33_m)*J_neg_5_3*(4.0/3.0)-F2_2*(F_21_33_m)*J_neg_5_3*(4.0/3.0)+(F_11_33_m)*(F_21_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_33_m)*(F_21_33_m)*2.0;
  ddw(4, 2) = C*(F3_1*J_neg_5_3*(F_sqnorm)*(2.0/3.0)-F1_3*(F_11_33_m)*J_neg_5_3*(4.0/3.0)-F2_2*(F_21_32_m)*J_neg_5_3*(4.0/3.0)+(F_11_33_m)*(F_21_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_33_m)*(F_21_32_m)*2.0+D*F3_1*(-J+1.0)*2.0;
  ddw(4, 3) = -C*(F2_1*(F_11_33_m)*J_neg_5_3*(4.0/3.0)-F2_2*(F_12_33_m)*J_neg_5_3*(4.0/3.0)+(F_11_33_m)*(F_12_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_33_m)*(F_12_33_m)*2.0;
  ddw(4, 4) = C*(J_neg_2_3*2.0+pow(F_11_33_m,2.0)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0)-F2_2*(F_11_33_m)*J_neg_5_3*(8.0/3.0))+D*pow(F_11_33_m,2.0)*2.0;
  ddw(4, 5) = -C*(F2_2*(F_11_32_m)*J_neg_5_3*(-4.0/3.0)+F2_3*(F_11_33_m)*J_neg_5_3*(4.0/3.0)+(F_11_32_m)*(F_11_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_32_m)*(F_11_33_m)*2.0;
  ddw(4, 6) = C*(F1_3*J_neg_5_3*(F_sqnorm)*(2.0/3.0)-F2_2*(F_12_23_m)*J_neg_5_3*(4.0/3.0)-F3_1*(F_11_33_m)*J_neg_5_3*(4.0/3.0)+(F_12_23_m)*(F_11_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_12_23_m)*(F_11_33_m)*2.0+D*F1_3*(-J+1.0)*2.0;
  ddw(4, 7) = -C*(F2_2*(F_11_23_m)*J_neg_5_3*(-4.0/3.0)+F3_2*(F_11_33_m)*J_neg_5_3*(4.0/3.0)+(F_11_23_m)*(F_11_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_23_m)*(F_11_33_m)*2.0;
  ddw(4, 8) = -C*(F1_1*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F2_2*(F_11_22_m)*J_neg_5_3*(4.0/3.0)+F3_3*(F_11_33_m)*J_neg_5_3*(4.0/3.0)-(F_11_22_m)*(F_11_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_22_m)*(F_11_33_m)*2.0-D*F1_1*(-J+1.0)*2.0;
  ddw(5, 0) = C*(F3_2*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F1_1*(F_11_32_m)*J_neg_5_3*(4.0/3.0)-F2_3*(F_22_33_m)*J_neg_5_3*(4.0/3.0)-(F_11_32_m)*(F_22_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_32_m)*(F_22_33_m)*2.0+D*F3_2*(-J+1.0)*2.0;
  ddw(5, 1) = C*(F3_1*J_neg_5_3*(F_sqnorm)*(-2.0/3.0)+F1_2*(F_11_32_m)*J_neg_5_3*(4.0/3.0)+F2_3*(F_21_33_m)*J_neg_5_3*(4.0/3.0)+(F_11_32_m)*(F_21_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_32_m)*(F_21_33_m)*2.0-D*F3_1*(-J+1.0)*2.0;
  ddw(5, 2) = -C*(F1_3*(F_11_32_m)*J_neg_5_3*(-4.0/3.0)+F2_3*(F_21_32_m)*J_neg_5_3*(4.0/3.0)+(F_11_32_m)*(F_21_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_32_m)*(F_21_32_m)*2.0;
  ddw(5, 3) = C*(F2_1*(F_11_32_m)*J_neg_5_3*(4.0/3.0)+F2_3*(F_12_33_m)*J_neg_5_3*(4.0/3.0)+(F_11_32_m)*(F_12_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_32_m)*(F_12_33_m)*2.0;
  ddw(5, 4) = -C*(F2_2*(F_11_32_m)*J_neg_5_3*(-4.0/3.0)+F2_3*(F_11_33_m)*J_neg_5_3*(4.0/3.0)+(F_11_32_m)*(F_11_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_32_m)*(F_11_33_m)*2.0;
  ddw(5, 5) = C*(J_neg_2_3*2.0+pow(F_11_32_m,2.0)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0)+F2_3*(F_11_32_m)*J_neg_5_3*(8.0/3.0))+D*pow(F_11_32_m,2.0)*2.0;
  ddw(5, 6) = -C*(F1_2*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F2_3*(F_12_23_m)*J_neg_5_3*(4.0/3.0)-F3_1*(F_11_32_m)*J_neg_5_3*(4.0/3.0)+(F_12_23_m)*(F_11_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_12_23_m)*(F_11_32_m)*2.0-D*F1_2*(-J+1.0)*2.0;
  ddw(5, 7) = C*(F1_1*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F2_3*(F_11_23_m)*J_neg_5_3*(4.0/3.0)+F3_2*(F_11_32_m)*J_neg_5_3*(4.0/3.0)+(F_11_23_m)*(F_11_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_23_m)*(F_11_32_m)*2.0+D*F1_1*(-J+1.0)*2.0;
  ddw(5, 8) = -C*(F2_3*(F_11_22_m)*J_neg_5_3*(4.0/3.0)-F3_3*(F_11_32_m)*J_neg_5_3*(4.0/3.0)+(F_11_22_m)*(F_11_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_22_m)*(F_11_32_m)*2.0;
  ddw(6, 0) = -C*(F1_1*(F_12_23_m)*J_neg_5_3*(4.0/3.0)+F3_1*(F_22_33_m)*J_neg_5_3*(4.0/3.0)-(F_12_23_m)*(F_22_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_12_23_m)*(F_22_33_m)*2.0;
  ddw(6, 1) = -C*(F2_3*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F1_2*(F_12_23_m)*J_neg_5_3*(4.0/3.0)-F3_1*(F_21_33_m)*J_neg_5_3*(4.0/3.0)+(F_12_23_m)*(F_21_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_12_23_m)*(F_21_33_m)*2.0-D*F2_3*(-J+1.0)*2.0;
  ddw(6, 2) = C*(F2_2*J_neg_5_3*(F_sqnorm)*(2.0/3.0)-F1_3*(F_12_23_m)*J_neg_5_3*(4.0/3.0)-F3_1*(F_21_32_m)*J_neg_5_3*(4.0/3.0)+(F_12_23_m)*(F_21_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_12_23_m)*(F_21_32_m)*2.0+D*F2_2*(-J+1.0)*2.0;
  ddw(6, 3) = -C*(F2_1*(F_12_23_m)*J_neg_5_3*(4.0/3.0)-F3_1*(F_12_33_m)*J_neg_5_3*(4.0/3.0)+(F_12_23_m)*(F_12_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_12_23_m)*(F_12_33_m)*2.0;
  ddw(6, 4) = C*(F1_3*J_neg_5_3*(F_sqnorm)*(2.0/3.0)-F2_2*(F_12_23_m)*J_neg_5_3*(4.0/3.0)-F3_1*(F_11_33_m)*J_neg_5_3*(4.0/3.0)+(F_12_23_m)*(F_11_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_12_23_m)*(F_11_33_m)*2.0+D*F1_3*(-J+1.0)*2.0;
  ddw(6, 5) = -C*(F1_2*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F2_3*(F_12_23_m)*J_neg_5_3*(4.0/3.0)-F3_1*(F_11_32_m)*J_neg_5_3*(4.0/3.0)+(F_12_23_m)*(F_11_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_12_23_m)*(F_11_32_m)*2.0-D*F1_2*(-J+1.0)*2.0;
  ddw(6, 6) = C*(J_neg_2_3*2.0+pow(F_12_23_m,2.0)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0)-F3_1*(F_12_23_m)*J_neg_5_3*(8.0/3.0))+D*pow(F_12_23_m,2.0)*2.0;
  ddw(6, 7) = -C*(F3_1*(F_11_23_m)*J_neg_5_3*(-4.0/3.0)+F3_2*(F_12_23_m)*J_neg_5_3*(4.0/3.0)+(F_11_23_m)*(F_12_23_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_23_m)*(F_12_23_m)*2.0;
  ddw(6, 8) = -C*(F3_1*(F_11_22_m)*J_neg_5_3*(4.0/3.0)+F3_3*(F_12_23_m)*J_neg_5_3*(4.0/3.0)-(F_11_22_m)*(F_12_23_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_22_m)*(F_12_23_m)*2.0;
  ddw(7, 0) = C*(F2_3*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F1_1*(F_11_23_m)*J_neg_5_3*(4.0/3.0)-F3_2*(F_22_33_m)*J_neg_5_3*(4.0/3.0)-(F_11_23_m)*(F_22_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_23_m)*(F_22_33_m)*2.0+D*F2_3*(-J+1.0)*2.0;
  ddw(7, 1) = C*(F1_2*(F_11_23_m)*J_neg_5_3*(4.0/3.0)+F3_2*(F_21_33_m)*J_neg_5_3*(4.0/3.0)+(F_11_23_m)*(F_21_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_23_m)*(F_21_33_m)*2.0;
  ddw(7, 2) = -C*(F2_1*J_neg_5_3*(F_sqnorm)*(2.0/3.0)-F1_3*(F_11_23_m)*J_neg_5_3*(4.0/3.0)+F3_2*(F_21_32_m)*J_neg_5_3*(4.0/3.0)+(F_11_23_m)*(F_21_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_23_m)*(F_21_32_m)*2.0-D*F2_1*(-J+1.0)*2.0;
  ddw(7, 3) = C*(F1_3*J_neg_5_3*(F_sqnorm)*(-2.0/3.0)+F2_1*(F_11_23_m)*J_neg_5_3*(4.0/3.0)+F3_2*(F_12_33_m)*J_neg_5_3*(4.0/3.0)+(F_11_23_m)*(F_12_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_23_m)*(F_12_33_m)*2.0-D*F1_3*(-J+1.0)*2.0;
  ddw(7, 4) = -C*(F2_2*(F_11_23_m)*J_neg_5_3*(-4.0/3.0)+F3_2*(F_11_33_m)*J_neg_5_3*(4.0/3.0)+(F_11_23_m)*(F_11_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_23_m)*(F_11_33_m)*2.0;
  ddw(7, 5) = C*(F1_1*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F2_3*(F_11_23_m)*J_neg_5_3*(4.0/3.0)+F3_2*(F_11_32_m)*J_neg_5_3*(4.0/3.0)+(F_11_23_m)*(F_11_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_23_m)*(F_11_32_m)*2.0+D*F1_1*(-J+1.0)*2.0;
  ddw(7, 6) = -C*(F3_1*(F_11_23_m)*J_neg_5_3*(-4.0/3.0)+F3_2*(F_12_23_m)*J_neg_5_3*(4.0/3.0)+(F_11_23_m)*(F_12_23_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_23_m)*(F_12_23_m)*2.0;
  ddw(7, 7) = C*(J_neg_2_3*2.0+pow(F_11_23_m,2.0)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0)+F3_2*(F_11_23_m)*J_neg_5_3*(8.0/3.0))+D*pow(F_11_23_m,2.0)*2.0;
  ddw(7, 8) = -C*(F3_2*(F_11_22_m)*J_neg_5_3*(4.0/3.0)-F3_3*(F_11_23_m)*J_neg_5_3*(4.0/3.0)+(F_11_22_m)*(F_11_23_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_22_m)*(F_11_23_m)*2.0;
  ddw(8, 0) = -C*(F2_2*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F1_1*(F_11_22_m)*J_neg_5_3*(4.0/3.0)+F3_3*(F_22_33_m)*J_neg_5_3*(4.0/3.0)-(F_11_22_m)*(F_22_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_22_m)*(F_22_33_m)*2.0-D*F2_2*(-J+1.0)*2.0;
  ddw(8, 1) = C*(F2_1*J_neg_5_3*(F_sqnorm)*(2.0/3.0)-F1_2*(F_11_22_m)*J_neg_5_3*(4.0/3.0)+F3_3*(F_21_33_m)*J_neg_5_3*(4.0/3.0)-(F_11_22_m)*(F_21_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_22_m)*(F_21_33_m)*2.0+D*F2_1*(-J+1.0)*2.0;
  ddw(8, 2) = -C*(F1_3*(F_11_22_m)*J_neg_5_3*(4.0/3.0)+F3_3*(F_21_32_m)*J_neg_5_3*(4.0/3.0)-(F_11_22_m)*(F_21_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_22_m)*(F_21_32_m)*2.0;
  ddw(8, 3) = C*(F1_2*J_neg_5_3*(F_sqnorm)*(2.0/3.0)-F2_1*(F_11_22_m)*J_neg_5_3*(4.0/3.0)+F3_3*(F_12_33_m)*J_neg_5_3*(4.0/3.0)-(F_11_22_m)*(F_12_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_22_m)*(F_12_33_m)*2.0+D*F1_2*(-J+1.0)*2.0;
  ddw(8, 4) = -C*(F1_1*J_neg_5_3*(F_sqnorm)*(2.0/3.0)+F2_2*(F_11_22_m)*J_neg_5_3*(4.0/3.0)+F3_3*(F_11_33_m)*J_neg_5_3*(4.0/3.0)-(F_11_22_m)*(F_11_33_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_22_m)*(F_11_33_m)*2.0-D*F1_1*(-J+1.0)*2.0;
  ddw(8, 5) = -C*(F2_3*(F_11_22_m)*J_neg_5_3*(4.0/3.0)-F3_3*(F_11_32_m)*J_neg_5_3*(4.0/3.0)+(F_11_22_m)*(F_11_32_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_22_m)*(F_11_32_m)*2.0;
  ddw(8, 6) = -C*(F3_1*(F_11_22_m)*J_neg_5_3*(4.0/3.0)+F3_3*(F_12_23_m)*J_neg_5_3*(4.0/3.0)-(F_11_22_m)*(F_12_23_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))+D*(F_11_22_m)*(F_12_23_m)*2.0;
  ddw(8, 7) = -C*(F3_2*(F_11_22_m)*J_neg_5_3*(4.0/3.0)-F3_3*(F_11_23_m)*J_neg_5_3*(4.0/3.0)+(F_11_22_m)*(F_11_23_m)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0))-D*(F_11_22_m)*(F_11_23_m)*2.0;
  ddw(8, 8) = C*(J_neg_2_3*2.0+pow(F_11_22_m,2.0)*J_neg_8_3*(F_sqnorm)*(1.0E+1/9.0)-F3_3*(F_11_22_m)*J_neg_5_3*(8.0/3.0))+D*pow(F_11_22_m,2.0)*2.0;
}