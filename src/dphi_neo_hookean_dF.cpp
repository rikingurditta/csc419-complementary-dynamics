#include <dphi_neo_hookean_dF.h>

void dpsi_neo_hookean_dF(Eigen::Vector9d &dw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
    // used matlab lol
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
    double F_sqnorm = F.squaredNorm();

    double F_11_22_m = F(0,0)*F(1,1)-F(0,1)*F(1,0);
    double F_11_32_m = F(0,0)*F(2,1)-F(0,1)*F(2,0);
    double F_21_32_m = F(1,0)*F(2,1)-F(1,1)*F(2,0);
    double F_11_23_m = F(0,0)*F(1,2)-F(0,2)*F(1,0);
    double F_12_23_m = F(0,1)*F(1,2)-F(0,2)*F(1,1);
    double F_11_33_m = F(0,0)*F(2,2)-F(0,2)*F(2,0);
    double F_21_33_m = F(1,0)*F(2,2)-F(1,2)*F(2,0);
    double F_12_33_m = F(0,1)*F(2,2)-F(0,2)*F(2,1);
    double F_22_33_m = F(1,1)*F(2,2)-F(1,2)*F(2,1);

    // F was flattened row-wise, as in lecture
    dw(0) = C*(F1_1*J_neg_2_3*2.0-(F_22_33_m)*J_neg_5_3*(F_sqnorm)*(2.0/3.0))-D*(F_22_33_m)*(-J+1.0)*2.0;
    dw(1) = C*(F1_2*J_neg_2_3*2.0+(F_21_33_m)*J_neg_5_3*(F_sqnorm)*(2.0/3.0))+D*(F_21_33_m)*(-J+1.0)*2.0;
    dw(2) = C*(F1_3*J_neg_2_3*2.0-(F_21_32_m)*J_neg_5_3*(F_sqnorm)*(2.0/3.0))-D*(F_21_32_m)*(-J+1.0)*2.0;
    dw(3) = C*(F2_1*J_neg_2_3*2.0+(F_12_33_m)*J_neg_5_3*(F_sqnorm)*(2.0/3.0))+D*(F_12_33_m)*(-J+1.0)*2.0;
    dw(4) = C*(F2_2*J_neg_2_3*2.0-(F_11_33_m)*J_neg_5_3*(F_sqnorm)*(2.0/3.0))-D*(F_11_33_m)*(-J+1.0)*2.0;
    dw(5) = C*(F2_3*J_neg_2_3*2.0+(F_11_32_m)*J_neg_5_3*(F_sqnorm)*(2.0/3.0))+D*(F_11_32_m)*(-J+1.0)*2.0;
    dw(6) = C*(F3_1*J_neg_2_3*2.0-(F_12_23_m)*J_neg_5_3*(F_sqnorm)*(2.0/3.0))-D*(F_12_23_m)*(-J+1.0)*2.0;
    dw(7) = C*(F3_2*J_neg_2_3*2.0+(F_11_23_m)*J_neg_5_3*(F_sqnorm)*(2.0/3.0))+D*(F_11_23_m)*(-J+1.0)*2.0;
    dw(8) = C*(F3_3*J_neg_2_3*2.0-(F_11_22_m)*J_neg_5_3*(F_sqnorm)*(2.0/3.0))-D*(F_11_22_m)*(-J+1.0)*2.0;
}