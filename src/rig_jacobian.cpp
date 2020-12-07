#include "rig_jacobian.h"


void rig_jacobian(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXd & W,
        Eigen::MatrixXd & J) {
  J.resize(V.rows() * 3, W.cols() * 12);
  for (int i = 0; i < V.rows(); i++) {
    for (int j = 0; j < W.cols(); j++) {
//      Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
      Eigen::MatrixXd horizontal(1, 4);
      horizontal.block<1, 3>(0, 0) = V.row(i);
      horizontal(3) = 1;
      Eigen::MatrixXd kronecker = Eigen::MatrixXd::Zero(3, 12);
      kronecker.block<1, 4>(0, 0) = horizontal;
      kronecker.block<1, 4>(1, 4) = horizontal;
      kronecker.block<1, 4>(2, 8) = horizontal;
      J.block<3, 12>(i * 3, j * 12) = kronecker;
    }
  }
}