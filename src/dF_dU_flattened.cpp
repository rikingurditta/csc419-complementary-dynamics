
void dF_dU_flattened(
  const Eigen::MatrixXd & V,
  const Eigen::RowVector4i & tet,
  Eigen::Matrix912d & B) {
  // TODO: maybe change this to transposed form so it doesn't need to be transposed later
  Eigen::Matrix3d Te = Eigen::Matrix3d::Zero();
  Te.col(0) = V.row(v1) - V.row(v0);
  Te.col(1) = V.row(v2) - V.row(v0);
  Te.col(2) = V.row(v3) - V.row(v0);
  Eigen::Matrix3d Te_inv = Te.inverse();
  Eigen::Matrix43d dphi;
  dphi.row(0) = -Eigen::RowVector3d::Ones() * Te_inv;
  dphi.block<3, 3>(1, 0) = Te_inv;

  B = Eigen::Matrix912d::Zero();
  B.block<3, 1>(0, 0) = dphi.block<1, 3>(0, 0);
  B.block<3, 1>(3, 1) = dphi.block<1, 3>(0, 0);
  B.block<3, 1>(6, 2) = dphi.block<1, 3>(0, 0);
  B.block<3, 1>(0, 3) = dphi.block<1, 3>(1, 0);
  B.block<3, 1>(3, 4) = dphi.block<1, 3>(1, 0);
  B.block<3, 1>(6, 5) = dphi.block<1, 3>(1, 0);
  B.block<3, 1>(0, 6) = dphi.block<1, 3>(2, 0);
  B.block<3, 1>(3, 7) = dphi.block<1, 3>(2, 0);
  B.block<3, 1>(6, 8) = dphi.block<1, 3>(2, 0);
  B.block<3, 1>(0, 9) = dphi.block<1, 3>(3, 0);
  B.block<3, 1>(3, 10) = dphi.block<1, 3>(3, 0);
  B.block<3, 1>(6, 11) = dphi.block<1, 3>(3, 0);
}