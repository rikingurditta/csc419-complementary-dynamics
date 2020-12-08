#include <Eigen/Core>

void dF_dU_flattened(
  const Eigen::MatrixXd & V,
  const Eigen::RowVector4i & tet,
  Eigen::MatrixXd & B);