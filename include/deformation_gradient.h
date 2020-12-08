#include <Eigen/Core>

void deformation_gradient(
  const Eigen::MatrixXd & V,
  const Eigen::RowVector4d & tet,
  const Eigen::VectorXd & Ur,
  const Eigen::VectorXd & Uc,
  Eigen::Matrix3d & F);