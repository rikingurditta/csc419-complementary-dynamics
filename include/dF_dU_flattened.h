#include <Eigen/Core>
#include <EigenTypes.h>

void dF_dU_flattened(
  const Eigen::MatrixXd & V,
  const Eigen::RowVector4i & tet,
  Eigen::Matrix912d & B);