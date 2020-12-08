#include <Eigen/Core>
#include <EigenTypes.h>

// Inputs:
//   V    #V by 3 list of rest (initial pose) mesh vertex positions
//   tet  list of vertices of tetrahedron
//   Ur   #V by 3 rig displacement at the current frame
//   Uc   #V by 3 complementary displacement at the current frame
// Outputs:
//   F    3 by 3 deformation gradient of tetrahedron
void deformation_gradient(
  const Eigen::MatrixXd & V,
  const Eigen::RowVector4i & tet,
  const Eigen::VectorXd & Ur,
  const Eigen::VectorXd & Uc,
  Eigen::Matrix3d & F);