#include <Eigen/Core>
#include <EigenTypes.h>

// Inputs:
//   V  #V by 3 list of rest (initial pose) mesh vertex positions
//   T  #T by 4 list tetrahedra indices into rows of V
//   Ur #V by 3 rig displacement at the current frame
//   Uc #V by 3 complementary displacement at the current frame
//   dt time step
// Outputs:
//   H  #V*3 by #V*3 hessian of total energy of the system
double hessian(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & T,
  const Eigen::VectorXd & Ur,
  const Eigen::VectorXd & Uc,
  const double dt,
  const double neohookean_C,
  const double neohookean_D,
  Eigen::SparseMatrixd H);