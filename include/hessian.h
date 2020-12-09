#include <Eigen/Core>
#include <EigenTypes.h>

// Inputs:
//   V  #V by 3 list of rest (initial pose) mesh vertex positions
//   T  #T by 4 list tetrahedra indices into rows of V
//   U  #V by 3 displacement at the current frame
//   dt time step
// Outputs:
//   H  #V*3 by #V*3 hessian of total energy of the system
void hessian(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & T,
  const Eigen::MatrixXd & U,
  const double dt,
  const double neohookean_C,
  const double neohookean_D,
  Eigen::SparseMatrixd H);