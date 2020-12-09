#include <Eigen/Core>

// Inputs:
//   V  #V by 3 list of rest (initial pose) mesh vertex positions
//   T  #T by 4 list tetrahedra indices into rows of V
//   U  #V by 3 displacement at the current frame
//   dt time step
// Outputs:
//   return total energy of the system
double energy(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & T,
  const Eigen::MatrixXd & U,
  double dt,
  double neohookean_C,
  double neohookean_D);