#include <Eigen/Dense>
#include <EigenTypes.h>

// Inputs:
//   V        #V by 3 list of rest (initial pose) mesh vertex positions
//   T        #T by 4 list tetrahedra indices into rows of V
//   density  density of material
// Outputs:
//   M        #V*3 by #V*3 sparse FEM mass matrix
void mass_matrix(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & T,
  double density,
  Eigen::SparseMatrixd & M);
