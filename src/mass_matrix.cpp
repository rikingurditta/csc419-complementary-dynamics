#include <mass_matrix.h>
#include <igl/volume.h>

void mass_matrix(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &T,
  double density,
  Eigen::SparseMatrixd &M) {
  M.setZero();
  M.resize(V.rows() * 3, V.rows() * 3);
  std::vector <Eigen::Triplet<double>> tl;
  tl.reserve(48 * T.rows());
  Eigen::VectorXd tet_volumes(T.rows());
  igl::volume(V, T, tet_volumes);
  for (int t = 0; t < T.rows(); t++) {
    // el_row, el_col loops: iterate over pairs of vertices of the tetrahedron
    for (int el_row = 0; el_row < 4; el_row++) {
      for (int el_col = 0; el_col < 4; el_col++) {
        // i loop: iterate over diagonal of 3x3 block for el_row, el_col in tetrahedron mass matrix
        // only need to iterate over diagonal because each 3x3 block is diagonal
        for (int i = 0; i < 3; i++) {
          if (el_row == el_col)
            tl.emplace_back(T(t, el_row) * 3 + i, T(t, el_col) * 3 + i, density * tet_volumes(t) / 10.);
          else
            tl.emplace_back(T(t, el_row) * 3 + i, T(t, el_col) * 3 + i, density * tet_volumes(t) / 20.);
        }
      }
    }
  }
  M.setFromTriplets(tl.begin(), tl.end());
}