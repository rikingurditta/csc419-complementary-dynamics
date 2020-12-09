#include <iostream>
#include <Eigen/Core>

#include <igl/readMESH.h>
#include <igl/boundary_facets.h>
#include <igl/bounding_box.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/parula.h>
#include <igl/lbs_matrix.h>

#include "rig_jacobian.h"
#include "mass_matrix.h"
#include "energy.h"
#include "gradient.h"
#include "hessian.h"

#include "complementary_displacement.h"


int main(int argc, char *argv[]) {
  // load tet mesh
  Eigen::MatrixXd V;
  Eigen::MatrixXi T;
  Eigen::MatrixXi F;
  igl::readMESH("../data/coarser_bunny.mesh", V, T, F);
  igl::boundary_facets(T, F);
  // boundary_facets returns inside out surface, so reverse direction of every face to turn the correct side outward
  for (int f = 0; f < F.rows(); f++) {
    int t = F(f, 0);
    F(f, 0) = F(f, 1);
    F(f, 1) = t;
  }

  double dt = 0.1;

  //material parameters
  double density = 0.1;
  double YM = 6e5; //young's modulus
  double mu = 0.4; //poissons ratio
  double neohookean_D = 0.5 * (YM * mu) / ((1.0 + mu) * (1.0 - 2.0 * mu));
  double neohookean_C = 0.5 * YM / (2.0 * (1.0 + mu));
  std::cout << "C: " << neohookean_C << " D: " << neohookean_D << "\n";

  Eigen::SparseMatrixd M;
  mass_matrix(V, T, density, M);
  auto curriedMass_matrix = [=](Eigen::SparseMatrixd &M) { mass_matrix(V, T, 1., M); };

  // choose two bones based on corners of bounding box
  // ideally bones would be chosen by animator in a 3d modelling program, but the transfer of data was too difficult,
  // so here we manually place bones
  int NUM_BONES = 2;
  Eigen::MatrixXd BV;
  Eigen::MatrixXd BF;
  igl::bounding_box(V, BV, BF);
  Eigen::RowVector3d lo = (BV.row(2) + BV.row(3)) / 2;
  Eigen::RowVector3d ro = (BV.row(4) + BV.row(5)) / 2;
  Eigen::RowVector3d l = (3 * lo + ro) / 4;
  Eigen::RowVector3d r = (3 * ro + lo) / 4;
  Eigen::MatrixXd bones(NUM_BONES, 3);
  bones.row(0) = l;
  bones.row(1) = r;

  // create bone transformations matrices for each frame of animation
  // again, ideally this would be created by animator in modelling program, but the data transfer was too difficult
  std::vector<Eigen::MatrixXd> T_list;
  for (int i = 0; i < 48; i++) {
    Eigen::Matrix34d T1 = Eigen::Matrix34d::Zero();
    T1.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
    Eigen::Vector3d d = Eigen::Vector3d::Zero();
    d << i, 0, 0;  // displacement increases with frame, creating translation animation
    T1.block<3, 1>(0, 3) = d * 2;
    Eigen::Matrix34d T2 = Eigen::Matrix34d::Zero();
    double angle = (i - 24) * 2. * M_PI / 180;  // angle increases with frame, creating rotation animation
    T2.block<3, 3>(0, 0) << 1, 0, 0,
      0, cos(angle), -sin(angle),
      0, sin(angle), cos(angle);
    // T2.block<3, 1>(0, 3) = d * 2;
    Eigen::Matrix<double, 8, 3> T_curr;
    T_curr.block<4, 3>(0, 0) = T1.transpose();
    T_curr.block<4, 3>(4, 0) = T2.transpose();
    T_list.emplace_back(T_curr);
  }

  // temporary: calculate weights based on inverse squared distance
  // in actual setup, weights are chosen by animator
  // we haven't figured out how to export weights from Blender
  Eigen::MatrixXd weights(V.rows(), 2);
  for (int i = 0; i < V.rows(); i++) {
    weights(i, 0) = 100 / (V.row(i) - l).squaredNorm();
    weights(i, 1) = 100 / (V.row(i) - r).squaredNorm();
    double s = weights.row(i).sum();
    weights.row(i) /= s;
  }

  // calculate skinning matrix
  Eigen::MatrixXd LBS;
  igl::lbs_matrix(V, weights, LBS);

  // calculate rig jacobian
  Eigen::MatrixXd J;
  rig_jacobian(V, weights, J);

  igl::opengl::glfw::Viewer viewer;

  // current vertex positions in animation
  Eigen::MatrixXd X = V;

  Eigen::MatrixXd Uc = Eigen::MatrixXd::Zero(V.rows(), 3);
  Eigen::MatrixXd Uclast = Uc;
  Eigen::MatrixXd duLast = Eigen::MatrixXd::Zero(V.rows(), 3);
  Eigen::VectorXd g = Eigen::VectorXd::Zero(V.rows() * 3);
  Eigen::SparseMatrixd H(V.rows() * 3, V.rows() * 3);

  bool oc = false;
  int frame = 0;
  viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer &) -> bool {
    if (viewer.core().is_animating) {
      // current transformation
      Eigen::MatrixXd T_curr = T_list[frame];
      // get rig deformation using linear blend skinning
      Eigen::MatrixXd Ur = LBS * T_curr - V;
      // get complementary displacement
      complementary_displacement(V, T, M, Ur, Uclast, X, duLast, J, g, H, dt, neohookean_C, neohookean_D, Uc);
      Uclast = Uc;
      std::cout << "|Uc|: " << Uc.squaredNorm() << "\n";
      X = V + Ur + Uc;
      if (oc) {
        // only visualize complementary displacement
        viewer.data().set_vertices(V + Uc);
      } else {
        viewer.data().set_vertices(X);
      }

      frame++;
      if (frame == T_list.size()) {
        frame = 0;
        viewer.core().is_animating = false;
      }
    }
    return false;
  };
  viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod) {
    if (key == ' ') {
      viewer.core().is_animating = !viewer.core().is_animating;
    }
    if (key == 'c') {  // this doesn't work
      oc = !oc;
    }
    return false;
  };

  viewer.data().set_mesh(X, F);
  viewer.core().is_animating = false;
  viewer.core().animation_max_fps = 16.;
  viewer.data().show_lines = false;
  viewer.data().show_overlay_depth = false;
  viewer.data().set_mesh(V, F);

  // colour based on bones - will need to change this when we have more than 2 bones
  Eigen::MatrixXd C;
  igl::parula(weights.col(0), true, C);
  viewer.data().set_colors(C);

  // visualize handles as points
  const Eigen::RowVector3d orange(1.0, 0.7, 0.2);
  viewer.data().set_points(bones, orange);

  viewer.launch();

  return 0;
}
