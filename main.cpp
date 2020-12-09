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

double t = 0;
double dt = 0.01;

void simulate() {

}

int main(int argc, char *argv[]) {
  Eigen::MatrixXd V;
  Eigen::MatrixXi T;
  Eigen::MatrixXi F;
  igl::readMESH("../data/coarse_bunny.mesh", V, T, F);
  igl::boundary_facets(T, F);
  // boundary_facets returns inside out surface, so reverse direction of every face
  for (int f = 0; f < F.rows(); f++) {
    int t = F(f, 0);
    F(f, 0) = F(f, 1);
    F(f, 1) = t;
  }

  //material parameters
  double density = 0.1;
  double YM = 6e5; //young's modulus
  double mu = 0.4; //poissons ratio
  double neohookean_D = 0.5 * (YM * mu) / ((1.0 + mu) * (1.0 - 2.0 * mu));
  double neohookean_C = 0.5 * YM / (2.0 * (1.0 + mu));

  Eigen::MatrixXd U = Eigen::MatrixXd::Ones(V.rows(), 3);
  std::cout << energy(V, T, U, 0.1, neohookean_C, neohookean_D) << "\n";
  auto curriedEnergy = [=](Eigen::MatrixXd Ut){ energy(V, T, Ut,0.1, neohookean_C, neohookean_D);};
  Eigen::VectorXd grad;
  gradient(V, T, U, 0.1, neohookean_C, neohookean_D, grad);
  auto curriedGradient = [=](Eigen::MatrixXd Ut, Eigen::VectorXd & grad){ gradient(V, T, Ut, 0.1, neohookean_C, neohookean_D, grad);};
  Eigen::SparseMatrixd hess;
  hessian(V, T, U, 0.1, neohookean_C, neohookean_D, hess);
  auto curriedHessian = [=](Eigen::MatrixXd Ut, Eigen::SparseMatrixd &hess) { hessian(V, T, Ut, 0.1, neohookean_C, neohookean_D, hess);};
  Eigen::SparseMatrixd M;
  mass_matrix(V, T, 1., M);
  auto curriedMass_matrix = [=](Eigen::SparseMatrixd &M) { mass_matrix(V, T, 1., M);};

  // choose two bones based on corners of bounding box
  // not final, just for now
  // in an actual setup, bones are chosen by animator
  Eigen::MatrixXd BV;
  Eigen::MatrixXd BF;
  igl::bounding_box(V, BV, BF);
  Eigen::RowVector3d l = (BV.row(2) + BV.row(3)) / 2;
  Eigen::RowVector3d r = (BV.row(4) + BV.row(5)) / 2;

  int NUM_BONES = 2;
  // create bones matrix
  Eigen::MatrixXd bones(NUM_BONES, 3);
  bones.row(0) = l;
  bones.row(1) = r;

  // create transformations matrix
  Eigen::MatrixXd bonesT(NUM_BONES, 12);
  // identity affine transformation is identity matrix and 0 vector translation
  Eigen::Matrix<double, 3, 4> affine_identity = Eigen::MatrixXd::Zero(3, 4);
  affine_identity.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
  // flatten matrix into row for storage
  // not sure what flattening is required in paper, not sure if it matters
  bonesT.row(0) = Eigen::Map<Eigen::Matrix<double, 12, 1>>(affine_identity.data());
  bonesT.row(1) = Eigen::Map<Eigen::Matrix<double, 12, 1>>(affine_identity.data());

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

  // calculate rig jacobian
  Eigen::MatrixXd J;
  rig_jacobian(V, weights, J);
  // std::cout << J << "\n";

  igl::opengl::glfw::Viewer viewer;

  double angle = 2. * M_PI / 180;
  Eigen::Matrix3d rot_angle;
  rot_angle << 1, 0, 0,
    0, cos(angle), -sin(angle),
    0, sin(angle), cos(angle);

  std::vector<Eigen::MatrixXd> T_list;
  for (int i = 0; i < 24; i++) {
    Eigen::Matrix34d T1 = Eigen::Matrix34d::Zero();
    T1.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
    Eigen::Vector3d d = Eigen::Vector3d::Zero();
    d << i, 0, 0;
    T1.block<3, 1>(0, 3) = d * 2;
    Eigen::Matrix34d T2 = Eigen::Matrix34d::Zero();
    T2.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
    for (int j = 0; j < i; j++) {
      T2.block<3, 3>(0, 0) *= rot_angle;
    }
//    T2.block<3, 1>(0, 3) = d * 2;
    Eigen::Matrix<double, 8, 3> T_curr;
    T_curr.block<4, 3>(0, 0) = T1.transpose();
    T_curr.block<4, 3>(4, 0) = T2.transpose();
    T_list.emplace_back(T_curr);
  }

  Eigen::MatrixXd LBS;
  igl::lbs_matrix(V, weights, LBS);

  Eigen::MatrixXd TU = V;

  int frame = 0;
  viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer &) -> bool {
    if (viewer.core().is_animating) {
      Eigen::MatrixXd T = T_list[frame];
      TU = LBS * T;
      frame++;

      viewer.data().set_vertices(TU);
      viewer.data().compute_normals();

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
    return false;
  };


  viewer.data().set_mesh(TU, F);

  viewer.core().is_animating = false;
  viewer.core().animation_max_fps = 24.;
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
