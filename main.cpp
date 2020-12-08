#include <iostream>
#include <Eigen/Core>

#include <igl/readMESH.h>
#include <igl/boundary_facets.h>
#include <igl/bounding_box.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/parula.h>

#include "rig_jacobian.h"

int main(int argc, char *argv[]) {
  Eigen::MatrixXd V;
  Eigen::MatrixXd T;
  Eigen::MatrixXi F;
  igl::readMESH("../data/coarse_bunny.mesh", V, T, F);
  igl::boundary_facets(T, F);
  // boundary_facets returns inside out surface, so reverse direction of every face
  for (int f = 0; f < F.rows(); f++) {
    int t = F(f, 0);
    F(f, 0) = F(f, 1);
    F(f, 1) = t;
  }

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
