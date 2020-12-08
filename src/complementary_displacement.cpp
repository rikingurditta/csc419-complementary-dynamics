#include "complementary_displacement.h"
#include "rig_jacobian.h"

template<typename ENERGY>
void complementary_displacement(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & T, 
  const Eigen::SparseMatrix<double> & M,
  const Eigen::MatrixXd & Ur,
  const Eigen::MatrixXd & UcLast
  const Eigen::MatrixXd & J,
  const Eigen::MatrixXd & g,
  const Eigen::SparseMatrix<double> & H,
  const double & dt,
  ENERGY & energy,
  Eigen::MatrixXd & Uc){
    Eigen::MatrixXd tmp_g;
    Eigen::SparseMatrix<double> tmp_h;
    newtons_method(UcLast, energy, ..., ..., 10 /* max steps for now i guess lol */,
		   tmp_g, tmp_h);
}


//Input:
//  x0 - initial point for newtons search
//  f(x) - function that evaluates and returns the cost function at x
//  g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//  H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//  max steps - the maximum newton iterations to take
//  tmp_g and tmp_H are scratch space to store gradients and hessians
//Output: 
//  x0 - update x0 to new value
template<typename Objective, typename Jacobian, typename Hessian>
double newtons_method(Eigen::MatrixXd &x0, Objective &f, Jacobian &g, Hessian &H, unsigned int maxSteps,
                      Eigen::MatrixXd &tmp_g, Eigen::SparseMatrix<double> &tmp_H) {
    // get initial gradient
    g(tmp_g, x0);
    int steps = 0;
    while (tmp_g.norm() >= 1e-8 and steps < maxSteps) {
        // get gradient and hessian of objective function
        g(tmp_g, x0);
        H(tmp_H, x0);
        // solve Hd = -g
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(tmp_H);
        Eigen::VectorXd d = solver.solve(-tmp_g);
        // line search for good alpha
        double alpha = 1.;
        while (f(x0 + alpha * d) > f(x0) + 1e-8 * alpha * tmp_g.dot(d) and alpha > 1e-8) {
            alpha = alpha * 0.5;  // scaling factor is 0.5
        }
        // update x0
        x0 = x0 + alpha * d;
        steps++;
    }
    return 0.0;
}
