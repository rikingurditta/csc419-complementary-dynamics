#include "complementary_displacement.h"
#include "mass_matrix.h"
//#include "gradient.h"
//#include "hessian.h"
#include "rig_jacobian.h"

template<typename ENERGY, typename Gradient, typename Hessian>
void complementary_displacement(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & T, 
  const Eigen::SparseMatrix<double> & M,
  const Eigen::MatrixXd & Ur,
  const Eigen::MatrixXd & UcLast,
  const Eigen::MatrixXd & ULast,
  const Eigen::MatrixXd & duLast,
  const Eigen::MatrixXd & J,
  const Eigen::MatrixXd & g,
  const Eigen::SparseMatrix<double> & H,
  const double & dt,
  ENERGY & energy,
  Gradient & gradient,
  Hessian & hessian,
  Eigen::MatrixXd & Uc){
    Eigen::MatrixXd tmp_g;
    Eigen::SparseMatrix<double> tmp_h;
    Uc = UcLast;
    newtons_method(Uc, energy, gradient, hessian, 10, dt, Ur, duLast, ULast, J, //placeholder value for the last derivative and last U value
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
void newtons_method(Eigen::MatrixXd &x0,
		    Objective &f,
		    Jacobian &g,
		    Hessian &H,
		    unsigned int maxSteps,
		    double dt,
		    Eigen::MatrixXd Ur,
		    Eigen::VectorXd lastDut,
		    Eigen::MatrixXd lastDu,
		    Eigen::MatrixXd & J,
		    Eigen::SparseMatrixd &M,
		    Eigen::MatrixXd rigJacobian,
		    Eigen::MatrixXd &tmp_g,
		    Eigen::SparseMatrix<double> &tmp_H) {
    // get initial gradient
    g(tmp_g, x0);
    int steps = 0;
    while (tmp_g.norm() >= 1e-8 and steps < maxSteps) {
        // get gradient and hessian of objective function
        g(tmp_g, x0);
        H(tmp_H, x0); // (would be K in the pseudo code)
	Eigen::MatrixXd Q = tmp_H + dt * dt * M;
	Eigen::VectorXd l = -tmp_g + 1 / dt * M * ((Ur - lastDu)/ dt - lastDut);
	Eigen::MatrixXd C = J.transpose() * M;
	Eigen::MatrixXd block = Eigen::MatrixXd::Zero(Q.cols() + C.rows(),
						      Q.cols() + C.rows());
	// construct the block matrix
	block.block(0, 0, Q.rows(), Q.cols()) = Q;
	block.block(0, Q.cols(), C.cols(), C.rows()) = C.transpose();
	block.block(Q.rows(), 0, C.rows(), C.cols()) = C;
        // solve Hd = -g
	Eigen::PartialPivLU<Eigen::MatrixXd> solver(block);
	Eigen::VectorXd sol = Eigen::VectorXd::Zero(block.cols());
	sol.head(l.size()) = l;
	Eigen::VectorXd tmp = solver.inverse() * sol;
	// ignore last value to get our x value which we store in d
	Eigen::VectorXd d = tmp.head(tmp.size() - 1);
        // line search for good alpha
        double alpha = 1.;
	// TODO: double check, will also need with the Uc and Ur thing
        while (f(x0 + alpha * d + Ur) > f(x0 + Ur) + 1e-8 * alpha * tmp_g.dot(d) and alpha > 1e-8) {
            alpha = alpha * 0.5;  // scaling factor is 0.5
        }
        // update x0
        x0 = x0 + alpha * d;
        steps++;
    }
}
