#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  x0 - initial point for newtons search
//  f(x) - function that evaluates and returns the cost function at x
//  g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//  H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//  max steps - the maximum newton iterations to take
//  tmp_g and tmp_H are scratch space to store gradients and hessians
//Output:
//  x0 - update x0 to new value
template <typename Objective, typename Jacobian, typename Hessian>
double newtons_method(Eigen::VectorXd& x0, Objective& f, Jacobian& g, Hessian& H, unsigned int maxSteps, Eigen::VectorXd& tmp_g, Eigen::SparseMatrixd& tmp_H)
{
    for (int i = 0; i < maxSteps; i++) {
        g(tmp_g, x0);

        if (tmp_g.norm() < 1e-3)
            return 0;

        double energy = f(x0);
        H(tmp_H, x0);

        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(tmp_H);
        Eigen::VectorXd d = solver.solve(-1. * tmp_g);
        double alpha = 1.0;
        double threshold = energy + 0.75 * d.transpose() * tmp_g;

        while (true) {
            if (f(x0 + 0.01 * alpha * d) <= threshold || alpha < 1e-6) {
                break;
            }
            else {
                alpha *= 0.01;
            }
        }

        x0 = x0 + alpha * d;
    }

    return 0.0;
}