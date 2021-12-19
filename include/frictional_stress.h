#include <Eigen/Dense>
#include <EigenTypes.h>
#include <newtons_method.h>

void frictional_stress(Eigen::MatrixXd& stress, double& mean_stress, double& shear_stress, 
                        Eigen::VectorXd q, Eigen::VectorXd qdot,
                        Eigen::Vector3d p0, Eigen::TensorP& P, Eigen::TensorXV& xv, 
                        Eigen::TensorYV& yv, Eigen::TensorZV& zv, const double dg);