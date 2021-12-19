#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/*

Perform PIC to transver grid velocities back onto the particles

*/

void particle_to_grid(Eigen::TensorXV& xv, Eigen::TensorYV& yv, Eigen::TensorZV& zv,
	Eigen::Vector3d p0, double dg, Eigen::VectorXd &q, Eigen::VectorXd &qdot);