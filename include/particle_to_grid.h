#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/*
	Maps particle velocities to grid positions.

	OUTPUT:
	xv: (X+1) * Y * Z tensor of x-component velocities
	yv: X * (Y+1) * Z tensor of y-component velocities
	zv: X * Y * (Z+1) tensor of z-component velocities

	THE FOLLOWING VARIABLES ARE USED FOR RESOLVING ALL POSITIONS WIHTIN THE GRID
	dg: distance between points on the same grid
	p0: corner of the pressure grid closest to the origin.

	q: (n*3)x1 generalized coordinates
	qdot: (n*3)x1 generalized velocity

	weight count tensors provided as scratch space since they 
	take a lot of space and we don't want to blow the stack.
*/
void particle_to_grid(Eigen::TensorXV& xv, Eigen::TensorYV& yv, Eigen::TensorZV& zv,
	Eigen::Vector3d p0, const double dg, Eigen::VectorXd q, Eigen::VectorXd qdot,
	Eigen::TensorXV& count_x, Eigen::TensorYV& count_y, Eigen::TensorZV& count_z);