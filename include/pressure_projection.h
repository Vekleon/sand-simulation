#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
#include <set>

/*
	Performs pressure projection

	xv: (X+1) * Y * Z tensor of x-component velocities
	yv: X * (Y+1) * Z tensor of y-component velocities
	zv: X * Y * (Z+1) tensor of z-component velocities
	pressure: (a-1)x(b-1) grid of pressures
	q: (n*3)x1 generalized coordinates
	qdot: (n*3)x1 generalized velocity
	P: ((a-1)*4)x((a-1)*4) matrix representing boundary conditions. The block at (a*4,b*4) represents the boundary conditions at pressure point (a,b)
	dx: distance along axis between positions on xv
	dy: distance along axis between positions on yv
	dz: distance along axis between positions on zv
*/
void pressure_projection(
	Eigen::TensorXV& xv, Eigen::TensorYV& yv, Eigen::TensorZV& zv, Eigen::TensorP& pressure,
	std::set<int>& pressureIndices, Eigen::VectorXd& q, Eigen::VectorXd& qdot, Eigen::TensorPB& P,
	const double dg, const double density, const double dt);