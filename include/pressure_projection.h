#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/*
	Performs pressure projection

	xv: axb grid of x velocities
	yv: axb grid of y velocities
	zv: axb grid of z velocities
	pressure: (a-1)x(b-1) grid of pressures
	q: (n*3)x1 generalized coordinates
	qdot: (n*3)x1 generalized velocity
	P: ((a-1)*4)x((a-1)*4) matrix representing boundary conditions. The block at (a*4,b*4) represents the boundary conditions at pressure point (a,b)
	dx: distance along axis between positions on xv
	dy: distance along axis between positions on yv
	dz: distance along axis between positions on zv
*/
void pressure_projection(
	Eigen::MatrixXd &xv, Eigen::MatrixXd& yv, Eigen::MatrixXd& zv, Eigen::MatrixXd &pressure,
	Eigen::VectorXd q, Eigen::VectorXd qdot, Eigen::MatrixXd P,
	const double dx, const double dy);