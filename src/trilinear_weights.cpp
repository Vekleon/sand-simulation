#include <trilinear_weights.h>

void trilinear_weights(std::array<double, 8> & weights, Eigen::Matrix<double, 8, 3> corners, Eigen::Vector3d particle) {

	/*
	
	Using the Wikipedia page for trilinear interpolation: https://en.wikipedia.org/wiki/Trilinear_interpolation#Method

	Fully expanding the term c and collecting all corner terms we get: 

	c = 
		c000*(xd*yd - yd - zd - xd + xd*zd + yd*zd - xd*yd*zd + 1) +
		c001*(zd - xd*zd - yd*zd + xd*yd*zd) + 
		c010*(yd - xd*yd - yd*zd + xd*yd*zd) + 
		c011*(yd*zd - xd*yd*zd) + 
		c100*(xd - xd*yd - xd*zd + xd*yd*zd) + 
		c101*(xd*zd - xd*yd*zd) + 
		c110*(xd*yd - xd*yd*zd) + 
		c111*xd*yd*zd

	*/

	// Getting the important corners
	const Eigen::RowVector3d c000 = corners.row(0);
	const Eigen::RowVector3d c001 = corners.row(1);
	const Eigen::RowVector3d c010 = corners.row(2);
	const Eigen::RowVector3d c100 = corners.row(4);

	const double xd = (particle(0) - c000(0)) / (c100(0) - c000(0));
	const double yd = (particle(1) - c000(1)) / (c010(1) - c000(1));
	const double zd = (particle(2) - c000(2)) / (c001(2) - c000(2));

	weights.at(0) = xd * yd - yd - zd - xd + xd * zd + yd * zd - xd * yd * zd + 1;
	weights.at(1) = zd - xd * zd - yd * zd + xd * yd * zd;
	weights.at(2) = yd - xd * yd - yd * zd + xd * yd * zd;
	weights.at(3) = yd * zd - xd * yd * zd;
	weights.at(4) = xd - xd * yd - xd * zd + xd * yd * zd;
	weights.at(5) = xd * zd - xd * yd * zd;
	weights.at(6) = xd * yd - xd * yd * zd;
	weights.at(7) = xd * yd * zd;
}