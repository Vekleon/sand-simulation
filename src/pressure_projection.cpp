#include <pressure_projection.h>

void pressure_projection(
	Eigen::TensorXV& xv, Eigen::TensorYV& yv, Eigen::TensorZV& zv, Eigen::MatrixXd& pressure,
	Eigen::VectorXd q, Eigen::VectorXd qdot, Eigen::TensorP P,
	const double dx, const double dy, const double dz, const double density){
    // TODO
	Eigen::RowVectorXd B(6);
	B << -1 / dx, 1 / dx, -1 / dy, 1 / dy, -1 / dz, 1 / dz;
	Eigen::MatrixXd D(6, 7);
	D.setZero();
	D.col(4) << 1 / dx, -1 / dx, 1 / dy, -1 / dy, 1 / dz, -1 / dz;
	D(0, 0) = -1 / dx; D(1, 1) = 1 / dx;
	D(2, 2) = -1 / dy; D(3, 4) = 1 / dy;
	D(4, 5) = -1 / dz; D(5, 6) = 1 / dz;



	int xDimension = P.size();
	int yDimension = P.at(0).rows();
	int zDimension = P.at(0).cols();
	// need to change bound limits because 4-7 pressure points can be calculated at a single time ??????
	for (int z = 0; z < zDimension; z++) {
		for (int y = 0; y < yDimension; y++) {
			for (int x = 0; x < xDimension; x++) {
				Eigen::VectorXd qj(6);
				qj << xv.at(x)(y, z), xv.at(x + 1)(y, z), yv.at(x)(y, z), yv.at(x)(y + 1, z), zv.at(x)(y, z), zv.at(x)(y, z + 1);				
			}
		}
	}
}