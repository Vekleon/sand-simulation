#include <mass_matrix_linear_tetrahedron.h>

void mass_matrix_linear_tetrahedron(Eigen::Matrix1212d& M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
	M.setZero();
	double sames = volume * density / 10.0;
	double diffs = volume * density / 20.0;
	Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

	for (int y = 0; y < 4; y++) {
		for (int x = 0; x < 4; x++) {
			if (x == y)
				M.block(3 * y, 3 * x, 3, 3) = I * sames;
			else
				M.block(3 * y, 3 * x, 3, 3) = I * diffs;
		}
	}
}