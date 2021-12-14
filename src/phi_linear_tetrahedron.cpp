#include <phi_linear_tetrahedron.h>

void phi_linear_tetrahedron(Eigen::Vector4d& phi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> x) {
	Eigen::Matrix3d T;
	Eigen::Vector3d x0 = V.row(element[0]).transpose();
	phi.resize(4, 1);

	for (int i = 1; i < 4; i++) {
		T.col(i - 1) = V.row(element[i]).transpose() - x0;
	}

	Eigen::Vector3d prePhi = T.fullPivHouseholderQr().solve(x - x0);
	double phi0 = 1 - prePhi(0) - prePhi(1) - prePhi(2);
	phi(0) = phi0;
	phi.segment<3>(1) = prePhi.transpose();
}