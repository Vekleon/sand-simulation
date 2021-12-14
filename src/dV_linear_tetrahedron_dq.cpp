#include <dV_linear_tetrahedron_dq.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>

void dV_linear_tetrahedron_dq(Eigen::Vector12d& dV, Eigen::Ref<const Eigen::VectorXd> q,
	Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
	double C, double D) {

	auto neohookean_linear_tet = [&](Eigen::Vector12d& dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
		Eigen::Matrix3d F;
		Eigen::Matrix43d phi;
		Eigen::Matrix34d T;
		Eigen::Vector9d psi;

		for (int i = 0; i < element.size(); i++)
			T.col(i) = q.segment<3>(3 * element(i));

		dphi_linear_tetrahedron_dX(phi, V, element, X);
		F = T * phi;

		Eigen::Matrix<double, 9, 12> B;
		B.setZero();
		for (int i = 0; i < 4; i++) {
			B.block(0, 3 * i, 3, 1) = phi.row(i).transpose();
			B.block(3, 3 * i + 1, 3, 1) = phi.row(i).transpose();
			B.block(6, 3 * i + 2, 3, 1) = phi.row(i).transpose();
		}

		dpsi_neo_hookean_dF(psi, F, C, D);
		dV = B.transpose() * psi;
	};

	quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);

}