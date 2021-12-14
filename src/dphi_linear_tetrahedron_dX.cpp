#include <dphi_linear_tetrahedron_dX.h>
#include <phi_linear_tetrahedron.h>
#include <iostream>
void dphi_linear_tetrahedron_dX(Eigen::Matrix43d& dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
	Eigen::Matrix3d T;
	Eigen::Vector3d x0 = V.row(element(0)).transpose();
	for (int i = 1; i < element.size(); i++)
		T.col(i - 1) = V.row(element(i)).transpose() - x0;

	Eigen::Matrix3d Inv = T.inverse();
	dphi.block<3, 3>(1, 0) = Inv;
	dphi.block<1, 3>(0, 0) = (-Inv.row(0).transpose() - Inv.row(1).transpose() - Inv.row(2).transpose());
}