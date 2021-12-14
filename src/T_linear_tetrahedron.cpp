#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>
void T_linear_tetrahedron(double& T, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
	Eigen::Matrix1212d M;
	mass_matrix_linear_tetrahedron(M, qdot, element, density, volume);

	Eigen::Vector12d q;
	for (int i = 0; i < 4; i++)
		q.segment<3>(3 * i) = qdot.segment<3>(3 * element(i));
	T = 0.5 * q.transpose() * M * q;
}