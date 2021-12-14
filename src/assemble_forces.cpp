#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd& f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot,
	Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
	double C, double D)
{
	f.resize(q.rows());
	f.setZero();
	for (int i = 0; i < T.rows(); i++) {
		Eigen::RowVector4i indices = T.row(i);
		Eigen::Vector12d tempF;
		dV_linear_tetrahedron_dq(tempF, q, V, indices, v0(i), C, D);

		for (int index = 0; index < 4; index++)
			f.segment<3>(3 * indices(index)) -= tempF.segment<3>(3 * index);
	}
};