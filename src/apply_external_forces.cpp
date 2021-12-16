#include <apply_external_forces.h>

void apply_external_forces(Eigen::VectorXd& qdot, const double dt) {
	const int n = qdot.size() / 3;
	for (int i = 0; i < n; i++) {
		qdot.segment<3>(i * 3) += dt * Eigen::Vector3d(0., -9.8, 0.);
	}
}