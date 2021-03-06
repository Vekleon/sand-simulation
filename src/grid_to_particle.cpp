#include <grid_to_particle.h>
#include <trilinear_weights.h>

void particle_to_grid(Eigen::TensorXV& xv, Eigen::TensorYV& yv, Eigen::TensorZV& zv,
	Eigen::Vector3d p0, double dg, Eigen::VectorXd& q, Eigen::VectorXd& qdot) {

	/*
	Each col of this matrix represents the offset from p0 to the point on a velocity grid which is closest to the origin.
	For example, the (0, 0, 0) point on the X-Velocity grid is equal to p0 - <0.5, 0.0, 0.0>.
	*/
	const Eigen::Matrix3d GRID_OFFSETS = 0.5 * Eigen::Matrix3d::Identity();
	const int n = q.size() / 3;

	Eigen::Matrix<double, 8, 3> corners_matrix;
	std::array<double, 8> weights;
	Eigen::Vector3i cell_idx;
	Eigen::VectorXd total_weights(qdot.size());
	total_weights.setConstant(1e-32);
	qdot.setZero();

	// For each particle P
	for (int pi = 0; pi < n; pi++) {

		Eigen::Vector3d particle_pos = q.segment<3>(pi * 3);

		// For each grid G in the staggered grid...
		for (int gi = 0; gi < 3; gi++) {

			// Index the current cell by identifying the corner closest to the origin
			roundVectorDown(cell_idx, (1. / dg) * (particle_pos - p0) + GRID_OFFSETS.col(gi));
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						corners_matrix.row(getCornerIndex(i, j, k)) = dg * (cell_idx.cast<double>() + Eigen::Vector3d(i, j, k) - GRID_OFFSETS.col(gi)) + p0;
					}
				}
			}

			// Resolve weights and apply sum, keeping track of weights
			trilinear_weights(weights, corners_matrix, particle_pos);
			Eigen::Vector3i corner_idx;
			for (int wi = 0; wi < 8; wi++) {
				getBinaryIndices(corner_idx, wi);
				switch (gi) {
					case 0:
						if (!validCoordinates(xv, corner_idx + cell_idx)) continue;
						qdot.coeffRef(3 * pi + gi) += weights.at(wi) * tensorAt(xv, corner_idx + cell_idx);
						break;
					case 1:
						if (!validCoordinates(yv, corner_idx + cell_idx)) continue;
						qdot.coeffRef(3 * pi + gi) += weights.at(wi) * tensorAt(yv, corner_idx + cell_idx);
						break;
					case 2:
						if (!validCoordinates(zv, corner_idx + cell_idx)) continue;
						qdot.coeffRef(3 * pi + gi) += weights.at(wi) * tensorAt(zv, corner_idx + cell_idx);
						break;
					default:
						assert(false);
						break;
				}
				total_weights.coeffRef(3 * pi + gi) += weights.at(wi);
			}
		}
	}

	qdot = qdot.cwiseQuotient(total_weights);
}