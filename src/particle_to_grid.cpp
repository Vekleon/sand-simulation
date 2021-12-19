#include <particle_to_grid.h>
#include <trilinear_weights.h>

void particle_to_grid(Eigen::TensorXV& xv, Eigen::TensorYV& yv, Eigen::TensorZV& zv,
	Eigen::Vector3d p0, double dg, Eigen::VectorXd q, Eigen::VectorXd qdot, 
	Eigen::TensorXV &count_x, Eigen::TensorYV &count_y, Eigen::TensorZV &count_z) {
    
	/*
	
	Initialize all the velocity grids to 0
	Create count grids parallel to each grid in the staggered grid, initialize them to zero
		...call them Count_X, Count_Y, Count_Z

	For each particle P:
		For each grid G in the staggered grid:
			Cs = Corners of cube which P.pos is in
			For each corner C in Cs:
				DIST <- distance from P.pos to C
				WEIGHT <- (1 - (DIST / dG))
				G at C += P.velocity[G] * WEIGHT
				Count_G at C += WEIGHT

	For each G:
		Everywhere that Count_G is nonzero, elementwise divide G by Count_G
	
	*/

	const int n = q.size() / 3;

	/*
	Each row of this matrix represents the offset from p0 to the point on a velocity grid which is closest to the origin.
	For example, the (0, 0, 0) point on the X-Velocity grid is equal to p0 - <0.5, 0.5, 0.0>.
	*/
	Eigen::Matrix3d GRID_OFFSETS;
	GRID_OFFSETS <<
		0.5, 0.5, 0.0,
		0.0, 0.5, 0.5,
		0.5, 0.0, 0.5;

	// We initialize the weight counts at a small but positive number 
	// so we don't have to check for divide-by-zero errors later
	for (int i = 0; i < xv.size(); i++) {
		xv.at(i).setZero();
		count_x.at(i).setConstant(1e-32);
	}
	for (int i = 0; i < yv.size(); i++) {
		yv.at(i).setZero();
		count_x.at(i).setConstant(1e-32);
	}
	for (int i = 0; i < zv.size(); i++) {
		zv.at(i).setZero();
		count_x.at(i).setConstant(1e-32);
	}

	Eigen::Matrix<double, 8, 3> corners_matrix;
	std::array<double, 8> weights;
	Eigen::Vector3i cell_idx;
	
	// For each particle P
	for (int pi = 0; pi < n; pi++) {

		Eigen::Vector3d particle_pos = q.segment<3>(pi * 3);

		// For each grid G in the staggered grid...
		for (int gi = 0; gi < 3; gi++) {

			// Index the current cell by identifying the corner closest to the origin
			roundVectorDown(cell_idx, (1. / dg) * (particle_pos - p0) - GRID_OFFSETS.row(gi).transpose());
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						corners_matrix.row(getCornerIndex(i, j, k)) = dg * (cell_idx.cast<double>() + Eigen::Vector3d(i, j, k) + GRID_OFFSETS.row(gi).transpose()) + p0;
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
						tensorAt(xv, corner_idx + cell_idx) += weights.at(wi) * qdot(3 * pi + gi);
						tensorAt(count_x, corner_idx + cell_idx) += weights.at(wi);
						break;
					case 1:
						tensorAt(yv, corner_idx + cell_idx) += weights.at(wi) * qdot(3 * pi + gi);
						tensorAt(count_y, corner_idx + cell_idx) += weights.at(wi);
						break;
					case 2:
						tensorAt(zv, corner_idx + cell_idx) += weights.at(wi) * qdot(3 * pi + gi);
						tensorAt(count_z, corner_idx + cell_idx) += weights.at(wi);
						break;
					default:
						assert(false);
						break;
				}
			}
		}
	}

	// Divide by weights now...might be very slow :(
	// This is where that small-but-positive initial weight comes in handy
	for (int i = 0; i < xv.size(); i++) {
		xv.at(i) = xv.at(i).cwiseQuotient(count_x.at(i));
	}
	for (int i = 0; i < yv.size(); i++) {
		yv.at(i) = yv.at(i).cwiseQuotient(count_y.at(i));
	}
	for (int i = 0; i < zv.size(); i++) {
		zv.at(i) = zv.at(i).cwiseQuotient(count_z.at(i));
	}
}