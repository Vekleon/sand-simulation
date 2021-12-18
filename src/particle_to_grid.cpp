#include <particle_to_grid.h>
#include "trilinear_weights.h"

// Return the greatest distance along any one axist from A to B. Unsigned.
double chebDist(Eigen::Vector3d A, Eigen::Vector3d B) {
	return (A - B).cwiseAbs().maxCoeff();
}

void roundVectorDown(Eigen::Vector3i& out, Eigen::Vector3d in) {
	for (int i = 0; i < 3; i++) out.coeffRef(i) = std::floor(in.coeff(i));
}

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

	// Prep work
	const int n = q.size() / 3;
	const Eigen::Vector3d ONE_HALF = Eigen::Vector3d::Constant(0.5);
	const Eigen::Vector3d X_GRID_OFFSET = Eigen::Vector3d(0, 0, 0.5);
	const Eigen::Vector3d Y_GRID_OFFSET = Eigen::Vector3d(0, 0.5, 0);
	const Eigen::Vector3d Z_GRID_OFFSET = Eigen::Vector3d(0.5, 0, 0);
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

		// X VELOCITY GRID
		// Index a cell of the X component velocity grid by indexing the corner closest to the origin
		roundVectorDown(cell_idx, (1. / dg) * (particle_pos - p0) - X_GRID_OFFSET);
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				for (int k = 0; k < 2; k++) {
					corners_matrix.row(getCornerIndex(i, j, k)) = dg * (cell_idx.cast<double>() + Eigen::Vector3d(i, j, k) + X_GRID_OFFSET) + p0;
				}
			}
		}
		trilinear_weights(weights, corners_matrix, particle_pos);
		Eigen::Vector3i corner_idx;
		for (int wi = 0; wi < 8; wi++) {
			int i, j, k;
			getBinaryIndices(corner_idx, wi);
			tensorCoeffRef(xv, corner_idx + cell_idx) += weights.at(wi) * qdot(3 * pi + 0);
			tensorCoeffRef(count_x, corner_idx + cell_idx) += weights.at(wi);
		}

		//for (int x_offset = 0; x_offset <= 1; x_offset++) {
		//	for (int z_offset = 0; z_offset <= 1; z_offset++) {
		//		
		//		// Current position in the x velocity grid
		//		Eigen::Vector3i corner_idx;
		//		corner_idx << cell_idx(0) + x_offset, cell_idx(1), cell_idx(2) + z_offset;

		//		// Coordinates of current corner
		//		Eigen::Vector3d corner_pos = dg * (corner_idx.cast<double>() - ONE_HALF);

		//		// Resolve weight for this addition
		//		const double dist = chebDist(corner_pos, particle_pos);
		//		const double weight = 1 - (dist / dg);

		//		// Add with weight and track running weight total
		//		xv.at(corner_idx(0))(corner_idx(1), corner_idx(2)) += qdot(pi * 3 + 0) * weight;
		//		count_x.at(corner_idx(0))(corner_idx(1), corner_idx(2)) += weight;
		//	}
		//}

		// TODO: Y AND Z VELOCITY GRIDS
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