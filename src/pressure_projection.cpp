#include <pressure_projection.h>
#include <math.h>

void pressure_projection(
	Eigen::TensorXV& xv, Eigen::TensorYV& yv, Eigen::TensorZV& zv, Eigen::TensorP& pressure,
	std::set<int>& pressure_indices, Eigen::VectorXd& q, Eigen::VectorXd& qdot, Eigen::TensorPB& P,
	const double dg, const double density, const double dt){
    // TODO

	const double dg_inv = 1. / dg;
	const double rho_over_dt = density / dt;
	Eigen::RowVectorXd B(6);
	B << -dg_inv, dg_inv, -dg_inv, dg_inv, -dg_inv, dg_inv;
	Eigen::MatrixXd D(6, 7);
	D.setZero();
	for (int i = 0; i < 6; i++) {
		D.coeffRef(i, i) = dg_inv;
		D.coeffRef(i, 6) = dg_inv;
		if (i % 2 == 0) {
			D.coeffRef(i, i) *= -1.;
		}
		else {
			D.coeffRef(i, 6) *= -1.;
		}
	}
	std::array<Eigen::Vector3i, 7> corner_coord_offsets;
	corner_coord_offsets.at(0) = Eigen::Vector3i(-1, 0, 0);
	corner_coord_offsets.at(1) = Eigen::Vector3i(+1, 0, 0);
	corner_coord_offsets.at(2) = Eigen::Vector3i(0, -1, 0);
	corner_coord_offsets.at(3) = Eigen::Vector3i(0, +1, 0);
	corner_coord_offsets.at(4) = Eigen::Vector3i(0, 0, -1);
	corner_coord_offsets.at(5) = Eigen::Vector3i(0, 0, +1);
	corner_coord_offsets.at(6) = Eigen::Vector3i::Zero();

	const int xDimension = pressure.size();
	const int yDimension = pressure.at(0).rows();
	const int zDimension = pressure.at(0).cols();
	Eigen::Matrix<double, 6, 1> qj;
	Eigen::VectorXd d(xDimension * yDimension * zDimension);
	Eigen::SparseMatrixd A;
	Eigen::Matrix<double, 6, 6> PTP;
	Eigen::Matrix<double, 1, 7> Aj;
	std::vector<Eigen::Triplet<double>> triplets;

	auto get_p_idx = [&xDimension, &yDimension, &zDimension](int x, int y, int z) {
		return get_cell_idx(x, y, z, xDimension, yDimension, zDimension);
	};

	// need to change bound limits because 4-7 pressure points can be calculated at a single time ??????
	for (int z = 0; z < zDimension; z++) {
		for (int y = 0; y < yDimension; y++) {
			for (int x = 0; x < xDimension; x++) {
				const int cur_row = get_p_idx(x, y, z);

				qj <<
					tensorAtOrZero(xv, x, y, z),
					tensorAtOrZero(xv, x + 1, y, z),
					tensorAtOrZero(yv, x, y, z),
					tensorAtOrZero(yv, x, y + 1, z),
					tensorAtOrZero(zv, x, y, z),
					tensorAtOrZero(zv, x, y, z + 1);
				d.coeffRef(cur_row) = rho_over_dt * B.dot(qj);

				// Boundary consitions
				PTP = P.at(x).at(y).at(z);
				// A copy of D but accounting for ghost pressures
				Eigen::Matrix<double, 6, 7> Dj = D;
				// a bit array to zero to figure out which columns we have to zero out
				// not sure if this would be better than using a vector 
				char zeros = 0;

				if(pressure_indices.find(get_p_idx(x - 1, y, z)) != pressure_indices.end()) zeros += 1;
				if(pressure_indices.find(get_p_idx(x + 1, y, z)) != pressure_indices.end()) zeros += 2;
				if(pressure_indices.find(get_p_idx(x, y - 1, z)) != pressure_indices.end()) zeros += 4;
				if(pressure_indices.find(get_p_idx(x, y + 1, z)) != pressure_indices.end()) zeros += 8;
				if(pressure_indices.find(get_p_idx(x, y, z - 1)) != pressure_indices.end()) zeros += 16;
				if(pressure_indices.find(get_p_idx(x, y, z + 1)) != pressure_indices.end()) zeros += 32;
				if(pressure_indices.find(get_p_idx(x, y, z)) != pressure_indices.end()) zeros += 64;

				for(int i = 1; i < 8; i++){
					int result = pow(2, i);
					if(result & zeros) Dj.col(i - 1).setZero();
				}

				Aj = B * PTP * Dj; // TODO: GHOST PRESSURES

				// Assemble the current row of A
				for (int i = 0; i < corner_coord_offsets.size(); i++) {
					int cur_x = x + corner_coord_offsets.at(i)(0);
					int cur_y = y + corner_coord_offsets.at(i)(1);
					int cur_z = z + corner_coord_offsets.at(i)(2);
					// It's possible to go out of bounds so let's not do that
					if (!validCoordinates(pressure, cur_x, cur_y, cur_z)) continue;
					int cur_col = get_p_idx(cur_x, cur_y, cur_z);
					triplets.push_back(Eigen::Triplet<double>(cur_row, cur_col, Aj(i)));
				}
			}
		}
	}

	A.resize(xDimension * yDimension * zDimension, xDimension * yDimension * zDimension);
	A.setFromTriplets(triplets.begin(), triplets.end());

	Eigen::ConjugateGradient<Eigen::SparseMatrixd> solver;
	solver.compute(A);
	assert(solver.info() == Eigen::Success);
	Eigen::VectorXd p = solver.solve(d);
	assert(solver.info() == Eigen::Success);

	// I think this is how we translate the p vector into the pressure grid, if we're staying consistent with how we assembled A and d
	for (int z = 0; z < zDimension; z++) {
		for (int y = 0; y < yDimension; y++) {
			for (int x = 0; x < xDimension; x++) {
				tensorAt(pressure, x, y, z) = p(get_p_idx(x, y, z));
			}
		}
	}

	// Compute v^{t+1} on the grid
	const double dt_over_rho = dt / density;
	for (int z = 0; z < zDimension + 1; z++) {
		for (int y = 0; y < yDimension + 1; y++) {
			for (int x = 0; x < xDimension + 1; x++) {
				// UPDATE X
				if (validCoordinates(xv, x, y, z)) {
					tensorAt(xv, x, y, z) -= dt_over_rho * dg_inv * (tensorAtOrZero(pressure, x, y, z) - tensorAtOrZero(pressure, x - 1, y, z));
				}
				// UPDATE Y
				if (validCoordinates(yv, x, y, z)) {
					tensorAt(yv, x, y, z) -= dt_over_rho * dg_inv * (tensorAtOrZero(pressure, x, y, z) - tensorAtOrZero(pressure, x, y - 1, z));
				}
				// UPDATE Z
				if (validCoordinates(zv, x, y, z)) {
					tensorAt(zv, x, y, z) -= dt_over_rho * dg_inv * (tensorAtOrZero(pressure, x, y, z) - tensorAtOrZero(pressure, x, y, z - 1));
				}
			}
		}
	}
}