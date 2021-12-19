#include <pressure_projection.h>

void pressure_projection(
	Eigen::TensorXV& xv, Eigen::TensorYV& yv, Eigen::TensorZV& zv, Eigen::TensorP& pressure,
	Eigen::VectorXd& q, Eigen::VectorXd& qdot, Eigen::TensorPB& P,
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

	const int xDimension = pressure.size();
	const int yDimension = pressure.at(0).rows();
	const int zDimension = pressure.at(0).cols();
	Eigen::Matrix<double, 6, 1> qj;
	Eigen::VectorXd d(xDimension * yDimension * zDimension);
	Eigen::SparseMatrixd A;
	Eigen::Matrix<double, 6, 6> PTP;
	Eigen::Matrix<double, 6, 7> Dj; // A copy of D but accounting for ghost pressures
	Eigen::Matrix<double, 1, 7> Aj;
	std::vector<Eigen::Triplet<double>> triplets;

	auto get_p_idx = [&xDimension, &yDimension, &zDimension](int x, int y, int z) {
		return get_cell_idx(x, y, z, xDimension, yDimension, zDimension);
	};

	// need to change bound limits because 4-7 pressure points can be calculated at a single time ??????
	for (int z = 0; z < zDimension; z++) {
		for (int y = 0; y < yDimension; y++) {
			for (int x = 0; x < xDimension; x++) {

				qj <<
					tensorAtOrZero(xv, x, y, z),
					tensorAtOrZero(xv, x + 1, y, z),
					tensorAtOrZero(yv, x, y, z),
					tensorAtOrZero(yv, x, y + 1, z),
					tensorAtOrZero(zv, x, y, z),
					tensorAtOrZero(zv, x, y, z + 1);
				d.coeffRef(get_p_idx(x, y, z)) = rho_over_dt * B.dot(qj);

				// Boundary consitions
				PTP = P.at(x).at(y).at(z);
				Aj = B * PTP * Dj; // TODO: GHOST PRESSURES

				// I think it's supposed to look something like this??? No idea what coefficient we put in here though.
				triplets.push_back(Eigen::Triplet<double>(get_p_idx(x - 1, y, z), get_p_idx(x, y, z), Aj(0)));
				triplets.push_back(Eigen::Triplet<double>(get_p_idx(x, y, z), get_p_idx(x + 1, y, z), Aj(1)));
				triplets.push_back(Eigen::Triplet<double>(get_p_idx(x, y - 1, z), get_p_idx(x, y, z), Aj(2)));
				triplets.push_back(Eigen::Triplet<double>(get_p_idx(x, y, z), get_p_idx(x, y + 1, z), Aj(3)));
				triplets.push_back(Eigen::Triplet<double>(get_p_idx(x, y, z - 1), get_p_idx(x, y, z), Aj(4)));
				triplets.push_back(Eigen::Triplet<double>(get_p_idx(x, y, z), get_p_idx(x, y, z + 1), Aj(5)));
				triplets.push_back(Eigen::Triplet<double>(get_p_idx(x, y, z), get_p_idx(x, y, z),     Aj(6)));
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
}