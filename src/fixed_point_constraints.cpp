#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
	P.resize(q_size - 3 * indices.size(), q_size);
	P.setZero();
	std::vector<Eigen::Triplet<double>> triplets;

	int row = 0;
	for (int i = 0; i < q_size; i += 3) {
		triplets.push_back(Eigen::Triplet<double>(row, i, 1.));
		triplets.push_back(Eigen::Triplet<double>(row + 1, i + 1, 1.));
		triplets.push_back(Eigen::Triplet<double>(row + 2, i + 2, 1.));
		row += 3;
	}

	P.setFromTriplets(triplets.begin(), triplets.end());
}