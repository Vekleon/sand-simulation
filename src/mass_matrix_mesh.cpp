#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>

void mass_matrix_mesh(Eigen::SparseMatrixd& M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0) {
	M.resize(qdot.size(), qdot.size());
	M.setZero();
	typedef Eigen::Triplet<double> Trip;
	std::vector<Trip> triplets;

	for (int i = 0; i < T.rows(); i++) {
		Eigen::Matrix1212d tempM;
		Eigen::RowVectorXi indices = T.row(i);
		mass_matrix_linear_tetrahedron(tempM, qdot, indices, density, v0(i));

		// Traverse the mass matrix
		for (int a = 0; a < 4; a++) {
			int aOffset = 3 * indices(a);
			for (int b = 0; b < 4; b++) {
				int bOffset = 3 * indices(b);
				for (int row = 0; row < 3; row++) {
					for (int col = 0; col < 3; col++) {
						triplets.push_back(Trip(aOffset + row, bOffset + col, tempM(3 * a + row, 3 * b + col)));
					}
				}
			}
		}
	}

	M.setFromTriplets(triplets.begin(), triplets.end());
}