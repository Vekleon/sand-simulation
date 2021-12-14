#include <assemble_stiffness.h>

void assemble_stiffness(Eigen::SparseMatrixd& K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot,
	Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
	double C, double D)
{
	typedef Eigen::Triplet<double> Trip;
	std::vector<Trip> triplets;
	K.resize(q.rows(), q.rows());
	K.setZero();
	for (int i = 0; i < T.rows(); i++) {
		Eigen::RowVectorXi indices = T.row(i);
		Eigen::Matrix1212d tempK;

		d2V_linear_tetrahedron_dq2(tempK, q, V, indices, v0(i), C, D);
		for (int a = 0; a < 4; a++) {
			int aOffset = 3 * indices(a);
			for (int b = 0; b < 4; b++) {
				int bOffset = 3 * indices(b);
				for (int row = 0; row < 3; row++) {
					for (int col = 0; col < 3; col++) {
						triplets.push_back(Trip(aOffset + row, bOffset + col, -tempK(3 * a + row, 3 * b + col)));
					}
				}
			}
		}
	}

	K.setFromTriplets(triplets.begin(), triplets.end());
};
