#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
#include <array>

/*

	Resolve the trilinear weights for the corners a cube given the point of interpolation.

	The indexing of the corner CXYZ follows the following binary scheme:
	idx = X*2^2 + Y*2^1 + Z*2^0
	Where X,Y,Z \in {0, 1}, indicating edge movements along each axis.
	For example, the index for the corner C011 is 3 in both the weights array and the corners matrix.

	C011- - - - - - -C111
	| \                 | \
	|   \               |   \
	|     \             |     \
	|       \           |       \
	|         \         |         \
	|           C010- - - - - - -C110
	|           |       |           |
	|           |       |           |
	|           |       |           |
	C001- - - - | - -C101           |
	  \         |         \         |
		\       |           \       |
		  \     |             \     |
			\   |               \   |
			  \ |                 \ |
				C000- - - - - - -C100

	In this diagram X goes left/right, Y goes up/down, Z goes forward/back.

	OUTPUT:
	weights - an 8-long array of weights.

	INPUTS:
	corners - 8x3 matrix of corner coordinates.
	particle - coordinates of interpolation point.

*/
void trilinear_weights(std::array<double, 8> & weights, Eigen::Matrix<double, 8, 3> corners, Eigen::Vector3d particle);


inline int getCornerIndex(int x, int y, int z) {
	return x * (1 << 2) + y * (1 << 1) + z * (1 << 0);
}

inline void getBinaryIndices(Eigen::Vector3i indices, int cornerIndex) {
	for (int i = 0; i < 3; i++) indices.coeffRef(i) = cornerIndex && 1 << i ? 1 : 0;
}

inline void roundVectorDown(Eigen::Vector3i& out, Eigen::Vector3d in) {
	for (int i = 0; i < 3; i++) out.coeffRef(i) = std::floor(in.coeff(i));
}