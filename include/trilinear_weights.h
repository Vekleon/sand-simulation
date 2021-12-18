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
	|           C001- - - - - - -C101
	|           |       |           |
	|           |       |           |
	|           |       |           |
	C010- - - - | - -C110           |
	  \         |         \         |
		\       |           \       |
		  \     |             \     |
			\   |               \   |
			  \ |                 \ |
				C000- - - - - - -C100

	In this diagram X goes left/right, Y goes forward/back, Z goes up/down.

	OUTPUT:
	weights - an 8-long array of weights.

	INPUTS:
	corners - 8x3 matrix of corner coordinates.
	particle - coordinates of interpolation point.

*/
void trilinear_weights(std::array<double, 8> & weights, Eigen::Matrix<float, 8, 3> corners, Eigen::Vector3d particle);