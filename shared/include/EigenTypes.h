#ifndef EIGENTYPES_H
#define EIGENTYPES_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>

#define TENSOR_P_X 100
#define TENSOR_P_Y 100
#define TENSOR_P_Z 100

#define TENSOR_XV_X (TENSOR_P_X + 1)
#define TENSOR_XV_Y TENSOR_P_Y
#define TENSOR_XV_Z TENSOR_P_Z

#define TENSOR_YV_X TENSOR_P_X
#define TENSOR_YV_Y (TENSOR_P_Y + 1)
#define TENSOR_YV_Z TENSOR_P_Z

#define TENSOR_ZV_X TENSOR_P_X
#define TENSOR_ZV_Y TENSOR_P_Y
#define TENSOR_ZV_Z (TENSOR_P_Z + 1)

namespace Eigen {


    //dense types
    using Vector4d = Eigen::Matrix<double,4,1>;
    using Vector6d = Eigen::Matrix<double, 6,1>;
    using Vector9d = Eigen::Matrix<double, 9, 1>;
    using Vector12d = Eigen::Matrix<double, 12,1>;
    
    using Matrix36d = Eigen::Matrix<double, 3,6>;
    using Matrix34d = Eigen::Matrix<double, 3,4>;
    using Matrix43d = Eigen::Matrix<double, 4,3>;
    using Matrix66d  = Eigen::Matrix<double, 6,6>;
    using Matrix99d = Eigen::Matrix<double, 9, 9>;
    using Matrix1212d = Eigen::Matrix<double, 12,12>;
    
    using Matrix44f = Eigen::Matrix<float, 4,4>;
    
    //sparse types
    using SparseMatrixd = Eigen::SparseMatrix<double>;

	// Forgive us for this we had no other choice lmao
	using TensorP = std::array<Eigen::Matrix<float, TENSOR_P_Y, TENSOR_P_Z>, TENSOR_P_X>;
	using TensorXV = std::array<Eigen::Matrix<float, TENSOR_XV_Y, TENSOR_XV_Z>, TENSOR_XV_X>;
	using TensorYV = std::array<Eigen::Matrix<float, TENSOR_YV_Y, TENSOR_YV_Z>, TENSOR_YV_X>;
	using TensorZV = std::array<Eigen::Matrix<float, TENSOR_ZV_Y, TENSOR_ZV_Z>, TENSOR_ZV_X>;

}

template <int X, int Y, int Z>
inline float& tensorCoeffRef(std::array<Eigen::Matrix<float, Y, Z>, X> tensor, Eigen::Vector3i idx) {
	int x = idx(0), y = idx(1), z = idx(2);
	assert(0 <= x && x < X);
	assert(0 <= y && y < Y);
	assert(0 <= z && z < Z);
	return tensor.at(x).coeffRef(y, z);
}


inline double stablePow(double a, double b) {        
    return static_cast<double> (std::pow(std::cbrt(static_cast<double>(a)),static_cast<double>(b)));
}

inline void inverse33(Eigen::Ref<Eigen::Matrix3d> result, Eigen::Ref<const Eigen::Matrix3d> A) {
    double determinant =    +A(0,0)*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
                            -A(0,1)*(A(1,0)*A(2,2)-A(1,2)*A(2,0))
                            +A(0,2)*(A(1,0)*A(2,1)-A(1,1)*A(2,0));
    double invdet = 1/determinant;
    result(0,0) =  (A(1,1)*A(2,2)-A(2,1)*A(1,2))*invdet;
    result(1,0) = -(A(0,1)*A(2,2)-A(0,2)*A(2,1))*invdet;
    result(2,0) =  (A(0,1)*A(1,2)-A(0,2)*A(1,1))*invdet;
    result(0,1) = -(A(1,0)*A(2,2)-A(1,2)*A(2,0))*invdet;
    result(1,1) =  (A(0,0)*A(2,2)-A(0,2)*A(2,0))*invdet;
    result(2,1) = -(A(0,0)*A(1,2)-A(1,0)*A(0,2))*invdet;
    result(0,2) =  (A(1,0)*A(2,1)-A(2,0)*A(1,1))*invdet;
    result(1,2) = -(A(0,0)*A(2,1)-A(2,0)*A(0,1))*invdet;
    result(2,2) =  (A(0,0)*A(1,1)-A(1,0)*A(0,1))*invdet;
    result.transposeInPlace();
}



#endif 
