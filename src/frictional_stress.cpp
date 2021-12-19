#include <frictional_stress.h>
#include <cmath>

void frictional_stress(Eigen::MatrixXd& stress, double& mean_stress, double& shear_stress, 
                        Eigen::VectorXd q, Eigen::VectorXd qdot,
                        Eigen::Vector3d p0, Eigen::TensorP& P, Eigen::TensorXV& xv, 
                        Eigen::TensorYV& yv, Eigen::TensorZV& zv, const double dg) {
    stress.resize(q.size() / 3, q.size() / 3);
    Eigen::MatrixXd GRID_OFFSETS(4, 3);
	GRID_OFFSETS <<
		0.5, 0.5, 0.0,
		0.0, 0.5, 0.5,
		0.5, 0.0, 0.5,
        0.5, 0.5, 0.5;
    
    // fuck...
    Eigen::MatrixXd D;

    int particles = q.size() / 3;
    for(int pi = 0; pi < particles; pi++){
        Eigen::Vector3d particle_pos = q.segment<3>(pi * 3);
        Eigen::Vector3i cell_index = ((1. / dg) * (particle_pos - p0 + GRID_OFFSETS.row(3).transpose()));
        
    }

    mean_stress = stress.trace() / 3;
    Eigen::MatrixXd delta;
    delta.resize(q.size() / 3, q.size() / 3);
    for(int row = 0; row < delta.rows(); row++){
        for(int col = 0; col < delta.cols(); col++){
            delta(row, col) = 0.1 * mean_stress;
        }
    }
    shear_stress = (stress - delta).norm() / (std::sqrt(2));
}