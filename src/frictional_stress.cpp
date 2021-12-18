#include <frictional_stress.h>
#include <cmath>

void frictional_stress(Eigen::MatrixXd frictions, Eigen::VectorXd& q, Eigen::VectorXd& qdot,
                        Eigen::Vector3d& p0, Eigen::TensorP& P, Eigen::TensorXV& xv, 
                        Eigen::TensorYV& yv, Eigen::TensorZV& zv, const double dg) {
    
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
}