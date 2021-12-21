#include <frictional_stress.h>
#include <math.h>

void frictional_stress(Eigen::TensorS& stress, Eigen::TensorRF rigid_map, 
                        Eigen::TensorP& P, Eigen::TensorXV& xv, Eigen::TensorYV& yv, 
                        Eigen::TensorZV& zv, const double dg, const double density, 
                        double dt, double friction_angle, const double cohesion) {
    // constructs the velocity divergence matrix.
    auto construct_u = [&xv, &yv, &zv](Eigen::Matrix3d& u, int x, int y, int z) {
        double x_base = xv.at(x)(y, z);
        double y_base = yv.at(x)(y, z);
        double z_base = zv.at(x)(y, z);
        u(0, 0) = xv.at(x + 1)(y, z) - x_base;
        u(0, 1) = xv.at(x)(y + 1, z) - x_base;
        u(0, 2) = xv.at(x)(y, z + 1) - x_base;
        u(1, 0) = yv.at(x + 1)(y, z) - y_base;
        u(1, 1) = yv.at(x)(y + 1, z) - y_base;
        u(1, 2) = yv.at(x)(y, z + 1) - y_base;
        u(2, 0) = zv.at(x + 1)(y, z) - z_base;
        u(2, 1) = zv.at(x)(y + 1, z) - z_base;
        u(2, 2) = zv.at(x)(y, z + 1) - z_base;
    };

    Eigen::Matrix3d D, u, rigid_stress;
    int x_dimension = P.size();
    int y_dimension = P.at(0).rows();
    int z_dimension = P.at(0).cols();

    // Calculate the friction value of every grid point.
    for(int z = 0; z < z_dimension; z++){
        for(int y = 0; y < y_dimension; y++){
            for(int x = 0; x < x_dimension; x++){
                construct_u(u, x, y , z);
                u /= dg;
                D = (u + u.transpose()) / 2;
                rigid_stress = density * D * dg * dg / dt;
                double mean = rigid_stress.trace() / 3;
                double shear = (rigid_stress - mean * Eigen::Matrix3d::Identity()).norm() / sqrt(2);
                
                double left = sqrt(3) * shear;
                double right = friction_angle * mean + cohesion;
                if(left < right) {
                    stress[x][y][z] = rigid_stress;
                    rigid_map[x][y][z] = true;
                } else {
                    double denominator = sqrt(1/3) * D.norm();
                    stress[x][y][z] = - friction_angle * P.at(x)(y, z) * D / denominator;
                    rigid_map[x][y][z] = false;
                }
            }
        }
    }
}