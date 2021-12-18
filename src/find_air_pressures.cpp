#include <find_air_pressures.h>
#include <vector>
#include <cmath>

void find_air_pressure(std::vector<Eigen::Vector3i>& pressureIndices, Eigen::Vector3d& P0, Eigen::VectorXd& q,
                        const double dg) {
    int particles = q.size() / 3;
    double radius = dg / 2;
    for(int qi = 0; qi < particles; qi++) {
        // NOTE: there are potential duplicates
        Eigen::Vector3d particle_pos = q.segment<3>(qi * 3);
        Eigen::Vector3i cell_idx = ((1. / dg) * (P0 - particle_pos)).cast<int>();
        cell_idx = cell_idx.cwiseAbs();
        pressureIndices.push_back(cell_idx);
    }
}