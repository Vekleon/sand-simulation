#include <find_air_pressures.h>
#include <vector>

void find_air_pressures(std::vector<Eigen::Vector3i>& pressureIndices, Eigen::Vector3d& P0, Eigen::VectorXd& q,
                        const double dg) {
    int particles = q.size() / 3;
    double radius = dg / 2;
    Eigen::Vector3d grid_offset;
    grid_offset << 0.5, 0.5, 0.5;

    for(int qi = 0; qi < particles; qi++) {
        // NOTE: there are potential duplicates
        Eigen::Vector3d particle_pos = q.segment<3>(qi * 3);
        
        // finding the proper pressure index through particle positioning
        Eigen::Vector3i cell_idx = ((1. / dg) * (particle_pos - P0 + grid_offset));
        cell_idx = cell_idx.cwiseAbs();

        Eigen::Vector3d pressure_pos = P0 + dg * cell_idx;
        double distance = (pressure_pos - particle_pos).norm();

        if(distance <= radius)
            pressureIndices.push_back(cell_idx);
    }
}