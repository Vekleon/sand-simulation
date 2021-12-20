#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
#include <set>

void find_connected_rigid(std::vector<std::set<int>>& connections, Eigen::TensorRF& rigid_mapping);