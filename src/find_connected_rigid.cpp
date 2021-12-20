#include <find_connected_rigid.h>

void traverse(std::set<int>& new_connection, Eigen::TensorRF& rigid_mapping, int x, int y, int z) {

}

void find_connected_rigid(std::vector<std::set<int>>& connections, Eigen::TensorRF& rigid_mapping) {
    std::set<int> added;
    int target;

    for(int z = 0; z < TENSOR_P_Z; z++) {
        for(int y = 0; y < TENSOR_P_Y; y++) {
            for(int x = 0; x < TENSOR_P_X; x++) {
                target = get_cell_idx(x, y, z, TENSOR_P_X, TENSOR_P_Y, TENSOR_P_Z);
                if(added.find(target) != added.end() || !rigid_mapping[x][y][z]) continue;
                std::set<int> new_connection {target};
                traverse(new_connection, rigid_mapping, x + 1, y , z);
                traverse(new_connection, rigid_mapping, x, y + 1, z);
                traverse(new_connection, rigid_mapping, x, y, z + 1);

                added.insert(new_connection.begin(), new_connection.end());
            }
        }
    }
}