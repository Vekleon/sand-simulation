#include <find_connected_rigid.h>

bool out_of_bounds(int x, int y, int z) {
    return x >= TENSOR_P_X || y >= TENSOR_P_Y || z >= TENSOR_P_Z;
}

void traverse(std::set<int>& explored, std::set<int>& new_connection, Eigen::TensorRF& rigid_mapping, int x, int y, int z) {
    int target = get_cell_idx(x, y, z, TENSOR_P_X, TENSOR_P_Y, TENSOR_P_Z);
    // check to see if the point is in bounds and hasn't been explored yet or if it's rigid
    if(out_of_bounds(x, y, z) || !rigid_mapping[x][y][z] || explored.find(target) != explored.end()) {
        explored.insert(target);
        return;
    }

    // add point to the connections and traverse all other neighboring points.
    new_connection.insert(target);
    explored.insert(target);
    traverse(explored, new_connection, rigid_mapping, x + 1, y , z);
    traverse(explored, new_connection, rigid_mapping, x - 1, y , z);
    traverse(explored, new_connection, rigid_mapping, x, y + 1, z);
    traverse(explored, new_connection, rigid_mapping, x, y - 1, z);
    traverse(explored, new_connection, rigid_mapping, x, y, z + 1);
    traverse(explored, new_connection, rigid_mapping, x, y, z - 1);
}

void find_connected_rigid(std::vector<std::set<int>>& connections, Eigen::TensorRF& rigid_mapping) {
    std::set<int> explored;
    int target;

    // Search throughout the grid and connect any rigid grid points found next to each other
    for(int z = 0; z < TENSOR_P_Z; z++) {
        for(int y = 0; y < TENSOR_P_Y; y++) {
            for(int x = 0; x < TENSOR_P_X; x++) {
                target = get_cell_idx(x, y, z, TENSOR_P_X, TENSOR_P_Y, TENSOR_P_Z);
                // checks to see if this point has been explor yet or if it's rigid
                if(!rigid_mapping[x][y][z] || explored.find(target) != explored.end()) {
                    explored.insert(target);
                    continue;
                }
                std::set<int> new_connection {target};
                traverse(explored, new_connection, rigid_mapping, x + 1, y , z);
                traverse(explored, new_connection, rigid_mapping, x - 1, y , z);
                traverse(explored, new_connection, rigid_mapping, x, y + 1, z);
                traverse(explored, new_connection, rigid_mapping, x, y - 1, z);
                traverse(explored, new_connection, rigid_mapping, x, y, z + 1);
                traverse(explored, new_connection, rigid_mapping, x, y, z - 1);

                explored.insert(new_connection.begin(), new_connection.end());
                connections.push_back(new_connection);
            }
        }
    }
}