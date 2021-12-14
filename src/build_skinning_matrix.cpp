#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>
#include <iostream>

bool inside(Eigen::Vector4d bary) {
    for (int i = 0; i < 4; i++)
        if (bary[i] > 1 || bary[i] < 0) return false;
    return true;
}

void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T,
                           Eigen::Ref<const Eigen::MatrixXd> V_skin)
{
    //typedef Eigen::Triplet<double> Trip;
    //std::vector<Trip> triplets;
    //N.resize(V_skin.rows(), V.rows());
    //for (int i = 0; i < V_skin.rows(); i++) {
    //    Eigen::Vector3d X = V_skin.row(i).transpose();
    //    for (int j = 0; j < T.rows(); j++) {
    //        Eigen::Vector4d barycentric;
    //        Eigen::RowVectorXi element = T.row(j);
    //        phi_linear_tetrahedron(barycentric, V, element, X);
    //        if (inside(barycentric)) {
    //            for(int k = 0; k < 4; k++)
    //                triplets.push_back(Trip(i, element[k], barycentric[k]));
    //            break;
    //        }
    //    }
    //}
    //N.setFromTriplets(triplets.begin(), triplets.end());
}