#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/*
    Finds the pressure grid indices which not in the liquid

    OUTPUT:
    pressureIndices: a list of indices into our pressure grid that are not submerged

    P0: the position of the first pressure point in our grid
    q: the generalized coordinates
    dg: the distance between points on the grid
*/
void find_air_pressure(Eigen::MatrixXi& pressureIndices, Eigen::Vector3d& P0, Eigen::VectorXd& q,
                        const double dg);