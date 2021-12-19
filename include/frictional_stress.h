#include <Eigen/Dense>
#include <EigenTypes.h>
#include <newtons_method.h>

/**
    Calculates the frictional stress matrix for the particles 

    OUTPUT:
    stress: the stress matrix to be calculated
    mean_stress: the mean stress of the stress matrix
    shear_stress: the shear stress of the system

    INPUT:
    q: the generalized coordinates of the particles in the system
    qdot: the velocity vectors of the particles
    p0: the corner coordinate of the pressure grid
    P: a tensure holding values for the particle grid
    xv: (X+1) * Y * Z tensor of x-component velocities
    yv: X * (Y+1) * Z tensor of y-component velocities
    zv: X * Y * (Z+1) tensor of z-component velocities
    dg: the length of each grid cell.
 */

void frictional_stress(Eigen::MatrixXd& stress, double& mean_stress, double& shear_stress, 
                        Eigen::VectorXd q, Eigen::VectorXd qdot,
                        Eigen::Vector3d p0, Eigen::TensorP& P, Eigen::TensorXV& xv, 
                        Eigen::TensorYV& yv, Eigen::TensorZV& zv, const double dg);