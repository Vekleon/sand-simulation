#include <Eigen/Dense>
#include <EigenTypes.h>
#include <newtons_method.h>

/**
    Calculates the frictional stress matrix for the particles 

    OUTPUT:
    stress: the stress tensor used to store calculated stress for each grid cell
    rigid_map: a tensor used to tell us whether a cell is rigid or not

    INPUT:
    P: a tensure holding values for the particle grid
    xv: (X+1) * Y * Z tensor of x-component velocities
    yv: X * (Y+1) * Z tensor of y-component velocities
    zv: X * Y * (Z+1) tensor of z-component velocities
    dg: the length of each grid cell.
    density: the density of our sand
    dt: the timestep size
    cohesion: the cohesion coefficient.
 */

void frictional_stress(Eigen::TensorS& stress, Eigen::TensorRF rigid_map, 
                        Eigen::TensorP& P, Eigen::TensorXV& xv, Eigen::TensorYV& yv, 
                        Eigen::TensorZV& zv, const double dg, const double density, 
                        double dt, double friction_angle, const double cohesion);
