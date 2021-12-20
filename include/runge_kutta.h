//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Runge-Kutta time integration
//  qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration

template<typename FORCE> 
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {
    Eigen::VectorXd qCopy = q;
    Eigen::VectorXd qdotCopy = qdot;
    
    Eigen::VectorXd k1Force, k2Force, k3Force, k4Force;
    Eigen::VectorXd totalV, totalA, v2, v3;

    force(k1Force, q, qdot);
    totalV = qdot;
    totalA = k1Force / mass;

    //calculate k2
    forward_euler(qCopy, qdotCopy, dt * 0.5, mass, force);
    force(k2Force, qCopy, qdotCopy);
    v2 = qdotCopy;
    totalV += qdotCopy;
    totalA += k2Force / mass;

    q += dt / 2 * totalV;
    qdot += dt / 2 * totalA;
}