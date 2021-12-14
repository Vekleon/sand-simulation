#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {
    double cLength = (q1 - q0).norm();
    double fTerm = (cLength - l0) * stiffness;
    Eigen::Vector3d normalization = (q1 - q0).normalized();

    // need to figure out if I have to reverse these values
    f.segment<3>(0) = -fTerm * normalization;
    f.segment<3>(3) = fTerm * normalization;
}