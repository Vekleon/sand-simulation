#include <iostream>
#include <thread>

#include <visualization.h>
#include <igl/edges.h>
#include <igl/edge_lengths.h>
#include <igl/readMESH.h>
#include <igl/boundary_facets.h>
#include <Eigen/Dense>

#include <init_state.h>
#include <find_min_vertices.h>
#include <fixed_point_constraints.h>
#include <linearly_implicit_euler.h>
#include <V_spring_particle_particle.h>

//Simulation State
bool simulating = true;
Eigen::VectorXd q;
Eigen::VectorXd qdot;
Eigen::SparseMatrixd M;
Eigen::SparseMatrixd P; //fixed point constraints 
Eigen::VectorXd x0; //fixed point constraints 

double t = 0; //simulation time 
double dt = 0.005; //time step
double k  = 1e5;
double k_selected = 1e5; //stiff spring for pulling on object
double m = 1.;

Eigen::MatrixXd V; //vertices of simulation mesh 
Eigen::MatrixXi T; //faces of simulation mesh
Eigen::MatrixXi F; //faces of simulation mesh
Eigen::MatrixXi E; //edges of simulation mesh (which will become springs)
Eigen::VectorXd l0; //original length of all edges in the mesh

//working memory for integrator
Eigen::VectorXd tmp_force;
Eigen::SparseMatrixd tmp_stiffness;

std::vector<unsigned int> fixed_point_indices;

void simulate() {

    while(simulating) {
        //TODO
    }
}

bool draw(igl::opengl::glfw::Viewer & viewer) {

    //    double V_spring, V_gravity, T_p, KE,PE;
    //    KE = 0;
    //    PE = 0;

    //    for(unsigned int p = 0; p < V.rows(); p++) {
    //        T_particle(T_p, (P.transpose()*qdot).segment<3>(3*p), m);
    //        V_gravity_particle(V_gravity, (P.transpose()*q+x0).segment<3>(3*p), m, Eigen::Vector3d(0., -9.8, 0.));

    //        PE += V_gravity;
    //        KE += T_p;
    //    }
    //    
    //    for(unsigned int ei = 0; ei < E.rows(); ei++) {
    //        V_spring_particle_particle(V_spring, (P.transpose()*q+x0).segment<3>(3*E(ei,0)), (P.transpose()*q+x0).segment<3>(3*E(ei,1)), l0(ei), k);

    //        PE += V_spring;
    //    }
    //    
    //    Visualize::add_energy(t, KE, PE);

    ////update vertex positions using simulation
    //Visualize::update_vertex_positions(0, P.transpose()*q + x0);

    return false;
}

int main(int argc, char **argv) {
    
    std::cout<<"sand-simulation\n";

    //load geometric data 
    igl::readMESH("../../data/coarse_bunny.mesh",V,T, F);
    igl::boundary_facets(T, F);
    F = F.rowwise().reverse().eval();
    igl::edges(T,E);
    igl::edge_lengths(V,E,l0);

    //setup simulation 
    init_state(q,qdot,V);
    
    //setup constraint matrix
    // find_min_vertices(fixed_point_indices, V, 3);
    // P.resize(q.rows(),q.rows());
    // P.setIdentity();
    // fixed_point_constraints(P, q.rows(), fixed_point_indices);

    x0 = q - P.transpose()*P*q; //vector x0 contains position of all fixed nodes, zero for everything else
    
    //correct M, q and qdot so they are the right size
    q = P*q;
    qdot = P*qdot;
    M = P*M*P.transpose();

    //run simulation in seperate thread to avoid slowing down the UI
    std::thread simulation_thread(simulate);
    simulation_thread.detach();

    std::function<Eigen::Ref<const Eigen::VectorXd> ()> test = [&]()->Eigen::Ref<const Eigen::VectorXd> { return P.transpose()*q+x0; };

    //setup libigl viewer and activate 
    Visualize::setup(q, qdot, true);
    //Visualize::add_object_to_scene(V,F, Eigen::RowVector3d(244,165,130)/255.);
    Visualize::viewer().callback_post_draw = &draw;
    
    Visualize::viewer().launch_init(true,false,"Mass-Spring Systems",0,0);
    Visualize::viewer().launch_rendering(true);
    simulating = false;
    Visualize::viewer().launch_shut();

}
