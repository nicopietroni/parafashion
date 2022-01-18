#include "param/param_utils.h"

#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <Eigen/LU>
#include <iostream>
#include <cmath>

#include <igl/lscm.h>

void procustes(const Eigen::MatrixXd& points1, // to
               const Eigen::MatrixXd& points2, // from
               Eigen::MatrixXd& R_est,
               Eigen::VectorXd& T_est){

    // https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
    // https://math.stackexchange.com/questions/849217/estimate-rotation-and-translation-from-two-sets-of-points-in-different-coordinat
    // https://en.wikipedia.org/wiki/Procrustes_analysis

    Eigen::MatrixXd points1t = points1.transpose();
    Eigen::MatrixXd points2t = points2.transpose();

    Eigen::VectorXd pb = points1t.rowwise().mean();
    Eigen::VectorXd qb = points2t.rowwise().mean();

    Eigen::MatrixXd X = (points1t.colwise() - pb);
    Eigen::MatrixXd Y = (points2t.colwise() - qb);
    Eigen::MatrixXd S = X * Y.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd sigma = Eigen::MatrixXd::Identity(svd.matrixU().cols(), svd.matrixV().cols());
    sigma(sigma.rows() - 1, sigma.cols() - 1) = -(svd.matrixV() * svd.matrixU().transpose()).determinant();
    R_est = svd.matrixV() * sigma * svd.matrixU().transpose();
    T_est = qb - R_est * pb;
}

// Move triangle from 3D (Z != 0) to 2D (Z = 0), with
// same lengths but arbitrary orientation
Eigen::MatrixXd move3Dto2D(const Eigen::MatrixXd& V_tri){
    Eigen::MatrixXd V_2d(3,3); // TODO Matrix3d
    // First, move 3D triangle to 2D plane
    // V_2d: put A in (0,0), B in (0, |AB|), and find C 
    double r0 = (V_tri.row(1) - V_tri.row(0)).norm();
    double r1 = (V_tri.row(2) - V_tri.row(0)).norm();
    double r2 = (V_tri.row(2) - V_tri.row(1)).norm();
    V_2d.row(0) = Eigen::RowVector3d(0, 0, 0);
    V_2d.row(1) = Eigen::RowVector3d(r0, 0, 0);
    double CAB_angle = std::acos((r0*r0 + r1*r1 - r2*r2)/(2*r0*r1));
    double l1 = r1 * std::cos(CAB_angle);
    double h = l1 * std::tan(CAB_angle);
    V_2d.row(2) = Eigen::RowVector3d(l1, h, 0);

    r0 = (V_2d.row(1) - V_2d.row(0)).norm();
    r1 = (V_2d.row(2) - V_2d.row(0)).norm();
    r2 = (V_2d.row(2) - V_2d.row(1)).norm();
    if ( std::fabs(r0 - (V_tri.row(1) - V_tri.row(0)).norm()) // just checking...
        +std::fabs(r1 - (V_tri.row(2) - V_tri.row(0)).norm())
        +std::fabs(r2 - (V_tri.row(2) - V_tri.row(1)).norm()) > 0.0001){
        std::cout << "ERROR, flat triangle is different:" << std::endl;
        std::cout << r0 << " vs " << (V_tri.row(1) - V_tri.row(0)).norm() << std::endl;
    }

    return V_2d;
}

// TODO Matrix3d ?
Eigen::MatrixXd makeTriPoints(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int f_id){
    Eigen::MatrixXd p(3,3);
    p.row(0) = V.row(F(f_id,0));
    p.row(1) = V.row(F(f_id,1));
    p.row(2) = V.row(F(f_id,2));
    return p;
};

void makeTriPoints(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int f_id, Eigen::MatrixXd& out){
    out.row(0) = V.row(F(f_id,0));
    out.row(1) = V.row(F(f_id,1));
    out.row(2) = V.row(F(f_id,2));
}

Eigen::VectorXd vertices2dToVector(const Eigen::MatrixXd& V){
    Eigen::VectorXd res(2 * V.rows()); // TODO faster?
    for (int i=0; i<V.rows(); i++){
        res(2 * i) = V(i,0);
        res(2 * i + 1) = V(i,1);
    }
    return res;
}

Eigen::Vector3d barycentricCoords(const Eigen::RowVector3d& p, const Eigen::RowVector3d& a, 
                                         const Eigen::RowVector3d& b, const Eigen::RowVector3d& c){
    Eigen::RowVector3d v0 = b - a;
    Eigen::RowVector3d v1 = c - a;
    Eigen::RowVector3d v2 = p - a;
    float d00 = v0.dot(v0);
    float d01 = v0.dot(v1);
    float d11 = v1.dot(v1);
    float d20 = v2.dot(v0);
    float d21 = v2.dot(v1);
    float denom = d00 * d11 - d01 * d01;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0f - v - w;
    return Eigen::Vector3d(u, v, w);
}

Eigen::Matrix3d computeRotation(const Eigen::RowVector3d& from,
                                const Eigen::RowVector3d& to){
    // There might already be something like this in Eigen? Couldn't find it
    // https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    Eigen::RowVector3d a = from.normalized();
    Eigen::RowVector3d b = to.normalized();
    if (a==b || a == -b) std::cout << "ERROR: case not handled in computeRotation" << std::endl;
    Eigen::RowVector3d v = a.cross(b);
    double s = v.norm();
    double c = a.dot(b);
    Eigen::Matrix3d vs;
    vs <<     0, -v[2],  v[1], 
           v[2],     0, -v[0], 
          -v[1],  v[0],     0;

    return Eigen::Matrix3d::Identity() + vs + vs * vs * 1.0 / (1.0+c);
}

#include <igl/arap.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>

Eigen::MatrixXd paramARAP(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F,
                          const Eigen::VectorXi& bnd){

    Eigen::MatrixXd V_2d, V_3db;
    V_3db = V_3d;
    double scale = V_3d.maxCoeff();
    V_3db = V_3db.array() / scale; 

    // Compute the initial solution for ARAP (harmonic parametrization)
    Eigen::MatrixXd bnd_uv;
    igl::map_vertices_to_circle(V_3db, bnd, bnd_uv);

    Eigen::MatrixXd initial_guess;
    igl::harmonic(V_3db, F, bnd, bnd_uv, 1, initial_guess);

    // Add dynamic regularization to avoid to specify boundary conditions
    igl::ARAPData arap_data;
    arap_data.with_dynamics = true;
    Eigen::VectorXi b  = Eigen::VectorXi::Zero(0);
    Eigen::MatrixXd bc = Eigen::MatrixXd::Zero(0,0);

    // Initialize ARAP
    arap_data.max_iter = 100;
    // 2 means that we're going to *solve* in 2d
    arap_precomputation(V_3db, F, 2, b, arap_data);

    // Solve arap using the harmonic map as initial guess
    V_2d = initial_guess;
    arap_solve(bc, arap_data, V_2d);


    Eigen::MatrixXd V_2db = Eigen::MatrixXd::Zero(V_2d.rows(), 3);
    V_2db.col(0) = V_2d.col(0);
    V_2db.col(1) = V_2d.col(1);
    
    return V_2db.array() * scale;
}


Eigen::MatrixXd paramLSCM(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F,
                          const Eigen::VectorXi& bnd){
    Eigen::VectorXi b(2,1);
    b(0) = bnd(0);
    b(1) = bnd(bnd.size()/2);
    Eigen::MatrixXd bc(2, 2);
    bc<<0,0,1,0;

    Eigen::MatrixXd V_2d;
    igl::lscm(V_3d,F,b,bc,V_2d);
    Eigen::MatrixXd V_2db = Eigen::MatrixXd::Zero(V_2d.rows(), 3);
    V_2db.col(0) = V_2d.col(0);
    V_2db.col(1) = V_2d.col(1);
    return V_2db;
}

// SCAF
/*#include <igl/triangle/scaf.h>
#include <igl/arap.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/MappingEnergyType.h>
#include <igl/doublearea.h>
#include <igl/PI.h>
#include <igl/flipped_triangles.h>
#include <igl/topological_hole_fill.h>*/

Eigen::MatrixXd paramSCAF(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F, const Eigen::VectorXi& bnd){

    /*int scaf_iterations = 10;

    Eigen::MatrixXd bnd_uv, uv_init;
    igl::triangle::SCAFData scaf_data;

    Eigen::VectorXd M;
    igl::doublearea(V_3d, F, M);*/
    
    /*std::vector<std::vector<int>> all_bnds;
    igl::boundary_loop(F, all_bnds);
    // Heuristic primary boundary choice: longest
    auto primary_bnd = std::max_element(all_bnds.begin(), all_bnds.end(), [](const std::vector<int> &a, const std::vector<int> &b) { return a.size()<b.size(); });

    Eigen::VectorXi bnd = Eigen::Map<Eigen::VectorXi>(primary_bnd->data(), primary_bnd->size());*/
    /*
    igl::map_vertices_to_circle(V_3d, bnd, bnd_uv);
    bnd_uv *= sqrt(M.sum() / (2 * igl::PI));

    if (bnd.rows() == V_3d.rows()){ // case: all vertex on boundary
        uv_init.resize(V_3d.rows(), 2);
        for (int i = 0; i < bnd.rows(); i++)
        uv_init.row(bnd(i)) = bnd_uv.row(i);
    }
    else{
        igl::harmonic(V_3d, F, bnd, bnd_uv, 1, uv_init);
        if (igl::flipped_triangles(uv_init, F).size() != 0)
        igl::harmonic(F, bnd, bnd_uv, 1, uv_init); // fallback uniform laplacian
    }
    
    

    Eigen::VectorXi b; Eigen::MatrixXd bc;

    igl::triangle::scaf_precompute(V_3d, F, uv_init, scaf_data, igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b, bc, 0);

    for (int i=0; i<scaf_iterations; i++){
        igl::triangle::scaf_solve(scaf_data, 1);
    }

    Eigen::MatrixXd V_2d(V_3d.rows(), 2);
    V_2d = scaf_data.w_uv.topRows(V_3d.rows());


    Eigen::MatrixXd V_2db = Eigen::MatrixXd::Zero(V_3d.rows(), 3);
    V_2db.col(0) = V_2d.col(0);
    V_2db.col(1) = V_2d.col(1);
    return V_2d;*/
    std::cout << "ERROR SCAF DISABLED" << std::endl;
    return V_3d;
}


