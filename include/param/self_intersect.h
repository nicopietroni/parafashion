#pragma once

#include <Eigen/Geometry>

// Simplified case where we assume non colinearity
// https://stackoverflow.com/questions/3838329/how-can-i-check-if-two-segments-intersect
bool ccw(const Eigen::RowVector2d& A,
         const Eigen::RowVector2d& B,
         const Eigen::RowVector2d& C){

    return (C(1)-A(1)) * (B(0)-A(0)) > (B(1)-A(1)) * (C(0)-A(0));
}

bool segmentIntersect(const Eigen::RowVector2d& A,
                      const Eigen::RowVector2d& B,
                      const Eigen::RowVector2d& C,
                      const Eigen::RowVector2d& D){
    
    return ccw(A, C, D) != ccw(B, C, D) && ccw(A, B, C) != ccw(A, B, D);                  
}


bool selfIntersect(const Eigen::MatrixXd& V_2d, const Eigen::VectorXi& bnd){
    double t,u;
    for (int i=0; i<bnd.rows(); i++){
        Eigen::RowVector2d p1 = V_2d.row(bnd(i));
        int next1 = (i + 1) %  bnd.rows();
        Eigen::RowVector2d p2 = V_2d.row(bnd(next1));
        for (int j=i+2; j<bnd.rows(); j++){
            if ((j+1) % bnd.rows() == i) continue; // don't test consecutive edges
            Eigen::RowVector2d q1 = V_2d.row(bnd(j));
            int next2 = (j + 1) %  bnd.rows();
            Eigen::RowVector2d q2 = V_2d.row(bnd(next2));
            if (segmentIntersect(p1, p2, q1, q2)){
                //std::cout << "intersect: " << i << " --- " << j << " /" << bnd.rows() << std::endl;
                return true;
            }
        }    
    }

    return false;
}