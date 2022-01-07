#pragma once

#include <Eigen/Core>

#include "param/bary_optimizer.h"
#include "param/param_utils.h"

class ClothParam {
public:
    ClothParam(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F,
               double max_stretch = 0.05)
               : F_(F), V_3d_(V_3d), max_stretch_(max_stretch){
        V_2d_ = paramARAP(V_3d_, F_);
        bo_.allocateMemory(F.rows(), V_3d.rows());
        bo_.measureScore(V_2d_, V_3d_, F_, stretch_u_, stretch_v_);
    };

    bool paramAttempt(int max_iter = 10){
        for (int current_iter = 0; current_iter < max_iter; current_iter++){
            bo_.measureScore(V_2d_, V_3d_, F_, stretch_u_, stretch_v_);
            if (constraintSatisfied()){
                return true;
            }
            V_2d_ = bo_.localGlobal(V_2d_, V_3d_, F_);
        }
        return constraintSatisfied();
    }

    void printStretchStats(){
        if (stretch_u_.rows() < 1) return;
        std::cout << "Stretch min -> max (avg): " << std::endl;
        printf("U: %f -> %f (%f)\n", stretch_u_.minCoeff(), stretch_u_.maxCoeff(), stretch_u_.mean());
        printf("V: %f -> %f (%f)\n", stretch_v_.minCoeff(), stretch_v_.maxCoeff(), stretch_v_.mean());
    };

    Eigen::MatrixXd getV2d(){return V_2d_;}

private:
    // Set during init
    const Eigen::MatrixXd V_3d_;
    const Eigen::MatrixXi F_;
    const double max_stretch_;

    Eigen::MatrixXd V_2d_;
    BaryOptimizer bo_;

    Eigen::VectorXd stretch_u_, stretch_v_;

    bool constraintSatisfied(){
        if (stretch_u_.rows() < 1) return false;
        return stretch_u_.minCoeff() > -max_stretch_ 
            && stretch_u_.maxCoeff() < max_stretch_ 
            && stretch_v_.minCoeff() > -max_stretch_ 
            && stretch_v_.maxCoeff() < max_stretch_;
    };
};