#include "param/cloth_param.h"
#include <igl/avg_edge_length.h>
#include "param/self_intersect.h"

//#define DEBUG_CLOTH_PARAM
#ifdef DEBUG_CLOTH_PARAM
#include <igl/writeOBJ.h>
#endif

#define CHECK_TOPOLOGY_PARAM

ClothParam::ClothParam(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F,
                       double max_stretch,
                       const std::vector<std::vector<std::pair<int, int>>>& dart_duplicates,
                       const std::vector<int>& dart_tips,
                       int seam_size,
                       CLOTH_INIT_TYPE init_type)
            : F_(F), V_3d_(V_3d), max_stretch_(max_stretch), seam_size_(seam_size), init_type_(init_type){
    
    igl::boundary_loop(F, bnd_);

    #ifdef CHECK_TOPOLOGY_PARAM
    std::vector<std::vector<int>> bnds;
    igl::boundary_loop(F, bnds);
    if (bnds.size() != 1) 
        std::cout << "ERROR: wrong topology, number of boundaries = " << bnds.size() << " (should be 1)" << std::endl;
    #endif

    if (init_type_ == CLOTH_INIT_ARAP){
        V_2d_ = paramARAP(V_3d_, F_, bnd_);
    }
    else if (init_type_ == CLOTH_INIT_LSCM){
        V_2d_ = paramLSCM(V_3d_, F_, bnd_);
    }
    else if (init_type_ == CLOTH_INIT_SCAF){
        V_2d_ = paramSCAF(V_3d_, F_, bnd_);
    }

    bool rescale_init = true;
    if (rescale_init){
        //V_2d_ *= (V_3d_.col(0).maxCoeff() - V_3d_.col(0).minCoeff()) / (V_2d_.col(0).maxCoeff() - V_2d_.col(0).minCoeff()); // doesn't work because axes are not the same!
        V_2d_ *= igl::avg_edge_length(V_3d_, F_) / igl::avg_edge_length(V_2d_, F_);
    }

    /*if (checkSelfIntersect()){
        //std::cout << "Self intersecting at init" << std::endl;
        // do something ?
    }*/
    setDartPairs(dart_duplicates, dart_tips);
    bo_.setSeamSize(seam_size_);
    bo_.allocateMemory(F.rows(), V_3d.rows());
    //bo_.measureScore(V_2d_, V_3d_, F_, stretch_u_, stretch_v_);
    
    std::vector<int> selec = autoSelect(V_3d, bnd_);
    bo_.setSelectedVertices(selec);

    Eigen::RowVector3d from = V_2d_.row(selec[0]) - V_2d_.row(selec[1]);
    Eigen::RowVector3d to(1.0, 0, 0);  
    Eigen::Matrix3d R = computeRotation(from, to);
    V_2d_ = (R * V_2d_.transpose()).transpose();
};

void ClothParam::paramIter(int n_iter){
    for (int i=0; i<n_iter; i++){
        V_2d_ = bo_.localGlobal(V_2d_, V_3d_, F_);
    }
}

bool ClothParam::paramAttempt(int max_iter){
    if (!topology_ok) {
        std::cout << "Topology check failed, aborting parameterization" << std::endl;
        V_2d_ = Eigen::MatrixXd::Zero(V_3d_.rows(), 3);
        stretch_u_ = Eigen::VectorXd::Constant(F_.rows(), 100.0);
        stretch_v_ = Eigen::VectorXd::Constant(F_.rows(), 100.0);
        return false;
    }
    for (int current_iter = 0; current_iter < max_iter; current_iter++){
        bo_.measureScore(V_2d_, V_3d_, F_, stretch_u_, stretch_v_);

        #ifdef DEBUG_CLOTH_PARAM
        if (stretch_u_.maxCoeff() < -0.5){ // not really supposed to happen if the initialization is ok
            igl::writeOBJ("../data/buggy/not_good.obj", V_3d_, F_);
            igl::writeOBJ("../data/buggy/not_good_uv.obj", V_2d_, F_);
        }
        if (std::isnan(stretch_u_.maxCoeff())){ // not really supposed to happen if the initialization is ok
            igl::writeOBJ("../data/buggy/nanned.obj", V_3d_, F_);
            igl::writeOBJ("../data/buggy/nanned_uv.obj", V_2d_, F_);
        }
        #endif

        /*if (constraintSatisfied()){
            return true;
        }*/
        V_2d_ = bo_.localGlobal(V_2d_, V_3d_, F_);
        
        if (enable_intersection_check_){
            if (checkSelfIntersect()){
                // interrupt process early
                // note: we do this after one iteration so we don't punish bad initialization
                return false;
            }
        }
    }
    return constraintSatisfied();
}

void ClothParam::printStretchStats() const {
    if (stretch_u_.rows() < 1) return;
    std::cout << "Stretch min -> max (avg): " << std::endl;
    printf("U: %f -> %f (%f)\n", stretch_u_.minCoeff(), stretch_u_.maxCoeff(), stretch_u_.mean());
    printf("V: %f -> %f (%f)\n", stretch_v_.minCoeff(), stretch_v_.maxCoeff(), stretch_v_.mean());
};


void ClothParam::getStretchStats(Eigen::VectorXd& stretch_u, Eigen::VectorXd& stretch_v) const {
    stretch_u = stretch_u_;
    stretch_v = stretch_v_;
}


void ClothParam::measureStretchStats(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d, 
                                const Eigen::MatrixXi& F, Eigen::VectorXd& stretch_u, 
                                Eigen::VectorXd& stretch_v){
    BaryOptimizer bo;
    bo.measureScore(V_2d, V_3d, F, stretch_u, stretch_v);
}

void ClothParam::setAlignmentVertexPair(int v1_id, int v2_id){
    bo_.setSelectedVertices({v1_id, v2_id});
}

void ClothParam::setDartPairs(const std::vector<std::vector<std::pair<int, int>>>& dart_duplicates,
                    const std::vector<int>& dart_tips){
    bo_.setUnorderedDarts(dart_duplicates, dart_tips);
}

bool ClothParam::checkSelfIntersect() const {
    return selfIntersect(V_2d_, bnd_);
}