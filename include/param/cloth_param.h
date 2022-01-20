#pragma once

#include <Eigen/Core>
#include <igl/boundary_loop.h>

#include "param/bary_optimizer.h"
#include "param/param_utils.h"
#include "param/auto_select.h"

enum CLOTH_INIT_TYPE {CLOTH_INIT_ARAP, CLOTH_INIT_LSCM, CLOTH_INIT_SCAF};

class ClothParam {
public:

    // ----- MAIN FUNCTIONS ----- //

    // Initialize and allocate memory
    ClothParam(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F,
               double max_stretch = 0.05,
               const std::vector<std::vector<std::pair<int, int>>>& dart_duplicates = {},
               const std::vector<int>& dart_tips = {},
               int seam_size = 0,
               CLOTH_INIT_TYPE = CLOTH_INIT_LSCM);

    // Tries to parameterize. Returns whether max_stretch 
    // was satisfied within max number of iterations,
    // AND solution doesn't self intersect.
    bool paramAttempt(int max_iter = 10);

    // ----- OPTIONAL ----- //

    // Tries to align V axis in parameterization with
    // direction v1 -> v2 
    void setAlignmentVertexPair(int v1_id, int v2_id);

    bool checkSelfIntersect() const;

    void setSeamTargets(const std::vector<Eigen::MatrixXd>& targets_p,
                        const std::vector<Eigen::VectorXi>& p_ids){
        bo_.setSeamTargets(targets_p, p_ids);
    }

    // ----- MISC. ----- //

    // Runs param for specified number of iterations, without intersection check or 
    void paramIter(int n_iter);

    void printStretchStats() const;
    void getStretchStats(Eigen::VectorXd& stretch_u, Eigen::VectorXd& stretch_v) const;

    static void measureStretchStats(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d, 
                                    const Eigen::MatrixXi& F, Eigen::VectorXd& stretch_u, 
                                    Eigen::VectorXd& stretch_v);

    Eigen::MatrixXd getV2d(){return V_2d_;}

    bool constraintSatisfied() const {
        if (stretch_u_.rows() < 1) return false;
        return stretch_u_.minCoeff() > -max_stretch_ 
            && stretch_u_.maxCoeff() < max_stretch_ 
            && stretch_v_.minCoeff() > -max_stretch_ 
            && stretch_v_.maxCoeff() < max_stretch_;
    };

    bool topology_ok = true;

    void disableIntersectionCheck(){enable_intersection_check_ = false;};
    void enableIntersectionCheck(){enable_intersection_check_ = true;};

    void loadConfig(std::string config_path);

private:
    // Set during init
    const Eigen::MatrixXd V_3d_;
    const Eigen::MatrixXi F_;
    const double max_stretch_;
    int seam_size_;
    CLOTH_INIT_TYPE init_type_;

    Eigen::VectorXi bnd_;
    Eigen::MatrixXd V_2d_;
    BaryOptimizer bo_;

    Eigen::VectorXd stretch_u_, stretch_v_;

    // config
    bool enable_intersection_check_ = true;

    // Sets pairs of vertex ids which should be symmetrical in a given dart.
    // first indexing: dart id
    // second indexing: pair id in dart    
    // pair: two vertex ids, which should be symmetrical w.r.t. dart
    // note: needs to be called before memory allocation!
    void setDartPairs(const std::vector<std::vector<std::pair<int, int>>>& dart_duplicates,
                      const std::vector<int>& dart_tips);
};