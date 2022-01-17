#pragma once

#include <vector>
#include <memory>

#include "param/cloth_param.h"

struct Seam {
    int patch1_id;
    int patch2_id;
    std::vector<std::pair<int, int>> corres;
    // where corres[i].first is in patch1_id
    //       corres[i].second is in patch2_id
    
    // (ordering of pairs not necessary)
    // (it is possible that patch1_id == patch2_id)
};

// note: this doesn't take care of packing! All patches are
// more or less centered around (0,0)
bool finalParamMultiPatch(const std::vector<Eigen::MatrixXd>& vec_V_3d, 
                          const std::vector<Eigen::MatrixXi>& vec_F,
                          const std::vector<std::vector<std::vector<std::pair<int, int>>>>& vec_dart_duplicates,
                          const std::vector<std::vector<int>>& vec_dart_tips,
                          const std::vector<Seam>& seams,
                          std::vector<Eigen::MatrixXd>& vec_V_2d,
                          CLOTH_INIT_TYPE init_type = CLOTH_INIT_ARAP){
    
    int n_patches = vec_V_3d.size();
    if (n_patches != vec_F.size()){
        std::cout << "ERROR: non matching sizes in final param" << std::endl;
    }

    std::vector<std::unique_ptr<ClothParam>> cloth_ps;
    for (int patch_id=0; patch_id<n_patches; patch_id++){
        std::unique_ptr<ClothParam> ptr(new ClothParam(vec_V_3d[patch_id], 
                           vec_F[patch_id], 0.0, 
                           vec_dart_duplicates[patch_id], 
                           vec_dart_tips[patch_id],
                           init_type));
        cloth_ps.push_back(std::move(ptr));
    }
    
    int max_iter = 20;
    for (int current_iter = 0; current_iter < max_iter; current_iter++){

        // TODO update seam info here ("local" step)

        for (int patch_id=0; patch_id<n_patches; patch_id++){

            cloth_ps[patch_id]->paramIter(1);

            #ifdef DEBUG_CLOTH_PARAM
            if (stretch_u_.maxCoeff() < -0.5){ // not really supposed to happen if the initialization is ok
                igl::writeOBJ("../data/buggy/final_not_good.obj", V_3d_, F_);
                igl::writeOBJ("../data/buggy/final_not_good_uv.obj", V_2d_, F_);
            }
            if (std::isnan(stretch_u_.maxCoeff())){ // not really supposed to happen if the initialization is ok
                igl::writeOBJ("../data/buggy/final_nanned.obj", V_3d_, F_);
                igl::writeOBJ("../data/buggy/final_nanned_uv.obj", V_2d_, F_);
            }
            #endif
        }
    }

    bool all_satisfied = true;
    vec_V_2d.resize(n_patches);
    for (int patch_id=0; patch_id<n_patches; patch_id++){
        vec_V_2d[patch_id] = cloth_ps[patch_id]->getV2d();
        if (cloth_ps[patch_id]->checkSelfIntersect() 
        || !cloth_ps[patch_id]->constraintSatisfied()){
            all_satisfied = false;
        }
    }


    std::cout << "TODO finalParamMultiPatch" << std::endl;
    return all_satisfied;
}