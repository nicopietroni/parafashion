
#include <vector>
#include <Eigen/Core>

#include "cloth_param.h"

// Flattens a patch, given multiple poses (several V_3d, but same F) 
Eigen::MatrixXd multiplePosesParam(const std::vector<Eigen::MatrixXd>& vec_V_3d, 
                                   const Eigen::MatrixXi& F,
                                   double stretch,
                                   const std::vector<std::vector<std::pair<int, int>>>& vec_dart_duplicates,
                                   const std::vector<int>& vec_dart_tips){

    ClothParam main_cp(vec_V_3d[0], F);
    std::vector<std::unique_ptr<ClothParam>> other_cps;
    for (int i=1; i<vec_V_3d.size(); i++){
        std::unique_ptr<ClothParam> ptr(new ClothParam(vec_V_3d[i], F, stretch, vec_dart_duplicates, vec_dart_tips));
        other_cps.push_back(std::move(ptr));
    }

    int n_iter = 10;
    for (int iter = 0; iter < n_iter; iter ++){
        main_cp.sendMainV2d(other_cps);
        for (int cp_id = 0; cp_id < other_cps.size(); cp_id ++){
            other_cps[cp_id]->paramAttempt(1); // Computes b (+ more stuff, it's an overkill but it's fine)
        }

        main_cp.stealSubBs(other_cps);
        main_cp.paramAttempt(1);
    }

    return main_cp.getV2d();
}