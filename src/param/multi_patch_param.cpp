
#define USE_JSON_CONFIG

#include "param/multi_patch_param.h"
#include "param/metrics.h"

//#define MULTI_PARAM_DEBUG


void computeTargetPositions(const Eigen::MatrixXd& V1,
                            const Eigen::MatrixXd& V2,
                            const Seam& seam,
                            Eigen::MatrixXd& targets_p,
                            Eigen::MatrixXd& targets_q){

    Eigen::MatrixXd points_p(seam.corres.size(), 2);
    Eigen::MatrixXd points_q(seam.corres.size(), 2);

    for (int i=0; i<seam.corres.size(); i++){
        points_p.row(i) = V1.row(seam.corres[i].first).leftCols(2);
        points_q.row(i) = V2.row(seam.corres[i].second).leftCols(2);
    }

    Eigen::MatrixXd R_est;
    Eigen::VectorXd T_est;
    procustes(points_p.leftCols(2), points_q.leftCols(2), R_est, T_est);

    Eigen::MatrixXd points_pt = points_p.transpose();
    Eigen::MatrixXd points_qt = points_q.transpose();

    Eigen::MatrixXd points_qpt = points_qt.colwise() - T_est;
    points_qpt = (R_est.transpose() * points_qpt);

    Eigen::MatrixXd points_ppt = (R_est * points_pt);
    points_ppt = points_ppt.colwise() + T_est;

    Eigen::MatrixXd points_pp = points_ppt.transpose();
    Eigen::MatrixXd points_qp = points_qpt.transpose(); 

    targets_p = (points_qp + points_p) / 2;
    targets_q = (points_pp + points_q) / 2;
}

void computeTargetPositions(const std::vector<std::unique_ptr<ClothParam>>& cloth_ps, 
                            const Seam& seam,
                            Eigen::MatrixXd& targets_p,
                            Eigen::MatrixXd& targets_q){
    computeTargetPositions(cloth_ps[seam.patch1_id]->getV2d(), 
                           cloth_ps[seam.patch2_id]->getV2d(), seam, targets_p, targets_q);
}

// note: this doesn't take care of packing! All patches are
// more or less centered around (0,0)
bool finalParamMultiPatch(const std::vector<Eigen::MatrixXd>& vec_V_3d, 
                          const std::vector<Eigen::MatrixXi>& vec_F,
                          const std::vector<std::vector<std::vector<std::pair<int, int>>>>& vec_dart_duplicates,
                          const std::vector<std::vector<int>>& vec_dart_tips,
                          const std::vector<Seam>& seams,
                          std::vector<Eigen::MatrixXd>& vec_V_2d,
                          CLOTH_INIT_TYPE init_type){

    #ifdef USE_JSON_CONFIG
    std::string config_path = "../configs/final_config.json";
    nlohmann::json config;
    std::ifstream input_config(config_path);
    input_config >> config;
    #endif
    
    int n_patches = vec_V_3d.size();
    if (n_patches != vec_F.size()){
        std::cout << "ERROR: non matching sizes in final param" << std::endl;
    }

    Eigen::VectorXi per_patch_seam_size = Eigen::VectorXi::Zero(n_patches); 
    for (Seam seam: seams){
        per_patch_seam_size[seam.patch1_id] += seam.corres.size();
        per_patch_seam_size[seam.patch2_id] += seam.corres.size();
    }

    #ifdef MULTI_PARAM_DEBUG
    std::cout << "per_patch_seam_size "<< per_patch_seam_size << std::endl; 
    #endif

    std::vector<std::unique_ptr<ClothParam>> cloth_ps;
    for (int patch_id=0; patch_id<n_patches; patch_id++){
        std::unique_ptr<ClothParam> ptr(new ClothParam(vec_V_3d[patch_id], 
                           vec_F[patch_id], 0.0, 
                           vec_dart_duplicates[patch_id], 
                           vec_dart_tips[patch_id],
                           per_patch_seam_size[patch_id],
                           init_type));
        
        #ifdef USE_JSON_CONFIG
        ptr->loadConfig(config_path);

        if (config["debug"]["save_init_params"]){
            std::cout << "Saving init params..." << std::endl;
            Eigen::MatrixXd V2_init = ptr->getV2d();
            std::string path = config["debug"]["save_init_paths"];
            igl::writeOBJ(path + "init_param_" + std::to_string(patch_id) + ".obj", V2_init, vec_F[patch_id]);
        }
        #endif

        cloth_ps.push_back(std::move(ptr));
    }

    int max_iter = 20;
    #ifdef USE_JSON_CONFIG
    max_iter = config["max_iter"];
    #endif

    for (int current_iter = 0; current_iter < max_iter; current_iter++){

        // first index is a patch, second is each of its seam in an arbitrary order
        std::vector<std::vector<Eigen::MatrixXd>> targets;
        std::vector<std::vector<Eigen::VectorXi>> v_ids;
        targets.resize(cloth_ps.size());
        v_ids.resize(cloth_ps.size());

        for (int seam = 0; seam < seams.size(); seam ++){
            #ifdef MULTI_PARAM_DEBUG
            std::cout << "seam: "<< seam << std::endl; 
            #endif
            Eigen::MatrixXd targets_p;
            Eigen::MatrixXd targets_q;
            computeTargetPositions(cloth_ps, seams[seam], targets_p, targets_q);

            int n_pairs = seams[seam].corres.size();
            Eigen::VectorXi p_ids(n_pairs), q_ids(n_pairs);

            for (int i=0; i<n_pairs; i++){
                p_ids(i) = seams[seam].corres[i].first;
                q_ids(i) = seams[seam].corres[i].second;
            }

            targets[seams[seam].patch1_id].push_back(targets_p);
            v_ids[seams[seam].patch1_id].push_back(p_ids);

            targets[seams[seam].patch2_id].push_back(targets_q);
            v_ids[seams[seam].patch2_id].push_back(q_ids);

            //cloth_ps[seams[seam].patch1_id]->setSeamTargets(targets_p, p_ids);
            //cloth_ps[seams[seam].patch2_id]->setSeamTargets(targets_q, q_ids);
        }

        #ifdef MULTI_PARAM_DEBUG
        std::cout << "targets.size() "<< targets.size() << std::endl; 
        #endif

        for (int patch_id=0; patch_id < cloth_ps.size(); patch_id ++){
            cloth_ps[patch_id]->setSeamTargets(targets[patch_id], v_ids[patch_id]);
        }

        for (int patch_id=0; patch_id<n_patches; patch_id++){
            #ifdef MULTI_PARAM_DEBUG
            std::cout << "patch_id "<< patch_id << std::endl; 
            #endif
            Eigen::MatrixXd prev = cloth_ps[patch_id]->getV2d();
            cloth_ps[patch_id]->paramIter(1);
            Eigen::MatrixXd next = cloth_ps[patch_id]->getV2d();
            if (std::isnan(next.maxCoeff()) || std::isnan(next.minCoeff())){
                cloth_ps[patch_id]->setV2d(prev);
                std::cout << "WARNING, nan, revert to previous solution" << std::endl;
            }
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

    // force seam reflectability in final solution (debug)
    bool force_reflectability = false;
    if (force_reflectability){
        for (int seam = 0; seam < seams.size(); seam ++){ 
            Eigen::MatrixXd targets_p, targets_q;
            computeTargetPositions(cloth_ps, seams[seam], targets_p, targets_q);

            int patch1 = seams[seam].patch1_id;
            int patch2 = seams[seam].patch2_id;
            for (int i=0; i<seams[seam].corres.size(); i++){
                vec_V_2d[patch1](seams[seam].corres[i].first, 0) = targets_p(i,0);
                vec_V_2d[patch1](seams[seam].corres[i].first, 1) = targets_p(i,1);
                vec_V_2d[patch1](seams[seam].corres[i].first, 2) = 0;
                vec_V_2d[patch2](seams[seam].corres[i].second, 0) = targets_q(i,0);
                vec_V_2d[patch2](seams[seam].corres[i].second, 1) = targets_q(i,1);
                vec_V_2d[patch2](seams[seam].corres[i].second, 2) = 0;
            }
        }
    }
    
    double length_error, reflec_error;
    measureSeamScore(vec_V_2d, seams, length_error, reflec_error);

    return all_satisfied;
}