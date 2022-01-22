#include "param/metrics.h"


// TODO this is a duplicate of BaryOptimizer::measureScore
void measureStretchScore(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d, 
                         const Eigen::MatrixXi& F, Eigen::VectorXd& stretch_u_vec,
                         Eigen::VectorXd& stretch_v_vec){

    double stretch_u_score = 0;
    double stretch_v_score = 0;
    double edges_score = 0;

    stretch_u_vec.resize(F.rows());
    stretch_v_vec.resize(F.rows());
    
    for (int f_id=0; f_id<F.rows(); f_id++) {
        // each triangle gives us 2 equations:
        // one target u and one target v

        Eigen::RowVectorXd D = (V_2d.row(F(f_id, 0)) + V_2d.row(F(f_id, 1)) + V_2d.row(F(f_id, 2)))/3.0; // centroid
        Eigen::Vector3d D_bary(0.333333, 0.333333, 0.333333);
        Eigen::RowVectorXd DU = D;
        DU(0) += 1.0;
        Eigen::Vector3d DU_bary = barycentricCoords(DU, V_2d.row(F(f_id, 0)), V_2d.row(F(f_id, 1)), V_2d.row(F(f_id, 2)));
        Eigen::RowVectorXd DV = D;
        DV(1) += 1.0;
        Eigen::Vector3d DV_bary = barycentricCoords(DV, V_2d.row(F(f_id, 0)), V_2d.row(F(f_id, 1)), V_2d.row(F(f_id, 2)));
        /*Eigen::RowVectorXd DUV = D; // PERF: deduce DUV eq based on first two?
        DUV(0) += 1.0;
        DUV(1) += 1.0;
        Eigen::Vector3d DUV_bary = barycentricCoords(DUV, V_2d.row(F(f_id, 0)), V_2d.row(F(f_id, 1)), V_2d.row(F(f_id, 2)));*/

        Eigen::RowVectorXd Dp = D_bary(0) * V_3d.row(F(f_id, 0)) + D_bary(1) * V_3d.row(F(f_id, 1)) + D_bary(2) * V_3d.row(F(f_id, 2));
        Eigen::RowVectorXd DUp = DU_bary(0) * V_3d.row(F(f_id, 0)) + DU_bary(1) * V_3d.row(F(f_id, 1)) + DU_bary(2) * V_3d.row(F(f_id, 2));
        Eigen::RowVectorXd DVp = DV_bary(0) * V_3d.row(F(f_id, 0)) + DV_bary(1) * V_3d.row(F(f_id, 1)) + DV_bary(2) * V_3d.row(F(f_id, 2));
        //Eigen::RowVectorXd DUVp = DUV_bary(0) * V_3d.row(F(f_id, 0)) + DUV_bary(1) * V_3d.row(F(f_id, 1)) + DUV_bary(2) * V_3d.row(F(f_id, 2));

        double target_u = (DUp - Dp).norm();
        double target_v = (DVp - Dp).norm();
        //double target_uv_u = (DUVp - Dp).dot(DUp - Dp)/(DUp - Dp).norm();
        //double target_uv_v = (DUVp - Dp).dot(DVp - Dp)/(DVp - Dp).norm();;

        double actual_u = V_2d(F(f_id, 0), 0) * (DU_bary(0) - D_bary(0))
                        + V_2d(F(f_id, 1), 0) * (DU_bary(1) - D_bary(1))
                        + V_2d(F(f_id, 2), 0) * (DU_bary(2) - D_bary(2));
        
        stretch_u_score += std::pow(target_u - actual_u, 2);

        stretch_u_vec(f_id) = actual_u / target_u;
        //stretch_u_score += actual_u / target_u;

        double actual_v = V_2d(F(f_id, 0), 1) * (DV_bary(0) - D_bary(0))
                        + V_2d(F(f_id, 1), 1) * (DV_bary(1) - D_bary(1))
                        + V_2d(F(f_id, 2), 1) * (DV_bary(2) - D_bary(2));
        
        stretch_v_score += std::pow(target_v - actual_v, 2);
        stretch_v_vec(f_id) = actual_v / target_v;
        //stretch_v_score += actual_v / target_v;
            

        // --- BELOW : edges eqs ---
        /*Eigen::MatrixXd V_tri_2d, V_tri_3d;
        makeTriPoints(V_2d, F, f_id, V_tri_2d);
        makeTriPoints(V_3d, F, f_id, V_tri_3d);
        V_tri_3d = move3Dto2D(V_tri_3d);
        Eigen::MatrixXd R_est;
        Eigen::VectorXd T_est;
        procustes(V_tri_2d, V_tri_3d, R_est, T_est);
        Eigen::MatrixXd p2 = V_tri_3d;

        Eigen::MatrixXd p2_rt, p2_r;
        Eigen::MatrixXd p2t = p2.transpose(); // TODO transposeInPlace ?
        p2_rt = p2t.colwise() - T_est;
        p2_rt = (R_est.transpose() * p2_rt);
        p2_r = p2_rt.transpose();

        std::vector<std::pair<int, int>> edges = {std::make_pair(0,1), std::make_pair(0,2), std::make_pair(1,2)}; 

        for (std::pair<int, int> edge : edges){
            Eigen::RowVectorXd Ap = p2_r.row(edge.first);
            Eigen::RowVectorXd Bp = p2_r.row(edge.second);

            double actual_u = V_2d(F(f_id, edge.second), 0) - V_2d(F(f_id, edge.first), 0);
            double target_u = (Bp - Ap)(0);
            edges_score += std::pow(target_u - actual_u, 2);

            double actual_v = V_2d(F(f_id, edge.second), 1) - V_2d(F(f_id, edge.first), 1);
            double target_v = (Bp - Ap)(1);
            edges_score += std::pow(target_v - actual_v, 2);
        }*/
    }

    stretch_u_score /= F.rows();
    stretch_v_score /= F.rows();
    edges_score /= F.rows();

    /*double selected_score = -1;
    if (selected_vs_.size() >= 2 && selected_vs_[0] >= 0 && selected_vs_[1] >= 0){
        // selected vertices should have = V
        int v0 = selected_vs_[0];
        int v1 = selected_vs_[1];
        selected_score = std::pow(V_2d(v0, 1) - V_2d(v1, 1), 2); 
    }
    else { // not enough vertices selected
    }*/

    // No dart score for now
    // equationsFromDarts(V_2d, F); // TODO ?

    /*std::cout << "Stretch U:\t" << stretch_u_score << std::endl;
    std::cout << "Stretch V:\t" << stretch_v_score << std::endl;
    std::cout << "Edges sc.:\t" << edges_score << std::endl;
    std::cout << "Selec sc.:\t" << selected_score << std::endl;*/


    // From centered around 1 to centered around 0:
    stretch_u_vec = stretch_u_vec.array() - 1.0;
    stretch_v_vec = stretch_v_vec.array() - 1.0;
}

void measureSeamScore(const std::vector<Eigen::MatrixXd>& vec_V_2d,
                      const std::vector<Seam>& seams,
                      double& length_error, double& reflec_error){
    
    std::cout << "SEAM METRICS, #seams: " << seams.size() << std::endl;

    std::cout << "Seam patches: " << std::endl;
    for (Seam s: seams){
        std::cout << "\t" << s.patch1_id << " \t" << s.patch2_id << std::endl;
    }

    length_error = 0;
    reflec_error = 0;


    std::cout << "Reflectability error: (m per point)" << std::endl;
    for (Seam s: seams){
        if (s.corres.size() == 0) continue;

        double dist_from_symmetric = 0.0;
        Eigen::MatrixXd targets_p, targets_q;
        computeTargetPositions(vec_V_2d[s.patch1_id], vec_V_2d[s.patch2_id], s, targets_p, targets_q);
        for (int i=0; i<s.corres.size(); i++){
            Eigen::RowVector3d p = vec_V_2d[s.patch1_id].row(s.corres[i].first);
            Eigen::RowVector3d q = vec_V_2d[s.patch2_id].row(s.corres[i].second);
            dist_from_symmetric += (targets_p.row(i) - p.leftCols(2)).norm(); 
            dist_from_symmetric += (targets_q.row(i) - q.leftCols(2)).norm(); 
        }
        double seam_error = dist_from_symmetric / static_cast<double>(s.corres.size() * 2);
        std::cout << "\t" << seam_error << std::endl;
        if (std::isinf(seam_error) || std::isnan(seam_error) || seam_error >= 10e4){
            std::cout << "\tdebug info:" << std::endl;
            std::cout << "\t targets_p " << targets_p << std::endl;
            std::cout << "\t ps: " << std::endl;
            for (int i=0; i<s.corres.size(); i++){
                Eigen::RowVector3d p = vec_V_2d[s.patch1_id].row(s.corres[i].first);
                std::cout << "\t p " << p << std::endl;
            }
            std::cout  << std::endl;
        }
        if (seam_error > reflec_error){
            reflec_error = seam_error; 
        }
    }

    std::cout << "Seam length error: (m)" << std::endl;
    for (Seam s: seams){
        double length1 = 0;
        double length2 = 0;
        Eigen::RowVector3d current1 = vec_V_2d[s.patch1_id].row(s.corres[0].first);
        Eigen::RowVector3d current2 = vec_V_2d[s.patch2_id].row(s.corres[0].second);
        for (int i=1; i<s.corres.size(); i++){
            Eigen::RowVector3d next1 = vec_V_2d[s.patch1_id].row(s.corres[i].first);
            Eigen::RowVector3d next2 = vec_V_2d[s.patch2_id].row(s.corres[i].second);
            length1 += (next1 - current1).norm();
            current1 = next1;
            length2 += (next2 - current2).norm();
            current2 = next2;
        }
        length_error += std::fabs(length1 - length2);
        std::cout << "\t" << std::fabs(length1 - length2) << std::endl;
    }


    std::cout << "Worst seam ref error: " << reflec_error << std::endl;
    std::cout << "Total seam length error: " << length_error << std::endl;
};