#pragma once
#include <Eigen/Core>
#include <vector>

std::vector<int> autoSelect(const Eigen::MatrixXd& V_3d, const Eigen::VectorXi& bnd){
    const int x_axis = 0;
    const int y_axis = 2;
    const int z_axis = 1;

    double avg_x = 0;
    double avg_y = 0;
    for (int i=0; i<bnd.rows(); i++){
        avg_x += V_3d(bnd(i), x_axis);
        avg_y += V_3d(bnd(i), y_axis);
    }
    avg_x /= bnd.rows();
    avg_y /= bnd.rows();


    std::vector<int> best_pair = {-1, -1};
    double best_score = -1;
    for (int i=0; i<bnd.rows(); i++){
        double avg_score_i = std::pow(V_3d(bnd(i), x_axis)  - avg_x, 2) 
                           + std::pow(V_3d(bnd(i), y_axis) - avg_y, 2);
        //if (avg_score_i > best_score && best_score != -1) continue; // score already too bad
        for (int j=0; j<bnd.rows(); j++){
            if (i == j) continue;
            double avg_score_j = std::pow(V_3d(bnd(j), x_axis)  - avg_x, 2) 
                               + std::pow(V_3d(bnd(j), y_axis)  - avg_y, 2);
            double delta_x = std::fabs(V_3d(bnd(i), x_axis) - V_3d(bnd(j), x_axis));
            double delta_y = std::fabs(V_3d(bnd(i), y_axis) - V_3d(bnd(j), y_axis));
            double delta_z = std::fabs(V_3d(bnd(i), z_axis) - V_3d(bnd(j), z_axis));
            //double delta_score = (1.0 + delta_x * delta_y) / (0.01 + delta_z);
            double delta_score = (1.0 + 1000.0 * delta_x + delta_y) / (0.0001 + delta_z);
            double total_score = avg_score_i + avg_score_j + 10.0 * delta_score;
            if (total_score < best_score || best_score == -1){
                best_score = total_score;
                best_pair = {bnd(i), bnd(j)};
            }
        }
    } 

    std::cout << best_pair[0] << " " << best_pair[1] << std::endl;
    return best_pair;
}