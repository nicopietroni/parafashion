#pragma once
#include <Eigen/Core>
#include <vector>

std::vector<int> autoSelect(const Eigen::MatrixXd& V_3d, const Eigen::VectorXi& bnd);