#pragma once
#include <Eigen/Core>
#include <iostream>

/* FROM 2 TO 1 
line3 = line2t.colwise() - T_est;
line3 = (R_est.transpose() * line3);
*/

/* FROM 1 TO 2 
line3 = (R_est * line3);
line3 = line3.colwise() + T_est;
*/

void procustes(const Eigen::MatrixXd& points1, // to
               const Eigen::MatrixXd& points2, // from
               Eigen::MatrixXd& R_est,
               Eigen::VectorXd& T_est);

// TODO Matrix3d ?
Eigen::MatrixXd makeTriPoints(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int f_id);
void makeTriPoints(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int f_id, Eigen::MatrixXd& out);


// Move triangle from 3D (Z != 0) to 2D (Z = 0), with
// same lengths but arbitrary orientation
Eigen::MatrixXd move3Dto2D(const Eigen::MatrixXd& V_tri);

Eigen::VectorXd vertices2dToVector(const Eigen::MatrixXd& V);

// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
// credits https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
Eigen::Vector3d barycentricCoords(const Eigen::RowVector3d& p, const Eigen::RowVector3d& a, 
                                         const Eigen::RowVector3d& b, const Eigen::RowVector3d& c);

Eigen::MatrixXd paramARAP(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F, const Eigen::VectorXi& bnd);
Eigen::MatrixXd paramLSCM(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F, const Eigen::VectorXi& bnd);
Eigen::MatrixXd paramSCAF(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F, const Eigen::VectorXi& bnd);

Eigen::Matrix3d computeRotation(const Eigen::RowVector3d& from,
                                const Eigen::RowVector3d& to);

Eigen::Matrix3d rotationVote(const Eigen::MatrixXd& V_3d,
                             const Eigen::MatrixXd& V_2d,
                             const Eigen::MatrixXi& F,
                             const Eigen::RowVector3d& target_3d,
                             const Eigen::RowVector3d& target_2d);