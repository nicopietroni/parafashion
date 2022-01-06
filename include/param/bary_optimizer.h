#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

#include "param/dart.h"

typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagonalMatrixXd;

// Interesting discussion on Eigen performance for LSCM:
// https://forum.kde.org/viewtopic.php?f=74&t=125165

// benchmarking eigen solvers for |Ax - b|^2
// https://stackoverflow.com/questions/42116271/best-way-of-solving-sparse-linear-systems-in-c-gpu-possible

// WEIGHTED SOLVE: https://math.stackexchange.com/questions/709602/when-solving-an-overdetermined-linear-system-is-it-possible-to-weight-the-influ
// https://forum.kde.org/viewtopic.php?f=74&t=110784

class BaryOptimizer {
public:
    BaryOptimizer(int n_faces, int n_vs);

    Eigen::MatrixXd localGlobal(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d, 
                                const Eigen::MatrixXi& F);

    void measureScore(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d, 
                      const Eigen::MatrixXi& F, Eigen::VectorXd& stretch_u_vec,
                      Eigen::VectorXd& stretch_v_vec);

    // ------------- CONFIG parameters ------------- // 
    // Penalize stretch on U/V axes 
    bool enable_stretch_eqs_ = true;
    double stretch_coeff_ = 10.0;

    // Place first vertex in (0,0)
    bool enable_set_seed_eqs_ = true;

    // Preserve mesh edges projections on U/V axes
    bool enable_edges_eqs_ = true;
    double edges_coeff_ = 1.0;

    // Align V axis with pair of selected vertices
    bool enable_selected_eqs_ = true;
    double selected_coeff_ = 20.0;

    // Encourage pairs of points on dart to be symmetric
    bool enable_dart_sym_eqs_ = true;
    double dart_sym_coeff_ = 1.0;

    // this one doesn't do what it's supposed to
    bool enable_angle_eqs_ = false;
    double angle_coeff_ = 0.7;
    // ------------- CONFIG end        ------------- // 

    void setSelectedVertices(std::vector<int> selected_vs) {selected_vs_ = selected_vs;};
    void setDarts(std::vector<SimpleDart> simple_darts);
    void setDarts(std::vector<std::vector<int>> ordered_cuts);

private:

    int next_equation_id_ = 0;
    int n_equations_ = 0;
    int n_triplets_ = 0;
    std::vector<int> selected_vs_;
    std::vector<SimpleDart> simple_darts_;

    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
    Eigen::VectorXd x;
    DiagonalMatrixXd W; 
    Eigen::MatrixXd V_tri_2d;
    Eigen::MatrixXd V_tri_3d;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    std::vector<Eigen::Triplet<double>> triplet_list; 
    
    void equationsFromTriangle(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d,
                               const Eigen::MatrixXi& F, int f_id);

    void equationsFromDarts(const Eigen::MatrixXd& V_2d,
                            const Eigen::MatrixXi& F);

    void makeSparseMatrix(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d,
                      const Eigen::MatrixXi& F);

    bool canUseSelectedEquation(){
        return enable_selected_eqs_ && selected_vs_.size() >= 2 && selected_vs_[0] >= 0 && selected_vs_[1] >= 0;
    }
    
    int edges_eq_time = 0;
    int stretch_shear_eq_time = 0;
    int triplets_time = 0;
    int darts_eq_time = 0;
    int seed_sel_eqs_time = 0;
};