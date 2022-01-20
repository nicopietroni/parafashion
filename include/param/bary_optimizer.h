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
    BaryOptimizer(){};
    void allocateMemory(int n_faces, int n_vs);

    Eigen::MatrixXd localGlobal(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d, 
                                const Eigen::MatrixXi& F, int max_iterations);
    Eigen::MatrixXd localGlobal(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d, 
                                const Eigen::MatrixXi& F);

    static void measureScore(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d, 
                      const Eigen::MatrixXi& F, Eigen::VectorXd& stretch_u_vec,
                      Eigen::VectorXd& stretch_v_vec);

    // ------------- CONFIG parameters ------------- // 
    // Penalize stretch on U/V axes 
    bool enable_stretch_eqs_ = true;
    double stretch_coeff_ = 10.0;

    // Place first vertex in (0,0)
    bool enable_set_seed_eqs_ = false;
    double seed_coeff_ = 0.0001;

    // Preserve mesh edges projections on U/V axes
    bool enable_edges_eqs_ = true;
    double edges_coeff_ = 1.0;

    // Align V axis with pair of selected vertices
    bool enable_selected_eqs_ = false;
    double selected_coeff_ = 20.0;

    // Try to align each triangle vertically
    bool enable_tri_align_eqs_ = true;
    double tri_align_coeff_ = 0.1;

    // Encourage pairs of points on dart to be symmetric
    bool enable_dart_sym_eqs_ = true;
    double dart_sym_coeff_ = 1.0;

    // Match target positions for seams
    bool enable_seam_eqs_ = true;
    double seam_coeff_ = 000.0;//10.0;

    // this one doesn't do what it's supposed to
    bool enable_angle_eqs_ = false;
    double angle_coeff_ = 0.7;

    void setCoeffs(double stretch_coeff, double edges_coeff, double selected_coeff, 
                   double tri_align_coeff, double dart_sym_coeff, double seam_coeff){
        stretch_coeff_ = stretch_coeff;
        edges_coeff_ = edges_coeff;
        selected_coeff_ = selected_coeff;
        tri_align_coeff_ = tri_align_coeff;
        dart_sym_coeff_ = dart_sym_coeff;
        seam_coeff_ = seam_coeff;
    }
    // ------------- CONFIG end        ------------- // 

    void setSelectedVertices(std::vector<int> selected_vs) {selected_vs_ = selected_vs;};
    std::vector<int> getSelectedVertices() const {return selected_vs_;}

    void setDarts(std::vector<SimpleDart> simple_darts);
    void setDarts(std::vector<std::vector<int>> ordered_cuts);
    void setUnorderedDarts(const std::vector<std::vector<std::pair<int, int>>>& dart_duplicates,
                           const std::vector<int>& dart_tips);

    void setSeamSize(int seam_size) {
        seam_size_ = seam_size;
    }
    void setSeamTargets(const std::vector<Eigen::MatrixXd>& targets_p,
                        const std::vector<Eigen::VectorXi>& p_ids){
        targets_p_seams_ = targets_p;
        p_ids_seams_ = p_ids;
    }

private:

    int next_equation_id_ = 0;
    int n_equations_ = 0;
    int n_triplets_ = 0;
    std::vector<int> selected_vs_;
    //std::vector<SimpleDart> simple_darts_;
    std::vector<UnorderedDart> unordered_darts_;
    int seam_size_;
    
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
    Eigen::VectorXd x;
    DiagonalMatrixXd W; 
    Eigen::MatrixXd V_tri_2d;
    Eigen::MatrixXd V_tri_3d;

    std::vector<Eigen::MatrixXd> targets_p_seams_;
    std::vector<Eigen::VectorXi> p_ids_seams_;

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
