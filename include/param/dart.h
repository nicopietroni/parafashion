#pragma once
#include <Eigen/Core>
#include <vector>
#include <iostream>

class Dart {
protected:
    const std::vector<int> v_ids_; // Ordered vertex ids
    // Eigen::RowVector3d symmetry_axis_; // Note: desired symmetry axis is in 3D rather than 2D because it doesn't cost much and it allows viz in 3D
    // TODO remove symmetry_axis: this class is position agnostic 

public:
    Dart() {};
    Dart(const std::vector<int>& selected
         //const Eigen::MatrixXd& F
         ) : v_ids_(selected) {
        
    };

    //Eigen::RowVector3d getSymmetryAxis() const {return symmetry_axis_;}
};

class SimpleDart : public Dart {
public:

    SimpleDart() {};
    SimpleDart(const std::vector<int>& selected)
        : Dart(selected) {
        if (selected.size() % 2 == 0) std::cout << "ERROR: simple darts should be odd sized" << std::endl;
    };


    int size() const {return v_ids_.size();}
    int tip() const {return v_ids_.size()/2;}
    int tip_id() const {return v_ids_[v_ids_.size()/2];}
    int start_id() const {return v_ids_[0];}
    int end() const {return v_ids_.size() - 1;}
    int end_id() const {return v_ids_[v_ids_.size() - 1];}
    int v_id(int i) const {return v_ids_[i];}

    int symmetric(int id1) const {
        if (id1 > tip()) return end() - id1;
        else if (id1 == tip()) return tip();
        else return end() - id1;
    }

    void print() const {
        std::cout << "Simple dart of size: " << v_ids_.size() << std::endl;
        std::cout << "╔ ";
        for (int i=0; i<tip(); i++){
            std::cout << v_ids_[i] << " ═ ";
        }
        std::cout << std::endl;
        std::cout << v_ids_[tip()] << std::endl;
        std::cout << "╚ ";
        for (int i=end(); i>tip(); i--){
            std::cout << v_ids_[i] << " ═ ";
        }
        std::cout << std::endl;
    }

    void xPattern(const Eigen::MatrixXd& V, 
                  Eigen::MatrixXd& begs, Eigen::MatrixXd& ends) const {
        if (size() < 3) return;
        begs.resize(size() - 3, 3);
        ends.resize(size() - 3, 3);
        begs.row(0) = V.row(v_ids_[0]);
        ends.row(0) = V.row(v_ids_[size() - 2]);
        for (int i=1; i<tip()-1; i++){
            begs.row(2*i - 1) = V.row(v_ids_[i]);
            ends.row(2*i - 1) = V.row(v_ids_[size() - i - 2]);
            begs.row(2*i + 0) = V.row(v_ids_[i]);
            ends.row(2*i + 0) = V.row(v_ids_[size() - i + 0]);
        }
        begs.row(size() - 4) = V.row(v_ids_[tip() - 1]);
        ends.row(size() - 4) = V.row(v_ids_[tip() + 2]);       
    }

    Eigen::RowVector3d computeSymmetryAxis(const Eigen::MatrixXd& V_2d) const {
        if (size() < 3) std::cout << "ERROR: too short to compute sym axis" << std::endl;
        Eigen::RowVector2d semi_point = V_2d.row(tip_id());

        //Eigen::RowVector3d axis = (symm_end - semi_point).normalized();
        //viewer.data().add_edges(semi_point, symm_end, Eigen::RowVector3d(1., 0., 0.));

        // Simple linear regression without the intercept term:
        // find optimal dart symmetry axis, using regression on middle points
        // between duplicated vertices, but constrained to go through
        // dart end
        // https://en.wikipedia.org/wiki/Simple_linear_regression
        Eigen::MatrixXd points(tip(), 2);
        for (int i=0; i<tip(); i++){
            points.row(i) = (V_2d.row(v_ids_[i]) + V_2d.row(v_ids_[symmetric(i)]))/2.0;
        }

        points = points.rowwise() - V_2d.row(tip_id());

        //Eigen::RowVectorXd means = points.colwise().mean();
        Eigen::RowVectorXd product_xy(points.rows());
        Eigen::RowVectorXd product_xx(points.rows());
        for (int i=0; i<points.rows(); i++) {
            product_xy(i) = points(i,0) * points(i,1); 
            product_xx(i) = points(i,0) * points(i,0); 
        }

        double beta = product_xy.mean() / product_xx.mean();
        Eigen::RowVector2d sym_2d = Eigen::RowVector2d(-1, -beta);
        Eigen::RowVector3d symmetry_axis;
        symmetry_axis(0) = sym_2d(0);
        symmetry_axis(1) = sym_2d(1);
        symmetry_axis(2) = 0;

        // orient towards dart

        Eigen::RowVector2d symm_end = (V_2d.row(start_id()) + V_2d.row(end_id()))/2.0;
        Eigen::RowVector2d out_axis = symm_end - semi_point;

        if (sym_2d.dot(out_axis) < 0) symmetry_axis = -symmetry_axis;

        symmetry_axis = symmetry_axis * out_axis.norm() / sym_2d.norm();
        return symmetry_axis;
    }

    /*Eigen::MatrixXd getMiddlePoints(const Eigen::MatrixXd& V) const {
        Eigen::MatrixXd points(tip(), 3);
        for (int i=0; i<tip(); i++){
            points.row(i) = (V.row(v_ids_[i]) + V.row(v_ids_[symmetric(i)]))/2.0;
        }
        return points;
    };*/

    Eigen::MatrixXd getSymmetricPoints(const Eigen::MatrixXd& V_2d,
                                       const Eigen::RowVector3d& sym_axis) const {
        int len = size();
        Eigen::MatrixXd points_2d(len, 2);
        for (int j=0; j<len; j++){
            int i = symmetric(j);
            Eigen::RowVector3d A, B, C, AD, CE;
            A = V_2d.row(tip_id());
            B = A + sym_axis;
            C = V_2d.row(v_ids_[i]);
            AD = (B - A) * (C - A).dot(B - A) / std::pow((B - A).norm(), 2);
            CE = 2 * (AD - (C - A));
            points_2d.row(j) = C + CE;
        }
        return points_2d;
    };

    Eigen::MatrixXd computeSymmetryTargets(const Eigen::MatrixXd& V_2d,
                                           const Eigen::RowVector3d& sym_axis) const {
        Eigen::MatrixXd sym_p = getSymmetricPoints(V_2d, sym_axis);

        Eigen::MatrixXd new_pos(size(), 2);
        for (int i=0; i<size(); i++){
            new_pos.row(i) = (V_2d.row(v_ids_[i]) + sym_p.row(i))/2.0;
        }

        return new_pos;
    }

    void snapSymmetric(Eigen::MatrixXd& V_2d,
                       const Eigen::RowVector3d& sym_axis) const {
        Eigen::MatrixXd new_pos = computeSymmetryTargets(V_2d, sym_axis);

        for (int i=0; i<size(); i++){
            V_2d.row(v_ids_[i]) = new_pos.row(i);
        }
    }
};

class UnorderedDart {
private:
    
public:
    const std::vector<std::pair<int, int>> pairs_;
    int tip_;

    UnorderedDart() {};
    UnorderedDart(const std::vector<std::pair<int, int>>& pairs, int tip)
        : pairs_(pairs), tip_(tip){};

    void print() const {
        std::cout << "Dart:" << std::endl;
        std::cout << "\t tip: " << tip_ << std::endl;
        std::cout << "\t pairs: " << std::endl;
        for (int i=0; i<pairs_.size(); i++){
            std::cout << "\t\t" << pairs_[i].first << " " << pairs_[i].second << std::endl;
        }
    }

    Eigen::RowVector2d computeSymmetryAxis(const Eigen::MatrixXd& V_2d) const {
        if (pairs_.size() < 1) std::cout << "ERROR: dart too short to compute sym axis" << std::endl;
        if (V_2d.cols() != 2) std::cout << "ERROR: expected 2D mesh in computeSymmetryAxis" << std::endl;
        Eigen::RowVector2d semi_point = V_2d.row(tip_);

        //Eigen::RowVector3d axis = (symm_end - semi_point).normalized();
        //viewer.data().add_edges(semi_point, symm_end, Eigen::RowVector3d(1., 0., 0.));

        // Simple linear regression without the intercept term:
        // find optimal dart symmetry axis, using regression on middle points
        // between duplicated vertices, but constrained to go through
        // dart end
        // https://en.wikipedia.org/wiki/Simple_linear_regression
        Eigen::MatrixXd points(pairs_.size(), 2);
        for (int i=0; i<pairs_.size(); i++){
            points.row(i) = (V_2d.row(pairs_[i].first) + V_2d.row(pairs_[i].second))/2.0;
        }

        points = points.rowwise() - V_2d.row(tip_);

        //Eigen::RowVectorXd means = points.colwise().mean();
        Eigen::RowVectorXd product_xy(points.rows());
        Eigen::RowVectorXd product_xx(points.rows());
        for (int i=0; i<points.rows(); i++) {
            product_xy(i) = points(i,0) * points(i,1); 
            product_xx(i) = points(i,0) * points(i,0); 
        }

        double beta = product_xy.mean() / product_xx.mean();
        Eigen::RowVector2d sym_2d = Eigen::RowVector2d(-1, -beta);
        Eigen::RowVector2d symmetry_axis;
        symmetry_axis(0) = sym_2d(0);
        symmetry_axis(1) = sym_2d(1);
        //symmetry_axis(2) = 0;

        // orient towards dart

        Eigen::RowVector2d symm_end = (V_2d.row(pairs_[0].first) + V_2d.row(pairs_[0].second))/2.0;
        Eigen::RowVector2d out_axis = symm_end - semi_point;

        if (sym_2d.dot(out_axis) < 0) symmetry_axis = -symmetry_axis;

        symmetry_axis = symmetry_axis * out_axis.norm() / sym_2d.norm();
        return symmetry_axis;
    }

    std::vector<std::pair<Eigen::RowVector2d, Eigen::RowVector2d>> getSymmetricPoints(const Eigen::MatrixXd& V_2d,
                                                                                      const Eigen::RowVector2d& sym_axis) const {

        
        std::vector<std::pair<Eigen::RowVector2d, Eigen::RowVector2d>> res;
        res.resize(pairs_.size());
        for (int i=0; i<pairs_.size(); i++){
            Eigen::RowVector2d A, B, C, AD, CE;
            A = V_2d.row(tip_);
            B = A + sym_axis;
            C = V_2d.row(pairs_[i].second);
            AD = (B - A) * (C - A).dot(B - A) / std::pow((B - A).norm(), 2);
            CE = 2 * (AD - (C - A));

            res[i].first = C + CE;

            A = V_2d.row(tip_);
            B = A + sym_axis;
            C = V_2d.row(pairs_[i].first);
            AD = (B - A) * (C - A).dot(B - A) / std::pow((B - A).norm(), 2);
            CE = 2 * (AD - (C - A));

            res[i].second = C + CE;
        }

        return res;
    };

    std::vector<std::pair<Eigen::RowVector2d, Eigen::RowVector2d>> computeSymmetryTargets(const Eigen::MatrixXd& V_2d,
                                                                                          const Eigen::RowVector2d& sym_axis) const {
        std::vector<std::pair<Eigen::RowVector2d, Eigen::RowVector2d>> sym_p = getSymmetricPoints(V_2d, sym_axis);
        std::vector<std::pair<Eigen::RowVector2d, Eigen::RowVector2d>> new_pos;
        new_pos.resize(pairs_.size());
        for (int i=0; i<pairs_.size(); i++){
            //new_pos.row(i) = (V_2d.row(v_ids_[i]) + sym_p.row(i))/2.0;
            new_pos[i].first = (V_2d.row(pairs_[i].first) + sym_p[i].first)/2.0;
            new_pos[i].second = (V_2d.row(pairs_[i].second) + sym_p[i].second)/2.0;
        }

        return new_pos;
    }

    void snapSymmetric(Eigen::MatrixXd& V_2d,
                       const Eigen::RowVector2d& sym_axis) const {

        if (V_2d.cols() > 2) std::cout << "Warning, replacing 3D rows with 2D" << std::endl;

        std::vector<std::pair<Eigen::RowVector2d, Eigen::RowVector2d>> new_pos = computeSymmetryTargets(V_2d, sym_axis);

        for (int i=0; i<pairs_.size(); i++){
            V_2d.row(pairs_[i].first) = new_pos[i].first;
            V_2d.row(pairs_[i].second) = new_pos[i].second;
        }

    }
};

class DiamondDart : Dart {

};