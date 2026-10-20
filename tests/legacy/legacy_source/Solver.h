/*
 * File:        Solver.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        01/13/2025
 *
 * Copyright (c) Xin Tao
 *
 */

#ifndef SOLVER_H
#define SOLVER_H

#include "common.h"
#include "Mesh.h"
#include "Equation.h"

class Solver {
  public:
    Solver(const Mesh& m_in, Equation* eqp);

    void update();
    double t() const { return istep_ * m.dt(); }
    const Xtensor2d& f() const { return f_; }
    double f(const Ind& ind) const { return f_(ind.i, ind.j); }

  private:
    const Mesh& m;

    Equation& eq;

    std::size_t istep_; // used to calcualte time = istep_ * dt()

    // M f = R
    Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> eigen_solver;

    SpMat M_;
    std::vector<T> M_coeffs_;

    Xtensor2d f_;
    Eigen::VectorXd R_;
    Eigen::VectorXd ftmp_; // used to store results from the Eigen Solver.

    xt::xtensor<Eigen::Matrix2d, 2> Lambda_;  // Lambda_(i,j) would be the Lambda Matrix at cell (i,j)
    const Eigen::Matrix2d& Lambda(int i, int j) const { return Lambda_(i,j); }
    const Eigen::Matrix2d& Lambda(const Ind& ind) const { return Lambda_(ind.i,ind.j); }
    void update_Lambda();

    //
    // build a lookup table for f at vertices. 
    //
    Xtensor2d vertex_f_;

    double vertex_f(std::size_t i, std::size_t j) const { return vertex_f_(i,j); }
    double vertex_f(const Ind& ind) const { return vertex_f_(ind.i, ind.j); }
    void update_vertex_f();

    void assemble();

    // to calculate coefficients alpha_sigma_i and a_sigma_i
    void a_sigma_func(const Ind& ind, int inbr, Vector2* a_sigma_i_p, double* a_sigmap);

    void calculate_mu(double a_sigma_K, double a_sigma_L, double* mu_Kp, double* mu_Lp){
      double denom = std::abs(a_sigma_K) + std::abs(a_sigma_L) + 2*gEPS;
      (*mu_Kp) = (std::abs(a_sigma_L)+gEPS)/denom;
      (*mu_Lp) = 1 - (*mu_Kp);
    }

    void init();

    // update coefficients in M for cell, and its neighbor: inbr 
    // note that coefficients for both cell Ind and its neighbor are updated.
    void apply_inner_face_pair(const Ind& cell, int inbr);

    void apply_dirichlet_face(const Ind& cell, int inbr); 

    void apply_boundary_faces(); 

    void fill_vertex_from_cells();
    void fill_vertex_from_bcs(); 

};

#endif /* SOLVER_H */
