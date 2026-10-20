/*
 * File:        Equation.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        01/13/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef EQUATION_H_
#define EQUATION_H_

#include "common.h"
#include "utils.h"
#include "Mesh.h"
#include "Parameters.h"
#include "BCTypes.h"
#include <cstddef>

//
// Definitions of the Equation are given here
// including D and boundary conditions.
// To handle different cases, we Equation as an abstract class.
//

class Equation{
  public:
    Equation(const Mesh& m) {

      const std::array<std::size_t, 2> shape = {m.nx(), m.ny()};

      G_.resize(shape);

      Dxx_.resize(shape);
      Dyy_.resize(shape);
      Dxy_.resize(shape);

      inv_tau_.resize(shape); // 1/tau: the equation has an extra -f/tau term, representing loss if tau > 0.
      inv_tau_.fill(0.0);    // by setting 1/tau = 0, this term does nothing.

    }

    double G(const Ind& ind) const { return G_(ind.i, ind.j); } 

    double Dxx(const Ind& ind) const { return Dxx_(ind.i, ind.j); }
    double Dyy(const Ind& ind) const { return Dyy_(ind.i, ind.j); } 
    double Dxy(const Ind& ind) const { return Dxy_(ind.i, ind.j); } 

    double inv_tau(const Ind& ind) const { return inv_tau_(ind.i, ind.j); } 

    virtual BCType bc_type(BoundaryID side) const = 0;

    // pure virtual functions to be defined by the CASE of interest. 
    virtual double init_f(const Ind& ind) const = 0;

    virtual void update(double t) = 0;

    // Dirichlet value at a boundary vertex on side
    virtual bool dirichlet_vertex_value(BoundaryID side,
                                        std::size_t i, std::size_t j, double t,
                                        double* out_value) const {
      (void)side; (void)i; (void)j; (void)t; (void)out_value;
      return false;
    }

  protected:
    Xtensor2d G_; 

    Xtensor2d Dxx_; 
    Xtensor2d Dyy_; 
    Xtensor2d Dxy_;

    Xtensor2d inv_tau_;  // 1/tau
};


#endif /* EQUATION_H */

