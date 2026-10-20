/*
 * File:        Albert_Young.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        01/13/2025 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#include "Albert_Young.h"


Albert_Young::Albert_Young(const Parameters& paras_in, const Mesh& mesh_in):Equation(mesh_in),paras(paras_in), m(mesh_in), io(paras_in){ // initialize G, D, and etc.
  init();
  constructD(paras, io);
}


double Albert_Young::init_f(const Ind& ind) const { 
  double x, y;
  x = m.x(ind.i);
  y = m.y(ind.j);

  return calculate_init_f(x, y);
}


void Albert_Young::init(){

  double x, y;

  for (std::size_t i = 0; i < m.nx(); i++){
    x = m.x(i);
    for (std::size_t j = 0; j < m.ny(); j++){
      y = m.y(j);
      G_(i,j) = calculate_G(x, y);
    }
  }
}

BCType Albert_Young::bc_type(BoundaryID side) const {
  switch (side) {
    case BoundaryID::XMIN:
      return BCType::Dirichlet;

    case BoundaryID::YMIN:
      return BCType::Dirichlet;

    case BoundaryID::YMAX:
      return BCType::Dirichlet;

    case BoundaryID::XMAX:
      return BCType::ZeroFlux;

    default:
      throw std::runtime_error("Albert_Young::bc_type: unknown BoundaryID");
  }
}


void Albert_Young::update(double t) {}

bool Albert_Young::dirichlet_vertex_value(BoundaryID side,
                                          std::size_t i, std::size_t j, double t,
                                          double* out_value) const {
  (void)t;  // remove if time-dependent
  //
  const double x = m.x_edge(i);
  const double y = m.y_edge(j);

  switch (side) {
    case BoundaryID::XMIN:
      *out_value = xmin(y);
      return true;

    case BoundaryID::YMIN:
      *out_value = ymin(x);
      return true;

    case BoundaryID::YMAX:
      *out_value = ymax(x);
      return true;

    case BoundaryID::XMAX:
      // right boundary zeroflux in Albert_Young
      return false;

    default:
      throw std::runtime_error("Albert_Young::dirichlet_value: unknown BoundaryID");
  }
}

void Albert_Young::locate(double alpha0, double logE, Loc* locp){
    std::size_t i0, j0;
    double wi, wj;

    double pos_x = (alpha0 - io.xmin_D()) / (io.xmax_D() - io.xmin_D()) * (io.nx_D() - 1);
    double pos_y = (logE - log(io.ymin_D())) / (log(io.ymax_D()) - log(io.ymin_D())) * (io.ny_D() - 1);

    i0 = floor(pos_x); 
    j0 = floor(pos_y); 

    calWeight(i0, wi, io.nx_D()-1, pos_x);
    calWeight(j0, wj, io.ny_D()-1, pos_y);

    locp->i0 = i0;
    locp->j0 = j0; 
    locp->wi = wi;
    locp->wj = wj; 
}

void Albert_Young::constructD(const Parameters& par, const Albert_Young_IO& io){

  const double denormalize_factor = gME * gME * gC * gC;
  const double second_to_day = 3600 * 24;
  double alpha0, logE, p;

  Loc loc; 

  for(std::size_t i = 0; i < m.nx(); i++){
    alpha0 = m.x(i);
    for(std::size_t j = 0; j < m.ny(); j++){
        logE = m.y(j);
        p = e2p(exp(logE), gE0);

        locate(alpha0, logE, &loc);  

        Dxx_(i,j) = interp2D(io.Dxx_raw, loc) * denormalize_factor * second_to_day / (p*p); 
        Dxy_(i,j) = interp2D(io.Dxy_raw, loc) * denormalize_factor * second_to_day * dlogE_dp(logE, gE0) / p;
        Dyy_(i,j) = interp2D(io.Dyy_raw, loc) * denormalize_factor * second_to_day * pow(dlogE_dp(logE, gE0),2); 
      }
  }
}







