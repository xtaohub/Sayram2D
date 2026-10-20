/*
 * File:        Albert_Young_IO.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        01/13/2025 
 * 
 * Copyright (c) Xin Tao 
 *
 */
#include "Albert_Young_IO.h"

Albert_Young_IO::Albert_Young_IO(const Parameters& paras_in) : paras(paras_in){

    read_D();
}

void Albert_Young_IO::read_D(){

    std::string dfile_base = "D/" + paras.dID() + ".h5";

    x_D = xt::load_hdf5<Xarray1d>(dfile_base, "/alpha0");
    x_D = x_D * gPI / 180; // convert to radial
    y_D = xt::load_hdf5<Xarray1d>(dfile_base, "/E");

    nx_D_ = x_D.size();
    ny_D_ = y_D.size();

    xmin_D_ = x_D[0];
    xmax_D_ = x_D[nx_D_-1];
    ymin_D_ = y_D[0];
    ymax_D_ = y_D[ny_D_-1];

    Dxx_raw = xt::load_hdf5<Xtensor2d>(dfile_base, "/Daa");
    Dxy_raw = xt::load_hdf5<Xtensor2d>(dfile_base, "/Dap");
    Dyy_raw = xt::load_hdf5<Xtensor2d>(dfile_base, "/Dpp");
  
}