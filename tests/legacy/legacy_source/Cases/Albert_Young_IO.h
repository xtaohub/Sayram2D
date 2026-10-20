/*
 * File:        Albert_Young_IO.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        01/13/2025 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef ALBERT_YOUNG_IO_H_
#define ALBERT_YOUNG_IO_H_


#include "common.h"
#include "Parameters.h"

class Albert_Young_IO{
  public:
    Albert_Young_IO(const Parameters& paras);

    Xarray1d x_D; // angle(°) convert to radial
    Xarray1d y_D; // Energy(keV)

    Xtensor2d Dxx_raw;
    Xtensor2d Dxy_raw;
    Xtensor2d Dyy_raw;

    std::size_t nx_D() const { return nx_D_; }
    std::size_t ny_D() const { return ny_D_; }

    double xmin_D() const { return xmin_D_; }
    double xmax_D() const { return xmax_D_; }
    double ymin_D() const { return ymin_D_; }
    double ymax_D() const { return ymax_D_; }

    void update(double t) {} // do nothing

  private:

    const Parameters& paras; 

    std::size_t nx_D_, ny_D_;

    double xmin_D_, xmax_D_, ymin_D_, ymax_D_;

    void read_D(); 
}; 



#endif /* ALBERT_YOUNG_IO_H */