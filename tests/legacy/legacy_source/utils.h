#ifndef UTILS_H_
#define UTILS_H_

#include <cmath>
#include "common.h"

inline double p2e(double p, double E0){ // convert momentum to energy
  return sqrt(p * p * gC * gC + E0 * E0) - E0; 
}


inline double e2p(double E, double E0){ // convert energy to momentum
  return sqrt(E * (E + 2 * E0)) / gC; 
}


inline double dlogE_dp(double logE, double E0){
  double E = exp(logE);
  return e2p(E, gE0) * gC * gC / (E * (E + E0));
}


// for common structured mesh cell interpolation
inline double interp2D(const Xtensor2d& raw, const Loc& loc) {
  int i0,j0;
  double wi,wj; 

  i0 = loc.i0;
  j0 = loc.j0;
  wi = loc.wi;
  wj = loc.wj; 

  return raw(i0,j0)*wi*wj + raw(i0+1,j0)*(1-wi)*wj + raw(i0+1,j0+1)*(1-wi)*(1-wj) + raw(i0,j0+1)*wi*(1-wj);
        
}


// for common weight calculation including range judge
inline void calWeight(std::size_t &i, double &w, std::size_t n, double pos){
  if (i >= 0 && i < n) {
    w = 1.0 - (pos - i); 
  } 
  else if (i < 0) {
    i = 0;
    w = 1.0;
  }
  else if (i >= n) {
    i = n - 1; 
    w = 0.0; 
  }
}

#endif
