#ifndef COMMON_H_
#define COMMON_H_

#define EIGEN_NO_DEBUG

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <limits>
#include "Eigen/Dense"
#include "Eigen/Sparse" 
#include "Eigen/Core"
#include "Eigen/SparseLU"
#include "xtensor/xarray.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xio.hpp"
#include "xtensor-io/xhighfive.hpp"

typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMat;
typedef Eigen::Triplet<double> T;
typedef Eigen::Vector2d Vector2;
typedef Eigen::Vector2d Point; 

typedef xt::xarray<double> Xarray1d;
typedef xt::xtensor<double, 2> Xtensor2d; 

struct Loc{
  int i0; 
  int j0; 
  double wi; 
  double wj; 
}; 

// constants
const double gEPS = std::numeric_limits<double>::epsilon();
const double gPI = 3.141592653589793238462;
const double gD2R = gPI / 180.0; // convert degree to radian
const double gC = 1;
const double gE0 = 0.511875; // MeV
const double gME = gE0 / (gC * gC); 
const double gRE = 6371000;

using namespace std;  //bad practice, used here for convenience

#endif
