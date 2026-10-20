/*
 * fin:        main.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 *
 * Date:        05/12/2024
 *
 * Copyright (c) Xin Tao
 *
 */

#include <fstream>
#include <iostream>
#include "Parameters.h"
#include "Mesh.h"
#include "Albert_Young.h"
#include "Solver.h"
#include <ctime>

static Grid2D make_uniform(const Parameters& p) {
  std::vector<double> xe(p.nalpha0() + 1);
  std::vector<double> ye(p.nE() + 1);

  const double dx =
      (p.alpha0_max() - p.alpha0_min()) / static_cast<double>(p.nalpha0());
  const double dy =
      (p.logEmax() - p.logEmin()) / static_cast<double>(p.nE());

  for (std::size_t i = 0; i <= p.nalpha0(); ++i) {
    xe[i] = p.alpha0_min() + dx * static_cast<double>(i);
  }
  for (std::size_t j = 0; j <= p.nE(); ++j) {
    ye[j] = p.logEmin() + dy * static_cast<double>(j);
  }

  return Grid2D(std::move(xe), std::move(ye));
}

int main(int argc, char** argv) {

  Parameters paras(argc,argv);

  Grid2D grid = make_uniform(paras);

  // Create mesh
  Mesh m(grid, paras.dt());

  // Create Equation object
  Albert_Young eq(paras, m);

  Solver solver(m, &eq);

  std::string filename;
  std::ofstream out;

  // Create output H5 file
  filename = paras.output_path() + "/" + paras.run_id() + "_data.h5";
  HighFive::File file(filename, HighFive::File::Overwrite);

  Eigen::Map<const Eigen::VectorXd> x_eig(m.x().data(), m.x().size());
  Eigen::Map<const Eigen::VectorXd> y_eig(m.y().data(), m.y().size());

  const Eigen::VectorXd alpha0 = x_eig / gPI * 180.0;
  const Eigen::VectorXd logEN  = y_eig.array() - std::log(gE0);

  xt::dump(file, "/alpha0", alpha0, xt::dump_mode::overwrite);
  xt::dump(file, "/logEN", logEN, xt::dump_mode::overwrite);

  // The timer
  clock_t start, end;
  double cpu_time;
  start = clock();

  xt::dump(file, "/f/0", solver.f(), xt::dump_mode::overwrite);

  // Time loop for solving
  for (int tstep = 1; tstep <= paras.nsteps(); ++tstep) {

    // Solve using FVM solver
    solver.update();

    if(tstep % paras.save_every_step() == 0){
      xt::dump(file, "/f/" + std::to_string(int((tstep) / paras.save_every_step())), solver.f(), xt::dump_mode::overwrite);
    }
  }

  Xarray1d t = xt::linspace(0.0, paras.T(), paras.nplots()+1);

  xt::dump(file, "/t", t, xt::dump_mode::overwrite);

  end = clock();
  cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "CPU time used " << cpu_time << " seconds" << std::endl;

  return 0;
}
