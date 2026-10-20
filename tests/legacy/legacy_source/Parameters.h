/*
 * File:        Parameters.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        01/13/2024
 *
 * Copyright (c) Xin Tao
 *
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <filesystem>
#include "common.h"
#include "Ini_reader.h"


// For your project, if you need extra 
// information to be read in from the input 
// file, you can define your own class which 
// inherits the Parameters class. 

class Parameters{
public:
  Parameters(int argc, char** argv);

  // ------- READ FROM INP FILE -------
  const std::string& inp_file() const { return inp_file_; }
  const std::string& run_id() const { return run_id_; }

  std::size_t nalpha0() const { return nalpha0_; }
  std::size_t nE() const { return nE_; }

  double alpha0_min() const { return alpha0_min_ * gPI / 180; }
  double alpha0_max() const { return alpha0_max_ * gPI / 180; }
  double Emin() const { return Emin_; }
  double Emax() const { return Emax_; }
  double logEmin() const { return logEmin_; }
  double logEmax() const { return logEmax_; }

  double T() const { return T_; }
  int nsteps() const { return nsteps_; }
  double dt() const { return T_ / nsteps_;}

  int nplots() const { return nplots_; }
  int save_every_step() const { return save_every_step_; }
  const std::string& output_path() const { return output_path_; }
  const std::string& dID() const { return dID_; }

private:
  std::string inp_file_;
  std::string run_id_;

  std::size_t nalpha0_;
  std::size_t nE_;

  double alpha0_min_;
  double alpha0_max_;
  double Emin_;
  double Emax_;
  double logEmin_;
  double logEmax_;

  double T_;
  double nsteps_;

  int nplots_;
  int save_every_step_;
  std::string output_path_;

  std::string dID_;

  void handle_main_input(int argc, char* argv[]);
  void read_inp_file();
};

#endif /* PARAMETERS_H_ */
