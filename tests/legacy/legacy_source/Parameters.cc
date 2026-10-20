#include "Parameters.h"

Parameters::Parameters(int argc, char** argv){
    handle_main_input(argc, argv);
    read_inp_file();

    output_path_ = "./output/" + run_id() + "/";
    std::filesystem::create_directories(output_path_);

    // copy the parameter file
    std::string paras_file = "./output/" + run_id() + "/" + run_id() + ".ini";
    std::string command = "cp " + inp_file() + " " +  paras_file;
    int result = system(command.c_str());
    if(result != 0){
        std::cerr << "Command failed with " << result << std::endl;
    }
}

void Parameters::handle_main_input(int argc, char* argv[]){
  switch (argc) {
  case 1:
    inp_file_ = "p.ini";
    break;

  case 2:
    inp_file_.assign(argv[1]);
    break;

  default:
    std::cerr << "NParas() > 2! This program takes at most one argument: the parameter file name." << std::endl;
    exit(1);
  }
}

void Parameters::read_inp_file(){

  Ini_reader ireader(inp_file());

  ireader.set_section("basic");

  ireader.read("run_id", &run_id_);
  ireader.read("nalpha0", &nalpha0_);
  ireader.read("nE", &nE_);

  ireader.read("alpha0min", &alpha0_min_);
  ireader.read("alpha0max", &alpha0_max_);
  ireader.read("Emin", &Emin_);
  ireader.read("Emax", &Emax_);

  logEmin_ = log(Emin_);
  logEmax_ = log(Emax_);

  ireader.read("T", &T_);
  ireader.read("nsteps", &nsteps_);

  ireader.set_section("diagnostics");

  ireader.read("nplots", &nplots_);
  save_every_step_ = nsteps_ / nplots_;
  nsteps_ = save_every_step_ * nplots_;

  ireader.set_section("diffusion_coefficients");
  ireader.read("dID", &dID_);

}

