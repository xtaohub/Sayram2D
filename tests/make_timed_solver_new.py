#!/usr/bin/env python3
"""Create a timing-instrumented copy of src/solver.cpp for benchmarking.

Only std::chrono instrumentation is inserted; numerical code is copied
verbatim. The generated file is a build artifact and is not committed.
"""

from pathlib import Path
import sys

src = Path("src/solver.cpp").read_text()

needle = '#include "sayram2d/solver.h"\n'
if needle not in src:
    raise SystemExit("solver.cpp: include not found")

instrumentation = r'''
#include <chrono>

double g_new_vertex_sec = 0.0;
double g_new_assemble_sec = 0.0;
double g_new_factorize_sec = 0.0;
double g_new_solve_sec = 0.0;
double g_new_scatter_sec = 0.0;

static double new_now_sec() {
  using clock = std::chrono::steady_clock;
  return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}
'''

src = src.replace(needle, needle + instrumentation, 1)

old = """  update_lambda(equation);
  update_vertex_f(equation, f_in);
  assemble(equation, f_in, dt);

  linear_solver_.factorize(matrix_);
  if (linear_solver_.info() != Eigen::Success) {
    throw std::runtime_error("Solver::solve: sparse factorization failed");
  }
  solution_ = linear_solver_.solve(rhs_);
  if (linear_solver_.info() != Eigen::Success) {
    throw std::runtime_error("Solver::solve: sparse solve failed");
  }

  // All reads of f_in are complete at this point, so f_out may alias f_in.
  for (std::size_t i = 0; i < mesh_.n1(); ++i) {
    for (std::size_t j = 0; j < mesh_.n2(); ++j) {
      f_out(i, j) = solution_(static_cast<int>(mesh_.cell_index({i, j})));
    }
  }
"""

new = """  const double t_vertex0 = new_now_sec();
  update_lambda(equation);
  update_vertex_f(equation, f_in);
  g_new_vertex_sec += new_now_sec() - t_vertex0;

  const double t_assemble0 = new_now_sec();
  assemble(equation, f_in, dt);
  g_new_assemble_sec += new_now_sec() - t_assemble0;

  const double t_factor0 = new_now_sec();
  linear_solver_.factorize(matrix_);
  g_new_factorize_sec += new_now_sec() - t_factor0;
  if (linear_solver_.info() != Eigen::Success) {
    throw std::runtime_error("Solver::solve: sparse factorization failed");
  }

  const double t_solve0 = new_now_sec();
  solution_ = linear_solver_.solve(rhs_);
  g_new_solve_sec += new_now_sec() - t_solve0;
  if (linear_solver_.info() != Eigen::Success) {
    throw std::runtime_error("Solver::solve: sparse solve failed");
  }

  const double t_scatter0 = new_now_sec();
  // All reads of f_in are complete at this point, so f_out may alias f_in.
  for (std::size_t i = 0; i < mesh_.n1(); ++i) {
    for (std::size_t j = 0; j < mesh_.n2(); ++j) {
      f_out(i, j) = solution_(static_cast<int>(mesh_.cell_index({i, j})));
    }
  }
  g_new_scatter_sec += new_now_sec() - t_scatter0;
"""

if old not in src:
    raise SystemExit("solver.cpp: solve block not found")
src = src.replace(old, new, 1)

out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/sayram-new-timed/Solver.cpp")
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(src)
print(out)
