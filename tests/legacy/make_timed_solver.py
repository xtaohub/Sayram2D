#!/usr/bin/env python3
"""Create a timing-instrumented copy of source/Solver.cc for benchmarking.

Only std::chrono instrumentation is inserted; the numerical code is copied
verbatim. The generated file is a build artifact and is not committed.
"""

from pathlib import Path
import sys

src = Path("source/Solver.cc").read_text()

needle = '#include "BCTypes.h"\n'
if needle not in src:
    raise SystemExit("Solver.cc: include line not found")

instrumentation = r'''
#include <chrono>

double g_legacy_assemble_sec = 0.0;
double g_legacy_factorize_sec = 0.0;
double g_legacy_solve_sec = 0.0;
double g_legacy_scatter_sec = 0.0;
double g_legacy_post_sec = 0.0;

static double legacy_now_sec() {
  using clock = std::chrono::steady_clock;
  return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}
'''

src = src.replace(needle, needle + instrumentation, 1)

old_assemble = """  R_.setZero();
  M_coeffs_.clear();

  assemble();

  M_.makeCompressed();
  eigen_solver.factorize(M_);
  ftmp_ = eigen_solver.solve(R_);
"""

new_assemble = """  R_.setZero();
  M_coeffs_.clear();

  const double t_assemble0 = legacy_now_sec();
  assemble();
  g_legacy_assemble_sec += legacy_now_sec() - t_assemble0;

  const double t_factor0 = legacy_now_sec();
  M_.makeCompressed();
  eigen_solver.factorize(M_);
  g_legacy_factorize_sec += legacy_now_sec() - t_factor0;

  const double t_solve0 = legacy_now_sec();
  ftmp_ = eigen_solver.solve(R_);
  g_legacy_solve_sec += legacy_now_sec() - t_solve0;
"""

if old_assemble not in src:
    raise SystemExit("Solver.cc: update assembly block not found")
src = src.replace(old_assemble, new_assemble, 1)

old_scatter = """  for (std::size_t i=0; i<m.nx(); ++i){
    for (std::size_t j=0; j<m.ny(); ++j) {
      f_(i,j) = ftmp_(m.flatten_cell_index({i, j}));
    }
  }

  istep_ += 1;
  eq.update(t());
  update_Lambda();
  update_vertex_f();
"""

new_scatter = """  const double t_scatter0 = legacy_now_sec();
  for (std::size_t i=0; i<m.nx(); ++i){
    for (std::size_t j=0; j<m.ny(); ++j) {
      f_(i,j) = ftmp_(m.flatten_cell_index({i, j}));
    }
  }
  g_legacy_scatter_sec += legacy_now_sec() - t_scatter0;

  const double t_post0 = legacy_now_sec();
  istep_ += 1;
  eq.update(t());
  update_Lambda();
  update_vertex_f();
  g_legacy_post_sec += legacy_now_sec() - t_post0;
"""

if old_scatter not in src:
    raise SystemExit("Solver.cc: update scatter block not found")
src = src.replace(old_scatter, new_scatter, 1)

out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/sayram-timed/Solver.cc")
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(src)
print(out)
