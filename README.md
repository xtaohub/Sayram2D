# Sayram2D

Sayram2D is a lightweight C++20 library for two-dimensional Fokker–Planck
diffusion on structured, orthogonal, possibly nonuniform rectangular grids.
It implements the positivity-preserving finite-volume NTPFA discretization
described in [`docs/sayram2d_algorithm.tex`](docs/sayram2d_algorithm.tex).

The redesign separates geometry, coefficients, boundary conditions, and the
numerical solve. The library core depends only on Eigen; it does not depend on
HDF5 or xtensor.

## Components

| Component | Responsibility |
|---|---|
| `sayram2d::Mesh` | Owns `q1`/`q2` edges, centers, and widths. Provides cell centers/measures, vertex positions, face geometry, neighbors, and boundary queries. No time step or physical coefficients. |
| `sayram2d::ArrayView2D<T>` | Minimal non-owning strided 2-D view. Supports row-major, column-major, and noncontiguous slices of larger arrays. |
| `sayram2d::BoundaryConditions` | Per-side `Dirichlet` or `ZeroFlux` configuration. Dirichlet data are specified at boundary vertices, as required by the discretization. |
| `sayram2d::Equation` | Non-owning views of `G`, `D11`, `D12`, `D22`, optional inverse loss time, and a non-owning reference to the caller-owned boundary configuration. Ordinary `Eigen::MatrixXd` coefficients are accepted directly and without copying. |
| `sayram2d::Solver` | Reusable sparse-matrix workspace and NTPFA assembly. One call to `solve` performs the documented single linear solve. |

## Build

Requirements:

- C++20 compiler
- Eigen 3.4 or later (header-only)
- GNU Make and `ar`
- HDF5 with C headers/libraries for the standalone application
  (not required by the numerical library)

The default Makefile assumes Eigen is under `/Users/xtao/local/include`.
Override `EIGEN_INCLUDE` and `CXX` as needed:

```bash
make EIGEN_INCLUDE=/path/to/eigen CXX=g++
```

Build the static library:

```bash
make
```

Run the unit and legacy-regression tests:

```bash
make test
```

Run the benchmark (the target creates a temporary timing-instrumented copy of
the solver; the numerical code is unchanged):

```bash
make benchmark-new
```

## Usage

```cpp
#include <Eigen/Dense>

#include "boundary_conditions.h"
#include "equation.h"
#include "mesh.h"
#include "solver.h"

using namespace sayram2d;

Mesh mesh({0.0, 0.3, 0.7, 1.0}, {0.0, 0.4, 1.0, 1.7});

const std::size_t n1 = mesh.n1();
const std::size_t n2 = mesh.n2();

Eigen::MatrixXd G = Eigen::MatrixXd::Constant(n1, n2, 1.0);
Eigen::MatrixXd D11 = Eigen::MatrixXd::Constant(n1, n2, 1.0);
Eigen::MatrixXd D12 = Eigen::MatrixXd::Zero(n1, n2);
Eigen::MatrixXd D22 = Eigen::MatrixXd::Constant(n1, n2, 1.0);
Eigen::MatrixXd inv_tau = Eigen::MatrixXd::Zero(n1, n2);

BoundaryConditions bc;
bc.set_dirichlet(Side::Q1Min, 0.5);   // constant value
bc.set_zero_flux(Side::Q1Max);
bc.set_zero_flux(Side::Q2Min);
bc.set_zero_flux(Side::Q2Max);

Equation equation(G, D11, D12, D22, bc);
// Optional; inverse loss defaults to zero.
equation.set_loss(inv_tau);

Solver solver(mesh);

Eigen::MatrixXd f_in = Eigen::MatrixXd::Ones(n1, n2);
Eigen::MatrixXd f_out(n1, n2);
solver.solve(equation, f_in, f_out, 0.01);
```

`Eigen::MatrixXd` is supported directly (zero-copy, column-major). Only
ordinary `MatrixXd` is accepted this way; `Eigen::Block`, `Eigen::Map`, other
Eigen expressions, vectors, and other containers require an explicit
`ArrayView2D`.

The existing view-based interface remains available:

```cpp
std::vector<double> G(n1 * n2, 1.0), D11(n1 * n2, 1.0);
std::vector<double> D12(n1 * n2, 0.0), D22(n1 * n2, 1.0);

Equation equation_view(make_view2d(G, n1, n2), make_view2d(D11, n1, n2),
                       make_view2d(D12, n1, n2), make_view2d(D22, n1, n2), bc);
equation_view.set_loss(0.0);          // optional; default is zero
```

`ArrayView2D` can also describe a strided slice:

```cpp
// Logical (n1,n2) slice of a row-major (n1,n2,2) array, fixed k=1.
auto view = make_view2d(data.data() + 1, n1, n2, n2 * 2, 2);
```

The optional `container_adapter.h` provides
`make_view2d(container)` for xtensor-like 2-D containers with `data()`,
`dimension()`, `shape(k)`, and `strides()`; the core does not include xtensor.

## Albert–Young case

The real Albert–Young diffusion tables are preprocessed once by
`tools/convert_albert_young.py` into cell-centered `G`, `D11`, `D12`, and `D22`
arrays. The C++ case does not interpolate or convert the raw coefficients.
The default converted data are in
`data/albert_young/ALBERT_YOUNG_80x80.h5`; see
`data/albert_young/README.md` for units, formulas, and the HDF5 format. The
legacy binary fixture remains for the regression harness.

`make test` runs the actual 80×80 Albert–Young case for 10 steps and compares
it with the legacy reference output. The observed agreement is
`max |delta f10| = 1.07e-14`. The Albert–Young driver uses `Eigen::MatrixXd`
coefficient and distribution arrays, passed directly to `Equation` and
`Solver::solve`.

## Example driver

[`src/main.cc`](src/main.cc) is a short, complete introduction to the API. It
builds a small nonuniform grid, uses `Eigen::MatrixXd` coefficients and fields,
sets mixed Dirichlet/ZeroFlux boundaries, constructs `Equation` and `Solver`,
advances five steps with `solve(equation, f_in, f_out, dt)`, and prints the
center, minimum, and maximum values.

```bash
make example
```

CMake target:

```bash
cmake -S . -B build-cmake
cmake --build build-cmake --target sayram2d_example
./build-cmake/sayram2d_example
```

## Ownership and lifetime

- `Mesh` owns its geometry and is immutable after construction. It may be
  shared read-only by multiple solver instances.
- `Solver` stores a reference to the mesh; the mesh must outlive the solver.
- `Equation` holds non-owning coefficient views and a `const`
  `BoundaryConditions&`. The caller must keep both the coefficient matrices and
  the `BoundaryConditions` object alive for as long as the `Equation` is used.
  Changes to coefficient values, boundary values, boundary types, or data spans
  are visible on subsequent solves. Construction from temporary coefficient
  matrices or a temporary `BoundaryConditions` is rejected at compile time.
- `ArrayView2D` never allocates or copies.
- `Solver::solve(eq, f, f, dt)` is supported for in-place solves; the old field
  is fully read before the output is written.

## Indexing and strides

- `ArrayView2D<T>::operator()(i, j)` addresses element
  `data()[i * stride1() + j * stride2()]`. The view stores extents and strides
  explicitly; it never infers layout from the container.
- The Albert–Young HDF5 input stores named datasets `/alpha0_edges`,
  `/alpha0_centers`, `/logE_edges`, `/logE_centers`, `/G`, `/D11`, `/D12`, and
  `/D22`. `G[i,j]` corresponds to `alpha0_centers[i]` and `logE_centers[j]`.
  The loader reads these datasets explicitly and fills Eigen column-major
  matrices by `(i,j)`.
- `Eigen::MatrixXd` is column-major and is mapped to
  `stride1 = 1`, `stride2 = rows` without copying.
- Higher-dimensional slices, `Eigen::Block`/`Eigen::Map`, and other containers
  must construct an `ArrayView2D` with explicit extents and strides.

## Standalone application

The standalone Albert–Young application is under `src/app/`. It uses the
redesigned library but keeps configuration and I/O outside it.

```bash
make app
./build/sayram2d/sayram2d_app [config.ini]
```

The default configuration path is `app.ini` in the current working directory.
A sample is provided at [`app.ini`](app.ini).

### Configuration

```ini
[run]
run_id = AlbertYoung
output_dir = output

[input]
data_file = data/albert_young/ALBERT_YOUNG_80x80.h5

[time]
dt = 0.002
nsteps = 10

[output]
n_frames = 11
```

- `run_id` is required and may contain letters, digits, `_`, and `-`.
- `output_dir` defaults to `output`.
- `data_file` is the HDF5 mesh/coefficient file described in
  [`data/albert_young/README.md`](data/albert_young/README.md).
- `dt` is in days and must be positive.
- `nsteps` must be positive; it is used exactly as requested.
- `n_frames` includes the initial and final snapshots, defaults to
  `nsteps + 1`, and must satisfy `nsteps % (n_frames - 1) == 0`.
- Relative paths resolve against the process working directory.
- Unknown sections/keys, duplicates, malformed numbers, missing required keys,
  and invalid ranges are rejected with file/line diagnostics.
- The Albert–Young initial distribution and boundary functions are built in;
  boundary conditions are not configured through the INI in this phase.

### Output

The app writes one HDF5 file:

```text
<output_dir>/<run_id>/<run_id>_data.h5
  /alpha0     pitch-angle cell centers, degrees
  /logEN      log(E/MeV) - log(0.511875)
  /t          saved times in days, including 0 and nsteps*dt
  /f/0        initial distribution
  /f/1        first saved distribution
  ...
```

Fields have shape `(n1, n2)` with `(i, j)` indexing. Eigen matrices are packed
row-major into an explicit buffer before HDF5 writes; `Eigen::MatrixXd::data()`
is never written directly. HDF5 is an application dependency only; `make lib`
and a CMake build with `-DSAYRAM2D_BUILD_APP=OFF` do not link it.

CMake:

```bash
cmake -S . -B build-cmake
cmake --build build-cmake
ctest --test-dir build-cmake --output-on-failure
./build-cmake/sayram2d_app app.ini
```

The Makefile detects Homebrew HDF5 with `brew --prefix hdf5`; override
`HDF5_INCLUDE`, `HDF5_LIBDIR`, or `HDF5_LIBS` for other installations.

## Quick case workflow

A minimal case Makefile is provided separately from the main build/tests
Makefile:

```bash
make -f Makefile.case          # build and run the 1-day Albert-Young case
make -f Makefile.case plot     # run and write figures/*.pdf and *.png
make -f Makefile.case help
```

It uses [`config/albert_young_1day.ini`](config/albert_young_1day.ini):
`dt = 0.002` days, `nsteps = 500`, and 11 output frames at 0, 0.1, ..., 1.0
day. The resulting HDF5 file is
`output/AlbertYoung1day/AlbertYoung1day_data.h5`.

The plot compares the app result with the committed legacy 1-day reference
[`tests/albert_young/legacy_reference_1day.h5`](tests/albert_young/legacy_reference_1day.h5)
at 0.5 and 2 MeV for 0.1 and 1 day. Outputs are
[`figures/albert_young_profiles_log.png`](figures/albert_young_profiles_log.png) and
[`figures/albert_young_difference.png`](figures/albert_young_difference.png),
with PDF versions alongside them. The profile plot uses a logarithmic
vertical axis for `f`.

There are two `main.cc` files by design:

- [`src/main.cc`](src/main.cc) is the small, dependency-free library example
  built as `sayram2d_example`; it demonstrates the API with a 4×2 grid.
- [`src/app/main.cc`](src/app/main.cc) is the full standalone application that
  reads INI/HDF5 input and writes HDF5 output.

The library example stays tiny and independent of HDF5; the app contains the
input/configuration/output workflow.

## I/O and configuration

The core library (`Mesh`, `Equation`, `Solver`, `BoundaryConditions`) has no
file-format or configuration dependency. I/O and configuration live in the
application and test tools:

- `src/app/config.cc` implements the strict application INI parser;
- `src/app/io.cc` reads the HDF5 mesh/coefficient file and writes HDF5 output
  using the native C API;
- `tools/convert_albert_young.py` reads the raw HDF5 tables and writes
  `data/albert_young/ALBERT_YOUNG_80x80.h5`; it can also write the legacy
  `S2DAYD1` regression binary;
- `tests/albert_young_case.h` validates the legacy binary used by the
  regression harness; the application loader validates the HDF5 input;
- `tests/albert_young_reference.cc` writes the comparison CSV with explicit
  metadata: `i` is the q1/alpha index, `j` is the q2/log-E index, and the CSV
  row order is `i*n2+j`;
- `tests/compare_legacy.py`, `tests/compare_albert_young.py`, and
  `tests/compare_app_hdf5.py` are the comparison consumers used by the test
  targets.

The retired `p.ini` workflow exists only in the frozen legacy source under
`tests/legacy/legacy_source/`; the standalone app uses the simplified `app.ini`
schema above.

## Concurrency

A `Solver` owns its sparse matrix, factorization, and numerical workspaces and
is not safe for concurrent `solve()` calls on the same instance. Use one
`Solver` per concurrent solve. There is no restriction on sharing an immutable
`Mesh`, read-only coefficient data, or separate `BoundaryConditions` objects
across solver instances.

## Boundary conditions

- `ZeroFlux` is the default on all four sides.
- `Dirichlet` values are specified at **boundary vertices**, not cell centers:
  - `Q1Min`/`Q1Max` require `n2 + 1` values or a scalar constant;
  - `Q2Min`/`Q2Max` require `n1 + 1` values or a scalar constant.
- Two Dirichlet sides meeting at a corner must agree; the solver checks this
  on every solve and throws on conflict. Dirichlet + ZeroFlux uses the
  Dirichlet value; ZeroFlux + ZeroFlux retains the along-boundary
  interpolation used by the discretization.
- Nonzero Neumann flux is not implemented; the only non-Dirichlet condition is
  homogeneous zero flux.

## Numerical method

The solver preserves the documented NTPFA formulation:

- one linear solve per time step; no nonlinear iteration;
- each interior face is processed once, with the paired entries
  `M_KK += A_K`, `M_KL -= A_L`, `M_LL += A_L`, `M_LK -= A_K`;
- positive inverse loss time is treated implicitly; negative inverse loss time
  is treated as an explicit source on the right-hand side, as documented;
- the matrix is generally nonsymmetric and is solved with `Eigen::SparseLU`.

`Lambda = G * D`, cell measures, face orientation, vertex interpolation, and
boundary-face flux formulas match the pre-redesign implementation. The
redesigned solver agrees with the saved legacy outputs for zero and positive
inverse loss to roundoff (maximum observed difference about `2e-15`).

## Repository layout

```
src/include/        public headers
src/                implementation and main.cc example
src/app/            standalone Albert-Young application
config/             example case configurations
scripts/            plotting utilities
figures/            generated comparison figures
tests/              unit tests, legacy/Albert-Young regression, benchmarks
  legacy/           frozen pre-redesign source and reference outputs
  albert_young/     legacy and new Albert-Young reference outputs
data/raw/           raw Albert-Young HDF5 input
data/albert_young/  preprocessed cell-centered G/D data
tools/              data-conversion tools
docs/               numerical algorithm and design notes
```

The legacy command-line driver and xtensor dependency are not part of the new
library. The frozen legacy source is retained under
`tests/legacy/legacy_source/` only so the baseline harness remains
reproducible.

## License and citation

This code is released under the MIT License.

If you use Sayram-2D in your research, please cite:

> Peng, P., Tao, X., Peng, Z., Jiang, Y., Gao, Z., Yang, D., et al. (2024).
> Modeling radiation belt dynamics using a positivity-preserving finite volume
> method on general meshes. *Journal of Geophysical Research: Space Physics*,
> 129, e2024JA032919. https://doi.org/10.1029/2024JA032919
