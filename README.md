# Sayram2D

Sayram2D is a C++20 code for two-dimensional Fokker–Planck diffusion on
structured, orthogonal, possibly nonuniform rectangular grids. It implements
the positivity-preserving finite-volume NTPFA discretization described in
[`docs/sayram2d_algorithm.tex`](docs/sayram2d_algorithm.tex).

The project provides two entry points:

- a standalone Albert–Young simulation application with INI configuration and
  HDF5 input/output;
- a reusable numerical library that can be embedded in other C++ callers.

The core library depends only on Eigen. HDF5 is used by the standalone
application, not by the library.

## Standalone application

### Build and run

```bash
make              # default: build ./sayram2d
# or: make app
./sayram2d [config.ini]
```

The command-line argument is an optional INI path. If omitted, the app uses
`app.ini` in the current working directory. There are no other command-line
options. Relative paths in the INI resolve against the working directory.
`make` writes the executable to `./sayram2d`, next to `app.ini`.

The compiler, Eigen, and HDF5 paths are the four variables at the top of the
[`Makefile`](Makefile). Edit them once for your machine if needed.

The standalone app requires the HDF5 C headers and library. Their paths are
the `HDF5_INCLUDE` and `HDF5_LIBDIR` variables at the top of the
[`Makefile`](Makefile); edit them if your installation is elsewhere.

With CMake:

```bash
cmake -S . -B build
cmake --build build
./sayram2d app.ini
```

The CMake target is `sayram2d_cli`; its output name is `sayram2d` and it is
placed at the repository root, next to `app.ini`.

### Configuration expected by the app

The app accepts this strict schema:

```ini
[run]
run_id = AlbertYoung
output_dir = output

[input]
data_file = data/albert_young/ALBERT_YOUNG_80x80.h5

[time]
dt = 0.002
nsteps = 500

[output]
n_frames = 11
```

| Section/key | Required | Meaning |
|---|---|---|
| `[run] run_id` | yes | Output subdirectory and file prefix. Letters, digits, `_`, and `-` only. |
| `[run] output_dir` | no | Parent output directory. Default: `output`. |
| `[input] data_file` | yes | HDF5 mesh and coefficient file. |
| `[time] dt` | yes | Time step in days; must be positive and finite. |
| `[time] nsteps` | yes | Positive number of steps; used exactly as requested. |
| `[output] n_frames` | no | Total saved frames, including initial and final. Default: `nsteps + 1`. Must satisfy `nsteps % (n_frames - 1) == 0`. |

The parser rejects missing required keys, unknown sections or keys, duplicate
keys, malformed numbers, non-finite values, nonpositive time settings, and
invalid frame counts. Diagnostics include the configuration file and line.

### Input expected by the app

`data_file` is a self-describing HDF5 file:

```text
/alpha0_edges    (n1+1)  alpha0 in radians
/alpha0_centers  (n1)    alpha0 in radians
/logE_edges      (n2+1)  log(E/MeV)
/logE_centers    (n2)    log(E/MeV)
/G               (n1,n2)  cell-centered
/D11             (n1,n2)
/D12             (n1,n2)
/D22             (n1,n2)
```

Root attributes document `format_version`, row-major `i*n2+j` ordering, and
the coordinate units. The provided file is
`data/albert_young/ALBERT_YOUNG_80x80.h5`, produced by
`scripts/convert_albert_young.py`.

The Albert–Young initial distribution and boundary functions are built into
the application; boundary conditions are not configured through the INI in
this version.

### Output written by the app

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

Fields have shape `(n1, n2)`. Eigen matrices are packed row-major into an
explicit buffer before HDF5 writes.

### Quick Albert–Young case

A separate minimal Makefile provides a one-command case workflow:

```bash
make -f Makefile.case          # build and run the 1-day case
make -f Makefile.case plot     # run and write figures/*.pdf and *.png
make -f Makefile.case help
```

It uses [`app.ini`](app.ini): `dt = 0.002` days, `nsteps = 500`, and
11 frames at 0, 0.1, ..., 1.0 day. The single result is
`output/AlbertYoung/AlbertYoung_data.h5`, containing both time slices.

The plot compares this result with the committed legacy 1-day reference
[`tests/albert_young/legacy_reference_1day.h5`](tests/albert_young/legacy_reference_1day.h5)
at 0.5 and 2 MeV for 0.1 and 1 day:

- [`figures/albert_young_profiles_log.png`](figures/albert_young_profiles_log.png)
  — profiles with a logarithmic `f` axis;
- [`figures/albert_young_difference.png`](figures/albert_young_difference.png)
  — signed differences.

PDF versions are alongside them.

The default entry point is [`src/main.cc`](src/main.cc), the standalone
application with configuration and HDF5 I/O. The dependency-free library
example is separate, at [`src/library/example.cc`](src/library/example.cc).

## Core library

### Components

| Component | Responsibility |
|---|---|
| `sayram2d::Mesh` | Owns `q1`/`q2` edges, centers, and widths. Provides cell centers/measures, vertex positions, face geometry, neighbors, and boundary queries. |
| `sayram2d::ArrayView2D<T>` | Minimal non-owning strided 2-D view. Supports row-major, column-major, and noncontiguous slices of larger arrays. |
| `sayram2d::BoundaryConditions` | Per-side `Dirichlet` or `ZeroFlux` configuration. Dirichlet data are specified at boundary vertices. |
| `sayram2d::Equation` | Non-owning views of `G`, `D11`, `D12`, `D22`, optional inverse loss, and a non-owning reference to the caller-owned boundary configuration. |
| `sayram2d::Solver` | Reusable sparse-matrix workspace and NTPFA assembly. One call to `solve` performs one documented linear solve. |

### Build the library

Requirements:

- C++20 compiler
- Eigen 3.4 or later (header-only)
- GNU Make
- `ar` only when building the static library
- HDF5 C headers/libraries only for the standalone application

The compiler, Eigen, and HDF5 paths are the four variables at the top of the
[`Makefile`](Makefile). Edit them once for your machine if needed.

Build the standalone app (the default target):

```bash
make
```

Build the optional static library:

```bash
make lib
```

The library itself does not require HDF5. `make lib` uses `ar` to archive the
compiled objects into `libsayram2d.a`; the example links that archive. The
default `make` builds the app directly from source, so it does not require
`ar`.

### Embedding example

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

BoundaryConditions bc;
bc.set_dirichlet(Side::Q1Min, 0.5);
bc.set_zero_flux(Side::Q1Max);
bc.set_zero_flux(Side::Q2Min);
bc.set_zero_flux(Side::Q2Max);

Equation equation(G, D11, D12, D22, bc);
// Optional; inverse loss defaults to zero.
equation.set_loss(0.0);

Solver solver(mesh);
Eigen::MatrixXd f_in = Eigen::MatrixXd::Ones(n1, n2);
Eigen::MatrixXd f_out(n1, n2);
solver.solve(equation, f_in, f_out, 0.01);
```

`Eigen::MatrixXd` is accepted directly, zero-copy, in column-major layout.
Only ordinary `MatrixXd` gets this overload; `Eigen::Block`, `Eigen::Map`,
other Eigen expressions, vectors, and other containers require an explicit
`ArrayView2D`.

The view-based interface remains available:

```cpp
std::vector<double> G(n1 * n2, 1.0), D11(n1 * n2, 1.0);
std::vector<double> D12(n1 * n2, 0.0), D22(n1 * n2, 1.0);

Equation equation_view(make_view2d(G, n1, n2), make_view2d(D11, n1, n2),
                       make_view2d(D12, n1, n2), make_view2d(D22, n1, n2), bc);
```

`ArrayView2D` can describe a strided slice:

```cpp
// Logical (n1,n2) slice of a row-major (n1,n2,2) array, fixed k=1.
auto view = make_view2d(data.data() + 1, n1, n2, n2 * 2, 2);
```

The optional `container_adapter.h` provides `make_view2d(container)` for
xtensor-like 2-D containers; the core does not include xtensor.

### Using the solver from a 3D code

If a 3D caller stores fields with the slice index fastest, for example
`volume[(i*n2 + j)*n3 + k]`, each 2D layer is an explicit strided view:

```cpp
const std::size_t n3 = 8;
std::vector<double> G(n1 * n2 * n3, 1.0);
std::vector<double> D11(n1 * n2 * n3, 1.0);
std::vector<double> D12(n1 * n2 * n3, 0.0);
std::vector<double> D22(n1 * n2 * n3, 1.0);
std::vector<double> f(n1 * n2 * n3);
std::vector<double> f_next(n1 * n2 * n3);

auto layer = [&](std::vector<double>& values, std::size_t k) {
  return make_view2d(values.data() + k, n1, n2,
                     static_cast<std::ptrdiff_t>(n2 * n3),
                     static_cast<std::ptrdiff_t>(n3));
};

Mesh mesh(q1_edges, q2_edges);   // shared q1/q2 grid
BoundaryConditions bc;           // shared or per-layer; caller-owned
Solver solver(mesh);

for (std::size_t k = 0; k < n3; ++k) {
  // Cheap: stores four views and the BC reference, checks extents, and does
  // not allocate a matrix or rebuild the sparse pattern.
  Equation equation(layer(G, k), layer(D11, k), layer(D12, k),
                    layer(D22, k), bc);
  solver.solve(equation, layer(f, k), layer(f_next, k), dt);
}
```

The view stores the explicit strides; no 3D container is copied or
interpreted by the library. `Solver` is reusable across sequential layers with
the same mesh. Constructing `Equation` inside the loop is nearly free: it only
stores non-owning views and a boundary reference, checks extents, and does no
matrix allocation or symbolic factorization. The solve dominates the cost. If
all layers use identical coefficients, the `Equation` can be hoisted out of the
loop. Use a separate `Solver` instance for concurrent solves.

### Example driver

[`src/library/example.cc`](src/library/example.cc) solves the uniform-grid
PPFV test problem from the paper:

- 50×50 grid on `[0,1]²`;
- `G = D11 = D22 = 1`, `D12 = 0.8`;
- initial `f = exp(-y/0.02) sin(πx/2)`;
- `f=0` at `x=0`, zero flux at `x=1`, the prescribed lower-energy profile at
  `y=0`, and `f=0` at `y=1`.

It advances to `t=0.01` with `dt=1e-4` using `std::vector` coefficient/field
storage and explicit `make_view2d(...)` views.

```bash
make example
```

CMake target: `sayram2d_example`.

### Ownership and lifetime

- `Mesh` owns its geometry and is immutable after construction. It may be
  shared read-only by multiple solver instances.
- `Solver` stores a reference to the mesh; the mesh must outlive the solver.
- `Equation` holds non-owning coefficient views and a `const`
  `BoundaryConditions&`. The caller must keep the coefficient matrices and
  the boundary object alive for as long as the equation is used.
- Coefficient values, boundary values, boundary types, and data spans may be
  changed in place between solves. Construction from temporary coefficient
  matrices or a temporary `BoundaryConditions` is rejected at compile time.
- `ArrayView2D` never allocates or copies.
- In-place solves `solver.solve(eq, f, f, dt)` are supported.

### Indexing and strides

- `ArrayView2D<T>::operator()(i, j)` addresses
  `data()[i * stride1() + j * stride2()]`.
- `Eigen::MatrixXd` is column-major and maps to `stride1 = 1`,
  `stride2 = rows` without copying.
- The application's HDF5 input uses named datasets and explicit `(i,j)`
  correspondence; file order is never inferred from Eigen layout.
- Higher-dimensional slices, `Eigen::Block`/`Eigen::Map`, and other containers
  must construct an `ArrayView2D` with explicit extents and strides.

### Boundary conditions

- `ZeroFlux` is the default on all four sides.
- Dirichlet data are specified at boundary vertices:
  - `Q1Min`/`Q1Max` require `n2 + 1` values or a scalar constant;
  - `Q2Min`/`Q2Max` require `n1 + 1` values or a scalar constant.
- Two Dirichlet sides meeting at a corner must agree; the solver checks this
  on every solve. Dirichlet + ZeroFlux uses the Dirichlet value.
- Nonzero Neumann flux is not implemented.

### Numerical method

- One linear solve per time step; no nonlinear iteration.
- Each interior face is processed once with paired matrix entries.
- Positive inverse loss is implicit; negative inverse loss is an explicit
  right-hand-side source.
- The matrix is generally nonsymmetric and is solved with `Eigen::SparseLU`.

The solver agrees with the saved legacy outputs for zero and positive inverse
loss to roundoff (maximum observed difference about `2e-15`).

### Concurrency

A `Solver` owns its sparse matrix, factorization, and numerical workspaces and
is not safe for concurrent `solve()` calls on the same instance. Use one
`Solver` per concurrent solve. Immutable meshes and read-only coefficient data
may be shared.

## Tests and benchmarks

```bash
make -f Makefile.tests
make -f Makefile.benchmark
```

CMake/CTest:

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

`make -f Makefile.tests` includes unit tests, legacy numerical comparisons,
the Albert–Young regression, and standalone-app/HDF5 checks. The root
`Makefile` only builds the library, example, and `sayram2d` app.

## Albert–Young data

The raw tables are interpolated once by `scripts/convert_albert_young.py` into
`data/albert_young/ALBERT_YOUNG_80x80.h5`. The app consumes that file directly
and does not interpolate coefficients. See
[`data/albert_young/README.md`](data/albert_young/README.md) for units,
formulas, and the HDF5 layout.

## I/O and configuration separation

The core library (`Mesh`, `Equation`, `Solver`, `BoundaryConditions`) has no
file-format or configuration dependency. The application owns:

- `src/config.cc` — strict INI parser;
- `src/io.cc` — HDF5 reader/writer using the native C API;
- `src/albert_young.cc` — built-in initial and boundary functions.

The converter and regression tools remain outside the library.

## Repository layout

```text
Makefile            library, example, and sayram2d app
Makefile.tests      unit/regression/application tests
Makefile.benchmark  timing benchmark
Makefile.case       one-command Albert-Young case
src/main.cc         standalone application entry point
src/config.*        application configuration
src/io.*            application HDF5 I/O
src/albert_young.*  application physics preset
src/library/        numerical library and its example
src/library/include public library headers
scripts/            data conversion and plotting
figures/            generated comparison figures
tests/              unit tests, legacy/Albert-Young regression, benchmarks
  legacy/           frozen pre-redesign source and reference outputs
  albert_young/     legacy and new Albert-Young reference outputs
data/albert_young/  raw and preprocessed Albert-Young data
docs/               numerical algorithm and design notes
```

## License and citation

This code is released under the MIT License.

If you use Sayram-2D in your research, please cite:

> Peng, P., Tao, X., Peng, Z., Jiang, Y., Gao, Z., Yang, D., et al. (2024).
> Modeling radiation belt dynamics using a positivity-preserving finite volume
> method on general meshes. *Journal of Geophysical Research: Space Physics*,
> 129, e2024JA032919. https://doi.org/10.1029/2024JA032919
