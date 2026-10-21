# Sayram2D

Sayram2D is a C++20 code for two-dimensional Fokker–Planck diffusion on
structured, orthogonal, possibly nonuniform rectangular grids. It implements
the positivity-preserving finite-volume NTPFA discretization described in
[`docs/sayram2d_algorithm.md`](docs/sayram2d_algorithm.md).

The project is organized around a reusable core library and a set of examples.
The core library depends only on Eigen. The standalone Albert–Young
application is one example; HDF5 is used by that example, not by the library.

## 1. Core library

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
compiled objects into `libsayram2d.a`; the examples link that archive. The
default `make` builds the app directly from source, so it does not require
`ar`.


### Components

| Component | Responsibility |
|---|---|
| `sayram2d::Mesh` | Owns `q1`/`q2` edges, centers, and widths. Provides cell centers/measures, vertex positions, face geometry, neighbors, and boundary queries. |
| `sayram2d::ArrayView2D<T>` | Minimal non-owning strided 2-D view. Supports row-major, column-major, and noncontiguous slices of larger arrays. |
| `sayram2d::BoundaryConditions` | Per-side `Dirichlet` or `ZeroFlux` configuration. Dirichlet data are specified at boundary vertices. |
| `sayram2d::Equation` | Non-owning views of `G`, `D11`, `D12`, `D22`, optional inverse loss, and a non-owning reference to the caller-owned boundary configuration. |
| `sayram2d::Solver` | Reusable sparse-matrix workspace and NTPFA assembly. One call to `solve` performs one documented linear solve. |

### Mesh

`Mesh` owns the two edge arrays and is immutable after construction. It checks
that edges are finite, strictly increasing, and contain at least one cell.

```cpp
std::vector<double> q1_edges{0.0, 0.3, 0.7, 1.0};
std::vector<double> q2_edges{0.0, 0.4, 1.0, 1.7};
Mesh mesh(std::move(q1_edges), std::move(q2_edges));

const std::size_t n1 = mesh.n1();
const std::size_t n2 = mesh.n2();
```

`n1` and `n2` are the numbers of cells in the two directions. Cell centers are
available through `q1_center(i)` and `q2_center(j)`. A mesh may be shared
read-only by multiple solvers.

### ArrayView2D

`ArrayView2D<T>` is the bridge between a 3D or 4D host array and the 2D solver.
It does **not** represent the whole 3D array. It represents one logical 2D plane
with extents `(n1, n2)` and two element strides.

#### Memory model

```cpp
view(i, j) == view.data()[i * view.stride1() + j * view.stride2()]
```

`Shape2` holds the two logical extents `n1` and `n2`; `Strides2` holds
`stride1` and `stride2`. The strides are measured in elements, not bytes. The
view owns no memory; it only stores the pointer, extents, and strides. `size()`
is `n1 * n2`, and `empty()` is true when either extent is zero.

A mutable `ArrayView2D<T>` can be converted to `ArrayView2D<const T>` without
copying. The reverse conversion is not allowed.

#### Constructing views

For a contiguous `std::vector`, the default is row-major:

```cpp
std::vector<double> G(n1 * n2, 1.0);  // n1*n2 values
const Shape2 shape{n1, n2};
auto G2 = make_view2d(G, shape);      // G2(i, j) == G[i*n2 + j]
```

A `std::vector` has no intrinsic row-major or column-major layout; it is a
linear buffer. Its length must equal `n1 * n2`; the caller chooses both the
logical shape and the index mapping. `ArrayLayout` selects that mapping:

- `ArrayLayout::row_major`: `G[i*n2 + j]`
- `ArrayLayout::column_major`: `G[i + j*n1]`

Fill the buffer in the same order you intend to read it. For row-major:

```cpp
for (std::size_t i = 0; i < n1; ++i) {
  for (std::size_t j = 0; j < n2; ++j) {
    G[i * n2 + j] = static_cast<double>(i + j);
  }
}
auto row_view = make_view2d(G, shape);
```

Calling `make_view2d(G, shape, ArrayLayout::column_major)` on the same
row-major-filled buffer is legal, but it maps `(i,j)` to `G[i + j*n1]`. It does
not transpose or rearrange the data; it simply applies a different indexing
convention. Use `column_major` only when the buffer was deliberately filled as
`G[i + j*n1]`.

Do not treat a 1D profile as a 2D field unless `n×1` or `1×n` is genuinely
what you mean.

For arbitrary layouts, use the shape/stride form:

```cpp
auto view = make_view2d(pointer, Shape2{n1, n2},
                        Strides2{stride1, stride2});
```

There is also a lower-level five-argument form:

```cpp
auto raw_view = make_view2d(pointer, n1, n2, stride1, stride2);
```

`Eigen::MatrixXd` has a direct convenience path: a full `MatrixXd` lvalue is
accepted by `Equation` and `Solver` without an explicit view, using Eigen's
column-major layout:

```cpp
Eigen::MatrixXd G(n1, n2);
Eigen::MatrixXd D11(n1, n2);
Eigen::MatrixXd D12(n1, n2);
Eigen::MatrixXd D22(n1, n2);

Equation equation(G, D11, D12, D22, bc);  // direct, zero-copy
Solver solver(mesh);
solver.solve(equation, f_in, f_out, dt);
```

`Eigen::Block`, `Eigen::Map`, and slices of larger arrays are not implicitly
accepted; construct an explicit `ArrayView2D` for them.

#### Common 3D layouts

For a 3D buffer `data[i][j][k]`, the correct offset and strides depend on the
storage order:

| Storage order | Fixed-`k` offset | `stride1` | `stride2` |
|---|---:|---:|---:|
| `k` fastest: `data[(i*n2 + j)*n3 + k]` | `k` | `n2*n3` | `n3` |
| `k` slowest: `data[(k*n1 + i)*n2 + j]` | `k*n1*n2` | `n2` | `1` |

The first layout is common for arrays that store several fields or energy
channels together. The second layout makes each 2D plane contiguous.

For a 4D array with the last two indices fixed, the same rule applies. If

```text
data[((i*n2 + j)*n3 + k)*n4 + l]
```

then fixing `k` and `l` gives:

```cpp
const std::ptrdiff_t offset = k*n4 + l;
const std::ptrdiff_t stride1 = n2*n3*n4;
const std::ptrdiff_t stride2 = n3*n4;
```

In general, remove the two fixed axes from the shape, keep the strides of the
two remaining axes, and add the fixed-index contributions to the data pointer.

#### Complete 3D example

The following uses the `k`-fastest layout:

```cpp
const std::size_t n3 = 8;
std::vector<double> G(n1 * n2 * n3, 1.0);
std::vector<double> D11(n1 * n2 * n3, 1.0);
std::vector<double> D12(n1 * n2 * n3, 0.0);
std::vector<double> D22(n1 * n2 * n3, 1.0);
std::vector<double> f(n1 * n2 * n3);
std::vector<double> f_next(n1 * n2 * n3);

const Shape2 shape{n1, n2};
const std::ptrdiff_t stride1 =
    static_cast<std::ptrdiff_t>(n2 * n3);
const std::ptrdiff_t stride2 = static_cast<std::ptrdiff_t>(n3);

auto layer = [&](std::vector<double>& values, std::size_t k) {
  return make_view2d(values.data() + k, shape,
                     Strides2{stride1, stride2});
};

Mesh mesh(q1_edges, q2_edges);   // shared q1/q2 grid
BoundaryConditions bc;           // caller-owned
Solver solver(mesh);

for (std::size_t k = 0; k < n3; ++k) {
  // Cheap: stores four views and the BC reference, checks extents, and does
  // not allocate a matrix or rebuild the sparse pattern.
  Equation equation(layer(G, k), layer(D11, k), layer(D12, k),
                    layer(D22, k), bc);
  solver.solve(equation, layer(f, k), layer(f_next, k), dt);
}
```

No 3D container is copied or interpreted by the library. The caller controls
the layout, and the library only sees four coefficient views and two field
views.

If all layers use identical coefficients, the `Equation` can be constructed
once outside the loop. `Solver` is reusable across sequential layers on the
same mesh. The solve dominates the cost; constructing `Equation` per layer is
nearly free.

#### Parallel 3D solves

The natural parallel unit is a 2D plane:

- each worker processes a different `k` (or `(k,l)` in 4D);
- each worker needs its own `Solver`; a `Solver` is not safe for concurrent
  `solve()` calls;
- coefficient planes may be shared read-only;
- each worker writes to a distinct output plane.

The caller is responsible for keeping the parent arrays and the
`BoundaryConditions` object alive for the duration of every solve.

### Equation

An `Equation` combines the four diffusion coefficients, an optional inverse
loss, and a boundary-condition object. It does not own any of them.

```cpp
Equation equation(G, D11, D12, D22, bc);

// Optional; inverse loss defaults to zero.
equation.set_loss(0.0);
```

The coefficients may be full `Eigen::MatrixXd` objects or explicit
`ArrayView2D` objects. A full `MatrixXd` lvalue is accepted directly and uses
Eigen's column-major layout. `Eigen::Block`, `Eigen::Map`, and slices of larger
arrays must be passed as `ArrayView2D`.

Coefficient values and boundary values may be changed in place between solves.
The caller must keep the coefficient buffers and `BoundaryConditions` object
alive for as long as the equation is used.

### Boundary conditions

`BoundaryConditions` is caller-owned and is referenced, not copied, by
`Equation`.

```cpp
BoundaryConditions bc;
bc.set_dirichlet(Side::Q1Min, 0.5);
bc.set_zero_flux(Side::Q1Max);
bc.set_zero_flux(Side::Q2Min);
bc.set_dirichlet(Side::Q2Max, 0.0);
bc.validate(mesh);  // optional explicit validation
```

- `ZeroFlux` is the default on all four sides.
- Dirichlet data are specified at boundary vertices:
  - `Q1Min`/`Q1Max` require `n2 + 1` values or a scalar constant;
  - `Q2Min`/`Q2Max` require `n1 + 1` values or a scalar constant.
- Two Dirichlet sides meeting at a corner must agree; the solver checks this
  on every solve. Dirichlet + ZeroFlux uses the Dirichlet value.
- Nonzero Neumann flux is not implemented.

### Solver

`Solver` owns the sparse matrix, factorization, and numerical workspaces. One
`solve` call performs one time step.

```cpp
Solver solver(mesh);
solver.solve(equation, f_in, f_out, dt);
```

`f_in` and `f_out` may be `Eigen::MatrixXd` objects or `ArrayView2D` views. The
same storage may be used for both arguments for an in-place solve. A `Solver`
is reusable across sequential solves on the same mesh, but it must not be used
concurrently from multiple threads.

The discretization, positivity-preserving flux construction, and matrix
assembly are described in
[`docs/sayram2d_algorithm.md`](docs/sayram2d_algorithm.md).

### Physics and unit helpers

[`src/include/physics.h`](src/include/physics.h) provides shared
constants and conversions:

- `physics::pi`, `physics::speed_of_light`, `physics::electron_rest_energy_mev`;
- degrees ↔ radians;
- eV/keV/MeV ↔ joules;
- kinetic energy ↔ momentum, in SI and in MeV/c.

The momentum/energy functions take the rest energy as an explicit parameter, so
the same helpers work for other particles:

```cpp
#include "physics.h"

const double rest_energy_j = sayram2d::physics::ev_to_joule(0.511e6);
const double momentum =
    sayram2d::physics::momentum_si(
        sayram2d::physics::mev_to_joule(2.0),  // kinetic energy, 2 MeV
        rest_energy_j);
```

`physics::momentum_mev_c(kinetic_mev)` is the MeV/c convenience used by the
Albert–Young model.

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

### Concurrency

A `Solver` owns its sparse matrix, factorization, and numerical workspaces and
is not safe for concurrent `solve()` calls on the same instance. Use one
`Solver` per concurrent solve. Immutable meshes and read-only coefficient data
may be shared.

## 2. Examples

The examples are intentionally separate from the core library. Examples 1 and
2 are small library demos; the Albert–Young example is the complete standalone
application.

### Application support layer

[`src/app/`](src/app/) contains generic application infrastructure:

- [`config.h`](src/app/config.h) — strict INI parsing;
- [`io.h`](src/app/io.h) — HDF5 input and output.

This code is not part of the numerical library and is not required by library
users. It is kept separate so that other realistic 2D applications can reuse
the configuration and I/O patterns. The methods
`HDF5Output::write_vector(values, name, units)` and
`HDF5Output::write_field(frame, field, name, units)` let callers choose dataset
names and units; the Albert–Young driver uses `alpha0`, `logE`, `t`, and
`PSD`.

### Example 1: vector-based PPFV

[`src/examples/example1.cc`](src/examples/example1.cc) solves the uniform-grid
PPFV test problem from the paper:

- 50×50 grid on `[0,1]²`;
- `G = D11 = D22 = 1`, `D12 = 0.8`;
- initial `f = exp(-y/0.02) sin(πx/2)`;
- `f=0` at `x=0`, zero flux at `x=1`, the prescribed lower-energy profile at
  `y=0`, and `f=0` at `y=1`.

It advances to `t=0.01` with `dt=1e-4` using `std::vector` coefficient/field
storage and explicit `make_view2d(...)` views.

```bash
make example1
```

It builds and runs the example, printing the final minimum and maximum values.
The executable is also available at `build/sayram2d/example1`.

With CMake:

```bash
cmake -S . -B build
cmake --build build --target sayram2d_example1
./build/sayram2d_example1
```

### Example 2: Eigen embedding

This is the simplest complete example of embedding Sayram2D in another
program. It uses the direct `Eigen::MatrixXd` interface; the view-based
interface is described under [ArrayView2D](#arrayview2d).

A complete, compilable version is
[`src/examples/example_eigen.cc`](src/examples/example_eigen.cc).

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

Build and run it with:

```bash
make example_eigen
```

### Albert–Young example

The standalone application is a worked Albert–Young example. It adds INI
configuration on top of the library, and HDF5 input/output for the coefficient
tables and saved distributions. It is not required to use the core library.

The INI file exists so that run parameters can change without recompiling:
time step, number of steps, saved-frame count, input data path, run ID, and
output directory. This lets one executable serve many runs, keeps case
definitions reproducible, and makes batch runs scriptable. The Albert–Young
initial distribution, boundary functions, and coefficient interpretation remain
compiled into the example.

The generic HDF5 I/O layer in [`src/app/io.h`](src/app/io.h) is separate from
the solver and reads/writes named mesh, coefficient, and field arrays. The
writer takes configurable dataset names; by default the Albert–Young layout is
`alpha0`/`logE`/`t`/`PSD`. The INI parser is similarly reusable.

Application support types (`Config`, `InputData`, `HDF5Output`) are in
`sayram2d::app`. The Albert–Young model and initial/boundary functions are in
`sayram2d::albert_young`.

#### Build and run

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

#### Configuration

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

#### Input

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

#### Output

The app writes one HDF5 file:

```text
<output_dir>/<run_id>/<run_id>_data.h5
  /alpha0     pitch-angle cell centers, units: deg
  /logE       log(E/MeV), units: log(MeV)
  /t          saved times, units: day, including 0 and nsteps*dt
  /PSD/0      initial distribution, units: arbitrary
  /PSD/1      first saved distribution, units: arbitrary
  ...
```

Each dataset carries a `units` attribute. Fields have shape `(n1, n2)`.
Eigen matrices are packed row-major into an explicit buffer before HDF5 writes.

#### Build, run, and plot

```bash
make
./sayram2d app.ini
python3 scripts/plot_albert_young.py
```

[`app.ini`](app.ini) uses `dt = 0.002` days, `nsteps = 500`, and 11 frames at
0, 0.1, ..., 1.0 day. The single result is
`output/AlbertYoung/AlbertYoung_data.h5`, containing both time slices.

The plot script compares this result with the committed legacy 1-day reference
[`tests/albert_young/legacy_reference_1day.h5`](tests/albert_young/legacy_reference_1day.h5)
at 0.5 and 2 MeV for 0.1 and 1 day:

- [`figures/albert_young_profiles_log.png`](figures/albert_young_profiles_log.png)
  — profiles with a logarithmic `f` axis;
- [`figures/albert_young_difference.png`](figures/albert_young_difference.png)
  — signed differences.

PDF versions are alongside them.

The application entry point is [`src/examples/albert_young/main.cc`](src/examples/albert_young/main.cc). The
dependency-free library examples are separate, under [`src/examples/`](src/examples/).

#### Data

The raw tables are interpolated once by `scripts/convert_albert_young.py` into
`data/albert_young/ALBERT_YOUNG_80x80.h5`. The app consumes that file directly
and does not interpolate coefficients. See
[`data/albert_young/README.md`](data/albert_young/README.md) for units,
formulas, and the HDF5 layout.

#### I/O and configuration separation

The core library (`Mesh`, `Equation`, `Solver`, `BoundaryConditions`) has no
file-format or configuration dependency. The application owns:

- `src/app/config.cc` — strict INI parser;
- `src/app/io.cc` — HDF5 reader/writer using the native C API;
- `src/examples/albert_young/albert_young.cc` — built-in initial and boundary functions.

The converter and regression tools remain outside the library.

## Tests

CMake/CTest is the test entry point:

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

The registered tests include unit tests, legacy numerical comparisons, the
Albert–Young regression, and standalone-app/HDF5 checks. The root `Makefile`
only builds the `sayram2d` app and the library.

## Repository layout

```text
Makefile            sayram2d app, library, and examples
src/*.cc            numerical library implementation
src/include/        public library headers
src/app/            generic application configuration and HDF5 I/O
src/examples/       example1.cc and example_eigen.cc
src/examples/albert_young/
                    Albert–Young driver and case model
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
