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
| `sayram2d::Array2D<T>` | Owning, contiguous, column-major 2-D array. Provides indexing, sizing, fill, iteration, and `view()`. |
| `sayram2d::ArrayView2D<T>` | Minimal non-owning strided 2-D view. Used by `Equation` and `Solver` and for slices of larger arrays. |
| `sayram2d::BoundaryConditions` | Per-side `Dirichlet` or `ZeroFlux` configuration. Dirichlet data are specified at boundary **vertices**. |
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

### Array2D

`Array2D<T>` owns contiguous, column-major 2D storage. It is the simple owning
type for coefficients and fields:

```cpp
Array2D<double> G(n1, n2, 1.0);
Array2D<double> D11(n1, n2, 1.0);
Array2D<double> D12(n1, n2, 0.8);
Array2D<double> D22(n1, n2, 1.0);

G(i, j) = ...;          // no layout decision at the call site
auto G_view = G.view(); // mutable view; a const Array2D gives a const view
```

`Array2D` provides `extent1()`, `extent2()`, `size()`, `operator()`, `fill()`,
and iteration. Its layout is always column-major:

```text
Array2D(i, j) == storage[i + j*n1]
```

This makes a fixed `q3` plane contiguous for a column-major `(q1,q2,q3)` array.
HDF5 files remain row-major at the serialization boundary; the I/O layer
converts explicitly.

### ArrayView2D

`ArrayView2D<T>` is the non-owning view used by `Equation` and `Solver`, and the
bridge between a 3D or 4D host array and the 2D solver. It does **not**
represent the whole 3D array. It represents one logical 2D plane with extents
`(n1, n2)` and two element strides.

Its default contiguous layout is **column-major**:

```text
stride1 = 1
stride2 = n1
```

`ArrayView2D` can also represent arbitrary strides, including row-major and
noncontiguous slices, but all contiguous views created by this library use the
column-major convention.

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

For owned data, use `Array2D` and call `view()`:

```cpp
Array2D<double> G(n1, n2, 1.0);
auto G_view = G.view();
```

For an external buffer, construct an explicit view with pointer, shape, and
strides:

```cpp
auto view = make_view2d(pointer, Shape2{n1, n2},
                        Strides2{stride1, stride2});
```

There is also a lower-level five-argument form:

```cpp
auto raw_view = make_view2d(pointer, n1, n2, stride1, stride2);
```

The three-argument pointer form creates a contiguous column-major view:

```cpp
auto contiguous_view = make_view2d(pointer, n1, n2);
```

### Equation

An `Equation` combines the four diffusion coefficients, an optional inverse
loss, and a boundary-condition object. It does not own any of them.

```cpp
Equation equation(G.view(), D11.view(), D12.view(), D22.view(), bc);

// Optional; inverse loss defaults to zero.
equation.set_loss(0.0);
```

The coefficients are `ArrayView2D<const double>` views. They are typically
obtained from `Array2D::view()` or from explicit slices of larger arrays.

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
- Dirichlet data are specified at boundary **vertices**:
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
solver.solve(equation, f_in.view(), f_out.view(), dt);
```

`f_in` and `f_out` are `ArrayView2D` views. The same storage may be used for
both arguments for an in-place solve. A `Solver` is reusable across sequential
solves on the same mesh, but it must not be used concurrently from multiple
threads.

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
- `Array2D` owns its storage; a view from `Array2D::view()` is valid only while
  that `Array2D` exists and has not been reallocated or reassigned.
- `ArrayView2D` never allocates or copies.
- In-place solves `solver.solve(eq, f.view(), f.view(), dt)` are supported.

### Concurrency

A `Solver` owns its sparse matrix, factorization, and numerical workspaces and
is not safe for concurrent `solve()` calls on the same instance. Use one
`Solver` per concurrent solve. Immutable meshes and read-only coefficient data
may be shared.

## 2. Examples

The examples are intentionally separate from the core library. Examples 1–2
are small library demos; the Albert–Young example is the complete standalone
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

### Example 1: time-independent PPFV

[`src/examples/example1.cc`](src/examples/example1.cc) solves the uniform-grid
PPFV test problem from the paper:

- 50×50 grid on `[0,1]²`;
- `G = D11 = D22 = 1`, `D12 = 0.8`;
- initial `f = exp(-y/0.02) sin(πx/2)`;
- `f=0` at `x=0`, zero flux at `x=1`, the prescribed lower-energy profile at
  `y=0`, and `f=0` at `y=1`.

It advances to `t=0.01` with `dt=1e-4` using `Array2D` coefficient/field storage
and `ArrayView2D` views inside `Equation` and `Solver`.

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

### Example 2: time-dependent coefficients and boundary conditions

[`src/examples/example2.cc`](src/examples/example2.cc) shows how to change the
diffusion coefficients and a Dirichlet boundary array during a simulation.

- The 50×50 grid advances for 100 steps with `dt = 0.01`.
- `D11`, `D12`, and `D22` decay as `exp(-t/tau)` with `tau = 10`.
- The `Q1Min` Dirichlet boundary is a 1D array over the `n2 + 1` boundary
  vertices, with a fixed spatial profile multiplied by the same decay factor.
- At `t = 1`, the decay factor is still about `0.905`, so neither the
  coefficients nor the boundary value collapse to zero.

The `Equation` is constructed once. Its coefficient views remain valid while the
underlying `Array2D` objects are updated in place, and its boundary reference
sees changes to the `q1_min_values` vector passed to `BoundaryConditions`. The
example prints `D11`, the midpoint Dirichlet value, and the field
minimum/maximum at `t = 0`, `0.5`, and `1`.

```bash
make example2
```

It builds and runs the example. The executable is also available at
`build/sayram2d/example2`.

With CMake:

```bash
cmake --build build --target sayram2d_example2
./build/sayram2d_example2
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
/q1_edges        (n1+1)  q1 = alpha0 in radians
/q2_edges        (n2+1)  q2 = log(E/MeV)
/G               (n1,n2)  cell-centered
/D11             (n1,n2)  cell-centered
/D12             (n1,n2)  cell-centered
/D22             (n1,n2)  cell-centered
```

`G`, `D11`, `D12`, and `D22` are defined at cell centers and have shape
`(n1, n2)`. They are not vertex-centered arrays and must not have shape
`(n1+1, n2+1)`. Cell-center coordinates are derived from the edge arrays by the
`Mesh` constructor.

HDF5 remains row-major at the file boundary: the reader converts these arrays
into column-major `Array2D` storage. Root attributes document `format_version`,
row-major `i*n2+j` ordering, and the coordinate units. The provided file is
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
Fields are packed row-major into an explicit buffer at the HDF5 boundary,
even though the in-memory `Array2D`/`ArrayView2D` layout is column-major. The
output root attribute records `index_order = "row-major i*n2+j"`.

#### Build, run, and plot

```bash
make
./sayram2d app.ini
python3 scripts/plot_albert_young.py
```

[`app.ini`](app.ini) uses `dt = 0.002` days, `nsteps = 500`, and 11 frames at
0, 0.1, ..., 1.0 day. The single result is
`output/AlbertYoung/AlbertYoung_data.h5`, containing both time slices.

When a legacy 1-day reference is available, the plot script compares results at
0.5 and 2 MeV for 0.1 and 1 day:

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

## 3. Using Sayram2D in 3D simulations

We demonstrate how Sayram2D can be used in 3D simulation codes. To make good
use of the column-major `ArrayView2D` and mirror the `Array2D` interface, define
a thin `Array3D` in the 3D host code. `Array3D` belongs to `sayram`, not to
`sayram2d`:

```cpp
// In the 3D code; not part of sayram2d.
template <typename T>
class Array3D {
 public:
  Array3D(std::size_t n1, std::size_t n2, std::size_t n3,
          T value = T{})
      : n1_(n1), n2_(n2), n3_(n3),
        values_(n1 * n2 * n3, value) {}

  T& operator()(std::size_t i, std::size_t j, std::size_t k) noexcept {
    return values_[i + n1_ * (j + n2_ * k)];
  }

  sayram2d::ArrayView2D<T> slice2d(std::size_t k) & noexcept {
    return {values_.data() + k * n1_ * n2_, n1_, n2_};
  }

  sayram2d::ArrayView2D<const T> slice2d(std::size_t k) const & noexcept {
    return {values_.data() + k * n1_ * n2_, n1_, n2_};
  }

  sayram2d::ArrayView2D<T> slice2d(std::size_t) && = delete;
  sayram2d::ArrayView2D<const T> slice2d(std::size_t) const && = delete;

 private:
  std::size_t n1_;
  std::size_t n2_;
  std::size_t n3_;
  std::vector<T> values_;
};
```

The 3D driver then becomes:

```cpp
Array3D<double> G(n1, n2, nL);
Array3D<double> D11(n1, n2, nL);
Array3D<double> D12(n1, n2, nL);
Array3D<double> D22(n1, n2, nL);
Array3D<double> f(n1, n2, nL);
Array3D<double> f_next(n1, n2, nL);

Mesh mesh(q1_edges, q2_edges);   // shared q1/q2 grid
BoundaryConditions bc;           // caller-owned

#pragma omp parallel
{
  Solver solver(mesh);           // one solver per worker

  #pragma omp for schedule(dynamic)
  for (std::size_t k = 0; k < nL; ++k) {
    Equation equation(G.slice2d(k), D11.slice2d(k),
                      D12.slice2d(k), D22.slice2d(k), bc);
    solver.solve(equation, f.slice2d(k), f_next.slice2d(k), dt);
  }
}
```

Each `slice2d` is O(1), returns a contiguous column-major view, and copies
nothing. `Equation` and `Solver` only see `ArrayView2D` objects.

Each worker needs its own `Solver`; coefficients may be shared read-only, and
each worker writes a distinct output plane. If all layers use identical
coefficients, the `Equation` can be constructed once outside the loop. Build
OpenMP code with the appropriate compiler and linker flags (for example
`-fopenmp` on GCC/Clang).

## 4. Python bindings

The optional Python extension exposes `Mesh`, `BoundaryConditions`, `Equation`,
`Solver`, and `Side`. It links the C++ library and requires pybind11.

Build it with:

```bash
python3 -m pip install pybind11
cmake -S . -B build-python \
  -DSAYRAM2D_BUILD_PYTHON=ON \
  -DSAYRAM2D_BUILD_APP=OFF \
  -DBUILD_TESTING=OFF \
  -Dpybind11_DIR="$(python3 -m pybind11 --cmakedir)"
cmake --build build-python
```

Run the time-dependent example:

```bash
PYTHONPATH=build-python/python \
  python3 python/example_time_dependent.py
```

The binding requires 2D NumPy arrays in Fortran order:

```python
G = np.asfortranarray(G)
assert G.flags.f_contiguous
```

Coefficient arrays are held by the `Equation` wrapper. Update them in place:

```python
D11[:, :] = decay
D12[:, :] = 0.8 * decay
D22[:, :] = decay
```

Do not reassign `D11 = ...` after constructing the equation; rebuild the
equation if you need new buffers. Boundary arrays passed to
`BoundaryConditions.set_dirichlet` are copied, so they can be updated and passed
again each time step.

`python/example_time_dependent.py` demonstrates the full workflow, including
time-dependent diffusion coefficients and a time-dependent `Q1Min` Dirichlet
array.

## Repository layout

```text
Makefile            sayram2d app, library, and examples
src/*.cc            numerical library implementation
src/include/        public library headers
src/app/            generic application configuration and HDF5 I/O
src/examples/       example1.cc and example2.cc
src/examples/albert_young/
                    Albert–Young driver and case model
scripts/            data conversion and plotting
figures/            generated comparison figures
data/albert_young/  raw and preprocessed Albert-Young data
docs/               numerical algorithm and design notes
python/             optional pybind11 extension and Python example
```

## License and citation

This code is released under the MIT License.

If you use Sayram-2D in your research, please cite:

> Peng, P., Tao, X., Peng, Z., Jiang, Y., Gao, Z., Yang, D., et al. (2024).
> Modeling radiation belt dynamics using a positivity-preserving finite volume
> method on general meshes. *Journal of Geophysical Research: Space Physics*,
> 129, e2024JA032919. https://doi.org/10.1029/2024JA032919
