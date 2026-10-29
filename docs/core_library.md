# 1. Core library

[Back to README](../README.md)

## Build the library

Requirements:

- C++20 compiler
- Eigen 3.4 or later (header-only)
- GNU Make
- `ar` only when building the static library with Make
- HDF5 C headers/libraries only for the standalone application

The compiler, Eigen, and HDF5 settings are the variables at the top of the
[`Makefile`](../Makefile). Edit them once for your machine if needed.

Build the standalone app (the default target):

```bash
make
```

Build the optional static library:

```bash
make lib
```

The library itself does not require HDF5. `make lib` uses `ar` to archive the
compiled objects into `lib/libsayram2d.a`; the examples link that archive. The
default `make` builds the app directly from source, so it does not require
`ar`.

`make clean` removes temporary Make build files. `make distclean` also removes
`lib/libsayram2d.a` and the `sayram2d` executable.

## Use the library

After `make lib`, an external program needs:

- the public headers in [`src/include/`](../src/include/) — add this directory to
  the include path;
- the static archive [`lib/libsayram2d.a`](../lib/libsayram2d.a) — link it;
- Eigen 3.4 or later — the public headers include Eigen.

The core library does not require HDF5; HDF5 is only needed for the
Albert–Young application. A minimal program is:

```cpp
#include <iostream>
#include <vector>

#include "array2d.h"
#include "boundary_conditions.h"
#include "equation.h"
#include "mesh.h"
#include "solver.h"

int main() {
  using namespace sayram2d;

  Mesh mesh(std::vector<double>{0.0, 0.5, 1.0},
            std::vector<double>{0.0, 0.5, 1.0});

  Array2D<double> G(2, 2, 1.0);
  Array2D<double> D11(2, 2, 1.0);
  Array2D<double> D12(2, 2, 0.0);
  Array2D<double> D22(2, 2, 1.0);
  Array2D<double> f(2, 2, 1.0);

  BoundaryConditions bc;
  bc.set_dirichlet(Side::Q1Min, 0.0);

  Equation equation(G.view(), D11.view(), D12.view(), D22.view(), bc);
  Solver solver(mesh);
  solver.solve(equation, f.view(), f.view(), 0.01);

  std::cout << "f(0,0) = " << f(0, 0) << "\n";
}
```

A minimal `Makefile` for this program, assuming Sayram2D is at
`$(SAYRAM2D)`:

```make
SAYRAM2D ?= /path/to/sayram2d
EIGEN_INCLUDE ?= /path/to/eigen

CXX ?= c++
CXXFLAGS ?= -std=c++20 -O2
CPPFLAGS := -I$(SAYRAM2D)/src/include -I$(EIGEN_INCLUDE)
LDFLAGS := -L$(SAYRAM2D)/lib
LDLIBS := -lsayram2d

minimal: minimal.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@
```

`EIGEN_INCLUDE` must be the directory that contains `Eigen/` (for example
`/usr/include/eigen3` or a local Eigen checkout). The archive can also be
passed directly instead of using `-L`/`-l`:

```make
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(SAYRAM2D)/lib/libsayram2d.a -o minimal
```

## Components

| Component | Responsibility |
|---|---|
| `sayram2d::Mesh` | Owns `q1`/`q2` edges, centers, and widths. Provides cell centers/measures, vertex positions, face geometry, neighbors, and boundary queries. |
| `sayram2d::Array2D<T>` | Owning, contiguous, column-major 2-D array. Provides indexing, sizing, fill, iteration, and `view()`. |
| `sayram2d::ArrayView2D<T>` | Minimal non-owning strided 2-D view. Used by `Equation` and `Solver` and for slices of larger arrays. |
| `sayram2d::BoundaryConditions` | Per-side `Dirichlet` or `ZeroFlux` configuration. Dirichlet data are specified at boundary vertices. |
| `sayram2d::Equation` | Non-owning views of `G`, `D11`, `D12`, `D22`, optional inverse loss, and a non-owning reference to the caller-owned boundary configuration. |
| `sayram2d::Solver` | Reusable sparse-matrix workspace and NTPFA assembly. One call to `solve` performs one documented linear solve. |

## Mesh

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

## Array2D

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

We chose column-major for `Array2D` because this makes a fixed `q3` plane
contiguous for a column-major `(q1,q2,q3)` array. HDF5 files remain row-major at
the serialization boundary; the I/O layer converts explicitly.

Examples of using `Array2D` are shown in
[`example1.cc`](../src/examples/example1.cc),
[`example2.cc`](../src/examples/example2.cc), and the
[Albert–Young application](../src/examples/albert_young/main.cc).

## ArrayView2D

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
noncontiguous slices.

**Important:** The view does not check whether its computed memory offsets
remain inside the underlying allocation. It does not know the allocation size,
and arbitrary strides are not bounded by `n1 * n2`; the caller must ensure that
the pointer and strides describe valid memory. This is why `ArrayView2D` is a
low-level type: for ordinary 2D applications, `Array2D` is the simplest owning
container.
Use `ArrayView2D` directly when you already have external storage, as the
[Example 1 row-major](../src/examples/example1_row.cc) and
[column-major](../src/examples/example1_col.cc) variants do, or design an
analogous `Array3D` container for 3D applications.

### Memory model

```cpp
view(i, j) == view.data()[i * view.stride1() + j * view.stride2()]
```

`Shape2` holds the two logical extents `n1` and `n2`; `Strides2` holds
`stride1` and `stride2`. The strides are measured in elements, not bytes. The
view owns no memory; it only stores the pointer, extents, and strides. `size()`
is `n1 * n2`, and `empty()` is true when either extent is zero.

A mutable `ArrayView2D<T>` can be converted to `ArrayView2D<const T>` without
copying. The reverse conversion is not allowed.

### Constructing views

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

This is equivalent to `Strides2{1, n1}`; it is not correct for row-major
storage.

### Constructing views from external data

The examples below use the five-argument form
`make_view2d(pointer, n1, n2, stride1, stride2)`, where `n1` and `n2` are
`std::size_t` and `stride1` and `stride2` are `std::ptrdiff_t`. If a stride is
computed from a `std::size_t` dimension, cast it to `std::ptrdiff_t` first.

In every case the view must map logical `(i, j)` to the host element for
`(q1=i, q2=j)`. The host index formula and the strides must agree.

#### 1D data

A flat 1D buffer as a 2D view, row-major (`data[i*n2 + j]`):

```cpp
std::vector<double> data(n1 * n2);
const std::ptrdiff_t s1 = static_cast<std::ptrdiff_t>(n2);
const std::ptrdiff_t s2 = 1;
auto row_view = make_view2d(data.data(), n1, n2, s1, s2);
```

The same buffer viewed column-major (`data[i + j*n1]`):

```cpp
std::vector<double> data(n1 * n2);
const std::ptrdiff_t s1 = 1;
const std::ptrdiff_t s2 = static_cast<std::ptrdiff_t>(n1);
auto col_view = make_view2d(data.data(), n1, n2, s1, s2);
```

Both views can be used, as long as they are applied consistently throughout
the simulation. The full row-major and column-major variants are implemented
in [example1_row.cc](../src/examples/example1_row.cc) and
[example1_col.cc](../src/examples/example1_col.cc).

#### 2D data

A C 2D array is already row-major and contiguous:

```cpp
constexpr std::size_t n1 = 3;
constexpr std::size_t n2 = 4;
double data[n1][n2] = {};

auto view = make_view2d(&data[0][0], n1, n2,
                        static_cast<std::ptrdiff_t>(n2), 1);
```

A `std::vector<std::vector<double>>` does not store its rows in one allocation
and cannot be represented by a single `ArrayView2D`.

Column-major 2D data is equally direct. For a contiguous Fortran array
`a(n1, n2)`, element `(i, j)` is at `a[i + j*n1]`:

```cpp
// Fortran: double precision a(n1, n2)
double* a = fortran_array;  // pointer to a(1, 1)
auto view = make_view2d(a, n1, n2, 1, n1);
```

#### 3D data

Let `a` be a contiguous 3D array with extents `n1`, `n2`, `n3`, and select
the plane `k = 3` (so `n3 > 3`). Both examples below use the same logical
extents `(n1, n2, n3)`; only the storage order differs. All indices are
zero-based, including those in the Fortran example.

A C-style row-major array stores element `(i, j, k)` at
`(i*n2 + j)*n3 + k`. The plane at `k = 3` therefore starts at `data + k` and
has strides `(n2*n3, n3)`:

```cpp
std::vector<double> a(n1 * n2 * n3);
double* data = a.data();
const std::size_t k = 3;

auto plane = make_view2d(data + k, n1, n2,
                         static_cast<std::ptrdiff_t>(n2 * n3),
                         static_cast<std::ptrdiff_t>(n3));
```

A Fortran array `a(n1, n2, n3)` is column-major and stores element
`(i, j, k)` at `i + j*n1 + k*n1*n2`. The plane at `k = 3` starts at
`data + k*n1*n2` and has strides `(1, n1)`:

```cpp
std::vector<double> a(n1 * n2 * n3);
double* data = a.data();
const std::size_t k = 3;

auto plane = make_view2d(data + k * n1 * n2, n1, n2, 1,
                         static_cast<std::ptrdiff_t>(n1));
```

In both cases `k = 3` selects the same logical plane.

#### 4D data

Let `a` be a contiguous 4D array with extents `n1`, `n2`, `n3`, `n4`, and
fix the last two indices at `k3 = 3` and `k4 = 2` (so `n3 > 3`, `n4 > 2`).
Both examples below use the same logical extents `(n1, n2, n3, n4)`; only the
storage order differs. All indices are zero-based, including those in the
Fortran example.

A C-style row-major array stores element `(i, j, k3, k4)` at
`((i*n2 + j)*n3 + k3)*n4 + k4`. The plane at `(3, 2)` therefore starts at
`data + 3*n4 + 2` and has strides `(n2*n3*n4, n3*n4)`:

```cpp
std::vector<double> a(n1 * n2 * n3 * n4);
double* data = a.data();
const std::size_t k3 = 3;
const std::size_t k4 = 2;

auto plane = make_view2d(data + k3*n4 + k4, n1, n2,
                         static_cast<std::ptrdiff_t>(n2 * n3 * n4),
                         static_cast<std::ptrdiff_t>(n3 * n4));
```

A Fortran array `a(n1, n2, n3, n4)` is column-major and stores element
`(i, j, k3, k4)` at `i + j*n1 + k3*n1*n2 + k4*n1*n2*n3`. The plane at
`(3, 2)` starts at `data + 3*n1*n2 + 2*n1*n2*n3` and has strides `(1, n1)`:

```cpp
std::vector<double> a(n1 * n2 * n3 * n4);
double* data = a.data();
const std::size_t k3 = 3;
const std::size_t k4 = 2;

auto plane = make_view2d(data + k3*n1*n2 + k4*n1*n2*n3, n1, n2, 1,
                         static_cast<std::ptrdiff_t>(n1));
```

In both cases the fixed indices select the same logical 2D plane with
extents `(n1, n2)`.

A wrong stride is not detected by `ArrayView2D`; it silently maps to the
wrong element, so the host index and the strides must be kept consistent.

## Equation

An `Equation` combines the four diffusion coefficients, an optional inverse
loss, and a boundary-condition object. It does not own any of them.

```cpp
Equation equation(G.view(), D11.view(), D12.view(), D22.view(), bc);

// Optional; inverse loss defaults to zero.
equation.set_loss(0.0);
```

The coefficients are `ArrayView2D<const double>` views. They are typically
obtained from `Array2D::view()` or from explicit slices of larger arrays. This
is intentional: coefficients and boundary values are time-dependent, may be
slices of a larger array, and should not be copied. Grouping them also keeps the
solver call simple: `solver.solve(equation, f_in, f_out, dt)`.

Coefficient values and boundary values may be changed in place between solves.
The caller must keep the coefficient buffers and `BoundaryConditions` object
alive for as long as the equation is used. An example is given in
[Example 2](examples.md#example-2-time-dependent-coefficients-and-boundary-conditions).

## Boundary conditions

`BoundaryConditions` is caller-owned and is referenced, not copied, by
`Equation`.

```cpp
BoundaryConditions bc;
std::vector<double> q2_max_values(mesh.n1() + 1, 0.5);

bc.set_dirichlet(Side::Q1Min, 0.5);               // scalar
bc.set_zero_flux(Side::Q1Max);
bc.set_zero_flux(Side::Q2Min);
bc.set_dirichlet(Side::Q2Max, q2_max_values);     // vertex array
bc.validate(mesh);  // optional explicit validation
```

- `ZeroFlux` is the default on all four sides.
- Negative Dirichlet values are rejected; zero is allowed.
- Dirichlet data are specified at boundary vertices:
  - `Q1Min`/`Q1Max` require `n2 + 1` values or a scalar constant;
  - `Q2Min`/`Q2Max` require `n1 + 1` values or a scalar constant.
- Two Dirichlet sides meeting at a corner must agree; the solver checks this
  on every solve. Dirichlet + ZeroFlux uses the Dirichlet value.
- Nonzero Neumann flux is not implemented.

## Solver

`Solver` owns the sparse matrix, factorization, and numerical workspaces. One
`solve` call performs one time step.

```cpp
Solver solver(mesh);
solver.solve(equation, f_in.view(), f_out.view(), dt);
```

`f_in` and `f_out` are `ArrayView2D` views. The same storage may be used for
both arguments for an in-place solve. `f_in` must be nonnegative; negative
values are rejected. A `Solver` is reusable across sequential solves on the
same mesh, but it must not be used concurrently from multiple threads.

The discretization, positivity-preserving flux construction, and matrix
assembly are described in
[`docs/sayram2d_algorithm.md`](sayram2d_algorithm.md).

## Physics and unit helpers

[`src/include/physics.h`](../src/include/physics.h) provides shared
constants and conversions:

- `physics::pi`, `physics::speed_of_light`, `physics::electron_rest_energy_mev`;
- degrees ↔ radians;
- eV/keV/MeV ↔ joules;
- kinetic energy ↔ momentum, with caller-supplied units and speed of light.

The momentum/energy functions take the rest energy and speed of light as
explicit parameters, so the same helpers work for other particles and unit
systems:

```cpp
#include "physics.h"

const double rest_energy_j = sayram2d::physics::ev_to_joule(0.511e6);
const double momentum =
    sayram2d::physics::momentum(
        sayram2d::physics::mev_to_joule(2.0),  // kinetic energy, 2 MeV
        rest_energy_j,
        sayram2d::physics::speed_of_light);
```

The arguments must use a consistent unit system. For SI units, use joules and
metres per second; for `MeV` and `MeV/c`, pass `c = 1`. The Albert–Young case
uses the latter convention.

[`src/include/radiation_belt.h`](../src/include/radiation_belt.h) provides optional
radiation-belt helpers for the `(alpha0, log E)` coordinates:

- `dloge_dp(log_e, rest_energy, speed_of_light)`;
- `G_alpha0_p(alpha0, p)`;
- `G_alpha0_loge(alpha0, log_e, rest_energy, speed_of_light)`.

These helpers are not used by `Equation` or `Solver`; callers still supply the
`G` field.

## Ownership and lifetime

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

## Concurrency

A `Solver` owns its sparse matrix, factorization, and numerical workspaces and
is not safe for concurrent `solve()` calls on the same instance. Use one
`Solver` per concurrent solve. Immutable meshes and read-only coefficient data
may be shared.

## Repository layout

```text
Makefile            sayram2d app, library, and examples
LICENSE             MIT license text and copyright notice
src/*.cc            numerical library implementation
src/include/        public library headers
src/app/            generic application configuration and HDF5 I/O
src/examples/       example1.cc, example1_row.cc, example1_col.cc, and example2.cc
src/examples/albert_young/
                    Albert–Young driver and case model
scripts/            data conversion and plotting
figures/            generated comparison figures
data/albert_young/  raw and preprocessed Albert-Young data
docs/               algorithm, library, and example documentation
lib/                generated static library (libsayram2d.a)
```
