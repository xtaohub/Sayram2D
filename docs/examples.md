# 2. Examples

[Back to README](../README.md)

The examples are intentionally separate from the core library. Example 1,
its row- and column-major storage variants, and Example 2 are small library
demos; the Albert–Young example is the complete standalone application.

## Application support layer

[`src/app/`](../src/app/) contains generic application infrastructure:

- [`config.h`](../src/app/config.h) — strict INI parsing;
- [`io.h`](../src/app/io.h) — HDF5 input and output.

This code is not part of the numerical library and is not required by library
users. It is kept separate so that other realistic 2D applications can reuse
the configuration and I/O patterns. The methods
`HDF5Output::write_vector(values, name, units)` and
`HDF5Output::write_field(frame, field, name, units)` let callers choose dataset
names and units; the Albert–Young driver uses `alpha0`, `logE`, `t`, and
`PSD`.

## Example 1: a time-independent problem

[`src/examples/example1.cc`](../src/examples/example1.cc) solves the uniform-grid
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

## Example 1 variants: external storage

Two additional programs keep the Example 1 structure but replace `Array2D`
with flat `std::vector<double>` buffers and explicit `ArrayView2D` views:

- [`src/examples/example1_row.cc`](../src/examples/example1_row.cc) — row-major
  storage, `data[i*n2 + j]`, view strides `(n2, 1)`;
- [`src/examples/example1_col.cc`](../src/examples/example1_col.cc) —
  column-major storage, `data[i + j*n1]`, view strides `(1, n1)`.

Both construct the views before `Equation` and `Solver`, then run the same
time-independent problem and print the same minima and maxima as Example 1.

```bash
make example1_row
make example1_col
```

The executables are also available at `build/sayram2d/example1_row` and
`build/sayram2d/example1_col`.

## Example 2: time-dependent coefficients and boundary conditions

[`src/examples/example2.cc`](../src/examples/example2.cc) shows how to change the
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

## Albert–Young example

The standalone application is a worked Albert–Young example. It adds INI
configuration on top of the library, and HDF5 input/output for the coefficient
tables and saved distributions. It is not required to use the core library.

The INI file exists so that run parameters can change without recompiling:
time step, number of steps, saved-frame count, input data path, run ID, and
output directory. This lets one executable serve many runs, keeps case
definitions reproducible, and makes batch runs scriptable. The Albert–Young
initial distribution, boundary functions, and coefficient interpretation remain
compiled into the example.

The generic HDF5 I/O layer in [`src/app/io.h`](../src/app/io.h) is separate from
the solver and reads/writes named mesh, coefficient, and field arrays. The
writer takes configurable dataset names; by default the Albert–Young layout is
`alpha0`/`logE`/`t`/`PSD`. The INI parser is similarly reusable.

Application support types (`Config`, `InputData`, `HDF5Output`) are in
`sayram2d::app`. The Albert–Young model and initial/boundary functions are in
`sayram2d::albert_young`.

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

The compiler, Eigen, and HDF5 settings are the variables at the top of the
[`Makefile`](../Makefile). Edit them once for your machine if needed.

The standalone app requires the HDF5 C headers and library. Their paths are
the `HDF5_INCLUDE` and `HDF5_LIBDIR` variables at the top of the
[`Makefile`](../Makefile); edit them if your installation is elsewhere.

### Configuration

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

### Input

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

### Output

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

### Build, run, and plot

```bash
make
./sayram2d app.ini
python3 scripts/plot_albert_young.py
```

[`app.ini`](../app.ini) uses `dt = 0.002` days, `nsteps = 500`, and 11 frames at
0, 0.1, ..., 1.0 day. The single result is
`output/AlbertYoung/AlbertYoung_data.h5`, containing both time slices.

The plot script compares results with the bundled legacy PSD frames at 0.1 and
1 day at 0.5 and 2 MeV:

- [`figures/albert_young_profiles_log.png`](../figures/albert_young_profiles_log.png)
  — profiles with a logarithmic `f` axis;
- [`figures/albert_young_difference.png`](../figures/albert_young_difference.png)
  — signed differences.

PDF versions are alongside them.

The application entry point is [`src/examples/albert_young/main.cc`](../src/examples/albert_young/main.cc). The
dependency-free library examples are separate, under [`src/examples/`](../src/examples/).

### Performance

The following single-threaded timing was measured on an Apple M2 Max
(12-core CPU, 64 GB memory) running macOS 26.6.2, using the Makefile build
with Apple clang 17.0.0 and `-O2`:

- grid: 80 × 80 (6400 cells), `alpha0` in [5, 90] degrees, log-spaced `E` in
  [0.2, 5] MeV;
- time step: `0.002` days;
- 500 steps, for a total of 1.0 day, with 11 saved frames.

| Measurement | Time |
|---|---:|
| Full application, median of 5 runs | 3.84 s |
| Solver loop only, median of 5 runs | 3.70 s |
| Configuration, setup, and I/O difference | ~0.1 s |

The solver dominates this case.

### Data

The raw tables are interpolated once by `scripts/convert_albert_young.py` into
`data/albert_young/ALBERT_YOUNG_80x80.h5`. The app consumes that file directly
and does not interpolate coefficients. See
[`data/albert_young/README.md`](../data/albert_young/README.md) for units,
formulas, and the HDF5 layout.

### I/O and configuration separation

The core library (`Mesh`, `Equation`, `Solver`, `BoundaryConditions`) has no
file-format or configuration dependency. The application owns:

- `src/app/config.cc` — strict INI parser;
- `src/app/io.cc` — HDF5 reader/writer using the native C API;
- `src/examples/albert_young/albert_young.cc` — built-in initial and boundary functions.

The converter and regression tools remain outside the library.
