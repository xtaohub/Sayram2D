# Albert–Young example

`main.cc` demonstrates the redesigned Sayram2D API with the real Albert–Young
80×80 setup:

- constructs `Mesh` from the q1/q2 edge arrays (the former `Grid2D` was merged
  into `Mesh` in Phase 1);
- loads already-converted cell-centered `G`, `D11`, `D12`, `D22` from
  `data/albert_young/ALBERT_YOUNG_80x80.bin`;
- configures the Albert–Young Dirichlet/ZeroFlux boundaries and initial field;
- constructs `Equation` and `Solver`;
- advances 10 steps with `solver.solve(equation, f_in, f_out, dt)`;
- updates `D11`/`D22` in place at step 5, without rebuilding the solver;
- writes the final field to a CSV file.

The example does not interpolate or convert the diffusion coefficients; that
is done once by `tools/convert_albert_young.py`.

## Build and run

Using the Makefile:

```bash
make example
```

The executable is `build/sayram2d/example`, and `make example` runs it with
the default data and writes `build/albert_young_example.csv`.

Run manually with custom paths:

```bash
./build/sayram2d/example \
  data/albert_young/ALBERT_YOUNG_80x80.bin \
  build/albert_young_example.csv
```

Using CMake:

```bash
cmake -S . -B build-cmake
cmake --build build-cmake --target sayram2d_albert_young
./build-cmake/sayram2d_albert_young \
  data/albert_young/ALBERT_YOUNG_80x80.bin \
  build/albert_young_example.csv
```

The expected output starts with:

```
step 5: updated D11 and D22 in place
wrote build/albert_young_example.csv
  cells=6400 steps=10 dt=0.002 days
  f range=[1.33227e-16, 3.72745]
```

The CSV columns are `i,j,alpha0_deg,E_MeV,f`.
