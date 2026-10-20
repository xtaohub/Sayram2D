# Albert–Young diffusion data

The raw input is the Albert–Young chorus diffusion table:

```
data/raw/AlbertYoung_chorus.h5
  /alpha0  (91,)    degrees
  /E       (49,)    MeV
  /Daa     (91,49)  raw pitch-angle diffusion
  /Dap     (91,49)  raw cross diffusion
  /Dpp     (91,49)  raw energy/momentum diffusion
```

`tools/convert_albert_young.py` interpolates these tables once onto a chosen
Sayram2D cell-center grid and writes the **cell-centered, simulation-ready**
arrays `G`, `D11`, `D12`, `D22` directly. The C++ cases and the standalone
application (`src/app/`) load those arrays; they do not interpolate or convert
units.

## Converted files

The standalone application reads the self-describing HDF5 file
`ALBERT_YOUNG_80x80.h5`. The legacy `S2DAYD1` binary
`ALBERT_YOUNG_80x80.bin` is retained only for the existing regression harness.

Both were generated with the default conversion settings:

```
alpha0 in [5, 90] degrees, n1 = 80
E in [0.2, 5] MeV, log-spaced, n2 = 80
```

Reproduce both with:

```bash
python3 tools/convert_albert_young.py \
  --input data/raw/AlbertYoung_chorus.h5 \
  --output data/albert_young/ALBERT_YOUNG_80x80.h5 \
  --binary-output data/albert_young/ALBERT_YOUNG_80x80.bin \
  --n1 80 --n2 80 --alpha-min 5 --alpha-max 90 --E-min 0.2 --E-max 5
```

### HDF5 input layout

```text
/alpha0_edges    float64[n1+1]   alpha0 in radians
/alpha0_centers  float64[n1]     alpha0 in radians
/logE_edges      float64[n2+1]   log(E/MeV)
/logE_centers    float64[n2]     log(E/MeV)
/G               float64[n1,n2]  cell-centered
/D11             float64[n1,n2]
/D12             float64[n1,n2]
/D22             float64[n1,n2]
```

Root attributes:

```text
format_version = 1
index_order = "row-major i*n2+j"
alpha0_units = "rad"
logE_units = "log(E/MeV)"
```

`G[i,j]`, `D11[i,j]`, etc. correspond to `alpha0_centers[i]` and
`logE_centers[j]`. The application reads the HDF5 datasets explicitly by name;
the on-disk order is documented rather than inferred from Eigen layout.

### Legacy regression binary

The S2DAYD1 binary has the same arrays in a positional layout:

```text
magic[8] = b"S2DAYD1"
int32 version = 1
int32 n1, int32 n2
float64 q1_edges[n1+1]
float64 q2_edges[n2+1]
float64 q1_centers[n1]
float64 q2_centers[n2]
float64 G[n1*n2]       row-major, index i*n2+j
float64 D11[n1*n2]
float64 D12[n1*n2]
float64 D22[n1*n2]
```

It is used by `tests/albert_young_case.h` so the numerical regression tests do
not require HDF5. The standalone app does not read this binary.

## Conversion formulas

With `E = exp(q2)` and `p = sqrt(E(E + 2 E0))`, `E0 = 0.511875 MeV`:

```
dlogE/dp = p / (E (E + E0))
denormalize = (E0 / c^2)^2 c^2       (c = 1)
D11 = Daa_raw * denormalize * 86400 / p^2
D12 = Dap_raw * denormalize * 86400 * (dlogE/dp) / p
D22 = Dpp_raw * denormalize * 86400 * (dlogE/dp)^2
G   = p^2 * (1.30 - 0.56 sin(alpha0)) * sin(alpha0) cos(alpha0)
             / (dlogE/dp)
```

The interpolation uses the same alpha/logE bilinear rule as the legacy
`Albert_Young::constructD`, with safe clamping at the raw-grid endpoints.
`D21 = D12`; the tensor is validated positive definite after conversion.

## Verification

`tests/albert_young_reference.cc` runs the real 80×80 case for 10 steps at
`dt = 0.002` days using the legacy binary fixture. The standalone app reads the
HDF5 file and writes a matching HDF5 result. Both results match the legacy
Albert–Young run saved in `tests/albert_young/legacy_reference.csv`:

- maximum `|delta f0| = 8.88e-16`;
- maximum `|delta f10| = 1.07e-14`;
- minimum `f10 = 1.33e-16` (positive).

Run it with:

```bash
make test
# or
make albert-young
```

The committed `tests/albert_young/legacy_reference_1day.h5` is the legacy
standalone application output for the same 80×80 case advanced to 1 day
(`dt = 0.002`, `nsteps = 500`, 11 frames at 0, 0.1, ..., 1.0 day). It is used
by `make -f Makefile.case plot` for the 0.5/2 MeV profile comparison.
