# Albert–Young diffusion data

The original data were provided by Jay M. Albert. The tables and formulas
below describe the converted, simulation-ready form used by Sayram2D.

The raw input is the Albert–Young chorus diffusion table:

```
data/albert_young/AlbertYoung_chorus.h5
  /alpha0  (91,)    degrees
  /E       (49,)    MeV
  /Daa     (91,49)  raw pitch-angle diffusion
  /Dap     (91,49)  raw cross diffusion
  /Dpp     (91,49)  raw energy/momentum diffusion
```

`scripts/convert_albert_young.py` interpolates these tables once onto a chosen
Sayram2D cell-center grid and writes the **cell-centered, simulation-ready**
arrays `G`, `D11`, `D12`, `D22` directly. The C++ cases and the standalone
application (`src/`) load those arrays; they do not interpolate or convert
units.

## Converted files

The standalone application reads the self-describing HDF5 file
`ALBERT_YOUNG_80x80.h5`. The `ALBERT_YOUNG_legacy_0.1d_1d.h5` file stores the
legacy PSD frames at 0.1 and 1 day used by the plotting comparison.

The converted HDF5 file was generated with the default settings:

```
alpha0 in [5, 90] degrees, n1 = 80
E in [0.2, 5] MeV, log-spaced, n2 = 80
```

Reproduce it with:

```bash
python3 scripts/convert_albert_young.py \
  --input data/albert_young/AlbertYoung_chorus.h5 \
  --output data/albert_young/ALBERT_YOUNG_80x80.h5 \
  --n1 80 --n2 80 --alpha-min 5 --alpha-max 90 --E-min 0.2 --E-max 5
```

### HDF5 input layout

```text
/q1_edges        float64[n1+1]   q1 = alpha0 in radians
/q2_edges        float64[n2+1]   q2 = log(E/MeV)
/G               float64[n1,n2]  cell-centered
/D11             float64[n1,n2]  cell-centered
/D12             float64[n1,n2]  cell-centered
/D22             float64[n1,n2]  cell-centered
```

Root attributes:

```text
format_version = 1
index_order = "row-major i*n2+j"
q1_units = "rad"
q2_units = "log(E/MeV)"
coefficient_location = "cell centers"
```

`G[i,j]`, `D11[i,j]`, etc. are defined at cell centers and have shape
`(n1, n2)`. They are not vertex-centered arrays and must not have shape
`(n1+1, n2+1)`. The HDF5 flat buffer uses row-major order for file
compatibility; the application reads the datasets explicitly by name.

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

The converted 80×80 case is validated against the legacy Albert–Young
reference:

- maximum `|delta f0| = 8.88e-16`;
- maximum `|delta f10| = 1.07e-14`;
- minimum `f10 = 1.33e-16` (positive).

The committed `ALBERT_YOUNG_legacy_0.1d_1d.h5` file contains the legacy PSD
frames at 0.1 and 1 day for the same 80×80 case advanced to 1 day
(`dt = 0.002`, `nsteps = 500`, 11 frames at 0, 0.1, ..., 1.0 day). It is used
by `python3 scripts/plot_albert_young.py` for the 0.5/2 MeV profile comparison.
It stores `alpha0`, `logEN`, the full time vector, and `/f/1` and `/f/10` in the
legacy row-major layout.
