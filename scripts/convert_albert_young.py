#!/usr/bin/env python3
"""Preprocess the Albert-Young diffusion tables for Sayram2D.

The raw HDF5 file stores Daa, Dap, Dpp on an (alpha0, E) grid in the original
normalized units. This tool interpolates them to a chosen Sayram2D cell-center
grid and applies the coordinate/momentum/unit conversions once, so the C++
case code can consume G, D11, D12, D22 directly.

The primary output is a self-describing HDF5 file:

  /q1_edges        float64[n1+1]   alpha0 in radians
  /q2_edges        float64[n2+1]   log(E/MeV)
  /G               float64[n1,n2]  cell-centered, row-major, index i*n2+j
  /D11             float64[n1,n2]
  /D12             float64[n1,n2]
  /D22             float64[n1,n2]

The HDF5 attributes document the index order and coordinate units. The
standalone application reads this HDF5 file. An optional legacy binary can
also be written for the pre-existing regression harness:

  --binary-output data/albert_young/ALBERT_YOUNG_80x80.bin

q1 is alpha0 in radians; q2 is log(E), with E in MeV.

Example:
  python3 scripts/convert_albert_young.py \
    --input data/albert_young/AlbertYoung_chorus.h5 \
    --output data/albert_young/ALBERT_YOUNG_80x80.h5 \
    --binary-output data/albert_young/ALBERT_YOUNG_80x80.bin \
    --n1 80 --n2 80 --alpha-min 5 --alpha-max 90 --E-min 0.2 --E-max 5
"""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import h5py
import numpy as np

PI = 3.141592653589793238462643383279502884
E0 = 0.511875  # electron rest energy, MeV
C = 1.0
SECONDS_PER_DAY = 86400.0


def e2p(energy: np.ndarray | float) -> np.ndarray | float:
    return np.sqrt(energy * (energy + 2.0 * E0)) / C


def dloge_dp(log_e: np.ndarray | float) -> np.ndarray | float:
    energy = np.exp(log_e)
    p = e2p(energy)
    return p * C * C / (energy * (energy + E0))


def interpolate(raw: np.ndarray, x: np.ndarray, y: np.ndarray,
                q1_centers: np.ndarray, q2_centers: np.ndarray) -> np.ndarray:
    """Safe bilinear interpolation matching the legacy locate/interp2D rule."""
    nx, ny = raw.shape
    result = np.empty((q1_centers.size, q2_centers.size), dtype=np.float64)

    pos_x = (q1_centers - x[0]) / (x[-1] - x[0]) * (nx - 1)
    pos_y = (q2_centers - np.log(y[0])) / (np.log(y[-1]) - np.log(y[0])) * (ny - 1)

    for i, px in enumerate(pos_x):
        if px <= 0.0:
            i0, wx = 0, 0.0
        elif px >= nx - 1:
            i0, wx = nx - 2, 1.0
        else:
            i0 = int(np.floor(px))
            wx = px - i0
        for j, py in enumerate(pos_y):
            if py <= 0.0:
                j0, wy = 0, 0.0
            elif py >= ny - 1:
                j0, wy = ny - 2, 1.0
            else:
                j0 = int(np.floor(py))
                wy = py - j0
            result[i, j] = (
                raw[i0, j0] * (1.0 - wx) * (1.0 - wy)
                + raw[i0 + 1, j0] * wx * (1.0 - wy)
                + raw[i0, j0 + 1] * (1.0 - wx) * wy
                + raw[i0 + 1, j0 + 1] * wx * wy
            )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="data/albert_young/AlbertYoung_chorus.h5")
    parser.add_argument("--output",
                        default="data/albert_young/ALBERT_YOUNG_80x80.h5",
                        help="HDF5 output file")
    parser.add_argument("--binary-output", default=None,
                        help="optional legacy binary output")
    parser.add_argument("--n1", type=int, default=80)
    parser.add_argument("--n2", type=int, default=80)
    parser.add_argument("--alpha-min", type=float, default=5.0,
                        help="degrees")
    parser.add_argument("--alpha-max", type=float, default=90.0,
                        help="degrees")
    parser.add_argument("--E-min", type=float, default=0.2, dest="e_min")
    parser.add_argument("--E-max", type=float, default=5.0, dest="e_max")
    args = parser.parse_args()

    if args.n1 < 1 or args.n2 < 1:
        raise SystemExit("n1 and n2 must be positive")
    if not (0.0 <= args.alpha_min < args.alpha_max <= 90.0):
        raise SystemExit("alpha range must satisfy 0 <= min < max <= 90")
    if not (0.0 < args.e_min < args.e_max):
        raise SystemExit("energy range must satisfy 0 < min < max")

    with h5py.File(args.input, "r") as handle:
        alpha0_deg = np.asarray(handle["alpha0"][:], dtype=np.float64)
        energy = np.asarray(handle["E"][:], dtype=np.float64)
        daa = np.asarray(handle["Daa"][:], dtype=np.float64)
        dap = np.asarray(handle["Dap"][:], dtype=np.float64)
        dpp = np.asarray(handle["Dpp"][:], dtype=np.float64)

    if daa.shape != (alpha0_deg.size, energy.size):
        raise SystemExit("raw D shape is inconsistent with alpha0/E")
    if np.any(~np.isfinite(alpha0_deg)) or np.any(~np.isfinite(energy)):
        raise SystemExit("raw coordinate arrays contain non-finite values")

    alpha0 = alpha0_deg * PI / 180.0
    q1_edges = np.linspace(args.alpha_min, args.alpha_max,
                           args.n1 + 1) * PI / 180.0
    q2_edges = np.linspace(np.log(args.e_min), np.log(args.e_max),
                           args.n2 + 1)
    q1_centers = 0.5 * (q1_edges[:-1] + q1_edges[1:])
    q2_centers = 0.5 * (q2_edges[:-1] + q2_edges[1:])

    daa_i = interpolate(daa, alpha0, energy, q1_centers, q2_centers)
    dap_i = interpolate(dap, alpha0, energy, q1_centers, q2_centers)
    dpp_i = interpolate(dpp, alpha0, energy, q1_centers, q2_centers)

    e = np.exp(q2_centers)[None, :]
    p = np.asarray(e2p(e))
    dloge = np.asarray(dloge_dp(q2_centers))[None, :]
    # G depends on the equatorial pitch angle alpha0, not local alpha.
    alpha0_centers = q1_centers[:, None]

    # Jay Albert's diffusion coefficients use a normalized format; this factor
    # restores the physical momentum scale before converting to D11/D12/D22.
    denormalize = (E0 / (C * C)) ** 2 * (C * C)
    G = (p ** 2 * (1.30 - 0.56 * np.sin(alpha0_centers))
         * np.sin(alpha0_centers) * np.cos(alpha0_centers) / dloge)
    D11 = daa_i * denormalize * SECONDS_PER_DAY / p ** 2
    D12 = dap_i * denormalize * SECONDS_PER_DAY * dloge / p
    D22 = dpp_i * denormalize * SECONDS_PER_DAY * dloge ** 2

    for name, values in (("G", G), ("D11", D11), ("D12", D12), ("D22", D22)):
        if np.any(~np.isfinite(values)):
            raise SystemExit(f"{name} contains non-finite values")
    det = D11 * D22 - D12 * D12
    if np.any(det <= 0.0):
        raise SystemExit("converted diffusion tensor is not positive definite")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as handle:
        handle.attrs["format_version"] = 1
        handle.attrs["index_order"] = "row-major i*n2+j"
        handle.attrs["q1_units"] = "rad"
        handle.attrs["q2_units"] = "log(E/MeV)"
        handle.attrs["coefficient_location"] = "cell centers"
        for name, array in (
            ("q1_edges", q1_edges),
            ("q2_edges", q2_edges),
            ("G", G),
            ("D11", D11),
            ("D12", D12),
            ("D22", D22),
        ):
            handle.create_dataset(name, data=np.asarray(array, dtype=np.float64))

    print(f"wrote {output}")

    if args.binary_output:
        binary_output = Path(args.binary_output)
        binary_output.parent.mkdir(parents=True, exist_ok=True)
        with binary_output.open("wb") as stream:
            stream.write(struct.pack("<8siii", b"S2DAYD1", 1, args.n1, args.n2))
            for array in (q1_edges, q2_edges, q1_centers, q2_centers,
                          G, D11, D12, D22):
                np.ascontiguousarray(array, dtype="<f8").tofile(stream)
        print(f"wrote {binary_output} (legacy regression fixture)")
    print(f"  n1={args.n1} n2={args.n2}")
    print(f"  alpha=[{args.alpha_min},{args.alpha_max}] deg, "
          f"E=[{args.e_min},{args.e_max}] MeV")
    print(f"  G range       [{G.min():.6e}, {G.max():.6e}]")
    print(f"  D11 range     [{D11.min():.6e}, {D11.max():.6e}]")
    print(f"  D12 range     [{D12.min():.6e}, {D12.max():.6e}]")
    print(f"  D22 range     [{D22.min():.6e}, {D22.max():.6e}]")


if __name__ == "__main__":
    main()
