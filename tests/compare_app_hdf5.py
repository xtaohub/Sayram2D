#!/usr/bin/env python3
"""Check the standalone app's HDF5 result against the legacy reference."""

from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np

H5_PATH = Path("build/app_test_output/AppTest/AppTest_data.h5")
LEGACY_PATH = Path("tests/albert_young/legacy_reference.csv")
ATOL = 1e-12
RTOL = 1e-12


def read_legacy(path: Path, n1: int, n2: int):
    alpha = np.empty(n1)
    logen = np.empty(n2)
    f0 = np.empty((n1, n2))
    f10 = np.empty((n1, n2))
    for line in path.read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("i,"):
            continue
        i, j, a, b, c, d = line.split(",")
        i, j = int(i), int(j)
        alpha[i] = float(a)
        logen[j] = float(b)
        f0[i, j] = float(c)
        f10[i, j] = float(d)
    return alpha, logen, f0, f10


def main() -> int:
    h5_path = Path(sys.argv[1]) if len(sys.argv) > 1 else H5_PATH
    legacy_path = Path(sys.argv[2]) if len(sys.argv) > 2 else LEGACY_PATH
    if not h5_path.exists():
        raise SystemExit(f"missing HDF5 output: {h5_path}")

    with h5py.File(h5_path, "r") as handle:
        for name in ("alpha0", "logEN", "t", "f"):
            if name not in handle:
                raise SystemExit(f"missing HDF5 dataset/group /{name}")
        alpha = handle["alpha0"][...]
        logen = handle["logEN"][...]
        times = handle["t"][...]
        frames = handle["f"]
        expected_frames = [str(i) for i in range(11)]
        if sorted(frames.keys(), key=int) != expected_frames:
            raise SystemExit("HDF5 /f frames are not 0..10")
        f0 = frames["0"][...]
        f10 = frames["10"][...]

    if alpha.shape != (80,) or logen.shape != (80,):
        raise SystemExit("coordinate datasets have wrong shapes")
    if times.shape != (11,):
        raise SystemExit("time dataset has wrong shape")
    if f0.shape != (80, 80) or f10.shape != (80, 80):
        raise SystemExit("field datasets have wrong shapes")
    if not np.isclose(times[0], 0.0) or not np.isclose(times[-1], 0.02):
        raise SystemExit(f"unexpected saved times: {times[0]}..{times[-1]}")

    ref_alpha, ref_logen, ref_f0, ref_f10 = read_legacy(legacy_path, 80, 80)
    alpha_error = np.max(np.abs(alpha - ref_alpha))
    logen_error = np.max(np.abs(logen - ref_logen))
    f0_error = np.max(np.abs(f0 - ref_f0))
    f10_error = np.max(np.abs(f10 - ref_f10))
    if alpha_error > ATOL:
        raise SystemExit(f"alpha0 mismatch: {alpha_error:.3e}")
    if logen_error > ATOL:
        raise SystemExit(f"logEN mismatch: {logen_error:.3e}")
    for name, error in (("f0", f0_error), ("f10", f10_error)):
        if error > ATOL + RTOL * max(1.0, np.max(np.abs(ref_f0 if name == "f0" else ref_f10))):
            raise SystemExit(f"{name} mismatch: {error:.3e}")
    if np.min(f10) < -1e-13:
        raise SystemExit(f"negative final distribution: {np.min(f10)}")

    print(
        "[PASS] standalone HDF5: "
        f"max|df0|={f0_error:.3e}, max|df10|={f10_error:.3e}, "
        f"max|dalpha|={alpha_error:.3e}, max|dlogEN|={logen_error:.3e}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
