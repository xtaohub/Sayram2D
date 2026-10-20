#!/usr/bin/env python3
"""Compare the redesigned solver against the committed legacy references.

Usage:
  python3 tests/compare_legacy.py build/new_reference tests/legacy/reference
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

CASES = [
    "uniform_diag_all_dirichlet",
    "nonuniform_diag_mixed",
    "nonuniform_offdiag_mixed",
    "uniform_offdiag_all_zeroflux",
    "nonuniform_offdiag_loss",
]

ATOL = 1e-12
RTOL = 1e-12


def read_csv(path: Path):
    meta: dict[str, str] = {}
    rows: list[list[float]] = []
    header: list[str] | None = None
    with path.open() as stream:
        for line in stream:
            if line.startswith("#"):
                key, _, value = line[1:].strip().partition("=")
                meta[key] = value
            elif line.startswith("i,"):
                header = line.strip().split(",")
            elif line.strip():
                rows.append([float(value) for value in line.strip().split(",")])
    if header is None:
        raise RuntimeError(f"{path}: missing header")
    return meta, header, rows


def fail(message: str) -> None:
    raise RuntimeError(message)


def compare_case(new_dir: Path, legacy_dir: Path, name: str) -> float:
    legacy_meta, legacy_header, legacy_rows = read_csv(
        legacy_dir / f"{name}.csv"
    )
    new_meta, new_header, new_rows = read_csv(new_dir / f"{name}.csv")

    for key in ("n1", "n2", "dt", "steps"):
        if legacy_meta[key] != new_meta[key]:
            fail(f"{name}: metadata {key} differs")
    if legacy_header != new_header:
        fail(f"{name}: CSV header differs")
    if len(legacy_rows) != len(new_rows):
        fail(f"{name}: row count differs")

    legacy_f = legacy_header.index("f0")
    max_diff = 0.0
    for row_a, row_b in zip(legacy_rows, new_rows):
        for index in range(legacy_f, len(legacy_header)):
            a = row_a[index]
            b = row_b[index]
            scale = max(1.0, abs(a), abs(b))
            diff = abs(a - b)
            max_diff = max(max_diff, diff)
            if diff > ATOL + RTOL * scale:
                fail(
                    f"{name}: f column {legacy_header[index]} differs at "
                    f"cell ({row_a[0]:.0f},{row_a[1]:.0f}): "
                    f"legacy={a:.17g} new={b:.17g}"
                )
    return max_diff


def main() -> int:
    new_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "build/new_reference")
    legacy_dir = Path(
        sys.argv[2] if len(sys.argv) > 2 else "tests/legacy/reference"
    )

    failures: list[str] = []
    for name in CASES:
        try:
            diff = compare_case(new_dir, legacy_dir, name)
            print(f"[PASS] {name}: max |delta f| = {diff:.3e}")
        except Exception as error:  # noqa: BLE001 - report all failures
            failures.append(str(error))
            print(f"[FAIL] {name}: {error}")

    # Conservation and positivity checks on the new outputs.
    try:
        meta, _, _ = read_csv(new_dir / "uniform_offdiag_all_zeroflux.csv")
        mass0 = float(meta["mass0"])
        massN = float(meta["massN"])
        relative = abs(massN - mass0) / abs(mass0)
        if relative > 1e-12:
            fail(f"ZeroFlux mass conservation failed: relative={relative:.3e}")
        if float(meta["minN"]) < -1e-13:
            fail("ZeroFlux solution became negative")
        print(f"[PASS] ZeroFlux conservation: relative mass change {relative:.3e}")
    except Exception as error:  # noqa: BLE001
        failures.append(str(error))
        print(f"[FAIL] ZeroFlux conservation: {error}")

    try:
        mixed_meta, _, _ = read_csv(new_dir / "nonuniform_offdiag_mixed.csv")
        source_meta, _, _ = read_csv(
            new_dir / "nonuniform_offdiag_source.csv"
        )
        if float(source_meta["massN"]) <= float(mixed_meta["massN"]):
            fail("negative inverse loss did not add mass")
        if float(source_meta["minN"]) < -1e-13:
            fail("source solution became negative")
        print(
            "[PASS] negative inverse loss source: "
            f"mass {mixed_meta['massN']} -> {source_meta['massN']}"
        )
    except Exception as error:  # noqa: BLE001
        failures.append(str(error))
        print(f"[FAIL] negative inverse loss source: {error}")

    # The legacy source reference is expected to differ: the legacy code
    # implemented negative tau through the diagonal, while the redesigned
    # solver uses the documented explicit RHS source. Report it, do not fail.
    try:
        legacy_source, _, _ = read_csv(
            legacy_dir / "nonuniform_offdiag_source_legacy.csv"
        )
        new_source, _, _ = read_csv(new_dir / "nonuniform_offdiag_source.csv")
        delta = abs(float(legacy_source["massN"]) - float(new_source["massN"]))
        print(
            f"[INFO] legacy-vs-new negative-tau mass delta: {delta:.3e} "
            "(expected; documented source treatment differs)"
        )
    except Exception as error:  # noqa: BLE001
        print(f"[INFO] legacy source comparison skipped: {error}")

    if failures:
        print(f"\n{len(failures)} check(s) failed.")
        return 1
    print("\nAll legacy comparisons passed.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
