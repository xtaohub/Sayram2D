#!/usr/bin/env python3
"""Compare the redesigned Albert-Young case against the legacy 80x80 run.

Usage:
  python3 tests/compare_albert_young.py \
    tests/albert_young/legacy_reference.csv \
    tests/albert_young/new_reference.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

ATOL = 1e-12
RTOL = 1e-12


def read(path: Path):
    meta: dict[str, str] = {}
    header: list[str] | None = None
    rows: dict[tuple[int, int], dict[str, float]] = {}
    with path.open() as stream:
        for line in stream:
            if line.startswith("#"):
                key, _, value = line[1:].strip().partition("=")
                meta[key] = value
            elif line.startswith("i,"):
                header = line.strip().split(",")
            elif line.strip():
                if header is None:
                    raise RuntimeError(f"{path}: data before header")
                values = line.strip().split(",")
                row = {name: float(value) for name, value in zip(header, values)}
                rows[(int(row["i"]), int(row["j"]))] = row
    if header is None:
        raise RuntimeError(f"{path}: missing header")
    return meta, header, rows


def main() -> int:
    legacy_path = Path(
        sys.argv[1] if len(sys.argv) > 1
        else "tests/albert_young/legacy_reference.csv"
    )
    new_path = Path(
        sys.argv[2] if len(sys.argv) > 2
        else "tests/albert_young/new_reference.csv"
    )

    legacy_meta, legacy_header, legacy = read(legacy_path)
    new_meta, new_header, new = read(new_path)
    if "f0" not in legacy_header or "f10" not in legacy_header:
        raise SystemExit("legacy comparison file lacks f0/f10")
    if "f0" not in new_header or "f10" not in new_header:
        raise SystemExit("new comparison file lacks f0/f10")
    expected_index_order = "row-major i*n2+j; i=q1(alpha0), j=q2(logE)"
    if new_meta.get("index_order") != expected_index_order:
        raise SystemExit("new comparison file has missing/incorrect index_order")
    if new_meta.get("n1") != "80" or new_meta.get("n2") != "80":
        raise SystemExit("new comparison file has unexpected dimensions")
    if len(legacy) != 6400 or len(new) != 6400:
        raise SystemExit("expected 80x80 = 6400 cells in both files")

    max0 = 0.0
    max10 = 0.0
    for key, old in legacy.items():
        if key not in new:
            raise SystemExit(f"missing new cell {key}")
        current = new[key]
        for name, max_value in (("f0", max0), ("f10", max10)):
            a = old[name]
            b = current[name]
            diff = abs(a - b)
            scale = max(1.0, abs(a), abs(b))
            if diff > ATOL + RTOL * scale:
                raise SystemExit(
                    f"cell {key} {name}: legacy={a:.17g} new={b:.17g}"
                )
            if name == "f0":
                max0 = max(max0, diff)
            else:
                max10 = max(max10, diff)

    min0 = min(row["f0"] for row in new.values())
    min10 = min(row["f10"] for row in new.values())
    if min10 < -1e-13:
        raise SystemExit(f"new Albert-Young solution is negative: {min10}")

    print(
        f"[PASS] Albert-Young 80x80: max |delta f0|={max0:.3e}, "
        f"max |delta f10|={max10:.3e}, min f0={min0:.3e}, min f10={min10:.3e}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
