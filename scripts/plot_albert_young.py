#!/usr/bin/env python3
"""Plot the Albert-Young profiles at 0.5 and 2 MeV, 0.1 and 1 day.

Compares the standalone application output with the legacy 1-day reference.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("build/mpl").resolve()))

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import tpy.plot as tpy_plot  # noqa: E402
from tpy.plot import j_figsize  # noqa: E402

_style = (Path(tpy_plot.__file__).resolve().parents[2] /
          "config/stylelib/jdefault.mplstyle")
if _style.exists():
    plt.style.use(str(_style))
else:
    plt.style.use("default")

E0 = 0.511875  # MeV
ENERGIES = (0.5, 2.0)
FRAMES = ((1, "0.1 day"), (10, "1 day"))


def load(path: Path):
    with h5py.File(path, "r") as handle:
        alpha = handle["alpha0"][...]
        log_en = handle["logEN"][...]
        fields = {index: handle[f"f/{index}"][...] for index, _ in FRAMES}
    return alpha, log_en, fields


def profile(field: np.ndarray, log_en: np.ndarray, energy: float) -> np.ndarray:
    target = np.log(energy / E0)
    return np.array([np.interp(target, log_en, row) for row in field])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--new", type=Path,
                        default=Path("output/AlbertYoung1day/AlbertYoung1day_data.h5"))
    parser.add_argument("--legacy", type=Path,
                        default=Path("tests/albert_young/legacy_reference_1day.h5"))
    parser.add_argument("--outdir", type=Path, default=Path("figures"))
    args = parser.parse_args()

    alpha_new, logen_new, new_fields = load(args.new)
    alpha_legacy, logen_legacy, legacy_fields = load(args.legacy)
    if not np.allclose(alpha_new, alpha_legacy):
        raise SystemExit("alpha0 grids differ between new and legacy outputs")
    if not np.allclose(logen_new, logen_legacy):
        raise SystemExit("logEN grids differ between new and legacy outputs")

    new_profiles = {
        frame: {energy: profile(new_fields[frame], logen_new, energy)
                for energy in ENERGIES}
        for frame, _ in FRAMES
    }
    legacy_profiles = {
        frame: {energy: profile(legacy_fields[frame], logen_legacy, energy)
                for energy in ENERGIES}
        for frame, _ in FRAMES
    }
    args.outdir.mkdir(parents=True, exist_ok=True)

    colors = {1: "tab:blue", 10: "tab:orange"}
    fig, axes = plt.subplots(
        1, 2, figsize=j_figsize(nrows=1, ncols=2), sharex=True
    )
    handles = []
    labels = []
    for ax, panel, energy in zip(axes, ("(a)", "(b)"), ENERGIES):
        for frame, label in FRAMES:
            # Legacy first, then new: the solid curve remains visible while
            # the dashed reference is still apparent at the line edges.
            legacy_line, = ax.plot(
                alpha_legacy, legacy_profiles[frame][energy],
                color=colors[frame], linestyle="--", linewidth=1.7,
                label=f"legacy, {label}")
            new_line, = ax.plot(
                alpha_new, new_profiles[frame][energy], color=colors[frame],
                linestyle="-", linewidth=1.0, label=f"new, {label}")
            if ax is axes[0]:
                handles.extend([new_line, legacy_line])
                labels.extend([f"new, {label}", f"legacy, {label}"])
        ax.set_title(f"{panel}  $E = {energy:g}$ MeV")
        ax.set_xlabel(r"$\alpha_0$ [deg]")
        ax.set_ylabel(r"$f$" if ax is axes[0] else "")
        ax.set_xlim(alpha_new[0], alpha_new[-1])
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.02), fontsize="small")
    fig.subplots_adjust(bottom=0.24, wspace=0.35)
    fig.savefig(args.outdir / "albert_young_profiles.pdf", bbox_inches="tight")
    fig.savefig(args.outdir / "albert_young_profiles.png", dpi=300,
                bbox_inches="tight")

    fig, axes = plt.subplots(
        1, 2, figsize=j_figsize(nrows=1, ncols=2), sharex=True
    )
    for ax, panel, energy in zip(axes, ("(a)", "(b)"), ENERGIES):
        differences = {
            frame: new_profiles[frame][energy] -
                   legacy_profiles[frame][energy]
            for frame, _ in FRAMES
        }
        max_abs = max(np.max(np.abs(values)) for values in differences.values())
        exponent = int(np.floor(np.log10(max_abs))) if max_abs > 0 else 0
        scale = 10.0 ** exponent
        for frame, label in FRAMES:
            ax.plot(alpha_new, differences[frame] / scale, color=colors[frame],
                    label=f"t = {label}")
        ax.axhline(0.0, color="0.3", linewidth=0.7)
        ax.set_title(f"{panel}  $E = {energy:g}$ MeV")
        ax.set_xlabel(r"$\alpha_0$ [deg]")
        ax.set_ylabel(f"$\\Delta f\\;[10^{{{exponent}}}]$")
        ax.set_xlim(alpha_new[0], alpha_new[-1])
    axes[0].legend(frameon=False, fontsize="small", loc="lower left")
    fig.subplots_adjust(wspace=0.42)
    fig.savefig(args.outdir / "albert_young_difference.pdf", bbox_inches="tight")
    fig.savefig(args.outdir / "albert_young_difference.png", dpi=300,
                bbox_inches="tight")

    for frame, label in FRAMES:
        diff = max(np.max(np.abs(new_profiles[frame][energy] -
                                 legacy_profiles[frame][energy]))
                   for energy in ENERGIES)
        print(f"max |new - legacy| at {label} (0.5/2 MeV): {diff:.6e}")

    print(f"wrote {args.outdir / 'albert_young_profiles.pdf'}")
    print(f"wrote {args.outdir / 'albert_young_difference.pdf'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
