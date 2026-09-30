"""Regenerate the Source 2 first-in-TPC Compton-angle PDF figure."""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
ANALYSIS_DIR = HERE.parent
PROJECT = ANALYSIS_DIR.parent
sys.path.insert(0, str(HERE))

from plot_theta_first_cdf import compute_first_in_tpc


def main():
    root_path = PROJECT / "result" / "gamma_source2_1e6_clean.root"
    output = ANALYSIS_DIR / "plot_s2" / "theta_first_inTPC.png"
    theta_first, multiplicity = compute_first_in_tpc(root_path)

    selections = [
        (theta_first[multiplicity == 1], "Single-scatter (only 1 in-TPC)"),
        (theta_first[multiplicity >= 2], "Multi-scatter (>=2 in-TPC, first one)"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.5), sharey=True)
    for ax, (values, title) in zip(axes, selections):
        ax.hist(
            values,
            bins=90,
            range=(0, 180),
            histtype="step",
            color="navy",
            lw=2,
            density=True,
            label=rf"Simulation (N={len(values):,}, med={np.median(values):.0f}$^\circ$)",
        )
        ax.set_xlabel(r"$\theta$ [deg]")
        ax.set_ylabel("Probability Density")
        ax.set_title(title)
        ax.legend(fontsize=10)
        ax.set_xlim(0, 180)

    fig.suptitle(
        r"Simulation: LXe TPC gamma deflection angle $\theta=\angle(p_{in},p_{out})$",
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved: {output}")


if __name__ == "__main__":
    main()
