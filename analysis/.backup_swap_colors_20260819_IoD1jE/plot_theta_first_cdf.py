"""Plot CDFs of the first in-TPC Compton deflection angle.

The event selections and angle definition match gamma_beam_analysis_R10.ipynb:
volume == XenonDetector, theta = angle(p_in, p_out), and the earliest valid
TPC Compton record in each event.
"""

from pathlib import Path

import awkward as ak
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot


HERE = Path(__file__).resolve().parent
ANALYSIS_DIR = HERE.parent if HERE.name == ".py" else HERE
PROJECT = ANALYSIS_DIR.parent
SOURCE_FILES = {
    "Source 1": PROJECT / "result" / "gamma_source1_1e6_clean.root",
    "Source 2": PROJECT / "result" / "gamma_source2_1e6_clean.root",
}
SOURCE_CONFIG = {
    "Source 1": {
        "color": "#173f7a",
        "display_label": "Source 1",
        "plot_dir": ANALYSIS_DIR / "plot_s1",
    },
    "Source 2": {
        "color": "#c62f48",
        "display_label": "Simulation",
        "plot_dir": ANALYSIS_DIR / "plot_s2",
    },
}


def compute_first_in_tpc(root_path, volume="XenonDetector"):
    """Return first-angle values and per-event in-TPC Compton multiplicity."""
    with uproot.open(root_path) as root_file:
        tree = root_file["mcTree"]
        arrays = tree.arrays(
            [
                "gammaCompVol",
                "gammaCompPrePx",
                "gammaCompPrePy",
                "gammaCompPrePz",
                "gammaCompPx",
                "gammaCompPy",
                "gammaCompPz",
                "gammaCompT",
            ],
            library="ak",
        )

    vol = arrays["gammaCompVol"]
    prex = arrays["gammaCompPrePx"]
    prey = arrays["gammaCompPrePy"]
    prez = arrays["gammaCompPrePz"]
    px = arrays["gammaCompPx"]
    py = arrays["gammaCompPy"]
    pz = arrays["gammaCompPz"]
    time = arrays["gammaCompT"]

    in_tpc = vol == volume
    counts = np.asarray(ak.sum(in_tpc, axis=1))

    dot = prex * px + prey * py + prez * pz
    norm_in = np.sqrt(prex**2 + prey**2 + prez**2)
    norm_out = np.sqrt(px**2 + py**2 + pz**2)
    denominator = norm_in * norm_out
    safe_denominator = ak.where(denominator > 0, denominator, 1.0)
    cosine = ak.where(denominator > 0, dot / safe_denominator, np.nan)
    theta = np.degrees(
        np.arccos(
            np.clip(
                np.asarray(ak.flatten(cosine)),
                -1,
                1,
            )
        )
    )

    records_per_event = np.asarray(ak.num(vol, axis=1))
    event_row = np.repeat(np.arange(len(records_per_event)), records_per_event)
    flat_time = np.asarray(ak.flatten(time))
    flat_in_tpc = np.asarray(ak.flatten(in_tpc))
    valid = flat_in_tpc & np.isfinite(theta)

    theta_all = theta[valid]
    event_all = event_row[valid]
    time_all = flat_time[valid]
    order = np.lexsort((time_all, event_all))
    sorted_event = event_all[order]
    first_mask = np.concatenate(([True], sorted_event[1:] != sorted_event[:-1]))
    theta_first = theta_all[order][first_mask]

    return theta_first, counts[counts > 0]


def ecdf(values):
    """Return x and F(x) for an empirical cumulative distribution."""
    x = np.sort(np.asarray(values, dtype=float))
    x = x[np.isfinite(x)]
    y = np.arange(1, len(x) + 1, dtype=float) / len(x)
    return x, y


def main():
    results = {}
    for label, root_path in SOURCE_FILES.items():
        print(f"reading {label}: {root_path}")
        results[label] = compute_first_in_tpc(root_path)
        SOURCE_CONFIG[label]["plot_dir"].mkdir(parents=True, exist_ok=True)

    selections = [
        ("all", "All events (first in-TPC)"),
        ("single", "Single-scatter (only 1 in-TPC)"),
        ("multi", r"Multi-scatter ($\geq$2 in-TPC, first one)"),
    ]
    for label in ("Source 1", "Source 2"):
        color = SOURCE_CONFIG[label]["color"]
        display_label = SOURCE_CONFIG[label]["display_label"]
        plot_dir = SOURCE_CONFIG[label]["plot_dir"]
        theta_first, multiplicity = results[label]

        fig, axes = plt.subplots(
            1, 3, figsize=(17, 5.5), sharex=True, sharey=True
        )
        for ax, (selection, title) in zip(axes, selections):
            theta_first, multiplicity = results[label]
            if selection == "all":
                values = theta_first
            elif selection == "single":
                values = theta_first[multiplicity == 1]
            else:
                values = theta_first[multiplicity >= 2]

            x, y = ecdf(values)
            median = np.median(values)
            ax.step(
                x,
                y,
                where="post",
                color=color,
                lw=2.2,
                label=f"{display_label} (N={len(values):,}, med={median:.0f}$^\\circ$)",
            )

            ax.axhline(
                0.5, color="#6b7280", lw=1.0, ls=(0, (4, 4)), alpha=0.65
            )
            ax.set_title(title, fontsize=12)
            ax.set_xlabel(r"$\theta$ [deg]")
            ax.set_xlim(0, 180)
            ax.set_ylim(0, 1.0)
            ax.set_xticks(np.arange(0, 181, 30))
            ax.set_yticks(np.arange(0, 1.01, 0.2))
            ax.grid(True, color="#d7dde3", lw=0.7, alpha=0.7)
            ax.legend(loc="lower right", fontsize=9.2, framealpha=0.95)

        axes[0].set_ylabel(r"Cumulative probability $P(\Theta \leq \theta)$")
        fig.suptitle(
            display_label
            + ": "
            + r"CDF of LXe TPC $\gamma$ deflection angle "
            + r"$\theta=\angle(p_{in},p_{out})$ "
            + r"(Compton, volume = XenonDetector)",
            fontsize=15,
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        output_cdf = plot_dir / "theta_first_inTPC_cdf.png"
        fig.savefig(output_cdf, dpi=180, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"saved: {output_cdf}")

        # A compact PDF + CDF view for the multi-first selection.
        fig, (ax_pdf, ax_cdf) = plt.subplots(
            1, 2, figsize=(13.5, 5.3), sharex=True
        )
        bins = np.linspace(0, 180, 91)
        values = theta_first[multiplicity >= 2]
        median = np.median(values)
        legend_label = (
            f"{display_label} (N={len(values):,}, med={median:.0f}$^\\circ$)"
        )

        ax_pdf.hist(
            values,
            bins=bins,
            density=True,
            histtype="step",
            color=color,
            lw=2.2,
            label=legend_label,
        )
        x, y = ecdf(values)
        ax_cdf.step(
            x,
            y,
            where="post",
            color=color,
            lw=2.2,
            label=legend_label,
        )

        for ax in (ax_pdf, ax_cdf):
            ax.set_xlim(0, 180)
            ax.set_xticks(np.arange(0, 181, 30))
            ax.set_xlabel(r"$\theta_{\rm first}$ [deg]")
            ax.grid(True, color="#d7dde3", lw=0.7, alpha=0.7)

        ax_pdf.set_title("Probability density (PDF)")
        ax_pdf.set_ylabel(r"Probability density [deg$^{-1}$]")
        ax_pdf.legend(loc="upper right", fontsize=9.4, framealpha=0.95)

        ax_cdf.set_title("Cumulative distribution (CDF)")
        ax_cdf.set_ylabel(r"Cumulative probability $P(\Theta \leq \theta)$")
        ax_cdf.set_ylim(0, 1.0)
        ax_cdf.set_yticks(np.arange(0, 1.01, 0.2))
        ax_cdf.axhline(
            0.5, color="#6b7280", lw=1.0, ls=(0, (4, 4)), alpha=0.65
        )
        ax_cdf.legend(loc="lower right", fontsize=9.4, framealpha=0.95)

        fig.suptitle(
            display_label
            + ": "
            + r"multi-first sample: first in-TPC Compton angle for events with "
            + r"$\geq 2$ in-TPC Compton interactions",
            fontsize=14.5,
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        output_multi = plot_dir / "theta_multifirst_inTPC_pdf_cdf.png"
        fig.savefig(
            output_multi, dpi=180, bbox_inches="tight", facecolor="white"
        )
        plt.close(fig)
        print(f"saved: {output_multi}")


if __name__ == "__main__":
    main()
