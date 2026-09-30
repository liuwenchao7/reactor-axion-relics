"""Plot truth-angle CDFs and reconstructed-cluster PDF/CDF figures.

The three-panel first-angle CDF remains an explicitly labelled ROOT-truth
diagnostic.  The multi/double PDF+CDF figures use reconstructed HDF5 clusters:

    theta_reco = angle(known_beam_axis, r2 - r1)

Here multi-scatter means at least two valid TPC ER clusters, while
double-scatter means exactly two.  The first two clusters are ordered by er_t.
For comparison, theta_doublefirst_inTPC_pdf_cdf.png deliberately restores the
pre-cluster truth definition: exactly two ROOT in-TPC Compton records and the
truth momentum angle of the earliest record.  The cluster-reconstructed double
sample is retained under theta_double_cluster_reco_inTPC_pdf_cdf.png.
"""

from pathlib import Path

import awkward as ak
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot


HERE = Path(__file__).resolve().parent
ANALYSIS_DIR = HERE.parent if HERE.name == "python_code" else HERE
PROJECT = ANALYSIS_DIR.parent
SOURCE_FILES = {
    "Source 1": PROJECT / "result" / "gamma_source1_1e6_clean.root",
    "Source 2": PROJECT / "result" / "gamma_source2_1e6_clean.root",
}
H5_FILES = {
    "Source 1": PROJECT / "result" / "gamma_source1_1e6_clean.h5",
    "Source 2": PROJECT / "result" / "gamma_source2_1e6_clean.h5",
}
BEAM_AXIS = np.array([-1.0, 0.0, 0.0])
SOURCE_CONFIG = {
    "Source 1": {
        "color": "crimson",
        "display_label": "Source 1",
        "plot_dir": ANALYSIS_DIR / "plot_s1",
    },
    "Source 2": {
        "color": "navy",
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


def angle_deg_batch(v1, v2):
    """Return row-wise angles in degrees for two (N, 3) arrays."""
    denominator = np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)
    theta = np.full(len(v1), np.nan, dtype=float)
    valid = denominator > 0
    cosine = np.einsum("ij,ij->i", v1[valid], v2[valid]) / denominator[valid]
    theta[valid] = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    return theta


def compute_cluster_theta_reco(
    h5_path,
    beam_axis=BEAM_AXIS,
    volume="XenonDetector",
    creator_process=None,
):
    """Return theta_reco and multiplicity for a selected TPC ER-cluster set."""
    with h5py.File(h5_path, "r") as h5_file:
        clusters = h5_file["clusters"].fields(
            (
                "eventId",
                "energy",
                "er_x",
                "er_y",
                "er_z",
                "er_t",
                "er_volume",
                "er_creatorProcess",
            )
        )[:]

    selected = (clusters["energy"] > 0) & (
        clusters["er_volume"] == volume.encode()
    )
    if creator_process is not None:
        creator_process = (
            creator_process.encode()
            if isinstance(creator_process, str)
            else creator_process
        )
        selected &= clusters["er_creatorProcess"] == creator_process
    event_id = clusters["eventId"][selected]
    position = np.column_stack(
        [clusters[name][selected] for name in ("er_x", "er_y", "er_z")]
    )
    cluster_time = clusters["er_t"][selected]
    finite = np.isfinite(position).all(axis=1) & np.isfinite(cluster_time)
    event_id = event_id[finite]
    position = position[finite]
    cluster_time = cluster_time[finite]

    # er_t determines the interaction order.  For exact time ties, choose the
    # cluster farther upstream along the known beam axis first.
    beam_axis = np.asarray(beam_axis, dtype=float)
    beam_axis /= np.linalg.norm(beam_axis)
    beam_coordinate = position @ beam_axis
    order = np.lexsort((beam_coordinate, cluster_time, event_id))
    event_id = event_id[order]
    position = position[order]
    _, first_index, multiplicity = np.unique(
        event_id, return_index=True, return_counts=True
    )
    eligible = multiplicity >= 2
    first_index = first_index[eligible]
    multiplicity = multiplicity[eligible]

    r1 = position[first_index]
    r2 = position[first_index + 1]
    outgoing = r2 - r1
    incoming = np.broadcast_to(beam_axis, outgoing.shape)
    theta_reco = angle_deg_batch(incoming, outgoing)
    valid = np.isfinite(theta_reco)
    return theta_reco[valid], multiplicity[valid]


def ecdf(values):
    """Return x and F(x) for an empirical cumulative distribution."""
    x = np.sort(np.asarray(values, dtype=float))
    x = x[np.isfinite(x)]
    y = np.arange(1, len(x) + 1, dtype=float) / len(x)
    return x, y


def plot_precluster_truth_double(
    values, color, display_label, plot_dir
):
    """Plot the pre-cluster double-first definition based on ROOT truth records."""
    median = np.median(values)
    legend_label = (
        f"{display_label} (N={len(values):,}, med={median:.1f}$^\\circ$)"
    )
    fig, (ax_pdf, ax_cdf) = plt.subplots(
        1, 2, figsize=(13.5, 5.3), sharex=True
    )
    bins = np.linspace(0, 180, 91)
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
    ax_cdf.set_ylabel(r"Cumulative probability $P(\Theta_{first}\leq\theta)$")
    ax_cdf.set_ylim(0, 1.0)
    ax_cdf.set_yticks(np.arange(0, 1.01, 0.2))
    ax_cdf.axhline(
        0.5, color="#6b7280", lw=1.0, ls=(0, (4, 4)), alpha=0.65
    )
    ax_cdf.legend(loc="lower right", fontsize=9.4, framealpha=0.95)
    fig.suptitle(
        display_label
        + ": Pre-cluster double-first truth diagnostic"
        + "\n"
        + r"Exactly 2 in-TPC Compton truth records; "
        + r"$\theta_{first}=\angle(p_{in},p_{out})$ for the earliest record",
        fontsize=13.5,
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    output = plot_dir / "theta_doublefirst_inTPC_pdf_cdf.png"
    fig.savefig(output, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(
        f"saved pre-cluster truth double: N={len(values):,}; "
        f"median theta_first={median:.3f} deg; {output}"
    )


def main():
    results = {}
    cluster_results = {}
    compton_cluster_results = {}
    for label, root_path in SOURCE_FILES.items():
        print(f"reading {label}: {root_path}")
        results[label] = compute_first_in_tpc(root_path)
        cluster_results[label] = compute_cluster_theta_reco(
            H5_FILES[label], BEAM_AXIS
        )
        compton_cluster_results[label] = compute_cluster_theta_reco(
            H5_FILES[label], BEAM_AXIS, creator_process="compt"
        )
        SOURCE_CONFIG[label]["plot_dir"].mkdir(parents=True, exist_ok=True)

    selections = [
        ("all", "All events (first in-TPC truth record)"),
        ("single", "1 in-TPC truth Compton record"),
        ("multi", r"$\geq$2 in-TPC truth Compton records (first)"),
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
            + r"ROOT-truth CDF of LXe TPC $\gamma$ deflection angle "
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

        # Restore the pre-cluster double-first product requested for comparison:
        # exactly two ROOT in-TPC Compton truth records, earliest truth angle.
        plot_precluster_truth_double(
            theta_first[multiplicity == 2], color, display_label, plot_dir
        )

        # Reconstructed-cluster PDF + CDF views.  The historical "multifirst"
        # filename is retained, although the plotted quantity is theta_reco.
        theta_reco, cluster_multiplicity = cluster_results[label]
        theta_reco_compt, compton_multiplicity = compton_cluster_results[label]
        cluster_selections = [
            (
                "multi",
                theta_reco,
                cluster_multiplicity >= 2,
                r"Multi-scatter ($N_{cluster}\geq2$)",
                "theta_multifirst_inTPC_pdf_cdf.png",
            ),
            (
                "double",
                theta_reco,
                cluster_multiplicity == 2,
                r"Double-scatter ($N_{cluster}=2$)",
                "theta_double_cluster_reco_inTPC_pdf_cdf.png",
            ),
            (
                "compton-cluster multi",
                theta_reco_compt,
                compton_multiplicity >= 2,
                r"Compton-like multi ($N_{cluster}^{compt}\geq2$; "
                r"er_creatorProcess = compt)",
                "theta_multifirst_compton_cluster_inTPC_pdf_cdf.png",
            ),
        ]
        for selection, selected_theta, mask, sample_title, output_name in cluster_selections:
            values = selected_theta[mask]
            median = np.median(values)
            legend_label = (
                f"{display_label} (N={len(values):,}, med={median:.1f}$^\\circ$)"
            )

            fig, (ax_pdf, ax_cdf) = plt.subplots(
                1, 2, figsize=(13.5, 5.3), sharex=True
            )
            bins = np.linspace(0, 180, 91)
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
                ax.set_xlabel(r"$\theta_{\rm reco}$ [deg]")
                ax.grid(True, color="#d7dde3", lw=0.7, alpha=0.7)

            ax_pdf.set_title("Probability density (PDF)")
            ax_pdf.set_ylabel(r"Probability density [deg$^{-1}$]")
            ax_pdf.legend(loc="upper right", fontsize=9.4, framealpha=0.95)

            ax_cdf.set_title("Cumulative distribution (CDF)")
            ax_cdf.set_ylabel(
                r"Cumulative probability $P(\Theta_{reco} \leq \theta)$"
            )
            ax_cdf.set_ylim(0, 1.0)
            ax_cdf.set_yticks(np.arange(0, 1.01, 0.2))
            ax_cdf.axhline(
                0.5,
                color="#6b7280",
                lw=1.0,
                ls=(0, (4, 4)),
                alpha=0.65,
            )
            ax_cdf.legend(loc="lower right", fontsize=9.4, framealpha=0.95)

            fig.suptitle(
                display_label
                + ": "
                + sample_title
                + "\n"
                + r"Reconstructed-cluster angle "
                + r"$\theta_{reco}=\angle(\hat d_{beam},r_2-r_1)$, "
                + r"$\hat d_{beam}=(-1,0,0)$",
                fontsize=13.5,
                y=0.995,
            )
            fig.tight_layout(rect=(0, 0, 1, 0.88))
            output = plot_dir / output_name
            fig.savefig(output, dpi=180, bbox_inches="tight", facecolor="white")
            plt.close(fig)
            print(
                f"{label} {selection}: N={len(values):,}; "
                f"median theta_reco={median:.3f} deg; saved: {output}"
            )


if __name__ == "__main__":
    main()
