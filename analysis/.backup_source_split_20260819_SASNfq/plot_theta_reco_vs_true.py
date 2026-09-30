"""Compare HDF5-cluster angular reconstruction with MC truth.

Definitions
-----------
theta_true = angle(p_in1, p_out1)
theta_reco = angle(r1 - nominal_source_position, r2 - r1)

The reconstruction side uses only ER clusters stored in HDF5. TPC clusters are
ordered by their reconstructed ``er_t`` value and the earliest two define r1
and r2. ROOT gammaComp truth is used only to provide theta_true for the same
event. TrackID and gammaComp vertices are never used for reconstruction or
event selection.
"""

from pathlib import Path

import awkward as ak
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.ticker import FuncFormatter, PercentFormatter
from scipy.optimize import curve_fit


ANALYSIS_DIR = Path(__file__).resolve().parent
PROJECT_DIR = ANALYSIS_DIR.parent
PLOT_DIR = ANALYSIS_DIR / "plot"

DATASETS = {
    "Source 1": {
        "root": PROJECT_DIR / "result" / "gamma_source1_1e6_clean.root",
        "h5": PROJECT_DIR / "result" / "gamma_source1_1e6_clean.h5",
        "source_position": np.array([300.0, 0.0, -75.5]),
    },
    "Source 2": {
        "root": PROJECT_DIR / "result" / "gamma_source2_1e6_clean.root",
        "h5": PROJECT_DIR / "result" / "gamma_source2_1e6_clean.h5",
        "source_position": np.array([2000.0, 0.0, -75.5]),
    },
}
COLORS = {"Source 1": "navy", "Source 2": "crimson"}
CORE_BIN_WIDTH = 0.25
CORE_FIT_LIMIT = 15.0


def angle_deg_batch(v1, v2):
    """Row-wise angle in degrees for two arrays with shape (N, 3)."""
    denominator = np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)
    theta = np.full(len(v1), np.nan, dtype=float)
    valid = denominator > 0
    cosine = np.einsum("ij,ij->i", v1[valid], v2[valid]) / denominator[valid]
    theta[valid] = np.degrees(np.arccos(np.clip(cosine, -1, 1)))
    return theta


def gaussian(x, amplitude, mean, sigma, offset):
    """Gaussian peak plus a constant local background."""
    return offset + amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def measure_core_width(difference, bin_width, fit_limit=CORE_FIT_LIMIT):
    """Fit a Gaussian plus constant background to the resolved central core."""
    edges = np.arange(
        -fit_limit - 0.5 * bin_width,
        fit_limit + bin_width,
        bin_width,
    )
    counts, edges = np.histogram(difference, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])

    peak_index = int(np.argmax(counts))
    peak = counts[peak_index]

    selected = counts > 0
    mode = centers[peak_index]
    robust_sigma = np.quantile(np.abs(difference - np.median(difference)), 0.68)
    initial_sigma = min(max(robust_sigma, bin_width), 0.8 * fit_limit)
    initial = (peak, mode, initial_sigma, np.median(counts))
    bounds = (
        (0.0, -fit_limit, bin_width / 10.0, 0.0),
        (np.inf, fit_limit, fit_limit, np.inf),
    )
    parameters, covariance = curve_fit(
        gaussian,
        centers[selected],
        counts[selected],
        p0=initial,
        bounds=bounds,
        sigma=np.sqrt(np.maximum(counts[selected], 1.0)),
        absolute_sigma=True,
        maxfev=100000,
    )
    errors = np.sqrt(np.diag(covariance))
    prediction = gaussian(centers[selected], *parameters)
    chi2 = np.sum(
        (counts[selected] - prediction) ** 2 / np.maximum(prediction, 1.0)
    )
    chi2_ndf = chi2 / max(selected.sum() - len(parameters), 1)
    fit_sigma = parameters[2]
    gaussian_fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0)) * fit_sigma
    return {
        "edges": edges,
        "centers": centers,
        "counts": counts,
        "parameters": parameters,
        "errors": errors,
        "chi2_ndf": chi2_ndf,
        "fit_sigma": fit_sigma,
        "gaussian_fwhm": gaussian_fwhm,
    }


def load_truth_first(root_path, row_to_event_id, volume="XenonDetector"):
    """Return an eventId-indexed first-TPC-Compton truth-angle lookup."""
    branches = [
        "nGammaComp",
        "gammaCompVol",
        "gammaCompPrePx",
        "gammaCompPrePy",
        "gammaCompPrePz",
        "gammaCompPx",
        "gammaCompPy",
        "gammaCompPz",
        "gammaCompT",
    ]
    with uproot.open(root_path) as root_file:
        arrays = root_file["mcTree"].arrays(branches, library="ak")

    counts_all = np.asarray(arrays["nGammaComp"])
    if len(counts_all) != len(row_to_event_id):
        raise ValueError("ROOT tree rows and HDF5 events rows do not match")
    event_row = np.repeat(np.arange(len(counts_all)), counts_all)
    in_tpc = np.asarray(ak.flatten(arrays["gammaCompVol"] == volume))
    time = np.asarray(ak.flatten(arrays["gammaCompT"]))[in_tpc]
    event_row = event_row[in_tpc]

    def vectors(names):
        return np.column_stack(
            [np.asarray(ak.flatten(arrays[name]))[in_tpc] for name in names]
        )

    momentum_in = vectors(
        ("gammaCompPrePx", "gammaCompPrePy", "gammaCompPrePz")
    )
    momentum_out = vectors(("gammaCompPx", "gammaCompPy", "gammaCompPz"))

    order = np.lexsort((time, event_row))
    event_row = event_row[order]
    momentum_in = momentum_in[order]
    momentum_out = momentum_out[order]
    _, first_index = np.unique(event_row, return_index=True)
    event_id = row_to_event_id[event_row[first_index]]
    theta_true = angle_deg_batch(
        momentum_in[first_index], momentum_out[first_index]
    )
    lookup = np.full(int(row_to_event_id.max()) + 1, np.nan)
    lookup[event_id] = theta_true
    return lookup


def load_cluster_pairs(h5_path, source_position, volume="XenonDetector"):
    """Build one r1/r2 pair per event from time-ordered HDF5 ER clusters."""
    with h5py.File(h5_path, "r") as h5_file:
        row_to_event_id = h5_file["events"]["eventId"][:]
        clusters = h5_file["clusters"].fields(
            ("eventId", "energy", "er_x", "er_y", "er_z", "er_t", "er_volume")
        )[:]

    in_tpc = (clusters["energy"] > 0) & (
        clusters["er_volume"] == volume.encode()
    )
    event_id = clusters["eventId"][in_tpc]
    position = np.column_stack(
        [clusters[name][in_tpc] for name in ("er_x", "er_y", "er_z")]
    )
    cluster_time = clusters["er_t"][in_tpc]
    finite = np.isfinite(position).all(axis=1) & np.isfinite(cluster_time)
    event_id = event_id[finite]
    position = position[finite]
    cluster_time = cluster_time[finite]

    # er_t is an HDF5 cluster property (seconds). Its absolute unit is
    # irrelevant for sorting. Source distance only breaks exact time ties.
    source_distance = np.linalg.norm(position - source_position, axis=1)
    order = np.lexsort((source_distance, cluster_time, event_id))
    event_id = event_id[order]
    position = position[order]
    _, first_index, multiplicity = np.unique(
        event_id, return_index=True, return_counts=True
    )
    first_index = first_index[multiplicity >= 2]
    second_index = first_index + 1
    pair_event_id = event_id[first_index]
    r1 = position[first_index]
    r2 = position[second_index]
    incoming = r1 - source_position
    outgoing = r2 - r1
    theta_reco = angle_deg_batch(incoming, outgoing)
    return row_to_event_id, pair_event_id, theta_reco


def compute_theta_reco_true(
    root_path, h5_path, source_position, volume="XenonDetector"
):
    """Match cluster-based theta_reco to first-TPC-Compton theta_true by eventId."""
    row_to_event_id, pair_event_id, theta_reco = load_cluster_pairs(
        h5_path, np.asarray(source_position, dtype=float), volume
    )
    theta_true_lookup = load_truth_first(root_path, row_to_event_id, volume)
    in_range = pair_event_id < len(theta_true_lookup)
    theta_true = np.full(len(pair_event_id), np.nan)
    theta_true[in_range] = theta_true_lookup[pair_event_id[in_range]]
    valid = np.isfinite(theta_true) & np.isfinite(theta_reco)
    return theta_true[valid], theta_reco[valid]


def plot_reco_vs_true(results):
    """Create the two-source HDF5-cluster/truth comparison."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    for ax, label in zip(axes, ("Source 1", "Source 2")):
        theta_true, theta_reco = results[label]
        valid = np.isfinite(theta_true) & np.isfinite(theta_reco)
        hist = ax.hist2d(
            theta_true[valid],
            theta_reco[valid],
            bins=(90, 90),
            range=((0, 180), (0, 180)),
            cmap="plasma",
        )
        ax.plot([0, 180], [0, 180], "w--", lw=1.5, alpha=0.8)
        ax.set_xlabel(r"$\theta_{true}$ [deg]")
        ax.set_ylabel(
            r"$\theta_{reco}=\angle(r_1-r_{source},r_2-r_1)$ [deg]"
        )
        ax.set_title(f"{label}: HDF5-cluster reconstruction (N={valid.sum():,})")
        fig.colorbar(hist[3], ax=ax, label="Counts")

    fig.tight_layout()
    output = PLOT_DIR / "theta_reco_vs_true.png"
    fig.savefig(output, dpi=150, bbox_inches="tight", facecolor="white")
    # Keep the historical alias synchronized so old links do not show stale data.
    historical_output = PLOT_DIR / "theta_geo_vs_out.png"
    fig.savefig(historical_output, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"saved: {output}")
    print(f"saved historical alias: {historical_output}")


def plot_difference(results):
    """Create a normalized 1D histogram of theta_reco - theta_true."""
    differences = {}
    for label, (theta_true, theta_reco) in results.items():
        valid = np.isfinite(theta_true) & np.isfinite(theta_reco)
        differences[label] = theta_reco[valid] - theta_true[valid]

    bins = np.linspace(-180, 180, 361)  # full physical range, 1 degree/bin

    for log_scale, filename in (
        (False, "theta_reco_minus_true_hist.png"),
        (True, "theta_reco_minus_true_hist_log.png"),
    ):
        fig, ax = plt.subplots(figsize=(10, 6.5))
        for label in ("Source 1", "Source 2"):
            difference = differences[label]
            mean = np.mean(difference)
            std = np.std(difference)
            color = COLORS[label]
            weights = np.ones(len(difference)) * 100.0 / len(difference)
            display_label = label.replace(" ", "")

            # # Colored mean +/- 1 std interval; the two sources nearly overlap.
            # ax.axvspan(mean - std, mean + std, color=color, alpha=0.08, zorder=0)
            # ax.axvline(mean - std, color=color, ls=":", lw=1.0, alpha=0.75)
            # ax.axvline(mean + std, color=color, ls=":", lw=1.0, alpha=0.75)
            ax.hist(
                difference,
                bins=bins,
                weights=weights,
                histtype="step",
                lw=2,
                color=color,
                label=(
                    f"{display_label} (mean={mean:.2f}\N{DEGREE SIGN})"
                    # f"std={std:.2f}\N{DEGREE SIGN})"
                ),
            )

        ax.axvline(0, color="0.35", ls="--", lw=1.2)
        ax.set_xlabel(r"$\Delta\theta=\theta_{reco}-\theta_{true}$ [deg]")
        ax.set_ylabel(r"Events per $1^\circ$ bin")
        title = "HDF5-cluster angular reconstruction difference"
        ax.set_title(title + (" (log y-axis)" if log_scale else ""))
        ax.set_xlim(-180, 180)
        ax.set_xticks(np.arange(-180, 181, 30))
        if log_scale:
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(
                FuncFormatter(lambda value, _: f"{value:g}%")
            )
        else:
            ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=9)
        fig.tight_layout()

        output = PLOT_DIR / filename
        fig.savefig(output, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"saved: {output}")


def get_differences(results):
    """Return finite theta_reco - theta_true arrays for both sources."""
    differences = {}
    for label, (theta_true, theta_reco) in results.items():
        valid = np.isfinite(theta_true) & np.isfinite(theta_reco)
        differences[label] = theta_reco[valid] - theta_true[valid]
    return differences


def plot_core_gaussian_fwhm(results):
    """Plot the resolved Source2 cluster residual and a Gaussian core fit."""
    differences = get_differences(results)
    label = "Source 2"
    difference = differences[label]
    result = measure_core_width(difference, CORE_BIN_WIDTH)
    measurements = {label: result}
    scale = 100.0 / len(difference)
    counts_percent = result["counts"] * scale
    fit_x = np.linspace(-CORE_FIT_LIMIT, CORE_FIT_LIMIT, 2001)
    fit_y = gaussian(fit_x, *result["parameters"]) * scale
    sigma = result["fit_sigma"]
    gaussian_fwhm = result["gaussian_fwhm"]
    chi2_ndf = result["chi2_ndf"]
    absolute = np.abs(difference)
    median_abs = np.median(absolute)
    containment68 = np.quantile(absolute, 0.68)
    containment90 = np.quantile(absolute, 0.90)

    fig, ax = plt.subplots(figsize=(8.6, 6.3))
    ax.stairs(
        counts_percent,
        result["edges"],
        color=COLORS[label],
        lw=2.0,
        label=f"HDF5 clusters ({CORE_BIN_WIDTH:.2f}\N{DEGREE SIGN} bins)",
    )
    ax.plot(
        fit_x,
        fit_y,
        color="darkorange",
        lw=2.0,
        label="Gaussian fit",
    )
    ax.axvline(0, color="0.45", ls=":", lw=1.0)
    ax.set_xlim(-CORE_FIT_LIMIT, CORE_FIT_LIMIT)
    ax.set_ylim(0, max(counts_percent.max(), fit_y.max()) * 1.13)
    ax.set_xlabel(r"$\Delta\theta=\theta_{reco}-\theta_{true}$ [deg]")
    ax.set_ylabel(rf"Events per ${CORE_BIN_WIDTH:.2f}^\circ$ bin [%]")
    ax.set_title("Source2 HDF5-cluster residual and Gaussian core fit")
    ax.grid(alpha=0.22)
    ax.text(
        0.96,
        0.95,
        (
            rf"Gaussian $\sigma_{{fit}}={sigma:.3f}^\circ$"
            "\n"
            rf"Gaussian FWHM $={gaussian_fwhm:.3f}^\circ$"
            "\n"
            rf"$\chi^2/\mathrm{{ndf}}={chi2_ndf:.1f}$ (poor fit)"
            "\n"
            rf"median $|\Delta\theta|={median_abs:.2f}^\circ$"
            "\n"
            rf"68% / 90%: ${containment68:.2f}^\circ$ / ${containment90:.2f}^\circ$"
        ),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.75", alpha=0.95),
    )
    ax.legend(loc="upper left", fontsize=9.5)

    print(
        f"{label} core ({CORE_BIN_WIDTH:.2f} deg bins): "
        f"Gaussian sigma={sigma:.6f} deg; chi2/ndf={chi2_ndf:.3f}; "
        f"Gaussian FWHM={gaussian_fwhm:.6f} deg; "
        f"median abs={median_abs:.6f} deg; "
        f"containment68/90={containment68:.6f}/{containment90:.6f} deg"
    )

    fig.text(
        0.5,
        0.01,
        "r1/r2: time-ordered HDF5 TPC ER clusters; incident direction: nominal source. "
        "No TrackID selection; ROOT supplies only the matched truth angle.",
        ha="center",
        fontsize=10,
        color="0.3",
    )
    fig.tight_layout(rect=(0, 0.055, 1, 1))
    output = PLOT_DIR / "theta_reco_minus_true_gaussian_fwhm.png"
    fig.savefig(output, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"saved: {output}")
    return measurements

def main():
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    results = {}
    for label, dataset in DATASETS.items():
        print(f"reading {label}: {dataset['root']} + {dataset['h5']}")
        theta_true, theta_reco = compute_theta_reco_true(
            dataset["root"], dataset["h5"], dataset["source_position"]
        )
        results[label] = (theta_true, theta_reco)

        valid = np.isfinite(theta_true) & np.isfinite(theta_reco)
        difference = theta_reco[valid] - theta_true[valid]
        print(
            f"  N={valid.sum():,}; median true/reco="
            f"{np.median(theta_true[valid]):.2f}/{np.median(theta_reco[valid]):.2f} deg; "
            f"mean signed diff={np.mean(difference):.3f} deg; "
            f"mean abs diff={np.mean(np.abs(difference)):.3f} deg; "
            f"median abs diff={np.median(np.abs(difference)):.3f} deg; "
            f"68% containment={np.quantile(np.abs(difference), 0.68):.3f} deg; "
            f"corr={np.corrcoef(theta_true[valid], theta_reco[valid])[0, 1]:.4f}"
        )

    plot_reco_vs_true(results)
    plot_difference(results)
    plot_core_gaussian_fwhm(results)


if __name__ == "__main__":
    main()
