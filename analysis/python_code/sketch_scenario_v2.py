"""Clear two-panel sketch of the R10_v2_clean reactor-gamma scenario.

Panel A shows the source and the clean LXe geometry in an x-z side view.
Panel B connects the Monte-Carlo deflection angle to the experimentally
reconstructable two-cluster geometry.  Dimensions are in centimetres.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Arc, Circle, Patch, Rectangle


# Geometry from config/cevns_clean.xml.  Child-volume z positions are relative
# to XenonDetector, so electrode coordinates must include XD_ZC.
XD_R = 10.2579
XD_H = 17.885
XD_ZC = -8.9425
VETO_R0 = 8.25
VETO_R1 = XD_R
VETO_H = 17.635
VETO_ZC = XD_ZC + 0.125
TPC_R = 8.25
SOURCE_Z = -7.55
SOURCE_R = 8.0

SOURCES = {
    "source1": {
        "x": 30.0,
        "display_x": 30.0,
        "cone_deg": 18.0,
        "output": "sketch_scenario_v2.png",
    },
    "source2": {
        "x": 200.0,
        "display_x": 30.0,
        "cone_deg": 0.1,
        "output": "sketch_scenario_v2_source2.png",
    },
}

COL = {
    "world": "#f5f7f9",
    "active": "#dceeff",
    "active_edge": "#2563a6",
    "veto": "#8fc1e8",
    "teflon": "#f5d77b",
    "copper": "#b87333",
    "source": "#24863b",
    "gamma": "#e58b00",
    "electron": "#c33d2e",
    "cluster2": "#7551b2",
    "reco": "#1565c0",
    "text": "#17212b",
    "muted": "#5d6975",
}


def draw_detector(ax, *, detailed=True):
    """Draw the LXe detector using global x-z coordinates."""
    bottom = XD_ZC - XD_H / 2
    veto_bottom = VETO_ZC - VETO_H / 2

    ax.add_patch(
        Rectangle(
            (-XD_R, bottom),
            2 * XD_R,
            XD_H,
            facecolor=COL["veto"],
            edgecolor=COL["active_edge"],
            lw=2.0,
            zorder=2,
        )
    )
    ax.add_patch(
        Rectangle(
            (-TPC_R, veto_bottom),
            2 * TPC_R,
            VETO_H,
            facecolor=COL["active"],
            edgecolor=COL["active_edge"],
            lw=1.5,
            zorder=3,
        )
    )

    if detailed:
        # Side Teflon: 8.25--8.50 cm; copper shaping rings: 8.50--9.00 cm.
        for sign in (-1, 1):
            x0 = sign * TPC_R if sign > 0 else sign * (TPC_R + 0.25)
            ax.add_patch(
                Rectangle(
                    (x0, veto_bottom),
                    0.25,
                    VETO_H,
                    facecolor=COL["teflon"],
                    edgecolor="#a98016",
                    lw=0.8,
                    zorder=4,
                )
            )
        for z_ring in np.linspace(veto_bottom + 1.6, -1.6, 6):
            for sign in (-1, 1):
                x0 = 8.50 if sign > 0 else -9.00
                ax.add_patch(
                    Rectangle(
                        (x0, z_ring - 0.22),
                        0.50,
                        0.44,
                        facecolor=COL["copper"],
                        edgecolor="#7a451d",
                        lw=0.5,
                        zorder=5,
                    )
                )

        # XML shift_z values are relative to the detector centre.
        electrodes = {
            "Gate": XD_ZC + 8.6925,
            "Cathode": XD_ZC - 6.6650,
        }
        for name, z_el in electrodes.items():
            ax.plot([-TPC_R, TPC_R], [z_el, z_el], color="#4b5563", lw=1.3, zorder=6)
            ax.text(
                -TPC_R + 0.35,
                z_el + (0.38 if name == "Gate" else -0.55),
                name,
                fontsize=8,
                color="#4b5563",
                ha="left",
                va="center",
                zorder=7,
            )


def _draw_gap_break(ax, x_break):
    """Mark the compressed source-detector distance used for source 2."""
    # Erase only short pieces of the three beam lines.  A full-height mask
    # would hide labels and the geometry legend.
    for z_mark in (SOURCE_Z - SOURCE_R, SOURCE_Z, SOURCE_Z + SOURCE_R):
        ax.plot(
            [x_break - 0.58, x_break + 0.58],
            [z_mark, z_mark],
            color=COL["world"],
            lw=4.2,
            solid_capstyle="butt",
            zorder=9,
        )
        ax.plot(
            [x_break - 0.38, x_break - 0.05],
            [z_mark - 0.24, z_mark + 0.24],
            color="#6b7280",
            lw=1.1,
            clip_on=False,
            zorder=11,
        )
        ax.plot(
            [x_break + 0.05, x_break + 0.38],
            [z_mark - 0.24, z_mark + 0.24],
            color="#6b7280",
            lw=1.1,
            clip_on=False,
            zorder=11,
        )


def draw_overview(ax, source_name):
    """Panel A: source, beam cone, and detector geometry."""
    source = SOURCES[source_name]
    source_x = source["display_x"]
    cone_deg = source["cone_deg"]
    ax.set_facecolor(COL["world"])
    draw_detector(ax, detailed=True)

    # Planar circular source appears as a vertical bar in this x-z projection.
    ax.add_patch(
        Rectangle(
            (source_x - 0.35, SOURCE_Z - SOURCE_R),
            0.70,
            2 * SOURCE_R,
            facecolor=COL["source"],
            edgecolor="#124c22",
            lw=1.5,
            zorder=7,
        )
    )
    ax.text(
        source_x - 0.3,
        1.35,
        f"planar $\\gamma$ source\n$({source['x']:.0f},\\,0,\\,-7.55)$ cm; $r=8$ cm",
        ha="right",
        va="bottom",
        fontsize=9.2,
        color="#165d27",
    )

    cone_x = -11.0
    if source_name == "source1":
        # Source 1 has a broad angular cone about -x.
        half_height = (source_x - cone_x) * np.tan(np.deg2rad(cone_deg))
        ray_starts = (SOURCE_Z, SOURCE_Z)
        ray_ends = (SOURCE_Z - half_height, SOURCE_Z + half_height)
    else:
        # Source 2 is an 8 cm radius plane emitting a near-parallel beam.
        # Show the band edges; the additional 0.1 degree divergence is tiny.
        divergence = (source_x - cone_x) * np.tan(np.deg2rad(cone_deg))
        ray_starts = (SOURCE_Z - SOURCE_R, SOURCE_Z + SOURCE_R)
        ray_ends = (SOURCE_Z - SOURCE_R - divergence, SOURCE_Z + SOURCE_R + divergence)

    for z_start, z_end in zip(ray_starts, ray_ends):
        ax.plot(
            [source_x, cone_x],
            [z_start, z_end],
            color=COL["source"],
            lw=1.3,
            ls=(0, (5, 4)),
            alpha=0.75,
            zorder=1,
        )
    ax.annotate(
        "",
        xy=(11.8, SOURCE_Z),
        xytext=(source_x - 2.5, SOURCE_Z),
        arrowprops=dict(arrowstyle="-|>", color=COL["source"], lw=2.4, mutation_scale=15),
        zorder=8,
    )
    ax.text(
        20.0,
        SOURCE_Z - 1.2,
        f"$1$--$8$ MeV $\\gamma$ beam\n$0$--${cone_deg:g}^\\circ$ about $-x$",
        ha="center",
        va="top",
        fontsize=9.4,
        color="#165d27",
    )

    if source_name == "source2":
        _draw_gap_break(ax, x_break=13.2)

    # A compact label avoids arrows crossing the detector internals.
    ax.text(
        0,
        -8.0,
        "active LXe TPC\n$r<8.25$ cm",
        ha="center",
        va="center",
        fontsize=10,
        weight="semibold",
        color="#194f86",
        zorder=9,
    )

    handles = [
        Patch(facecolor=COL["active"], edgecolor=COL["active_edge"], label="active LXe TPC"),
        Patch(facecolor=COL["veto"], edgecolor=COL["active_edge"], label="LXe veto shell"),
        Patch(facecolor=COL["teflon"], edgecolor="#a98016", label="side Teflon"),
        Patch(facecolor=COL["copper"], edgecolor="#7a451d", label="shaping rings"),
    ]
    ax.legend(
        handles=handles,
        loc="lower right",
        ncol=2,
        fontsize=8.2,
        frameon=True,
        framealpha=0.96,
        borderpad=0.7,
        columnspacing=1.2,
        handlelength=1.5,
    )

    ax.text(
        0.02,
        0.975,
        "World: G4_AIR \nClean geometry: no lead castle or PMT array",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        color=COL["muted"],
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#c8d0d8", alpha=0.95),
    )

    ax.set_xlim(-13.0, 34.0)
    ax.set_ylim(-20.5, 4.2)
    ax.set_aspect("equal", adjustable="box")
    if source_name == "source2":
        ax.set_xticks([-10, 0, 10, 13.2, 30])
        ax.set_xticklabels(["−10", "0", "10", "…", "200"])
        ax.set_xlabel("$x$ [cm]")
    else:
        ax.set_xlabel("$x$ [cm]")
    ax.set_ylabel("$z$ [cm]")
    ax.set_title("Beam and clean detector geometry", loc="left", weight="bold", pad=10)


def _angle_deg(vector):
    return np.degrees(np.arctan2(vector[1], vector[0])) % 360


def draw_event(ax):
    """Panel B: one connected two-site gamma event and the angle definitions."""
    ax.set_facecolor("#fbfcfd")
    draw_detector(ax, detailed=False)

    r0 = np.array([9.5, -6.0])
    r1 = np.array([3.2, -6.0])
    r2 = np.array([-2.0, -12.8])
    electron_end = np.array([6.0, -10.1])

    arrow_kw = dict(arrowstyle="-|>", color=COL["gamma"], lw=3.0, mutation_scale=17)
    ax.annotate("", xy=r1, xytext=r0, arrowprops=arrow_kw, zorder=9)
    ax.annotate("", xy=r2, xytext=r1, arrowprops=arrow_kw, zorder=9)

    # p_in 的虚线延长线
    unit_in = (r1 - r0) / np.linalg.norm(r1 - r0)

    # 从圆点边缘开始，避免虚线穿过第一散射点
    extension_start = r1 + 0.48 * unit_in
    extension_end = r1 + 4.5 * unit_in

    ax.plot(
        [extension_start[0], extension_end[0]],
        [extension_start[1], extension_end[1]],
        color=COL["gamma"],
        lw=1.8,
        ls=(0, (4, 3)),
        alpha=0.8,
        zorder=8,
    )

    ax.annotate(
        "",
        xy=electron_end,
        xytext=r1,
        arrowprops=dict(arrowstyle="-|>", color=COL["electron"], lw=2.2, mutation_scale=15),
        zorder=10,
    )

    # A blue dashed segment emphasizes the experimentally inferred direction.
    unit_out = (r2 - r1) / np.linalg.norm(r2 - r1)
    normal = np.array([-unit_out[1], unit_out[0]])
    offset = 0.22 * normal
    ax.plot(
        [r1[0] + offset[0], r2[0] + offset[0]],
        [r1[1] + offset[1], r2[1] + offset[1]],
        color=COL["reco"],
        lw=1.8,
        ls=(0, (4, 3)),
        zorder=10,
    )

    for centre, radius, color in [(r1, 0.4, COL["electron"]), (r2, 0.4, COL["cluster2"])]:
        ax.add_patch(Circle(centre, radius, facecolor=color, edgecolor="white", lw=1.5, zorder=12))

    vin = r1 - r0
    vout = r2 - r1
    a1 = _angle_deg(vin)
    a2 = _angle_deg(vout)
    if a2 < a1:
        a2 += 360
    ax.add_patch(
        Arc(
            r1,
            5.2,
            5.2,
            theta1=a1,
            theta2=a2,
            color="#3d4650",
            lw=2.0,
            zorder=13,
        )
    )
    mid = np.deg2rad((a1 + a2) / 2)
    theta_pos = r1 + 3.4 * np.array([np.cos(mid), np.sin(mid)])
    ax.text(
        theta_pos[0],
        theta_pos[1],
        r"$\theta$",
        fontsize=13,
        weight="bold",
        ha="center",
        va="center",
        color=COL["text"],
        zorder=14,
    )

    ax.text(6.7, -4.45, r"$\vec p_{\rm in}$", fontsize=11, color="#9a5c00", ha="center")
    ax.text(1.4, -10.0, r"$\vec p_{\rm out}$", fontsize=11, color="#9a5c00", ha="left")
    ax.text(6.2, -10.35, "recoil $e^-$", fontsize=9.5, color=COL["electron"], ha="left")

    ax.annotate(
        "first Compton\n    cluster",
        xy=r1,
        xytext=(-1.55, -2.25),
        fontsize=9.5,
        color=COL["electron"],
        weight="semibold",
        arrowprops=dict(arrowstyle="->", color=COL["electron"], lw=1.1),
        ha="left",
    )
    ax.annotate(
        "second cluster\n    Compton or\n    photoabsorption",
        xy=r2,
        xytext=(-7.85, -8.85),
        fontsize=8.9,
        color=COL["cluster2"],
        weight="semibold",
        arrowprops=dict(arrowstyle="->", color=COL["cluster2"], lw=1.1),
        ha="left",
    )

    ax.text(
        -7.7,
        -17.1,
        r"MC truth:  $\theta=\angle(\vec p_{\rm in},\vec p_{\rm out})$" "\n"
        r"Reconstruction:  $\theta=\angle(\vec p_{\rm in},\vec r_2-\vec r_1)$" "\n"
        r"Requires $\geq 2$ spatially resolved TPC energy clusters.",
        fontsize=9.3,
        ha="left",
        va="bottom",
        color=COL["text"],
        linespacing=1.35,
        bbox=dict(boxstyle="round,pad=0.55", fc="#eef5fc", ec=COL["reco"], lw=1.3),
        zorder=15,
    )

    guide = [
        Line2D([0], [0], color=COL["gamma"], lw=3, label="$\\gamma$ track"),
        Line2D([0], [0], color=COL["reco"], lw=1.8, ls=(0, (4, 3)), label="$\\vec r_2-\\vec r_1$"),
        Line2D([0], [0], color=COL["electron"], lw=2.2, label="recoil electron"),
    ]
    ax.legend(handles=guide, loc="upper left", fontsize=8.5, framealpha=0.96)

    # Equal padding around the true detector bounds keeps the frame balanced.
    frame_pad = 0.85
    ax.set_xlim(-XD_R - frame_pad, XD_R + frame_pad)
    ax.set_ylim(XD_ZC - XD_H / 2 - frame_pad, XD_ZC + XD_H / 2 + frame_pad)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("$x$ [cm]")
    ax.set_ylabel("$z$ [cm]")
    ax.set_title("Schematic of interactions", loc="left", weight="bold", pad=10)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        choices=tuple(SOURCES),
        default="source2",
        help="source configuration to draw (default: source2)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    source = SOURCES[args.source]
    out = Path(__file__).with_name(source["output"])
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "axes.edgecolor": "#4b5563",
            "axes.linewidth": 0.9,
            "xtick.color": "#4b5563",
            "ytick.color": "#4b5563",
        }
    )

    fig, (ax_overview, ax_event) = plt.subplots(
        1,
        2,
        figsize=(14.2, 7.3),
        gridspec_kw={"width_ratios": [1.34, 1.0]},
    )
    draw_overview(ax_overview, args.source)
    draw_event(ax_event)

    fig.suptitle(
        "RELICS reactor-$\\gamma$ interaction scenario",
        fontsize=18,
        weight="bold",
        color=COL["text"],
        y=0.985,
    )
    fig.text(
        0.5,
        0.945,
        f"R10  ·  $x$–$z$ side view",
        # f"R10  ·  {args.source}  ·  $x$–$z$ side view",
        ha="center",
        va="top",
        fontsize=10.5,
        color=COL["muted"],
    )
    # fig.text(
    #     0.5,
    #     0.018,
    #     "Trajectories and interaction positions in the right panel are schematic.",
    #     ha="center",
    #     va="bottom",
    #     fontsize=8.7,
    #     color=COL["muted"],
    # )

    fig.subplots_adjust(left=0.065, right=0.985, bottom=0.085, top=0.89, wspace=0.18)
    fig.savefig(out, dpi=220, bbox_inches="tight", facecolor="white")
    print(f"saved: {out}")


if __name__ == "__main__":
    main()

    
