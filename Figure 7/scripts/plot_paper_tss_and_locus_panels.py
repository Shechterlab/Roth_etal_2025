#!/usr/bin/env python3
"""Paper-oriented side-by-side TSS profiles and a focused H3C1/H4C1 panel."""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
import numpy as np
import pyBigWig


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "paper_tss_and_locus_panels"
OUT.mkdir(parents=True, exist_ok=True)

ARIAL = Path("/mnt/c/Windows/Fonts/arial.ttf")
ARIAL_BOLD = Path("/mnt/c/Windows/Fonts/arialbd.ttf")
for font_path in (ARIAL, ARIAL_BOLD):
    if font_path.exists():
        font_manager.fontManager.addfont(str(font_path))
plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8,
        "axes.linewidth": 0.7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

GRAY = "#353535"
LIGHT_GRAY = "#a0a0a0"
GREEN = "#07834d"
LOCUS_GRAY = "#555555"


def save_all(fig, base):
    fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight", pad_inches=0.025)
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.025)
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.025)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Panel 1: expressed and lower/non-expressed profiles side by side
# ---------------------------------------------------------------------------

summary_path = (
    ROOT / "results/histone_gene_stacked_metaplots/histone_genes_metaplot_summary.csv"
)
profiles = {}
with summary_path.open() as handle:
    for row in csv.DictReader(handle):
        key = (row["gene_set"], row["track"])
        profiles.setdefault(key, []).append((int(row["bin"]), float(row["mean"])))
profiles = {
    key: np.array([value for _, value in sorted(rows)])
    for key, rows in profiles.items()
}

rows = [
    ("H3, DMSO 24 h", "H3, DMSO 24 h", GRAY, "H3/H4"),
    ("H4, DMSO 24 h", "H4, DMSO 24 h", GRAY, "H3/H4"),
    ("H4, GSK591 24 h", "H4, GSK591 24 h", GREEN, "H3/H4"),
    ("HINFP (U2OS)", "HINFP (U2OS)", GRAY, "U2OS"),
    ("NPAT (U2OS)", "NPAT (U2OS)", GRAY, "U2OS"),
    ("IgG (U2OS)", "IgG (U2OS)", LIGHT_GRAY, "U2OS"),
]

group_max = {}
for group in {row[3] for row in rows}:
    values = np.concatenate(
        [
            profiles[(gene_set, key)]
            for gene_set in ("expressed", "nonexpressed")
            for _, key, _, row_group in rows
            if row_group == group
        ]
    )
    group_max[group] = float(np.ceil(values.max() * 1.08 * 2) / 2)

fig, axes = plt.subplots(
    len(rows),
    2,
    figsize=(5.25, 5.45),
    sharex=False,
    gridspec_kw={"hspace": 0.76, "wspace": 0.17},
)
column_titles = ["Expressed\n(n = 56)", "Lower/non-expressed\n(n = 60)"]
x = np.arange(80)

for col, (gene_set, title) in enumerate(
    zip(("expressed", "nonexpressed"), column_titles)
):
    axes[0, col].set_title(title, fontsize=10, fontweight="bold", pad=6)
    for row_idx, (display, key, color, group) in enumerate(rows):
        ax = axes[row_idx, col]
        y = profiles[(gene_set, key)]
        ax.plot(x, y, color=color, linewidth=1.55, solid_capstyle="round")
        ax.axvline(40, color="#8a8a8a", linewidth=0.6, linestyle=(0, (2, 2)))
        ax.set_ylim(0, group_max[group])
        ax.set_xlim(0, 79)
        ax.set_xticks([0, 40, 79])
        ax.set_xticklabels(["−400", "TSS", "+400"], fontsize=7.1)
        ax.tick_params(axis="x", bottom=True, labelbottom=True, length=2.2, pad=1.5)
        ax.spines[["top", "right"]].set_visible(False)
        if col == 0:
            ax.set_yticks([0, group_max[group]])
            ax.set_yticklabels(["0", f"{group_max[group]:g}"], fontsize=7.1)
            ax.tick_params(axis="y", left=True, length=2.2, pad=1.5)
            ax.set_ylabel(display, rotation=0, ha="right", va="center", fontsize=8.4, labelpad=20)
        else:
            ax.spines["left"].set_visible(False)
            ax.tick_params(axis="y", left=False, labelleft=False)
        if row_idx == len(rows) - 1:
            ax.set_xlabel("Distance from TSS (bp)", fontsize=7.8, labelpad=2)

# Group-scale labels and a subtle divider between A549 and public U2OS tracks.
axes[0, 0].text(
    -0.27, 1.02, "CPM", transform=axes[0, 0].transAxes,
    ha="left", va="bottom", fontsize=7.2
)
axes[3, 0].text(
    -0.27, 1.02, "a.u.", transform=axes[3, 0].transAxes,
    ha="left", va="bottom", fontsize=7.2
)
fig.subplots_adjust(left=0.37, right=0.995, top=0.955, bottom=0.075)
save_all(fig, OUT / "expressed_vs_nonexpressed_TSS_profiles_side_by_side")


# ---------------------------------------------------------------------------
# Panel 2: focused H3C1/H4C1 locus, no treatment comparison
# ---------------------------------------------------------------------------

CHROM = "chr6"
START = 26_019_400
END = 26_022_650
BIN = 5
X = np.arange(START, END, BIN)

def bw_values(path):
    with pyBigWig.open(str(path)) as bw:
        values = bw.stats(CHROM, START, END, nBins=len(X), type="mean")
    return np.nan_to_num(np.asarray(values, dtype=float), nan=0.0)


core_locus_tracks = [
    ("NPAT (U2OS)", ROOT / "old/data/bw/NPAT_U2OS_REP1.mLb.clN.bigWig", LOCUS_GRAY, "NPAT"),
    ("H3, DMSO 24 h", ROOT / "bw/averaged/H3_DMSO_24hr_mean.bigWig", LOCUS_GRAY, "H3/H4"),
    ("H4, DMSO 24 h", ROOT / "bw/averaged/H4_DMSO_24hr_mean.bigWig", LOCUS_GRAY, "H3/H4"),
    ("PRO-seq, control", ROOT / "bw/CONTROL_REP1_T1.sorted.bothStrands.bigWig", LOCUS_GRAY, "PROseq"),
]

treatment_locus_tracks = core_locus_tracks + [
    (
        "PRO-seq, GSK591 2 d",
        ROOT / "bw/GSK591-2days_REP1_T1.sorted.bothStrands.bigWig",
        GREEN,
        "PROseq",
    )
]


def render_locus(locus_tracks, basename):
    locus_values = {label: bw_values(path) for label, path, _, _ in locus_tracks}
    locus_max = {}
    for group in {row[3] for row in locus_tracks}:
        vals = np.concatenate(
            [locus_values[label] for label, _, _, row_group in locus_tracks if row_group == group]
        )
        locus_max[group] = float(np.ceil(np.quantile(vals, 0.995) * 1.08 * 2) / 2)

    fig = plt.figure(figsize=(5.25, 0.72 * len(locus_tracks) + 1.15))
    gs = fig.add_gridspec(
        len(locus_tracks) + 1,
        1,
        height_ratios=[1] * len(locus_tracks) + [0.60],
        hspace=0.33,
    )
    locus_axes = []
    for idx, (label, _, color, group) in enumerate(locus_tracks):
        ax = fig.add_subplot(gs[idx, 0], sharex=locus_axes[0] if locus_axes else None)
        locus_axes.append(ax)
        y = locus_values[label]
        min_y = min(0, float(np.quantile(y, 0.005)) * 1.05)
        ax.fill_between(X, 0, y, color=color, alpha=1, linewidth=0)
        ax.plot(X, y, color=color, linewidth=1.05)
        ax.axhline(0, color="#777777", linewidth=0.35)
        ax.set_ylim(min_y, locus_max[group])
        ax.set_yticks([0, locus_max[group]])
        ax.set_yticklabels(["0", f"{locus_max[group]:g}"], fontsize=7)
        ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=8.4, labelpad=18)
        ax.spines[["top", "right", "bottom"]].set_visible(False)
        ax.tick_params(axis="x", bottom=False, labelbottom=False)
        ax.tick_params(axis="y", length=2, pad=1.5)

    gene_ax = fig.add_subplot(gs[-1, 0], sharex=locus_axes[0])
    gene_ax.set_ylim(0, 1)
    gene_ax.spines[["top", "right", "left"]].set_visible(False)
    gene_ax.tick_params(axis="y", left=False, labelleft=False)
    gene_ax.set_xticks([26_020_000, 26_021_000, 26_022_000])
    gene_ax.set_xticklabels(["26.020", "26.021", "26.022"], fontsize=7.5)
    gene_ax.set_xlabel("hg38 chr6 position (Mb)", fontsize=8, labelpad=2)

    gene_models = [
        (26_020_451, 26_020_958, "H3C1"),
        (26_021_649, 26_022_050, "H4C1"),
    ]
    for start, end, name in gene_models:
        y = 0.55
        height = 0.18
        gene_ax.add_patch(
            Rectangle(
                (start, y),
                end - start,
                height,
                facecolor="#9a9a9a",
                edgecolor="none",
            )
        )
        for chevron_x in np.arange(start + 80, end - 20, 105):
            gene_ax.text(
                chevron_x,
                y + height / 2,
                "›",
                ha="center",
                va="center",
                fontsize=7.2,
                fontweight="bold",
                color="white",
            )
        gene_ax.text(
            (start + end) / 2,
            y - 0.06,
            name,
            ha="center",
            va="top",
            fontsize=8.5,
        )

    fig.suptitle("H3C1–H4C1 histone-gene locus", fontsize=11, fontweight="bold", y=0.99)
    fig.subplots_adjust(left=0.30, right=0.995, top=0.91, bottom=0.12)
    save_all(fig, OUT / basename)


render_locus(
    treatment_locus_tracks,
    "H3_H4_NPAT_PROseq_H3C1_H4C1_locus",
)
render_locus(
    core_locus_tracks,
    "H3_H4_NPAT_PROseq_control_only_H3C1_H4C1_locus",
)
render_locus(
    treatment_locus_tracks,
    "H3_H4_NPAT_PROseq_control_GSK591_H3C1_H4C1_locus",
)
