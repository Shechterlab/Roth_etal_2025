#!/usr/bin/env python3
"""Compact, non-overlaid metaplots across PRO-seq-expressed histone genes."""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pyBigWig


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "histone_gene_stacked_metaplots"
OUT.mkdir(parents=True, exist_ok=True)

BED = (
    ROOT
    / "paper_github_upload/input_metadata/expressed_histone_genes/"
    "histone_body_HistoneDB20_PROseq_CONTROL_mean_ge100.bed"
)
ALL_BED = ROOT / "paper_github_upload/bed/histone_genes_HDB20.bed"

ARIAL = Path("/mnt/c/Windows/Fonts/arial.ttf")
ARIAL_BOLD = Path("/mnt/c/Windows/Fonts/arialbd.ttf")
for font_path in (ARIAL, ARIAL_BOLD):
    if font_path.exists():
        font_manager.fontManager.addfont(str(font_path))
plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 9,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

TSS_FLANK_BP = 400
N_TSS_BINS = 80
N_BOOT = 1000
RNG = np.random.default_rng(20260730)

CORE_TRACKS = [
    ("H3, DMSO 24 h", ROOT / "bw/averaged/H3_DMSO_24hr_mean.bigWig", "#3a3a3a", "H3/H4"),
    ("H4, DMSO 24 h", ROOT / "bw/averaged/H4_DMSO_24hr_mean.bigWig", "#3a3a3a", "H3/H4"),
    ("H4, GSK591 24 h", ROOT / "bw/averaged/H4_GSK591_24hr_mean.bigWig", "#07834d", "H3/H4"),
    ("HINFP (U2OS)", ROOT / "old/data/bw/HINFP_U2OS_REP1.mLb.clN.bigWig", "#3a3a3a", "U2OS"),
    ("NPAT (U2OS)", ROOT / "old/data/bw/NPAT_U2OS_REP1.mLb.clN.bigWig", "#3a3a3a", "U2OS"),
    ("IgG (U2OS)", ROOT / "old/data/bw/rabIgG_U2OS_REP1.mLb.clN.bigWig", "#3a3a3a", "U2OS"),
]

S5PH_TRACKS = [
    ("Pol II S5ph, DMSO 24 h", ROOT / "bw/averaged/PolIIS5ph_DMSO_24hr_mean.bigWig", "#3a3a3a", "S5ph"),
    ("Pol II S5ph, GSK591 24 h", ROOT / "bw/averaged/PolIIS5ph_GSK591_24hr_mean.bigWig", "#07834d", "S5ph"),
]


def read_genes():
    rows = []
    with BED.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 6:
                continue
            symbol = fields[3].split("|")[1] if "|" in fields[3] else fields[3]
            rows.append(
                {
                    "chrom": fields[0],
                    "start": int(fields[1]),
                    "end": int(fields[2]),
                    "symbol": symbol,
                    "strand": fields[5],
                }
            )
    return rows


def read_all_genes():
    rows = []
    with ALL_BED.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 6:
                continue
            rows.append(
                {
                    "chrom": fields[0],
                    "start": int(fields[1]),
                    "end": int(fields[2]),
                    "symbol": fields[3],
                    "strand": fields[5],
                }
            )
    return rows


def safe_stats(bw, chrom, start, end, bins):
    if start < 0 or end <= start:
        return np.zeros(bins)
    try:
        vals = bw.stats(chrom, start, end, nBins=bins, type="mean")
    except RuntimeError:
        vals = [0.0] * bins
    return np.nan_to_num(np.asarray(vals, dtype=float), nan=0.0)


def gene_profile(bw, gene):
    tss = gene["start"] if gene["strand"] == "+" else gene["end"]
    profile = safe_stats(
        bw,
        gene["chrom"],
        tss - TSS_FLANK_BP,
        tss + TSS_FLANK_BP,
        N_TSS_BINS,
    )
    if gene["strand"] == "+":
        return profile
    return profile[::-1]


def summarize(matrix):
    mean = matrix.mean(axis=0)
    boot = np.empty((N_BOOT, matrix.shape[1]), dtype=float)
    for i in range(N_BOOT):
        boot[i] = matrix[RNG.integers(0, matrix.shape[0], matrix.shape[0])].mean(axis=0)
    low, high = np.quantile(boot, [0.025, 0.975], axis=0)
    return mean, low, high


expressed_genes = read_genes()
expressed_symbols = {gene["symbol"] for gene in expressed_genes}
all_genes = read_all_genes()
nonexpressed_genes = [gene for gene in all_genes if gene["symbol"] not in expressed_symbols]
all_tracks = CORE_TRACKS + S5PH_TRACKS
matrices = {}
summaries = {}
for label, path, _, _ in all_tracks:
    with pyBigWig.open(str(path)) as bw:
        matrices[("expressed", label)] = np.vstack([gene_profile(bw, gene) for gene in expressed_genes])
        matrices[("nonexpressed", label)] = np.vstack([gene_profile(bw, gene) for gene in nonexpressed_genes])
    summaries[("expressed", label)] = summarize(matrices[("expressed", label)])
    summaries[("nonexpressed", label)] = summarize(matrices[("nonexpressed", label)])

GLOBAL_GROUP_YMAX = {}
for group in {row[3] for row in all_tracks}:
    highs = np.concatenate(
        [
            summaries[(gene_set, row[0])][0]
            for gene_set in ("expressed", "nonexpressed")
            for row in all_tracks
            if row[3] == group
        ]
    )
    GLOBAL_GROUP_YMAX[group] = max(0.25, float(np.quantile(highs, 0.995)) * 1.06)


def make_figure(tracks, basename, gene_set, title):
    n = len(tracks)
    n_genes = len(expressed_genes) if gene_set == "expressed" else len(nonexpressed_genes)
    group_ymax = GLOBAL_GROUP_YMAX
    fig, axes = plt.subplots(n, 1, figsize=(3.9, 0.58 * n + 0.76), sharex=True)
    if n == 1:
        axes = [axes]
    x = np.arange(N_TSS_BINS)
    groups_labeled = set()
    for ax, (label, _, color, group) in zip(axes, tracks):
        mean, low, high = summaries[(gene_set, label)]
        ax.plot(x, mean, color=color, linewidth=1.35)
        ax.set_ylim(0, group_ymax[group])
        ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=8.7, labelpad=17)
        ax.set_yticks([0, group_ymax[group]])
        ax.set_yticklabels(["0", f"{group_ymax[group]:.1f}"], fontsize=7.2)
        ax.tick_params(axis="y", left=True, labelleft=True, length=2.2, width=0.7, pad=2)
        if group not in groups_labeled:
            ax.text(
                -0.015, 1.13, "CPM", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7.2
            )
            groups_labeled.add(group)
        ax.axvline(N_TSS_BINS / 2, color="#777777", linewidth=0.65, linestyle=(0, (2, 2)))
        ax.spines[["top", "right", "bottom"]].set_visible(False)
        ax.tick_params(axis="x", bottom=False, labelbottom=False)
    axes[-1].spines["bottom"].set_visible(True)
    axes[-1].tick_params(axis="x", bottom=True, labelbottom=True, length=3)
    axes[-1].set_xticks([0, N_TSS_BINS / 2, N_TSS_BINS - 1])
    axes[-1].set_xticklabels([f"−{TSS_FLANK_BP}", "TSS", f"+{TSS_FLANK_BP}"], fontsize=8.2)
    axes[-1].set_xlabel("Distance from TSS (bp)", fontsize=8.5, labelpad=2)
    fig.suptitle(
        f"{title} (n={n_genes})",
        fontsize=10.5,
        fontweight="bold",
        y=0.992,
    )
    fig.tight_layout(pad=0.32, h_pad=0.01)
    base = OUT / basename
    fig.savefig(base.with_suffix(".png"), dpi=450, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


make_figure(
    CORE_TRACKS,
    "histone_genes_stacked_metaplot_core",
    "expressed",
    "Expressed histone genes",
)
make_figure(
    CORE_TRACKS[:3] + S5PH_TRACKS + CORE_TRACKS[3:],
    "histone_genes_stacked_metaplot_with_S5ph",
    "expressed",
    "Expressed histone genes",
)
make_figure(
    CORE_TRACKS,
    "nonexpressed_histone_genes_stacked_metaplot_core",
    "nonexpressed",
    "Non-expressed histone genes",
)

with (OUT / "histone_genes_metaplot_summary.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["gene_set", "track", "segment", "bin", "mean", "ci_low", "ci_high", "n_genes"])
    for gene_set, gene_rows in (("expressed", expressed_genes), ("nonexpressed", nonexpressed_genes)):
        for label, _, _, _ in all_tracks:
            mean, low, high = summaries[(gene_set, label)]
            for i, (m, lo, hi) in enumerate(zip(mean, low, high)):
                writer.writerow([gene_set, label, "TSS_pm400", i, m, lo, hi, len(gene_rows)])
