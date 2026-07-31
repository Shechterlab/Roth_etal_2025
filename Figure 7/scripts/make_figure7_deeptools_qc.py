#!/usr/bin/env python3
"""Build optional PRO-seq-ranked deepTools heatmaps for Figure 7 QC.

The exact PRO-seq ordering is stored in the small rank manifest committed with
the paper code, so the full featureCounts matrix is not required to reproduce
the heatmap order.
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import shutil
import subprocess


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=root)
    parser.add_argument(
        "--averaged-bigwig-dir",
        type=Path,
        default=root.parent / "bw" / "averaged",
    )
    parser.add_argument(
        "--public-bigwig-dir",
        type=Path,
        default=root.parent / "old" / "data" / "bw",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
    )
    return parser.parse_args()


def require_executable(name: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        raise RuntimeError(f"{name} is required but was not found on PATH")
    return executable


def run(args: list[object], env: dict[str, str]) -> None:
    command = [str(value) for value in args]
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True, env=env)


def write_bed(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w") as handle:
        for row in rows:
            handle.write(
                f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['symbol']}"
                f"\t{row['control_mean']}\t{row['strand']}\n"
            )


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    averaged_dir = args.averaged_bigwig_dir.resolve()
    public_dir = args.public_bigwig_dir.resolve()
    out = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else root / "results" / "figure7_deeptools_qc"
    )
    bed_dir = out / "bed"
    matrix_dir = out / "matrices"
    figure_dir = out / "figures"
    cache_dir = out / "cache"
    mpl_dir = out / "mplconfig"
    for directory in (bed_dir, matrix_dir, figure_dir, cache_dir, mpl_dir):
        directory.mkdir(parents=True, exist_ok=True)

    compute_matrix = require_executable("computeMatrix")
    plot_heatmap = require_executable("plotHeatmap")
    pdftoppm = require_executable("pdftoppm")

    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str(mpl_dir)
    env["XDG_CACHE_HOME"] = str(cache_dir)

    manifest = root / "input_metadata" / "figure7" / "histone_genes_PROseq_rank_manifest.csv"
    if not manifest.exists():
        raise FileNotFoundError(manifest)
    with manifest.open() as handle:
        rows = list(csv.DictReader(handle))
    expressed = [row for row in rows if row["expression_group"].startswith("Expressed")]
    lower = [row for row in rows if row["expression_group"].startswith("Lower")]
    expressed.sort(key=lambda row: int(row["rank_within_group"]))
    lower.sort(key=lambda row: int(row["rank_within_group"]))

    expressed_bed = bed_dir / "histone_genes_expressed_PROseq_desc.bed"
    lower_bed = bed_dir / "histone_genes_lower_nonexpressed_PROseq_desc.bed"
    write_bed(expressed_bed, expressed)
    write_bed(lower_bed, lower)

    comparisons = {
        "H3_H4": {
            "tracks": [
                averaged_dir / "H3_DMSO_24hr_mean.bigWig",
                averaged_dir / "H4_DMSO_24hr_mean.bigWig",
                averaged_dir / "H4_GSK591_24hr_mean.bigWig",
            ],
            "labels": ["H3", "H4", "H4 + GSK591"],
            "colors": ["#555555", "#555555", "#07834d"],
        },
        "HLB_factors": {
            "tracks": [
                public_dir / "HINFP_U2OS_REP1.mLb.clN.bigWig",
                public_dir / "NPAT_U2OS_REP1.mLb.clN.bigWig",
                public_dir / "rabIgG_U2OS_REP1.mLb.clN.bigWig",
            ],
            "labels": ["HINFP", "NPAT", "IgG"],
            "colors": ["#555555", "#555555", "#9a9a9a"],
        },
    }

    for name, config in comparisons.items():
        tracks = [Path(path) for path in config["tracks"]]
        missing = [path for path in tracks if not path.exists()]
        if missing:
            raise FileNotFoundError("\n".join(str(path) for path in missing))
        matrix = matrix_dir / f"{name}_histone_TSS_pm400_PROseq_ranked.matrix.gz"
        sorted_bed = bed_dir / f"{name}_matrix_regions.bed"
        run(
            [
                compute_matrix,
                "reference-point",
                "--referencePoint",
                "TSS",
                "-b",
                "400",
                "-a",
                "400",
                "--binSize",
                "10",
                "-R",
                expressed_bed,
                lower_bed,
                "-S",
                *tracks,
                "--samplesLabel",
                *config["labels"],
                "--missingDataAsZero",
                "--sortRegions",
                "keep",
                "--outFileSortedRegions",
                sorted_bed,
                "-p",
                "max/2",
                "-o",
                matrix,
            ],
            env,
        )
        pdf = figure_dir / f"{name}_histone_TSS_pm400_profile_heatmap.pdf"
        run(
            [
                plot_heatmap,
                "-m",
                matrix,
                "-out",
                pdf,
                "--sortRegions",
                "no",
                "--averageTypeSummaryPlot",
                "mean",
                "--plotType",
                "lines",
                "--plotTitle",
                "PRO-seq-ranked histone genes",
                "--regionsLabel",
                f"Expressed (n={len(expressed)})",
                f"Lower/non-expressed (n={len(lower)})",
                "--samplesLabel",
                *config["labels"],
                "--colorList",
                *[f"white,{color}" for color in config["colors"]],
                "--refPointLabel",
                "TSS",
                "--yAxisLabel",
                "Mean signal",
                "--heatmapHeight",
                "10",
                "--heatmapWidth",
                "4.5",
                "--whatToShow",
                "plot and heatmap",
                "--legendLocation",
                "best",
            ],
            env,
        )
        run(
            [pdftoppm, "-singlefile", "-png", "-r", "220", pdf, figure_dir / pdf.stem],
            env,
        )

    with (out / "README.md").open("w") as handle:
        handle.write(
            "# Figure 7 deepTools QC\n\n"
            "These heatmaps use the same ±400 bp TSS window and the exact "
            "PRO-seq rank stored in `input_metadata/figure7/"
            "histone_genes_PROseq_rank_manifest.csv`. They are supporting QC "
            "for the compact Python metaplots, not a separate biological test.\n"
        )
    print(f"Wrote deepTools QC panels to {out}")


if __name__ == "__main__":
    main()
