#!/usr/bin/env python3

import argparse
import csv
import math
import os
import re

import pyBigWig


ANTIBODIES = ("H3", "H4", "PolIIS5ph")
TREATMENTS = ("DMSO", "GSK591", "LLY283")
TIMES = ("3hr", "12hr", "24hr", "48hr")
TRACK_DRUGS = ("GSK591", "LLY283")
TRACK_TIMES = ("24hr", "48hr")
PC = 1e-3
BIGWIG_DIR_OVERRIDE = None


def read_bed(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            rows.append({
                "chr": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "symbol": fields[3],
                "strand": fields[5] if len(fields) > 5 else ".",
            })
    return rows


def region_id(region_set, row, start=None, end=None):
    start = row["start"] if start is None else start
    end = row["end"] if end is None else end
    return f"{region_set}|{row['symbol']}|{row['chr']}:{start}-{end}:{row['strand']}"


def tss_window(row, flank=1000):
    tss = row["start"] if row["strand"] != "-" else row["end"]
    start = max(0, tss - flank)
    end = tss + flank
    return start, end


def read_histone_symbols(path):
    with open(path, newline="") as handle:
        return {row["HGNC symbol"] for row in csv.DictReader(handle)}


def bw_path(root, antibody, treatment, time):
    if BIGWIG_DIR_OVERRIDE is not None:
        bw_dir = BIGWIG_DIR_OVERRIDE
    else:
        local_dir = os.path.join(root, "bw", "averaged")
        parent_dir = os.path.join(os.path.dirname(root), "bw", "averaged")
        bw_dir = local_dir if os.path.isdir(local_dir) else parent_dir
    path = os.path.join(bw_dir, f"{antibody}_{treatment}_{time}_mean.bigWig")
    return path if os.path.exists(path) else None


def bw_stat(bw, chrom, start, end, stat_type):
    try:
        if stat_type == "mean":
            mean_value = bw.stats(chrom, start, end, type="mean")[0]
            coverage = bw.stats(chrom, start, end, type="coverage")[0]
            if mean_value is None or coverage is None:
                return 0.0
            value = mean_value * coverage
        else:
            value = bw.stats(chrom, start, end, type=stat_type)[0]
    except RuntimeError:
        value = None
    return float(value) if value is not None and math.isfinite(value) else 0.0


def write_signal_tables(root, allgenes, histone_symbols):
    stats_dir = os.path.join(root, "results", "nfcore_CPM_h4_h3_histone_gene_enrichment", "stats")
    os.makedirs(stats_dir, exist_ok=True)
    region_defs = {
        "allgenes_body": [(row, row["start"], row["end"]) for row in allgenes],
        "allgenes_TSS1000": [(row, *tss_window(row)) for row in allgenes],
    }
    for region_set, regions in region_defs.items():
        for antibody in ("H3", "H4"):
            path = bw_path(root, antibody, "DMSO", "24hr")
            if path is None:
                raise FileNotFoundError(f"Missing averaged bigWig for {antibody} DMSO 24hr")
            out_path = os.path.join(stats_dir, f"{region_set}_{antibody}_DMSO24h_nfcore_CPM_peaklike_signal.csv")
            bw = pyBigWig.open(path)
            try:
                with open(out_path, "w", newline="") as handle:
                    writer = csv.writer(handle)
                    writer.writerow(["region_set", "antibody", "symbol", "region_id", "chr", "start", "end", "strand", "is_expressed_histone", "mean_signal", "max_signal"])
                    for row, start, end in regions:
                        writer.writerow([
                            region_set,
                            antibody,
                            row["symbol"],
                            region_id(region_set, row, start, end),
                            row["chr"],
                            start,
                            end,
                            row["strand"],
                            str(row["symbol"] in histone_symbols).upper(),
                            bw_stat(bw, row["chr"], start, end, "mean"),
                            bw_stat(bw, row["chr"], start, end, "max"),
                        ])
            finally:
                bw.close()


def write_condition_mean(root, allgenes):
    out_dir = os.path.join(root, "results", "nfcore_CPM_body_occupancy_ratio", "stats")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "condition_mean_signal_by_gene.csv")
    bws = {}
    for antibody in ANTIBODIES:
        for treatment in TREATMENTS:
            for time in TIMES:
                path = bw_path(root, antibody, treatment, time)
                if path is not None:
                    bws[(antibody, treatment, time)] = pyBigWig.open(path)
    try:
        with open(out_path, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["region_set", "region_id", "symbol", "chr", "start", "end", "strand", "antibody", "treatment", "time", "mean_signal", "n_reps"])
            for row in allgenes:
                rid = region_id("allgenes_body", row)
                for (antibody, treatment, time), bw in bws.items():
                    value = bw_stat(bw, row["chr"], row["start"], row["end"], "mean")
                    if math.isfinite(value):
                        writer.writerow([
                            "allgenes_body",
                            rid,
                            row["symbol"],
                            row["chr"],
                            row["start"],
                            row["end"],
                            row["strand"],
                            antibody,
                            treatment,
                            time.replace("hr", "h"),
                            value,
                            1,
                        ])
    finally:
        for bw in bws.values():
            bw.close()


def write_h4_official_inputs(root, histone_rows):
    out_dir = os.path.join(root, "results", "nfcore_CPM_h4_official_figure", "data")
    os.makedirs(out_dir, exist_ok=True)
    rows = []
    for antibody in ("H4", "PolIIS5ph"):
        for drug in TRACK_DRUGS:
            for time in TRACK_TIMES:
                dmso_path = bw_path(root, antibody, "DMSO", time)
                drug_path = bw_path(root, antibody, drug, time)
                if dmso_path is None or drug_path is None:
                    continue
                dmso_bw = pyBigWig.open(dmso_path)
                drug_bw = pyBigWig.open(drug_path)
                try:
                    for row in histone_rows:
                        dmso = bw_stat(dmso_bw, row["chr"], row["start"], row["end"], "mean")
                        drug_signal = bw_stat(drug_bw, row["chr"], row["start"], row["end"], "mean")
                        logfc = math.log2((drug_signal + PC) / (dmso + PC)) if math.isfinite(dmso) and math.isfinite(drug_signal) else math.nan
                        rows.append({
                            "antibody": antibody,
                            "drug": drug,
                            "symbol": row["symbol"],
                            "region_id": region_id("histone_body", row),
                            "chr": row["chr"],
                            "start": row["start"],
                            "end": row["end"],
                            "strand": row["strand"],
                            "time": time.replace("hr", "h"),
                            "dmso_mean_signal": dmso,
                            "drug_mean_signal": drug_signal,
                            "log2FC_drug_vs_DMSO": logfc,
                        })
                finally:
                    dmso_bw.close()
                    drug_bw.close()
    out_path = os.path.join(out_dir, "nfcore_CPM_H4_PolIIS5ph_drug_24h_48h_expressed_histone_log2FC.csv")
    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    gsk_rows = [
        {
            "symbol": row["symbol"],
            "region_id": row["region_id"],
            "chr": row["chr"],
            "start": row["start"],
            "end": row["end"],
            "strand": row["strand"],
            "time": row["time"],
            "dmso_mean_signal": row["dmso_mean_signal"],
            "gsk591_mean_signal": row["drug_mean_signal"],
            "log2FC_GSK591_vs_DMSO": row["log2FC_drug_vs_DMSO"],
        }
        for row in rows
        if row["antibody"] == "H4" and row["drug"] == "GSK591"
    ]
    with open(os.path.join(out_dir, "nfcore_CPM_H4_GSK591_24h_48h_expressed_histone_log2FC.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(gsk_rows[0].keys()))
        writer.writeheader()
        writer.writerows(gsk_rows)


def binned_rows(chrom, start, end, step, samples):
    open_map = {label: pyBigWig.open(path) for label, path in samples.items()}
    try:
        for label, bw in open_map.items():
            pos = start
            while pos < end:
                b_end = min(end, pos + step)
                yield [chrom, pos, b_end, (pos + b_end) / 2.0, label, bw_stat(bw, chrom, pos, b_end, "mean")]
                pos = b_end
    finally:
        for bw in open_map.values():
            bw.close()


def write_track_tables(root):
    out_dir = os.path.join(root, "results", "nfcore_CPM_h4_official_figure", "data")
    for drug in TRACK_DRUGS:
        samples = {
            "H3 DMSO 24h": bw_path(root, "H3", "DMSO", "24hr"),
            "H4 DMSO 24h": bw_path(root, "H4", "DMSO", "24hr"),
            f"H4 {drug} 24h": bw_path(root, "H4", drug, "24hr"),
        }
        if any(path is None for path in samples.values()):
            continue
        for filename, start, end, step in [
            (f"hist1_chr6p_nfcore_CPM_tracks_{drug}.csv", 0, 58_000_000, 100_000),
            (f"hist1_chr6_zoom_nfcore_CPM_tracks_{drug}.csv", 26_000_000, 26_285_000, 250),
        ]:
            with open(os.path.join(out_dir, filename), "w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["chr", "start", "end", "mid", "track", "signal"])
                writer.writerows(binned_rows("chr6", start, end, step, samples))
    for src, dst in [
        ("hist1_chr6p_nfcore_CPM_tracks_GSK591.csv", "hist1_chr6p_nfcore_CPM_tracks.csv"),
        ("hist1_chr6_zoom_nfcore_CPM_tracks_GSK591.csv", "hist1_chr6_zoom_nfcore_CPM_tracks.csv"),
    ]:
        src_path = os.path.join(out_dir, src)
        if os.path.exists(src_path):
            with open(src_path) as inp, open(os.path.join(out_dir, dst), "w") as out:
                out.write(inp.read())


def update_readme(root):
    path = os.path.join(root, "results", "averaged_bigwig_source_README.md")
    with open(path, "w") as handle:
        handle.write(
            "# Averaged BigWig Source\n\n"
            "Paper-facing signal tables and figures in this folder were regenerated from `bw/averaged/*_mean.bigWig` files.\n\n"
            "PolIIS5ph averaged bigWigs were available for 12 h, 24 h, and 48 h. No PolIIS5ph 3 h averaged file was present, so PolIIS5ph 3 h is omitted from averaged-source time-course outputs.\n"
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root_dir")
    parser.add_argument(
        "--bigwig-dir",
        default=None,
        help="Optional directory containing averaged *_mean.bigWig inputs.",
    )
    parser.add_argument(
        "--figure7-only",
        action="store_true",
        help="Skip the unrelated full treatment/time-course gene-body table.",
    )
    args = parser.parse_args()
    root = args.root_dir
    global BIGWIG_DIR_OVERRIDE
    BIGWIG_DIR_OVERRIDE = os.path.abspath(args.bigwig_dir) if args.bigwig_dir else None
    allgenes = read_bed(os.path.join(root, "bed", "A549_expressed-mRNA-avgTPMgreaterThan10.bed"))
    expressed_histone_bed = os.path.join(root, "input_metadata", "expressed_histone_genes", "histone_body_HistoneDB20_PROseq_CONTROL_mean_ge100.bed")
    histone_rows = read_bed(expressed_histone_bed if os.path.exists(expressed_histone_bed) else os.path.join(root, "bed", "histone_genes_HDB20.bed"))
    histone_symbols = read_histone_symbols(os.path.join(root, "histone_genes", "human_histones.csv"))
    write_signal_tables(root, allgenes, histone_symbols)
    if not args.figure7_only:
        write_condition_mean(root, allgenes)
    write_h4_official_inputs(root, histone_rows)
    write_track_tables(root)
    update_readme(root)
    print(f"Wrote averaged-source signal tables under {root}/results")


if __name__ == "__main__":
    main()
