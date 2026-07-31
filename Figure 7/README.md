# PRMT5 H3/H4 CUT&Tag manuscript figures

This contains the final figure files, the numerical tables needed to inspect the plotted values, and unchanged copies of the scripts used to generate them.

## Included panels

| Manuscript panel | Files | Description |
|---|---|---|
| Figure 7a | `figures/Figure7a_H3C1_H4C1_browser.*` | H3C1–H4C1 locus browser with NPAT, H3 and H4 CUT&Tag, and control/GSK591 PRO-seq |
| Figure 7b | `figures/Figure7b_H3_H4_CUTTag.*` | H3/H4 CUT&Tag signal around histone-gene transcription start sites |
| Figure 7c | `figures/Figure7c_histone_gene_metaplot.*` | Stacked metaplot for the 56 expressed histone genes |
| Supplemental Figure 7a | `figures/FigureS7a_H4_occupancy_and_log2FC.*` | H4 CUT&Tag occupancy and expressed-histone-gene log2 fold change after GSK591 |

Each figure is supplied as PDF and high-resolution PNG; SVG is also included where produced by the original plotting workflow.

## Data

- `data/H3C1_H4C1_locus_binned_signal.csv`: binned values for the Figure 7a locus tracks.
- `data/histone_genes_metaplot_summary.csv`: binned TSS-profile summary used for Figure 7c. This curated copy contains no Pol II S5ph rows.
- `data/H4_GSK591_24h_48h_expressed_histone_log2FC.csv`: H4-only values underlying the log2 fold-change portion of Supplemental Figure 7a.
- `data/H4_GSK591_24h_48h_summary.csv`: H4-only summary statistics.
- `data/histone_genes_expressed_PROseq_desc.bed`: expressed histone-gene regions used for the CUT&Tag TSS panel.

## Scripts

The scripts in `scripts/` are copied unchanged from the analysis workspace. Source bigWigs are not included because of their size. The committed figure and data files allow the displayed results to be inspected without those source tracks.

## Exclusions

Pol II S5ph data, Pol II S5ph figures and statistics, lower/non-expressed-gene exploratory panels, deepTools QC panels not used in the paper, caches, temporary plotting files, and legacy results are intentionally omitted.
