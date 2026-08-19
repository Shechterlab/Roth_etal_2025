# PRMT5 activity sustains histone production to maintain genome integrity

This repository contains the analysis scripts, FIJI/ImageJ macros, processed data, and selected figure outputs associated with:

> Roth JS, DeAngelo JD, Young DL, Maron MI, et al. **PRMT5 activity sustains histone production to maintain genome integrity.** *Nature Communications* (provisionally accepted).

The study uses time-resolved PRO-seq, quantitative proteomics, CUT&RUN, CUT&Tag, and fluorescence imaging to show that PRMT5 activity supports replication-dependent histone transcription and histone protein production. Disrupting PRMT5 reduces histone supply and produces replication-associated nuclear abnormalities. The manuscript further identifies enrichment of histone H4 and H4R3me2s at histone locus bodies.

## Repository organization

| Location | Contents |
|---|---|
| `Figure 1/` | FIJI macros for image preparation, nuclear/NPAT quantification, and line-profile analysis used for cellular morphology and immunofluorescence panels |
| `Figure 2/` | PRO-seq processing and analysis, differential expression, TPM calculation, heatmaps, permutation testing, chromosome plots, deepTools metagene plots, histone-gene annotations, and selected outputs |
| `Figure 3/` | FIJI macros for image preparation and NPAT foci quantification following HINFP knockdown or PRMT5 inhibition |
| `Figure 4/` | FIJI macros for image preparation and NPAT foci quantification in p53-knockdown experiments |
| `Figure 6/` | SILAC histone-PTM analysis, subcellular-fraction TMT proteomics, FIJI image-processing macros, and fluorescence line-profile analysis |
| `Figure 7/` | H3/H4 CUT&Tag quantification and plotting, deepTools QC, figure-ready data and outputs, and FIJI macros for H4/H4R3me2s imaging analyses |
| `scripts/proseq_analysis/` | Archival copies of the core PRO-seq count-matrix and DESeq2 helper scripts also present under `Figure 2/proseq_analysis/` |

There is no figure-specific computational directory for Figure 5; that figure is based primarily on biochemical, SILAC, and flow-cytometry experiments described in the manuscript. Relevant SILAC analysis files are stored under `Figure 6/SILAC/`.

## Script guide

### Figures 1, 3, 4, 6, and 7: fluorescence imaging

The figure directories reuse several FIJI/ImageJ macros:

- `ImageProcessingMacrosForDisplay/Macro1_Level-Images.txt` through `Macro5_SubsetCrop.txt` form a sequential workflow for selecting display levels, applying the same settings across an experiment, saving individual channels, applying lookup-table colors, and generating consistent crops.
- `ImmunofluorescenceQuantifications/Macro_NPAT-CountPerNuclei.txt` segments nuclei and counts NPAT foci per nucleus.
- `ImmunofluorescenceQuantifications/Macro_QuantifyNPAT.txt` measures NPAT foci properties, including count, area, and intensity.
- `Figure 7/ImmunofluorescenceQuantifications/Macro_Quantify-RegionOfInterestInNPAT-Puncta.txt` measures fluorescence signal within NPAT-positive regions for the H4 and H4R3me2s analyses.
- `ProfilePlotIntensityFromImages/Macro_ProfilePlot.txt` exports line-profile measurements from microscopy images. `ProfilePlotIntensityFromImages/ProfilePlot.Rmd` combines the exported CSV files and creates normalized channel-intensity plots.

These macros are archived from the original analysis environment. Before running them, review the expected channel order, file naming, pixel scale, thresholds, and input/output directories. Some macros call the FIJI command `Fresh Start` and therefore assume the same FIJI environment used for the study.

### Figure 2: PRO-seq

#### Count generation and differential analysis

- `proseq_analysis/20250608_make_hg38_refseq_NOtss_gtf.R` creates hg38 transcription-start-site intervals and prepares a TSS-excluded annotation used for gene-body counting. The documented workflow also uses `bedtools subtract`.
- `proseq_analysis/PROseq_deseq2_genebody.R` accepts a sample sheet, results directory, and GTF; runs `Rsubread::featureCounts`; filters low-count genes; performs DESeq2 analysis; and produces sample-distance, PCA, volcano, and comparison outputs.
- `proseq_analysis/merge_featureCounts_tables.R` merges individual featureCounts tables into one count matrix.
- `proseq_analysis/addGeneName.R` maps Ensembl identifiers to human gene symbols with `org.Hs.eg.db`.
- `proseq_analysis/proSeq_DESeq2.R` is a compact DESeq2 template for a user-supplied count matrix and sample metadata.
- `proseq_analysis/proseq_merged_DESeq2.R` runs pairwise treatment-versus-control DESeq2 comparisons from the merged count matrix.
- `proseq_analysis/20250608_PROseq_TPM_calculation.R` counts reads over transcript features, calculates transcript-length-normalized TPM values, filters low-expression genes, and writes the TPM matrix used by downstream analyses.
- `proseq_analysis/GSK591_PROseq_DE-comparisons.Rmd` performs the integrated time-course analysis, including p-value diagnostics, volcano plots, overlaps, pairwise comparisons, TPM heatmaps, PCA, and gene-ontology enrichment. `comparisons.txt` identifies the DESeq2 result files to load.
- `proseq_analysis/ENSG_to_bed.R` queries Ensembl BioMart and converts Ensembl gene identifiers to strand-aware BED6 records.

The `PROseq_DESeq2/` directory contains selected differential-expression tables used by the notebooks. The top-level `scripts/proseq_analysis/` directory contains duplicate archival copies of several core helpers; use the copies under `Figure 2/proseq_analysis/` when following the figure-organized layout.

#### Figure-specific visualization

- `Heatmap/DataViz_JSRe0230_AnalysisOfPROSeqNewDeseq2_20250619.Rmd` assembles early and extended GSK591 time-course results and produces the PRO-seq heatmap used in Figure 2.
- `RandomPermutationTesting/DataViz_JSRe0115_PROSeq-RandomPermutationTest.Rmd` compares histone-gene log2 fold changes with the rest of the transcriptome using random permutations across treatment times.
- `KaryoploteR/KaryoploteR.Rmd` maps regulated histone genes to the major chromosome 1 and chromosome 6 histone clusters.
- `deeptools/deeptools_histoneGenes_Proseq9.sh` creates scale-regions PRO-seq and reference-point CUT&RUN matrices/heatmaps for histone-gene classes.
- `deeptools/deeptools_allExpressedGenes_Proseq.sh` performs the analogous analysis for expressed protein-coding genes.
- `proseq_analysis/generate_bigWigs-1bp.sh` creates strand-specific, CPM-normalized, 1-bp bigWig tracks from BAM files with `bamCoverage`.

The `histone_genes/` and `deeptools/` directories include the BED/GTF gene sets needed by these analyses. Large source BAM and bigWig files are not stored in this repository.

### Figure 6: quantitative proteomics

- `SILAC/SILAC_HistonePTM_Analysis_PRMT5i.Rmd` analyzes heavy/light histone peptide measurements after release into S phase with or without GSK591. It calculates fraction heavy, PTM occupancy, and normalized abundance and exports summary tables, heatmaps, and stacked bar plots. The source workbooks, exported CSV files, plots, and rendered HTML notebook are included.
- `TMT_Proteome/TMT_analysis_DEP_byFraction_GSK591_geometricmeanInternalNorm_SpikeIn.Rmd` analyzes cytoplasmic/hypotonic, nucleoplasmic, and chromatin TMT fractions. It applies internal-reference scaling, constructs DEP objects, evaluates differential abundance, and produces PCA, scatter, violin/permutation, time-course, and volcano plots. The experimental design, input workbook, selected normalized tables, rendered HTML notebook, and plot outputs are included.

### Figure 7: H3/H4 CUT&Tag

See `Figure 7/README.md` for a panel-level inventory. The principal scripts are:

- `scripts/quantify_averaged_bigwigs_for_paper.py` quantifies averaged bigWig signal over gene sets and the H3C1-H4C1 locus and writes figure-ready tables.
- `scripts/plot_nfcore_CPM_h4_official_figure.R` creates the official H4 CUT&Tag occupancy and fold-change plots from the quantified tables.
- `scripts/plot_histone_gene_stacked_metaplots.py` generates stacked H3/H4 profiles across expressed histone genes.
- `scripts/plot_paper_tss_and_locus_panels.py` generates transcription-start-site profiles and the focused H3C1-H4C1 locus panel.
- `scripts/make_figure7_deeptools_qc.py` creates optional PRO-seq-ranked deepTools heatmaps for quality control.

The repository includes the numerical CSV tables and final PDF/PNG/SVG panels needed to inspect the plotted values. Source bigWigs are omitted because of their size.

## Data availability

Raw sequencing data are available from NCBI GEO:

- PRO-seq: [GSE275217](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275217) and [GSE275220](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275220)
- CUT&RUN: [GSE301721](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301721)
- CUT&Tag: [GSE334611](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334611)

Raw total-cell and fractionated LC-MS/MS data are available from ProteomeXchange/PRIDE under [PXD065294](https://www.ebi.ac.uk/pride/archive/projects/PXD065294).

Processed tables required for selected analyses are included alongside their scripts. Large alignment and signal-track files must be downloaded from the public repositories or regenerated from the raw data.

## Software requirements

The manuscript reports the following principal software versions: FIJI/ImageJ2 2.14.0/1.54f, R 4.4.0, PEPPRO 1.0, deepTools 3.5.1, IGV 2.14.1, nf-core/cutandrun 3.2.1-3.2.2, and nf-core/chipseq 2.1.0. The archived Slurm scripts activate a deepTools 3.5.6 environment.

Primary dependencies include:

- **R/Bioconductor:** DESeq2, Rsubread, data.table, tidyverse, biomaRt, org.Hs.eg.db, rtracklayer, GenomicRanges, karyoploteR, DEP, SummarizedExperiment, ComplexHeatmap, clusterProfiler, pheatmap, ggplot2, and readxl.
- **Python:** NumPy, Matplotlib, and pyBigWig.
- **Command line:** deepTools (`bamCoverage`, `computeMatrix`, and `plotHeatmap`), bedtools, Bowtie2, and a Unix-like shell; the supplied deepTools scripts contain Slurm directives.
- **Imaging:** FIJI/ImageJ2 with support for the commands used by the archived macros.

This repository does not provide a single locked software environment. Consult the `library(...)` and `import` statements in each script for its complete dependency list.

## Running the analyses

Clone the repository:

```bash
git clone https://github.com/Shechterlab/Roth_etal_2025.git
cd Roth_etal_2025
```

Then select the figure-specific workflow and inspect its script before execution. In particular:

1. Replace workstation- or cluster-specific paths with paths on your system.
2. Confirm sample names, comparisons, treatment labels, genome annotations, and image channel order.
3. Download or regenerate any raw BAM/bigWig files that are not included.
4. Create the required R, Python, FIJI, or deepTools environment.
5. Run notebooks from the directory containing their inputs, or update their input/output parameters first.

Many files are direct archival copies of the scripts used during analysis and contain absolute paths from the original workstations or HPC cluster. They are provided to document the published analysis rather than as a portable, one-command workflow. Rendered HTML notebooks, processed tables, and final figure files are included where available to aid inspection.

## Citation

Please cite the accepted article when the final Nature Communications citation and DOI become available. Until then, cite the manuscript as:

> Roth JS, DeAngelo JD, Young DL, Maron MI, et al. PRMT5 activity sustains histone production to maintain genome integrity. *Nature Communications*. Provisionally accepted.

## Contact and reuse

For questions about the study or repository, contact the corresponding author, David Shechter, at `david.shechter @ einsteinmed.edu`.
