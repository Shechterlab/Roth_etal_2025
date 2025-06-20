#!/bin/bash
#SBATCH --job-name=histone_deeptools_scaleregions
#SBATCH --output=histone_deeptools_scaleregions_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --get-user-env

# Activate deeptools 3.5.6 conda environment
source /gs/gsfs0/hpc01/rhel8/apps/conda3/bin/activate /gs/gsfs0/users/dshecht1/.conda/envs/deeptools

# Working directory
cd /gs/gsfs0/home/dshecht1/analysis/proseq/PRMTi/deeptools

# Create output directories
mkdir -p deeptools_output-mat deeptools_output-pdf

# Parameters
THREADS=30
BINSIZE=10
UPSTREAM=3000
DOWNSTREAM=3000


# BED files (correct order)
BED_FILES=(
    A549_expressed-mRNA-avgTPMgreaterThan10.bed
)

# Region labels for plotting (same order as BED files)
BED_LABELS=(
    "All Expressed Genes (TPM > 10)"
)

#########################################################
# PROSEQ DATASET
#########################################################
#note misnaming of CONTROL_REP1 is actually REP2 bc of T2

# BigWig input files
PROSEQ_FILES=(
    2020-07_bw/CONTROL_REP1_T2.sorted.COMBINEDPOSITIVE.bigWig
    2020-07_bw/GSK591-15min_REP2_T1.sorted.COMBINEDPOSITIVE.bigWig
    2020-07_bw/GSK591-3hr_REP2_T1.sorted.COMBINEDPOSITIVE.bigWig
    2020-01_bw/CONTROL_REP2_T1.sorted.COMBINEDPOSITIVE.bigWig
    2020-01_bw/GSK591-2days_REP2_T1.sorted.COMBINEDPOSITIVE.bigWig
    2020-01_bw/GSK591-4days_REP2_T1.sorted.COMBINEDPOSITIVE.bigWig
    2020-01_bw/GSK591-7days_REP2_T1.sorted.COMBINEDPOSITIVE.bigWig
    2020-01_bw/MS023-7days_REP2_T1.sorted.COMBINEDPOSITIVE.bigWig
)

# Sample labels for plotting
PROSEQ_LABELS=(
    CONTROL-202007
    GSK591-15min
    GSK591-3hr
    CONTROL-202001
    GSK591-2days
    GSK591-4days
    GSK591-7days
    MS023-7days
)

# ComputeMatrix
computeMatrix scale-regions \
    -R "${BED_FILES[@]}" \
    -S "${PROSEQ_FILES[@]}" \
    -b ${UPSTREAM} -a ${DOWNSTREAM} \
    --skipZeros \
    --binSize ${BINSIZE} \
    --numberOfProcessors ${THREADS} \
    --outFileName deeptools_output-mat/matrix_proseq_scaled-allExpressed.gz \
    --outFileSortedRegions deeptools_output-mat/sorted_regions_proseq_scaled-allExpressed.bed \
    --missingDataAsZero \
    --verbose

# PlotHeatmap with correct regionsLabel
plotHeatmap \
    -m deeptools_output-mat/matrix_proseq_scaled-allExpressed.gz \
    -out deeptools_output-pdf/heatmap_proseq_scaled-allExpressed.pdf \
    --colorMap Blues \
    --samplesLabel "${PROSEQ_LABELS[@]}" \
    --regionsLabel "${BED_LABELS[@]}" \
    --sortUsing mean \
    --sortRegions descend \
    --heatmapHeight 5 \
    --heatmapWidth 3 \
    --zMin 0

#########################################################
# CUT&RUN DATASET
#########################################################

CUTRUN_FILES=(
    03_bed_to_bigwig/D0_K4me3_R1.bigWig
    03_bed_to_bigwig/D4_K4me3_R1.bigWig
    03_bed_to_bigwig/D0_K27me3_R1.bigWig
    03_bed_to_bigwig/D4_K27me3_R1.bigWig
)

CUTRUN_LABELS=(
    D0-K4me3
    D4-K4me3
    D0-K27me3
    D4-K27me3
)

# ComputeMatrix (NO regionsLabel)
computeMatrix reference-point \
    --referencePoint TSS \
    -R "${BED_FILES[@]}" \
    -S "${CUTRUN_FILES[@]}" \
	-b ${UPSTREAM} -a ${DOWNSTREAM} \
    --skipZeros \
    --binSize ${BINSIZE} \
    --numberOfProcessors ${THREADS} \
    --outFileName deeptools_output-mat/matrix_cutrun_scaled-allExpressed.gz \
    --outFileSortedRegions deeptools_output-mat/sorted_regions_cutrun_scaled-allExpressed.bed \
    --missingDataAsZero \
    --verbose

# PlotHeatmap
plotHeatmap \
    -m deeptools_output-mat/matrix_cutrun_scaled-allExpressed.gz \
    -out deeptools_output-pdf/heatmap_cutrun_scaled-allExpressed.pdf \
    --colorMap Greens \
    --samplesLabel "${CUTRUN_LABELS[@]}" \
    --regionsLabel "${BED_LABELS[@]}" \
    --sortUsing mean \
    --sortRegions descend \
    --heatmapHeight 5 \
    --heatmapWidth 3 \
    --zMin 0

echo "=== ALL PROCESSING COMPLETE ==="
