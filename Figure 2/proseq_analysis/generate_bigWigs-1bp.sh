#!/bin/bash
# David Shechter 2024-12-06
# USE as follows: sbatch --array=0-<N-1> generate_bigwig_array.sh --bam-dir /path/to/your/bam/files
#SBATCH -p normal                # Partition
#SBATCH -N 1                     # Number of nodes
#SBATCH -n 8                     # Number of CPUs per task
#SBATCH --time=04:00:00          # Runtime (adjust as needed)
#SBATCH --mem=16G                # Memory per task
#SBATCH -o logs/deeptools_%A_%a.out  # Standard output log (%A = job ID, %a = task ID)
#SBATCH -e logs/deeptools_%A_%a.err  # Standard error log
#SBATCH --array=0-99             # Array range (adjust based on number of BAM files)

# Load necessary modules
source /gs/gsfs0/hpc01/rhel8/apps/conda3/bin/activate
conda activate deeptools || { echo "Error: Could not activate deeptools environment."; exit 1; }

# Function to display usage
usage() {
    echo "Usage: $0 --bam-dir BAM_DIR --output-dir OUTPUT_DIR"
    echo "  --bam-dir      Directory containing BAM files"
    echo "  --output-dir   Directory to save output BigWig files"
    exit 1
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --bam-dir)
            BAM_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [[ -z "$BAM_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Both --bam-dir and --output-dir must be specified."
    usage
fi

# Create output and log directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Get list of BAM files
BAM_FILES=($(find "$BAM_DIR" -name "*.bam" | sort))
NUM_BAMS=${#BAM_FILES[@]}

# Ensure the task index is within bounds
if [[ $SLURM_ARRAY_TASK_ID -ge $NUM_BAMS ]]; then
    echo "Task ID ${SLURM_ARRAY_TASK_ID} exceeds number of BAM files (${NUM_BAMS}). Exiting."
    exit 1
fi

# Get the BAM file corresponding to this array task
BAM_FILE=${BAM_FILES[$SLURM_ARRAY_TASK_ID]}
BASE_NAME=$(basename "$BAM_FILE" .bam)

# Check if the BAM file has an associated .bai index file
if [[ ! -f "${BAM_FILE}.bai" ]]; then
    echo "Error: Index file (.bai) not found for ${BAM_FILE}. Skipping."
    exit 1
fi

echo "Processing $BAM_FILE..."

# Generate BigWig for plus strand
bamCoverage --bam "$BAM_FILE" \
            --outFileName "$OUTPUT_DIR/${BASE_NAME}.plus.bigWig" \
            --binSize 1 \
            --normalizeUsing CPM \
            --filterRNAstrand forward \
            --effectiveGenomeSize 2913022398 \
            --numberOfProcessors 8 || { echo "Error processing plus strand for $BAM_FILE"; exit 1; }

# Generate BigWig for minus strand
bamCoverage --bam "$BAM_FILE" \
            --outFileName "$OUTPUT_DIR/${BASE_NAME}.minus.bigWig" \
            --binSize 1 \
            --normalizeUsing CPM \
            --filterRNAstrand reverse \
            --effectiveGenomeSize 2913022398 \
            --numberOfProcessors 8 || { echo "Error processing minus strand for $BAM_FILE"; exit 1; }

# Combine plus and minus strand BigWig files into a single file
bigwigCompare --bigwig1 "$OUTPUT_DIR/${BASE_NAME}.plus.bigWig" \
              --bigwig2 "$OUTPUT_DIR/${BASE_NAME}.minus.bigWig" \
              --operation subtract \
              --outFileName "$OUTPUT_DIR/${BASE_NAME}.bothStrands.bigWig" \
              --binSize 1 \
              --numberOfProcessors 8 || { echo "Error combining strands for $BAM_FILE"; exit 1; }

echo "Finished processing $BAM_FILE."

# Deactivate the conda environment
conda deactivate
