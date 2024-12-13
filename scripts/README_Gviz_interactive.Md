# Gviz Dynamic Visualization Workflow

## Overview

This repository contains an R Markdown script for dynamically generating genomic plots using **Gviz**. It processes a CSV file as input to define track data and produces a publication-quality plot. The workflow supports **DataTrack**, **AnnotationTrack**, and **BiomartGeneRegionTrack**, with customization options for color, stacking, height, and genomic regions.

### Features
- **Dynamic Input**: Tracks and plot parameters are read from a CSV file for easy configuration.
- **Customizable Outputs**: Includes options for track-specific appearance (e.g., color, stacking, height).
- **Flexible Execution**: Designed to run in RStudio or on High-Performance Computing (HPC) systems.
- **Reproducibility**: All outputs are saved to a user-specified directory outside of the Git repository.

---

## Requirements

### Software
- R (≥ 4.0)
- Required R packages:
  - `Gviz`
  - `tidyverse`
  - `fs`

### Input CSV File
The CSV file should include the following columns:

| Column Name         | Description                                                                                             | Required | Default (if missing)          |
|---------------------|---------------------------------------------------------------------------------------------------------|----------|--------------------------------|
| `file_path`         | Path to the data file (e.g., BAM, BigWig, BED).                                                        | Yes      | N/A                            |
| `track_name`        | A user-friendly name for the track to display.                                                         | Yes      | N/A                            |
| `type`              | Track type (`DataTrack`, `AnnotationTrack`, or `BiomartGeneRegionTrack`).                              | Yes      | N/A                            |
| `genome`            | Genome assembly (e.g., `hg19`, `dm6`).                                                                 | Yes      | N/A                            |
| `chromosome`        | Chromosome for the genomic region (required only in the first row).                                     | No       | N/A                            |
| `start`             | Start coordinate for the region (required only in the first row).                                       | No       | N/A                            |
| `end`               | End coordinate for the region (required only in the first row).                                         | No       | N/A                            |
| `color`             | Track color.                                                                                           | No       | Default Gviz colors            |
| `height`            | Height of the track in the plot.                                                                        | No       | Default Gviz height            |
| `stacking`          | Stacking mode (`squish`, `dense`, or `full`).                                                           | No       | `squish`                       |
| `ylim`              | Y-axis limits, specified as R code (e.g., `c(0, 100)`).                                                 | No       | Default Gviz limits            |
| `other_annotations` | Any additional metadata or annotations (e.g., manual labels).                                           | No       | N/A                            |

### Example CSV File
```csv
file_path,track_name,type,genome,chromosome,start,end,color,height,stacking,ylim,other_annotations
/path/to/track1.bw,Example BigWig,DataTrack,dm6,chr3L,100000,200000,blue,3,squish,c(0,100),coverage
/path/to/track2.bed,Example BED,AnnotationTrack,dm6,,,,green,2,full,,annotations
/path/to/track3.bam,Example BAM,DataTrack,dm6,,,,red,4,,,
```

### Usage

Running the Workflow
Clone the repository:
Use the following command to clone the repository to your local machine:
git clone https://github.com/your-repo-name.git
cd your-repo-name

* Prepare the input CSV file:

Ensure the CSV file is formatted as described above.
Save the file in a directory accessible to R.

*Run the R Markdown Script:
Open the R Markdown file (gviz_workflow.Rmd) in RStudio. Specify the input CSV file and output directory in the YAML header or dynamically using rmarkdown::render().

*YAML Parameters
Set the following parameters in the YAML header of the R Markdown script:

`input_csv: Path to your input CSV file.`
`output_dir: Directory where outputs will be saved.`
Example YAML configuration:
``` params:
  input_csv: "path_to_your_input_file.csv"
  output_dir: "path_to_output_directory" ```
  
  Alternatively, you can specify these parameters when rendering the script programmatically:
`rmarkdown::render("gviz_workflow.Rmd", params = list(input_csv = "path/to/input.csv", output_dir = "path/to/output_directory"), output_file = "report.html")`

Inspect Outputs:
The plot will be saved as a PNG file in the specified output directory.
Additional logs or metadata can also be saved alongside the plot.
