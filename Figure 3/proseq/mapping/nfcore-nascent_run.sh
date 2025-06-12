#!bin/bash
#David Shechter 2023-11-22
#script to run nf-core/nascent on proseq data
#use refgenie builds in shared folder

#Load conda envrionment and modules
.  /gs/gsfs0/hpc01/rhel8/apps/conda3/bin/activate
conda activate zlib_nextflow_rnaseq
module load singularity
conda activate nextflow

module load singularity
conda activate nextflow

nextflow run nf-core/nascent -r 2.2.0 --input `pwd`/samplesheet.csv \
-profile singularity \
--outdir `pwd`/nfcore_out \
--assay_type PROseq \
--aligner bwa \
--fasta /gs/gsfs0/users/shechter-lab/data/NGS/stds/refgenie/alias/hg38/fasta/default/hg38.fa \
--gtf /gs/gsfs0/users/shechter-lab/data/NGS/stds/refgenie/alias/hg38/ensembl_gtf/default/hg38.gtf.gz \
--bwa_index /gs/gsfs0/users/shechter-lab/data/NGS/stds/refgenie/alias/hg38/bwa_index/default/ \
-resume


conda deactivate
module unload singularity
conda deactivate

###END NEXTFLOW###
