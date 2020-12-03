#!/bin/bash
#PBS -P DementiaCFDNA
#PBS -N bismark_align_mbias
#PBS -l select=1:ncpus=24:mem=96GB
#PBS -l walltime=24:00:00

# load modules
module load bismark
module load perl
module load bowtie2
module load python/3.8.2

## specify variables
data_dir=/project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/results_demultiplex_trim
code_dir=/project/RDS-FMH-DementiaCFDNA-RW/local_lib/code/scimet_scripts
BISGENOME=/project/RDS-FMH-DementiaCFDNA-RW/local_lib/genomes/normalized_hg38_GRCm39_pUC19_Lambda # reference genome directory (should have run bismark_genome_prepation before alignment)
NUM_CORES=1
FASTQ=${data_dir}/temp_demultiplex_R1.L00.1_trimmed.fq.gz

cd ${data_dir}

# # bismark genome preparation - if not done
# bismark_genome_preparation --bowtie2 --verbose ${BISGENOME}

# run bismark align and methylation extraction
bash ${code_dir}/bismark_align_mbias.sh $FASTQ $NUM_CORES $BISGENOME