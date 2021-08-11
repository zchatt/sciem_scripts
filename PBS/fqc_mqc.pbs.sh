#!/bin/bash
#PBS -P FFbigdata
#PBS -N fqc_mqc
#PBS -l select=1:ncpus=8:mem=16GB
#PBS -l walltime=2:00:00
#PBS -M zacchatt@gmail.com

module load python
module load fastqc
pip install multiqc

indir=/project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/fastq/results_demultiplex_trim
cd $indir

# fastqc/ multiqc
#mkdir fastQC
#fastqc *passed_reads*
#mv *passed_reads_fastqc* fastQC
#cd fastQC
#multiqc . -o .

mkdir fastQC_posttrim
fastqc *trimmed*
mv *trimmed_fastqc* fastQC_posttrim
cd fastQC_posttrim
multiqc . -o .