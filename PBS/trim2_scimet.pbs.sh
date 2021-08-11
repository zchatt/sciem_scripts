#!/bin/bash
#PBS -P DementiaCFDNA
#PBS -N trim2_scimet
#PBS -l select=1:ncpus=8:mem=24GB
#PBS -l walltime=10:00:00
#PBS -M plam7692@uni.sydney.edu.au

# load modules
module load cutadapt

# set directories and variables
data_dir=/project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/fastq/results_demultiplex_trim
cd $data_dir

# run over all trimmed .fastq files
for FQ in *_trimmed.fq.gz;do
	cutadapt --anywhere=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --anywhere=GAAGAGCACACGTCTGAACTC --anywhere=ATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGG --minimum-length=20 --times=2 -o ${FQ%%_trimmed.fq.gz}_trimmed2.fq.gz $FQ
done
