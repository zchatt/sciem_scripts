#!/bin/bash
#PBS -P FFbigdata
#PBS -N split_i7
#PBS -l select=1:ncpus=8:mem=16GB
#PBS -l walltime=12:00:00
#PBS -M zacchatt@gmail.com

module load python 

indir=/project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/fastq/results_demultiplex_trim
read_file=Undetermined.merged_demultiplex_R1.L00.1.fq.gz
barcode_file=/project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/barcodes_scimet.txt
code_location=/project/RDS-FMH-DementiaCFDNA-RW/local_lib/git_repo/scimet_scripts

cd $indir

python ${code_location}/python/split_fq_i7.py $read_file $barcode_file