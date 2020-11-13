#!/bin/bash
#PBS -P DementiaCFDNA
#PBS -N bam_to_cgmap
#PBS -l select=1:ncpus=8:mem=24GB
#PBS -l walltime=12:00:00

module load cgmaptools
module load parallel

ref_genome=/project/RDS-FMH-DementiaCFDNA-RW/local_lib/genomes/temp2.fa
cd /project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/results_demultiplex_trim/cell_cluster

# make list of .bam files
ls *bam > bam_list
# sort bam file
parallel --xapply -j+0 --eta samtools sort {1} -o {1.}_sorted.bam :::: bam_list
# convert bam to cgmap
parallel --xapply -j+0 --eta cgmaptools convert bam2cgmap -b {1.}_sorted.bam -g $ref_genome -o {1.} :::: bam_list


