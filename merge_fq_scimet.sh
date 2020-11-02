## We only had 4 lanes. This script combines R1, R2, I1, and I2 from all 4 lanes into 1 file for each of the files. Note: I2 will be split further into I2 and I3.
for n in 1 2 3 4; do cat /project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/Undetermined_S0_L00${n}_I1_001.fastq.gz; done > Undetermined.merged.I1.fq.gz
for n in 1 2 3 4; do cat /project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/Undetermined_S0_L00${n}_R1_001.fastq.gz; done > Undetermined.merged.R1.fq.gz
for n in 1 2 3 4; do cat /project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/Undetermined_S0_L00${n}_R2_001.fastq.gz; done > Undetermined.merged.R2.fq.gz
for n in 1 2 3 4; do cat /project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/Undetermined_S0_L00${n}_I2_001.fastq.gz; done > Undetermined.merged.I2_I3.fq.gz
