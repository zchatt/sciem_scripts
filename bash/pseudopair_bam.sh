#!/bin/bash

## run example
# pseudopair_bam.sh sample_R1.bam sample_R1.bam

BAM1=$1
BAM2=$2

samtools view $BAM1 > ${BAM1%%.bam}.sam
samtools view $BAM2 > ${BAM2%%.bam}.sam
 
# join by read names
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.3,2.4,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16 <(sort -k1 ${BAM1%%.bam}_R1.sam) <(sort -k1 ${BAM2%%.bam}_R2.sam) > ${BAM1%%.bam}_tmp

# Falsify FLAG and calculate TLEN
awk '{print $1"\t""83""\t"$3"\t"$4"\t"$5"\t"$6"\t"$8"\t"$8"\t"($8 + 142 - $4)"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16'} ${BAM1%%.bam}_tmp > ${BAM1%%.bam}_tmp2

# combine header of R1 sam file
cat <(samtools view -H $BAM1) ${BAM1%%.bam}_tmp2 > ${BAM1%%_R1.bam}_pseudopaired.sam

# cleanup tmp files
rm ${BAM1%%.bam}_tmp ${BAM1%%.bam}_tmp2

