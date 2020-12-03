#!/bin/bash
FASTQ=${1}
NUM_CORES=${2}
BISGENOME="${3}"
outDir=$(pwd)/"alignment_results"; mkdir -p $outDir
##########################################################################################################################
##########################################################################################################################
### Alignment using bismark
#########################################################################################################################
#########################################################################################################################
bismark --fastq --non_directional --un --phred33-quals --score_min L,0,-0.2 \
  --multicore ${NUM_CORES} \
  --output_dir ${outDir} \
  --temp_dir ${outDir} \
  --bowtie2 ${BISGENOME} ${FASTQ}

##########################################################################################################################
##########################################################################################################################
### DNA methylation extraction
#########################################################################################################################
#########################################################################################################################
###Reset input file to the aligned .bam file

aligned=$outDir/$(basename ${FASTQ} .fq.gz)_bismark_bt2.bam
bismark_methylation_extractor -s \
  --comprehensive \
  --cytosine_report \
  --bedGraph \
  --CX \
  --genome_folder $BISGENOME \
  --multicore ${NUM_CORES} \
  --output ${outDir} \
  --gzip $aligned

## m-bias plot
# python ${code_dir}/m_bias_plot.py $(basename ${aligned} .bam).M-bias.txt

if [ $? -eq 0 ]
then
  echo "Successfully completed"
  exit 0
else
  echo "Script failed" >&2
  exit 1
fi
