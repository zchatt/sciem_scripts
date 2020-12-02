#!/bin/bash
## to be run with .pbs script
# load modules
module load bismark
module load perl
module load bowtie2
module load python/3.8.2

# ## specify variables
# # project directory where demultiplexed .fastq file is contained
# data_dir=/project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/results_demultiplex_trim
# code_dir=/project/RDS-FMH-DementiaCFDNA-RW/local_lib/code/scimet_scripts
# cd ${data_dir}
#
# # reference genome directory (should have run bismark_genome_prepation before alignment)
# BISGENOME=/project/RDS-FMH-DementiaCFDNA-RW/local_lib/genomes/normalized_hg38_pUC19_Lambda
#
# # # bismark genome preparation
# # /project/RDS-FMH-DementiaCFDNA-RW/Praves/git_repo/Bismark/bismark_genome_preparation --bowtie2 --verbose ${BISGENOME}
#
# # sample name
# SAMPLE_NAME=${data_dir}/temp_demultiplex_R1.L00.1_trimmed.fq.gz
#
# # number of cores to use
# NUM_CORES=16
#
# bash ${code_dir}/bismark_align_mbias.sh $SAMPLE_NAME $NUM_CORES $BISGENOME

SAMPLE_NAME=${1}
NUM_CORES=${2}
BISGENOME="${3}"
outDir=$(pwd)/"alignment_results"; mkdir -p $outDir
##########################################################################################################################
##########################################################################################################################
### Alignment of WGBS to HG19+lambda+pUC19 genome
#########################################################################################################################
#########################################################################################################################
# genome preparation should be performed once
## Specify path to bowtie
# bowtie2_path=$(which bowtie2)
## bismark genome preparation
# /project/RDS-FMH-DementiaCFDNA-RW/Praves/git_repo/Bismark/bismark_genome_preparation --bowtie2 --verbose ${BISGENOME}

bismark --fastq --non_directional --un --phred33-quals --score_min L,0,-0.2 \
  --multicore ${NUM_CORES} \
  --output_dir ${outDir} \
  --temp_dir ${outDir} \
  --bowtie2 ${BISGENOME} ${SAMPLE_NAME}

##########################################################################################################################
##########################################################################################################################
### DNA methylation extraction
#########################################################################################################################
#########################################################################################################################
###Reset input file to the aligned .bam file

aligned=$outDir/$(basename ${SAMPLE_NAME} .fq.gz)_bismark_bt2.bam
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
