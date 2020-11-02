#!/bin/bash
# #PBS -P DementiaCFDNA
# #PBS -N demultiplex_trim_scimet
# #PBS -l select=1:ncpus=32:mem=60GB
# #PBS -l walltime=10:00:00
# #PBS -M plam7692@uni.sydney.edu.au
#
# load modules
# module load trimgalore
# module load perl
# module load python/3.8.2
# module load fastqc
#
# # set directories
# data_dir=/project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet
# cd ${data_dir}
# code_dir=/project/RDS-FMH-DementiaCFDNA-RW/local_lib/code/scimet_scripts
#
# # set variables
# SAMPLE_NAME=Undetermined.merged
#
# bash ${code_dir}/demultiplex_trim_scimet.sh $SAMPLE_NAME $data_dir $code_dir
#

# preprocessing steps: combine reads and indexes from lanes 1-4
# ## We only had 4 lanes. This script combines R1, R2, I1, and I2 from all 4 lanes into 1 file for each of the files. Note: I2 will be split further into I2 and I3.
# for n in 1 2 3 4; do cat /project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/Undetermined_S0_L00${n}_I1_001.fastq.gz; done > Undetermined.merged.I1.fq.gz
# for n in 1 2 3 4; do cat /project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/Undetermined_S0_L00${n}_R1_001.fastq.gz; done > Undetermined.merged.R1.fq.gz
# for n in 1 2 3 4; do cat /project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/Undetermined_S0_L00${n}_R2_001.fastq.gz; done > Undetermined.merged.R2.fq.gz
# for n in 1 2 3 4; do cat /project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/Undetermined_S0_L00${n}_I2_001.fastq.gz; done > Undetermined.merged.I2_I3.fq.gz
#
# split I2_I3 into I2 and I3
# zcat ${data_dir}/100000_random_i2_i3.fq.gz | awk 'NR%4==1{print$0} NR%4==2{print substr($1,1,11)} NR%4==3{print$0} NR%4==0{print substr($1,1,11)}' | gzip > 100000_random_i2.fq.gz
# zcat ${data_dir}/100000_random_i2_i3.fq.gz | awk 'NR%4==1{print$0} NR%4==2{print substr($1,12,9)} NR%4==3{print$0} NR%4==0{print substr($1,12,9)}' | gzip > 100000_random_i3.fq.gz

## check if all arguments are passed
if [ -z $1 ]; then
    echo "Need to submit name of sample" && exit
fi

if [ -z $2 ]; then
    echo "Need the path to folder of data" && exit
fi

if [ -z $3 ]; then
    echo "Need the path to folder of code" && exit
fi

## Define file names and project folder
SAMPLE_NAME=${1}

data_dir=${2}
echo "data_dir = $data_dir"

code_dir=${3}
echo "code_dir = $code_dir"

outDir=${data_dir}"/results_demultiplex_trim"
echo $outDir && mkdir -p $outDir

fastQCDir=${outDir}"/fastQC"
echo "FastQCDir is $fastQCDir" && mkdir -p $fastQCDir

fastQC_trim_Dir=${outDir}"/fastQC_post_trim"
echo "FastQCDir_trim is $fastQC_trim_Dir" && mkdir -p $fastQC_trim_Dir

trim_Dir=${outDir}"/fastq_trim"
echo "TrimmedDir is $trim_Dir" && mkdir -p ${trim_Dir}

READ_1=${data_dir}/${SAMPLE_NAME}.R1.fq.gz
READ_2=${data_dir}/${SAMPLE_NAME}.R2.fq.gz
INDEX_1=${data_dir}/${SAMPLE_NAME}.I1.fq.gz
INDEX_2=${data_dir}/${SAMPLE_NAME}.I2.fq.gz
INDEX_3=${data_dir}/${SAMPLE_NAME}.I3.fq.gz

##########################################################################################################################
##########################################################################################################################
### Check if files are there
##########################################################################################################################
##########################################################################################################################
echo "Read_1 is $READ_1"
echo "Read_2 is $READ_2"
echo "Index_1 is $INDEX_1"
echo "Index_2 is $INDEX_2"
echo "Index_3 is $INDEX_3"


files=( "${READ_1}" "${READ_2}" "${INDEX_1}" "${INDEX_2}" "${INDEX_3}" )

for file in ${files[@]}; do
  if [ ! -f $file ];
  then
      find $file -type f
      if [ ! -f $file ];
      then
          echo "Can't find $file fastq file"
      else
          echo "Found" $file
      fi
  else
      echo "Found" $file
  fi
done

##########################################################################################################################
##########################################################################################################################
### demultiplex read_1 with hamming distance 2
##########################################################################################################################
##########################################################################################################################

if [ ! -f ${outDir}"/demultiplex_R1.L00.1.fq.gz" ];
then
  echo "Demultiplexing with HD 2"

  echo "perl ${code_dir}/demultiplex_scimet_prav_modified.pl \
  $READ_1:$INDEX_1:$INDEX_2:$INDEX_3 \
  ${code_dir}/barcodes_scimet.txt ${outDir}/demultiplex_R1"
  echo "[TIME: demultiplex_R1]"
  time perl ${code_dir}/demultiplex_scimet_prav_modified.pl \
  $READ_1:$INDEX_1:$INDEX_2:$INDEX_3 \
  ${code_dir}/barcodes_scimet.txt ${outDir}/demultiplex_R1
else
  echo "Found" ${outDir}"/demultplex_R1.L00.1.fq.gz"
fi

## change file name of demultiplexed R1
mv ${outDir}"/demultiplex_R1.L00.1.fq.gz" ${outDir}"/temp_demultiplex_R1.L00.1.fq.gz"
mv ${outDir}"/demultiplex_R1.L00.fail.1.fq.gz" ${outDir}"/temp_demultiplex_R1.L00.fail.1.fq.gz"

##########################################################################################################################
##########################################################################################################################
### demultiplex read_2 with hamming distance 2
##########################################################################################################################
##########################################################################################################################

if [ ! -f ${outDir}"/demultiplex_R2.L00.1.fq.gz" ];
then
  echo "Demultiplexing with HD 2"

  echo "perl ${code_dir}/demultiplex_scimet_prav_modified.pl \
  $READ_2:$INDEX_1:$INDEX_2:$INDEX_3 \
  ${code_dir}/barcodes_scimet.txt ${outDir}/demultiplex_R2"
  echo "[TIME: demultiplex_R1]"
  time perl ${code_dir}/demultiplex_scimet_prav_modified.pl \
  $READ_2:$INDEX_1:$INDEX_2:$INDEX_3 \
  ${code_dir}/barcodes_scimet.txt ${outDir}/demultiplex_R2
else
  echo "Found" ${outDir}"/demultplex_R2.L00.1.fq.gz"
fi

## merge passed R1 & R2 and failed R1 and R2
cat ${outDir}"/temp_demultiplex_R1.L00.1.fq.gz" ${outDir}"/demultiplex_R2.L00.1.fq.gz" > ${outDir}"/demultiplex_R1_R2.L00.1.fq.gz"
cat ${outDir}"/temp_demultiplex_R1.L00.fail.1.fq.gz" ${outDir}"/demultiplex_R2.L00.fail.1.fq.gz" > ${outDir}"/demultiplex_R1_R2.L00.fail.1.fq.gz"

## change file name variables
demultiplex_hd_pass=${outDir}/demultiplex_R1_R2.L00.1.fq.gz
demultiplex_hd_fail=${outDir}/demultiplex_R1_R2.L00.fail.1.fq.gz

##########################################################################################################################
##########################################################################################################################
### FastQC preTrim
##########################################################################################################################
##########################################################################################################################

if [ ! -f $fastQCDir/$(basename $demultiplex_hd_pass .fq.gz)_fastqc.zip ];
then
  echo "Running fastQC preTrim"

  echo "fastqc \
  --outdir=$fastQCDir \
  $demultiplex_hd_pass"
  echo "[TIME: fastqc postTrim]"
  time fastqc \
  --outdir=$fastQCDir \
  $demultiplex_hd_pass
else
  echo "Found" $fastQCDir/$(basename $demultiplex_hd_pass .fq.gz)_fastqc.zip
fi

##########################################################################################################################
##########################################################################################################################
### Trimming
##########################################################################################################################
##########################################################################################################################

if [ ! -f $trim_Dir/$(basename $demultiplex_hd_pass .fq.gz)_trimmed.fq.gz ];
then
  echo "Running trim galore"

  echo " trim_galore --quality 30 --phred33 -a AGATCGGAAGAGC --stringency 1 -e 0.1 \
  --gzip --length 20 --max_n 10 --output_dir $trim_Dir $demultiplex_hd_pass"
  echo "[TIME: trim galore]"
  time trim_galore --quality 30 --phred33 -a AGATCGGAAGAGC --stringency 1 -e 0.1 \
  --gzip --length 20 --max_n 10 --output_dir $trim_Dir $demultiplex_hd_pass
else
  echo "Found" $trim_Dir/$(basename $demultiplex_hd_pass .fq.gz)_trimmed.fq.gz
fi

###Reset file name to the trimmed fastq's

demultiplex_hd_pass_trimmed=$trim_Dir/$(basename $demultiplex_hd_pass .fq.gz)_trimmed.fq.gz

##########################################################################################################################
##########################################################################################################################
### FastQC postTrim
##########################################################################################################################
##########################################################################################################################

if [ ! -f $fastQC_trim_Dir/$(basename $demultiplex_hd_pass_trimmed .fq.gz)_fastqc.zip ];
then
  echo "Running fastQC postTrim"

  echo "fastqc --outdir=$fastQC_trim_Dir $demultiplex_hd_pass_trimmed"
  echo "[TIME: fastqc postTrim]"
  time fastqc \
  --outdir=$fastQC_trim_Dir \
  $demultiplex_hd_pass_trimmed
else
  echo "Found" $fastQC_trim_Dir/$(basename $demultiplex_hd_pass_trimmed .fq.gz)_fastqc.zip
fi

if [ $? -eq 0 ]
then
  echo "Successfully completed"
  exit 0
else
  echo "Script failed" >&2
  exit 1
fi
