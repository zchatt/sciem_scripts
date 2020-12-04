#!/bin/bash
# The script performs demultiplexing and trimming with fastqc QC
## check if all arguments are passed
if [ -z $1 ]; then
    echo "Need to submit name of sample" && exit
fi

if [ -z $2 ]; then
    echo "Need the path to data" && exit
fi

if [ -z $3 ]; then
    echo "Need the path to code" && exit
fi

if [ -z $4 ]; then
    echo "Need barcodes file" && exit
fi


## Define file names and project folder
SAMPLE_NAME=${1}

data_dir=${2}
echo "data_dir = $data_dir"

code_dir=${3}
echo "code_dir = $code_dir"

barcodes_scimet=${4}
echo "barcodes_scimet = $barcodes_scimet"

outDir=${data_dir}"/results_demultiplex_trim"
echo $outDir && mkdir -p $outDir

fastQCDir=${outDir}"/fastQC"
echo "FastQCDir is $fastQCDir" && mkdir -p $fastQCDir

fastQC_trim_Dir=${outDir}"/fastQC_post_trim"
echo "FastQCDir_trim is $fastQC_trim_Dir" && mkdir -p $fastQC_trim_Dir

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

if [ ! -f ${outDir}"/${SAMPLE_NAME}_demultiplex_R1.fq.gz" ];
then
  echo "Demultiplexing R1 with HD==2"

  echo "perl ${code_dir}/perl/demultiplex_scimet_r1.pl \
  $READ_1:$INDEX_1:$INDEX_2:$INDEX_3 \
  $barcodes_scimet ${outDir}/${SAMPLE_NAME}_demultiplex_R1"
  echo "[TIME: demultiplex_R1]"
  time perl ${code_dir}/perl/demultiplex_scimet_r1.pl \
  $READ_1:$INDEX_1:$INDEX_2:$INDEX_3 \
  $barcodes_scimet ${outDir}/${SAMPLE_NAME}_demultiplex_R1
else
  echo "Found" ${outDir}"/${SAMPLE_NAME}_demultiplex_R1.fq.gz"
fi

##########################################################################################################################
##########################################################################################################################
### demultiplex read_2 with hamming distance 2
##########################################################################################################################
##########################################################################################################################

if [ ! -f ${outDir}"/${SAMPLE_NAME}_demultiplex_R2.fq.gz" ];
then
  echo "Demultiplexing R2 with HD==2"

  echo "perl ${code_dir}/perl/demultiplex_scimet_r2.pl \
  $READ_2:$INDEX_1:$INDEX_2:$INDEX_3 \
  $barcodes_scimet ${outDir}/${SAMPLE_NAME}_demultiplex_R2"
  echo "[TIME: demultiplex_R2]"
  time perl ${code_dir}/perl/demultiplex_scimet_r2.pl \
  $READ_2:$INDEX_1:$INDEX_2:$INDEX_3 \
  $barcodes_scimet ${outDir}/${SAMPLE_NAME}_demultiplex_R2
else
  echo "Found" ${outDir}"/${SAMPLE_NAME}_demultiplex_R2.fq.gz"
fi

## assign variable names
demultiplex_hd_pass_r1=${outDir}/${SAMPLE_NAME}_demultiplex_R1.L00.1.fq.gz
demultiplex_hd_pass_r2=${outDir}/${SAMPLE_NAME}_demultiplex_R2.L00.1.fq.gz

demultiplex_hd_fail_r1=${outDir}/${SAMPLE_NAME}_demultiplex_R1.L00.fail.1.fq.gz
demultiplex_hd_fail_r2=${outDir}/${SAMPLE_NAME}_demultiplex_R2.L00.fail.1.fq.gz

##########################################################################################################################
##########################################################################################################################
### FastQC preTrim
##########################################################################################################################
##########################################################################################################################

## Read 1
if [ ! -f $fastQCDir/$(basename $demultiplex_hd_pass_r1 .fq.gz)_fastqc.zip ];
then
  echo "Running fastQC preTrim"

  echo "fastqc \
  --outdir=$fastQCDir \
  $demultiplex_hd_pass_r1"
  echo "[TIME: fastqc postTrim]"
  time fastqc \
  --outdir=$fastQCDir \
  $demultiplex_hd_pass_r1
else
  echo "Found" $fastQCDir/$(basename $demultiplex_hd_pass_r1 .fq.gz)_fastqc.zip
fi

## Read 2
if [ ! -f $fastQCDir/$(basename $demultiplex_hd_pass_r2 .fq.gz)_fastqc.zip ];
then
  echo "Running fastQC preTrim"

  echo "fastqc \
  --outdir=$fastQCDir \
  $demultiplex_hd_pass_r2"
  echo "[TIME: fastqc postTrim]"
  time fastqc \
  --outdir=$fastQCDir \
  $demultiplex_hd_pass_r2
else
  echo "Found" $fastQCDir/$(basename $demultiplex_hd_pass_r2 .fq.gz)_fastqc.zip
fi

##########################################################################################################################
##########################################################################################################################
### Trimming
##########################################################################################################################
##########################################################################################################################

## Read 1
if [ ! -f $outDir/$(basename $demultiplex_hd_pass_r1 .fq.gz)_trimmed.fq.gz ];
then
  echo "Running trim galore"

  echo " trim_galore --quality 20 --phred33 --illumina --stringency 3 -e 0.1 \
  --gzip --length 20 --output_dir $outDir $demultiplex_hd_pass_r1"
  echo "[TIME: trim galore]"
  time trim_galore --quality 20 --phred33 --illumina --stringency 3 -e 0.1 \
  --gzip --length 20 --output_dir $outDir $demultiplex_hd_pass_r1
else
  echo "Found" $outDir/$(basename $demultiplex_hd_pass_r1 .fq.gz)_trimmed.fq.gz
fi

## Read 2
if [ ! -f $outDir/$(basename $demultiplex_hd_pass_r2 .fq.gz)_trimmed.fq.gz ];
then
  echo "Running trim galore"

  echo " trim_galore --quality 20 --phred33 --illumina --stringency 3 -e 0.1 \
  --gzip --length 20 --output_dir $outDir $demultiplex_hd_pass_r2"
  echo "[TIME: trim galore]"
  time trim_galore --quality 20 --phred33 --illumina --stringency 3 -e 0.1 \
  --gzip --length 20 --output_dir $outDir $demultiplex_hd_pass_r2
else
  echo "Found" $outDir/$(basename $demultiplex_hd_pass_r2 .fq.gz)_trimmed.fq.gz
fi

## assign variable names
demultiplex_hd_pass_r1_trimmed=$outDir/$(basename $demultiplex_hd_pass_r1 .fq.gz)_trimmed.fq.gz
demultiplex_hd_pass_r2_trimmed=$outDir/$(basename $demultiplex_hd_pass_r2 .fq.gz)_trimmed.fq.gz

##########################################################################################################################
##########################################################################################################################
### FastQC postTrim
##########################################################################################################################
##########################################################################################################################

## Read 1
if [ ! -f $fastQC_trim_Dir/$(basename $demultiplex_hd_pass_r1_trimmed .fq.gz)_fastqc.zip ];
then
  echo "Running fastQC postTrim"

  echo "fastqc --outdir=$fastQC_trim_Dir $demultiplex_hd_pass_r1_trimmed"
  echo "[TIME: fastqc postTrim]"
  time fastqc \
  --outdir=$fastQC_trim_Dir \
  $demultiplex_hd_pass_r1_trimmed
else
  echo "Found" $fastQC_trim_Dir/$(basename $demultiplex_hd_pass_r1_trimmed .fq.gz)_fastqc.zip
fi

## Read 2
if [ ! -f $fastQC_trim_Dir/$(basename $demultiplex_hd_pass_r2_trimmed .fq.gz)_fastqc.zip ];
then
  echo "Running fastQC postTrim"

  echo "fastqc --outdir=$fastQC_trim_Dir $demultiplex_hd_pass_r2_trimmed"
  echo "[TIME: fastqc postTrim]"
  time fastqc \
  --outdir=$fastQC_trim_Dir \
  $demultiplex_hd_pass_r2_trimmed
else
  echo "Found" $fastQC_trim_Dir/$(basename $demultiplex_hd_pass_r2_trimmed .fq.gz)_fastqc.zip
fi

