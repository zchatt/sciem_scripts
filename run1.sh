#!/bin/bash
#PBS -P FFbigdata
#PBS -N run1
#PBS -l select=1:ncpus=20:mem=96GB
#PBS -l walltime=20:00:00
#PBS -M zacchatt@gmail.com

module load bcl2fastq2
module load trimgalore
module load perl
module load python/2.7.9
module load fastqc
module load bbmap
module load cutadapt
module load bowtie2
module load samtools
module load fastp

# run using bcl2fastq v2.19.0.316 software
tar=/project/RDS-FMH-FFEPIGENETICS-RW/scwgbs/210804_A00152_0453_AHHGWYDRXY.tar
topdir=/project/RDS-FMH-FFEPIGENETICS-RW/scwgbs
runfolder_dir=/project/RDS-FMH-FFEPIGENETICS-RW/scwgbs/210804_A00152_0453_AHHGWYDRXY
input_dir=${runfolder_dir}/Data/Intensities/BaseCalls
code_dir=/project/RDS-FMH-DementiaCFDNA-RW/local_lib/git_repo/scimet_scripts
barcodes_scimet=${topdir}/barcodes_scimet_10_11_10_NovaSeq.txt
SAMPLE_NAME=Undetermined.merged
scbsmap_location=/project/RDS-FMH-DementiaCFDNA-RW/local_lib/scBS-map/
bs_seeker2_location=/project/RDS-FMH-DementiaCFDNA-RW/local_lib/BSseeker2-2.1.1/
ref_genome=/project/RDS-FMH-DementiaCFDNA-RW/local_lib/BSseeker2-2.1.1/bs_utils/reference_genomes/normalized_hg38_pUC19_Lambda.fa
export PATH=/project/RDS-FMH-DementiaCFDNA-RW/local_lib/BSseeker2-2.1.1/:$PATH

# ## 0. untar
# cd $topdir
# tar -xvf $tar

# 1. Convert bcl to fastq 
cd $runfolder_dir

bcl2fastq \
-i $input_dir \
--create-fastq-for-index-reads \
## --with-failed-reads \ # cannot be applied to cbcl file from novaseq
--use-bases-mask Y*,I10,I21,Y*

## 2a. Combine all lanes of run, rename and split I2 to I2 & I3. 
## We only had 4 lanes. This script combines R1, R2, I1, and I2 from all 4 lanes into 1 file for each of the files. Note: I2 will be split further into I2 and I3.
for n in 1 2 3 4; do cat ${runfolder_dir}/Data/Intensities/BaseCalls/Undetermined_S0_L00${n}_I1_001.fastq.gz; done > Undetermined.merged.I1.fq.gz
for n in 1 2 3 4; do cat ${runfolder_dir}/Data/Intensities/BaseCalls/Undetermined_S0_L00${n}_R1_001.fastq.gz; done > Undetermined.merged.R1.fq.gz
for n in 1 2 3 4; do cat ${runfolder_dir}/Data/Intensities/BaseCalls/Undetermined_S0_L00${n}_R2_001.fastq.gz; done > Undetermined.merged.R2.fq.gz
for n in 1 2 3 4; do cat ${runfolder_dir}/Data/Intensities/BaseCalls/Undetermined_S0_L00${n}_I2_001.fastq.gz; done > Undetermined.merged.I2_I3.fq.gz

## 2b. Index_2 (I2) 
#The I2 index contains 2 barcode sequences concataneted together into a single file that is split into I2 and I3 indexes for downstream processing. Note the difference in NovaSeq and NextSeq indexing when splitting
## split I2 into I2 and I3
zcat ${runfolder_dir}/Undetermined.merged.I2_I3.fq.gz | awk 'NR%4==1{print$0} NR%4==2{print substr($1,11,11)} NR%4==3{print$0} NR%4==0{print substr($1,11,11)}' | gzip > Undetermined.merged.I2.fq.gz
zcat ${runfolder_dir}/Undetermined.merged.I2_I3.fq.gz | awk 'NR%4==1{print$0} NR%4==2{print substr($1,1,10)} NR%4==3{print$0} NR%4==0{print substr($1,1,10)}' | gzip > Undetermined.merged.I3.fq.gz

# 3. Demultiplexing R1 and R2, QC and Trimming
bash ${code_dir}/bash/demultiplex_trim_scimet.sh $SAMPLE_NAME $runfolder_dir $code_dir $barcodes_scimet 

## 4a. Trim 2 - linear primer sequences 
cd ${runfolder_dir}/results_demultiplex_trim

for FQ in *_trimmed.fq.gz;do
echo $FQ
cutadapt --anywhere=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
 --anywhere=GAAGAGCACACGTCTGAACTC \
 --anywhere=ATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGG \
 --minimum-length=20 --times=2 -o ${FQ%%_trimmed.fq.gz}_trimmed2.fq.gz $FQ
done

# ## 4b. Trim 3 - hard trim to 60bp for R2 dud to intesity ddop-off and overcalling of "G"
cd ${runfolder_dir}/results_demultiplex_trim
fastp --max_len1 60 -i Undetermined.merged_demultiplex_R2.L00.1_trimmed2.fq.gz -o Undetermined.merged_demultiplex_R2.L00.1_trimmed3.fq.gz

## 5. scBS-MAP alignment of R1
cd ${runfolder_dir}/results_demultiplex_trim
FQ=Undetermined.merged_demultiplex_R1.L00.1_trimmed2.fq.gz
perl $scbsmap_location/scBS-map.pl -l 10 -p 20 -n 10 -f $FQ -g $ref_genome -o ${FQ%%_trimmed2.fq.gz}.bam
FQ=Undetermined.merged_demultiplex_R2.L00.1_trimmed3.fq.gz
perl $scbsmap_location/scBS-map.pl -l 10 -p 20 -n 10 -f $FQ -g $ref_genome -o ${FQ%%_trimmed3.fq.gz}.bam

## 6. split and QC of R1
mkdir ${runfolder_dir}/fq_split
cd ${runfolder_dir}/fq_split
demuxbyname.sh in=${runfolder_dir}/results_demultiplex_trim/Undetermined.merged_demultiplex_R1.L00.1.fq.gz out=%_R1.fq length=31 prefixmode=t
demuxbyname.sh in=${runfolder_dir}/results_demultiplex_trim/Undetermined.merged_demultiplex_R1.L00.1_trimmed.fq.gz out=%_trimmed_R1.fq length=31 prefixmode=t
demuxbyname.sh in=${runfolder_dir}/results_demultiplex_trim/Undetermined.merged_demultiplex_R1.L00.1_trimmed2.fq.gz out=%_trimmed2_R1.fq length=31 prefixmode=t
gzip *.fq
fastqc *fq.gz
#multiqc . -o .





#### intermediate testing scripts ####
# ### evaluating R2 trimming
# module load fastp

# cd /project/RDS-FMH-FFEPIGENETICS-RW/scwgbs/210804_A00152_0453_AHHGWYDRXY/results_demultiplex_trim
# zcat Undetermined.merged_demultiplex_R2.L00.1_trimmed2.fq.gz | head -4000000 | gzip > r2_trim_analysis/1M_R2.L00.1_trimmed2.fq.gz

# cd r2_trim_analysis
# fastqc 1M_R2.L00.1_trimmed2.fq.gz

# # trimming of polyG
# fastp -g --poly_g_min_len 3 -i 1M_R2.L00.1_trimmed2.fq.gz -o 1M_R2.L00.1_trimmed3.fq.gz

# # trimming to 60bp R2
# fastp --max_len1 60 -i 1M_R2.L00.1_trimmed2.fq.gz -o 1M_R2.L00.1_trimmed3.fq.gz

