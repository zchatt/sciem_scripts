# scimet_scripts
_scimet_scripts_ is a collection of scripts used for the analysis of single-cell combinatorial indexing Whole Genome Bisulfite Sequencing (WGBS) data.

## Instructions

### Software Required

[bedtools](https://bedtools.readthedocs.io/en/latest/),
[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (2.2.5),
[parallel](https://www.gnu.org/software/parallel/) (20160222),
[picard](https://broadinstitute.github.io/picard/) (2.7.1),
[bs_seeker2](http://pellegrini-legacy.mcdb.ucla.edu/bs_seeker2/)
[python](https://www.python.org/) (2.7.9),
[R](https://www.r-project.org/) (3.6.3),
[samtools](http://www.htslib.org/) (1.9),
[trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (0.36)

### Inputs
code_location=/project/RDS-FMH-DementiaCFDNA-RW/local_lib/git_repo/scimet_scripts
data_location=/project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet

## Pipeline

## Step-by-step

## 1. Convert bcl to fastq 
files using bcl2fastq2 software (installed on Artemis). This generates .fastq files from 4 lanes of the sequencing flowcell. The files generated are: Read_1, Index_1, Index_2, Read_2.  

	${code_location}/PBS/bcl2fastq_scimet.pbs

## 2a. Combine all lanes of run, rename and split I2 to I2 & I3. 
This script merges R1, R2, I1, and I2 (which contains I2 and I3) from 4 lanes of fastq files.
	
	${code_location}/bash/merge_fq_scimet.sh 

## 2b. Index_2 (I2) 
The I2 index contains 2 barcode sequences concataneted together into a single file that is split into I2 and I3 indexes for downstream processing.

	${code_location}/PBS/split_i2.pbs

## 3. Demultiplexing R1 and R2, QC and Trimming
a) To assign individual reads to a particular barcode/ single-cell, index sequences (I1, I2 and I3) are compared to a list of expected barcodes. A hamming distance <3 algorithm is used to for index match. A concatenated barcode consisting of I1,I2 and I3 is then written to each read name for post-trimming and/or post-alignment demultiplexing. Each R1 and R2 are processed indepedntly. The following scripts are modified from CPT_HiSeq_fq_to_methCPT_split_hamming_PerLane_SE.pl as generously supplied by Andrew Adey.
b) fastqc on pre-trimmed .fastq files
c) Trimming - In addition to trimming the adaper sequences, this will trim low quality bases (denoted by the phred score of 20). Also, any read sequences resulting in a length of 20 bp or lower will not be included in the output file.  
d) fastqc on post-trimmed .fastq files

	${code_location}/PBS/demultiplex_trim_scimet.pbs

## 4. Split fastq by experiment
	FASTQ=${data_dir}/Undetermined.merged_demultiplex_R1.L00.1.fq.gz
	python ${code_location}/python/split_fq_i7.py $FASTQ ${code_location}/data/barcodes_scimet.txt

## 5. fastqc/multiqc per experiment

	${code_location}/PBS/fqc_mqc.pbs

## 6. Trimming per experiment

	${code_location}/PBS/trim_scimet.pbs

## 7. Trimming 2 
	${code_location}/PBS/trim2_scimet.pbs

## 4. Hybrid genome generation
The hybrid hg38_GRCm39_pUC19_Lambda.fa is located /project/RDS-FMH-DementiaCFDNA-RW/local_lib/genomes/normalized_hg38_GRCm39_pUC19_Lambda.fa. The hybrid hg38_GRCm39_pUC19_Lambda.fa genome was created using the following scripts;
	
	${code_location}/PBS/make_genome.pbs

## 5. scBS-MAP alignment.
	
	${code_location}/PBS/scbsmap.pbs

## 6. Annotate .bam with cell and experiment relevent info columns

	for bam_in in *trimmed2.bam;do
	echo $bam_in 
	samtools view -h $bam_in | awk '{$5 != 0; print $0}' | awk 'BEGIN{OFS="\t"}{split($1, a, ":"); print $0,"CB:Z:"a[1]}' | awk -v name=$bam_in 'BEGIN{OFS="\t"}{split($1, a, ":"); print $0,"EX:Z:"name}' | samtools view -bS > ${bam_in%%.bam}_CB.bam
	done

## 7. Split aligned .bam file to single-cells and/ or experiments. Note; needs to be just one .bam file as python will create all cluster.bam subs to spill into.

	samtools merge merged_trimmed2.bam *_trimmed2_CB.bam
	qsub -v bam_in=merged_trimmed2.bam split_bam.pbs

## 8. Bam to CGmap - sorts bam files and extracts DNA methylation information into ATCGmap and CGmap format.

	${code_location}/PBS/bam_to_cgmap.pbs

## 9a. Quantify cross-over. Note 1) cluster, 2) origin, 3) not-origin

	rm tmp
	while read x y z; do
	comm -23 <(samtools view realign_${z}/cluster${x}.bam | awk '{print $1}' | sort) <(samtools view realign_${y}/cluster${x}.bam | awk '{print $1}' | sort) | wc -l >> tmp 
	done < cluster_origin.txt
	mv tmp cross_over.txt

## 9b. select thresholded nuclei; 100 reads/ <10% cross-over
	
	# determine high quality nuclei selection threshold & plotting
	${code_location}/R/crossover_scwgbs.R

	# create soft link with threshold nuclei
	cd  /project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/fastq/results_demultiplex_trim/nuclei_select
	while read x y z; do
	ln -s ../realign_${y}/cluster${x}.CGmap.gz .
	done < nuclei_select.txt

## 10. Global DNA Methyltion statistics/ cell

	for cgmap in *CGmap.gz;do
	echo ${cgmap%%.CGmap.gz}
	cgmaptools mstat -i ${cgmap%%.CGmap.gz}.CGmap.gz -c 1 -p ${cgmap%%.CGmap.gz} > ${cgmap%%.CGmap.gz}.mstat.data
	done

	rm mean_mC
	sed -n '1 p' cluster749.mstat.data >> mean_mC
	for cgmap in *.mstat.data;do
	awk -v a=${cgmap%%.mstat.data} '{print $0,a}' $cgmap | sed -n '2 p' >> mean_mC
	done

## 11. Plotting DNA methylation density accross gene bodies per experiment
	
	${code_location}/PBS/merge_cgmap.pbs
	${code_location}/bash/


## Plotting functions 
	# Tracking reads
	${code_location}/R/read_tracking_plots

	# M-bias
	${code_location}/R/m_bias_scbsmap.R

	# tSNE of global DNA methylation
	${code_location}/R/NMF_tSNE.R

	# multiqc
	pip install multiqc
	multiqc . -o multiqc_out









