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
code_location=/project/RDS-FMH-DementiaCFDNA-RW/local_lib/code/scimet_scripts
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

## 3a. Demultiplexing R1 and R2.
To assign individual reads to a particular barcode/ single-cell, index sequences (I1, I2 and I3) are compared to a list of expected barcodes. A hamming distance <3 algorithm is used to for index match. A concatenated barcode consisting of I1,I2 and I3 is then written to each read name for post-trimming and/or post-alignment demultiplexing. Each R1 and R2 are processed indepedntly. The following scripts are modified from CPT_HiSeq_fq_to_methCPT_split_hamming_PerLane_SE.pl as generously supplied by Andrew Adey.
## 3b. FastQC pre-trim
## 3c. Trimming
In addition to trimming the adaper sequences, this will trim low quality bases (denoted by the phred score of 20). Also, any read sequences resulting in a length of 20 bp or lower will not be included in the output file.  
## 3d. FastQC post-trim

	${code_location}/PBS/demultiplex_trim_scimet.pbs

## 4. Hybrid genome generation
The hybrid hg38_GRCm39_pUC19_Lambda.fa is located /project/RDS-FMH-DementiaCFDNA-RW/local_lib/genomes/normalized_hg38_GRCm39_pUC19_Lambda.fa. The hybrid hg38_GRCm39_pUC19_Lambda.fa genome was created using the following scripts;
	
	${code_location}/make_genome.pbs

## 5. Alignment using scBS-MAP. 
	
	${code_location}/scbsmap.pbs

## 8. Split aligned .bam file to single-cells and/ or experiments
a) Select unique reads and create cell-barcode col "CB" 

	bam_in=demultiplex_illumina_trimmed_R1.bam

	samtools view -h $bam_in | awk '{$5 != 0; print $0}' | awk 'BEGIN{OFS="\t"}{split($1, a, ":"); print $0,"CB:Z:"a[1]}' | samtools view -bS > tmp.bam

b) Split .bam file based on "Cluster". The "Cluster" can be set to per-cell, per-experiment or any other grouping.
	
	cluster_file=/project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/clusters.csv

	module load python
	${code_location}/split_bam.py $cluster_file tmp.bam

## In Development - Additional steps to add into pip
	bbmap (installed on Artemis) - why???
	bamtools (installed on Artemis)  - handling bam files post alignment
	CGmap tools (not installed on Artemis) - post alignment bisulfite seq handling
	Clustering (Can be done using SciKit-learn - not installed on Artemis)

## 9. Collect bam stats for each experiment

	module load samtools
	cd /project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/results_demultiplex_trim/cell_cluster
	cluster_file=/project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/clusters.csv

	# format/ remove header
	head -1 $cluster_file > header
	awk -F "\"*,\"*" 'FNR > 1 {print $0}' $cluster_file > tmp

	# collect stats 1)#reads, 2)#reads hg38,3) #reads GRCm39, 4) #reads pUC19, 5) #reads Lambda
	rm summary_data.csv
	while read -r line ; do
	cell=$(echo $line | cut -d "," -f 2)

	reads1=$(samtools view -c cluster${cell}.bam)
	reads2=$(samtools view cluster${cell}.bam | awk '{split($3,a,"r"); print a[1]}' | sort | uniq -c | grep "ch" | awk '{print $1}')
	reads3=$(samtools view cluster${cell}.bam | awk '{split($3,a,"_"); print a[1]}' | sort | uniq -c | grep "NC" | awk '{print $1}')
	reads4=$(samtools view cluster${cell}.bam | awk '{split($3,a,"_"); print a[1]}' | sort | uniq -c | grep "pUC19" | awk '{print $1}')
	reads5=$(samtools view cluster${cell}.bam | awk '{split($3,a,"_"); print a[1]}' | sort | uniq -c | grep "Lambda" | awk '{print $1}')

	echo "$cell,$reads1,$reads2,$reads3,$reads4,$reads5"
	done < tmp > summary_data.csv

# plotting functions with R

	${code_location}/low_level_scwgbs.R

## multiqc
	pip install multiqc
	multiqc . -o multiqc_out

## By cell
	# Mean unique reads per cell  - summarise for each experiment
	# Collision rate; 1) % unique aligned mouse/ human genome for each cell 2) set threshold for annotating mouse OR human (% of cells distinctly mouse or human) 3) within distinct mouse or human % of reads uniquely aligned to other.

## By experiment
	# Mean alignment of each experiment
	# Number of barcodes identified compared to expected & % useable reads (uniquely aligned reads/all reads assigned to a barcode)

## 10. Bam to CGmap - sorts bam files and extracts DNA methylation information into ATCGmap and CGmap format.

	${code_location}/bam_to_cgmap.pbs

## 1. CGmap coverage cgmaptools oac


## QuickRun - A single bash script was written to combine all the steps (1-Nth) above to this point. 
	
	${code_location}/demultiplex_trim_scimet.sh










