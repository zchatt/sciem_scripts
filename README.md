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

	${code_location}/bcl2fastq_scimet.pbs

## 2a. Combine all lanes of run, rename and split I2 to I2 & I3. 
This script merges R1, R2, I1, and I2 (which contains I2 and I3) from 4 lanes of fastq files.
	
	${code_location}/merge_fq_scimet.sh 

## 2b. Index_2 (I2) 
contains 2 barcode sequences concataneted together into a single file. Therefore, I2 was split into I2 and I3.

	${code_location}/split_i2.pbs

## 3. Demultiplexing R1, then demultiplex R2 by rename.
To assign individual reads to a particular barcode, reads must be first assessed based on the index sequences corresponding to those specific reads. This is accomplished in a process called demultiplexing. Here, all the index sequences (I1, I2, I3) are used to assign a particular read to a prticular cell type. For this, Hamming Distance algorithm is used to account for sequencing errors. This code needs barcode file with i5 index 1 base short i.e 9 bases. Also, i7 tags from only i7_1 - i7_9 were used. So, only those tags put on the barcode file. This was generated with the script "barcodes_scimet_from_csv.py".After Read_1 is demultiplexed, the reads that pass/fail are separated into two files i.e. pass and fail. The names for both the files are changed to temp_{previous_name}. Then, Read_2 is demultiplexed following the same steps. Then, the reads that passed from both R1 and R2 are concatenated into a single R1_R2 file (similar for the failed reads). This R1_R2 file is then passed to fastqc.  The following is modified from modified from CPT_HiSeq_fq_to_methCPT_split_hamming_PerLane_SE.pl as supplied by Andrew Adey.

	#${code_location}/demultiplex_scimet_prav_modified.pl (modified from CPT_HiSeq_fq_to_methCPT_split_hamming_PerLane_SE.pl)
	${code_location}/demultiplex_scimet_r1.pl 
	${code_location}/demultiplex_scimet_r2.pl

## 4. FastQC pre-trim
	
	fastqc --outdir=$fastQCDir
  	$demultiplex_hd_pass -> here, demultiplex_hd_pass is the concatenated file containing R1 and R2 that passed with HD 2.  

## 5. Trimming
In addition to trimming the adaper sequences, this will trim low quality bases (denoted by the phred score of 30). Also, any read sequences resulting in a length of 20 bp or lower will not be included in the output file.  
	
	trim_galore --quality 30 --phred33 -a AGATCGGAAGAGC --stringency 1 -e 0.1 \
  --gzip --length 20 --max_n 10 --output_dir $trim_Dir $demultiplex_hd_pass

## 6. FastQC post-trim

	fastqc --outdir=$fastQC_trim_Dir \

	$demultiplex_hd_pass_trimmed -> here, demultiplexed_hd_pass_trimmed is the file obtained after trimming $demultiplex_hd_pass file. 

## 7a. Hybrid genome
The hybrid hg38_GRCm39_pUC19_Lambda.fa is located /project/RDS-FMH-DementiaCFDNA-RW/local_lib/genomes/normalized_hg38_GRCm39_pUC19_Lambda.fa. The hybrid hg38_GRCm39_pUC19_Lambda.fa genome was created using the following scripts;
	
	${code_location}/make_genome.pbs

## 7b. Alignment to human-mouse-pUC19-Lambda hybrid genome using scBS-MAP. 
	
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










