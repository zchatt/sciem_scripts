# scimet_scripts
_scimet_scripts_ is s collection of scripts used for the analysis of single-cell combinatorial indexing Whole Genome Bisulfite Sequencing (WGBS) data.

## Instructions

### Software Required

[bedtools](https://bedtools.readthedocs.io/en/latest/),
[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (2.2.5),
[parallel](https://www.gnu.org/software/parallel/) (20160222),
[picard](https://broadinstitute.github.io/picard/) (2.7.1),
[](http://pellegrini-legacy.mcdb.ucla.edu/bs_seeker2/)
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

	${code_location}/bcl2fastq_scimet.sub


## 2a. Combine all lanes of run, rename and split I2 to I2 & I3. 
This script merges R1, R2, I1, and I2 (which contains I2 and I3) from 4 lanes of fastq files.
	
	${code_location}/merge_fq_scimet.sh 

## 2b. Index_2 (I2) 
contains 2 barcode sequences concataneted together into a single file. Therefore, I2 was split into I2 and I3.

	${code_location}/split_i2.pbs

## 3. Demultiplexing R1, then demultiplex R2 by rename.
To assign individual reads to a particular barcode, reads must be first assessed based on the index sequences corresponding to those specific reads. This is accomplished in a process called demultiplexing. Here, all the index sequences (I1, I2, I3) are used to assign a particular read to a prticular cell type. For this, Hamming Distance algorithm is used to account for sequencing errors. This code needs barcode file with i5 index 1 base short i.e 9 bases. Also, i7 tags from only i7_1 - i7_9 were used. So, only those tags put on the barcode file. This was generated with the script "barcodes_scimet_from_csv.py".After Read_1 is demultiplexed, the reads that pass/fail are separated into two files i.e. pass and fail. The names for both the files are changed to temp_{previous_name}. Then, Read_2 is demultiplexed following the same steps. Then, the reads that passed from both R1 and R2 are concatenated into a single R1_R2 file (similar for the failed reads). This R1_R2 file is then passed to fastqc.  The following is modified from modified from CPT_HiSeq_fq_to_methCPT_split_hamming_PerLane_SE.pl as supplied by Andrew Adey.

	${code_location}/demultiplex_scimet_prav_modified.pl (modified from CPT_HiSeq_fq_to_methCPT_split_hamming_PerLane_SE.pl)

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
	
	cat /project/RDS-SMS-FFbigdata-RW/local_lib/genomes/hg38/hg38.fa /project/RDS-SMS-FFbigdata-RW/local_lib/genomes/GRCm39/GCF_000001635.27_GRCm39_genomic.fna /project/RDS-SMS-FFbigdata-RW/local_lib/genomes/pUC19/pUC19_addgene.fa /project/RDS-SMS-FFbigdata-RW/local_lib/genomes/lambda/lambda.fa > hg38_GRCm39_pUC19_Lambda.fa

	picard NormalizeFasta I=hg38_GRCm39_pUC19_Lambda.fa O=normalized_hg38_GRCm39_pUC19_Lambda.fa

	samtools faidx normalized_hg38_GRCm39_pUC19_Lambda.fa

## In Development - Additional steps to add into pip
## bbmap (installed on Artemis) - why???
## scBS-MAP (Dependencies: SAMtools, bowtie2, BS-Seeker2 (This requires bowtie2, python2.6 or higher, pysam)). From this, install BS-Seeker and pysam. Needed fro alignment
## bamtools (installed on Artemis)  - handling bam files post alignment
## CGmap tools (not installed on Artemis) - post alignment bisulfite seq handling
## Clustering (Can be done using SciKit-learn - not installed on Artemis)

## 7b. Alignment to human-mouse-pUC19-Lambda hybrid genome using scBS-MAP. 
	
	export PATH=//project/RDS-FMH-DementiaCFDNA-RW/local_lib/BSseeker2-2.1.1/:$PATH

## 8. Split alignment to single-cells


## QuickRun - A single bash script was written to combine all the steps (1-Nth) above to this point. 
	
	${code_location}/demultiplex_trim_scimet.sh
