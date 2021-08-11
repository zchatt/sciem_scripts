# script for plotting density around all gene bodies in hg38 and mm10
module load bedtools
module load cgmaptools

# obtain TSS and TES for protin coding genes annotated to hg38 
cd /project/RDS-SMS-FFbigdata-RW/local_lib/genomes/hg38/annotation
wgett http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37gencode.v37.annotation.gtf.gz 
gunzip gencode.v37.annotation.gtf
grep 'gene_type "protein_coding"' gencode.v37.annotation.gtf | awk '($3=="gene") {printf("%s\t%s\t%s\t%s\n",$1,int($4)-1,$5,$7);}' | sort -T . -k1,1 -k2,2n > gencode_hg38_proteincoding_genes.bed

# create fragmented regions for all bed regions
cat gencode_hg38_proteincoding_genes.bed | cgmaptools bed2fragreg -n 30 -F 1000,1000,1000,1000,1000,1000,1000,1000,1000,1000 -T 1000,1000,1000,1000,1000,1000,1000,1000,1000,1000 > fragreg.bed

# obtain TSS and TES for protein coding genes annotated to mm10/GRCm39
cd /project/RDS-SMS-FFbigdata-RW/local_lib/genomes/GRCm39/annotation

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/gencode.vM26.annotation.gtf.gz
gunzip gencode.vM26.annotation.gtf.gz
grep 'gene_type "protein_coding"' gencode.vM26.annotation.gtf | awk '($3=="gene") {printf("%s\t%s\t%s\t%s\n",$1,int($4)-1,$5,$7);}' | sort -T . -k1,1 -k2,2n > gencode_GRCm39_proteincoding_genes.bed

# change chromosome names to refseq accession to match ref genome used for alignment
sed -i 's/chr1/NC_000067.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr2/NC_000068.8/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr3/NC_000069.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr4/NC_000070.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr5/NC_000071.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr6/NC_000072.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr7/NC_000073.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr8/NC_000074.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr9/NC_000075.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr10/NC_000076.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr11/NC_000077.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr12/NC_000078.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr13/NC_000079.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr14/NC_000080.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr15/NC_000081.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr16/NC_000082.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr17/NC_000083.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr18/NC_000084.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chr19/NC_000085.7/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chrX/NC_000086.8/g' gencode_GRCm39_proteincoding_genes.bed
sed -i 's/chrY/NC_000087.8/g' gencode_GRCm39_proteincoding_genes.bed

# create fragmented regions for all bed regions
cat gencode_GRCm39_proteincoding_genes.bed  | cgmaptools bed2fragreg -n 30 -F 1000,1000,1000,1000,1000,1000,1000,1000,1000,1000 -T 1000,1000,1000,1000,1000,1000,1000,1000,1000,1000 > fragreg.bed

# combine fragreg.bed from human and mouse
cd /project/RDS-FMH-DementiaCFDNA-RW/Epigenetics/scimet/fastq/results_demultiplex_trim/nuclei_select

cat /project/RDS-SMS-FFbigdata-RW/local_lib/genomes/hg38/annotation/fragreg.bed \
<(awk 'NR>1' /project/RDS-SMS-FFbigdata-RW/local_lib/genomes/GRCm39/annotation/fragreg.bed) > fragreg_hg38GRCm39_genes.bed

# calculate methylation in FraGmented regions
gunzip -c merge_sciEMMG.CGmap.gz | cgmaptools mfg -r fragreg_hg38GRCm39_genes.bed -c 1 -x CG > merge_sciEMMG.CG_genes.mfg
gunzip -c merge_sciEMN9.CGmap.gz | cgmaptools mfg -r fragreg_hg38GRCm39_genes.bed -c 1 -x CG > merge_sciEMN9.CG_genes.mfg
gunzip -c merge_sciMETMG.CGmap.gz | cgmaptools mfg -r fragreg_hg38GRCm39_genes.bed -c 1 -x CG > merge_sciMETMG.CG_genes.mfg
gunzip -c merge_sciMETN9.CGmap.gz | cgmaptools mfg -r fragreg_hg38GRCm39_genes.bed -c 1 -x CG > merge_sciMETN9.CG_genes.mfg

gunzip -c merge_sciEM.CGmap.gz | cgmaptools mfg -r fragreg_hg38GRCm39_genes.bed -c 1 -x CG > merge_sciEM_total.CG_genes.mfg
gunzip -c merge_sciMET.CGmap.gz | cgmaptools mfg -r fragreg_hg38GRCm39_genes.bed -c 1 -x CG > merge_sciMET_total.CG_genes.mfg

# Plot distributions accross bed regions
(head -1 merge_sciEM_total.CG_genes.mfg | gawk '{$1="Sample"; print $0;}';
 for F in *total.CG_genes.mfg; do
   awk -v SampleName=$(echo $F | sed s/.mfg//g) '/total_ave_mC/{$1=SampleName; print $0;}' $F
 done
) > mfg_merge.xls

cgmaptools fragreg -i mfg_merge.xls -o merge.fragreg.pdf -f pdf







awk -v SampleName='echo $F' $F | sed s/.mfg//g`

awk -v SampleName=$(echo $F | sed s/.mfg//g) '/total_ave_mC/{$1=SampleName; print $0;}' $F


awk -v SampleName='echo $F | sed s/.mfg//g' '/total_ave_mC/{$1=SampleName; print $0;}' $F

