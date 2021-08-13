library(R.utils)
library(miscTools)

# read in experiment sheet produced from experiment_barcode_combine.R
#dat <- read.delim("/Users/zacc/USyd/scWGBS/sci_experiment.txt", sep='\t', header=T)
dat <- read.delim("/project/RDS-FMH-FFEPIGENETICS-RW/scwgbs/210804_A00152_0453_AHHGWYDRXY/sci_experiment.txt", sep='\t', header=T)

# directory of the .fastq files with names same as index_sequence_R1.fq
#fq_indir="/Users/zacc/USyd/scWGBS/bioinformatics/statcollect_test"
fq_indir="/project/RDS-FMH-FFEPIGENETICS-RW/scwgbs/210804_A00152_0453_AHHGWYDRXY/fq_split"
setwd(fq_indir)

############################
### read count on .fastq ###
############################
file_name <- c("_R1.fq")

col_name <- paste0("reads",file_name)
file_names <- paste0(dat$index_sequence,file_name)[paste0(dat$index_sequence,file_name) %in% list.files()]
file_position <- which(paste0(dat$index_sequence,file_name) %in% list.files())
dat[,col_name] <- rep("NA",nrow(dat))
for (i in 1:length(file_position)){
  print(paste0(round(i/length(file_position) * 100), "% done"))
  dat[,col_name][file_position[i]] <- as.numeric(countLines(file_names[i]))/ 4
}
file_names <- paste0(dat$index_sequence,"_R1.fq")[paste0(dat$index_sequence,"_R1.fq") %in% list.files()]
file_position <- which(paste0(dat$index_sequence,"_R1.fq") %in% list.files())

#################################
### read count on aligned.bam ### 
#################################
# # _R1.bam
# file_name <- c("_R1.bam")
# 
# col_name <- paste0("reads",file_name)
# file_names <- paste0(dat$index_sequence,file_name)[paste0(dat$index_sequence,file_name) %in% list.files()]
# file_position <- which(paste0(dat$index_sequence,file_name) %in% list.files())
# dat[,col_name] <- rep("NA",nrow(dat))
# for (i in 1:length(file_position)){
#   print(paste0(round(i/length(file_position) * 100), "% done"))
#   dat[,col_name][file_position[i]] <- as.numeric(countLines(file_names[i]))
# }
# 
# # _R1.end2end.bam
# file_name <- c("_R1.end2end.bam")
# 
# col_name <- paste0("reads",file_name)
# file_names <- paste0(dat$index_sequence,file_name)[paste0(dat$index_sequence,file_name) %in% list.files()]
# file_position <- which(paste0(dat$index_sequence,file_name) %in% list.files())
# dat[,col_name] <- rep("NA",nrow(dat))
# for (i in 1:length(file_position)){
#   print(paste0(round(i/length(file_position) * 100), "% done"))
#   dat[,col_name][file_position[i]] <- as.numeric(countLines(file_names[i]))
# }
# 
# # _R1.local.bam
# file_name <- c("_R1.local.bam")
# 
# col_name <- paste0("reads",file_name)
# file_names <- paste0(dat$index_sequence,file_name)[paste0(dat$index_sequence,file_name) %in% list.files()]
# file_position <- which(paste0(dat$index_sequence,file_name) %in% list.files())
# dat[,col_name] <- rep("NA",nrow(dat))
# for (i in 1:length(file_position)){
#   print(paste0(round(i/length(file_position) * 100), "% done"))
#   dat[,col_name][file_position[i]] <- as.numeric(countLines(file_names[i]))
# }

# "scBS-map.report"
file_name <- c("_R1.scBS-map.report")

col_name <- paste0("reads",file_name)
file_names <- paste0(dat$index_sequence,file_name)[paste0(dat$index_sequence,file_name) %in% list.files()]
tmp <- read.delim(file_names[1], sep='\t',header=F,skip=3)
col_names<- paste0(as.character(t(data.frame(do.call('rbind', strsplit(as.character(tmp$V1),': ',fixed=TRUE)))$X1[1:14]))," ",file_name)
file_position <- which(paste0(dat$index_sequence,file_name) %in% list.files())
dat <- cbind(dat,matrix("NA",nrow(dat),length(col_names)))
colnames(dat) <- c(colnames(dat)[1:(ncol(dat)-length(col_names))],col_names)

for (i in 1:length(file_position)){
  print(paste0(round(i/length(file_position) * 100), "% done"))
  tmp <- read.delim(file_names[i], sep='\t',header=F,skip=3)
  tmp2 <- cbind(dat[file_position[i], which(!colnames(dat) %in% col_names)],
                t(as.data.frame(t(gsub("%","",t(data.frame(do.call('rbind', strsplit(as.character(tmp$V1),': ',fixed=TRUE)))$X2[1:14]))))))
  colnames(tmp2) <- colnames(dat)
  dat[file_position[i],] <- tmp2
  }

# samtools stats on pseudopaired.bam

# write to file
write.table(dat,file = "sci_experiment_collected.txt", sep='\t', header=T)