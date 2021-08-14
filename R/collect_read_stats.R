#!/usr/bin/env Rscript
### collect expeirment stats and bind to experiment sheet from experiment_barcode_combine.R
# local data for testing
# dat <- read.delim("/Users/zacc/USyd/scWGBS/sci_experiment.txt", sep='\t', header=T)
# trim_fq="/Users/zacc/USyd/scWGBS/bioinformatics/statcollect_test"

# experiment sheet produced from experiment_barcode_combine.R
dat <- read.delim("/project/RDS-FMH-FFEPIGENETICS-RW/scwgbs/210804_A00152_0453_AHHGWYDRXY/sci_experiment.txt", sep='\t', header=T)

# output dir
outdir="/project/RDS-FMH-FFEPIGENETICS-RW/scwgbs/210804_A00152_0453_AHHGWYDRXY"
  
# samtools pseudo stats dir
pseudostats="/project/RDS-FMH-FFEPIGENETICS-RW/scwgbs/210804_A00152_0453_AHHGWYDRXY/fq_split"

# scBS-map.report dir
scbsmap="/project/RDS-FMH-FFEPIGENETICS-RW/scwgbs/210804_A00152_0453_AHHGWYDRXY/fq_split"

# fastq dirs
trim_fq="/project/RDS-FMH-FFEPIGENETICS-RW/scwgbs/210804_A00152_0453_AHHGWYDRXY/fq_split"
trim1_fq="/project/RDS-FMH-FFEPIGENETICS-RW/scwgbs/210804_A00152_0453_AHHGWYDRXY/trim1_fq_split"
raw_fq="/project/RDS-FMH-FFEPIGENETICS-RW/scwgbs/210804_A00152_0453_AHHGWYDRXY/demult_fq_split"

############################
### read count on .fastq ###
############################
# trimmed .fastq
setwd(trim_fq)
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

# trim1 .fastq
setwd(trim1_fq)
file_name <- c("_trimmed_R1.fq")

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

# raw demultiplexed .fastq
setwd(raw_fq)
file_name <- c("_demultiplex_R1.fq")

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

# "scBS-map.report" R1
setwd(scbsmap)
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

# "scBS-map.report" R2
setwd(scbsmap)
file_name <- c("_R2.scBS-map.report")

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
# note 1. that pseudopaired.bam files are created using pseudopair_bam.sh run on R1 & R2 that have been aligned (each SE) by scbsmap 
# note 2. stats are taken using samtools stats command
setwd(pseudostats)
file_name <- c("_pseudopaired.stats")

col_name <- paste0("reads",file_name)
file_names <- paste0(dat$index_sequence,file_name)[paste0(dat$index_sequence,file_name) %in% list.files()]
system(paste("cat",file_names[1],"| grep ^SN | cut -f 2- > tmp"))
tmp <- read.delim("tmp", sep='\t',header=F)
tmp <- tmp[grep(":",tmp$V1),]
col_names<- paste0(as.character(t(data.frame(do.call('rbind', strsplit(as.character(tmp$V1),':',fixed=TRUE)))[,1]))," ",file_name)
file_position <- which(paste0(dat$index_sequence,file_name) %in% list.files())
dat <- cbind(dat,matrix("NA",nrow(dat),length(col_names)))
colnames(dat) <- c(colnames(dat)[1:(ncol(dat)-length(col_names))],col_names)

for (i in 1:length(file_position)){
  print(paste0(round(i/length(file_position) * 100), "% done"))
  
  system(paste("cat",file_names[i],"| grep ^SN | cut -f 2- > tmp"))
  tmp <- read.delim("tmp", sep='\t',header=F)
  tmp <- tmp[grep(":",tmp$V1),]
  tmp2 <- cbind(dat[file_position[i], which(!colnames(dat) %in% col_names)],
                t(as.character(tmp$V2)))
  colnames(tmp2) <- colnames(dat)
  dat[file_position[i],] <- tmp2
}

# write to file
write.table(dat,file = paste0(outdir,"/sci_experiment_collected.txt"), sep='\t', header=T)