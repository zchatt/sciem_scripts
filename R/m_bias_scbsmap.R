library(ggplot2)
library(gridExtra)
library(reshape)
library(zoo)

setwd("/Users/zacc/USyd/scWGBS/results/low_level_seq")

# BSseeker DNA methylation output
#X: methylated CG
#x: un-methylated CG
#Y: methylated CHG
#y: un-methylated CHG
#Z: methylated CHH
#z: un-methylated CHH

# scbsmap .bam files were firstly converted to .sam (samtools view)

#read in entire SAM file
sam_files <- list.files(pattern=".sam$")
for ( s1 in 1:length(sam_files)){
sam <- read.delim(file=sam_files[s1],header=F,nrows=10000000)

# remove duplicate sequences
sam <- sam[!duplicated(sam$V10),]

# remove lambda and pUC19
sam <- sam[sam$V3 != "Lambda_NEB",]
sam <- sam[sam$V3 != "pUC19",]

# convert meth calls to dataframe
m1 <- gsub("XM:Z:","",sam$V15)
n_reads <- c(1000000)
if (length(m1) < n_reads) {
  n_reads = length(m1)
}
m1<-m1[1:n_reads]
tmp <- data.frame(do.call(rbind, strsplit(m1, "", fixed=TRUE)))

# break into 1K chunks
res<-tmp[1:2,]
n_split <- 1000
for (z in 1:(n_reads/n_split)){
  z1 <- z*n_split - n_split + 1
  z2 <- (z+1)*n_split - n_split
  print(z1)
  tmpe<-tmp[z1:z2,]
  for (i in 1:nrow(tmpe)){
    ni <- c(nchar(m1[(z*n_split)+i-n_split])):ncol(tmpe)
    tmpe[i,ni]<-rep(NA,length(ni))
  }
  res<-rbind(res,tmpe)
}
tmp <- res[-c(1,2),]

fac_cols <- sapply(tmp, is.factor)                        
tmp[fac_cols] <- lapply(tmp[fac_cols], as.character) 
tmp <- data.frame(lapply(tmp, function(x){gsub("-", NA, x)}))

# CG
tmp2 <- data.frame(lapply(tmp, function(x){gsub("x",0, x)}))
tmp2 <- data.frame(lapply(tmp2, function(x){gsub("X",1, x)}))
tmp2 <- data.frame(lapply(tmp2, function(x){gsub("y",NA, x)}))
tmp2 <- data.frame(lapply(tmp2, function(x){gsub("Y",NA, x)}))
tmp2 <- data.frame(lapply(tmp2, function(x){gsub("z",NA, x)}))
tmp2 <- data.frame(lapply(tmp2, function(x){as.numeric(gsub("Z",NA, x))}))

cg <- colMeans(tmp2,na.rm=T)*100

# CHG
tmp2 <- data.frame(lapply(tmp, function(x){gsub("x",NA, x)}))
tmp2 <- data.frame(lapply(tmp2, function(x){gsub("X",NA, x)}))
tmp2 <- data.frame(lapply(tmp2, function(x){gsub("y",0, x)}))
tmp2 <- data.frame(lapply(tmp2, function(x){gsub("Y",1, x)}))
tmp2 <- data.frame(lapply(tmp2, function(x){gsub("z",NA, x)}))
tmp2 <- data.frame(lapply(tmp2, function(x){as.numeric(gsub("Z",NA, x))}))

chg <- colMeans(tmp2,na.rm=T)*100

# CHH
tmp2 <- data.frame(lapply(tmp, function(x){gsub("x",NA, x)}))
tmp2 <- data.frame(lapply(tmp2, function(x){gsub("X",NA, x)}))
tmp2 <- data.frame(lapply(tmp2, function(x){gsub("y",NA, x)}))
tmp2 <- data.frame(lapply(tmp2, function(x){gsub("Y",NA, x)}))
tmp2 <- data.frame(lapply(tmp2, function(x){gsub("z",0, x)}))
tmp2 <- data.frame(lapply(tmp2, function(x){as.numeric(gsub("Z",1, x))}))

chh <- colMeans(tmp2,na.rm=T)*100

# plot m-bias
dplot <- melt(cbind(cg,chg,chh))
dplot$X1 <- as.numeric(gsub("X","",dplot$X1))
colnames(dplot)<-c("position","context","meth")

p<-ggplot(dplot, aes(x=position, y=meth, group=context,colour = context)) +
  geom_line(aes(y=rollmean(meth, 7, na.pad=TRUE))) +
  xlab("Position") + ylab("DNA methylation %") + ylim(0,100) +
  theme_bw() + ggtitle(paste0(sam_files[s1],n_reads," reads"))

# print to file
pdf(file=paste0("mbias_scbsmap_",sam_files[s1],".pdf"))
grid.arrange(p,nrow=2, ncol=2)
dev.off()

write.table(dplot,file=paste0("mbias_scbsmap_demultiplexed",sam_files[s1],".txt"),sep='\t')
}
