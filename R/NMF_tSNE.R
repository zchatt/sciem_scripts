library(Rtsne)
#install.packages("/Users/zacc/Downloads/NNLM_0.4.2.tar.gz", repos = NULL, type="source")
library(NNLM)
library(ggplot2)
library(openxlsx)
library(gridExtra)
library(PerformanceAnalytics)

setwd("/Users/zacc/USyd/scWGBS/results/nuclei_select")

# 10 is an example number, need to select the optimal number in practice...
nnmf_k <- 10

# dataset: Rows are features and columns are samples/cells...
dataset <- read.delim(file="/Users/zacc/USyd/scWGBS/results/nuclei_selectmean_mC", sep='\t', header=T)
row.names(dataset) <- as.character(unlist(as.data.frame(strsplit(dataset$CW," "))[2,])) 
dataset$CW <- as.character(unlist(as.data.frame(strsplit(dataset$CW," "))[1,])) 

# sample_label: group information of samples...
dat <- read.xlsx("/Users/zacc/USyd/scWGBS/results/low_level_seq/clusters_280121.xlsx")
dat$library <- gsub("_.*","",dat$Experiment)
row.names(dat) <- paste0("cluster",dat$Cluster)
dat$log10_reads_origin <- log10(dat$reads_origin)

dat <- merge(dat,dataset,by="row.names")
dat <- dat[dat$library %in% c("sciMETN9","sciEMMG"),]

# create data matrix
d2 <- t(as.matrix(dat[,c("CG","CA","CC","CT")]))

# get list of labels you want to overlay
lab_list <- c("library","log10_reads_origin","origin","CA","CG","CH")

for (i in 1:length(lab_list)){
sample_label <- dat[,lab_list[i]]
    
set.seed(100)
#data_nmf <- nnmf(d2, nnmf_k, verbose = FALSE)                     
data_nmf <- nnmf(d2, nnmf_k, verbose = FALSE,check.k = FALSE)                     
#tsne_out <- Rtsne(t(data_nmf$H), theta = 0.5, perplexity = 30, check_duplicates = FALSE)

data_nmf_clean <-t(data_nmf$H)
sample_label <- sample_label[complete.cases(data_nmf_clean)]
data_nmf_clean <- data_nmf_clean[complete.cases(data_nmf_clean),]

tsne_out <- Rtsne(data_nmf_clean, theta = 0.5, perplexity = 30, check_duplicates = FALSE)
tsne_xy  <- as.data.frame(tsne_out$Y)
plotfile <- cbind(c(1:nrow(data_nmf_clean)), tsne_xy, sample_label)

colnames(plotfile) <- c('sample_id','tSNE1','tSNE2','sample_label')
assign(paste0("g",i),ggplot(plotfile, aes(x = tSNE1, y = tSNE2, color = sample_label)) + geom_point(cex=2))

}

# write plots to file
pdf("tSNE_Ccontext_scwgbs.pdf")
grid.arrange(g1,g2,g3,g4,g5,g6,nrow=3, ncol=2)
dev.off()


## check correlations of meth and reads
d2 <- dat[,c("reads_total","reads_hg38","reads_GRCm39","reads_Lambda","reads_pUC19","reads_realigned_onlyhg38",
             "reads_realigned_onlymm10","reads_origin","crossover_reads","cross_over_prc","log10_reads_origin",
             "C","CG","CHG","CHH","CA","CC","CT","CH")]

chart.Correlation(d2, histogram=TRUE, pch=19)
