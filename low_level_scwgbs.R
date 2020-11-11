library(devtools)
library(ggplot2)
library(easyGgplot2)

setwd("/Users/zacc/USyd/scWGBS/results/low_level_seq")

# read in bam summary data
dat<-read.csv(file="clusters.csv")

# summarise for experiment
dat$method<-gsub("_.*","",dat$Experiment)

# barcodes detected
dat$detected<-dat$reads_total > 0

g1<-ggplot(dat, aes(fill=detected, y=rep(1,nrow(dat)), x=method)) + 
  geom_bar(position="stack", stat="identity") + ylab("barcodes detected")

# get cells with reads
dplot<-dat[dat$reads_total > 0,]

# plot reads/cell
g2<-ggplot(data=dplot, aes(x=method,y=log(reads_total+1), fill=method)) + 
  geom_point(position = position_jitterdodge(),alpha=0.5, aes(colour = method)) + 
  geom_boxplot(alpha=0.6,outlier.alpha=0,) + ylab("log(reads/cell + 1)") + theme_bw() + xlab("") + 
  theme(legend.position="none")

# plot hg38 aligned reads/cell
g3<-ggplot(data=dat, aes(x=method,y=reads_hg38/reads_total, fill=method)) + 
  geom_point(position = position_jitterdodge(),alpha=0.5, aes(colour = method)) + 
  geom_boxplot(alpha=0.6,outlier.alpha=0,) + ylab("reads hg38 (% total))") + theme_bw() + xlab("") + 
  theme(legend.position="none")

# plot GRCm39 aligned reads/cell
g4<-ggplot(data=dat, aes(x=method,y=reads_GRCm39/reads_total, fill=method)) + 
  geom_point(position = position_jitterdodge(),alpha=0.5, aes(colour = method)) + 
  geom_boxplot(alpha=0.6,outlier.alpha=0,) + ylab("reads GRCm39 (% total)") + theme_bw() + xlab("") + 
  theme(legend.position="none")

# reads mapped human v mouse
g5<-ggplot(data=dat, aes(x=log(reads_GRCm39+1),y=log(reads_hg38 +1), fill=method)) + 
  geom_point(alpha=0.5, aes(colour = method)) + 
  ylab("log(hg38 reads/cell + 1)") + theme_bw() + xlab("log(GRCm39 reads/cell + 1)") + 
  theme(legend.position="none")

# reads mapping to experimental controls
g6<-ggplot(data=dat, aes(x=lambda,y=log(reads_Lambda+1), fill=method)) + 
  geom_point(position = position_jitterdodge(),alpha=0.5, aes(colour = method)) + 
  geom_boxplot(alpha=0.6,outlier.alpha=0,) + ylab("log(Lambda reads/cell + 1)") + theme_bw() + xlab("") + 
  theme(legend.position="none")

g7<-ggplot(data=dat, aes(x=pUC19,y=log(reads_pUC19+1), fill=method)) + 
  geom_point(position = position_jitterdodge(),alpha=0.5, aes(colour = method)) + 
  geom_boxplot(alpha=0.6,outlier.alpha=0,) + ylab("log(pUC19 reads/cell + 1)") + theme_bw() + xlab("") + 
  theme(legend.position="none")

# write plots to file
pdf("plot_adj_fqnorm_unn_posnegbox_all_breachers.pdf")
grid.arrange(g1,g2,g3,g4,g5,g6,g7,nrow=2, ncol=2)
dev.off()
