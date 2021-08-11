library(devtools)
library(ggplot2)
library(dplyr)
library(easyGgplot2)
library(factoextra)
library(openxlsx)
library(inflection)

setwd("/Users/zacc/USyd/scWGBS/results/low_level_seq")

dat <- read.xlsx("clusters_280121.xlsx")
dat$library <- gsub("_.*","",dat$Experiment)

# plot total reads and percentage of cross-over reads
dplot <- dat[dat$library == "sciEMMG",]
ggplot(data=dplot, aes(y=cross_over_prc,x=log10(reads_origin), fill=library)) + 
  geom_point(alpha=0.5, aes(colour = library)) + 
  ylab("cross-over (%)") + theme_bw() + xlab("log10(reads)") + geom_density_2d() + ylim(0,100)


dplot <- dat[dat$library == "sciMETN9",]
g1 <- ggplot(data=dplot, aes(y=cross_over_prc,x=log10(reads_origin), fill=library)) + 
  geom_point(alpha=0.5,color = "red") + geom_density_2d(color = "black") + ggtitle("sciMETN9") +
  ylab("cross-over (%)") + theme_bw() + xlab("log10(reads)") + ylim(0,100) + theme(legend.position="none") 

dplot <- dat[dat$library == "sciMETMG",]
g2 <- ggplot(data=dplot, aes(y=cross_over_prc,x=log10(reads_origin), fill=library)) + 
  geom_point(alpha=0.5,color = "red") + geom_density_2d(color = "black") + ggtitle("sciMETMG") +
  ylab("cross-over (%)") + theme_bw() + xlab("log10(reads)") + ylim(0,100) + theme(legend.position="none") 

dplot <- dat[dat$library == "sciEMN9",]
g3 <- ggplot(data=dplot, aes(y=cross_over_prc,x=log10(reads_origin), fill=library)) + 
  geom_point(alpha=0.5,color = "red") + geom_density_2d(color = "black") + ggtitle("sciEMN9") +
  ylab("cross-over (%)") + theme_bw() + xlab("log10(reads)") + ylim(0,100) + theme(legend.position="none") 

dplot <- dat[dat$library == "sciEMMG",]
g4 <- ggplot(data=dplot, aes(y=cross_over_prc,x=log10(reads_origin), fill=library)) + 
  geom_point(alpha=0.5,color = "red") + geom_density_2d(color = "black") + ggtitle("sciEMMG") +
  ylab("cross-over (%)") + theme_bw() + xlab("log10(reads)") + ylim(0,100) + theme(legend.position="none") 

# save plots
pdf(file="crossover_scwgbs.pdf")
grid.arrange(g1,g2,g3,g4,nrow=2, ncol=2)
dev.off()

# Histogram with density plot
ggplot(dat, aes(x=cross_over_prc)) + 
  geom_histogram(aes(y=..density..),binwidth=0.4, colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") + theme_bw() + xlab("cross-over (%)") +
  geom_vline(aes(xintercept=mean(weight)),
             color="blue", linetype="dashed", size=1)

# Find Extremum Distance Estimator and use as cutoff for allowable cross-over %
EME <- bede(dat$cross_over_prc[!is.na(dat$cross_over_prc)],c(1:length(dat$cross_over_prc[!is.na(dat$cross_over_prc)])), 0)$iplast

# Histogram with density plot
g1 <- ggplot(dat, aes(x=cross_over_prc)) + 
  geom_histogram(aes(y=..density..),binwidth=0.4, colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") + theme_bw() + xlab("cross-over (%)") +
  geom_vline(aes(xintercept=EME), color="black", linetype="dashed", size=0.5)

# save plots
pdf(file="crossover_density_scwgbs.pdf")
grid.arrange(g1,nrow=2, ncol=2)
dev.off()

# plot reads/cell
g2<-ggplot(data=dat, aes(x=library,y=cross_over_prc, fill=library)) + 
  geom_point(position = position_jitterdodge(),alpha=0.5, aes(colour = library)) + 
  geom_boxplot(alpha=0.6,outlier.alpha=0,) + ylab("cross-over (%)") + theme_bw() + xlab("") + 
  theme(legend.position="none")

pdf(file="boxplot_crossover_scwgbs.pdf")
grid.arrange(g2,nrow=2, ncol=2)
dev.off()


# select nuclei that pass EME threshold
dat2 <- dat[which(dat$cross_over_prc < EME),]

# Find Extremum Distance Estimator and use as cutoff for allowable cross-over %
EME2 <- bede(log10(dat2$reads_origin),c(1:length(dat2$reads_origin)), 0)$iplast

g1 <- ggplot(dat2, aes(x=log10(reads_origin))) + 
  geom_histogram(aes(y=..density..),binwidth=0.2, colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") + theme_bw() + xlab("log10(reads)")

# save plots
pdf(file="reads_origin_passEME_density_scwgbs.pdf")
grid.arrange(g1,nrow=2, ncol=2)
dev.off()

# select nulcei with <17% cross-over and > 10 reads fro clustering
R10_CO17 <- dat$cross_over_prc < 17 & dat$reads_origin > 10
R10_CO17[is.na(R10_CO17)] <- "FALSE"
dat$R10_CO17 <- R10_CO17

nuclei_select <- cbind(dat$Cluster,dat$origin,R10_CO17)
nuclei_select <- nuclei_select[nuclei_select[,3] == "TRUE",]
write.table(nuclei_select, file="nuclei_select.txt",sep='\t', col.names=F, row.names=F, quote = F)

# barplot of human mouse of good quality nuclei
dat2 <- dat[which(dat$R10_CO17 == "TRUE"),]

ggplot(data=dat2, aes(x=library, y=origin, fill=origin)) +
  geom_bar(stat="identity")

g1 <- ggplot(dat2, aes(fill=origin, y=rep(1,nrow(dat2)), x=library)) + 
  geom_bar(position="stack", stat="identity") + ylab("nuclei") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

g2 <- ggplot(dat2, aes(fill=origin, y=rep(1,nrow(dat2)), x=library)) + 
  geom_bar(position="fill", stat="identity") + ylab("fraction of nuclei") + theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# save plots
pdf(file="barplots_origin_passEMEperlib_scwgbs.pdf")
grid.arrange(g1,g2,nrow=2, ncol=2)
dev.off()

# reads v cross over good quality nuclei
dplot <- dat[dat$R10_CO17 == "TRUE",]
ggplot(data=dplot, aes(y=cross_over_prc,x=log10(reads_origin), fill=library)) + 
  geom_point(alpha=0.5, aes(colour = library)) + 
  ylab("cross-over (%)") + theme_bw() + xlab("log10(reads)") + geom_density_2d() + ylim(0,25)

# plot reads/cell
g2<-ggplot(data=dplot, aes(x=library,y=cross_over_prc, fill=library)) + 
  geom_point(position = position_jitterdodge(),alpha=0.5, aes(colour = library)) + 
  geom_boxplot(alpha=0.6,outlier.alpha=0,) + ylab("cross-over (%)") + theme_bw() + xlab("") + 
  theme(legend.position="none")

t.test(dplot$cross_over_prc[dplot$library %in% c("sciEMMG","sciEMN9")],
       dplot$cross_over_prc[dplot$library %in% c("sciMETMG","sciMETN9")])

mean(dplot$cross_over_prc[dplot$library %in% c("sciEMMG")])
