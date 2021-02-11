library(openxlsx)
library(ggplot2)
library(gridExtra)

setwd("/Users/zacc/USyd/scWGBS/results")

# read in sequencing read tracking information
dat <- read.xlsx("sci_read_tracking.xlsx", sheet ="read_tacking_new")

# rearrange data
m <- melt(dat)
m$library <- as.factor(m$library)

# plotting 
for(i in levels(m$library)) {
  p <- ggplot(subset(m, library==i), aes(variable, value,  fill = variable)) + 
    facet_wrap(~ library) +
    geom_bar(stat="identity", show_guide=FALSE)
  
  ggsave(paste0("figure_",i,".pdf"), p)
}


id1<- c("reads_demultiplexed")
g1 <- ggplot(subset(m, library==id1), aes(variable, value,  fill = variable)) + 
  facet_wrap(~ library) + xlab("") + ylab ("reads demultiplexed") +
  geom_bar(stat="identity", show_guide=FALSE) + theme_bw() + theme_minimal()

id1<- c("reads_demultiplexed_asprc_loading")
g2 <- ggplot(subset(m, library==id1), aes(variable, value,  fill = variable)) + 
  facet_wrap(~ library) + xlab("") + ylab ("reads (% loading)") +
  geom_bar(stat="identity", show_guide=FALSE) + theme_bw() + theme_minimal()

id1<- c("reads_posttrim")
g3 <- ggplot(subset(m, library==id1), aes(variable, value,  fill = variable)) + 
  facet_wrap(~ library) + xlab("") + ylab ("reads post-trim") +
  geom_bar(stat="identity", show_guide=FALSE) + theme_bw() + theme_minimal()

id1<- c("reads_posttrim_asprc_demultiplexed")
g4 <- ggplot(subset(m, library==id1), aes(variable, value,  fill = variable)) + 
  facet_wrap(~ library) + xlab("") + ylab ("reads post-trim (% demultiplexed)") +
  geom_bar(stat="identity", show_guide=FALSE) + theme_bw() + theme_minimal()

# write plots to file
pdf("barplots_scwgbs_pre_trim2.pdf")
grid.arrange(g1,g2,g3,g4,nrow=2, ncol=2)
dev.off()


id1<- c("reads_posttrim_prcwith_GAAGAGCACACGTCTGAACTCC")
g1 <- ggplot(subset(m, library==id1), aes(variable, value,  fill = variable)) + 
  facet_wrap(~ library) + xlab("") + ylab ("reads w adapter (% of post-trim)") +
  geom_bar(stat="identity", show_guide=FALSE) + theme_bw() + theme_minimal()

id1<- c("reads_posttrim2_asprc_demultiplexed")
g2 <- ggplot(subset(m, library==id1), aes(variable, value,  fill = variable)) + 
  facet_wrap(~ library) + xlab("") + ylab ("reads post-trim2 (% demultiplexed)") +
  geom_bar(stat="identity", show_guide=FALSE) + theme_bw() + theme_minimal()

# write plots to file
pdf("barplots_scwgbs_trimtrim2.pdf")
grid.arrange(g1,g2,nrow=2, ncol=2)
dev.off()


id1<- c("reads_mapped_total_posttrim2")
g1 <- ggplot(subset(m, library==id1), aes(variable, value,  fill = variable)) + 
  facet_wrap(~ library) + xlab("") + ylab ("reads mapped") +
  geom_bar(stat="identity", show_guide=FALSE) + theme_bw() + theme_minimal()

id1<- c("reads_posttrim2_asprc_demultiplexed")
g2 <- ggplot(subset(m, library==id1), aes(variable, value,  fill = variable)) + 
  facet_wrap(~ library) + xlab("") + ylab ("reads mapped (% demultiplexed)") +
  geom_bar(stat="identity", show_guide=FALSE) + theme_bw() + theme_minimal()

id1<- c("reads_mapped_noncontrolgenomes_asprc_trimmed2")
g3 <- ggplot(subset(m, library==id1), aes(variable, value,  fill = variable)) + 
  facet_wrap(~ library) + xlab("") + ylab ("reads mapped non-ctr (% post-trim2)") +
  geom_bar(stat="identity", show_guide=FALSE) + theme_bw() + theme_minimal()


# write plots to file
pdf("barplots_scwgbs_align.pdf")
grid.arrange(g1,g2,g3,nrow=2, ncol=2)
dev.off()


