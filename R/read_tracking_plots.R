library(openxlsx)
library(ggplot2)
library(gridExtra)

setwd("/Users/zacc/USyd/scWGBS/results")

# read in sequencing read tracking information
dat <- read.xlsx("sci_read_tracking.xlsx", sheet ="read_tacking_new")
colnames(dat) <-c("library","sciMET","sciMET(mg)","sciEM(n9)","sciEM(mg)")

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

# plots for sciMET, sciEMn9 and sciEMmg
# rearrange data
m <- melt(dat[,-3])
m$cols1 <- as.character(m$variable)
m$cols1[m$cols1 == "sciMET" ] <- c("grey")
m$cols1[m$cols1 == "sciEM(n9)" ] <- c("cyan3")
m$cols1[m$cols1 == "sciEM(mg)" ] <- c("dodgerblue")
m$library <- as.factor(m$library)
m$variable <- factor(m$variable, levels=c("sciMET","sciEM(n9)","sciEM(mg)"))

id1<- c("reads_demultiplexed")
g1 <- ggplot(subset(m, library==id1), aes(variable, value,  fill = variable)) + 
  facet_wrap(~ library) + xlab("") + ylab ("Reads demultiplexed (M)") + scale_fill_manual(values=c("grey","cyan3","dodgerblue")) +
  geom_bar(stat="identity", show_guide=FALSE) + theme_bw() + theme_minimal() + scale_x_discrete(guide = guide_axis(angle = 45))

id1<- c("reads_demultiplexed_asprc_loading")
g2 <- ggplot(subset(m, library==id1), aes(variable, value,  fill = variable)) + 
  facet_wrap(~ library) + xlab("") + ylab ("Reads demultiplexed (% of loading)") + scale_fill_manual(values=c("grey","cyan3","dodgerblue")) +
  geom_bar(stat="identity", show_guide=FALSE) + theme_bw() + theme_minimal() + scale_x_discrete(guide = guide_axis(angle = 45))

dplot<-dat[dat$library  %in% c("plot_reads_trim_demult_prc","plot_reads_trim2_demult_prc",
                               "plot_reads_unmapped_demult_prc","plot_reads_mapped_demult_prc"),c(1,2,4,5)]
m <- melt(dplot)
m$library <- factor(m$library, levels=c("plot_reads_trim_demult_prc","plot_reads_trim2_demult_prc",
                                        "plot_reads_unmapped_demult_prc","plot_reads_mapped_demult_prc"))

g3 <- ggplot(m, aes(x=variable, y=value,fill=library)) + 
  geom_col() + scale_fill_grey()+ ylab("Reads (% demultiplexed)") + xlab("") +
  scale_y_continuous(label = scales::percent) + theme_bw() + theme_minimal() + theme(legend.position = "top")

# write plots to file
pdf("barplots_scwgbs_n3.pdf")
grid.arrange(g1,g2,nrow=2, ncol=4)
grid.arrange(g3,nrow=2, ncol=2)
dev.off()
