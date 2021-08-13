library(readxl)

#setwd
setwd("/Users/zacc/USyd/scWGBS")

# read in barcodes
# olig <- read_excel("/Users/zacc/USyd/scWGBS/scimet_oligos_idt.xlsx", sheet = "oligos")
olig <- read.delim("/Users/zacc/github_repo/scimet_scripts/data/barcodes_scimet_10_11_10_NovaSeq.txt", sep='\t', header=F)
colnames(olig) <- c("Name","index_number","Index")

# read in experiment layout
edat <- read_excel("/Users/zacc/USyd/scWGBS/Lab_130421/sci_130421_explayouts.xlsx", sheet = "experiment_file")

# make all unique combinations of i7, Tn5 + i5
# Novaseq indexing combinations = i7[10bp,I1],Tn5[11bp,I2],i5[10bp,I3]
tn5 <- olig[grep("sciMET_Tn5", olig$Name),]
i5 <- olig[grep("sciMET_i5", olig$Name),]
i7 <- olig[grep("sciMET_i7", olig$Name),]

res <- expand.grid(i7$Index,tn5$Index,i5$Index)
index_sequence <- paste0(res$Var1,res$Var2,res$Var3)
res2 <- expand.grid(i7$Name,tn5$Name,i5$Name)

index_df <- as.data.frame(cbind(res,res2,index_sequence))
colnames(index_df) <- c("i7_sequence","tn5_sequence","i5_sequence",
                        "i7_name","tn5_name","i5_name",
                        "index_sequence")
index_df$sciMET_i5.i7 <- paste0(index_df$i5_name,sep=".", index_df$i7_name)

# merge experiment data
tn5 <- edat[,grep("Tn5", colnames(edat))]
id2 <- merge(index_df,tn5,by.x="tn5_name", by.y="sciMET_Tn5")

i5i7 <- edat[,grep("i5.i7", colnames(edat))]
id3 <- merge(id2,i5i7,by="sciMET_i5.i7")

# select barcodes used within experiment
#index_used <- index_df[index_df$i5_index %in% paste0("sciMET_i5_",c(1:8)) &
#                       index_df$index1_name %in% paste0("sciMET_i7_",c(1:10)), ]
# remove unused barcodes
#index_final <- index_used[-which(index_used$i5_index %in% paste0("sciMET_i5_",c(3:8)) &
#                           index_used$index1_name %in% paste0("sciMET_i7_",c(10))), ]

# write to files
index_final <- id3
write.table(index_final,file="sci_experiment.txt",sep='\t', quote=FALSE, row.names=FALSE)
