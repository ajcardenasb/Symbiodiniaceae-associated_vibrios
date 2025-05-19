setwd("~/Documents/Bioinformatics_scripts/R_scripts/Fractionation_optimization/")
library(dplyr)
library(stringr)


##################################################################
#####Identifying and removing contaminant ASVs normalized data#####
###################################################################

asv=read.table("inputs/FractionationOptimization_ASV_table.txt", header = TRUE, row.names = 1, sep = " ")[,1:50]
tax=read.table("inputs/FractionationOptimization_ASV_table.txt", header = TRUE, row.names = 1, sep = " ")[,51:58]
met=read.table("Inputs/FractionationOptimization_metafile.txt", header = TRUE,row.names = 1, sep = "\t")

# all samples have more than 10000 reads
options(scipen=5)
hist(colSums(asv),  breaks = 50, labels = F)
axis(side=1, at=seq(0,5000, 1000))
colSums(asv)
names(asv)

colnames(asv[, colSums(asv) < 5000]) # 13 samples under 5000, mainly water samples
asv.o=asv[, colSums(asv) > 1000]
#message(ncol(asv.o)," samples with > 5000 reads were retained out of ", ncol(asv), " total samples")

#Identify and removing contaminant ASVs raw data
asv.r=as.data.frame(sweep(asv.o,2,colSums(asv.o),`/`))

asv.r$average_negatives=rowMeans(asv.r[,colnames(asv.r) %in% c("Neg.1","Neg.2")])
asv.r$average_samples=rowMeans(asv.r[, !colnames(asv.r) %in% c("Neg.1","Neg.2")])
asv.r$contaFactor=(asv.r$average_negatives/asv.r$average_samples)
Conta=subset(asv.r, asv.r$contaFactor > 10 )[,51:ncol(asv.r)]
Conta$Family=tax$Family[match(rownames(Conta), rownames(asv))]
message("Number of total ASVs: ", nrow(asv))
message("Number of identified contaminant ASVs removed from the analysis: ", length(rownames(Conta)), "\n", Conta$Family[1],"\n", Conta$Family[2],"\n", Conta$Family[3],"\n", Conta$Family[4],"\n", Conta$Family[5])

#remove any chloroplast or mitochobdria
unwant_tax=tax %>% filter_all(any_vars(str_detect(., 'Mitochondria|Chloroplast')))
colnames(asv.o)
asv.NoMito_noNegs=subset(asv.o, !rownames(asv.o) %in% rownames(unwant_tax))[,!colnames(asv.o)%in% c("Neg.1","Neg.2")]

## remove ASVs with presence considered as counts > 50 reads and in at least in 5 samples
asv_binary=as.data.frame(lapply(asv.NoMito_noNegs, function(x) ifelse(x >= 10, 1, 0))) 
rownames(asv_binary)=rownames(asv.NoMito_noNegs)
asv_noRares=asv_binary %>% filter(rowSums(.) >= 3)

final_asv=subset(asv.NoMito_noNegs, rownames(asv.NoMito_noNegs) %in% rownames(asv_noRares))
dim(asv_noRares)
colnames(final_asv)
las_col_stats=as.data.frame(colSums(final_asv))

# Export normalized and raw ASV tables
final=merge(asv_noF1s, tax[,-1], by="row.names")
colnames(final)
write.table(final, "Outputs/FractionationOptimization_QCfiltered_V3.txt",  quote = FALSE, row.names=F, sep = "\t") 
message("Number of ASVs used in the analysis: ", length(rownames(final)))

# exporting sequences
library(seqinr)
asv_seq=read.table("inputs/FractionationOptimization_ASV_table.txt", header = TRUE, row.names = 1, sep = " ")
names(asv_seq)

anc=read.table("outputs/ANCOM_final", header = T)
vibrios=subset(anc, Family == "Vibrionaceae" & Diff_more_abundant  == "F2")

vibrios_asv_seq = subset(asv_seq, rownames(asv_seq) %in% vibrios$ASV | rownames(asv_seq) %in% c("ASV0011", "ASV0014", "ASV0015")) 
write.fasta(as.list(vibrios_asv_seq$Sequence), rownames(vibrios_asv_seq), "outputs/vibriosASV_outgroup.fasta", open = "w", nbchar = 60, as.string = FALSE)



