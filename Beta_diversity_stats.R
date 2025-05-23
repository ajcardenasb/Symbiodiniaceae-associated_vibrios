
library(vegan)
library(compositions)
library(pairwiseAdonis)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Fractionation_optimization/")

asv=read.table("outputs/FractionationOptimization_QCfiltered_V3.txt", header = TRUE, row.names = 1, sep = "\t")[,1:36]
tax=read.table("outputs/FractionationOptimization_QCfiltered_V3.txt", header = TRUE, row.names = 1, sep = "\t")[,37:43]
met=read.table("Inputs/FractionationOptimization_metafileV3.txt", header = TRUE,row.names = 1, sep = "\t")
met$label=rownames(met)
met$Group=paste(met$SampleType, met$NucleicAcid, sep="-")
met$Group=gsub("FRE", "Fresh", met$Group) 
met$Group=gsub("FRO", "Frozen", met$Group) 
names(asv)



#dist=vegdist(t(asv_clr), method="euclidean", upper = T)
#dist_df=as.matrix(dist)

######## Stats on community composition ############
asv_clr=apply(asv,2,clr)
asv.n=as.data.frame(t(asv_clr))

#asv.n=as.data.frame(t(sweep(asv,2,colSums(asv),`/`)))

asv.n$SampleType=met$SampleType[match(rownames(asv.n), rownames(met))]
asv.n$Fraction=met$Fraction[match(rownames(asv.n), rownames(met))]
asv.n$Protocol=met$Protocol[match(rownames(asv.n), rownames(met))]
asv.n$NucleicAcid=met$NucleicAcid[match(rownames(asv.n), rownames(met))]
asv.n$Group=paste(asv.n$SampleType,"-",asv.n$NucleicAcid)
asv.n$Group2=paste(asv.n$Fraction,"-",asv.n$Protocol)

## Model across unfractionated samples
asv.f1=subset(asv.n, Fraction == "F1")
adonis2(asv.f1[,1:1150] ~ asv.f1$NucleicAcid * asv.f1$SampleType, method = "euclidian", by = "terms")
adonis2(asv.f1[,1:1150] ~ asv.f1$SampleType, method = "euclidean" )
adonis2(asv.f1[,1:1150] ~ asv.f1$Group , method = "euclidean" )

### Model to test for fraction across groups of samples
adonis2(asv.n[,1:1150] ~ asv.n$Fraction + asv.n$NucleicAcid + asv.n$SampleType + asv.n$Protocol, method = "euclidian", by = "terms")

#Sub-setting per group to test for specific differences F1 vs F2
rna.fro=subset(asv.n, NucleicAcid == "RNA" & SampleType == "Frozen")
#pairwise.adonis(rna.fro[,1:1150], rna.fro$Fraction,  sim.method = "euclidian", p.adjust.m = "fdr", perm = 999) 
adonis2(rna.fro[,1:1150] ~ rna.fro$Fraction, method = "euclidian")

#Sub-setting per group to test for specific differences P1 vs P2
rna.fro=subset(asv.n, NucleicAcid == "RNA" & SampleType == "Fresh" & !Protocol == "Unfractionated")
adonis2(rna.fro[,1:1150] ~ rna.fro$Protocol, method = "euclidian")
