library(phyloseq)
library(ggplot2)
library(reshape2)
library(patchwork)
library(microbiome)
library(scales)
library(plyr)
library(gridExtra)
library(vegan)
library(compositions)
library(GUniFrac)


setwd("~/Documents/Bioinformatics_scripts/R_scripts/Fractionation_optimization/")

asv=read.table("outputs/FractionationOptimization_QCfiltered_V3.txt", header = TRUE, row.names = 1, sep = "\t")[,1:36]
tax=read.table("outputs/FractionationOptimization_QCfiltered_V3.txt", header = TRUE, row.names = 1, sep = "\t")[,37:43]
met=read.table("Inputs/FractionationOptimization_metafileV3.txt", header = TRUE,row.names = 1, sep = "\t")
met$label=rownames(met)
met$Group=paste(met$SampleType, met$NucleicAcid, sep="-")
met$Group=gsub("FRE", "Fresh", met$Group) 
met$Group=gsub("FRO", "Frozen", met$Group) 
met$Group=factor(met$Group, levels = c("Fresh-DNA","Frozen-DNA", "Fresh-RNA", "Frozen-RNA" ))


P_fractions=c("#DEBA6F", "#823038")
P_groups=c( "#ECA1A5", "#D3565D","#60d5b9",  "#2c967e")

##################### Figure 1A: Ordination plot ######################

otu.t= otu_table(asv, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(met))
tax.t= tax_table(as.matrix(tax))

phy.all= phyloseq(otu.t, tax.t,  sam.t)
phy.t=microbiome::transform(phy.all, transform = "clr", target = "OTU", shift = 0, scale = 1)
PCOA_all = ordinate(phy.t, method = "RDA", distance = "euclidean")
pdf("./outputs/Fractionation_ordination.pdf",  width = 8, height =3, pointsize = 12) 
plot_ordination(phy.t,PCOA_all, color = "Group", shape = "Protocol") + geom_point(size = 3, alpha = 1) +
  scale_colour_manual(values=P_groups)  + theme_bw()+  theme( legend.position = "right", legend.title=element_blank()) #legend.position = c(0.80,0.25)
dev.off()


##################### Figure 1B: Dissimilarity ######################

##### dissimilarities ####
asv_clr=apply(asv,2,clr)
dist=vegdist(t(asv_clr), method="euclidean", upper = T)
dist_df=as.matrix(dist)
long_dist = reshape2::melt(dist_df, value.name = c("Distance"))#### Melting Data into Long Format
long_dist$Fraction1=met$Fraction[match(long_dist$Var1, rownames(met))]
long_dist$Fraction2=met$Fraction[match(long_dist$Var2, rownames(met))]
long_dist$Group1=met$Group[match(long_dist$Var1, rownames(met))]
long_dist$Group2=met$Group[match(long_dist$Var2, rownames(met))]
long_dist$Protocol1=met$Protocol[match(long_dist$Var1, rownames(met))]
long_dist$Protocol2=met$Protocol[match(long_dist$Var2, rownames(met))]
long_dist$Comparison=paste(long_dist$Protocol1, long_dist$Protocol2)
long_s=subset(long_dist, Fraction1 == "F1" & Fraction2 == "F2" & Group1==Group2 & !Comparison %in% c("Unfractionated Unfractionated", "P1 P2" ,"P2 P1" ) )
long_s$Group1=factor(long_s$Group1, levels = c("Fresh-DNA","Frozen-DNA", "Fresh-RNA", "Frozen-RNA" ))

dissi_plot=ggplot(long_s, aes(x=Protocol2, y=Distance, fill=Group1)) + 
  geom_boxplot()  +
  labs(y = "Microbiome dissimilarity\nbetween fractions", x= "") + theme_minimal() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_grid(~Group1) +  scale_fill_manual(values=P_groups) 

##### alpha-diversity ####
cnts=t(asv[,!colnames(asv) == "Z.F2.D.P2.C2"])
min(rowSums(cnts)) # determine sample with lowest counts
asv.rar=Rarefy(cnts, min(rowSums(cnts)))$otu.tab.rff
alpha=as.data.frame(t(estimateR(asv.rar)))
alpha$Shannon=diversity(asv.rar, index = "shannon")
alpha$Group=met$Group[match(rownames(alpha), rownames(met))]
alpha$Fraction=met$Fraction[match(rownames(alpha), rownames(met))]
alpha$Protocol=met$Protocol[match(rownames(alpha), rownames(met))]
alpha$Protocol=gsub("Unfractionated", "F1", alpha$Protocol)
alpha$Protocol=factor(alpha$Protocol, levels = c("F1", "P1", "P2"))
alpha$Group=factor(alpha$Group, levels = c("Fresh-DNA","Frozen-DNA", "Fresh-RNA", "Frozen-RNA" ))

alpha_plot=ggplot(alpha, aes(x=Protocol, y=S.chao1, fill= Fraction)) + 
  geom_boxplot()  + theme_minimal() + theme(legend.position = c(0.05,0.75)) + 
  labs(y = "Estimated bacterial\nrichness (Chao1)", x= "")  +
  scale_fill_manual(values=P_fractions) + facet_grid(~Group) 

pdf("Outputs/fractionation_dissimilaritis_alpha.pdf", height = 4, width = 8, pointsize = 12)
dissi_plot+alpha_plot + plot_annotation(tag_levels = 'A')
dev.off()
