library(phyloseq)
library(ggplot2)
library(reshape2)
library(patchwork)
library(microbiome)
library(ANCOMBC)
library(scales)
library(plyr)
library(gridExtra)


setwd("~/Documents/Bioinformatics_scripts/R_scripts/Fractionation_optimization/")

asv=read.table("outputs/FractionationOptimization_QCfiltered_V3.txt", header = TRUE, row.names = 1, sep = "\t")[,1:36]
tax=read.table("outputs/FractionationOptimization_QCfiltered_V3.txt", header = TRUE, row.names = 1, sep = "\t")[,37:43]
map=read.table("Inputs/FractionationOptimization_metafileV3.txt", header = TRUE,row.names = 1, sep = "\t")
map$Template <- paste(map$NucleicAcid,map$SampleType, sep="-")
map$Method <- paste(map$Protocol,map$Fraction,sep="-")

tax$Phylum=ifelse(is.na(tax$Phylum), paste("Unclassified ",tax$Kingdom, sep = ""), as.character(tax$Phylum) )
tax$Class=ifelse(is.na(tax$Class), as.character(tax$Phylum), as.character(tax$Class) )
tax$Order=ifelse(is.na(tax$Order), as.character(tax$Class), as.character(tax$Order) )
tax$Family=ifelse(is.na(tax$Family), as.character(tax$Order), as.character(tax$Family) )
tax$Genus=ifelse(is.na(tax$Genus), as.character(tax$Family), as.character(tax$Genus) )

tax$Species=rownames(tax)
map$Template <- paste(map$NucleicAcid,map$SampleType, sep="-")

otu.t= otu_table(asv, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t= tax_table(as.matrix(tax))

phy.all= phyloseq(otu.t, tax.t,  sam.t)

##### ASV level across each treatment group #####

### DNA-fresh P1
dna_fre_p1=subset_samples(phy.all, Template == "DNA-Fresh" & !Method == "P2-F2")
res1=ancombc2(data=dna_fre_p1,fix_formula="Fraction" ,p_adj_method = "fdr", group = "Method", tax_level = "Species" , global = NULL)
res1_df=res1$res
res1_sig=(res1_df[res1_df[,13]=="TRUE",])[,c(1,3,5,7,9,11)] # If this command produces a table with 0 rows, it means there are no significant genera
colnames(res1_sig)=c("ASV","LFC","se_LFC","W",	"pval",	"qval")
res1_sig$Diff_more_abundant=ifelse( res1_sig$W < 0 , "F1", "F2")
res1_sig$Genus=tax$Genus[match(res1_sig$ASV, rownames(tax))]
res1_sig$Family=tax$Family[match(res1_sig$ASV, rownames(tax))]
res1_sig$Comparison = "dna_fre_p1"

### DNA-fresh P2
dna_fre_p2=subset_samples(phy.all, Template == "DNA-Fresh" & !Method == "P1-F2")
res2=ancombc2(data=dna_fre_p2,fix_formula="Fraction" ,p_adj_method = "fdr", group = "Method", tax_level = "Species" , global = NULL)
res2_df=res2$res
res2_sig=(res2_df[res2_df[,13]=="TRUE",])[,c(1,3,5,7,9,11)] # If this command produces a table with 0 rows, it means there are no significant genera
colnames(res2_sig)=c("ASV","LFC","se_LFC","W",	"pval",	"qval")
res2_sig$Diff_more_abundant=ifelse( res2_sig$W < 0 , "F1", "F2")
res2_sig$Genus=tax$Genus[match(res2_sig$ASV, rownames(tax))]
res2_sig$Family=tax$Family[match(res2_sig$ASV, rownames(tax))]
res2_sig$Comparison = "dna_fre_p2"
  
### DNA-Frozen P1
dna_fro_p1=subset_samples(phy.all, Template == "DNA-Frozen" & !Method == "P2-F2")
res3=ancombc2(data=dna_fro_p1,fix_formula="Fraction" ,p_adj_method = "fdr", group = "Method", tax_level = "Species" , global = NULL)
res3_df=res3$res
res3_sig=(res3_df[res3_df[,13]=="TRUE",])[,c(1,3,5,7,9,11)] # If this command produces a table with 0 rows, it means there are no significant genera
colnames(res3_sig)=c("ASV","LFC","se_LFC","W",	"pval",	"qval")
res3_sig$Diff_more_abundant=ifelse( res3_sig$W < 0 , "F1", "F2")
res3_sig$Genus=tax$Genus[match(res3_sig$ASV, rownames(tax))]
res3_sig$Family=tax$Family[match(res3_sig$ASV, rownames(tax))]
res3_sig$Comparison = "dna_fro_p1"

### DNA-Frozen P2
dna_fro_p2=subset_samples(phy.all, Template == "DNA-Frozen" & !Method == "P1-F2")
res4=ancombc2(data=dna_fro_p2,fix_formula="Fraction" ,p_adj_method = "fdr", group = "Method", tax_level = "Species" , global = NULL)
res4_df=res4$res
res4_sig=(res4_df[res4_df[,13]=="TRUE",])[,c(1,3,5,7,9,11)] # If this command produces a table with 0 rows, it means there are no significant genera
colnames(res4_sig)=c("ASV","LFC","se_LFC","W",	"pval",	"qval")
res4_sig$Diff_more_abundant=ifelse( res4_sig$W < 0 , "F1", "F2")
res4_sig$Genus=tax$Genus[match(res4_sig$ASV, rownames(tax))]
res4_sig$Family=tax$Family[match(res4_sig$ASV, rownames(tax))]
res4_sig$Comparison = "dna_fro_p2"


### RNA-fresh P1
rna_fre_p1=subset_samples(phy.all, Template == "RNA-Fresh" & !Method == "P2-F2")
res5=ancombc2(data=rna_fre_p1,fix_formula="Fraction" ,p_adj_method = "fdr", group = "Method", tax_level = "Species" , global = NULL)
res5_df=res5$res
res5_sig=(res5_df[res5_df[,13]=="TRUE",])[,c(1,3,5,7,9,11)] # If this command produces a table with 0 rows, it means there are no significant genera
colnames(res5_sig)=c("ASV","LFC","se_LFC","W",	"pval",	"qval")
res5_sig$Diff_more_abundant=ifelse( res5_sig$W < 0 , "F1", "F2")
res5_sig$Genus=tax$Genus[match(res5_sig$ASV, rownames(tax))]
res5_sig$Family=tax$Family[match(res5_sig$ASV, rownames(tax))]
res5_sig$Comparison = "rna_fre_p1"

### RNA-fresh P2
rna_fre_p2=subset_samples(phy.all, Template == "RNA-Fresh" & !Method == "P1-F2")
res6=ancombc2(data=rna_fre_p2,fix_formula="Fraction" ,p_adj_method = "fdr", group = "Method", tax_level = "Species" , global = NULL)
res6_df=res6$res
res6_sig=(res6_df[res6_df[,13]=="TRUE",])[,c(1,3,5,7,9,11)] # If this command produces a table with 0 rows, it means there are no significant genera
colnames(res6_sig)=c("ASV","LFC","se_LFC","W",	"pval",	"qval")
res6_sig$Diff_more_abundant=ifelse( res6_sig$W < 0 , "F1", "F2")
res6_sig$Genus=tax$Genus[match(res6_sig$ASV, rownames(tax))]
res6_sig$Family=tax$Family[match(res6_sig$ASV, rownames(tax))]
res6_sig$Comparison = "rna_fre_p2"

### RNA-Frozen P1
rna_fro_p1=subset_samples(phy.all, Template == "RNA-Frozen" & !Method == "P2-F2")
res7=ancombc2(data=rna_fro_p1,fix_formula="Fraction" ,p_adj_method = "fdr", group = "Method", tax_level = "Species" , global = NULL)
res7_df=res7$res
res7_sig=(res7_df[res7_df[,13]=="TRUE",])[,c(1,3,5,7,9,11)] # If this command produces a table with 0 rows, it means there are no significant genera
colnames(res7_sig)=c("ASV","LFC","se_LFC","W",	"pval",	"qval")
res7_sig$Diff_more_abundant=ifelse( res7_sig$W < 0 , "F1", "F2")
res7_sig$Genus=tax$Genus[match(res7_sig$ASV, rownames(tax))]
res7_sig$Family=tax$Family[match(res7_sig$ASV, rownames(tax))]
res7_sig$Comparison = "rna_fro_p1"

### RNA-Frozen P2
rna_fro_p2=subset_samples(phy.all, Template == "RNA-Frozen" & !Method == "P1-F2")
res8=ancombc2(data=rna_fro_p2,fix_formula="Fraction" ,p_adj_method = "fdr", group = "Method", tax_level = "Species" , global = NULL)
res8_df=res8$res
res8_sig=(res8_df[res8_df[,13]=="TRUE",])[,c(1,3,5,7,9,11)] # If this command produces a table with 0 rows, it means there are no significant genera
colnames(res8_sig)=c("ASV","LFC","se_LFC","W",	"pval",	"qval")
res8_sig$Diff_more_abundant=ifelse( res8_sig$W < 0 , "F1", "F2")
res8_sig$Genus=tax$Genus[match(res8_sig$ASV, rownames(tax))]
res8_sig$Family=tax$Family[match(res8_sig$ASV, rownames(tax))]
res8_sig$Comparison = "rna_fro_p2"

### after completing all comparisons you can merge them into one file
all=rbind(res1_sig, res2_sig,res3_sig, res4_sig,res5_sig, res6_sig,res7_sig, res8_sig ) 
#write.table(all, "outputs/ANCOM_final",quote = F, row.names = F)


##### ASV level overall #####
res9=ancombc2(data=phy.all,fix_formula="Fraction",p_adj_method = "fdr", tax_level = "Species")
res9_df=res9$res
res9_sig=(res9_df[res9_df[,13]=="TRUE",])[,c(1,3,5,7,9,11)]
colnames(res9_sig)=c("ASV","LFC","se_LFC","W",	"pval",	"qval")
res9_sig$Diff_more_abundant=ifelse( res9_sig$W < 0 , "F1", "F2")
res9_sig$Genus=tax$Genus[match(res9_sig$ASV, rownames(tax))]
res9_sig$Family=tax$Family[match(res9_sig$ASV, rownames(tax))]
write.table(res9_sig, "outputs/ANCOM_final_overall_ASV",quote = F, row.names = F)

##### genus level overall #####
res10=ancombc2(data=phy.all,fix_formula="Fraction",p_adj_method = "fdr", tax_level = "Genus")
res10_df=res10$res
res10_sig=(res10_df[res10_df[,13]=="TRUE",])[,c(1,3,5,7,9,11)]
colnames(res10_sig)=c("ASV","LFC","se_LFC","W",	"pval",	"qval")
res10_sig$Diff_more_abundant=ifelse( res10_sig$W < 0 , "F1", "F2")
res10_sig$Genus=tax$Genus[match(res10_sig$ASV, rownames(tax))]
res10_sig$Family=tax$Family[match(res10_sig$ASV, rownames(tax))]
write.table(res10_sig, "outputs/ANCOM_final_overall_genus",quote = F, row.names = F)
