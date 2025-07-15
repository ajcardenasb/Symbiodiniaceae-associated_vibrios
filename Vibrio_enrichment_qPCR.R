setwd("~/Documents/Bioinformatics_scripts/R_scripts/Fractionation_optimization/qPCR/Vibrios/")
library(stringr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(patchwork)

met1=read.table("metadata_plate3_cDNA.txt", header = T, sep = "\t")
met1=met1[!(met1$Sample == ""), ]
met2=read.table("metadata_plate4_cDNA.txt", header = T, sep = "\t")
met2=met2[!(met2$Sample == ""), ]
met2$Sample=gsub("Mst_", "Mste_", met2$Sample)
dat1=read.table("results_plate3_cDNA.txt", header = T, sep = "\t")
dat2=read.table("results_plate4_cDNA.txt", header = T, sep = "\t")
dat1$Position=gsub(" ", "", dat1$Position)
dat2$Position=gsub(" ", "", dat2$Position)
met3=read.table("metadata_Anov_Vib.txt", header = T, sep = "\t")
met3=met3[!(met3$Sample == ""), ]
dat3=read.table("results_Anov_Vib.txt", header = T, sep = "\t")
dat3$Position=gsub(" ", "", dat3$Position)
met4=read.table("../Fractionation_efficiency_P4_metadata.txt", header = T)
met4=met4[!(met4$Sample == ""), ]
dat4=read.table("../Fractionation_efficiency_P4_results.txt", header = T, sep = "\t")
dat4$Position=gsub(" ", "", dat4$Position)


#create a single table with sample and cq information
met1$Cq=dat1$Cq[match(met1$Well,dat1$Position)]
met2$Cq=dat2$Cq[match(met2$Well,dat2$Position)]
met3$Cq=dat3$Cq[match(met3$Well,dat3$Position)]
met4$Cq=dat4$Cq[match(met4$Well,dat4$Position)]

all_merged=rbind(met1, met2, met3, met4) 
all_merged$Sample=gsub("FRO", "Anov", all_merged$Sample)
all_merged=all_merged[!(all_merged$Sample == "blank"), ]

#remove missing values and blanks, calculate mean between technical replicates
all_merged_NoMissing<- all_merged[!(all_merged$Sample == "" | all_merged$Cq == "-"), ]
all_merged_NoMissing$Cq=as.numeric(as.character(all_merged_NoMissing$Cq))
replicate_Means=all_merged_NoMissing %>% 
  group_by(Sample, Primer) %>% 
  dplyr::summarise(Mean_Cq=mean(Cq, na.rm = TRUE)) 

# add metadata info and prepare table for ratio calculations 
meta_info=as.data.frame(str_split_fixed(replicate_Means$Sample, '_', 4))
colnames(meta_info)=c("Species", "Fraction", "Template", "Replicate")
final=cbind(replicate_Means, meta_info)
final_m=reshape2::dcast(final, Species+Fraction+Template+Replicate~Primer, value.var = "Mean_Cq", mean)

## delta delta
final_m$Delta_AlgaeCoral=final_m$algae-final_m$coral
final_m$FoldChange_AlgaeCoral=2^(-final_m$Delta_AlgaeCoral) #Calculate ΔCt for each sample
final_m$Delta_Vibrios=final_m$vibrio-final_m$bacteria
final_m$FoldChange_Vibrios=2^(-final_m$Delta_Vibrios) #Calculate ΔCt for each sample
final_m=final_m[!(final_m$Species == "Mste"), ]
final_m$Species=factor(final_m$Species, levels = c("Anov", "Spis", "Pdam")) #c("Anov", "Spis", "Pdam", "Mste")


#Now, compute the fold increase in algal abundance relative to unfractionated:
final_fractionation=reshape2::dcast(final_m, Species+Template+Replicate~Fraction, value.var = "FoldChange_AlgaeCoral", mean)
final_fractionation$Fold_algal_enrichment=final_fractionation$F2/final_fractionation$F1 # Convert ΔCt to relative abundance (fold enrichment)

final_fractionation_long=reshape2::melt(final_fractionation[,-c(4:5)], id = c("Species", "Template", "Replicate"), value.name = "AlgalEnrichment")
#pdf("../outputs/qPCT_algal_enrichment.pdf", width = 5, height = 3, pointsize = 12)
A_plot=ggplot(final_fractionation_long, aes(x=Species, y=AlgalEnrichment)) + 
  geom_boxplot() +  labs(y="Fold Algal enrichment", x = "") + theme_bw() #+ facet_wrap(~Species+Template, scale="free", ncol =4) 
#dev.off()

### same but with Vibrio enrichment
final_vibrios=reshape2::dcast(final_m, Species+Template+Replicate~Fraction, value.var = "FoldChange_Vibrios", mean)
final_vibrios$Fold_Vibrio_enrichment=final_vibrios$F2/final_vibrios$F1 # Convert ΔCt to relative abundance (fold enrichment)


final_vibrios_long=reshape2::melt(final_vibrios[,-c(4:5)], id = c("Species", "Template", "Replicate"), value.name = "VibrioEnrichment")

#pdf("../outputs/qPCT_algal_enrichment.pdf", width = 5, height = 3, pointsize = 12)
B_plot=ggplot(final_vibrios_long, aes(x=Species, y=VibrioEnrichment)) + 
  geom_boxplot() +  labs(y="Fold Vibrio enrichment", x = "") + theme_bw() #+ facet_wrap(~Species+Template, scale="free", ncol =4) 
#dev.off()

A_plot+B_plot

### correlation?
final_fractionation_long$VibrioEnrichment = final_vibrios_long$VibrioEnrichment
ggplot(final_fractionation_long, aes(x=VibrioEnrichment, y=AlgalEnrichment)) + 
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") #+ facet_wrap(~Species)#+  labs(y="Fold Algal enrichment", x = "") + theme_bw() #+ facet_wrap(~Species+Template, scale="free", ncol =4) 

library(ggpubr)

#assess outliers Cook's Distance (influential point detection)
model <- lm(AlgalEnrichment ~ VibrioEnrichment, data = final_fractionation_long)
plot(cooks.distance(model), type = "h")
abline(h = 4/(nrow(final_fractionation_long)-2), col = "red", lty = 2)
which.max(final_fractionation_long$VibrioEnrichment)
final_fractionation_long=final_fractionation_long[-(which.max(final_fractionation_long$VibrioEnrichment)), ]

# Create the scatter plot with correlation
C_plot=ggscatter(final_fractionation_long, x = "VibrioEnrichment", y = "AlgalEnrichment",
          add = "reg.line",              # Adds regression line
          conf.int = TRUE,               # Adds confidence interval
          cor.coef = TRUE,               # Adds correlation coefficient
          cor.method = "pearson",        # Method: pearson, spearman, kendall
          cor.coeff.args = list(label.x = 2.5, label.y = -15)) +  # Adjust label position
  xlab("Vibrio Enrichment") +
  ylab("Algal Enrichment") +
  theme_minimal()

pdf("../../outputs/qPCT_vibrio_enrichment.pdf", width = 8, height = 3, pointsize = 12)
A_plot+B_plot+C_plot
dev.off()
