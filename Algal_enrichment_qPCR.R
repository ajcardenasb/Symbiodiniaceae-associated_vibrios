setwd("~/Documents/Bioinformatics_scripts/R_scripts/Fractionation_optimization/qPCR")
library(stringr)
library(reshape2)
library(dplyr)
library(ggplot2)

P1_met=read.table("Fractionation_efficiency_P1_metadata.txt", header = T)
P2_met=read.table("Fractionation_efficiency_P2_metadata.txt", header = T)
P3_met=read.table("Fractionation_efficiency_P3_metadata.txt", header = T)
P4_met=read.table("Fractionation_efficiency_P4_metadata.txt", header = T)
P1_dat=read.table("Fractionation_efficiency_P1_results.txt", header = T, sep = "\t")
P2_dat=read.table("Fractionation_efficiency_P2_results.txt", header = T, sep = "\t")
P3_dat=read.table("Fractionation_efficiency_P3_results.txt", header = T, sep = "\t")
P4_dat=read.table("Fractionation_efficiency_P4_results.txt", header = T, sep = "\t")

#create a single table with sample and cq information
P1_met$Cq=P1_dat$Cq[match(P1_met$Well,P1_dat$Position)]
P2_met$Cq=P2_dat$Cq[match(P2_met$Well,P2_dat$Position)]
P3_met$Cq=P3_dat$Cq[match(P3_met$Well,P3_dat$Position)]
P4_met$Cq=P4_dat$Cq[match(P4_met$Well,P4_dat$Position)]
all_merged=rbind(P1_met,P2_met, P3_met, P4_met) 
all_merged$Sample=gsub("C", "_C", all_merged$Sample)


#remove missing values and blanks, calculate mean between technical replicates
all_merged_NoMissing=subset(all_merged, !Cq == "-    ")
all_merged_NoMissing$Cq=as.numeric(all_merged_NoMissing$Cq)
replicate_Means=all_merged_NoMissing %>% filter(!Sample == "blank") %>%
 group_by(Sample, Primer) %>% dplyr::summarise(Mean_Cq=mean(Cq)) 
            
# add metadata info and prepare table for ratio calculations 
meta_info=as.data.frame(str_split_fixed(replicate_Means$Sample, '_', 5))
colnames(meta_info)=c("SampleType", "Fraction", "Template", "Protocol", "Replicate")
meta_info$Protocol=ifelse(meta_info$Fraction == "F1", "Unfractionated", as.character(meta_info$Protocol))
final=cbind(replicate_Means, meta_info) %>% filter(Primer %in% c("coral", "algae"))
final_m=reshape2::dcast(final, SampleType+Template+Protocol+Replicate~Primer, value.var = "Mean_Cq", mean)

test=subset(final_m, SampleType == "FRE" & Template == "DNA" & Replicate == "C3")
            
## delta delta
final_m$Delta1=final_m$algae-final_m$coral
final_m$FoldChange=2^(-final_m$Delta1) #Calculate ΔCt for each sample
#Now, compute the fold increase in algal abundance relative to unfractionated:
final_m2=reshape2::dcast(final_m, SampleType+Template+Replicate~Protocol, value.var = "FoldChange", mean)
final_m2$FC_P1=final_m2$P1/final_m2$Unfractionated # Convert ΔCt to relative abundance (fold enrichment)
final_m2$FC_P2=final_m2$P2/final_m2$Unfractionated # Convert ΔCt to relative abundance (fold enrichment)
final_m2_long=reshape2::melt(final_m2[,-c(4:6)], id = c("SampleType", "Template", "Replicate"))
#aesthetics 
final_m2_long$variable=gsub("FC_", "", final_m2_long$variable)
final_m2_long$Group=paste(final_m2_long$SampleType, " - ", final_m2_long$Template)
final_m2_long$Group=gsub("FRE", "Fresh", final_m2_long$Group) 
final_m2_long$Group=gsub("FRO", "Frozen", final_m2_long$Group) 

pdf("../outputs/qPCT_algal_enrichment.pdf", width = 6, height = 2.5, pointsize = 12)
ggplot(final_m2_long, aes(x=variable, y=value)) + 
  geom_boxplot() +  labs(y="Fold algal enrichment", x = "Centrifugation protocol") + theme_bw() + facet_wrap(~Group, scale="free", ncol =4) 
dev.off()


