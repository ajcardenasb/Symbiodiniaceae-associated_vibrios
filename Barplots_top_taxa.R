library(reshape2)
library(ggplot2)
library(scales)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Fractionation_optimization/")


asv=read.table("outputs/FractionationOptimization_QCfiltered_V3.txt", header = TRUE, row.names = 1, sep = "\t")[,1:36]
tax=read.table("outputs/FractionationOptimization_QCfiltered_V3.txt", header = TRUE, row.names = 1, sep = "\t")[,37:43]
met=read.table("Inputs/FractionationOptimization_metafileV3.txt", header = TRUE,row.names = 1, sep = "\t")
asv$ASV=rownames(asv)
met$Group=paste(met$SampleType, met$NucleicAcid, sep="-")
met$Group=gsub("FRE", "Fresh", met$Group) 
met$Group=gsub("FRO", "Frozen", met$Group) 
names(asv)


asv$ASV=ifelse(asv$ASV %in% asv$ASV[1:20], as.character(asv$ASV), "Others")
asv.gg=aggregate(asv[, 1:36], by = list(asv$ASV), FUN =  sum)
asv.l=melt(asv.gg, id.vars=c("Group.1"), variable.name = "Taxa", value.name = "Abundance")
colnames(asv.l)=c("Taxa","Sample","Abundance")

asv.l$Group=met$Group[match(asv.l$Sample, rownames(met))]
asv.l$Fraction=met$Fraction[match(asv.l$Sample, rownames(met))]
asv.l$Protocol=met$Protocol[match(asv.l$Sample, rownames(met))]
asv.l$Colony=met$Replicate[match(asv.l$Sample, rownames(met))]
asv.l$Taxa2=paste(tax$Genus)[match(asv.l$Taxa, rownames(tax))]
#asv.l$Taxa2=paste(tax$Genus, tax$Family, sep = " | ")[match(asv.l$Taxa, rownames(tax))]
asv.l$Label3=paste(asv.l$Taxa, " (", asv.l$Taxa2, ")", sep = "")

final_asv=asv.l%>% group_by( Colony, Fraction, Group, Protocol, Label3) %>% dplyr::summarise(Abundance=sum(Abundance))
final_asv$Protocol=gsub("Unfractionated", "F1", final_asv$Protocol)
final_asv$Protocol=factor(final_asv$Protocol, levels = c("F1", "P1", "P2"))

P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")

pdf("Outputs/fractionation_barplots_nolegend.pdf", height = 3, width = 6.5, pointsize = 12)
ggplot() + geom_bar(aes(y = Abundance, x = Colony, fill = Label3), 
                             data = final_asv, stat="identity", position = "fill")  +
  labs( y= "Percentage of 16S rRNA sequences", x="") + 
  scale_fill_manual("Taxa",values=P21)  +   theme_classic() + 
  theme( legend.key.size = unit(0.4, "cm"),legend.key.width = unit(0.4,"cm"), legend.text = element_text(size=8),
         legend.position = 'none', axis.text.x=element_text(angle=90,hjust=1),) + 
  guides(fill=guide_legend(ncol=4,  byrow=F)) + facet_grid(~Group+Protocol) + scale_y_continuous( expand = c(0, 0)) 
dev.off()


#legend on the right
ggplot() + geom_bar(aes(y = Abundance, x = Colony, fill = Label3), 
                    data = final_asv, stat="identity", position = "fill")  +
  labs( y= "Percentage of 16S rRNA sequences", x="") + 
  scale_fill_manual("Taxa",values=P21)  +   theme_classic() + 
  theme( legend.key.size = unit(0.4, "cm"),legend.key.width = unit(0.4,"cm"), legend.text = element_text(size=8),
         legend.position = 'right', axis.text.x=element_text(angle=90,hjust=1),) + 
  guides(fill=guide_legend(ncol=1)) + facet_grid(~Group+Protocol) + scale_y_continuous( expand = c(0, 0)) 

#legend at the bottom
pdf("Outputs/fractionation_barplots_legend.pdf", height = 6, width = 10, pointsize = 12)
ggplot() + geom_bar(aes(y = Abundance, x = Colony, fill = Label3), 
                    data = final_asv, stat="identity", position = "fill")  +
  labs( y= "Percentage of 16S rRNA sequences", x="") + 
  scale_fill_manual("Taxa",values=P21)  +   theme_classic() + 
  theme( legend.key.size = unit(0.4, "cm"),legend.key.width = unit(0.4,"cm"), legend.text = element_text(size=8),
         legend.position = 'bottom', axis.text.x=element_text(angle=90,hjust=1),) + 
  guides(fill=guide_legend(ncol=3,  byrow=F)) + facet_grid(~Group+Protocol) + scale_y_continuous( expand = c(0, 0)) 
dev.off()
