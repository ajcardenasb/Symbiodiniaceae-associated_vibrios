setwd("~/Documents/Bioinformatics_scripts/R_scripts/Fractionation_optimization/")
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(tibble)
library(ggplot2)
library(dplyr)

viFa=read.table("inputs/virulence_ggtree", header = T, sep = "\t", row.names = 1)
viFa[is.na(viFa)] <- 0
tree <- read.tree("outputs/vibrio_phylogenomics_no_references.nwk") 

# Genome names in gene_matrix need to be in same order as tree tip labels
viFa <- viFa[tree$tip.label, ]
viFa[] <- lapply(viFa, function(x) as.factor(as.character(x)))

p <- ggtree(tree, branch.length = "none") +  # branch.length = "none" removes variability in branch lengths — all tips end at the same vertical point.
  geom_tiplab(size = 3, align = TRUE, linetype = "dotted") #align = TRUE in geom_tiplab() makes labels extend to the same x-axis. Optional: linetype = "dotted" adds guide lines between tip and label.


pdf(file = "./outputs/Virulence_ggtree.pdf",  width = 7, height = 5, pointsize = 12) 
gheatmap(p, viFa,
         offset = 0.05,
         width = 2,
         font.size = 3,
         colnames_angle = 90,
         hjust = 0,
         colnames_position = "top") +
  scale_fill_manual(values = c("0" = "#F7F7F7",  # light gray for absence
                               "1" = "#0E6C8C")) +  # blue for presence
  theme(legend.position = "none") +
  coord_cartesian(clip = 'off')
dev.off()

