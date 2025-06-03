setwd("~/Documents/Bioinformatics_scripts/R_scripts/Fractionation_optimization/")
library(ggplot2)
library(dplyr)

#anc=read.table("outputs/ANCOM_results.txt", header = T)
anc=read.table("outputs/ANCOM_final", header = T, sep = "\t")
tax=read.table("inputs/FractionationOptimization_ASV_table.txt", header = TRUE, row.names = 1, sep = " ")[,51:58]
anc$Genus=tax$Genus[match(anc$ASV,rownames(tax) )]
anc$Family=tax$Family[match(anc$ASV,rownames(tax) )]


#### barplots with DE numbers ####
sum1= anc %>% group_by( Diff_more_abundant, Comparison) %>% tally() #%>% filter(Diff_more_abundant=="F2")

sum1_totals <- sum1 %>%
  group_by(Comparison) %>%
  dplyr::summarise(total = sum(n)) %>%
  arrange(total)

# Reorder Comparison levels based on total
sum1$Comparison <- factor(sum1$Comparison, levels = sum1_totals$Comparison)

sum1$n =ifelse(sum1$Diff_more_abundant == "F1", sum1$n*-1, sum1$n)


# Plot
pdf("Outputs/DAA_Number_barplots.pdf", height = 4, width = 8, pointsize = 12)
ggplot(sum1, aes(x = Comparison, y = n, fill = Diff_more_abundant)) +
  geom_bar(stat = "identity", position = "identity", width = 0.6) +
  scale_y_continuous(name = "Number of Differentially Abundant ASVs") +
  scale_fill_manual(values = c("F1" = "#DEBA6F", "F2" = "#823038")) +
  coord_flip() +
  labs(x = "Comparison", fill = "More abundant In") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )
dev.off()

#### Figure Common and unique patterns  ####
# tax$Phylum=ifelse(is.na(tax$Phylum), paste("Unclassified ",tax$Kingdom, sep = ""), as.character(tax$Phylum) )
# tax$Class=ifelse(is.na(tax$Class), as.character(tax$Phylum), as.character(tax$Class) )
# tax$Order=ifelse(is.na(tax$Order), as.character(tax$Class), as.character(tax$Order) )
# tax$Family=ifelse(is.na(tax$Family), as.character(tax$Order), as.character(tax$Family) )
# tax$Genus=ifelse(is.na(tax$Genus), as.character(tax$Family), as.character(tax$Genus) )

# Summarize and calculate percentages

anc_full_tax=merge(anc, tax, by.x = "ASV", by.y =0)
anc_full_tax %>% group_by( ASV, Comparison) %>% tally() 

summary_ASVs= anc_full_tax %>%
  group_by(ASV) %>%  filter(Diff_more_abundant=="F2")%>%
 dplyr::summarise(n_comparisons = n_distinct(Comparison), percent_comparisons = (n_comparisons / 8) * 100,
                  .groups = "drop") %>%
  arrange(desc(percent_comparisons)) %>%
  mutate(ASV = factor(ASV, levels = ASV[order(percent_comparisons)]))  # preserve sort order

summary_ASVs$Genus=tax$Genus[match(summary_ASVs$ASV,rownames(tax),)]

top_ASVs = summary_ASVs[c(1:12, 15:17, 19:23),]
top_ASVs$ASV=paste(top_ASVs$ASV," (",top_ASVs$Genus, ")", sep = "")
top_ASVs$ASV <- factor(top_ASVs$ASV, levels = top_ASVs$ASV[order(top_ASVs$percent_comparisons)])


# Bar plot showing top ASVs across comparisons

P20=c("#1f77b4","#ff7f0f","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#16bece", "#669900","#99cc33","#ccee66","#006699","#3399cc","#990066","#cc3399","#ff6600","#ff9900","#BABCCA", "#D6D6D6")

pdf("Outputs/ASV_percentage_comparisons.pdf", height = 3, width = 6, pointsize = 12)
ggplot(top_ASVs, aes(x = ASV, y = percent_comparisons)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "",
    x = "",
    y = "Percent of Comparisons",
    fill = "Genus"
  ) +   theme_minimal() #+ scale_fill_manual(values=P20) +
dev.off()


# Count total number of appearances per genus/family
genus_counts <- summary_ASVs %>%
  group_by(Genus) %>%
  dplyr::summarise(total_ASV_hits = sum(n_comparisons)) %>%
  arrange(desc(total_ASV_hits))

summary_genus= anc_full_tax %>%
  group_by(Genus.y) %>%  filter(Diff_more_abundant=="F2") %>%
  dplyr::summarise(n_comparisons = n_distinct(Comparison), percent_comparisons = (n_comparisons / 8) * 100,
                   .groups = "drop") %>%
  arrange(desc(percent_comparisons)) %>%
  mutate(Genus = factor(Genus.y, levels = Genus.y[order(percent_comparisons)]))  # preserve sort order

summary_family= anc_full_tax %>%
  group_by(Family.y) %>%  filter(Diff_more_abundant=="F2") %>%
  dplyr::summarise(n_comparisons = n_distinct(Comparison), percent_comparisons = (n_comparisons / 8) * 100,
                   .groups = "drop") %>%
  arrange(desc(percent_comparisons)) %>%
  mutate(Family = factor(Family.y, levels = Family.y[order(percent_comparisons)]))  # preserve sort order


# plot 

top_families = summary_family[c(1:20),]

pdf("Outputs/Family_percentage_comparisons_F2.pdf", height = 3, width = 6, pointsize = 12)
ggplot(top_families, aes(x = Family, y = percent_comparisons)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "",
    x = "",
    y = "Percent of Comparisons",
    fill = "Genus"
  ) +   theme_minimal() #+ scale_fill_manual(values=P20) +
dev.off()

# Spider / Radar plot at the genus level F2
library(dplyr)
library(tidyr)
library(fmsb)
library(tibble)

anc$Comparison2=gsub("_p1|_p2", "", anc$Comparison)

# Step 1: Count ASVs per Genus and Comparison
asv_counts <- anc %>% filter(Diff_more_abundant=="F2") %>%
  group_by(Family, Comparison2) %>%
  dplyr::summarise(n = n(), .groups = "drop")

# Step 2: Identify Top 10 Genera with most ASVs across all comparisons
top_family <- asv_counts %>%
  group_by(Family) %>%
  dplyr::summarise(total = sum(n)) %>%
  arrange(desc(total)) %>%
  slice(c(1:2, 4:10, 13)) %>% # slice_head(n = 10)
  pull(Family)

test=asv_counts %>%
  group_by(Family) %>%
  dplyr::summarise(total = sum(n)) %>%
  arrange(desc(total)) 

# Step 3: Filter to keep only top genera and reshape data
radar_data <- asv_counts %>%
  filter(Family %in% top_family) %>%
  pivot_wider(names_from = Family, values_from = n, values_fill = 0)

radar_data$Comparison2=gsub("fre", "Fresh",radar_data$Comparison2) 
radar_data$Comparison2=gsub("fro", "Frozen",radar_data$Comparison2) 
radar_data$Comparison2=gsub("dna", "DNA",radar_data$Comparison2) 
radar_data$Comparison2=gsub("rna", "RNA",radar_data$Comparison2) 

# Step 4: fmsb requires max and min rows
radar_ready <- radar_data %>%
 column_to_rownames("Comparison2") %>%
  as.data.frame()

radar_ready <- radar_ready[ , top_family] # Ensure only top genera columns
radar_ready <- rbind(apply(radar_ready, 2, max),  # max
                     rep(0, ncol(radar_ready)),    # min
                     radar_ready)



# step 5 plot
colors_border <- c("#D3565D","#2c967e", "#ECA1A5",  "#60d5b9") # Customize as needed

pdf("Outputs/spider_plot_F2.pdf", height = 5, width = 6, pointsize = 12)
radar_plot=radarchart(radar_ready,
           axistype = 1,
           # Custom polygon
           pcol = colors_border,
           plwd = 2,
           plty = 1,
           # Custom the grid
           cglcol = "grey", cglty = 1, axislabcol = "grey", caxislabels = seq(0, max(radar_ready), 2), cglwd = 0.8,
           # Custom labels
           vlcex = 0.8
)

legend(x = "topright", legend = rownames(radar_ready)[3:nrow(radar_ready)],
       bty = "n", pch = 20, col = colors_border, text.col = "black", cex = 0.9, pt.cex = 1.5)
dev.off()

######## lollipop plot for the overall comparison

# Summarize number of ASVs and median logFC per genus
anc_overall=read.table("outputs/ANCOM_final_overall_ASV", header = T, sep = "\t", fill = T)
anc_overall$LFC <- as.numeric(anc_overall$LFC)

summary_table_lolli <- anc_overall %>%
  group_by(Genus, Family, Diff_more_abundant) %>%
  dplyr::summarise(
    n_ASVs = n(),
    median_logFC = median(LFC, na.rm = TRUE)) %>%
  arrange(desc(n_ASVs))

summary_table_lolli$Taxa=ifelse(is.na(summary_table_lolli$Genus), 
                                 paste("Unclassified",summary_table_lolli$Family, sep= " " ),
                                 paste(summary_table_lolli$Genus, " (",summary_table_lolli$Family,")", sep = "" ) )

#top_lolli=summary_table_lolli[c(1:2, 5:19, 21:23, 26:32, 34,37:38 ),]
top_lolli=summary_table_lolli[c(2:31),]

# Plot
pdf("Outputs/lollipop_plot_overall_F2.pdf", height = 8.6, width = 7, pointsize = 12)
ggplot(top_lolli, aes(x = reorder(Taxa, median_logFC), y = median_logFC)) +
  geom_segment(aes(xend = Taxa, y = 0, yend = median_logFC), color = "grey") +
  geom_point(aes(size = n_ASVs), color = "black") +
  geom_text(aes(label = n_ASVs), vjust = 0.5, color = "white", size = 3.5) +
  scale_size(range = c(5, 12)) +  # Adjust bubble sizes
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none", axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14))+
  labs(
    x = "",
    y = "Median log fold change"  )
dev.off()


#### Vibrios only ####

#DA_vibrios_treatment=anc %>% filter(Diff_more_abundant == "F2" & Family == "Vibrionaceae") 
DA_vibrios_treatment=anc %>% filter(Diff_more_abundant == "F2" & Genus == "Vibrio") 
#DA_vibrios_overall=anc_overall %>% filter(Diff_more_abundant == "F2" & Genus == "Vibrio") 
#intersect(DA_vibrios_treatment$ASV, DA_vibrios_overall$ASV)

# exporting sequences
library(seqinr)
asv_seq=read.table("inputs/FractionationOptimization_ASV_table.txt", header = TRUE, row.names = 1, sep = " ")
names(asv_seq)

# anc=read.table("outputs/ANCOM_final", header = T, sep = "\t")
# vibrios=subset(anc, Genus == "Vibrio" & Diff_more_abundant  == "F2")
# 
# vibrios_asv_seq = subset(asv_seq, rownames(asv_seq) %in% vibrios$ASV | rownames(asv_seq) %in% c("ASV0011", "ASV0014", "ASV0015")) 
# write.fasta(as.list(vibrios_asv_seq$Sequence), rownames(vibrios_asv_seq), "outputs/vibriosASV_outgroup.fasta", open = "w", nbchar = 60, as.string = FALSE)

vibrios_asv=subset(tax, Genus == "Vibrio")
vibrios_asv_seq = subset(asv_seq, rownames(asv_seq) %in% rownames(vibrios_asv) | rownames(asv_seq) %in% c("ASV0011", "ASV0014", "ASV0015")) 
write.fasta(as.list(vibrios_asv_seq$Sequence), rownames(vibrios_asv_seq), "outputs/all_vibriosASV_outgroup.fasta", open = "w", nbchar = 60, as.string = FALSE)

