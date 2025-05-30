##### Phylo tree
library(msa)
library(Biostrings)
library(DECIPHER)

# Input Genome sequences
genomic_seqs <- readDNAStringSet("~/Documents/Projects_AU/Fractionation_project/ML_ASV_tree/Vibrio_Genomes_Final11_16S.fasta")
asv_seqs <- readDNAStringSet("~/Documents/Projects_AU/Fractionation_project/ML_ASV_tree/all_vibriosASV_outgroup.fasta")
reference_seqs<-readDNAStringSet("~/Documents/Projects_AU/Fractionation_project/ML_ASV_tree/Vibrio_references.fasta")

# Combine them
all_seqs <- c(asv_seqs, reference_seqs,genomic_seqs) #c(genomic_seqs, asv_seqs, reference_seqs)
all_seqs <- c(reference_seqs,genomic_seqs) #c(genomic_seqs, asv_seqs, reference_seqs)

# Alignment using MUSCLE
#alignment <- msa(all_seqs, method = "Muscle") # lacks memory to complete the task
aligned_set <- msaConvert(alignment, type = "Biostrings::DNAStringSet") # Convert to DNAStringSet or AAStringSet first and then save as fasta
seqinr_alignment <- msaConvert(alignment, type = "seqinr::alignment")# Convert to seqinr format
aligned_seqs <- sapply(seqinr_alignment$seq, paste0, collapse = "") Convert each aligned sequence to a single character string
names(aligned_seqs) <- seqinr_alignment$nam
fasta_lines <- unlist(mapply(function(name, seq) c(paste0(">", name), seq), names(aligned_seqs), aligned_seqs))
writeLines(fasta_lines, con = "outputs/vibrio_aligned.fasta")


# Alignment using DECIPHER
alignment <- AlignSeqs(all_seqs , processors = 4) 
BrowseSeqs(alignment) #visualize the alignment

alignment_matrix <- as.matrix(alignment)
gap_fraction <- colSums(alignment_matrix == "-") / nrow(alignment_matrix)
plot(gap_fraction, type = "h", ylab = "Gap frequency", main = "Gap percentage per column")# Visualize gap profile
gap_fraction <- colSums(alignment_matrix == "-") / nrow(alignment_matrix)# Calculate gap frequency per column
threshold <- 0.5 # Define threshold — e.g., remove columns with >50% gaps
start_col <- which(gap_fraction < threshold)[1]# Find left trim position
end_col <- tail(which(gap_fraction < threshold), 1)# Find right trim position
trimmed_matrix <- alignment_matrix[, start_col:end_col]# Trim only ends
trimmed_seqs <- DNAStringSet(apply(trimmed_matrix, 1, paste0, collapse = ""))# Rebuild DNAStringSet
names(trimmed_seqs) <- rownames(trimmed_matrix)
BrowseSeqs(trimmed_seqs) #visualize the trimmed alignment
writeXStringSet(trimmed_seqs, "~/Documents/Projects_AU/Fractionation_project/ML_ASV_tree/Vibrios_alignment_trimmed.fasta") #export the trimmed alignment

# ML tree
library(phangorn)
library(ape)

alignment_matrix <- as.matrix(trimmed_seqs)# Assume you have: trimmed_seqs (DNAStringSet)
alignment_phydat <- phyDat(alignment_matrix, type = "DNA")# Convert to phyDat object
dm <- dist.ml(alignment_phydat)# Calculate distance matrix
nj_tree <- NJ(dm)# Build a neighbor-joining tree as a starting point
plot(nj_tree, main = "Starting NJ Tree")# Plot it (optional)
fit <- pml(nj_tree, data = alignment_phydat)# Create a pml object (phylogenetic ML model)
fit_ml <- optim.pml(fit, model = "GTR", optInv = TRUE, optGamma = TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))# Optimize the tree with ML
plot(fit_ml$tree, main = "Maximum Likelihood Tree") # Plot the Final ML Tree
tree_rooted <- root(fit_ml$tree, outgroup =  c("ASV0011", "ASV0014", "ASV0015"), resolve.root = TRUE)# Optional: root or ladderize
plot(tree_rooted)
# bs <- bootstrap.pml(fit, bs = 100, optNni = TRUE)
# plotBS(fit_ml$tree, bs, p = 50, main = "ML Tree with Bootstrap")
write.tree(tree_rooted, file = "~/Documents/Projects_AU/Fractionation_project/ML_ASV_tree/Vibrio_ASV+genome_ML_tree_rooted.newick")# Save to Newick file
