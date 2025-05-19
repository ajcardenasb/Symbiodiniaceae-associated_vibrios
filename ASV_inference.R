setwd("~/Documents/Projects_AU/Symbiodiniaceae_microbiome/Optimization_protocol/Kendalls_data/")
library(dada2)

# Set your working directory

### ### ### ### ### 
### Read files  ### 
### ### ### ### ###  
Sys.setenv('R_MAX_VSIZE'=32000000000)
cat("Reading files")
path <- "."
list.files(path) # prints file names

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

cat("Processing",length(sample.names),"samples:", sample.names)


#Inspect read quality profiles
plotQualityProfile(fnFs[15:17]) # neg controls are bad but all other samples are fine
plotQualityProfile(fnRs[15:17]) # neg controls are bad but all other samples are fine

### ### ### ### ### ### ### ### 
### step 1. Filter and trim ### 
### ### ### ### ### ### ### ### 

cat("Filtering and trimming")
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
### set parameters  maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE) 
head(out) # shows first 10 lines of output


### ### ### ### ### ### ### 
### step 2. Dereplicate ### 
### ### ### ### ### ### ### 
derepF1 <- derepFastq(filtFs, verbose=TRUE)
derepR1 <- derepFastq(filtRs, verbose=TRUE)


### ### ### ### ### ### ### ### ###
### step 3. Learn error rates ### 
### ### ### ### ### ### ### ### ###

cat("Learning error rates")
errF <- learnErrors(derepF1, multithread=TRUE) #1119735570 total bases in 536961 reads from 7 samples will be used for learning the error rates.
errR <- learnErrors(derepR1, multithread=TRUE) #120810308 total bases in 536961 reads from 7 samples will be used for learning the error rates.


### ### ### ### ### ### ### ### 
### step 4. Sample inference ### 
### ### ### ### ### ### ### ### 

#apply the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFs <- dada(derepF1, err=errF, multithread=TRUE)
dadaRs <- dada(derepR1, err=errR, multithread=TRUE)
#Inspecting the returned dada-class object:
#dadaFs[[1]]


### ### ### ### ### ### ### ### ### 
### step 5. Merge paired reads ### 
### ### ### ### ### ### ### ### ### 
mergers <- mergePairs(dadaFs, derepF1, dadaRs, derepR1, verbose=TRUE)


### ### ### ### ### ### ### ### ### ###
### step 6. Construct sequence table ### 
### ### ### ### ### ### ### ### ### ###
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


### ### ### ### ### ### ### ### 
### step 7. Remove chimeras ### 
### ### ### ### ### ###  ### ### 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) # Identified 16005 bimeras out of 47542 input sequences.

## calculate the  frequency of chimeric sequences 
sum(seqtab.nochim)/sum(seqtab) #0.9159348

### ### ### ### ### ### ### ## ### ###
### Track reads through the pipeline ### 
### ### ### ### ### ### ### ## ### ### 

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

### ### ### ### ### ### ### ### 
### step 8. Assign taxonomy ### 
### ### ### ### ### ### ### ### 
taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/Bioinformatics_scripts/SILVA-files/silva_nr_v138_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "~/Documents/Bioinformatics_scripts/SILVA-files/silva_species_assignment_v138.fa.gz")

#taxa <- assignTaxonomy(seqtab.nochim, "/home/cardena/databases/Ribosomal/SILVA_138/dada2/silva_nr_v138_train_set.fa", multithread=TRUE)
#taxa <- addSpecies(taxa, "/home/cardena/databases/Ribosomal/SILVA_138/dada2/silva_species_assignment_v138.fa")

#inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


### ### ### ### ### ###
### Export table  #### 
### ### ### ### ### ### 

samples.out <- rownames(seqtab.nochim)

asv.1=t(seqtab.nochim)
asv.1=cbind(asv.1, "sum"=rowSums(asv.1)) 
asv.2=merge(asv.1, taxa, by="row.names")
colnames(asv.2)[1]="Sequence"
asv.3=asv.2[,c(2:ncol(asv.2),1)]
asv.final=asv.3[order(-asv.3$sum),]
rownames(asv.final) = sprintf("ASV%04d", 1:nrow(asv.final))
write.table(asv.final, "FractionationOptimization_ASV_table.txt",  quote = FALSE)
write.table(track, "FractionationOptimization_ASVs_stats.txt", quote = FALSE)
dim(asv.final)
names(asv.final)

##
