# Repository - Symbiodiniaceae-associated_vibrios

This repository contains the scripts used to analyze data and create figures for the manuscript "Vibrio spp. maintain close spatial association with the coral endosymbiotic algae in healthy corals"

Raw sequencing data are deposited in the NCBI Sequence Read Archive (SRA) under BioProject [PRJNA1265228](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1265228)

## Workflow

### Metabaracoding - 16S rRNA data analysis
1. Amplicon Sequence variance (ASV) were inferred using [dada2](https://github.com/benjjneb/dada2) using the script `/Symbiodiniaceae-associated_vibrios/ASV_inference.R`
2. Quality checks (i.e., removal of putatively contaminant ASVs and removal of samples with < 1000 reads) were done using the script `/Symbiodiniaceae-associated_vibrios/QC.R`
3. The script `/Symbiodiniaceae-associated_vibrios/Alpha_Beta_diversity.R` includes alpha diversity estimates done using [Vegan](https://github.com/vegandevs/vegan), as well as ordination plots to represent beta diversity.
4. Statistical comparisons of beta diversity using PERMANOVAs were done using the script `/Symbiodiniaceae-associated_vibrios/Beta_diversity_stats.R`
5. Statistical comparisons of alpha diversity using linear models were done using the script `Alpha_diversity_stats.R`
6. Bar plots of most abundant bacterial taxa were created using the script`/Symbiodiniaceae-associated_vibrios/Barplots_top_taxa.R`
7. Differentially abundance analysis was done using [ANCOMBC](https://github.com/FrederickHuangLin/ANCOMBC) following the script `/Symbiodiniaceae-associated_vibrios/Differential_Abundance_Analysis.R`
8. Results from the differentially abundance analysis were plotted using the script `/Symbiodiniaceae-associated_vibrios/Plotting_ANCOM-BC_results.R`

### qPCR assessment of algal and Vibrio enrichment
1. Algal enrichment was calculated usign the script `/Symbiodiniaceae-associated_vibrios/Algal_enrichment_qPCR.R`

### Vibrio genomics
#### Assembly, ORF prediction, gene quantification and annotation.
1. Genome assembly, assesment, taxonomy placement, phylogenomic analysis, and annotation were done using the script `/Symbiodiniaceae-associated_vibrios/Create vibrio_genomics.sh.`
