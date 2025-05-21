# Repository - Symbiodiniaceae-associated_vibrios

This repository contains the scripts used to analyze data and create figures for the manuscript "Vibrio spp. maintain close spatial association with the coral endosymbiotic algae in healthy corals"

Raw sequencing data are deposited in the NCBI Sequence Read Archive (SRA) under BioProject [PRJNA1265228](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1265228)

## Workflow

### Metabaracoding - 16S rRNA data analysis
1. Amplicon Sequence variance (ASV) were inferred using [dada2](https://github.com/benjjneb/dada2) using the script `/Symbiodiniaceae-associated_vibrios/ASV_inference.R`
2. Quality checks (i.e., removal of putatively contaminant ASVs and removal of samples with < 1000 reads) were done using the script `/Symbiodiniaceae-associated_vibrios/QC.R`
3. Bar plots of most abundant bacterial taxa were created using the script`/Symbiodiniaceae-associated_vibrios/endoliths_16S_barplots.R`
4. Ordination plots and PERMANOVAs were done using [Vegan](https://github.com/vegandevs/vegan) with the script `/Symbiodiniaceae-associated_vibrios/endoliths_16S_ordination_permanovas.R`
5. Alpha diversity estimates and statistical comparisons were done using [Vegan](https://github.com/vegandevs/vegan) following the script `/Symbiodiniaceae-associated_vibrios/endoliths_16S_alphaDiversity.R`
6. Differentially abundance analysis was done using [ANCOMBC](https://github.com/FrederickHuangLin/ANCOMBC) following the script `/Symbiodiniaceae-associated_vibrios/endoliths_16S_ANCOMBC.R`
7. Venn diagrams were plotted to compare metagenomic with metabarcoding results using the script `/Symbiodiniaceae-associated_vibrios/endoliths_16S_vs_metaG_vennDiagram.R`


### qPCR assessment of algal and Vibrio enrichment
1. Algal enrichment was calculated usign the script `/Symbiodiniaceae-associated_vibrios/Algal_enrichment_qPCR.R`

### Vibrio genomics
#### Assembly, ORF prediction, gene quantification and annotation.
1. Assembly, ORF prediction and estimation of gene counts were done using the script `/Symbiodiniaceae-associated_vibrios/endoliths_metaG_assembly.sh`
2. Assess contig length cutoffs to keep for downstream analysis 'ORFs_master_to_taxa_KOs' and `/Symbiodiniaceae-associated_vibrios/Mapping_contigLength_Profiles.R`
3. Create KO and KEGG module count matrix removing contigs < 500 bp `/Symbiodiniaceae-associated_vibrios/ORF500_to_taxa&KOs.R`
