#!/bin/bash

#BSUB -J final11

#BSUB -q normal

#BSUB -e %J_error.txt

source activate base
conda activate /project/Fractionation/anaconda_environments/flye

main_path="/project/Fractionation/data/Fractionation_vibrio/final11/"
genome_list="/project/Fractionation/data/Fractionation_vibrio/final11/genome_list"

## assemblies
cat "$genome_list" | \
while read genome; do
flye --nano-corr ${main_path}/raws/${genome}.fastq --out-dir ${main_path}/assemblies/${genome}_flye_output
done

## renaming files
for i in ${main_path}/assemblies/*_flye_output/assembly.fasta; do
  # Extract base folder name
  folder=$(basename $(dirname "$i"))
  
  # Extract the prefix (everything before _flye_output)
  prefix=${folder%%_flye_output}
  
  # Copy and rename
  cp "$i" "${main_path}/assemblies/${prefix}.fasta"
done

## quast
conda activate /project/Fractionation/anaconda_environments/quast
quast  ${main_path}/assemblies/*.fasta -o ${main_path}/quast/

#checkM
conda activate /project/Fractionation/anaconda_environments/checkM
checkm lineage_wf -t 24 -x fasta ${main_path}/assemblies/ ${main_path}/checkm/

## taxonomy
conda activate /project/Fractionation/anaconda_environments/gtdbtk
gtdbtk classify_wf --genome_dir ${main_path}/assemblies/  --out_dir ${main_path}/taxonomy/ --extension fasta

## phylogenomics
conda activate /project/Fractionation/anaconda_environments/gtdbtk
gtdbtk de_novo_wf --genome_dir /project/Fractionation/data/Fractionation_vibrio/final11/phylogenomics/  --extension fasta --out_dir /project/Fractionation/data/Fractionation_vibrio/final11/phylogenomics/test3 --bacteria --outgroup_taxon g__Photobacterium --skip_gtdb_refs  --custom_taxonomy_file /project/Fractionation/data/Fractionation_vibrio/final11/phylogenomics/tax_list

## annotation
conda activate /project/Fractionation/anaconda_environments/prokka

cat "$genome_list" | while read -r genome; do
  prokka \
    --prefix "$genome" \
    --force \
    --outdir "${main_path}/annotations/${genome}_prokka_output" \
    --locustag "$genome" \
    --cpus 64 \
    --rawproduct \
    "${main_path}/assemblies/${genome}.fasta"
done



## retrieval of 16S rRNA sequences

source activate /project/Fractionation/anaconda_environments/barrnap
export EMBOSS_ACDROOT=/project/Fractionation/anaconda_environments/barrnap/share/EMBOSS/acd
export EMBOSS_DATA=/project/Fractionation/anaconda_environments/barrnap/share/EMBOSS/data

cat "$genome_list" | while read -r genome; do
  barrnap -o "${main_path}/16S_sequences/${genome}_rrna.fasta" < "${main_path}/assemblies/${genome}.fasta" > "${main_path}/16S_sequences/${genome}_rrna.gff"

  grep ">16S" -A 1 "${main_path}/16S_sequences/${genome}_rrna.fasta" > "${main_path}/16S_sequences/${genome}_16S.fasta"

  mafft --auto "${main_path}/16S_sequences/${genome}_16S.fasta" > "${main_path}/16S_sequences/${genome}_16S_aligned.fasta"

  cons -sequence "${main_path}/16S_sequences/${genome}_16S_aligned.fasta" \
       -outseq "${main_path}/16S_sequences/${genome}_16S_consensus.fasta" \
       -setcase 1
done


## renaming headers  ### needs checking

cat "$genome_list" | while read -r genome; do
  input_file="${genome}_16S_consensus.fasta"
  output_file="${genome}_16S_consensus_renamed.fasta"

  if [[ -f "$input_file" ]]; then
    sed "1s/^>.*/>${genome}/" "$input_file" > "$output_file"
  else
    echo "Warning: File $input_file not found."
  fi
done

### comparison between 16S rRNA genomic and ASVs

source activate /project/Fractionation/anaconda_environments/blast

cutadapt -g AGGATTAGATACCCTGGTA -a CRRCACGAGCTGACGAC -o Vibrio_Genomes_trimmed_16S.fasta Vibrio_Genomes_Final11_16S.fasta

makeblastdb -in Vibrio_Genomes_trimmed_16S.fasta -dbtype nucl -out isolate_16S_db

blastn -query all_vibriosASV_outgroup.fasta -db isolate_16S_db -out ASV_vs_isolates_besthit.tsv -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_target_seqs 1 -max_hsps 1

