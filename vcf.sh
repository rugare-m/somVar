#!/bin/bash

# Job name:
#SBATCH --job-name=fastq
#SBATCH --mem=80G
# Wall clock limit:
#SBATCH --time=24:00:30

# Variables
accession="SRR10556218"
threshold=0.7
#options: high_depth_model.pkl or low_depth_model.pkl
model=high_depth_model.pkl


# Paths
cwd=$(pwd)
fastq="$cwd/fastq"
ref_dir="$cwd/fasta"

# Files
compressed_fasta="$ref_dir/hg38.fa.gz"
dbSNP="$fastq/common_all_20180418.vcf.gz"
COSMIC="$fastq/CosmicCodingMutsV98.vcf.gz"

high_confidence_vcf="${accession}.hc.vcf.gz"

# "$fastq/"

time python "$fastq/indexer.py" -d "$cwd/fastq"

time python "$fastq/annotate_serum.py" -f "$compressed_fasta" -d "$cwd"
"$fastq/serum_filter.sh" -c "$dbSNP" -d "$COSMIC"
time python "$fastq/merge.py" -f "$compressed_fasta" -d "$cwd/"
time python "$fastq/dataframe.py" -f "$compressed_fasta" -v1 COSMIC_gatkbcftools.$accession.vcf.gz -v2 COSMIC_gatkfbayes.$accession.vcf.gz -v3 COSMIC_lofreq.$accession.vcf.gz -v4 COSMIC_mutect2.$accession.vcf.gz -vx $accession.dedup.snps.vcf.gz
time python "$fastq/predictor.py" -o "$high_confidence_vcf" -x "$threshold" -d $accession.dedup.snps.vcf.gz -m "$model"

# clean up
rm -f *.bai COSMOS_* COSMIC_* gatkbcftools*vcf* gatkfbayes*vcf* lofreq.*vcf* mutect2*vcf* bcftools*vcf* *.bqsr.bam* fbayes.*vcf* *.dedup.snps.vcf.gz* *_vcfs.list


