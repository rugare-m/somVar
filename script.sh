#!/bin/bash

# Job name:
#SBATCH --job-name=ml_pipeline
#SBATCH --mem=80G
# Wall clock limit:
#SBATCH --time=24:00:30

# Paths
ref_dir="/users/rugarem/volatile/chapter3/data/reference"
cwd="/users/rugarem/volatile/chapter3/nextflow"

# Files
compressed_fasta="$ref_dir/hg38.fa.gz"
uncompressed_fasta="$ref_dir/uhg38.fa"
dbSNP="common_all_20180418.vcf.gz"
COSMIC="CosmicCodingMutsV98.vcf.gz"

# Variables
accession="SRR10556218"
high_confidence_vcf=$accession.high_conf.vcf.gz

# Recommended threshold for high depth model is 0.7 and 0.75 for low depth model
threshold=0.7
model=high_depth_model.pkl


time python bwa.py -f "$ref_dir/hg38" -d "$cwd" -t 60
time python indexer.py -d "$cwd"
time python gatk.py -f "$compressed_fasta" -k "1000G_omni2.5.hg38.vcf.gz"
time python bcftools.py -f "$compressed_fasta" -d "$cwd"
time python freebayes.py -f "$uncompressed_fasta" -d "$cwd"
time python lofreq.py -f "$compressed_fasta" -d "$cwd"
time python mutect2.py -R "$compressed_fasta" -d "$cwd"
time python annotate_serum.py -f "$compressed_fasta" -d "$cwd"
./serum_filter.sh -c "$dbSNP" -d "$COSMIC"
time python merge.py -f "$compressed_fasta" -d "$cwd/"
time python dataframe.py -f "$compressed_fasta" -b COSMIC_gatkbcftools.$accession.vcf.gz -r COSMIC_gatkfbayes.$accession.vcf.gz -l COSMIC_lofreq.$accession.vcf.gz -m COSMIC_mutect2.$accession.vcf.gz -x $accession.dedup.snps.vcf.gz
time python predictor.py -o "$high_confidence_vcf" -x "$threshold" -d $accession.dedup.snps.vcf.gz -m "$model"

# clean up
rm -f *.bai COSMOS_* COSMIC_* gatkbcftools*vcf* gatkfbayes*vcf* lofreq.*vcf* mutect2*vcf* bcftools*vcf* *.bqsr.bam* fbayes.*vcf* *.dedup.snps.vcf.gz* *_vcfs.list


