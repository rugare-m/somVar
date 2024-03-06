#!/bin/bash

# Job name:
#SBATCH --job-name=axolotl
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
axolotl="$cwd/axolotl"
ref_dir="$cwd/FASTA"

# Files
compressed_fasta="$ref_dir/hg38.fa.gz"
uncompressed_fasta="$ref_dir/uhg38.fa"
dbSNP="$axolotl/common_all_20180418.vcf.gz"
COSMIC="$axolotl/CosmicCodingMutsV98.vcf.gz"
high_confidence_vcf="${accession}_hc.vcf.gz"

# "$axolotl/"
time python "$axolotl/bwa.py" -f "$ref_dir/hg38" -d "$cwd" -t 60
time python "$axolotl/indexer.py" -d "$cwd"
time python "$axolotl/gatk.py" -f "$compressed_fasta" -k "$axolotl/1000G_omni2.5.hg38.vcf.gz"

# Run commands in parallel
(time python "$axolotl/bcftools.py" -f "$compressed_fasta" -d "$cwd") &
(time python "$axolotl/freebayes.py" -f "$uncompressed_fasta" -d "$cwd") &
(time python "$axolotl/lofreq.py" -f "$compressed_fasta" -d "$cwd") &
(time python "$axolotl/mutect2.py" -R "$compressed_fasta" -d "$cwd") &
wait

time python "$axolotl/annotate_serum.py" -f "$compressed_fasta" -d "$cwd"
"$axolotl/serum_filter.sh" -c "$axolotl/$dbSNP" -d "$axolotl/$COSMIC"
time python "$axolotl/merge.py" -f "$compressed_fasta" -d "$cwd/"
time python "$axolotl/dataframe.py" -f "$compressed_fasta" -b COSMIC_gatkbcftools.$accession.vcf.gz -r COSMIC_gatkfbayes.$accession.vcf.gz -l COSMIC_lofreq.$accession.vcf.gz -m COSMIC_mutect2.$accession.vcf.gz -x $accession.dedup.snps.vcf.gz
time python "$axolotl/predictor.py" -o "$high_confidence_vcf" -x "$threshold" -d $accession.dedup.snps.vcf.gz -m "$axolotl/$model"

# clean up
rm -f *.bai COSMOS_* COSMIC_* gatkbcftools*vcf* gatkfbayes*vcf* lofreq.*vcf* mutect2*vcf* bcftools*vcf* *.bqsr.bam* fbayes.*vcf* *.dedup.snps.vcf.gz* *_vcfs.list


