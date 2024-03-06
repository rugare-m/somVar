#!/bin/bash

# Job name:
#SBATCH --job-name=callers
#SBATCH --mem=32G
# Wall clock limit:
#SBATCH --time=24:00:30



#time python bwa.py  -f "/users/rugarem/volatile/chapter3/data/reference/hg38" -d "/users/rugarem/volatile/chapter3/nextflow"  -t 60
#time python indexer.py -d "/users/rugarem/volatile/chapter3/nextflow"
#time python gatk.py -f "/users/rugarem/volatile/chapter3/data/reference/hg38.fa.gz" -k "1000G_omni2.5.hg38.vcf.gz"
#time python bcftools.py  -f "/users/rugarem/volatile/chapter3/data/reference/hg38.fa.gz" -d "/users/rugarem/volatile/chapter3/nextflow"
#time python freebayes.py -f "/users/rugarem/volatile/chapter3/data/reference/uhg38.fa" -d "/users/rugarem/volatile/chapter3/nextflow"
#time python lofreq.py    -f "/users/rugarem/volatile/chapter3/data/reference/hg38.fa.gz" -d "/users/rugarem/volatile/chapter3/nextflow"
#time python mutect2.py   -R "/users/rugarem/volatile/chapter3/data/reference/hg38.fa.gz" -d "/users/rugarem/volatile/chapter3/nextflow"
#time python annotate_serum.py -f "/users/rugarem/volatile/chapter3/data/reference/hg38.fa.gz" -d "/users/rugarem/volatile/chapter3/nextflow"
#./serum_filter.sh -c common_all_20180418.vcf.gz -d CosmicCodingMutsV98.vcf.gz
#time python merge.py -f "/users/rugarem/volatile/chapter3/data/reference/hg38.fa.gz" -d "/users/rugarem/volatile/chapter3/nextflow/"
time python dataframe.py -f "/users/rugarem/volatile/chapter3/data/reference/hg38.fa.gz" -b COSMIC_gatkbcftools.SRR10556218.vcf.gz -r COSMIC_gatkfbayes.SRR10556218.vcf.gz -l COSMIC_lofreq.SRR10556218.vcf.gz -m COSMIC_mutect2.SRR10556218.vcf.gz -x SRR10556218.dedup.snps.vcf.gz
time python predictor.py -o SRR10556218_predictions.vcf.gz -x 0.7 -d SRR10556218.dedup.snps.vcf.gz -m high_depth_model.pkl
