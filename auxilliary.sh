#!/bin/bash

# Job name:
#SBATCH --job-name=fastq
#SBATCH --mem=80G
# Wall clock limit:
#SBATCH --time=24:00:30


cwd=$(pwd)

#reference genome
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
gunzip hg38.fa.gz
bgzip hg38.fa
gatk CreateSequenceDictionary -R hg38.fa.gz -O "$cwd/hg38.dict"
samtools faidx "$cwd/hg38.fa.gz"
bwa-mem2 index -p hg38 "$cwd/hg38.fa.gz"
mv hg38* "$cwd/fasta/"

cp "$cwd/fasta/hg38.fa.gz" "$cwd/fasta/uhg38.fa.gz"
gunzip "$cwd/fasta/uhg38.fa.gz"
samtools faidx "$cwd/fasta/uhg38.fa"

#Omni100G
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz
mv 1000G_omni2.5.hg38.vcf.gz "$cwd/fastq"

#dbSNP
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz
mv common_all_20180418.vcf.gz "$cwd/fastq"

# Example FASTQ
prefetch SRR10556218
fasterq-dump SRR10556218
rm "$cwd/SRR10556218.fastq"
rm -r "$cwd/SRR10556218"
bgzip "$cwd/SRR10556218_1.fastq"
bgzip "$cwd/SRR10556218_2.fastq"
