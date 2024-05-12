
### Background
This Nextflow pipeline outputs a set of high confidence ctDNA somatic variants, given paired FASTQ files as input. The pipeline maps reads to the reference using bwa2, calls variants with bcftools, FreeBayes, LoFreq and Mutect2, merges the output FASTQs, and finally predicts the high confidence variants using a Random Forest model. The tool was fitted on WES data and assumes WES data is being passed in. There are two models to choose from, a low depth model for data sequenced at ~10X, and a high depth model for data sequenced at ~200X. 


### Requirements

<details>
<summary>Reference Genomes</summary>

 - 1x GRCh38 reference genome, compressed with bgzip.This genome must be indexed with bwa2, and samtools faidx & we'll need a dictionary created with gatk CreateSequenceDictionary. GATK, samtools and bwa2 will be installed in the virtual environment further on!

 - 1x uncompressed GRCh38 reference genome. The uncompressed genome will also need to be indexed with gatk CreateSequenceDictionary dict and samtools faidx
</details>


<details>
<summary>Databases</summary>
 - 1000G_omni2.5.hg38.vcf.gz
 - CosmicCodingMutsV98.vcf.gz
 - dbSNP common variants vcf
</details>


### Usage
First we need to setup the virtual environment with all the requirements:

```bash
conda env create -f somvar.yml
```
