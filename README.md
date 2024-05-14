
### Background
This Nextflow pipeline outputs a set of high confidence ctDNA somatic variants, given paired FASTQ files as input. The pipeline maps reads to the reference using bwa2, calls variants with bcftools, FreeBayes, LoFreq and Mutect2, merges the output VCF files, and finally predicts the high confidence variants using a Random Forest model. The tool was fitted on WES data and assumes WES data is being passed in. There are two models to choose from, a low depth model for data sequenced at ~10X, and a high depth model for data sequenced at ~200X. 

### Requirements

The latest version of conda must be installed - tested version is conda 23.3.1

<details>
<summary>Reference Genomes</summary>
It's recommended to store the reference files in the 'references' directory

 - 1x GRCh38 reference genome, compressed with bgzip.This genome must be indexed with bwa2, and samtools faidx & we'll need a dict created with gatk CreateSequenceDictionary. GATK, samtools and bwa2 will be installed in the virtual environment further on! 

 - 1x uncompressed GRCh38 reference genome. The uncompressed genome will also need to be indexed with gatk CreateSequenceDictionary dict and samtools faidx
</details>

<details>
<summary>Databases</summary>
 
The following VCF files **must** be stored in the 'databases' directory
 
 - 1000G_omni2.5.hg38.vcf.gz

 - CosmicCodingMutsV98.vcf.gz

 - dbSNP common variants vcf
   
</details>

### Usage
First we need to setup the virtual environment with all the dependencies:

```bash
conda env create -f somvar.yml
```


#### Next, we'll have to modify the params.json file
genome = path to reference genome

hg38 = path to reference genome bwa index prefix, eg "path/to/reference/bwa/index/hg38

known = path to 1000G_omni2.5.hg38.vcf.gz

cosmic = path to CosmicCodingMutsV98.vcf.gz

dbsnp = path to dbSNP common variants

model = path to the model for predicting high confidence variants  ["high_depth_model.pkl" OR "low_depth_model.pkl"]

threshold = probability threshold to predict high confidence variants - default is 0.75

output = path to name of output file 

### Running the pipeline

The pipeline expects FASTQ files compressed with bgzip in the somVar directory. 

The pipeline expects FASTQs with the naming in this format ACCESSION_1/2.fastq.gz. 

You can find the path to the virtual environment with:

```bash
conda env list
```
Then:
```bash
nextflow run main.nf -resume -with-conda /path/to/virtual/env/somvar -params-file params.json --f1 SRR000000_1.fastq.gz --f2 SRR000000_2.fastq.gz
```
