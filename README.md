
### Background
This Nextflow pipeline outputs a set of high confidence ctDNA somatic variants, give paired FASTQ files as input. 


<details>
<summary>Click to expand dropdown</summary>

- Option 1
- Option 2
- Option 3

</details>


### Research Poster Abstract
Tumour variants can predict response to treatment and inform personalized treatment for cancer patients. Tissue biopsies provide the gold standard for the identification of tumour variants; however, these are typically difficult to obtain and require an invasive procedure. Circulating tumour DNA (ctDNA) is a promising, minimally invasive cancer biomarker that can be used to inform treatment of cancer patients. ctDNA is released from tumour cells into the bloodstream. ctDNA is easily accessible and contains tumour variants. Detecting real ctDNA variants with Next Generation Sequencing (NGS) technology can be a challenge due to the low abundance of ctDNA in the pool of cell free DNA in the bloodstream. Rule-based filtering strategies either remove a substantial number of true positive ctDNA variants along with false variant calls or retains an implausibly large number of total variants. Machine Learning (ML) enables identification of complex, non-linear patterns which may improve ability to distinguish between real low frequency ctDNA variants, and false positive calls arising from sequencing errors. The aim of this study is to develop a machine learning model for detecting high confidence ctDNA somatic variants in the absence of a matched tissue sample. We used a ctDNA dataset (SRP073475) sequenced with matched tissue samples (1) to train and test models. Samples in this dataset were whole exome sequenced at an average depth of 70x. Model validation included a holdout sample from the train/test datasets, and a high depth targeted sequencing dataset. 

### Pipeline for developing and evaluating ML models 
<img width="450" alt="Screenshot 2024-04-15 at 14 41 09" src="https://github.com/rugare-m/Predicting-High-Confidence-ctDNA-Somatic-Variants-with-an-Ensemble-Machine-Learning-Model/assets/88198662/361c440b-3a69-465d-afd4-dcd00618322d">
