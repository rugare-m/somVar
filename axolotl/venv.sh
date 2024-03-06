# Create Conda virtual environment
conda create --name axolotl python=3.9.13

# Activate  virtual environment
conda activate axolotl

# Install the required packages
conda install bwa-mem2=2.2.1 gatk4=4.3.0.0 bcftools=1.15 samtools=1.15 freebayes=1.3.6 lofreq=2.1.5 sra-tools=3.0.5
pip install joblib==1.3.2 pandas==2.2.1 scikit-learn==1.3.2 pysam==0.20.0 biopython==1.79 cyvcf2==0.30.18 pyfastx==1.1.0

#pip install pandas==2.2.1
#pip install scikit-learn==1.3.2
#pip install pysam==0.20.0
#pip install biopython==1.79
#pip install cyvcf2==0.30.18
#pip install pyfastx==1.1.0
