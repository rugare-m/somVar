#!/bin/bash

#SBATCH --job-name=tissue_filter           # Job name
#SBATCH --mem=50gb
#SBATCH --time=48:05:00              # Time limit hrs:min:sec
#SBATCH --output=parallel_%j.log     # Standard output and error log


#bcftools filter
for vcf in /users/rugarem/volatile/chapter3/data/erve/tissue/fastq/gatk*bcftools*vcf.gz
do
tabix -f $vcf
#filtering for DP and QUAL
bcftools filter -O z -o `dirname $vcf`/filter_`basename $vcf` -i 'INFO/DP>70 & QUAL>=50' $vcf
bcftools filter -O z -o `dirname $vcf`/lowerx_`basename $vcf` -i 'FORMAT/AF < 0.4 | FORMAT/AF > 0.6' `dirname $vcf`/filter_`basename $vcf` # include less than 45% 
bcftools filter -O z -o `dirname $vcf`/upper_`basename $vcf` -i 'FORMAT/AF < 0.90' `dirname $vcf`/lowerx_`basename $vcf` #exclude 1 ish germlines
#annotate with dbsnp (all common SNPs)
java -Xmx8G -jar SnpSift.jar annotate -id /users/rugarem/volatile/barkla/snpEff/databases/common_all_20180418.vcf.gz `dirname $vcf`/upper_`basename $vcf` >  `dirname $vcf`/annotated_`basename $vcf .gz`
bgzip `dirname $vcf`/annotated_`basename $vcf .gz`
#filter out all common snps
bcftools filter -O z -o `dirname $vcf`/dbsnp_`basename $vcf` -e 'ID~"rs"' `dirname $vcf`/annotated_`basename $vcf`
#Annotate for COSMIC presence
java -Xmx8G -jar SnpSift.jar annotate -id -noInfo /users/rugarem/volatile/barkla/snpEff/databases/CosmicCodingMutsV98.vcf.gz `dirname $vcf`/dbsnp_`basename $vcf` > `dirname $vcf`/COSMIC_`basename $vcf .gz`
bgzip `dirname $vcf`/COSMIC_`basename $vcf .gz`
#keep COSMIC variants
bcftools view -O z -i 'ID~"COSV"' `dirname $vcf`/COSMIC_`basename $vcf` > `dirname $vcf`/complete_`basename $vcf`
#delete files
rm `dirname $vcf`/annotated_`basename $vcf`
rm `dirname $vcf`/filter_`basename $vcf`
rm `dirname $vcf`/dbsnp_`basename $vcf`
#rm `dirname $vcf`/COSMIC_`basename $vcf`
rm `dirname $vcf`/upper_`basename $vcf`
#rm `dirname $vcf`/lower_`basename $vcf`
rm `dirname $vcf`/lowerx_`basename $vcf`
done

#freebayes 
for vcf in /users/rugarem/volatile/chapter3/data/erve/tissue/fastq/gatk*fbayes*.vcf.gz
do
tabix -f $vcf
#filtering for DP and QUAL
bcftools filter -O z -o `dirname $vcf`/filter_`basename $vcf` -i 'INFO/DP>70 & QUAL>=50' $vcf
bcftools filter -O z -o `dirname $vcf`/lowerx_`basename $vcf` -i 'FORMAT/AF < 0.4 | FORMAT/AF > 0.6' `dirname $vcf`/filter_`basename $vcf` # include less than 45% 
bcftools filter -O z -o `dirname $vcf`/upper_`basename $vcf` -i 'FORMAT/AF < 0.90' `dirname $vcf`/lowerx_`basename $vcf` #exclude 1 ish germlines
#annotate with dbsnp (all common SNPs)
java -Xmx8G -jar SnpSift.jar annotate -id /users/rugarem/volatile/barkla/snpEff/databases/common_all_20180418.vcf.gz `dirname $vcf`/upper_`basename $vcf` >  `dirname $vcf`/annotated_`basename $vcf .gz`
bgzip `dirname $vcf`/annotated_`basename $vcf .gz`
#filter out all common snps
bcftools filter -O z -o `dirname $vcf`/dbsnp_`basename $vcf` -e 'ID~"rs"' `dirname $vcf`/annotated_`basename $vcf`
#Annotate for COSMIC presence
java -Xmx8G -jar SnpSift.jar annotate -id -noInfo /users/rugarem/volatile/barkla/snpEff/databases/CosmicCodingMutsV98.vcf.gz `dirname $vcf`/dbsnp_`basename $vcf` > `dirname $vcf`/COSMIC_`basename $vcf .gz`
bgzip `dirname $vcf`/COSMIC_`basename $vcf .gz`
#keep COSMIC variants
bcftools view -O z -i 'ID~"COSV"' `dirname $vcf`/COSMIC_`basename $vcf` > `dirname $vcf`/complete_`basename $vcf`
#delete files
rm `dirname $vcf`/annotated_`basename $vcf`
rm `dirname $vcf`/filter_`basename $vcf`
rm `dirname $vcf`/dbsnp_`basename $vcf`
#rm `dirname $vcf`/COSMIC_`basename $vcf`
rm `dirname $vcf`/upper_`basename $vcf`
#rm `dirname $vcf`/lower_`basename $vcf`
rm `dirname $vcf`/lowerx_`basename $vcf`
done

#mutect2
for vcf in /users/rugarem/volatile/chapter3/data/erve/tissue/fastq/mutect2*.vcf.gz
do
tabix -f $vcf
#filtering for DP and QUAL
bcftools filter -O z -o `dirname $vcf`/filter_`basename $vcf` -i 'INFO/DP>70' $vcf
bcftools filter -O z -o `dirname $vcf`/lowerx_`basename $vcf` -i 'FORMAT/AF < 0.4 | FORMAT/AF > 0.6' `dirname $vcf`/filter_`basename $vcf` # include less than 45% 
bcftools filter -O z -o `dirname $vcf`/upper_`basename $vcf` -i 'FORMAT/AF < 0.90' `dirname $vcf`/lowerx_`basename $vcf` #exclude 1 ish germlines
#annotate with dbsnp (5% + in all dbsnp populations)
java -Xmx8G -jar SnpSift.jar annotate -id /users/rugarem/volatile/barkla/snpEff/databases/common_all_20180418.vcf.gz `dirname $vcf`/upper_`basename $vcf` >  `dirname $vcf`/annotated_`basename $vcf .gz`
bgzip `dirname $vcf`/annotated_`basename $vcf .gz`
#filter out all common snps
bcftools filter -O z -o `dirname $vcf`/dbsnp_`basename $vcf` -e 'ID~"rs"' `dirname $vcf`/annotated_`basename $vcf`
#Annotate for COSMIC presence
java -Xmx8G -jar SnpSift.jar annotate -id -noInfo /users/rugarem/volatile/barkla/snpEff/databases/CosmicCodingMutsV98.vcf.gz `dirname $vcf`/dbsnp_`basename $vcf` > `dirname $vcf`/COSMIC_`basename $vcf`
#keep COSMIC variants
bcftools view -O z -i 'ID~"COSV"' `dirname $vcf`/COSMIC_`basename $vcf` > `dirname $vcf`/complete_`basename $vcf`
rm `dirname $vcf`/annotated_`basename $vcf`
rm `dirname $vcf`/filter_`basename $vcf`
rm `dirname $vcf`/dbsnp_`basename $vcf`
#rm `dirname $vcf`/COSMIC_`basename $vcf`
rm `dirname $vcf`/upper_`basename $vcf`
#rm `dirname $vcf`/lower_`basename $vcf`
rm `dirname $vcf`/lowerx_`basename $vcf`
done

for vcf in /users/rugarem/volatile/chapter3/data/erve/tissue/fastq/lofreq*vcf.gz
do
tabix -f $vcf
#filtering for DP and QUAL
bcftools filter -O z -o `dirname $vcf`/filter_`basename $vcf` -i 'INFO/DP> 70 & QUAL>=50' $vcf
bcftools filter -O z -o `dirname $vcf`/lowerx_`basename $vcf` -i 'INFO/AF < 0.4 | INFO/AF > 0.6' `dirname $vcf`/filter_`basename $vcf` # remove 0.5 germliines 
bcftools filter -O z -o `dirname $vcf`/upper_`basename $vcf` -i 'INFO/AF < 0.90' `dirname $vcf`/lowerx_`basename $vcf` #exclude 1 ish germlines

#annotate with dbsnp (all common SNPs)
java -Xmx8G -jar SnpSift.jar annotate -id /users/rugarem/volatile/barkla/snpEff/databases/common_all_20180418.vcf.gz `dirname $vcf`/upper_`basename $vcf` >  `dirname $vcf`/annotated_`basename $vcf .gz`
bgzip `dirname $vcf`/annotated_`basename $vcf .gz`

#filter out all common snps
bcftools filter -O z -o `dirname $vcf`/dbsnp_`basename $vcf` -e 'ID~"rs"' `dirname $vcf`/annotated_`basename $vcf` 

#Annotate for COSMIC presence
java -Xmx8G -jar SnpSift.jar annotate -id -noInfo /users/rugarem/volatile/barkla/snpEff/databases/CosmicCodingMutsV98.vcf.gz `dirname $vcf`/dbsnp_`basename $vcf` > `dirname $vcf`/COSMO_`basename $vcf .gz`
bgzip `dirname $vcf`/COSMO_`basename $vcf .gz`
#keep COSMIC variants
bcftools view -O z -i 'ID~"COSV"' `dirname $vcf`/COSMO_`basename $vcf` > `dirname $vcf`/comp_`basename $vcf`

#add synthetic sample name
/users/rugarem/volatile/chapter3/data/dietz/serum/bwa2/lofreq_gt.sh -i `dirname $vcf`/comp_`basename $vcf` -g 1/1 -n 20 -o `dirname $vcf`/complete_`basename $vcf`

#delete files
rm `dirname $vcf`/annotated_`basename $vcf`
rm `dirname $vcf`/filter_`basename $vcf`
rm `dirname $vcf`/upper_`basename $vcf`
rm `dirname $vcf`/lowerx_`basename $vcf`
#rm `dirname $vcf`/comp_`basename $vcf`
rm `dirname $vcf`/COSMO_`basename $vcf`
rm `dirname $vcf`/dbsnp_`basename $vcf`

done


