#!/bin/bash

while getopts ":c:d:" opt; do
  case $opt in
    c)
      common_vcf="$OPTARG"
      ;;
    d)
      cosmic_vcf="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Check if the required options are provided
if [ -z "$common_vcf" ] || [ -z "$cosmic_vcf" ]; then
  echo "Usage: $0 -c <common_vcf> -d <cosmic_vcf>"
  exit 1
fi


cwd=$(pwd)


for vcf in "$cwd"/*vcf.gz
do
tabix -f $vcf
#filtering for DP and QUAL
bcftools filter -O z -o `dirname $vcf`/filter_`basename $vcf` -i 'INFO/DP>0 & QUAL>=0' $vcf
bcftools filter -O z -o `dirname $vcf`/lowerx_`basename $vcf` -i 'FORMAT/AF < 0.4 | FORMAT/AF > 0.6' `dirname $vcf`/filter_`basename $vcf` # include less than 45% 
bcftools filter -O z -o `dirname $vcf`/upper_`basename $vcf` -i 'FORMAT/AF < 0.90' `dirname $vcf`/lowerx_`basename $vcf` #exclude 1 ish germlines

#annotate with dbsnp (all common SNPs)
SnpSift annotate -v -id "$common_vcf" `dirname $vcf`/upper_`basename $vcf` >  `dirname $vcf`/annotated_`basename $vcf .gz`
bgzip `dirname $vcf`/annotated_`basename $vcf .gz`

#Annotate for COSMIC presence
SnpSift annotate -v -id -noInfo "$cosmic_vcf" `dirname $vcf`/annotated_`basename $vcf` > `dirname $vcf`/COSMIC_`basename $vcf .gz`
bgzip `dirname $vcf`/COSMIC_`basename $vcf .gz`

rm `dirname $vcf`/annotated_`basename $vcf`
rm `dirname $vcf`/filter_`basename $vcf`
rm `dirname $vcf`/upper_`basename $vcf`
rm `dirname $vcf`/lowerx_`basename $vcf`
done
