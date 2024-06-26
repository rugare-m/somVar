#!/bin/bash

while getopts ":c:d:b:l:f:m:p:" opt; do
  case $opt in
    c)
      common_vcf="$OPTARG"
      ;;
    d)
      cosmic_vcf="$OPTARG"
      ;;
    b)
      bcftools="$OPTARG"
      ;;
    l)
      lofreq="$OPTARG"
      ;;
    f)
      freebayes="$OPTARG"
      ;;
    m)
      mutect2="$OPTARG"
      ;;
    p)
      path="$OPTARG"
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




#lofreq
for vcf in $lofreq
do
tabix -f $vcf
#filtering for DP and QUAL
bcftools filter -O z -o `dirname $vcf`/filter_`basename $vcf` -i 'INFO/DP> 0 & QUAL>=0' $vcf
bcftools filter -O z -o `dirname $vcf`/lowerx_`basename $vcf` -i 'INFO/AF < 0.4 | INFO/AF > 0.6' `dirname $vcf`/filter_`basename $vcf` # remove 0.5 germliines 
bcftools filter -O z -o `dirname $vcf`/upper_`basename $vcf` -i 'INFO/AF < 0.90' `dirname $vcf`/lowerx_`basename $vcf` #exclude 1 ish germlines

#annotate with dbsnp (all common SNPs)
SnpSift annotate -id "$common_vcf" `dirname $vcf`/upper_`basename $vcf` >  `dirname $vcf`/annotated_`basename $vcf .gz`
bgzip -f `dirname $vcf`/annotated_`basename $vcf .gz`

#Annotate for COSMIC presence
SnpSift annotate -id -noInfo "$cosmic_vcf" `dirname $vcf`/annotated_`basename $vcf` > `dirname $vcf`/COSMO_`basename $vcf .gz`
bgzip -f `dirname $vcf`/COSMO_`basename $vcf .gz`

#add synthetic sample name
$path/lofreq_gt.sh -i `dirname $vcf`/COSMO_`basename $vcf` -g 1/1 -n 20 -o `dirname $vcf`/COSMOS_`basename $vcf`

#remove contig names in vcf header
gzip -cd `dirname $vcf`/COSMOS_`basename $vcf` | grep -v '^##contig=' > `dirname $vcf`/COSMIC_`basename $vcf .gz`
bgzip -f `dirname $vcf`/COSMIC_`basename $vcf .gz`

#delete files
rm `dirname $vcf`/annotated_`basename $vcf`
rm `dirname $vcf`/filter_`basename $vcf`
rm `dirname $vcf`/upper_`basename $vcf`
rm `dirname $vcf`/lowerx_`basename $vcf`
rm `dirname $vcf`/COSMO_`basename $vcf`
rm `dirname $vcf`/COSMOS_`basename $vcf`

done

#bcftools filter
for vcf in $bcftools
do
tabix -f $vcf
#filtering for DP and QUAL
bcftools filter -O z -o `dirname $vcf`/filter_`basename $vcf` -i 'INFO/DP>0 & QUAL>=0' $vcf
bcftools filter -O z -o `dirname $vcf`/lowerx_`basename $vcf` -i 'FORMAT/AF < 0.4 | FORMAT/AF > 0.6' `dirname $vcf`/filter_`basename $vcf` # include less than 45% 
bcftools filter -O z -o `dirname $vcf`/upper_`basename $vcf` -i 'FORMAT/AF < 0.90' `dirname $vcf`/lowerx_`basename $vcf` #exclude 1 ish germlines

#annotate with dbsnp (all common SNPs)
SnpSift annotate -id "$common_vcf" `dirname $vcf`/upper_`basename $vcf` >  `dirname $vcf`/annotated_`basename $vcf .gz`
bgzip -f `dirname $vcf`/annotated_`basename $vcf .gz`

#Annotate for COSMIC presence
SnpSift annotate -id -noInfo "$cosmic_vcf" `dirname $vcf`/annotated_`basename $vcf` > `dirname $vcf`/COSMIC_`basename $vcf .gz`
bgzip -f `dirname $vcf`/COSMIC_`basename $vcf .gz`

rm `dirname $vcf`/annotated_`basename $vcf`
rm `dirname $vcf`/filter_`basename $vcf`
rm `dirname $vcf`/upper_`basename $vcf`
rm `dirname $vcf`/lowerx_`basename $vcf`
done

#freebayes 
for vcf in $freebayes
do
tabix -f $vcf
#filtering for DP and QUAL
bcftools filter -O z -o `dirname $vcf`/filter_`basename $vcf` -i 'INFO/DP>0 & QUAL>=0' $vcf
bcftools filter -O z -o `dirname $vcf`/lowerx_`basename $vcf` -i 'FORMAT/AF < 0.4 | FORMAT/AF > 0.6' `dirname $vcf`/filter_`basename $vcf` # include less than 45% 
bcftools filter -O z -o `dirname $vcf`/upper_`basename $vcf` -i 'FORMAT/AF < 0.90' `dirname $vcf`/lowerx_`basename $vcf` #exclude 1 ish germlines

#annotate with dbsnp (all common SNPs)
SnpSift annotate -id "$common_vcf" `dirname $vcf`/upper_`basename $vcf` >  `dirname $vcf`/annotated_`basename $vcf .gz`
bgzip -f `dirname $vcf`/annotated_`basename $vcf .gz`

#Annotate for COSMIC presence
SnpSift annotate -id -noInfo "$cosmic_vcf" `dirname $vcf`/annotated_`basename $vcf` > `dirname $vcf`/COSMIC_`basename $vcf .gz`
bgzip -f `dirname $vcf`/COSMIC_`basename $vcf .gz`

#delete files
rm `dirname $vcf`/annotated_`basename $vcf`
rm `dirname $vcf`/filter_`basename $vcf`
rm `dirname $vcf`/upper_`basename $vcf`
rm `dirname $vcf`/lowerx_`basename $vcf`
done


#mutect2
for vcf in $mutect2
do
tabix -f $vcf
#filtering for DP and QUAL
bcftools filter -O z -o `dirname $vcf`/filter_`basename $vcf` -i 'INFO/DP>0' $vcf
bcftools filter -O z -o `dirname $vcf`/lowerx_`basename $vcf` -i 'FORMAT/AF < 0.4 | FORMAT/AF > 0.6' `dirname $vcf`/filter_`basename $vcf` # include less than 45% 
bcftools filter -O z -o `dirname $vcf`/upper_`basename $vcf` -i 'FORMAT/AF < 0.90' `dirname $vcf`/lowerx_`basename $vcf` #exclude 1 ish germlines


#annotate with dbsnp (5% + in all dbsnp populations)
SnpSift annotate -id "$common_vcf" `dirname $vcf`/upper_`basename $vcf` >  `dirname $vcf`/annotated_`basename $vcf .gz`
bgzip -f `dirname $vcf`/annotated_`basename $vcf .gz`

#Annotate for COSMIC presence
SnpSift annotate -id -noInfo "$cosmic_vcf" `dirname $vcf`/annotated_`basename $vcf` > `dirname $vcf`/COSMIC_`basename $vcf .gz`
bgzip -f `dirname $vcf`/COSMIC_`basename $vcf .gz`

rm `dirname $vcf`/annotated_`basename $vcf`
rm `dirname $vcf`/filter_`basename $vcf`
rm `dirname $vcf`/upper_`basename $vcf`
rm `dirname $vcf`/lowerx_`basename $vcf`
done

