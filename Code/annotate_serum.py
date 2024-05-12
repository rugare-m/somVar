import subprocess
import os
import glob
from multiprocessing import Pool
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Variant Annotation Pipeline")
    parser.add_argument('-f', '--reference_file', required=True, help="Reference genome in fasta format")
    parser.add_argument('-v', '--vcf', required=True, help="VCF file")
    parser.add_argument('-b', '--bam', required=True, help="Single BAM file")

    return parser.parse_args()

def process_vcf_files(vcf, reference, bam):
    tabix = ["tabix", "-p", "vcf", "-f", vcf]
    subprocess.run(tabix)
    
    print(f"File {vcf} indexed")
    # Extract the sample name from the vcf file name
    sample_name = os.path.basename(vcf).split('.')[0]+"."+os.path.basename(vcf).split('.')[1]
    output_file = f"gatk{sample_name}.vcf.gz"
    
    bam = bam
       
    annotate = ["gatk", "VariantAnnotator", "-R", reference, "-I", bam, "-V", vcf, "-O", output_file, "-A", "AlleleFraction"]
    subprocess.run(annotate)

args = parse_arguments()

vcf = args.vcf
reference = args.reference_file
bam = args.bam

process_vcf_files(vcf, reference, bam)
