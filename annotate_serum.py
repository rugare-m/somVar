import subprocess
import os
import glob
from multiprocessing import Pool
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Variant Annotation Pipeline")
    parser.add_argument('-f', '--reference_file', required=True, help="Reference genome in fasta format")
    parser.add_argument('-d', '--dirs', nargs='+', required=True, help="List of input directories containing VCF files")

    return parser.parse_args()

def process_vcf_files(vcf_files, output_dir, reference):
    for vcf_file in vcf_files:
        tabix = ["tabix", "-p", "vcf", "-f", vcf_file]
        subprocess.run(tabix)
    
        print(f"File {vcf_file} indexed")
        # Extract the sample name from the vcf file name
        sample_name = os.path.basename(vcf_file).split('.')[0]+"."+os.path.basename(vcf_file).split('.')[1]
        output_file = f"gatk{sample_name}.vcf.gz"
        
        out_path = os.path.join(output_dir, output_file)
        bam = os.path.basename(vcf_file).split('.')[1] + ".bqsr.bam"
        bam_path = os.path.join(output_dir, bam)
        
        annotate = ["gatk", "VariantAnnotator", "-R", reference, "-I", bam_path, "-V", vcf_file, "-O", out_path, "-A", "AlleleFraction"]
        subprocess.run(annotate)

if __name__ == "__main__":
    args = parse_arguments()
    reference_file = args.reference_file
    dirs = args.dirs

    for directory in dirs:
        directory_path = directory
        vcf_files_bwa2 = glob.glob(os.path.join(directory_path, "bcftools*.vcf.gz"))
        vcf_files_fbayes = glob.glob(os.path.join(directory_path, "fbayes*.vcf.gz"))
        
        # Create a Pool for parallel processing
        with Pool(processes=2) as pool:
            pool.starmap(process_vcf_files, [(vcf_files_bwa2, directory_path, reference_file), (vcf_files_fbayes, directory_path, reference_file)])

