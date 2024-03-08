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

def process_vcf_files(chunk, output_dir, reference):
    for vcf_file in chunk:
        tabix = ["tabix", "-p", "vcf", "-f", vcf_file]
        subprocess.run(tabix)
    
        print(f"File {vcf_file} indexed")
        # Extract the sample name from the vcf file name
        sample_name = os.path.basename(vcf_file).split('.')[0]+"."+os.path.basename(vcf_file).split('.')[1]
        output_file = f"gatk{sample_name}.vcf.gz"
        
        out_path = os.path.join(output_dir, output_file)
        bam = os.path.basename(vcf_file).split('.')[1] + ".bqsr.bam"
        bam_path = os.path.join(output_dir, bam)
        
        annotate = ["gatk", "VariantAnnotator", "-R", reference, "-I", bam_path, "-V", vcf_file, "-O", out_path, "-A", "q"]
        subprocess.run(annotate)

if __name__ == "__main__":
    args = parse_arguments()
    reference_file = args.reference_file
    dirs = args.dirs

    for directory in dirs:
        directory_path = directory
        vcf_files = glob.glob(os.path.join(directory_path, "*.vcf.gz"))

        # Check if there are at least 4 VCF files in the directory
        if len(vcf_files) < 4:
            print(f"Error: Directory {directory_path} does not contain enough VCF files for parallel processing (at least 4 required).")
            continue

        # Split the list of VCF files into chunks of 4
        vcf_chunks = [vcf_files[i:i + 4] for i in range(0, len(vcf_files), 4)]

        # Create a Pool for parallel processing
        with Pool() as pool:
            pool.starmap(process_vcf_files, [(chunk, directory_path, reference_file) for chunk in vcf_chunks])


