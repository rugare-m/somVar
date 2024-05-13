import argparse
import os
import glob
import subprocess
import concurrent.futures

def parse_arguments():
    parser = argparse.ArgumentParser(description="GATK Best Practices Pipeline")
    parser.add_argument('-f', '--reference_file', required=True, help="Reference genome in fasta format")
    parser.add_argument('-k', '--known_sites', required=True, help="Known sites for base recalibration in VCF format")

    return parser.parse_args()

args = parse_arguments()

reference_file = args.reference_file
known_sites = args.known_sites

# Index BAM files
bam_files = glob.glob("./*.bam") 

def index_bam(bam_file):
    print(f" Processing file: {bam_file} ", flush=True)
    # Run samtools index command
    command = f"samtools index {bam_file}"
    index_file = f"{bam_file}.bai"
    if os.path.exists(index_file):
        print(f" Index file already exists: {index_file} ", flush=True)
        return
    else:
        subprocess.run(command, shell=True)
        print(f" Finished indexing file: {bam_file} ", flush=True)

with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
    executor.map(index_bam, bam_files)

# GATK Best Practices
# 1. Remove duplicates
bam_files = glob.glob("./*.bam") 

def remove_duplicates(bam_file):
    print(f" Processing file: {bam_file} ", flush=True)
    output_file = bam_file.replace(".sort.bam", ".markdup.bam")
    metrics_file = bam_file.replace(".sort.bam", ".markdup.metrics.txt")
    command = f"gatk MarkDuplicates -I {bam_file} -O {output_file} -M {metrics_file} --REMOVE_DUPLICATES false"
    subprocess.run(command, shell=True)
    print(f" Finished processing file: {bam_file} ", flush=True)
    os.remove(bam_file)
    os.remove(f"{bam_file}.bai")
    os.remove(metrics_file)
    print(f" Deleted input BAM file: {bam_file} ", flush=True)

with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
    executor.map(remove_duplicates, bam_files)

# 2. Add read groups
bam_files = glob.glob("./*.markdup.bam")  

def add_read_groups(bam_file):
    print(f" Processing file: {bam_file} ", flush=True)
    output_file = bam_file.replace(".markdup.bam", ".rg.bam")
    command = f"gatk AddOrReplaceReadGroups -I {bam_file} -O {output_file} -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM 20"
    subprocess.run(command, shell=True)
    print(f" Finished processing file: {bam_file} ", flush=True)
    os.remove(bam_file)
    print(f" Deleted input BAM file: {bam_file} ", flush=True)

with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
    executor.map(add_read_groups, bam_files)


