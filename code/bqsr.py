import argparse
import os
import glob
import subprocess
import concurrent.futures

def parse_arguments():
    parser = argparse.ArgumentParser(description="GATK Best Practices Pipeline")
    parser.add_argument('-f', '--reference_file', required=True, help="Reference genome in fasta format")
    parser.add_argument('-k', '--known_sites', required=True, help="Known sites for base recalibration in VCF format")
    parser.add_argument('-b', '--bam_files', required=True, help="BAM files to process")

    return parser.parse_args()

args = parse_arguments()

reference_file = args.reference_file
known_sites = args.known_sites


# 3. Base recalibration
bam_files = args.bam_files.split(",")

def base_recalibration(bam_file):
    print(f" Processing file: {bam_file} ", flush=True)
    output_file = bam_file.replace(".rg.bam", ".table")
    command = f"gatk BaseRecalibrator -I {bam_file} -O {output_file} -R {reference_file} --known-sites {known_sites}"
    subprocess.run(command, shell=True)   
    print(f" Finished generating recalibration file: {bam_file} ", flush=True)

with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
    executor.map(base_recalibration, bam_files)

## 4. Apply recalibration
#def apply_recalibration(bam_file):
#    print(f" Applying recalibration to: {bam_file} ", flush=True)
#    output_bam = bam_file.replace(".rg.bam", ".bqsr.bam")
#    recalibration_file = bam_file.replace(".rg.bam", ".table")
#    command = f"gatk ApplyBQSR -R {reference_file} -I {bam_file} --bqsr-recal-file {recalibration_file} -O {output_bam}"
#    subprocess.run(command, shell=True)
#    os.remove(bam_file) 
#    os.remove(recalibration_file)    
#    print(" Processing bam file complete! ", flush=True)
#    subprocess.run(f"samtools index {output_bam}", shell=True)
# index bam files


#with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
#        executor.map(apply_recalibration, bam_files)
