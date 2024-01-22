#works - parallel processing, 3 subprocesses at a time
import os
import glob
import subprocess
import concurrent.futures

bam_files = glob.glob("/users/rugarem/volatile/chapter3/data/erve/*/fastq/*.bam") 

def index_bam(bam_file):
    print(f" Processing file: {bam_file} ", flush=True)
    index_file = bam_file + ".bai"  # Index file name for the BAM file
    # Run samtools index command
    command = f"samtools index {bam_file}"
    subprocess.run(command, shell=True)
    print(f" Finished indexing file: {bam_file} ", flush=True)

# Index bam files
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
    executor.map(index_bam, bam_files)

bam_files = glob.glob("/users/rugarem/volatile/chapter3/data/erve/*/fastq/*.bam")  
reference_file = "/users/rugarem/volatile/chapter3/data/reference/hg38.fa.gz"

def remove_duplicates(bam_file):
    print(f" Processing file: {bam_file} ", flush=True)
    output_file = bam_file.replace(".sort.bam", ".markdup.bam")
    metrics_file = bam_file.replace(".sort.bam", ".markdup.metrics.txt")
    command = f"gatk MarkDuplicates -I {bam_file} -O {output_file} -M {metrics_file}"
    subprocess.run(command, shell=True)
    print(f" Finished processing file: {bam_file} ", flush=True)
    os.remove(bam_file)
    os.remove(bam_file + ".bai") 
    os.remove(metrics_file)  # Delete the input BAM file
    print(f" Deleted input BAM file: {bam_file} ", flush=True)

# gatk best practices
##1, remove duplicates
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
    executor.map(remove_duplicates, bam_files)

def add_read_groups(bam_file):
    print(f" Processing file: {bam_file} ", flush=True)
    output_file = bam_file.replace(".markdup.bam", ".rg.bam")
    command = f"gatk AddOrReplaceReadGroups -I {bam_file} -O {output_file} -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM 20"
    subprocess.run(command, shell=True)
    print(f" Finished processing file: {bam_file} ", flush=True)
    os.remove(bam_file)  # Delete the input BAM file
    print(f" Deleted input BAM file: {bam_file} ", flush=True)

#2 add read groups
bam_files = glob.glob("/users/rugarem/volatile/chapter3/data/erve/*/fastq/*.markdup.bam")  
with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
    executor.map(add_read_groups, bam_files)

def base_recalibration(bam_file):
    known_sites = "/users/rugarem/volatile/chapter3/data/1000G_omni2.5.hg38.vcf.gz"
    
    print(f" Processing file: {bam_file} ", flush=True)
    output_file = bam_file.replace(".rg.bam", ".table")
    command = f"gatk BaseRecalibrator -I {bam_file} -O {output_file} -R {reference_file} --known-sites {known_sites}"
    subprocess.run(command, shell=True)   
    print(f" Finished generating recalibration file: {bam_file} ", flush=True)
    
#3, base recalibration
bam_files = glob.glob("/users/rugarem/volatile/chapter3/data/erve/*/fastq/*.rg.bam")
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
    executor.map(base_recalibration, bam_files)

def apply_recalibration(bam_file):
    
#4, apply recalibration
    print(f" Applying recalibration to: {bam_file} ", flush=True)
    output_bam = bam_file.replace(".rg.bam", ".bqsr.bam")
    recalibration_file = bam_file.replace(".rg.bam", ".table")
    command = f"gatk ApplyBQSR -R {reference_file} -I {bam_file} --bqsr-recal-file {recalibration_file} -O {output_bam}"
    subprocess.run(command, shell=True)
    
# deleteing bam file and calibration file...
    
    os.remove(bam_file) 
    os.remove(recalibration_file)
    
# processing bam file complete!
    
    print (" processing bam file complete! ", flush=True)

bam_files = glob.glob("/users/rugarem/volatile/chapter3/data/erve/*/fastq/*.rg.bam")
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        executor.map(apply_recalibration, bam_files)
