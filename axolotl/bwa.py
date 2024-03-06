import argparse
import concurrent.futures
import subprocess
import glob
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description="Map reads to reference and convert SAM to BAM.")
    parser.add_argument('-d', '--input_dir', default="./", help="Input directory for fastq files")
    parser.add_argument('-f', '--reference_fa', required=True, help="Reference genome in fasta format")
    parser.add_argument('-t', '--threads', default="40", help="Number of threads to use")

    return parser.parse_args()

args = parse_arguments()

# Map reads to reference
input_dir = args.input_dir
input_pattern = "*_1.fastq.gz"
reference_fa = args.reference_fa
threads = args.threads

fastq_files = glob.glob(f"{input_dir}/{input_pattern}")

def map_reads(fastq):
    fastq_prefix = fastq[:-11]  # Remove "_1.fastq.gz" from the file name
    output_sam = f"{fastq_prefix}.sam"

    bwamem = ["bwa-mem2", "mem", "-t", threads, reference_fa, fastq, f"{fastq_prefix}_2.fastq.gz"]

    print(f"Mapping reads from {fastq} to reference...", flush=True)
    with open(output_sam, "w") as outfile:
        subprocess.run(bwamem, check=True, stdout=outfile)
    print(f"Mapping of reads from {fastq} completed.", flush=True)

    # Delete fastq files
    #print(f"Deleting {fastq} and {fastq_prefix}_2.fastq.gz", flush=True)
    #os.remove(fastq)
    #os.remove(f"{fastq_prefix}_2.fastq.gz")

with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
    executor.map(map_reads, fastq_files)

print("All reads mapped to the reference.", flush=True)

# Convert sam to bam
sam_files = glob.glob("./*.sam")
print(sam_files, flush=True)

def process_sam(sam_file):
    bam_file = os.path.splitext(sam_file)[0] + ".bam"
    sorted_bam_file = os.path.splitext(sam_file)[0] + ".sort.bam"

    # Convert SAM to BAM
    samtools_convert_cmd = ["samtools", "view", "-@", threads, "-b", "-o", bam_file, sam_file]
    print(f"Converting {sam_file} to BAM...", flush=True)
    subprocess.run(samtools_convert_cmd, check=True)
    print(f"Conversion of {sam_file} to BAM completed.", flush=True)

    # Sort the BAM file
    samtools_sort_cmd = ["samtools", "sort", "-o", sorted_bam_file, bam_file]
    print(f"Sorting BAM file {bam_file}...", flush=True)
    subprocess.run(samtools_sort_cmd, check=True)
    print(f"Sorting of {bam_file} completed.", flush=True)

    # Delete the unsorted BAM file
    os.remove(bam_file)
    print(f"Deleted unsorted BAM file: {bam_file}", flush=True)

    # Delete the SAM file
    os.remove(sam_file)
    print(f"Deleted SAM file: {sam_file}", flush=True)

with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
    executor.map(process_sam, sam_files)

print("All SAM files processed.", flush=True)
