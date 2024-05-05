import argparse
import concurrent.futures
import subprocess
import os

# Parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Map paired-end reads to reference genome, convert SAM to BAM, and sort BAM files.")
    parser.add_argument('-r', '--reference_fa', required=True, help="Path to reference genome in FASTA format")
    parser.add_argument('-t', '--threads', default=40, type=int, help="Number of threads to use")
    parser.add_argument('-f1', '--fastq1', required=True, help="Path to the first FASTQ file of the paired-end reads")
    parser.add_argument('-f2', '--fastq2', required=True, help="Path to the second FASTQ file of the paired-end reads")

    return parser.parse_args()

args = parse_arguments()

# Extract prefix from the first FASTQ filename to use as output prefix
output_prefix = os.path.basename(args.fastq1).rsplit("_1.fastq.gz", 1)[0]

# Function to map reads to reference genome
def map_reads():
    output_sam = f"{output_prefix}.sam"

    bwamem = ["bwa-mem2", "mem", "-t", str(args.threads), args.reference_fa, args.fastq1, args.fastq2]

    print(f"Mapping reads from {args.fastq1} and {args.fastq2} to reference...", flush=True)
    with open(output_sam, "w") as outfile:
        subprocess.run(bwamem, check=True, stdout=outfile)
    print(f"Mapping completed for {args.fastq1} and {args.fastq2}.", flush=True)

    return output_sam

# Function to convert SAM to BAM and sort BAM files
def process_sam(sam_file):
    bam_file = os.path.splitext(sam_file)[0] + ".bam"
    sorted_bam_file = os.path.splitext(sam_file)[0] + ".sort.bam"

    # Convert SAM to BAM
    samtools_convert_cmd = ["samtools", "view", "-@", str(args.threads), "-b", "-o", bam_file, sam_file]
    subprocess.run(samtools_convert_cmd, check=True)

    # Sort BAM
    samtools_sort_cmd = ["samtools", "sort", "-@", str(args.threads), "-o", sorted_bam_file, bam_file]
    subprocess.run(samtools_sort_cmd, check=True)

    # Delete unsorted BAM and SAM files to save space
    os.remove(bam_file)  # Delete unsorted BAM
    os.remove(sam_file)  # Delete SAM

    return sorted_bam_file

# Map reads to generate SAM files
sam_file = map_reads()

# Process SAM to create sorted BAM
sorted_bam_file = process_sam(sam_file)

print("All operations completed successfully.", flush=True)
