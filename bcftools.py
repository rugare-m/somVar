import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="BCFTools Variant Calling Pipeline")
    parser.add_argument('-f', '--reference_file', required=True, help="Reference genome in fasta format")
    parser.add_argument('-b', '--bam', required=True, help="Single BAM file")

    return parser.parse_args()

args = parse_arguments()

reference_file = args.reference_file

def process_file(alignments_file):
    sample_name = os.path.basename(alignments_file).split(".")[0]
    output_directory = os.path.dirname(alignments_file)
    output_file = os.path.join(output_directory, f"bcftools.{sample_name}.vcf.gz")
    log_file = os.path.join(output_directory, f"bcftools.{sample_name}.log")  # Log file name for caller output

    command1 = ["bcftools", "mpileup", "-Ob", "-f", reference_file, alignments_file]
    command2 = ["bcftools", "call", "-mv", "-Oz", "-o", output_file]

    with open(log_file, "w") as log:
        process1 = subprocess.Popen(command1, stdout=subprocess.PIPE, stderr=log)
        process2 = subprocess.Popen(command2, stdin=process1.stdout, stderr=log)
        process1.stdout.close()
        process2.communicate()

bam = args.bam

process_file(bam)