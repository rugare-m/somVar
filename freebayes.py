import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="FreeBayes Variant Calling Pipeline")
    parser.add_argument('-f', '--reference_file', required=True, help="Reference genome in fasta format")
    parser.add_argument('-b', '--bam', required=True, help="Single BAM file")

    return parser.parse_args()

args = parse_arguments()

reference_file = args.reference_file

def process_file(alignments_file):
    sample_name = os.path.basename(alignments_file).split(".")[0]
    output_directory = os.path.dirname(alignments_file)
    output_file = os.path.join(output_directory, f"fbayes.{sample_name}.vcf")
    log_file = os.path.join(output_directory, f"fbayes.{sample_name}.log")  # Log file name for caller output

    fbayes_command = ["freebayes", "-f", reference_file, alignments_file]

    with open(log_file, "w") as log:
        process = subprocess.Popen(fbayes_command, stdout=open(output_file, "w"), stderr=log)
        process.communicate()

    process2 = subprocess.Popen(["bgzip", output_file], stdout=subprocess.PIPE)
    process2.communicate() 

bam = args.bam

process_file(bam)

