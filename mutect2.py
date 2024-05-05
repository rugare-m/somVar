import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Mutect2 Variant Calling Pipeline")
    parser.add_argument('-R', '--reference_file', required=True, help="Reference genome in fasta format")
    parser.add_argument('-b', '--bam', required=True, help="Single BAM file")

    return parser.parse_args()

args = parse_arguments()

reference_file = args.reference_file
bam = args.bam

def process_file(alignments_file):
    sample_name = os.path.basename(alignments_file).split(".")[0]
    output_directory = os.path.dirname(alignments_file)
    output_file = os.path.join(output_directory, f"mutect2.{sample_name}.vcf")
    #filtered = os.path.join(output_directory, f"filtered_mutect2.{sample_name}.vcf")
    log_file = os.path.join(output_directory, f"mutect2.{sample_name}.log") 

    mutect2_command = ["gatk", "Mutect2", "-R", reference_file, "-I", alignments_file, "-O", output_file]
    compress = ["bgzip", output_file]
    #mutect2_filter = ["gatk", "FilterMutectCalls", "-R", reference_file, "-V", output_file, "-O", filtered]
    
    with open(log_file, "w") as log:
        process = subprocess.run(mutect2_command, stdout=subprocess.PIPE, stderr=log)
        #process2 = subprocess.run(mutect2_filter, stdout=subprocess.PIPE, stderr=log)

    subprocess.run(compress, stdout=subprocess.PIPE)


process_file(bam)


