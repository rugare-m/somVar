import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="FreeBayes Variant Calling Pipeline")
    parser.add_argument('-f', '--reference_file', required=True, help="Reference genome in fasta format")
    parser.add_argument('-d', '--dirs', nargs='+', required=True, help="List of input directories containing BAM files")

    return parser.parse_args()

args = parse_arguments()

reference_file = args.reference_file
dirs = args.dirs

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
dirs = args.dirs

for dir in dirs:
    bam_directory = dir
    bam_files = [os.path.join(bam_directory, f) for f in os.listdir(bam_directory) if f.endswith(".bqsr.bam")]
    
    with ProcessPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(process_file, bam_file) for bam_file in bam_files]
        for future in as_completed(futures):
            pass

    logs_directory = "logs"
    if not os.path.exists(logs_directory):
        os.makedirs(logs_directory)

    for log_file in os.listdir(bam_directory):
        if log_file.endswith(".log"):
            source_path = os.path.join(bam_directory, log_file)
            dest_path = os.path.join(logs_directory, log_file)
            os.rename(source_path, dest_path)



