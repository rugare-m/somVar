import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="BCFTools Variant Calling Pipeline")
    parser.add_argument('-f', '--reference_file', required=True, help="Reference genome in fasta format")
    parser.add_argument('-d', '--dirs', nargs='+', required=True, help="List of input directories containing BAM files")

    return parser.parse_args()

args = parse_arguments()

reference_file = args.reference_file
dirs = args.dirs

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

dirs = args.dirs

for dir in dirs:
    bam_directory = dir
    bam_files = [os.path.join(bam_directory, f) for f in os.listdir(bam_directory) if f.endswith(".bqsr.bam")]
    
    for bam_file in bam_files:
        subprocess.run(["samtools", "index", bam_file])

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

            