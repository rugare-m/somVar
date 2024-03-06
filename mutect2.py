import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Mutect2 Variant Calling Pipeline")
    parser.add_argument('-R', '--reference_file', required=True, help="Reference genome in fasta format")
    parser.add_argument('-d', '--dirs', nargs='+', required=True, help="List of input directories containing BAM files")

    return parser.parse_args()

args = parse_arguments()

reference_file = args.reference_file
dirs = args.dirs

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

