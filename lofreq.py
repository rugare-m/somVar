import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="LoFreq Variant Calling Pipeline")
    parser.add_argument('-f', '--reference_file', required=True, help="Reference genome in fasta format")
    parser.add_argument('-b', '--bam', required=True, help="Single BAM file")

    return parser.parse_args()

args = parse_arguments()

reference_file = args.reference_file
bam = args.bam

def process_file(alignments_file):
    subprocess.run(["samtools", "index", bam])
    sample_name = os.path.basename(alignments_file).split(".")[0]
    output_directory = os.path.dirname(alignments_file)
    output_file = os.path.join(output_directory, f"lofreq.{sample_name}.vcf.gz")
    with_gt = os.path.join(output_directory, f"lofreq_gt.{sample_name}.vcf.gz")
    log_file = os.path.join(output_directory, f"lofreq.{sample_name}.log")  # Log file name for caller output

    lofreq_command = ["lofreq", "call-parallel","--pp-threads", "8", "-f", reference_file, "-o", output_file, alignments_file]

    with open(log_file, "w") as log:
        process = subprocess.Popen(lofreq_command, stdout=subprocess.PIPE, stderr=log)
        process.communicate()
    
    gt = ["lofreq_gt.sh", "-i", output_file, "-g", 
      "1/1", "-n", "20", "-o", with_gt]
    #subprocess.run(gt)

    #clean up
    #os.remove(output_file)


process_file(bam)





