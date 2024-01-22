#works
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

def process_file(alignments_file):
    sample_name = os.path.basename(alignments_file).split(".")[0]
    output_directory = os.path.dirname(alignments_file)
    reference_file = "/users/rugarem/volatile/chapter3/data/reference/uhg38.fa"
    #reference_file = "/users/rugarem/volatile/chapter3/data/reference/hg38.fa.gz"
    output_file = os.path.join(output_directory, f"fbayes.{sample_name}.vcf")
    log_file = os.path.join(output_directory, f"fbayes.{sample_name}.log")  # Log file name for caller output

    fbayes_command = ["freebayes", "-f", reference_file, alignments_file]

    with open(log_file, "w") as log:
        process = subprocess.Popen(fbayes_command, stdout=open(output_file, "w"), stderr=log)
        process.communicate()
    
    process2 = subprocess.Popen(["bgzip", output_file], stdout=subprocess.PIPE)

dirs = ["/users/rugarem/volatile/chapter3/data/erve/tissue/fastq/"]
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

    for log_file in os.listdir(dir):
        if log_file.endswith(".log"):
            source_path = os.path.join(dir, log_file)
            dest_path = os.path.join(logs_directory, log_file)
            os.rename(source_path, dest_path)
            
        