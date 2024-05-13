import os
import subprocess

def index_vcf_files():
    # Get the directory path where the script is located
    script_directory = os.path.dirname(__file__)
    # Set the directory path to the parent directory of the script's directory
    directory = os.path.abspath(os.path.join(script_directory, os.pardir))

    # Append 'databases' to the directory path
    directory = os.path.join(directory, 'databases')

    # Get a list of all VCF files in the directory
    vcf_files = [file for file in os.listdir(directory) if file.endswith('.vcf.gz')]

    # Iterate over each VCF file and index it using tabix
    for vcf_file in vcf_files:
        vcf_path = os.path.join(directory, vcf_file)
        subprocess.run(['tabix','-f', '-p', 'vcf', vcf_path], check=True)
        print(f"Indexing completed for {vcf_file}")


# Example usage
index_vcf_files()
