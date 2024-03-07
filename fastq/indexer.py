import os
import argparse
import subprocess

def index_vcf_files(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".vcf.gz"):
                file_path = os.path.join(root, file)
                print("Indexing " + file_path, flush=True)
                subprocess.run(["tabix", "-p", "vcf", file_path])

def main():

    parser = argparse.ArgumentParser(description="Index VCF files in a specified directory.")

    parser.add_argument("-d","--directory", required=True, help="The directory containing VCF files.")
    args = parser.parse_args()
    index_vcf_files(args.directory)

if __name__ == "__main__":
    main()
