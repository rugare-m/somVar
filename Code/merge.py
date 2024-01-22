import subprocess
import os
import glob

dirs = ["SRR3401417/", "SRR3401416/", "SRR3401415/", "SRR3401407/", "SRR3401418/"]
dirs = ["fastq/"]
for dir in dirs:
    directory_path = f"/users/rugarem/volatile/chapter3/data/erve/plasma/{dir}" 
    vcf_files = glob.glob(os.path.join(directory_path, "COSMIC_*.vcf.gz"))
    #vcf_files = glob.glob(os.path.join(directory_path, "complete_*.vcf.gz"))
    for vcf_file in vcf_files:
        print(f"Processing file: {vcf_file}", flush=True)

        tabix = ["tabix", "-f", vcf_file]
        subprocess.run(tabix)

        print("File indexed", flush=True)

        sort_out = vcf_file.replace(".vcf.gz", ".sort.vcf.gz")
        sort = ["bcftools", "sort", "-O", "z", "-o", sort_out, vcf_file]
        subprocess.run(sort)

        print("File sorted", flush=True)

        norm_out = vcf_file.replace(".vcf.gz", ".norm.vcf.gz")
        norm = ["bcftools", "norm", "-O","z", "-m", "-both", "-o", norm_out, sort_out]
        subprocess.run(norm)

        print("File normalized", flush=True)
    
        index = ["gatk", "IndexFeatureFile", "-I", norm_out]
        #tabix = ["tabix", "-p","-f", "vcf", norm_out]
        subprocess.run(index)
        #subprocess.run(tabix)

        print("File indexed and tabixed", flush=True)

        sample_id = vcf_file.split('.')[-3]


        annotate_out = vcf_file.replace(".vcf.gz", ".annotated.vcf.gz")

        #index BAM if not already indexed
        if directory_path+f"{sample_id}.bqsr.bam.bai" in os.listdir(directory_path):
            print("BAM file already indexed")
        else:
            index = ["samtools", "index", directory_path+f"{sample_id}.bqsr.bam"]
            subprocess.run(index)
            

        if directory_path+f"{sample_id}.bqsr.bam.bai" in os.listdir(directory_path):
            print("BAM file already indexed")
        else:
            subprocess.run(index)

        #annotate vcf 
        annotator = [
               "gatk", 
               "VariantAnnotator", 
               "-R", "/users/rugarem/volatile/chapter3/data/reference/hg38.fa.gz", 
               "-I", directory_path+f"{sample_id}.bqsr.bam", 
               "-V",norm_out, 
               "-O", annotate_out, 
               "-A", "FisherStrand",
               "-A", "MappingQuality",
               "-A", "BaseQuality", 
               "-A", "FragmentLength", 
               "-A", "ReadPosition",
               "-A", "AlleleFraction",
               "-A", "DepthPerAlleleBySample"
               ]

        subprocess.run(annotator)

        index = ["bcftools", "index", annotate_out]
        tabix = ["tabix", "-f", annotate_out]
        subprocess.run(index)
        subprocess.run(tabix)
    
    
    
    
    
    print("Starting merging process...")
    vcf_files = glob.glob(os.path.join(directory_path, "*.annotated.vcf.gz"))

    sample_id_groups = {}
    for vcf_file in vcf_files:
        sample_id = vcf_file.split('.')[-4]
        if sample_id not in sample_id_groups:
            sample_id_groups[sample_id] = []
        sample_id_groups[sample_id].append(vcf_file)
    
    
    for sample_id, vcf_files in sample_id_groups.items():
        print(f"Processing files for Sample ID: {sample_id}")

        # Write the file names to a text file
        output_file = directory_path+f"{sample_id}_vcfs.list"
        with open(output_file, 'w') as file:
            for vcf_file in vcf_files:
                file.write(f"{vcf_file}\n")

        gatk_merge = [
            'gatk', 'MergeVcfs',
            '-R', '/users/rugarem/volatile/chapter3/data/reference/hg38.fa.gz',
            '-I', directory_path+f"{sample_id}_vcfs.list",
            '-O', directory_path+f"{sample_id}.merged.vcf.gz"
        ]

        subprocess.run(gatk_merge)

        print(f"Completed merging files for Sample ID: {sample_id}", flush=True)


        snps = [
            "bcftools", "view", "-v", "snps", "-O", "z",
            "-o", f"{directory_path}{sample_id}.merged.snps.vcf.gz",
            f"{directory_path}{sample_id}.merged.vcf.gz"
        ]
        print ("Starting SNP extraction -  command: ",snps, flush=True)
        subprocess.run(snps)
        
        dups = ["bcftools", "norm", "-d", "all", "-o", 
                f"{directory_path}{sample_id}.dedup.snps.vcf.gz", 
                f"{directory_path}{sample_id}.merged.snps.vcf.gz"]
        
        print ("Removing duplicate SNPs -  command: ",dups, flush=True)
        subprocess.run(dups)


        # Clean up
        rm_patterns = [directory_path+"*tbi", 
                       directory_path+"*sort*", 
                       directory_path+"*norm*",
                       directory_path+"*merged*"]
        rm_command = "rm -f {}".format(" ".join(rm_patterns))

        subprocess.run(rm_command, shell=True)
        print("Cleanup complete", flush=True)

    #delete  vcf files
    rm_patterns = [directory_path+"*.merged.vcf.gz"]
    rm_command = "rm -f {}".format(" ".join(rm_patterns))
    subprocess.run(rm_command, shell=True)
    _ = [print(f"Deleting file: {file}") or os.remove(file) for file in vcf_files]

