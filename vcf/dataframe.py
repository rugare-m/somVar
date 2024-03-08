import pysam
import pandas as pd
from Bio import SeqIO
import pyfastx
from cyvcf2 import VCF
import subprocess
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description="VCF Processing Pipeline")
    parser.add_argument('-f', '--reference_genome_path', required=True, help="Path to the reference genome in fasta format")
    parser.add_argument('-v1', '--vcf1', required=True, help="Path to the vcf 1")
    parser.add_argument('-v2', '--vcf2', required=True, help="Path to the vcf 2")
    parser.add_argument('-v3', '--vcf3', required=True, help="Path to the vcf 3 ")
    parser.add_argument('-v4', '--vcf4', required=True, help="Path to the vcf 4 ")
    parser.add_argument('-vx', '--merged', required=True, help="Path to the merged vcf file")
    return parser.parse_args()

# Function to extract the variants from a VCF file in format of chrom:pos:ref>alt
def extract_variants_from_vcf(vcf_path):
    
    #Open VCF
    vcf_file = pysam.VariantFile(vcf_path)

    # Initialize a set to store the variants
    variants_set = set()

    # Iterate through the variants in the VCF file
    for record in vcf_file:
        chrom = record.chrom
        pos = record.pos
        ref = record.ref
        alts = record.alts

        # Create a string in the format of chrom:pos:ref>alt
        variant = f"{record.chrom}:{record.pos}:{record.ref}>{record.alts[0]}"
        # Add the variant to the set
        variants_set.add(variant)
        


    vcf_file.close()
    #Return the set of variants
    return variants_set



# Path to VCFs
args = parse_arguments()
caller1 = args.bcftools
caller2 = args.freebayes
caller3 = args.lofreq
caller4 = args.mutect2

merged = args.merged






vcf_files = [caller1, caller2, caller3, caller4, merged]

for vcf_file in vcf_files:
    subprocess.run(["bcftools", "index", "-f", vcf_file])
    
    
    
# get the set of variants for each VCF file
print("Extracting variants from VCF files...", flush=True)
merge_set = extract_variants_from_vcf(merged)
c1_set = extract_variants_from_vcf(caller1)
c2_set = extract_variants_from_vcf(caller2)
c3_set = extract_variants_from_vcf(caller3)
c4_set = extract_variants_from_vcf(caller4)
print("Done extracting variants.", flush=True)

#dictionary to store the variants and their presence in each VCF file
variants_dict = {}

# Iterate through the variants in merged VCF
for variant in merge_set:
    variants_dict[variant] = {
        "Caller 1":int(variant in c1_set),
        "Caller 2":int(variant in c2_set),
        "Caller 3":int(variant in c3_set),
        "Caller 4":int(variant in c4_set),
    }
        
print("Done creating dictionary.", flush=True)


#Parse the VCF file and extract the features - input vcf is merged vcf
def parse_vcf(vcf_file):
    variant_info = {}

    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            variant = f"{record.chrom}:{record.pos}:{record.ref}>{record.alts[0]}"
            
            # alt median base quality
            if 'MBQ' in record.info:
                mbq = record.info['MBQ'][1]  # Extract the value from the tuple
            else:
                mbq = "NA"
            
            
            rsr = record.samples[0]['AD'][0]
            rsa = record.samples[0]['AD'][1]
            # depth
            dp = rsa + rsr
                
            # Ref minus Alt fragment length    
            if 'MFRL' in record.info:
                ref_fragment = record.info['MFRL'][0]
                alt_fragment = record.info['MFRL'][1]
                delta = alt_fragment - ref_fragment
            else:
                delta = "NA"
            
            #strand bias
            sb = record.info.get('FS', None)
            #mapping quality
            mapq = record.info.get('MMQ', None)
            mapq = mapq[1] if mapq else None  #unpack tuple > get second alt value
            #distance from end of read
            mpos = record.info.get('MPOS', None) if 'MPOS' in record.info else None
            mpos = mpos[0] if mpos else None  #unpack tuple if it exists
            #allele frequency
            af = rsa/dp if dp else None
            #dbSNP
            dbsnp =  1 if record.id and 'rs' in record.id else 0
            #COSMIC
            cosmic = 1 if record.id and 'COSV' in record.id else 0

            variant_info[variant] = {"Median ALT base quality": mbq,
                                     "Reads Supporting Reference": rsr,
                                     "Reads Supporting Alternate": rsa,
                                     "Read Depth": dp,
                                     "Median Alt Fragment Length": alt_fragment,
                                     "Strand Bias": sb,
                                     "Mapping Quality": mapq,
                                     "Distance From End of Read": mpos,
                                     "Allele Frequency": af,
                                     "dbSNP": dbsnp,
                                     "COSMIC": cosmic,
                                    }

    return variant_info

#extract features from merged vcf file
vcf_file_path = merged
result = parse_vcf(vcf_file_path)

# Merge the two dictionaries - results contains features and variants_dict contains presence in each VCF
merged_dict = {}
for key in result:
    if key in variants_dict:
        merged_dict[key] = {**result[key], **variants_dict[key]}

#print (merged_dict)

# Print the len of the merged dictionary - confirm that the number of variants is the same
#print ("length of merged dict: ",len(merged_dict))
#print ("length of features dict: ",len(result))
#print ("length of concordance dict: ",len(variants_dict))


def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    total_bases = len(sequence)
    
    if total_bases == 0:
        return 0.0
    gc_content = (gc_count / total_bases) * 100
    return gc_content

def find_homopolymers(sequence):
    homopolymers = []
    current_homopolymer_length = 1

    # iterate through the sequence starting from the second nucleotide
    for i in range(1, len(sequence)):
        # check if the current nucleotide is the same as the previous one
        if sequence[i] == sequence[i - 1]:
            current_homopolymer_length += 1
        else:
            # if a different nucleotide is encountered, check the length of the homopolymer
            if current_homopolymer_length >= 4:
                homopolymers.append((sequence[i - 1], i - current_homopolymer_length, current_homopolymer_length))
            # reset the homopolymer length for the new nucleotide
            current_homopolymer_length = 1

    # check for homopolymers at the end of the sequence
    if current_homopolymer_length >= 4:
        homopolymers.append((sequence[-1], len(sequence) - current_homopolymer_length, current_homopolymer_length))

    return homopolymers

def get_flanking_sequences(fasta, sequence_id, position, flank_length):
    left_flank = ""
    right_flank = ""

    sequence = fasta[sequence_id]

    if position >= flank_length and position + flank_length < len(sequence):
        lpos = position - 1
        rpos = position
        left_flank = str(sequence[lpos - flank_length:lpos]).upper()
        right_flank = str(sequence[rpos:rpos + flank_length]).upper()

    return left_flank, right_flank

# Provide file paths directly
vcf_file_path = merged

args = parse_arguments()
reference_genome_path = args.reference_genome_path
flank_length = 20
snp = 'REF'  

hg38 = pyfastx.Fasta(reference_genome_path)

variant_cont = {}
for variant in VCF(vcf_file_path):
    print(f"\nProcessing variant: {variant.CHROM}:{variant.POS}:{variant.REF}>{variant.ALT[0]}")
    flank_length = int(flank_length)
    gapp = variant.POS + 1
    gapn = variant.POS - 1
    intervalp = (gapp, variant.POS + flank_length)
    intervaln = (variant.POS - flank_length, variant.POS)
    fasta = []

    left_flank, right_flank = get_flanking_sequences(
        hg38, variant.CHROM, variant.POS, flank_length
    )

    seq = left_flank + variant.REF + right_flank

    # GC content
    gc = calculate_gc_content(seq)
    print(f"GC Content: {gc}")

    # WHR from > A machine learning model to determine the accuracy of variant calls in capture based next generation sequencing
    homopolymers = find_homopolymers(seq)
    sum = 0
    total = 0

    for nucleotide, start_index, length in homopolymers:
        sum += (length ** 2)
        total += 1

    if sum == 0:
        whr = 0
    else:
        whr = (sum / total)

    print(f"WHR: {whr}")

    variant_key = f"{variant.CHROM}:{variant.POS}:{variant.REF}>{variant.ALT[0]}"
    variant_cont[variant_key] = {"WHR": whr, "GC %": gc}

len_variants = len(variant_cont)
print(f"\nTotal variants processed: {len_variants}")


def merge_dictionaries(variant_cont, merged_dict):
    dfx = {}

    # Iterate over keys in both dictionaries
    for key in set(variant_cont.keys()) | set(merged_dict.keys()):
        # Merge inner dictionaries if the key is common
        if key in variant_cont and key in merged_dict:
            dfx[key] = {**variant_cont[key], **merged_dict[key]}
        elif key in variant_cont:
            dfx[key] = variant_cont[key]
        elif key in merged_dict:
            dfx[key] = merged_dict[key]

    return dfx

result = merge_dictionaries(variant_cont, merged_dict)

# Convert the dictionary to a Pandas DataFrame and save it to a CSV file
df = pd.DataFrame.from_dict(result, orient='index')
#df.to_csv('all_feature_data.csv', index_label='Variant')
df.to_csv('dataframe.csv', index_label='Variant')

