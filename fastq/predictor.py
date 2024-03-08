import joblib
import pandas as pd
import gzip
import subprocess
from sklearn.preprocessing import StandardScaler
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter VCF file based on predicted variants")
    parser.add_argument('-o', '--output_vcf', required=True, help="Output VCF file path (.gz)")
    parser.add_argument('-x', '--threshold', required=True, help="Custom threshold for the model")
    parser.add_argument('-m', '--model', choices=['low_depth_model.pkl', 'high_depth_model.pkl'], required=True, help="Select the model (low or high depth)")
    parser.add_argument('-d', '--merged', required=True, help="Path to merged vcf file")
    return parser.parse_args()

print("Loading data...")
dataframe = pd.read_csv('dataframe.csv')

# Keep SNPs only
print("Filtering SNPs...", flush=True)
snp_pattern = r'^chr\d+:\d+:[ACGT]>[ACGT]$'
snp_mask = dataframe['Variant'].str.match(snp_pattern)
snp_val_df = dataframe[snp_mask]
snp_val_df.reset_index(drop=True, inplace=True)

# Drop rows with missing values and 'Read Depth' 1 or lower
print("Cleaning data...")
snp_val_df = snp_val_df.dropna()
snp_val_df = snp_val_df[snp_val_df['Read Depth'] > 1]

# Scale after dropping READ DEPTH 1 or lower
print("Scaling data...", flush=True)
column_to_scale = 'Read Depth'
scaler = StandardScaler()
snp_val_df.loc[:, [column_to_scale]] = scaler.fit_transform(snp_val_df[[column_to_scale]].values.reshape(-1, 1))

# Split the data into features and targets
X_val = snp_val_df.drop(['Variant','Reads Supporting Reference','Reads Supporting Alternate'], axis=1)

# Load the trained model and custom threshold from the file
args = parse_arguments()
print("Loading trained model and threshold...", flush=True)
model_filename = args.model
loaded_data = joblib.load(model_filename)
loaded_model = loaded_data['model']
threshold = args.threshold
custom_threshold = float(threshold)

# Make predictions on the new data using the custom threshold
print("Making predictions...", flush=True)
new_data = X_val
predictions = (loaded_model.predict_proba(new_data)[:, 1] > custom_threshold).astype(int)
snp_val_df['Predicted'] = predictions

# Filter df based on the 'Predicted' column
predicted_variants = snp_val_df[snp_val_df['Predicted'] == 1]
variant_identifiers = predicted_variants['Variant'].str.split(':', expand=True).iloc[:, :4].apply(':'.join, axis=1).tolist()

if __name__ == "__main__":
    args = parse_arguments()
    output_vcf = args.output_vcf
    merged = args.merged

    print("Filtering VCF...")
    # Filter VCF lines based on variant identifiers
    def filter_vcf(vcf_path, output_vcf, variant_identifiers):
        with gzip.open(vcf_path, 'rt') as vcf_file, gzip.open(output_vcf, 'wt') as output_file:
            for line in vcf_file:
                if line.startswith('#'):
                    # Write header lines to the output VCF file
                    output_file.write(line)
                else:
                    # Extract variant identifier from the VCF line
                    vcf_fields = line.split('\t')
                    vcf_variant_identifier = '{}:{}:{}>{}'.format(vcf_fields[0], vcf_fields[1], vcf_fields[3], vcf_fields[4])
                    # Check if the variant identifier is in the list of predicted variants
                    if vcf_variant_identifier in variant_identifiers:
                        output_file.write(line)

    filter_vcf(merged, output_vcf, variant_identifiers)

    # Remove duplicate variants 
    #subprocess.run(['bcftools', 'norm', '-d', 'all', '-o', output_vcf.replace('.vcf.gz', 'hc.vcf.gz'), output_vcf])
    subprocess.run(['bcftools', 'norm', '-d', 'all', '-o', output_vcf.replace('.dedup.snps.vcf.gz', 'hc.vcf.gz'), output_vcf])
    # Delete the original vcf file
    subprocess.run(['rm', output_vcf])
    print("Process completed successfully.", flush=True)
#end
