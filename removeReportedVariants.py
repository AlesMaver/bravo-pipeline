###########################
# Modify the SGP VCF file to remove reported genotypes
###########################
from Bio import bgzf
import argparse
import pandas as pd

argParser = argparse.ArgumentParser()
argParser.add_argument("-i", "--input_file", help="Input VCF file")
argParser.add_argument("-o", "--output_file", help="Output VCF file")
argParser.add_argument("-v", "--reported_variants", help="Reported variants file")

args = argParser.parse_args()

input_file = args.input_file
output_file = args.output_file
df = pd.read_csv(args.reported_variants)

header_row = ""
empty_genotype_string = "./.:.:.:.:."
variants_to_remove_hg38 = df["variant_unique_id_hg38"].tolist()

with bgzf.open(input_file, 'rt') as input_f, bgzf.open(output_file, 'wt') as output_f:
    # for _ in range(500000):
    #     _discard = input_f.readline()  # Read and discard the line

    for line in input_f:
        current_variant = ""
        row = line.strip().split('\t')
        
        # Get the header row with sample names
        if line.startswith('#CHROM'):
            header_row = row

        if len(row)>1:
            current_variant = "-".join([row[0], row[1], row[3], row[4]])

            if current_variant in variants_to_remove_hg38:
                # print("FOUND A MATCH - will remove the genotypes!")
                matched_samples = df[df["variant_unique_id_hg38"]==current_variant]["patient_identifier"].tolist()
                # print(matched_samples)
                genotypes_to_clear = [index for index, item in enumerate(header_row) if item in matched_samples]
                if len(genotypes_to_clear)>0:
                    for genotype_index in genotypes_to_clear:
                        # print(row[genotype_index])
                        row[genotype_index] = empty_genotype_string

        _ = output_f.write('\t'.join(row) + '\n') # Redirect output to _ to prevent printing to console
