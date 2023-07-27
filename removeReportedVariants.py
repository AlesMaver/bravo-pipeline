###########################
# Modify the SGP VCF file to remove reported genotypes
###########################
import gzip
import time
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

chunk_size = 10000
limit = 100000000
header_row = ""
empty_genotype_string = "./.:.:.:.:."
variants_to_remove_hg38 = df["variant_unique_id_hg38"].tolist()
start_time = time.time()

with gzip.open(input_file, 'rt') as input_f, gzip.open(output_file, 'wt') as output_f:
    # for _ in range(500000):
    #     _discard = input_f.readline()  # Read and discard the line

    processed_rows = 0
    chunked_lines = []

    for line in input_f:
        current_variant = ""
        row = line.strip().split('\t')
        
        # Get the header row with sample names
        if line.startswith('#CHROM'):
                header_row = row

        if len(row)>1:
            current_variant = "-".join([row[0], row[1], row[3], row[4]])

            if current_variant in variants_to_remove_hg38:
                print("FOUND A MATCH - will remove the genotypes!")
                matched_samples = df[df["variant_unique_id_hg38"]==current_variant]["patient_identifier"].tolist()
                print(matched_samples)
                genotypes_to_clear = [index for index, item in enumerate(header_row) if item in matched_samples]
                if len(genotypes_to_clear)>0:
                    for genotype_index in genotypes_to_clear:
                        print(row[genotype_index])
                        row[genotype_index] = empty_genotype_string

        output_line = '\t'.join(row) + '\n'
        chunked_lines.append(output_line)
        #_ = output_f.write(output_line) # Redirect output to _ to prevent printing to console

        processed_rows += 1
        if processed_rows % chunk_size == 0:
            print("Working on row: " + str(processed_rows))
            
            # Write data in chunks
            print("Writing data for chunk...")
            _ = output_f.writelines(chunked_lines) # Redirect output to _ to prevent printing to console
            chunked_lines = []

            print("Continue searching for matches...")

            # Track the timing per chunk
            end_time = time.time()
            elapsed_time = end_time - start_time
            anticipated_time = ((33000000-processed_rows)/chunk_size * elapsed_time)/3600
            anticipated_time = ((34000000-processed_rows) * (elapsed_time/processed_rows))/3600
            print(f"Elapsed time for {processed_rows} records: {elapsed_time:.4f} seconds")
            print(f"Elapsed time for {chunk_size} records: {anticipated_time:.1f} hours")
            #start_time = time.time()

        if processed_rows >= limit:
            break
    
    # Write the final chunk
    print("Final chunk...")
    _ = output_f.writelines(chunked_lines) # Redirect output to _ to prevent printing to console
    chunked_lines = []
