#!/usr/bin/env python
# coding: utf-8

import os
import glob
import pandas as pd
import time

# Path to your featureCounts output folder
path = "/Users/aakash/Downloads/rnaseq/quants"

# Find all featureCount files
files = glob.glob(os.path.join(path, "*_featurecounts.txt"))
print("Files found:", files)

all_counts = []

for file in files:
    start_time = time.time()
    
    # Read file (skip comment lines starting with #)
    df = pd.read_csv(file, sep="\t", comment="#")
    
    # Extract sample name from filename
    sample_name = os.path.basename(file).replace("_featurecounts.txt", "")
    
    # Keep only Geneid and last column (counts)
    df = df[["Geneid", df.columns[-1]]]
    df.rename(columns={df.columns[-1]: sample_name}, inplace=True)
    
    all_counts.append(df)
    
    elapsed = (time.time() - start_time) / 60  # minutes
    print(f"Completed {sample_name} | Rows: {df.shape[0]} | Time: {elapsed:.2f} min")

# Merge all dataframes
counts_matrix = all_counts[0]
for df in all_counts[1:]:
    counts_matrix = counts_matrix.merge(df, on="Geneid", how="outer")

# Save output
output_file = os.path.join(path, "counts_matrix.csv")
counts_matrix.to_csv(output_file, index=False)

print("\nAll files processed!")
print("Merged matrix shape:", counts_matrix.shape)
print("Saved to:", output_file)

