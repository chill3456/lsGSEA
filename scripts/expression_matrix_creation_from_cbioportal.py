#!python

import os
import sys
import pandas as pd

def transform_gene_counts_file(file_path):
    # Open the input file for reading
    with open(file_path, "r") as input_fh:
        # Read the content of the file
        lines = input_fh.readlines()

    # Skip the first, third, fourth, fifth, and sixth rows
    lines = lines[1:2] + lines[6:]

    # Create a new list to store modified rows
    modified_rows = []

    # Process each line
    for line in lines:
        # Split the line by tab delimiter
        columns = line.strip().split("\t")

        # Keep only the required columns (0, 1, 2, 6)
        modified_row = "\t".join(columns[0:3] + [columns[6]])

        # Add the modified row to the list
        modified_rows.append(modified_row)

    # Construct the output file path
    output_file = file_path.replace(".rna_seq.augmented_star_gene_counts.tsv", "_modified.tsv")

    # Write the modified content to the output file
    with open(output_file, "w") as output_fh:
        output_fh.write("\n".join(modified_rows))

    print(f"File {file_path} transformed successfully. Output written to {output_file}.")
    return output_file

def transform_gene_counts_file_in_directory(root_directory):
    # Walk through all subdirectories and files in the root directory
    for dirpath, dirnames, filenames in os.walk(root_directory):
        # Process each file in the current directory
        for filename in filenames:
            if filename.endswith(".rna_seq.augmented_star_gene_counts.tsv"):
                file_path = os.path.join(dirpath, filename)
                transform_gene_counts_file(file_path)

def combine_modified_files(root_directory):
    # Create an empty DataFrame to store the combined data
    combined_data = pd.DataFrame()

    # Track if the columns have been added to the combined_data DataFrame
    columns_added = False

    # Walk through all subdirectories and files in the root directory
    for dirpath, dirnames, filenames in os.walk(root_directory):
        # Process each file in the current directory
        for filename in filenames:
            if filename.endswith("_modified.tsv"):
                file_path = os.path.join(dirpath, filename)
                df = pd.read_csv(file_path, sep="\t")

                # Get the file name without the extension and "_modified" suffix
                file_name = os.path.splitext(filename)[0].replace("_modified", "")

                # If columns have not been added yet, add only the desired columns
                if not columns_added:
                    combined_data = df[["gene_id", "gene_name", "gene_type"]]
                    columns_added = True

                # Rename the "tpm_unstranded" column to the modified file name
                df = df.rename(columns={"tpm_unstranded": file_name})

                # Merge the current file's data with the combined data based on "gene_id"
                combined_data = pd.merge(combined_data, df[["gene_id", file_name]], on="gene_id", how="outer")

    # Save the combined data to a file
    output_file = os.path.join(root_directory, "combined_gene_counts.tsv")
    combined_data.to_csv(output_file, sep="\t", index=False)

    print(f"Combined file created successfully. Output written to {output_file}.")



if __name__ == "__main__":
    # Check if the root directory path argument is provided
    if len(sys.argv) < 2:
        print("Please provide the root directory path as an argument.")
        sys.exit(1)

    # Get the root directory path from the command line argument
    root_directory = sys.argv[1]

    # Combine modified files in the root directory
    combine_modified_files(root_directory)
