"""
RNA Analysis: AL Bulges Separation
==================================

Description:
------------
This script analyzes RNA bulge separation data stored in files with two columns:
1. Distance (integer)
2. Sequence (string)

It performs the following:
1. Reads all files from the specified directory (`results/stemloops/AL_bulges_separation`).
2. Parses filenames to extract AL length (XXX) and bulge size (YYY).
3. Computes:
   - Occurrence (number of lines in each file).
   - Mean distance and standard deviation (from the 1st column).
   - Aggregated sequences.
4. Generates:
   - Heatmaps for occurrence (absolute and relative frequencies).
   - Heatmap for mean distances.
   - RNA conservation logos for all sequences.

Input:
------
Text files named in the format: `AL_XXXnt_YYYbulge.txt`
  - XXX: Length of AL
  - YYY: Bulge size

Output:
-------
1. Heatmaps (`.png`): Occurrence (absolute and relative) and mean distances.
2. Sequence logos (`.png`): For all sequences.
3. Statistics (`_stats` text files): For occurrences, mean distances, and standard deviations.

Dependencies:
-------------
- Python libraries: `os`, `pandas`, `logomaker`, `matplotlib`, `seaborn`, `numpy`.

Author: [Your Name]
Date: [Current Date]
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict
import logomaker

# Function to parse filenames and extract metadata
def parse_filename(filename):
    """
    Extract AL length (XXX) and bulge size (YYY) from filename.
    Args:
        filename (str): Filename in the format AL_XXXnt_YYYbulge.txt.
    Returns:
        tuple: (XXX, YYY) as integers.
    """
    filename = filename.replace("AL_", "").replace("bulge.txt", "")
    XXX, YYY = filename.split("nt_")
    return int(XXX), int(YYY)

# Function to process all files in the directory
def process_bulge_files(folder):
    """
    Process files in the specified folder to extract data for analysis.
    Args:
        folder (str): Path to the folder containing RNA data files.
    Returns:
        dict: Processed data containing occurrences, distances, and sequences.
    """
    data = defaultdict(lambda: defaultdict(list))

    for dirpath, _, filenames in os.walk(folder):
        for filename in filenames:
            if not filename.endswith(".txt"):
                continue

            filepath = os.path.join(dirpath, filename)
            print(f"Processing file: {filepath}")

            # Parse filename to get XXX and YYY
            XXX, YYY = parse_filename(filename)

            # Read file and populate data
            with open(filepath, "r") as f:
                for line in f:
                    fields = line.strip().split("\t")
                    if len(fields) != 2:
                        print(f"Skipping malformed line: {line.strip()}")
                        continue

                    distance = float(fields[0])
                    sequence = fields[1]
                    key = f"{YYY}"

                    # Populate data
                    data[XXX][key].append((distance, sequence))

    # Convert data into matrices and collect sequences
    rows = sorted(data.keys())
    cols = sorted({key for subdict in data.values() for key in subdict.keys()})
    occurrences = pd.DataFrame(index=rows, columns=cols).fillna(0)
    mean_distances = pd.DataFrame(index=rows, columns=cols).fillna(0)
    most_frequent_values = pd.DataFrame(index=rows, columns=cols).fillna(0)
    std_distances = pd.DataFrame(index=rows, columns=cols).fillna(0)
    sequences = defaultdict(lambda: defaultdict(list))

    for row in rows:
        for col in cols:
            if col in data[row]:
                distances = [entry[0] for entry in data[row][col]]
                seq_list = [entry[1] for entry in data[row][col]]
                occurrences.at[row, col] = len(distances)

                #report mean distances or closer-to-mean value from the list
                mean_dist = sum(distances) / len(distances) if distances else 0
                closest_to_mean = min(distances, key=lambda x: abs(x - mean_dist)) if distances else None
                mean_distances.at[row, col] = closest_to_mean
                #mean_distances.at[row, col] = sum(distances) / len(distances) if distances else 0
                std_distances.at[row, col] = (pd.Series(distances).std() if len(distances) > 1 else 0)
                
                # Find the value with the highest frequency (mode)
                if distances:
                    modes = pd.Series(distances).mode()
                    most_frequent_value = modes.iloc[0]  # Pick the first mode if multiple
                else:
                    most_frequent_value = None

                most_frequent_values.at[row, col] = most_frequent_value


                sequences[row][col] = seq_list

    return {
        "occurrences": occurrences,
        "mean_distances": mean_distances,
        "std_distances": std_distances,
        "most_freq": most_frequent_values,
        "sequences": sequences,
    }

# (Reuse `sort_data_matrix`, `save_matrices_to_text`, `create_heatmaps`, and `generate_rna_logo` from the original script.)

# Sort the matrix by XXX and YYYxZZZ
def sort_data_matrix(data):
    """
    Sorts data matrices (occurrences, distances) by XXX values (row index).
    Args:
        data (dict): Dictionary containing the data matrices.
    Returns:
        dict: Sorted data matrices.
    """
    # Extract occurrences and distances
    occurrences = data["occurrences"]
    mean_distances = data["mean_distances"]
    std_distances = data["std_distances"]
    most_frequent_values = data["most_freq"]

    # Sort by XXX (row index)
    occurrences = occurrences.sort_index(ascending=False)
    mean_distances = mean_distances.sort_index(ascending=False)
    std_distances = std_distances.sort_index(ascending=False)
    most_frequent_values = most_frequent_values.sort_index(ascending=False)

    # Sort by YYYxZZZ (column index)
    occurrences = occurrences.sort_index(axis=1, ascending=True)
    mean_distances = mean_distances.sort_index(axis=1, ascending=True)
    std_distances = std_distances.sort_index(axis=1, ascending=True)
    most_frequent_values = most_frequent_values.sort_index(axis=1, ascending=True)

    # Return sorted matrices
    return {
        "occurrences": occurrences,
        "mean_distances": mean_distances,
        "std_distances": std_distances,
        "most_freq": most_frequent_values,
    }


# save matrix to text
def save_matrices_to_text(data, output_folder):
    """
    Save sorted matrices to text files.
    Args:
        data (dict): Dictionary containing occurrences, mean distances, and std deviations.
        output_folder (str): Path to save the text files.
    """
    os.makedirs(output_folder, exist_ok=True)

    for matrix_name, matrix_data in data.items():
        filename = os.path.join(output_folder, f"{matrix_name}.txt")
        matrix_data.to_csv(filename, sep="\t", float_format="%.3f")
        print(f"Matrix {matrix_name} saved to {filename}")


# Function to create heatmaps
def create_heatmaps(data, output_folder):
    """
    Generate heatmaps for occurrences, relative occurrences, mean distances, and standard deviations.
    Args:
        data (dict): Dictionary containing occurrences, distances, and sequences.
        output_folder (str): Path to save the heatmap plots.
    """
    os.makedirs(output_folder, exist_ok=True)

    # Sort matrices by XXX and YYYxZZZ
    sorted_data = sort_data_matrix(data)
    occurrences = sorted_data["occurrences"]
    mean_distances = sorted_data["mean_distances"]
    std_distances = sorted_data["std_distances"]
    most_frequent_values = sorted_data["most_freq"]


    # Save matrices as text files
    save_matrices_to_text(sorted_data, os.path.join(output_folder, "text_matrices"))
    
    # Calculate relative occurrences
    total_occurrences = occurrences.sum().sum()
    relative_occurrences = occurrences / total_occurrences

    # Replace zero with NaN
    occurrences.replace(0, float('nan'), inplace=True)
    relative_occurrences.replace(0, float('nan'), inplace=True)
    mean_distances.replace(0, float('nan'), inplace=True)
    std_distances.replace(0, float('nan'), inplace=True)
    most_frequent_values.replace(0, float('nan'), inplace=True)



    # Plot heatmaps
    plt.figure(figsize=(12, 6))
    sns.heatmap(occurrences, annot=True, fmt=".0f", cmap="Blues", cbar_kws={'label': 'Occurrences'}, annot_kws={"rotation": 90})
    plt.title("Absolute Occurrences Heatmap")
    plt.savefig(os.path.join(output_folder, "absolute_occurrences.png"))
    plt.close()

    plt.figure(figsize=(12, 6))
    sns.heatmap(relative_occurrences, annot=True, fmt=".2f", cmap="Purples", cbar_kws={'label': 'Relative Occurrences'}, annot_kws={"rotation": 90})
    plt.title("Relative Occurrences Heatmap")
    plt.savefig(os.path.join(output_folder, "relative_occurrences.png"))
    plt.close()

    plt.figure(figsize=(12, 6))
    sns.heatmap(mean_distances, annot=True, fmt=".0f", cmap="Greens", cbar_kws={'label': 'Mean Distance'}, annot_kws={"rotation": 90})
    plt.title("Mean Distances Heatmap")
    plt.savefig(os.path.join(output_folder, "mean_distances.png"))
    plt.close()

    plt.figure(figsize=(12, 6))
    sns.heatmap(std_distances, annot=True, fmt=".2f", cmap="Oranges", cbar_kws={'label': 'Standard Deviation'}, annot_kws={"rotation": 90})
    plt.title("Standard Deviation of Distances Heatmap")
    plt.savefig(os.path.join(output_folder, "std_distances.png"))
    plt.close()

    plt.figure(figsize=(12, 6))
    sns.heatmap(most_frequent_values, annot=True, fmt=".0f", cmap="Greens", cbar_kws={'label': 'Occurrences'}, annot_kws={"rotation": 90})
    plt.title("Absolute Occurrences Heatmap")
    plt.savefig(os.path.join(output_folder, "most_frequent_distance.png"))
    plt.close()







# Function to generate RNA logo (reusing from the main script)
def generate_rna_logo(sequences, output_file, title="RNA Conservation Logo", reverse=False):
    if reverse:
        sequences = [seq[::-1] for seq in sequences]

    data = pd.DataFrame([list(seq) for seq in sequences])
    counts_df = data.apply(lambda x: x.value_counts(normalize=True)).fillna(0).T

    num_columns = counts_df.shape[0]
    width = max(2, num_columns * 2)
    height = 6

    fig, ax = plt.subplots(figsize=(width, height), dpi=300)
    logo = logomaker.Logo(counts_df, ax=ax, color_scheme="classic")
    logo.style_xticks(rotation=90, anchor=0)
    ax.set_title(title, fontsize=18)
    ax.set_ylabel("Frequency", fontsize=16)
    ax.set_xlabel("Position", fontsize=16)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)

    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)



# Function to generate RNA sequence logos
def generate_sequence_logos(data, output_folder):
    """
    Generate RNA conservation logos for sequences in the dataset.
    Args:
        data (dict): Dictionary containing sequences for analysis.
        output_folder (str): Path to save the logo plots.
    """
    os.makedirs(output_folder, exist_ok=True)

    sequences = data.get("sequences", {})
    for XXX, length_data in sequences.items():
        for key, seq_list in length_data.items():
            if seq_list:
                # Generate sequence logo
                output_file = os.path.join(output_folder, f"logo_{XXX}_{key}.png")
                generate_rna_logo(
                    seq_list,
                    output_file,
                    title=f"Sequence Logo (XXX: {XXX}, {key})"
                )

# Main function for analysis
def analyze_bulge_separation(folder, output_folder):
    """
    Main function to analyze AL bulge separation data.
    Args:
        folder (str): Input directory containing data files.
        output_folder (str): Directory to save output results.
    """
    print(f"Analyzing folder: {folder}")
    data = process_bulge_files(folder)

    # Create heatmaps
    create_heatmaps(data, os.path.join(output_folder, "heatmaps"))

    # Generate sequence logos
    generate_sequence_logos(data, os.path.join(output_folder, "logos"))

# Run the analysis
if __name__ == "__main__":
    input_folder = "results/stemloops/AL_bulges_separation"
    output_folder = input_folder
    analyze_bulge_separation(input_folder, output_folder)

