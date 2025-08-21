"""
RNA Sequence Data Analysis Script
=================================

Description:
------------
This script performs a comprehensive analysis of RNA sequence data stored in text files. It generates various visualizations and statistics to study the properties of RNA sequences.

Key Features:
-------------
1. **Load and parse data**: Reads RNA sequence files with varying numbers of columns (1, 2, or 3).
2. **Length distribution plots**: Generates bar plots for length distributions (absolute counts and relative frequencies).
3. **RNA conservation logos**: Creates RNA sequence logos for visualization of sequence conservation.
4. **Output statistics**: Saves statistics (counts and relative frequencies) in sorted text files.

Input:
------
The script expects folders containing `.txt` files with the following formats:
- **1-column files**: Single sequences.
- **2-column files**: 3' sequences and 5' sequences (e.g., symmetric/asymmetric loops).
- **3-column files**: Sequence length, 3' sequences, and 5' sequences (e.g., fully base-paired stems).

Output:
-------
1. Bar plots (`.png`): Visualizing length distributions.
2. Sequence conservation logos (`.png`): For single sequences, 3' sequences, and 5' sequences.
3. Statistics (`_stats` text files): Counts and relative frequencies of sequence lengths.

Dependencies:
-------------
- Python libraries: `os`, `pandas`, `logomaker`, `matplotlib`, `seaborn`, `collections`.

Execution:
----------
1. Modify `folders_to_analyze` to specify the directories and expected column numbers.
2. Run the script. Outputs are saved in the respective directories.

Author: [Your Name]
Date: [Current Date]
"""

import os
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

# Function to generate RNA Logo using logomaker
def generate_rna_logo(sequences, output_file, title="RNA Conservation Logo", reverse=False):
    """
    Generate an RNA conservation logo from a list of RNA sequences using logomaker.
    Args:
        sequences (list): List of RNA sequences (strings).
        output_file (str): Path to save the RNA logo as an image.
        title (str): Title for the logo.
        reverse (bool): Reverse sequences for visualization.
    """
    if reverse:
        # Reverse sequences if specified
        sequences = [seq[::-1] for seq in sequences]

    # Convert sequences to a DataFrame where each row represents one sequence
    data = pd.DataFrame([list(seq) for seq in sequences])

    # Create a frequency matrix (proportional counts for each nucleotide at each position)
    counts_df = data.apply(lambda x: x.value_counts(normalize=True)).fillna(0).T

    # Adjust plot width dynamically based on the number of positions in the sequences
    num_columns = counts_df.shape[0]
    width = max(2, num_columns * 2)  # Ensure sufficient width for long sequences
    height = 6

    # Set up the figure
    fig, ax = plt.subplots(figsize=(width, height), dpi=300)

    # Generate the logo
    logo = logomaker.Logo(counts_df, ax=ax, color_scheme="classic")

    # Style adjustments
    logo.style_xticks(rotation=90, anchor=0)
    ax.set_title(title, fontsize=18)
    ax.set_ylabel("Frequency", fontsize=16)
    ax.set_xlabel("Position", fontsize=16)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)

    # Save the logo
    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)  # Close the figure to prevent memory issues
    print(f"RNA logo saved to {output_file}")

# Function to load data from files in a directory
def load_files_from_dir(folder, num_columns):
    """
    Load RNA sequence data from files in a specified directory.
    Args:
        folder (str): Path to the directory containing files.
        num_columns (int): Number of expected columns in the files.
    Returns:
        dict: A dictionary where keys are lengths and values are sequence data.
        dict: A dictionary mapping filenames to their lengths.
    """
    data = defaultdict(list)
    file_labels = {}

    for dirpath, _, filenames in os.walk(folder):
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            if filename.endswith('.txt'):
                print(f"Processing file: {filepath}")
                with open(filepath, "r") as f:
                    seq_lengths = None
                    for line in f:
                        fields = line.strip().split("\t")
                        if len(fields) == num_columns:
                            # Process based on the number of columns
                            if num_columns == 3:
                                length = int(fields[0])
                                seq_3p = fields[1]
                                seq_5p = fields[2]
                                data[length].append((seq_3p, seq_5p))
                            elif num_columns == 2:
                                seq_3p = fields[0]
                                seq_5p = fields[1]
                                seq_lengths = f"{len(seq_3p)}x{len(seq_5p)}"
                                data[seq_lengths].append((seq_3p, seq_5p))
                            elif num_columns == 1:
                                seq = fields[0]
                                data[len(seq)].append(seq)
                        else:
                            print(f"Skipping malformed line: {line.strip()}")
                    if seq_lengths:
                        file_labels[filename] = seq_lengths

    return data, file_labels

# Function to plot the distribution of lengths
def plot_length_distribution(data, output_file, relative=False):
    """
    Plot the distribution of lengths as a bar plot and save it to a file.
    Additionally, save statistics (counts and relative frequencies) to a .txt file.
    
    Args:
        data (dict): Dictionary with lengths as keys and sequences as values.
        output_file (str): Path to save the plot.
        relative (bool): Whether to plot relative proportions.
    """
    if not data:
        print("No data found for length distribution.")
        return

    # Prepare data for plotting and statistics
    length_labels = list(data.keys())
    counts = [len(seqs) for seqs in data.values()]
    total_count = sum(counts)
    relative_frequencies = [count / total_count for count in counts]

    # Create a DataFrame and sort by count
    df = pd.DataFrame({
        'Length': length_labels,
        'Count': counts,
        'Relative Frequency': relative_frequencies
    }).sort_values(by='Count')

    # Save statistics to a text file
    stat_file = os.path.splitext(output_file)[0] + "_stats"
    with open(stat_file, "w") as f:
        f.write(f"Statistics for {output_file}\n\n")
        f.write(df.to_string(index=False))
    print(f"Statistics saved to {stat_file}")
    
    # sort by Length to plot
    df=df.sort_values(by='Length')
    # Plot the distribution
    plt.figure(figsize=(9, 6), dpi=300)
    ax = sns.barplot(x='Length', y='Relative Frequency' if relative else 'Count', data=df, edgecolor="black")

    # Customize fonts and add black outlines to bars
    ax.set_xlabel("Length", fontsize=16)
    ax.set_ylabel("Proportion" if relative else "Count", fontsize=16)
    ax.set_title("Length Distribution" + (" (Relative)" if relative else ""), fontsize=18)
    ax.tick_params(axis='x', labelsize=14)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.tick_params(axis='y', labelsize=14)

    for bar in ax.patches:
        bar.set_edgecolor('black')
        bar.set_linewidth(1.5)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Barplot saved to {output_file}")

# Main function to handle all analyses
def analyze_folder(folder, num_columns, is_reverse=False):
    """
    Perform comprehensive analysis of a folder with RNA sequence data.
    Args:
        folder (str): Path to the folder containing RNA data files.
        num_columns (int): Number of columns in the files.
        is_reverse (bool): Whether to reverse sequences for the logo.
    """
    print(f"Analyzing folder: {folder}")
    data, file_labels = load_files_from_dir(folder, num_columns)

    if not data:
        print(f"No valid data found in {folder}. Skipping.")
        return

    # Plot length distributions
    plot_length_distribution(data, os.path.join(folder, "length_distribution.png"))
    plot_length_distribution(data, os.path.join(folder, "length_distribution_relative.png"), relative=True)

    # Generate RNA Logos
    for length, seqs in data.items():
        if num_columns == 1:
            # Single-column sequences
            generate_rna_logo(seqs, os.path.join(folder, f"logo_length_{length}.png"), title=f"Bulge")
        else:
            # Two- or three-column sequences
            seq_3p = [seq[0] for seq in seqs]
            seq_5p = [seq[1] for seq in seqs]
            generate_rna_logo(seq_3p, os.path.join(folder, f"3p_logo_length_{length}.png"), title=f"3' sequence (Length {length})")
            generate_rna_logo(seq_5p, os.path.join(folder, f"5p_logo_length_{length}_reversed.png"), title=f"5' sequence (Length {length})", reverse=is_reverse)

# Define folders and expected columns
folders_to_analyze = {
    "results/stems/fully_BP_stems/": 3,
    "results/stems/simmetric_internal_loop/": 2,
    "results/stems/asimmetric_internal_loop/": 2,
    "results/stems/bulges/": 1,
    "results/stemloops/fully_BP_stemLoops/": 3,
    "results/stemloops/stemloop_simmetric_internal_loop": 2,
    "results/stemloops/asimmetric_internal_loop": 2,
    "results/stemloops/stemloop_bulges": 1,
    "results/stemloops/all_apical_loops": 1,
}

# Run analysis for each folder
if __name__ == "__main__":
    for folder, num_columns in folders_to_analyze.items():
        analyze_folder(folder, num_columns, is_reverse=(num_columns > 1))

