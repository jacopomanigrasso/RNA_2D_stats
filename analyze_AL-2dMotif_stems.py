import os
import re
import pandas as pd
import logomaker
import matplotlib.pyplot as plt


# Function to create RNA sequence logo
def generate_rna_logo(sequences, output_file, title="RNA Conservation Logo", reverse=False):
    """
    Generate an RNA conservation logo.
    Args:
        sequences (list): List of RNA sequences.
        output_file (str): Path to save the logo.
        title (str): Title of the logo.
        reverse (bool): Whether to reverse the sequences for visualization.
    """
    if reverse:
        print("Reversing sequences for 5' logo generation...")
        [seq[::-1] for seq in sequences]  #if you want to reverse the 5'end sequence for plotting purposes, do: sequences = [seq[::-1] for seq in sequences]

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


# Function to generate logos for each file
def generate_logos_for_files(input_folder, output_folder):
    """
    Process files and generate sequence logos for 3' and 5' sequences.
    Args:
        input_folder (str): Folder containing input files.
        output_folder (str): Folder to save generated logos.
    """
    # Create output folder if it doesn't exist
    logo_folder = os.path.join(output_folder, "stem_logos")
    os.makedirs(logo_folder, exist_ok=True)

    # Regex patterns for matching filenames
    file_pattern_1 = re.compile(r"^AL_\d+nt_\d+x\d+IL_\d+ntStem\.txt$")
    file_pattern_2 = re.compile(r"^AL_\d+nt_\d+bulge_\d+ntStem\.txt$")

    # Process each file in the folder
    for dirpath, _, filenames in os.walk(input_folder):
        for filename in filenames:
            if file_pattern_1.match(filename) or file_pattern_2.match(filename):
                filepath = os.path.join(dirpath, filename)
                print(f"Processing file: {filepath}")

                # Read the file
                sequences_3p = []
                sequences_5p = []
                with open(filepath, "r") as f:
                    for line in f:
                        fields = line.strip().split("\t")
                        if len(fields) != 3:
                            print(f"Skipping malformed line: {line.strip()}")
                            continue
                        _, seq_3p, seq_5p = fields
                        sequences_3p.append(seq_3p)
                        sequences_5p.append(seq_5p)

                #Reverse 5' sequences for plotting purposes
                sequences_5p = [seq[::-1] for seq in sequences_5p]

                # Generate logos for 3' and 5' sequences
                if sequences_3p:
                    output_file_3p = os.path.join(logo_folder, f"{filename}_3p_logo.png")
                    generate_rna_logo(
                        sequences_3p,
                        output_file_3p,
                        title=f"3' Sequence Logo for {filename}"
                    )

                if sequences_5p:
                    output_file_5p = os.path.join(logo_folder, f"{filename}_5p_logo.png")
                    generate_rna_logo(
                        sequences_5p,
                        output_file_5p,
                        title=f"5' Sequence Logo for {filename}",
                        reverse=False,
                    )


# Main function for analysis
def main(input_folder, output_folder):
    """
    Main function to analyze RNA data and generate sequence logos.
    Args:
        input_folder (str): Input directory containing RNA data files.
        output_folder (str): Directory to save output results.
    """
    print(f"Analyzing folder: {input_folder}")
    generate_logos_for_files(input_folder, output_folder)


# Run the analysis
if __name__ == "__main__":
    input_folder = "results/stemloops/AL_2dMotif_stems/selection_bugles"  # Replace with your input folder path
    output_folder = input_folder  # Replace with your output folder path
    main(input_folder, output_folder)

