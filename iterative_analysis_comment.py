# RNA Structural Element Analysis Script
# ========================================
# This script processes RNA sequences and their dot-bracket notations to identify and analyze 
# structural elements such as stems, stem-loops, bulges, and interior loops. It organizes the results 
# into a categorized output directory structure. The script is designed to handle multiple RNA sequences 
# from a CSV file and provide detailed information about their secondary structure features.
#
# Key functionalities:
# 1. **Identification of RNA structural elements**:
#    - Detect stems and stem-loops based on sequence structure and defined constraints.
#    - Differentiate between fully base-paired structures and those containing bulges or interior loops.
#
# 2. **Analysis of stems**:
#    - Classify stems into fully base-paired and non-fully base-paired (with bulges or interior loops).
#    - Extract information about:
#        - Length of the stem.
#        - Sequences of the left and right arms.
#    - For stems with bulges or interior loops:
#        - Calculate their lengths.
#        - Extract the sequences involved.
#    - Classify bulges and interior loops as symmetric or asymmetric, storing this classification in separate files.
#
# 3. **Analysis of stem-loops**:
#    - Identify stem-loops with fully base-paired stems and loops.
#    - Extract information about:
#        - Length of the stem and loop.
#        - Sequences of the left arm, right arm, and the loop.
#    - For stem-loops with bulges or interior loops:
#        - Identify and classify each structural feature.
#        - Extract sequences and calculate distances from the apical loop.
#
# 4. **Proximity analysis**:
#    - For stem-loops with both bulges and interior loops:
#        - Merge the structural features into a unified list.
#        - Sort the features by proximity to the apical loop.
#        - Identify the closest structural feature to the apical loop.
#        - Record detailed information about this feature, including:
#            - Distance from the apical loop.
#            - Sequence details of the bulge or interior loop.
#
# 5. **Output**:
#    - Organizes results into a `results` folder with two main subdirectories:
#        - `stems`: Contains information about fully base-paired stems, bulges, symmetric, and asymmetric loops.
#        - `stemloops`: Contains information about fully base-paired stem-loops, apical loops, and proximity analyses.
#    - Outputs are stored in text files, with naming conventions indicating the structural feature lengths and types.
#
# Input:
# - A CSV file where:
#    - The first column contains RNA sequences.
#    - The second column contains dot-bracket notations corresponding to the sequences.
#
# Output:
# - Categorized text files within a "results" folder, storing sequences and metadata for identified RNA structural features.
#
#
# Abbreviations and info:
# AL = apical loop
# IL = internal loop = Regions where unpaired bases occur on both sides of the stem.
# Bulge = Regions where unpaired bases occur on only one side of the stem.


# Read sequences and dot-brackets from file
input_file = "input_example_file.txt"  # Replace with your actual file name


import os
import csv
import math
from functions import find_stem, find_stem_loop, find_bulge_interiorLoop  # Import necessary functions

# Define output folder structure
base_dir = "results"
stems_dir = os.path.join(base_dir, "stems")
stemloops_dir = os.path.join(base_dir, "stemloops")

# Ensure directories exist
os.makedirs(stems_dir, exist_ok=True)
os.makedirs(stemloops_dir, exist_ok=True)
os.makedirs(os.path.join(stems_dir, "fully_BP_stems"), exist_ok=True)
os.makedirs(os.path.join(stems_dir, "bulges"), exist_ok=True)
os.makedirs(os.path.join(stems_dir, "simmetric_internal_loop"), exist_ok=True)
os.makedirs(os.path.join(stems_dir, "asimmetric_internal_loop"), exist_ok=True)
os.makedirs(os.path.join(stemloops_dir, "fully_BP_stemLoops"), exist_ok=True)
os.makedirs(os.path.join(stemloops_dir, "fully_BP_stemLoops/ALvsSTEMlength"), exist_ok=True)
os.makedirs(os.path.join(stemloops_dir, "all_apical_loops"), exist_ok=True)
os.makedirs(os.path.join(stemloops_dir, "stemloop_bulges"), exist_ok=True)
os.makedirs(os.path.join(stemloops_dir, "stemloop_simmetric_internal_loop"), exist_ok=True)
os.makedirs(os.path.join(stemloops_dir, "asimmetric_internal_loop"), exist_ok=True)
os.makedirs(os.path.join(stemloops_dir, "AL_bulges_separation"), exist_ok=True)
os.makedirs(os.path.join(stemloops_dir, "AL_internal_loop_separation"), exist_ok=True)

os.makedirs(os.path.join(stemloops_dir, "AL_2dMotif_stems"), exist_ok=True)

# Function to write results to a file
def append_to_file(folder, filename, data):
    """
    Append a given string to a file in the specified folder.
    If the file does not exist, it will be created.
    """
    with open(os.path.join(folder, filename), "a") as f:
        f.write(data + "\n")

# Open the input file and process each sequence and dot-bracket pair
with open(input_file, "r") as file:
    reader = csv.reader(file)
    for row in reader:
        seq = row[0].strip()       # First column: Sequence
        dot = row[1].strip()       # Second column: Dot-bracket notation

        # Step 1: Identify all stems with at most 6 nt gaps
        all_2d = find_stem(dot, max_stem_gap=6, min_stem_len=4)

        # Step 2: Separate stem-loops and stems
        stemLoops = find_stem_loop(dot, max_loop_len=10, max_stem_gap=6, min_stem_len=4)
        stems = [stem for stem in all_2d if stem not in stemLoops]

        # Step 3: Analyze stems
        for stem in stems:
            # a. Check if the stem has bulges or interior loops
            bulges, interiorLoops = find_bulge_interiorLoop(dot, stem)

            # b. Handle fully base-paired stems (no bulges or interior loops)
            if len(interiorLoops) == 0 and len(bulges) == 0:
                # Calculate length of the stem
                length = stem[1] - stem[0] + 1
                # Extract left and right sequences
                left_seq = seq[stem[0]-1:stem[1]]
                right_seq = seq[stem[2]-1:stem[3]]
                # Write results to the "fully_BP_stems" directory
                filename = f"stem_{length}nts.txt"
                data = f"{length}\t{left_seq}\t{right_seq}"
                append_to_file(os.path.join(stems_dir, "fully_BP_stems"), filename, data)

            # c. Handle stems with bulges
            for bulge in bulges:
                # Measure length and extract the sequence of the interior loop
                il_len = bulge[1] - bulge[0] + 1
                il_seq = seq[bulge[0]-1:bulge[1]]
                # Write results to the "bulges" directory
                filename = f"bulges_{il_len}nt.txt"
                append_to_file(os.path.join(stems_dir, "bulges"), filename, il_seq)

            # d. Handle stems with interiorLoops
            for interiorLoop in interiorLoops:
                # Measure lengths of the left and right sides of the bulge
                ls_len = interiorLoop[1] - interiorLoop[0] + 1
                rs_len = interiorLoop[3] - interiorLoop[2] + 1
                # Extract left and right sequences
                left_seq = seq[interiorLoop[0]-1:interiorLoop[1]]
                right_seq = seq[interiorLoop[2]-1:interiorLoop[3]]
                # Write to appropriate directory based on symmetry
                if ls_len == rs_len:
                    filename = f"IL_{ls_len}x{ls_len}.txt"
                    append_to_file(os.path.join(stems_dir, "simmetric_internal_loop"), filename, f"{left_seq}\t{right_seq}")
                else:
                    filename = f"IL_{ls_len}x{rs_len}.txt"
                    append_to_file(os.path.join(stems_dir, "asimmetric_internal_loop"), filename, f"{left_seq}\t{right_seq}")

            # e. other sections - in case



        # Step 4: Analyze stem-loops
        for stemLoop in stemLoops:
            # a. Check for bulges and interior loops
            bulges, interiorLoops = find_bulge_interiorLoop(dot, stemLoop)

            # b. Fully base-paired stem-loops
            if len(interiorLoops) == 0 and len(bulges) == 0:
                # Calculate stem and loop lengths
                stem_len = stemLoop[1] - stemLoop[0] + 1
                loop_len = stemLoop[2] - stemLoop[1] - 1
                # Extract left, right, and loop sequences
                left_seq = seq[stemLoop[0]-1:stemLoop[1]]
                right_seq = seq[stemLoop[2]-1:stemLoop[3]]
                loop_seq = seq[stemLoop[1]:stemLoop[2]-1]
                # Write results to "fully_BP_stemLoops" and "all_apical_loops"
                filename = f"SL_{stem_len}nts.txt"
                data = f"{stem_len}\t{left_seq}\t{right_seq}"
                append_to_file(os.path.join(stemloops_dir, "fully_BP_stemLoops"), filename, data)             
                filename = f"AL_{loop_len}nt.txt"
                append_to_file(os.path.join(stemloops_dir, "all_apical_loops"), filename, loop_seq)

                #Write AL vs Stem length files
                filename = f"AL_{loop_len}nt_{stem_len}x{stem_len}stem.txt"
                data = f"{loop_len}\t{left_seq}\t{right_seq}"
                append_to_file(os.path.join(stemloops_dir, "fully_BP_stemLoops/ALvsSTEMlength"), filename, data)

            # c. Stem-loops with bulges
            for bulge in bulges:
                il_len = bulge[1] - bulge[0] + 1
                il_seq = seq[bulge[0]-1:bulge[1]]
                loop_len = stemLoop[2] - stemLoop[1] - 1
                loop_seq = seq[stemLoop[1]:stemLoop[2]-1]
                filename = f"bulges_{il_len}nt_{loop_len}SL.txt"
                append_to_file(os.path.join(stemloops_dir, "stemloop_bulges"), filename, il_seq)
                filename = f"AL_{loop_len}nt.txt"
                append_to_file(os.path.join(stemloops_dir, "all_apical_loops"), filename, loop_seq)

            # d. Stem-loops with interiorLoop
            for interiorLoop in interiorLoops:
                loop_len = stemLoop[2] - stemLoop[1] - 1
                loop_seq = seq[stemLoop[1]:stemLoop[2]-1]
                filename = f"AL_{loop_len}nt.txt"
                append_to_file(os.path.join(stemloops_dir, "all_apical_loops"), filename, loop_seq)
                #
                ls_len = interiorLoop[1] - interiorLoop[0] + 1
                rs_len = interiorLoop[3] - interiorLoop[2] + 1
                left_seq = seq[interiorLoop[0]-1:interiorLoop[1]]
                right_seq = seq[interiorLoop[2]-1:interiorLoop[3]]
                if ls_len == rs_len:
                    filename = f"IL_{ls_len}x{ls_len}_{loop_len}SL.txt"
                    append_to_file(os.path.join(stemloops_dir, "stemloop_simmetric_internal_loop"), filename, f"{left_seq}\t{right_seq}")
                else:
                    filename = f"IL_{ls_len}x{rs_len}_{loop_len}SL.txt"
                    append_to_file(os.path.join(stemloops_dir, "asimmetric_internal_loop"), filename, f"{left_seq}\t{right_seq}")


            # e. If bulges or interior loops exist, analyze them for proximity to the apical loop
            if len(bulges) > 0 or len(interiorLoops) > 0:
                # e1. Merge bulges and interior loops into a single list and sort
                il_bulge = bulges + interiorLoops
                il_bulge.sort()

                # e2. Compute half of the stemLoop length
                half_sL = math.floor((stemLoop[3] - stemLoop[0]) / 2)

                # e3-e4. Find the closest element to the apical loop
                min_p = float('inf')
                closest_element = None

                for element in il_bulge:
                    if len(element) == 4:  # It's a bulge
                        if element[0] < stemLoop[1]:  # Bulge on the left
                            p = abs(stemLoop[1] - element[1])
                        else:  # Bulge on the right
                            p = abs(stemLoop[2] - element[2])
                    elif len(element) == 2:  # It's an interior loop
                        p = abs(stemLoop[1] - element[1])

                    # If p > half_sL, recalculate p from the right side
                    if p > half_sL:
                        p = abs(stemLoop[2] - element[1])

                    # Update minimum p and closest element
                    if p < min_p:
                        min_p = p
                        closest_element = element

                # e5. Process the closest element
                if len(closest_element) == 2:  # Interior loop
                    il_len = closest_element[1] - closest_element[0] + 1
                    il_seq = seq[closest_element[0] - 1: closest_element[1]]
                    loop_len = stemLoop[2] - stemLoop[1] - 1
                    filename = f"AL_{loop_len}nt_{il_len}bulge.txt"
                    data = f"{min_p}\t{il_seq}"
                    append_to_file(os.path.join(stemloops_dir, "AL_bulges_separation"), filename, data)

                    #get the sequence of the stem between the closest element and the apical loop
                    if closest_element[0] > stemLoop[1]: # it's a 5' bulge
                        right_seq = seq[stemLoop[2] -1 : closest_element[1] - 1]
                        start_3p_stem= stemLoop[1] - len(right_seq)
                        stem_len = len(right_seq) 
                        left_seq = seq[start_3p_stem : stemLoop[1]]
                        
                        filename=f"AL_{loop_len}nt_{il_len}bulge_{stem_len}ntStem.txt"
                        data = f"{stem_len}\t{left_seq}\t{right_seq}"
                        append_to_file(os.path.join(stemloops_dir, "AL_2dMotif_stems"), filename, data)
                    else:
                        left_seq = seq[closest_element[0]: stemLoop[1]]
                        end_5p_stem = stemLoop[2] + len(left_seq)
                        right_seq = seq[stemLoop[2] -1: end_5p_stem - 1]
                        stem_len = len(left_seq)

                        filename=f"AL_{loop_len}nt_{il_len}bulge_{stem_len}ntStem.txt"
                        data = f"{stem_len}\t{left_seq}\t{right_seq}"
                        append_to_file(os.path.join(stemloops_dir, "AL_2dMotif_stems"), filename, data)



                elif len(closest_element) == 4:  # Bulge
                    ls_len = closest_element[1] - closest_element[0] + 1
                    rs_len = closest_element[3] - closest_element[2] + 1
                    left_seq = seq[closest_element[0] - 1: closest_element[1]]
                    right_seq = seq[closest_element[2] - 1: closest_element[3]]
                    loop_len = stemLoop[2] - stemLoop[1] - 1
                    filename = f"AL_{loop_len}nt_{ls_len}x{rs_len}IL.txt"
                    data = f"{min_p}\t{left_seq}\t{right_seq}"
                    append_to_file(os.path.join(stemloops_dir, "AL_internal_loop_separation"), filename, data)

                    #get the sequence of the stem between the closest element and the apical loop

                    left_seq = seq[closest_element[1]: stemLoop[1]]
                    right_seq = seq[stemLoop[2] - 1: closest_element[2] -1]
                    len_left_seq = len(left_seq)
                    filename = f"AL_{loop_len}nt_{ls_len}x{rs_len}IL_{len_left_seq}ntStem.txt"
                    data = f"{len_left_seq}\t{left_seq}\t{right_seq}"
                    append_to_file(os.path.join(stemloops_dir, "AL_2dMotif_stems"), filename, data)



