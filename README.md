# RNA_2D_stats

A toolkit for analyzing RNA secondary structure motifs and generating statistical insights from RNA sequence-structure data.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Input Format](#input-format)
- [Output Structure](#output-structure)
- [Configuration](#configuration)
- [Citation](#citation)
- [Support](#support)

## Overview

RNA_2D_stats is a Python-based analysis pipeline designed to identify, classify, and analyze RNA secondary structure elements from sequence and dot-bracket notation data. The toolkit provides statistical analysis and visualization of structural motifs including stems, stem-loops, bulges, and interior loops, with analysis modules for proximity relationships and conservation patterns.

## Features

### Core Analysis
- **Structural Element Detection**: Automatically identifies stems, stem-loops, bulges, and interior loops
- **Motif Classification**: Categorizes structures as fully base-paired, symmetric, or asymmetric
- **Statistical Analysis**: Generates comprehensive statistics on motif occurrence and properties
- **Sequence Conservation Analysis**: Creates RNA sequence logos
- **Proximity Analysis**: Analyzes spatial relationships between structural elements

### Additional Analysis Modules
- **Apical Loop-Internal Loop Analysis**: Analysis of distance relationships between apical loops and internal loops with heatmap visualizations
- **Apical Loop-Bulge Analysis**: Distance analysis between apical loops and bulges with occurrence frequency mapping
- **Stem Analysis**: Analysis of stem sequences connecting apical loops to structural motifs with 3' and 5' end conservation
- **2D Motif-Stem Analysis**: Comprehensive analysis of stems connecting apical loops to nearest 2D motifs with sequence logo generation

### Visualization and Stats Modules
- **Length Distribution Plots**: Both absolute counts and relative frequencies
- **Conservation Logos**: RNA sequence logos for structural elements
- **Heatmaps**: Occurrence and distance relationship visualizations
- **Statistical Summaries**: Detailed numerical analysis outputs

## Input Format

The tool expects a tab-separated text file with two columns:

- **Column 1**: RNA sequence (A, U, G, C nucleotides)
- **Column 2**: Dot-bracket notation representing secondary structure

### Example:
```
UCUAGGUGAUUUCUGUGAAAUCGAGCCCACUUGAUUGUUUCUGUGAAACACUCUA	....(((((((((...))))))..))).....((.((((((...)))))).))..
GGCUUAUCAAGAGAGGUGGAGGGACUGGCCCGAUGAAACCCGGAGCCGGUAAGC	((((((((((((((((((((((((((((((((((((((((((((((((((((
```

## Installation

### Prerequisites
- Python 3.9+
- Conda or Miniconda

### Setup Environment

1. Clone the repository:
```bash
git clone https://github.com/yourusername/RNA_2D_stats.git
cd RNA_2D_stats
```

2. Create and activate the conda environment:
```bash
conda env create -f environment.yml
conda activate RNA_2D_stats
```

## Usage

### Quick Start

1. **Prepare your input file**: Ensure your RNA data is in the two-column tab-separated format described above.

2. **Run the main analysis pipeline**:
```bash
python iterative_analysis_comment.py
```
*Note: Edit line 63 in the script to specify your input file name*

3. **Generate statistics and visualizations**:
```bash
python get_stats_final.py
```

### Detailed Workflow

#### Step 1: Core Structure Analysis
The main analysis script (`iterative_analysis_comment.py`) processes your input file and:
- Identifies all structural elements (stems, stem-loops, bulges, interior loops)
- Classifies elements by type and symmetry
- Organizes results into a structured directory hierarchy

#### Step 2: Statistical Analysis
The statistics script (`get_stats_final.py`) generates:
- Length distribution plots: Both absolute counts and relative frequencies
- RNA sequence logos: Conservation analysis for each structural element type
- Statistical summaries: Detailed counts and frequencies saved as text files

#### Step 3: Proximity Analysis
Run analysis modules for detailed motif relationships:

**Apical Loop-Internal Loop Analysis**:
```bash
python analyze_AL-IL.py
```
- Analyzes distance relationships between apical loops and internal loops
- Generates occurrence and distance heatmaps
- Creates sequence conservation logos for 3' and 5' ends

**Apical Loop-Bulge Analysis**:
```bash
python analyze_AL-bulges.py
```
- Examines proximity between apical loops and bulges
- Produces occurrence frequency matrices
- Generates statistical summaries and sequence logos

**Stem Sequence Analysis**:
```bash
python analyze_AL-stems.py
```
- Analyzes stem sequences connecting apical loops to motifs
- Creates conservation logos for stem sequences
- Provides length-based statistical analysis

**2D Motif-Stem Analysis**:
```bash
python analyze_AL-2dMotif_stems.py
```
- Analysis of stems connecting apical loops to nearest 2D motifs
- Generates both 3' and 5' sequence conservation logos
- Provides structural relationship mapping

## Output Structure

The analysis creates a `results/` directory with the following organization:

```
results/
├── stems/
│   ├── fully_BP_stems/           # Fully base-paired stems
│   ├── bulges/                   # Single-sided unpaired regions
│   ├── simmetric_internal_loop/  # Symmetric interior loops
│   └── asimmetric_internal_loop/ # Asymmetric interior loops
└── stemloops/
    ├── fully_BP_stemLoops/       # Fully base-paired stem-loops
    │   └── ALvsSTEMlengh/        # Apical loop vs stem length analysis
    ├── all_apical_loops/         # Apical loop sequences
    ├── stemloop_bulges/          # Stem-loops with bulges
    ├── AL_bulges_separation/     # Apical loop-bulge proximity analysis
    ├── AL_internal_loop_separation/ # Apical loop-interior loop proximity
    └── AL_2dMotif_stems/         # Stems connecting apical loops to nearest motifs
```

### Additional Analysis Outputs

Each analysis module generates:

**Heatmaps**:
- `occurrence_heatmap.png`: Absolute occurrence frequencies
- `occurrence_relative_heatmap.png`: Relative occurrence frequencies
- `mean_distances_heatmap.png`: Average distance relationships
- `most_frequent_values_heatmap.png`: Mode distance values

**Sequence Logos**:
- `all_sequences_logo.png`: Overall sequence conservation
- `3p_logo_.png`: 3' end sequence conservation
- `5p_logo_.png`: 5' end sequence conservation
- `stem_logos/`: Directory containing motif-specific logos

**Statistical Files**:
- `*_stats`: Detailed numerical summaries
- `occurrences_stats`: Occurrence frequency data
- `mean_distances_stats`: Distance relationship statistics
- `std_distances_stats`: Standard deviation analysis

## Output Files

### Core Analysis
- **Length distribution plots**: `length_distribution.png` and `length_distribution_relative.png`
- **Sequence logos**: `logo_length{X}.png`, `3p_logo_length_{X}.png`, `5p_logo_length_{X}.png`
- **Statistics files**: `*_stats` files containing detailed numerical summaries
- **Structural data**: Tab-separated files with sequences and metadata for each motif type

### Advanced Analysis
- **Proximity heatmaps**: Occurrence and distance relationship visualizations
- **Conservation logos**: Sequence conservation analysis for structural elements
- **Statistical matrices**: Comprehensive numerical analysis of motif relationships

## Key Functions

The `functions.py` module provides core utilities:

- `dot2ct()`: Convert dot-bracket to CT list format
- `dot2bpmap()`: Convert dot-bracket to base-pair mapping
- `ct2dot()`: Convert CT list back to dot-bracket notation
- `find_stem()`: Identify stem structures with configurable parameters
- `find_stem_loop()`: Detect stem-loop motifs
- `find_bulge_interiorLoop()`: Locate bulges and interior loops within stems

## Configuration

Key parameters can be adjusted in the analysis scripts:

- **Maximum stem gap**: Default 6 nucleotides
- **Minimum stem length**: Default 4-5 nucleotides
- **Maximum loop length**: Default 10 nucleotides for stem-loops
- **Minimum stem length for logos**: Configurable in visualization scripts
- **Distance analysis parameters**: Customizable in proximity analysis modules

## Citation

If you use RNA_2D_stats in your research, please cite the appropriate manuscript (citation details to be added upon publication).

## Support

For questions, bug reports, or suggestions, please contact Jacopo Manigrasso at: [jacopo.manigrasso@gmail.com](mailto:jacopo.manigrasso@gmail.com)


Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
