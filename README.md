# Protein Amino Acid Topology Analysis and Visualization

This repository provides tools for analyzing amino acid topology from protein structures and visualizing the results using PyMOL.

## Overview

The project consists of two main components:

1. **AA_TopoAttr.R** - An R script that calculates various amino acid topology features from PDB files
2. **pymol_visualization.py** - A PyMOL script for visualizing the calculated topology properties

## Features

### Topology Analysis (AA_TopoAttr.R)
- Calculates contact flexibility (CF) metrics 
- Computes local density (LD) measures
- Generates quantile-scaled (QS) versions of all metrics
- Calculates CF-Force Polarity Index (CF-FPI) metrics
- Outputs tab-separated files with comprehensive topology attributes

### Visualization (pymol_visualization.py)
- Creates protein structure visualizations
- Colors protein chains based on selected topology properties
- Uses Red-Blue (RdBu) color scheme for intuitive value representation
- Customizable cartoon representation with smooth rendering
- High-resolution output (2400×2400 pixels)

## Prerequisites

### For Topology Analysis
- R (version 4.0+)
- Required R packages:
  - bio3d
  - dplyr
  - stringr

### For Visualization
- PyMOL (version 2.0+)
- Python (with PyMOL module)

## Installation

1. Install R dependencies:
```R
install.packages(c("bio3d", "dplyr", "stringr"))
```

2. Ensure PyMOL is installed and accessible from command line.

## Usage

### Step 1: Calculate Topology Attributes

Run the R script with a PDB file:

```bash
Rscript AA_TopoAttr.R <input_pdb_file> <output_directory>
```

**Example:**
```bash
Rscript AA_TopoAttr.R pymol_demo/1ycs.pdb pymol_demo/
```

This generates a topology attribute file: `AA_TopoAttr_1ycs.txt`

### Step 2: Visualize in PyMOL

1. Open PyMOL
2. Load and run the visualization script:
```python
run pymol_visualization.py
```

### Customizing Visualization

Edit the user-configurable parameters in `pymol_visualization.py`:

```python
# ===================== User-configurable parameters =====================
gene_name = 'TP53'
pdb_id = '1ycs'
chain_id = 'A'
plot_var = 'CF10QS'  # Options: CF10, CF10QS, LD15, LD15QS, CF10QS_FPI

project_dir = 'your/project/directory/'

MIN_VALUE = -1
MAX_VALUE = 1  
# =======================================================================
```

## Output Files

### Topology Analysis Output
- `AA_TopoAttr_<pdb_id>.txt` - Tab-separated file containing:
  - Basic information (chain, position, amino acid)
  - Contact flexibility metrics (CF10, CF10QS)
  - Local density metrics (LD15, LD15QS) 
  - CF-Force Polarity Index (CF10QS_FPI)

### Visualization Output
- `<pdb_id>_<property>_<chain>.png` - High-quality protein visualization

## Example

![Example Visualization](pymol_demo/1ycs_CF10QS_A.png)
*Figure 1: Visualization of CF10QS topology property for chain A of 1ycs PDB structure. Red regions indicate lower contact flexibility, Blue regions indicate higher contact flexibility.*


## License

This project is licensed under the MIT License - see the LICENSE file for details.
