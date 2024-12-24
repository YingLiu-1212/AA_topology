# AA_topology

## Usage

### Dependencies

Before running this script, please ensure that the following R packages are installed:

- `bio3d`: For reading and analyzing PDB files.
- `dplyr`: For data manipulation.

If these packages are not yet installed, you can use the following commands to install them:

```R
install.packages("bio3d", dependencies=TRUE)
install.packages("dplyr")
```

### Running the Script

The script receives the PDB file path and chain ID through command-line arguments. For example:

```bash
Rscript AA_topology.R path/to/your/file.pdb A
```

Where `path/to/your/file.pdb` is the path to the PDB file, and `A` is the chain ID you are interested in.

