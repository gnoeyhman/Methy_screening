# Methy_screening
A simple python script for screening methyl group for target ligand

## Description

Use the package manager pip to install Biopython

`pip install biopython`


## Usage
Before you run the script, you need to prepare the apo structure of designed protein and selected conformation for target ligand

put the pdbs in input directory

```
indir = "~/DIRECTORY/input"
outdir = "~/DIRECTORY/output"

apo_str = "apo.pdb"
lig_molecule = "lig.pdb"
```

log_file = "log_bump.txt"
vdw_cut_off = 0.3

select_cut_off = 3.8


run the python code in PyMol session

```python
run bump_cal_6.py
```
