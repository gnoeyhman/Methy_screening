# Methy_screening
----
A simple python script for screening methyl group for target ligand

## Description

Use the package manager pip to install Biopython

`pip install biopython`


## Usage
Before you run the script, you need to prepare the apo structure of designed protein and selected conformation for target ligand

put the pdbs in input directory

    indir = "/Users/gnoeyhman/computation/DeGrado_Lab/carbon_screening/03_04_2024_final/input_2"
    outdir = "/Users/gnoeyhman/computation/DeGrado_Lab/carbon_screening/03_04_2024_final/output"

    log_file = "log_bump.txt"
    vdw_cut_off = 0.3

    select_cut_off = 3.8

    apo_str = "apo_2.pdb"
    ref_molecule = "bas.pdb"

run the python code in PyMol session

```python
run bump_cal_6.py
```
