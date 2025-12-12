# Aromatic-Sulfur Interactions Analysis

Python pipeline to analyze cysteine–aromatic interactions in protein structures.  
This repository implements the workflow used to replicate the analysis shown in **Figure 4** of:

> J. Phys. Chem. Lett. 2021, 12, 2865−2870.

The code identifies aromatic ring centers and cysteine sulfur atoms, calculates distances and angles between them, checks for steric obstruction, and categorizes interactions into **sigma**, **mid**, and **pi** types.

---

## Features

- Downloads PDB structures automatically from RCSB.
- Extracts all atom coordinates, focusing on chain A.
- Computes aromatic ring centers and plane equations for HIS, PHE, TYR, and the benzene ring of TRP.
- Calculates vector from cysteine sulfur (SG) to aromatic ring center.
- Computes the angle between this vector and the aromatic ring plane.
- Filters interactions based on distance thresholds (3–10.5 Å).
- Checks for steric obstruction by other atoms near the line connecting sulfur to ring center.
- Categorizes angles into:
  - **sigma**: 0–20°  
  - **mid**: 20–70°  
  - **pi**: 70–90°  

---

## Installation

Clone the repository:

```bash
git clone https://github.com/jakedixon18/aromatic-sulfur-interactions.git
cd aromatic-sulfur-interactions
```
Make sure you have the following dependencies installed:

Python 3.8+

NumPy

Install via pip if necessary:
```bash
pip install numpy
```

---

## Usage
```python
from aromatic_sulfur_analysis import run_protein_structure_solver

# Example: analyze PDB structure 4F5S
final_output = run_protein_structure_solver("4F5S")

# Print first 10 entries
for row in final_output[:10]:
    print(row)
```

---

## Output
```markdown
The pipeline produces a list of tuples, each containing:

cys_resname – residue name of the cysteine

cys_resnum – residue number of the cysteine

aro_resname – aromatic residue name

aro_resnum – aromatic residue number

distance – distance from sulfur to aromatic center (Å)

angle_type – "sigma", "mid", or "pi"

obstruction – "obstructed" or "un-obstructed"

Example:

('CYS', 42, 'PHE', 113, 7.2, 'mid', 'un-obstructed')
```

---

## References

J. Phys. Chem. Lett. 2021, 12, 2865−2870
DOI: https://doi.org/10.1021/acs.jpclett.1c00415
