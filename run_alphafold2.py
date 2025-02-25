#!/usr/bin/env python3
# Purpose: Run AlphaFold2 Multimer via ColabFold to generate structural models
# Inputs: 
#   - Input FASTA (e.g., preprocess/mature_${seq}.fasta)
#   - Output PDB path (e.g., preprocess/mature_${seq}.pdb)
#   - Optional: --disulfide_constraints flag
# Output: PDB file with structural model
# Notes: Optimized for multi-chain ILP prediction, including prepropeptides and propeptides
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

import subprocess
import sys
import glob
import os

input_fasta, output_pdb = sys.argv[1:3]
# Build ColabFold command with multimer model and recycling for accuracy
cmd = ["colabfold_batch", input_fasta, os.path.dirname(output_pdb), "--num-recycle", "3", "--model-type", "alphafold2_multimer_v3"]
if "--disulfide_constraints" in sys.argv:
    cmd.append("--disulfide")  # Add disulfide constraints if specified
subprocess.run(cmd, check=True)
# Move the predicted PDB to the specified output path
for f in glob.glob(f"{os.path.dirname(output_pdb)}/predict_*.pdb"):
    os.rename(f, output_pdb)
