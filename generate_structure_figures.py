#!/usr/bin/env python3
# Purpose: Generate 3D structure figures for TM-align references and predicted ILPs
# Inputs: output figures dir, preprocess/, analysis/pdbs/
# Outputs: output/figures/ref_*.png, output/figures/*_*.png
# Notes: Uses PyMOL for rendering
# Author: Jorge L. Pérez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

import glob
import os
import subprocess
import sys
import yaml

output_dir, preprocess_dir, pdbs_dir = sys.argv[1:4]
with open("config.yaml") as f:
    config = yaml.safe_load(f)

os.makedirs(output_dir, exist_ok=True)

# Reference structures
for ref_name, ref_path in config["tmalign_refs"].items():
    if os.path.exists(ref_path):
        subprocess.run(["pymol", "-c", "-d", f"load {ref_path}; zoom; ray; save {output_dir}/ref_{ref_name}.png; quit"])

# Predicted ILPs
for type in ["prepro", "pro", "mature"]:
    for pdb in glob.glob(f"{preprocess_dir}/{type}_*.pdb") + glob.glob(f"{pdbs_dir}/{type}_*.pdb"):
        base = os.path.basename(pdb).replace(".pdb", "")
        subprocess.run(["pymol", "-c", "-d", f"load {pdb}; zoom; ray; save {output_dir}/{base}.png; quit"])
