"""
Wrapper for mp2rage_genUniDen.py
"""
from mp2rage_genUniDen import mp2rage_genUniDen
import os

# Access inputs/outputs via the injected `snakemake` object (no import needed)
uni_path = str(snakemake.input.uni)
inv1_path = str(snakemake.input.inv1)
inv2_path = str(snakemake.input.inv2)
output_path = str(snakemake.output.uni_den)


print(f"UNI path: {uni_path} (exists: {os.path.exists(uni_path)})")
print(f"INV1 path: {inv1_path} (exists: {os.path.exists(inv1_path)})")
print(f"INV2 path: {inv2_path} (exists: {os.path.exists(inv2_path)})")

# Call your function with string paths
mp2rage_genUniDen(uni_path, inv1_path, inv2_path, output_path, 6)