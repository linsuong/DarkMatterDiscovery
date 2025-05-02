import os
import re

# === Configuration ===
input_dir = "your_directory_here"      # <-- Replace with your folder path
output_file = "cross_sections.dat"

data = []

# Loop through all .txt files in the specified directory
for filename in os.listdir(input_dir):
    if filename.endswith(".txt"):
        match = re.search(r"Mh(\d+)DMP(\d+)", filename)
        if match:
            Mh = int(match.group(1))
            DMP = int(match.group(2))
            with open(os.path.join(input_dir, filename), "r") as f:
                for line in f:
                    if "Total" in line:
                        parts = line.split()
                        sigma = float(parts[1])
                        data.append((Mh, DMP, sigma))
                        break

# Write to .dat file
with open(output_file, "w") as out:
    out.write("# Mh\tDMP\tCrossSection_fb\n")
    for Mh, DMP, sigma in sorted(data):
        out.write(f"{Mh}\t{DMP}\t{sigma:.6f}\n")

print(f"Data written to {output_file}")