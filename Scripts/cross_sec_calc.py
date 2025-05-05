import os
import re
import pandas as pd

# === Configuration ===
input_dir = "Work/batch_results_DM0_1/pp-h+h-_10_step_increment/html/runs"      # <-- Replace with your folder path
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

df = pd.read_csv(output_file, sep=r'\s+', low_memory= False)
print(df)

def calcBr(df, widths, total_width, process):
    total_br = 0
    for width in widths:
        br_col = f'BR_{width}'
        df[br_col] = df[width] / df[total_width]
        total_br += df[br_col]
    df[f'totalBr_{process}'] = total_br

#process pp->h1h2
widthsh1h2 = ['~h2->e1,E1,~h1', '~h2->e2,E2,~h1', '~h2->e3,E3,~h1', 'widthsTotalh1h2']

#process pp->h+h-
widthshphm = ['~h-->N1,e1,~h1', '~h-->N2,e2,~h1 ', '~h-->N3,e3,~h1 ', 
          '~h+->n1,E1,~h1 ', '~h+->n2,E2,~h1', '~h+->n3,E3,~h1', 'widthsTotalhphm']

brhphm = calcBr(df, widthshphm, 'widthsTotalhphm', 'hphm')

df['total_CS_hphm'] = df['CS_hphm'] * df['Br_hphm']

#total cross section:
df['total_CS'] = (2 * df['total_CS_hphm']) + df['total_CS_h1h2']