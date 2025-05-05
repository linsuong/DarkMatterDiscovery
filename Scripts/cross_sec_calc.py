import os
import re
import pandas as pd
import numpy as np

# === Configuration ===
input_dir = "Work/batch_results/pp-h1h2_10_step_increment/html/runs"  # Replace with your actual path
output_file = "cs_widths.dat"

data = []

# Decay width channels to search for
width_labels = {
    '~h2->e1,E1,~h1': '~h2->e1,E1,~h1',
    '~h2->e2,E2,~h1': '~h2->e2,E2,~h1',
    '~h2->e3,E3,~h1': '~h2->e3,E3,~h1',
    '~h-->N1,e1,~h1': '~h-->N1,e1,~h1',
    '~h-->N2,e2,~h1': '~h-->N2,e2,~h1',
    '~h-->N3,e3,~h1': '~h-->N3,e3,~h1',
    '~h+->n1,E1,~h1': '~h+->n1,E1,~h1',
    '~h+->n2,E2,~h1': '~h+->n2,E2,~h1',
    '~h+->n3,E3,~h1': '~h+->n3,E3,~h1'
}

# Loop through .txt files
for filename in os.listdir(input_dir):
    if filename.endswith(".txt"):
        match = re.search(r"Mh(\d+)DMP(\d+)", filename)
        if match:
            Mh = int(match.group(1))
            DMP = int(match.group(2))
            filepath = os.path.join(input_dir, filename)
            with open(filepath, "r") as f:
                lines = f.readlines()

            # Initialize values
            sigma = np.nan
            widths = {key: np.nan for key in width_labels}

            # === Parse cross section ===
            for line in lines:
                if line.strip().startswith("Total"):
                    parts = line.split()
                    try:
                        sigma = float(parts[1])
                    except (IndexError, ValueError):
                        pass
                    break

            # === Parse decay widths ===
            for line in lines:
                for label in width_labels:
                    if label in line:
                        parts = line.strip().split()
                        try:
                            width_value = float(parts[1])
                            widths[label] = width_value
                        except (IndexError, ValueError):
                            pass

            # Combine data
            row = {'Mh': Mh, 'DMP': DMP, 'CrossSection_fb': sigma}
            row.update(widths)
            data.append(row)

# Convert to DataFrame
df = pd.DataFrame(data)

# Save to .dat file
df.to_csv(output_file, sep="\t", index=False, float_format="%.6e", na_rep="NaN")
print(f"Written to {output_file}")


'''
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
df['total_CS'] = (2 * df['total_CS_hphm']) + df['total_CS_h1h2']'''