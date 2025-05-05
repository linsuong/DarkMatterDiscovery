import os
import re
import pandas as pd
import numpy as np

# === Configuration ===
"""
input_dirs = [
            'Work/batch_results/DM0_10/pp-h1h2_10_step_increment/html/runs/',
             'Work/batch_results/DM0_10/pp-h1h2_mh1_1/html/runs/',
             'Work/batch_results/DM0_10/pp-h+h-_10_step_increment/html/runs/',
             'Work/batch_results/DM0_10/pp-h+h-_mh1_1/html/runs/'
             ]
"""

input_dirs = [
            'Work/batch_results/DM0_1/pp-h1h2_10_step_increment/html/runs/',
             'Work/batch_results/DM0_1/pp-h1h2_mh1_1/html/runs/',
             'Work/batch_results/DM0_1/pp-hphm_10_step_increment/html/runs/',
             'Work/batch_results/DM0_1/pp-hphm_mh1_1/html/runs/'
             ]

output_file = "cross_sections/DM_1/cs_widths.dat"

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

for input_dir in input_dirs:
    for filename in os.listdir(input_dir):
        if not filename.endswith(".txt"):
            continue

        Mh, DMP = None, None

        # Try full match: Mh10DMP20
        match = re.search(r"Mh1(\d+)DMP(\d+)", filename)
        if match:
            Mh = int(match.group(1))
            DMP = int(match.group(2))
        else:
            # Try fallback: DMP10 → assume Mh = 1
            match_simple = re.search(r"DMP(\d+)", filename)
            if match_simple:
                Mh = 1
                DMP = int(match_simple.group(1))
            else:
                continue  # skip if pattern not matched

        filepath = os.path.join(input_dir, filename)
        with open(filepath, "r") as f:
            lines = f.readlines()

        # === Initialize ===
        sigma = np.nan
        widths = {key: np.nan for key in width_labels}

        # Parse total cross section
        for line in lines:
            if line.strip().startswith("Total"):
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        sigma = float(parts[1])
                    except ValueError:
                        pass
                break

        # Parse decay widths
        for line in lines:
            for label in width_labels:
                if label in line:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            widths[label] = float(parts[1])
                        except ValueError:
                            pass

        # Combine data
        row = {'MD1': Mh, 'DMP': DMP, 'DM3': 1, 'CrossSection_fb': sigma}
        row.update(widths)
        data.append(row)

# === Convert to DataFrame and Save ===
df = pd.DataFrame(data)
os.makedirs(os.path.dirname(output_file), exist_ok=True)
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