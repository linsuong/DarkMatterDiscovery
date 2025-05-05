import pandas as pd
import numpy as np

# --- File paths ---
cs_files = {
    "DM_1": "cross_sections/DM_1/cs_widths.dat",
    "DM_10": "cross_sections/DM_10/cs_widths.dat",
    "DM_100": "cross_sections/DM_100/cs_widths.dat"
}

branching_files = {
    "Z": 'Work/Z_fixed/scan2.dat',
    "W+": 'Work/W+_fixed/scan2.dat',
    "W-": 'Work/W-_fixed/scan2.dat'
}

# --- Valid scanned mass values ---
valid_MD1 = [1, 10, 20, 30, 40, 50, 60, 70, 80]
valid_DMP = [5, 20, 40, 60, 80, 120]
valid_DM3 = [1, 10, 100]

# --- Load cross section data and combine ---
cs_dfs = []
for label, path in cs_files.items():
    df = pd.read_csv(path, sep=r'\s+', low_memory=False)
    df["DM_label"] = label
    cs_dfs.append(df)

cs_all = pd.concat(cs_dfs, ignore_index=True)

# --- Filter to keep only scanned points ---
cs_all = cs_all[
    cs_all['MD1'].isin(valid_MD1) &
    cs_all['DMP'].isin(valid_DMP) &
    cs_all['DM3'].isin(valid_DM3)
].copy()


# Keep only relevant columns and drop duplicates
cs_all = cs_all.drop_duplicates()

# --- Load branching ratio data ---
dfZ = pd.read_csv(branching_files["Z"], sep=r'\s+', low_memory=False)
dfWm = pd.read_csv(branching_files["W-"], sep=r'\s+', low_memory=False)
dfWp = pd.read_csv(branching_files["W+"], sep=r'\s+', low_memory=False)

# --- Compute total BRs ---
dfZ['Br(h2->e+e-h1)_total'] = dfZ["Br(Z->e+e-)"] * dfZ["Br(h2->Zh1)"] + dfZ['Br(h2->e+e-h1)']
dfZ['Br(h2->mu+mu-h1)_total'] = dfZ["Br(Z->mu+mu-)"] * dfZ["Br(h2->Zh1)"] + dfZ['Br(h2->mu+mu-h1)']
dfZ['Br(h2->tau+tau-h1)_total'] = dfZ["Br(Z->tau+tau-)"] * dfZ["Br(h2->Zh1)"] + dfZ['Br(h2->tau+tau-h1)']
dfZ['Br(h2->n+n-h1)_total'] = dfZ["Br(Z->n+n-)"] * dfZ["Br(h2->Zh1)"] + dfZ['Br(h2->n+n-h1)']

dfWm['Br(h-->e-nu)_total'] = dfWm["Br(W-->e-nu)"] * dfWm["Br(h-->W-h1)"] + dfWm['Br(h-->e-nu)']
dfWm['Br(h-->mu-nu)_total'] = dfWm["Br(W-->mu-nu)"] * dfWm["Br(h-->W-h1)"] + dfWm['Br(h-->mu-nu)']
dfWm['Br(h-->tau-nu)_total'] = dfWm["Br(W-->tau-nu)"] * dfWm["Br(h-->W-h1)"] + dfWm['Br(h-->tau-nu)']

dfWp['Br(h+->e+nu)_total'] = dfWp["Br(W+->e+nu)"] * dfWp["Br(h+->W+h1)"] + dfWp['Br(h+->e+nu)']
dfWp['Br(h+->mu+nu)_total'] = dfWp["Br(W+->mu+nu)"] * dfWp["Br(h+->W+h1)"] + dfWp['Br(h+->mu+nu)']
dfWp['Br(h+->tau+nu)_total'] = dfWp["Br(W+->tau+nu)"] * dfWp["Br(h+->W+h1)"] + dfWp['Br(h+->tau+nu)']

# --- Keep only relevant columns and remove duplicates ---
#dfZ = dfZ[['MD1', 'DMP', 'DM3', 'Br(h2->e+e-h1)_total', 'Br(h2->mu+mu-h1)_total', 'Br(h2->tau+tau-h1)_total', 'Br(h2->n+n-h1)_total']].drop_duplicates()
#dfWm = dfWm[['MD1', 'DMP', 'DM3', 'Br(h-->e-nu)_total', 'Br(h-->mu-nu)_total', 'Br(h-->tau-nu)_total']].drop_duplicates()
#dfWp = dfWp[['MD1', 'DMP', 'DM3', 'Br(h+->e+nu)_total', 'Br(h+->mu+nu)_total', 'Br(h+->tau+nu)_total']].drop_duplicates()

# --- Merge all data ---
merged = cs_all.merge(dfZ, on=["MD1", "DMP", "DM3"], how="left") \
               .merge(dfWm, on=["MD1", "DMP", "DM3"], how="left") \
               .merge(dfWp, on=["MD1", "DMP", "DM3"], how="left")

# Optional: fill NaNs with 0 or keep as NaN
merged = merged.fillna(0)

columns_to_remove = [
    # Original columns to remove
    '~h2->e1,E1,~h1', '~h2->e2,E2,~h1', '~h2->e3,E3,~h1',
    '~h-->N1,e1,~h1', '~h-->N2,e2,~h1', '~h-->N3,e3,~h1',
    '~h+->n1,E1,~h1', '~h+->n2,E2,~h1', '~h+->n3,E3,~h1',
    'DM_label', 'MD2_x',
    'MD2', 'Br(h+->e+nu)', 'Br(h+->mu+nu)', 'Br(h+->tau+nu)',
    'Br(h+->W+h1)', 'Br(W+->e+nu)', 'Br(W+->mu+nu)', 'Br(W+->tau+nu)',
    'Br(h2->e+e-h1)', 'Br(h2->mu+mu-h1)', 'Br(h2->tau+tau-h1)',
    'Br(h2->n+n-h1)', 'Br(h2->Zh1)', 'Br(Z->e+e-)', 'Br(Z->mu+mu-)',
    'Br(Z->tau+tau-)', 'Br(Z->n+n-)',
    'MD2_y', 'Br(h-->e-nu)', 'Br(h-->mu-nu)', 'Br(h-->tau-nu)',
    'Br(h-->W-h1)', 'Br(W-->e-nu)', 'Br(W-->mu-nu)', 'Br(W-->tau-nu)', "Br(h2->n+n-h1)_total"
]
# Drop the unwanted columns
merged = merged.drop(columns=columns_to_remove, errors='ignore')  # 'errors=ignore' prevents errors if a column doesn't exist

# --- Sort by DM3 and save ---
merged = merged.sort_values(by='DM3').reset_index(drop=True)
merged.to_csv("cross_sections/combined.dat", index=False, sep='\t')

print("✅ Combined file saved as cross_sections/combined.dat (unwanted columns removed)")

print(merged)