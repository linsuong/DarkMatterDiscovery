import pandas as pd
import numpy as np

# Load datasets
file = 'cross_sections/combined.dat'
limits = 'cross_sections/limits_extracted_full.dat'

df = pd.read_csv(file, sep=r'\s+', low_memory=False)
df_limit = pd.read_csv(limits, sep=r'\s+', low_memory=False)

# Define valid mass points
valid_MD1 = [1, 10, 20, 30, 40, 50, 60, 70, 80]
valid_DMP = [5, 20, 40, 60, 80, 120]
valid_DM3 = [1, 10, 100]

# Filter df_limit to only keep valid points (optional, if needed)
df_limit = df_limit[
    df_limit['m_D1'].isin(valid_MD1) &
    df_limit['delta_m_plus'].isin(valid_DMP) &
    df_limit['delta_m_0'].isin(valid_DM3)
]

# Rename columns in df_limit to match df's columns
df_limit_renamed = df_limit.rename(columns={
    'm_D1': 'MD1',
    'delta_m_plus': 'DMP',
    'delta_m_0': 'DM3'
})

# Merge the sigma_2l_A column into df
df_merged = df.merge(
    df_limit_renamed[['MD1', 'DMP', 'DM3', 'sigma_2l_A', "err_2l_A"]],
    on=['MD1', 'DMP', 'DM3'],
    how='left'  # keeps all rows from df, even if no match in df_limit
)

# Calculate total cross section
br_columns = [
    'Br(h2->e+e-h1)_total',
    'Br(h2->mu+mu-h1)_total',
    'Br(h2->tau+tau-h1)_total',
    'Br(h-->e-nu)_total',
    'Br(h-->mu-nu)_total',
    'Br(h-->tau-nu)_total',
    'Br(h+->e+nu)_total',
    'Br(h+->mu+nu)_total',
    'Br(h+->tau+nu)_total'
]

df_merged['total_br'] = df[br_columns].sum(axis=1)

df_merged['total_CS'] = df_merged['total_br'] * df_merged['CrossSection_fb']
# Filter to keep only rows where sigma_2l_A is not NaN
df_filtered = df_merged[df_merged['sigma_2l_A'].notna()]
df_filtered = df_filtered.drop_duplicates(subset=['MD1', 'DMP', 'DM3'])

L0 = 39.5
L_new = 3000

df_filtered['limit_CS'] = df_filtered['sigma_2l_A'] * np.sqrt(L0/L_new)
df_filtered['limit_CS_err'] = df_filtered['err_2l_A'] * np.sqrt(L0/L_new)
df_filtered = df_filtered.sort_values(by=['DM3', 'MD1', 'DMP'])
df_filtered['Cross old limit'] = df_filtered.apply(lambda row: 'Yes' if row['total_CS'] > row['sigma_2l_A'] else 'No', axis=1)
df_filtered['Cross new limit'] = df_filtered.apply(lambda row: 'Yes' if row['total_CS'] > row['limit_CS'] else 'No', axis=1)

# Select and output the desired columns
output = ['MD1', 'DMP', 'DM3', 'total_CS', 'limit_CS', 'Cross new limit', 'sigma_2l_A', 'Cross old limit']
df_filtered[output].to_csv('cross_sections/filtered_data.dat', sep='\t', index=False, float_format='%.6e')

print("Filtered data saved to 'cross_sections/filtered_data.dat'")
print(f"Original rows: {len(df_merged)}, Filtered rows: {len(df_filtered)}")

latex_table = df_filtered[output].to_latex(
    index=False,
    float_format="%.3e",
    column_format="r" * len(output),  # right-align columns
    caption="Projected cross section limits and total signal cross sections for selected mass points.",
    label="tab:cross_section_limits"
)

with open('cross_sections/filtered_data_table.tex', 'w') as f:
    f.write(latex_table)

print("LaTeX table saved to 'cross_sections/filtered_data_table.tex'")