import pandas as pd 
import numpy as np 

file = "cross_sections/DM_100/cs_widths.dat"
output_file = "cross_sections/DM_100/cs_calculated.dat"

df = pd.read_csv(file, sep=r'\s+', low_memory= False)

decay_cols = [
    '~h2->e1,E1,~h1', '~h2->e2,E2,~h1', '~h2->e3,E3,~h1',
    '~h-->N1,e1,~h1', '~h-->N2,e2,~h1', '~h-->N3,e3,~h1',
    '~h+->n1,E1,~h1', '~h+->n2,E2,~h1', '~h+->n3,E3,~h1'
]

df['totalWidth'] = df[decay_cols].sum(axis=1)

decayh2 = ['~h2->e1,E1,~h1', '~h2->e2,E2,~h1', '~h2->e3,E3,~h1']
decayhp = ['~h-->N1,e1,~h1', '~h-->N2,e2,~h1', '~h-->N3,e3,~h1']
decayhm = ['~h+->n1,E1,~h1', '~h+->n2,E2,~h1', '~h+->n3,E3,~h1']

# Calculate branching ratios (BR) for each process
df['brh2'] = (df['~h2->e1,E1,~h1'] / df['totalWidth'] + 
              df['~h2->e2,E2,~h1'] / df['totalWidth'] + 
              df['~h2->e3,E3,~h1'] / df['totalWidth'])

df['brhp'] = (df['~h-->N1,e1,~h1'] / df['totalWidth'] + 
              df['~h-->N2,e2,~h1'] / df['totalWidth'] + 
              df['~h-->N3,e3,~h1'] / df['totalWidth'])

df['brhm'] = (df['~h+->n1,E1,~h1'] / df['totalWidth'] + 
              df['~h+->n2,E2,~h1'] / df['totalWidth'] + 
              df['~h+->n3,E3,~h1'] / df['totalWidth'])

'''
# H2 contribution (where the decay widths are not NaN)
mask_h2 = df[decayh2].notna().all(axis=1)
df['brh2_weighted'] = np.where(mask_h2, df['brh2'] * df['CrossSection_fb'], np.nan)

# HP contribution (where the decay widths are not NaN)
mask_hp = df[decayhp].notna().all(axis=1)
df['brhp_weighted'] = np.where(mask_hp, df['brhp'] * df['CrossSection_fb'], np.nan)

# HM contribution (where the decay widths are not NaN)
mask_hm = df[decayhm].notna().all(axis=1)
df['brhm_weighted'] = np.where(mask_hm, df['brhm'] * df['brhp'] * df['CrossSection_fb'], np.nan)
'''

df['cs_total'] = df[['CrossSection_fb', 'brh2', 'brhp', 'brhm']].prod(axis=1)

#brhphm = calcBr(df, widthshphm, 'widthsTotalhphm', 'hphm')


df.to_csv(output_file, sep="\t", index=False, float_format="%.6e", na_rep="NaN")

print(f"Data written to {output_file}")

"""
df['total_CS_hphm'] = df['CS_hphm'] * df['Br_hphm']

#total cross section:
df['total_CS'] = (2 * df['total_CS_hphm']) + df['total_CS_h1h2']

def calcBr(df, widths, total_width, process):
    total_br = 0
    for width in widths:
        br_col = f'BR_{width}'
        df[br_col] = df[width] / df[total_width]
        total_br += df[br_col]
    df[f'totalBr_{process}'] = total_br"""