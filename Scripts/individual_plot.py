import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cuts
from matplotlib.colors import LogNorm
import itertools

plt.close('all')
plt.rcParams['font.serif'] = ['Times New Roman'] 

file_path = 'scans/5-D_scans/combined.dat'

df = pd.read_csv(file_path, sep=r'\s+', low_memory=False)

params_dict = {
    'MD1': r'$m_{h1}$',
    'MD2': r'$m_{h2}$',
    'MDP': r'$m_{h_{\pm}}$',
    'Omegah2': r'$\Omega h_2$',
    'l345': r'$\lambda_{345}$',
    'DM2': r'$\Delta m_1$',  # mass diff mh2 - mh1
    'DMP': r'$\Delta m_+$',  # mass diff mh+ - mh2
    'DM3': r'$\Delta m_+$'   # mass diff mh2 - mh+
}

df_cut = cuts.cuts(df, cut1=True, cut2=True, cut3=True, cut4=True, 
                   cut5=True, cut6=True, cut7=True, cut8=True)

# Ensure 'DMP' and 'MD1' columns exist and avoid errors
if 'DMP' in df_cut.columns and 'MD1' in df_cut.columns:
    # Create the new column with the ratio DMP/MD1
    df_cut.loc[:, 'values_divide'] = df_cut['DMP'] / df_cut['MD1']
else:
    print("Error: 'DMP' or 'MD1' column is not present in the dataframe")

# Ensure there are no NaN values in the new column
df_cut.dropna(subset=['values_divide'], inplace=True)

# Plotting
plt.scatter(df_cut['MD1'], df_cut['values_divide'], c='red', label=r'$\frac{DMP}{MD1}$')
plt.axhline(y=0.1, color='blue', linestyle='--', label=r'Expected ratio $\frac{DMP}{MD1} = 0.1$')

# Optional: plot DMP vs MD1 as well to compare
plt.scatter(df_cut['MD1'], df_cut['DMP'], c='black', label='DMP vs MD1')

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$m_{h1}$')
plt.ylabel(r'$DMP / m_{h1}$')
plt.title('Scatter plot of ratio DMP/MD1 with expected value of 0.1')
plt.legend()
plt.show()