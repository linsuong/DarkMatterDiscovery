import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
import re

path = "Work/CalcHEP_Scan_Z_MD1_100GeV/scan2.dat"

df = pd.read_csv(path, sep=r'\s+', low_memory= False)

#df['Br(h2->e+e-h1)_total'] = (df["Br(Z->e+e-)"] * df["Br(h2->Zh1)"]) + df['Br(h2->e+e-h1)']

df['Br(h2->e+e-h1)_total'] = (df["Br(Z->e+e-)"] * df["Br(h2->Zh1)"]) + df['Br(h2->e+e-h1)']

cutBr = df['Br(h2->e+e-h1)_total'] < 0.01

df_cut = df[cutBr]

print(df_cut)
print(df_cut['MD1'], df_cut['DM3'], df_cut['DMP'], df_cut['Br(h2->e+e-h1)'])

#plt.scatter(df['DMP'], df['Br(h2->e+e-h1)_total'],)
#plt.show()