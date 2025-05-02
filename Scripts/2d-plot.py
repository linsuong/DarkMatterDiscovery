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


plt.scatter(df_cut['MD1'], df_cut['l345'], s = 0.5, c = 'black')
#plt.xlim(50, 100)
plt.show()