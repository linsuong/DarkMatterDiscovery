import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
import re

def sanitize_filename(name):
    return re.sub(r'[^a-zA-Z0-9_\-]', '_', name)

path = "Work/CalcHEP_Scan_Z_MD1_100GeV/scan2.dat"

df = pd.read_csv(path, sep=r'\s+', low_memory= False)

params_dict = {
    'MD1': r'$m_{h1}$',
    'MD2': r'$m_{h2}$',
    'MDP': r'$m_{h_{\pm}}$',
    'DM2' : r'$\Delta m^2$', #mass diff mh2 - mh1
    'DMP' : r'$\Delta m^+$', #mass diff mh+ - mh1
    'DM3' : r'$\Delta m^3$',
    
    #for decay via Z boson:
    "Br(h2->e+e-h1)" 	: r'Br($h_2 \rightarrow e^- e^+ h_1$)',
    "Br(h2->mu+mu-h1)" : r'Br($h_2 \rightarrow \mu^- \mu^+ h_1$)',
    "Br(h2->tau+tau-h1)" : r'Br($h_2 \rightarrow \tau^-  \tau^+ h_1$)',	
    "Br(h2->n+n-h1)" : r'Br($h_2 \rightarrow n n h_1$)',
    
    #for decay via W- boson:
    "Br(h-->e-nu)" : r'Br($h_- \rightarrow e^- \bar{\nu}_e h_1$)',
    "Br(h-->mu-nu)" : r'Br($h_- \rightarrow \mu^- \bar{\nu}_{\mu} h_1$)',
    "Br(h-->tau-nu)" : r'Br($h_- \rightarrow \tau^- \bar{\nu}_{\tau} h_1$)',
    
    #for decay via W+ boson:
    "Br(h+->e-nu)" : r'Br($h_+ \rightarrow e^+ \nu_e h_1$)',
    "Br(h+->mu-nu)" : r'Br($h_+ \rightarrow \mu^+ \nu_{\mu} h_1$)',
    "Br(h+->tau-nu)" : r'Br($h_+ \rightarrow \tau^+ \nu_{\tau} h_1$)'
}


if path == "Work/CalcHEP_Scan_Z/scan2.dat" or path == "Work/CalcHEP_Scan_Z_MD1_100GeV/scan2.dat":
    
    branching_ratios = ['Br(h2->e+e-h1)', 'Br(h2->mu+mu-h1)', 'Br(h2->tau+tau-h1)', 'Br(h2->n+n-h1)']

if path == "Work/CalcHEP_Scan_W-/scan2.dat":
    branching_ratios = ['Br(h-->e-nu)', 'Br(h-->mu-nu)', 'Br(h-->tau-nu)']
    
if path == "Work/CalcHEP_Scan_W+/scan2.dat":    
    branching_ratios = ['Br(h+->e+nu)', 'Br(h+->mu+nu)', 'Br(h+->tau+nu)']

def plotfig(df1, df2, df3, xlog=False, ylog=False, label_dict=params_dict):
    label1 = label_dict.get(df1, df1)
    label2 = label_dict.get(df2, df2)
    label3 = label_dict.get(df3, df3)
    label5 = label_dict.get('MD1', 'MD1')
    
    if xlog:
        plt.xscale('log')
    if ylog:
        plt.yscale('log')

    vmin = df[df3].min()
    vmax = df[df3].max()
    
    sc = plt.scatter(df[df1], df[df2], c=df[df3], rasterized=True, s=1,
                     cmap='plasma_r', norm=Normalize(vmin=vmin, vmax=vmax))
    
    cbar = plt.colorbar(sc)
    cbar.set_label(label3, fontsize=15)
            
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.xlabel(label1, fontsize=15)
    plt.ylabel(label2, fontsize=15)
    plt.tick_params(axis='both', labelsize=15)
    plt.title(f'{label2} against {label1}, scaled by {label3} (1 to 3), {label5} = 100 GeV', fontsize=15)

    # Sanitize filename
    safe_label1 = sanitize_filename(label1)
    safe_label3 = sanitize_filename(label3)
    safe_label4 = sanitize_filename(label2)
    
    plt.savefig(f'branching_ratio_plots/1_to_3/MD1_100GeV/{safe_label4}{safe_label1}_scale_{safe_label3}_MD1_100GeV.jpg', bbox_inches='tight')
    plt.close()


for i in range(len(branching_ratios)):
    print(i)
    plotfig('DMP', branching_ratios[i], 'DM3')
    plotfig('DM3', branching_ratios[i], 'DMP')
