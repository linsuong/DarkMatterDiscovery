import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
import re

def sanitize_filename(name):
    return re.sub(r'[^a-zA-Z0-9_\-]', '_', name)

path = "Work/CalcHEP_Scan_W/scan2.dat"

df = pd.read_csv(path, sep=r'\s+', low_memory= False)

params_dict = {
    'MD1': r'$m_{h1}$',
    'MD2': r'$m_{h2}$',
    'MDP': r'$m_{h_{\pm}}$',
    'DM2' : r'$\Delta m_1$', #mass diff mh2 - mh1
    'DMP' : r'$\Delta m_+$', #mass diff mh+ - mh1
    'DM3' : r'$\Delta m_2$',
    
    #for decay via Z boson:
    "Br(h2->e+e-h1)_total" 	: r'Br($h_2 \rightarrow e^- e^+ h_1$)',
    "Br(h2->mu+mu-h1)_total" : r'Br($h_2 \rightarrow \mu^- \mu^+ h_1$)',
    "Br(h2->tau+tau-h1)_total" : r'Br($h_2 \rightarrow \tau^-  \tau^+ h_1$)',	
    "Br(h2->n+n-h1)_total" : r'Br($h_2 \rightarrow n n h_1$)',
    
    #for decay via W- boson:
    "Br(h-->e-nu)_total" : r'Br($h_- \rightarrow e^- \bar{\nu}_e h_1$)',
    "Br(h-->mu-nu)_total" : r'Br($h_- \rightarrow \mu^- \bar{\nu}_{\mu} h_1$)',
    "Br(h-->tau-nu)_total" : r'Br($h_- \rightarrow \tau^- \bar{\nu}_{\tau} h_1$)',
    
    #for decay via W+ boson:
    "Br(h+->e-nu)_total" : r'Br($h_+ \rightarrow e^+ \nu_e h_1$)',
    "Br(h+->mu-nu)_total" : r'Br($h_+ \rightarrow \mu^+ \nu_{\mu} h_1$)',
    "Br(h+->tau-nu)_total" : r'Br($h_+ \rightarrow \tau^+ \nu_{\tau} h_1$)'
}

"""
"Br(h2->Zh1)" : r'Br(h_2 \rightarrow h_1 Z)',
    
"Br(Z->e+e-)" 	: r'Br($Z \rightarrow e^- e^\+$)',
"Br(Z->mu+mu-)" : r'Br($Z \rightarrow \mu^- \mu\+$)',
"Br(Z->tau+tau-)" : r'Br($Z \rightarrow \tau^- \tau^\+$)',
"Br(Z->nn)" 	: r'Br($Z \rightarrow n n$)',

"Br(h-->W-h1)" : r'Br($h_- \rightarrow W^- h_1$)',	

"Br(W-->e-nu)" : r'Br($W^- \rightarrow e^- \nu_e$)',
"Br(W-->mu-nu)" : r'Br($W^- \rightarrow \mu^- \nu_{\mu}$)',
"Br(W-->tau-nu)" : r'Br($W^- \rightarrow \tau^- \nu_{\tau}$)'
"""

savefile_dict = {
    #for decay via Z boson:
    "Br(h2->e+e-h1)_total" 	: r'Br(h2->ee)',
    "Br(h2->mu+mu-h1)_total" : r'Br(h2->mm)',
    "Br(h2->tau+tau-h1)_total" : r'Br(h2->tt)',	
    "Br(h2->n+n-h1)_total" : r'Br(h2->nnh1)',
    
    #for decay via W- boson:
    "Br(h-->e-nu)_total" : r'Br(h-->en)',
    "Br(h-->mu-nu)_total" : r'Br(h-->mm)',
    "Br(h-->tau-nu)_total" : r'Br(h-->tn)',	
    
    #for decay via W+ boson:
    "Br(h-->e-nu)_total" : r'Br(h+->en)',
    "Br(h-->mu-nu)_total" : r'Br(h+->mm)',
    "Br(h-->tau-nu)_total" : r'Br(h+->tn)'
}

if path == "Work/CalcHEP_Scan_Z/scan2.dat":
    # Total BRs for h2 decays (direct + via Z)
    df['Br(h2->e+e-h1)_total'] = (df["Br(Z->e+e-h1)"] * df["Br(h2->Zh1)"]) + df['Br(h2->e+e-h1)']
    df['Br(h2->mu+mu-h1)_total'] = (df["Br(Z->mu+mu-h1)"] * df["Br(h2->Zh1)"]) + df['Br(h2->mu+mu-h1)']
    df['Br(h2->tau+tau-h1)_total'] = (df["Br(Z->tau+tau-h1)"] * df["Br(h2->Zh1)"]) + df['Br(h2->tau+tau-h1)']
    df['Br(h2->n+n-h1)_total'] = (df["Br(Z->n+n-h1)"] * df["Br(h2->Zh1)"]) + df['Br(h2->n+n-h1)']
    
    branching_ratios = ['Br(h2->e+e-h1)_total', 'Br(h2->mu+mu-h1)_total', 'Br(h2->tau+tau-h1)_total', 'Br(h2->n+n-h1)_total']

if path == "Work/CalcHEP_Scan_W-/scan2.dat":
    # Total BRs for h- decays (direct + via W)
    df['Br(h-->e-nu)_total'] = (df["Br(W-->e-nu)"] * df["Br(h-->W-h1)"]) + df['Br(h-->e-nu)']
    df['Br(h-->mu-nu)_total'] = (df["Br(W-->mu-nu)"] * df["Br(h-->W-h1)"]) + df['Br(h-->mu-nu)']
    df['Br(h-->tau-nu)_total'] = (df["Br(W-->tau-nu)"] * df["Br(h-->W-h1)"]) + df['Br(h-->tau-nu)']
    
    branching_ratios = ['Br(h-->e-nu)_total', 'Br(h-->mu-nu)_total', 'Br(h-->tau-nu)_total']
    
if path == "Work/CalcHEP_Scan_W-/scan2.dat":
    # Total BRs for h+ decays (direct + via W)
    df['Br(h+->e+nu)_total'] = (df["Br(W+->e+nu)"] * df["Br(h+->W+h1)"]) + df['Br(h+->e+nu)']
    df['Br(h+->mu+nu)_total'] = (df["Br(W+->mu+nu)"] * df["Br(h+->W+h1)"]) + df['Br(h+->mu+nu)']
    df['Br(h+->tau+nu)_total'] = (df["Br(W+->tau+nu)"] * df["Br(h+->W+h1)"]) + df['Br(h+->tau+nu)']
    
    branching_ratios = ['Br(h+->e+nu)_total', 'Br(h+->mu+nu)_total', 'Br(h+->tau+nu)_total']

def plotfig(df1, df2, df3, xlog=False, ylog=False, label_dict=params_dict):
    label1 = label_dict.get(df1, df1)
    label2 = label_dict.get(df2, df2)
    label3 = label_dict.get(df3, df3)
    label4 = savefile_dict.get(df2, df2)
    
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
    plt.title(f'{label2} against {label1}, scaled by {label3}', fontsize=15)

    # Sanitize filename
    safe_label1 = sanitize_filename(label1)
    safe_label3 = sanitize_filename(label3)
    safe_label4 = sanitize_filename(label4)
    
    plt.savefig(f'branching_ratio_plots/images/{safe_label4}{safe_label1}_scale_{safe_label3}.jpg', bbox_inches='tight')
    plt.close()


for i in range(len(branching_ratios)):
    print(i)
    plotfig('DMP', branching_ratios[i], 'DM3')
    plotfig('DM3', branching_ratios[i], 'DMP')
