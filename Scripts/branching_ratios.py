import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 

path = "Work/CalcHEP_Scan_W/scan2.dat"

df = pd.read_csv(path, sep=r'\s+', low_memory= False)

params_dict = {
    'MD1': r'$m_{h1}$',
    'MD2': r'$m_{h2}$',
    'MDP': r'$m_{h_{\pm}}$',
    'DM2' : r'$\Delta m_1$', #mass diff mh2 - mh1
    'DM3' : r'$\Delta m_+$', #mass diff mh+ - mh2
    "Br(h1->e nu)" 	: r'Br($h_1 \rightarrow e \nu h_1$)',
    "Br(h1->mu nu)" : r'Br($h_1 \rightarrow \mu \nu_{\mu} h_1$)',
    "Br(h1->tau nu)" : r'Br($h_1 \rightarrow \tau \nu_{\tau} h_1$)',	 
    "Br(h2->e nu)" 	: r'Br($h_2 \rightarrow e \nu h_1$)',
    "Br(h2->mu nu)" : r'Br($h_2 \rightarrow \mu \nu_{\mu} h_1$)',
    "Br(h2->tau nu)" : r'Br($h_2 \rightarrow \tau \nu_{\tau} h_1$)'
}

def plotfig(df1, df2, xlog = False, ylog = False ):
    label1 = label_dict.get(df1, df1)
    label2 = label_dict.get(df2, df2)
    
    plt.figure(figsize=(8, 6))
    if xlog == True:
        plt.xscale('log')
    
    if ylog == True:
        plt.yscale('log')
        
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.xlabel(label1, fontsize=12)
    plt.ylabel(label2, fontsize=12)
    
    plt.scatter(dataframe[df1], dataframe[df2])

cutBr = ((df["Br(h1->e nu)"] > 0) & (df["Br(h1->mu nu)"] > 0) & (df["Br(h1->tau nu)"] > 0) & 
        (df["Br(h2->e nu)"] > 0) &  (df["Br(h2->mu nu)"] > 0) & (df["Br(h2->tau nu)"] > 0)
        )

dataframe_cut = df[cutBr]

plotfig('MD1', 'DM2')