import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cuts
from matplotlib.colors import LogNorm
import itertools

plt.close('all')

plt.rcParams['font.serif'] = ['Times New Roman'] 

file_path = 'scans/5-D_scans/combined.dat'

df = pd.read_csv(file_path, sep=r'\s+', low_memory= False)

params_dict = {
    'MD1': r'$m_{h1}$',
    'MD2': r'$m_{h2}$',
    'MDP': r'$m_{h_{\pm}}$',
    'Omegah2': r'$\Omega h_2$',
    'l345' : r'$\lambda_{345}$',
    'DM2' : r'$\Delta m_1$', #mass diff mh2 - mh1
    'DMP' : r'$\Delta m_+$' #mass diff mh+ - mh2
}

cuts_dict = {
    'cut1' : '+ Theory Constraints',
    'cut2' : '+ LEP',
    'cut3' : '+ DM Direct Detection',
    'cut4' : '+ CMB',
    'cut5' : r'+ $Br(H \rightarrow \mathrm{inv}) < 0.145$',
    'cut6' : '+ LZ 2024',
    'cut7' : '+ EWPT',
    'cut8' : '+ Relic Density'
    
}

def plot_cuts_grid(df, xvar, yvar, scalevar, scale=True, xlog=True, ylog=True,
                   label_dict=params_dict, cuts_func=cuts.cuts, cuts_dict=cuts_dict):
    
    # Create subplots with manual spacing control
    fig, axes = plt.subplots(2, 4, figsize=(16,8), sharex=True, sharey=True, constrained_layout=True)
    #plt.subplots_adjust(hspace=0.3)
    axes = axes.flatten()

    cut_kwargs = {}  # Cumulative cuts to apply

    for i, (cut_flag, cut_label) in enumerate(cuts_dict.items()):
        ax = axes[i]
        
        # Add current cut to the argument list
        cut_kwargs[cut_flag] = True
        
        # Apply cumulative cuts
        colorbar_range = cuts_func(df.copy(), cut1=True)
        filtered_df = cuts_func(df.copy(), **cut_kwargs)

        if xlog:
            ax.set_xscale('log')
        if ylog:
            ax.set_yscale('log')

        label1 = label_dict.get(xvar, xvar)
        label2 = label_dict.get(yvar, yvar)
        label3 = label_dict.get(scalevar, scalevar)

        # Main scatter plot
        if scale == True:
            sc = ax.scatter(
                filtered_df[xvar], filtered_df[yvar],
                c=filtered_df[scalevar],
                s=1, cmap='plasma',
                norm=LogNorm(vmin=1e-5, vmax=colorbar_range[scalevar].max()),
                rasterized=True
            )
        else:
            sc = ax.scatter(filtered_df[xvar], filtered_df[yvar], color='red', s=1)

        ax.set_title(cut_label, fontsize=25)
        ax.grid(True, linestyle='--', alpha=0.5)
        if i >= 4:
            ax.set_xlabel(label1, fontsize=25)
        else:
            ax.set_xlabel("")
            ax.tick_params(labelbottom=False)
        if i % 4 == 0:
            ax.set_ylabel(label2, fontsize=25)
        ax.tick_params(axis='both', labelsize=20)
        
        plt.ylim(-1.5, 1.5)

    # Single colorbar on the far right (works with constrained_layout)
    if scale:
        cbar = fig.colorbar(sc, ax=axes.ravel().tolist(), location='right')
        cbar.set_label(label3, fontsize=25)
        cbar.ax.tick_params(labelsize=20)

    # Title for the full figure
    fig.suptitle(f'Cumulative Cuts on {label2} against {label1}, scaled by {label3}', fontsize=28)
    
    #fig.savefig(f"big_plots/{yvar}_against_{xvar}_{scalevar}.pdf", bbox_inches='tight', dpi=150)
    print(f'figure {yvar}_against_{xvar} saved')        
    #plt.close('all')           
    plt.show()

"""
pairs = [
    ['MD1', 'MD2'], ['MD1', 'MDP'], ['MD1', 'DMP'], ['MD1', 'DM2'], ['MD1', 'DM3'],
    ['MD2', 'MD1'], ['MD2', 'MDP'], ['MD2', 'DMP'], ['MD2', 'DM2'], ['MD2', 'DM3'],
    ['MDP', 'MD1'], ['MDP', 'MD2'], ['MDP', 'DMP'], ['MDP', 'DM2'], ['MDP', 'DM3'],
    ['MD1', 'MD2'], ['MDP', 'MD1'], ['DMP', 'MD1'], ['DM2', 'MD1'], ['DM3', 'MD1'],
    ['MD2', 'MD1'], ['MDP', 'MD2'], ['DMP', 'MD2'], ['DM2', 'MD2'], ['DM3', 'MD2'],
    ['MD1', 'MDP'], ['MD2', 'MDP'], ['DMP', 'MDP'], ['DM2', 'MDP'], ['DM3', 'MDP'],
    ['MD1', 'l345'], ['MD2', 'l345'], ['MDP', 'l345'],
    ['DMP', 'l345'], ['DM2', 'l345'], ['DM3', 'l345']
]"""

elements = ['MD1', 'MD2', 'MDP', 'DMP', 'DM2', 'DM3', 'l345']
pairs = [list(p) for p in itertools.permutations(elements, 2)]

#for x, y in pairs:
    #plot_cuts_grid(df, xvar=x, yvar=y, scalevar= 'Omegah2')

plot_cuts_grid(df, xvar='MD1', yvar='l345',ylog = False, scalevar= 'Omegah2')