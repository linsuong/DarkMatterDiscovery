import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cuts
from matplotlib.colors import LogNorm

file_path = 'scans/5-D_scans/combined.dat'

df = pd.read_csv(file_path, sep=r'\s+', low_memory= False)

params_dict = {
    'MD1': r'$m_{h1}$',
    'MD2': r'$m_{h2}$',
    'MDP': r'$m_{h_{\pm}}$',
    'Omegah2': r'$\Omega h_2$',
    'l345' : r'$\lambda_{345}$',
    'DM2' : r'$\Delta m_1$', #mass diff mh2 - mh1
    'DM3' : r'$\Delta m_+$' #mass diff mh+ - mh2
}

cuts_dict = {
    'cut1' : '+ Vacuum stability',
    'cut2' : '+ LEP',
    'cut3' : '+ DM Direct Detection',
    'cut4' : '+ CMB',
    'cut5' : r'+ $Br(H \rightarrow \mathrm{invisible}) < 0.145$',
    'cut6' : '+ LZ 2024',
    'cut7' : '+ EWPT',
    'cut8' : '+ Relic Density'
    
}

def plot_cuts_grid(df, xvar, yvar, scalevar, scale=True, xlog=True, ylog=True,
                   label_dict=params_dict, cuts_func=cuts.cuts, cuts_dict=cuts_dict):
    
    # Create subplots with auto-managed spacing
    fig, axes = plt.subplots(2, 4, figsize=(22, 10), sharex=True, sharey=True, constrained_layout=True)
    axes = axes.flatten()

    cut_kwargs = {}  # Cumulative cuts to apply

    for i, (cut_flag, cut_label) in enumerate(cuts_dict.items()):
        ax = axes[i]
        
        # Add current cut to the argument list
        cut_kwargs[cut_flag] = True
        
        # Apply cumulative cuts
        colorbar_range = cuts_func(df.copy(), cut1=True, cut2=True)
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

        ax.set_title(cut_label, fontsize=15)
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.set_xlabel(label1, fontsize=15)
        ax.set_xlabel(label1, fontsize=15)
        if i % 4 == 0:
            ax.set_ylabel(label2, fontsize=15)

    # Single colorbar on the far right (works with constrained_layout)
    if scale:
        cbar = fig.colorbar(sc, ax=axes.ravel().tolist(), location='right', shrink=0.8)
        cbar.set_label(label3, fontsize=15)

    # Title for the full figure
    fig.suptitle(f'Cumulative Cuts on {label2} against {label1}, scaled by {label3}', fontsize=20)

    # Show plot
    #plt.show()
    
    #fig.savefig(f"big_plots_(low_dpi)/{xvar}_against_{yvar}_{scalevar}.pdf", bbox_inches='tight', dpi = 150)
    #print('figure saved')                    
    plt.show()


pairs = [
    ['MD1', 'MD2'], ['MD1', 'MDP'], ['MD1', 'DMP'], ['MD1', 'DM2'], ['MD1', 'DM3'],
    ['MD2', 'MD1'], ['MD2', 'MDP'], ['MD2', 'DMP'], ['MD2', 'DM2'], ['MD2', 'DM3'],
    ['MDP', 'MD1'], ['MDP', 'MD2'], ['MDP', 'DMP'], ['MDP', 'DM2'], ['MDP', 'DM3'],
    ['MD1', 'MD2'], ['MDP', 'MD1'], ['DMP', 'MD1'], ['DM2', 'MD1'], ['DM3', 'MD1'],
    ['MD2', 'MD1'], ['MDP', 'MD2'], ['DMP', 'MD2'], ['DM2', 'MD2'], ['DM3', 'MD2'],
    ['MD1', 'MDP'], ['MD2', 'MDP'], ['DMP', 'MDP'], ['DM2', 'MDP'], ['DM3', 'MDP'],
    ['MD1', 'l345'], ['MD2', 'l345'], ['MDP', 'l345'],
    ['DMP', 'l345'], ['DM2', 'l345'], ['DM3', 'l345']
]
for x, y in pairs:
    plot_cuts_grid(df, xvar=x, yvar=y, scalevar= 'Omegah2')


#plot_cuts_grid(df, xvar='MD1', yvar='l345', scalevar= 'Omegah2')