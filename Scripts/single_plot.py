import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import itertools
import pandas as pd
import cuts

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


def plot_cut(df, xvar, yvar, scalevar, xlog=True, ylog=True,
                   label_dict=params_dict, cuts_func=cuts.cuts):

    # Apply cuts
    filtered_df = cuts_func(df, cut1 = True,cut2  = True, cut3  = True, cut4  = True, cut5  = True, cut6  = True, cut7  = True, cut8  = True)

    # Axis labels
    label1 = label_dict.get(xvar, xvar)
    label2 = label_dict.get(yvar, yvar)
    label3 = label_dict.get(scalevar, scalevar)

    # Start plotting
    plt.figure(figsize=(8, 6))
    sc = plt.scatter(
        filtered_df[xvar], filtered_df[yvar],
        c=filtered_df[scalevar],
        s=1.5, cmap='plasma',
        norm=LogNorm(vmin=filtered_df[scalevar].min(), vmax=filtered_df[scalevar].max()),
        rasterized=True
    )

    plt.xlabel(label1, fontsize=18)
    plt.ylabel(label2, fontsize=18)

    if xlog:
        plt.xscale('log')
    if ylog:
        plt.yscale('log')

    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()

    plt.xlabel(label1, fontsize=18)
    plt.ylabel(label2, fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.colorbar(sc, label=label3).ax.tick_params(labelsize=18)
    
    if xvar == 'MD1':
        plt.xlim = (40, 70)
        
    if yvar == 'MD1':
        plt.ylim = (40, 70)

    print(f'figure {yvar}_against_{xvar} saved')
    plt.savefig(f'single_plots/{yvar}_against_{xvar}_{scalevar}.pdf',  bbox_inches='tight', dpi=80)
    plt.close("all")



elements = ['MD1', 'MD2', 'MDP', 'DMP', 'DM2', 'DM3', 'l345']
pairs = [list(p) for p in itertools.permutations(elements, 2)]

for x, y in pairs:
    plot_cut(df, x, y, 'Omegah2', xlog = False, ylog = False)