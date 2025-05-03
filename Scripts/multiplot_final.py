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
    'DMP' : r'$\Delta m_+$', #mass diff mh+ - mh2
<<<<<<< HEAD
    'DM3' : r'$\Delta m_3$' #mass diff mh2 - mh+
=======
    'DM3' : r'$\Delta m_+$' #mass diff mh2 - mh+
>>>>>>> 744fba5 (batch_file reorganise, new dataset)
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
    axes = axes.flatten()

    cut_kwargs = {}  # Cumulative cuts to apply
    ylims_row2 = None  # Initialize ylims for row 2

    for i, (cut_flag, cut_label) in enumerate(cuts_dict.items()):
        ax = axes[i]
        
        # Add current cut to the argument list
        cut_kwargs[cut_flag] = True
        
        # Apply cumulative cuts
        filtered_df = cuts_func(df.copy(), **cut_kwargs)

        if xlog:
            ax.set_xscale('log')
            
        if ylog:
            if yvar != 'l345':
                ax.set_yscale('log')

        label1 = label_dict.get(xvar, xvar)
        label2 = label_dict.get(yvar, yvar)
        label3 = label_dict.get(scalevar, scalevar)

        # Main scatter plot
        if scale:
            sc = ax.scatter(
                filtered_df[xvar], filtered_df[yvar],
                c=filtered_df[scalevar],
                s=2, cmap='plasma',
                norm=LogNorm(vmin=1e-5, vmax=filtered_df[scalevar].max()),
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
        
        """
        if i == 0:
            ydata = filtered_df[yvar]
            ydata = ydata[np.isfinite(ydata)]  # Clean NaNs and infs
            if ylog and yvar != 'l345':
                ydata = ydata[ydata > 0]  # Remove zeros or negatives for log scale
            ylims_row1 = (ydata.min(), ydata.max())
            print("Captured manual y-limits for row 2:", ylims_row2)
        # Capture y-limits for the second row
        
        if i == 4:
            ydata = filtered_df[yvar]
            ydata = ydata[np.isfinite(ydata)]  # Clean NaNs and infs
            if ylog and yvar != 'l345':
                ydata = ydata[ydata > 0]  # Remove zeros or negatives for log scale
            ylims_row2 = (ydata.min(), ydata.max())
            print("Captured manual y-limits for row 2:", ylims_row2)
        
        if  1<= i <= 3:
            # Apply captured y-limits to the second row plots
            if ylims_row2:
                axes[i].set_ylim(ylims_row1)
        
        if 5 <= i <= 7:
            # Apply captured y-limits to the second row plots
            if ylims_row2:
                axes[i].set_ylim(ylims_row2)
        
        """
    # Single colorbar on the far right (works with constrained_layout)
    if scale:
        cbar = fig.colorbar(sc, ax=axes.ravel().tolist(), location='right')
        cbar.set_label(label3, fontsize=25)
        cbar.ax.tick_params(labelsize=20)

    # Title for the full figure
    fig.suptitle(f'Cumulative Cuts on {label2} against {label1}, scaled by {label3}', fontsize=28)
    
    fig.savefig(f"big_plots/{yvar}_against_{xvar}_{scalevar}.pdf", bbox_inches='tight', dpi=150)
    print(f'figure {yvar}_against_{xvar} saved')        
    #plt.close('all')           
    #plt.show()
 
#plot_cuts_grid(df, xvar='MD1', yvar='l345', scalevar= 'Omegah2', ylog =  False)

elements = ['MD1', 'MD2', 'MDP', 'DMP', 'DM2', 'DM3', 'l345']
pairs = [list(p) for p in itertools.permutations(elements, 2)]

for x, y in pairs:
    plot_cuts_grid(df, xvar=x, yvar=y, scalevar= 'Omegah2')

#plot_cuts_grid(df, xvar='MD1', yvar='l345',ylog = False, scalevar= 'Omegah2')