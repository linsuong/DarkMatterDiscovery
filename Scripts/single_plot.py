import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import itertools
import pandas as pd
import cuts

plt.close('all')

plt.rcParams['font.serif'] = ['Times New Roman'] 

file_path = 'scan_shrink.dat'

df = pd.read_csv(file_path, sep=r'\s+', low_memory= False)

params_dict = {
    'MD1': r'$m_{h1}$',
    'MD2': r'$m_{h2}$',
    'MDP': r'$m_{h_{\pm}}$',
    'Omegah2': r'$\Omega h_2$',
    'l345' : r'$\lambda_{345}$',
    'DM2' : r'$\Delta m_1$', #mass diff mh2 - mh1
    'DMP' : r'$\Delta m_+$', #mass diff mh+ - mh2
    'DM3' : r'$\Delta m_3$' #mass diff mh2 - mh+
}


def plot_cut(df, xvar, yvar, scalevar, xlog=True, ylog=True, scale = True,
                   limits = False, limits2 = False, label_dict=params_dict, cuts_func=cuts.cuts):

    # Apply cuts
    filtered_df = cuts_func(df, cut1 = True,cut2  = True, cut3  = True, cut4  = True, 
                            cut5  = True, cut6  = True, cut7  = True, cut8  = True)

    # Axis labels
    label1 = label_dict.get(xvar, xvar)
    label2 = label_dict.get(yvar, yvar)
    label3 = label_dict.get(scalevar, scalevar)

    # Start plotting
    plt.figure(figsize=(8, 6))
    if scale:
        sc = plt.scatter(
            filtered_df[xvar], filtered_df[yvar],
            c=filtered_df[scalevar],
            s=10, cmap='plasma',
            norm=LogNorm(vmin=filtered_df[scalevar].min(), vmax=filtered_df[scalevar].max()),
            rasterized=True
        )
        
        cbar = plt.colorbar(sc)
        cbar.set_label(label3, fontsize=20)
        cbar.ax.tick_params(labelsize=20)

    else: 
        sc = plt.scatter(
            filtered_df[xvar], filtered_df[yvar], 
            s = 10, c = 'red', rasterized=True)

    if xlog:
        plt.xscale('log')
    if ylog:
        plt.yscale('log')

    plt.grid(which='major', linestyle='--', linewidth=0.8, alpha=0.5)
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.4)
    
    plt.tight_layout()

    plt.xlabel(label1, fontsize=25)
    plt.ylabel(label2, fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    
    if limits:
        if xvar in ['MD1', 'MD2']:
            plt.xlim(0, 300)

        if yvar in ['MD1', 'MD2']:
            plt.ylim(0, 300)
    
        if yvar in ['l345'] and xvar in ['MD1']: #noted
            plt.ylim(-0.0001, 0.0001)
            plt.xlim(60, 80)
            
        if yvar in ['DMP'] and xvar in ['MD1']: #noted
            plt.ylim(bottom = 8)
            plt.xlim(60, 80)
    
        if yvar in ['DM3'] and xvar in ['MD2']: #noted
            plt.xlim(400, 1000)
            plt.ylim(-12, 2)
        
        if yvar in ['MD2'] and xvar in ['MD1']: #noted
            #plt.ylim(bottom = 400)
            plt.xlim(40, 100)  
        
        #plt.show()
<<<<<<< HEAD
        plt.savefig(f'single_plots_(shrink)/{yvar}_against_{xvar}_{scalevar}.pdf',  bbox_inches='tight', dpi=80)
    
    else:   
=======
        plt.savefig(f'single_plots_(small_range)/{yvar}_against_{xvar}_{scalevar}.pdf',  bbox_inches='tight', dpi=80)

    if limits2:
        if xvar in ['MD1', 'MD2']:
            plt.xlim(left = 425)

        if yvar in ['MD1', 'MD2']:
            plt.ylim(bottom = 425)
            
        if yvar in ['MD1'] and xvar in ['DM2']:
            plt.xlim(-5, 5)
            
        if yvar in ['MD1'] and xvar in ['DMP']:
            plt.xlim(-5, 10)
            
        plt.savefig(f'single_plots_(big_range)/{yvar}_against_{xvar}_{scalevar}.pdf',  bbox_inches='tight', dpi=80)
            
    if limits == False and limits2 == False:   
>>>>>>> 27a133560f1962cfaf939c5978d5121ea7f58624
        #plt.show() 
        plt.savefig(f'single_plots/{yvar}_against_{xvar}_{scalevar}.pdf',  bbox_inches='tight', dpi=80)
    
    
    print(f'figure {yvar}_against_{xvar} saved')
    plt.close("all")

#plot_cut(df, 'MD1', 'MD2', 'Omegah2', xlog = False, ylog = False, scale = True, limits = False)


elements = ['MD1', 'MD2', 'MDP', 'DMP', 'DM2', 'DM3', 'l345']
pairs = [list(p) for p in itertools.permutations(elements, 2)]

for x, y in pairs:
    #plot_cut(df, x, y, 'Omegah2', xlog = False, ylog = False, scale = False, limits = True)
    plot_cut(df, x, y, 'Omegah2', xlog = False, ylog = False, scale = False, limits2 = True)
