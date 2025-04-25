import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cuts
from matplotlib.colors import LogNorm

file_path = 'scans\\5-D_scans\\combined.dat'

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
    'cut1' : 'Vacuum stability',
    'cut2' : '+ LEP',
    'cut3' : '+ EWPT',
    'cut4' : '+ Relic Density',
    'cut5' : '+ DM DD',
    'cut6' : '+ CMB',
    'cut7' : '+ Branching Ratio',
    'cut8' : '+ SI'
}


def plotfig(dataframe, df1, df2, omegah2bar = False, xlog = True, ylog = True, savefig = False, ext = None, label_dict = params_dict):
    
    dataframe
    label1 = label_dict.get(df1, df1)
    label2 = label_dict.get(df2, df2)
    
    plt.figure(figsize=(8, 6))
    if xlog == True:
        plt.xscale('log')
    
    if ylog == True:
        plt.yscale('log')
        
    if omegah2bar == True:
        sc = plt.scatter(dataframe[df1], dataframe[df2], c=dataframe['Omegah2'], rasterized=True, s=1,
                        cmap='plasma', norm=LogNorm(vmin=10e-6, vmax=dataframe['Omegah2'].max()))

        plt.title(f'Plot of {label1} against {label2}, coloured by $\\Omega h_2$')
        
        cbar = plt.colorbar(sc)
        cbar.set_label('$\\Omega h_2$', fontsize=15)
    
    else:
        sc = plt.scatter(dataframe[df1], dataframe[df2], c = 'red', s = 1)
        plt.title(f'Plot of {label1} against {label2}')

    if omegah2bar & savefig:
        plt.savefig(f"plots/plot_{df1}_{df2}(no_grad){ext}.pdf", format='pdf') #pdf is full quality
    
    if omegah2bar == False & savefig == True:
        plt.savefig(f"plots/plot_{df1}_{df2}{ext}.pdf", format='pdf') #pdf is full quality

    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.xlabel(label1, fontsize=15)
    plt.ylabel(label2, fontsize=15)
    
def constraintplot(df1, df2, label_dict = params_dict):
    label1 = label_dict.get(df1, df1)
    label2 = label_dict.get(df2, df2)
    
    plt.figure(figsize=(8, 6))
    
    lz5tmedian_df = dataframe[dataframe['expName'] == 'LZ5Tmedian']

    notappl_df = dataframe[dataframe['expName'] == 'NotAppl']

    xenon1t_df = dataframe[dataframe['expName'] == 'XENON1T_2018']
    
    plt.scatter(notappl_df[df1], notappl_df[df2], c = 'blue', s = 0.8, label = 'NotAppl')
    plt.scatter(lz5tmedian_df[df1], lz5tmedian_df[df2], c = 'red', s = 0.8, label = 'LZ5Tmedian')
    plt.scatter(xenon1t_df[df1], xenon1t_df[df2], c = 'green', s = 0.8, label = 'XENON1T')
    #plt.scatter(notappl_df[df1], notappl_df[df2], c = 'blue', s = 0.8, label = 'NotAppl')

    #plt.ylim(-2, 2)
    plt.xscale('log')
    plt.xlabel(label1)
    plt.ylabel(label2)
    plt.legend()
    #plt.savefig('run/run_Dec2/allowedParamsExps1.pdf')
    plt.show()


"""
def multi_plot_old(df, x_axis, y_axis, ylog = None, xlog = None, savefig = True):
    df_f = cuts(df, cut1 = True)
    plotfig(df_f, x_axis, y_axis, omegah2bar = True, xlog = xlog, ylog = ylog, savefig = savefig)

    df_f = cuts(df, cut1 = True, cut2= True)
    plotfig(df_f, x_axis, y_axis, omegah2bar = True, xlog = xlog, ylog = ylog, savefig = savefig)

    df_f = cuts(df, cut1 = True, cut2= True, cut3 = True)
    plotfig(df_f, x_axis, y_axis, omegah2bar = True, xlog = xlog, ylog = ylog, savefig = savefig)

    #df_f = cuts(df, cut1 = True, cut2= True, cut3_strict = True)
    
    df_f = cuts(df, cut1 = True, cut2= True, cut3 = True, cut4 = True)
    plotfig(df_f, x_axis, y_axis, omegah2bar = True, xlog = xlog, ylog = ylog, savefig = savefig)

    df_f = cuts(df, cut1 = True, cut2= True, cut3 = True, cut4 = True, cut5 = True)
    plotfig(df_f, x_axis, y_axis, omegah2bar = True, xlog = xlog, ylog = ylog, savefig = savefig)
   
    plt.show()
    
def multi_plot(df, x_axis, y_axis, xlog=None, ylog=None, omegah2bar=True, cut_flags=None):
    '
    Plot data with different combinations of cuts dynamically.

    Parameters:
    - df: DataFrame
    - x_axis: column name for x-axis
    - y_axis: column name for y-axis
    - ylog_values: List of booleans to toggle y-log scale for each plot
    - omegah2bar: Whether to color by Omega h2
    - cut_flags: A list of dictionaries with cut booleans for each plot
    '
    if cut_flags is None:
        # Default list of cuts if none are provided
        cut_flags = [
            {'cut1': True},  # Cut 1 only
            {'cut1': True, 'cut2': True},  # Cut 1 + Cut 2
            {'cut1': True, 'cut2': True, 'cut3': True},  # Cut 1 + Cut 2 + Cut 3
            {'cut1': True, 'cut2': True, 'cut3': True, 'cut4': True},  # Add Cut 4
            {'cut1': True, 'cut2': True, 'cut3': True, 'cut4': True, 'cut5': True}  # Add Cut 5
        ]

    # Loop over each set of cuts and corresponding ylog value
    for i, flags in enumerate(cut_flags):
        print(f"Applying cuts: {flags}")
        df_f = cuts(df, **flags)  # Unpack cuts dynamically
        plotfig(df_f, x_axis, y_axis, omegah2bar=omegah2bar, xlog=xlog, ylog=ylog, savefig=False)

    plt.show()
    
multi_plot_old(df, 'MD1', 'l345', xlog = True, ylog = False)

"""
