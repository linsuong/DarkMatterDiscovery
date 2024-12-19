import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

file_path = '/Users/linusong/Repositories/DarkMatterDiscovery/scans/5-D_scans/combined.dat.gz'

df = pd.read_csv(file_path, sep=r'\s+', compression= 'gzip', low_memory= False)

params_dict = {
    'MD1': r'$m_{h1}$',
    'MD2': r'$m_{h2}$',
    'MDP': r'$m_{h_{\pm}}$',
    'Omegah2': r'$\Omega h_2$',
    'l345' : r'$\lambda_{345}$',
    'DM2' : r'$\Delta m_1$', #mass diff mh2 - mh1
    'DM3' : r'$\Delta m_+$' #mass diff mh+ - mh2
}

def calculate_ST(dataframe):
    """
    Add S and T parameter columns to the dataframe based on mass differences.
    """
    v = 246  # Higgs vev in GeV

    # Calculate S and T based on mass splittings
    x1 = dataframe['MD1'] / dataframe['MDP']
    x2 = dataframe['MD2'] / dataframe['MDP']
    log_term = np.log(dataframe['MDP'] / dataframe['MD1'])
    
    # Approximate expressions for S and T (from paper)
    dataframe['S'] = (1 / (72 * np.pi)) * ((x2**2 - x1**2) * log_term)
    dataframe['T'] = (1 / (32 * np.pi**2 * v**2)) * (dataframe['MDP']**2 - dataframe['MD2']**2)
    
    return dataframe

def chi_squared(S, T, S0=0.06, T0=0.1, sigma_S=0.09, sigma_T=0.07, rho=0.91):
    """
    Compute the chi-squared value for the S and T parameters.
    """
    # Covariance matrix and its inverse
    cov = np.array([[sigma_S**2, rho * sigma_S * sigma_T],
                    [rho * sigma_S * sigma_T, sigma_T**2]])
    inv_cov = np.linalg.inv(cov)

    # Delta vector
    delta = np.array([S - S0, T - T0])
    
    # Compute chi-squared
    chi2 = delta.T @ inv_cov @ delta
    return chi2

def cuts(dataframe, cut1=False, cut2=False, cut3=False, cut3_strict = False, cut4=False, cut5=False):
    """
    Applies constraints to the dataframe and returns a filtered dataframe.
    """
    if cut3 == True and cut3_strict == True:
        raise Exception('please choose for a strict or relaxed bound on relic density')
    MW = 80.377
    MZ = 91.19
    l1 = 0.129
    v = 246

    dataframe['MDP'] = dataframe['DMP'] + dataframe['MD1']
    dataframe['DM2'] = dataframe['DM3'] + dataframe['DMP']
    dataframe['MD2'] = dataframe['DM2'] + dataframe['MD1']
    
    #columns for vacuum stability 
    dataframe['R'] = dataframe['l345'] /(2 * np.sqrt(l1))
    
    #columns for LEP constraints:
    dataframe['MD1+MD2'] = dataframe['MD1'] + dataframe['MD2']
    dataframe['MD1+MDP'] = dataframe['MD1'] + dataframe['MDP']
    dataframe['MD2+MDP'] = dataframe['MD2'] + dataframe['MDP']

    cutl345 = dataframe['l345'] > -np.inf
    cutMD1 = dataframe['MD1'] > -np.inf
    cutMD2 = dataframe['MD2'] > -np.inf
    cutDM2 = dataframe['DM2'] > -np.inf
    cutOM = dataframe['Omegah2'] > -np.inf
    cutDD = dataframe['PvalDD'] > -np.inf
    cutCMB = dataframe['CMB_ID'] > -np.inf
    
    cutMD1MD2 = True
    cutMD2MDP = True
    cutMDP = True
    cutEWPT = True

    if cut1:
        #vs1: mh1^2 > 0 for |R| < 1
        vs1 = (dataframe['R'].abs() < 1) & (dataframe['MD1'] ** 2 > 0)
        
        #vs2: mh1^2 > threshold for R > 1
        threshold = (dataframe['R'] - 1) * np.sqrt(l1) * (v ** 2) 
        vs2 = (dataframe['MD1']**2 > threshold) & (dataframe['R'] > 1)
        
        cutMD1_1 = vs1 | vs2
        
        cutl345 = (
            (dataframe['l345'] < 2 * (((df['MD1'] ** 2)/(v ** 2)) + np.sqrt(l1))) & 
            (dataframe['l345'] < ((16/3) * np.pi) - l1) & 
            (dataframe['l345'] < 4 * np.pi - ((3/2) * l1)) &
            (dataframe['l345'] > -1.47) &
            (dataframe['l345'] < 4 * np.pi)
        )
        
    if cut2:
        # LEP constraints
        cutMD1MD2 = (dataframe['MD1'] + dataframe['MD2']) > MZ
        cutMD1MDP = (dataframe['MD1'] + dataframe['MDP']) > MW
        cutMD2MDP = (dataframe['MD2'] + dataframe['MDP']) > MW
        cutMDP = (2 * dataframe['MDP'] > MZ) & (dataframe['MDP'] > 70)

        # EWPT constraints
        x1 = dataframe['MD1'] / dataframe['MDP']
        x2 = dataframe['MD2'] / dataframe['MDP']
        log_term = np.log(dataframe['MDP'] / dataframe['MD1'])
        
        # Approximate expressions for S and T (from paper)
        dataframe['S'] = (1 / (72 * np.pi)) * ((x2**2 - x1**2) * log_term)
        dataframe['T'] = (1 / (32 * np.pi**2 * v**2)) * (dataframe['MDP']**2 - dataframe['MD2']**2)

        chi2_threshold = 5.99  # 95% confidence level
        dataframe['chi2'] = dataframe.apply(lambda row: chi_squared(row['S'], row['T']), axis=1)
        cutEWPT = dataframe['chi2'] < chi2_threshold

        #cutMD1 = dataframe['MD1'] > 80
        #cutMD2 = dataframe['MD2'] > 100
        #cutDM2 = dataframe['DM2'] < 8

    if cut3:
        cutOM = dataframe['Omegah2'] < 0.12024
        #cutOM = (df['Omegah2'] > 0.10) & (df['Omegah2'] < 0.12024) #strict bound of Omegah2
        
    if cut3_strict:
        cutOM = (df['Omegah2'] > 0.10737) & (df['Omegah2'] < 0.13123) #strict bound of Omegah2

    if cut4:
        cutDD = dataframe['PvalDD'] > 0.1

    if cut5:
        cutCMB = dataframe['CMB_ID'] < 1

    # Combine all cuts
    cut_tot = cutMD1_1 & cutMD1 & cutMD2 & cutDM2 & cutl345 & cutEWPT & cutOM & cutDD & cutCMB & cutMD1MD2 & cutMD2MDP & cutMDP

    # Apply the combined cuts
    dataframe_cut = dataframe[cut_tot]

    return dataframe_cut

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
        cbar.set_label('$\\Omega h_2$', fontsize=12)
    
    else:
        sc = plt.scatter(dataframe[df1], dataframe[df2], c = 'red', s = 1)
        plt.title(f'Plot of {label1} against {label2}')

    if omegah2bar & savefig:
        plt.savefig(f"plots/plot_{df1}_{df2}(no_grad){ext}.pdf", format='pdf') #pdf is full quality
    
    if omegah2bar == False & savefig == True:
        plt.savefig(f"plots/plot_{df1}_{df2}{ext}.pdf", format='pdf') #pdf is full quality

    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.xlabel(label1, fontsize=12)
    plt.ylabel(label2, fontsize=12)
    
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
    
'''
cut 1: constraints from model
cut 2: constraints from LEP 
cut 3: constraints from relic density
cut 4: DM DD constraints
cut 5: CMB constraints
'''

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
    """
    Plot data with different combinations of cuts dynamically.

    Parameters:
    - df: DataFrame
    - x_axis: column name for x-axis
    - y_axis: column name for y-axis
    - ylog_values: List of booleans to toggle y-log scale for each plot
    - omegah2bar: Whether to color by Omega h2
    - cut_flags: A list of dictionaries with cut booleans for each plot
    """
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