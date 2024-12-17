import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

file_path = '/Users/linusong/Repositories/DarkMatterDiscovery/scans/5-D_scans/Dec_16_combined.dat'

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

def cuts(dataframe, cut1=False, cut2=False, cut3=False, cut4=False, cut5=False):
    """
    Applies constraints to the dataframe and returns a filtered dataframe.
    """

    dataframe['MD2'] = dataframe['DMP'] + dataframe['DM3'] + dataframe['MD1']
    dataframe['MDP'] = dataframe['DMP'] + dataframe['MD1']
    dataframe['DM2'] = dataframe['MD2'] - dataframe['MD1']

    cutl345 = dataframe['l345'] > -np.inf
    cutMD1 = dataframe['MD1'] > -np.inf
    cutMD2 = dataframe['MD2'] > -np.inf
    cutDM2 = dataframe['DM2'] > -np.inf
    cutOM = dataframe['Omegah2'] > -np.inf
    cutDD = dataframe['PvalDD'] > -np.inf
    cutCMB = dataframe['CMB_ID'] > -np.inf

    if cut1:
        cutl345 = dataframe['l345'] > -1.47

    if cut2:
        cutMD1 = dataframe['MD1'] < 80
        cutMD2 = dataframe['MD2'] < 100
        cutDM2 = dataframe['DM2'] < 8

    if cut3:
        cutOM = (dataframe['Omegah2'] < 0.12024)
        #cutOM = (df['Omegah2'] > 0.10) & (df['Omegah2'] < 0.12024) #strict bound of Omegah2


    if cut4:
        cutDD = dataframe['PvalDD'] > 0.1

    if cut5:
        cutCMB = dataframe['CMB_ID'] < 1

    # Combine all cuts
    cut_tot = cutMD1 & cutMD2 & cutDM2 & cutl345 & cutOM & cutDD & cutCMB

    # Apply the combined cuts
    dataframe_cut = dataframe[cut_tot]

    return dataframe_cut


'''
cutDD=(df['PvalDD'] > 0.1)
cutOM=(df['Omegah2'] < 0.12024)
#cutOM = (df['Omegah2'] > 0.10) & (df['Omegah2'] < 0.12024) #strict bound of Omegah2
#cutCMB=(df['CMB_ID'] <1) 


#cut_tot=(cutDD & cutOM)
cut_tot=(cutOM)
df_f = df[cut_tot]
'''

def plotfig(dataframe, df1, df2, omegah2bar = False, xlog = True, ylog = True, savefig = False, label_dict = params_dict):
    
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
                        cmap='viridis', norm=LogNorm(vmin=dataframe['Omegah2'].min(), vmax=dataframe['Omegah2'].max()))

        plt.title(f'Plot of {label1} against {label2}, coloured by $\\Omega h_2$')
        
        cbar = plt.colorbar(sc)
        cbar.set_label('$\\Omega h_2$', fontsize=12)
    
    else:
        sc = plt.scatter(dataframe[df1], dataframe[df2], c = 'red', s = 1)
        plt.title(f'Plot of {label1} against {label2}')

    if savefig:
        plt.savefig(f"plots/plot_{df1}_{df2}(no_grad).pdf", format='pdf') #pdf is full quality

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
    

df_f = cuts(df, cut2 = False, cut3 = True, cut4 = True)
plotfig(df_f, 'MD1', 'l345', omegah2bar= True, ylog = False, savefig = False)

df_f = cuts(df, cut3 = False, cut4 = True)
plotfig(df_f, 'MD1', 'l345', omegah2bar= True, ylog = False, savefig = False)

plt.show()