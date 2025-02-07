import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

file_path = 'scans/5-D_scans/combined.dat'

df = pd.read_csv(file_path, sep=r'\s+', low_memory= False)
df['brH_DMDM'] = pd.to_numeric(df['brH_DMDM'], errors='coerce')


params_dict = {
    'MD1': r'$m_{h1}$',
    'MD2': r'$m_{h2}$',
    'MDP': r'$m_{h_{\pm}}$',
    'Omegah2': r'$\Omega h_2$',
    'l345' : r'$\lambda_{345}$',
    'DMP' : r'$\Delta m^+$', #mass diff mhp - mh1
    'DM2' : r'$\Delta m^0$', #mass diff mh2 - mh1
    'DM3' : r'$\Delta m^1$' #mass diff mh2 - mh+
}

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
    dataframe['MD2'] = dataframe['DM3'] + dataframe['DMP'] + dataframe['MD1']
    dataframe['DM2'] = dataframe['DM3'] + dataframe['DMP']
    
    #columns for vacuum stability 
    dataframe['R'] = dataframe['l345'] /(2 * np.sqrt(l1))
    
    #columns for LEP constraints:
    dataframe['MD1+MD2'] = dataframe['MD1'] + dataframe['MD2']
    dataframe['MD1+MDP'] = dataframe['MD1'] + dataframe['MDP']
    dataframe['MD2+MDP'] = dataframe['MD2'] + dataframe['MDP']

    cutl345 = dataframe['l345'] > -np.inf
    cutOM = dataframe['Omegah2'] > -np.inf
    cutDD = dataframe['PvalDD'] > -np.inf
    cutCMB = dataframe['CMB_ID'] > -np.inf
    cutBr = dataframe['brH_DMDM'] > -np.inf
    
    cutMD1 = True 
    cutMDP = True
    cutl345 = True
    cutMass = True
    cutLEP = True
    cutLEP2 = True
    

    if cut1:
        #vs1: mh1^2 > 0 for |R| < 1
        vs1 = (dataframe['R'].abs() < 1) & (dataframe['MD1'] ** 2 > 0)
        
        #vs2: mh1^2 > threshold for R > 1
        threshold = (dataframe['R'] - 1) * np.sqrt(l1) * (v ** 2) 
        vs2 = (dataframe['MD1']**2 > threshold) & (dataframe['R'] > 1)
        
        cutMD1 = vs1 | vs2
        
        cutl345 = (
            (dataframe['l345'] < 2 * (((dataframe['MD1'] ** 2)/(v ** 2)) + np.sqrt(l1))) & 
            #(dataframe['l345'] < ((16/3) * np.pi) - l1) & 
            #(dataframe['l345'] < 4 * np.pi - ((3/2) * l1)) &
            (dataframe['l345'] > -1.47) &
            (dataframe['l345'] < 8 * np.pi)
        )
        
        cutMass = ((dataframe['MD1']<1000) & (dataframe['MD2']<1000) & (dataframe['MDP']<1000) &
                    (dataframe['MD1']>10) & (dataframe['MD2']>10) & (dataframe['MDP']>10))
        
    if cut2:
        #cutLEP2 = (dataframe['MD1'] > 80) & (dataframe['MD2'] > 100) & (dataframe['DM2'] < 8) 
        cutMDP = (dataframe['MDP'] > 70)
        cutLEP = ((dataframe['MD1+MD2'] > MZ) & (dataframe['MD1+MDP'] > MW) & 
                    (dataframe['MD2+MDP'] > MW) & (2 * dataframe['MDP'] > MZ))
        
        
        
    if cut3:
        
        cutOM = dataframe['Omegah2'] < 0.12024
        #cutOM = (df['Omegah2'] > 0.10) & (df['Omegah2'] < 0.12024) #strict bound of Omegah2
        
    if cut3_strict:
        cutOM = (df['Omegah2'] > 0.10737) & (df['Omegah2'] < 0.13123) #strict bound of Omegah2

    if cut4:
        cutDD = dataframe['PvalDD'] > 0.1

        cutCMB = dataframe['CMB_ID'] < 1
        
        cutBr = dataframe['brH_DMDM'] < 0.145

    # Combine all cuts
    cut_tot = (
                cutMD1 & cutl345 & cutMass & 
                cutLEP & cutLEP2 & cutMDP &
                cutOM & 
                cutDD & cutCMB & cutBr
                )

    # Apply the combined cuts
    dataframe_cut = dataframe[cut_tot]

    return dataframe_cut

def plotfig(dataframe, df1, df2, omegah2bar = False, colbar = True, xlog = True, ylog = True, savefig = False, label_dict = params_dict, info = None, label_y = True):

    dataframe
    label1 = label_dict.get(df1, df1)
    label2 = label_dict.get(df2, df2)

    if xlog:
        plt.xscale('log')
    if ylog:
        plt.yscale('log')

    if omegah2bar:
        sc = plt.scatter(dataframe[df1], dataframe[df2], c=dataframe['Omegah2'], rasterized=True, s=1,
                        cmap='plasma_r', norm=LogNorm(vmin=1e-6, vmax=dataframe['Omegah2'].max()))
        if colbar:
            cbar = plt.colorbar(sc)
            cbar.set_label('$\\Omega h_2$', fontsize=13)
    else:
        plt.scatter(dataframe[df1], dataframe[df2], c='red', s=1)

    plt.xlabel(label1, fontsize=13)
    
    if label_y == True:
        plt.ylabel(label2, fontsize=13)
    
    plt.tick_params(axis='both', labelsize=13)
    
    plt.grid(True, linestyle='--', alpha=0.7)
    
    plt.title(f'Plot of {label1} against {label2}', fontsize = 13)
    #plt.xlim()
    #plt.ylim()
    
    
    if savefig == True:   
        if omegah2bar == False:
            plt.savefig(f"plots/plot_{label1}_{label2}{info}(no_grad).pdf", format='pdf') #pdf is full quality
        
        if omegah2bar == True:
            plt.savefig(f"plots/plot_{label1}_{label2}{info}.pdf", format='pdf') #pdf is full quality 

def multiplotfig(dataframe, df1, df2, ax, omegah2bar = False, colbar = True, xlog = True, ylog = True, 
                 savefig = False, label_dict = params_dict, info = None, label_y = True):
    label1 = label_dict.get(df1, df1)
    label2 = label_dict.get(df2, df2)

    if xlog:
        ax.set_xscale('log')
    if ylog:
        ax.set_yscale('log')

    if omegah2bar:
        sc = ax.scatter(dataframe[df1], dataframe[df2], c=dataframe['Omegah2'], rasterized=True, s=1,
                        cmap='plasma_r', norm=LogNorm(vmin=1e-6, vmax=dataframe['Omegah2'].max()))
        if colbar:
            cbar = plt.colorbar(sc, ax=ax)
            cbar.ax.tick_params(labelsize=30)
            cbar.set_label('$\\Omega h_2$', fontsize=30)
    else:
        ax.scatter(dataframe[df1], dataframe[df2], c='red', s=1)

    ax.set_xlabel(label1, fontsize=30)
    
    if label_y == True:
        ax.set_ylabel(label2, fontsize=30)
    
    ax.tick_params(axis='both', labelsize=30)
    
    ax.grid(True, linestyle='--', alpha=0.7)
    
    
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
    #plt.show()


pairs = [
    ['MD1', 'MD2'], ['MD1', 'MDP'], ['MD1', 'DMP'], ['MD1', 'DM2'], ['MD1', 'DM3'],
    ['MD2', 'MD1'], ['MD2', 'MDP'], ['MD2', 'DMP'], ['MD2', 'DM2'], ['MD2', 'DM3'],
    ['MDP', 'MD1'], ['MDP', 'MD2'], ['MDP', 'DMP'], ['MDP', 'DM2'], ['MDP', 'DM3'],
]
pairs2 = [
    ['MD1', 'l345'], ['MD2', 'l345'], ['MDP', 'l345'],
    ['DMP', 'l345'], ['DM2', 'l345'], ['DM3', 'l345']
]

def oneplots():
    # Define your cuts and plotting loop
    for cut_number in range(1, 5):  # Iterating over cut levels (cut1, cut12, cut123, cut1234)
        cut_args = {f'cut{i+1}': True for i in range(cut_number)}  # Dynamically set cut arguments

        for x, y in pairs:
            df_f = cuts(df, **cut_args)  # Apply cuts dynamically
            info_label = ''.join([f'cut{i+1}' for i in range(cut_number)])  # Generate info string

            # Call plotfig with dynamic x, y, and info
            #print(info_label)
            #plotfig(df_f, x, y, omegah2bar=True, ylog= False, savefig=False, info=info_label)
            
            plotfig(df_f, x, y, omegah2bar=True, colbar = False, ylog= True, savefig=False, info=info_label)
                
            plt.savefig(f'plots/{x}_{y}_{info_label}.pdf')
            plt.close()

    for cut_number in range(1, 5):  # Iterating over cut levels (cut1, cut12, cut123, cut1234)
        cut_args = {f'cut{i+1}': True for i in range(cut_number)}  # Dynamically set cut arguments

        for x, y in pairs2:
            df_f = cuts(df, **cut_args)  # Apply cuts dynamically
            info_label = ''.join([f'cut{i+1}' for i in range(cut_number)])  # Generate info string

            # Call plotfig with dynamic x, y, and info
            #print(info_label)
            #plotfig(df_f, x, y, omegah2bar=True, ylog= False, savefig=False, info=info_label)
            
            plotfig(df_f, x, y, omegah2bar=True, colbar = True, ylog= False, savefig=False, info=info_label)
                
            plt.savefig(f'plots/{x}_{y}_{info_label}.pdf')
            plt.close()
            

def oneplots_nobar():
    # Define your cuts and plotting loop
    for cut_number in range(1, 5):  # Iterating over cut levels (cut1, cut12, cut123, cut1234)
        cut_args = {f'cut{i+1}': True for i in range(cut_number)}  # Dynamically set cut arguments

        for x, y in pairs:
            df_f = cuts(df, **cut_args)  # Apply cuts dynamically
            info_label = ''.join([f'cut{i+1}' for i in range(cut_number)])  # Generate info string

            # Call plotfig with dynamic x, y, and info
            #print(info_label)
            #plotfig(df_f, x, y, omegah2bar=True, ylog= False, savefig=False, info=info_label)
            
            plotfig(df_f, x, y, omegah2bar=True, colbar = False, ylog= True, savefig=False, info=info_label)
                
            plt.savefig(f'plots1/{x}_{y}_{info_label}.pdf')
            plt.close()

    for cut_number in range(1, 5):  # Iterating over cut levels (cut1, cut12, cut123, cut1234)
        cut_args = {f'cut{i+1}': True for i in range(cut_number)}  # Dynamically set cut arguments

        for x, y in pairs2:
            df_f = cuts(df, **cut_args)  # Apply cuts dynamically
            info_label = ''.join([f'cut{i+1}' for i in range(cut_number)])  # Generate info string

            # Call plotfig with dynamic x, y, and info
            #print(info_label)
            #plotfig(df_f, x, y, omegah2bar=True, ylog= False, savefig=False, info=info_label)
            
            plotfig(df_f, x, y, omegah2bar=True, colbar = False, ylog= False, savefig=False, info=info_label)
                
            plt.savefig(f'plots1/{x}_{y}_{info_label}.pdf')
            plt.close()

def four_plot():
    for x, y in pairs:
        fig, axes = plt.subplots(1, 4, figsize=(28, 7), constrained_layout=True)  # Use axes array
        ylog = True
        for cut_number in range(1, 5):
            cut_args = {f'cut{i+1}': True for i in range(cut_number)}
            df_f = cuts(df, **cut_args)
            print(x, y, cut_args)
            print(np.shape(df_f))
            
            if cut_number > 1 and cut_number < 4:
                axes[cut_number - 1].set_ylabel("")
                axes[cut_number - 1].tick_params(axis='y', left=False, labelleft=False)

                multiplotfig(df_f, x, y, axes[cut_number - 1], omegah2bar=True, ylog=ylog, colbar= False, label_y = False)

            elif cut_number == 4:
                axes[cut_number - 1].set_ylabel("")
                axes[cut_number - 1].tick_params(axis='y', left=False, labelleft=False)

                multiplotfig(df_f, x, y, axes[cut_number - 1], omegah2bar=True, ylog=ylog, colbar= True, label_y = False)
        
            elif cut_number == 1:
                multiplotfig(df_f, x, y, axes[cut_number - 1], omegah2bar=True, ylog=ylog, colbar= False)

        plt.suptitle(f'{params_dict.get(y, y)} against {params_dict.get(x, x)}', fontsize=35)
        plt.savefig(f'4plot/{y}_{x}.pdf')
        print(f'saved {y}_{x} plot')
        
    for x, y in pairs2:
        fig, axes = plt.subplots(1, 4, figsize=(28, 7), constrained_layout=True)  # Use axes array
        ylog = False
        for cut_number in range(1, 5):
            cut_args = {f'cut{i+1}': True for i in range(cut_number)}
            df_f = cuts(df, **cut_args)
            print(x, y, cut_args)
            print(np.shape(df_f))
            
            if cut_number > 1 and cut_number < 4:
                axes[cut_number - 1].set_ylabel("")
                axes[cut_number - 1].tick_params(axis='y', left=False, labelleft=False)

                multiplotfig(df_f, x, y, axes[cut_number - 1], omegah2bar=True, ylog=ylog, colbar= False, label_y = False)

            elif cut_number == 4:
                axes[cut_number - 1].set_ylabel("")
                axes[cut_number - 1].tick_params(axis='y', left=False, labelleft=False)

                multiplotfig(df_f, x, y, axes[cut_number - 1], omegah2bar=True, ylog=ylog, colbar= True, label_y = False)
        
            elif cut_number == 1:
                multiplotfig(df_f, x, y, axes[cut_number - 1], omegah2bar=True, ylog=ylog, colbar= False)

        plt.suptitle(f'{params_dict.get(y, y)} against {params_dict.get(x, x)}', fontsize=35)
        plt.savefig(f'4plot/{y}_{x}.pdf')
        print(f'saved {y}_{x} plot')
        
        #plt.show()

#oneplots()
#oneplots_nobar()
four_plot()

#df_cut = cuts(df, cut1= True, cut2= True, cut3= True, cut4= True)
#plotfig(df_cut, 'MD1', 'MD2', omegah2bar= True, colbar=True)
#plt.show()