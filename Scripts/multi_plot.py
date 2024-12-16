import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

# Specify the path to your file
file_path = 'run/5-D_scans/combined.dat'  # Replace with your actual file name

# Read the file into a pandas DataFrame; you can use .gz format to save space
df = pd.read_csv(file_path, sep=r'\s+')

# Display the DataFrame
#print("DataFrame:")
print(df['Omegah2'])
print(df.head())

# Example: Access a specific column as a vector
#vector = df['Omegah2'].values
#print("\nOmegah2 column as a vector:")
#print(vector)

# Filter the DataFrame to include rows where PvalDD > 0.1 

cutDD=(df['PvalDD'] > 0.1)
cutOM=(df['Omegah2'] < 0.12024)
#cutOM = (df['Omegah2'] > 0.11) & (df['Omegah2'] < 0.12024) #strict bound of Omegah2
cutCMB=(df['CMB_ID'] <1) 
df['MD2'] = df['DMP'] + df['DM3'] + df['MD1']
df['MDP'] = df['DMP'] + df['MD1']

cut_tot=(cutDD & cutOM)

df_f = df[cut_tot] #make dataframe smaller, quicker to deal with

params_dict = {
    'DM1': r'$m_{h1}$',
    'DM2': r'$m_{h2}$',
    'DMP': r'$m_{h_{\\pm}}$',
    'Omegah2': r'$\\Omega h_2$'
}

# Create a scatter plot of MD1 vs DMP, coloured by Omegah2
def plotfig(df1, df2, omegah2bar = False, xlog = True, ylog = True, savefig = False, label_dict = params_dict):
    label1 = label_dict.get(df1, df1)
    label2 = label_dict.get(df2, df2)
    
    plt.figure(figsize=(8, 6))

    if omegah2bar == True:
        sc = plt.scatter(df_f[df1], df_f[df2], c=df_f['Omegah2'], rasterized=True, s=1,
                        cmap='viridis', norm=LogNorm(vmin=df_f['Omegah2'].min(), vmax=df_f['Omegah2'].max()))

        plt.title(f'Plot of {label1} against {label2}, coloured by $\\Omega h_2$')
        
        cbar = plt.colorbar(sc)
        cbar.set_label('$\\Omega h_2$', fontsize=12)
    
    else:
        sc = plt.scatter(df_f[df1], df_f[df2], c = 'red', s = 1)
        plt.title(f'Plot of {label1} against {label2}')

        
    if xlog == True:
        plt.xscale('log')
    
    if ylog == True:
        plt.yscale('log')
        
    if savefig:
        plt.savefig(f"run/run_Dec2/plot_{df1}_{df2}(no_grad).pdf", format='pdf') #pdf is full quality

    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.xlabel(label1, fontsize=12)
    plt.ylabel(label2, fontsize=12)
    plt.show() 
    
def constraintplot(df1, df2, label_dict = params_dict):
    label1 = label_dict.get(df1, df1)
    label2 = label_dict.get(df2, df2)
    
    plt.figure(figsize=(8, 6))
    
    lz5tmedian_df = df_f[df_f['expName'] == 'LZ5Tmedian']

    notappl_df = df_f[df_f['expName'] == 'NotAppl']

    xenon1t_df = df_f[df_f['expName'] == 'XENON1T_2018']
    
    plt.scatter(notappl_df[df1], notappl_df[df2], c = 'blue', s = 0.8, label = 'NotAppl')
    plt.scatter(lz5tmedian_df[df1], lz5tmedian_df[df2], c = 'red', s = 0.8, label = 'LZ5Tmedian')
    plt.scatter(xenon1t_df[df1], xenon1t_df[df2], c = 'green', s = 0.8, label = 'XENON1T')
    #plt.scatter(notappl_df[df1], notappl_df[df2], c = 'blue', s = 0.8, label = 'NotAppl')

    plt.ylim(-2, 2)
    plt.xscale('log')
    plt.xlabel(label1)
    plt.ylabel(label2)
    plt.legend()
    #plt.savefig('run/run_Dec2/allowedParamsExps1.pdf')
    plt.show()
    

plotfig('MD2', 'MDP', omegah2bar= True, xlog = False, ylog = False)
