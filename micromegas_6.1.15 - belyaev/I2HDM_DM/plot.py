import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

# Specify the path to your file
file_path = '~/Documents/DMDiscovery/I2HDM_DM/scan.dat.gz'  # Replace with your actual file name

# Read the file into a pandas DataFrame; you can use .gz format to save space
df = pd.read_csv(file_path, sep=r'\s+', compression='gzip')

# Display the DataFrame
#print("DataFrame:")
print(df['Omegah2'])
print(df.head())

# Example: Access a specific column as a vector
#vector = df['Omegah2'].values
#print("\nOmegah2 column as a vector:")
#print(vector)

# Filter the DataFrame to include rows where PvalDD > 0.1 

cutDM1_test=(df['MD1'] <20)
cutDD=(df['PvalDD'] > 0.1)
cutOM=(df['Omegah2'] < 1.)
cutCMB=(df['CMB_ID'] <1)

cut_tot=(cutDD & cutOM & cutCMB)

df_f = df[cut_tot]

df_test=df[cut_tot & cutDM1_test]
print(df_test)
#df_f = df[ (df['PvalDD'] > 0.1)]



print(df_test[['MD1','DMP','DM3', 'l345', 'Omegah2', 'protonSI','PvalDD','CMB_ID','brH_DMDM']])
#print(df_f['PvalDD'] )

# Create a scatter plot of MD1 vs DMP, coloured by Omegah2
plt.figure(figsize=(8, 6))

#sc = plt.scatter(df_f['MD1'], df_f['DMP'], c=df_f['Omegah2'],
#                 cmap='viridis', norm=LogNorm(vmin=df_f['Omegah2'].min(), vmax=df_f['Omegah2'].max()))

sc = plt.scatter(df_f['MD1'], df_f['l345'], c=df_f['Omegah2'], rasterized=True,
                 cmap='viridis', norm=LogNorm(vmin=df_f['Omegah2'].min(), vmax=df_f['Omegah2'].max()))

plt.xscale('log')
#plt.yscale('log')

# Add colour bar
cbar = plt.colorbar(sc)
cbar.set_label('Omegah2', fontsize=12)

# Add labels and title
plt.xlabel('MD1', fontsize=12)
plt.ylabel('$\\lambda_{345}$', fontsize=12)
plt.title('Scatter Plot of MD1 vs DMP Coloured by $\\Omega h_2$', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)


plt.tight_layout()

# Save  the plot to pdf format
plt.savefig("plot.pdf", format='pdf')

# Show the plot
plt.show() 
