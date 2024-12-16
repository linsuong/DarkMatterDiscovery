import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

# Specify the path to your file
file_path = 'run/run_Dec11/scan.dat'  # Replace with your actual file name

# Read the file into a pandas DataFrame; you can use .gz format to save space
df = pd.read_csv(file_path, sep=r'\s+')

# Display the DataFrame
#print("DataFrame:")
#print(df['Omegah2'])
#print(df.head())

# Example: Access a specific column as a vector
#vector = df['Omegah2'].values
#print("\nOmegah2 column as a vector:")
#print(vector)

# Filter the DataFrame to include rows where PvalDD > 0.1 

#cutDM1=(df['MD1']) 
#cutDD=(df['PvalDD'] > 0.1)
#cutOM = (df['Omegah2'] > 0.11) & (df['Omegah2'] < 0.12024) #strict bound of Omegah2
#cutOM = df['Omegah2'] < 0.12024 #upper limit only
#cutCMB=(df['CMB_ID'] < 1) 

#cut_tot=(cutOM)

#df_f = df[cut_tot] #make dataframe smaller, quicker to deal with

#print(df_f)
#print(df_f['MD1'])
#print(df_f['Omegah2'])
# Create a scatter plot of MD1 vs DMP, coloured by Omegah2
#plt.figure(figsize=(8, 6))

#sc = plt.scatter(df_f['MD1'], df_f['DMP'], c=df_f['Omegah2'],
#                 cmap='viridis', norm=LogNorm(vmin=df_f['Omegah2'].min(), vmax=df_f['Omegah2'].max()))

sc = plt.scatter( df['MD1'], df['Omegah2'],rasterized=True, s = 0.5)
plt.axhline(y = 0.12024, color = 'r', linestyle = '-') 

#12 projections - this is one of them.
plt.xscale('log')
plt.yscale('log')

# Add colour bar
#cbar = plt.colorbar(sc)
#cbar.set_label('$\\Omega h_2$', fontsize=12)

# Add labels and title
plt.xlim(10)
plt.xlabel('$m_{h_1}$', fontsize=12)
plt.ylabel('$\\Omega h^2$', fontsize=12)
plt.title('Plot where $m_{h_2} = m_{h_\pm} = m_{h_1}$ + 1 GeV', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()

# Save  the plot to pdf format
#plt.savefig("run/run_Dec2/plot_MD1_l345.pdf", format='pdf') #pdf is full quality

# Show the plot
plt.show() 
