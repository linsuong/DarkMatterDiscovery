
"""
    import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

# Specify the path to your file
files = [
         'run/Omegah2_MD1_l345_+1/scan.dat', 
         'run/Omegah2_MD1_l345_+1/scan_1.dat', 
         'run/Omegah2_MD1_l345_+1/scan_2.dat', 
         'run/Omegah2_MD1_l345_+1/scan_3.dat',
         'run/Omegah2_MD1_l345_+1/scan_4.dat',
         'run/Omegah2_MD1_l345_+1/scan_5.dat',
         'run/Omegah2_MD1_l345_+1/scan_6.dat',
         'run/Omegah2_MD1_l345_+1/scan_7.dat'
         ]

#list of filepaths

# Read the file into a pandas DataFrame; you can use .gz format to save space
for i in range(len(files)):
    df = pd.read_csv(files[i], sep=r'\s+')
    if df['l345'][1] > 0:   
        plt.plot( df['MD1'], df['Omegah2'],'-',linewidth = 1, rasterized=True, label =f'$\lambda_{{345}} = ${df['l345'][1]}' )
        
    else:
        plt.plot( df['MD1'], df['Omegah2'],'--',linewidth = 1, rasterized=True, label =f'$\lambda_{{345}} = ${df['l345'][1]}' )
    
plt.axhline(y = 0.12024, color = 'r', linestyle = 'dashdot', label = '$\Omega h^2 = 0.12024$') 
#plt.gca().add_patch(Rectangle((0, 10e-7), 45, 10e1, color='red', alpha=0.2))
plt.xscale('log')
plt.yscale('log')
plt.xlim(10, 10e2)
plt.ylim(10e-7, 10e1)
plt.xlabel('$m_{h_1}$ (GeV)', fontsize=12)
plt.ylabel('Relic Density, $\\Omega h^2$', fontsize=12)
plt.title('Plot where $m_{h_2} = m_{h_\pm} = m_{h_1}$ + 1 GeV', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.legend(loc = 'upper right', prop={'size': 8}, ncol = 2)
plt.savefig("run/Omegah2_MD1_l345_+1/plot_MD1_l345.pdf", format='pdf')
plt.show() 
    
    
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# List of file paths
files = [
    'run/Omegah2_MD1_l345_+1/scan.dat',
    'run/Omegah2_MD1_l345_+1/scan_1.dat',
    'run/Omegah2_MD1_l345_+1/scan_2.dat',
    'run/Omegah2_MD1_l345_+1/scan_3.dat',
    'run/Omegah2_MD1_l345_+1/scan_4.dat',
    'run/Omegah2_MD1_l345_+1/scan_5.dat',
    'run/Omegah2_MD1_l345_+1/scan_6.dat',
    'run/Omegah2_MD1_l345_+1/scan_7.dat'
]

# Desired order for legend
desired_order = [1, 0.1, 0.01, 0.001, -1, -0.1, -0.01, -0.001]

# Prepare lists to store plot elements
lines = []
labels = []

# Read the file and plot data
for filepath in files:
    df = pd.read_csv(filepath, sep=r'\s+')
    l345_value = df['l345'][1]
    
    if l345_value > 0:   
        line, = plt.plot(df['MD1'], df['Omegah2'], '-', linewidth=1, rasterized=True,
                         label=f'$\lambda_{{345}} = {l345_value}$')
    else:
        line, = plt.plot(df['MD1'], df['Omegah2'], '--', linewidth=1, rasterized=True,
                         label=f'$\lambda_{{345}} = {l345_value}$')
    
    # Store the line and label
    lines.append(line)
    labels.append(l345_value)

# Sort lines and labels based on the desired order
sorted_indices = [labels.index(val) for val in desired_order if val in labels]
sorted_lines = [lines[i] for i in sorted_indices]
sorted_labels = [f'$\lambda_{{345}} = {desired_order[i]}$' for i in range(len(sorted_indices))]

# Reference line for Omega h^2
plt.axhline(y=0.12024, color='r', linestyle='dashdot', label='$\Omega h^2 = 0.12024$')

# Formatting the plot
plt.xscale('log')
plt.yscale('log')
plt.xlim(10, 10e2)
plt.ylim(10e-7, 10e1)
plt.xlabel('$m_{h_1}$ (GeV)', fontsize=12)
plt.ylabel('Relic Density, $\\Omega h^2$', fontsize=12)
plt.title('Plot where $m_{h_2} = m_{h_\pm} = m_{h_1}$ + 1 GeV', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()

# Add sorted legend
plt.legend(sorted_lines, sorted_labels, loc='upper right', prop={'size': 8}, ncol=2)

# Save and show the plot
plt.savefig("run/Omegah2_MD1_l345_+1/plot_MD1_l345.pdf", format='pdf')
plt.show()
