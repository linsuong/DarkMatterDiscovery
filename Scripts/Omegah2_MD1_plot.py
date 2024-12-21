
"""
    import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

# Specify the path to your file
files = [
         'scans/1-D_scans/Omegah2_MD1_l345_+1/scan.dat', 
         'scans/1-D_scans/Omegah2_MD1_l345_+1/scan_1.dat', 
         'scans/1-D_scans/Omegah2_MD1_l345_+1/scan_2.dat', 
         'scans/1-D_scans/Omegah2_MD1_l345_+1/scan_3.dat',
         'scans/1-D_scans/Omegah2_MD1_l345_+1/scan_4.dat',
         'scans/1-D_scans/Omegah2_MD1_l345_+1/scan_5.dat',
         'scans/1-D_scans/Omegah2_MD1_l345_+1/scan_6.dat',
         'scans/1-D_scans/Omegah2_MD1_l345_+1/scan_7.dat'
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
plt.savefig("scans/1-D_scans/Omegah2_MD1_l345_+1/plot_MD1_l345.pdf", format='pdf')
plt.show() 
    
    
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

files = [
    'scans/1-D_scans/Omegah2_MD1_l345_+1/scan.dat',
    'scans/1-D_scans/Omegah2_MD1_l345_+1/scan_1.dat',
    'scans/1-D_scans/Omegah2_MD1_l345_+1/scan_2.dat',
    'scans/1-D_scans/Omegah2_MD1_l345_+1/scan_3.dat',
    'scans/1-D_scans/Omegah2_MD1_l345_+1/scan_4.dat',
    'scans/1-D_scans/Omegah2_MD1_l345_+1/scan_5.dat',
    'scans/1-D_scans/Omegah2_MD1_l345_+1/scan_6.dat',
    'scans/1-D_scans/Omegah2_MD1_l345_+1/scan_7.dat'
]

desired_order = [1, 0.1, 0.01, 0.001, -1, -0.1, -0.01, -0.001]

lines = []
labels = []

# Read the file and plot data
for filepath in files:
    df = pd.read_csv(filepath, sep=r'\s+')
    l345_value = df['l345'][1]
    
    if l345_value > 0:   
        line, = plt.plot(df['MD1'], df['Omegah2'], '-', linewidth=1, rasterized=False,
                         label=f'$\lambda_{{345}} = {l345_value}$')
    else:
        line, = plt.plot(df['MD1'], df['Omegah2'], '--', linewidth=1, rasterized=False,
                         label=f'$\lambda_{{345}} = {l345_value}$')
    
    lines.append(line)
    labels.append(l345_value)

sorted_indices = [labels.index(val) for val in desired_order if val in labels]
sorted_lines = [lines[i] for i in sorted_indices]
sorted_labels = [f'$\lambda_{{345}} = {desired_order[i]}$' for i in range(len(sorted_indices))]

plt.text(x = 400, y = 0.17, s = '$\Omega h^2 \\approx 0.119$', color = 'red')
plt.axhline(y=0.11933, color='r', linestyle='dashdot', label='$\Omega h^2 = 0.12024$')

plt.xscale('log')
plt.yscale('log')
plt.xlim(10, 10e2)
plt.ylim(10e-7, 10e1)
plt.xlabel('$m_{h_1}$ (GeV)', fontsize=12)
plt.ylabel('Relic Density, $\\Omega h^2$', fontsize=12)
plt.title('Plot where $m_{h_2} = m_{h_\pm} = m_{h_1}$ + 1 GeV', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.gca().add_patch(Rectangle((0, 10e-7), 45, 10e1, color='red', alpha=0.2))
plt.tight_layout()

plt.legend(sorted_lines, sorted_labels, loc='upper right', prop={'size': 8}, ncol=2)

#plt.savefig("plots/plot_MD1_l345+1.pdf", format='pdf')
plt.show()
