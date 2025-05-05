import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data
df = pd.read_csv('cross_sections/filtered_data.dat', sep='\s+')

# Create figure
plt.figure(figsize=(15, 5))

# Plot for each DM3 value
for i, dm3 in enumerate([1, 10, 100]):
    plt.subplot(1, 3, i+1)
    subset = df[df['DM3'] == dm3].sort_values('DMP')
    
    # Plot limits as scatter with error bars
    plt.errorbar(subset['DMP'], subset['sigma_2l_A'],
                yerr=subset['err_2l_A'],
                fmt='r_', capsize=3, label='Limit')
    
    # Plot theory points - color by allowed/excluded
    for _, row in subset.iterrows():
        plt.scatter(row['DMP'], row['total_CS'], color='green', s=40)
    
    # Basic formatting
    plt.title(f'DM3 = {dm3} GeV')
    plt.xlabel('DMP (GeV)')
    plt.ylabel('Cross Section (fb)')
    plt.yscale('log')
    plt.grid(True, which='both', linestyle=':', alpha=0.3)
    if i == 0:
        plt.legend()

plt.tight_layout()
plt.savefig('scatter_limit_plot.png', dpi=120)
plt.show()