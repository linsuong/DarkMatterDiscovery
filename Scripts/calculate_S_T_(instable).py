from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd

file_path = 'scans/5-D_scans/combined.dat'

data = pd.read_csv(file_path, sep=r'\s+', low_memory= False)

alpha = 1/137 #fine structure constant,  = e**2/(4 * pi * epsilon_0 * h_bar * c)
nu = 246 #VEV

def f_a(x):
    return -5 + 12 * np.log(x)

def f_b(x):
    return 3 - 4 * np.log(x)

def f_c(x, y):
    mask = np.isclose(x, y, rtol=1e-10)
    result = np.zeros_like(x)
    result[~mask] = ((x[~mask] + y[~mask]) / 2) - ((x[~mask] * y[~mask]) / (x[~mask] - y[~mask])) * np.log(x[~mask] / y[~mask])
    return result

data['MDP'] = data['DMP'] + data['MD1']
data['MD2'] = data['DM3'] + data['DMP'] + data['MD1']
data['DM2'] = data['DM3'] + data['DMP']
    
data['x1'] = data['MD1']/data['MDP']
data['x2'] = data['MD2']/data['MDP']
    
data['S'] = (1/(72 * np.pi * ((data['x2']**2 - data['x1']**2) ** 3))) * ((data['x2'] ** 6) * f_a(data['x2']) - 
                                                                            ((data['x1'] ** 6) * f_a(data['x1'])) + 
                                                                            (9 * (data['x2'] ** 2) * (data['x1']** 2) * 
                                                                            ((data['x2'] ** 2)) * f_b(data['x2']) - 
                                                                            (data['x1'] ** 2) * f_b(data['x1'])
                                                                            )
                                                                            )

data['T'] = (1/(32 * (np.pi ** 2) * alpha * (nu ** 2))) * (f_c(data['MDP'] ** 2, data['MD2'] ** 2) + 
                                                            f_c(data['MDP'] ** 2, data['MD1'] ** 2) -
                                                            f_c(data['MD2'] ** 2, data['MD1'] ** 2)
                                                            )
"""
version = 'all'


if version == 'less':
    cutT = (data['T'] < (T_central + T_error))
    cutS = (data['S'] < (S_central + S_error))

if version == 'more':
    cutT = (data['T'] > (T_central - T_error))
    cutS = (data['S'] > (S_central - S_error))

if version == 'all':
    cutT = (data['T'] > (T_central - T_error)) & (data['T'] < (T_central + T_error))
    cutS = (data['S'] > (S_central - S_error)) & (data['S'] < (S_central + S_error))

cuts = cutT & cutS
data = data[cuts]
    
cuts = cutT & cutS

data = data[cuts]

print(np.shape(data))
#plt.scatter(data['S'], data['T'], s = 0.1)

#plt.close('all')"""
S_central, S_error = 0.06, 0.09
T_central, T_error = 0.1, 0.07

plt.axvline(x = S_central - S_error, color = "black", lw = 0.1)
plt.axvline(x = S_central + S_error, color = "black", lw = 0.1)
plt.axhline(y = T_central + T_error, color = "black", lw = 0.1)
plt.axhline(y = T_central - T_error, color = "black", lw = 0.1)
 
plt.scatter(data['S'], data['T'], s=0.5)

plt.xlim(-0.2, 0.3)
plt.ylim(-0.1, 0.3)
plt.xlabel('S')
plt.ylabel('T')
plt.legend()
plt.show()
"""
plt.axvline(x = S_central - S_error, color = "black", lw = 0.1)
plt.axvline(x = S_central + S_error, color = "black", lw = 0.1)
plt.axhline(y = T_central + T_error, color = "black", lw = 0.1)
plt.axhline(y = T_central - T_error, color = "black", lw = 0.1)

sc = plt.scatter(data['S'], data['T'], s=0.5, c=data['DM4'], rasterized=True,
             cmap='plasma', norm=LogNorm(vmin=10e-2, vmax=data['DM4'].max()))

cbar = plt.colorbar(sc)
cbar.set_label('$m_{h_\pm} - m_{h_2}$', fontsize=12)

plt.scatter(data['S'], data['T'], s=0.5)
plt.xlim(-0.2, 0.3)
plt.ylim(-0.1, 0.3)
plt.title('S and T plot, unstable')
plt.xlabel('S')
plt.ylabel('T')
plt.legend()
plt.show()"""