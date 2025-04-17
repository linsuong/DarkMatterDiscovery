from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
import matplotlib.colors as colors

file_path = 'scans/5-D_scans/combined.dat'

data = pd.read_csv(file_path, sep=r'\s+', low_memory= False)

data['MDP'] = data['DMP'] + data['MD1']
data['MD2'] = data['DM3'] + data['DMP'] + data['MD1']
data['DM2'] = data['DM3'] + data['DMP']
data['DM4'] = data['MDP'] - data['MD2']
    
data['x1'] = data['MD1']/data['MDP']
data['x2'] = data['MD2']/data['MDP']

"""
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
"""
data["S"] = ((1/(72*np.pi)) * (1/((((data["Mh2"]/data["Mh+"])**2) - ((data["MD1"]/data["Mh+"])**2))**3))
            * ((((data["Mh2"]/data["Mh+"])**6) * fa((data["Mh2"]/data["Mh+"]))) - (((data["MD1"]/data["Mh+"])**6)
            * fa((data["MD1"]/data["Mh+"]))) + (9 * ((data["Mh2"]/data["Mh+"])**2) * ((data["MD1"]/data["Mh+"])**2 )
            * ((((data["Mh2"]/data["Mh+"])**2) * fb((data["Mh2"]/data["Mh+"]))) - (((data["MD1"]/data["Mh+"])**2)
                                                                        * fb((data["MD1"]/data["Mh+"])))))))
data["T"] = ((1/(32*(np.pi**2)*alpha*(nu**2)))
            * (f_c((data["Mh+"]**2),(data["Mh2"]**2)) + f_c((data["Mh+"]**2),(data["MD1"]**2)) - f_c((data["Mh2"]**2),(data["MD1"]**2))))

"""

# Constants
alpha = 1 / 137  # Fine structure constant
nu = 246  # VEV in GeV

# Define helper functions
def f_a(x):
    return -5 + 12 * np.log(x)

def f_a_1(x):
    return 12/x

def f_a_2(x):
    return -12/(x ** 2)

def f_a_3(x):
    return 24/(x ** 3)

def f_b(x):
    return 3 - 4 * np.log(x)

def f_b_1(x):
    return -4/x

def f_b_2(x):
    return 4/(x ** 2)

def f_b_3(x):
    return -8/(x ** 3)

def f_c(x, y):
    mask = np.isclose(x, y, rtol=1e-10)
    result = np.zeros_like(x)
    result[~mask] = ((x[~mask] + y[~mask]) / 2) - ((x[~mask] * y[~mask]) / (x[~mask] - y[~mask])) * np.log(x[~mask] / y[~mask])
    
    return result

def compute_S(x1, x2):
    mask = np.isclose(x1, x2, rtol=1e-10)
    
    S = np.zeros_like(x1)
    
    # Case where x1 ≠ x2
    denominator = 72 * np.pi * ((x2**2 - x1**2) ** 3)
    numerator = (x2**6) * f_a(x2) - (x1**6) * f_a(x1) + (9 * (x2**2) * (x1**2)) * ((x2**2) * f_b(x2) - (x1**2) * f_b(x1))
    S[~mask] = numerator[~mask] / denominator[~mask]
    
    # Case where x1 ≈ x2 (use the limit)
    if np.any(mask):
        x = x1[mask]
        #S_limit = (1 / (24 * np.pi)) * (-5 + 12 * np.log(x) + 3 * x - 4 * x * np.log(x))
        
        S_limit = (1/(72 * 48 * np.pi * x)) * ((120 * (x ** 3) * f_a(x)) + (90 * (x ** 4) * f_a_1(x)) + (18 * (x ** 5) * f_a_2(x)) + ((x ** 6) * f_a_3(x)) + (216 * (x ** 3) * f_b(x)) + (324 * (x ** 4) * f_b_1(x)) + (108 * (x ** 5) * f_b_2(x)) + (9 * (x ** 4) * f_b_3(x)))
        
        S[mask] = S_limit
    
    return S


def compute_T(MDP, MD1, MD2):
    return (1 / (32 * (np.pi ** 2) * alpha * (nu ** 2))) * (
        f_c(MDP**2, MD2**2) 
        + f_c(MDP**2, MD1**2) 
        - f_c(MD2**2, MD1**2)
    )

# Apply to DataFrame
data['S'] = compute_S(data['x1'], data['x2'])
data['T'] = compute_T(data['MDP'], data['MD1'], data['MD2'])

version = 'all'

S_central, S_error = 0.06, 0.09
T_central, T_error = 0.1, 0.07

"""
if version == 'less':
    cutT = (data['T'] < (T_central + T_error))
    cutS = (data['S'] < (S_central + S_error))

if version == 'more':
    cutT = (data['T'] > (T_central - T_error))
    cutS = (data['S'] > (S_central - S_error))

if version == 'all':
    cutT = (data['T'] > (T_central - T_error)) & (data['T'] < (T_central + T_error))
    cutS = (data['S'] > (S_central - S_error)) & (data['S'] < (S_central + S_error))

if version == 'none':
    cutT = True
    cutS = True
    
cuts = cutT & cutS
data = data[cuts]
    
cuts = cutT & cutS

data = data[cuts]
"""
print(np.shape(data))
#plt.scatter(data['S'], data['T'], s = 0.1)

#plt.close('all')

plt.axvline(x = S_central - S_error, color = "black", lw = 0.1)
plt.axvline(x = S_central + S_error, color = "black", lw = 0.1)
plt.axhline(y = T_central + T_error, color = "black", lw = 0.1)
plt.axhline(y = T_central - T_error, color = "black", lw = 0.1)

norm = colors.SymLogNorm(linthresh=1e-3, vmin=data['DM4'].min(), vmax=data['DM4'].max())

sc = plt.scatter(data['S'], data['T'], s=0.5, c=data['DM4'], rasterized=True,
             cmap='plasma', norm=norm)

#vmin=data['DM4'].min(), vmax=data['DM4'].max()
cbar = plt.colorbar(sc)
cbar.set_label('$m_{h_\pm} - m_{h_2}$', fontsize=12)
        
plt.xlim(-0.2, 0.3)
plt.ylim(-0.1, 0.3)
plt.xlabel('S')
plt.ylabel('T')
plt.title('Plot of S against T', fontsize = 15)
plt.legend()
plt.show()
