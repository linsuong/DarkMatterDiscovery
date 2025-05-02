import numpy as np
import pandas as pd
import yaml

def cuts(dataframe, cut1=False, cut2=False, cut3=False, cut8_strict = False, cut4=False, cut5=False, cut6 = False, cut7 = False, cut8 = False):
    """
    Applies constraints to the dataframe and returns a filtered dataframe.
    """
    if cut3 == True and cut8_strict == True:
        raise Exception('please choose for a strict or relaxed bound on relic density')
    
    MW = 80.377
    MZ = 91.19
    l1 = 0.129
    v = 246
        
    yaml_file = "Data/SI_WS2022+WS2024.yaml"
    
    with open(yaml_file, 'r') as file:
        LZ = yaml.safe_load(file)
        
    conversion_factor = 1e36  # Convert from cm^-2 to pb
    y_data = {}
    # have MD1, DM3 and DMP, so convert to other values such as MDP, MD2 and DM2.
    dataframe['MDP'] = dataframe['DMP'] + dataframe['MD1']
    dataframe['MD2'] = dataframe['DM3'] + dataframe['DMP'] + dataframe['MD1']
    dataframe['DM2'] = dataframe['DM3'] + dataframe['DMP']
    
    #columns for vacuum stability 
    dataframe['R'] = dataframe['l345'] /(2 * np.sqrt(l1))
    
    #columns for LEP constraints:
    dataframe['MD1+MD2'] = dataframe['MD1'] + dataframe['MD2']
    dataframe['MD1+MDP'] = dataframe['MD1'] + dataframe['MDP']
    dataframe['MD2+MDP'] = dataframe['MD2'] + dataframe['MDP']

    #columns for EWPT constraints:
    dataframe['x1'] = dataframe['MD1']/dataframe['MDP']
    dataframe['x2'] = dataframe['MD2']/dataframe['MDP']
    
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
    dataframe['S'] = compute_S(dataframe['x1'], dataframe['x2'])
    dataframe['T'] = compute_T(dataframe['MDP'], dataframe['MD1'], dataframe['MD2'])

    dataframe['brH_DMDM'] = pd.to_numeric(dataframe['brH_DMDM'], errors='coerce')
    
    cutl345 = dataframe['l345'] > -np.inf
    cutOM = dataframe['Omegah2'] > -np.inf
    cutDD = dataframe['PvalDD'] > -np.inf
    cutCMB = dataframe['CMB_ID'] > -np.inf
    cutBr = dataframe['brH_DMDM'] > -np.inf
    
    cutMD1 = cutl345 = cutMass = cutLEP = cutMDP = cutOM = cutDD = cutCMB = cutBr = cutLZ = cutT = cutS = True
    
    if cut1:
        #vs1: mh1^2 > 0 for |R| < 1
        vs1 = (dataframe['R'].abs() < 1) & (dataframe['MD1'] ** 2 > 0)
        
        #vs2: mh1^2 > threshold for R > 1
        threshold = (dataframe['R'] - 1) * np.sqrt(l1) * (v ** 2) 
        vs2 = (dataframe['MD1']**2 > threshold) & (dataframe['R'] > 1)
        
        cutMD1 = vs1 | vs2
        
        cutl345 = (
            (dataframe['l345'] < 2 * (((dataframe['MD1'] ** 2)/(v ** 2)) + np.sqrt(l1))) & 
            (dataframe['l345'] < ((16/3) * np.pi) - l1) & 
            (dataframe['l345'] < 4 * np.pi - ((3/2) * l1)) &
            (dataframe['l345'] > -1.47) 
            #(dataframe['l345'] < 8 * np.pi)
        )
        
        cutMass = ((dataframe['MD1']<1000) & (dataframe['MD2']<1000) & (dataframe['MDP']<1000) &
                    (dataframe['MD1']>10) & (dataframe['MD2']>10) & (dataframe['MDP']>10))
        
    if cut2:
        # --- 1. Universal low-mass exclusion (Eq. 31) ---
        # Exclude ALL points with Mh1,Mh2 < 45 GeV OR Mh+ < 70 GeV
        cutLEP_universal = ~( 
            ((dataframe['MD1'] < 45) & (dataframe['MD2'] < 45)) | (dataframe['MDP'] < 70)
        )        
        # --- 2. Kinematic LEP cuts (Eq. 19) ---
        # Forbid Z/W decays to inert scalars
        cutLEP_kinematic = (
            (dataframe['MD1+MD2'] > MZ) & 
            (dataframe['MD1+MDP'] > MW) & 
            (dataframe['MD2+MDP'] > MW) & 
            (2 * dataframe['MDP'] > MZ)
        )
        
        # --- 3. Di-lepton LEP cut (Eq. 23) ---
        # Exclude Mh1 < 80 GeV + Mh2 < 100 GeV + ΔM > 8 GeV
        cutLEP_dilepton = ~(
            (dataframe['MD1'] < 80) & 
            (dataframe['MD2'] < 100) & 
            (dataframe['MD2'] - dataframe['MD1'] > 8)
        )
        
        # --- Combine all LEP cuts ---
        cutLEP = cutLEP_universal & cutLEP_kinematic & cutLEP_dilepton
        
        # --- Additional charged scalar mass cut ---
        cutMDP = (dataframe['MDP'] > 70)  # Redundant if already in cutLEP_universal
        
    if cut3:
        cutDD = dataframe['PvalDD'] > 0.1

    if cut4:
        cutCMB = dataframe['CMB_ID'] < 1
    
    if cut5:
        cutBr = dataframe['brH_DMDM'] < 0.145
        
    if cut6:
        if 'independent_variables' in LZ:
            for var in LZ['independent_variables']:
                if var['header']['name'] == 'mass':
                    x_values = [point['value'] for point in var['values']]
                    break
        
            if 'dependent_variables' in LZ:
                for var in LZ['dependent_variables']:
                    name = var['header']['name']
                    y_values = [point['value'] * conversion_factor for point in var['values']]
                    y_data[name] = y_values
                        
        cutLZ=(dataframe['protonSI'] < np.interp(dataframe['MD1'], x_values, y_data["limit"]))    
                     
    if cut7:
        #cutT = (dataframe['T'] > (0.1 - 0.07)) & (dataframe['T'] < (0.17))
        #cutS = (dataframe['S'] > (-0.03)) & (dataframe ['S'] < (0.06 + 0.09))
        
        cutT = (dataframe['T'] > (0.04 - 0.08)) & (dataframe['T'] < (0.04 + 0.08))
        cutS = (dataframe['S'] > (0.08 - 0.07)) & (dataframe ['S'] < (0.08 + 0.07))
        
    if cut8:
        #cutOM = dataframe['Omegah2'] < 0.12024
        #cutOM = (dataframe['Omegah2'] > 0.10) & (dataframe['Omegah2'] < 0.12024) #strict bound of Omegah2 (10%)
        cutOM = (dataframe['Omegah2'] > 0.084168) & (dataframe['Omegah2'] < 0.156312) #30% bound
        
    if cut8_strict:
        cutOM = (dataframe['Omegah2'] > 0.10737) & (dataframe['Omegah2'] < 0.13123) #strict bound of Omegah2

    # Combine all cuts
    cut_tot = cutMD1 & cutl345 & cutMass & cutLEP & cutMDP
    cut_tot &= cutT & cutS
    cut_tot &=  cutOM & cutDD & cutCMB & cutBr & cutLZ 
    #& cutT & cutS
    # Apply the combined cuts
    dataframe_cut = dataframe[cut_tot]

    return dataframe_cut
