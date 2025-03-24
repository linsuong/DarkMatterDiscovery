import numpy as np
import pandas as pd
import yaml

def cuts(dataframe, cut1=False, cut2=False, cut3=False, cut3_strict = False, cut4=False, cut5=False, cut6 = False, cut7 = False, cut8 = False):
    """
    Applies constraints to the dataframe and returns a filtered dataframe.
    """
    if cut3 == True and cut3_strict == True:
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
    
    def f_a(x):
        return -5 + 12 * np.log(x)

    def f_b(x):
        return 3 - 4 * np.log(x)

    def f_c(x, y):
        mask = np.isclose(x, y, rtol=1e-10)
        result = np.zeros_like(x)
        result[~mask] = ((x[~mask] + y[~mask]) / 2) - ((x[~mask] * y[~mask]) / (x[~mask] - y[~mask])) * np.log(x[~mask] / y[~mask])
        return result
    
    alpha = 1/137
    nu = 246    
        
    dataframe['S'] = (1/(72 * np.pi * ((dataframe['x2']**2 - dataframe['x1']**2) ** 3))) * ((dataframe['x2'] ** 6) * f_a(dataframe['x2']) - 
                                                                                ((dataframe['x1'] ** 6) * f_a(dataframe['x1'])) + 
                                                                                (9 * ((dataframe['x2'] * dataframe['x1']) ** 2) * 
                                                                                ((dataframe['x2'] ** 2)) * f_b(dataframe['x2']) - 
                                                                                (dataframe['x1'] ** 2) * f_b(dataframe['x1']))
                                                                                )

    dataframe['T'] = (1/(32 * (np.pi ** 2) * alpha * (nu ** 2))) * (f_c(dataframe['MDP'] ** 2, dataframe['MD2'] ** 2) + 
                                                                f_c(dataframe['MDP'] ** 2, dataframe['MD1'] ** 2) -
                                                                f_c(dataframe['MD2'] ** 2, dataframe['MD1'] ** 2)
                                                                )
    
    
    
    cutl345 = dataframe['l345'] > -np.inf
    cutOM = dataframe['Omegah2'] > -np.inf
    cutDD = dataframe['PvalDD'] > -np.inf
    cutCMB = dataframe['CMB_ID'] > -np.inf
    cutBr = dataframe['brH_DMDM'] > -np.inf
    
    cutMD1 = cutMDP = cutl345 = cutMass = cutLEP = cutLEP2 = cutLZ = True
    
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
        cutLEP2 = (dataframe['MD1'] > 80) & (dataframe['MD2'] > 100) & (dataframe['DM2'] < 8) 
        cutMDP = (dataframe['MDP'] > 70)
        cutLEP = ((dataframe['MD1+MD2'] > MZ) & (dataframe['MD1+MDP'] > MW) & 
                    (dataframe['MD2+MDP'] > MW) & (2 * dataframe['MDP'] > MZ))
    
    if cut8:
        cutT = (dataframe['T'] > (0.1 - 0.07)) & (dataframe['T'] < (0.1 + 0.07))
        cutS = (dataframe['S'] > (0.06 - 0.09)) & (dataframe ['S'] < (0.06 + 0.09))
        
        cutT = (dataframe['T'] > (0.1 - 0.07)) & (dataframe['T'] < (0.1 + 0.07))
        cutS = (dataframe['S'] > (0.06 - 0.09)) & (dataframe ['S'] < (0.06 + 0.09))


    if cut3:
        cutOM = dataframe['Omegah2'] < 0.12024
        #cutOM = (df['Omegah2'] > 0.10) & (df['Omegah2'] < 0.12024) #strict bound of Omegah2
        
    if cut3_strict:
        cutOM = (dataframe['Omegah2'] > 0.10737) & (dataframe['Omegah2'] < 0.13123) #strict bound of Omegah2

    if cut4:
        cutDD = dataframe['PvalDD'] > 0.1

    if cut5:
        cutCMB = dataframe['CMB_ID'] < 1
    
    if cut6:
        cutBr = dataframe['brH_DMDM'] < 0.145
        
    if cut7:
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
                        
        cutLZ=(dataframe['protonSI'] > np.interp(dataframe['MD1'], x_values, y_data["limit"]))       
                     
    # Combine all cuts
    cut_tot = cutMD1 & cutl345 & cutMass & cutLEP & cutLEP2 & cutMDP
    cut_tot &=  cutOM & cutDD & cutCMB & cutBr & cutLZ 
    #& cutT & cutS
    # Apply the combined cuts
    dataframe_cut = dataframe[cut_tot]

    return dataframe_cut
