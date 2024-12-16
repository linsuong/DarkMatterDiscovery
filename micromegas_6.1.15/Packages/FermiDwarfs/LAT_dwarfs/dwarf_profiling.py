from __future__ import division
from optparse import OptionParser
import math
import numpy as np
from astropy.io import fits
import os
import os.path
from scipy import stats
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.integrate import trapezoid as trapz
from scipy.signal import savgol_filter
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
import fileinput
import multiprocessing as mp
from multiprocessing import Pool
from functools import partial
from iminuit import Minuit
from matplotlib import rcParams
# this is an absurd command to make room for xlabel
rcParams.update({'figure.autolayout': True}) 


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++ some global declarations +++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# global parameters
rad = 0.5 # deg
Area = np.pi*rad**2
DOmega = Area*(np.pi/180.)**2 # Delta Omega of regions

# original energy binning of the default LAT data
NbinsMAX = 24
Elist0 = np.logspace(np.log10(0.5),np.log10(500.),NbinsMAX+1)

# Set path to dwarf working directory
file_path = "./data/misc/"


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++ declaration of functions        ++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def kde_skl(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scikit-learn"""
    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(x)
    # score_samples() returns the log-likelihood of the samples
    log_pdf = kde_skl.score_samples(x_grid[:, np.newaxis])
    return np.exp(log_pdf)

def optbandwidth(lista,bwmin,bwmax,Np):
    grid = GridSearchCV(KernelDensity(), {'bandwidth': np.linspace(bwmin, bwmax,Np)},cv=10)
    grid.fit(lista)
    string=str(grid.best_params_)
    trash,val = string.split()
    return val[:-1]

# reads integrated dNdE
# make it flexible by integration in the code to express the possibility of a user-defined ' of Ebins -> user provides input spectrum as dN/dE [impose units?]
# Need m_DM, <sv> for reference
# needed function to read the input table and perform the integration depending on the chosen Nebins (interpolation to make it easier)
def DMflux(sv, mDM, Ngammabin, log10_Jfac):
    #log10_Jfac=[18.2, 17.4, 17.6, 17.87, 19.0, 18.69, 18.03, 16.9, 18.2, 17.8, 17.47, 17.48, 16.3, 16.4, 17.6, 18.9, 18.40, 18.86, 18.09, 18.6, 17.9, 19.4, 18.41, 18.9, 17.9]
    #log10_Jfac=[0,0,0,17.87,0,18.69,18.03,0,0,0,17.47,17.48,0,0,0,0,18.40,18.86,18.09,0,0,0,18.41,0,0] # Justin's J
    flux = 1/(4.*np.pi)*sv/(mDM**2)*(10**log10_Jfac)*Ngammabin  #### Majorana/Dirac prefactor covered by micrOMEGA
    return flux

def f_dNdE_renormalized(E, E_samples, f_dNdE_initial):
    res = np.where(E_samples.max() < E, 1e-180, f_dNdE_initial(E) / 2.30258509299404568402 / E)
    return res ###E in GeV
    
def read_DM_spectrum_PPPC(M_DM, annihil_channel, is_wEW = True):
    assert annihil_channel in ["bb", "WW", "tautau"]
    
    if is_wEW:
        PPPC_column = {"bb": 13,
                       "WW": 17,
                        "tautau": 10,
        }
    else:
        PPPC_column = {"bb": 13,
                       "WW": 9,
                        "tautau": 10,
        }
    
    annihil_channel = PPPC_column[annihil_channel]

    #### Read PPPC table with electroweak corrections
    dNdE_wEW = np.loadtxt(file_path + "AtProduction_gammas.dat")
    masses_wEW = sorted(list(set(dNdE_wEW[:, 0])))
    ncol_E = int(len(dNdE_wEW) / len(masses_wEW))
    size_row = len(dNdE_wEW[0])
    dNdE_wEW = dNdE_wEW.reshape((len(masses_wEW), ncol_E, size_row))
    
    #### Read PPPC table without electroweak corrections
    dNdE_woEW = np.loadtxt(file_path + "AtProductionNoEW_gammas.dat")
    masses_woEW = sorted(list(set(dNdE_woEW[:, 0])))
    ncol_E = int(len(dNdE_woEW) / len(masses_woEW))
    size_row = len(dNdE_woEW[0])
    dNdE_woEW = dNdE_woEW.reshape((len(masses_woEW), ncol_E, size_row))
    
    if is_wEW:
        try:
            get_table_at_MDM = masses_wEW.index(M_DM)
            energies = np.power(10.,dNdE_wEW[get_table_at_MDM, :, 1]) * M_DM
            f_dNdE = interp1d(energies, dNdE_wEW[get_table_at_MDM, :, annihil_channel], fill_value = 'extrapolate')
            f_dNdE_final = lambda E: f_dNdE_renormalized(E, energies, f_dNdE)
            return f_dNdE_final
        except (RuntimeError, TypeError, NameError, ValueError):
            print("Requested dark matter mass not in the original PPPC table!")
    else:
        try:
            get_table_at_MDM = masses_woEW.index(M_DM)
            energies = np.power(10.,dNdE_woEW[get_table_at_MDM, :, 1]) * M_DM
            f_dNdE = interp1d(energies, dNdE_woEW[get_table_at_MDM, :, annihil_channel], fill_value = 'extrapolate')
            f_dNdE_final = lambda E: f_dNdE_renormalized(E, energies, f_dNdE)
            return f_dNdE_final
        except (RuntimeError, TypeError, NameError, ValueError):
            print("Requested dark matter mass not in the original PPPC table!")
    

def integrate_DM_dNdE(m_DM, spectrum_file = None, LAT_ebins = [0.5, 0.67, 0.89, 1.19, 1.58, 2.81, 500]):
    #### Expected file format:
    #### M_DM [GeV] ---- E [GeV] ---- dN/dE [GeV^-1]
    if spectrum_file:
        try:
            dNdE_data = np.loadtxt(spectrum_file, comments = '#')
        except ValueError:
            print("File not found or incorrect data table format. The user may add comments to the table. They have to be preceded by '#'.")
        
        selected_E = []
        selected_dNdE = []
        for M, E, dNdE in zip(dNdE_data[:, 0], dNdE_data[:, 1], dNdE_data[:, 2]):
            if M == m_DM:
                selected_E.append(E)
                selected_dNdE.append(dNdE)
        if len(selected_E) == 0:
            raise ValueError("Called dark matter mass is not part of the input differential spectrum table!")
        
        f_dNdE = interp1d(selected_E, selected_dNdE, kind = 'slinear', bounds_error = False, fill_value = 0.0)
    else:
        f_dNdE = read_DM_spectrum_PPPC(m_DM, "bb")
    
    DM_integrated = np.zeros(len(LAT_ebins) - 1)
    for k, (E1, E2) in enumerate(zip(LAT_ebins, LAT_ebins[1:])):
        sample_dNdE = np.logspace(np.log10(E1), np.log10(E2), 500)
        DM_integrated[k] += trapz(f_dNdE(sample_dNdE), sample_dNdE)
    return DM_integrated

def DMcount(sv, mDM, Ngammabin, Jmean_log10, exp,i,b):
    return float(DMflux(sv, mDM, Ngammabin,Jmean_log10[i]))*float(exp[i][b])

def getdSphs_data_from_table(dSphs_table, N_data_cols):
    dsph_data = np.loadtxt(dSphs_table, usecols = np.arange(1, N_data_cols))
    dshp_names = np.loadtxt(dSphs_table, dtype = "S", usecols = (0, ))
    dshp_names = list(map(lambda x: x.decode('UTF-8'), dshp_names))
    dshp_dict = {}
    for name, data_array in zip(dshp_names, dsph_data):
        dsph_properties = {
            "POS": (data_array[0], data_array[1]),    ###Galactic coordinates: LON/LAT [deg/deg]
            "JFACTOR": (data_array[2], data_array[3]),###J_factors as log_10(J) and scatter
            "EXP": data_array[4:]                     ###Fermi-LAT exposure at given dwarf position (uses the original data binning, rebinned during the profiling)
        }
        dshp_dict[name] = dsph_properties
    return dshp_dict

def binnedE(arr):
    b = [0,1,2,3,4,6,len(arr)]
    conc=[]
    for i in range(len(b)-1):
        conc.append(arr[b[i]])
    conc.append(arr[-1])
    return conc


def binned2(arr):
    b = [0,1,2,3,4,6,len(arr[0])]
    arr = np.array(arr)
    conc = arr[:,b[0]].reshape((-1,1))
    for i in range(1,len(b)-1):
        if(b[i]==b[i+1]-1):
            conc=np.concatenate((conc, arr[:,b[i]].reshape((-1,1))),axis=1 )
        else:
            bin=np.mean(arr[:,b[i]:b[i+1]],axis=1)
            conc=np.concatenate((conc,bin.reshape((-1,1))), axis=1)
    return conc
    

# this is the part of likelihood only with Poisson 
def logL_P(logJ,lnB,sv, mDM, Ngammabin,i,b,edata,lnN,log10Jmean):
    n = np.round(np.exp(lnN))
    n = n.astype('int')
    lamb = ((10**logJ)/(10**log10Jmean[i]))*DMcount(sv, mDM, Ngammabin,log10Jmean, edata,i,b) + np.exp(lnB)
    res = n[i]*math.log(lamb) - lamb - math.log(math.factorial(n[i]))
    return -2*res

# piece of likelihood on J
def logL_J(logJ,i,log10Jmean,dJ):
    denJ = np.log(10)*np.sqrt(2.*np.pi)*dJ[i]
    res = - (logJ-log10Jmean[i])**2/(2*dJ[i]**2) - math.log(denJ)
    return -2*res

# empirical distribution of J
def logL_J_emp(lista,x,bandwidth=0.2, **kwargs):
    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(lista)
    logpdf=kde_skl.score_samples(x.reshape(-1,1))        
    return -2.*logpdf[0]



def logL_B(lny,x,lnyi,xi,sigma,sigmaY):
    N = len(xi)
    X = np.array(x*len(xi)).reshape((-1,2))
    lnY = np.array([lny]*len(xi))
    AuxD = X-xi
    D =  AuxD[:,0]**2 + AuxD[:,1]**2 
    Dy = (lnY - lnyi)**2
    expX = np.longdouble(np.exp(-D/(2*sigma**2)))
    expY = np.longdouble( np.exp(-Dy/(2*sigmaY**2)))
    res = np.sum(expX*expY) /( (2*np.pi)**(3/2)*(sigma**2)*sigmaY*N  )
    return -2*np.log(res)


def Bratios(lnB):
    Nbins=len(lnB)
    lnB=np.array(lnB)
    ratios=[]
    for b in range(Nbins):
        ratios.append(lnB[b,:]/lnB[0,:])
    return ratios

def scanner(scan_params, dwarf_data, ylist, spectrum, log10J, dJ, lnC, expo, dSph_tag, is_empirical):
    n_scan, mDM = scan_params
    dNdE_integrated = integrate_DM_dNdE(mDM, spectrum_file = spectrum)
    Nbins=6
    j = list(dwarf_data.keys()).index(dSph_tag)
    Ngammabins = dNdE_integrated
    
    os.system('rm profiling_caseJ/sv_caseJ_'+ dSph_tag +'_'+str(n_scan)+'.dat')
    for sigv in ylist:
        fileout = open('profiling_caseJ/sv_caseJ_'+ dSph_tag +'_'+str(n_scan)+'.dat',"a+")
        #creating array of functions
        arrlogL = []
        for b in range(Nbins):
            arrlogL.append( lambda logJ, b=b: logL_P(logJ, dwarf_data[dSph_tag]["LOG_BKG"][1][b], sigv, mDM, Ngammabins[b],j,b,expo, lnC[b],log10J))
        
        P = lambda logJ: sum(arrlogL[i](logJ,i) for i in range(len(arrlogL)))
        if is_empirical:
            logL_totE = lambda logJ: P(logJ) + logL_J_emp(dwarf_data[dSph_tag]["J_EMP"],np.array(logJ).reshape(-1,1), bandwidth=2*dwarf_data[dSph_tag]["OPTBW"])
        else:
            logL_totE = lambda logJ: P(logJ) + logL_J(logJ, list(dwarf_data.keys()).index(dSph_tag), log10J, dJ) ###For Fermi J-factors
        m=Minuit(logL_totE, logJ =19., error_logJ=2., limit_logJ=(10.,22.),print_level=0, errordef =1)
        m.migrad()

        fileout.write(" ".join([str(sigv),str(m.values["logJ"]), str(m.errors["logJ"]),str(m.fval),"\n"]))
        fileout.close()


def scanner_stack(scan_params, ylist, dwarf_data, dlist, spectrum, expo, lnC, log10J, dJ, is_empirical):
    n_scan, mDM = scan_params
    dNdE_integrated = integrate_DM_dNdE(mDM, spectrum_file = spectrum)
    Nbins=6
    Ngammabins = dNdE_integrated
    dwarfs_emp = ["Draco", "Carina", "Fornax", "LeoI", "LeoII", "Sculptor", "Sextans", "UrsaMinor"]

    os.system('rm profiling_caseJ/sv_caseJ_' + "+".join(dlist) + '_' + str(n_scan) +'.dat')
    for sigv in ylist:
        fileout = open('profiling_caseJ/sv_caseJ_' + "+".join(dlist) + '_' + str(n_scan) +'.dat',"a+")
        #creating array of functions
        arrlogL2 = []
        for dSph_tag in dlist:
            arrlogL1 = []
            for b in range(Nbins):
                arrlogL1.append(lambda logJ_j, j = dSph_tag, b=b: logL_P(logJ_j,dwarf_data[j]["LOG_BKG"][1][b], sigv,mDM, Ngammabins[b], list(dwarf_data.keys()).index(j), b,expo, lnC[b],log10J))
            P_j = lambda logJ_j, j = dSph_tag: sum(arrlogL1[i](logJ_j, j ,i) for i in range(Nbins))
            if is_empirical and dSph_tag in dwarfs_emp:
                logtot_j_emp = lambda logJ_j, j = dSph_tag: P_j(logJ_j, j) + logL_J_emp(dwarf_data[j]["J_EMP"], np.array(logJ_j).reshape(-1,1),bandwidth=2 * dwarf_data[j]["OPTBW"])
                arrlogL2.append(lambda logJ_j, j = dSph_tag: logtot_j_emp(logJ_j,j))
            else:
                logtot_j = lambda logJ_j, j = dSph_tag: P_j(logJ_j, j) + logL_J(logJ_j, list(dwarf_data.keys()).index(j), log10J, dJ) ###For Fermi J-factors
                arrlogL2.append(lambda logJ_j, j = dSph_tag: logtot_j(logJ_j,j))
    
        def logL_tot(*args):
            args = args[0]
            s = sum([arrlogL2[k](var, dlist[k]) for k, var in enumerate(args)])
            return s
        
        kw = {}
        param_keys = []
        for dSph_tag in dlist:
            kw["J_"+dSph_tag] = 18.0
            kw["error_J_"+dSph_tag] = 2.0
            kw["limit_J_"+dSph_tag] = (10.,22.)
            param_keys.append("J_"+dSph_tag)
            
        m = Minuit(lambda x: logL_tot(x), use_array_call = True, forced_parameters = param_keys, print_level=0, pedantic = False, errordef=1, **kw)
        m.migrad()

        fit_parameter_ = []
        for key, val in m.values.items():
            fit_parameter_.append(str(val))
            fit_parameter_.append(str(m.errors[key]))
        fileout.write(" ".join([str(sigv)] + fit_parameter_ + [str(m.fval), "\n"]))

def scanner_stackB(scan_params, ylist, dwarf_data, dlist, spectrum, log10J, dJ, ratios, expo, lnC, lnCV0, xD, xV, is_empirical):
    n_scan, mDM = scan_params
    dNdE_integrated = integrate_DM_dNdE(mDM, spectrum_file = spectrum)
    Nbins=6
    Ngammabins = dNdE_integrated
    os.system('rm profiling_caseJB/sv_caseJB_' + "+".join(dlist) + '_' + str(n_scan) + '.dat')
    dwarfs_emp = ["Draco", "Carina", "Fornax", "LeoI", "LeoII", "Sculptor", "Sextans", "UrsaMinor"]

    # bkg PDF parameters
    for k in dwarf_data.keys():
        dwarf_data[k]["sigma"] = 1.4
        dwarf_data[k]["sigmaY"] = 0.15
    
    for sigv in ylist:
        fileout = open('profiling_caseJB/sv_caseJB_' + "+".join(dlist) + '_' + str(n_scan) + '.dat',"a+")
        #creating array of functions
        arrlogL2 = []
        for dSph_tag in dlist:
            arrlogL1 = []
            d_pos = lambda j: [xD[list(dwarf_data.keys()).index(j),0],xD[list(dwarf_data.keys()).index(j),1]]
            for b in range(Nbins):
                arrlogL1.append(lambda logJ, lnB, j = dSph_tag, b=b: logL_P(logJ, 
                                                                            lnB * ratios[b][list(dwarf_data.keys()).index(j)],
                                                                            sigv, 
                                                                            mDM, 
                                                                            Ngammabins[b], 
                                                                            list(dwarf_data.keys()).index(j),
                                                                            b,
                                                                            expo,
                                                                            lnC[b],
                                                                            log10J))
            P = lambda logJ, lnB, j = dSph_tag: sum(arrlogL1[i](logJ,lnB,j,i) for i in range(Nbins))
            if is_empirical and dSph_tag in dwarfs_emp:
                PJ_emp = lambda logJ, lnB, j = dSph_tag: P(logJ, lnB, j) + logL_J_emp(dwarf_data[j]["J_EMP"], np.array(logJ).reshape(-1,1), bandwidth=2 * dwarf_data[j]["OPTBW"])
                PJB_emp = lambda logJ, lnB, j = dSph_tag: PJ_emp(logJ, lnB, j) + logL_B(lnB, d_pos(j), lnCV0, xV,
                                                             dwarf_data[j]["sigma"], dwarf_data[j]["sigmaY"])
                arrlogL2.append(lambda logJ, lnB, j = dSph_tag: PJB_emp(logJ, lnB, j))
            else:
                PJ = lambda logJ, lnB, j = dSph_tag: P(logJ, lnB, j) + logL_J(logJ, list(dwarf_data.keys()).index(j), log10J, dJ) # for Fermi J-factors
                PJB = lambda logJ, lnB, j = dSph_tag: PJ(logJ, lnB, j) + logL_B(lnB, d_pos(j), lnCV0, xV,
                                                             dwarf_data[j]["sigma"], dwarf_data[j]["sigmaY"])
                arrlogL2.append(lambda logJ, lnB, j = dSph_tag: PJB(logJ, lnB, j))
            
            
        def logL_tot(*args):
            args = args[0]
            s = sum([arrlogL2[k](args[2*k], args[2*k+1], dSph) for k, dSph in enumerate(dlist)])
            return s
        
        kw = {}
        param_keys = []
        for dSph_tag in dlist:
            kw["J_"+dSph_tag] = 18.0
            kw["error_J_"+dSph_tag] = 2.0
            kw["limit_J_"+dSph_tag] = (10.,22.)
            param_keys.append("J_"+dSph_tag)
            kw["B_"+dSph_tag] = 4.0
            kw["error_B_"+dSph_tag] = 0.1
            kw["limit_B_"+dSph_tag] = (0.,7.)
            param_keys.append("B_"+dSph_tag)
        
        m = Minuit(lambda x: logL_tot(x), use_array_call = True, forced_parameters = param_keys, print_level=0, pedantic = False, errordef=1, **kw)
        m.migrad()

        fit_parameter_ = []
        for key in param_keys:
            fit_parameter_.append(str(m.values[key]))
            fit_parameter_.append(str(m.errors[key]))
        fileout.write(" ".join([str(sigv)] + fit_parameter_ + [str(m.fval), "\n"]))
                                        
def profiling(dwarf_data, case = "J", dwarf_list = ["Draco", ], emp_jlist = "./Jfactors", is_empirical = True, bkg_path = "./bkg_data", svmin=1e-27, svmax = 1e-22, ncpu=4, filescan='misc/scan_mV_EV.dat', spectrum = None):
    ###################################
    ### Loading parameters and files ##
    ###################################
    
    #   Jfactors mean, and error from Read 2019
    log10Jfacs, deltaJ = np.array([dSph["JFACTOR"][0] for k, dSph in dwarf_data.items()]), np.array([dSph["JFACTOR"][1] for k, dSph in dwarf_data.items()])
    
    #   optimum bandwidths obtained by companion code
    optbws = [0.011045,0.011045,0.011045,0.011045,0.011045,0.011045,0.021090,0.011045]
    
    #   Jfactors from Justin's files
    dwarfs_JJustin = ["Draco", "Carina", "Fornax", "LeoI", "LeoII", "Sculptor", "Sextans", "UMi"]

    for k, dSph_name in enumerate(dwarfs_JJustin):
        tmp_data = np.log10(np.loadtxt(emp_jlist + "/Jfac_" + dSph_name + ".txt"))
        if dSph_name == "UMi":
            dSph_name = "UrsaMinor"
        dwarf_data[dSph_name]["J_EMP"] = np.array(tmp_data).reshape(-1,1)
        dwarf_data[dSph_name]["OPTBW"] = optbws[k]
            
    #   bckg estimation @ dSphs
    os.system('rm fcbinds.dat')
    os.system('ls -v ' + bkg_path + '/lnbckg_dSphs_bin_*.dat > fcbinds.dat')
    fcbins = np.loadtxt('fcbinds.dat',dtype='str')
    lncounts_meas = []
    lnbckg_est = []
    for f in fcbins:
        vals = np.loadtxt(f)
        lncounts_meas.append(vals[:,0])
        lnbckg_est.append(vals[:,1])
        
    for k, dSph in enumerate(dwarf_data.keys()):
        dwarf_data[dSph]["LOG_BKG"] = (np.array(lncounts_meas)[:, k], np.array(lnbckg_est)[:, k])
        
    # binned dShps exposure
    exp_dSphs = np.array([dSph["EXP"] for key, dSph in dwarf_data.items()])
    binexp_dSphs=binned2(np.array(exp_dSphs))
    
    # dSphs positions
    lat_dSphs, lon_dSphs = np.array([dSph["POS"][1] for key, dSph in dwarf_data.items()]), np.array([dSph["POS"][0] for key, dSph in dwarf_data.items()])
    x_dSphs = np.concatenate((lon_dSphs.reshape((-1,1)),lat_dSphs.reshape((-1,1))),axis=1)
    
    # counts from voids, for PDF construction
    x_voids = np.loadtxt('./data/misc/voids_all_25dSphs.dat')
    count_voids = np.loadtxt('./data/misc/Newbins_counts_voids.dat')
    
    # remove all voids with zero counts in any of the macro energy bins
    reduced_voids = [(valid, void) for valid, void in zip(count_voids, x_voids) if valid.all() > 0]
    x_voids = np.array([x[1] for x in reduced_voids])
    lncountB0 = np.array([np.log(x[0][0]) for x in reduced_voids])
    
    # read the DM masses to sample from an external file provided by the filescan parser option
    Nsamples = np.loadtxt(filescan, usecols=(0, ), unpack=True)
    
    Nsamples_tt = Nsamples.T if np.shape(Nsamples) != np.shape(1.0) else np.array([[Nsamples]])
    svlist=np.logspace(np.log10(svmin),np.log10(svmax),10)

    # +++++++++ profiling according to the case ++++++++++++++
    if case == "J" and len(dwarf_list) == 1:
        scanner_1arg = partial(scanner,
                               dwarf_data = dwarf_data, 
                               ylist=svlist,
                               spectrum = spectrum,
                               log10J=log10Jfacs,
                               dJ = deltaJ,
                               lnC = lncounts_meas,
                               expo=binexp_dSphs,
                               dSph_tag = dwarf_list[0],
                               is_empirical = is_empirical,
        ) # works only for 1 dwarf in list
        
        p=Pool(ncpu)
        p.map(scanner_1arg, zip(range(len(Nsamples_tt)), Nsamples_tt))

    elif case == "J" and len(dwarf_list) > 1:

        scanner_1arg = partial(scanner_stack, 
                               dwarf_data = dwarf_data, 
                               ylist = svlist, 
                               dlist= dwarf_list, 
                               spectrum = spectrum,
                               expo = binexp_dSphs,
                               lnC = lncounts_meas,
                               log10J = log10Jfacs, 
                               dJ= deltaJ,
                               is_empirical = is_empirical,
        )

        p=Pool(ncpu)
        p.map(scanner_1arg, zip(range(len(Nsamples_tt)), Nsamples_tt))
        
    elif case == "JB":
      
        scanner_1arg = partial(scanner_stackB, 
                               ylist = svlist,
                               dwarf_data = dwarf_data, 
                               dlist = dwarf_list,
                               ratios=Bratios(lnbckg_est),
                               spectrum = spectrum,
                               log10J = log10Jfacs,
                               dJ= deltaJ,
                               expo=binexp_dSphs,
                               lnC=lncounts_meas,
                               lnCV0=lncountB0,
                               xD=x_dSphs,
                               xV=x_voids,
                               is_empirical = is_empirical,
        )

        p=Pool(ncpu)
        p.map(scanner_1arg, zip(range(len(Nsamples_tt)), Nsamples_tt))

def main():
    parser = OptionParser()
    
    ############################################
    # case J: profiling over J (depending on listed dwarfs either single or stacked)
    # case JB: profiling over J & b (stacked)
    ############################################
    parser.add_option("-c", "--case", dest="case", help="Case: J - profiling only the J-factor; JB - combined profiling of dSph J-factor and background contribution", metavar="MODE", type='str', default = "J")
    
    # Summary table of all dwarf spheroidal galaxies accessible with this code, the ordering of the dwarfs is essential for the script to function properly.
    # The table cotains, besides the names, the position of the dwarf in Galactic coordinates, its J-factor and uncertainty (if available) from Table 1 in arXiv:
    # 1611.03184 as well as the exposure of the Fermi LAT at the position of the dwarf with respect to the reference LAT data set.
    parser.add_option("-l", "--dictionary", dest="dwarf_dict", help="dSph dict", metavar="DICT", type='str', default = "default_dwarf_summary_table.dat")

    # Sample dSph for either single case or combined limits by giving the name of the respective dwarfs listed in the dwarf data table. Dwarfs needs have to be joined with a '+' without addtional spaces.
    parser.add_option("-d", "--dwarfs", dest="dwarf_list", help="dSph list", metavar="DWARF", type='str', default = "Draco")
    
    # The user may provide their own assessment of the J-factor and statistical scatter of a dwarf from a series of trials. This script comes with a number of dwarfs, which has been analysed in this respect. The profiling over J-factor uncertainties may include these values (1) or not (0). In the latter case, the values in the Fermi LAT publication arXiv:1611.03184 are used. Dwarfs for which these custom J-factors exist may be mixed with those for which analogous data are not available.
    parser.add_option("-j", "--jfactors", dest="j_type", help="Case: 1 -- use data-driven J-factors for those dwarfs for which they are available; 0 -- use the values of Table 1 in arXiv:1611.03184 of the Fermi-LAT collaboration", metavar="EMP", type=int, default = 1) 
    
    parser.add_option("--svmin", dest="svmin", help="svmin for scan", metavar="SVMIN", type='float', default =1e-27)
    parser.add_option("--svmax", dest="svmax", help="svmax for scan", metavar="SVMAX", type='float', default =1e-22)
    parser.add_option("--filescan", dest="filescan", help="Sample points file", metavar="SCAN", type='str', default = "misc/scan_DM_mass.dat")
    parser.add_option("--emp_J", dest="j_path", help="relative path to the folder with data-driven J-factor values per dwarf", metavar="JPATH", type='str', default = './Jfactors')
    parser.add_option("--lnB", dest="bkg_path", help="relative path to the folder with measured and predicted background events at a chosen dwarf position", metavar="BKGPATH", type='str', default = './bkg_data')
    parser.add_option("--spectrum", dest="spectrum", help="differential DM gamma-ray spectrum file\nenergies expected in GeV", metavar="DNDE", type='str', default=None)

    (options, args) = parser.parse_args()
    
    if options.case is None:
        parser.print_help()
        quit()

    num_cores = mp.cpu_count()
    num_cores -= 2 
    
    dwarf_table = options.dwarf_dict
    #### The second argument is the total number of columns in the dwarf data file. It needs to be updated by hand.
    dwarf_data = getdSphs_data_from_table(dwarf_table, 29)
    listd = options.dwarf_list.split('+')
        
    profiling(dwarf_data, case = options.case, dwarf_list = listd, emp_jlist = options.j_path, bkg_path = options.bkg_path, is_empirical = bool(options.j_type), svmin = options.svmin, svmax = options.svmax, ncpu=num_cores, filescan=options.filescan, spectrum = options.spectrum)

if __name__ == '__main__':
    main()


