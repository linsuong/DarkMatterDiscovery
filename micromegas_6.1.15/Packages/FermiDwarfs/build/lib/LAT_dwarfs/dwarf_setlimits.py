from __future__ import division
from optparse import OptionParser
import math
import numpy as np
import os
import os.path
from scipy import interpolate
import matplotlib.pyplot as plt
import multiprocessing as mp
from multiprocessing import Pool
from functools import partial
from matplotlib import rcParams

def limits(case="J", ndwarf = ["Draco",], deltaEXCL = 3.84, filescan="misc/scan_mV_EV.dat", root = "profiling_case2"):
    
    # extract mV, EV values
    Nsamples = np.loadtxt(filescan, usecols=(0, ), unpack=True)
    Nsamples_tt = Nsamples.T if np.shape(Nsamples) != np.shape(1.0) else np.array([Nsamples])
    nscan = range(len(Nsamples_tt))
    
    if case == "J":
        if len(ndwarf) == 1:
            ndwarf_mod = ndwarf
            for i, nd in enumerate(ndwarf_mod): # run over dwarf
                
                os.system('rm excl_limits_case{}/sv_limits_case{}_'.format(case, case) + str(nd) + '.dat')
                fileout = open('excl_limits_case{}/sv_limits_case{}_'.format(case, case) + str(nd) + '.dat',"a+")

                for ns in nscan:
                    file = "sv_case{}_".format(case) + str(nd)+"_"+str(ns)+".dat"

                    arrays = np.loadtxt(root+file)
                    svarr = arrays[:, 0]
                    logL = arrays[:, -1]
                    
                    sv_new = np.logspace(math.log10(min(svarr)),math.log10(max(svarr)),10000)
                    svlogL_interp = interpolate.interp1d(svarr,logL,kind='cubic')
                    logL_new =svlogL_interp(sv_new)
                    
                    logLmin = min(logL_new)
                    kmin = np.argmin(logL_new)
                    
                    svEXCL =sv_new[kmin]

                    for k in range(kmin,len(sv_new)):
                        if(logL_new[k] - logLmin > deltaEXCL):
                                svEXCL = sv_new[k]
                                print(' sv excl = ',svEXCL)
                                break
                    mDM = Nsamples_tt[ns]
                    fileout.write(" ".join([str(mDM), str(svEXCL),"\n"]) )
        else:
            os.system('rm excl_limits_case{}/sv_limits_caseJ_'.format(case) + '+'.join(ndwarf) + '.dat')
            fileout = open('excl_limits_case{}/sv_limits_caseJ_'.format(case) + '+'.join(ndwarf) + '.dat', "a+")
        
            for ns in nscan:
                file = 'sv_caseJ_' + "+".join(ndwarf) + '_' + str(ns) +'.dat'
            
                arrays = np.loadtxt(root+file)
                svarr = arrays[:,0]
                logL = arrays[:,-1]

                sv_new = np.logspace(math.log10(min(svarr)),math.log10(max(svarr)),10000)
                svlogL_interp = interpolate.interp1d(svarr,logL,kind='cubic')
                logL_new =svlogL_interp(sv_new)
                    
                logLmin = min(logL_new)
                kmin = np.argmin(logL_new)

    #            plt.loglog(sv_new, logL_new)
    #            plt.loglog(svarr, logL, ".r")
    #            plt.axvline(x=sv_new[kmin])
                svEXCL =sv_new[len(sv_new)-1] # for some points in the scan we may not find the UL

                for k in range(kmin,len(sv_new)):
                    if(logL_new[k] - logLmin > deltaEXCL):
                        svEXCL = sv_new[k]
                        break

                # mV, EV corresponding values
                mDM = Nsamples_tt[ns]

                fileout.write(" ".join([str(mDM), str(svEXCL),"\n"]) )
            
    elif case == "JB":
        os.system('rm excl_limits_case{}/sv_limits_caseJB_'.format(case) + '+'.join(ndwarf) + '.dat')
        fileout = open('excl_limits_case{}/sv_limits_caseJB_'.format(case) + '+'.join(ndwarf) + '.dat', "+a")
    
        for ns in nscan:
            file = 'sv_caseJB_' + "+".join(ndwarf) + '_' + str(ns) + '.dat'
        
            arrays = np.loadtxt(root+file)
            svarr = arrays[:,0]
            logL = arrays[:,-1]

            sv_new = np.logspace(math.log10(min(svarr)),math.log10(max(svarr)),10000)
            svlogL_interp = interpolate.interp1d(svarr,logL,kind='cubic')
            logL_new =svlogL_interp(sv_new)
                
            logLmin = min(logL_new)
            kmin = np.argmin(logL_new)

#            plt.loglog(sv_new, logL_new)
#            plt.loglog(svarr, logL, ".r")
#            plt.axvline(x=sv_new[kmin])
            svEXCL =sv_new[len(sv_new)-1] # for some points in the scan we may not find the UL

            for k in range(kmin,len(sv_new)):
                if(logL_new[k] - logLmin > deltaEXCL):
                    svEXCL = sv_new[k]
                    break

            # mV, EV corresponding values
            mDM = Nsamples_tt[ns]

            fileout.write(" ".join([str(mDM), str(svEXCL),"\n"]) )

#            plt.axvline(x=svEXCL)
#            plt.title(str(ns))
#            plt.show()
    
    else:
        quit()


def interactive():
    parser = OptionParser()
    
    ############################################
    # case J: profiling over J (depending on listed dwarfs either single or stacked)
    # case JB: profiling over J & b (stacked)
    ############################################
    parser.add_option("-c", "--case", dest="case", help="Case: J - profiling only the J-factor; JB - combined profiling of dSph J-factor and background contribution", metavar="MODE", type='str', default = "J")
    
    # Sample dSph for either single case or combined limits
    parser.add_option("-d", "--dwarfs", dest="dwarf_list", help="dSph list", metavar="DWARF", type='str', default = "Draco")
    parser.add_option("--deltaEXCL", dest="deltaEXCL", help="delta in LogL for exclusion limit calculation", metavar="deltaEXCL", type='float', default =3.84)
    parser.add_option("--filescan", dest="filescan", help="Sample points list", metavar="SCAN", type='str', default = "misc/scan_DM_mass.dat")
    parser.add_option("--root", dest="root", help="Path for profiling files", metavar="ROOT", type='str', default = "profiling_case2/")
    
    (options, args) = parser.parse_args()
    
    if options.case is None:
        parser.print_help()
        quit()
    
    listd = options.dwarf_list.split('+')

    limits(case= options.case, ndwarf = listd, deltaEXCL = options.deltaEXCL, filescan = options.filescan, root = options.root)

if __name__ == '__main__':
    interactive()




