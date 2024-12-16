from LAT_dwarfs import dwarf_profiling, dwarf_setlimits
import sys, os
import numpy as np
import multiprocessing as mp

def main():
    spectrum_file_micrOMEGA = sys.argv[1]
    spectrum_data = np.loadtxt(spectrum_file_micrOMEGA)
    formatted_spectrum = np.stack((np.ones(len(spectrum_data[:, 0])), spectrum_data[:, 0], spectrum_data[:, 1]), axis =1)
    np.savetxt("./micrOMEGA_spectrum_formatted.dat", formatted_spectrum)
    
    case = "JB"
    
    dwarf_table = "./data/default_dwarf_summary_table.dat"
    dwarf_data = dwarf_profiling.getdSphs_data_from_table(dwarf_table, 29)
    dwarf_list = "Draco+Sculptor+LeoII+UrsaMinor".split('+')
    
    os.system("mkdir ./profiling_case{}".format(case))
    os.system("mkdir ./excl_limits_case{}".format(case))
    
    dwarf_profiling.profiling(dwarf_data, case = case, dwarf_list = dwarf_list, emp_jlist = "./data/Jfactors", bkg_path = "./data/bkg_data", is_empirical = True, svmin = 1e-2, svmax = 1e2, ncpu=mp.cpu_count()-2, filescan="./data/misc/scan_test_spec.txt", spectrum = "./micrOMEGA_spectrum_formatted.dat")
    dwarf_setlimits.limits(case= case, ndwarf = dwarf_list, deltaEXCL = 3.84, filescan = "./data/misc/scan_test_spec.txt", root = "./profiling_case{}/".format(case))
    res = np.loadtxt('./excl_limits_case{}/sv_limits_case{}_'.format(case, case) + '+'.join(dwarf_list) + '.dat')
    if res[1] > 1.0:
        print("Point in DM parameter space NOT excluded by Fermi-LAT gamma-ray data from dwarf spheroidal galaxies!")
    else:
        print("Point in DM parameter space excluded by Fermi-LAT gamma-ray data from dwarf spheroidal galaxies!")
    return int(res[1] <= 1.0)

if __name__ == "__main__":
    main()
