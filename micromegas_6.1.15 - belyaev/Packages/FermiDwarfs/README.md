Python code to derive data-driven upper limits on the thermally averaged, velocity-weighted pair-annihilation cross-section (velocity-independent; $s$-wave) of a user-defined particle dark matter model using the expected differential gamma-ray spectrum of pair-annihilation events (provided by the user) as well as 10 years of Fermi-LAT data from observations of the Milky Way´s dwarf spheroidal galaxies.

## Prerequisites

Python 3.6 or higher and the following packages:

    - numpy 
    - scipy
    - astropy
    - scikit-learn
    - iminuit (version < 2.0)

## Installation (general; not for micrOMEGA users)

A full version of this project can be installed as follows:

```sh
  $ git clone https://gitlab.in2p3.fr/christopher.eckner/mlfermilatdwarfs.git
  $ cd mlfermilatdwarfs
  $ pip install .
``` 

Note that the code is designed to be run via the command line interface as it requires parser arguments. However, each routine of the project maybe run on its own after the installation.

## Current test version (for micrOMEGA users)

The version that comes as a test case for interfacing it with micrOMEGA is a slimmer package which offers reduced verbosity and an easier handling adapted to the needs of a micrOMEGA user. The installation works as above but using the .zip package.

```sh
  $ cd ID_MLFermiDwarfs
  $ pip install .
``` 

The user needs to call the python script "run_ID_MLFermiDwarf.py" in the command line of the terminal and provide as an argument in the command line the file where the spectrum of the DM model is stored. For example:

```sh
  $ python3 run_ID_MLFermiDwarf.py idm.txt
``` 
The result will be a prompted sentence in the terminal stating that this point is either excluded or not. The script itself returns a Boolean value (0: not excluded, 1: excluded). The file "idm.txt" is an example of a micrOMEGA output file compatible with the code structure.
    
## License

This project is licensed under a MIT License - see the ``LICENSE`` file.

## Contact

Email to: calore [at] lapth.cnrs.fr / serpico [at] lapth.cnrs.fr / eckner [at] lapth.cnrs.fr
