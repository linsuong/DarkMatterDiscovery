from setuptools import setup, find_packages
import LAT_dwarfs

setup(
    name='LAT_dwarfs',
    version=LAT_dwarfs.__version__,
    include_package_data=True,
    description='Python code to derive data-driven upper limits on the thermally averaged, velocity-weighted pair-annihilation cross-section (velocity-independent) of a user-defined particle dark matter model using the expected differential gamma-ray spectrum of an annihilation event (provided by the user) as well as 10 years of Fermi-LAT data from observations of the Milky Way´s dwarf spheroidal galaxies.',

    # The project's main homepage.
    #url='https://',

    # Author details
    author='Francesca Calore, Pasquale Serpico, Bryan Zaldivar',
    author_email='calore@lapth.cnrs.fr, serpico@lapth.cnrs.fr, b.zaldivar.m@csic.es',

    # Choose your license
    license='Open Source. MIT license. See LICENSE file.',

    # See https://PyPI.python.org/PyPI?%3Aaction=list_classifiers
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],

    # What does your project relate to?
    keywords=['gamma-ray astronomy', 'particle physics', 'indirect dark matter searches'],

    packages=find_packages(),

    install_requires=[
        'numpy >= 1.19',
        'scipy >= 1.6.2',
        'iminuit <= 1.5.4',
        'astropy >= 4.1',
        'scikit-learn >= 0.24'
    ],


)
