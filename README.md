# README of the RELXILL_NK model

RELXILL_NK is a relativistic reflection model that is based on [RELXILL](http://www.sternwarte.uni-erlangen.de/~dauser/research/relxill/)
and extends it to non-Kerr spacetimes.
Please follow the instructions below carefully in order to ensure a properly working version of the model. This model is designed and tested to work within the X-ray spectral fitting software XSPEC. 1.6.3 is the current version of the model.

If you are using the RELXILL_NK model in your work please cite the following:

* Bambi et al., _Testing the Kerr black hole hypothesis using X-ray reflection spectroscopy_, [The Astrophysical Journal, 842, 76 (2017)](https://doi.org/10.3847/1538-4357/aa74c0)
* Abdikamalov et al., _Public Release of RELXILL_NK: A Relativistic Refletion Model for Testing Einstein's Gravity_, [The Astrophysical Journal, 878, 91 (2019)](https://doi.org/10.3847/1538-4357/ab1f89)

## Getting started

1. Installing the model:

  - Run following to execute the compile script:  
        chmod u+x compile_relxill.sh

  - Execute  
        ./compile_relxill.sh

2. Setting up the model environment:

  - Most importantly the model needs to know where the pre-calculated
    tables are located. By default it assumed that they are in the
    current working directory.

  - The recommended approach is to set the environment variable
    "RELXILL_TABLE_PATH" to the path where the tables are stored.  
    That means for a csh shell putting a line like the following in
    the .cshrc:  
        setenv RELXILL_TABLE_PATH /home/user/data/relline_tables/  
    or for a bash shell in .bashrc:  
        export RELXILL_TABLE_PATH='/home/user/data/relline_tables/'

3. Loading the model in XSPEC:

   From the directory where the model has been installed, the model can
   simply be loaded inside XSPEC:  
       lmod relxill_nk .

4. For questions and bug reports, please contact <relxill_nk@fudan.edu.cn>

## Usage instructions

* To work with the model, one would require certain FITS files.
    - Xillver table (which is the same as for relxilll and can be downloaded from the [webpage of relxill model](https://www.sternwarte.uni-erlangen.de/~dauser/research/relxill/index.html)).
    - FITS file of the transfer functions of the Johannsen metric are available here: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13906295.svg)](https://doi.org/10.5281/zenodo.13906295)
    - The models with lamppost/ring-like/disk-like coronae geometry require another FITS file (currently can be provided upon request, we will provide a link to download at a later stage).
 
 * The parameter def_par_type specifies the spacetime metric. One should freeze the value of def_par_type to the one associated to the non-Kerr metric of interest.  
     For example,  
     def_par_type = 1 specifies Johanssen metric with non-zero $\alpha_{13}$  
     def_par_type = 11 specifies KRZ metric with non-zero $\delta_{1}$  
 
 * The parameter mdot_type specifies the accretion rate in Eddington units and controls the thickness of the accretion disk according to the model proposed by [Taylor and Reynolds](https://iopscience.iop.org/article/10.3847/1538-4357/aaad63). NOTE: The FITS file for disk with finite thickness are only available for the Johanssen metric with non-zero $\alpha_{13}$.  
     mdot_type = 0 is infinitesimally thin disk  
     mdot_type = 1 is disk of finite thickness which has a mass accretion rate of 5% of the Eddington limit  
     mdot_type = 2 is disk of finite thickness which has a mass accretion rate of 10% of the Eddington limit  
     mdot_type = 3 is disk of finite thickness which has a mass accretion rate of 20% of the Eddington limit  
     mdot_type = 4 is disk of finite thickness which has a mass accretion rate of 30% of the Eddington limit  
 
 * The range of deformation parameters are not uniform and depend on the spin of the black hole, since XSPEC requires a uniform range of the parameter values we scale the deformation parameter values in the range [-1, 1]. We provide pyhton scripts (unscale.py and unscale_batch.py) to obtain the actual (unscaled) values of the deformation parameters after the fitting.
