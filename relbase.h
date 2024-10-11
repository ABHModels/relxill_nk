/*
   This file is part of the RELXILL model code.

   RELXILL is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   RELXILL is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.
   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.

    Copyright 2019 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef RELBASE_H_
#define RELBASE_H_

#define _GNU_SOURCE

#include "common.h"

#include "fitsio.h"

#include "relutility.h"
#include "reltable.h"
#include "rellp.h"
#include "xilltable.h"
#include "relcache.h"


/*********** DEFINE STATEMENTS *********/

#define PARAM_DEFAULT 0.0

#define version_major 1
#define version_minor 3
#define version_build 5
#define version_dev ""
#define version_major_nk 1
#define version_minor_nk 6
#define version_build_nk 3
#define version_dev_nk ""

/** path to all RELXILL tables */
#define RELXILL_TABLE_PATH "/media/linux/relxill_nk_v1.6.3"

/** dimensions of the RELLINE table */
#define RELTABLE_NA 25
#define RELTABLE_NR 100
#define RELTABLE_NG 40
#define RELTABLE_MAX_R 1000.0
#define RELTABLE_FILENAME "rel_table_v0.5a.fits"
#define RELTABLE_NMU0 30


/** dimensions of the LP table */
#define LPTABLE_NA 20
#define LPTABLE_NH 250
#define LPTABLE_FILENAME "rel_lp_table_v0.5b.fits"
#define LPTABLE_NR 100

#define RELRCTABLE_NA 19
#define RELRCTABLE_NH 34
// #define RELRCTABLE_NHR 15
// #define RELRCTABLE_NHR 244
#define RELRCTABLE_NHR 122
// #define RELRCTABLE_NH 100
#define RELRCTABLE_NR 100


/** parameters for interpolation an interagration **/
#define N_FRAD 1000      // values of radial bins (from rmin to rmax)
#define N_ZONES 10       // number of radial zones (as each zone is convolved with the input spectrum N_ZONES < N_FRAD)
#define N_ZONES_MAX 50  // maximal number of radial zones

/** parameters for the convolution **/
#define N_ENER_CONV  4096  // number of bins for the convolution, not that it needs to follow 2^N because of the FFT
#define EMIN_RELXILL 0.00035  // minimal energy of the convolution (in keV)
#define EMAX_RELXILL 2000.0 // minimal energy of the convolution (in keV)
#define EMIN_XILLVER 0.01
#define EMAX_XILLVER 1000.0

/** minimal and maximal energy for reflection strength calculation **/
#define RSTRENGTH_EMIN 20.0
#define RSTRENGTH_EMAX 40.0

/** file to (additionally) store the version number **/

// non-Kerr definations

// #define version_major_nk 1
// #define version_minor_nk 4
// #define version_build_nk 0
// #define version_dev_nk ""

// relbase.h
#define RELTABLE_NDEFPAR 30 							// nonkerr modifications

#define NK_FILENAME_1 "Trf_Johannsen_a13.fits"
#define NK_FILENAME_2 "Trf_Johannsen_a22.fits"
#define NK_FILENAME_3 "Trf_Johannsen_e3.fits"
#define NK_FILENAME_11 "Trf_KRZ_d1.fits"
#define NK_FILENAME_12 "Trf_KRZ_d2.fits"
#define NK_FILENAME_13 "Trf_KRZ_d3.fits"
#define NK_FILENAME_14 "Trf_KRZ_d4.fits"
#define NK_FILENAME_15 "Trf_KRZ_d5.fits"
#define NK_FILENAME_16 "Trf_KRZ_d6.fits"
// #define NK_FILENAME_4 "Trf_KRZ_d1_2.76e+00.fits"
// #define NK_FILENAME_5 "Trf_conformal_gravity.fits"
// #define NK_FILENAME_6 "Trf_KRZ_d2.fits"
#define NK_FILENAME_1_1 "Trf_a13_T_5.fits"
#define NK_FILENAME_1_2 "Trf_a13_T_10.fits"
#define NK_FILENAME_1_3 "Trf_a13_T_20.fits"
#define NK_FILENAME_1_4 "Trf_a13_T_30.fits"
#define LPNKTABLE_FILENAME_1 "rellp_nk_a13.fits"
#define LPNKTABLE_FILENAME_2 "rellp_nk_a22.fits"
// #define RELRCTABLE_FILENAME "Ring_corona_a13_new.fits"
// #define RELRCTABLE_FILENAME "Ring_corona_a13_2.fits"
#define RELRCTABLE_FILENAME "Ring_corona_a13new_del02_norm.fits"

#define RELNKTABLE_NA 30
#define RELNKTABLE_NMU0 22
#define RELNKTABLE_NR 100
#define RELNKTABLE_NG 40
#define RELNKTABLE_NG_OLD 20


/****** TYPE DEFINITIONS ******/

typedef struct{
	int save_g_ind;


	double cache_rad_relb_fun;
	double cache_val_relb_func[2];

	double cache_bin_ener;
	int cached_relbf;


	double re;
	double gmax;
	double gmin;
	double del_g;
	double emis;

	int limb_law;

	int ng;
	double** trff;
	double** cosne;
	double* gstar;
} str_relb_func;

/****** FUNCTION DEFINITIONS ******/

char* get_filename(int* status);  // non-Kerr mods
char* get_defparcolname(int* status);
char* get_lp_filename(int* status);
char* get_rc_filename(int* status);
/* the relbase function calculating the basic relativistic line shape for a given parameter setup*/
rel_spec* relbase(double* ener, const int n_ener,relParam* param, xillTable* xill_tab, int* status);

/** calculate the relline profile(s) for all given zones **/
void relline_profile(rel_spec* spec, relSysPar* sysPar, int* status);

void save_relline_profile(rel_spec* spec);

void free_cached_tables(void );

relSysPar* new_relSysPar(int nr, int ng, int* status);

void free_relSysPar(relSysPar* sysPar);

void relxill_kernel(double* ener_inp, double* spec_inp, int n_ener_inp, xillParam* xill_param, relParam* rel_param, int* status );

void relconv_kernel(double* ener_inp, double* spec_inp, int n_ener_inp, relParam* rel_param, int* status );

relSysPar* get_system_parameters(relParam* param, int* status);


/** function adding a primary component with the proper norm to the flux **/
void add_primary_component(double* ener, int n_ener, double* flu, relParam* rel_param, xillParam* xill_param, int* status);

void free_rel_spec(rel_spec* spec);
rel_spec* new_rel_spec(int nzones, const int n_ener, int*status);

/** caching routines **/
void init_specCache(specCache** spec, int* status);
void free_specCache(void);
void free_fft_cache(double*** sp,int n1, int n2);
void free_out_spec(out_spec* spec);
out_spec* init_out_spec(int n_ener, double* ener, int* status);

int redo_xillver_calc(relParam* rel_param, xillParam* xill_param, relParam* ca_rel, xillParam* ca_xill);
int redo_relbase_calc(relParam* rel_param, relParam* ca_rel_param);

void set_cached_rel_param(relParam* par, relParam** ca_rel_param, int* status);

int comp_xill_param(xillParam* cpar, xillParam* par);


#endif /* RELBASE_H_ */
