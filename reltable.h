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

    Copyright 2017 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef RELTABLE_H_
#define RELTABLE_H_

#include "relbase.h"
#include "relutility.h"



/** a single element in the RELLINE table array */
typedef struct{
	float* r;
	float* gmin;
	float* gmax;
	float** trff1;
	float** trff2;
	float** cosne1;
	float** cosne2;
}relDat;

/** the RELLINE table structure */
// typedef struct{
// 	float* a; // spin
// 	int n_a;

// 	float* mu0; // inclination
// 	int n_mu0;

// 	relDat*** arr; // relline data array

// 	// dimensions of relline array
// 	int n_r;
// 	int n_g;

// }relTable;

typedef struct{
	float* a; // spin
	int n_a;

	float** def_par; // deformation parameter  // need to be modified yet 
	int n_def_par; 

	float* mu0; // inclination
	int n_mu0;

	relDat**** arr; // relline data array

	// dimensions of relline array
	int n_r;
	int n_g;

}relnkTable; // non-Kerr modification



/** the LAMP POST single data structure */
typedef struct{
	float* h; // height

	float* rad; // radius

	float** intens;
	float** del;
	float** del_inc;

}lpDat;

/** the LAMP POST table structure */
// typedef struct{
// 	float* a; // spin
// 	int n_a;
// 	int n_h;
// 	int n_rad;
// 	lpDat** dat;
// }lpTable;

typedef struct{
	float* a; // spin
	int n_a;

	float** def_par; // deformation parameter  // need to be modified yet 
	int n_def_par; 

	int n_h;
	int n_rad;
	lpDat*** dat;

}lpnkTable;

typedef struct {
	float * rad;
	float * intens;
} rcDat;


typedef struct {
	float * a;
	int n_a;
	
	float ** def_par;
	int n_dp;
	
	float ** height;
	int n_h;
	
	float *r_ring;
	int n_r;

	int n_intens;
	float ****rad;
	float *****intens;
	
	// rcDat ***** dat;

} rcTable;


// rcTable* new_rcTable(int n_a, int n_dp, int n_h, int n_r, int *status);
rcTable* new_rcTable(int n_a, int n_dp, int n_h, int n_r, int n_intens, int* status);

void free_rcTable(rcTable * tab);

/* create a new LP table */
// lpTable* new_lpTable(int n_a, int n_h, int n_intens, int* status);
lpnkTable* new_lpnkTable(int n_a, int n_defpar, int n_h, int n_rad, int* status);

/* destroy the LP table structure */
// void free_lpTable(lpTable* tab);
void free_lpnkTable(lpnkTable* tab);

/* destroy the LP dat structure */
void free_lpDat(lpDat* dat, int nh);

/* create a new relline table structure */
// relTable* new_relTable(int n_a, int n_mu0, int n_r, int n_g, int* status);
relnkTable* new_relnkTable(int n_a, int n_mu0, int n_r, int n_g, int n_def_par, int* status);

/* destroy the relline table structure */
// void free_relTable(relTable* tab);
void free_relnkTable(relnkTable* tab);

/* routine to read the RELLINE table */
// void read_relline_table(char* filename, relTable** tab, int* status);
void read_rellinenk_table(char* filename, relnkTable** inp_tab, int* status);

/* routine to read the LP table */
// void read_lp_table(char* filename, lpTable** inp_tab, int* status);
void read_lp_nk_table(char* filename, lpnkTable** inp_tab, int* status);

void check_rcTable_cache(char* fname, rcTable* tab, int* ind, int* status);

fitsfile* open_rc_tab(char* filename, int* status);

void init_rc_table(char* filename, rcTable** inp_tab, int* status);

#endif /* RELTABLE_H_ */
