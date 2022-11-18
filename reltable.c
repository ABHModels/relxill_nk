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

#include "reltable.h"
#include "time.h"

extern int loaded_file_def_par;

static relDat* new_relDat(int nr, int ng, int* status){
	relDat* dat = (relDat*) malloc (sizeof(relDat));
	CHECK_MALLOC_RET_STATUS(dat,status,dat);

	dat->r = (float*) malloc( sizeof(float) * nr);
	CHECK_MALLOC_RET_STATUS(dat->r,status,dat);
	dat->gmin = (float*) malloc( sizeof(float) * nr);
	CHECK_MALLOC_RET_STATUS(dat->gmin,status,dat);
	dat->gmax = (float*) malloc( sizeof(float) * nr);
	CHECK_MALLOC_RET_STATUS(dat->gmax,status,dat);

	dat->cosne1 = (float**) malloc( sizeof(float*) * nr);
	CHECK_MALLOC_RET_STATUS(dat->cosne1,status,dat);
	dat->cosne2 = (float**) malloc( sizeof(float*) * nr);
	CHECK_MALLOC_RET_STATUS(dat->cosne2,status,dat);
	dat->trff1 = (float**) malloc( sizeof(float*) * nr);
	CHECK_MALLOC_RET_STATUS(dat->trff1,status,dat);
	dat->trff2 = (float**) malloc( sizeof(float*) * nr);
	CHECK_MALLOC_RET_STATUS(dat->trff2,status,dat);

	int ii;
	for (ii=0; ii<nr; ii++){
		dat->cosne1[ii] = (float*) malloc( sizeof(float) * ng);
		CHECK_MALLOC_RET_STATUS(dat->cosne1,status,dat);
		dat->cosne2[ii] = (float*) malloc( sizeof(float) * ng);
		CHECK_MALLOC_RET_STATUS(dat->cosne2,status,dat);
		dat->trff1[ii] = (float*) malloc( sizeof(float) * ng);
		CHECK_MALLOC_RET_STATUS(dat->trff1,status,dat);
		dat->trff2[ii] = (float*) malloc( sizeof(float) * ng);
		CHECK_MALLOC_RET_STATUS(dat->trff2,status,dat);
	}
	return dat;
}

/* create a new LP table structure*/
static lpDat* new_lpDat(int n_h, int n_rad, int* status){
	lpDat* dat = (lpDat*) malloc(sizeof(lpDat)*n_h);
	CHECK_MALLOC_RET_STATUS(dat,status,NULL);

	dat->h = (float*) malloc (sizeof(float)*n_h);
	CHECK_MALLOC_RET_STATUS(dat->h,status,dat);
	dat->rad = (float*) malloc (sizeof(float)*n_rad);
	CHECK_MALLOC_RET_STATUS(dat->rad,status,dat);

	dat->intens = (float**) malloc (sizeof(float*)*n_h);
	CHECK_MALLOC_RET_STATUS(dat->intens,status,dat);
	dat->del = (float**) malloc (sizeof(float*)*n_h);
	CHECK_MALLOC_RET_STATUS(dat->del,status,dat);
	dat->del_inc = (float**) malloc (sizeof(float*)*n_h);
	CHECK_MALLOC_RET_STATUS(dat->del_inc,status,dat);

	int ii;
	for (ii=0;ii<n_h;ii++){
		dat->intens[ii] = (float*) malloc (sizeof(float)*n_rad);
		CHECK_MALLOC_RET_STATUS(dat->intens[ii],status,dat);
		dat->del[ii] = (float*) malloc (sizeof(float)*n_rad);
		CHECK_MALLOC_RET_STATUS(dat->del[ii],status,dat);
		dat->del_inc[ii] = (float*) malloc (sizeof(float)*n_rad);
		CHECK_MALLOC_RET_STATUS(dat->del_inc[ii],status,dat);
	}
	return dat;
}

/** get a new and empty rel table (structure will be allocated)  */
// relTable* new_relTable(int n_a, int n_mu0, int n_r, int n_g, int* status){
// 	relTable* tab = (relTable*) malloc (sizeof(relTable));
// 	CHECK_MALLOC_RET_STATUS(tab,status,tab);

// 	// we know which dimensions the table should have
// 	tab->n_a   =  n_a;
// 	tab->n_mu0 =  n_mu0;
// 	tab->n_r   =  n_r;
// 	tab->n_g   =  n_g;

// 	tab->a = NULL;
// 	tab->mu0 = NULL;

// 	tab->arr=NULL;

// 	tab->arr = (relDat***) malloc (sizeof(relDat**)*tab->n_a);
// 	CHECK_MALLOC_RET_STATUS(tab->arr,status,tab);

// 	int ii; int jj;
// 	for (ii=0; ii<tab->n_a; ii++){

// 		tab->arr[ii] = (relDat**) malloc (sizeof(relDat*)*tab->n_mu0);
// 		CHECK_MALLOC_RET_STATUS(tab->arr[ii],status,tab);

// 		for (jj=0; jj<tab->n_mu0; jj++){
// 			tab->arr[ii][jj] = NULL;
// 			CHECK_STATUS_RET(*status,tab);
// 		}
// 	}

// 	return tab;
// }

// relTable* new_relTable(int n_a, int n_mu0, int n_r, int n_g, int* status)
relnkTable* new_relnkTable(int n_a, int n_mu0, int n_r, int n_g, int n_def_par, int* status){
	relnkTable* tab = (relnkTable*) malloc (sizeof(relnkTable));
	CHECK_MALLOC_RET_STATUS(tab,status,tab);

	// we know which dimensions the table should have
	tab->n_a   =  n_a;
	tab->n_mu0 =  n_mu0;
	tab->n_r   =  n_r;
	tab->n_g   =  n_g;
	tab->n_def_par = n_def_par;

	tab->a = NULL;
	tab->mu0 = NULL;
	tab->def_par = NULL;

	tab->arr=NULL;

	tab->arr = (relDat****) malloc (sizeof(relDat***)*tab->n_a);
	CHECK_MALLOC_RET_STATUS(tab->arr,status,tab);

	int ii; int jj; int kk;
	for (ii=0; ii<tab->n_a; ii++){

		tab->arr[ii] = (relDat***) malloc (sizeof(relDat**)*tab->n_def_par);
		CHECK_MALLOC_RET_STATUS(tab->arr[ii],status,tab);

		for (jj=0; jj<tab->n_def_par; jj++){
			tab->arr[ii][jj] = (relDat**) malloc (sizeof(relDat*)*tab->n_mu0);
			CHECK_MALLOC_RET_STATUS(tab->arr[ii][jj],status,tab);

			for(kk = 0; kk < tab->n_mu0; kk++){
				tab->arr[ii][jj][kk] = NULL;
				CHECK_STATUS_RET(*status,tab);
			}
		}

	}

	// printf(" new NK object is initialized\n");
	return tab;
}


/** read one axis of the rel table from the FITS file   */
static void get_reltable_axis(int nrows, float** val, char* extname, char* colname, fitsfile* fptr, int* status){

	int extver = 0;
	fits_movnam_hdu(fptr, BINARY_TBL, extname, extver ,status);
	if (*status!=EXIT_SUCCESS){
		printf(" *** error moving to extension %s\n",extname);
		return;
	}

	// get the column id-number
	int colnum;
	if(fits_get_colnum(fptr, CASEINSEN, colname, &colnum, status)) return;

	// get the number of rows
	long n;
	if (fits_get_num_rows(fptr, &n, status)) return;

	if (nrows != n){
		RELXILL_ERROR("wrong dimension of at least one axis given in the rel_table",status);
	}

	// allocate memory for the array
	*val=(float*)malloc(n*sizeof(float));
	CHECK_MALLOC_VOID_STATUS(*val,status);

    int anynul=0;
    double nullval=0.0;
    LONGLONG nelem = (LONGLONG) n;
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nelem ,&nullval,*val, &anynul, status);

	return;
}

/**     */
static void load_single_relDat_2dcol(fitsfile* fptr, float** val,int n1, int n2, int colnum, int* status){

    int anynul=0;
    double nullval=0.0;

    assert(val!=NULL);

    LONGLONG nelem = (LONGLONG) n2;

    int ii;
    for (ii=0; ii<n1;ii++){
        if(fits_read_col(fptr, TFLOAT, colnum, ii+1, 1, nelem ,&nullval,val[ii], &anynul, status)) return;

    }
}


/** load one single data extension from the relline table   */
static relDat* load_single_relDat(fitsfile* fptr, char* extname, int nhdu, int* status){

	// int extver = 0;
	// fits_movnam_hdu(fptr, BINARY_TBL, extname, extver ,status);
	int exttype;
	fits_movabs_hdu(fptr,nhdu,&exttype, status);
	if (*status!=EXIT_SUCCESS){
		printf(" *** error moving to extension %s\n",extname);
		return NULL;
	}

	// get the column id-number
	int colnum_r,colnum_gmin,colnum_gmax;
	int colnum_trff1, colnum_trff2;
	int colnum_cosne1, colnum_cosne2;
	if(fits_get_colnum(fptr, CASEINSEN, "r", &colnum_r, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "gmin", &colnum_gmin, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "gmax", &colnum_gmax, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "trff1", &colnum_trff1, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "trff2", &colnum_trff2, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "cosne1", &colnum_cosne1, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "cosne2", &colnum_cosne2, status)) return NULL;

	// check the number of rows (need to coincide with RELTABLE_NR
	long n;
	if (fits_get_num_rows(fptr, &n, status)) return NULL;
	if (n != RELNKTABLE_NR){
		RELXILL_ERROR("inconsistent number of rows in rel table",status);
		printf("    -> expecting %i, but found %ld in extensions %s",RELNKTABLE_NR,n,extname);
		return NULL;
	}

	// allocate the memory for the table
	relDat* dat = new_relDat(RELNKTABLE_NR,RELNKTABLE_NG,status);
	CHECK_STATUS_RET(*status,NULL);

	// now load the table column by column

    // (1) start with the 1D columns
    int anynul=0;
    double nullval=0.0;
    LONGLONG nelem = (LONGLONG) RELNKTABLE_NR;
    fits_read_col(fptr, TFLOAT, colnum_r, 1, 1, nelem ,&nullval,dat->r, &anynul, status);
    CHECK_STATUS_RET(*status,dat);
    fits_read_col(fptr, TFLOAT, colnum_gmin, 1, 1, nelem ,&nullval,dat->gmin, &anynul, status);
    CHECK_STATUS_RET(*status,dat);
    fits_read_col(fptr, TFLOAT, colnum_gmax, 1, 1, nelem ,&nullval,dat->gmax, &anynul, status);
    CHECK_STATUS_RET(*status,dat);

    // (2) and finally the 2D columns
    load_single_relDat_2dcol(fptr,dat->trff2,RELNKTABLE_NR,RELNKTABLE_NG,colnum_trff2,status);
    CHECK_STATUS_RET(*status,dat);
    load_single_relDat_2dcol(fptr,dat->trff1,RELNKTABLE_NR,RELNKTABLE_NG,colnum_trff1,status);
    CHECK_STATUS_RET(*status,dat);
    load_single_relDat_2dcol(fptr,dat->cosne1,RELNKTABLE_NR,RELNKTABLE_NG,colnum_cosne1,status);
    CHECK_STATUS_RET(*status,dat);
    load_single_relDat_2dcol(fptr,dat->cosne2,RELNKTABLE_NR,RELNKTABLE_NG,colnum_cosne2,status);
    CHECK_STATUS_RET(*status,dat);

    return dat;

}

/** load the complete relline table */
// void read_relline_table(char* filename, relTable** inp_tab, int* status){

// 	relTable* tab = (*inp_tab);
// 	fitsfile *fptr=NULL;

// 	char* fullfilename=NULL;
// 	char* extname=NULL;

// 	do{ // Errot handling loop
// 		if (tab != NULL){
// 			RELXILL_ERROR("relline table already loaded",status);
// 			break;
// 		}

// 		tab = new_relTable(RELTABLE_NA,RELTABLE_NMU0,RELTABLE_NR,RELTABLE_NG,status);
// 		CHECK_STATUS_BREAK(*status);

// 		// should be set by previous routine
// 		assert(tab!=NULL);
// 		assert(tab->arr!=NULL);

// 		// get the full filename
// 		if (asprintf(&fullfilename, "%s/%s", get_relxill_table_path() ,filename) == -1){
// 			RELXILL_ERROR("failed to construct full path the rel table",status);
// 			break;
// 		}

// 		// open the file
// 		if (fits_open_table(&fptr, fullfilename, READONLY, status)) {
// 			CHECK_RELXILL_ERROR("opening of the rel table failed",status);
// 			printf("    either the full path given (%s) is wrong \n",fullfilename);
// 			printf("    or you need to download the table ** %s **  from \n",filename);
// 			printf("    http://www.sternwarte.uni-erlangen.de/research/relxill/ \n");
// 			break;
// 		}

// 		// first read the axes of the table
// 		get_reltable_axis(tab->n_a, &(tab->a), "a", "a", fptr, status);
// 		CHECK_RELXILL_ERROR("reading of spin axis failed",status);

// 		get_reltable_axis(tab->n_mu0, &(tab->mu0), "mu0", "mu0", fptr, status);
// 		CHECK_RELXILL_ERROR("reading of mu0 axis failed",status);

// 		//now load the full table (need to go through all extensions)
// 		int ii; int jj;
// 		for (ii=0; ii<tab->n_a; ii++){
// 			for (jj=0; jj<tab->n_mu0; jj++){

// 				if (asprintf(&extname, "%i_%i", ii+1,jj+1) == -1){
// 					RELXILL_ERROR("failed to construct full path the rel table",status);
// 					break;
// 				}

// 				assert(tab->arr[ii][jj]==NULL);
// 				int nhdu = (ii)*tab->n_mu0+jj+4;
// 				tab->arr[ii][jj] = load_single_relDat(fptr, extname, nhdu, status);
// 				free(extname);
// 				if (*status!=EXIT_SUCCESS){
// 					RELXILL_ERROR("failed to load data from the rel table into memory",status);
// 					break;
// 				}
// 			}
// 		}

// 	} while(0);

// 	if (*status==EXIT_SUCCESS){
// 		// assigne the value
// 		(*inp_tab) = tab;
// 	} else {
// 		free_relTable(tab);
// 	}
// 	free(fullfilename);

// 	if (fptr!=NULL) {fits_close_file(fptr,status);}

// 	return;
// }


void read_rellinenk_table(char* filename, relnkTable** inp_tab, int* status){

	relnkTable* tab = (*inp_tab);
	fitsfile *fptr=NULL;

	char* fullfilename=NULL;
	char* extname=NULL;
	do{ // Errot handling loop
		// if (tab != NULL && !check_caching_defpartype()){
		if (tab != NULL){// && is_file_loaded()){
			RELXILL_ERROR("relline_nk table already loaded",status);
			break;
		}

		// if (tab != NULL){
		// 	free_relnkTable(tab);
		// 	tab = NULL;
		// }
		tab = new_relnkTable(RELNKTABLE_NA,RELNKTABLE_NMU0,RELNKTABLE_NR,RELNKTABLE_NG,RELTABLE_NDEFPAR,status);
		CHECK_STATUS_BREAK(*status);

		// should be set by previous routine
		assert(tab!=NULL);
		assert(tab->arr!=NULL);

		// get the full filename
		if (asprintf(&fullfilename, "%s/%s", get_relxill_table_path() , filename) == -1){
			RELXILL_ERROR("failed to construct full paththe rel table",status);
			break;
		}

		// open the file
		if (fits_open_table(&fptr, fullfilename, READONLY, status)) {
			CHECK_RELXILL_ERROR("opening of the rel table failed",status);
			printf("    either the full path given (%s) is wrong \n",fullfilename);
			printf("    or you need to download the table ** %s **  from \n",filename);
			printf("    http://www.sternwarte.uni-erlangen.de/research/relxill/ \n");
			break;
		}

		int _nhdu_ = 2;
		int exttype;
		fits_movabs_hdu(fptr, _nhdu_, &exttype, status);
		if (*status!=EXIT_SUCCESS){
			printf(" *** error moving to spin and defpar table %s\n",extname);
			return;
		}
		int colnum_a, colnum_defpar;
		long _n_;
		if(fits_get_colnum(fptr, CASEINSEN, "a", &colnum_a, status)) { printf(" error at fits get colnum a"); break; }
		// printf(" colnum_a -- %d\n", colnum_a);

		if(fits_get_colnum(fptr, CASEINSEN, get_defparcolname(status), &colnum_defpar, status)) { printf(" error at fits get colnum defpar"); break; }				// *** need to be checked more
		// printf(" colnum_defpar -- %d\n", colnum_defpar);

		if (fits_get_num_rows(fptr, &_n_, status)) return;

		if (_n_ != RELNKTABLE_NA){
			// printf(" %d %d\n", _n_, RELTABLE_NA);
			RELXILL_ERROR("inconsistent number of rows in spin table",status);
			return;
		}

		tab->a = (float*) malloc(tab->n_a * sizeof(float));

		int anynul=0;
		double nullval=0.0;
		LONGLONG nelem = (LONGLONG) _n_;
		fits_read_col(fptr, TFLOAT, colnum_a, 1, 1, nelem ,&nullval,tab->a, &anynul, status);

		if (*status!=EXIT_SUCCESS){
			printf(" *** error reading spin tables %s\n",extname);
			return;
		}

		tab->def_par = (float**) malloc( sizeof(float*) * tab->n_a);
		int ii;
		for (ii = 0; ii < tab->n_a; ii++){
			tab->def_par[ii] = (float*) malloc( sizeof(float) * tab->n_def_par);
		}

		anynul=0;
		nullval=0.0;
   		nelem = (LONGLONG) RELTABLE_NDEFPAR;
   		float vall[RELTABLE_NDEFPAR];

		for (ii=0; ii<RELNKTABLE_NA;ii++){
	        	if(fits_read_col(fptr, TFLOAT, colnum_defpar, ii+1, 1, nelem ,&nullval, tab->def_par[ii], &anynul, status)) return;
	        	// cpy(tab->defpar[ii], vall, RELTABLE_NDEFPAR);
		}

		if (*status!=EXIT_SUCCESS){
			printf(" *** error reading defpar table %s\n",extname);
			return;
		}

	    _nhdu_ = _nhdu_ + 1;
		fits_movabs_hdu(fptr, _nhdu_, &exttype, status);
		if (*status!=EXIT_SUCCESS){
			printf(" *** error moving to mu0 table %s\n",extname);
			return;
		}
		int colnum_mu0;
		if(fits_get_colnum(fptr, CASEINSEN, "mu0", &colnum_mu0, status)) { printf(" error at fits get colnum mu0"); break; }
		
		if (fits_get_num_rows(fptr, &_n_, status)) return;
		if (_n_ != RELNKTABLE_NMU0){
			RELXILL_ERROR("inconsistent number of rows in mu0 table",status);
			return;
		}

		tab->mu0 = (float*) malloc(tab->n_mu0 * sizeof(float));
	
		anynul=0;
		nullval=0.0;
		nelem = (LONGLONG) _n_;
		//fits_read_col(fptr, TFLOAT, colnum, 1, 1, nelem ,&nullval,*val, &anynul, status);
		fits_read_col(fptr, TFLOAT, colnum_mu0, 1, 1, nelem ,&nullval,tab->mu0, &anynul, status);

		if (*status!=EXIT_SUCCESS){
			printf(" *** error reading mu0 table %s\n",extname);
			return;
		}
		//now load the full table (need to go through all extensions)
		int jj, kk, nhdu = 4;
		for (ii=0; ii<tab->n_a; ii++){
			for (jj=0; jj < tab->n_def_par; jj++){    								// ????????
				 for(kk = 0; kk < tab->n_mu0; kk++) {

					if (asprintf(&extname, "%i_%i_%i", ii+1,jj+1,kk+1) == -1){ 									// 		?????? extname  - where is it used? *** found
						RELXILL_ERROR("failed to construct full path the rel nk table",status);
						break;
					}

					assert(tab->arr[ii][jj][kk]==NULL); 												// nonkerr modifications
					//int nhdu = (ii)*tab->n_defpar + (jj)*tab->n_mu0 + kk + 4; 							// nonkerr modifications
					tab->arr[ii][jj][kk] = load_single_relDat(fptr, extname, nhdu, status);
					nhdu = nhdu + 1;
					free(extname);
					if (*status!=EXIT_SUCCESS){
						RELXILL_ERROR("failed to load data from the rel table into memory",status);
						break;
					}	
				}
			}
		}

	} while(0);
	if (*status==EXIT_SUCCESS){
		// assigne the value
		(*inp_tab) = tab;
	} else {
		free_relnkTable(tab);
	}
	free(fullfilename);

	if (fptr!=NULL) {fits_close_file(fptr,status);}

	return;
}





// lpDat* load_single_lpDat(fitsfile* fptr, int n_h, int n_rad, int rownum, int* status){
// lpDat* load_single_lpDat(fitsfile* fptr, int n_h, int i, int j, int n_rad, int rownum, lpDat** lptab, double defpar[], int nhdu, int* status){
// 	lpDat* dat = new_lpDat(n_h,n_rad,status);
// 	CHECK_MALLOC_RET_STATUS(dat,status,NULL);

// 	int colnum_r;
// 	int colnum_hgrid;
// 	int colnum_h;
// 	int colnum_del;
// 	int colnum_del_inc;


// 	char* colname_h=NULL;
// 	char* colname_del=NULL;
// 	char* colname_del_inc=NULL;

// 	LONGLONG nelem;
//     int anynul=0;
//     double nullval=0.0;


// 	if(fits_get_colnum(fptr, CASEINSEN, "r", &colnum_r, status)) return NULL;
// 	if(fits_get_colnum(fptr, CASEINSEN, "hgrid", &colnum_hgrid, status)) return NULL;


// 	nelem = (LONGLONG) n_h;
//     fits_read_col(fptr, TFLOAT, colnum_hgrid, rownum, 1, nelem ,&nullval,dat->h, &anynul, status);

//     nelem = (LONGLONG) n_rad;
//     fits_read_col(fptr, TFLOAT, colnum_r, rownum, 1, nelem ,&nullval,dat->rad, &anynul, status);
//     CHECK_STATUS_RET(*status,NULL);


// 	int ii;
//     nelem = (LONGLONG) n_rad;
// 	for (ii=0; ii<n_h; ii++){

// 		if (asprintf(&colname_h, "h%i", ii+1) == -1){
// 			RELXILL_ERROR("failed to construct colname of the lp table",status);
// 			return NULL;
// 		}
// 		if(fits_get_colnum(fptr, CASEINSEN, colname_h, &colnum_h, status)) return NULL;
// 		free(colname_h);

// 		if (asprintf(&colname_del, "del%i", ii+1) == -1){
// 			RELXILL_ERROR("failed to construct colname of the lp table",status);
// 			return NULL;
// 		}
// 		if(fits_get_colnum(fptr, CASEINSEN, colname_del, &colnum_del, status)) return NULL;
// 		free(colname_del);

// 		if (asprintf(&colname_del_inc, "del_inc%i", ii+1) == -1){
// 			RELXILL_ERROR("failed to construct colname of the lp table",status);
// 			return NULL;
// 		}
// 		if(fits_get_colnum(fptr, CASEINSEN, colname_del_inc, &colnum_del_inc, status)) return NULL;
// 		free(colname_del_inc);

// 	    fits_read_col(fptr, TFLOAT, colnum_h, rownum, 1, nelem ,&nullval,dat->intens[ii], &anynul, status);
// 	    fits_read_col(fptr, TFLOAT, colnum_del, rownum, 1, nelem ,&nullval,dat->del[ii], &anynul, status);
// 	    fits_read_col(fptr, TFLOAT, colnum_del_inc, rownum, 1, nelem ,&nullval,dat->del_inc[ii], &anynul, status);

// 	    CHECK_STATUS_RET(*status,NULL);
// 	}


// 	return dat;
// }



double h_a13[20] = {2.0, 2.7, 3.0, 3.0, 3.01, 3.0, 2.99, 2.97, 2.96, 2.9, 2.74, 2.6, 2.51, 2.465, 2.37, 2.315, 2.256, 2.14, 2.045, 1.901};
double h_k[20] = {1.1589, 1.89, 2.13397, 2.24539, 2.27, 2.23798, 2.19106, 2.11397, 2.04, 1.95895, 1.86547, 1.75591, 1.69944, 1.6485, 1.5887, 1.53764, 1.47306, 1.4, 1.30573, 1.15537};
double hg[30];
int ci=-1;

lpDat* load_single_lp_nk_Dat(fitsfile* fptr, int n_h, int n_rad, int i, int j, int rownum, double defpar[], int nhdu, int* status){
	lpDat* dat = new_lpDat(n_h,n_rad,status);
	CHECK_MALLOC_RET_STATUS(dat,status,NULL);

	int exttype;
	fits_movabs_hdu(fptr,nhdu,&exttype, status);
	if (*status!=EXIT_SUCCESS){
		printf(" *** error moving to extension %d\n",nhdu);
		return NULL;
	}
	// int iii;

	int colnum_r;
	int colnum_h;
	int colnum_del;
	int colnum_del_inc;

	char* colname_h=NULL;
	char* colname_del=NULL;
	char* colname_del_inc=NULL;

	LONGLONG nelem;
    int anynul=0;
    double nullval=0.0;

	int k, zero;
	double dpdel, hdel, shift, shift2;
    	
	if(ci != i) {

		if (loaded_file_def_par == 1) {
			dpdel = (defpar[RELTABLE_NDEFPAR-1] - defpar[0]) / 28.0;
			shift = dpdel * (int)(defpar[0] / dpdel);
				hdel = (h_a13[i]-h_k[i]) / 29.0;
			shift2 = (defpar[0] - shift) * hdel / dpdel;
			zero = 0; 
			for(k=0; k<RELTABLE_NDEFPAR; k++) 
				if(defpar[k] == 0.0) {zero = k-1; k=RELTABLE_NDEFPAR;}
			hg[0] = h_a13[i];
			for(k=0; k<zero;k++) 
				hg[k+1] = h_a13[i]+(shift2 - 1.5*(int)k*hdel);
			hg[zero] = h_k[i];
			for(k=zero+1; k < RELTABLE_NDEFPAR; k++)
				// hg[k] = 1.1*h_k[i];
				hg[k] = 0.3 * (double)(k - zero) / (30.0 - zero) * h_k[i] + h_k[i];
		}
		if (loaded_file_def_par == 2) {
			for(k=0;k<RELTABLE_NDEFPAR;k++) {
				hg[k] = 1.1*h_k[i];
				if(defpar[k] == 0.0) { hg[k] = h_k[i]; }
			}
		}
		ci = i;
	}
	

	for (k = 0; k<250; k++) {
		dat->h[k] = ((float)(k) / 249.0) * ((float)(k) / 249.0) * (500.0 - (float)hg[j]) + (float)hg[j];
	}


	if(fits_get_colnum(fptr, CASEINSEN, "r", &colnum_r, status)) return NULL;

    nelem = (LONGLONG) n_rad;
    fits_read_col(fptr, TFLOAT, colnum_r, rownum, 1, nelem ,&nullval,dat->rad, &anynul, status);
    CHECK_STATUS_RET(*status,NULL);


	int ii;
    nelem = (LONGLONG) n_rad;
	for (ii=0; ii<n_h; ii++){

		if (asprintf(&colname_h, "h%i", ii+1) == -1){
			RELXILL_ERROR("failed to construct colname of the lp table",status);
			return NULL;
		}
		if(fits_get_colnum(fptr, CASEINSEN, colname_h, &colnum_h, status)) return NULL;
		free(colname_h);

		if (asprintf(&colname_del, "del%i", ii+1) == -1){
			RELXILL_ERROR("failed to construct colname of the lp table",status);
			return NULL;
		}
		if(fits_get_colnum(fptr, CASEINSEN, colname_del, &colnum_del, status)) return NULL;
		free(colname_del);

		if (asprintf(&colname_del_inc, "del_inc%i", ii+1) == -1){
			RELXILL_ERROR("failed to construct colname of the lp table",status);
			return NULL;
		}
		if(fits_get_colnum(fptr, CASEINSEN, colname_del_inc, &colnum_del_inc, status)) return NULL;
		free(colname_del_inc);

	    fits_read_col(fptr, TFLOAT, colnum_h, rownum, 1, nelem ,&nullval,dat->intens[ii], &anynul, status);
	    fits_read_col(fptr, TFLOAT, colnum_del, rownum, 1, nelem ,&nullval,dat->del[ii], &anynul, status);
	    fits_read_col(fptr, TFLOAT, colnum_del_inc, rownum, 1, nelem ,&nullval,dat->del_inc[ii], &anynul, status);

	    CHECK_STATUS_RET(*status,NULL);
	}


	return dat;
}





/** load the complete relline table */
// void read_lp_table(char* filename, lpTable** inp_tab, int* status){

// 	lpTable* tab = (*inp_tab);
// 	fitsfile *fptr=NULL;

// 	char* fullfilename=NULL;

// 	do{ // Errot handling loop
// 		if (tab != NULL){
// 			RELXILL_ERROR("relline LP table already loaded",status);
// 			break;
// 		}

// 		tab = new_lpTable(LPTABLE_NA,LPTABLE_NH,LPTABLE_NR,status);
// 		CHECK_STATUS_BREAK(*status);

// 		// should be set by previous routine
// 		assert(tab!=NULL);
// 		assert(tab->dat!=NULL);

// 		// get the full filename
// 		if (asprintf(&fullfilename, "%s/%s", get_relxill_table_path(),filename) == -1){
// 			RELXILL_ERROR("failed to construct full path the lp table",status);
// 			break;
// 		}

// 		// open the file
// 		if (fits_open_table(&fptr, fullfilename, READONLY, status)) {
// 			CHECK_RELXILL_ERROR("opening of the lp table failed",status);
// 			printf("    full path given: %s \n",fullfilename);
// 			break;
// 		}

// 		// first read the spin axes of the table AND move to the correct extension
// 		get_reltable_axis(tab->n_a, &(tab->a), "I_h", "a", fptr, status);
// 		CHECK_RELXILL_ERROR("reading of spin axis failed",status);


// 		//now load the full table (need to go through all extensions)
// 		int rownum=-1;
// 		int ii;
// 		for (ii=0; ii<tab->n_a; ii++){

// 			rownum = ii+1; // number of the row we read in the fits table

// 			assert(tab->dat[ii]==NULL);
// 			tab->dat[ii] = load_single_lpDat(fptr,tab->n_h, tab->n_rad, rownum, status);
// 			if (*status!=EXIT_SUCCESS){
// 				RELXILL_ERROR("failed to load data from the lp table into memory",status);
// 				break;
// 			}
// 		}

// 	} while(0);

// 	if (*status==EXIT_SUCCESS){
// 		// assigne the value
// 		(*inp_tab) = tab;
// 	} else {
// 		free_lpTable(tab);
// 	}
// 	free(fullfilename);

// 	if (fptr!=NULL) {fits_close_file(fptr,status);}

// 	return;
// }


void read_lp_nk_table(char* filename, lpnkTable** inp_tab, int* status){

	// printf(" lp 1.1\n");
	lpnkTable* tab = (*inp_tab);
	fitsfile *fptr=NULL;
	char* extname=NULL;

	char* fullfilename=NULL;

	do{ // Errot handling loop
		if (tab != NULL){   // needs modification, new caching proecedure for lamppost model
			RELXILL_ERROR("relline LP table already loaded",status);
			break;
		}

		// printf(" lp 1.2\n");

		tab = new_lpnkTable(LPTABLE_NA,RELTABLE_NDEFPAR,LPTABLE_NH,LPTABLE_NR,status);
		CHECK_STATUS_BREAK(*status);

		// should be set by previous routine
		assert(tab!=NULL);
		assert(tab->dat!=NULL);

		// get the full filename
		if (asprintf(&fullfilename, "%s/%s", get_relxill_table_path(), filename) == -1){
			RELXILL_ERROR("failed to construct full path the lp table",status);
			break;
		}

		// open the file
		if (fits_open_table(&fptr, fullfilename, READONLY, status)) {
			CHECK_RELXILL_ERROR("opening of the lp table failed",status);
			printf("    full path given: %s \n",fullfilename);
			break;
		}

		// printf(" lp 1.3\n");

		int _nhdu_ = 2;
		int exttype;
		fits_movabs_hdu(fptr, _nhdu_, &exttype, status);
		if (*status!=EXIT_SUCCESS){
			printf(" *** error moving to spin table %s\n",extname);
			return;
		}
		int colnum_a, colnum_defpar;
		long _n_;
		if(fits_get_colnum(fptr, CASEINSEN, "a", &colnum_a, status)) { printf(" Lamppost error at fits get colnum a"); break; }

		if(fits_get_colnum(fptr, CASEINSEN, get_defparcolname(status), &colnum_defpar, status)) { printf(" Lamppost error at fits get colnum defpar"); break; }				// *** need to be checked more
		// printf(" colnum_defpar -- %d\n", colnum_defpar);

		if (fits_get_num_rows(fptr, &_n_, status)) return;

		if (_n_ != LPTABLE_NA){
			RELXILL_ERROR("inconsistent number of rows in spin table",status);
			return;
		}

		// printf(" length of a, defpar tables -- %d\n", (int)_n_);

		tab->a = (float*) malloc(tab->n_a * sizeof(float));

		int anynul=0;
		double nullval=0.0;
		LONGLONG nelem = (LONGLONG) _n_;
		fits_read_col(fptr, TFLOAT, colnum_a, 1, 1, nelem ,&nullval,tab->a, &anynul, status);

		if (*status!=EXIT_SUCCESS){
			printf(" *** error reading spin tables %s\n",extname);
			return;
		}
		int ii;

    	tab->def_par = (float**) malloc( sizeof(float*) * tab->n_a);
		for (ii = 0; ii < tab->n_a; ii++){
			tab->def_par[ii] = (float*) malloc( sizeof(float) * tab->n_def_par);
		}

		anynul=0;
		nullval=0.0;
   		nelem = (LONGLONG) RELTABLE_NDEFPAR;
   		float vall[RELTABLE_NDEFPAR];

		for (ii=0; ii<LPTABLE_NA;ii++){
	        	if(fits_read_col(fptr, TFLOAT, colnum_defpar, ii+1, 1, nelem ,&nullval, tab->def_par[ii], &anynul, status)) return;
	        	// cpy(tab->defpar[ii], vall, RELTABLE_NDEFPAR);
		}

		if (*status!=EXIT_SUCCESS){
			printf(" *** error reading defpar table %s\n",extname);
			return;
		}

		//now load the full table (need to go through all extensions)
    	anynul=0;
    	nullval=0.0;
		int rownum=-1;
		int jj;

		rownum=-1;
		for (ii=0; ii<tab->n_a; ii++){
			for(jj=0; jj<tab->n_def_par; jj++) {

				rownum = ii+1; // number of the row we read in the fits table
				_nhdu_++;
			
				// tab->dat[ii][jj] = load_single_lpDat(fptr,tab->n_h, tab->n_rad, rownum, status);
				// tab->dat[ii][jj] = new_lpDat(LPTABLE_NH,LPTABLE_NR,status);
				tab->dat[ii][jj] = load_single_lp_nk_Dat(fptr, tab->n_h, tab->n_rad, ii, jj, 1, tab->def_par[ii], _nhdu_, status);

				if (*status!=EXIT_SUCCESS){
					RELXILL_ERROR("failed to load data from the lp table into memory",status);
					break;
				}
			}
		}

	} while(0);

	// printf(" lp 1.6\n");

	if (*status==EXIT_SUCCESS){
		// assigne the value
		(*inp_tab) = tab;
	} else {
		free_lpnkTable(tab);
	}
	free(fullfilename);

	if (fptr!=NULL) {fits_close_file(fptr,status);}

	return;
}




static void free_relDat(relDat* dat, int nr){
	if (dat!=NULL){
		int ii;
		for (ii=0; ii<nr; ii++){
			if (dat->cosne1 !=NULL) free(dat->cosne1[ii]);
			if (dat->cosne2 !=NULL) free(dat->cosne2[ii]);
			if (dat->trff1 !=NULL) free(dat->trff1[ii]);
			if (dat->trff2 !=NULL) free(dat->trff2[ii]);
		}
		free(dat->cosne1);
		free(dat->cosne2);
		free(dat->trff1);
		free(dat->trff2);

		free(dat->r);
		free(dat->gmin);
		free(dat->gmax);
	}
}

// void free_relTable(relTable* tab){
// 	if(tab!=NULL){
// 		if (tab->arr!=NULL){
// 			int ii;
// 			for (ii=0; ii<tab->n_a; ii++){
// 				if (tab->arr[ii] !=NULL){
// 					int jj;
// 					for (jj=0; jj<tab->n_mu0; jj++){
// 						free_relDat(tab->arr[ii][jj],tab->n_r);
// 						free(tab->arr[ii][jj]);
// 					}
// 					free(tab->arr[ii]);
// 				}
// 			}
// 			free(tab->arr);
// 		}
// 		free(tab->a);
// 		free(tab->mu0);
// 		free(tab);
// 	}
// }
void free_relnkTable(relnkTable* tab){
	if(tab!=NULL){
		if (tab->arr!=NULL){
			int ii;
			for (ii=0; ii<tab->n_a; ii++){
				if (tab->arr[ii] !=NULL){
					int jj;
					for (jj=0; jj<tab->n_def_par; jj++){
						if(tab->arr[ii][jj] != NULL){
							int kk;
							for(kk=0; kk<tab->n_mu0; kk++){
								free_relDat(tab->arr[ii][jj][kk],tab->n_r);
								free(tab->arr[ii][jj][kk]);
							}
							free(tab->arr[ii][jj]);
						}
					}
					free(tab->arr[ii]);
				}
			}
			free(tab->arr);
		}
		free(tab->a);
		free(tab->mu0);
		free(tab);
	}
}



// lpTable* new_lpTable(int n_a,int n_h, int n_rad, int* status){
// 	lpTable* tab = (lpTable*) malloc (sizeof(lpTable));
// 	CHECK_MALLOC_RET_STATUS(tab,status,NULL);

// 	tab->n_a = n_a;
// 	tab->n_h = n_h;
// 	tab->n_rad = n_rad;

// 	tab->a = NULL;

// 	tab->dat = (lpDat**) malloc (sizeof(lpDat*)*tab->n_a);
// 	CHECK_MALLOC_RET_STATUS(tab->dat,status,tab);

// 	int ii;
// 	for (ii=0; ii<tab->n_a; ii++){
// 		tab->dat[ii] = NULL;
// 	}
// 	return tab;
// }

lpnkTable* new_lpnkTable(int n_a, int n_def_par, int n_h, int n_rad, int* status){
	lpnkTable* tab = (lpnkTable*) malloc (sizeof(lpnkTable));
	CHECK_MALLOC_RET_STATUS(tab,status,NULL);

	tab->n_a = n_a;
	tab->n_h = n_h;
	tab->n_rad = n_rad;
	tab->n_def_par = n_def_par;

	tab->a = NULL;
	tab->def_par = NULL;

	tab->dat = (lpDat***) malloc (sizeof(lpDat*)*tab->n_a);
	CHECK_MALLOC_RET_STATUS(tab->dat,status,tab);

	int ii, jj;
	for (ii=0; ii<tab->n_a; ii++){
		tab->dat[ii] = (lpDat**) malloc (sizeof(lpDat*)*tab->n_def_par);
		CHECK_MALLOC_RET_STATUS(tab->dat[ii],status,tab);

		for(jj=0; jj<tab->n_def_par; jj++) {
			tab->dat[ii][jj] = NULL;
		}
	}


	return tab;
}



/* destroy the LP table structure */
void free_lpDat(lpDat* dat,int nh){
	if (dat!=NULL){
		int ii;
		for (ii=0;ii<nh;ii++){
			if (dat->del!=NULL) free(dat->del[ii]);
			if (dat->del_inc!=NULL) free(dat->del_inc[ii]);
			if (dat->intens!=NULL) free(dat->intens[ii]);
		}
		free(dat->del);
		free(dat->del_inc);
		free(dat->intens);

		free(dat->h);
		free(dat->rad);

	}
}

/* destroy the LP table structure */
// void free_lpTable(lpTable* tab){
// 	if (tab!=NULL){
// 		if (tab->dat!=NULL){
// 			int ii;
// 			for (ii=0;ii<tab->n_a;ii++){
// 				free_lpDat(tab->dat[ii],tab->n_h);
// 				free(tab->dat[ii]);
// 			}
// 			free(tab->dat);
// 		}
// 		free(tab->a);
// 		free(tab);
// 	}
// }

void free_lpnkTable(lpnkTable* tab){
	if (tab!=NULL){
		if (tab->dat!=NULL){
			int ii;
			for (ii=0;ii<tab->n_a;ii++){
				if(tab->dat[ii] != NULL) {
					int jj;
					for(jj=0; jj<tab->n_def_par; jj++) {
						free_lpDat(tab->dat[ii][jj],tab->n_h);
						free(tab->dat[ii][jj]);
					} 
				free(tab->dat[ii]);
				}
			}
			free(tab->dat);
		}
		free(tab->a);
		free(tab);
	}
}


rcTable* new_rcTable(int n_a, int n_dp, int n_h, int n_r, int n_intens, int* status) {
	rcTable* tab = (rcTable*) malloc (sizeof(rcTable));
	CHECK_MALLOC_RET_STATUS(tab,status,NULL);

	tab->n_a = n_a;
	tab->n_h = n_h;
	tab->n_r = n_r;
	tab->n_dp = n_dp;
	tab->n_intens = n_intens;

	tab->a = NULL;
	tab->def_par = NULL;
	tab->r_ring = NULL;
	tab->height = NULL;


	tab->rad = (float****) malloc(sizeof(float***)*tab->n_a);
	CHECK_MALLOC_RET_STATUS(tab->rad,status,tab);

	int ii, jj, kk, ll;
	for (ii = 0; ii < tab->n_a; ii++){
		tab->rad[ii] = (float***) malloc(sizeof(float**)*tab->n_dp);
		CHECK_MALLOC_RET_STATUS(tab->rad[ii],status,tab);
		for (jj = 0; jj < tab->n_dp; jj++) {
			tab->rad[ii][jj] = (float**) malloc(sizeof(float*)*tab->n_h);
			CHECK_MALLOC_RET_STATUS(tab->rad[ii][jj],status,tab);
			for (kk = 0; kk < tab->n_h; kk++) {
				tab->rad[ii][jj][kk] = NULL;
			}
		}
	}


	
	tab->intens = (float*****) malloc(sizeof(float****)*tab->n_a);
	CHECK_MALLOC_RET_STATUS(tab->intens,status,tab);

	for (ii = 0; ii < tab->n_a; ii++){
		tab->intens[ii] = (float****) malloc(sizeof(float***)*tab->n_dp);
		CHECK_MALLOC_RET_STATUS(tab->intens[ii],status,tab);
		for (jj = 0; jj < tab->n_dp; jj++) {
			tab->intens[ii][jj] = (float***) malloc(sizeof(float**)*tab->n_h);
			CHECK_MALLOC_RET_STATUS(tab->intens[ii][jj],status,tab);
			for (kk = 0; kk < tab->n_h; kk++) {
				tab->intens[ii][jj][kk] = (float**) malloc(sizeof(float*)*tab->n_r);
				CHECK_MALLOC_RET_STATUS(tab->intens[ii][jj][kk],status,tab);
				for (ll = 0; ll < tab->n_r; ll++) {
					tab->intens[ii][jj][kk][ll] = NULL;
					// tab->dat[ii][jj][kk][ll][mm] = NULL;			
				}
			}
		}
	}

	return tab;
}




// static fitsfile* open_xillver_tab(char* filename, int* status){
fitsfile* open_rc_tab(char* filename, int* status){	
	// printf(" rc file open start\n");
	fitsfile* fptr = NULL;
	char* fullfilename=NULL;
	// get the full filename
	if (asprintf(&fullfilename, "%s/%s", get_relxill_table_path(),filename) == -1){
		RELXILL_ERROR("failed to construct full path the rel table",status);
		return NULL;
	}

	// open the file
	if (fits_open_table(&fptr, fullfilename, READONLY, status)) {
		CHECK_RELXILL_ERROR("opening of the table containing the emissivity profile of the ring-like coronae failed",status);
		printf("    either the full path given (%s) is wrong \n",fullfilename);
		printf("    or you need to download the table ** %s **  from \n",filename);
		printf("    http://www.physics.fudan.edu.cn/tps/people/bambi/Site/RELXILL_NK.html \n");
		return NULL;
	}
	// printf(" rc file open end\n");

	free(fullfilename);

	return fptr;
}




// load_single_rcDat(fname, fptr, tab, ind[0], ind[1]+jj, ind[3]+kk, ll, status);
void load_single_rcDat(char* fname, fitsfile** fptr, rcTable* tab, int ii, int jj, int kk, int nhdu, int* status){

	// printf("%s %d %d %d %d\n", fname, ii, jj, kk, *status);

	assert(tab->rad[ii][jj][kk]==NULL);
	// assert(tab->intens[ii][jj][kk][ll]==NULL);

	// open the fits file if not already open
	if (*fptr==NULL){
		// printf("check open\n");
		*fptr = open_rc_tab(fname,status);
		CHECK_STATUS_VOID(*status);
	}

	// rcDat* dat = new_rcDat(RELRCTABLE_NR, status);
	// CHECK_STATUS_VOID(*status);

	// int rownum = ((ii * tab->n_a + jj) * tab->n_dp + kk) * tab->n_h + ll + 3;
	int rownum = (ii * tab->n_a + jj) * tab->n_dp + kk + 3;
	rownum = nhdu;

	// int _nhdu_ = 2;
	int exttype;
	fits_movabs_hdu(*fptr, rownum, &exttype, status);
	CHECK_STATUS_VOID(*status);

	char * colname_r=NULL, *colname_intens=NULL;
	int colnum_r, colnum_intens;
	long _n_;
	int anynul=0;
	double nullval=0.0;
	LONGLONG nelem = (LONGLONG) _n_;

	float* arr = (float*) malloc (tab->n_intens*sizeof(float));
	CHECK_MALLOC_VOID_STATUS(arr,status);


	if (asprintf(&colname_r, "r_%d", rownum-3) == -1){
		RELXILL_ERROR("failed to construct full path to row #",status);
		return;
	}
	if(fits_get_colnum(*fptr, CASEINSEN, colname_r, &colnum_r, status)) { printf(" *** error *** failed getting column number of %s", colname_r); return; }
	// printf(" colnum_r -- %d, %s\n", colnum_r, colname_r);
	if (fits_get_num_rows(*fptr, &_n_, status)) return;
	if (_n_ != RELRCTABLE_NR){
		// printf(" %d %d\n", _n_, RELTABLE_NA);
		RELXILL_ERROR("inconsistent number of rows in ring-like coronae intensity table",status);
		return;
	}
	anynul=0;
	nullval=0.0;
	nelem = (LONGLONG) _n_;
	fits_read_col(*fptr, TFLOAT, colnum_r, 1, 1, nelem , &nullval, arr, &anynul, status);
	CHECK_STATUS_VOID(*status);

	tab->rad[ii][jj][kk] = arr;
	// free(arr);
	// printf(" rad[%d][%d][%d][0] = %f, rad[%d][%d][%d][99] = %f\n", ii, jj, kk, tab->rad[ii][jj][kk][0], ii, jj, kk, tab->rad[ii][jj][kk][99]);

	int ll;
	for (ll = 0; ll < tab->n_r; ll++) {
		assert(tab->intens[ii][jj][kk][ll]==NULL);

		arr = (float*) malloc (tab->n_intens*sizeof(float));
		CHECK_MALLOC_VOID_STATUS(arr,status);

		if (asprintf(&colname_intens, "intensity_%d_%d", rownum-3, ll) == -1){
			RELXILL_ERROR("failed to construct full path to row #",status);
			return;
		}

		if(fits_get_colnum(*fptr, CASEINSEN, colname_intens, &colnum_intens, status)) { printf(" *** error *** failed getting column number of %s", colname_intens); return; }
		// printf(" colnum_defpar -- %d\n", colnum_defpar);

		if (fits_get_num_rows(*fptr, &_n_, status)) return;
		if (_n_ != RELRCTABLE_NR){
			// printf(" %d %d\n", _n_, RELTABLE_NA);
			RELXILL_ERROR("inconsistent number of rows in ring-like coronae intensity table",status);
			return;
		}

		anynul=0;
		nullval=0.0;
		fits_read_col(*fptr, TFLOAT, colnum_intens, 1, 1, nelem , &nullval, arr, &anynul, status);
		CHECK_STATUS_VOID(*status);

		tab->intens[ii][jj][kk][ll] = arr;

		free(colname_intens);

	}
}






// void check_rcTable_cache(char* fname, rcTable* tab, int* ind, int* status) {
// 	// =2=  check if the necessary spectra are loaded (we only open the file once)
// 	fitsfile* fptr = NULL;
// 	int ii, jj, kk, ll;

//     // ind[0] = ind_a;
//     // ind[1] = ind_dp1;
//     // ind[2] = ind_dp2;
//     // ind[3] = ind_h1;
//     // ind[4] = ind_h2;
//     // ind[5] = ind_r;

// 	for (jj = 0; jj < 2; jj++){
// 		for (kk = 0; kk < 2; kk++){
// 			if (tab->rad[ind[0]][ind[1] + jj][ind[3] + kk] == NULL) {
// 				load_single_rcDat(fname, &fptr, tab, ind[0], ind[1]+jj, ind[3]+kk, status);
// 				CHECK_STATUS_VOID(*status); 
// 			}
// 			if (tab->rad[ind[0]+1][ind[2] + jj][ind[4] + kk] == NULL) {
// 				// load_single_rcDat(fname, &fptr, tab, ind[0]+ii, ind[2+jj], ind[1]+kk, status);
// 				load_single_rcDat(fname, &fptr, tab, ind[0]+1, ind[2]+jj, ind[4]+kk, status);
// 				CHECK_STATUS_VOID(*status);
// 				// for (ll = 0; ll < tab->n_r; ll++){
// 				// 	load_single_rcDat(fname, fptr, tab, ind[0], ind[1]+jj, ind[3]+kk, ll, status);
// 				// 		// load_single_rcDat(filename, fptr, tab, ii, jj, kk, ll, status);
// 				// }

// 				// for (ll = 0; ll < tab->n_r; ll++){
// 				// 	load_single_rcDat(fname, fptr, tab, ind[0]+1, ind[2]+jj, ind[4]+kk, ll, status);

// 					// if (tab->dat[ind[0] + ii][ind[2 + jj]][ind[1] + kk][ll] == NULL) {
// 					// 	load_single_rcDat(fname, &fptr, tab, ind[0]+ii, ind[2+jj], ind[1]+kk, ll, status);
// 					// }
// 					// for (mm=0;mm<tab->n_incl; mm++){
// 					// 	if (tab->dat[ind[0]+ii][ind[1]+jj][ind[2]+kk][ind[3]+ll][mm] == NULL){
// 					// 		load_single_spec(fname, &fptr, tab, ind[0]+ii,ind[1]+jj,ind[2]+kk,ind[3]+ll,mm,status);
// 					// 		CHECK_STATUS_VOID(*status);
// 					// 		// todo: check if we can remove part of the cache
// 					// 	}

// 					// }
// 				// }
// 			}
// 		}
// 	}
// 	if (fptr != NULL) {
// 		if (fits_close_file(fptr,status)) {
// 			RELXILL_ERROR(" *** error closing FITS file", status);
// 		}
// 	}
// 	return ;

// }




// static void free_rcDat(rcDat* dat) {//, int nr){
// 	if (dat!=NULL){
// 		free(dat->rad);
// 		free(dat->intens);
// 	}
// }



void free_rcTable(rcTable* tab){
	if (tab!=NULL){
		free(tab->a);
		// free(tab->height);
		int ii;
		if (tab->def_par != NULL) {
			for (ii = 0; ii < tab->n_a; ii++)
				free(tab->def_par[ii]);
		}
		if (tab->height != NULL) {
			for (ii = 0; ii < tab->n_a; ii++)
				free(tab->height[ii]);
		}

		int jj; int kk; int ll;
		if (tab->rad != NULL) {
			for (ii = 0; ii < tab->n_a; ii++) {
				if (tab->rad[ii] != NULL) {
					for (jj = 0; jj < tab->n_dp; jj++) {
						if (tab->rad[ii][jj] != NULL) {
							for (kk = 0; kk < tab->n_h; jj++) {
								free(tab->rad[ii][jj][kk]);
							}
							free(tab->rad[ii][jj]);
						}
					}
					free(tab->rad[ii]);
				}
			}
			free(tab->rad);
		}
		

		if (tab->intens != NULL) {
			for (ii = 0; ii < tab->n_a; ii ++) {
				if (tab->intens[ii] != 0) {
					for (jj = 0; jj < tab->n_dp; jj++) {
						if (tab->intens[ii][jj] != NULL) {
							for (kk = 0; kk < tab->n_h; kk++) {
								if (tab->intens[ii][jj][kk] != NULL) {
									for (ll = 0; ll < tab->n_r; ll++) {
										free(tab->intens[ii][jj][kk][ll]);
									}
									free(tab->intens[ii][jj][kk]);
								}
							}
							free(tab->intens[ii][jj]);
						}
					}
					free(tab->intens[ii]);	
				}
			}
			free(tab->intens);
		}
		free(tab);
	}
}


/** load the complete relline table */
void init_rc_table(char* filename, rcTable** inp_tab, int* status){
	// printf(" rc 20\n");
	rcTable* tab = (*inp_tab);
	fitsfile* fptr = NULL;
	char *fullfilename;
	// print_version_number(status);
	CHECK_STATUS_VOID(*status);
	// printf(" rc 21\n");


	do{ // Errot handling loop

		fptr = open_rc_tab(filename,status);
		CHECK_STATUS_BREAK(*status);
		// printf(" rc 22\n");

		assert (tab == NULL);

		tab = new_rcTable(RELRCTABLE_NA, RELTABLE_NDEFPAR, RELRCTABLE_NH, RELRCTABLE_NHR, RELRCTABLE_NR, status);
		// printf(" n_a=%d, n_dp=%d, n_h=%d, n_r=%d\n", tab->n_a, tab->n_dp, tab->n_h, tab->n_r);
		// printf(" rc 23\n");

		int _nhdu_ = 2;
		int exttype;
		fits_movabs_hdu(fptr, _nhdu_, &exttype, status);
		CHECK_STATUS_BREAK(*status);
		// printf(" rc 24\n");
		int colnum_a, colnum_defpar;
		long _n_;
		if(fits_get_colnum(fptr, CASEINSEN, "a", &colnum_a, status)) { printf(" error at fits get colnum a"); break; }
		// printf(" colnum_a -- %d\n", colnum_a);

		if(fits_get_colnum(fptr, CASEINSEN, get_defparcolname(status), &colnum_defpar, status)) { printf(" error at fits get colnum defpar"); break; }				// *** need to be checked more
		// printf(" colnum_defpar -- %d\n", colnum_defpar);

		if (fits_get_num_rows(fptr, &_n_, status)) return;

		if (_n_ != RELRCTABLE_NA){
			// printf(" %d %d\n", _n_, RELTABLE_NA);
			RELXILL_ERROR("inconsistent number of rows in spin table",status);
			return;
		}
		// printf(" rc 25\n");

		tab->a = (float*) malloc(tab->n_a * sizeof(float));

		int anynul=0;
		double nullval=0.0;
		LONGLONG nelem = (LONGLONG) _n_;
		fits_read_col(fptr, TFLOAT, colnum_a, 1, 1, nelem ,&nullval,tab->a, &anynul, status);
		CHECK_STATUS_BREAK(*status);

		// printf(" rc 26\n");

		tab->def_par = (float**) malloc( sizeof(float*) * tab->n_a);
		
		int ii, jj, kk, ll;
		for (ii = 0; ii < tab->n_a; ii++){
			tab->def_par[ii] = (float*) malloc( sizeof(float) * tab->n_dp);
		}
		// printf(" rc 27\n");
		anynul=0;
		nullval=0.0;
   		nelem = (LONGLONG) RELTABLE_NDEFPAR;
   		// float vall[RELTABLE_NDEFPAR];

		for (ii=0; ii<tab->n_a;ii++){
			// printf(" %d\n", ii);
        	if(fits_read_col(fptr, TFLOAT, colnum_defpar, ii+1, 1, nelem ,&nullval, tab->def_par[ii], &anynul, status)) return;
        	CHECK_STATUS_BREAK(*status);
        	// printf(" a[%d]=%f, d[%d][0]=%f, d[%d][29]=%f\n", ii, tab->a[ii], ii, tab->def_par[ii][0], ii, tab->def_par[ii][29]);
		}

		tab->height = (float**) malloc(sizeof(float*)*tab->n_a);
		CHECK_MALLOC_VOID_STATUS(tab->height, status);

		double r_h;
		for (ii = 0; ii < tab->n_a; ii++) {
			tab->height[ii] = (float*) malloc(sizeof(float)*tab->n_h);
			CHECK_MALLOC_VOID_STATUS(tab->height, status);

			r_h = 1.75 * (1.0 + sqrt(1 - tab->a[ii] * tab->a[ii]));
			for (kk = 0; kk < tab->n_h; kk++) {
				// tab->height[ii][kk] = (float) kk / 249.0 * (500.0 - r_h) + r_h;
				tab->height[ii][kk] = pow((float) kk / 249.0, 2) * (500.0 - r_h) + r_h;
			}
		}


		tab->r_ring = (float*) malloc(sizeof(float)*tab->n_r);
		CHECK_MALLOC_VOID_STATUS(tab->r_ring, status);

		r_h = 0.5;
		for (ii = 0; ii < tab->n_r; ii++) {
			tab->r_ring[ii] = r_h + 0.2 * (float)ii;
			// printf(" r[%d]=%f\n", ii, tab->r_ring[ii]);
			// r_h += 25. / 15.;
		}



		for (ii = 0; ii < tab->n_a; ii++) {
			for (jj = 0; jj < tab->n_dp; jj++) {
				for (kk = 0; kk < tab->n_h; kk++) {
					// for (ll = 0; ll < tab->n_r; ll++) {
						_nhdu_++;
						load_single_rcDat(filename, &fptr, tab, ii, jj, kk, _nhdu_, status);
						// load_single_rcDat(fname, &fptr, tab, ind[0], ind[1]+jj, ind[3]+kk, status);
						if (*status!=EXIT_SUCCESS){
							// assigne the value
							printf(" error reading the ring coronae fits file at ii=%d, jj=%d, kk=%d\n", ii, jj, kk);
						} 

						// CHECK_STATUS_BREAK(*status);
					// }
				}
			}
		}

		// should be set by previous routine
		assert(tab!=NULL);
		assert(tab->a!=NULL);
		assert(tab->def_par!=NULL);
		assert(tab->height!=NULL);
		assert(tab->r_ring!=NULL);

	} while(0);

	if (*status==EXIT_SUCCESS){
		// assigne the value
		// printf("success\n");
		(*inp_tab) = tab;
	} else {
		free_rcTable(tab);
		// printf("error\n");
	}

	if (fptr != NULL) {
		if (fits_close_file(fptr,status)) {
			RELXILL_ERROR(" *** error closing the ring coronae FITS file", status);
		}
	}

	return;
}
