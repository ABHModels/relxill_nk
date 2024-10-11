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

#include "rellp.h"
#include "relrc.h"

// lpTable* cached_lp_table = NULL;
lpnkTable* cached_lp_table = NULL;

int loaded_lp_nk_file = -1, loaded_mdot_type = -1;  // non-Kerr mods

extern double R_ISCO;
// extern double SumEmis;

/** routine for the broken power law emissivity **/
static void get_emis_bkn(double* emis, double* re,int nr,
		double index1, double index2, double rbr){

	// TODO: make it independent of rbr (maybe 1 at r=1rg?)
	double alpha;

	int ii;
	for (ii=0; ii<nr; ii++){
		alpha = index1;
		if (re[ii] > rbr){
			alpha = index2;
		}
		emis[ii] = pow(re[ii]/rbr,-alpha);
	}

	return;
}

/** routine for the twice broken power law emissivity **/
static void get_emis_bkn2(double* emis, double* re,int nr,
		double index1, double index2, double index3, double rbr1, double rbr2){

	// TODO: make it independent of rbr (maybe 1 at r=1rg?)
	double alpha;

	int ii;
	for (ii=0; ii<nr; ii++){
		alpha = index1;
		if (re[ii] > rbr2){
			alpha = index3;
			emis[ii] = pow(re[ii]/rbr2,-alpha) * pow(rbr2/rbr1, -index2);
		}
		else if (re[ii] > rbr1) {
			alpha = index2;
			emis[ii] = pow(re[ii]/rbr1,-alpha);
		}
		else {
			alpha = index1;
			emis[ii] = pow(re[ii]/rbr1,-alpha);
		}
		// SumEmis += emis[ii]/re[ii];
	}

	return;
}



// static void get_emis_jet_point_source(relParam* param, double* emis, double* del_emit, double* del_inc,
// 		double* re, int n_r, lpnkTable* tab, int ind_a, double ifac_a, int* status){

// 	double defparval = param->def_par_unscaled;


// 	lpDat* dat[2];
// 	int ind_h[2];
// 	double ifac_h[2];

// 	// we set the
// 	int ii;
// 	for (ii=0; ii<2; ii++){
// 		dat[ii] = tab->dat[ind_a+ii];
// 		ind_h[ii] = binary_search_float(dat[ii]->h,tab->n_h,param->height);
// 		ifac_h[ii]   = (param->height-dat[ii]->h[ind_h[ii]])/
// 					   (dat[ii]->h[ind_h[ii]+1]-dat[ii]->h[ind_h[ii]]);

// 		// make sure the incident angle is defined as positive value (otherwise the interpolation
// 		// will create problems / jumps )
// 		int jj; int kk;
// 		for (jj=0; jj<tab->n_h; jj++){
// 			for (kk=0; kk<tab->n_rad; kk++){
// 				dat[ii]->del_inc[jj][kk] = fabs(dat[ii]->del_inc[jj][kk]);
// 				dat[ii]->del[jj][kk] = fabs(dat[ii]->del[jj][kk]);
// 			}
// 		}

// 	}


// 	double jet_rad[tab->n_rad];
// 	double jet_emis[tab->n_rad];
// 	double jet_del[tab->n_rad];
// 	double jet_del_inc[tab->n_rad];

// 	// interpolate everything for the given a-h values on the original grid
// 	for (ii=0; ii<tab->n_rad; ii++){
// 	     // #1: intensity
// 	     jet_emis[ii]  =
// 	          (1.0-ifac_a)*(1.0-ifac_h[0])*dat[0]->intens[ind_h[0]][ii]
// 	          +(1.0-ifac_a)*(ifac_h[0])*dat[0]->intens[ind_h[0]+1][ii]
// 	          +(ifac_a)*(1.0-ifac_h[1])*dat[1]->intens[ind_h[1]][ii]
// 	          +(ifac_a)*(ifac_h[1])*dat[1]->intens[ind_h[1]+1][ii];

// 	     jet_del[ii]  =
// 	          (1.0-ifac_a)*(1.0-ifac_h[0])*dat[0]->del[ind_h[0]][ii]
// 	          +(1.0-ifac_a)*(ifac_h[0])*dat[0]->del[ind_h[0]+1][ii]
// 	          +(ifac_a)*(1.0-ifac_h[1])*dat[1]->del[ind_h[1]][ii]
// 	          +(ifac_a)*(ifac_h[1])*dat[1]->del[ind_h[1]+1][ii];

// 	     jet_del_inc[ii]  =
// 	          (1.0-ifac_a)*(1.0-ifac_h[0])*dat[0]->del_inc[ind_h[0]][ii]
// 	          +(1.0-ifac_a)*(ifac_h[0])*dat[0]->del_inc[ind_h[0]+1][ii]
// 	          +(ifac_a)*(1.0-ifac_h[1])*dat[1]->del_inc[ind_h[1]][ii]
// 	          +(ifac_a)*(ifac_h[1])*dat[1]->del_inc[ind_h[1]+1][ii];

// 	     // #2: r-grid
// 	     jet_rad[ii] = interp_lin_1d(ifac_a,dat[0]->rad[ii],dat[1]->rad[ii]);
// 	 }

// 	// and now rebin it to the given radial grid (TODO: check validity as we differ from the original code here)
// 	double inter_r;

// 	// get the extent of the disk (indices are defined such that tab->r[ind+1] <= r < tab->r[ind]
// //	int ind_rmin = binary_search(jet_rad,tab->n_rad,param->rin);
// 	int ind_rmin = binary_search(jet_rad,tab->n_rad,re[n_r-1]);
// 	assert(ind_rmin>0);
// 	int kk=ind_rmin;
// 	for (ii=n_r-1 ; ii>=0 ;ii--){
// 		while((re[ii] >= jet_rad[kk+1])){
// 			kk++;
// 			if (kk>=tab->n_rad-1) { //TODO: construct table such that we don't need this?
// 				if ( re[ii]-RELTABLE_MAX_R <= 1e-6){
// 					kk=tab->n_rad-2;
// 					break;
// 				} else {
// 					RELXILL_ERROR("interpolation of rel_table on fine radial grid failed due to corrupted grid",status);
// 					printf("   --> radius %.4e ABOVE the maximal possible radius of %.4e \n",
// 							re[ii], RELTABLE_MAX_R);
// 					CHECK_STATUS_VOID(*status);
// 				}
// 			}
// 		}


// 		// for larger angles logarithmic interpolation works slightly better
// 		if (jet_del[kk]/M_PI*180.0<=75.0){
// 			inter_r = (re[ii]-jet_rad[kk])/(jet_rad[kk+1]-jet_rad[kk]);
// 		} else {
// 			inter_r = (log(re[ii])-log(jet_rad[kk]))/
// 					(log(jet_rad[kk+1])-log(jet_rad[kk]));
// 		}

// 		//  log grid for the intensity (due to the function profile)
// 		emis[ii] = interp_log_1d(inter_r, jet_emis[kk], jet_emis[kk+1]);

// 		del_emit[ii] = interp_lin_1d(inter_r, jet_del[kk], jet_del[kk+1]);
// 		del_inc[ii] = interp_lin_1d(inter_r, jet_del_inc[kk], jet_del_inc[kk+1]);

// 	    /** multiply by the additional factor gi^gamma (see Dauser et al., 2013) **/
// 		 emis[ii] *= pow(gi_potential_lp(re[ii],param->a,param->height,param->beta,del_emit[ii]), param->gamma);

// 	     // take the beaming of the jet into account (see Dauser et al., 2013)
// 	     if (param->beta > 1e-6) {
// 	        emis[ii] *= pow(doppler_factor(del_emit[ii],param->beta),2);
// 	     }
// 	}

// }


static void get_emis_jet_point_source(relParam* param, double* emis, double* del_emit, double* del_inc,
		double* re, int n_r, lpnkTable* tab, int ind_a, double ifac_a, int* status){

	// printf(" lamppost interpolation\n");
	double defparval = param->def_par_unscaled;
	// printf(" defparval=%f\n", defparval);
	int ii;
	int idefl=-1, idefr=-1, idl1, idl2, idr1, idr2;
	int idefl1, idefr1, idefl2, idefr2;



    for (ii = 0; ii < RELTABLE_NDEFPAR; ii++) {
		//printf(" ** defpar[%d][%d], defpar[%d][%d] --  %f, %f\n", ind_a, ii, ind_a+1, ii, reltab->def_par[ind_a][ii], reltab->def_par[ind_a+1][ii]);
		if (tab->def_par[ind_a][ii] < defparval) {
			idefl = ii;
        }
		if (tab->def_par[ind_a + 1][ii] < defparval) {
			idefr = ii;
        }
	}
    
    if(idefl>=0 && idefl<(RELTABLE_NDEFPAR-1) && idefr>=0 && idefr<(RELTABLE_NDEFPAR-1)) {
        idl1 = 0; idl2 = 0; idr1 = 0; idr2 = 0;
    }
    else if(idefl<0 && idefr>=0 && idefr <(RELTABLE_NDEFPAR-1)){
        idl1 = 1; idl2 = 0; idr1 = 0; idr2 = 0;
    }
    else if(idefl>=0 && idefl<(RELTABLE_NDEFPAR-1) && idefr<0){
        //printf("Here\n");
        idl1 = 0; idl2 = 0; idr1 = 1; idr2 = 0;
    }
    else if(idefl==(RELTABLE_NDEFPAR-1) && idefr>=0 && idefr<(RELTABLE_NDEFPAR-1)){
        idl1 = 0; idl2 = -1; idr1 = 0; idr2 = 0;
    }
    else if(idefl>=0 && idefl<(RELTABLE_NDEFPAR-1) && idefr==(RELTABLE_NDEFPAR-1)){
        idl1 = 0; idl2 = 0; idr1 = 0; idr2 = -1;
    }
    else{
        printf(" Lamppost error: defpar outside the scope of the grid. idefl = %d, idefr = %d\n",idefl,idefr);
        printf(" Note, for now we have a different range of deformation parameter for LP_NK\n Please, choose a different value of the deformation parameter...\n");
        exit(1);
    }

    idefl1 = idefl + idl1;
    idefl2 = idefl + 1 + idl2;
    idefr1 = idefr + idr1;
    idefr2 = idefr + 1 + idr2;

    double tifl, tifr;


	if(idefl1 != idefl2)
		tifl = (defparval - tab->def_par[ind_a][idefl1]) / (tab->def_par[ind_a][idefl2] - tab->def_par[ind_a][idefl1]);
	else tifl = 0.;
	if(idefr1 != idefr2)
		tifr = (defparval - tab->def_par[ind_a+1][idefr1]) / (tab->def_par[ind_a+1][idefr2] - tab->def_par[ind_a+1][idefr1]);
	else tifr = 0.;

	lpDat* dat[4];
	int ind_h[4];
	double ifac_h[4];

	dat[0] = tab->dat[ind_a][idefl1];
	dat[1] = tab->dat[ind_a][idefl2];
	dat[2] = tab->dat[ind_a+1][idefr1];
	dat[3] = tab->dat[ind_a+1][idefr2];

	// we set the
	// int ii;
	for (ii=0; ii<4; ii++){
		// dat[ii] = tab->dat[ind_a+ii];
		ind_h[ii] = binary_search_float(dat[ii]->h,tab->n_h,param->height);
		ifac_h[ii]   = (param->height-dat[ii]->h[ind_h[ii]])/
					   (dat[ii]->h[ind_h[ii]+1]-dat[ii]->h[ind_h[ii]]);
		// printf(" ih[%d]=%f, ifh[%d]=%f\n", ind_h[ii], ind_h[ii], ii, ifac_h[ii]);
		// in_h[ii] = interp_lin_1d(ifac_h[ii], dat[ii]->h[ind_h[ii]], dat[ii]->h[ind_h[ii]+1]);
	}

	double jet_rad[tab->n_rad];
	double jet_emis[tab->n_rad];
	double jet_del[tab->n_rad];
	double jet_del_inc[tab->n_rad];

	double te1, te2, te3, te4;
	double td1, td2, td3, td4;
	double tdi1, tdi2, tdi3, tdi4;
	// double tr1, tr2, tr3, tr4;
	double te, td, tdi, tr, te12, te34, td12, td34, tdi12, tdi34, tr12, tr34;
	// double tifl, tifr;


	// interpolate everything for the given a-h values on the original grid
	for (ii=0; ii<tab->n_rad; ii++){

		te1 = interp_lin_1d(ifac_h[0], dat[0]->intens[ind_h[0]][ii], dat[0]->intens[ind_h[0]+1][ii]);
		te2 = interp_lin_1d(ifac_h[1], dat[1]->intens[ind_h[1]][ii], dat[1]->intens[ind_h[1]+1][ii]);
		te3 = interp_lin_1d(ifac_h[2], dat[2]->intens[ind_h[2]][ii], dat[2]->intens[ind_h[2]+1][ii]);
		te4 = interp_lin_1d(ifac_h[3], dat[3]->intens[ind_h[3]][ii], dat[3]->intens[ind_h[3]+1][ii]);

		td1 = interp_lin_1d(ifac_h[0], dat[0]->del[ind_h[0]][ii], dat[0]->del[ind_h[0]+1][ii]);
		td2 = interp_lin_1d(ifac_h[1], dat[1]->del[ind_h[1]][ii], dat[1]->del[ind_h[1]+1][ii]);
		td3 = interp_lin_1d(ifac_h[2], dat[2]->del[ind_h[2]][ii], dat[2]->del[ind_h[2]+1][ii]);
		td4 = interp_lin_1d(ifac_h[3], dat[3]->del[ind_h[3]][ii], dat[3]->del[ind_h[3]+1][ii]);

		tdi1 = interp_lin_1d(ifac_h[0], dat[0]->del_inc[ind_h[0]][ii], dat[0]->del_inc[ind_h[0]+1][ii]);
		tdi2 = interp_lin_1d(ifac_h[1], dat[1]->del_inc[ind_h[1]][ii], dat[1]->del_inc[ind_h[1]+1][ii]);
		tdi3 = interp_lin_1d(ifac_h[2], dat[2]->del_inc[ind_h[2]][ii], dat[2]->del_inc[ind_h[2]+1][ii]);
		tdi4 = interp_lin_1d(ifac_h[3], dat[3]->del_inc[ind_h[3]][ii], dat[3]->del_inc[ind_h[3]+1][ii]);

		// printf(" %.10f %.10f %.10f %.10f %f %f %f %f %f %f %f %f\n", te1, te2, te3, te4, td1, td2, td3, td4, tdi1, tdi2, tdi3, tdi4);

		te12 = tifl * (te2 - te1) + te1;
		td12 = tifl * (td2 - td1) + td1;
		tdi12 = tifl * (tdi2 - tdi1) + tdi1;

		te34 = tifr * (te4 - te3) + te3;
		td34 = tifr * (td4 - td3) + td3;
		tdi34 = tifr * (tdi4 - tdi3) + tdi3;

		te = ifac_a * (te34 - te12) + te12;
		td = ifac_a * (td34 - td12) + td12;
		tdi = ifac_a * (tdi34 - tdi12) + tdi12;

		tr12 = tifl * (dat[1]->rad[ii] - dat[0]->rad[ii]) + dat[0]->rad[ii];
		tr34 = tifr * (dat[3]->rad[ii] - dat[2]->rad[ii]) + dat[2]->rad[ii];
		// printf(" %f %f, %f %f %f %f\n", tr12, tr34, dat[0]->rad[ii], dat[1]->rad[ii], dat[2]->rad[ii], dat[3]->rad[ii]);

		tr = ifac_a * (tr34 - tr12) + tr12;

		jet_emis[ii] = te;
		jet_del[ii] = td;
		jet_del_inc[ii] = tdi;
		jet_rad[ii] = tr;
		// printf(" int[%d]=%.10f, del[%d]=%f, deli[%d]=%f, rad[%d]=%f\n\n", ii, te, ii, td, ii, tdi, ii, tr);
	 }

	// and now rebin it to the given radial grid (TODO: check validity as we differ from the original code here)
	double inter_r;

	// get the extent of the disk (indices are defined such that tab->r[ind+1] <= r < tab->r[ind]
	// printf(" param_rin=%f\n", param->rin);
	int ind_rmin = binary_search(jet_rad,tab->n_rad,param->rinp);
	// printf(" rin=%f, ind_rmin=%d\n\n", param->rin, ind_rmin);
	if(ind_rmin == 0) ind_rmin = 1;
	assert(ind_rmin>0);
	int kk=ind_rmin;
	for (ii=n_r-1 ; ii>=0 ;ii--){
		while((re[ii] >= jet_rad[kk+1])){
			kk++;
			if (kk>=tab->n_rad-1) { //TODO: construct table such that we don't need this?
				if ( re[ii]-RELTABLE_MAX_R <= 1e-6){
					kk=tab->n_rad-2;
					break;
				} else {
					RELXILL_ERROR("interpolation of rel_table on fine radial grid failed due to corrupted grid",status);
					printf("   --> radius %.4e ABOVE the maximal possible radius of %.4e \n",
							re[ii], RELTABLE_MAX_R);
					CHECK_STATUS_VOID(*status);
				}
			}
		}


		// for larger angles logarithmic interpolation works slightly better
		if (jet_del[kk]/M_PI*180.0<=75.0){
			inter_r = (re[ii]-jet_rad[kk])/(jet_rad[kk+1]-jet_rad[kk]);
		} else {
			inter_r = (log(re[ii])-log(jet_rad[kk]))/
					(log(jet_rad[kk+1])-log(jet_rad[kk]));
		}

		//  log grid for the intensity (due to the function profile)
		emis[ii] = interp_log_1d(inter_r, jet_emis[kk], jet_emis[kk+1]);

		if (isnan(emis[ii])) 
			emis[ii] = interp_lin_1d(inter_r, jet_emis[kk], jet_emis[kk+1]);
		
		
		del_emit[ii] = interp_lin_1d(inter_r, jet_del[kk], jet_del[kk+1]);
		del_inc[ii] = interp_lin_1d(inter_r, jet_del_inc[kk], jet_del_inc[kk+1]);

	    // multiply by the additional factor gi^gamma (see Dauser et al., 2013)
		// emis[ii] *= pow(gi_potential_lp_nk(re[ii],param->a,param->height,param->beta,del_emit[ii], param->def_par_unscaled, param->def_par_type), param->gamma);
		// emis[ii] *= pow(gi_potential(re[ii], param), param->gamma);
		emis[ii] *= pow(gfactor(re[ii], param), param->gamma);
		// printf(" emis gfac1=%f, gfac2=%f\n", gi_potential(re[ii], param), gfactor(re[ii], param));


	    // take the beaming of the jet into account (see Dauser et al., 2013)
	    if (param->beta > 1e-6) {
	       emis[ii] *= pow(doppler_factor(del_emit[ii],param->beta),2);
	    }
	}

}





/** routine to calculate the emissivity in the lamp post geometry**/
// void get_emis_jet(relParam* param, double* emis, double* del_emit, double* del_inc,
// 		double* re, int n_r, int* status){

// 	CHECK_STATUS_VOID(*status);

// 	if (cached_lp_table==NULL){
// 		read_lp_table(LPTABLE_FILENAME, &cached_lp_table,status);
// 		CHECK_STATUS_VOID(*status);
// 	}

// 	// calculate the emissivity
// 	assert(emis!=NULL);

// 	int ind_a   = binary_search_float(cached_lp_table->a,cached_lp_table->n_a,param->a);
// 	double ifac_a   = (param->a-cached_lp_table->a[ind_a])/
// 				   (cached_lp_table->a[ind_a+1]-cached_lp_table->a[ind_a]);


// 	// TODO: make it possible for an extended jet
// 	get_emis_jet_point_source(param, emis, del_emit, del_inc, re, n_r,
// 			cached_lp_table, ind_a, ifac_a,status);

// 	return;
// }


/** routine to calculate the emissivity in the lamp post geometry**/ // non-Kerr mods
void get_emis_jet(relParam* param, double* emis, double* del_emit, double* del_inc,
		double* re, int n_r, int* status){
	// printf(" get_emis_jet start\n");
	// printf(" lp1\n");
	// printf(" loaded_lp_nk_file - %d\n", loaded_lp_nk_file);

	if (loaded_lp_nk_file != param->def_par_type || loaded_mdot_type != param->mdot_type) {
		free_lpnkTable(cached_lp_table);
		cached_lp_table = NULL;
	}

	if (cached_lp_table==NULL){// || check_caching_defpartype()){
		// read_lp_table(LPTABLE_FILENAME, &cached_lp_table,status);
		// printf(" lamppost file read\n");
		read_lp_nk_table(get_lp_filename(status), &cached_lp_table, status);
		loaded_lp_nk_file = param->def_par_type;
		loaded_mdot_type = param->mdot_type;
		CHECK_STATUS_VOID(*status);
	}
	// printf(" lp2\n");
	// calculate the emissivity
	assert(emis!=NULL);

	int ind_a   = binary_search_float(cached_lp_table->a,cached_lp_table->n_a,param->a);
	double ifac_a   = (param->a-cached_lp_table->a[ind_a])/
				   (cached_lp_table->a[ind_a+1]-cached_lp_table->a[ind_a]);

	// TODO: make it possible for an extended jet
	get_emis_jet_point_source(param, emis, del_emit, del_inc, re, n_r,
			cached_lp_table, ind_a, ifac_a,status);
	// printf(" get_emis_jet end\n");
	return;
}



// calculate the angles of emission from the primary source to git Rin and Rout
void get_ad_del_lim(relParam* param, relSysPar* sysPar, int* status) {
	int nr_lim = 2;
	double del_emit[2];
	double del_dummy[2];
	double emis_dummy[2];
	double rad[2] = { AD_ROUT_MAX, R_ISCO };  // non-Kerr mods
	// get the primary source emission angle for the simulated inner and out edge of the disk
	get_emis_jet(param, emis_dummy, del_emit, del_dummy, rad, nr_lim, status);
	sysPar->del_ad_rmax = del_emit[0];
	sysPar->del_ad_risco = del_emit[1];

	if (*status!=EXIT_SUCCESS){
		printf(" *** failed calculating the primary source photon emission angles for Rin=%.2e and Rmax=%.2e",rad[1], rad[0]);
	}
}


void calc_emis_profile(double* emis, double* del_emit, double* del_inc, double* re, int nr, relParam* param, int* status){

	CHECK_STATUS_VOID(*status);

	double invalid_angle = -1.0;

	/**  *** Broken Power Law Emissivity ***  **/
	if (param->emis_type==EMIS_TYPE_BKN){

		// if (param->emis3 != PARAM_DEFAULT && param->rbr2 != PARAM_DEFAULT)
			get_emis_bkn2(emis, re, nr,
				param->emis1,param->emis2,param->emis3,param->rbr,param->rbr2);
		// else
		// 	get_emis_bkn(emis, re, nr,
		// 		param->emis1,param->emis2,param->rbr);
		// set the angles in this case to a default value
		int ii;
		for (ii=0; ii<nr; ii++){
			del_emit[ii] = invalid_angle;
			del_inc[ii] = invalid_angle;
		}

	/**  *** Lamp Post Emissivity ***  **/
	} else if (param->emis_type==EMIS_TYPE_LP){


		// if (redo_get_emis_lp(param, cached_param, sysPar)){  // XXX Currently we always re-do it if we get here
		get_emis_jet(param,emis, del_emit, del_inc,
				re, nr, status);
	
	} else if (param->emis_type == EMIS_TYPE_RING) {

		// printf(" ring calc start\n");

		get_emis_ring_corona(param, emis, re, nr, status);

		int ii;
		for (ii=0; ii<nr; ii++){
			del_emit[ii] = invalid_angle;
			del_inc[ii] = invalid_angle;
		}
	} else if (param->emis_type == EMIS_TYPE_DISK) {

		// printf(" disk calc start\n");

		// printf(" 1 r_ring - %f\n", param->r_ring);

		get_emis_disk_corona(param, emis, re, nr, status);

		// printf(" 2 r_ring - %f\n", param->r_ring);
		int ii;
		for (ii=0; ii<nr; ii++){
			del_emit[ii] = invalid_angle;
			del_inc[ii] = invalid_angle;
		}
	} else {

		RELXILL_ERROR(" calculation of emissivity profile not possible \n",status);
		printf("   -> emis_type=%i not known \n",param->emis_type);
		return;
	}


	if (*status!=EXIT_SUCCESS){
		RELXILL_ERROR("calculating the emissivity profile failed due to previous error", status);
	}

		// int ii;
		// for (ii=0; ii<nr; ii++){
		// 	if (isnan(emis[ii])) {
		// 		printf(" emis %d\n", emis[ii]);
		// 	}
		// 	// if (isnan(del_emit[ii])) {
		// 	// 	printf(" del_emit %d\n", del_emit[ii]);
		// 	// }
		// 	// if (isnan(del_inc[ii])) {
		// 	// 	printf(" del_inc %d\n", del_inc[ii]);
		// 	// }
		// }


	return;
}

void free_cached_lpTable(void){
	free_lpnkTable(cached_lp_table);
}
