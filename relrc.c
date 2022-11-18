

#include "relrc.h"

// double d2[9];
rcTable *cached_rc_tab = NULL;
int loaded_mdot_type_rc = -1, loaded_rc_tab = -1;



static void interp_rc_table(double * emis, double * re, int n_r, rcTable * tab, relParam * param, int * status) {
	// int ind_a = ind[0];
	// int ind_dp1 = ind[1];
	// int ind_dp2 = ind[2];
	// int ind_h1 = ind[3];
	// int ind_h2 = ind[4];
	// int ind_r = ind[5];

	int ind_a   = binary_search_float(tab->a,tab->n_a, param->a);
	int ind_r	= binary_search_float(tab->r_ring,tab->n_r, param->r_ring);
	int ind_h1 	= binary_search_float(tab->height[ind_a],tab->n_h, param->height);
	int ind_h2 	= binary_search_float(tab->height[ind_a+1],tab->n_h, param->height);
	int ind_dp1 = binary_search_float(tab->def_par[ind_a], tab->n_dp, param->def_par_unscaled);
	int ind_dp2 = binary_search_float(tab->def_par[ind_a+1], tab->n_dp, param->def_par_unscaled);


	double ifac_a   = (param->a - tab->a[ind_a])/
				   (tab->a[ind_a+1] - tab->a[ind_a]);
	double ifac_h1 = (param->height - tab->height[ind_a][ind_h1])/
				   (tab->height[ind_a][ind_h1+1] - tab->height[ind_a][ind_h1]);
	double ifac_h2 = (param->height - tab->height[ind_a+1][ind_h2])/
				   (tab->height[ind_a+1][ind_h2+1]-tab->height[ind_a+1][ind_h2]);
	double ifac_r = (param->r_ring - tab->r_ring[ind_r])/
				   (tab->r_ring[ind_r+1] - tab->r_ring[ind_r]);
	double ifac_dp1 = (param->def_par_unscaled - tab->def_par[ind_a][ind_dp1])/
				   (tab->def_par[ind_a][ind_dp1+1] - tab->def_par[ind_a][ind_dp1]);
	double ifac_dp2 = (param->def_par_unscaled - tab->def_par[ind_a+1][ind_dp2])/
				   (tab->def_par[ind_a+1][ind_dp2+1] - tab->def_par[ind_a+1][ind_dp2]);


	// double re2[2], intens2[2];
	double intens2[2];
	// double * re_ring = (double*) malloc( sizeof(double) * tab->n_intens);
	// double * intens_ring = (double*) malloc( sizeof(double) * tab->n_intens);
	double re_ring[tab->n_intens];
	double intens_ring[tab->n_intens];
	int ii, jj;
	double t1, t2, t3, t4, t5, t6;
	for (ii = 0; ii < tab->n_intens; ii++) {

		t1 = ifac_h1 * (tab->rad[ind_a][ind_dp1][ind_h1+1][ii] - tab->rad[ind_a][ind_dp1][ind_h1][ii]) + tab->rad[ind_a][ind_dp1][ind_h1][ii];
		t2 = ifac_h1 * (tab->rad[ind_a][ind_dp1+1][ind_h1+1][ii] - tab->rad[ind_a][ind_dp1+1][ind_h1][ii]) + tab->rad[ind_a][ind_dp1+1][ind_h1][ii];

		t3 = ifac_h2 * (tab->rad[ind_a+1][ind_dp2][ind_h2+1][ii] - tab->rad[ind_a+1][ind_dp2][ind_h2][ii]) + tab->rad[ind_a+1][ind_dp2][ind_h2][ii];
		t4 = ifac_h2 * (tab->rad[ind_a+1][ind_dp2+1][ind_h2+1][ii] - tab->rad[ind_a+1][ind_dp2+1][ind_h2][ii]) + tab->rad[ind_a+1][ind_dp2+1][ind_h2][ii];

		t5 = ifac_dp1 * (t2 - t1) + t1;
		t6 = ifac_dp2 * (t4 - t3) + t3;

		re_ring[ii] = ifac_a * (t6 - t5) + t5;

			// if (isnan(re_ring[ii])) {
			// 	printf(" intens radius %d ",re_ring[ii]);
			// }


		for (jj = 0; jj < 2; jj++) {
			
			// t1 = ifac_h1 * (tab->dat[ind_a][ind_dp1][ind_h1+1][ind_r+jj]->rad[ii] - tab->dat[ind_a][ind_dp1][ind_h1][ind_r+jj]->rad[ii]) + tab->dat[ind_a][ind_dp1][ind_h1][ind_r+jj]->rad[ii];
			// t2 = ifac_h1 * (tab->dat[ind_a][ind_dp1+1][ind_h1+1][ind_r+jj]->rad[ii] - tab->dat[ind_a][ind_dp1+1][ind_h1][ind_r+jj]->rad[ii]) + tab->dat[ind_a][ind_dp1+1][ind_h1][ind_r+jj]->rad[ii];

			// t3 = ifac_h2 * (tab->dat[ind_a+1][ind_dp2][ind_h2+1][ind_r+jj]->rad[ii] - tab->dat[ind_a+1][ind_dp2][ind_h2][ind_r+jj]->rad[ii]) + tab->dat[ind_a+1][ind_dp2][ind_h2][ind_r+jj]->rad[ii];
			// t4 = ifac_h2 * (tab->dat[ind_a+1][ind_dp2+1][ind_h2+1][ind_r+jj]->rad[ii] - tab->dat[ind_a+1][ind_dp2+1][ind_h2][ind_r+jj]->rad[ii]) + tab->dat[ind_a+1][ind_dp2+1][ind_h2][ind_r+jj]->rad[ii];

			// t5 = ifac_dp1 * (t2 - t1) + t1;
			// t6 = ifac_dp2 * (t4 - t3) + t3;

			// re2[jj] = ifac_a * (t6 - t5) + t5;

			t1 = ifac_h1 * (fabs(tab->intens[ind_a][ind_dp1][ind_h1+1][ind_r+jj][ii]) - fabs(tab->intens[ind_a][ind_dp1][ind_h1][ind_r+jj][ii])) + fabs(tab->intens[ind_a][ind_dp1][ind_h1][ind_r+jj][ii]);
			t2 = ifac_h1 * (fabs(tab->intens[ind_a][ind_dp1+1][ind_h1+1][ind_r+jj][ii]) - fabs(tab->intens[ind_a][ind_dp1+1][ind_h1][ind_r+jj][ii])) + fabs(tab->intens[ind_a][ind_dp1+1][ind_h1][ind_r+jj][ii]);

			t3 = ifac_h2 * (fabs(tab->intens[ind_a+1][ind_dp2][ind_h2+1][ind_r+jj][ii]) - fabs(tab->intens[ind_a+1][ind_dp2][ind_h2][ind_r+jj][ii])) + fabs(tab->intens[ind_a+1][ind_dp2][ind_h2][ind_r+jj][ii]);
			t4 = ifac_h2 * (fabs(tab->intens[ind_a+1][ind_dp2+1][ind_h2+1][ind_r+jj][ii]) - fabs(tab->intens[ind_a+1][ind_dp2+1][ind_h2][ind_r+jj][ii])) + fabs(tab->intens[ind_a+1][ind_dp2+1][ind_h2][ind_r+jj][ii]);

			t5 = ifac_dp1 * (t2 - t1) + t1;
			t6 = ifac_dp2 * (t4 - t3) + t3;

			intens2[jj] = ifac_a * (t6 - t5) + t5;
 
		}
		// interp_lin_1d(ifac_h[0], dat[0]->intens[ind_h[0]][ii], dat[0]->intens[ind_h[0]+1][ii]);

		intens_ring[ii] = interp_lin_1d(ifac_r, intens2[0], intens2[1]);
		// printf(" re[%d]=%f, intens=%f\n", ii, re_ring[ii], intens_ring[ii]);
		// int ii;
		// for (ii=0; ii<nr; ii++){
			// if (isnan(intens_ring[ii])) {
			// 	printf(" intens %d ",intens_ring[ii]);
			// }


		// re_ring[ii] = interp_lin_1d(ifac_r, re2[0], re2[1]);
	}
	// printf("\n");


	double inter_r;

	// get the extent of the disk (indices are defined such that tab->r[ind+1] <= r < tab->r[ind]
	// printf(" param_rin=%f\n", param->rin);
	int ind_rmin = binary_search(re_ring,tab->n_intens,param->rinp);
	// printf(" rin=%f, ind_rmin=%d\n\n", param->rin, ind_rmin);
	if(ind_rmin == 0) ind_rmin = 1;
	assert(ind_rmin>0);
	int kk=ind_rmin;
	for (ii=n_r-1 ; ii>=0 ;ii--){
		while((re[ii] >= re_ring[kk+1])){
			kk++;
			if (kk>=tab->n_intens-1) { //TODO: construct table such that we don't need this?
				if ( re[ii]-RELTABLE_MAX_R <= 1e-6){
					kk=tab->n_intens-2;
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
		// if (jet_del[kk]/M_PI*180.0<=75.0){
		inter_r = (re[ii]-re_ring[kk])/(re_ring[kk+1]-re_ring[kk]);
		// } else {
		// 	inter_r = (log(re[ii])-log(re_ring[kk]))/
		// 			(log(re_ring[kk+1])-log(re_ring[kk]));
		// }
		// if (isnan(inter_r)) 
		// 	printf(" %d, r=%f, potential frac=%f\n", ii, re[ii], inter_r);


		//  log grid for the intensity (due to the function profile)
		// emis[ii] = interp_log_1d(inter_r, intens_ring[kk], intens_ring[kk+1]);
		emis[ii] = interp_lin_1d(inter_r, intens_ring[kk], intens_ring[kk+1]);
		// if (isnan(emis[ii])) 
		// 	printf(" %d, r=%f, frac=%f, re[%d]=%f, re[%d]=%f, potential emis=%f\n", ii, re[ii], inter_r, kk, re_ring[kk], kk+1, re_ring[kk+1], emis[ii]);
		
		// del_emit[ii] = interp_lin_1d(inter_r, jet_del[kk], jet_del[kk+1]);
		// del_inc[ii] = interp_lin_1d(inter_r, jet_del_inc[kk], jet_del_inc[kk+1]);

	     // multiply by the additional factor gi^gamma (see Dauser et al., 2013)

		 // emis[ii] *= pow(gi_potential_lp_nk(re[ii],param->a,param->height,param->beta,del_emit[ii], param->def_par_unscaled, param->def_par_type), param->gamma);
		// emis[ii] *= pow(gi_potential(re[ii], param), param->gamma);
		// double g2 = gi_potential(re[ii], param);
		double g1 = gfactor(re[ii], param);
		if (isnan(g1)) {
			printf(" !!! NaN value is encountered in redshift function, r[%d]=%f, g1=%f\n", ii, re[ii], g1);
		}

		emis[ii] *= pow(g1, param->gamma);
		// printf(" %f %f\n", g2, g1);
		// emis[ii] *= pow(gfactor(re[ii], param), param->gamma);
		if (isnan(emis[ii])) 
			printf(" !!! NaN value is encountered, r[%d]=%f, emissivity=%f\n", ii, re[ii], emis[ii]);
		// printf(" emis r[%d]=%f, gfac1=%f, gfac2=%f\n", ii, re[ii], g2, g1);
	     // take the beaming of the jet into account (see Dauser et al., 2013)
	     // if (param->beta > 1e-6) {
	     //    emis[ii] *= pow(doppler_factor(del_emit[ii],param->beta),2);
	     // }
	}

	// free(re_ring);
	// free(intens_ring);
	// free(ind);

}

void get_emis_ring_corona(relParam* param, double* emis,
		double* re, int n_r, int* status){

// xill_spec* get_xillver_spectra(xillParam* param, int* status){
	// printf(" rc 1\n");
	// rcTable* tab = NULL;

	if (loaded_rc_tab != param->def_par_type || loaded_mdot_type_rc != param->mdot_type) {
		free_rcTable(cached_rc_tab);
		cached_rc_tab = NULL;
	}


	if (cached_rc_tab == NULL) {
		char* fname = get_rc_filename(status);
		CHECK_STATUS_VOID(*status);
		// printf(" rc 2\n");
		init_rc_table(fname, &cached_rc_tab, status);
		CHECK_STATUS_VOID(*status);
		// printf(" rc 3\n");
		loaded_rc_tab = param->def_par_type;
		loaded_mdot_type_rc = param->mdot_type;
	}
	// init_ring_grid(param, &tab, status);
	// CHECK_STATUS_VOID(*status);


	assert(cached_rc_tab!=NULL);
	// assert(fname!=NULL);
	// printf(" rc 4\n");

	// =1=  get the inidices
	// int* ind = get_xill_indices(param, tab, status);
	// int* ind = get_rc_indices(param, cached_rc_tab, status);
	// CHECK_STATUS_VOID(*status);
	// printf(" rc 5\n");

	// =2=  check if the necessary spectra are loaded (we only open the file once)
	// check_xillTable_cache(fname, tab, ind, status);
	// check_rcTable_cache(fname, tab, ind, status);
	// printf(" rc 6\n");

	// =3= interpolate values

	interp_rc_table(emis, re, n_r, cached_rc_tab, param, status);
	assert(emis!=NULL);
	// printf(" rc 7\n");

	// xill_spec* spec = interp_xill_table(tab,param,ind,status);
	// free(ind);
	// return spec;
	return;
}




void get_emis_disk_corona(relParam* param, double* emis,
		double* re, int n_r, int* status){

// xill_spec* get_xillver_spectra(xillParam* param, int* status){
	// printf(" disk 1\n");
	// rcTable* tab = NULL;

	if (loaded_rc_tab != param->def_par_type || loaded_mdot_type_rc != param->mdot_type) {
		free_rcTable(cached_rc_tab);
		cached_rc_tab = NULL;
	}


	if (cached_rc_tab == NULL) {
		char* fname = get_rc_filename(status);
		CHECK_STATUS_VOID(*status);
		// printf(" disk 2\n");
		init_rc_table(fname, &cached_rc_tab, status);
		CHECK_STATUS_VOID(*status);
		// printf(" disk 3\n");
		loaded_rc_tab = param->def_par_type;
		loaded_mdot_type_rc = param->mdot_type;
	}
	// init_ring_grid(param, &tab, status);
	// CHECK_STATUS_VOID(*status);


	assert(cached_rc_tab!=NULL);
	// assert(fname!=NULL);
	// printf(" disk 4\n");

	// =1=  get the inidices
	// int* ind = get_xill_indices(param, tab, status);
	// int* ind = get_rc_indices(param, cached_rc_tab, status);
	// CHECK_STATUS_VOID(*status);
	// printf(" rc 5\n");

	// =2=  check if the necessary spectra are loaded (we only open the file once)
	// check_xillTable_cache(fname, tab, ind, status);
	// check_rcTable_cache(fname, tab, ind, status);
	// printf(" rc 6\n");

	// =3= interpolate values

	double rs[15], emis_part[N_FRAD];
	int ii, jj, count = 2;
	double rin, rout, r_ring;
	rout = param->r_ring_out;
	rin = param->r_ring;
	// r_ring = param->r_ring;
	interp_rc_table(emis_part, re, n_r, cached_rc_tab, param, status);
	for (jj = 0; jj < N_FRAD; jj++)
		emis[jj] = emis_part[jj];
	param->r_ring = rout;
	interp_rc_table(emis_part, re, n_r, cached_rc_tab, param, status);
	for (jj = 0; jj < N_FRAD; jj++)
		emis[jj] += emis_part[jj];
	// param->r_ring = rin;

	int ind_r1 = binary_search_float(cached_rc_tab->r_ring, cached_rc_tab->n_r, rin);
	int ind_r2 = binary_search_float(cached_rc_tab->r_ring, cached_rc_tab->n_r, rout);

	for (ii = ind_r1+1; ii <= ind_r2; ii++) {
		param->r_ring = cached_rc_tab->r_ring[ii];
		interp_rc_table(emis_part, re, n_r, cached_rc_tab, param, status);
		for (jj = 0; jj < N_FRAD; jj++)
			emis[jj] += emis_part[jj];
		count += 1;
	}
	param->r_ring = rin;
	for (jj = 0; jj < N_FRAD; jj++)
		emis[jj] /= (float)count;

	// double ifac_r = (param->r_ring - cached_rc_tab->r_ring[ind_r])/
	// 			   (cached_rc_tab->r_ring[ind_r+1] - cached_rc_tab->r_ring[ind_r]);	

	// for (ii = 0; ii < 15; ii++) {
	// 	r_ring = (double) ii * (rout - rin) / 14.0 + rin;
	// 	param->r_ring = r_ring;
	// 	// printf(" r_ring - %f\n", r_ring);
	// 	interp_rc_table(emis_part, re, n_r, cached_rc_tab, param, status);
	// 	if (r_ring == rin)
	// 		for (jj = 0; jj < N_FRAD; jj++)
	// 			emis[jj] = emis_part[jj];
	// 	else
	// 		for (jj = 0; jj < N_FRAD; jj++)
	// 			emis[jj] += emis_part[jj];
	// }
	// param->r_ring = rin;

	// interp_rc_table(emis, re, n_r, cached_rc_tab, param, status);
	assert(emis!=NULL);
	// printf(" rc 7\n");

	// xill_spec* spec = interp_xill_table(tab,param,ind,status);
	// free(ind);
	// return spec;
	return;
}
