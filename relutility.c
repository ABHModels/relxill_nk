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

#include "relutility.h"


/** linear interpolation in 1 dimension **/
double interp_lin_1d(double ifac_r, double rlo, double rhi){
	return ifac_r*rhi + (1.0-ifac_r)*rlo;
}

double interp_log_1d(double ifac_r, double rlo, double rhi){
	return exp(ifac_r*log(rhi) + (1.0-ifac_r)*log(rlo));
}

/** linear interpolation in 2 dimensions **/
double interp_lin_2d(double ifac1, double ifac2, double r11, double r12, double r21, double r22){
	return (1.0 - ifac1) * (1.0 - ifac2) * r11 +
		   (ifac1)       * (1.0 - ifac2) * r12 +
		   (1.0 - ifac1) * (ifac2)       * r21 +
		   (ifac1)       * (ifac2)       * r22;
}

double interp_lin_2d_float(double ifac1, double ifac2, float r11, float r12, float r21, float r22){
	return (1.0 - ifac1) * (1.0 - ifac2) * r11 +
		   (ifac1)       * (1.0 - ifac2) * r12 +
		   (1.0 - ifac1) * (ifac2)       * r21 +
		   (ifac1)       * (ifac2)       * r22;
}

// non-Kerr modifications
double interp_lin_3d(double v1, double v2, double v3, double v4, double v5, double v6, double v7, double v8, double d[]) {
	// double d, double d1, double d2, double d3, double d4, double d5, double d6, double d7, double d8){
	return (v1*d[1] + v2*d[2] + v3*d[3] + v4*d[4] + v5*d[5] + v6*d[6] + v7*d[7] + v8*d[8]) / d[0];
	// return res;
}


void relxill_error(const char* const func, const char* const msg, int* status){
	*status = EXIT_FAILURE;
	printf(" *** error in relxill (%s): %s!\n", func, msg);
}


void relxill_warning(const char* const msg){
	printf(" *** warning from relxill: %s!\n", msg);
}


int is_xill_model(int model_type){
	if ((model_type == MOD_TYPE_XILLVERDENS) || (model_type == MOD_TYPE_XILLVER)){
		return 1;
	} else {
		return 0;
	}
}

int is_xill_model2(int model_type){
	if ((model_type == MOD_TYPE_XILLVERDENS) || (model_type == MOD_TYPE_XILLVER) || (model_type == MOD_TYPE_RELXILLION_NK2) || (model_type == MOD_TYPE_RELXILL2)){
	// if ((model_type == MOD_TYPE_XILLVERDENS) || (model_type == MOD_TYPE_XILLVER)){
		return 1;
	} else {
		return 0;
	}
}

int is_densgrad_model(int rmodel_type, int xmodel_type, int ion_grad_type ){
	if (xmodel_type == MOD_TYPE_RELXILLDENS && rmodel_type == MOD_TYPE_RELXILLDENSGRAD_NK && ion_grad_type!=ION_GRAD_TYPE_CONST){
	// if ((model_type == MOD_TYPE_XILLVERDENS) || (model_type == MOD_TYPE_XILLVER)){
		return 1;
	} else {
		return 0;
	}
}



// ion-gradient model, which is not set to constant ionization
//  - in case ion_grad_type=constant, it is working as a normal model
int is_iongrad_model(int model_type, int ion_grad_type ){
	// if ( (model_type == MOD_TYPE_RELXILLLPION) &&  ( ion_grad_type!=ION_GRAD_TYPE_CONST ) ) {
	if ( (model_type == MOD_TYPE_RELXILLLPION || MOD_TYPE_RELXILLION_NK || MOD_TYPE_RELXILLION_NK2 || MOD_TYPE_RELXILLDENSGRAD_NK) &&  ( ion_grad_type!=ION_GRAD_TYPE_CONST ) ) {
		return 1;
	} else {
		return 0;
	}
}

// #define MOD_TYPE_RELXILLLPION -21
// #define MOD_TYPE_RELXILLION_NK -30
// #define MOD_TYPE_RELXILLION_NK2 -31


int is_relxillion_model(int model_type, int ion_grad_type ){  // relxill_nk variable ionization
	if ( (model_type == MOD_TYPE_RELXILLION_NK || model_type == MOD_TYPE_RELXILLION_NK2) &&  ( ion_grad_type!=ION_GRAD_TYPE_CONST ) ) {
		return 1;
	} else {
		return 0;
	}
}


/** calculate the gravitational redshift **/
// double grav_redshift(relParam* param){
// 	if (param->emis_type==EMIS_TYPE_LP){
// 		return 1.0 / sqrt( 1.0 - 2*param->height /
// 				(param->height*param->height + param->a*param->a))-1.0;
// 	} else {
// 		// important: without a geometrical assumption no grav. redshift can be calculated
// 		return 0.0;
// 	}
// }

// /** calculate the gravitational redshift **/
double grav_redshift(relParam* param){
	if (param->emis_type==EMIS_TYPE_LP || param->emis_type==EMIS_TYPE_RING || param->emis_type==EMIS_TYPE_DISK){
		// double tmp = 1.0 / sqrt( 1.0 - 2*param->height /
		// 		(param->height*param->height + param->a*param->a))-1.0;
		double tmp2 = 1.0 / get_g_inf(param) - 1.0;
		// printf(" grav_redshift %f %f\n", tmp, tmp2);
		return tmp2;

	} else {
		// important: without a geometrical assumption no grav. redshift can be calculated
		return 0.0;
	}
}


static double relat_abberation(double del,double beta) {
	return acos( (cos(del)-beta) / (1-beta*cos(del)));
}


/** calculate the reflection fraction as defined in Dauser+2016 **/
lpReflFrac* calc_refl_frac(relSysPar* sysPar, relParam* param, int* status){

	// in case there is no relativity information, the refl_frac is 1
	if (param==NULL){
		printf(" *** Warning: can not calculate reflection fraction as no relat. parameters are given \n");
		return NULL;
	}

	/** important: get the radial values for which the RELLINE is calculated
	 *             should be Rin=r_isco & Rout=1000rg  **/

	// get the angle emitted in the rest-frame of the primary source, which hits the inner and outer edge of the disk
	double del_bh  = sysPar->del_emit[inv_binary_search(sysPar->re, sysPar->nr, param->rinp)];
	double del_ad = sysPar->del_emit[inv_binary_search(sysPar->re, sysPar->nr, param->rout)];

	/** calculate the coordinate transformation / relat abberation
	 *   - an observer on the accretion disk sees the rays from
	 *     del_bh up to del_ad
	 *   - for the reflection fraction we therefore need to convert from
	 *     the moving source (which the disk observer sees) into the
	 *     local frame
	 *   -> therefore we need to calculate the abberation fro -beta
	 */
	if (param->beta>1e-6) {
	     del_bh = relat_abberation(del_bh, -1.*param->beta);
	     del_ad = relat_abberation(del_ad, -1.*param->beta);
	}

	lpReflFrac* str = (lpReflFrac*) malloc (sizeof(lpReflFrac));
	CHECK_MALLOC_RET_STATUS(str,status,NULL);

	str->f_bh  = 0.5*(1.0 - cos(del_bh));
	str->f_ad  = 0.5*(cos(del_bh) - cos(del_ad));
	/** photons are not allowed to cross the disk
	 *  (so they only reach infinity if they don't hit the disk plane) */
	str->f_inf = 0.5*(1.0 + cos(sysPar->del_ad_rmax));

	/** fraction of photons which would hit the maximally
	 *  simulated accretion disk. Do not change this with BETA (!) */
	str->f_ad_norm = 0.5*(cos(sysPar->del_ad_risco) - cos(sysPar->del_ad_rmax));

	// photons are not allowed to crosstalk the disk plane
	if (str->f_inf > 0.5){
	  str->f_inf = 0.5;
	}
	
	str->refl_frac = str->f_ad/str->f_inf;
	str->refl_frac_norm = str->f_ad_norm/str->f_inf;

	return str;
}

void check_relxill_error(const char* const func, const char* const msg, int* status){
	if (*status!=EXIT_SUCCESS){
		*status = EXIT_FAILURE;
		printf(" *** error in relxill (%s): %s!\n", func, msg);
	}
}

void get_version_number(char** vstr, int* status){

	if (strcmp(version_dev,"")==0){
		if (asprintf(vstr, "%i.%i.%i", version_major, version_minor, version_build) == -1){
			RELXILL_ERROR("failed to get version number",status);
		}
	} else {
		if (asprintf(vstr, "%i.%i.%i%s", version_major, version_minor, version_build, version_dev) == -1){
			RELXILL_ERROR("failed to get version number",status);
		}
	}
}

void get_version_number_nk(char** vstr, int* status){

	if (strcmp(version_dev,"")==0){
		if (asprintf(vstr, "%i.%i.%i", version_major_nk, version_minor_nk, version_build_nk) == -1){
			RELXILL_ERROR("failed to get version number",status);
		}
	} else {
		if (asprintf(vstr, "%i.%i.%i%s", version_major_nk, version_minor_nk, version_build_nk, version_dev_nk) == -1){
			RELXILL_ERROR("failed to get version number",status);
		}
	}
}


/**  FLOAT search for value "val" in array "arr" (sorted ASCENDING!) with length n and
 	 return bin k for which arr[k]<=val<arr[k+1] **/
int binary_search_float(float* arr,int n,float val){

	int klo=0;
	int khi=n-1;
	int k=-1;
	while ( (khi-klo) > 1 ){
		k=(khi+klo)/2;
		if(arr[k]>val){
			khi=k;
		} else {
			klo=k;
		}
	}
	return klo;
}

/**  search for value "val" in array "arr" (sorted ASCENDING!) with length n and
 	 return bin k for which arr[k]<=val<arr[k+1] **/
int binary_search(double* arr,int n,double val){

	int klo=0;
	int khi=n-1;
	int k=-1;
	while ( (khi-klo) > 1 ){
		k=(khi+klo)/2;
		if(arr[k]>val){
			khi=k;
		} else {
			klo=k;
		}
	}
	return klo;
}


/**  FLOAT search for value "val" in array "arr" (sorted DESCENDING!) with length n and
 	 return bin k for which arr[k]<=val<arr[k+1] **/
int inv_binary_search_float(float* arr,int n,float val){

	int klo=0;
	int khi=n-1;
	int k=-1;
	while ( (khi-klo) > 1 ){
		k=(khi+klo)/2;
		if(arr[k]<val){
			khi=k;
		} else {
			klo=k;
		}
	}
	return klo;
}


/**  search for value "val" in array "arr" (sorted DESCENDING!) with length n and
 	 return bin k for which arr[k]<=val<arr[k+1] **/
int inv_binary_search(double* arr,int n,double val){

	int klo=0;
	int khi=n-1;
	int k=-1;
	while ( (khi-klo) > 1 ){
		k=(khi+klo)/2;
		if(arr[k]<val){
			khi=k;
		} else {
			klo=k;
		}
	}
	return klo;
}

/** test if it is a relxill flavour model **/
int is_relxill_model(int model_type){
	if (model_type < 0 && model_type != MOD_TYPE_RELXILLOLD) {
		return 1;
	} else {
		return 0;
	}
}

int fifty_zones(int model_type) {
	// if (model_type < 0 && model_type != MOD_TYPE_RELXILLOLD) {
	if (model_type == MOD_TYPE_RELXILL || model_type == MOD_TYPE_RELXILLDENS || model_type == MOD_TYPE_RELXILLION_NK || model_type == MOD_TYPE_RELXILLION_NK2) {
		return 1;
	} else {
		return 0;
	}	
}

// #define MOD_TYPE_RELXILL -1
// #define MOD_TYPE_RELXILLDENS -10
// #define MOD_TYPE_RELXILLION_NK -30
// #define MOD_TYPE_RELXILLION_NK2 -31

/** trapez integration around a single bin
 *  caveat: only returns half of the full 2*PI*r*dr due to computational speed**/
double trapez_integ_single(double* re, int ii, int nr){
	double dr;
	// dr is defined such that the full disk is covered once, with NO overlapping bins
	if (ii==0){
		dr = 0.5*(re[ii] - re[ii+1]) ;
	} else if (ii==nr-1){
        dr =  0.5*(re[ii-1] - re[ii]) ;
	} else {
        dr = 0.5*( re[ii-1] - re[ii+1]);
	}
	return re[ii]*dr*M_PI;
}


/** convert gstar to energy */
double gstar2ener(double g, double gmin, double gmax, double ener){
	return (g*(gmax-gmin) + gmin)*ener;
}

/** get a radial grid on the accretion disk in order to calculate a relline for each zone **/
void get_rzone_grid(relParam * param, double rmin, double rmax, double* rgrid, int nzones, double h, int* status){
// void get_rzone_grid(double rmin, double rmax, double* rgrid, int nzones, double h, int* status){

	if (nzones==1){
		rgrid[0] = rmin;
		rgrid[1] = rmax;
		// printf(" rzone_grid 1\n");
	} else if (h < 1e-8) {
		get_log_grid(rgrid,nzones+1,rmin,rmax);
		// double rlo = rmin;
		// double rhi = rmax;
		// // // printf("ionization rgrid\n");
		// int ii;
	 //    for (ii=0; ii<nzones+1; ii++) {
		// 	rgrid[ii] = 1.0*(ii) / (nzones) * ( 1.0/rhi - 1.0/rlo) + 1.0/rlo;
		// 	rgrid[ii] = fabs(1.0/rgrid[ii]);

		// 	// rgrid[ii] = pow(1.0*(ii) / (int)(nzones), 1.5) * (rhi - rlo) + rlo;
	 //    	// printf("%f\n", rgrid[ii]);

		// }
		// printf(" rzone_grid 2\n");

	} else {

		double r_transition = rmin;
		int indr = 0;

		// if h > rmin we choose a log grid for r<h
		if (h > rmin){

			r_transition = h;

			get_log_grid(rgrid,nzones+1,rmin,rmax);
			indr = binary_search(rgrid,nzones+1,r_transition);

			r_transition = rgrid[indr];

		}

		if (indr < nzones ){

			double rlo = r_transition;
			double rhi = rmax; // rgrid[nzones];
			// add 1/r for larger radii
            int ii;
            for (ii=indr; ii<nzones+1; ii++){
				rgrid[ii] = 1.0*(ii-indr) / (nzones-indr) * ( 1.0/rhi - 1.0/rlo) + 1.0/rlo;
				rgrid[ii] = fabs(1.0/rgrid[ii]);
			}

		}
		// printf(" rzone_grid 3\n");

	}
	return;
}

void get_relxillion_rzone_grid(double rmin, double rmax, double* rgrid, int nzones, int* status){ // relxill_nk variable ionization
	// printf(" xillion grid\n");
	// get_log_grid(rgrid,nzones+1,rmin,rmax);
	double rlo = rmin;
	double rhi = rmax;
	// printf("ionization rgrid\n");
		int ii;
	    for (ii=0; ii<nzones+1; ii++){
			rgrid[ii] = 1.0*(ii) / (int)(nzones) * ( 1.0/rhi - 1.0/rlo) + 1.0/rlo;
			rgrid[ii] = fabs(1.0/rgrid[ii]);

			// rgrid[ii] = pow(1.0*(ii) / (int)(nzones), 1.5) * (rhi - rlo) + rlo;
	    	// printf("%f\n", rgrid[ii]);

		}
	return;
}

/** get the relxill table path (dynamically from env variable)  **/
char* get_relxill_table_path( void ){
	char* path;
	path = getenv("RELXILL_TABLE_PATH");
	if (path!=NULL){
		return path;
	} else {
		return RELXILL_TABLE_PATH;
	}
}


/** check if we are currently debugging the model **/
int is_debug_run( void ){
	// return 1;
	char* env;
	env = getenv("DEBUG_RELXILL");
	if (env != NULL){
		int debug = atof(env);
		if (debug == 1){
			return 1;
		}
	}
	return 0;
}


/** check if we should return the relline/relconv physical norm from ENV **/
int do_not_normalize_relline( void ){
	char* env;
	env = getenv("RELLINE_PHYSICAL_NORM");
	if (env != NULL){
		int phys_norm = atof(env);
		if (phys_norm == 1){
			return 1;
		}
	}
	return 0;
}


/* get a logarithmic grid from emin to emax with n_ener bins  */
void get_log_grid(double* ener, int n_ener, double emin, double emax){
	int ii;
	for (ii=0; ii<n_ener; ii++){
		ener[ii] = 1.0*ii / (n_ener-1) * ( log(emax) - log(emin)) + log(emin);
		ener[ii] = exp(ener[ii]);
	}
}

/* get a logarithmic grid from emin to emax with n_ener bins  */
void get_rgrid(double* ener, int n_ener, double emin, double emax){
	int ii;
	for (ii=0; ii<n_ener; ii++){
		ener[ii] = 1.0*ii / (n_ener-1) * ( 1.0/emax - 1.0/emin) + 1.0/emin;
		ener[ii] = fabs(1.0/ener[ii]);
	}
}


/* get a logarithmic grid from emin to emax with n_ener bins  */
void get_lin_grid(double* ener, int n_ener, double emin, double emax){
	int ii;
	for (ii=0; ii<n_ener; ii++){
		ener[ii] = 1.0*ii / (n_ener-1) * ( emax - emin) + emin;
	}
}


/* get RMS (ISCO) for the Kerr Case */
double kerr_rms(double a){
	//	 accounts for negative spin
  double sign = 1.0;
  if (a<0) {
     sign = -1.0;
  }

  double Z1 = 1.0+pow(1.0-a*a,1.0/3.0)*(pow(1.0+a,1.0/3.0)+pow(1.0-a,1.0/3.0));
  double Z2=sqrt((3.0*a*a)+(Z1*Z1));

  return 3.0+Z2-sign*sqrt((3.0-Z1)*(3.0+Z1+(2*Z2)));
}

/* get the rplus value (size if the black hole event horizon */
double kerr_rplus(double a){
	return 1 + sqrt (1 - a*a);
}

/** calculate the doppler factor for a moving primary source **/
double doppler_factor(double del, double bet) {
	return sqrt(1.0 - bet*bet) / (1.0 + bet*cos(del));
}


/** calculates g = E/E_i in the lamp post geometry (see, e.g., 27 in Dauser et al., 2013, MNRAS) **/
double gi_potential_lp(double r, double a, double h, double bet, double del){

	/** ! calculates g = E/E_i in the lamp post geometry
	  ! (see, e.g., page 48, Diploma Thesis, Thomas Dauser) **/
	double ut_d = ((r*sqrt(r)+a)/(sqrt(r)*sqrt(r*r -3*r + 2*a*sqrt(r))));
	double ut_h = sqrt((h*h + a*a)/(h*h - 2*h + a*a));

	double gi = ut_d/ut_h;

	// check if we need to calculate the additional factor for the velocity
	if (fabs(bet) < 1e-6){
		return gi;
	}

	double gam = 1.0/sqrt(1.0-bet*bet);

	// get the sign for the equation
	double sign = 1.0;
	if (del > M_PI/2) {
		sign = -1.0;
	}

	double delta_eq = h*h - 2*h + a*a;
	double q2 = (pow(sin(del),2))*(  pow((h*h + a*a),2) / delta_eq ) - a*a;

	double beta_fac = sqrt(  pow((h*h + a*a ),2) - delta_eq*(q2 + a*a) );
	beta_fac = gam*(1.0 + sign*beta_fac / (h*h + a*2) *bet);

	return gi / beta_fac;
}

double get_g_inf(relParam * param) {
	
	double epsi3 = 0, a13 = 0, a22 = 0, a52 = 0, spin, r0, th0, gtt_src;

	if (param->def_par_type == 1)
	{
		a13 = param->def_par_unscaled;
		// a22 = 0.0;
	}
	if (param->def_par_type == 2)
	{
		a22 = param->def_par_unscaled;
		// a13 = 0.0;
	}

	if (param->emis_type == EMIS_TYPE_LP) {
		r0 = param->height;
		th0 = 0.0;
	}
	else if (param->emis_type == EMIS_TYPE_RING) {
		r0 = sqrt(param->height * param->height + param->r_ring * param->r_ring);
		th0 = atan(param->r_ring / param->height);	
	}
	else if (param->emis_type == EMIS_TYPE_DISK) {
		double rring2 = (param->r_ring + param->r_ring_out) / 2.0;
		r0 = sqrt(param->height * param->height + rring2 * rring2);
		th0 = atan(rring2 / param->height);	
	}
	else {
		printf(" The redshift factor cannot be calculated for the given emission type\n");
		exit(1);
	}
	spin = param->a;

	double t2 = pow(r0,3);
	double t3 = pow(spin,2);
	double t11 = pow(r0,2);
	double t14 = a22 + t11;
	double t15 = sin(th0);
	double t16 = pow(t15,2);
	
	gtt_src = (r0*(pow(t14,2)*t16*t3 - pow(r0,4)*((-2 + r0)*r0 + t3))*(epsi3 + t2 + r0*t3*pow(cos(th0),2)))/pow(-(r0*t14*t16*t3) + (a13 + t2)*(t11 + t3),2);	

	return sqrt(fabs(gtt_src));
}

/** calculates g = E/E_i in the lamp post geometry (see, e.g., Abdikamalov et al., 2019, ApJ) **/
double gi_potential_lp_nk(double r, double a, double h, double bet, double del, double def_par, int def_par_type) {
	double Omega, gi_lp;
	double th, g_tt, g_tp, g_pp;
	double dgttdr, dgtpdr, dgrrdr, dgppdr, dgththdr;
	double spin = a;
	double a13 = 0.0, a22 = 0.0, a52 = 0.0, eps3 = 0.0;
	double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, /*t16*/ t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, /*t34*/ t35, t36, t37, t38, t39, t40, /*t41*/ t42, t43, t44, t45, t46, t47, t48;
	double r0, th0;
	double g_tt_src;

	if (def_par_type == 1)
	{
		a13 = def_par;
		// a22 = 0.0;
	}
	if (def_par_type == 2)
	{
		a22 = def_par;
		// a13 = 0.0;
	}

	r0 = h;
	th0 = 0.0;
	g_tt_src = ((-pow(r0,2) - eps3/r0 - pow(a,2)*pow(cos(th0),2))*(pow(a,2) - 2*r0 + pow(r0,2) - pow(a,2)*pow(1 + a22/pow(r0,2),2)*pow(sin(th0),2)))/
		pow((pow(a,2) + pow(r0,2))*(1 + a13/pow(r0,3)) - pow(a,2)*(1 + a22/pow(r0,2))*pow(sin(th0),2),2);

	spin = a;
	th = M_PI/2.0;
	
	t1 = cos(th);
	t2 = r * r;
	t3 = pow(t2, 0.2e1);
	t4 = t2 * t3;
	t5 = r * t3;
	t6 = r * t2;
	t7 = spin * spin;
	t8 = t7 * pow(t1, 0.2e1);
	t9 = (t8 + t2) * r + eps3;
	t10 = a22 + t2;
	t11 = sin(th);
	t12 = t7 + t2;
	t13 = pow(t10, 0.2e1);
	t14 = pow(t11, 0.2e1);
	t15 = t7 * t13 * t14;
	t17 = a13 + t6;
	t18 = t7 * r;
	t19 = t18 * t10 * t14;
	t20 = t12 * t17;
	t21 = t20 - t19;
	t22 = 2.0 * r;
	t23 = -t22 + t12;
	t17 = pow(t12, 0.2e1) * pow(t17, 0.2e1) - t7 * t4 * t23 * t14;
	t24 = a22 + t7;
	t25 = r * a22;
	t26 = t7 * a13;
	t4 = 2.0 * t4 + t2 * (a13 * t24 + (a22 * t7 + (a13 + t25) * r) * r) + t26 * a22;
	t27 = a52 + t2;
	t28 = t23 * t3 - t15;
	t29 = t8 / 0.3e1 + t2;
	t30 = 0.2e1 / 0.3e1 * t7;
	t31 = 0.3e1 / 0.5e1 * t7;
	t32 = (t31 + t2) * r + 0.2e1 / 0.5e1 * a13;
	t31 = r * t32 - t31 * (a22 / 0.3e1 + t2) * t14;
	t33 = 0.1e1 / t21;
	t35 = pow(t33, 0.2e1);
	t36 = t33 * t35;
	t37 = 0.3e1 * t29;
	t38 = t35 * t14;
	t39 = t7 * a52;
	t40 = a52 + t7;
	t42 = 0.1e1 / r;
	t43 = t9 * t42;
	t44 = -2.0 * t11 * t35 * t1 * (-t43 * t17 + (t5 * t9 * t23 + t17) * t14 * t7) + 4.0 * t1 * t7 * t10 * t11 * t14 * t9 * t17 * t36;
	t21 = 0.1e1 / t21;
	t45 = 0.1e1 / t9;
	t46 = 0.1e1 / t27;
	t47 = 0.1e1 / t23;
	t48 = pow(t21, 0.2e1);
	t19 = 2.0 * t4 * t1 * spin * t11 * (t18 * t14 * (t19 + a22 * eps3 - ((-(a22 - t7) * r + a13 - eps3) * r - t8 * t10) * r - t26) + t20 * t9) * t21 * t48;
	t21 = 2.0 * t7 * t1;
	t5 = r * t9 * (-t12 * t3 + 2.0 * t5 + t15) * t48;
	t6 = (-2.0 * (t3 * (-t40 + r) + t8 * ((t2 * (r - 0.1e1) + a52) * r - t39)) * r + t3 * (a52 * (-6.0) - 0.3e1 * eps3) + (-t2 * t40 + t39) * eps3 + 4.0 * (eps3 + t39) * t6) * pow(t47, 0.2e1) * pow(t46, 0.2e1);
	
	g_tt = t5;
	g_pp = t43 * t17 * t14 * t48;
	g_tp = -t4 * t9 * spin * t14 * t48;
	
	dgttdr = t35 * (t28 * (t9 * (0.10e2 * r * t31 * t33 - 0.1e1) - t37 * r) + (-6.0) * (t2 * ((-0.5e1 / 0.3e1 + r) * r + t30) - t30 * t10 * t14) * t2 * t9);
	dgtpdr = t38 * spin * ((0.20e2 * t9 * t31 * t33 + t29 * (-6.0)) * t4 / 0.2e1 - 0.12e2 * r * ((t24 / 0.6e1 + t2 / 0.3e1) * a13 + t3 + t25 * (0.5e1 / 0.12e2 * t2 + t7 / 0.4e1)) * t9);
	dgrrdr = t6;
	dgththdr = -eps3 * pow(t42, 0.2e1) + t22;
	dgppdr = 0.10e2 * t38 * t9 * (-t31 * t17 * t33 * t42 + t20 * t32 - 0.4e1 / 0.5e1 * t3 * ((-0.7e1 / 0.4e1 + r) * r + 0.3e1 / 0.4e1 * t7) * t7 * t14) + t38 * t17 * t42 * (-t43 + t37);
	

	Omega = (-dgtpdr + sqrt(dgtpdr*dgtpdr - dgttdr*dgppdr)) / dgppdr;
	gi_lp = sqrt(-g_tt_src / (-g_tt - 2.0 * g_tp * Omega - g_pp * Omega * Omega));

	return gi_lp;
}

double gfactor(double r, relParam * param) {

	double epsi3 = 0, a13 = 0, a22 = 0, a52 = 0;
	double gtt_s, gtt, gpp, gtp, dgttdr, dgppdr, dgtpdr;
	double spin = param->a;
	double r0, th, Omega, gfac, denom, num;

	if (param->def_par_type == 1)
	{
		a13 = param->def_par_unscaled;
		// a22 = 0.0;
	}
	if (param->def_par_type == 2)
	{
		a22 = param->def_par_unscaled;
		// a13 = 0.0;
	}

	if (param->emis_type == EMIS_TYPE_LP) {
		r0 = param->height;
		th = 0.0;
	}
	else if (param->emis_type == EMIS_TYPE_RING || param->emis_type == EMIS_TYPE_DISK) {
		r0 = sqrt(param->height * param->height + param->r_ring * param->r_ring);
		th = atan(param->r_ring / param->height);	
	}
	else {
		printf(" The redshift factor cannot be calculated for the given emission type\n");
		exit(1);
	}


	double t2 = pow(r0,3);
	double t3 = pow(spin,2);
	double t11 = pow(r0,2);
	double t14 = a22 + t11;
	double t15 = sin(th);
	double t16 = pow(t15,2);
	gtt_s = (r0*(pow(t14,2)*t16*t3 - pow(r0,4)*((-2 + r0)*r0 + t3))*(epsi3 + t2 + r0*t3*pow(cos(th),2)))/pow(-(r0*t14*t16*t3) + (a13 + t2)*(t11 + t3),2);

	// printf(" new g_tt_src=%f ", gtt_s);
	th = M_PI / 2.0;
	t2 = pow(r,3);
	t3 = pow(spin,2);
	t11 = pow(r,2);
	t14 = a22 + t11;
	t15 = sin(th);
	t16 = pow(t15,2);
	double t5 = cos(th);
	double t6 = pow(t5,2);
	double t7 = r*t3*t6;
	double t8 = epsi3 + t2 + t7;
	double t9 = a13 + t2;
	double t12 = t11 + t3;
	double t13 = t12*t9;
	double t17 = -(r*t14*t16*t3);
	double t18 = t13 + t17;
	double t19 = pow(t18,-2);
	double t21 = -2 + r;
	double t22 = r*t21;
	double t23 = t22 + t3;
	double t29 = 1/r;
	double t20 = pow(r,4);
	double t38 = 2*r;
	double t24 = -(t20*t23);
	double t25 = pow(t14,2);
	double t26 = t16*t25*t3;
	double t27 = t24 + t26;
	double t45 = t3*t6;
	double t30 = pow(t9,2);
	double t32 = pow(t12,2);
	double t34 = pow(r,6);
	double t58 = -2 + t38;
	double t64 = 2*r*t9;
	double t65 = 3*t11*t12;
	double t66 = -2*t11*t16*t3;
	double t67 = -(t14*t16*t3);
	double t68 = t64 + t65 + t66 + t67;
	double t69 = pow(t18,-3);
	double t33 = t30*t32;
	double t35 = -(t16*t23*t3*t34);
	double t36 = t33 + t35;
	double t71 = 3*t11;
	double t72 = t45 + t71;
	double t51 = pow(r,-2);
	double t39 = -t11;
	double t40 = -t3;
	double t41 = pow(r,-5);
	double t42 = t12*t14*t41*t9;
	double t43 = t38 + t39 + t40 + t42;
	double t47 = pow(r,-3);
	double t48 = a13*t47;
	double t49 = 1 + t48;
	double t50 = t12*t49;
	double t52 = a22*t51;
	double t53 = 1 + t52;
	double t54 = -(t16*t3*t53);
	double t55 = t50 + t54;
	double t56 = pow(t55,-2);
	double t91 = pow(r,-4);
	double t44 = epsi3*t29;
	double t46 = t11 + t44 + t45;
	gtt = r*t19*t27*t8;
	gpp = t16*t19*t29*t36*t8;
	gtp = -(spin*t16*t43*t46*t56);
	dgttdr = r*t19*t27*t72 + t19*t27*t8 + r*t19*(-4*t2*t23 + 4*r*t14*t16*t3 - t20*t58)*t8 - 2*r*t27*t68*t69*t8;
	dgppdr = t16*t19*t29*t36*t72 - t16*t19*t36*t51*t8 - 2*t16*t29*t36*t68*t69*t8 + t16*t19*t29*t8*(-6*pow(r,5)*t16*t23*t3 + 4*r*t12*t30 - t16*t3*t34*t58 + 6*t11*t32*t9);
	dgtpdr = -(spin*t16*t43*(t38 - epsi3*t51)*t56) + (2*spin*t16*t43*t46*(2*a22*t16*t3*t47 + 2*r*t49 - 3*a13*t12*t91))/pow(t55,3) - spin*t16*t46*t56*(2 - 2*r + 3*t12*t14*t47 - (5*t12*t14*t9)/pow(r,6) + 2*t12*t9*t91 + 2*t14*t9*t91);

	Omega = (-dgtpdr + sqrt(dgtpdr * dgtpdr - dgttdr * dgppdr)) / dgppdr;
	num = sqrt(-gtt_s);
	denom = sqrt(-gtt - 2.*gtp*Omega - gpp*Omega*Omega); //t-component of 4-velocity
	// printf(" denom=%f \n", gtt + 2.*gtp*Omega + gpp*Omega*Omega);
	if (isnan(Omega)) 
		printf(" !!! Nan value in the redshift calculation 1\n");
	if (isnan(num)) 
		printf(" !!! Nan value in the redshift calculation 2\n");
	if (isnan(denom)) 
		printf(" !!! Nan value in the redshift calculation 3\n");
	// gfac = num / denom;
	// gfac = sqrt(gtt_s / (gtt + 2.*gtp*Omega + gpp*Omega*Omega));
	// printf(" %f\n", sqrt(gtt_s / (gtt + 2.*gtp*Omega + gpp*Omega*Omega)));

	// Omega = (-dgtpdr + sqrt(dgtpdr*dgtpdr - dgttdr*dgppdr)) / dgppdr;
	// gi_lp = sqrt(-g_tt_src / (-g_tt - 2.0 * g_tp * Omega - g_pp * Omega * Omega));
	

	return sqrt(gtt_s / (gtt + 2.*gtp*Omega + gpp*Omega*Omega));
}


double gi_potential(double r, relParam * param) {

	double epsi3 = 0, a13 = 0, a22 = 0, a52 = 0;
	double gtt, gtp, gpp, dgttdr, dgtpdr, dgppdr, r0, th0, gtt_src, gi_lp, Omega, spin, th = M_PI/2;

	if (param->def_par_type == 1)
	{
		a13 = param->def_par_unscaled;
		// a22 = 0.0;
	}
	if (param->def_par_type == 2)
	{
		a22 = param->def_par_unscaled;
		// a13 = 0.0;
	}

	if (param->emis_type == EMIS_TYPE_LP) {
		r0 = param->height;
		th0 = 0.0;
	}
	else if (param->emis_type == EMIS_TYPE_RING || param->emis_type == EMIS_TYPE_DISK) {
		r0 = sqrt(param->height * param->height + param->r_ring * param->r_ring);
		th0 = atan(param->r_ring / param->height);	
	}
	else {
		printf(" The redshift factor cannot be calculated for the given emission type\n");
		exit(1);
	}
	spin = param->a;

	double t2 = pow(r0,3);
	double t3 = pow(spin,2);
	double t11 = pow(r0,2);
	double t14 = a22 + t11;
	double t15 = sin(th0);
	double t16 = pow(t15,2);
	
	gtt_src = (r0*(pow(t14,2)*t16*t3 - pow(r0,4)*((-2 + r0)*r0 + t3))*(epsi3 + t2 + r0*t3*pow(cos(th0),2)))/pow(-(r0*t14*t16*t3) + (a13 + t2)*(t11 + t3),2);

	// printf(" old g_tt_src=%f ", gtt_src);

	// double spin = param->a, th = M_PI/2.0;
	t2 = pow(r,3);
	t3 = pow(spin,2);
	t11 = pow(r,2);
	t14 = a22 + t11;
	t15 = sin(th);
	t16 = pow(t15,2);
	double t5 = cos(th);
	double t6 = pow(t5,2);
	double t7 = r*t3*t6;
	double t8 = epsi3 + t2 + t7;
	double t9 = a13 + t2;
	double t12 = t11 + t3;
	double t13 = t12*t9;
	double t17 = -(r*t14*t16*t3);
	double t18 = t13 + t17;
	double t19 = pow(t18,-2);
	double t21 = -2 + r;
	double t22 = r*t21;
	double t23 = t22 + t3;
	double t29 = 1/r;
	double t20 = pow(r,4);
	double t38 = 2*r;
	double t24 = -(t20*t23);
	double t25 = pow(t14,2);
	double t26 = t16*t25*t3;
	double t27 = t24 + t26;
	double t45 = t3*t6;
	double t30 = pow(t9,2);
	double t32 = pow(t12,2);
	double t34 = pow(r,6);
	double t58 = -2 + t38;
	double t64 = 2*r*t9;
	double t65 = 3*t11*t12;
	double t66 = -2*t11*t16*t3;
	double t67 = -(t14*t16*t3);
	double t68 = t64 + t65 + t66 + t67;
	double t69 = pow(t18,-3);
	double t33 = t30*t32;
	double t35 = -(t16*t23*t3*t34);
	double t36 = t33 + t35;
	double t71 = 3*t11;
	double t72 = t45 + t71;
	double t51 = pow(r,-2);
	double t39 = -t11;
	double t40 = -t3;
	double t41 = pow(r,-5);
	double t42 = t12*t14*t41*t9;
	double t43 = t38 + t39 + t40 + t42;
	double t47 = pow(r,-3);
	double t48 = a13*t47;
	double t49 = 1 + t48;
	double t50 = t12*t49;
	double t52 = a22*t51;
	double t53 = 1 + t52;
	double t54 = -(t16*t3*t53);
	double t55 = t50 + t54;
	double t56 = pow(t55,-2);
	double t91 = pow(r,-4);
	double t44 = epsi3*t29;
	double t46 = t11 + t44 + t45;
	gtt = r*t19*t27*t8;
	gpp = t16*t19*t29*t36*t8;
	gtp = -(spin*t16*t43*t46*t56);
	dgttdr = r*t19*t27*t72 + t19*t27*t8 + r*t19*(-4*t2*t23 + 4*r*t14*t16*t3 - t20*t58)*t8 - 2*r*t27*t68*t69*t8;
	dgppdr = t16*t19*t29*t36*t72 - t16*t19*t36*t51*t8 - 2*t16*t29*t36*t68*t69*t8 + t16*t19*t29*t8*(-6*pow(r,5)*t16*t23*t3 + 4*r*t12*t30 - t16*t3*t34*t58 + 6*t11*t32*t9);
	dgtpdr = -(spin*t16*t43*(t38 - epsi3*t51)*t56) + (2*spin*t16*t43*t46*(2*a22*t16*t3*t47 + 2*r*t49 - 3*a13*t12*t91))/pow(t55,3) - spin*t16*t46*t56*(2 - 2*r + 3*t12*t14*t47 - (5*t12*t14*t9)/pow(r,6) + 2*t12*t9*t91 + 2*t14*t9*t91);

	Omega = (-dgtpdr + sqrt(dgtpdr*dgtpdr - dgttdr*dgppdr)) / dgppdr;
	gi_lp = sqrt(-gtt_src / (-gtt - 2.0 * gtp * Omega - gpp * Omega * Omega));
	// printf(" denom=%f \n", gtt + 2.0 * gtp * Omega + gpp * Omega * Omega);
	return gi_lp;
}




/** print the xillver spectrum   **/
void save_xillver_spectrum(double* ener, double* flu, int n_ener, char* fname){

	FILE* fp =  fopen ( fname,"w+" );
	int ii;
	for (ii=0; ii<n_ener; ii++){
		fprintf(fp, " %e \t %e \t %e \n",ener[ii],ener[ii+1],flu[ii]);
	}
	if (fclose(fp)) exit(1);
}


/* A simple implementation of the FFT taken from http://paulbourke.net/miscellaneous/dft/
   (uses the Radix-2 Cooley-Tukey algorithm)

   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform
*/
void FFT_R2CT(short int dir,long m,double *x,double *y){

   long n,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   n = 1;
   for (i=0;i<m;i++){
      n *= 2;
   }

   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0;
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0;
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1;
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1)
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<n;i++) {
         x[i] /= n;
         y[i] /= n;
      }
   }

   return;
}



// /** get the number of zones on which we calculate the relline-spectrum **/
// int get_num_zones(int model_type, int emis_type, int ion_grad_type){

// 	char* env;
// 	env = getenv("RELXILL_NUM_RZONES");
// 	int env_n_zones = 0;
// 	if (env != NULL){
// 		env_n_zones = atof(env);
// 	}


// 	// set the number of zones in radial direction (1 for relline/conv model, N_ZONES for xill models)
// 	if (is_iongrad_model(model_type, ion_grad_type)){
// 		if (env != NULL){
// 			if ( (env_n_zones > 9 ) && (env_n_zones <= N_ZONES_MAX) ){
// 				return env_n_zones;
// 			} else {
// 				printf(" *** warning: value of %i for RELXILL_NUM_RZONES not within required interval of [10,%i] \n",env_n_zones,N_ZONES_MAX);
// 			}
// 		}
// 		return N_ZONES_MAX;
// 	} else if (is_relxill_model(model_type) && (emis_type==EMIS_TYPE_LP)){

// 		if (env != NULL){
// 			if ( (env_n_zones > 0 ) && (env_n_zones <= N_ZONES_MAX) ){
// 				return env_n_zones;
// 			} else {
// 				printf(" *** warning: value of %i for RELXILL_NUM_RZONES not within required interval of [1,%i] \n",env_n_zones,N_ZONES_MAX);
// 			}
// 		}
// 		// return 1;
// 		return N_ZONES;
// 	} else if (is_relxill_model(model_type) && (emis_type==EMIS_TYPE_BKN) ){
// 	// } else if (is_relxillion_model(model_type, ion_grad_type) && emis_type==EMIS_TYPE_BKN) { // relxill_nk variable ionization
// 		if (env != NULL){
// 			if ( (env_n_zones > 9 ) && (env_n_zones <= N_ZONES_MAX) ){
// 				return env_n_zones;
// 			} else {
// 				printf(" *** warning: value of %i for RELXILL_NUM_RZONES not within required interval of [10,%i] \n",env_n_zones,N_ZONES_MAX);
// 			}
// 		}
// 		return N_ZONES_MAX;
// 	} else if (is_relxill_model(model_type) && (emis_type==EMIS_TYPE_RING)) {
// 		if (env != NULL){
// 			if ( (env_n_zones > 0 ) && (env_n_zones <= N_ZONES_MAX) ){
// 				return env_n_zones;
// 			} else {
// 				printf(" *** warning: value of %i for RELXILL_NUM_RZONES not within required interval of [1,%i] \n",env_n_zones,N_ZONES_MAX);
// 			}
// 		}
// 		// return 1;
// 		return N_ZONES;
// 	} else if (is_relxill_model(model_type) && (emis_type==EMIS_TYPE_DISK)) {
// 		if (env != NULL){
// 			if ( (env_n_zones > 0 ) && (env_n_zones <= N_ZONES_MAX) ){
// 				return env_n_zones;
// 			} else {
// 				printf(" *** warning: value of %i for RELXILL_NUM_RZONES not within required interval of [1,%i] \n",env_n_zones,N_ZONES_MAX);
// 			}
// 		}
// 		// return 1;
// 		return N_ZONES;
// 	} else if (model_type == MOD_TYPE_RELXILLOLD || model_type > 0){ 
// 		return 1;
// 	} else {
// 		// return 1;
// 		return N_ZONES_MAX;
// 	}

// }


/** get the number of zones on which we calculate the relline-spectrum **/
int get_num_zones(int model_type, int emis_type, int ion_grad_type){

	char* env;
	env = getenv("RELXILL_NUM_RZONES");
	int env_n_zones = 0;
	if (env != NULL){
		env_n_zones = atof(env);
	}


	// set the number of zones in radial direction (1 for relline/conv model, N_ZONES for xill models)
	if (is_iongrad_model(model_type, ion_grad_type) || fifty_zones(model_type)){
		if (env != NULL){
			if ( (env_n_zones > 9 ) && (env_n_zones <= N_ZONES_MAX) ){
				return env_n_zones;
			} else {
				printf(" *** warning: value of %i for RELXILL_NUM_RZONES not within required interval of [10,%i] \n",env_n_zones,N_ZONES_MAX);
			}
		}
		return N_ZONES_MAX;
	} else if (is_relxill_model(model_type) && (emis_type==EMIS_TYPE_LP || emis_type==EMIS_TYPE_RING || emis_type==EMIS_TYPE_DISK)){

		if (env != NULL){
			if ( (env_n_zones > 0 ) && (env_n_zones <= N_ZONES_MAX) ){
				return env_n_zones;
			} else {
				printf(" *** warning: value of %i for RELXILL_NUM_RZONES not within required interval of [1,%i] \n",env_n_zones,N_ZONES_MAX);
			}
		}
		return N_ZONES;
	} else {
		return 1;
	}

}


/** rebin spectrum to a given energy grid
 *  length of ener is nbins+1       **/

void rebin_spectrum(double* ener, double* flu, int nbins, double* ener0, double* flu0, int nbins0){

	int ii; int jj;
	int imin = 0;
	int imax = 0;

	for (ii=0; ii<nbins; ii++){

		flu[ii] = 0.0;

		/* check of the bin is outside the given energy range */
		if ( (ener0[0] <= ener[ii+1]) && (ener0[nbins0] >= ener[ii]) ){

			/* need to make sure we are in the correct bin */
			while ( ener0[imin]<=ener[ii] && imin<=nbins0){
				imin++;
			}
			// need to set it back, as we just crossed to the next bin
			if (imin>0){
				imin--;
			}
			while ( (ener0[imax]<=ener[ii+1] && imax<nbins0)){
				imax++;
			}
			if (imax>0){
				imax--;
			}

			double elo = ener[ii];
			double ehi = ener[ii+1];
			if (elo < ener0[imin]) elo=ener0[imin];
			if (ehi > ener0[imax+1]) ehi=ener0[imax+1];

			if (imax==imin){
				flu[ii] = (ehi-elo) / (ener0[imin+1] - ener0[imin]) * flu0[imin];
			} else {

				double dmin=(ener0[imin+1]-elo)/(ener0[imin+1]-ener0[imin]);
				double dmax=(ehi-ener0[imax])/(ener0[imax+1]-ener0[imax]);

				flu[ii] += flu0[imin]*dmin + flu0[imax]*dmax;

				for (jj=imin+1; jj <= imax-1; jj++) {
					flu[ii] += flu0[jj];
				}

			}

		}

	}
}

int do_renorm_model(relParam* rel_param){


	int renorm = 0;

	if ( is_relxill_model(rel_param->model_type) ){
			if (rel_param->emis_type == EMIS_TYPE_LP || rel_param->emis_type == EMIS_TYPE_RING || rel_param->emis_type == EMIS_TYPE_DISK){
				renorm=0;
			} else {
				renorm=1;
			}
	} else {
		if ( do_not_normalize_relline() ){
			renorm=0;
		} else {
			renorm=1;
		}
	}

	return renorm;
}

void get_nthcomp_param( double* nthcomp_param, double gam, double kte, double z){
	  nthcomp_param[0] = gam;
	  nthcomp_param[1] = kte;
	  nthcomp_param[2] = 0.05; // ktbb
	  nthcomp_param[3] = 1.0;  // inp_type
	  nthcomp_param[4] = z;

}




/*** we calculate the disk density from  Shakura & Sunyaev (1973)
 *    - for zone A as describe in their publication,  formula (Eq 2.11)
 *    - only the radial dependence is picked up here  (viscosity alpha=const.)
 *    - normalized such that dens(rms) = 1
 *                                                               ***/
double density_ss73_zone_a(double radius, double rms){
	return pow((radius/rms),(3./2)) * pow((1-sqrt(rms/radius)),-2);
}


// calculate the log(xi) for given density and emissivity
static double cal_lxi(double dens, double emis){
	return log10(4.0*M_PI*emis/dens);
}

// determine the radius of maximal ionization
static double cal_lxi_max_ss73(double* re, double* emis, int nr, double rin, int* status){

	CHECK_STATUS_RET(*status,0.0);

	double rad_max_lxi  = pow((11./9.),2) * rin;  // we use the same definition as Adam with r_peak = (11/9)^2 rin to be consistent (does not matter much)

	// radial AD grid is sorted descending (!)
	int kk = inv_binary_search(re, nr, rad_max_lxi);
	double interp = (rad_max_lxi - re[kk+1]) / ( re[kk] - re[kk+1] );

	double emis_max_lxi = interp_lin_1d(interp,emis[kk+1],emis[kk]);


	double lxi_max =  cal_lxi( density_ss73_zone_a(rad_max_lxi,rin) , emis_max_lxi);

	return lxi_max;
}

static void save_ion_profile(ion_grad* ion){

	FILE* fp =  fopen ( "test_ion_grad_relxill.dat","w+" );
	int ii;
	for (ii=0; ii<ion->nbins; ii++){
		fprintf(fp, " %e \t %e \t %e \n",ion->r[ii],ion->r[ii+1],ion->lxi[ii]);
	}
	if (fclose(fp)) exit(1);

}

/** *** set log(xi) to obey the limits of the xillver table: TODO: check if we need to adjust the normalization as well  ***
 *  NOTE: with correctly set xpsec/isis limits, it is only possible to reach the lower boundary       **/
static void lxi_set_to_xillver_bounds(double* pt_lxi){
	/**  TODO: Need to define this globally **/
	double xlxi_tab_min = 0.0;
	double xlxi_tab_max = 4.7;

	// #1: set the value of xi to the lowest value of the table
    if (*pt_lxi < xlxi_tab_min){
    	*pt_lxi = xlxi_tab_min;
    } else if (*pt_lxi > xlxi_tab_max){
    //	#2: high ionization: we approximately assume such a highly ionized disk acts as a mirror
    	*pt_lxi = xlxi_tab_max;
    }

}

ion_grad* calc_ion_gradient(relParam* rel_param, double xlxi0, double xindex, int type, double* rgrid, int n, int* status) {

	CHECK_STATUS_RET(*status,NULL);

	ion_grad* ion = new_ion_grad(rgrid, n, status);
	CHECK_STATUS_RET(*status,NULL);

	double rmean[n];
	double del_inc[n];
	int ii;
	for (ii=0; ii<n; ii++){
		rmean[ii] = 0.5*(rgrid[ii] + rgrid[ii+1]);
	}

	if (type==ION_GRAD_TYPE_PL){
		for (ii=0; ii<n; ii++){
		     ion->lxi[ii] = (exp(xlxi0))* pow((rmean[ii]/rmean[0]),-1.0*xindex);  // TODO: check if we need to subtract xlxi_tab_min here
		     ion->lxi[ii] = log(ion->lxi[ii]);


		     lxi_set_to_xillver_bounds(&(ion->lxi[ii]));

		}

	} else if (type==ION_GRAD_TYPE_ALPHA){
		double dens[n];
		double rin = rgrid[0];

		// TODO: use a better approach to not linearly interpolate but rather average over the profile?
		double emis_zones[n];

		// we need the emissivity profile (should be cached, so no extra effort required here)
		relSysPar* sysPar = get_system_parameters(rel_param,status);
		assert(sysPar->emis!=NULL);
		inv_rebin_mean(sysPar->re, sysPar->emis, sysPar->nr, rmean, emis_zones, n,  status);
		inv_rebin_mean(sysPar->re, sysPar->del_inc, sysPar->nr, rmean, del_inc, n,  status);
		inv_rebin_mean(sysPar->re, sysPar->del_emit, sysPar->nr, rmean, ion->del_emit, n,  status);

		// calculate the maximal ionization assuming r^-3 and SS73 alpha disk
		double lxi_max = cal_lxi_max_ss73(sysPar->re, sysPar->emis, sysPar->nr, rin, status);

		// the maximal ionization is given as input parameter, so we need to normalize our calculation by this value
		double fac_lxi_norm = xlxi0 - lxi_max; // subtraction instead of division because of the log

		/** calculate the density for a  stress-free inner boundary condition, i.e., R0=rin in SS73)  **/
		for (ii=0; ii<n; ii++){
			dens[ii] = density_ss73_zone_a(rmean[ii],rin);

			// now we can use the emissivity to calculate the ionization
			ion->lxi[ii] = cal_lxi(dens[ii], emis_zones[ii]) + fac_lxi_norm;

			ion->lxi[ii] += log10(cos(M_PI/4)/cos(del_inc[ii]));

		    lxi_set_to_xillver_bounds(&(ion->lxi[ii]));
		}


	} else if (type==ION_GRAD_TYPE_CONST) {  // should not happen, as this will be approximated by 1 zone (but just in case we get here...)
		for (ii=0; ii<n; ii++){
		     ion->lxi[ii] = xlxi0;
		}
	} else {
		printf(" *** ionization type with number %i not implemented \n",type);
		printf("     choose either %i for the PL, %i for the ALPHA-disk, or %i for constant\n",
				ION_GRAD_TYPE_PL,ION_GRAD_TYPE_ALPHA,ION_GRAD_TYPE_CONST);
		RELXILL_ERROR("unknown ionization gradient type", status);
	}



	if (is_debug_run()){
		save_ion_profile(ion);
	}

	if (*status!=EXIT_SUCCESS){
		RELXILL_ERROR("calculating the ionization gradient failed due to previous error", status);
	}

	return ion;
}


ion_grad* new_ion_grad(double* r, int n, int* status){

	ion_grad* ion = (ion_grad*) malloc(  sizeof(ion_grad));
	CHECK_MALLOC_RET_STATUS(ion,status,NULL);

	ion->r = (double*) malloc ( (n+1)*sizeof(double));
	CHECK_MALLOC_RET_STATUS(ion->r,status,NULL);
	ion->lxi = (double*) malloc ( (n)*sizeof(double));
	CHECK_MALLOC_RET_STATUS(ion->lxi,status,NULL);
	ion->fx = (double*) malloc ( (n)*sizeof(double));
	CHECK_MALLOC_RET_STATUS(ion->fx,status,NULL);
	ion->del_emit = (double*) malloc ( (n)*sizeof(double));
	CHECK_MALLOC_RET_STATUS(ion->del_emit,status,NULL);

	ion->nbins = n;

	int ii;
	for (ii=0; ii<n; ii++){
		ion->r[ii] = r[ii];
		ion->lxi[ii] = 0.0;
		ion->fx[ii] = 0.0;
		ion->del_emit[ii] = M_PI/4.; // assume default 45 deg (xillver assumption), only used if beta>0
	}
	// radius goes to n+1
	ion->r[n] = r[n];

	return ion;
}

void free_ion_grad(ion_grad* ion){

	if (ion!=NULL){
		if (ion->r!=NULL){
			free(ion->r);
		}
		if (ion->lxi!=NULL){
			free(ion->lxi);
		}
		if (ion->fx!=NULL){
			free(ion->fx);
		}
		free(ion);
	}
}

// for x0 descending and xn ascending, calculate the mean at xn
void inv_rebin_mean(double* x0, double* y0, int n0, double*  xn, double* yn, int nn, int* status){
	if (xn[0]>xn[nn-1] || x0[nn-1]>x0[0]){
		RELXILL_ERROR(" *** grid in wrong order",status);
		return;
	}
	if (xn[0]<x0[n0-1] || xn[nn-1]>x0[0]){
		RELXILL_ERROR(" *** new grid is larger than the input grid",status);
		return;
	}

	int in = nn-1; // traverse new array backwards
	int ii;
	for (ii=0; ii<n0-2; ii++){  // only go to the second to last bin

		if (x0[ii]>xn[in] && x0[ii+1]<=xn[in]){

			double ifac_r = (xn[in]-x0[ii+1]) /  (x0[ii] - x0[ii+1] );
			yn[in] = interp_lin_1d(ifac_r, y0[ii+1], y0[ii]);

			in--;
			if (in<0){ // we can stop if once we are below zero as the array is filled
				break;
			}
		}

	}

}
