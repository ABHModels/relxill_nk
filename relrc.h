

#ifndef RELRC_H_
#define RELRC_H_

#include "relbase.h"
#include "relutility.h"

void get_emis_ring_corona(relParam* param, double* emis, double* re, int n_r, int* status);


void get_emis_disk_corona(relParam* param, double* emis,
		double* re, int n_r, int* status);


#endif /* RELRC_H_ */
