/************************************************************************

************************************************************************/
#ifndef SOI_GLOBALS_H	    
#define SOI_GLOBALS_H	    
#include <stdlib.h>

char *input_file = NULL;
char *led_file = NULL;
char *out_amp_file = NULL;
char *out_data_file = NULL;
char *deskew = NULL;
char *iqflip = NULL;
char *off_vid = NULL;
char *srm = NULL;
char *ref_file = NULL;
char *orbdir = NULL;
char *lookdir = NULL;

int debug_flag = 0;
int bytes_per_line = 0;
int good_bytes = 0;
int first_line = 0;
int num_patches = 0;
int first_sample = 0;
int num_valid_az = 0;
int st_rng_bin = 0;
int num_rng_bins = 0;
int nextend = 0;
int nlooks = 0;
int xshift = 0;
int yshift = 0;
int fdc_ystrt = 0;
int fdc_strt = 0;

/*New parameters 4/23/97 -EJP */
int rec_start = 0;
int rec_stop = 0;
/* End new parameters 4/23/97 -EJP */ 

/* New parameters 4/23/97 -DTS */
int SC_identity = 0;	/* (1)-ERS1 (2)-ERS2 (3)-Radarsat (4)-Envisat (5)-ALOS  (6)-  (7)-TSX (8)-CSK (9)-RS2 (10)-S1A*/
int ref_identity = 0;	/* (1)-ERS1 (2)-ERS2 (3)-Radarsat (4)-Envisat (5)-ALOS  (6)-  (7)-TSX (8)-CSK (9)-RS2 (10)-S1A*/
double SC_clock_start = 0;	/* YYDDD.DDDD */
double SC_clock_stop = 0;	/* YYDDD.DDDD */
double icu_start = 0;       /* onboard clock counter */
double clock_start = 0;     /* DDD.DDDDDDDD  clock without year has more precision */
double clock_stop = 0;      /* DDD.DDDDDDDD  clock without year has more precision */
/* End new parameters 4/23/97 -DTS */

double caltone = 0;
double RE = 0;		/* Local Earth radius */
double raa = 0;             /* ellipsoid semi-major axis - added by RJM */
double rcc = 0;             /* ellipsoid semi-minor axis - added by RJM */
double vel1 = 0;		/* Equivalent SC velocity */
double ht1 = 0;		/* (SC_radius - RE) center of frame*/
double ht0 = 0;		/* (SC_radius - RE) start of frame */
double htf = 0;		/* (SC_radius - RE) end of frame */
double near_range = 0;
double far_range = 0;
double prf1 = 0;
double xmi1 = 0;
double xmq1 = 0;
double az_res = 0;
double fs = 0;
double slope = 0;
double pulsedur = 0;
double lambda = 0;
double rhww = 0;
double pctbw = 0;
double pctbwaz = 0;
double fd1 = 0;
double fdd1 = 0;
double fddd1 = 0;
double sub_int_r = 0;
double sub_int_a = 0;
double stretch_r = 0;
double stretch_a = 0;
double a_stretch_r = 0;
double a_stretch_a = 0;

/* New parameters 8/28/97 -DTS */
double baseline_start = 0;
double baseline_center = 0;
double baseline_end = 0;
double alpha_start = 0;
double alpha_center = 0;
double alpha_end = 0;
/* End new parameters 8/28/97 -DTS */
double bparaa = 0;               /* parallel baseline - added by RJM */
double bperpp = 0;               /* perpendicular baseline - added by RJM */

/* New parameters 4/26/06 */
int nrows = 0;
int num_lines = 0;

/* New parameters 09/18/08 */
double TEC_start = 0;
double TEC_end = 0;
#endif /* SOI_GLOBALS_H	*/
