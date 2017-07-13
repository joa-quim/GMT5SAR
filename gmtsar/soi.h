/************************************************************************
* soi.h is the include file for the esarp SAR processor.		*
************************************************************************/
/************************************************************************
* Creator: Evelyn J. Price	(Scripps Institution of Oceanography)	*
* Date   : 11/18/96							*
************************************************************************/
/************************************************************************
* Modification History							*
*									*
* Date									*
*									*
*  4/23/97- 	added parameters for orbit calculations: x_target,      *
*		y_target,z_target,baseline,alpha,sc_identity,		*
*		ref_identity,SC_clock_start,SC_clock_stop,              *
*		clock_start,clock_stop   				*
*		-DTS							*
*									*
* 4/23/97-	added parameters: rec_start, rec_stop			*
*		-EJP							*
*									*
* 8/28/97-	added parameters baseline_start baseline_end		*
*		alpha_start alpha_end					*
*									*
* 9/12/97	added clipi2 function to clip to short int		*
*									*
* 4/26/06	added nrows, num_lines					*
************************************************************************/
#ifndef SOI_H	    
#define SOI_H	    
#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <time.h>
#include <stdlib.h>
#define SOL 299792456.0
#define PI 3.1415926535897932
#define PI2 6.2831853071795864
#ifndef I2MAX
#	define I2MAX 32767
#endif
#define I2SCALE 4.e6
#define TRUE 1
#define FALSE 0
#define RW 0666
#define MULT_FACT 1000.0 
#define sgn(A) ((A) >= 0.0 ? 1.0 : -1.0)
#define clipi2(A) ( ((A) > I2MAX) ? I2MAX : (((A) < -I2MAX) ? -I2MAX : A) )
#include "sfd_complex.h"
#include "declspec.h"

EXTERN_MSC char *input_file;
EXTERN_MSC char *led_file;
EXTERN_MSC char *out_amp_file;
EXTERN_MSC char *out_data_file;
EXTERN_MSC char *deskew;
EXTERN_MSC char *iqflip;
EXTERN_MSC char *off_vid;
EXTERN_MSC char *srm;
EXTERN_MSC char *ref_file;
EXTERN_MSC char *orbdir;
EXTERN_MSC char *lookdir;

EXTERN_MSC int debug_flag;
EXTERN_MSC int bytes_per_line;
EXTERN_MSC int good_bytes;
EXTERN_MSC int first_line;
EXTERN_MSC int num_patches;
EXTERN_MSC int first_sample;
EXTERN_MSC int num_valid_az;
EXTERN_MSC int st_rng_bin;
EXTERN_MSC int num_rng_bins;
EXTERN_MSC int nextend;
EXTERN_MSC int nlooks;
EXTERN_MSC int xshift;
EXTERN_MSC int yshift;
EXTERN_MSC int fdc_ystrt;
EXTERN_MSC int fdc_strt;

/*New parameters 4/23/97 -EJP */
EXTERN_MSC int rec_start;
EXTERN_MSC int rec_stop;
/* End new parameters 4/23/97 -EJP */ 

/* New parameters 4/23/97 -DTS */
EXTERN_MSC int SC_identity;	/* (1)-ERS1 (2)-ERS2 (3)-Radarsat (4)-Envisat (5)-ALOS  (6)-Envisat_SLC  (7)-TSX (8)-CSK (9)-RS2 (10)-S1A*/
EXTERN_MSC int ref_identity;	/* (1)-ERS1 (2)-ERS2 (3)-Radarsat (4)-Envisat (5)-ALOS  (6)-Envisat_SLC  (7)-TSX (8)-CSK (9)-RS2 (10)-S1A*/
EXTERN_MSC double SC_clock_start;	/* YYDDD.DDDD */
EXTERN_MSC double SC_clock_stop;	/* YYDDD.DDDD */
EXTERN_MSC double icu_start;       /* onboard clock counter */
EXTERN_MSC double clock_start;     /* DDD.DDDDDDDD  clock without year has more precision */
EXTERN_MSC double clock_stop;      /* DDD.DDDDDDDD  clock without year has more precision */
/* End new parameters 4/23/97 -DTS */

EXTERN_MSC double caltone;
EXTERN_MSC double RE;		/* Local Earth radius */
EXTERN_MSC double raa;             /* ellipsoid semi-major axis - added by RJM */
EXTERN_MSC double rcc;             /* ellipsoid semi-minor axis - added by RJM */
EXTERN_MSC double vel1;		/* Equivalent SC velocity */
EXTERN_MSC double ht1;		/* (SC_radius - RE) center of frame*/
EXTERN_MSC double ht0;		/* (SC_radius - RE) start of frame */
EXTERN_MSC double htf;		/* (SC_radius - RE) end of frame */
EXTERN_MSC double near_range;
EXTERN_MSC double far_range;
EXTERN_MSC double prf1;
EXTERN_MSC double xmi1;
EXTERN_MSC double xmq1;
EXTERN_MSC double az_res;
EXTERN_MSC double fs;
EXTERN_MSC double slope;
EXTERN_MSC double pulsedur;
EXTERN_MSC double lambda;
EXTERN_MSC double rhww;
EXTERN_MSC double pctbw;
EXTERN_MSC double pctbwaz;
EXTERN_MSC double fd1;
EXTERN_MSC double fdd1;
EXTERN_MSC double fddd1;
EXTERN_MSC double sub_int_r;
EXTERN_MSC double sub_int_a;
EXTERN_MSC double stretch_r;
EXTERN_MSC double stretch_a;
EXTERN_MSC double a_stretch_r;
EXTERN_MSC double a_stretch_a;

/* New parameters 8/28/97 -DTS */
EXTERN_MSC double baseline_start;
EXTERN_MSC double baseline_center;
EXTERN_MSC double baseline_end;
EXTERN_MSC double alpha_start;
EXTERN_MSC double alpha_center;
EXTERN_MSC double alpha_end;
/* End new parameters 8/28/97 -DTS */
EXTERN_MSC double bparaa;               /* parallel baseline - added by RJM */
EXTERN_MSC double bperpp;               /* perpendicular baseline - added by RJM */

/* New parameters 4/26/06 */
EXTERN_MSC int nrows;
EXTERN_MSC int num_lines;

/* New parameters 09/18/08 */
EXTERN_MSC double TEC_start;
EXTERN_MSC double TEC_end;
#endif /* SOI_H	*/
