/* taken from soi.h */
#ifndef IMAGE_SIO_H
#define IMAGE_SIO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#	include "../../../unistd.h"
#else
#	include <unistd.h>
#endif
#include <math.h>
#include <time.h>
#include "../../../gmtsar/PRM.h"
#include "../../../gmtsar/sfd_complex.h"

#define SOL 299792456.0
#define PI 3.1415926535897932
#define PI2 6.2831853071795864
#define I2MAX 32767.0
#define I2SCALE 4.e6
#define TRUE 1
#define FALSE 0
#define RW 0666
#define MULT_FACT 1000.0
#define sgn(A) ((A) >= 0.0 ? 1.0 : -1.0)
#define clipi2(A) ( ((A) > I2MAX) ? I2MAX : (((A) < -I2MAX) ? -I2MAX : A) )
#define nint(x) (int)rint(x)
#define ERS1 1
#define ERS2 2
#define RSAT 3
#define ENVS 4
#define ALOS 5

#define EXIT_FLAG 1
#define paka(p) {perror((p)); exit(EXIT_FLAG);}
#define MALLOC(p,s) if (((p) = malloc(s)) == NULL) {paka("error: malloc()  ");}

#define NULL_DATA 15
#define NULL_INT -99999
#define NULL_DOUBLE -99999.9999
#define NULL_CHAR  "XXXXXXXX"


/*
offset_video 		off_vid		
chirp_ext 		nextend
-------------------------------
scnd_rng_mig 		srm
Flip_iq 		iqflip
reference_file		ref_file
rng_spec_wgt 		rhww
rm_rng_band 		pctbw
rm_az_band 		pctbwaz
rng_samp_rate		fs
good_bytes_per_line 	good_bytes
earth_radius		RE
SC_vel			vel
SC_height		ht
SC_height_start		ht_start
SC_height_end		ht_end
PRF			prf
I_mean			xmi
Q_mean			xmq
pulse_dur		pulsedur
radar_wavelength	lambda
rng_spec_wgt		rhww

*/
int	verbose;	/* controls minimal level of output 	*/ 
int	debug; 		/* more output 				*/
int	roi; 		/* more output 				*/
int	swap; 		/* whether to swap bytes 		*/
int	quad_pol; 	/* quad polarization data 		*/
int	ALOS_format; 	/* AUIG:  ALOS_format = 0  		*/
			/* ERSDAC:  ALOS_format = 1  		*/
			/* ALOS2:   ALOS_format =2              */
int	force_slope;	/* whether to set the slope	 	*/
int	dopp;		/* whether to calculate doppler 	*/
int	quiet_flag;	/* reduce output			*/
int	SAR_mode;	/* 0 => high-res 			*/	
			/* 1 => wide obs 			*/
			/* 2 => polarimetry 			*/
			/* from ALOS Product Format 3-2		*/
int	prefix_off;	/* offset needed for ALOS-2 prefix size */
double	forced_slope;	/* value to set chirp_slope to          */
double  tbias;          /* time bias for clock bias             */
double  rbias;          /* range bias for near range corr       */
double  slc_fact;       /* factor to convert float to int slc   */
#endif	/* IMAGE_SIO_H */
