#ifndef LIB_FUNCTIONS_H
#define LIB_FUNCTIONS_H

/* Shut up MSVC " possible loss of data" warnings */ 
#		pragma warning( disable : 4244 )
#		pragma warning( disable : 4273 )	/* the bloody "inconsistent dll linkage" -- DRASTIC */
#		pragma warning( disable : 4305 )

#include <math.h>
#include <stdio.h>
/* include files to define sarleader structure */
#include "data_ALOS.h"
#include "data_ALOSE.h"

#include "orbit_ALOS.h"
#include "sarleader_ALOS.h"
#include "sarleader_fdr.h"
#include "image_sio.h"
#include "../../../gmtsar/declspec.h"
#ifndef M_PI
#	define M_PI   3.14159265358979323846
#endif

/* function prototypes 				*/
EXTERN_MSC void calc_height_velocity(struct ALOS_ORB *, struct PRM *, double, double, double *, double *, double *, double *, double *);
EXTERN_MSC void calc_dop(struct PRM *);
void cfft1d_(int *, fcomplex *, int *);
EXTERN_MSC void read_data(fcomplex *, unsigned char *, int, struct PRM *);
EXTERN_MSC void null_sio_struct(struct PRM *);
EXTERN_MSC void get_sio_struct(FILE *, struct PRM *);
EXTERN_MSC void put_sio_struct(struct PRM, FILE *);
EXTERN_MSC void get_string(char *, char *, char *, char *);
EXTERN_MSC void get_int(char *, char *, char *, int *);
EXTERN_MSC void get_double(char *, char *, char *, double *);
EXTERN_MSC void hermite_c(double *, double *, double *, int, int, double, double *, int *);
EXTERN_MSC void interpolate_ALOS_orbit_slow(struct ALOS_ORB *, double, double *, double *, double *, int *);
EXTERN_MSC void interpolate_ALOS_orbit(struct ALOS_ORB *, double *, double *, double *, double, double *, double *, double *, int *);
EXTERN_MSC void get_orbit_info(struct ALOS_ORB *, struct SAR_info);
EXTERN_MSC void get_attitude_info(struct ALOS_ATT *, int, struct SAR_info);
EXTERN_MSC void print_binary_position(struct sarleader_binary *, int, FILE *, FILE *);
EXTERN_MSC int  is_big_endian_(void);
EXTERN_MSC int  is_big_endian__(void);
EXTERN_MSC void die (char *, char *);
EXTERN_MSC void cross3_(double *, double *, double *);
EXTERN_MSC void get_seconds(struct PRM, double *, double *);
EXTERN_MSC void plh2xyz(double *, double *, double, double);
EXTERN_MSC void xyz2plh(double *, double *, double, double);
EXTERN_MSC void polyfit(double *, double *, double *, int *, int *);
void gauss_jordan(double **, double *, double *, int *);
EXTERN_MSC int  find_fft_length(int);
EXTERN_MSC void rng_expand(fcomplex *, int,  fcomplex *, int );
EXTERN_MSC void rng_compress(fcomplex *cin, int nffti, fcomplex *cout, int nffto);
EXTERN_MSC void rng_filter(fcomplex *cin, int nffti, fcomplex *cout);

EXTERN_MSC void ALOS_ldr_orbit(struct ALOS_ORB *, struct PRM *);
EXTERN_MSC void read_ALOS_sarleader(FILE *, struct PRM *, struct ALOS_ORB *);
EXTERN_MSC void ALOS_ldr_prm(struct SAR_info, struct PRM *);
EXTERN_MSC void set_ALOS_defaults(struct PRM *prm);
EXTERN_MSC void print_ALOS_defaults(struct PRM *prm);
EXTERN_MSC long read_ALOS_data_SLC (FILE *imagefile, FILE *outfile, struct PRM *prm, long *byte_offset);
void swap_ALOS_data_info(struct sardata_info *sdr);
EXTERN_MSC void swap32(char *in, char *out, int n);

/* ffpack functions */
#ifndef MAXFAC
#	define MAXFAC 13    /* maximum number of factors in factorization of n */
#endif
#ifdef DOUBLE
#define Treal double
#else
#define Treal float
#endif

EXTERN_MSC void cfftf(int n, Treal c[], Treal wsave[]);
EXTERN_MSC void cfftb(int n, Treal c[], Treal wsave[]);
EXTERN_MSC void cffti(int n, Treal wsave[]);
EXTERN_MSC void rfftb1(int n, Treal c[], Treal ch[], const Treal wa[], const int ifac[MAXFAC+2]);
EXTERN_MSC void rfftf(int n, Treal r[], Treal wsave[]);
EXTERN_MSC void rfftb(int n, Treal r[], Treal wsave[]);
EXTERN_MSC void rffti(int n, Treal wsave[]);


#endif  /* LIB_FUNCTIONS_H */
