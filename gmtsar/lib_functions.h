/*	$Id: lib_functions.h 39 2013-04-07 00:49:34Z pwessel $	*/
/* include files to define sarleader structure */
#ifndef LIB_FUNCTIONS_H
#define LIB_FUNCTIONS_H

/* Shut up MSVC " possible loss of data" warnings */ 
#		pragma warning( disable : 4244 )
#		pragma warning( disable : 4273 )	/* the bloody "inconsistent dll linkage" -- DRASTIC */
#		pragma warning( disable : 4305 )

#include "sarleader_ALOS.h"
#include "sarleader_fdr.h"
#include "PRM.h"
#include "xcorr.h"
#include "sfd_complex.h"
#include "declspec.h"
#include "orbit_ALOS.h"
#include "gmt.h"

#ifndef M_PI
#define M_PI          3.14159265358979323846
#endif

/* function prototypes 				*/
EXTERN_MSC void null_sio_struct(struct PRM *);
EXTERN_MSC void get_sio_struct(FILE *, struct PRM *);
EXTERN_MSC void put_sio_struct(struct PRM, FILE *);
EXTERN_MSC void get_string(char *, char *, char *, char *);
EXTERN_MSC void get_int(char *, char *, char *, int *);
EXTERN_MSC void get_double(char *, char *, char *, double *);
EXTERN_MSC void ALOS_ldr_prm(struct SAR_info, struct PRM *);
EXTERN_MSC int is_big_endian_(void);
EXTERN_MSC int is_big_endian__(void);
EXTERN_MSC void die (char *, char *);
EXTERN_MSC int get_prm(struct PRM *p, char *filename);
EXTERN_MSC void cross3(double *, double *, double *);
EXTERN_MSC void get_seconds(struct PRM, double *, double *);
EXTERN_MSC void plh2xyz(double *, double *, double, double);
EXTERN_MSC void xyz2plh(double *, double *, double, double);
EXTERN_MSC void find_unit_vector(double *, double *);
EXTERN_MSC double	find_length(double *);
EXTERN_MSC double find_distance(double *, double *);
EXTERN_MSC double find_distance3(double, double, double, double, double, double);
EXTERN_MSC int geo2latlon(double *, double *, struct PRM);
EXTERN_MSC void geoxyz(double, double, double, double *, double *);
EXTERN_MSC int spline_(int *istart, int *nn, double *x, double *u, double *s, double *a);
EXTERN_MSC int evals_(int *istart, double *y, int *nn, double *x, double *u, double *s, double *eval);
EXTERN_MSC int find_fft_length(int n);
EXTERN_MSC void aastretch (fcomplex **fdata, int ipatch, int nrows, int num_valid_az, int num_rng_bins, float coef);
EXTERN_MSC void acpatch (void *API, fcomplex **data, int nrows, double delr, double fd, double fdd, double fddd);
EXTERN_MSC void conv2d (float *rdat, int *ni, int *nj, float *filt, int *nif, int *njf, float *fdat, int *ic, int *jc, float *rnorm);
EXTERN_MSC void do_freq_corr (void *API, struct xcorr xc, int iloc);
EXTERN_MSC void do_time_corr(struct xcorr xc, int iloc);
EXTERN_MSC double calc_time_corr(struct xcorr xc, int ioff, int joff);
EXTERN_MSC int fft_bins (int num);
EXTERN_MSC void fft_interpolate_1d (void *API, struct FCOMPLEX *in, int N, struct FCOMPLEX *out, int ifactor);
EXTERN_MSC void fft_interpolate_2d(void *API, struct FCOMPLEX *in, int N1, int M1, struct FCOMPLEX *out, int N,  int M, int ifactor);
EXTERN_MSC void print_prm_params(struct PRM p1, struct PRM p2);
EXTERN_MSC void fix_prm_params(struct PRM *p, char *s);
EXTERN_MSC void get_locations(struct xcorr *xc);
EXTERN_MSC void get_params(FILE *fh);
EXTERN_MSC void do_highres_corr(void *API, struct xcorr xc, int iloc);
EXTERN_MSC void intp_coef(int nfilter, float *xintp);
EXTERN_MSC void print_params(struct xcorr *xc);
EXTERN_MSC void set_defaults(struct xcorr *xc);
EXTERN_MSC void parse_command_line(int na, char **a, struct xcorr *xc, int *nfiles, int *input_flag, char *USAGE);
EXTERN_MSC void handle_prm(char **argv, struct xcorr *xc, int nfiles);
EXTERN_MSC void print_results(struct xcorr xc, int iloc);
EXTERN_MSC void print_complex(struct FCOMPLEX *a, int ny, int nx, int real_flag);
EXTERN_MSC void print_float(float *a, int ny, int nx);
EXTERN_MSC void print_double(double *a, int ny, int nx);
EXTERN_MSC void print_int(int *a, int ny, int nx);
EXTERN_MSC void radopp (double *fd, double *fdd, double *fddd, double r, double del);
EXTERN_MSC void read_complex_short(FILE *f, int *d, int iy, int jx, int npx, int npy, int nx);
EXTERN_MSC void read_real_float(FILE *f, int *d, int iy, int jx, int npx, int npy, int nx);
EXTERN_MSC void read_data(struct xcorr xc);
EXTERN_MSC void read_complex_short2float(FILE *f, float *d, int iy, int jx, int npx, int npy, int nx);
EXTERN_MSC void read_optional_args(void *API, int argc, char **argv, struct PRM *tp, int *topoflag, struct PRM *mp, int *modelflag);
EXTERN_MSC void read_xcorr_data(struct xcorr xc, int iloc);
EXTERN_MSC void rmpatch (fcomplex **data, int nrows, double delr, double fd, double fdd, double fddd);
EXTERN_MSC void rng_cmp (void *API, int ranfft, fcomplex *data, fcomplex *ref);
EXTERN_MSC void rng_ref (void *API, int ranfft, float delr, fcomplex *ref1);
EXTERN_MSC void shift (void *API, int ranfft, fcomplex *data, double shift);
EXTERN_MSC int trans_col (void *API, int xnum, int ynum, fcomplex **data);
EXTERN_MSC int read_SLC_short2float(FILE *SLCfile, char *name, short *sdata, fcomplex *cdata, int xdim, int psize, double dfact);
EXTERN_MSC int read_SLC_short2double(FILE *SLCfile, char *name, short *sdata, dcomplex *cdata, int xdim, int psize, double dfact);
EXTERN_MSC void handle_input(char *, struct xcorr *);
EXTERN_MSC void read_params(struct xcorr *, FILE *);
EXTERN_MSC void make_mask(struct xcorr);
EXTERN_MSC void do_highres(struct xcorr, int);
EXTERN_MSC void allocate_arrays(struct xcorr *);

EXTERN_MSC void llt2rat_sub(char * filename, double *target_llt, double *target_rat );
EXTERN_MSC void read_orb(FILE *, struct PRM *, struct ALOS_ORB *);
EXTERN_MSC void write_orb(FILE *, struct ALOS_ORB *);
EXTERN_MSC void set_prm_defaults(struct PRM *prm);
EXTERN_MSC void hermite_c(double *x, double *y, double *z, int nmax, int nval, double xp, double *yp, int *ir);
EXTERN_MSC void polyfit(double *T, double *Y, double *C, int *Mp, int *Np);
EXTERN_MSC void ldr_orbit(struct ALOS_ORB *, struct PRM *);
EXTERN_MSC void calc_dop(struct PRM *);
EXTERN_MSC void conv2d (float *rdat, int *ni, int *nj, float *filt, int *nif, int *njf, float *fdat, int *ic, int *jc, float *rnorm);

/* sbas functions */ 
EXTERN_MSC int parse_command_ts(int64_t agc, char **agv, float *sf ,double *wl, double *theta, double *rng, int64_t *flag_rms, int64_t *flag_dem, int64_t *atm);
EXTERN_MSC int init_array_ts(double *G, double *Gs, float *res, float *dem, float *disp, int64_t n, int64_t m, int64_t xdim, int64_t ydim, int64_t N, int64_t S);
EXTERN_MSC int read_table_data_ts(void *API,FILE *infile,FILE *datefile,char **gfile,char **cfile,int64_t *H,float *bperp,int64_t *flag,float *var,float *phi,int64_t S,int64_t N,int64_t xdim,int64_t ydim,struct GMT_GRID **Out, int64_t *L, double *time);
EXTERN_MSC int init_G_ts(double *G, double *Gs, int64_t N, int64_t S, int64_t m, int64_t n, int64_t *L, int64_t *H, double *time, float sf, float *bperp, double scale);
EXTERN_MSC int64_t lsqlin_sov_ts(int64_t xdim, int64_t ydim, float *disp, float *vel, int64_t *flag, double *d, double *ds,double *time, double *G, double *Gs, double *A, float *var, float *phi, int64_t N, int64_t S,int64_t m, int64_t n, double *work,int64_t lwork, int64_t flag_dem, float *dem, int64_t flag_rms, float *res, int64_t *jpvt, double wl, double *atm_rms);
EXTERN_MSC int write_output_ts(void *API, struct GMT_GRID *Out,int64_t agc,char **agv, int64_t xdim, int64_t ydim, int64_t S, int64_t flag_rms, int64_t flag_dem, float *disp, float *vel, float *res, float *dem, float *screen, double wl, int64_t n_atm, int64_t *L);
EXTERN_MSC int free_memory_ts(int64_t N,float *phi,float *var,char **gfile,char **cfile,float *disp,double *G,double *A,double *Gs,int64_t *H,double *d,double *ds,int64_t *L,float *res,float *vel,double *time,int64_t *flag,float *bperp,float *dem,double *work,int64_t *jpvt,int64_t *hit);
EXTERN_MSC int allocate_memory_ts(int64_t **jpvt,double **work,double **d,double **ds,float **bperp,char ***gfile,char ***cfile,int64_t **L,double **time,int64_t **H,double **G,double **A,double **Gs,int64_t **flag,float **dem,float **res,float **vel,float **phi,float **var,float **disp,int64_t n,int64_t m,int64_t lwork,int64_t ldb,int64_t N,int64_t S,int64_t xdim,int64_t ydim, int64_t **hit);
EXTERN_MSC int remove_ts(float *phi, float *ts, int64_t xdim, int64_t ydim, int64_t N, int64_t S, int64_t *H, int64_t *L);
EXTERN_MSC int sum_intfs(float *phi, int64_t *mark, float *screen, int64_t xdim, int64_t ydim, int64_t N);
EXTERN_MSC int connect(int64_t *L, int64_t *H, double *time, int64_t *hit, int64_t *mark, int64_t N, int64_t S, int64_t n, int64_t mode);
EXTERN_MSC double compute_noise(float *screen,int64_t xdim, int64_t ydim);
EXTERN_MSC int rank_double(double *nums, int64_t *seq, int64_t n);
EXTERN_MSC int apply_screen(float *screen, float *phi, int64_t xdim, int64_t ydim, int64_t N, int64_t *mark);

#endif /* LIB_FUNCTIONS_H */
