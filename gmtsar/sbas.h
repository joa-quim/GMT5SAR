 /*      $Id: sbas.h 39 2016-06-18 03/16/24 Xiaohua Xu $  */

/* sbas functions */
EXTERN_MSC int parse_command_ts(int , char **, float *, double *, double *, double *, int *, int *); 
EXTERN_MSC int init_array_ts(double *, double *, float *, float *, float *, int , int , int , int , int, int);
EXTERN_MSC int read_table_data_ts(void *,FILE *,FILE *,char **,char **,int *,float *,int *,float *,float *,int,int,int,int,struct  GMT_GRID **,int *,double *); 
EXTERN_MSC int init_G_ts(double *, double *, int, int, int, int, int *, int *, double *, float, float *, double);
EXTERN_MSC int lsqlin_sov_ts(int,int,float *,float *,int *,double *,double *,double *,double *,double *,double *,float *,float *,int,int,int,int,double *,int,int,float *,int,float *,int *,double);
EXTERN_MSC int write_output_ts(void *API, struct GMT_GRID *Out,int agc,char **agv, int xdim, int ydim, int S, int flag_rms, int flag_dem, float *disp, float *vel, float *res, float *dem, float *screen, double wl);
EXTERN_MSC int free_memory_ts(int N,float *phi,float *var,char **gfile,char **cfile,float *disp,double *G,double *A,double *Gs,int *H,double *d,double *ds,int *L,float *res,float *vel,double *time,int *flag,float *bperp,float *dem,double *work,int *jpvt,int *hit);
EXTERN_MSC int allocate_memory_ts(int **jpvt,double **work,double **d,double **ds,float **bperp,char ***gfile,char ***cfile,int **L,double **time,int **H,double **G,double **A,double **Gs,int **flag,float **dem,float **res,float **vel,float **phi,float **var,float **disp,int n,int m,int lwork,int ldb,int N,int S,int xdim,int ydim, int **hit);

