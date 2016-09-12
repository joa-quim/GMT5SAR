 /*      $Id: sbas.h 39 2016-06-18 03/16/24 Xiaohua Xu $  */

/* sbas functions */
EXTERN_MSC int parse_command_ts(int , char **, float *, double *, double *, double *, int *, int *); 
EXTERN_MSC int allocate_memory_ts(int **,double **,double **,double **,float **,char ***,char ***,int **,double **,int **,double **,double **,double **,int **,float **,float **,float **,float **,float **,float **,int,int,int,int,int,int,int,int);
EXTERN_MSC int init_array_ts(double *, double *, float *, float *, float *, int , int , int , int , int, int);
EXTERN_MSC int read_table_data_ts(void *,FILE *,FILE *,char **,char **,int *,float *,int *,float *,float *,int,int,int,int,struct  GMT_GRID **,int *,double *); 
EXTERN_MSC int init_G_ts(double *, double *, int, int, int, int, int *, int *, double *, float, float *, double);
EXTERN_MSC int lsqlin_sov_ts(int,int,float *,float *,int *,double *,double *,double *,double *,double *,double *,float *,float *,int,int,int,int,double *,int,int,float *,int,float *,int *,double);
EXTERN_MSC int write_output_ts(void *,struct GMT_GRID *,int,char **,int,int,int,int,int,float *,float *,float *,float *); 
EXTERN_MSC int free_memory_ts(int ,float *,float *,char **,char **,float *,double *,double *,double *,int *,double *,double *,int *,float *,float *,double *,int *,float *,float *,double *,int *);
