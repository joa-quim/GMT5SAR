/*************************************************************************

  snaphu header file
  Written by Curtis W. Chen
  Copyright 2002 Board of Trustees, Leland Stanford Jr. University
  Please see the supporting documentation for terms of use.
  No warranty.

*************************************************************************/

#ifndef SNAPHU_H
#define SNAPHU_H

/* Shut up MSVC " possible loss of data" warnings */ 
#		pragma warning( disable : 4244 )
#		pragma warning( disable : 4305 )

/**********************/
/* defined constants  */
/**********************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#ifdef _WIN32
#	include "../../unistd.h"
#else
#	include <unistd.h>
#endif
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#ifdef _WIN32
#	include <windows.h>
#	include <time.h>
#	include <process.h>
#	include <signal.h>
#	define pid_t int
#else
#	include <sys/wait.h>
#	include <sys/time.h>
#	include <sys/resource.h>
#endif

#define PROGRAMNAME          "snaphu"
#define VERSION              "1.4.2"
#ifdef PI
#undef PI
#endif
#define PI                   3.14159265358979323846
#define TWOPI                6.28318530717958647692
#define SQRTHALF             0.70710678118654752440
#define MAXSTRLEN            512
#define MAXTMPSTRLEN         1024
#define MAXLINELEN           2048
#define TRUE                 1
#define FALSE                0
#define LARGESHORT           32000
#define LARGELONG            2000000000
#define LARGELONGLONG        9000000000000000000
#define LARGEFLOAT           1.0e35
#define VERYFAR              LARGELONG
#define GROUNDROW            -2
#define GROUNDCOL            -2
#define MAXGROUPBASE         LARGELONG
#define ONTREE               1
#define INBUCKET             2
#define NOTINBUCKET          3
#define POSINCR              0
#define NEGINCR              1
#define NOCOSTSHELF          -LARGESHORT
#define MINSCALARCOST        1
#define INITARRSIZE          500
#define NEWNODEBAGSTEP       500
#define CANDIDATEBAGSTEP     500
#define NEGBUCKETFRACTION    1.0
#define POSBUCKETFRACTION    1.0
#define CLIPFACTOR           0.6666666667
#define DEF_OUTFILE          "snaphu.out"
#define DEF_SYSCONFFILE      ""     /* "/usr/local/snaphu/snaphu.conf" */
#define DEF_WEIGHTFILE       ""     /* "snaphu.weight" */
#define DEF_AMPFILE          ""     /* "snaphu.amp" */
#define DEF_AMPFILE2         ""     /* "snaphu.amp" */
#define DEF_MAGFILE          ""     /* "snaphu.mag" */
#define DEF_CORRFILE         ""     /* "snaphu.corr" */
#define DEF_ESTFILE          ""     /* "snaphu.est" */
#define DEF_COSTINFILE       ""
#define DEF_INITFILE         ""
#define DEF_FLOWFILE         ""
#define DEF_EIFILE           ""
#define DEF_ROWCOSTFILE      ""
#define DEF_COLCOSTFILE      ""
#define DEF_MSTROWCOSTFILE   ""
#define DEF_MSTCOLCOSTFILE   ""
#define DEF_MSTCOSTSFILE     ""
#define DEF_CORRDUMPFILE     ""
#define DEF_RAWCORRDUMPFILE  ""
#define DEF_CONNCOMPFILE     ""
#define DEF_COSTOUTFILE      ""
#define DEF_LOGFILE          ""
#define MAXITERATION         5000
#define NEGSHORTRANGE        SHRT_MIN
#define POSSHORTRANGE        SHRT_MAX
#define MAXRES               SCHAR_MAX
#define MINRES               SCHAR_MIN
#define PROBCOSTP            (-16)
#define NULLFILE             "/dev/null"
#define DEF_ERRORSTREAM      stderr
#define DEF_OUTPUTSTREAM     stdout
#define DEF_VERBOSESTREAM    NULL
#define DEF_COUNTERSTREAM    NULL
#define DEF_INITONLY         FALSE
#define DEF_INITMETHOD       MSTINIT
#define DEF_UNWRAPPED        FALSE
#define DEF_REGROWCONNCOMPS  FALSE
#define DEF_EVAL             FALSE
#define DEF_WEIGHT           1
#define DEF_COSTMODE         TOPO
#define DEF_VERBOSE          FALSE
#define DEF_AMPLITUDE        TRUE
#define AUTOCALCSTATMAX      0
#define USEMAXCYCLEFRACTION  (-123)
#define COMPLEX_DATA         1         /* file format */
#define FLOAT_DATA           2         /* file format */
#define ALT_LINE_DATA        3         /* file format */
#define ALT_SAMPLE_DATA      4         /* file format */
#define ABNORMAL_EXIT        1         /* exit code */
#define NORMAL_EXIT          0         /* exit code */
#define DUMP_PATH            "/tmp/"   /* default location for writing dumps */
#define NARMS                8         /* number of arms for Despeckle() */
#define ARMLEN               5         /* length of arms for Despeckle() */
#define KEDGE                5         /* length of edge detection window */
#define ARCUBOUND            200       /* capacities for cs2 */
#define MSTINIT              1         /* initialization method */
#define MCFINIT              2         /* initialization method */
#define BIGGESTDZRHOMAX      10000.0
#define SECONDSPERPIXEL      0.000001  /* for delay between thread creations */
#define MAXTHREADS           64
#define TMPTILEDIRROOT       "snaphu_tiles_"
#define TILEDIRMODE          511
#define TMPTILEROOT          "tmptile_"
#define TMPTILECOSTSUFFIX    "cost_"
#define TMPTILEOUTFORMAT     ALT_LINE_DATA
#define REGIONSUFFIX         "_regions"
#define LOGFILEROOT          "tmptilelog_"
#define RIGHT                1
#define DOWN                 2
#define LEFT                 3
#define UP                   4
#define TILEDPSICOLFACTOR    0.8
#define ZEROCOSTARC          -LARGELONG
#define PINGPONG             2
#define SINGLEANTTRANSMIT    1
#define NOSTATCOSTS          0
#define TOPO                 1
#define DEFO                 2
#define SMOOTH               3


/* SAR and geometry parameter defaults */

#define DEF_ORBITRADIUS      7153000.0
#define DEF_ALTITUDE         0.0
#define DEF_EARTHRADIUS      6378000.0
#define DEF_BASELINE         150.0
#define DEF_BASELINEANGLE    (1.25*PI)
#define DEF_BPERP            0
#define DEF_TRANSMITMODE     PINGPONG
#define DEF_NLOOKSRANGE      1
#define DEF_NLOOKSAZ         5
#define DEF_NLOOKSOTHER      1
#define DEF_NCORRLOOKS       23.8
#define DEF_NCORRLOOKSRANGE  3  
#define DEF_NCORRLOOKSAZ     15
#define DEF_NEARRANGE        831000.0
#define DEF_DR               8.0
#define DEF_DA               20.0 
#define DEF_RANGERES         10.0
#define DEF_AZRES            6.0
#define DEF_LAMBDA           0.0565647


/* scattering model defaults */

#define DEF_KDS              0.02
#define DEF_SPECULAREXP      8.0
#define DEF_DZRCRITFACTOR    2.0
#define DEF_SHADOW           FALSE
#define DEF_DZEIMIN          -4.0
#define DEF_LAYWIDTH         16 
#define DEF_LAYMINEI         1.25
#define DEF_SLOPERATIOFACTOR 1.18
#define DEF_SIGSQEI          100.0


/* decorrelation model parameters */

#define DEF_DRHO             0.005
#define DEF_RHOSCONST1       1.3
#define DEF_RHOSCONST2       0.14
#define DEF_CSTD1            0.4
#define DEF_CSTD2            0.35
#define DEF_CSTD3            0.06
#define DEF_DEFAULTCORR      0.01
#define DEF_RHOMINFACTOR     1.3


/* pdf model parameters */

#define DEF_DZLAYPEAK        -2.0
#define DEF_AZDZFACTOR       0.99
#define DEF_DZEIFACTOR       4.0 
#define DEF_DZEIWEIGHT       0.5 
#define DEF_DZLAYFACTOR      1.0
#define DEF_LAYCONST         0.9
#define DEF_LAYFALLOFFCONST  2.0
#define DEF_SIGSQSHORTMIN    1
#define DEF_SIGSQLAYFACTOR   0.1


/* deformation mode parameters */

#define DEF_DEFOAZDZFACTOR   1.0
#define DEF_DEFOTHRESHFACTOR 1.2
#define DEF_DEFOMAX          1.2
#define DEF_SIGSQCORR        0.05
#define DEF_DEFOLAYCONST     0.9


/* algorithm parameters */

#define DEF_FLIPPHASESIGN    FALSE
#define DEF_MAXFLOW          4
#define DEF_KROWEI           65
#define DEF_KCOLEI           257
#define DEF_KPARDPSI         7
#define DEF_KPERPDPSI        7
#define DEF_THRESHOLD        0.001
#define DEF_INITDZR          2048.0
#define DEF_INITDZSTEP       100.0
#define DEF_MAXCOST          1000.0
#define DEF_COSTSCALE        100.0 
#define DEF_COSTSCALEAMBIGHT 80.0 
#define DEF_DNOMINCANGLE     0.01
#define DEF_SRCROW           -1
#define DEF_SRCCOL           -1
#define DEF_P                PROBCOSTP
#define DEF_NSHORTCYCLE      200
#define DEF_MAXNEWNODECONST  0.0008
#define DEF_MAXCYCLEFRACTION 0.00001
#define DEF_SOURCEMODE       0
#define DEF_MAXNFLOWCYCLES   USEMAXCYCLEFRACTION
#define DEF_INITMAXFLOW      9999
#define INITMAXCOSTINCR      200
#define NOSTATINITMAXFLOW    15
#define DEF_ARCMAXFLOWCONST  3
#define DEF_DUMPALL          FALSE
#define DUMP_INITFILE        "snaphu.init"
#define DUMP_FLOWFILE        "snaphu.flow"
#define DUMP_EIFILE          "snaphu.ei"
#define DUMP_ROWCOSTFILE     "snaphu.rowcost"
#define DUMP_COLCOSTFILE     "snaphu.colcost"
#define DUMP_MSTROWCOSTFILE  "snaphu.mstrowcost"
#define DUMP_MSTCOLCOSTFILE  "snaphu.mstcolcost"
#define DUMP_MSTCOSTSFILE    "snaphu.mstcosts"
#define DUMP_CORRDUMPFILE    "snaphu.corr"
#define DUMP_RAWCORRDUMPFILE "snaphu.rawcorr"
#define INCRCOSTFILEPOS      "snaphu.incrcostpos"
#define INCRCOSTFILENEG      "snaphu.incrcostneg"
#define DEF_CS2SCALEFACTOR   8


/* default tile parameters */

#define DEF_NTILEROW         1
#define DEF_NTILECOL         1
#define DEF_ROWOVRLP         0
#define DEF_COLOVRLP         0
#define DEF_PIECEFIRSTROW    1
#define DEF_PIECEFIRSTCOL    1
#define DEF_PIECENROW        0
#define DEF_PIECENCOL        0
#define DEF_TILECOSTTHRESH   500
#define DEF_MINREGIONSIZE    100
#define DEF_NTHREADS         1
#define DEF_SCNDRYARCFLOWMAX 8
#define DEF_TILEEDGEWEIGHT   2.5
#define DEF_ASSEMBLEONLY     FALSE
#define DEF_RMTMPTILE        FALSE


/* default connected component parameters */
#define DEF_MINCONNCOMPFRAC  0.01
#define DEF_CONNCOMPTHRESH   300
#define DEF_MAXNCOMPS        32


/* default file formats */

#define DEF_INFILEFORMAT              COMPLEX_DATA
#define DEF_UNWRAPPEDINFILEFORMAT     ALT_LINE_DATA
#define DEF_MAGFILEFORMAT             FLOAT_DATA
#define DEF_OUTFILEFORMAT             ALT_LINE_DATA
#define DEF_CORRFILEFORMAT            ALT_LINE_DATA
#define DEF_ESTFILEFORMAT             ALT_LINE_DATA
#define DEF_AMPFILEFORMAT             ALT_SAMPLE_DATA

/* command-line usage help strings */

#define OPTIONSHELPFULL\
 "usage:  snaphu [options] infile linelength [options]\n"\
 "options:\n"\
 "  -t              use topography mode costs (default)\n"\
 "  -d              use deformation mode costs\n"\
 "  -s              use smooth-solution mode costs\n"\
 "  -f <filename>   read configuration parameters from file\n"\
 "  -o <filename>   write output to file\n"\
 "  -a <filename>   read amplitude data from file\n"\
 "  -A <filename>   read power data from file\n"\
 "  -m <filename>   read interferogram magnitude data from file\n"\
 "  -c <filename>   read correlation data from file\n"\
 "  -e <filename>   read coarse unwrapped-phase estimate from file\n"\
 "  -w <filename>   read scalar weights from file\n"\
 "  -b <decimal>    perpendicular baseline (meters, topo mode only)\n"\
 "  -p <decimal>    Lp-norm parameter p\n"\
 "  -i              do initialization and exit\n"\
 "  -n              do not use statistical costs (with -p or -i)\n"\
 "  -u              infile is already unwrapped; initialization not needed\n"\
 "  -q              quantify cost of unwrapped input file then exit\n"\
 "  -g <filename>   grow connected components mask and write to file\n"\
 "  -G <filename>   grow connected components mask for unwrapped input\n"\
 "  -l <filename>   log runtime parameters to file\n"\
 "  -v              give verbose output\n"\
 "  --mst           use MST algorithm for initialization (default)\n"\
 "  --mcf           use MCF algorithm for initialization\n"\
 "  --aa <filename1> <filename2>    read amplitude from next two files\n"\
 "  --AA <filename1> <filename2>    read power from next two files\n"\
 "  --costinfile <filename>         read statistical costs from file\n"\
 "  --costoutfile <filename>        write statistical costs to file\n"\
 "  --tile <nrow> <ncol> <rowovrlp> <colovrlp>  unwrap as nrow x ncol tiles\n"\
 "  --nproc <integer>               number of processors used in tile mode\n"\
 "  --assemble <dirname>            assemble unwrapped tiles in dir\n"\
 "  --piece <firstrow> <firstcol> <nrow> <ncol>  unwrap subset of image\n" \
 "  --debug, --dumpall              dump all intermediate data arrays\n"\
 "  --copyright, --info             print copyright and bug report info\n"\
 "  -h, --help                      print this help text\n"\
 "\n"

#define OPTIONSHELPBRIEF\
 "usage:  snaphu [options] infile linelength [options]\n"\
 "most common options:\n"\
 "  -t              use topography mode costs (default)\n"\
 "  -d              use deformation mode costs\n"\
 "  -s              use smooth-solution mode costs\n"\
 "  -f <filename>   read configuration parameters from file\n"\
 "  -o <filename>   write output to file\n"\
 "  -a <filename>   read amplitude data from file\n"\
 "  -c <filename>   read correlation data from file\n"\
 "  -b <decimal>    perpendicular baseline (meters)\n"\
 "  -i              do initialization and exit\n"\
 "  -l <filename>   log runtime parameters to file\n"\
 "  -v              give verbose output\n"\
 "  --mst           use MST algorithm for initialization (default)\n"\
 "  --mcf           use MCF algorithm for initialization\n"\
 "\n"\
 "type snaphu -h for a complete list of options\n"\
 "\n"

#define COPYRIGHT\
 "Copyright 2002 Board of Trustees, Leland Stanford Jr. University\n"\
 "\n"\
 "Except as noted below, permission to use, copy, modify, and\n"\
 "distribute, this software and its documentation for any purpose is\n"\
 "hereby granted without fee, provided that the above copyright notice\n"\
 "appear in all copies and that both that copyright notice and this\n"\
 "permission notice appear in supporting documentation, and that the\n"\
 "name of the copyright holders be used in advertising or publicity\n"\
 "pertaining to distribution of the software with specific, written\n"\
 "prior permission, and that no fee is charged for further distribution\n"\
 "of this software, or any modifications thereof.  The copyright holder\n"\
 "makes no representations about the suitability of this software for\n"\
 "any purpose.  It is provided \"as is\" without express or implied\n"\
 "warranty.\n"\
 "\n"\
 "THE COPYRIGHT HOLDER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS\n"\
 "SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND\n"\
 "FITNESS, IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY\n"\
 "SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER\n"\
 "RESULTING FROM LOSS OF USE, DATA, PROFITS, QPA OR GPA, WHETHER IN AN\n"\
 "ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT\n"\
 "OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.\n"\
 "\n"\
 "The parts of this software derived from the CS2 minimum cost flow\n"\
 "solver written by A. V. Goldberg and B. Cherkassky are governed by the\n"\
 "terms of the copyright holder of that software.  Permission has been\n"\
 "granted to use and distrubute that software for strictly noncommercial\n"\
 "purposes as part of this package, provided that the following\n"\
 "copyright notice from the original distribution and URL accompany the\n"\
 "software:\n"\
 "\n"\
 "  COPYRIGHT C 1995 IG Systems, Inc.  Permission to use for\n"\
 "  evaluation purposes is granted provided that proper\n"\
 "  acknowledgments are given.  For a commercial licence, contact\n"\
 "  igsys@eclipse.net.\n"\
 "\n"\
 "  This software comes with NO WARRANTY, expressed or implied. By way\n"\
 "  of example, but not limitation, we make no representations of\n"\
 "  warranties of merchantability or fitness for any particular\n"\
 "  purpose or that the use of the software components or\n"\
 "  documentation will not infringe any patents, copyrights,\n"\
 "  trademarks, or other rights.\n"\
 "\n"\
 "  http://www.igsystems.com/cs2\n"\
 "\n"\
 "\n"\
 "Send snaphu bug reports to Curtis W. Chen (curtis@nova.stanford.edu).\n"\
 "\n"


/********************/
/* type definitions */
/********************/

/* node data structure */
typedef struct nodeST{
  short row,col;                /* row, col of this node */
  uint64_t level;          /* tree level */
  struct nodeST *next;          /* ptr to next node in thread or bucket */
  struct nodeST *prev;          /* ptr to previous node in thread or bucket */
  struct nodeST *pred;          /* parent node in tree */
  int64_t group;                   /* for marking label */
  int64_t incost,outcost;          /* costs to, from root of tree */
}nodeT;


/* arc cost data structure */
typedef struct costST{
  short offset;                 /* offset of wrapped phase gradient from 0 */
  short sigsq;                  /* variance due to decorrelation */
  short dzmax;                  /* largest discontinuity on shelf */
  short laycost;                /* cost of layover discontinuity shelf */
}costT;


/* arc cost data structure for smooth costs */
typedef struct smoothcostST{
  short offset;                 /* offset of wrapped phase gradient from 0 */
  short sigsq;                  /* variance due to decorrelation */
}smoothcostT;


/* incremental cost data structure */
typedef struct incrcostST{
  short poscost;                /* cost for positive flow increment */
  short negcost;                /* cost for negative flow increment */
}incrcostT;


/* arc candidate data structure */
typedef struct candidateST{
  nodeT *from, *to;             /* endpoints of candidate arc */
  int64_t violation;               /* magnitude of arc violation */
  short arcrow,arccol;          /* indexes into arc arrays */
  signed char arcdir;           /* direction of arc (1=fwd, -1=rev) */
}candidateT;


/* bucket data structure */
typedef struct bucketST{
  int64_t size;                    /* number of buckets in list */
  int64_t curr;                    /* current bucket index */
  int64_t maxind;                  /* maximum bucket index */
  int64_t minind;                  /* smallest (possibly negative) bucket index */
  nodeT **bucket;               /* array of first nodes in each bucket */
  nodeT **bucketbase;           /* real base of bucket array */
  signed char wrapped;          /* flag denoting wrapped circular buckets */
}bucketT;


/* secondary arc data structure */
typedef struct scndryarcST{
  short arcrow;                 /* row of arc in secondary network array */
  short arccol;                 /* col of arc in secondary network array */
  nodeT *from;                  /* secondary node at tail of arc */
  nodeT *to;                    /* secondary node at head of arc */
  signed char fromdir;          /* direction from which arc enters head */
}scndryarcT;


/* supplementary data structure for secondary nodes */
typedef struct nodesuppST{
  short row;                    /* row of node in primary network problem */
  short col;                    /* col of node in primary network problem */
  nodeT **neighbornodes;        /* pointers to neighboring secondary nodes */
  scndryarcT **outarcs;         /* pointers to secondary arcs to neighbors */
  short noutarcs;               /* number of arcs from this node */
}nodesuppT;


/* run-time parameter data structure */
typedef struct paramST{

  /* SAR and geometry parameters */
  double orbitradius;     /* radius of platform orbit (meters) */
  double altitude;        /* SAR altitude (meters) */
  double earthradius;     /* radius of earth (meters) */
  double bperp;           /* nominal perpendiuclar baseline (meters) */
  signed char transmitmode; /* transmit mode (PINGPONG or SINGLEANTTRANSMIT) */
  double baseline;        /* baseline length (meters, always postive) */
  double baselineangle;   /* baseline angle above horizontal (rad) */
  int64_t nlooksrange;       /* number of looks in range for input data */ 
  int64_t nlooksaz;          /* number of looks in azimuth for input data */ 
  int64_t nlooksother;       /* number of nonspatial looks for input data */ 
  double ncorrlooks;      /* number of independent looks in correlation est */
  int64_t ncorrlooksrange;   /* number of looks in range for correlation */ 
  int64_t ncorrlooksaz;      /* number of looks in azimuth for correlation */ 
  double nearrange;       /* slant range to near part of swath (meters) */
  double dr;              /* range bin spacing (meters) */
  double da;              /* azimuth bin spacing (meters) */
  double rangeres;        /* range resolution (meters) */
  double azres;           /* azimuth resolution (meters) */
  double lambda;          /* wavelength (meters) */

  /* scattering model parameters */
  double kds;             /* ratio of diffuse to specular scattering */
  double specularexp;     /* power specular scattering component */
  double dzrcritfactor;   /* fudge factor for linearizing scattering model */
  signed char shadow;     /* allow discontinuities from shadowing */
  double dzeimin;         /* lower limit for backslopes (if shadow = FALSE) */
  int64_t laywidth;          /* width of window for summing layover brightness */
  double layminei;        /* threshold brightness for assuming layover */
  double sloperatiofactor;/* fudge factor for linearized scattering slopes */
  double sigsqei;         /* variance (dz, meters) due to uncertainty in EI */

  /* decorrelation model parameters */
  double drho;            /* step size of correlation-slope lookup table */
  double rhosconst1,rhosconst2;/* for calculating rho0 in biased rho */
  double cstd1,cstd2,cstd3;/* for calculating correlation power given nlooks */
  double defaultcorr;     /* default correlation if no correlation file */
  double rhominfactor;    /* threshold for setting unbiased correlation to 0 */

  /* pdf model parameters */
  double dzlaypeak;       /* range pdf peak for no discontinuity when bright */
  double azdzfactor;      /* fraction of dz in azimuth vs. rnage */
  double dzeifactor;      /* nonlayover dz scale factor */
  double dzeiweight;      /* weight to give dz expected from intensity */
  double dzlayfactor;     /* layover regime dz scale factor */
  double layconst;        /* normalized constant pdf of layover edge */
  double layfalloffconst; /* factor of sigsq for layover cost increase */
  int64_t sigsqshortmin;     /* min short value for costT variance */
  double sigsqlayfactor;  /* fration of ambiguityheight^2 for layover sigma */

  /* deformation mode parameters */
  double defoazdzfactor;  /* scale for azimuth ledge in defo cost function */
  double defothreshfactor;/* factor of rho0 for discontinuity threshold */
  double defomax;         /* max discontinuity (cycles) from deformation */
  double sigsqcorr;       /* variance in measured correlation */
  double defolayconst;    /* layconst for deformation mode */

  /* algorithm parameters */
  signed char eval;       /* evaluate unwrapped input file if TRUE */
  signed char unwrapped;  /* input file is unwrapped if TRUE */
  signed char regrowconncomps;/* grow connected components and exit if TRUE */
  signed char initonly;   /* exit after initialization if TRUE */
  signed char initmethod; /* MST or MCF initialization */
  signed char costmode;   /* statistical cost mode */
  signed char dumpall;    /* dump intermediate files */
  signed char verbose;    /* print verbose output */
  signed char amplitude;  /* intensity data is amplitude, not power */
  signed char havemagnitude; /* flag to create correlation from other inputs */
  signed char flipphasesign; /* flag to flip phase and flow array signs */
  int64_t initmaxflow;       /* maximum flow for initialization */
  int64_t arcmaxflowconst;   /* units of flow past dzmax to use for initmaxflow */
  int64_t maxflow;           /* max flow for tree solve looping */
  int64_t krowei, kcolei;    /* size of boxcar averaging window for mean ei */
  int64_t kpardpsi;          /* length of boxcar for mean wrapped gradient */
  int64_t kperpdpsi;         /* width of boxcar for mean wrapped gradient */
  double threshold;       /* thershold for numerical dzrcrit calculation */
  double initdzr;         /* initial dzr for numerical dzrcrit calc. (m) */
  double initdzstep;      /* initial stepsize for spatial decor slope calc. */
  double maxcost;         /* min and max float values for cost arrays */
  double costscale;       /* scale factor for discretizing to integer costs */
  double costscaleambight;/* ambiguity height for auto costs caling */
  double dnomincangle;    /* step size for range-varying param lookup table */
  int64_t srcrow,srccol;     /* source node location */
  double p;               /* power for Lp-norm solution (less than 0 is MAP) */
  int64_t nshortcycle;       /* number of points for one cycle in short int dz */
  double maxnewnodeconst; /* number of nodes added to tree on each iteration */
  int64_t maxnflowcycles;    /* max number of cycles to consider nflow done */
  double maxcyclefraction;/* ratio of max cycles to pixels */
  int64_t sourcemode;        /* 0, -1, or 1, determines how tree root is chosen */
  int64_t cs2scalefactor;    /* scale factor for cs2 initialization (eg, 3-30) */

  /* tiling parameters */
  int64_t ntilerow;          /* number of tiles in azimuth */
  int64_t ntilecol;          /* number of tiles in range */
  int64_t rowovrlp;          /* pixels of overlap between row tiles */
  int64_t colovrlp;          /* pixels of overlap between column tiles */
  int64_t piecefirstrow;     /* first row (indexed from 1) for piece mode */
  int64_t piecefirstcol;     /* first column (indexed from 1) for piece mode */
  int64_t piecenrow;         /* number of rows for piece mode */
  int64_t piecencol;         /* number of rows for piece mode */
  int64_t tilecostthresh;    /* maximum cost within single reliable tile region */
  int64_t minregionsize;     /* minimum number of pixels in a region */
  int64_t nthreads;          /* number of parallel processes to run */
  int64_t scndryarcflowmax;  /* max flow increment for which to keep cost data */
  double tileedgeweight;  /* weight applied to tile-edge secondary arc costs */
  signed char assembleonly; /* flag for assemble-only (no unwrap) mode */
  signed char rmtmptile;  /* flag for removing temporary tile files */
  char tiledir[MAXSTRLEN];/* directory for temporary tile files */

  /* connected component parameters */
  double minconncompfrac; /* min fraction of pixels in connected component */
  int64_t conncompthresh;    /* cost threshold for connected component */
  int64_t maxncomps;         /* max number of connected components */

  
}paramT;


/* input file name data structure */
typedef struct infileST{
  char infile[MAXSTRLEN];             /* input interferogram */
  char magfile[MAXSTRLEN];            /* interferogram magnitude (optional) */
  char ampfile[MAXSTRLEN];            /* image amplitude or power file */
  char ampfile2[MAXSTRLEN];           /* second amplitude or power file */
  char weightfile[MAXSTRLEN];         /* arc weights */
  char corrfile[MAXSTRLEN];           /* correlation file */
  char estfile[MAXSTRLEN];            /* unwrapped estimate */
  char costinfile[MAXSTRLEN];         /* file from which cost data is read */
  signed char infileformat;           /* input file format */
  signed char unwrappedinfileformat;  /* input file format if unwrapped */
  signed char magfileformat;          /* interferogram magnitude file format */
  signed char corrfileformat;         /* correlation file format */
  signed char weightfileformat;       /* weight file format */
  signed char ampfileformat;          /* amplitude file format */
  signed char estfileformat;          /* unwrapped-estimate file format */
}infileT;


/* output file name data structure */
typedef struct outfileST{
  char outfile[MAXSTRLEN];            /* unwrapped output */
  char initfile[MAXSTRLEN];           /* unwrapped initialization */
  char flowfile[MAXSTRLEN];           /* flows of unwrapped solution */
  char eifile[MAXSTRLEN];             /* despckled, normalized intensity */
  char rowcostfile[MAXSTRLEN];        /* statistical azimuth cost array */
  char colcostfile[MAXSTRLEN];        /* statistical range cost array */
  char mstrowcostfile[MAXSTRLEN];     /* scalar initialization azimuth costs */
  char mstcolcostfile[MAXSTRLEN];     /* scalar initialization range costs */
  char mstcostsfile[MAXSTRLEN];       /* scalar initialization costs (all) */
  char corrdumpfile[MAXSTRLEN];       /* correlation coefficient magnitude */
  char rawcorrdumpfile[MAXSTRLEN];    /* correlation coefficient magnitude */
  char conncompfile[MAXSTRLEN];       /* connected component map or mask */
  char costoutfile[MAXSTRLEN];        /* file to which cost data is written */
  char logfile[MAXSTRLEN];            /* file to which parmeters are logged */
  signed char outfileformat;          /* output file format */
}outfileT;


/* tile parameter data structure */
typedef struct tileparamST{
  int64_t firstcol;          /* first column of tile to process (index from 0) */
  int64_t ncol;              /* number of columns in tile to process */
  int64_t firstrow;          /* first row of tile to process (index from 0) */
  int64_t nrow;              /* number of rows in tile to process */
}tileparamT;


/* type for total cost of solution (may overflow int64_t) */
typedef double totalcostT;            /* typedef int64_t int64_t totalcostT; */
#define INITTOTALCOST LARGEFLOAT      /* #define INITTOTALCOST LARGELONGLONG */



/***********************/
/* function prototypes */
/***********************/

/* functions in snaphu.c */

void Unwrap(infileT *infiles, outfileT *outfiles, paramT *params, 
	    int64_t linelen, int64_t nlines);
void UnwrapTile(infileT *infiles, outfileT *outfiles, paramT *params, 
		tileparamT *tileparams, int64_t nlines, int64_t linelen);


/* functions in snaphu_tile.c */

void SetupTile(int64_t nlines, int64_t linelen, paramT *params, 
	       tileparamT *tileparams, outfileT *outfiles, 
	       outfileT *tileoutfiles, int64_t tilerow, int64_t tilecol);
void GrowRegions(void **costs, short **flows, int64_t nrow, int64_t ncol, 
		 incrcostT **incrcosts, outfileT *outfiles, paramT *params);
void GrowConnCompsMask(void **costs, short **flows, int64_t nrow, int64_t ncol, 
		       incrcostT **incrcosts, outfileT *outfiles, 
		       paramT *params);
int64_t ThickenCosts(incrcostT **incrcosts, int64_t nrow, int64_t ncol);
nodeT *RegionsNeighborNode(nodeT *node1, int64_t *arcnumptr, nodeT **nodes, 
			   int64_t *arcrowptr, int64_t *arccolptr, 
			   int64_t nrow, int64_t ncol);
void ClearBuckets(bucketT *bkts);
void MergeRegions(nodeT **nodes, nodeT *source, int64_t *regionsizes, 
		  int64_t closestregion, int64_t nrow, int64_t ncol);
void RenumberRegion(nodeT **nodes, nodeT *source, int64_t newnum, 
		    int64_t nrow, int64_t ncol);
void AssembleTiles(outfileT *outfiles, paramT *params, 
		   int64_t nlines, int64_t linelen);
void ReadNextRegion(int64_t tilerow, int64_t tilecol, int64_t nlines, int64_t linelen,
		    outfileT *outfiles, paramT *params, 
		    short ***nextregionsptr, float ***nextunwphaseptr,
		    void ***nextcostsptr, 
		    int64_t *nextnrowptr, int64_t *nextncolptr);
void SetTileReadParams(tileparamT *tileparams, int64_t nexttilenlines, 
		       int64_t nexttilelinelen, int64_t tilerow, int64_t tilecol, 
		       int64_t nlines, int64_t linelen, paramT *params);
void ReadEdgesAboveAndBelow(int64_t tilerow, int64_t tilecol, int64_t nlines, 
			    int64_t linelen, paramT *params, outfileT *outfiles, 
			    short *regionsabove, short *regionsbelow,
			    float *unwphaseabove, float *unwphasebelow,
			    void *costsabove, void *costsbelow);
void TraceRegions(short **regions, short **nextregions, short **lastregions, 
		  short *regionsabove, short *regionsbelow, float **unwphase, 
		  float **nextunwphase, float **lastunwphase, 
		  float *unwphaseabove, float *unwphasebelow, void **costs, 
		  void **nextcosts, void **lastcosts, void *costsabove, 
		  void *costsbelow, int64_t prevnrow, int64_t prevncol, int64_t tilerow,
		  int64_t tilecol, int64_t nrow, int64_t ncol, nodeT **scndrynodes,
		  nodesuppT **nodesupp, scndryarcT **scndryarcs, 
		  int64_t ***scndrycosts, short *nscndrynodes, 
		  short *nscndryarcs, int64_t *totarclens, short **bulkoffsets, 
		  paramT *params);
int64_t FindNumPathsOut(nodeT *from, paramT *params, int64_t tilerow, int64_t tilecol, 
		     int64_t nnrow, int64_t nncol, short **regions, 
		     short **nextregions, short **lastregions,
		     short *regionsabove, short *regionsbelow, int64_t prevncol);
void RegionTraceCheckNeighbors(nodeT *from, nodeT **nextnodeptr, 
			       nodeT **primarynodes, short **regions, 
			       short **nextregions, short **lastregions, 
			       short *regionsabove, short *regionsbelow,  
			       int64_t tilerow, int64_t tilecol, int64_t nnrow, 
			       int64_t nncol, nodeT **scndrynodes, 
			       nodesuppT **nodesupp, scndryarcT **scndryarcs, 
			       int64_t *nnewnodesptr, int64_t *nnewarcsptr, 
			       int64_t flowmax, int64_t nrow, int64_t ncol, 
			       int64_t prevnrow, int64_t prevncol, paramT *params, 
			       void **costs, void **rightedgecosts, 
			       void **loweredgecosts, void **leftedgecosts, 
			       void **upperedgecosts, short **flows, 
			       short **rightedgeflows, short **loweredgeflows,
			       short **leftedgeflows, short **upperedgeflows,
			       int64_t ***scndrycosts, 
			       nodeT ***updatednontilenodesptr, 
			       int64_t *nupdatednontilenodesptr, 
			       int64_t *updatednontilenodesizeptr,
			       short **inontilenodeoutarcptr, 
			       int64_t *totarclenptr);
void SetUpperEdge(int64_t ncol, int64_t tilerow, int64_t tilecol, void **voidcosts, 
		  void *voidcostsabove, float **unwphase, 
		  float *unwphaseabove, void **voidupperedgecosts, 
		  short **upperedgeflows, paramT *params, short **bulkoffsets);
void SetLowerEdge(int64_t nrow, int64_t ncol, int64_t tilerow, int64_t tilecol, 
		  void **voidcosts, void *voidcostsbelow, 
		  float **unwphase, float *unwphasebelow, 
		  void **voidloweredgecosts, short **loweredgeflows, 
		  paramT *params, short **bulkoffsets);
void SetLeftEdge(int64_t nrow, int64_t prevncol, int64_t tilerow, int64_t tilecol, 
		 void **voidcosts, void **voidlastcosts, float **unwphase, 
		 float **lastunwphase, void **voidleftedgecosts, 
		 short **leftedgeflows, paramT *params, short **bulkoffsets);
void SetRightEdge(int64_t nrow, int64_t ncol, int64_t tilerow, int64_t tilecol, 
		  void **voidcosts, void **voidnextcosts, 
		  float **unwphase, float **nextunwphase, 
		  void **voidrightedgecosts, short **rightedgeflows, 
		  paramT *params, short **bulkoffsets);
void TraceSecondaryArc(nodeT *primaryhead, nodeT **scndrynodes, 
		       nodesuppT **nodesupp, scndryarcT **scndryarcs, 
		       int64_t ***scndrycosts, int64_t *nnewnodesptr, 
		       int64_t *nnewarcsptr, int64_t tilerow, int64_t tilecol, 
		       int64_t flowmax, int64_t nrow, int64_t ncol, 
		       int64_t prevnrow, int64_t prevncol, paramT *params, 
		       void **tilecosts, void **rightedgecosts, 
		       void **loweredgecosts, void **leftedgecosts,
		       void **upperedgecosts, short **tileflows, 
		       short **rightedgeflows, short **loweredgeflows, 
		       short **leftedgeflows, short **upperedgeflows,
		       nodeT ***updatednontilenodesptr, 
		       int64_t *nupdatednontilenodesptr, 
		       int64_t *updatednontilenodesizeptr,
		       short **inontilenodeoutarcptr, int64_t *totarclenptr);
nodeT *FindScndryNode(nodeT **scndrynodes, nodesuppT **nodesupp, 
		      int64_t tilenum, int64_t primaryrow, int64_t primarycol);
void IntegrateSecondaryFlows(int64_t linelen, int64_t nlines, nodeT **scndrynodes, 
			     nodesuppT **nodesupp, scndryarcT **scndryarcs, 
			     short *nscndryarcs, short **scndryflows, 
			     short **bulkoffsets, outfileT *outfiles, 
			     paramT *params);
void ParseSecondaryFlows(int64_t tilenum, short *nscndryarcs, short **tileflows, 
			 short **regions, short **scndryflows, 
			 nodesuppT **nodesupp, scndryarcT **scndryarcs, 
			 int64_t nrow, int64_t ncol, int64_t ntilerow, int64_t ntilecol,
			 paramT *params);


/* functions in snaphu_solver.c */

int64_t TreeSolve(nodeT **nodes, nodesuppT **nodesupp, nodeT *ground, 
	       nodeT *source, candidateT **candidatelistptr, 
	       candidateT **candidatebagptr, int64_t *candidatelistsizeptr,
	       int64_t *candidatebagsizeptr, bucketT *bkts, short **flows, 
	       void **costs, incrcostT **incrcosts, nodeT ***apexes, 
	       signed char **iscandidate, int64_t ngroundarcs, int64_t nflow, 
	       float **mag, float **wrappedphase, char *outfile, 
	       int64_t nnoderow, short *nnodesperrow, int64_t narcrow, 
	       short *narcsperrow, int64_t nrow, int64_t ncol,
	       outfileT *outfiles, paramT *params);
void AddNewNode(nodeT *from, nodeT *to, int64_t arcdir, bucketT *bkts, 
		int64_t nflow, incrcostT **incrcosts, int64_t arcrow, int64_t arccol, 
		paramT *params);
void CheckArcReducedCost(nodeT *from, nodeT *to, nodeT *apex, 
			 int64_t arcrow, int64_t arccol, int64_t arcdir, 
			 int64_t nflow, nodeT **nodes, nodeT *ground, 
			 candidateT **candidatebagptr, 
			 int64_t *candidatebagnextptr, 
			 int64_t *candidatebagsizeptr, incrcostT **incrcosts, 
			 signed char **iscandidate, paramT *params);
int64_t InitTree(nodeT *source, nodeT **nodes, nodesuppT **nodesupp, 
	      nodeT *ground, int64_t ngroundarcs, bucketT *bkts, int64_t nflow, 
	      incrcostT **incrcosts, nodeT ***apexes, 
	      signed char **iscandidate, int64_t nnoderow, short *nnodesperrow, 
	      int64_t narcrow, short *narcsperrow, int64_t nrow, int64_t ncol, 
	      paramT *params);
nodeT *FindApex(nodeT *from, nodeT *to);
int CandidateCompare(const void *c1, const void *c2);
nodeT *NeighborNodeGrid(nodeT *node1, int64_t arcnum, int64_t *upperarcnumptr,
			nodeT **nodes, nodeT *ground, int64_t *arcrowptr, 
			int64_t *arccolptr, int64_t *arcdirptr, int64_t nrow, 
			int64_t ncol, nodesuppT **nodesupp);
nodeT *NeighborNodeNonGrid(nodeT *node1, int64_t arcnum, int64_t *upperarcnumptr,
			   nodeT **nodes, nodeT *ground, int64_t *arcrowptr, 
			   int64_t *arccolptr, int64_t *arcdirptr, int64_t nrow, 
			   int64_t ncol, nodesuppT **nodesupp);
void GetArcGrid(nodeT *from, nodeT *to, int64_t *arcrow, int64_t *arccol, 
		int64_t *arcdir, int64_t nrow, int64_t ncol, nodesuppT **nodesupp);
void GetArcNonGrid(nodeT *from, nodeT *to, int64_t *arcrow, int64_t *arccol, 
		   int64_t *arcdir, int64_t nrow, int64_t ncol, nodesuppT **nodesupp);
void NonDegenUpdateChildren(nodeT *startnode, nodeT *lastnode, 
			    nodeT *nextonpath, int64_t dgroup, 
			    int64_t ngroundarcs, int64_t nflow, nodeT **nodes,
			    nodesuppT **nodesupp, nodeT *ground, 
			    nodeT ***apexes, incrcostT **incrcosts, 
			    int64_t nrow, int64_t ncol, paramT *params);
void InitNetwork(short **flows, int64_t *ngroundarcsptr, int64_t *ncycleptr, 
		 int64_t *nflowdoneptr, int64_t *mostflowptr, int64_t *nflowptr, 
		 int64_t *candidatebagsizeptr, candidateT **candidatebagptr, 
		 int64_t *candidatelistsizeptr, candidateT **candidatelistptr, 
		 signed char ***iscandidateptr, nodeT ****apexesptr, 
		 bucketT **bktsptr, int64_t *iincrcostfileptr, 
		 incrcostT ***incrcostsptr, nodeT ***nodesptr, nodeT *ground, 
		 int64_t *nnoderowptr, short **nnodesperrowptr, int64_t *narcrowptr,
		 short **narcsperrowptr, int64_t nrow, int64_t ncol, 
		 signed char *notfirstloopptr, totalcostT *totalcostptr,
		 paramT *params);
void InitNodeNums(int64_t nrow, int64_t ncol, nodeT **nodes, nodeT *ground);
void InitBuckets(bucketT *bkts, nodeT *source, int64_t nbuckets);
void InitNodes(int64_t nrow, int64_t ncol, nodeT **nodes, nodeT *ground);
void BucketInsert(nodeT *node, int64_t ind, bucketT *bkts);
void BucketRemove(nodeT *node, int64_t ind, bucketT *bkts);
nodeT *ClosestNode(bucketT *bkts);
nodeT *ClosestNodeCircular(bucketT *bkts);
nodeT *MinOutCostNode(bucketT *bkts);
nodeT *SelectSource(nodeT **nodes, nodeT *ground, int64_t nflow, 
		    short **flows, int64_t ngroundarcs, 
		    int64_t nrow, int64_t ncol, paramT *params);
short GetCost(incrcostT **incrcosts, int64_t arcrow, int64_t arccol, 
	      int64_t arcdir);
int64_t ReCalcCost(void **costs, incrcostT **incrcosts, int64_t flow, 
		int64_t arcrow, int64_t arccol, int64_t nflow, int64_t nrow, 
		paramT *params);
void SetupIncrFlowCosts(void **costs, incrcostT **incrcosts, short **flows,
			int64_t nflow, int64_t nrow, int64_t narcrow, 
			short *narcsperrow, paramT *params);
totalcostT EvaluateTotalCost(void **costs, short **flows, int64_t nrow, int64_t ncol,
			     short *narcsperrow,paramT *params);
void MSTInitFlows(float **wrappedphase, short ***flowsptr, 
		  short **mstcosts, int64_t nrow, int64_t ncol, 
		  nodeT ***nodes, nodeT *ground, int64_t maxflow);
void SolveMST(nodeT **nodes, nodeT *source, nodeT *ground, 
	      bucketT *bkts, short **mstcosts, signed char **residue, 
	      signed char **arcstatus, int64_t nrow, int64_t ncol);
int64_t DischargeTree(nodeT *source, short **mstcosts, short **flows,
		   signed char **residue, signed char **arcstatus, 
		   nodeT **nodes, nodeT *ground, int64_t nrow, int64_t ncol);
signed char ClipFlow(signed char **residue, short **flows, 
		     short **mstcosts, int64_t nrow, int64_t ncol, 
		     int64_t maxflow);
void MCFInitFlows(float **wrappedphase, short ***flowsptr, short **mstcosts, 
		  int64_t nrow, int64_t ncol, int64_t cs2scalefactor);


/* functions in snaphu_cost.c */

void BuildCostArrays(void ***costsptr, short ***mstcostsptr, 
		     float **mag, float **wrappedphase, 
		     float **unwrappedest, int64_t linelen, int64_t nlines, 
		     int64_t nrow, int64_t ncol, paramT *params, 
		     tileparamT *tileparams, infileT *infiles, 
		     outfileT *outfiles);
void **BuildStatCostsTopo(float **wrappedphase, float **mag, 
			  float **unwrappedest, float **pwr, 
			  float **corr, short **rowweight, short **colweight,
			  int64_t nrow, int64_t ncol, tileparamT *tileparams, 
			  outfileT *outfiles, paramT *params);
void **BuildStatCostsDefo(float **wrappedphase, float **mag, 
			  float **unwrappedest, float **corr, 
			  short **rowweight, short **colweight,
			  int64_t nrow, int64_t ncol, tileparamT *tileparams, 
			  outfileT *outfiles, paramT *params);
void **BuildStatCostsSmooth(float **wrappedphase, float **mag, 
			    float **unwrappedest, float **corr, 
			    short **rowweight, short **colweight,
			    int64_t nrow, int64_t ncol, tileparamT *tileparams, 
			    outfileT *outfiles, paramT *params);
void GetIntensityAndCorrelation(float **mag, float **wrappedphase, 
				float ***pwrptr, float ***corrptr, 
				infileT *infiles, int64_t linelen, int64_t nlines,
				int64_t nrow, int64_t ncol, outfileT *outfiles, 
				paramT *params, tileparamT *tileparams);
void RemoveMean(float **ei, int64_t nrow, int64_t ncol, 
                int64_t krowei, int64_t kcolei);
float *BuildDZRCritLookupTable(double *nominc0ptr, double *dnomincptr, 
			       int64_t *tablesizeptr, tileparamT *tileparams,
			       paramT *params);
double SolveDZRCrit(double sinnomincangle, double cosnomincangle, 
		    paramT *params, double threshold);
void SolveEIModelParams(double *slope1ptr, double *slope2ptr, 
			double *const1ptr, double *const2ptr, 
			double dzrcrit, double dzr0, double sinnomincangle, 
			double cosnomincangle, paramT *params);
double EIofDZR(double dzr, double sinnomincangle, double cosnomincangle,
	       paramT *params);
float **BuildDZRhoMaxLookupTable(double nominc0, double dnominc, 
				 int64_t nominctablesize, double rhomin, 
				 double drho, int64_t nrho, paramT *params);
double CalcDZRhoMax(double rho, double nominc, paramT *params, 
		    double threshold);
void CalcCostTopo(void **costs, int64_t flow, int64_t arcrow, int64_t arccol, 
		  int64_t nflow, int64_t nrow, paramT *params, 
		  int64_t *poscostptr, int64_t *negcostptr);
void CalcCostDefo(void **costs, int64_t flow, int64_t arcrow, int64_t arccol, 
		  int64_t nflow, int64_t nrow, paramT *params, 
		  int64_t *poscostptr, int64_t *negcostptr);
void CalcCostSmooth(void **costs, int64_t flow, int64_t arcrow, int64_t arccol, 
		    int64_t nflow, int64_t nrow, paramT *params, 
		    int64_t *poscostptr, int64_t *negcostptr);
void CalcCostL0(void **costs, int64_t flow, int64_t arcrow, int64_t arccol, 
		int64_t nflow, int64_t nrow, paramT *params, 
		int64_t *poscostptr, int64_t *negcostptr);
void CalcCostL1(void **costs, int64_t flow, int64_t arcrow, int64_t arccol, 
		int64_t nflow, int64_t nrow, paramT *params, 
		int64_t *poscostptr, int64_t *negcostptr);
void CalcCostL2(void **costs, int64_t flow, int64_t arcrow, int64_t arccol, 
		int64_t nflow, int64_t nrow, paramT *params, 
		int64_t *poscostptr, int64_t *negcostptr);
void CalcCostLP(void **costs, int64_t flow, int64_t arcrow, int64_t arccol, 
		int64_t nflow, int64_t nrow, paramT *params, 
		int64_t *poscostptr, int64_t *negcostptr);
void CalcCostNonGrid(void **costs, int64_t flow, int64_t arcrow, int64_t arccol, 
		     int64_t nflow, int64_t nrow, paramT *params, 
		     int64_t *poscostptr, int64_t *negcostptr);
int64_t EvalCostTopo(void **costs, short **flows, int64_t arcrow, int64_t arccol,
		  int64_t nrow, paramT *params);
int64_t EvalCostDefo(void **costs, short **flows, int64_t arcrow, int64_t arccol,
		  int64_t nrow, paramT *params);
int64_t EvalCostSmooth(void **costs, short **flows, int64_t arcrow, int64_t arccol,
		    int64_t nrow, paramT *params);
int64_t EvalCostL0(void **costs, short **flows, int64_t arcrow, int64_t arccol,
		int64_t nrow, paramT *params);
int64_t EvalCostL1(void **costs, short **flows, int64_t arcrow, int64_t arccol,
		int64_t nrow, paramT *params);
int64_t EvalCostL2(void **costs, short **flows, int64_t arcrow, int64_t arccol,
		int64_t nrow, paramT *params);
int64_t EvalCostLP(void **costs, short **flows, int64_t arcrow, int64_t arccol,
		int64_t nrow, paramT *params);
int64_t EvalCostNonGrid(void **costs, short **flows, int64_t arcrow, int64_t arccol, 
		     int64_t nrow, paramT *params);
void CalcInitMaxFlow(paramT *params, void **costs, int64_t nrow, int64_t ncol);


/* functions in snaphu_util.c */

int IsTrue(char *str);
int IsFalse(char *str);
signed char SetBooleanSignedChar(signed char *boolptr, char *str);
double ModDiff(double f1, double f2);
void WrapPhase(float **wrappedphase, int64_t nrow, int64_t ncol);
void CalcWrappedRangeDiffs(float **dpsi, float **avgdpsi, float **wrappedphase,
			   int64_t kperpdpsi, int64_t kpardpsi,
			   int64_t nrow, int64_t ncol);
void CalcWrappedAzDiffs(float **dpsi, float **avgdpsi, float **wrappedphase,
			int64_t kperpdpsi, int64_t kpardpsi, int64_t nrow, int64_t ncol);
void CycleResidue(float **phase, signed char **residue, 
		  int nrow, int ncol);
void CalcFlow(float **phase, short ***flowsptr, int64_t nrow, int64_t ncol);
void IntegratePhase(float **psi, float **phi, short **flows,
		    int64_t nrow, int64_t ncol);
float **ExtractFlow(float **unwrappedphase, short ***flowsptr, 
		    int64_t nrow, int64_t ncol);
void FlipPhaseArraySign(float **arr, paramT *params, int64_t nrow, int64_t ncol);
void FlipFlowArraySign(short **arr, paramT *params, int64_t nrow, int64_t ncol);
void **Get2DMem(int nrow, int ncol, int psize, size_t size);
void **Get2DRowColMem(int64_t nrow, int64_t ncol, int psize, size_t size);
void **Get2DRowColZeroMem(int64_t nrow, int64_t ncol, int psize, size_t size);
void *MAlloc(size_t size);
void *CAlloc(size_t nitems, size_t size);
void *ReAlloc(void *ptr, size_t size);
void Free2DArray(void **array, unsigned int nrow);
void Set2DShortArray(short **arr, int64_t nrow, int64_t ncol, int64_t value);
signed char ValidDataArray(float **arr, int64_t nrow, int64_t ncol);
signed char IsFinite(double d);
int64_t LRound(double a);
int64_t Short2DRowColAbsMax(short **arr, int64_t nrow, int64_t ncol);
float LinInterp1D(float *arr, double index, int64_t nelem);
float LinInterp2D(float **arr, double rowind, double colind , 
                  int64_t nrow, int64_t ncol);
void Despeckle(float **mag, float ***ei, int64_t nrow, int64_t ncol);
float **MirrorPad(float **array1, int64_t nrow, int64_t ncol, int64_t krow, int64_t kcol);
void BoxCarAvg(float **avgarr, float **padarr, int64_t nrow, int64_t ncol, 
	       int64_t krow, int64_t kcol);
char *StrNCopy(char *dest, const char *src, size_t n);
void FlattenWrappedPhase(float **wrappedphase, float **unwrappedest, 
			 int64_t nrow, int64_t ncol);
void Add2DFloatArrays(float **arr1, float **arr2, int64_t nrow, int64_t ncol);
int StringToDouble(char *str, double *d);
int StringToLong(char *str, int64_t *l);
void CatchSignals(void (*SigHandler)(int));
void SetDump(int signum);
void KillChildrenExit(int signum);
void SignalExit(int signum);
void StartTimers(time_t *tstart, double *cputimestart);
void DisplayElapsedTime(time_t tstart, double cputimestart);
int LongCompare(const void *c1, const void *c2);

/* functions in snaphu_io.c */

void SetDefaults(infileT *infiles, outfileT *outfiles, paramT *params);
void ProcessArgs(int argc, char *argv[], infileT *infiles, outfileT *outfiles,
		 int64_t *ncolptr, paramT *params);
void CheckParams(infileT *infiles, outfileT *outfiles, 
		 int64_t linelen, int64_t nlines, paramT *params);
void ReadConfigFile(char *conffile, infileT *infiles, outfileT *outfiles,
		    int64_t *ncolptr, paramT *params);
void WriteConfigLogFile(int argc, char *argv[], infileT *infiles, 
			outfileT *outfiles, int64_t linelen, paramT *params);
void LogStringParam(FILE *fp, char *key, char *value);
void LogBoolParam(FILE *fp, char *key, signed char boolvalue);
void LogFileFormat(FILE *fp, char *key, signed char fileformat);
int64_t GetNLines(infileT *infiles, int64_t linelen);
void WriteOutputFile(float **mag, float **unwrappedphase, char *outfile, 
		     outfileT *outfiles, int64_t nrow, int64_t ncol);
FILE *OpenOutputFile(char *outfile, char *realoutfile);
void WriteAltLineFile(float **mag, float **phase, char *outfile, 
		      int64_t nrow, int64_t ncol);
void WriteAltSampFile(float **arr1, float **arr2, char *outfile, 
		      int64_t nrow, int64_t ncol);
void Write2DArray(void **array, char *filename, int64_t nrow, int64_t ncol, 
		  size_t size);
void Write2DRowColArray(void **array, char *filename, int64_t nrow, 
			int64_t ncol, size_t size);
void ReadInputFile(infileT *infiles, float ***magptr, float ***wrappedphaseptr,
		   short ***flowsptr, int64_t linelen, int64_t nlines, 
		   paramT *params, tileparamT *tileparams);
void ReadMagnitude(float **mag, infileT *infiles, int64_t linelen, int64_t nlines, 
		   tileparamT *tileparams);
void ReadUnwrappedEstimateFile(float ***unwrappedestptr, infileT *infiles, 
			       int64_t linelen, int64_t nlines, 
			       paramT *params, tileparamT *tileparams);
void ReadWeightsFile(short ***weightsptr,char *weightfile, 
		     int64_t linelen, int64_t nlines, tileparamT *tileparams);
void ReadIntensity(float ***pwrptr, float ***pwr1ptr, float ***pwr2ptr, 
		   infileT *infiles, int64_t linelen, int64_t nlines, 
		   paramT *params, tileparamT *tileparams);
void ReadCorrelation(float ***corrptr, infileT *infiles,
		     int64_t linelen, int64_t nlines, tileparamT *tileparams);
void ReadAltLineFile(float ***mag, float ***phase, char *alfile, 
		     int64_t linelen, int64_t nlines, tileparamT *tileparams);
void ReadAltLineFilePhase(float ***phase, char *alfile, 
			  int64_t linelen, int64_t nlines, tileparamT *tileparams);
void ReadComplexFile(float ***mag, float ***phase, char *rifile, 
		     int64_t linelen, int64_t nlines, tileparamT *tileparams);
void Read2DArray(void ***arr, char *infile, int64_t linelen, int64_t nlines, 
		 tileparamT *tileparams, size_t elptrsize, size_t elsize);
void ReadAltSampFile(float ***arr1, float ***arr2, char *infile,
		     int64_t linelen, int64_t nlines, tileparamT *tileparams);
void Read2DRowColFile(void ***arr, char *filename, int64_t linelen, int64_t nlines, 
		      tileparamT *tileparams, size_t size);
void Read2DRowColFileRows(void ***arr, char *filename, int64_t linelen, 
			  int64_t nlines, tileparamT *tileparams, size_t size);
void SetDumpAll(outfileT *outfiles, paramT *params);
void SetStreamPointers(void);
void SetVerboseOut(paramT *params);
void ChildResetStreamPointers(pid_t pid, int64_t tilerow, int64_t tilecol,
			      paramT *params);
void DumpIncrCostFiles(incrcostT **incrcosts, int64_t iincrcostfile, 
		       int64_t nflow, int64_t nrow, int64_t ncol);
void MakeTileDir(paramT *params, outfileT *outfiles);
void ParseFilename(char *filename, char *path, char *basename);


/* functions in snaphu_cs2.c  */

void SolveCS2(signed char **residue, short **mstcosts, int64_t nrow, int64_t ncol, 
	      int64_t cs2scalefactor, short ***flowsptr);



/*******************************************/
/* global (external) variable declarations */
/*******************************************/

/* flags used for signal handling */
extern char dumpresults_global;
extern char requestedstop_global;

/* ouput stream pointers */
/* sp0=error messages, sp1=status output, sp2=verbose, sp3=verbose counter */
extern FILE *sp0, *sp1, *sp2, *sp3;

/* node pointer for marking arc not on tree in apex array */
/* this should be treat as a constant */
extern nodeT NONTREEARC[1];

/* pointers to functions which calculate arc costs */
extern void (*CalcCost)();
extern int64_t (*EvalCost)();

/* pointers to functions for tailoring network solver to specific topologies */
extern nodeT *(*NeighborNode)();
extern void (*GetArc)();

#endif	/* end of snaphu.h */




