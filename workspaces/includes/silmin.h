#ifndef _Silmin_h
#define _Silmin_h

/*
MELTS Source Code: RCS $Log: silmin.h,v $
MELTS Source Code: RCS Revision 1.8  2007/09/13 16:12:02  ghiorso
MELTS Source Code: RCS (1) Revised standard state liquid properties.
MELTS Source Code: RCS (2) Revised standard state solid properties (removed non-Berman) Cp, and
MELTS Source Code: RCS     removed Saxena EOS treatment.  All EOS parameterizations are Vinet.
MELTS Source Code: RCS     Updated K, K', alpha to conform to Knittle (1995) and Fei (1995)
MELTS Source Code: RCS     except where refitted Berman (1988) makes more sense.
MELTS Source Code: RCS (3) Updated code to allow for fusion entropies of liquid components to
MELTS Source Code: RCS     be adjusted (fusion enthalpies are dependent).
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.7  2007/08/26 21:38:31  ghiorso
MELTS Source Code: RCS Normalized residuals (for xMELTS calibration) to the number of atoms in the
MELTS Source Code: RCS endmember mineral formula.  Revised residual-statistics.out file accordingly.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.6  2007/08/23 16:09:39  ghiorso
MELTS Source Code: RCS Database updates.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.5  2007/06/16 01:01:55  ghiorso
MELTS Source Code: RCS Revised EOS regression to have K', K'', and K''' as parameters.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.4  2007/06/08 17:25:43  ghiorso
MELTS Source Code: RCS Added code to allow regression of Ghiorso EOS parameters
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.3  2006/10/20 00:59:22  ghiorso
MELTS Source Code: RCS (1) Made initial modifications for thread safe code.
MELTS Source Code: RCS (2) Added support for XML I/O in batch mode
MELTS Source Code: RCS (3) Added support for Melts-batch listener for eventual integration into VIGMCS
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2006/08/17 16:47:19  ghiorso
MELTS Source Code: RCS Made modifications to protect strings.  These modifications allow removal
MELTS Source Code: RCS of the flag -fwritable-strings during gcc compilation.  This brings the
MELTS Source Code: RCS code up to gcc 4.x standards.
MELTS Source Code: RCS
MELTS Source Code: RCS Other minor rearrangements and cleanup.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2006/08/15 16:57:36  ghiorso
MELTS Source Code: RCS xMELTS gcc 3.x sources
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.3  2005/01/21 18:18:22  cvsaccount
MELTS Source Code: RCS
MELTS Source Code: RCS Added data structures and code to implement coordination number transformations
MELTS Source Code: RCS in the liquid phase EOS model.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2004/09/15 02:47:56  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2004/01/02 19:21:49  cvsaccount
MELTS Source Code: RCS CTserver University of Chicago
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.7  2003/09/27 15:35:22  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.6  2003/05/01 17:33:54  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.5  2002/07/26 01:01:37  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.4  2002/06/28 03:33:12  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.3  2002/01/21 19:48:03  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2002/01/10 02:28:04  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:48:56  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
 * Revision 3.10  1997/06/21  22:49:27  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.9  1997/05/03  20:23:06  ghiorso
 * *** empty log message ***
 *
 * Revision 3.8  1997/03/27  17:03:11  ghiorso
 * *** empty log message ***
 *
 * Revision 3.7  1996/09/24  20:33:21  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 3.6  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 3.5  1995/11/23  22:37:42  ghiorso
 * Final implementation of subsolidus fO2 buffering.
 *
 * Revision 3.4  1995/11/01  22:48:53  ghiorso
 * Implementation of subsolidus options after Asimow.
 * Additional implementation of nepheline solid solutions.
 *
 * Revision 3.3  1995/09/04  23:29:01  ghiorso
 * Added new function copySilminStateStructure, for copy the silminState
 * structure. Called in create_xy_plot.padb.c for allocating history records.
 *
 * Revision 3.2  1995/09/04  19:59:58  ghiorso
 * Update to allow display of bulk composition (in grams) in the text entry
 * fields of the main silmin display. Liquid composition is no longer
 * display here, and is available only through the popup selection.
 *
 * Revision 3.1  1995/08/18  19:13:43  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      SILMIN include file (file: SILMIN.H)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  September 10, 1990   Original Version
**      V1.1    Mark S. Ghiorso  April 27, 1991
**              Enlargement of solid structures and addition of public
**              functions for minerals
**      V1.2-1  Mark S. Ghiorso  September 5, 1991
**              (1) reorganized SilminState and SilminHistory structures
**              (2) added reference to checkStateAgainstInterface() function
**      V1.2-2  Mark S. Ghiorso  September 7, 1991
**              (1) added declaration of getlog10fo2
**      V1.2-3  Mark S. Ghiorso  September 9, 1991
**              Added macros for silminStateChange function
**      V1.2-4  Mark S. Ghiorso  September 14, 1991
**              Changed parameter list for con*** routine declarations and
**              altered solids[] structure for corresponding declaration
**              of (*convert)
**      V1.2-5  Mark S. Ghiorso  September 14, 1991
**              Added declaration of new update* procedures
**              Added fracMass to silminState structure
**                               September 20, 1991
**              (1) Removed wkarea extern declaration
**      V1.2-6  Mark S. Ghiorso  September 23, 1991
**              Altered declaration of getAffinityAndComposition
**      V1.2-7  Mark S. Ghiorso  September 24, 1991
**              Changed parameter list for con*** and test*** routine
**              declarations and altered solids[] structure for corresponding
**              declaration of (*convert) and (*test)
**      V1.2-8  Mark S. Ghiorso  October 12, 1991
**              Added checForCoexistingSolids function and solidCoexist
**              element to silminState structure
**      V1.2-9  Mark S. Ghiorso  October 15, 1991
**              (1) Altered definition of silminState->solidComp to be a
**                  two dimensional array.
**              (2) same for silminState->solidDelta.
**              (3) changed silminState->solidCoexist from pointer to double
**                  to int *(silminState->nSolidCoexist), to count the
**                  number of columns (immiscible solids) in solidComp and
**                  solidDelta arrays
**      V2.1-10 Mark S. Ghiorso  October 24, 1991
**              Added declarations of *Pyx routines for pyroxene solid
**              solutions
**      V2.1-11 Mark S. Ghiorso  November 23, 1991
**              (1) Removed references to debug information on interface
**                  display
**      V2.1-12 Mark S. Ghiorso  December 10, 1991
**              (1) Removed reference to magma_padb
**      V2.1-13 Mark S. Ghiorso  December 14, 1991
**              (1) Altered silminState structure to remove magma entries
**                  and redefine assimilant composition arrays
**      V2.1-14 Mark S. Ghiorso  January 2, 1992
**              (1) modified silminState structure for fractionated solids
**      V2.1-15 Mark S. Ghiorso  January 15, 1992
**              (1) redefined argument list to updateAssimilantPADB()
**      V2.1-16 Mark S. Ghiorso  January 21, 1992
**              (1) added definition of visLiq() for viscosity of the
**                  liquid
**      V2.1-17 Mark S. Ghiorso  February 18, 1992
**              (1) Minor changes for ANSI compliance
**      V2.1-18 Mark S. Ghiorso  April 7, 1992
**              Added *Bio and *Grn solid routines
**              Mark S. Ghiorso  April 8, 1992
**              Added *Mel solid routines
**      V2.1-19 Mark S. Ghiorso  April 16, 1992
**              Added *Cum solid routine
**      V2.2-1  Mark S. Ghiorso  April 20, 1992
**              (1) Added external declaration of muO2Liq() function
**                               April 28, 1992
**              (2) Added typedef and extern declaration of constraint struct
**      V2.2-2  Mark S. Ghiorso  May 1, 1992
**              (1) Added reference enthalpy and volume entries to
**                  silminState structure
**      V2.2-3  Mark S. Ghiorso  May 4, 1992
**              (1) Added *cylSolids to silminState structure
**      V2.3-1  Mark S. Ghiorso  July 11, 1992
**              (1) Removed solid solution function declarations
**              (2) Added (*display) function to solids structure
**      V2.3-2  Mark S. Ghiorso  June 14, 1993
**              Added new macros for fo2 path constraints
**      V2.4-1  Mark S. Ghiorso  July 3, 1993
**              (1) Altered modelParameter[] structure to change value
**                  member to enthalpy member.
**              (2) Altered modelParameter[] structure to add entropy
**                  member
**      V2.4-2  Mark S. Ghiorso July 13, 1993
**              (1) Altered modelParameter[] structure to change active
**                  member to activeH member
**              (2) Altered modelParameter[] structure to add activeS
**                  member
**      V2.4-3  Mark S. Ghiorso July 19, 1993
**              (1) Altered modelParameter[] structure to add activeV
**                  and volume members
**      V3.0-1  Mark S. Ghiorso May 10, 1994
**              (1) altered solid structure declaration to reflect entries
**                  for new isenthalpic, isentropic, isochoric derivatives
**              (2) made other corrections to include isentropic constraints
**                              May 27, 1994
**              (3) Altered constraint structure so that tDelta and pDelta
**                  are replaced by T and P.
**                              May 28, 1994
**              (4) Added fo2 field to constraint structure
**                              June 10, 1994
**              (5) Added new members (dsp*) to silminState structure
**                              June 13, 1994
**              (6) Added external references to (SILMIN_SUPPORT.C):
**                  void correctTforChangeInEnthalpy(void)
**                  void correctTforChangeInEntropy(void)
**                  void correctPforChangeInVolume(void)
**                              June 17, 1994
**              (7) Added declaration of getdlog10fo2dt(), getdlog10fo2dp(),
**                  getd2log10fo2dt2() and getd2log10fo2dp2().
**    V3.1-1  Paul D. Asimow  July 31, 1995
**            (1) Added d3rdm3 option to (*solids.convert)[]
**            (2) Added d3gdr3 option to (*solids.gmix)[]
**            (3) Added **solidDelta to Constraints struct
**--
*/

#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef MINGW
#include <windows.h>
#endif

#include "mthread.h"

/*
 *==============================================================================
 * Error handling
 */
#if defined(MINGW) &&  !defined(__USING_SJLJ_EXCEPTIONS__)
#define USESEH 1
#else
#define USESJLJ 1
#endif

/*
 *==============================================================================
 * Numerical constants
 */
#ifdef MINGW
#ifndef DBL_EPSILON
#define DBL_EPSILON __DBL_EPSILON__
#endif
#endif

#define TAU   DBL_EPSILON /* machine precision               */
#define BIG   DBL_MAX     /* maximum double precision number */
#ifndef TRUE
#define TRUE  1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/*
 *==============================================================================
 * Macros
 */
#ifdef ABS
#undef ABS
#endif
#define ABS(x)   ((x) < 0 ? -(x) : (x))
#ifdef MAX
#undef MAX
#endif
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/*
 *==============================================================================
 * Public liquid solution functions:
 */

void conLiq  (int inpMask, int outMask, double t, double p, double *o,
                            double *m, double *r, double *x, double **dm, double ***dm2,
                            double *logfo2);
int testLiq (int mask, double t, double p, int na, int nr, char **names,
             char **formulas, double *r, double *m);
void dispLiq (int mask, double t, double p, double *x, char **formula);
void setModeToMixingLiq(int flag);
void actLiq  (int mask, double t, double p, double *x, double *a,
                            double *mu, double **dx, double **dw);
void gmixLiq (int mask, double t, double P, double *x, double *gmix,
                            double *dx, double **dx2);
void hmixLiq (int mask, double t, double P, double *x, double *hmix, double *dw);
void smixLiq (int mask, double t, double P, double *x, double *smix,
                            double *dx, double **dx2, double *dw);
void cpmixLiq(int mask, double t, double P, double *x, double *cpmix,
                            double *dt, double *dx);
void vmixLiq (int mask, double t, double P, double *x, double *vmix,
                            double *dx, double **dx2, double *dt, double *dp, double *dt2,
                            double *dtdp, double *dp2, double *dxdt, double *dxdp, double *dw);
void muO2Liq(int mask, double t, double p, double *m, double *muO2, double *dm,
                            double *dt, double *dp, double **d2m, double *d2mt, double *d2mp,
                            double *d2t2, double *d2tp, double *d2p2);
void visLiq  (int mask, double t, double P, double *x, double *viscosity);

#ifdef PHMELTS_ADJUSTMENTS
void muH2OLiq(int mask, double t, double p, double *m, double *muH2O, double *dm,
	                     double *dt, double *dp, double **d2m, double *d2mt, double *d2mp,
                             double *d2t2, double *d2tp, double *d2p2);
#endif

/*
 *==============================================================================
 * Argument BITMASKs for public solid solution functions:
 */
#define FIRST       00000001 /* octal for binary 00000000000000000001 */
#define SECOND      00000002 /* octal for binary 00000000000000000010 */
#define THIRD       00000004 /* octal for binary 00000000000000000100 */
#define FOURTH      00000010 /* octal for binary 00000000000000001000 */
#define FIFTH       00000020 /* octal for binary 00000000000000010000 */
#define SIXTH       00000040 /* octal for binary 00000000000000100000 */
#define SEVENTH     00000100 /* octal for binary 00000000000001000000 */
#define EIGHTH      00000200 /* octal for binary 00000000000010000000 */
#define NINTH       00000400 /* octal for binary 00000000000100000000 */
#define TENTH       00001000 /* octal for binary 00000000001000000000 */
#define ELEVENTH    00002000 /* octal for binary 00000000010000000000 */
#define TWELFTH     00004000 /* octal for binary 00000000100000000000 */
#define THIRTEENTH  00010000 /* octal for binary 00000001000000000000 */
#define FOURTEENTH  00020000 /* octal for binary 00000010000000000000 */
#define FIFTEENTH   00040000 /* octal for binary 00000100000000000000 */
#define SIXTEENTH   00100000 /* octal for binary 00001000000000000000 */
#define SEVENTEENTH 00200000 /* octal for binary 00010000000000000000 */
#define EIGHTEENTH  00400000 /* octal for binary 00100000000000000000 */
#define NINETEENTH  01000000 /* octal for binary 01000000000000000000 */
#define TWENTIETH   02000000 /* octal for binary 10000000000000000000 */

/*
 *==============================================================================
 * Algorithmic constants
 */
#define MASSIN  1.0e-5   /* initial mass of included phases           */
#define MASSOUT 1.0e-7   /* minimal mass of included phases           */
#define SCALET  1000.0   /* scaling factor for T in isenthalpic calcs */
#define SCALEP  1000.0   /* scaling factor for P in isochoric calcs   */
#define ITERMX  100      /* maximum number of quadratic interations   */

#define R       8.3143   /* gas constant (J/K)                        */
#define TR      298.15   /* reference temperature (K)                 */
#define PR      1.0      /* reference pressure (bars)                 */

/*
 *==============================================================================
 * Thermodynamic data structures
 *   (1) These are type definitions, and are incorporated into actual
 *       data structures which characterize each phase
 */
#define EOS_BERMAN  0
#define EOS_GHIORSO 1
#define EOS_KRESS   2
#define EOS_VINET   3
#define EOS_SAXENA  4

typedef union _EOS {
    struct { double a0, a1, a2, a3, b0, b1, b2, b3, dKdP, d2KdTdP; }  Saxena;
    struct { double alpha, c, dcdt, cp, mw, d2vdp2, d3vdp3, d4vdp4; } Ghiorso;
    struct { double v1, v2, v3, v4; }                                 Berman;
    struct { double dvdt, dvdp, d2vdtp, d2vdp2; }                     Kress;
    struct { double alpha, K, Kp;   }                                 Vinet;
} EOS;

#define CP_BERMAN 0
#define CP_SAXENA 1

typedef union _CP {
    struct { double k0, k1, k2, k3, Tt, deltah, l1, l2; } Berman;
    struct { double a, b, c, d, e, g, h; }                Saxena;
} CP;

typedef struct _thermoRef {
    double h, s, v;        /* enthalpy (J), entropy (J/K), volume (J/bar)   */
    int    cp_type;        /* Cp formulation for solid                      */
    CP     cp;             /* Cp  parameters                                */
    int    eos_type;       /* EOS formulation for solid                     */
    EOS    eos;            /* EOS parameters                                */
} ThermoRef;

typedef struct _thermoLiq {
    double v;              /* volume of liquid (at Tr)                      */
    int    eos_type;       /* EOS formulation for liquid                    */
    EOS    eos;            /* EOS parameters                                */
    double tfus, sfus, cp; /* T fusion (K), S fusion (J/K), liquid Cp (J/K) */
    double tglass;         /* T glass transition (K)                        */
} ThermoLiq;

typedef struct _thermoData {
    double g, h, s, v, cp; /* gibbs energy (J), enthalpy (J), entropy (J/K) */
                         /* volume (J/bar) and heat capacity (J/K)        */
    double dcpdt, dvdt, dvdp;       /* various T, P derivatives */
    double d2vdt2, d2vdp2, d2vdtdp;
} ThermoData;

typedef struct _modelParameters {
    const char *label;   /* label for the model parameter                      */
    double     enthalpy; /* value of the model parameter (J)		     */
    double     entropy;  /* value of the model parameter (J/K)		     */
    double     volume;   /* value of the model parameter (J/bar)  	     */
    int        activeH;  /* active enthalpy parameter (TRUE), inactive (FALSE) */
    int        activeS;  /* active entropy  parameter (TRUE), inactive (FALSE) */
    int        activeV;  /* active volume   parameter (TRUE), inactive (FALSE) */
    int        activeF;  /* active fusion   parameter (TRUE), inactive (FALSE) */
} ModelParameters;

extern
ModelParameters *modelParameters; /* data structure
                                     containing nls*(nls-1)/2+nls entries */

typedef struct _eosModelParameters {
    const char *label;     /* label for the model parameter                    */
    double     Kp;         /* value of Kp   for oxide                          */
    double     Kpp;        /* value of Kpp  for oxide                          */
    double     Kppp;       /* value of Kppp for oxide                          */
    int        activeKp;   /* active Kp   parameter (TRUE), inactive (FALSE)   */
    int        activeKpp;  /* active Kpp  parameter (TRUE), inactive (FALSE)   */
    int        activeKppp; /* active Kppp parameter (TRUE), inactive (FALSE)   */
    double     v2;         /* internal storage                                 */
    double     v3;         /* internal storage                                 */
    double     v4;         /* internal storage                                 */
} EosModelParameters;

typedef struct _ssCnModelParameters {
    const char *label;     /* label for the model parameter                    */
    double     enthalpy;   /* Ghiorso-Kress, enthalpy CN[*] - CN[ref]          */
    double     entropy;    /* Ghiorso-Kress, entropy  CN[*] - CN[ref]          */
} SsCnModelParameters;

typedef struct _wCnModelParameters {
    const char *label;     /* label for the model parameter                    */
    double     W[3];       /* Delta W_CN[ref][*] - W_CN[ref][ref]              */
                         /* Delta W_CN[*][ref] - W_CN[ref][ref]              */
		         /* Delta W_CN[*][*]   - W_CN[ref][ref]              */
} WCnModelParameters;

extern double fCN[]; /* configurational collapse coefficiets              */

extern
EosModelParameters eosModelParameters[];     /* data structure containing
                                                               nc entries */
extern
SsCnModelParameters ssCnModelParameters[];   /* data structure containing
                                                               nc entries */

extern
WCnModelParameters wCnModelParameters[];     /* data structure containing
                                                                                                        [nc*(nc-1)/2 entries] */

extern int nCN;                    /* total number of coordination states */

/*
 *==============================================================================
 * Bulk System, Liquid, and Solids data
 */
#define FEO   1
#define FE2O3 2
#define OTHER 0

typedef struct _bulkSystem {
    const char *label;     /* label for oxide compositional variables            */
    int        type;       /* FEO, FE2O3, OTHER oxide                            */
    double     coeff;      /* coefficient if OTHER, for ferrous/ferric calc      */
    double     mw;         /* molecular weights of oxides                        */
    double     *oxToLiq;   /* pointer to an array of length [nc] which converts
                                                        moles of oxides to moles of liquid components      */
    int        *oxToElm;   /* pointer to an array of length [106] which converts
                                                        moles of oxides to moles of elements               */
    double     gk_v;       /* Ghiorso-Kress, v-bar-sub-i                         */
    double     gk_dvdt;    /* Ghiorso-Kress, dvdt-bar-sub-i                      */
    double     gk_c;       /* Ghiorso-Kress, c-sub-i                             */
    double     gk_cXal2o3; /* Ghiorso-Kress, c-sub-i cross coeff with Al2O3      */
    double     gk_dcdt;    /* Ghiorso-Kress, dcdt-bar-sub-i                      */
    double     gk_cp;      /* Ghiorso-Kress, cp model values                     */
    double     gk_d2vdp2;  /* Ghiorso-Kress, d2vdp2-bar-sub-i                    */
    double     gk_d3vdp3;  /* Ghiorso-Kress, d3vdp3-bar-sub-i                    */
    double     gk_d4vdp4;  /* Ghiorso-Kress, d4vdp4-bar-sub-i                    */
} BulkSystem;
extern BulkSystem bulkSystem[];
extern int nc;

typedef struct _liquid {
    const char *label;   /* label for the liquid component                    */
    double     *liqToOx; /* pointer to an array of length [nc] which converts
                                                    component moles to moles of oxides                */
    ThermoRef  ref;      /* reference properties (solid) at TR, PR            */
    ThermoLiq  liq;      /* reference properties (liquid)                     */
    ThermoData fus;      /* reference properties (liquid) at T fusion, PR     */
    ThermoData cur;      /* current thermodynamic data at T, P                */
    double     tfus;     /* Reference value of .liq.tfus (calibration)        */
} Liquid;
extern Liquid *liquid;
extern int nlc;
extern int nls;

extern int liqERRstate;

#define PHASE     1
#define COMPONENT 0

typedef struct _solids {
    const char *label;   /* label for solid phase or solid phase components   */
    int        type;     /* PHASE or COMPONENT                                */
    const char *formula; /* character string formula                          */
    int        inclInClb;/* TRUE/FALSE - include phase in calibration         */
    int        inStdSet; /* TRUE/FALSE - member of standard set of solids     */
    double     *solToOx; /* pointer to an array of length [nc] which converts
                                                    COMPONENT moles to moles of oxides                */
    double     *solToLiq;/* pointer to an array of length [nc] which converts
                                                    COMPONENT moles to moles of liquid components     */
    double     mw;       /* molecular weight pure phase or component          */
    double     nAtoms;   /* number of atoms in the formula unit of the phase  */
    /* Defined if type is COMPONENT, else empty                               */
    ThermoRef  ref;      /* ref thermodynamic data, constants at TR, PR       */
    ThermoData cur;      /* current thermodynamic data                        */
    int        na;       /* if type == PHASE; number of endmember components  */
    int        nr;       /* if type == PHASE; number of independent variables */
    /* Defined if type is PHASE, else pointers to void. common arguments:
     mask  -  bitwise mask for selecting input/output
     t     -  Temperature (K)
     p     -  Pressure (bars)
     *x    -  (pointer to x[]) Array of independent compositional variables
    */
    /* returns TRUE if values are correct/within bounds, else returns FALSE    */
    int (*test) (int mask, double t, double p,
     int    na,       /* # of components in solution     BINARY MASK: 000001 */
     int    nr,       /* # of indep compos variables     BINARY MASK: 000010 */
     char **names,    /* names compon, expected order    BINARY MASK: 000100 */
     char **formulas, /* form of compon, expected order  BINARY MASK: 001000 */
     double *r,       /* indep compositional variables   BINARY MASK: 010000 */
     double *m        /* moles of endmember components   BINARY MASK: 100000 */
    );           /* r[] and m[] are tested for bound constraints               */
    void (*convert) (int inpMask, int outMask, double t, double p,
     double *e,     /* moles of elements               BINARY MASK: 00000001 */
     double *m,     /* moles of endmember components   BINARY MASK: 00000010 */
     double *r,     /* indep compositional variables   BINARY MASK: 00000100 */
     double *x,     /* mole fractions of endmember cmp BINARY MASK: 00001000 */
     double **dm,   /* matrix[i][j]: dr[i]/dm[j]       BINARY MASK: 00010000 */
     double ***d2m, /* cube[i][j][k]: d2r[i]/dm[j]dm[k]BINARY MASK: 00100000 */
     double **dr,   /* matrix[i][j]: dx[i]/dr[j]       BINARY MASK: 01000000 */
     double ****d3m /* 4d[i][j][k][l]: d3r[i]/dm[j]dm[k]dm[l] MASK: 10000000 */
    );
    void (*activity) (int mask, double t, double p, double *x,
     double *a,   /* (pointer to a[]) activities           BINARY MASK: 0001 */
     double *mu,  /* (pointer to mu[]) chemical potentials BINARY MASK: 0010 */
     double **dx  /* (pointer to dx[][]) d(a[])/d(x[])     BINARY MASK: 0100 */
    );            /* exclusion applied to activities if:     BINARY MASK: 1000 */
    void (*gmix) (int mask, double t, double p, double *x,
     double *gmix,  /* Gibbs energy of mixing              BINARY MASK: 0001 */
     double *dx,    /* (pointer to dx[]) d(g)/d(x[])       BINARY MASK: 0010 */
     double **dx2,  /* (pointer to dx2[][]) d2(g)/d(x[])2  BINARY MASK: 0100 */
     double ***dx3  /* (pointer to dx3[][][]) d3(g)/d(x[])3BINARY MASK: 1000 */
    );
    void (*hmix) (int mask, double t, double p, double *x,
     double *hmix  /* Enthalpy of mixing                      BINARY MASK: 1 */
    );
    void (*smix) (int mask, double t, double p, double *x,
     double *smix,  /* Entropy of mixing                    BINARY MASK: 001 */
     double *dx,    /* (pointer to dx[]) d(s)/d(x[])        BINARY MASK: 010 */
     double **dx2   /* (pointer to dx2[][]) d2(s)/d(x[])2   BINARY MASK: 100 */
    );
    void (*cpmix) (int mask, double t, double p, double *x,
     double *cpmix, /* Heat capacity of mixing              BINARY MASK: 001 */
     double *dt,    /* d(cp)/d(t)                           BINARY MASK: 010 */
     double *dx     /* d(cp)/d(x[])                         BINARY MASK: 100 */
    );
    void (*vmix) (int mask, double t, double p, double *x,
     double *vmix, /* Volume of mixing               BINARY MASK: 0000000001 */
     double *dx,   /* (pointer to dx[]) d(v)/d(x[])  BINARY MASK: 0000000010 */
     double **dx2, /* (point dx2[][]) d(v)/d(x[])2   BINARY MASK: 0000000100 */
     double *dt,   /* d(v)/d(t)                      BINARY MASK: 0000001000 */
     double *dp,   /* d(v)/d(p)                      BINARY MASK: 0000010000 */
     double *dt2,  /* d2(v)/d(t)2                    BINARY MASK: 0000100000 */
     double *dtdp, /* d2(v)/d(t)d(p)                 BINARY MASK: 0001000000 */
     double *dp2,  /* d2(v)/d(p)2                    BINARY MASK: 0010000000 */
     double *dxdt, /* d2(v)/d(x[])d(t)               BINARY MASK: 0100000000 */
     double *dxdp  /* d2(v)/d(x[])d(p)               BINARY MASK: 1000000000 */
    );
    void (*display) (int mask, double t, double p, double *x,
     char **formula /* Mineral formula for interface display  BINARY MASK: 1 */
    );
} Solids;
extern Solids *solids;
extern int npc;

typedef struct _oxygen {
    double     *liqToOx; /* pointer to an array of length [nc] which converts
                                                    moles of liquid components to moles of oxygen      */
    double     *solToOx; /* pointer to an array of length [npc] which converts
                                                    moles of solid phase/components to moles of oxygen */
    ThermoRef  ref;      /* reference thermodynamic data, constants at TR, PR  */
    ThermoData cur;      /* current thermodynamic data                         */
} Oxygen;
extern Oxygen oxygen;

/*
 *=============================================================================
 * Status of calculation and intensive and extensive variables for the system
 */

#define MODE_xMELTS           0
#define MODE__MELTS           1
#define MODE_pMELTS           2
#define MODE__MELTSandCO2     3
#define MODE__MELTSandCO2_H2O 4
#define MODE_DEFAULT          MODE_xMELTS

extern int calculationMode;

#define FO2_NONE    0
#define FO2_HM      1
#define FO2_NNO     2
#define FO2_QFM     3
#define FO2_COH     4
#define FO2_IW      5
#define FO2_QFM_P3  6
#define FO2_QFM_P2  7
#define FO2_QFM_P1  8
#define FO2_QFM_M1  9
#define FO2_QFM_M2 10
#define FO2_QFM_M3 11
#define FO2_QFM_M4 12
#define FO2_QFM_M5 13
#define FO2_QFM_M6 14
#define FO2_QFM_M7 15
#define FO2_QFM_M8 16
#define FO2_QFM_M9 17
#define FO2_QFM_P0_5 18
#define FO2_QFM_P1_5 19

#ifdef PHMELTS_ADJUSTMENTS
#define FO2_ABS      -1
#endif

#ifdef PHMELTS_ADJUSTMENTS
#define TEXT_NONE      0
#define TEXT_TABLE     1
#define TEXT_BOTH_0    2
#define TEXT_BOTH_1    3
#define TEXT_ALPHA_0   4
#define TEXT_ALPHA_1   5
#endif

typedef struct _silminState {
    double  *bulkComp;      /* current bulk composition (moles of oxides)      */
    double  *dspBulkComp;   /* displayed bulk composition (grams of oxides)    */
    ThermoData bulkTD;      /* current thermodynamic properties of the system  */

    double  **liquidComp;   /* current liquid composition (moles of liq comp)  */
    int     nLiquidCoexist; /* number of coexisting liquids                    */
    double  **liquidDelta;  /* correction to liquidComp[] from quadratic min   */
    double  liquidMass;     /* current mass of liquids (grams)                 */
    ThermoData liquidTD;    /* current thermodynamic properties of the liquids */
#ifdef PHMELTS_ADJUSTMENTS
    int     incLiquids;     /* number of allowed liquids (default nlc)         */
    int     cylLiquids;     /* current liquid phases suppressed due to cycling */
#else
    int     multipleLiqs;   /* current value of toggle mode (TRUE/FALSE)       */
#endif

    double  **solidComp;    /* current solid composition (moles of endmembers) */
    int     *nSolidCoexist; /* number of coexisting solids (cols of solidComp) */
    double  **solidDelta;   /* correction to solidComp[] from quadratic min    */
    double  solidMass;      /* current mass of solids (grams)                  */
    int     *incSolids;     /* current solid phases to be allowed to precip    */
    int     *cylSolids;     /* current solid phases suppressed due to cycling  */
    ThermoData solidTD;     /* current thermodynamic properties of the solids  */

    /* Note: in Mark's file dspTstart etc. are listed as in C                  */
    /* His version uses dspTstart etc. in K                                    */
    /* alphaMELTS version has dspTstart etc. in displayed units (i.e. C)       */
#ifdef PHMELTS_ADJUSTMENTS
    double  dspTstart;      /* displayed initial temperature (C)               */
    double  dspTstop;       /* displayed final temperature (C)                 */
    double  dspTinc;        /* displayed temperature increment (C)             */
#else
    double  dspTstart;      /* displayed initial temperature (K)               */
    double  dspTstop;       /* displayed final temperature (K)                 */
    double  dspTinc;        /* displayed temperature increment (K)             */
#endif
    double  dspHstop;       /* displayed final enthalpy (J)                    */
    double  dspHinc;        /* displayed enthalpy increment (J)                */
    double  dspSstop;       /* displayed final entropy (J/K)                   */
    double  dspSinc;        /* displayed entropy increment (J/K)               */
    double  T;              /* current temperature (K)                         */

    double  dspPstart;      /* displayed initial pressure (bars) [may be 0]    */
    double  dspPstop;       /* displayed final pressure (bars)                 */
    double  dspPinc;        /* displayed pressure increment (bars)             */
    double  dspVstop;       /* displayed final volume (J/bar)                  */
    double  dspVinc;        /* displayed volume increment (J/bar)              */
    double  dspDPDt;        /* displayed dPdT gradient (bars/K)                */
    double  dspDPDH;        /* displayed dPdH gradient (bars/J)                */
    double  dspDPDS;        /* displayed dPdS gradient (bars-K/J)              */
#ifdef PHMELTS_ADJUSTMENTS
    /* seems to be some inconsistency in use of dspDVDt; use dspDTDV instead   */
    double  dspDTDV;        /* displayed dTdV gradient (J/bar-K)               */
#else
    double  dspDVDt;        /* displayed dVdT gradient (J/bar-K)               */
#endif
    double  P;              /* current pressure (bars) [rounded to >= 1]       */

    double  fo2;            /* current value of fo2 (numeric, base 10 log)     */
    int     fo2Path;        /* current value of fo2 path (i.e. FO2_NONE, etc   */
    double  fo2Delta;       /* offset from fo2Path                             */
    double  oxygen;         /* reference value of O2 content in the system     */

#ifdef PHMELTS_ADJUSTMENTS
    int     fo2Alt;         /* same behaviour as ALPHAMELTS_ALTERNATIVE_FO2    */
    int     fo2Liq;         /* stores fo2Path like ALPHAMELTS_LIQUID_FO2       */
    int     fo2Sol;         /* stores fo2Path for subsolidus start             */
    int     fo2Iter;        /* iterate fo2 on-off so get true equilibrium      */
    int     H2Obuffer;      /* current value of aH2O buffer mode (TRUE/FALSE)  */
    double  aH2O;           /* reference water activity for buffering option   */
    double  refMass;        /* initial mass - used for scaling open system(s)  */
#endif

    int     fractionateSol; /* current value of fractionation mode (TRUE/FALSE)*/
    int     fractionateLiq; /* current value of fractionation mode (TRUE/FALSE)*/
    int     fractionateFlu; /* current value of fractionation mode (TRUE/FALSE)*/
    double  **fracSComp;    /* -> current fractionated compos (moles of end.)  */
    int     *nFracCoexist;  /* number of coexisting frac solids                */
    double  *fracLComp;     /* -> current fractionated compos (moles of end.)  */
    double  fracMass;       /* current mass of fractionated liquid+solids      */

    int     isenthalpic;    /* current value of isenthalpic mode (TRUE/FALSE)  */
    double  refEnthalpy;    /* reference enthalpy of the system                */
    int     isentropic;     /* current value of isentropic mode (TRUE/FALSE)   */
    double  refEntropy;     /* reference entropy of the system                 */
    double  tDelta;         /* correction to temperature from quadratic min    */

    int     isochoric;      /* current value of isochoric mode (TRUE/FALSE)    */
    double  refVolume;      /* reference volume of the system                  */
    double  pDelta;         /* correction to pressure from quadratic min       */

    int     assimilate;     /* current value of assimilation mode (TRUE/FALSE) */
    double  **dspAssimComp; /* -> displayed assimilant compos (units and X)    */
    int     *nDspAssimComp; /* -> number of each type of phase in previous     */
    int     dspAssimUnits;  /* -> displayed  assim units (vol % or wt %)       */
    double  dspAssimMass;   /* -> displayed total mass of assimilant (grams)   */
    double  dspAssimT;      /* -> displayed temperature of assimilant          */
    int     dspAssimInc;    /* -> displayed total number of assim increments   */
    double  dspAssimLiqM;   /* -> displayed total mass of Liq assimilant (gms) */
    double  **assimComp;    /* -> current assimilant compos (moles of end.)    */
    int     *nAssimComp;    /* -> number of each type of phase in previous     */
    double  assimMass;      /* -> actual mass of assimilant added (grams)      */
    double  assimT;         /* -> assimilant temperature (C)                   */
    int     assimInc;       /* -> actual number of assim increments            */
    ThermoData assimTD;     /* -> thermodynamic properties of assimilant       */

    // Put in silminInputData ??
#ifdef ALPHAMELTS_UPDATE_SYSTEM
    char    *PTfile;        /* filename for P (or V) - T (or H or S) pairs     */
#endif
#ifdef PHMELTS_ADJUSTMENTS
    int     txtOutput;      /* current value of text output options            */
#else
    int     plotState;      /* current value of user configurable plot state   */
#endif
    double  *ySol;          /* array output from evaluateSaturationState       */
    double  *yLiq;          /* array output from evaluateSaturationState       */


    //int     nMajors;
    //int     nTraces;
#ifdef PHMELTS_ADJUSTMENTS
    double  minF;
    double  *fracLiquids;   /* fractionation coefficients for liquids          */
    double  minLiqPhi;
    double  *fracFluids;    /* fractionation coefficients for fluids           */
    double  minFluPhi;
    double  *fracSolids;    /* fractionation coefficients for solids           */
    double  maxF;
    double  fracOut;        /* default fractionation coefficient for all       */
#endif

} SilminState;

extern SilminState *silminState;

//#define SILMIN_STATE_CHANGE_NONE        0000000 /* octal mask   */
/*#define SILMIN_STATE_CHANGE_FATAL_ERROR 0000001
#define SILMIN_STATE_CHANGE_BULK        0000002
#define SILMIN_STATE_CHANGE_INC_SOLIDS  0000004
#define SILMIN_STATE_CHANGE_T           0000010
#define SILMIN_STATE_CHANGE_P           0000020
#define SILMIN_STATE_CHANGE_FO2PATH     0000040
#define SILMIN_STATE_CHANGE_FRAC_SOL    0000100
#define SILMIN_STATE_CHANGE_ISENTHALPIC 0000200
#define SILMIN_STATE_CHANGE_ISENTROPIC  0000400
#define SILMIN_STATE_CHANGE_ISOCHORIC   0001000
#define SILMIN_STATE_CHANGE_ASSIM_COMP  0002000
#define SILMIN_STATE_CHANGE_ASSIM_MASS  0004000
#define SILMIN_STATE_CHANGE_ASSIM_T     0010000
#define SILMIN_STATE_CHANGE_ASSIM_INC   0020000
#define SILMIN_STATE_CHANGE_PLOT        0040000
#define SILMIN_STATE_CHANGE_MUL_LIQUIDS 0100000
#define SILMIN_STATE_CHANGE_FRAC_LIQ    0200000
#define SILMIN_STATE_CHANGE_FRAC_FLU    0400000*/

typedef struct _silminInputData {
    char *name;        /* name of file which contains current input parameters */
    char *title;       /* title information                                    */
} SilminInputData;

extern SilminInputData silminInputData;

typedef struct _silminHistory *SilminHistoryPtr;
typedef struct _silminHistory {
    SilminState      *state;                   /* current state of the system */
    SilminHistoryPtr next;
} SilminHistory;

extern SilminHistory *silminHistory;

typedef struct _constraints {
    double *lambda;       /* Array of length nc+3 of Lagrange multipliers       */
    double *lambdaO2;     /* Lagrange multiplier for f O2 constraints           */
#ifdef PHMELTS_ADJUSTMENTS
    double lambdaH2O;     /* Lagrange multiplier for a H2O constraint       */
#endif
    double lambdaH;       /* Lagrange multiplier for isenthalpic constraint     */
    double lambdaS;       /* Lagrange multiplier for isentropic constraint      */
    double lambdaV;       /* Lagrange multiplier for isochoric constraint       */
    double **liquidDelta; /* Correction vector of liquid moles to maintain feas */
    double **solidDelta;  /* Correction vector of solid moles to maintain feasi */
    double T;             /* Consistent Temperature of the system               */
    double P;             /* Consistent Pressure of the system                  */
    double fo2;           /* Constraint log 10 fO2 of the system                */
#ifdef PHMELTS_ADJUSTMENTS
    double aH2O;          /* Constraint aH2O of the system           */
#endif
} Constraints;

extern Constraints *constraints;

extern int quad_tol_modifier;

/*
 *==============================================================================
 * Externally defined support functions:
 */

#define HESSIAN_TYPE_NORMAL 1
#define HESSIAN_TYPE_ONE    2

int         addOrDropLiquid(double *deltaBulkComp);
Constraints *allocConstraintsPointer(void);
SilminState *allocSilminStatePointer(void);
int         checkForCoexistingLiquids(void);
int         checkForCoexistingSolids(void);
int         checkStateAgainstInterface(void);
SilminState *copySilminStateStructure(SilminState *pOld, SilminState *pNew);
void        correctTforChangeInEnthalpy(void);
void        correctTforChangeInEntropy(void);
void        correctPforChangeInVolume(void);
#ifdef PHMELTS_ADJUSTMENTS
void        correctTforChangeInVolume(void);
void        correctXforChangeInBulkComp(void);
int         pdaNorm(void);
#endif
void        destroyConstraintsStructure(void *p);
void        destroySilminStateStructure(void *p);
int         evaluateSaturationState(double *rSol, double *rLiq);
double      formulaToMwStoich(char *formula, double *stoich);
int         getAffinityAndComposition(double t, double p, int index, int *zeroX,
                            double *muMinusMu0, double *affinity, double *indepVar);
void        getEqualityConstraints(int *conRows, int *conCols, double ***cMatrixPt,
                            double **hVectorPt, double **dVectorPt, double **yVectorPt);
double      getlog10fo2(double t, double p, int buffer);
double      getdlog10fo2dt(double t, double p, int buffer);
double      getdlog10fo2dp(double t, double p, int buffer);
double      getd2log10fo2dt2(double t, double p, int buffer);
double      getd2log10fo2dp2(double t, double p, int buffer);
int         getProjGradientAndHessian(int conRows, int conCols, double ***eMatrixPt,
                            double ***bMatrixPt, double **cMatrix, double *hVector, double *dVector,
                            double *yVector);
void        gibbs(double t, double p, char *name, ThermoRef *phase,
                            ThermoLiq *liquid, ThermoData *fusion, ThermoData *result);
void        InitComputeDataStruct(void);
void        intenToExtenGradient(double pMix, double *dpMix, int nr,  double *dp,
                            int na, double mTotal, double **drdm);
void        intenToExtenHessian(double pMix, double *dpMix, double **d2pMix,
                            int nr, double **d2p, int na, double mTotal, double **drdm,
                            double ***d2rdm2);
double      linearSearch(double lambda, int *notcomp);
int         spinodeTest(void);
int         subsolidusmuO2(int mask, double *muO2, double *dm, double *dt, double *dp,
                            double **d2m, double *d2mt, double *d2mp, double *d2t2, double *d2tp,
                            double *d2p2);
void        updateAssimilantPADB(char *member);
void        updateBulkADB(void);
void        updateSolidADB(double *rSol, double *rLiq);
void        updateCompADB(void);
void        updateTpPADB(int member);
void        updateUserGraphGW(void);

#ifdef PHMELTS_ADJUSTMENTS
double getLiquidProperties(SilminState *state, int ns, ThermoData *value,
            double *viscosity, int derivatives);
double getSolidProperties(SilminState *state, int index, int ns, ThermoData *value,
			  int derivatives);
double getSystemProperties(SilminState *state, int derivatives);
double getBulkSolidProperties(SilminState *state, int derivatives);
double getBulkLiquidProperties(SilminState *state, int derivatives);
#endif

#ifdef BATCH_VERSION

int liquidus(void);
int findWetLiquidus(void);
int putOutputDataToFile(char *);
int putSequenceDataToXmlFile(int);
int silmin(void);

#define ASSIM_PADB_INDEX_MASS        0  /* + npc + nc */
#define ASSIM_PADB_INDEX_T           1  /* + npc + nc */
#define ASSIM_PADB_INDEX_INCREMENT   2  /* + npc + nc */
#define ASSIM_PADB_INDEX_LIQUID_MASS 3  /* + npc + nc */
#define ASSIM_PADB_UNITS_VOLUME      0
#define ASSIM_PADB_UNITS_WEIGHT      1

#endif

#endif /* _Silmin_h */
