#ifndef FRAMEWORK_VERSION
//static const char compileDate[] = __DATE__;
//static const char compileTime[] = __TIME__;
#endif /* NOT FRAMEWORK_VERSION */
const char *liquid_ver(void) { return "$Id: liquid.c,v 1.42 2009/05/14 04:24:00 ghiorso Exp $"; }

/*
MELTS Source Code: RCS $Log: liquid.c,v $
MELTS Source Code: RCS Revision 1.30  2008/05/24 21:32:43  ghiorso
MELTS Source Code: RCS Added NaAlSiO4, KAlSiO4, and Ca(1/2)AlSiO4 species.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.29  2008/05/23 17:07:31  ghiorso
MELTS Source Code: RCS Added more hydroxyl species.
MELTS Source Code: RCS Improved preclb interface for larger number of species.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.28  2008/05/20 17:32:21  ghiorso
MELTS Source Code: RCS Added FeSiO3 and FeAlO2.5 species to balance Mg-equivalents.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.27  2008/05/03 18:16:16  ghiorso
MELTS Source Code: RCS Revised Fe2O3 and FeO1.3 EOS coefficients.
MELTS Source Code: RCS Corrected Ferric/ferrous calculation in liquid.c
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.26  2008/05/02 19:03:52  ghiorso
MELTS Source Code: RCS Revised liquid speciation model.
MELTS Source Code: RCS Created new test routine for homogeneous equilibrium fO2 at P.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.25  2008/01/06 22:35:59  ghiorso
MELTS Source Code: RCS Updated param and liq data structures to modify EOS parameters for FeO1.3
MELTS Source Code: RCS and to compute f O2 and ferric/ferrous from KC '89 + EOS integral
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.24  2007/12/22 22:43:30  ghiorso
MELTS Source Code: RCS Fixed error in BM integration in Gibbs.c
MELTS Source Code: RCS Updated param_struct_data.h file for AGU 2007 xMELTS parameters
MELTS Source Code: RCS Added support for status file production in MELTS-batch (XML)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.23  2007/11/28 01:05:37  ghiorso
MELTS Source Code: RCS Reverted back to version:
MELTS Source Code: RCS 1.21 liquid.c
MELTS Source Code: RCS 1.13 liq_struct_data.h
MELTS Source Code: RCS 1.33 param_struct_data.h
MELTS Source Code: RCS i.e. MgSiO3 and CaSiO3 species only.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.21  2007/11/23 17:47:50  ghiorso
MELTS Source Code: RCS Added the species CaSiO3.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.20  2007/11/22 04:08:13  ghiorso
MELTS Source Code: RCS Corrected infinite loop error in order() in albite.c
MELTS Source Code: RCS Removed arbitrary volume corrections in sol_struct_data.h
MELTS Source Code: RCS Turned on non-quadrilateral cpx endmembers for regression.
MELTS Source Code: RCS Added MgSiO3 species to liquid model.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.19  2007/10/18 00:33:37  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.18  2007/10/18 00:01:42  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.17  2007/10/17 16:38:20  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.16  2007/10/15 17:43:40  ghiorso
MELTS Source Code: RCS Improved convergence criteria, roundoff errors and numerical stability of
MELTS Source Code: RCS the ordering routines in liquid.c
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.15  2007/10/09 01:08:59  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.14  2007/10/06 16:47:56  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.13  2007/10/04 17:39:21  ghiorso
MELTS Source Code: RCS Relaxed convergence criteria in liquid ordering routine.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.12  2007/10/04 00:46:05  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.11  2007/10/03 21:33:48  ghiorso
MELTS Source Code: RCS Updated liquid eos thermodynamics.
MELTS Source Code: RCS Added regression of ferric/ferrous parameters from data file.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.10  2007/08/30 18:26:10  ghiorso
MELTS Source Code: RCS Revised high-pressure log fO2 calculations (buffers and liquid)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.9  2007/07/09 21:33:41  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.8  2007/06/16 01:01:55  ghiorso
MELTS Source Code: RCS Revised EOS regression to have K', K'', and K''' as parameters.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.7  2007/06/13 15:12:47  ghiorso
MELTS Source Code: RCS Revised strategy for EOS parameter fitting.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.6  2007/06/09 20:30:04  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.5  2007/06/08 17:25:42  ghiorso
MELTS Source Code: RCS Added code to allow regression of Ghiorso EOS parameters
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.4  2007/02/13 21:48:29  ghiorso
MELTS Source Code: RCS Modifications to read XML database files for LEPER calibration.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.3  2006/10/20 00:59:22  ghiorso
MELTS Source Code: RCS (1) Made initial modifications for thread safe code.
MELTS Source Code: RCS (2) Added support for XML I/O in batch mode
MELTS Source Code: RCS (3) Added support for Melts-batch listener for eventual integration into VIGMCS
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2006/08/17 16:47:18  ghiorso
MELTS Source Code: RCS Made modifications to protect strings.  These modifications allow removal
MELTS Source Code: RCS of the flag -fwritable-strings during gcc compilation.  This brings the
MELTS Source Code: RCS code up to gcc 4.x standards.
MELTS Source Code: RCS
MELTS Source Code: RCS Other minor rearrangements and cleanup.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2006/08/15 16:57:36  ghiorso
MELTS Source Code: RCS xMELTS gcc 3.x sources
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.16  2005/02/12 18:42:34  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.15  2005/02/11 03:27:10  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.14  2005/01/25 03:25:03  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.13  2005/01/24 03:38:04  cvsaccount
MELTS Source Code: RCS
MELTS Source Code: RCS Added new files and modifications to perform builds for MgO-SiO2 system
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.12  2005/01/23 20:13:27  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.11  2005/01/21 18:18:22  cvsaccount
MELTS Source Code: RCS
MELTS Source Code: RCS Added data structures and code to implement coordination number transformations
MELTS Source Code: RCS in the liquid phase EOS model.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.10  2005/01/05 17:27:26  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.9  2004/12/11 22:19:43  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.8  2004/10/03 22:40:08  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.7  2004/09/27 18:24:41  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.6  2004/09/25 19:01:50  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.5  2004/09/24 18:26:35  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.4  2004/09/24 02:31:21  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.3  2004/09/23 21:18:06  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2004/09/23 20:07:04  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2004/01/02 19:21:49  cvsaccount
MELTS Source Code: RCS CTserver University of Chicago
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.11  2003/09/30 17:36:38  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.10  2003/09/27 15:35:22  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.9  2002/08/03 00:21:40  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.8  2002/07/31 20:42:27  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.7  2002/07/29 20:44:03  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.6  2002/07/26 23:55:04  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.5  2002/07/26 01:01:37  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.4  2002/04/10 20:39:27  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.3  2002/04/10 00:42:21  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2002/04/06 00:51:39  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:58:12  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Routines to compute liquid solution properties
**      (file: LIQUID.C)
**
**--
*/

#ifdef DEBUG
#undef DEBUG
#endif

#define DO_NOT_USE_GHIORSO_KRESS_MODEL
#define USE_KRESS_CARMICHAEL_FO2

#include "melts_gsl.h"

#include "silmin.h"
#include "mthread.h"
#ifdef USESEH
#include <windows.h>
void raise_sigabrt(DWORD dwType);
extern int doInterrupt;
#else
#include <signal.h>
#endif

#include "param_struct_data.h"

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))

// from olivine.c
#define USE_SVD 0
#define MAX_ITER 200  /* Maximum number of iterations allowed in order */

enum {
    ERR_NONE,
    ERR_DGDR_GMAP_1,
    ERR_DGDR_GMAP_2,
    ERR_DGDT_GMAP,
    ERR_D2GDRDT_GMAP_1,
    ERR_D2GDRDT_GMAP_2,
    ERR_D3GDRDT2_GMAP_1,
    ERR_D3GDRDT2_GMAP_2,
    ERR_D2GDT2_GMAP,
    ERR_D3GDT3_GMAP,
    ERR_D2GDR2_GMAP,
    ERR_D3GDR2DT_GMAP,
    ERR_D3GDR3_GMAP,
    ERR_AB_ZERO,
    ERR_A_ZERO,
    ERR_B_ZERO,
    ERR_SUM_ZERO,
    ERR_SUM_LT_ZERO,
    ERR_SUM_GT_ZERO
};

#define printERR(name, num, description, value) \
    { fprintf(stderr, "*-->Exception in %s (%d) (file liquid.c). Bad argument to log(%s = %13.6e)\n", name, num, description, value); liqERRstate = num; return 0.0;}

/* This global variable is only accessed from preclb.  The construction here
   is not thread safe.                                                       */

int liqERRstate = ERR_NONE;

/*
 *===========================================================================
 * Backward compatible liquid functions defined in liquid_v34.c
*/

void conLiq_v34  (int inpMask, int outMask, double t, double p, double *o, double *m, double *r, double *x, double **dm, double ***dm2, double *logfo2);
int  testLiq_v34 (int mask, double t, double p, int na, int nr, char **names, char **formulas, double *r, double *m);
void dispLiq_v34 (int mask, double t, double p, double *x, char **formula);
void actLiq_v34  (int mask, double t, double p, double *x, double *a,     double *mu,   double **dx);
void gmixLiq_v34 (int mask, double t, double p, double *x, double *gmix,  double *dx,   double **dx2);
void hmixLiq_v34 (int mask, double t, double p, double *x, double *hmix);
void smixLiq_v34 (int mask, double t, double p, double *x, double *smix,  double *dx,   double **dx2);
void cpmixLiq_v34(int mask, double t, double p, double *x, double *cpmix, double *dt,   double *dx);
void vmixLiq_v34 (int mask, double t, double p, double *x, double *vmix,  double *dx,   double **dx2, double *dt, double *dp,   double *dt2,  double *dtdp,
                                                           double *dp2,   double *dxdt, double *dxdp);
void muO2Liq_v34 (int mask, double t, double p, double *m, double *muO2,  double *dm,   double *dt,   double *dp, double **d2m, double *d2mt, double *d2mp,
                                                           double *d2t2,  double *d2tp, double *d2p2);
void visLiq_v34  (int mask, double t, double p, double *x, double *viscosity);

void conLiq_CO2  (int inpMask, int outMask, double t, double p, double *o, double *m, double *r, double *x, double **dm, double ***dm2, double *logfo2);
int  testLiq_CO2 (int mask, double t, double p, int na, int nr, char **names, char **formulas, double *r, double *m);
void dispLiq_CO2 (int mask, double t, double p, double *x, char **formula);
void actLiq_CO2  (int mask, double t, double p, double *x, double *a,     double *mu,   double **dx);
void gmixLiq_CO2 (int mask, double t, double p, double *x, double *gmix,  double *dx,   double **dx2);
void hmixLiq_CO2 (int mask, double t, double p, double *x, double *hmix);
void smixLiq_CO2 (int mask, double t, double p, double *x, double *smix,  double *dx,   double **dx2);
void cpmixLiq_CO2(int mask, double t, double p, double *x, double *cpmix, double *dt,   double *dx);
void vmixLiq_CO2 (int mask, double t, double p, double *x, double *vmix,  double *dx,   double **dx2, double *dt, double *dp,   double *dt2,  double *dtdp,
                                    double *dp2,   double *dxdt, double *dxdp);
void muO2Liq_CO2 (int mask, double t, double p, double *m, double *muO2,  double *dm,   double *dt,   double *dp, double **d2m, double *d2mt, double *d2mp,
                                    double *d2t2,  double *d2tp, double *d2p2);
void visLiq_CO2  (int mask, double t, double p, double *x, double *viscosity);

void conLiq_CO2_H2O  (int inpMask, int outMask, double t, double p, double *o, double *m, double *r, double *x, double **dm, double ***dm2, double *logfo2);
int  testLiq_CO2_H2O (int mask, double t, double p, int na, int nr, char **names, char **formulas, double *r, double *m);
void dispLiq_CO2_H2O (int mask, double t, double p, double *x, char **formula);
void actLiq_CO2_H2O  (int mask, double t, double p, double *x, double *a,     double *mu,   double **dx);
void gmixLiq_CO2_H2O (int mask, double t, double p, double *x, double *gmix,  double *dx,   double **dx2);
void hmixLiq_CO2_H2O (int mask, double t, double p, double *x, double *hmix);
void smixLiq_CO2_H2O (int mask, double t, double p, double *x, double *smix,  double *dx,   double **dx2);
void cpmixLiq_CO2_H2O(int mask, double t, double p, double *x, double *cpmix, double *dt,   double *dx);
void vmixLiq_CO2_H2O (int mask, double t, double p, double *x, double *vmix,  double *dx,   double **dx2, double *dt, double *dp,   double *dt2,  double *dtdp,
                                    double *dp2,   double *dxdt, double *dxdp);
void muO2Liq_CO2_H2O (int mask, double t, double p, double *m, double *muO2,  double *dm,   double *dt,   double *dp, double **d2m, double *d2mt, double *d2mp,
                                    double *d2t2,  double *d2tp, double *d2p2);
void visLiq_CO2_H2O  (int mask, double t, double p, double *x, double *viscosity);

/*#ifdef PHMELTS_ADJUSTMENTS
void conLiq_Alt  (int inpMask, int outMask, double t, double p, double *o, double *m, double *r, double *x, double **dm, double ***dm2, double *logfo2);
void muO2Liq_Alt (int mask, double t, double p, double *m, double *muO2,  double *dm,   double *dt,   double *dp, double **d2m, double *d2mt, double *d2mp,
                                                    double *d2t2,  double *d2tp, double *d2p2);
        #endif*/

/*
 *=============================================================================
 * Private functions and globals:
*/

#define NA 19                    /* Number of liquid components                   */
#define NS  1                    /* Number of ordering parameters for liquid species      */
#define NY  0                    /* Number of ordering parameters for coordination states */

static const int iOxAl2O3      =  2; /* Index of Al2O3 in bulksystem[] structure array    */
static const int iOxFe2O3      =  3; /* Index of Fe2O3 in bulksystem[] structure array    */
static const int iOxFeO        =  5; /* Index of FeO in bulksystem[] structure array      */
//static const int iOxCaO        = 10; /* Index of CaO in bulksystem[] structure array      */
//static const int iOxNa2O       = 11; /* Index of Na2O in bulksystem[] structure array     */
//static const int iOxK2O        = 12; /* Index of K2O in bulksystem[] structure array      */
static const int iOxFeO1_3     = 19; /* Index of FeO1.3 in bulksystem[] structure array   */
//static const int iCmpAl2O3     =  1; /* Index of Al2O3 in r[] array                       */
//static const int iCmpFe2SiO5   = -1; /* Index of Fe2SiO5 in r[] array                     */
//static const int iCmpFe2SiO4   =  4; /* Index of Fe2SiO4 in r[] array                     */
static const int iCmpFe2SiO4_6 = -1; /* Index of Fe2SiO4.6 in s[] array                   */
//static const int iCmpFe2AlO4_5 = -1; /* Index of Fe2AlO4.5 in s[] array                   */
//static const int iCmpFe2AlO3_5 = -1; /* Index of Fe2AlO3.5 in s[] array                   */
static const int iCmpFe2AlO4_1 = -1; /* Index of Fe2AlO4.1 in s[] array                   */

#define NT (NS+NY)               /* Number of ordering parameters                         */
#define NR (NA-1)                /* Number of independent mole fraction variables         */
#define NW ((NA+NS)*(NA+NS-1)/2) /* Number of regular solution interaction parameters     */
#define NV (NR+NS)               /* Number of independent variables in the model          */
#define NP (NA+NS+NW)            /* Number of model parameters                    */

static int nH2O, nCO2;

static double gT, dgdrT[NR], d2gdr2T[NR][NR], d3gdr3T[NR][NR][NR];

static void ternaryH2OCO2terms(int mask, double *r) {
    static const double WHCX[NA] = {
                                                                        0.0*3.0, //  0    SiO2
                                                                        0.0*3.0, //  1  0 TiO2
                                                                        0.0*3.0, //  2  1 Al2O3
                                                                        0.0*3.0, //  3  2 FE2O3
                                                                        0.0*3.0, //  4  3 MgCr2O3
                                                                        0.0*3.0, //  5  4 Fe2SiO4
                                                                        0.0*3.0, //  6  5 MnSi1/2O2
                                                                        0.0*3.0, //  7  6 Mg2SiO4
                                                                        0.0*3.0, //  8  7 NiSi1/2O2
                                                                        0.0*3.0, //  9  8 CoSi1/2O2
                                                                        0.0*3.0, // 10  9 CaSiO3
                                                                        0.0*3.0, // 11 10 Na2SiO3
                                                                        0.0*3.0, // 12 11 KAlSiO4
                                                                        0.0*3.0, // 13 12 Ca3(PO4)2
                                                                        0.0*3.0, // 14 13 CO2
                                                                        0.0*3.0, // 15 14 SO3
                                                                        0.0*3.0, // 16 15 Cl2O
                                                                        0.0*3.0, // 17 16 F2O
                                                                        0.0*3.0  // 18 17 H2O
                             };
    static double x[NA];
    int i;

    for (i=0, x[0]=1.0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }

    if (mask & FIRST) {
        for (i=0, gT=0.0; i<NA; i++) gT += WHCX[i]*x[nH2O]*x[nCO2]*x[i];
    }

    if (mask & SECOND) {
        int j;
        for (i=0; i<NR; i++) {
            dgdrT[i] = 0.0;
            switch (i) {
                case 13:
                {
                    dgdrT[nCO2-1] = 0.0;
                    for (j=1; j<NA; j++) dgdrT[nCO2-1] += x[nH2O]*(WHCX[j]-WHCX[0])*x[j];
                    dgdrT[nCO2-1] += WHCX[0]*x[nH2O] + (WHCX[nCO2]-WHCX[0])*x[nH2O]*x[nCO2];
                    break;
                }
                case 17:
                {
                    dgdrT[nH2O-1] = 0.0;
                    for (j=1; j<NA; j++) dgdrT[nH2O-1] += x[nCO2]*(WHCX[j]-WHCX[0])*x[j];
                    dgdrT[nH2O-1] += WHCX[0]*x[nCO2] + (WHCX[nH2O]-WHCX[0])*x[nH2O]*x[nCO2];
                    break;
                }
                default:
                {
                    dgdrT[i] = (WHCX[i+1]-WHCX[0])*x[nH2O]*x[nCO2];
                }
            }
        }
    }

    if (mask & THIRD) {
        int j;
        for (i=0; i<NR; i++) for (j=0; j<NR; j++) d2gdr2T[i][j] = 0.0;
        d2gdr2T[nH2O-1][nH2O-1] = 2.0*(WHCX[nH2O]-WHCX[0])*x[nCO2];
        d2gdr2T[nCO2-1][nCO2-1] = 2.0*(WHCX[nCO2]-WHCX[0])*x[nH2O];

        d2gdr2T[nH2O-1][nCO2-1] = WHCX[0]*(x[0]-x[nCO2]-x[nH2O]);
            for (i=1; i<NA; i++) d2gdr2T[nH2O-1][nCO2-1] += WHCX[i]*x[i];
            d2gdr2T[nH2O-1][nCO2-1] += WHCX[nCO2]*x[nCO2] + WHCX[nH2O]*x[nH2O];
            d2gdr2T[nCO2-1][nH2O-1] = d2gdr2T[nH2O-1][nCO2-1];

        for (i=1; i<NA; i++) {
            if (i != nCO2 && i != nH2O) {
                d2gdr2T[i-1][nH2O-1] = (WHCX[i]-WHCX[0])*x[nCO2];
                d2gdr2T[nH2O-1][i-1] = d2gdr2T[i-1][nH2O-1];

                d2gdr2T[i-1][nCO2-1] = (WHCX[i]-WHCX[0])*x[nH2O];
                d2gdr2T[nCO2-1][i-1] = d2gdr2T[i-1][nCO2-1];
            }
        }
    }

    if (mask & FOURTH) {
        int j, k;
        for (i=0; i<NR; i++) for (j=0; j<NR; j++) for (k=0; k<NR; k++) d3gdr3T[i][j][k] = 0.0;
        d3gdr3T[nH2O-1][nH2O-1][nCO2-1] = 2.0*(WHCX[nH2O]-WHCX[0]);
            d3gdr3T[nH2O-1][nCO2-1][nH2O-1] = d3gdr3T[nH2O-1][nH2O-1][nCO2-1];
            d3gdr3T[nCO2-1][nH2O-1][nH2O-1] = d3gdr3T[nH2O-1][nH2O-1][nCO2-1];

        d3gdr3T[nH2O-1][nCO2-1][nCO2-1] = 2.0*(WHCX[nCO2]-WHCX[0]);
            d3gdr3T[nCO2-1][nH2O-1][nCO2-1] = d3gdr3T[nH2O-1][nCO2-1][nCO2-1];
            d3gdr3T[nCO2-1][nCO2-1][nH2O-1] = d3gdr3T[nH2O-1][nCO2-1][nCO2-1];

        for (i=1; i<NA; i++) {
            if (i != nCO2 && i != nH2O) {
                d3gdr3T[i-1][nH2O-1][nCO2-1] = WHCX[i] - WHCX[0];
                d3gdr3T[i-1][nCO2-1][nH2O-1] = d3gdr3T[i-1][nH2O-1][nCO2-1];
                d3gdr3T[nH2O-1][i-1][nCO2-1] = d3gdr3T[i-1][nH2O-1][nCO2-1];
                d3gdr3T[nCO2-1][i-1][nH2O-1] = d3gdr3T[i-1][nH2O-1][nCO2-1];
                d3gdr3T[nH2O-1][nCO2-1][i-1] = d3gdr3T[i-1][nH2O-1][nCO2-1];
                d3gdr3T[nCO2-1][nH2O-1][i-1] = d3gdr3T[i-1][nH2O-1][nCO2-1];
            }
        }
    }

}

/* The statics from here to ... */

static int convergedInOrder;

static double *xSpecies;         /* Mole fractions of endmember species                                                  */
static double **dxSpeciesdr;     /* d(xSpecies)/dr                                                                       */
static double **dxSpeciesds;     /* d(xSpecies)/ds                                                                       */
static double ***d2xSpeciesdrds; /* d2(xSpecies)/drds                                                                    */
static double nSpecies;          /* Total moles of all solution species relative to 1 mole of basis species              */
static double *dnSpeciesds;      /* d(nSpecies)/ds                                                                       */
static double **d2nSpeciesds2;   /* d2(nSpecies)/ds2                                                                     */
static double ***d3nSpeciesds3;  /* d3(nSpecies)/ds3                                                                     */

static MTHREAD_MUTEX_T global_data_mutex = MTHREAD_MUTEX_INITIALIZER;

/* ... here are dealt with by creating a mutex lock on the code that requires access to these quantities */

/***********************************************/
/* Statics for class initialization structures */
/***********************************************/

static MTHREAD_ONCE_T initThreadBlock = MTHREAD_ONCE_INIT;

static void initializeLiquid(void);

static void threadInit(void) {
    initializeLiquid();
}

/* The statics from here ... */

static int NE;   /* Number of liquid endmembers (species) */

static double Gconst,       *gr,       *gs,       **grr,        **grs,       **gss;
static double Hconst,       *hr,       *hs,       **hrr,        **hrs,       **hss;
static double Sconst,       *sr,       *ss,       **srr,        **srs,       **sss;
static double Vconst,       *vr,       *vs,       **vrr,        **vrs,       **vss;
static double CPconst,      *cpr,      *cps;
static double DCPDTconst,   *dcpdtr,   *dcpdts;
static double DVDTconst,    *dvdtr,    *dvdts;
static double DVDPconst,    *dvdpr,    *dvdps;
static double D2VDT2const,  *d2vdt2r,  *d2vdt2s;
static double D2VDTDPconst, *d2vdtdpr, *d2vdtdps;
static double D2VDP2const,  *d2vdp2r,  *d2vdp2s;

static double **taylorCoeff;  /* Taylor Expansion coefficients: [endmember species | W(i,j)][g0, gr, gs, grr, grs, gss] */
static double **rsEndmembers; /* r and s coefficients for species endmembers                                            */

/* ... to here are set in initializeLiquid() */

/* This is adapted from olivine.c */
static MTHREAD_ONCE_T initThreadOBlock = MTHREAD_ONCE_INIT;

static MTHREAD_KEY_T tOldKey;
static MTHREAD_KEY_T pOldKey;
static MTHREAD_KEY_T rOldKey;
static MTHREAD_KEY_T sOldKey;
static MTHREAD_KEY_T ptToD2gds2Key;
static MTHREAD_KEY_T d2gds2Key;
static MTHREAD_KEY_T ptToVKey;
static MTHREAD_KEY_T ptToSKey;
static MTHREAD_KEY_T indexD2gds2Key;

static void freeNSvector(void *NSarray) {
    gsl_vector_free((gsl_vector *) NSarray);
}

static void freeNSmatrix(void *NSarray) {
    gsl_matrix_free((gsl_matrix *) NSarray);
}

static void freeIndexD2gds2(void *indexD2gds2) {
    gsl_permutation_free((gsl_permutation *) indexD2gds2);
}

static void threadOInit(void) {
    MTHREAD_KEY_CREATE(&tOldKey,       free);
    MTHREAD_KEY_CREATE(&pOldKey,       free);
    MTHREAD_KEY_CREATE(&rOldKey,       freeNSvector);
    MTHREAD_KEY_CREATE(&sOldKey,       freeNSvector);
    MTHREAD_KEY_CREATE(&ptToD2gds2Key, freeNSmatrix);
    MTHREAD_KEY_CREATE(&d2gds2Key,     free);
    MTHREAD_KEY_CREATE(&ptToVKey,      freeNSmatrix);
    MTHREAD_KEY_CREATE(&ptToSKey,      freeNSvector);
    MTHREAD_KEY_CREATE(&indexD2gds2Key, freeIndexD2gds2);
}

static double getTOld() {
    double *tOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    tOldPt = (double *) MTHREAD_GETSPECIFIC(tOldKey);
    if (tOldPt == NULL) {
        tOldPt  = (double *) malloc(sizeof(double));
        *tOldPt = -9999.0;
        MTHREAD_SETSPECIFIC(tOldKey, (void *) tOldPt);
    }
    return *tOldPt;
}

static void setTOld(double tOld) {
    double *tOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    tOldPt = (double *) MTHREAD_GETSPECIFIC(tOldKey);
    if (tOldPt == NULL) {
        tOldPt  = (double *) malloc(sizeof(double));
        *tOldPt = -9999.0;
        MTHREAD_SETSPECIFIC(tOldKey, (void *) tOldPt);
    }
    *tOldPt = tOld;
}

static double getPOld() {
    double *pOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    pOldPt = (double *) MTHREAD_GETSPECIFIC(pOldKey);
    if (pOldPt == NULL) {
        pOldPt  = (double *) malloc(sizeof(double));
        *pOldPt = -9999.0;
        MTHREAD_SETSPECIFIC(pOldKey, (void *) pOldPt);
    }
    return *pOldPt;
}

static void setPOld(double pOld) {
    double *pOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    pOldPt = (double *) MTHREAD_GETSPECIFIC(pOldKey);
    if (pOldPt == NULL) {
        pOldPt  = (double *) malloc(sizeof(double));
        *pOldPt = -9999.0;
        MTHREAD_SETSPECIFIC(pOldKey, (void *) pOldPt);
    }
    *pOldPt = pOld;
}

static double *getROld() {
    gsl_vector *rOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    rOldPt = (gsl_vector *) MTHREAD_GETSPECIFIC(rOldKey);
    if (rOldPt == NULL) {
        rOldPt = gsl_vector_alloc((size_t) NR);
        gsl_vector_set_all(rOldPt, -9999.0);
        MTHREAD_SETSPECIFIC(rOldKey, (void *) rOldPt);
    }
    return rOldPt->data;
}

static double *getSOld() {
    gsl_vector *sOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    sOldPt = (gsl_vector *) MTHREAD_GETSPECIFIC(sOldKey);
    if (sOldPt == NULL) {
        sOldPt = gsl_vector_alloc((size_t) NT);
        gsl_vector_set_all(sOldPt, 0.0);
        MTHREAD_SETSPECIFIC(sOldKey, (void *) sOldPt);
    }
    return sOldPt->data;
}

static gsl_matrix *getPtToD2gds2() {
    gsl_matrix *ptToD2gds2Pt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    ptToD2gds2Pt = (gsl_matrix *) MTHREAD_GETSPECIFIC(ptToD2gds2Key);
    if (ptToD2gds2Pt == NULL) {
        ptToD2gds2Pt  = gsl_matrix_alloc((size_t) NT, (size_t) NT);
        gsl_matrix_set_zero(ptToD2gds2Pt);
        MTHREAD_SETSPECIFIC(ptToD2gds2Key, (void *) ptToD2gds2Pt);
    }
    return ptToD2gds2Pt;
}

static double **getD2gds2() {
    double **d2gds2Pt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    d2gds2Pt = (double **) MTHREAD_GETSPECIFIC(d2gds2Key);
    if (d2gds2Pt == NULL) {
        int i;
        gsl_matrix *ptToD2gds2Pt = getPtToD2gds2();
        d2gds2Pt  = (double **) malloc((size_t) NT*sizeof(double *));
        for (i=0; i<NT; i++) d2gds2Pt[i] = &(ptToD2gds2Pt->data[i*NT]);
        MTHREAD_SETSPECIFIC(d2gds2Key, (void *) d2gds2Pt);
    }
    return d2gds2Pt;
}

static gsl_matrix *getPtToV() {
    gsl_matrix *ptToVPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    ptToVPt = (gsl_matrix *) MTHREAD_GETSPECIFIC(ptToVKey);
    if (ptToVPt == NULL) {
        ptToVPt  = gsl_matrix_alloc((size_t) NT, (size_t) NT);
        gsl_matrix_set_zero(ptToVPt);
        MTHREAD_SETSPECIFIC(ptToVKey, (void *) ptToVPt);
    }
    return ptToVPt;
}

static gsl_vector *getPtToS() {
    gsl_vector *ptToSPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    ptToSPt = (gsl_vector *) MTHREAD_GETSPECIFIC(ptToSKey);
    if (ptToSPt == NULL) {
        ptToSPt  = gsl_vector_alloc((size_t) NT);
        gsl_vector_set_zero(ptToSPt);
        MTHREAD_SETSPECIFIC(ptToSKey, (void *) ptToSPt);
    }
    return ptToSPt;
}

static gsl_permutation *getIndexD2gds2() {
    gsl_permutation *indexD2gds2Pt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    indexD2gds2Pt = (gsl_permutation *) MTHREAD_GETSPECIFIC(indexD2gds2Key);
    if (indexD2gds2Pt == NULL) {
        indexD2gds2Pt = gsl_permutation_alloc((size_t) NT);
        gsl_permutation_init(indexD2gds2Pt);
        MTHREAD_SETSPECIFIC(indexD2gds2Key, (void *) indexD2gds2Pt);
    }
    return indexD2gds2Pt;
}

static double **matrix_alloc(int n1, int n2) {
    int i;
    double *m0 = (double *) malloc((size_t) n1*n2*sizeof(double));
    double **m = (double **) malloc((size_t) n1*sizeof(double *));
    for (i=0; i<n1; i++) m[i] = &m0[i*n2];
    return m;
}

static void matrix_free(double **m, int n1, int n2) {
    free(m[0]);
    free(m);
}

static double *vector_alloc(int n) {
    double *v = (double *) malloc((size_t) n*sizeof(double));
    return v;
}

static void vector_free(double *v, int n) {
    free(v);
}

static void initializeTaylor(gsl_vector *array) {
    gsl_vector_set_zero(array);
}

static void loadTaylor(gsl_vector *vTaylor, double coeff, double *p) {
    int i, j, n;
    double array[NP];
    gsl_vector_view temp;

    array[0] = coeff;
    for (i=0, n=1; i<NV; i++, n++) array[n] = (p[i] != 0.0) ? coeff*p[i] : 0.0;
    for (i=0; i<NV; i++) for (j=i; j<NV; j++, n++)
        array[n] =  (p[i] != 0.0 && p[j] != 0.0) ? coeff*p[i]*p[j] : 0.0;

    temp = gsl_vector_view_array(array, (size_t) NP);
    gsl_vector_add(vTaylor, &temp.vector);

}

static void initializeLiquid(void) {
    int i, j, k, m, n;
    static int contentsOK = FALSE;
    double *temp, *coeff;
    gsl_matrix *mTaylor, *bTaylor;

#ifdef DEBUG
    printf("Entering function initializeLiquid...\n");
#endif

    /* Initialize global constants */
    NE = nls; /* Number of liquid endmembers (species) */
    if ( (NE == 0) || (NA != nlc) || (NS != (nls-nlc)) ) {
        if (NE == 0)         printf("Error in initializeLiquid() [liquid.c].  Number of species is zero.\n");
        if (NA != nlc)       printf("Error in initializeLiquid() [liquid.c].  Constant NA is not equal to nlc (%d).\n", nlc);
        if (NS != (nls-nlc)) printf("Error in initializeLiquid() [liquid.c].  Constant NS is not equal to nls-nlc (%d).\n", nls-nlc);
        exit(0);
    }

    nH2O = -1;
    for (i=0; i<NE; i++) if ((strcmp(liquid[i].label, "H2O") == 0) || (strcmp(liquid[i].label, "h2o") == 0)) { nH2O = i; break; }
    nCO2 = -1;
    for (i=0; i<NE; i++) if ((strcmp(liquid[i].label, "CO2") == 0) || (strcmp(liquid[i].label, "co2") == 0)) { nCO2 = i; break; }

    /* Static global storage for endmember species mole fractions */
    xSpecies       = vector_alloc(NE);
    dxSpeciesdr    = matrix_alloc(NE, NR);
    dxSpeciesds    = matrix_alloc(NE, NS);
    d2xSpeciesdrds = (double ***) malloc((size_t) NE*sizeof(double **));
                   for (i=0; i<NE; i++) d2xSpeciesdrds[i] = matrix_alloc(NR, NS);
    dnSpeciesds    = vector_alloc(NS);
    d2nSpeciesds2  = matrix_alloc(NS, NS);
    d3nSpeciesds3  = (double ***) malloc((size_t) NS*sizeof(double **));
                   for (i=0; i<NS; i++) d3nSpeciesds3[i] = matrix_alloc(NS, NS);

    /* Allocate static storage for the derivatives of g */

    gr = vector_alloc(NR);
    gs = vector_alloc(NS);
    grr = matrix_alloc(NR, NR);
    grs = matrix_alloc(NR, NS);
    gss = matrix_alloc(NS, NS);
    hr = vector_alloc(NR);
    hs = vector_alloc(NS);
    hrr = matrix_alloc(NR, NR);
    hrs = matrix_alloc(NR, NS);
    hss = matrix_alloc(NS, NS);
    sr = vector_alloc(NR);
    ss = vector_alloc(NS);
    srr = matrix_alloc(NR, NR);
    srs = matrix_alloc(NR, NS);
    sss = matrix_alloc(NS, NS);
    vr = vector_alloc(NR);
    vs = vector_alloc(NS);
    vrr = matrix_alloc(NR, NR);
    vrs = matrix_alloc(NR, NS);
    vss = matrix_alloc(NS, NS);
    cpr = vector_alloc(NR);
    cps = vector_alloc(NS);
    dcpdtr = vector_alloc(NR);
    dcpdts = vector_alloc(NS);
    dvdtr = vector_alloc(NR);
    dvdts = vector_alloc(NS);
    dvdpr = vector_alloc(NR);
    dvdps = vector_alloc(NS);
    d2vdt2r = vector_alloc(NR);
    d2vdt2s = vector_alloc(NS);
    d2vdtdpr = vector_alloc(NR);
    d2vdtdps = vector_alloc(NS);
    d2vdp2r = vector_alloc(NR);
    d2vdp2s = vector_alloc(NS);

#ifdef DEBUG
    printf("...Constructing r,s --> Endmember species matrix.\n");
#endif
    /* r and s coefficients for endmember species */
    rsEndmembers = matrix_alloc (NE, NV);
    coeff        = vector_alloc (NE);
    for (i=0; i<NE; i++) for (j=0; j<NV; j++) rsEndmembers[i][j] = 0.0;
    for (i=0; i<NA; i++) coeff[i] = 1.0;
    /* --> basis species     */
    for (i=1; i<NA; i++) rsEndmembers[i][i-1] = 1.0;
    /* --> dependent species */
    for (i=NA; i<NE; i++) {
        for (k=0, coeff[i]=0.0; k<nc; k++) coeff[i] += (bulkSystem[k].oxToLiq)[0]*(liquid[i].liqToOx)[k];
        for (j=1; j<NA; j++) for (k=0, rsEndmembers[i][j-1]=0.0; k<nc; k++)
            rsEndmembers[i][j-1] += (bulkSystem[k].oxToLiq)[j]*(liquid[i].liqToOx)[k];
        for (j=0; j<NR; j++) coeff[i] += rsEndmembers[i][j];
        if (coeff[i] != 0.0) for (j=0; j<NR; j++) rsEndmembers[i][j] /= coeff[i];
        rsEndmembers[i][i-NA+NR] = 1.0;  /* ordering parameter */
    }
#ifdef DEBUG
    printf("...rsEndmembers matrix:\n");
    for (i=0; i<NE; i++) {
        printf("coeff[%-15.15s] = %10.4f\n", liquid[i].label, coeff[i]);
        for (j=0; j<NV; j++) printf("%4.1f ", rsEndmembers[i][j]);
        printf("\n");
    }
#endif

#ifdef DEBUG
    printf("...Zeroing dxdr, dxds, and d2xdrds matrices.\n");
#endif
    for (i=0; i<NE; i++) {
        for (j=0; j<NR; j++) {
            dxSpeciesdr[i][j] = 0.0;
            for (k=0; k<NS; k++) d2xSpeciesdrds[i][j][k] = 0.0;
        }
        for (j=0; j<NS; j++) dxSpeciesds[i][j] = 0.0;
    }

    /********************************************
    * Create the liquid model Taylor expansion *
    ********************************************/

    if (!contentsOK) {

        gsl_vector *S;
        gsl_matrix *V;
        gsl_matrix_view A;
        gsl_vector_view x, b;


#ifdef DEBUG
        printf("...Constructing Taylor expansion coefficients of liquid model.\n");
#endif
        mTaylor = gsl_matrix_alloc (NP, NP);

        /* First-order Taylor terms */
        for (i=0; i<NE; i++) {
            x = gsl_matrix_row(mTaylor, i);
            initializeTaylor(&x.vector);
            loadTaylor (&x.vector, coeff[i], &rsEndmembers[i][0]);
        }

        /* Second-order Taylor terms: binary interaction parameters
        join A-B ==>  4 ( G(A/2+B/2) - G(A)/2 - G(B)/2 )        */
        m = NE; n = NP; //1 + NV + NV*NV;
        temp = vector_alloc(NE);
        for (i=0; i<NE; i++) {
            for (j=i+1; j<NE; j++) {
                x = gsl_matrix_row(mTaylor, m);
                initializeTaylor(&x.vector);
                for (k=0; k<NV; k++) temp[k] = (rsEndmembers[i][k] + rsEndmembers[j][k])/2.0;
                loadTaylor (&x.vector, (double)  4.0, temp);
                loadTaylor (&x.vector, (double) -2.0, &rsEndmembers[i][0]);
                loadTaylor (&x.vector, (double) -2.0, &rsEndmembers[j][0]);
                m++;
            }
        }
        vector_free(temp,  NE);
        vector_free(coeff, NE);

        // SVD not implemented for M < N in GSL
        if (m < n) gsl_matrix_transpose(mTaylor);

        /* Fill in uninitialized entries */
        bTaylor    = gsl_matrix_alloc((size_t) NP, (size_t) NP);

        taylorCoeff = matrix_alloc(NP, NP);
        S = gsl_vector_alloc((size_t) NP);
        V = gsl_matrix_alloc((size_t) NP, (size_t) NP);
        A = gsl_matrix_submatrix(mTaylor, (size_t) 0, (size_t) 0, (size_t) NP, (size_t) NP);
        x = gsl_vector_view_array(taylorCoeff[0], (size_t) NP);

#ifdef DEBUG
        printf("...Performing singular value analysis of coefficient matrix...\n");
#endif
        gsl_linalg_SV_decomp(&A.matrix, V, S, &x.vector); // first row of taylorCoeff used as working space
        for (i=0, j=0; i<NP; i++) if (fabs(gsl_vector_get(S, i)) < DBL_EPSILON) { gsl_vector_set(S, i, 0.0); j++; }
#ifdef DEBUG
        if (NP-j < n) printf("...Problem is rank deficient! rank = %d\n", NP-j);
        printf("...Performing back-substitution phase...\n");
#endif

        gsl_matrix_set_identity(bTaylor);
        for (i=0; i<NP; i++) {
            b = gsl_matrix_row(bTaylor, i);
            x = gsl_vector_view_array(taylorCoeff[i], (size_t) NP);
            // SVD not implemented for M < N in GSL
            if (m >= n) gsl_linalg_SV_solve(&A.matrix, V, S, &b.vector, &x.vector);
            else gsl_linalg_SV_solve(V, &A.matrix, S, &b.vector, &x.vector);
        }

        for (i=0; i<NP; i++) for (j=0; j<NP; j++)
            if (fabs(taylorCoeff[i][j]) < sqrt(DBL_EPSILON)) taylorCoeff[i][j] = 0.0;

        gsl_matrix_free(V); gsl_vector_free(S);
        gsl_matrix_free(bTaylor); gsl_matrix_free(mTaylor);
        contentsOK = TRUE;

    }

#ifdef DEBUG
    printf("...Exiting function initializeLiquid.\n");
#endif
}

static int rANDsTOx (double r[NR], double s[NT]) {
    static double tolerance;
    double rSum, coeff, dcoeffds[NS], denom, dDenomds[NS];
    int i, j, k, okay = TRUE;
    /* const double y = 0.3; */ /* Fe2SiO(4+2*0.3) or Fe2AlO(3.5+2*0.3) */

    for (i=0, rSum=0.0; i<NR; i++) rSum += r[i];

    coeff = 1.0;
        for (i=0; i<NS; i++) dcoeffds[i] = 0.0;

    /* xSpecies */
    xSpecies[ 0] = 1.0 - rSum*coeff;                   /* SiO2  */
    for (i=0; i<NR; i++) xSpecies[ i+1] = r[i]*coeff;  /* basis */
    for (i=0; i<NS; i++) xSpecies[NA+i] = s[i];        /* depen */

    xSpecies[ 0] += s[0]; // special case SiO2
    xSpecies[10] -= s[0]; // special case CaSiO3
    xSpecies[14] -= s[0]; // special case CO2

    /* Catch bad input data */
    for (i=0;  i<NE; i++) okay &= (xSpecies[i] >= 0.0);
    if (!okay) return okay;

    for (i=NS; i<NT; i++) okay &= ((s[i] > (10.0*DBL_EPSILON)) && (s[i] < (1.0-10.0*DBL_EPSILON)));
    if (okay && (NT > NS)) {
        double yIV = 1.0;
        for (i=NS; i<NT; i++) yIV -= s[i];
        okay &= ((yIV > (10.0*DBL_EPSILON)) && (yIV < (1.0-10.0*DBL_EPSILON)));
    }
    if (!okay) return okay;

    /* Correct roundoff problems - removed check on 4/10/02 when MgO species was included */
    if (tolerance == 0.0) tolerance = pow(DBL_EPSILON, (double) (2.0/3.0));
    /*  for (i=0; i<(NA+NS); i++) if (fabs(xSpecies[i]) < tolerance) xSpecies[i] = 0.0; */

    /* d xSpecies / dr */
    for (i=0; i<NR; i++) {
        dxSpeciesdr[  0][i] = - coeff;  /* SiO2  */
        dxSpeciesdr[i+1][i] =   coeff;  /* other */
    }

    /* d xSpecies / ds */
    for (i=0; i<NS; i++) {
        dxSpeciesds[   0][i] = - rSum*dcoeffds[i]; /* SiO2  */
        for (j=0; j<NR; j++) dxSpeciesds[ j+1][i] =   r[j]*dcoeffds[i]; /* basis */
        dxSpeciesds[NA+i][i] = 1.0;                /* depen */
    }

    dxSpeciesds[ 0][0] += 1.0; // special case SiO2
    dxSpeciesds[10][0] -= 1.0; // special case CaCO3
    dxSpeciesds[14][0] -= 1.0; // special case CO2

    /* d2 xSpecies / dr ds */
    for (i=0; i<NR; i++) {
        for (j=0; j<NS; j++) {
            d2xSpeciesdrds[  0][i][j] = -dcoeffds[j];  /* SiO2  */
            d2xSpeciesdrds[i+1][i][j] =  dcoeffds[j];  /* other */
        }
    }

    /* Total moles of species relative to 1 mole of basis components */
    denom = 12.0;                                    /* Special case */
    dDenomds[0] = 0.0;  // CaCO3

    nSpecies = 12.0/denom;
    for (i=0; i<NS; i++) {
        dnSpeciesds[i] = -12.0*dDenomds[i]/(denom*denom);
        for (j=0; j<NS; j++) {
            d2nSpeciesds2[i][j] = 24.0*dDenomds[i]*dDenomds[j]/(denom*denom*denom);
            for (k=0; k<NS; k++) d3nSpeciesds3[i][j][k] = -72.0*dDenomds[i]*dDenomds[j]*dDenomds[k]/(denom*denom*denom*denom);
        }
    }

    return okay;
}

#define WH(k) (modelParameters[k].enthalpy)
#define WS(k) (modelParameters[k].entropy)
#define WV(k) (modelParameters[k].volume)
#define W(k)  (modelParameters[k].enthalpy - t*modelParameters[k].entropy  + (p-1.0)*modelParameters[k].volume)

#define G(i)       ((liquid[i].cur).g + modelParameters[NW+i].enthalpy - t*modelParameters[NW+i].entropy + (p-1.0)*modelParameters[NW+i].volume)
#define H(i)       ((liquid[i].cur).h + modelParameters[NW+i].enthalpy)
#define S(i)       ((liquid[i].cur).s + modelParameters[NW+i].entropy)
#define V(i)       ((liquid[i].cur).v + modelParameters[NW+i].volume)
#define CP(i)      ((liquid[i].cur).cp)
#define DCPDT(i)   ((liquid[i].cur).dcpdt)
#define DVDT(i)    ((liquid[i].cur).dvdt)
#define DVDP(i)    ((liquid[i].cur).dvdp)
#define D2VDT2(i)  ((liquid[i].cur).d2vdt2)
#define D2VDTDP(i) ((liquid[i].cur).d2vdtdp)
#define D2VDP2(i)  ((liquid[i].cur).d2vdp2)

static void loadTaylorCoefficients(double t, double p)
{
    int i, j, k, l, m, n;

#ifdef DEBUG
    printf("Call to loadTaylorCoefficients ...\n");
#endif

    for (i=0; i<NE; i++)  gibbs(t,   p, (char *) liquid[i].label, &(liquid[i].ref), &(liquid[i].liq), &(liquid[i].fus), &(liquid[i].cur));

    Gconst       = 0.0;  Hconst       = 0.0;  Sconst       = 0.0;  Vconst       = 0.0;
    CPconst      = 0.0;  DCPDTconst   = 0.0;  DVDTconst    = 0.0;  DVDPconst    = 0.0;
    D2VDT2const  = 0.0;  D2VDTDPconst = 0.0;  D2VDP2const  = 0.0;

    memset(gr      , '\0', (size_t) NR*sizeof(double));
    memset(hr      , '\0', (size_t) NR*sizeof(double));
    memset(sr      , '\0', (size_t) NR*sizeof(double));
    memset(vr      , '\0', (size_t) NR*sizeof(double));
    memset(cpr     , '\0', (size_t) NR*sizeof(double));
    memset(dcpdtr  , '\0', (size_t) NR*sizeof(double));
    memset(dvdtr   , '\0', (size_t) NR*sizeof(double));
    memset(dvdpr   , '\0', (size_t) NR*sizeof(double));
    memset(d2vdt2r , '\0', (size_t) NR*sizeof(double));
    memset(d2vdtdpr, '\0', (size_t) NR*sizeof(double));
    memset(d2vdp2r , '\0', (size_t) NR*sizeof(double));
    for (j=0; j<NR; j++) {
        memset(grr[j], '\0', (size_t) NR*sizeof(double));
        memset(hrr[j], '\0', (size_t) NR*sizeof(double));
        memset(srr[j], '\0', (size_t) NR*sizeof(double));
        memset(vrr[j], '\0', (size_t) NR*sizeof(double));
        memset(grs[j], '\0', (size_t) NS*sizeof(double));
        memset(hrs[j], '\0', (size_t) NS*sizeof(double));
        memset(srs[j], '\0', (size_t) NS*sizeof(double));
        memset(vrs[j], '\0', (size_t) NS*sizeof(double));
    }
    memset(gs      , '\0', (size_t) NS*sizeof(double));
    memset(hs      , '\0', (size_t) NS*sizeof(double));
    memset(ss      , '\0', (size_t) NS*sizeof(double));
    memset(vs      , '\0', (size_t) NS*sizeof(double));
    memset(cps     , '\0', (size_t) NS*sizeof(double));
    memset(dcpdts  , '\0', (size_t) NS*sizeof(double));
    memset(dvdts   , '\0', (size_t) NS*sizeof(double));
    memset(dvdps   , '\0', (size_t) NS*sizeof(double));
    memset(d2vdt2s , '\0', (size_t) NS*sizeof(double));
    memset(d2vdtdps, '\0', (size_t) NS*sizeof(double));
    memset(d2vdp2s , '\0', (size_t) NS*sizeof(double));
    for (j=0; j<NS; j++) {
        memset(gss[j], '\0', (size_t) NS*sizeof(double));
        memset(hss[j], '\0', (size_t) NS*sizeof(double));
        memset(sss[j], '\0', (size_t) NS*sizeof(double));
        memset(vss[j], '\0', (size_t) NS*sizeof(double));
    }

    for (i=0; i<NE; i++) {
        Gconst       += taylorCoeff[i][0]*G(i);
        Hconst       += taylorCoeff[i][0]*H(i);
        Sconst       += taylorCoeff[i][0]*S(i);
        Vconst       += taylorCoeff[i][0]*V(i);
        CPconst      += taylorCoeff[i][0]*CP(i);
        DCPDTconst   += taylorCoeff[i][0]*DCPDT(i);
        DVDTconst    += taylorCoeff[i][0]*DVDT(i);
        DVDPconst    += taylorCoeff[i][0]*DVDP(i);
        D2VDT2const  += taylorCoeff[i][0]*D2VDT2(i);
        D2VDTDPconst += taylorCoeff[i][0]*D2VDTDP(i);
        D2VDP2const  += taylorCoeff[i][0]*D2VDP2(i);
        for (j=0; j<NR; j++) {
            gr[j]       += taylorCoeff[i][1+j]*G(i);
            hr[j]       += taylorCoeff[i][1+j]*H(i);
            sr[j]       += taylorCoeff[i][1+j]*S(i);
            vr[j]       += taylorCoeff[i][1+j]*V(i);
            cpr[j]      += taylorCoeff[i][1+j]*CP(i);
            dcpdtr[j]   += taylorCoeff[i][1+j]*DCPDT(i);
            dvdtr[j]    += taylorCoeff[i][1+j]*DVDT(i);
            dvdpr[j]    += taylorCoeff[i][1+j]*DVDP(i);
            d2vdt2r[j]  += taylorCoeff[i][1+j]*D2VDT2(i);
            d2vdtdpr[j] += taylorCoeff[i][1+j]*D2VDTDP(i);
            d2vdp2r[j]  += taylorCoeff[i][1+j]*D2VDP2(i);
        }
        for (j=0; j<NS; j++) {
            gs[j]       += taylorCoeff[i][1+NR+j]*G(i);
            hs[j]       += taylorCoeff[i][1+NR+j]*H(i);
            ss[j]       += taylorCoeff[i][1+NR+j]*S(i);
            vs[j]       += taylorCoeff[i][1+NR+j]*V(i);
            cps[j]      += taylorCoeff[i][1+NR+j]*CP(i);
            dcpdts[j]   += taylorCoeff[i][1+NR+j]*DCPDT(i);
            dvdts[j]    += taylorCoeff[i][1+NR+j]*DVDT(i);
            dvdps[j]    += taylorCoeff[i][1+NR+j]*DVDP(i);
            d2vdt2s[j]  += taylorCoeff[i][1+NR+j]*D2VDT2(i);
            d2vdtdps[j] += taylorCoeff[i][1+NR+j]*D2VDTDP(i);
            d2vdp2s[j]  += taylorCoeff[i][1+NR+j]*D2VDP2(i);
        }
    }
    /* Code below is optimized for speed of execution ... */
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        register double w  = W(n);
        if (w != 0.0) {
            Gconst += taylorCoeff[n+NE][0]*w;
            for (j=0, m=0; j<NR; j++) {
                gr[j] += taylorCoeff[n+NE][1+j]*w;
                for (k=j; k<NR; k++, m++) grr[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*w;
                for (k=0; k<NS; k++, m++) grs[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*w;
            }
            for (j=0; j<NS; j++) {
                gs[j] += taylorCoeff[n+NE][1+NR+j]*w;
                for (k=j; k<NS; k++, m++) if (taylorCoeff[n+NE][1+NR+NS+m] != 0.0) gss[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*w;
            }
        }
    }
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        register double wh = WH(n);
        if (wh != 0.0) {
        Hconst += taylorCoeff[n+NE][0]*wh;
        for (j=0, m=0; j<NR; j++) {
            hr[j] += taylorCoeff[n+NE][1+j]*wh;
            for (k=j; k<NR; k++, m++) hrr[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*wh;
            for (k=0; k<NS; k++, m++) hrs[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*wh;
        }
        for (j=0; j<NS; j++) {
            hs[j] += taylorCoeff[n+NE][1+NR+j]*wh;
            for (k=j; k<NS; k++, m++) if (taylorCoeff[n+NE][1+NR+NS+m] != 0.0) hss[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*wh;
        }
        }
    }
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        register double ws = WS(n);
        if (ws != 0.0) {
        Sconst += taylorCoeff[n+NE][0]*ws;
        for (j=0, m=0; j<NR; j++) {
            sr[j] += taylorCoeff[n+NE][1+j]*ws;
            for (k=j; k<NR; k++, m++) srr[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*ws;
            for (k=0; k<NS; k++, m++) srs[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*ws;
        }
        for (j=0; j<NS; j++) {
            ss[j] += taylorCoeff[n+NE][1+NR+j]*ws;
            for (k=j; k<NS; k++, m++) if (taylorCoeff[n+NE][1+NR+NS+m] != 0.0) sss[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*ws;
        }
        }
    }
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        register double wv = WV(n);
        if (wv != 0.0) {
        Vconst += taylorCoeff[n+NE][0]*wv;
        for (j=0, m=0; j<NR; j++) {
            vr[j] += taylorCoeff[n+NE][1+j]*wv;
            for (k=j; k<NR; k++, m++) vrr[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*wv;
            for (k=0; k<NS; k++, m++) vrs[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*wv;
        }
        for (j=0; j<NS; j++) {
            vs[j] += taylorCoeff[n+NE][1+NR+j]*wv;
            for (k=j; k<NS; k++, m++) if (taylorCoeff[n+NE][1+NR+NS+m] != 0.0) vss[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*wv;
        }
        }
    }
    /* ... end speed optimization. */

    for (j=0; j<NR; j++) for (k=j+1; k<NR; k++) {
        grr[k][j]       = grr[j][k];
        hrr[k][j]       = hrr[j][k];
        srr[k][j]       = srr[j][k];
        vrr[k][j]       = vrr[j][k];
    }

    for (j=0; j<NS; j++) for (k=j+1; k<NS; k++) {
        gss[k][j]       = gss[j][k];
        hss[k][j]       = hss[j][k];
        sss[k][j]       = sss[j][k];
        vss[k][j]       = vss[j][k];
    }

/* ============================================================================== */
/* Do this to insure that mixing properties are reference to a T,P standard State */
/* ============================================================================== */

#ifdef USE_GHIORSO_KRESS_MODEL
    for (i=0; i<NE; i++)  gibbs(t, p, (char *) liquid[i].label, &(liquid[i].ref), &(liquid[i].liq), &(liquid[i].fus), &(liquid[i].cur));
#endif /* USE_GHIORSO_KRESS_MODEL */

}

/* ================================================ */
/* The fill (some thermodynamic property) routines. */
/* ================================================ */

#ifdef USE_GHIORSO_KRESS_MODEL

/* -----> These routines are coded by Maple */


/* ***************************************************************************************
 * These functions return integral properties of V-liquid, as described by the Ghiorso EOS
 *****************************************************************************************/

static double integralV_GK(double r[NR], double s[NT], double t, double p) {
    const double pr      = 1.0;
    const double tr      = 1673.15;
    double m[NA], mOx[NA+1], mOxTot, v, dvdt, cRef, c, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum;
    double gInt;
    int i, j;

    if (fabs(p-pr) < 10.0*DBL_EPSILON) return (double) 0.0;

    /* Convert input composition (r) to liquid moles (m)  */
    for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* Compute moles and total moles of oxides */
    for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
    if (mOxTot == 0.0) return 0.0;

    /* Deal with the special case of FeO1.3 */
    mOx[NA] = 0.0;
    if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
        const double y = 0.3;
        mOx[iOxFeO1_3] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
        }
    }

    for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
        v        += mOx[i]*bulkSystem[i].gk_v;
        dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
        cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
        dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
        cp      += mOx[i]*bulkSystem[i].gk_cp;
        d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
        d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
        d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        mw      += mOx[i]*bulkSystem[i].mw;
    }
    if (v == 0.0) return 0.0;

    alpha   = dvdt/v;
    v0      = v*exp(alpha*(t-tr));
    v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
    c      = cRef + (t-tr)*dcdt;
    v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
    v2      = d2vdp2;
    a      = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;
    b      = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;
    sum      = a*a - 4.0*b;
    gInt    = 0.0;

    if ((a == 0.0) && (b == 0.0)) {
        gInt      += v0*(p-pr) + v1*(p-pr)*(p-pr)/2.0 + v2*(p-pr)*(p-pr)*(p-pr)/6.0;

    } else if ((a != 0.0) && (b == 0.0)) {
        gInt      += (v0 - v2/(2.0*a*a))*(p-pr) + (v1 + v2/(2.0*a))*(p-pr)*(p-pr)/2.0 + v2*log(1.0+a*(p-pr))/(2.0*a*a*a);

    } else if ((a == 0.0) && (b != 0.0)) {
        gInt      += (v0 + v2/(2.0*b))*(p-pr) + v1*log(1.0 + b*(p-pr)*(p-pr))/(2.0*b);
        gInt      += (b > 0.0) ? -v2*atan(sqrt(b)*(p-pr))/(2.0*b*sqrt(b)) : -v2*log((1.0+sqrt(-b)*(p-pr))/(1.0-sqrt(-b)*(p-pr)))/(4.0*b*sqrt(-b));

    } else if (sum > 0.0) {
        double x = sqrt(sum);
        double y = (a + 2.0*b*(p-pr))/x;
        double z = a/x;
        double PcA = (2.0*b*pr - a + x)/(2.0*b);
        double PcB = (2.0*b*pr - a - x)/(2.0*b);
        double arg = (1.0 + y)*(1.0 - z)/((1.0 - y)*(1.0 + z));

        if (((pr < PcA) && (PcA < p)) || ((pr < PcB) && (PcB < p))) printf("v %g,alpha %g, v1Ref %g, v1 %g\nc %g, v2 %g, a %g, b %g\n", v, alpha, v1Ref, v1, c, v2, a, b);
        if (((pr < PcA) && (PcA < p)) || ((pr < PcB) && (PcB < p))) printERR("integralV_GK", ERR_SUM_GT_ZERO, "1.0+a*(p-pr)+b*(p-pr)*(p-pr)", 1.0+a*(p-pr)+b*(p-pr)*(p-pr))
        if (arg <= 0.0)                        printERR("integralV_GK", ERR_SUM_GT_ZERO, "arg", arg)

        gInt      += (v0 + a*v1/b + v2/(2.0*b))*(p-pr);
        gInt      += (v1*(1.0-a*a/b) - v2*a/(2.0*b))*log(1.0+a*(p-pr)+b*(p-pr)*(p-pr))/(2.0*b);
        gInt      += (v1*a*(3.0-a*a/b) + v2*(1.0 - a*a/(2.0*b)))*log(arg)/(2.0*b*x);

    } else if (sum == 0.0) {
        gInt      += (v0 + 4.0*v1/a + 2.0*v2/(a*a))*(p-pr);
        gInt      += -8.0*(v2/a + v1)/(a*a*(2.0+a*(p-pr))) + 4.0*(v2/a + v1)/(a*a);
        gInt      += -4.0*(3.0*v1 + 2.0*v2/a)*log(1.0 + a*(p-pr)/2.0)/(a*a);

    } else if(sum < 0.0) {
        double x = sqrt(-sum);
        double y = (a + 2.0*b*(p-pr))/x;
        double z = a/x;

        gInt      += (v0 + a*v1/b + v2/(2.0*b))*(p-pr);
        gInt      += (v1*(1.0-a*a/b) - v2*a/(2.0*b))*log(1.0+a*(p-pr)+b*(p-pr)*(p-pr))/(2.0*b);
        gInt      += (v1*a*(3.0-a*a/b) + v2*(1.0 - a*a/(2.0*b)))*atan((z-y)/(1.0+z*y))/(b*x);

    }

    return gInt;
}

static double dIntegralV_GKdT(double r[NR], double s[NT], double t, double p) {
    const double pr      = 1.0;
    const double tr      = 1673.15;
    double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum;
    double dv1dt, dgIntdt;
    int i, j;

    if (fabs(p-pr) < 10.0*DBL_EPSILON) return (double) 0.0;

    /* Convert input composition (r) to liquid moles (m)  */
    for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* Compute moles and total moles of oxides */
    for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
    if (mOxTot == 0.0) return 0.0;

    /* Deal with the special case of FeO1.3 */
    mOx[NA] = 0.0;
    if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
        const double y = 0.3;
        mOx[iOxFeO1_3] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
        }
    }

    for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
        v        += mOx[i]*bulkSystem[i].gk_v;
        dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
        cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
        dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
        cp      += mOx[i]*bulkSystem[i].gk_cp;
        d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
        d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
        d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        mw      += mOx[i]*bulkSystem[i].mw;
    }
    if (v == 0.0) return 0.0;

    alpha   = dvdt/v;
    v0      = v*exp(alpha*(t-tr));
    v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
    c      = cRef + (t-tr)*dcdt;
    v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
    v2      = d2vdp2;
    a      = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;
    b      = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;
    sum      = a*a - 4.0*b;
    dv1dt   = - 2.0*v0*v0*alpha*(1000.0/(mw*c*c) + t*alpha*alpha/(cp)) - v0*v0*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp));
    dgIntdt = 0.0;

    if ((a == 0.0) && (b == 0.0)) {
        dgIntdt   += alpha*v0*(p-pr) + dv1dt*(p-pr)*(p-pr)/2.0;

    } else if ((a != 0.0) && (b == 0.0)) {
        dgIntdt   += alpha*v0*(p-pr) + dv1dt*(p-pr)*(p-pr)/2.0;

    } else if ((a == 0.0) && (b != 0.0)) {
        dgIntdt   += alpha*v0*(p-pr) + dv1dt*log(1.0 + b*(p-pr)*(p-pr))/(2.0*b);

    } else if (sum > 0.0) {
        double x = sqrt(sum);
        double y = (a + 2.0*b*(p-pr))/x;
        double z = a/x;
        double PcA = (2.0*b*pr - a + x)/(2.0*b);
        double PcB = (2.0*b*pr - a - x)/(2.0*b);
        double arg = (1.0 + y)*(1.0 - z)/((1.0 - y)*(1.0 + z));

        if (((pr < PcA) && (PcA < p)) || ((pr < PcB) && (PcB < p))) printERR("dIntegralV_GKdT", ERR_SUM_GT_ZERO, "1.0+a*(p-pr)+b*(p-pr)*(p-pr)", 1.0+a*(p-pr)+b*(p-pr)*(p-pr))
        if (arg <= 0.0)                        printERR("dIntegralV_GKdT", ERR_SUM_GT_ZERO, "arg", arg)

        dgIntdt   += (alpha*v0 + a*dv1dt/b)*(p-pr);
        dgIntdt   += dv1dt*(1.0-a*a/b)*log(1.0+a*(p-pr)+b*(p-pr)*(p-pr))/(2.0*b);
        dgIntdt   += dv1dt*a*(3.0-a*a/b)*log(arg)/(2.0*b*x);

    } else if (sum == 0.0) {
        dgIntdt   += (alpha*v0 + 4.0*dv1dt/a)*(p-pr);
        dgIntdt   += -8.0*dv1dt/(a*a*(2.0+a*(p-pr))) + 4.0*dv1dt/(a*a);
        dgIntdt   += -4.0*3.0*dv1dt*log(1.0 + a*(p-pr)/2.0)/(a*a);

    } else if(sum < 0.0) {
        double x = sqrt(-sum);
        double y = (a + 2.0*b*(p-pr))/x;
        double z = a/x;

        dgIntdt   += (alpha*v0 + a*dv1dt/b)*(p-pr);
        dgIntdt   += dv1dt*(1.0-a*a/b)*log(1.0+a*(p-pr)+b*(p-pr)*(p-pr))/(2.0*b);
        dgIntdt   += dv1dt*a*(3.0-a*a/b)*atan((z-y)/(1.0+z*y))/(b*x);

    }

    return dgIntdt;
}

static void dIntegralV_GKdr(double r[NR], double s[NT], double t, double p, double *dr) {
    const double pr     = 1.0;
    const double tr     = 1673.15;
    double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum;
    double dvdr[NR], d2vdrdt[NR], dcRefdr[NR], d2cdrdt[NR], dmwdr[NR], dcpdr[NR], d3vdrdp2[NR], d4vdrdp3[NR], d5vdrdp4[NR], dmOxTotdr[NR], denom;
    int i, j;

    for (i=0; i<NR; i++) dr[i] = 0.0;
    if (fabs(p-pr) < 10.0*DBL_EPSILON) return;

    /* Convert input composition (r) to liquid moles (m)  */
    for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* Compute moles and total moles of oxides */
    for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
    if (mOxTot == 0.0) return;

    /* Deal with the special case of FeO1.3 */
    mOx[NA] = 0.0;
    if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
        const double y = 0.3;
        mOx[iOxFeO1_3] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
        }
    }

    for (i=0; i<NR; i++) {
        dvdr[i]    = 0.0; d2vdrdt[i]  = 0.0; dcRefdr[i]  = 0.0; d2cdrdt[i]  = 0.0; dmwdr[i] = 0.0;
        dcpdr[i]   = 0.0; d3vdrdp2[i] = 0.0; d4vdrdp3[i] = 0.0; d5vdrdp4[i] = 0.0;
        for (j=0, dmOxTotdr[i]=0.0; j<nc; j++) dmOxTotdr[i] += (liquid[i+1].liqToOx)[j] - (liquid[0].liqToOx)[j];
    }

    for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
        v        += mOx[i]*bulkSystem[i].gk_v;
        dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
        cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
        dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
        cp      += mOx[i]*bulkSystem[i].gk_cp;
        d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
        d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
        d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        mw      += mOx[i]*bulkSystem[i].mw;
        for (j=0; j<NR; j++) {
            dvdr[j]       += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_v;
            d2vdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dvdt;
            dmwdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].mw;
            dcpdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_cp;
            dcRefdr[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c/mOxTot
              - mOx[i]*bulkSystem[i].gk_c*dmOxTotdr[j]/(mOxTot*mOxTot);
            dcRefdr[j]   += (iOxAl2O3 != -1) ? ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
              - 2.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot) : 0.0;
            d2cdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt/mOxTot
              - mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotdr[j]/(mOxTot*mOxTot);
            d3vdrdp2[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d4vdrdp3[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d5vdrdp4[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        }
        if (iCmpAl2O3 != -1) dcRefdr[iCmpAl2O3] += ((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*mOx[i]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot);
    }
    if (v == 0.0) return;

    alpha   = dvdt/v;
    v0      = v*exp(alpha*(t-tr));
    v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
    c      = cRef + (t-tr)*dcdt;
    v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
    v2      = d2vdp2;
    denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
    a      = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
    b      = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;
    sum      = a*a - 4.0*b;

    for (i=0; i<NR; i++) {
        double dalphadr   = d2vdrdt[i]/v - dvdt*dvdr[i]/(v*v);
        double dv0dr      = (dvdr[i] + v*dalphadr*(t-tr))*exp(alpha*(t-tr));
        double dv1Refdr   = -2.0*v*dvdr[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
        -v*v*(-1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
                                + 2.0*tr*alpha*dalphadr/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp));
        double dcdr       = dcRefdr[i] + (t-tr)*d2cdrdt[i];
        double dv1dr      = -2.0*v0*dv0dr*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
        -v0*v0*(-1000.0*(c*dmwdr[i]+2.0*mw*dcdr)/(mw*mw*c*c*c)
              + 2.0*t*alpha*dalphadr/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp));
        double dv2dr      = d3vdrdp2[i];
        double ddenomdr   = 2.0*dv1Refdr*d3vdp3+2.0*v1Ref*d4vdrdp3[i]-6.0*d2vdp2*d3vdrdp2[i];
        double dadr       = (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdr*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)/denom
                    - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdr;
        double dbdr       = (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])/denom
                    - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdr;
        double dgIntdr    = 0.0;

        if ((a == 0.0) && (b == 0.0)) {
            dgIntdr  = dv0dr*(p-pr) + dv1dr*(p-pr)*(p-pr)/2.0 + dv2dr*(p-pr)*(p-pr)*(p-pr)/6.0;

        } else if ((a != 0.0) && (b == 0.0)) {
            dgIntdr  = (dv0dr - dv2dr/(2.0*a*a) + v2*dadr/(a*a*a))*(p-pr) + (dv1dr + dv2dr/(2.0*a) - v2*dadr/(2.0*a*a))*(p-pr)*(p-pr)/2.0;
            dgIntdr += (dv2dr*log(1.0+a*(p-pr)) + v2*dadr*(p-pr)/(1.0+a*(p-pr)))/(2.0*a*a*a) - 3.0*dadr*v2*log(1.0+a*(p-pr))/(2.0*a*a*a*a);

        } else if ((a == 0.0) && (b != 0.0)) {
            dgIntdr  = (dv0dr + dv2dr/(2.0*b) - v2*dbdr/(2.0*b*b))*(p-pr);
            dgIntdr += (dv1dr*log(1.0 + b*(p-pr)*(p-pr)) + v1*dbdr*(p-pr)*(p-pr)/(1.0 + b*(p-pr)*(p-pr)))/(2.0*b) - dbdr*v1*log(1.0 + b*(p-pr)*(p-pr))/(2.0*b*b);

            printf("*-->Exception in dIntegralV_GKdr (liquid.c). a is zero, b is greater than zero.\n"); liqERRstate = ERR_A_ZERO;
            dgIntdr  = 0.0;

        } else if (sum > 0.0) {
            dgIntdr  = dgdrGMAP(p/10000.0, pr/10000.0,
                                                    v0,       v1*10000.0,    v2*10000.0*10000.0,    a*10000.0,    b*10000.0*10000.0,
        dv0dr, dv1dr*10000.0, dv2dr*10000.0*10000.0, dadr*10000.0, dbdr*10000.0*10000.0)*10000.0;

        } else if (sum == 0.0) {
            printf("*-->Exception in dIntegralV_GKdr (liquid.c). a*a-4*b is equal to zero.\n"); liqERRstate = ERR_SUM_ZERO;
            dgIntdr  = 0.0;

        } else if(sum < 0.0) {
            dgIntdr  = dgdrLMAP(p/10000.0, pr/10000.0,
                                                    v0,       v1*10000.0,    v2*10000.0*10000.0,    a*10000.0,    b*10000.0*10000.0,
        dv0dr, dv1dr*10000.0, dv2dr*10000.0*10000.0, dadr*10000.0, dbdr*10000.0*10000.0)*10000.0;

        }

        dr[i] = dgIntdr;
    }
}

static void dIntegralV_GKds(double r[NR], double s[NT], double t, double p, double *ds) {
    const double pr = 1.0;
    const double tr = 1673.15;
    const double y  = 0.3;
    double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum;
    double dvds[NS], d2vdsdt[NS], dcRefds[NS], d2cdsdt[NS], dmwds[NS], dcpds[NS], d3vdsdp2[NS], d4vdsdp3[NS], d5vdsdp4[NS], denom;
    int i, j;

    for (i=0; i<NS; i++) ds[i] = 0.0;
    if (fabs(p-pr) < 10.0*DBL_EPSILON) return;

    /* Convert input composition (r) to liquid moles (m)  */
    for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* Compute moles and total moles of oxides */
    for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
    if (mOxTot == 0.0) return;

    /* Deal with the special case of FeO1.3 */
    mOx[NA] = 0.0;
    if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
        const double y = 0.3;
        mOx[iOxFeO1_3] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
        }
    }

    for (i=0; i<NS; i++) {
        dvds[i]    = 0.0; d2vdsdt[i]  = 0.0; dcRefds[i]  = 0.0; d2cdsdt[i]  = 0.0; dmwds[i] = 0.0;
        dcpds[i]   = 0.0; d3vdsdp2[i] = 0.0; d4vdsdp3[i] = 0.0; d5vdsdp4[i] = 0.0;
    }

    for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
        double dmOxds[NS], dmOxTotds[NS];
        v        += mOx[i]*bulkSystem[i].gk_v;
        dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
        cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
        dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
        cp      += mOx[i]*bulkSystem[i].gk_cp;
        d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
        d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
        d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        mw      += mOx[i]*bulkSystem[i].mw;
        for (j=0; j<NS; j++) { dmOxds[j] = 0.0; dmOxTotds[j] = 0.0; }
        if (iCmpFe2SiO4_6 != -1) {
            for (j=0; j<NS; j++) {
                if      (iOxFe2O3  == i) dmOxds[j] += -2.0*y*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
                else if (iOxFeO       == i) dmOxds[j] += -2.0*(1.0-2.0*y)*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
                else if (iOxFeO1_3 == i) dmOxds[j] += 2.0*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
                dmOxTotds[j] += 2.0*y*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
            }
            if      (iOxFe2O3  == i) dmOxds[iCmpFe2SiO4_6] += -2.0*y*nSpecies;
            else if (iOxFeO     == i) dmOxds[iCmpFe2SiO4_6] += -2.0*(1.0-2.0*y)*nSpecies;
            else if (iOxFeO1_3 == i) dmOxds[iCmpFe2SiO4_6] += 2.0*nSpecies;
            dmOxTotds[iCmpFe2SiO4_6] += 2.0*y*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            for (j=0; j<NS; j++) {
                if      (iOxFe2O3  == i) dmOxds[j] += -2.0*y*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
                else if (iOxFeO       == i) dmOxds[j] += -2.0*(1.0-2.0*y)*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
                else if (iOxFeO1_3 == i) dmOxds[j] += 2.0*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
                dmOxTotds[j] += 2.0*y*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
            }
            if      (iOxFe2O3  == i) dmOxds[iCmpFe2AlO4_1] += -2.0*y*nSpecies;
            else if (iOxFeO     == i) dmOxds[iCmpFe2AlO4_1] += -2.0*(1.0-2.0*y)*nSpecies;
            else if (iOxFeO1_3 == i) dmOxds[iCmpFe2AlO4_1] += 2.0*nSpecies;
            dmOxTotds[iCmpFe2AlO4_1] += 2.0*y*nSpecies;
        }
        for (j=0; j<NS; j++) {
            dvds[j]       += dmOxds[j]*bulkSystem[i].gk_v;
            d2vdsdt[j]   += dmOxds[j]*bulkSystem[i].gk_dvdt;
            dmwds[j]     += dmOxds[j]*bulkSystem[i].mw;
            dcpds[j]     += dmOxds[j]*bulkSystem[i].gk_cp;
            dcRefds[j]   += dmOxds[j]*bulkSystem[i].gk_c/mOxTot - mOx[i]*bulkSystem[i].gk_c*dmOxTotds[j]/(mOxTot*mOxTot);
            dcRefds[j]   += (iOxAl2O3 != -1) ? dmOxds[j]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
              - 2.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotds[j]/(mOxTot*mOxTot*mOxTot) : 0.0;
            d2cdsdt[j]   += dmOxds[j]*bulkSystem[i].gk_dcdt/mOxTot - mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotds[j]/(mOxTot*mOxTot);
            d3vdsdp2[j]  += dmOxds[j]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d4vdsdp3[j]  += dmOxds[j]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d5vdsdp4[j]  += dmOxds[j]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        }
    }
    if (v == 0.0) return;

    alpha   = dvdt/v;
    v0      = v*exp(alpha*(t-tr));
    v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
    c      = cRef + (t-tr)*dcdt;
    v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
    v2      = d2vdp2;
    denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
    a      = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
    b      = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;
    sum      = a*a - 4.0*b;

    for (i=0; i<NS; i++) {
        double dalphads   = d2vdsdt[i]/v - dvdt*dvds[i]/(v*v);
        double dv0ds      = (dvds[i] + v*dalphads*(t-tr))*exp(alpha*(t-tr));
        double dv1Refds   = -2.0*v*dvds[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
        -v*v*(-1000.0*(cRef*dmwds[i]+2.0*mw*dcRefds[i])/(mw*mw*cRef*cRef*cRef)
                                + 2.0*tr*alpha*dalphads/(cp) - tr*alpha*alpha*dcpds[i]/(cp*cp));
        double dcds       = dcRefds[i] + (t-tr)*d2cdsdt[i];
        double dv1ds      = -2.0*v0*dv0ds*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
        -v0*v0*(-1000.0*(c*dmwds[i]+2.0*mw*dcds)/(mw*mw*c*c*c)
              + 2.0*t*alpha*dalphads/(cp) - t*alpha*alpha*dcpds[i]/(cp*cp));
        double dv2ds      = d3vdsdp2[i];
        double ddenomds   = 2.0*dv1Refds*d3vdp3+2.0*v1Ref*d4vdsdp3[i]-6.0*d2vdp2*d3vdsdp2[i];
        double dads       = (d3vdsdp2[i]*d3vdp3+d2vdp2*d4vdsdp3[i]-dv1Refds*d4vdp4/2.0-v1Ref*d5vdsdp4[i]/2.0)/denom
                    - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomds;
        double dbds       = (d3vdsdp2[i]*d4vdp4/4.0+d2vdp2*d5vdsdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdsdp3[i])/denom
                    - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomds;
        double dgIntds    = 0.0;

        if ((a == 0.0) && (b == 0.0)) {
            dgIntds  = dv0ds*(p-pr) + dv1ds*(p-pr)*(p-pr)/2.0 + dv2ds*(p-pr)*(p-pr)*(p-pr)/6.0;

        } else if ((a != 0.0) && (b == 0.0)) {
            dgIntds  = (dv0ds - dv2ds/(2.0*a*a) + v2*dads/(a*a*a))*(p-pr) + (dv1ds + dv2ds/(2.0*a) - v2*dads/(2.0*a*a))*(p-pr)*(p-pr)/2.0;
            dgIntds += (dv2ds*log(1.0+a*(p-pr)) + v2*dads*(p-pr)/(1.0+a*(p-pr)))/(2.0*a*a*a) - 3.0*dads*v2*log(1.0+a*(p-pr))/(2.0*a*a*a*a);

        } else if ((a == 0.0) && (b != 0.0)) {
            dgIntds  = (dv0ds + dv2ds/(2.0*b) - v2*dbds/(2.0*b*b))*(p-pr);
            dgIntds += (dv1ds*log(1.0 + b*(p-pr)*(p-pr)) + v1*dbds*(p-pr)*(p-pr)/(1.0 + b*(p-pr)*(p-pr)))/(2.0*b) - dbds*v1*log(1.0 + b*(p-pr)*(p-pr))/(2.0*b*b);

            printf("*-->Exception in dIntegralV_GKds (liquid.c). a is zero, b is greater than zero.\n"); liqERRstate = ERR_A_ZERO;
            dgIntds  = 0.0;

        } else if (sum > 0.0) {
            eosIntegralBranch = GMAPeosBRANCH;
            dgIntds  = dgdrGMAP(p/10000.0, pr/10000.0,
                                                    v0,       v1*10000.0,    v2*10000.0*10000.0,    a*10000.0,    b*10000.0*10000.0,
                                                    dv0ds, dv1ds*10000.0, dv2ds*10000.0*10000.0, dads*10000.0, dbds*10000.0*10000.0)*10000.0;
            eosIntDGDS[i] = dgIntds;

        } else if (sum == 0.0) {
            printf("*-->Exception in dIntegralV_GKds (liquid.c). a*a-4*b is equal to zero.\n"); liqERRstate = ERR_SUM_ZERO;
            dgIntds  = 0.0;

        } else if(sum < 0.0) {
            eosIntegralBranch = LMAPeosBRANCH;
            dgIntds  = dgdsLMAP(p/10000.0, pr/10000.0,
                                                    v0,       v1*10000.0,    v2*10000.0*10000.0,    a*10000.0,    b*10000.0*10000.0,
        dv0ds, dv1ds*10000.0, dv2ds*10000.0*10000.0, dads*10000.0, dbds*10000.0*10000.0)*10000.0;
            eosIntDGDS[i] = dgIntds;

        }

        ds[i] = dgIntds;
    }
}

static void d2IntegralV_GKdr2(double r[NR], double s[NT], double t, double p, double d2r[NR][NR]) {
    const double pr      = 1.0;
    const double tr      = 1673.15;
    double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum;
    double dvdr[NR], d2vdrdt[NR], dcRefdr[NR], d2cdrdt[NR], dmwdr[NR], dcpdr[NR], d3vdrdp2[NR], d4vdrdp3[NR], d5vdrdp4[NR], dmOxTotdr[NR], denom;
    double d2cRefdr2[NR][NR], d3cdr2dt[NR][NR];
    int i, j, k;

    for (i=0; i<NR; i++) for (j=0; j<NR; j++) d2r[i][j] = 0.0;
    if (fabs(p-pr) < 10.0*DBL_EPSILON) return;

    /* Convert input composition (r) to liquid moles (m)  */
    for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* Compute moles and total moles of oxides */
    for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
    if (mOxTot == 0.0) return;

    /* Deal with the special case of FeO1.3 */
    mOx[NA] = 0.0;
    if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
        const double y = 0.3;
        mOx[iOxFeO1_3] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
        }
    }

    for (i=0; i<NR; i++) {
        dvdr[i]    = 0.0; d2vdrdt[i]  = 0.0; dcRefdr[i]  = 0.0; d2cdrdt[i]  = 0.0; dmwdr[i] = 0.0;
        dcpdr[i]   = 0.0; d3vdrdp2[i] = 0.0; d4vdrdp3[i] = 0.0; d5vdrdp4[i] = 0.0;
        for (j=0, dmOxTotdr[i]=0.0; j<nc; j++) dmOxTotdr[i] += (liquid[i+1].liqToOx)[j] - (liquid[0].liqToOx)[j];
        for (j=i; j<NR; j++) { d2cRefdr2[i][j] = 0.0; d3cdr2dt[i][j] = 0.0; }
    }

    for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
        v        += mOx[i]*bulkSystem[i].gk_v;
        dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
        cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
        dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
        cp      += mOx[i]*bulkSystem[i].gk_cp;
        d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
        d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
        d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        mw      += mOx[i]*bulkSystem[i].mw;

        for (j=0; j<NR; j++) {
            dvdr[j]       += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_v;
            d2vdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dvdt;
            dmwdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].mw;
            dcpdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_cp;
            dcRefdr[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c/mOxTot
              - mOx[i]*bulkSystem[i].gk_c*dmOxTotdr[j]/(mOxTot*mOxTot);
            dcRefdr[j]   += (iOxAl2O3 != -1) ? ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
              - 2.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot) : 0.0;
            d2cdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt/mOxTot
              - mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotdr[j]/(mOxTot*mOxTot);
            d3vdrdp2[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d4vdrdp3[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d5vdrdp4[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        }
        if (iCmpAl2O3 != -1) dcRefdr[iCmpAl2O3] += ((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*mOx[i]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot);

        for (j=0; j<NR; j++) {
            for (k=j; k<NR; k++) {
                d2cRefdr2[j][k] += -((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c*dmOxTotdr[k]/(mOxTot*mOxTot)
             -  ((liquid[k+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c*dmOxTotdr[j]/(mOxTot*mOxTot)
             + 2.0*mOx[i]*bulkSystem[i].gk_c*dmOxTotdr[j]*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot);
                d2cRefdr2[j][k] += (iOxAl2O3 != -1) ?
                                    - 2.0*((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot)
             - 2.0*((liquid[k+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot)
        + 6.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot*mOxTot)
             : 0.0;
                d3cdr2dt[j][k]  += -((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt*dmOxTotdr[k]/(mOxTot*mOxTot)
             -  ((liquid[k+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt*dmOxTotdr[j]/(mOxTot*mOxTot)
             + 2.0*mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotdr[j]*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot);
            }
            if (iCmpAl2O3 != -1) d2cRefdr2[(j<iCmpAl2O3) ? j : iCmpAl2O3][(j>iCmpAl2O3) ? j : iCmpAl2O3] +=
                ((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
                - 2.0*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*mOx[i]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot);
        }
        if (iCmpAl2O3 != -1) d2cRefdr2[iCmpAl2O3][iCmpAl2O3] +=
            ((liquid[iCmpAl2O3+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
            - 2.0*mOx[i]*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*bulkSystem[i].gk_cXal2o3*dmOxTotdr[iCmpAl2O3]/(mOxTot*mOxTot*mOxTot);

    }
    if (v == 0.0) return;

    alpha   = dvdt/v;
    v0      = v*exp(alpha*(t-tr));
    v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
    c      = cRef + (t-tr)*dcdt;
    v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
    v2      = d2vdp2;
    denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
    a      = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
    b      = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;
    sum      = a*a - 4.0*b;

    for (i=0; i<NR; i++) {
        double dalphadri   = d2vdrdt[i]/v - dvdt*dvdr[i]/(v*v);
        double dv0dri      = (dvdr[i] + v*dalphadri*(t-tr))*exp(alpha*(t-tr));
        double dv1Refdri   = -2.0*v*dvdr[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
            -v*v*(- 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
                   + 2.0*tr*alpha*dalphadri/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp));
        double dcdri       = dcRefdr[i] + (t-tr)*d2cdrdt[i];
        double dv1dri      = -2.0*v0*dv0dri*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
            -v0*v0*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)/(mw*mw*c*c*c)
                        + 2.0*t*alpha*dalphadri/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp));
        double ddenomdri   = 2.0*dv1Refdri*d3vdp3+2.0*v1Ref*d4vdrdp3[i]-6.0*d2vdp2*d3vdrdp2[i];
        double dadri       = (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdri*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)/denom
             - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdri;
        double dbdri       = (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])/denom
             - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdri;

        for (j=i; j<NR; j++) {
            double dalphadrj   = d2vdrdt[j]/v - dvdt*dvdr[j]/(v*v);
            double dv0drj     = (dvdr[j] + v*dalphadrj*(t-tr))*exp(alpha*(t-tr));
            double dv1Refdrj   = -2.0*v*dvdr[j]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                -v*v*(- 1000.0*(cRef*dmwdr[j]+2.0*mw*dcRefdr[j])/(mw*mw*cRef*cRef*cRef)
           + 2.0*tr*alpha*dalphadrj/(cp) - tr*alpha*alpha*dcpdr[j]/(cp*cp));
            double dcdrj     = dcRefdr[j] + (t-tr)*d2cdrdt[j];
            double dv1drj     = -2.0*v0*dv0drj*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                -v0*v0*(- 1000.0*(c*dmwdr[j]+2.0*mw*dcdrj)/(mw*mw*c*c*c)
             + 2.0*t*alpha*dalphadrj/(cp) - t*alpha*alpha*dcpdr[j]/(cp*cp));
            double ddenomdrj   = 2.0*dv1Refdrj*d3vdp3+2.0*v1Ref*d4vdrdp3[j]-6.0*d2vdp2*d3vdrdp2[j];
            double dadrj     = (d3vdrdp2[j]*d3vdp3+d2vdp2*d4vdrdp3[j]-dv1Refdrj*d4vdp4/2.0-v1Ref*d5vdrdp4[j]/2.0)/denom
            - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdrj;
            double dbdrj     = (d3vdrdp2[j]*d4vdp4/4.0+d2vdp2*d5vdrdp4[j]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[j])/denom
             - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdrj;

            double d2alphadr2   = -d2vdrdt[i]*dvdr[j]/(v*v)-d2vdrdt[j]*dvdr[i]/(v*v)+2.0*dvdt*dvdr[i]*dvdr[j]/(v*v*v);
            double d2v0dr2      = dvdr[i]*dalphadrj*(t-tr)*exp(alpha*(t-tr))
                        + dvdr[j]*dalphadri*(t-tr)*exp(alpha*(t-tr))
                        + v*d2alphadr2*(t-tr)*exp(alpha*(t-tr))
                        + v*dalphadri*pow(t-tr,2.0)*dalphadrj*exp(alpha*(t-tr));
            double d2v1Refdr2   = -2.0*dvdr[j]*dvdr[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                        - 2.0*v*dvdr[i]*(- 1000.0*(cRef*dmwdr[j]+2.0*mw*dcRefdr[j])/(mw*mw*cRef*cRef*cRef)
                   + 2.0*tr*alpha*dalphadrj/(cp) - tr*alpha*alpha*dcpdr[j]/(cp*cp))
          - 2.0*v*dvdr[j]*(- 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
                    + 2.0*tr*alpha*dalphadri/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp))
               -v*v*(- 1000.0*(dcRefdr[j]*dmwdr[i]+2.0*dmwdr[j]*dcRefdr[i]+2.0*mw*d2cRefdr2[i][j])/(mw*mw*cRef*cRef*cRef)
                        + 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])*(2.0*dmwdr[j]*cRef+3.0*mw*dcRefdr[j])/(mw*mw*mw*cRef*cRef*cRef*cRef)
                        + 2.0*tr*alpha*d2alphadr2/(cp)
                        + 2.0*tr*dalphadrj*dalphadri/(cp) - 2.0*tr*alpha*dalphadri*dcpdr[j]/(cp*cp)
                        - 2.0*tr*alpha*dalphadrj*dcpdr[i]/(cp*cp) + 2.0*tr*alpha*alpha*dcpdr[i]*dcpdr[j]/(cp*cp*cp));
            double d2cdr2     = d2cRefdr2[i][j] + (t-tr)*d3cdr2dt[i][j];
            double d2v1dr2     = -2.0*dv0drj*dv0dri*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                        - 2.0*v0*d2v0dr2*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                        - 2.0*v0*dv0dri*(- 1000.0*(c*dmwdr[j]+2.0*mw*dcdrj)/(mw*mw*c*c*c)
                   + 2.0*t*alpha*dalphadrj/(cp) - t*alpha*alpha*dcpdr[j]/(cp*cp))
          - 2.0*v0*dv0drj*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)/(mw*mw*c*c*c)
                    + 2.0*t*alpha*dalphadri/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp))
               -v0*v0*(- 1000.0*(dcdrj*dmwdr[i]+2.0*dmwdr[j]*dcdri+2.0*mw*d2cdr2)/(mw*mw*c*c*c)
                            + 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)*(2.0*dmwdr[j]*c+3.0*mw*dcdrj)/(mw*mw*mw*c*c*c*c)
                            + 2.0*t*alpha*d2alphadr2/(cp)
                            + 2.0*t*dalphadrj*dalphadri/(cp) - 2.0*t*alpha*dalphadri*dcpdr[j]/(cp*cp)
                            - 2.0*t*alpha*dalphadrj*dcpdr[i]/(cp*cp) + 2.0*t*alpha*alpha*dcpdr[i]*dcpdr[j]/(cp*cp*cp));
            double d2denomdr2   = 2.0*d2v1Refdr2*d3vdp3 + 2.0*dv1Refdri*d4vdrdp3[j] + 2.0*dv1Refdrj*d4vdrdp3[i] - 6.0*d3vdrdp2[j]*d3vdrdp2[i];
            double d2adr2     = (d3vdrdp2[i]*d4vdrdp3[j]+d3vdrdp2[j]*d4vdrdp3[i]-d2v1Refdr2*d4vdp4/2.0-dv1Refdri*d5vdrdp4[j]/2.0-dv1Refdrj*d5vdrdp4[i]/2.0)/denom
             - (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdri*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)*ddenomdrj/(denom*denom)
            - (d3vdrdp2[j]*d3vdp3+d2vdp2*d4vdrdp3[j]-dv1Refdrj*d4vdp4/2.0-v1Ref*d5vdrdp4[j]/2.0)*ddenomdri/(denom*denom)
            + 2.0*(d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*ddenomdri*ddenomdrj/(denom*denom*denom)
            - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*d2denomdr2/(denom*denom);
            double d2bdr2     = (d3vdrdp2[i]*d5vdrdp4[j]/4.0+d3vdrdp2[j]*d5vdrdp4[i]/4.0-2.0/3.0*d4vdrdp3[j]*d4vdrdp3[i])/denom
            - (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])*ddenomdrj/(denom*denom)
            - (d3vdrdp2[j]*d4vdp4/4.0+d2vdp2*d5vdrdp4[j]/4.0-2.0*d4vdrdp3[j]*d3vdp3/3.0)*ddenomdri/(denom*denom)
             + 2.0*(d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)*ddenomdri*ddenomdrj/(denom*denom*denom)
             - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)*d2denomdr2/(denom*denom);
            double d2gIntdr2   = 0.0;

            if ((a == 0.0) && (b == 0.0)) {
    d2gIntdr2  = d2v0dr2*(p-pr) + d2v1dr2*(p-pr)*(p-pr)/2.0;

            } else if ((a != 0.0) && (b == 0.0)) {
    printf("*-->Exception in d2IntegralV_GKdr2 (liquid.c). a is not equal to zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
    d2gIntdr2  = 0.0;

            } else if ((a == 0.0) && (b != 0.0)) {
    printf("*-->Exception in d2IntegralV_GKdr2 (liquid.c). a is zero, b is not equal to zero.\n"); liqERRstate = ERR_A_ZERO;
    d2gIntdr2  = 0.0;

            } else if (sum > 0.0) {
    d2gIntdr2  = d2gdr2GMAP(p/10000.0, pr/10000.0,
                          v0,           v1*10000.0,      d2vdp2*10000.0*10000.0,      a*10000.0,      b*10000.0*10000.0,
        dv0dri,   dv1dri*10000.0, d3vdrdp2[i]*10000.0*10000.0,  dadri*10000.0,  dbdri*10000.0*10000.0,
            dv0drj,   dv1drj*10000.0, d3vdrdp2[j]*10000.0*10000.0,  dadrj*10000.0,  dbdrj*10000.0*10000.0,
              d2v0dr2, d2v1dr2*10000.0,                         0.0, d2adr2*10000.0, d2bdr2*10000.0*10000.0)*10000.0;

            } else if (sum == 0.0) {
    printf("*-->Exception in d2IntegralV_GKdr2 (liquid.c). a*a-4*b is equal to zero.\n"); liqERRstate = ERR_SUM_ZERO;
    d2gIntdr2  = 0.0;

            } else if(sum < 0.0) {
    d2gIntdr2  = d2gdr2LMAP(p/10000.0, pr/10000.0,
                          v0,           v1*10000.0,      d2vdp2*10000.0*10000.0,      a*10000.0,      b*10000.0*10000.0,
        dv0dri,   dv1dri*10000.0, d3vdrdp2[i]*10000.0*10000.0,  dadri*10000.0,  dbdri*10000.0*10000.0,
            dv0drj,   dv1drj*10000.0, d3vdrdp2[j]*10000.0*10000.0,  dadrj*10000.0,  dbdrj*10000.0*10000.0,
              d2v0dr2, d2v1dr2*10000.0,                         0.0, d2adr2*10000.0, d2bdr2*10000.0*10000.0)*10000.0;

            }

            d2r[i][j] += d2gIntdr2;
            if (i != j) d2r[j][i] += d2gIntdr2;

        }
    }

}

static void d2IntegralV_GKdrds(double r[NR], double s[NT], double t, double p, double drds[NR][NS]) {
    const double pr = 1.0;
    const double tr = 1673.15;
    const double y  = 0.3;
    double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum;
    double dvdr[NR], d2vdrdt[NR], dcRefdr[NR], d2cdrdt[NR], dmwdr[NR], dcpdr[NR], d3vdrdp2[NR], d4vdrdp3[NR], d5vdrdp4[NR], dmOxTotdr[NR], denom;
    double dvds[NS], d2vdsdt[NS], dcRefds[NS], d2cdsdt[NS], dmwds[NS], dcpds[NS], d3vdsdp2[NS], d4vdsdp3[NS], d5vdsdp4[NS], dmOxTotds[NS];
    double d2vdrds[NR][NS], d3vdrdsdt[NR][NS], d2cRefdrds[NR][NS], d3cdrdsdt[NR][NS], d2mwdrds[NR][NS], d2cpdrds[NR][NS], d4vdrdsdp2[NR][NS],
         d5vdrdsdp3[NR][NS], d6vdrdsdp4[NR][NS], d2mOxTotdrds[NR][NS];
    int i, j, k;

    for (i=0; i<NR; i++) for (j=0; j<NS; j++) drds[i][j] = 0.0;
    if (fabs(p-pr) < 10.0*DBL_EPSILON) return;

    /* Convert input composition (r) to liquid moles (m)  */
    for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* Compute moles and total moles of oxides */
    for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
    if (mOxTot == 0.0) return;

    /* Deal with the special case of FeO1.3 */
    mOx[NA] = 0.0;
    if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
        const double y = 0.3;
        mOx[iOxFeO1_3] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
        }
    }

    for (i=0; i<NR; i++) {
        dvdr[i]    = 0.0; d2vdrdt[i]  = 0.0; dcRefdr[i]  = 0.0; d2cdrdt[i]  = 0.0; dmwdr[i] = 0.0;
        dcpdr[i]   = 0.0; d3vdrdp2[i] = 0.0; d4vdrdp3[i] = 0.0; d5vdrdp4[i] = 0.0;
        for (j=0, dmOxTotdr[i]=0.0; j<nc; j++) dmOxTotdr[i] += (liquid[i+1].liqToOx)[j] - (liquid[0].liqToOx)[j];
    }
    for (i=0; i<NS; i++) {
        dvds[i]    = 0.0; d2vdsdt[i]  = 0.0; dcRefds[i]  = 0.0; d2cdsdt[i]  = 0.0; dmwds[i] = 0.0;
        dcpds[i]   = 0.0; d3vdsdp2[i] = 0.0; d4vdsdp3[i] = 0.0; d5vdsdp4[i] = 0.0; dmOxTotds[i] = 0.0;
    }
    if (iCmpFe2SiO4_6 != -1) {
        for (i=0; i<NS; i++) dmOxTotds[i] += 2.0*y*s[iCmpFe2SiO4_6]*dnSpeciesds[i];
        dmOxTotds[iCmpFe2SiO4_6] += 2.0*y*nSpecies;
    }
    if (iCmpFe2AlO4_1 != -1) {
        for (i=0; i<NS; i++) dmOxTotds[i] += 2.0*y*s[iCmpFe2AlO4_1]*dnSpeciesds[i];
        dmOxTotds[iCmpFe2AlO4_1] += 2.0*y*nSpecies;
    }
    for (j=0; j<NR; j++) for (k=0; k<NS; k++) {
        d2mOxTotdrds[j][k] = 0.0; d2vdrds[j][k]      = 0.0; d3vdrdsdt[j][k]    = 0.0; d2cRefdrds[j][k]   = 0.0; d3cdrdsdt[j][k]    = 0.0;
        d2mwdrds[j][k]     = 0.0; d2cpdrds[j][k]     = 0.0; d4vdrdsdp2[j][k]   = 0.0; d5vdrdsdp3[j][k]   = 0.0; d6vdrdsdp4[j][k]   = 0.0;
    }

    for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
        double dmOxds[NS];
        v        += mOx[i]*bulkSystem[i].gk_v;
        dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
        cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
        dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
        cp      += mOx[i]*bulkSystem[i].gk_cp;
        d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
        d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
        d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        mw      += mOx[i]*bulkSystem[i].mw;

        for (j=0; j<NS; j++) dmOxds[j] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            for (j=0; j<NS; j++) {
                if      (iOxFe2O3  == i) dmOxds[j] += -2.0*y*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
                else if (iOxFeO       == i) dmOxds[j] += -2.0*(1.0-2.0*y)*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
                else if (iOxFeO1_3 == i) dmOxds[j] += 2.0*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
            }
            if      (iOxFe2O3  == i) dmOxds[iCmpFe2SiO4_6] += -2.0*y*nSpecies;
            else if (iOxFeO     == i) dmOxds[iCmpFe2SiO4_6] += -2.0*(1.0-2.0*y)*nSpecies;
            else if (iOxFeO1_3 == i) dmOxds[iCmpFe2SiO4_6] += 2.0*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            for (j=0; j<NS; j++) {
                if      (iOxFe2O3  == i) dmOxds[j] += -2.0*y*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
                else if (iOxFeO       == i) dmOxds[j] += -2.0*(1.0-2.0*y)*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
                else if (iOxFeO1_3 == i) dmOxds[j] += 2.0*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
            }
            if      (iOxFe2O3  == i) dmOxds[iCmpFe2AlO4_1] += -2.0*y*nSpecies;
            else if (iOxFeO     == i) dmOxds[iCmpFe2AlO4_1] += -2.0*(1.0-2.0*y)*nSpecies;
            else if (iOxFeO1_3 == i) dmOxds[iCmpFe2AlO4_1] += 2.0*nSpecies;
        }
        for (j=0; j<NS; j++) {
            dvds[j]       += dmOxds[j]*bulkSystem[i].gk_v;
            d2vdsdt[j]   += dmOxds[j]*bulkSystem[i].gk_dvdt;
            dmwds[j]     += dmOxds[j]*bulkSystem[i].mw;
            dcpds[j]     += dmOxds[j]*bulkSystem[i].gk_cp;
            dcRefds[j]   += dmOxds[j]*bulkSystem[i].gk_c/mOxTot - mOx[i]*bulkSystem[i].gk_c*dmOxTotds[j]/(mOxTot*mOxTot);
            dcRefds[j]   += (iOxAl2O3 != -1) ? dmOxds[j]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
              - 2.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotds[j]/(mOxTot*mOxTot*mOxTot) : 0.0;
            d2cdsdt[j]   += dmOxds[j]*bulkSystem[i].gk_dcdt/mOxTot - mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotds[j]/(mOxTot*mOxTot);
            d3vdsdp2[j]  += dmOxds[j]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d4vdsdp3[j]  += dmOxds[j]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d5vdsdp4[j]  += dmOxds[j]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        }

        for (j=0; j<NR; j++) {
            dvdr[j]       += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_v;
            d2vdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dvdt;
            dmwdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].mw;
            dcpdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_cp;
            dcRefdr[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c/mOxTot
              - mOx[i]*bulkSystem[i].gk_c*dmOxTotdr[j]/(mOxTot*mOxTot);
            dcRefdr[j]   += (iOxAl2O3 != -1) ? ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
              - 2.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot) : 0.0;
            d2cdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt/mOxTot
              - mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotdr[j]/(mOxTot*mOxTot);
            d3vdrdp2[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d4vdrdp3[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d5vdrdp4[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        }
        if (iCmpAl2O3 != -1) dcRefdr[iCmpAl2O3] += ((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*mOx[i]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot);
        for (j=0; j<NR; j++) for (k=0; k<NS; k++) {
            d2cRefdrds[j][k] += -((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c*dmOxTotds[k]/(mOxTot*mOxTot)
                  - (dmOxds[k]*bulkSystem[i].gk_c*dmOxTotdr[j] + mOx[i]*bulkSystem[i].gk_c*d2mOxTotdrds[j][k])/(mOxTot*mOxTot)
                  + 2.0*mOx[i]*bulkSystem[i].gk_c*dmOxTotdr[j]*dmOxTotds[k]/(mOxTot*mOxTot*mOxTot);
            d2cRefdrds[j][k] += (iOxAl2O3 != -1) ? -2.0*((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotds[k]/(mOxTot*mOxTot*mOxTot)
                  - 2.0*(dmOxds[k]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j] + mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*d2mOxTotdrds[j][k])/(mOxTot*mOxTot*mOxTot)
      + 6.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]*dmOxTotds[k]/(mOxTot*mOxTot*mOxTot*mOxTot)
      : 0.0;
            d3cdrdsdt[j][k]  += -((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt*dmOxTotds[k]/(mOxTot*mOxTot)
                  - (dmOxds[k]*bulkSystem[i].gk_dcdt*dmOxTotdr[j] + mOx[i]*bulkSystem[i].gk_dcdt*d2mOxTotdrds[j][k])/(mOxTot*mOxTot)
      + 2.0*mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotdr[j]*dmOxTotds[k]/(mOxTot*mOxTot*mOxTot);
        }
        if (iCmpAl2O3 != -1) for (k=0; k<NS; k++) {
            d2cRefdrds[iCmpAl2O3][k] += ((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*dmOxds[k]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
                                                                - 2.0*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*mOx[i]*bulkSystem[i].gk_cXal2o3*dmOxTotds[k]/(mOxTot*mOxTot*mOxTot);
        }
    }
    if (v == 0.0) return;

    alpha   = dvdt/v;
    v0      = v*exp(alpha*(t-tr));
    v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
    c      = cRef + (t-tr)*dcdt;
    v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
    v2      = d2vdp2;
    denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
    a      = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
    b      = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;
    sum      = a*a - 4.0*b;

    for (i=0; i<NR; i++) {
        double dalphadr   = d2vdrdt[i]/v - dvdt*dvdr[i]/(v*v);
        double dv0dr      = (dvdr[i] + v*dalphadr*(t-tr))*exp(alpha*(t-tr));
        double dv1Refdr   = -2.0*v*dvdr[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
        -v*v*(-1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
                                + 2.0*tr*alpha*dalphadr/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp));
        double dcdr       = dcRefdr[i] + (t-tr)*d2cdrdt[i];
        double dv1dr      = -2.0*v0*dv0dr*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
        -v0*v0*(-1000.0*(c*dmwdr[i]+2.0*mw*dcdr)/(mw*mw*c*c*c)
              + 2.0*t*alpha*dalphadr/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp));
        double dv2dr      = d3vdrdp2[i];
        double ddenomdr   = 2.0*dv1Refdr*d3vdp3+2.0*v1Ref*d4vdrdp3[i]-6.0*d2vdp2*d3vdrdp2[i];
        double dadr       = (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdr*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)/denom
                    - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdr;
        double dbdr       = (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])/denom
                    - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdr;
        for (j=0; j<NS; j++) {
            double dalphads    = d2vdsdt[j]/v - dvdt*dvds[j]/(v*v);
            double dv0ds    = (dvds[j] + v*dalphads*(t-tr))*exp(alpha*(t-tr));
            double dv1Refds    = -2.0*v*dvds[j]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                    -v*v*(-1000.0*(cRef*dmwds[j]+2.0*mw*dcRefds[j])/(mw*mw*cRef*cRef*cRef)
             + 2.0*tr*alpha*dalphads/(cp) - tr*alpha*alpha*dcpds[j]/(cp*cp));
            double dcds    = dcRefds[j] + (t-tr)*d2cdsdt[j];
            double dv1ds    = -2.0*v0*dv0ds*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                    -v0*v0*(-1000.0*(c*dmwds[j]+2.0*mw*dcds)/(mw*mw*c*c*c)
               + 2.0*t*alpha*dalphads/(cp) - t*alpha*alpha*dcpds[j]/(cp*cp));
            double dv2ds    = d3vdsdp2[j];
            double ddenomds    = 2.0*dv1Refds*d3vdp3+2.0*v1Ref*d4vdsdp3[j]-6.0*d2vdp2*d3vdsdp2[j];
            double dads    = (d3vdsdp2[j]*d3vdp3+d2vdp2*d4vdsdp3[j]-dv1Refds*d4vdp4/2.0-v1Ref*d5vdsdp4[j]/2.0)/denom
                - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomds;
            double dbds    = (d3vdsdp2[j]*d4vdp4/4.0+d2vdp2*d5vdsdp4[j]/4.0-2.0/3.0*d3vdp3*d4vdsdp3[j])/denom
                - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomds;

            double d2alphadrds = d3vdrdsdt[i][j]/v - d2vdrdt[i]*dvds[j]/(v*v) - (d2vdsdt[j]*dvdr[i]+dvdt*d2vdrds[i][j])/(v*v)
                         + 2.0*dvdt*dvdr[i]*dvds[j]/(v*v*v);
            double d2v0drds     = (d2vdrds[i][j] + dvds[j]*dalphadr*(t-tr) + v*d2alphadrds*(t-tr))*exp(alpha*(t-tr))
                         + (dvdr[i] + v*dalphadr*(t-tr))*dalphads*(t-tr)*exp(alpha*(t-tr));
            double d2v1Refdrds = -2.0*(dvds[j]*dvdr[i]+v*d2vdrds[i][j])*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                         - 2.0*v*dvdr[i]*(-1000.0*(cRef*dmwds[j]+2.0*mw*dcRefds[j])/(mw*mw*cRef*cRef*cRef)
              + 2.0*tr*alpha*dalphads/(cp) - tr*alpha*alpha*dcpds[j]/(cp*cp))
                    - 2.0*v*dvds[j]*(-1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
              + 2.0*tr*alpha*dalphadr/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp))
        - v*v*(-1000.0*(dcRefds[j]*dmwdr[i]+cRef*d2mwdrds[i][j]+2.0*dmwds[j]*dcRefdr[i]+2.0*mw*d2cRefdrds[i][j])
                      /(mw*mw*cRef*cRef*cRef)
              +1000.0*(2.0*dmwds[j]*cRef+3.0*mw*dcRefds[j])*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])
                        /(mw*mw*mw*cRef*cRef*cRef*cRef)
              + 2.0*tr*(dalphads*dalphadr+alpha*d2alphadrds)/(cp)
        - 2.0*tr*alpha*dalphadr*dcpds[j]/(cp*cp)
        - tr*(2.0*alpha*dalphads*dcpdr[i]+alpha*alpha*d2cpdrds[i][j])/(cp*cp)
        + 2.0*tr*alpha*alpha*dcpdr[i]*dcpds[j]/(cp*cp*cp)
        );
            double d2cdrds    = d2cRefdrds[i][j] + (t-tr)*d3cdrdsdt[i][j];
            double d2v1drds    = -2.0*(dv0ds*dv0dr+v0*d2v0drds)*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                                                - 2.0*v0*dv0dr*(-(dmwds[j]*c+2.0*mw*dcds)*1000.0/(mw*mw*c*c*c) + 2.0*t*alpha*dalphads/(cp)
                      - t*alpha*alpha*dcpds[j]/(cp*cp))
            - 2.0*v0*dv0ds*(-1000.0*(c*dmwdr[i]+2.0*mw*dcdr)/(mw*mw*c*c*c)
                      + 2.0*t*alpha*dalphadr/(cp)
                        - t*alpha*alpha*dcpdr[i]/(cp*cp)
                        )
      - v0*v0*(-1000.0*(dcds*dmwdr[i]+c*d2mwdrds[i][j]+2.0*dmwds[j]*dcdr+2.0*mw*d2cdrds)/(mw*mw*c*c*c)
                        +1000.0*(2.0*dmwds[j]*c+3.0*mw*dcds)*(c*dmwdr[i]+2.0*mw*dcdr)/(mw*mw*mw*c*c*c*c)
                + 2.0*t*(dalphads*dalphadr+alpha*d2alphadrds)/(cp)
            - 2.0*t*alpha*dalphadr*dcpds[j]/(cp*cp)
            - t*(2.0*alpha*dalphads*dcpdr[i]+alpha*alpha*d2cpdrds[i][j])/(cp*cp)
            + 2.0*t*alpha*alpha*dcpdr[i]*dcpds[j]/(cp*cp*cp)

            );
            double d2v2drds    = d4vdrdsdp2[i][j];
            double d2denomdrds= 2.0*d2v1Refdrds*d3vdp3 + 2.0*dv1Refdr*d4vdsdp3[j]
                                                + 2.0*dv1Refds*d4vdrdp3[i] + 2.0*v1Ref*d5vdrdsdp3[i][j]
      - 6.0*d3vdsdp2[j]*d3vdrdp2[i] - 6.0*d2vdp2*d4vdrdsdp2[i][j];
            double d2adrds    = (  d4vdrdsdp2[i][j]*d3vdp3 + d3vdrdp2[i]*d4vdsdp3[j]
                           + d3vdsdp2[j]*d4vdrdp3[i] + d2vdp2*d5vdrdsdp3[i][j]
                           - d2v1Refdrds*d4vdp4/2.0 - dv1Refdr*d5vdsdp4[j]/2.0
            - dv1Refds*d5vdrdp4[i]/2.0 - v1Ref*d6vdrdsdp4[i][j]/2.0
        )/denom
                                                - (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdr*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)*ddenomds/(denom*denom)
            - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*d2denomdrds/(denom*denom)
      - (d3vdsdp2[j]*d3vdp3 + d2vdp2*d4vdsdp3[j] - dv1Refds*d4vdp4/2.0 - v1Ref*d5vdsdp4[j]/2.0)*ddenomdr/(denom*denom)
            + 2.0*(d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*ddenomdr*ddenomds/(denom*denom*denom);


            double d2bdrds    = (d4vdrdsdp2[i][j]*d4vdp4/4.0 + d3vdrdp2[i]*d5vdsdp4[j]/4.0
                           + d3vdsdp2[j]*d5vdrdp4[i]/4.0 + d2vdp2*d6vdrdsdp4[i][j]/4.0
            - 2.0/3.0*d4vdsdp3[j]*d4vdrdp3[i] - 2.0/3.0*d3vdp3*d5vdrdsdp3[i][j]
                           )/denom
                                                - (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])*ddenomds/(denom*denom)
                        - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)*d2denomdrds/(denom*denom)
      - (d3vdsdp2[j]*d4vdp4/4.0+d2vdp2*d5vdsdp4[j]/4.0-d4vdsdp3[j]*d3vdp3/3.0-d3vdp3*d4vdsdp3[j]/3.0)*ddenomdr/(denom*denom)
      + 2.0*(d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)*ddenomdr*ddenomds/(denom*denom*denom);

            double d2gIntdrds = 0.0;

            if ((a == 0.0) && (b == 0.0)) {
            d2gIntdrds = d2v0drds*(p-pr) + d2v1drds*(p-pr)*(p-pr)/2.0 + d2v2drds*(p-pr)*(p-pr)*(p-pr)/6.0;

            } else if ((a != 0.0) && (b == 0.0)) {
            d2gIntdrds = (d2v0drds - d2v2drds/(2.0*a*a) + dv2dr*dads/(a*a*a) + (dv2ds*dadr + v2*d2adrds)/(a*a*a) - 3.0*v2*dadr*dads/(a*a*a*a)
         )*(p-pr)
                        + (d2v1drds + d2v2drds/(2.0*a) - dv2dr*dads/(2.0*a*a) - (dv2ds*dadr + v2*d2adrds)/(2.0*a*a) + v2*dadr*dads/(a*a*a)
         )*(p-pr)*(p-pr)/2.0;
            d2gIntdrds += (d2v2drds*log(1.0+a*(p-pr)) + dv2dr*dads*(p-pr)/(1.0+a*(p-pr))
                                + (dv2ds*dadr + v2*d2adrds)*(p-pr)/(1.0+a*(p-pr))
           - v2*dadr*dads*(p-pr)*(p-pr)/((1.0+a*(p-pr))*(1.0+a*(p-pr)))
                )/(2.0*a*a*a)
            - 3.0*(dv2dr*log(1.0+a*(p-pr)) + v2*dadr*(p-pr)/(1.0+a*(p-pr)))*dads/(2.0*a*a*a*a)
              - 3.0*(d2adrds*v2*log(1.0+a*(p-pr)) + dadr*dv2ds*log(1.0+a*(p-pr)) + dadr*v2*dads*(p-pr)/(1.0+a*(p-pr)))/(2.0*a*a*a*a)
            + 12.0*dadr*v2*log(1.0+a*(p-pr))*dads/(2.0*a*a*a*a*a);

            } else if ((a == 0.0) && (b != 0.0)) {
            printf("*-->Exception in d2IntegralV_GKdrds (liquid.c). a is zero, b is greater than zero.\n"); liqERRstate = ERR_A_ZERO;
            d2gIntdrds = 0.0;

            } else if (sum > 0.0) {
            d2gIntdrds = d2gdr2GMAP(p/10000.0, pr/10000.0,
                          v0,             v1*10000.0,       v2*10000.0*10000.0,       a*10000.0,       b*10000.0*10000.0,
                          dv0dr,       dv1dr*10000.0,    dv2dr*10000.0*10000.0,    dadr*10000.0,    dbdr*10000.0*10000.0,
                                                                dv0ds,       dv1ds*10000.0,    dv2ds*10000.0*10000.0,    dads*10000.0,    dbds*10000.0*10000.0,
              d2v0drds, d2v1drds*10000.0, d2v2drds*10000.0*10000.0, d2adrds*10000.0, d2bdrds*10000.0*10000.0)*10000.0;

            } else if (sum == 0.0) {
            printf("*-->Exception in d2IntegralV_GKdrds (liquid.c). a*a-4*b is equal to zero.\n"); liqERRstate = ERR_SUM_ZERO;
            d2gIntdrds = 0.0;

            } else if(sum < 0.0) {
            d2gIntdrds = d2gdr2LMAP(p/10000.0, pr/10000.0,
                          v0,             v1*10000.0,       v2*10000.0*10000.0,       a*10000.0,       b*10000.0*10000.0,
                          dv0dr,       dv1dr*10000.0,    dv2dr*10000.0*10000.0,    dadr*10000.0,    dbdr*10000.0*10000.0,
                                                                dv0ds,       dv1ds*10000.0,    dv2ds*10000.0*10000.0,    dads*10000.0,    dbds*10000.0*10000.0,
              d2v0drds, d2v1drds*10000.0, d2v2drds*10000.0*10000.0, d2adrds*10000.0, d2bdrds*10000.0*10000.0)*10000.0;

            }

            drds[i][j] = d2gIntdrds;
        } /* end loop on j [NS] */
    } /* end loop on i [NR] */
}

static void d2IntegralV_GKdrdt(double r[NR], double s[NT], double t, double p, double *drdt) {
    const double pr     = 1.0;
    const double tr     = 1673.15;
    double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum;
    double dvdr[NR], d2vdrdt[NR], dcRefdr[NR], d2cdrdt[NR], dmwdr[NR], dcpdr[NR], d3vdrdp2[NR], d4vdrdp3[NR], d5vdrdp4[NR], dmOxTotdr[NR], denom;
    int i, j;

    for (i=0; i<NR; i++) drdt[i] = 0.0;
    if (fabs(p-pr) < 10.0*DBL_EPSILON) return;

    /* Convert input composition (r) to liquid moles (m)  */
    for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* Compute moles and total moles of oxides */
    for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
    if (mOxTot == 0.0) return;

    /* Deal with the special case of FeO1.3 */
    mOx[NA] = 0.0;
    if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
        const double y = 0.3;
        mOx[iOxFeO1_3] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
        }
    }

    for (i=0; i<NR; i++) {
        dvdr[i]    = 0.0; d2vdrdt[i]  = 0.0; dcRefdr[i]  = 0.0; d2cdrdt[i]  = 0.0; dmwdr[i] = 0.0;
        dcpdr[i]   = 0.0; d3vdrdp2[i] = 0.0; d4vdrdp3[i] = 0.0; d5vdrdp4[i] = 0.0;
        for (j=0, dmOxTotdr[i]=0.0; j<nc; j++) dmOxTotdr[i] += (liquid[i+1].liqToOx)[j] - (liquid[0].liqToOx)[j];
    }

    for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
        v        += mOx[i]*bulkSystem[i].gk_v;
        dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
        cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
        dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
        cp      += mOx[i]*bulkSystem[i].gk_cp;
        d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
        d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
        d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        mw      += mOx[i]*bulkSystem[i].mw;

        for (j=0; j<NR; j++) {
            dvdr[j]       += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_v;
            d2vdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dvdt;
            dmwdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].mw;
            dcpdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_cp;
            dcRefdr[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c/mOxTot
              - mOx[i]*bulkSystem[i].gk_c*dmOxTotdr[j]/(mOxTot*mOxTot);
            dcRefdr[j]   += (iOxAl2O3 != -1) ? ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
              - 2.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot) : 0.0;
            d2cdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt/mOxTot
              - mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotdr[j]/(mOxTot*mOxTot);
            d3vdrdp2[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d4vdrdp3[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d5vdrdp4[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        }
        if (iCmpAl2O3 != -1) dcRefdr[iCmpAl2O3] += ((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*mOx[i]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot);
    }
    if (v == 0.0) return;

    alpha   = dvdt/v;
    v0      = v*exp(alpha*(t-tr));
    v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
    c      = cRef + (t-tr)*dcdt;
    v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
    v2      = d2vdp2;
    denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
    a      = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
    b      = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;
    sum      = a*a - 4.0*b;

    for (i=0; i<NR; i++) {
        double dalphadr = d2vdrdt[i]/v - dvdt*dvdr[i]/(v*v);
        double dv0dr    = (dvdr[i] + v*dalphadr*(t-tr))*exp(alpha*(t-tr));
        double dv1Refdr = -2.0*v*dvdr[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                    -v*v*(- 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
                            + 2.0*tr*alpha*dalphadr/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp));
        double dcdr     = dcRefdr[i] + (t-tr)*d2cdrdt[i];
        double dv1dr    = -2.0*v0*dv0dr*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                    -v0*v0*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdr)/(mw*mw*c*c*c)
                                + 2.0*t*alpha*dalphadr/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp));
        double dv2dr    = d3vdrdp2[i];
        double ddenomdr = 2.0*dv1Refdr*d3vdp3+2.0*v1Ref*d4vdrdp3[i]-6.0*d2vdp2*d3vdrdp2[i];
        double dadr     = (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdr*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)/denom
                - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdr;
        double dbdr     = (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])/denom
                - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdr;
        double dv0dt    = v*alpha*exp(alpha*(t-tr));
        double dv1dt    = - 2.0*v0*v0*alpha*(1000.0/(mw*c*c) + t*alpha*alpha/(cp)) - v0*v0*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp));

        double d2v0drdt = v*dalphadr*exp(alpha*(t-tr)) + (dvdr[i] + v*dalphadr*(t-tr))*alpha*exp(alpha*(t-tr));
        double d2v1drdt = -2.0*(dv0dt*dv0dr + v0*d2v0drdt)*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                - 2.0*v0*dv0dr*(- 2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp))
                - 2.0*v0*dv0dt*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdr)/(mw*mw*c*c*c)
                  + 2.0*t*alpha*dalphadr/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp))
              - v0*v0*(- 1000.0*(dcdt*dmwdr[i]+2.0*mw*d2cdrdt[i])/(mw*mw*c*c*c) + 3000.0*(c*dmwdr[i]+2.0*mw*dcdr)*dcdt/(mw*mw*c*c*c*c)
                 + 2.0*alpha*dalphadr/(cp) - alpha*alpha*dcpdr[i]/(cp*cp));
        double d2gIntdrdt = 0.0;

        if ((a == 0.0) && (b == 0.0)) {
            d2gIntdrdt       = d2v0drdt*(p-pr) + d2v1drdt*(p-pr)*(p-pr)/2.0;

        } else if ((a != 0.0) && (b == 0.0)) {
            printf("*-->Exception in d2IntegralV_GKdrdt (liquid.c). a is greater than zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
            d2gIntdrdt       = 0.0;

        } else if ((a == 0.0) && (b != 0.0)) {
            printf("*-->Exception in d2IntegralV_GKdrdt (liquid.c). a is zero, b is greater than zero.\n"); liqERRstate = ERR_A_ZERO;
            d2gIntdrdt       = 0.0;

        } else if (sum > 0.0) {
            d2gIntdrdt       = d2gdrdtGMAP(p/10000.0, pr/10000.0,
                                     v0,             v1*10000.0,    v2*10000.0*10000.0,    a*10000.0,    b*10000.0*10000.0,
             dv0dr,       dv1dr*10000.0, dv2dr*10000.0*10000.0, dadr*10000.0, dbdr*10000.0*10000.0,
             dv0dt,       dv1dt*10000.0,
             d2v0drdt, d2v1drdt*10000.0)*10000.0;

        } else if (sum == 0.0) {
            printf("*-->Exception in d2IntegralV_GKdrdt (liquid.c). a*a-4*b is equal to zero.\n"); liqERRstate = ERR_SUM_ZERO;
            d2gIntdrdt       = 0.0;

        } else if(sum < 0.0) {
            d2gIntdrdt       = d2gdrdtLMAP(p/10000.0, pr/10000.0,
                                     v0,             v1*10000.0,    v2*10000.0*10000.0,    a*10000.0,    b*10000.0*10000.0,
             dv0dr,       dv1dr*10000.0, dv2dr*10000.0*10000.0, dadr*10000.0, dbdr*10000.0*10000.0,
             dv0dt,       dv1dt*10000.0,
             d2v0drdt, d2v1drdt*10000.0)*10000.0;

        }

        drdt[i] += d2gIntdrdt;
    }
}


static void d2IntegralV_GKdrdp(double r[NR], double s[NT], double t, double p, double *drdp) {
    const double pr     = 1.0;
    const double tr     = 1673.15;
    double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum;
    double dvdr[NR], d2vdrdt[NR], dcRefdr[NR], d2cdrdt[NR], dmwdr[NR], dcpdr[NR], d3vdrdp2[NR], d4vdrdp3[NR], d5vdrdp4[NR], dmOxTotdr[NR], denom;
    int i, j;

    for (i=0; i<NR; i++) drdp[i] = 0.0;

    /* Convert input composition (r) to liquid moles (m)  */
    for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* Compute moles and total moles of oxides */
    for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
    if (mOxTot == 0.0) return;

    /* Deal with the special case of FeO1.3 */
    mOx[NA] = 0.0;
    if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
        const double y = 0.3;
        mOx[iOxFeO1_3] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
        }
    }

    for (i=0; i<NR; i++) {
        dvdr[i]    = 0.0; d2vdrdt[i]  = 0.0; dcRefdr[i]  = 0.0; d2cdrdt[i]  = 0.0; dmwdr[i] = 0.0;
        dcpdr[i]   = 0.0; d3vdrdp2[i] = 0.0; d4vdrdp3[i] = 0.0; d5vdrdp4[i] = 0.0;
        for (j=0, dmOxTotdr[i]=0.0; j<nc; j++) dmOxTotdr[i] += (liquid[i+1].liqToOx)[j] - (liquid[0].liqToOx)[j];
    }

    for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
        v        += mOx[i]*bulkSystem[i].gk_v;
        dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
        cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
        dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
        cp      += mOx[i]*bulkSystem[i].gk_cp;
        d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
        d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
        d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        mw      += mOx[i]*bulkSystem[i].mw;

        for (j=0; j<NR; j++) {
            dvdr[j]       += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_v;
            d2vdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dvdt;
            dmwdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].mw;
            dcpdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_cp;
            dcRefdr[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c/mOxTot
              - mOx[i]*bulkSystem[i].gk_c*dmOxTotdr[j]/(mOxTot*mOxTot);
            dcRefdr[j]   += (iOxAl2O3 != -1) ? ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
              - 2.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot) : 0.0;
            d2cdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt/mOxTot
              - mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotdr[j]/(mOxTot*mOxTot);
            d3vdrdp2[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d4vdrdp3[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d5vdrdp4[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        }
        if (iCmpAl2O3 != -1) dcRefdr[iCmpAl2O3] += ((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*mOx[i]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot);
    }
    if (v == 0.0) return;

    alpha   = dvdt/v;
    v0      = v*exp(alpha*(t-tr));
    v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
    c      = cRef + (t-tr)*dcdt;
    v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
    v2      = d2vdp2;
    denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
    a      = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
    b      = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;
    sum      = a*a - 4.0*b;

    for (i=0; i<NR; i++) {
        double dalphadr   = d2vdrdt[i]/v - dvdt*dvdr[i]/(v*v);
        double dv0dr      = (dvdr[i] + v*dalphadr*(t-tr))*exp(alpha*(t-tr));
        double dv1Refdr   = -2.0*v*dvdr[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
        -v*v*(- 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
                                + 2.0*tr*alpha*dalphadr/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp));
        double dcdr       = dcRefdr[i] + (t-tr)*d2cdrdt[i];
        double dv1dr      = -2.0*v0*dv0dr*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
        -v0*v0*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdr)/(mw*mw*c*c*c)
              + 2.0*t*alpha*dalphadr/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp));
        double dv2dr      = d3vdrdp2[i];
        double ddenomdr   = 2.0*dv1Refdr*d3vdp3+2.0*v1Ref*d4vdrdp3[i]-6.0*d2vdp2*d3vdrdp2[i];
        double dadr       = (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdr*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)/denom
                    - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdr;
        double dbdr       = (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])/denom
                    - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdr;
        double d2gdrdp;

        if ((a == 0.0) && (b == 0.0)) {
            d2gdrdp  = dv1dr + dv2dr*(p-pr);

        } else if ((a != 0.0) && (b == 0.0)) {
            printf("*-->Exception in d2IntegralV_GKdrdp (liquid.c). a is greater than zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
            d2gdrdp  = 0.0;

        } else if ((a == 0.0) && (b != 0.0)) {
            printf("*-->Exception in d2IntegralV_GKdrdp (liquid.c). a is zero, b is greater than zero.\n"); liqERRstate = ERR_A_ZERO;
            d2gdrdp  = 0.0;

        } else {
            d2gdrdp  = d2gdrdpMAP(p, pr, v0, v1, v2, a, b, dv0dr, dv1dr, dv2dr, dadr, dbdr);

        }

        drdp[i] += d2gdrdp;
    }
}

static void d2IntegralV_GKds2(double r[NR], double s[NT], double t, double p, double d2s[NS][NS]) {
    const double pr = 1.0;
    const double tr = 1673.15;
    const double y  = 0.3;
    double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum;
    double dvds[NS], d2vdsdt[NS], dcRefds[NS], d2cdsdt[NS], dmwds[NS], dcpds[NS], d3vdsdp2[NS], d4vdsdp3[NS], d5vdsdp4[NS], denom;
    double d2vds2[NS][NS], d3vds2dt[NS][NS], d2cRefds2[NS][NS], d3cds2dt[NS][NS], d2mwds2[NS][NS], d2cpds2[NS][NS], d4vds2dp2[NS][NS],
         d5vds2dp3[NS][NS], d6vds2dp4[NS][NS];
    int i, j, k;

    for (i=0; i<NS; i++) for (j=0; j<NS; j++) d2s[i][j] = 0.0;
    if (fabs(p-pr) < 10.0*DBL_EPSILON) return;

    /* Convert input composition (r) to liquid moles (m)  */
    for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* Compute moles and total moles of oxides */
    for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
    if (mOxTot == 0.0) return;

    /* Deal with the special case of FeO1.3 */
    mOx[NA] = 0.0;
    if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
        const double y = 0.3;
        mOx[iOxFeO1_3] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
        }
    }

    for (i=0; i<NS; i++) {
        dvds[i]    = 0.0; d2vdsdt[i]  = 0.0; dcRefds[i]  = 0.0; d2cdsdt[i]  = 0.0; dmwds[i] = 0.0;
        dcpds[i]   = 0.0; d3vdsdp2[i] = 0.0; d4vdsdp3[i] = 0.0; d5vdsdp4[i] = 0.0;
        for (j=0; j<NS; j++) {
            d2vds2[i][j]  = 0.0; d3vds2dt[i][j]  = 0.0; d2cRefds2[i][j] = 0.0; d3cds2dt[i][j]  = 0.0; d2mwds2[i][j] = 0.0;
            d2cpds2[i][j] = 0.0; d4vds2dp2[i][j] = 0.0; d5vds2dp3[i][j] = 0.0; d6vds2dp4[i][j] = 0.0;
        }
    }

    for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
        double dmOxds[NS], dmOxTotds[NS], d2mOxds2[NS][NS], d2mOxTotds2[NS][NS];
        v        += mOx[i]*bulkSystem[i].gk_v;
        dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
        cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
        dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
        cp      += mOx[i]*bulkSystem[i].gk_cp;
        d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
        d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
        d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        mw      += mOx[i]*bulkSystem[i].mw;

        for (j=0; j<NS; j++) { dmOxds[j] = 0.0; dmOxTotds[j] = 0.0; }
        for (j=0; j<NS; j++) for (k=0; k<NS; k++) { d2mOxds2[j][k] = 0.0; d2mOxTotds2[j][k] = 0.0; }
        if (iCmpFe2SiO4_6 != -1) {

            for (j=0; j<NS; j++) {
                if      (iOxFe2O3  == i) dmOxds[j] += -2.0*y*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
                else if (iOxFeO       == i) dmOxds[j] += -2.0*(1.0-2.0*y)*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
                else if (iOxFeO1_3 == i) dmOxds[j] += 2.0*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
  dmOxTotds[j] += 2.0*y*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
  for (k=0; k<NS; k++) {
                    if      (iOxFe2O3  == i) d2mOxds2[j][k] += -2.0*y*s[iCmpFe2SiO4_6]*d2nSpeciesds2[j][k];
                    else if (iOxFeO    == i) d2mOxds2[j][k] += -2.0*(1.0-2.0*y)*s[iCmpFe2SiO4_6]*d2nSpeciesds2[j][k];
                    else if (iOxFeO1_3 == i) d2mOxds2[j][k] += 2.0*s[iCmpFe2SiO4_6]*d2nSpeciesds2[j][k];
    d2mOxTotds2[j][k] += 2.0*y*s[iCmpFe2SiO4_6]*d2nSpeciesds2[j][k];
  }
            }

            if      (iOxFe2O3  == i) dmOxds[iCmpFe2SiO4_6] += -2.0*y*nSpecies;
            else if (iOxFeO     == i) dmOxds[iCmpFe2SiO4_6] += -2.0*(1.0-2.0*y)*nSpecies;
            else if (iOxFeO1_3 == i) dmOxds[iCmpFe2SiO4_6] += 2.0*nSpecies;
            dmOxTotds[iCmpFe2SiO4_6] += 2.0*y*nSpecies;
            for (j=0; j<NS; j++) {
                if (j != iCmpFe2SiO4_6) {
                    if      (iOxFe2O3  == i) { d2mOxds2[iCmpFe2SiO4_6][j] += -2.0*y*dnSpeciesds[j];        d2mOxds2[j][iCmpFe2SiO4_6] += -2.0*y*dnSpeciesds[j];       }
                    else if (iOxFeO    == i) { d2mOxds2[iCmpFe2SiO4_6][j] += -2.0*(1.0-2.0*y)*dnSpeciesds[j]; d2mOxds2[j][iCmpFe2SiO4_6] += -2.0*(1.0-2.0*y)*dnSpeciesds[j]; }
                    else if (iOxFeO1_3 == i) { d2mOxds2[iCmpFe2SiO4_6][j] += 2.0*dnSpeciesds[j];            d2mOxds2[j][iCmpFe2SiO4_6] += 2.0*dnSpeciesds[j];              }
                    d2mOxTotds2[iCmpFe2SiO4_6][j] += 2.0*y*dnSpeciesds[j];
                    d2mOxTotds2[j][iCmpFe2SiO4_6] += 2.0*y*dnSpeciesds[j];
  }
            }
            if      (iOxFe2O3  == i) d2mOxds2[iCmpFe2SiO4_6][iCmpFe2SiO4_6] += -4.0*y*dnSpeciesds[iCmpFe2SiO4_6];
            else if (iOxFeO     == i) d2mOxds2[iCmpFe2SiO4_6][iCmpFe2SiO4_6] += -4.0*(1.0-2.0*y)*dnSpeciesds[iCmpFe2SiO4_6];
            else if (iOxFeO1_3 == i) d2mOxds2[iCmpFe2SiO4_6][iCmpFe2SiO4_6] += 4.0*dnSpeciesds[iCmpFe2SiO4_6];
            d2mOxTotds2[iCmpFe2SiO4_6][iCmpFe2SiO4_6] += 4.0*y*dnSpeciesds[iCmpFe2SiO4_6];

        }
        if (iCmpFe2AlO4_1 != -1) {

            for (j=0; j<NS; j++) {
                if      (iOxFe2O3  == i) dmOxds[j] += -2.0*y*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
                else if (iOxFeO       == i) dmOxds[j] += -2.0*(1.0-2.0*y)*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
                else if (iOxFeO1_3 == i) dmOxds[j] += 2.0*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
  dmOxTotds[j] += 2.0*y*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
  for (k=0; k<NS; k++) {
                    if      (iOxFe2O3  == i) d2mOxds2[j][k] += -2.0*y*s[iCmpFe2AlO4_1]*d2nSpeciesds2[j][k];
                    else if (iOxFeO    == i) d2mOxds2[j][k] += -2.0*(1.0-2.0*y)*s[iCmpFe2AlO4_1]*d2nSpeciesds2[j][k];
                    else if (iOxFeO1_3 == i) d2mOxds2[j][k] += 2.0*s[iCmpFe2AlO4_1]*d2nSpeciesds2[j][k];
    d2mOxTotds2[j][k] += 2.0*y*s[iCmpFe2AlO4_1]*d2nSpeciesds2[j][k];
  }
            }

            if      (iOxFe2O3  == i) dmOxds[iCmpFe2AlO4_1] += -2.0*y*nSpecies;
            else if (iOxFeO     == i) dmOxds[iCmpFe2AlO4_1] += -2.0*(1.0-2.0*y)*nSpecies;
            else if (iOxFeO1_3 == i) dmOxds[iCmpFe2AlO4_1] += 2.0*nSpecies;
            dmOxTotds[iCmpFe2AlO4_1] += 2.0*y*nSpecies;
            for (j=0; j<NS; j++) {
                if (j != iCmpFe2AlO4_1) {
                    if      (iOxFe2O3  == i) { d2mOxds2[iCmpFe2AlO4_1][j] += -2.0*y*dnSpeciesds[j];        d2mOxds2[j][iCmpFe2AlO4_1] += -2.0*y*dnSpeciesds[j];       }
                    else if (iOxFeO    == i) { d2mOxds2[iCmpFe2AlO4_1][j] += -2.0*(1.0-2.0*y)*dnSpeciesds[j]; d2mOxds2[j][iCmpFe2AlO4_1] += -2.0*(1.0-2.0*y)*dnSpeciesds[j]; }
                    else if (iOxFeO1_3 == i) { d2mOxds2[iCmpFe2AlO4_1][j] += 2.0*dnSpeciesds[j];            d2mOxds2[j][iCmpFe2AlO4_1] += 2.0*dnSpeciesds[j];              }
                    d2mOxTotds2[iCmpFe2AlO4_1][j] += 2.0*y*dnSpeciesds[j];
                    d2mOxTotds2[j][iCmpFe2AlO4_1] += 2.0*y*dnSpeciesds[j];
  }
            }
            if      (iOxFe2O3  == i) d2mOxds2[iCmpFe2AlO4_1][iCmpFe2AlO4_1] += -4.0*y*dnSpeciesds[iCmpFe2AlO4_1];
            else if (iOxFeO     == i) d2mOxds2[iCmpFe2AlO4_1][iCmpFe2AlO4_1] += -4.0*(1.0-2.0*y)*dnSpeciesds[iCmpFe2AlO4_1];
            else if (iOxFeO1_3 == i) d2mOxds2[iCmpFe2AlO4_1][iCmpFe2AlO4_1] += 4.0*dnSpeciesds[iCmpFe2AlO4_1];
            d2mOxTotds2[iCmpFe2AlO4_1][iCmpFe2AlO4_1] += 4.0*y*dnSpeciesds[iCmpFe2AlO4_1];

        }
        for (j=0; j<NS; j++) {
            dvds[j]       += dmOxds[j]*bulkSystem[i].gk_v;
            d2vdsdt[j]   += dmOxds[j]*bulkSystem[i].gk_dvdt;
            dmwds[j]     += dmOxds[j]*bulkSystem[i].mw;
            dcpds[j]     += dmOxds[j]*bulkSystem[i].gk_cp;
            dcRefds[j]   += dmOxds[j]*bulkSystem[i].gk_c/mOxTot - mOx[i]*bulkSystem[i].gk_c*dmOxTotds[j]/(mOxTot*mOxTot);
            dcRefds[j]   += (iOxAl2O3 != -1) ? dmOxds[j]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
              - 2.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotds[j]/(mOxTot*mOxTot*mOxTot) : 0.0;
            d2cdsdt[j]   += dmOxds[j]*bulkSystem[i].gk_dcdt/mOxTot - mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotds[j]/(mOxTot*mOxTot);
            d3vdsdp2[j]  += dmOxds[j]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d4vdsdp3[j]  += dmOxds[j]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d5vdsdp4[j]  += dmOxds[j]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);

            for (k=0; k<NS; k++) {
                d2vds2[j][k]     += d2mOxds2[j][k]*bulkSystem[i].gk_v;
                d3vds2dt[j][k]   += d2mOxds2[j][k]*bulkSystem[i].gk_dvdt;
                d2mwds2[j][k]    += d2mOxds2[j][k]*bulkSystem[i].mw;
                d2cpds2[j][k]    += d2mOxds2[j][k]*bulkSystem[i].gk_cp;
                d2cRefds2[j][k]  += d2mOxds2[j][k]*bulkSystem[i].gk_c/mOxTot
                    - dmOxds[j]*bulkSystem[i].gk_c*dmOxTotds[k]/(mOxTot*mOxTot)
                    - dmOxds[k]*bulkSystem[i].gk_c*dmOxTotds[j]/(mOxTot*mOxTot)
        - mOx[i]*bulkSystem[i].gk_c*d2mOxTotds2[j][k]/(mOxTot*mOxTot)
        + 2.0*mOx[i]*bulkSystem[i].gk_c*dmOxTotds[j]*dmOxTotds[k]/(mOxTot*mOxTot*mOxTot);
                d2cRefds2[j][k]  += (iOxAl2O3 != -1) ?
                      d2mOxds2[j][k]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
        - 2.0*dmOxds[j]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotds[k]/(mOxTot*mOxTot*mOxTot)
                    - 2.0*dmOxds[k]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotds[j]/(mOxTot*mOxTot*mOxTot)
        - 2.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*d2mOxTotds2[j][k]/(mOxTot*mOxTot*mOxTot)
        + 6.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotds[j]*dmOxTotds[k]/(mOxTot*mOxTot*mOxTot*mOxTot)
        : 0.0;
                d3cds2dt[j][k]   += d2mOxds2[j][k]*bulkSystem[i].gk_dcdt/mOxTot
                    - dmOxds[j]*bulkSystem[i].gk_dcdt*dmOxTotds[k]/(mOxTot*mOxTot)
                    - dmOxds[k]*bulkSystem[i].gk_dcdt*dmOxTotds[j]/(mOxTot*mOxTot)
        - mOx[i]*bulkSystem[i].gk_dcdt*d2mOxTotds2[j][k]/(mOxTot*mOxTot)
        + 2.0*mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotds[j]*dmOxTotds[k]/(mOxTot*mOxTot*mOxTot);
                d4vds2dp2[j][k]  += d2mOxds2[j][k]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
                d5vds2dp3[j][k]  += d2mOxds2[j][k]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
                d6vds2dp4[j][k]  += d2mOxds2[j][k]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
            }

        }
    }
    if (v == 0.0) return;

    alpha   = dvdt/v;
    v0      = v*exp(alpha*(t-tr));
    v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
    c      = cRef + (t-tr)*dcdt;
    v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
    v2      = d2vdp2;
    denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
    a      = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
    b      = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;
    sum      = a*a - 4.0*b;

    for (i=0; i<NS; i++) {
        double dalphadsi  = d2vdsdt[i]/v - dvdt*dvds[i]/(v*v);
        double dv0dsi     = (dvds[i] + v*dalphadsi*(t-tr))*exp(alpha*(t-tr));
        double dv1Refdsi  = -2.0*v*dvds[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
             -v*v*(-1000.0*(cRef*dmwds[i]+2.0*mw*dcRefds[i])/(mw*mw*cRef*cRef*cRef)
                                            + 2.0*tr*alpha*dalphadsi/(cp) - tr*alpha*alpha*dcpds[i]/(cp*cp));
        double dcdsi      = dcRefds[i] + (t-tr)*d2cdsdt[i];
        double dv1dsi     = -2.0*v0*dv0dsi*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
             -v0*v0*(-1000.0*(c*dmwds[i]+2.0*mw*dcdsi)/(mw*mw*c*c*c)
                                                + 2.0*t*alpha*dalphadsi/(cp) - t*alpha*alpha*dcpds[i]/(cp*cp));
        double ddenomdsi  = 2.0*dv1Refdsi*d3vdp3+2.0*v1Ref*d4vdsdp3[i]-6.0*d2vdp2*d3vdsdp2[i];
        double dadsi      = (d3vdsdp2[i]*d3vdp3+d2vdp2*d4vdsdp3[i]-dv1Refdsi*d4vdp4/2.0-v1Ref*d5vdsdp4[i]/2.0)/denom
                    - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdsi;
        double dbdsi      = (d3vdsdp2[i]*d4vdp4/4.0+d2vdp2*d5vdsdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdsdp3[i])/denom
                    - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdsi;

        for (j=i; j<NS; j++) {
            double dalphadsj   = d2vdsdt[j]/v - dvdt*dvds[j]/(v*v);
            double dv0dsj     = (dvds[j] + v*dalphadsj*(t-tr))*exp(alpha*(t-tr));
            double dv1Refdsj   = -2.0*v*dvds[j]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                -v*v*(- 1000.0*(cRef*dmwds[j]+2.0*mw*dcRefds[j])/(mw*mw*cRef*cRef*cRef)
           + 2.0*tr*alpha*dalphadsj/(cp) - tr*alpha*alpha*dcpds[j]/(cp*cp));
            double dcdsj     = dcRefds[j] + (t-tr)*d2cdsdt[j];
            double dv1dsj     = -2.0*v0*dv0dsj*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                -v0*v0*(- 1000.0*(c*dmwds[j]+2.0*mw*dcdsj)/(mw*mw*c*c*c)
             + 2.0*t*alpha*dalphadsj/(cp) - t*alpha*alpha*dcpds[j]/(cp*cp));
            double ddenomdsj   = 2.0*dv1Refdsj*d3vdp3+2.0*v1Ref*d4vdsdp3[j]-6.0*d2vdp2*d3vdsdp2[j];
            double dadsj     = (d3vdsdp2[j]*d3vdp3+d2vdp2*d4vdsdp3[j]-dv1Refdsj*d4vdp4/2.0-v1Ref*d5vdsdp4[j]/2.0)/denom
            - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdsj;
            double dbdsj     = (d3vdsdp2[j]*d4vdp4/4.0+d2vdp2*d5vdsdp4[j]/4.0-2.0/3.0*d3vdp3*d4vdsdp3[j])/denom
               - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdsj;

            double d2alphads2   = d3vds2dt[i][j]/v - d2vdsdt[i]*dvds[j]/(v*v) - d2vdsdt[j]*dvds[i]/(v*v) - dvdt*d2vds2[i][j]/(v*v) + 2.0*dvdt*dvds[i]*dvds[j]/(v*v*v);
            double d2v0ds2      = d2vds2[i][j]*exp(alpha*(t-tr))
                                                    + dvds[i]*dalphadsj*(t-tr)*exp(alpha*(t-tr))
                        + dvds[j]*dalphadsi*(t-tr)*exp(alpha*(t-tr))
                        + v*d2alphads2*(t-tr)*exp(alpha*(t-tr))
                        + v*dalphadsi*pow(t-tr,2.0)*dalphadsj*exp(alpha*(t-tr));
            double d2v1Refds2   = -2.0*v*d2vds2[i][j]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                                                    - 2.0*dvds[j]*dvds[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                        - 2.0*v*dvds[i]*(- 1000.0*(cRef*dmwds[j]+2.0*mw*dcRefds[j])/(mw*mw*cRef*cRef*cRef)
                   + 2.0*tr*alpha*dalphadsj/(cp) - tr*alpha*alpha*dcpds[j]/(cp*cp))
          - 2.0*v*dvds[j]*(- 1000.0*(cRef*dmwds[i]+2.0*mw*dcRefds[i])/(mw*mw*cRef*cRef*cRef)
                    + 2.0*tr*alpha*dalphadsi/(cp) - tr*alpha*alpha*dcpds[i]/(cp*cp))
               -v*v*(- 1000.0*(dcRefds[j]*dmwds[i]+cRef*d2mwds2[i][j]+2.0*dmwds[j]*dcRefds[i]+2.0*mw*d2cRefds2[i][j])/(mw*mw*cRef*cRef*cRef)
                        + 1000.0*(cRef*dmwds[i]+2.0*mw*dcRefds[i])*(2.0*dmwds[j]*cRef+3.0*mw*dcRefds[j])/(mw*mw*mw*cRef*cRef*cRef*cRef)
                        + 2.0*tr*alpha*d2alphads2/(cp)
                        + 2.0*tr*dalphadsj*dalphadsi/(cp) - 2.0*tr*alpha*dalphadsi*dcpds[j]/(cp*cp)
                        - 2.0*tr*alpha*dalphadsj*dcpds[i]/(cp*cp) + 2.0*tr*alpha*alpha*dcpds[i]*dcpds[j]/(cp*cp*cp));
            double d2cds2     = d2cRefds2[i][j] + (t-tr)*d3cds2dt[i][j];
            double d2v1ds2     = -2.0*dv0dsj*dv0dsi*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                        - 2.0*v0*d2v0ds2*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                        - 2.0*v0*dv0dsi*(- 1000.0*(c*dmwds[j]+2.0*mw*dcdsj)/(mw*mw*c*c*c)
                   + 2.0*t*alpha*dalphadsj/(cp) - t*alpha*alpha*dcpds[j]/(cp*cp))
          - 2.0*v0*dv0dsj*(- 1000.0*(c*dmwds[i]+2.0*mw*dcdsi)/(mw*mw*c*c*c)
                    + 2.0*t*alpha*dalphadsi/(cp) - t*alpha*alpha*dcpds[i]/(cp*cp))
               -v0*v0*(- 1000.0*(dcdsj*dmwds[i]+c*d2mwds2[i][j]+2.0*dmwds[j]*dcdsi+2.0*mw*d2cds2)/(mw*mw*c*c*c)
                            + 1000.0*(c*dmwds[i]+2.0*mw*dcdsi)*(2.0*dmwds[j]*c+3.0*mw*dcdsj)/(mw*mw*mw*c*c*c*c)
                            + 2.0*t*alpha*d2alphads2/(cp)
                            + 2.0*t*dalphadsj*dalphadsi/(cp) - 2.0*t*alpha*dalphadsi*dcpds[j]/(cp*cp)
                            - 2.0*t*alpha*dalphadsj*dcpds[i]/(cp*cp) + 2.0*t*alpha*alpha*dcpds[i]*dcpds[j]/(cp*cp*cp));
            double d2denomds2   = 2.0*d2v1Refds2*d3vdp3+2.0*dv1Refdsi*d4vdsdp3[j]+2.0*dv1Refdsj*d4vdsdp3[i]+2.0*v1Ref*d5vds2dp3[i][j]-6.0*d3vdsdp2[j]*d3vdsdp2[i]-6.0*d2vdp2*d4vds2dp2[i][j];
            double d2ads2     = (d4vds2dp2[i][j]*d3vdp3+d3vdsdp2[i]*d4vdsdp3[j]+d3vdsdp2[j]*d4vdsdp3[i]+d2vdp2*d5vds2dp3[i][j]
                         - d2v1Refds2*d4vdp4/2.0-dv1Refdsi*d5vdsdp4[j]/2.0-dv1Refdsj*d5vdsdp4[i]/2.0-v1Ref*d6vds2dp4[i][j]/2.0)/denom
             - (d3vdsdp2[i]*d3vdp3+d2vdp2*d4vdsdp3[i]-dv1Refdsi*d4vdp4/2.0-v1Ref*d5vdsdp4[i]/2.0)*ddenomdsj/(denom*denom)
            - (d3vdsdp2[j]*d3vdp3+d2vdp2*d4vdsdp3[j]-dv1Refdsj*d4vdp4/2.0-v1Ref*d5vdsdp4[j]/2.0)*ddenomdsi/(denom*denom)
            + 2.0*(d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*ddenomdsi*ddenomdsj/(denom*denom*denom)
            - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*d2denomds2/(denom*denom);
            double d2bds2     = (d4vds2dp2[i][j]*d4vdp4/4.0+d3vdsdp2[i]*d5vdsdp4[j]/4.0+d3vdsdp2[j]*d5vdsdp4[i]/4.0+d2vdp2*d6vds2dp4[i][j]/4.0
                         - 2.0/3.0*d4vdsdp3[j]*d4vdsdp3[i]-2.0/3.0*d3vdp3*d5vds2dp3[i][j])/denom
            - (d3vdsdp2[i]*d4vdp4/4.0+d2vdp2*d5vdsdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdsdp3[i])*ddenomdsj/(denom*denom)
            - (d3vdsdp2[j]*d4vdp4/4.0+d2vdp2*d5vdsdp4[j]/4.0-2.0*d4vdsdp3[j]*d3vdp3/3.0)*ddenomdsi/(denom*denom)
             + 2.0*(d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)*ddenomdsi*ddenomdsj/(denom*denom*denom)
             - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)*d2denomds2/(denom*denom);
            double d2gIntds2   = 0.0;

            if ((a == 0.0) && (b == 0.0)) {
    d2gIntds2  = d2v0ds2*(p-pr) + d2v1ds2*(p-pr)*(p-pr)/2.0;

            } else if ((a != 0.0) && (b == 0.0)) {
    printf("*-->Exception in d2IntegralV_GKds2 (liquid.c). a is not equal to zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
    d2gIntds2  = 0.0;

            } else if ((a == 0.0) && (b != 0.0)) {
    printf("*-->Exception in d2IntegralV_GKds2 (liquid.c). a is zero, b is not equal to zero.\n"); liqERRstate = ERR_A_ZERO;
    d2gIntds2  = 0.0;

            } else if (sum > 0.0) {
    d2gIntds2  = d2gdr2GMAP(p/10000.0, pr/10000.0,
                          v0,           v1*10000.0,          d2vdp2*10000.0*10000.0,      a*10000.0,      b*10000.0*10000.0,
                          dv0dsi,   dv1dsi*10000.0,     d3vdsdp2[i]*10000.0*10000.0,  dadsi*10000.0,  dbdsi*10000.0*10000.0,
            dv0dsj,   dv1dsj*10000.0,     d3vdsdp2[j]*10000.0*10000.0,  dadsj*10000.0,  dbdsj*10000.0*10000.0,
              d2v0ds2, d2v1ds2*10000.0, d4vds2dp2[i][j]*10000.0*10000.0, d2ads2*10000.0, d2bds2*10000.0*10000.0)*10000.0;

            } else if (sum == 0.0) {
    printf("*-->Exception in d2IntegralV_GKds2 (liquid.c). a*a-4*b is equal to zero.\n"); liqERRstate = ERR_SUM_ZERO;
    d2gIntds2  = 0.0;

            } else if(sum < 0.0) {
    d2gIntds2  = d2gdr2LMAP(p/10000.0, pr/10000.0,
                          v0,           v1*10000.0,          d2vdp2*10000.0*10000.0,      a*10000.0,      b*10000.0*10000.0,
                          dv0dsi,   dv1dsi*10000.0,     d3vdsdp2[i]*10000.0*10000.0,  dadsi*10000.0,  dbdsi*10000.0*10000.0,
            dv0dsj,   dv1dsj*10000.0,     d3vdsdp2[j]*10000.0*10000.0,  dadsj*10000.0,  dbdsj*10000.0*10000.0,
              d2v0ds2, d2v1ds2*10000.0, d4vds2dp2[i][j]*10000.0*10000.0, d2ads2*10000.0, d2bds2*10000.0*10000.0)*10000.0;

            }

            d2s[i][j] += d2gIntds2;
            if (i != j) d2s[j][i] += d2gIntds2;

        }
    }
}

static void d2IntegralV_GKdsdt(double r[NR], double s[NT], double t, double p, double dsdt[NS]) {
    const double pr = 1.0;
    const double tr = 1673.15;
    const double y  = 0.3;
    double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum;
    double dvds[NS], d2vdsdt[NS], dcRefds[NS], d2cdsdt[NS], dmwds[NS], dcpds[NS], d3vdsdp2[NS], d4vdsdp3[NS], d5vdsdp4[NS], denom;
    int i, j;

    for (i=0; i<NS; i++) dsdt[i] = 0.0;
    if (fabs(p-pr) < 10.0*DBL_EPSILON) return;

    /* Convert input composition (r) to liquid moles (m)  */
    for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* Compute moles and total moles of oxides */
    for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
    if (mOxTot == 0.0) return;

    /* Deal with the special case of FeO1.3 */
    mOx[NA] = 0.0;
    if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
        const double y = 0.3;
        mOx[iOxFeO1_3] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
        }
    }

    for (i=0; i<NS; i++) {
        dvds[i]    = 0.0; d2vdsdt[i]  = 0.0; dcRefds[i]  = 0.0; d2cdsdt[i]  = 0.0; dmwds[i] = 0.0;
        dcpds[i]   = 0.0; d3vdsdp2[i] = 0.0; d4vdsdp3[i] = 0.0; d5vdsdp4[i] = 0.0;
    }

    for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
        double dmOxds[NS], dmOxTotds[NS];
        v        += mOx[i]*bulkSystem[i].gk_v;
        dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
        cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
        dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
        cp      += mOx[i]*bulkSystem[i].gk_cp;
        d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
        d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
        d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        mw      += mOx[i]*bulkSystem[i].mw;
        for (j=0; j<NS; j++) { dmOxds[j] = 0.0; dmOxTotds[j] = 0.0; }
        if (iCmpFe2SiO4_6 != -1) {
            for (j=0; j<NS; j++) {
                if      (iOxFe2O3  == i) dmOxds[j] += -2.0*y*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
                else if (iOxFeO       == i) dmOxds[j] += -2.0*(1.0-2.0*y)*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
                else if (iOxFeO1_3 == i) dmOxds[j] += 2.0*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
  dmOxTotds[j] += 2.0*y*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
            }
            if      (iOxFe2O3  == i) dmOxds[iCmpFe2SiO4_6] += -2.0*y*nSpecies;
            else if (iOxFeO     == i) dmOxds[iCmpFe2SiO4_6] += -2.0*(1.0-2.0*y)*nSpecies;
            else if (iOxFeO1_3 == i) dmOxds[iCmpFe2SiO4_6] += 2.0*nSpecies;
            dmOxTotds[iCmpFe2SiO4_6] += 2.0*y*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            for (j=0; j<NS; j++) {
                if      (iOxFe2O3  == i) dmOxds[j] += -2.0*y*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
                else if (iOxFeO       == i) dmOxds[j] += -2.0*(1.0-2.0*y)*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
                else if (iOxFeO1_3 == i) dmOxds[j] += 2.0*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
  dmOxTotds[j] += 2.0*y*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
            }
            if      (iOxFe2O3  == i) dmOxds[iCmpFe2AlO4_1] += -2.0*y*nSpecies;
            else if (iOxFeO     == i) dmOxds[iCmpFe2AlO4_1] += -2.0*(1.0-2.0*y)*nSpecies;
            else if (iOxFeO1_3 == i) dmOxds[iCmpFe2AlO4_1] += 2.0*nSpecies;
            dmOxTotds[iCmpFe2AlO4_1] += 2.0*y*nSpecies;
        }
        for (j=0; j<NS; j++) {
            dvds[j]       += dmOxds[j]*bulkSystem[i].gk_v;
            d2vdsdt[j]   += dmOxds[j]*bulkSystem[i].gk_dvdt;
            dmwds[j]     += dmOxds[j]*bulkSystem[i].mw;
            dcpds[j]     += dmOxds[j]*bulkSystem[i].gk_cp;
            dcRefds[j]   += dmOxds[j]*bulkSystem[i].gk_c/mOxTot - mOx[i]*bulkSystem[i].gk_c*dmOxTotds[j]/(mOxTot*mOxTot);
            dcRefds[j]   += (iOxAl2O3 != -1) ? dmOxds[j]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
              - 2.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotds[j]/(mOxTot*mOxTot*mOxTot) : 0.0;
            d2cdsdt[j]   += dmOxds[j]*bulkSystem[i].gk_dcdt/mOxTot - mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotds[j]/(mOxTot*mOxTot);
            d3vdsdp2[j]  += dmOxds[j]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d4vdsdp3[j]  += dmOxds[j]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d5vdsdp4[j]  += dmOxds[j]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        }
    }
    if (v == 0.0) return;

    alpha   = dvdt/v;
    v0      = v*exp(alpha*(t-tr));
    v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
    c      = cRef + (t-tr)*dcdt;
    v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
    v2      = d2vdp2;
    denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
    a      = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
    b      = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;
    sum      = a*a - 4.0*b;

    for (i=0; i<NS; i++) {
        double dalphads = d2vdsdt[i]/v - dvdt*dvds[i]/(v*v);
        double dv0ds    = (dvds[i] + v*dalphads*(t-tr))*exp(alpha*(t-tr));
        double dv1Refds = -2.0*v*dvds[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                    -v*v*(-1000.0*(cRef*dmwds[i]+2.0*mw*dcRefds[i])/(mw*mw*cRef*cRef*cRef)
                      + 2.0*tr*alpha*dalphads/(cp) - tr*alpha*alpha*dcpds[i]/(cp*cp));
        double dcds     = dcRefds[i] + (t-tr)*d2cdsdt[i];
        double dv1ds    = -2.0*v0*dv0ds*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                    -v0*v0*(-1000.0*(c*dmwds[i]+2.0*mw*dcds)/(mw*mw*c*c*c)
                        + 2.0*t*alpha*dalphads/(cp) - t*alpha*alpha*dcpds[i]/(cp*cp));
        double dv2ds    = d3vdsdp2[i];
        double ddenomds = 2.0*dv1Refds*d3vdp3+2.0*v1Ref*d4vdsdp3[i]-6.0*d2vdp2*d3vdsdp2[i];
        double dads     = (d3vdsdp2[i]*d3vdp3+d2vdp2*d4vdsdp3[i]-dv1Refds*d4vdp4/2.0-v1Ref*d5vdsdp4[i]/2.0)/denom
                - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomds;
        double dbds     = (d3vdsdp2[i]*d4vdp4/4.0+d2vdp2*d5vdsdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdsdp3[i])/denom
                - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomds;
        double dv0dt    = v*alpha*exp(alpha*(t-tr));
        double dv1dt    = - 2.0*v0*v0*alpha*(1000.0/(mw*c*c) + t*alpha*alpha/(cp)) - v0*v0*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp));

        double d2v0dsdt = v*dalphads*exp(alpha*(t-tr)) + (dvds[i] + v*dalphads*(t-tr))*alpha*exp(alpha*(t-tr));
        double d2v1dsdt = -2.0*(dv0dt*dv0ds + v0*d2v0dsdt)*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                - 2.0*v0*dv0ds*(- 2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp))
                - 2.0*v0*dv0dt*(- 1000.0*(c*dmwds[i]+2.0*mw*dcds)/(mw*mw*c*c*c)
                  + 2.0*t*alpha*dalphads/(cp) - t*alpha*alpha*dcpds[i]/(cp*cp))
              - v0*v0*(- 1000.0*(dcdt*dmwds[i]+2.0*mw*d2cdsdt[i])/(mw*mw*c*c*c) + 3000.0*(c*dmwds[i]+2.0*mw*dcds)*dcdt/(mw*mw*c*c*c*c)
                 + 2.0*alpha*dalphads/(cp) - alpha*alpha*dcpds[i]/(cp*cp));
        double d2gIntdsdt = 0.0;

        if ((a == 0.0) && (b == 0.0)) {
            d2gIntdsdt       = d2v0dsdt*(p-pr) + d2v1dsdt*(p-pr)*(p-pr)/2.0;

        } else if ((a != 0.0) && (b == 0.0)) {
            printf("*-->Exception in d2IntegralV_GKdsdt (liquid.c). a is greater than zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
            d2gIntdsdt       = 0.0;

        } else if ((a == 0.0) && (b != 0.0)) {
            printf("*-->Exception in d2IntegralV_GKdsdt (liquid.c). a is zero, b is greater than zero.\n"); liqERRstate = ERR_A_ZERO;
            d2gIntdsdt       = 0.0;

        } else if (sum > 0.0) {
            d2gIntdsdt       = d2gdrdtGMAP(p/10000.0, pr/10000.0,
                                     v0,             v1*10000.0,    v2*10000.0*10000.0,    a*10000.0,    b*10000.0*10000.0,
             dv0ds,       dv1ds*10000.0, dv2ds*10000.0*10000.0, dads*10000.0, dbds*10000.0*10000.0,
             dv0dt,       dv1dt*10000.0,
             d2v0dsdt, d2v1dsdt*10000.0)*10000.0;

        } else if (sum == 0.0) {
            printf("*-->Exception in d2IntegralV_GKdsdt (liquid.c). a*a-4*b is equal to zero.\n"); liqERRstate = ERR_SUM_ZERO;
            d2gIntdsdt       = 0.0;

        } else if(sum < 0.0) {
            d2gIntdsdt       = d2gdrdtLMAP(p/10000.0, pr/10000.0,
                                     v0,             v1*10000.0,    v2*10000.0*10000.0,    a*10000.0,    b*10000.0*10000.0,
             dv0ds,       dv1ds*10000.0, dv2ds*10000.0*10000.0, dads*10000.0, dbds*10000.0*10000.0,
             dv0dt,       dv1dt*10000.0,
             d2v0dsdt, d2v1dsdt*10000.0)*10000.0;

        }

        dsdt[i] += d2gIntdsdt;

    }
}

static double d2IntegralV_GKdT2(double r[NR], double s[NT], double t, double p)  {
    const double pr      = 1.0;
    const double tr      = 1673.15;
    double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum, denom;
    double dv0dt, dv1dt, d2v0dt2, d2v1dt2, d2gIntdt2;
    int i, j;

    if (fabs(p-pr) < 10.0*DBL_EPSILON) return (double) 0.0;

    /* Convert input composition (r) to liquid moles (m)  */
    for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* Compute moles and total moles of oxides */
    for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
    if (mOxTot == 0.0) return 0.0;

    /* Deal with the special case of FeO1.3 */
    mOx[NA] = 0.0;
    if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
        const double y = 0.3;
        mOx[iOxFeO1_3] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
        }
    }

    for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
        v        += mOx[i]*bulkSystem[i].gk_v;
        dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
        cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
        dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
        cp      += mOx[i]*bulkSystem[i].gk_cp;
        d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
        d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
        d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        mw      += mOx[i]*bulkSystem[i].mw;
    }
    if (v == 0.0) return 0.0;

    alpha   = dvdt/v;
    v0      = v*exp(alpha*(t-tr));
    v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
    c      = cRef + (t-tr)*dcdt;
    v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
    v2      = d2vdp2;
    denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
    a      = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
    b      = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;
    sum      = a*a - 4.0*b;

    dv0dt   = alpha*v*exp(alpha*(t-tr));
    dv1dt   = - 2.0*v0*dv0dt*(1000.0/(mw*c*c) + t*alpha*alpha/(cp)) - v0*v0*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp));
    d2v0dt2 = alpha*alpha*v*exp(alpha*(t-tr));
    d2v1dt2 = - 2.0*(dv0dt*dv0dt+v0*d2v0dt2)*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
        - 4.0*v0*dv0dt*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp))
                        - v0*v0*(6000.0*dcdt*dcdt/(mw*c*c*c*c));
    d2gIntdt2 = 0.0;

    if ((a == 0.0) && (b == 0.0)) {
        d2gIntdt2 += alpha*alpha*v0*(p-pr) + d2v1dt2*(p-pr)*(p-pr)/2.0;

    } else if ((a != 0.0) && (b == 0.0)) {
        printf("*-->Exception in d2IntegralV_GKdT2 (liquid.c). a is greater than zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
        d2gIntdt2        = 0.0;

    } else if ((a == 0.0) && (b != 0.0)) {
        printf("*-->Exception in d2IntegralV_GKdT2 (liquid.c). a is zero, b is greater than zero.\n"); liqERRstate = ERR_A_ZERO;
        d2gIntdt2        = 0.0;

    } else if (sum > 0.0) {
        d2gIntdt2 = d2gdt2GMAP(p/10000.0, pr/10000.0,
                           v0,           v1*10000.0, v2*10000.0*10000.0, a*10000.0, b*10000.0*10000.0,
            dv0dt,     dv1dt*10000.0,
            d2v0dt2, d2v1dt2*10000.0)*10000.0;

    } else if (sum == 0.0) {
        printf("*-->Exception in d2IntegralV_GKdT2 (liquid.c). a*a-4*b is equal to zero.\n"); liqERRstate = ERR_SUM_ZERO;
        d2gIntdt2        = 0.0;

    } else if(sum < 0.0) {
        d2gIntdt2 = d2gdt2LMAP(p/10000.0, pr/10000.0,
                           v0,           v1*10000.0, v2*10000.0*10000.0, a*10000.0, b*10000.0*10000.0,
            dv0dt,     dv1dt*10000.0,
            d2v0dt2, d2v1dt2*10000.0)*10000.0;

    }

    return d2gIntdt2;
}

static double d2IntegralV_GKdTdP(double r[NR], double s[NT], double t, double p) {
    const double pr      = 1.0;
    const double tr      = 1673.15;
    double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum;
    double dv0dt, dv1dt, d2gIntdtdp;
    int i, j;

    /* Convert input composition (r) to liquid moles (m)  */
    for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* Compute moles and total moles of oxides */
    for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
    if (mOxTot == 0.0) return 0.0;

    /* Deal with the special case of FeO1.3 */
    mOx[NA] = 0.0;
    if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
        const double y = 0.3;
        mOx[iOxFeO1_3] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
        }
    }

    for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
        v        += mOx[i]*bulkSystem[i].gk_v;
        dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
        cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
        dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
        cp      += mOx[i]*bulkSystem[i].gk_cp;
        d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
        d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
        d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        mw      += mOx[i]*bulkSystem[i].mw;
    }
    if (v == 0.0) return 0.0;

    alpha   = dvdt/v;
    v0      = v*exp(alpha*(t-tr));
    v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
    c      = cRef + (t-tr)*dcdt;
    v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
    v2      = d2vdp2;
    a      = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;
    b      = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;
    sum      = a*a - 4.0*b;
    dv0dt   = alpha*v0;
    dv1dt   = - 2.0*v0*v0*alpha*(1000.0/(mw*c*c) + t*alpha*alpha/(cp)) - v0*v0*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp));

    if ((a == 0.0) && (b == 0.0)) {
        d2gIntdtdp = alpha*v0 + dv1dt*(p-pr);

    } else if ((a != 0.0) && (b == 0.0)) {
        printf("*-->Exception in d2IntegralV_GKdTdP (liquid.c). a is greater than zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
        d2gIntdtdp = 0.0;

    } else if ((a == 0.0) && (b != 0.0)) {
        printf("*-->Exception in d2IntegralV_GKdTdP (liquid.c). a is zero, b is greater than zero.\n"); liqERRstate = ERR_A_ZERO;
        d2gIntdtdp = 0.0;

    } else {
        d2gIntdtdp = d2gdtdpMAP(p, pr, v0, v1, v2, a, b, dv0dt, dv1dt);

    }
    return d2gIntdtdp;
}

static double d2IntegralV_GKdP2(double r[NR], double s[NT], double t, double p) {
    const double pr      = 1.0;
    const double tr      = 1673.15;
    double m[NA], mOx[NA+1], mOxTot, v, dvdt, cRef, c, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum;
    double d2gdp2;
    int i, j;

    /* Convert input composition (r) to liquid moles (m)  */
    for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* Compute moles and total moles of oxides */
    for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
    if (mOxTot == 0.0) return 0.0;

    /* Deal with the special case of FeO1.3 */
    mOx[NA] = 0.0;
    if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
        const double y = 0.3;
        mOx[iOxFeO1_3] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
        }
    }

    for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
        v        += mOx[i]*bulkSystem[i].gk_v;
        dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
        cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
        dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
        cp      += mOx[i]*bulkSystem[i].gk_cp;
        d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
        d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
        d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        mw      += mOx[i]*bulkSystem[i].mw;
    }
    if (v == 0.0) return 0.0;

    alpha   = dvdt/v;
    v0      = v*exp(alpha*(t-tr));
    v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
    c      = cRef + (t-tr)*dcdt;
    v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
    v2      = d2vdp2;
    a      = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;
    b      = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;
    sum      = a*a - 4.0*b;

    if ((a == 0.0) && (b == 0.0)) {
        d2gdp2 = v1 + v2*(p-pr);

    } else if ((a != 0.0) && (b == 0.0)) {
        printf("*-->Exception in d2IntegralV_GKdP2 (liquid.c). a is greater than zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
        d2gdp2 = 0.0;

    } else if ((a == 0.0) && (b != 0.0)) {
        printf("*-->Exception in d2IntegralV_GKdP2 (liquid.c). a is zero, b is greater than zero.\n"); liqERRstate = ERR_A_ZERO;
        d2gdp2 = 0.0;

    } else {
        d2gdp2 = d2gdp2MAP(p, pr, v0, v1, v2, a, b);

    }

    return d2gdp2;
}

#endif /* USE_GHIORSO_KRESS_MODEL */

static double fillG (double r[NR], double s[NT], double t, double p) {
    double result, config;
    int i, j;

    /* Taylor expansion and standard state terms */
    result = Gconst;
    for (i=0; i<NR; i++) {
        result += gr[i]*r[i];
        for (j=i; j<NR; j++) result += grr[i][j]*r[i]*r[j];
        for (j=0; j<NS; j++) result += grs[i][j]*r[i]*s[j];
    }
    for (i=0; i<NS; i++) {
        result += gs[i]*s[i];
        for (j=i; j<NS; j++) result += gss[i][j]*s[i]*s[j];
    }

    /* Configurational entropy terms */
    for (i=0, config=0.0; i<NE; i++) if (xSpecies[i] > 0.0) config += xSpecies[i]*log(xSpecies[i]);
    if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0) config += xSpecies[nH2O]*log(xSpecies[nH2O]) + (1.0-xSpecies[nH2O])*log(1.0-xSpecies[nH2O]);
    result += R*t*nSpecies*config;

    ternaryH2OCO2terms(FIRST, r);
    result += gT;

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double yIV    = 1.0;
        double gIntIV = integralV_GK(r, s, t, p);
        double gInt;

        for (i=0; i<NY; i++) yIV -= (s[NS+i] > 0.0) ? s[NS+i] : DBL_EPSILON;
        gInt = gIntIV + R*t*yIV*log(yIV);

        for (i=0; i<NY; i++) {
            double f = fCN[i];
            double deltaG_GK, y;

            y         = (s[NS+i] > 0.0) ? s[NS+i] : DBL_EPSILON;
            deltaG_GK = ssCnModelParameters[i*nc+0].enthalpy - t*ssCnModelParameters[i*nc+0].entropy;
            for (j=0; j<NR; j++) deltaG_GK += (  (ssCnModelParameters[i*nc+j+1].enthalpy - t*ssCnModelParameters[i*nc+j+1].entropy)
                                                                    - (ssCnModelParameters[i*nc  +0].enthalpy - t*ssCnModelParameters[i*nc  +0].entropy))*r[j];

            gInt += R*t*y*log(y) + y*deltaG_GK + (f-1.0)*y*gIntIV;
        }

        result += gInt;
    }
#endif /* USE_GHIORSO_KRESS_MODEL */

    return result;
}

static double fillS (double r[NR], double s[NT], double t, double p) {
    double result, config;
    int i, j;

    /* Taylor expansion and standard state terms */
    result = Sconst;
    for (i=0; i<NR; i++) {
        result += sr[i]*r[i];
        for (j=i; j<NR; j++) result += srr[i][j]*r[i]*r[j];
        for (j=0; j<NS; j++) result += srs[i][j]*r[i]*s[j];
    }
    for (i=0; i<NS; i++) {
        result += ss[i]*s[i];
        for (j=i; j<NS; j++) result += sss[i][j]*s[i]*s[j];
    }

    /* Configurational entropy terms */
    for (i=0, config=0.0; i<NE; i++) if (xSpecies[i] > 0.0) config += xSpecies[i]*log(xSpecies[i]);
    if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0) config += xSpecies[nH2O]*log(xSpecies[nH2O]) + (1.0-xSpecies[nH2O])*log(1.0-xSpecies[nH2O]);
    result -= R*nSpecies*config;

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double yIV       = 1.0;
        double dgIntIVdt = dIntegralV_GKdT(r, s, t, p);
        double dgIntdt;

        for (i=0; i<NY; i++) yIV -= (s[NS+i] > 0.0) ? s[NS+i] : DBL_EPSILON;
        dgIntdt = dgIntIVdt + R*yIV*log(yIV);

        for (i=0; i<NY; i++) {
            double f = fCN[i];
            double dDeltaG_GKdT, y;

            y            = (s[NS+i] > 0.0) ? s[NS+i] : DBL_EPSILON;
            dDeltaG_GKdT = - ssCnModelParameters[i*nc+0].entropy;
            for (j=0; j<NR; j++) dDeltaG_GKdT += (-ssCnModelParameters[i*nc+j+1].entropy + ssCnModelParameters[i*nc+0].entropy)*r[j];

            dgIntdt += R*y*log(y) + y*dDeltaG_GKdT + (f-1.0)*y*dgIntIVdt;
        }

        result += -dgIntdt;
    }

/*  Maple-derived equivalent functions:
dgdtLMAP(p, pr, v0, v1, v2, a, b, dv0dt, dv1dt)
dgdtGMAP(p, pr, v0, v1, v2, a, b, dv0dt, dv1dt)
*/

#endif /* USE_GHIORSO_KRESS_MODEL */

    return result;
}

#ifdef USE_GHIORSO_KRESS_MODEL

static double fillV (double r[NR], double s[NT], double t, double p) {
    const double pr      = 1.0;
    const double tr      = 1673.15;
    double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, coeffCN;
    int i, j;
    double result = 0.0;

    /* Convert input composition (r) to liquid moles (m)  */
    for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* Compute moles and total moles of oxides */
    for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
    if (mOxTot == 0.0) return result;

    /* Deal with the special case of FeO1.3 */
    mOx[NA] = 0.0;
    if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
        const double y = 0.3;
        mOx[iOxFeO1_3] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
        }
    }

    for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
        v        += mOx[i]*bulkSystem[i].gk_v;
        dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
        cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
        dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
        cp      += mOx[i]*bulkSystem[i].gk_cp;
        d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
        d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
        d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        mw      += mOx[i]*bulkSystem[i].mw;
    }
    if (v == 0.0) return result;

    alpha   = dvdt/v;
    v0      = v*exp(alpha*(t-tr));
    v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
    c      = cRef + (t-tr)*dcdt;
    v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
    v2      = d2vdp2;
    a      = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;
    b      = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;

    if (1.0 + a*(p-pr) + b*(p-pr)*(p-pr) == 0.0) return result;

    for (i=0, coeffCN=1.0; i<NY; i++) coeffCN += (fCN[i]-1.0)*s[NS+i];

    result = coeffCN*((v0 + (v1+v0*a)*(p-pr) + (v2/2.0+v1*a+v0*b)*(p-pr)*(p-pr))/(1.0 + a*(p-pr) + b*(p-pr)*(p-pr)));
    return result;
}

#else

static double fillV (double r[NR], double s[NT], double t, double p) {
    double result;
    int i, j;

    /* Taylor expansion and standard state terms */
    result = Vconst;
    for (i=0; i<NR; i++) {
        result += vr[i]*r[i];
        for (j=i; j<NR; j++) result += vrr[i][j]*r[i]*r[j];
        for (j=0; j<NS; j++) result += vrs[i][j]*r[i]*s[j];
    }
    for (i=0; i<NS; i++) {
        result += vs[i]*s[i];
        for (j=i; j<NS; j++) result += vss[i][j]*s[i]*s[j];
    }

    return result;
}

#endif /* USE_GHIORSO_KRESS_MODEL */


static void fillDGDR (double r[NR], double s[NT], double t, double p, double *result) {
    int i, j;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) {
        result[i] = gr[i] + grr[i][i]*r[i];
        for (j=0; j<NR; j++) result[i] += grr[i][j]*r[j];
        for (j=0; j<NS; j++) result[i] += grs[i][j]*s[j];
    }

    /* Configurational entropy terms */
    for (j=0; j<NR; j++) {
        double config = 0.0;
        for (i=0; i<NE; i++) if (xSpecies[i] > 0.0 && dxSpeciesdr[i][j] != 0.0) config += dxSpeciesdr[i][j]*(1.0 + log(xSpecies[i]));
        if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0 && dxSpeciesdr[nH2O][j] != 0.0)
            config += dxSpeciesdr[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]));
        result[j] += R*t*nSpecies*config;
    }

    ternaryH2OCO2terms(SECOND, r);
    for (j=0; j<NR; j++) result[j] += dgdrT[j];

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double dgIntIVdr[NR];

        dIntegralV_GKdr(r, s, t, p, dgIntIVdr);

        for (i=0; i<NR; i++) {
            result[i] += dgIntIVdr[i];
            for (j=0; j<NY; j++) {
                double y = (s[NS+j] > 0.0) ? s[NS+j] : DBL_EPSILON;
                double f = fCN[j];

                result[i] += y*( (ssCnModelParameters[j*nc+i+1].enthalpy - t*ssCnModelParameters[j*nc+i+1].entropy)
                  -(ssCnModelParameters[j*nc  +0].enthalpy - t*ssCnModelParameters[j*nc  +0].entropy))
                        + (f-1.0)*y*dgIntIVdr[i];
            }
        }

    }
#endif /* USE_GHIORSO_KRESS_MODEL */

}

static void fillDGDS (double r[NR], double s[NT], double t, double p, double *result) {
    int i, j;

    memset(result, '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) {
        result[i] = gs[i] + gss[i][i]*s[i];
        for (j=0; j<NS; j++) result[i] += gss[i][j]*s[j];
        for (j=0; j<NR; j++) result[i] += grs[j][i]*r[j];
    }

    /* Configurational entropy terms */
    for (j=0; j<NS; j++) {
        double config = 0.0;
        for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
            config += dnSpeciesds[j]*xSpecies[i]*log(xSpecies[i]) + nSpecies*dxSpeciesds[i][j]*(1.0 + log(xSpecies[i]));
        if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
            config += dnSpeciesds[j]*(xSpecies[nH2O]*log(xSpecies[nH2O]) + (1.0-xSpecies[nH2O])*log(1.0-xSpecies[nH2O]))
                            + nSpecies*dxSpeciesds[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]));
        result[j] += R*t*config;
    }

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double yIV    = 1.0;
        double gIntIV = integralV_GK(r, s, t, p);
        double dgIntIVds[NS];

        dIntegralV_GKds(r, s, t, p, dgIntIVds);

        for (i=0; i<NS; i++) {
            result[i] += dgIntIVds[i];
            for (j=0; j<NY; j++) {
                double y = (s[NS+j] > 0.0) ? s[NS+j] : DBL_EPSILON;
                double f = fCN[j];
                result[i] += y*( (ssCnModelParameters[j*nc+i+1].enthalpy - t*ssCnModelParameters[j*nc+i+1].entropy)
                  -(ssCnModelParameters[j*nc  +0].enthalpy - t*ssCnModelParameters[j*nc  +0].entropy))
                        + (f-1.0)*y*dgIntIVds[i];
            }
        }

        for (i=0; i<NY; i++) {
            double f =fCN[i];
            double deltaG_GK, y;

            y         = (s[NS+i] > 0.0) ? s[NS+i] : DBL_EPSILON;
            deltaG_GK = ssCnModelParameters[i*nc+0].enthalpy - t*ssCnModelParameters[i*nc+0].entropy;

            for (j=0; j<NR; j++) deltaG_GK += ( (ssCnModelParameters[i*nc+j+1].enthalpy - t*ssCnModelParameters[i*nc+j+1].entropy)
                                                                    -(ssCnModelParameters[i*nc  +0].enthalpy - t*ssCnModelParameters[i*nc  +0].entropy))*r[j];

            yIV  -= y;
            result[NS+i] += R*t*log(y) + deltaG_GK + (f-1.0)*gIntIV;
        }

        for (i=0; i<NY; i++) result[NS+i] -= R*t*log(yIV);

    }
#endif /* USE_GHIORSO_KRESS_MODEL */

}

static void fillDGDW (double r[NR], double s[NT], double t, double p, double *result) {
    int i, j, k, l, m, n;

    /*******************************
    * Parameters: NW WH(), NE H() *
    *             NW WS(), NE S() *
    *             NW WV(), NE V() *
    *******************************/
    memset(result, '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
    * NW W parameters solution are first *
    **************************************/
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        result[     n] +=         taylorCoeff[n+NE][0];
        result[  NP+n] +=      -t*taylorCoeff[n+NE][0];
#ifndef USE_GHIORSO_KRESS_MODEL
        result[2*NP+n] += (p-1.0)*taylorCoeff[n+NE][0];
#endif /* USE_GHIORSO_KRESS_MODEL */
        m = 0;
        for (j=0, m=0; j<NR; j++) {
            result[     n] +=         taylorCoeff[n+NE][1+j]*r[j];
            result[  NP+n] +=      -t*taylorCoeff[n+NE][1+j]*r[j];
#ifndef USE_GHIORSO_KRESS_MODEL
            result[2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+j]*r[j];
#endif /* USE_GHIORSO_KRESS_MODEL */
            for (k=j; k<NR; k++, m++) {
                result[     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*r[j]*r[k];
                result[  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*r[j]*r[k];
#ifndef USE_GHIORSO_KRESS_MODEL
                result[2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*r[j]*r[k];
#endif /* USE_GHIORSO_KRESS_MODEL */
            }
            for (k=0; k<NS; k++, m++) {
                result[     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*r[j]*s[k];
                result[  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*r[j]*s[k];
#ifndef USE_GHIORSO_KRESS_MODEL
                result[2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*r[j]*s[k];
#endif /* USE_GHIORSO_KRESS_MODEL */
            }
        }
        for (j=0; j<NS; j++) {
            result[     n] +=         taylorCoeff[n+NE][1+NR+j]*s[j];
            result[  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+j]*s[j];
#ifndef USE_GHIORSO_KRESS_MODEL
            result[2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+j]*s[j];
#endif /* USE_GHIORSO_KRESS_MODEL */
            for (k=j; k<NS; k++, m++) {
                result[     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*s[j]*s[k];
                result[  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*s[j]*s[k];
#ifndef USE_GHIORSO_KRESS_MODEL
                result[2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*s[j]*s[k];
#endif /* USE_GHIORSO_KRESS_MODEL */
            }
        }
    }

    /**************************************
    * NE standard state terms are second *
    **************************************/
    for (i=0; i<NE; i++, n++) {
        result[        n] +=         taylorCoeff[i][0];
        result[  NP+n] +=      -t*taylorCoeff[i][0];
#ifndef USE_GHIORSO_KRESS_MODEL
        result[2*NP+n] += (p-1.0)*taylorCoeff[i][0];
#endif /* USE_GHIORSO_KRESS_MODEL */
        for (j=0; j<NR; j++) {
            result[        n] +=           taylorCoeff[i][1+j]*r[j];
            result[  NP+n] +=      -t*taylorCoeff[i][1+j]*r[j];
#ifndef USE_GHIORSO_KRESS_MODEL
            result[2*NP+n] += (p-1.0)*taylorCoeff[i][1+j]*r[j];
#endif /* USE_GHIORSO_KRESS_MODEL */
        }
        for (j=0; j<NS; j++) {
            result[        n] +=         taylorCoeff[i][1+NR+j]*s[j];
            result[  NP+n] +=      -t*taylorCoeff[i][1+NR+j]*s[j];
#ifndef USE_GHIORSO_KRESS_MODEL
            result[2*NP+n] += (p-1.0)*taylorCoeff[i][1+NR+j]*s[j];
#endif /* USE_GHIORSO_KRESS_MODEL */
        }
    }
}

static void fillD2GDR2 (double r[NR], double s[NT], double t, double p, double result[NR][NR]) {
    int i, j, k;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) {
        result[i][i] = 2.0*grr[i][i];
        for (j=i+1; j<NR; j++) {
            result[i][j] = grr[i][j];
            result[j][i] = grr[i][j];
        }
    }

    /* Configurational entropy terms */
    for (j=0; j<NR; j++) {
        for (k=j; k<NR; k++) {
            double config = 0.0;
            for (i=0; i<NE; i++) if (xSpecies[i] > 0.0 && dxSpeciesdr[i][j] != 0.0 && dxSpeciesdr[i][k] != 0.0)
                config += dxSpeciesdr[i][j]*dxSpeciesdr[i][k]/xSpecies[i];
            if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0 && dxSpeciesdr[nH2O][j] != 0.0 && dxSpeciesdr[nH2O][k] != 0.0)
                config += dxSpeciesdr[nH2O][j]*dxSpeciesdr[nH2O][k]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]));
            result[j][k] += R*t*nSpecies*config;
            result[k][j]  = result[j][k];
        }
    }

    ternaryH2OCO2terms(THIRD, r);
    for (j=0; j<NR; j++) for (i=j; i<NR; i++) { result[j][i] += d2gdr2T[j][i]; result[i][j] = result[j][i]; }

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double d2gIntIVdr2[NR][NR];

        d2IntegralV_GKdr2(r, s, t, p, d2gIntIVdr2);

        for (i=0; i<NR; i++) {
            for (j=i; j<NR; j++) {

                result[i][j] += d2gIntIVdr2[i][j];
  if (i != j) result[j][i] += d2gIntIVdr2[i][j];

  for (k=0; k<NY; k++) {
                    double y = (s[NS+k] > 0.0) ? s[NS+k] : DBL_EPSILON;
                    double f = fCN[k];
    result[i][j] += (f-1.0)*y*d2gIntIVdr2[i][j];
    if (i != j) result[j][i] += (f-1.0)*y*d2gIntIVdr2[i][j];
  }

            }
        }
    }
#endif /* USE_GHIORSO_KRESS_MODEL */

}

static void fillD2GDRDS (double r[NR], double s[NT], double t, double p, double result[NR][NT]) {
    int i, j, k;

    for (j=0; j<NR; j++) memset(result[j], '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) for (j=0; j<NS; j++) result[i][j] = grs[i][j];

    /* Configurational entopy terms */
    for (j=0; j<NR; j++) {
        for (k=0; k<NS; k++) {
            double config = 0.0;
            for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
                config += nSpecies*(d2xSpeciesdrds[i][j][k]*log(xSpecies[i]) + dxSpeciesdr[i][j]*dxSpeciesds[i][k]/xSpecies[i])
                    + dnSpeciesds[k]*(dxSpeciesdr[i][j]*(1.0 + log(xSpecies[i])));
            if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
                config += dnSpeciesds[k]*dxSpeciesdr[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                    + nSpecies*(d2xSpeciesdrds[nH2O][j][k]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                    + dxSpeciesds[nH2O][k]*dxSpeciesdr[nH2O][j]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O])));
            result[j][k] += R*t*config;
        }
    }

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double dgIntIVdr[NR];
        double d2gIntIVdrds[NR][NS];

        d2IntegralV_GKdrds(r, s, t, p, d2gIntIVdrds);

        for (i=0; i<NR; i++) {
            for (j=0; j<NS; j++) {

                result[i][j] += d2gIntIVdrds[i][j];

  for (k=0; k<NY; k++) {
                    double y = (s[NS+k] > 0.0) ? s[NS+k] : DBL_EPSILON;
                    double f = fCN[k];
    result[i][j] += (f-1.0)*y*d2gIntIVdrds[i][j];
  }

            }
        }

        dIntegralV_GKdr(r, s, t, p, dgIntIVdr);

        for (i=0; i<NR; i++) {
            for (j=0; j<NY; j++) {
                double f = fCN[j];

                result[i][NS+j] += ( (ssCnModelParameters[j*nc+i+1].enthalpy - t*ssCnModelParameters[j*nc+i+1].entropy)
                      -(ssCnModelParameters[j*nc  +0].enthalpy - t*ssCnModelParameters[j*nc  +0].entropy))
                                    + (f-1.0)*dgIntIVdr[i];
            }
        }

    }
#endif /* USE_GHIORSO_KRESS_MODEL */
}

static void fillD2GDRDT (double r[NR], double s[NT], double t, double p, double *result) {
    int i, j;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) {
        result[i] = -sr[i] + -srr[i][i]*r[i];
        for (j=0; j<NR; j++) result[i] += -srr[i][j]*r[j];
        for (j=0; j<NS; j++) result[i] += -srs[i][j]*s[j];
    }

    /* Configurational entropy terms */
    for (j=0; j<NR; j++) {
        double config = 0.0;
        for (i=0; i<NE; i++) if (xSpecies[i] > 0.0 && dxSpeciesdr[i][j] != 0.0) config += dxSpeciesdr[i][j]*(1.0 + log(xSpecies[i]));
        if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0 && dxSpeciesdr[nH2O][j] != 0.0)
            config += dxSpeciesdr[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]));
        result[j] += R*nSpecies*config;
    }

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double d2gIntIVdrdt[NR];

        d2IntegralV_GKdrdt(r, s, t, p, d2gIntIVdrdt);

        for (i=0; i<NR; i++) {
            result[i] += d2gIntIVdrdt[i];

            for (j=0; j<NY; j++) {
                double y = (s[NS+j] > 0.0) ? s[NS+j] : DBL_EPSILON;
                double f = fCN[j];

                result[i] += y*(-ssCnModelParameters[j*nc+i+1].entropy + ssCnModelParameters[j*nc  +0].entropy) + (f-1.0)*y*d2gIntIVdrdt[i];
            }

        }
    }
#endif /* USE_GHIORSO_KRESS_MODEL */

}

static void fillD2GDRDP (double r[NR], double s[NT], double t, double p, double *result) {
    int i, j;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) {
#ifdef USE_GHIORSO_KRESS_MODEL
        result[i] = 0.0;
#else
        result[i] = vr[i] + vrr[i][i]*r[i];
        for (j=0; j<NR; j++) result[i] += vrr[i][j]*r[j];
        for (j=0; j<NS; j++) result[i] += vrs[i][j]*s[j];
#endif /* USE_GHIORSO_KRESS_MODEL */
    }

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double coeffCN, d2gIntIVdrdp[NR];

        d2IntegralV_GKdrdp(r, s, t, p, d2gIntIVdrdp);
        for (j=0, coeffCN=1.0; j<NY; j++) coeffCN += (fCN[j]-1.0)*s[NS+j];

        for (i=0; i<NR; i++) result[i] += coeffCN*d2gIntIVdrdp[i];
    }
#endif /* USE_GHIORSO_KRESS_MODEL */
}

static void fillD2GDRDW (double r[NR], double s[NT], double t, double p, double result[NR][3*NP]) {
    int i, j, k, l, m, n, ii;

    /*******************************
    * Parameters: NW WH(), NE H() *
    *             NW WS(), NE S() *
    *             NW WV(), NE V() *
    *******************************/
    for (ii=0; ii<NR; ii++) memset(result[ii], '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
    * NW W parameters solution are first *
    **************************************/
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        for (ii=0; ii<NR; ii++) {
            result[ii][     n] +=        taylorCoeff[n+NE][1+ii];
            result[ii][  NP+n] +=     -t*taylorCoeff[n+NE][1+ii];
#ifndef USE_GHIORSO_KRESS_MODEL
            result[ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+ii];
#endif /* USE_GHIORSO_KRESS_MODEL */
            m = 0;
            for (j=0, m=0; j<NR; j++) {
                for (k=j; k<NR; k++, m++) {
                    if (j == ii) {
                        result[ii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*s[k];
                        result[ii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*s[k];
#ifndef USE_GHIORSO_KRESS_MODEL
                        result[ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*s[k];
#endif /* USE_GHIORSO_KRESS_MODEL */
                    }
                    if (k == ii) {
                        result[ii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*r[j];
                        result[ii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*r[j];
#ifndef USE_GHIORSO_KRESS_MODEL
                        result[ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*r[j];
#endif /* USE_GHIORSO_KRESS_MODEL */
                    }
                }
                for (k=0; k<NS; k++, m++) {
                    if (j == ii) {
                        result[ii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*s[k];
                        result[ii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*s[k];
#ifndef USE_GHIORSO_KRESS_MODEL
                        result[ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*s[k];
#endif /* USE_GHIORSO_KRESS_MODEL */
                    }
                }
            }
        }
    }

    /**************************************
   * NE standard state terms are second *
   **************************************/
    for (i=0; i<NE; i++, n++) {
        for (ii=0; ii<NR; ii++) {
            result[ii][     n] +=         taylorCoeff[i][1+ii];
            result[ii][  NP+n] +=      -t*taylorCoeff[i][1+ii];
#ifndef USE_GHIORSO_KRESS_MODEL
            result[ii][2*NP+n] += (p-1.0)*taylorCoeff[i][1+ii];
#endif /* USE_GHIORSO_KRESS_MODEL */
        }
    }
}

static void fillD2GDS2New (double r[NR], double s[NT], double t, double p, double **result) {
    int i, j, k;

    for (j=0; j<NT; j++) memset(result[j], '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) {
        result[i][i] = 2.0*gss[i][i];
        for (j=i+1; j<NS; j++) {
            result[i][j] = gss[i][j];
            result[j][i] = gss[i][j];
        }
    }

    /* Configurational entropy terms */
    for (j=0; j<NS; j++) {
        for (k=j; k<NS; k++) {
            double config = 0.0;
            for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
                config += nSpecies*dxSpeciesds[i][j]*dxSpeciesds[i][k]/xSpecies[i]
                + dnSpeciesds[k]*dxSpeciesds[i][j]*(1.0 + log(xSpecies[i]))
                + dnSpeciesds[j]*dxSpeciesds[i][k]*(1.0 + log(xSpecies[i]))
                + d2nSpeciesds2[j][k]*xSpecies[i]*log(xSpecies[i]);
            if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
                config += d2nSpeciesds2[j][k]*(xSpecies[nH2O]*log(xSpecies[nH2O]) + (1.0-xSpecies[nH2O])*log(1.0-xSpecies[nH2O]))
                + dnSpeciesds[j]*dxSpeciesds[nH2O][k]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                + dnSpeciesds[k]*dxSpeciesds[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                + nSpecies*dxSpeciesds[nH2O][j]*dxSpeciesds[nH2O][k]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]));
            result[j][k] += R*t*config;
            result[k][j]  = result[j][k];
        }
    }
}

static void fillD2GDS2 (double r[NR], double s[NT], double t, double p, double result[NT][NT]) {
    int i, j, k;

    for (j=0; j<NT; j++) memset(result[j], '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) {
        result[i][i] = 2.0*gss[i][i];
        for (j=i+1; j<NS; j++) {
            result[i][j] = gss[i][j];
            result[j][i] = gss[i][j];
        }
    }

    /* Configurational entropy terms */
    for (j=0; j<NS; j++) {
        for (k=j; k<NS; k++) {
            double config = 0.0;
            for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
                config += nSpecies*dxSpeciesds[i][j]*dxSpeciesds[i][k]/xSpecies[i]
                + dnSpeciesds[k]*dxSpeciesds[i][j]*(1.0 + log(xSpecies[i]))
                + dnSpeciesds[j]*dxSpeciesds[i][k]*(1.0 + log(xSpecies[i]))
                + d2nSpeciesds2[j][k]*xSpecies[i]*log(xSpecies[i]);
            if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
                config += d2nSpeciesds2[j][k]*(xSpecies[nH2O]*log(xSpecies[nH2O]) + (1.0-xSpecies[nH2O])*log(1.0-xSpecies[nH2O]))
                + dnSpeciesds[j]*dxSpeciesds[nH2O][k]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                + dnSpeciesds[k]*dxSpeciesds[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                + nSpecies*dxSpeciesds[nH2O][j]*dxSpeciesds[nH2O][k]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]));
            result[j][k] += R*t*config;
            result[k][j]  = result[j][k];
        }
    }

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double yIV = 1.0;
        double d2gIntIVds2[NS][NS];

        d2IntegralV_GKds2(r, s, t, p, d2gIntIVds2);

        for (i=0; i<NS; i++) {
            for (j=i; j<NS; j++) {

                result[i][j] += d2gIntIVds2[i][j];
                if (i != j) result[j][i] += d2gIntIVds2[i][j];

                for (k=0; k<NY; k++) {
                    double y = (s[NS+k] > 0.0) ? s[NS+k] : DBL_EPSILON;
                    double f = fCN[k];
                    result[i][j] += (f-1.0)*y*d2gIntIVds2[i][j];
                    if (i != j) result[j][i] += (f-1.0)*y*d2gIntIVds2[i][j];
                }
            }
        }

        for (i=0; i<NY; i++) yIV -= (s[NS+i] > 0.0) ? s[NS+i] : DBL_EPSILON;

        for (i=0; i<NY; i++) {
            double y = (s[NS+i] > 0.0) ? s[NS+i] : DBL_EPSILON;
            result[NS+i][NS+i] = R*t*(1.0/y + 1.0/yIV);
            for (j=i+1; j<NY; j++) { result[NS+i][NS+j] = R*t/yIV; result[NS+j][NS+i] = result[NS+i][NS+j]; }
        }

    }
#endif /* USE_GHIORSO_KRESS_MODEL */

}

static void fillD2GDSDT (double r[NR], double s[NT], double t, double p, double *result) {
    int i, j;

    memset(result, '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) {
        result[i] = -ss[i] - sss[i][i]*s[i];
        for (j=0; j<NS; j++) result[i] += -sss[i][j]*s[j];
        for (j=0; j<NR; j++) result[i] += -srs[j][i]*r[j];
    }

    /* Configurational entropy terms */
    for (j=0; j<NS; j++) {
        double config = 0.0;
        for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
            config += dnSpeciesds[j]*xSpecies[i]*log(xSpecies[i]) + nSpecies*dxSpeciesds[i][j]*(1.0 + log(xSpecies[i]));
        if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
            config += dnSpeciesds[j]*(xSpecies[nH2O]*log(xSpecies[nH2O]) + (1.0-xSpecies[nH2O])*log(1.0-xSpecies[nH2O]))
                            + nSpecies*dxSpeciesds[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]));
        result[j] += R*config;
    }

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double yIV       = 1.0;
        double dgIntIVdt = dIntegralV_GKdT(r, s, t, p);
        double d2gIntIVdsdt[NS];

        d2IntegralV_GKdsdt(r, s, t, p, d2gIntIVdsdt);

        for (i=0; i<NS; i++) {
            result[i] += d2gIntIVdsdt[i];

            for (j=0; j<NY; j++) {
                double y = (s[NS+j] > 0.0) ? s[NS+j] : DBL_EPSILON;
                double f = fCN[j];

                result[i] += y*(-ssCnModelParameters[j*nc+i+1].entropy + ssCnModelParameters[j*nc  +0].entropy) + (f-1.0)*y*d2gIntIVdsdt[i];
            }

        }

        for (i=0; i<NY; i++) yIV -= (s[NS+i] > 0.0) ? s[NS+i] : DBL_EPSILON;

        for (i=0; i<NY; i++) {
            double f = fCN[i];
            double dDeltaG_GKdT, y;

            y            = (s[NS+i] > 0.0) ? s[NS+i] : DBL_EPSILON;
            dDeltaG_GKdT = - ssCnModelParameters[i*nc+0].entropy;
            for (j=0; j<NR; j++) dDeltaG_GKdT += (-ssCnModelParameters[i*nc+j+1].entropy + ssCnModelParameters[i*nc+0].entropy)*r[j];

            result[NS+i] += R*(log(y)-log(yIV)) + dDeltaG_GKdT + (f-1.0)*dgIntIVdt;
        }
    }
#endif /* USE_GHIORSO_KRESS_MODEL */

}

#ifdef USE_GHIORSO_KRESS_MODEL

static void fillD2GDSDP (double r[NR], double s[NT], double t, double p, double *result) {
    const double pr      = 1.0;
    const double tr      = 1673.15;
    const double y  = 0.3;
    double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b;
    double dvds[NS], d2vdsdt[NS], dcRefds[NS], d2cdsdt[NS], dmwds[NS], dcpds[NS], d3vdsdp2[NS], d4vdsdp3[NS], d5vdsdp4[NS], denom;
    int i, j;

    memset(result, '\0', (size_t) NT*sizeof(double));

    /* Convert input composition (r) to liquid moles (m)  */
    for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* Compute moles and total moles of oxides */
    for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
    if (mOxTot == 0.0) return;

    /* Deal with the special case of FeO1.3 */
    mOx[NA] = 0.0;
    if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
        const double y = 0.3;
        mOx[iOxFeO1_3] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
        }
    }

    for (i=0; i<NS; i++) {
        dvds[i]    = 0.0; d2vdsdt[i]  = 0.0; dcRefds[i]  = 0.0; d2cdsdt[i]  = 0.0; dmwds[i] = 0.0;
        dcpds[i]   = 0.0; d3vdsdp2[i] = 0.0; d4vdsdp3[i] = 0.0; d5vdsdp4[i] = 0.0;
    }

    for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
        double dmOxds[NS], dmOxTotds[NS];
        v        += mOx[i]*bulkSystem[i].gk_v;
        dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
        cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
        dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
        cp      += mOx[i]*bulkSystem[i].gk_cp;
        d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
        d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
        d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        mw      += mOx[i]*bulkSystem[i].mw;
        for (j=0; j<NS; j++) { dmOxds[j] = 0.0; dmOxTotds[j] = 0.0; }
        if (iCmpFe2SiO4_6 != -1) {
            for (j=0; j<NS; j++) {
                if      (iOxFe2O3  == i) dmOxds[j] += -2.0*y*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
                else if (iOxFeO       == i) dmOxds[j] += -2.0*(1.0-2.0*y)*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
                else if (iOxFeO1_3 == i) dmOxds[j] += 2.0*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
                dmOxTotds[j] += 2.0*y*s[iCmpFe2SiO4_6]*dnSpeciesds[j];
            }
            if      (iOxFe2O3  == i) dmOxds[iCmpFe2SiO4_6] += -2.0*y*nSpecies;
            else if (iOxFeO     == i) dmOxds[iCmpFe2SiO4_6] += -2.0*(1.0-2.0*y)*nSpecies;
            else if (iOxFeO1_3 == i) dmOxds[iCmpFe2SiO4_6] += 2.0*nSpecies;
            dmOxTotds[iCmpFe2SiO4_6] += 2.0*y*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            for (j=0; j<NS; j++) {
                if      (iOxFe2O3  == i) dmOxds[j] += -2.0*y*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
                else if (iOxFeO       == i) dmOxds[j] += -2.0*(1.0-2.0*y)*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
                else if (iOxFeO1_3 == i) dmOxds[j] += 2.0*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
                dmOxTotds[j] += 2.0*y*s[iCmpFe2AlO4_1]*dnSpeciesds[j];
            }
            if      (iOxFe2O3  == i) dmOxds[iCmpFe2AlO4_1] += -2.0*y*nSpecies;
            else if (iOxFeO     == i) dmOxds[iCmpFe2AlO4_1] += -2.0*(1.0-2.0*y)*nSpecies;
            else if (iOxFeO1_3 == i) dmOxds[iCmpFe2AlO4_1] += 2.0*nSpecies;
            dmOxTotds[iCmpFe2AlO4_1] += 2.0*y*nSpecies;
        }
        for (j=0; j<NS; j++) {
            dvds[j]       += dmOxds[j]*bulkSystem[i].gk_v;
            d2vdsdt[j]   += dmOxds[j]*bulkSystem[i].gk_dvdt;
            dmwds[j]     += dmOxds[j]*bulkSystem[i].mw;
            dcpds[j]     += dmOxds[j]*bulkSystem[i].gk_cp;
            dcRefds[j]   += dmOxds[j]*bulkSystem[i].gk_c/mOxTot - mOx[i]*bulkSystem[i].gk_c*dmOxTotds[j]/(mOxTot*mOxTot);
            dcRefds[j]   += (iOxAl2O3 != -1) ? dmOxds[j]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
              - 2.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotds[j]/(mOxTot*mOxTot*mOxTot) : 0.0;
            d2cdsdt[j]   += dmOxds[j]*bulkSystem[i].gk_dcdt/mOxTot - mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotds[j]/(mOxTot*mOxTot);
            d3vdsdp2[j]  += dmOxds[j]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d4vdsdp3[j]  += dmOxds[j]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d5vdsdp4[j]  += dmOxds[j]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        }
    }
    if (v == 0.0) return;

    alpha   = dvdt/v;
    v0      = v*exp(alpha*(t-tr));
    v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
    c      = cRef + (t-tr)*dcdt;
    v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
    v2      = d2vdp2;
    denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
    a      = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
    b      = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;

    if (1.0 + a*(p-pr) + b*(p-pr)*(p-pr) == 0.0) return;


    for (i=0; i<NS; i++) {
        double dalphads = d2vdsdt[i]/v - dvdt*dvds[i]/(v*v);
        double dv0ds    = (dvds[i] + v*dalphads*(t-tr))*exp(alpha*(t-tr));
        double dv1Refds = -2.0*v*dvds[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                    -v*v*(-1000.0*(cRef*dmwds[i]+2.0*mw*dcRefds[i])/(mw*mw*cRef*cRef*cRef)
                      + 2.0*tr*alpha*dalphads/(cp) - tr*alpha*alpha*dcpds[i]/(cp*cp));
        double dcds     = dcRefds[i] + (t-tr)*d2cdsdt[i];
        double dv1ds    = -2.0*v0*dv0ds*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                    -v0*v0*(-1000.0*(c*dmwds[i]+2.0*mw*dcds)/(mw*mw*c*c*c)
                        + 2.0*t*alpha*dalphads/(cp) - t*alpha*alpha*dcpds[i]/(cp*cp));
        double dv2ds    = d3vdsdp2[i];
        double ddenomds = 2.0*dv1Refds*d3vdp3+2.0*v1Ref*d4vdsdp3[i]-6.0*d2vdp2*d3vdsdp2[i];
        double dads     = (d3vdsdp2[i]*d3vdp3+d2vdp2*d4vdsdp3[i]-dv1Refds*d4vdp4/2.0-v1Ref*d5vdsdp4[i]/2.0)/denom
                - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomds;
        double dbds     = (d3vdsdp2[i]*d4vdp4/4.0+d2vdp2*d5vdsdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdsdp3[i])/denom
                - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomds;

        result[i] += (dv0ds + (dv1ds+dv0ds*a+v0*dads)*(p-pr) + (dv2ds/2.0+dv1ds*a+v1*dads+dv0ds*b+v0*dbds)*(p-pr)*(p-pr))
                 /(1.0 + a*(p-pr) + b*(p-pr)*(p-pr))
               - (v0 + (v1+v0*a)*(p-pr) + (v2/2.0+v1*a+v0*b)*(p-pr)*(p-pr))*(dads*(p-pr) + dbds*(p-pr)*(p-pr))
                    /((1.0 + a*(p-pr) + b*(p-pr)*(p-pr))*(1.0 + a*(p-pr) + b*(p-pr)*(p-pr)));
    }

    for (i=0; i<NY; i++) result[NS+i] = (fCN[i]-1.0)*((v0 + (v1+v0*a)*(p-pr) + (v2/2.0+v1*a+v0*b)*(p-pr)*(p-pr))/(1.0 + a*(p-pr) + b*(p-pr)*(p-pr)));

}

#else

static void fillD2GDSDP (double r[NR], double s[NT], double t, double p, double *result) {
    int i, j;

    memset(result, '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) {
        result[i] = vs[i] + vss[i][i]*s[i];
        for (j=0; j<NS; j++) result[i] += vss[i][j]*s[j];
        for (j=0; j<NR; j++) result[i] += vrs[j][i]*r[j];
    }
}

#endif /* USE_GHIORSO_KRESS_MODEL */

static void fillD2GDSDW (double r[NR], double s[NT], double t, double p, double result[NT][3*NP]) {
    int i, j, k, l, m, n, ii;

    /*******************************
   * Parameters: NW WH(), NE H() *
   *             NW WS(), NE S() *
   *             NW WV(), NE V() *
   *******************************/
    for (ii=0; ii<NT; ii++) memset(result[ii], '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
   * NW W parameters solution are first *
   **************************************/
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        for (ii=0; ii<NS; ii++) {
            for (j=0, m=0; j<NR; j++) {
                m += NR - j;
                m += ii;
                result[ii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*r[j];
                result[ii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*r[j];
#ifndef USE_GHIORSO_KRESS_MODEL
                result[ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*r[j];
#endif /* USE_GHIORSO_KRESS_MODEL */
                m += NS - ii;
            }
            result[ii][     n] +=         taylorCoeff[n+NE][1+NR+ii];
            result[ii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+ii];
#ifndef USE_GHIORSO_KRESS_MODEL
            result[ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+ii];
#endif /* USE_GHIORSO_KRESS_MODEL */

            for (k=ii; k<NS; k++) {
                m = ii*NS+(k+1)-(ii+1)*(ii+2)/2+(ii+1)-1+NR*(NR-1)/2+NR+NR*NS;
                if (taylorCoeff[n+NE][1+NR+NS+m] != 0.0) {
                    result[ii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*s[k];
                    result[ii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*s[k];
#ifndef USE_GHIORSO_KRESS_MODEL
                    result[ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*s[k];
#endif /* USE_GHIORSO_KRESS_MODEL */
                }
            }

            for (j=0; j<=ii; j++) {
                m = j*NS+(ii+1)-(j+1)*(j+2)/2+(j+1)-1+NR*(NR-1)/2+NR+NR*NS;
                if (taylorCoeff[n+NE][1+NR+NS+m] != 0.0) {
                    result[ii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*s[j];
                    result[ii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*s[j];
#ifndef USE_GHIORSO_KRESS_MODEL
                    result[ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*s[j];
#endif /* USE_GHIORSO_KRESS_MODEL */
                }
            }

        }
    }

    /**************************************
    * NE standard state terms are second *
    **************************************/
    for (i=0; i<NE; i++, n++) {
        for (ii=0; ii<NS; ii++) {
            result[ii][     n] +=         taylorCoeff[i][1+NR+ii];
            result[ii][  NP+n] +=      -t*taylorCoeff[i][1+NR+ii];
#ifndef USE_GHIORSO_KRESS_MODEL
            result[ii][2*NP+n] += (p-1.0)*taylorCoeff[i][1+NR+ii];
#endif /* USE_GHIORSO_KRESS_MODEL */
        }
    }
}

static double fillD2GDT2 (double r[NR], double s[NT], double t, double p) {
    double result;
    int i;

    /* Taylor expansion and standard state terms */
    result = CPconst;
    for (i=0; i<NR; i++) result += cpr[i]*r[i];
    for (i=0; i<NS; i++) result += cps[i]*s[i];
    result /= -t;

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double coeffCN;
        double d2gIntIVdt2 = d2IntegralV_GKdT2(r, s, t, p);

        for (i=0, coeffCN=1.0; i<NY; i++) coeffCN += (fCN[i]-1.0)*s[NS+i];
        result += coeffCN*d2gIntIVdt2;
    }
#endif /* USE_GHIORSO_KRESS_MODEL */

    return result;
}

static double fillD2GDTDP (double r[NR], double s[NT], double t, double p) {
    double result;
    int i;

    /* Taylor expansion and standard state terms */
    result = DVDTconst;
    for (i=0; i<NR; i++) result += dvdtr[i]*r[i];
    for (i=0; i<NS; i++) result += dvdts[i]*s[i];

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double coeffCN;
        double d2gIntIVdtdp = d2IntegralV_GKdTdP(r, s, t, p);

        for (i=0, coeffCN=1.0; i<NY; i++) coeffCN += (fCN[i]-1.0)*s[NS+i];
        result += coeffCN*d2gIntIVdtdp;
    }
#endif /* USE_GHIORSO_KRESS_MODEL */

    return result;
}

static void fillD2GDTDW (double r[NR], double s[NT], double t, double p, double result[3*NP]) {
    int i, j, k, l, m, n;

    /*******************************
    * Parameters: NW WH(), NE H() *
    *             NW WS(), NE S() *
    *             NW WV(), NE V() *
    *******************************/
    memset(result, '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
    * NW W parameters solution are first *
    **************************************/
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        result[NP+n] += -taylorCoeff[n+NE][0];
        m = 0;
        for (j=0, m=0; j<NR; j++) {
            result[NP+n] += -taylorCoeff[n+NE][1+j]*r[j];
            for (k=j; k<NR; k++, m++) result[NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*r[j]*r[k];
            for (k=0; k<NS; k++, m++) result[NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*r[j]*s[k];
        }
        for (j=0; j<NS; j++) {
            result[NP+n] += -taylorCoeff[n+NE][1+NR+j]*s[j];
            for (k=j; k<NS; k++, m++) result[NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*s[j]*s[k];
        }
    }

    /**************************************
    * NE standard state terms are second *
    **************************************/
    for (i=0; i<NE; i++, n++) {
        result[NP+n] += -taylorCoeff[i][0];
        for (j=0; j<NR; j++) result[NP+n] += -taylorCoeff[i][1+j]*r[j];
        for (j=0; j<NS; j++) result[NP+n] += -taylorCoeff[i][1+NR+j]*s[j];
    }
}

static void fillD2GDPDW (double r[NR], double s[NT], double t, double p, double *result) {

    /*******************************
    * Parameters: NW WH(), NE H() *
    *             NW WS(), NE S() *
    *             NW WV(), NE V() *
    *******************************/
    memset(result, '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
   * NW W parameters solution are first *
   **************************************/
#ifndef USE_GHIORSO_KRESS_MODEL
{
    int i, j, k, l, m, n;
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        result[2*NP+n] += taylorCoeff[n+NE][0];
        m = 0;
        for (j=0, m=0; j<NR; j++) {
            result[2*NP+n] += taylorCoeff[n+NE][1+j]*r[j];
            for (k=j; k<NR; k++, m++) result[2*NP+n] += taylorCoeff[n+NE][1+NR+NS+m]*r[j]*r[k];
            for (k=0; k<NS; k++, m++) result[2*NP+n] += taylorCoeff[n+NE][1+NR+NS+m]*r[j]*s[k];
        }
        for (j=0; j<NS; j++) {
            result[2*NP+n] += taylorCoeff[n+NE][1+NR+j]*s[j];
            for (k=j; k<NS; k++, m++) result[2*NP+n] += taylorCoeff[n+NE][1+NR+NS+m]*s[j]*s[k];
        }
    }

    /**************************************
    * NE standard state terms are second *
    **************************************/
    for (i=0; i<NE; i++, n++) {
        result[2*NP+n] += taylorCoeff[i][0];
        for (j=0; j<NR; j++) result[2*NP+n] += taylorCoeff[i][1+j]*r[j];
        for (j=0; j<NS; j++) result[2*NP+n] += taylorCoeff[i][1+NR+j]*s[j];
    }
}
#endif /* USE_GHIORSO_KRESS_MODEL */
}

static double fillD2GDP2 (double r[NR], double s[NT], double t, double p) {
    double result;
    int i;

    /* Taylor expansion and standard state terms */
    result = DVDPconst;
    for (i=0; i<NR; i++) result += dvdpr[i]*r[i];
    for (i=0; i<NS; i++) result += dvdps[i]*s[i];

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double coeffCN;
        double d2gIntIVdp2 = d2IntegralV_GKdP2(r, s, t, p);

        for (i=0, coeffCN=1.0; i<NY; i++) coeffCN += (fCN[i]-1.0)*s[NS+i];
        result += coeffCN*d2gIntIVdp2;
    }
#endif /* USE_GHIORSO_KRESS_MODEL */

    return result;
}

static void fillD3GDR3 (double r[NR], double s[NT], double t, double p, double result[NR][NR][NR]) {
    int i, j, k, l;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) for (j=0; j<NR; j++) memset(result[i][j], '\0', (size_t) NR*sizeof(double));

    /* Configurational entropy terms */
    for (j=0; j<NR; j++) {
        for (k=j; k<NR; k++) {
            for (l=k; l<NR; l++) {
                double config = 0.0;
                for (i=0; i<NE; i++) if (xSpecies[i] > 0.0 && dxSpeciesdr[i][j] != 0.0 && dxSpeciesdr[i][k] != 0.0 && dxSpeciesdr[i][l] != 0.0)
                    config += -dxSpeciesdr[i][j]*dxSpeciesdr[i][k]*dxSpeciesdr[i][l]/(xSpecies[i]*xSpecies[i]);
                if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0 && dxSpeciesdr[nH2O][j] != 0.0 && dxSpeciesdr[nH2O][k] != 0.0 && dxSpeciesdr[i][l] != 0.0)
                    config += -dxSpeciesdr[nH2O][j]*dxSpeciesdr[nH2O][k]*dxSpeciesdr[nH2O][l]
                        *(1.0/(xSpecies[nH2O]*xSpecies[nH2O]) - 1.0/((1.0-xSpecies[nH2O])*(1.0-xSpecies[nH2O])));
                result[j][k][l] += R*t*nSpecies*config;
                result[k][j][l]  = result[j][k][l];
                result[l][j][k]  = result[j][k][l];
                result[l][k][j]  = result[j][k][l];
                result[j][l][k]  = result[j][k][l];
                result[k][l][j]  = result[j][k][l];
            }
        }
    }

    ternaryH2OCO2terms(FOURTH, r);
    for (j=0; j<NR; j++) {
        for (k=j; k<NR; k++) {
            for (l=k; l<NR; l++) {
                result[j][k][l] += d3gdr3T[j][k][l];
                result[k][j][l]  = result[j][k][l];
                result[l][j][k]  = result[j][k][l];
                result[l][k][j]  = result[j][k][l];
                result[j][l][k]  = result[j][k][l];
                result[k][l][j]  = result[j][k][l];
            }
        }
    }

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        const double pr         = 1.0;
        const double tr         = 1673.15;
        double m[NA], mOx[NA], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum, coeffCN;
        double dvdr[NR], d2vdrdt[NR], dcRefdr[NR], d2cdrdt[NR], dmwdr[NR], dcpdr[NR], d3vdrdp2[NR], d4vdrdp3[NR], d5vdrdp4[NR], dmOxTotdr[NR], denom;
        double d2cRefdr2[NR][NR], d3cdr2dt[NR][NR], d3cRefdr3[NR][NR][NR], d4cdr3dt[NR][NR][NR];

        if (fabs(p-pr) < 10.0*DBL_EPSILON) return;

        /* Convert input composition (r) to liquid moles (m)  */
        for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

        /* Compute moles and total moles of oxides */
        for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
        if (mOxTot == 0.0) return;

        for (i=0; i<NR; i++) {
            dvdr[i]    = 0.0; d2vdrdt[i]  = 0.0; dcRefdr[i]  = 0.0; d2cdrdt[i]  = 0.0; dmwdr[i] = 0.0;
            dcpdr[i]   = 0.0; d3vdrdp2[i] = 0.0; d4vdrdp3[i] = 0.0; d5vdrdp4[i] = 0.0;
            for (j=0, dmOxTotdr[i]=0.0; j<nc; j++) dmOxTotdr[i] += (liquid[i+1].liqToOx)[j] - (liquid[0].liqToOx)[j];
            for (j=i; j<NR; j++) {
                d2cRefdr2[i][j] = 0.0; d3cdr2dt[i][j] = 0.0;
  for (k=j; k<NR; k++) { d3cRefdr3[i][j][k] = 0.0; d4cdr3dt[i][j][k]  = 0.0; }
            }
        }

        for (i=0, coeffCN=1.0; i<NY; i++) coeffCN += (fCN[i]-1.0)*s[NS+i];

        for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<nc; i++) {
            v       += mOx[i]*bulkSystem[i].gk_v;
            dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
            cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
            dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
            cp      += mOx[i]*bulkSystem[i].gk_cp;
            d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
            mw      += mOx[i]*bulkSystem[i].mw;

            for (j=0; j<NR; j++) {
                dvdr[j]      += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_v;
                d2vdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dvdt;
                dmwdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].mw;
                dcpdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_cp;
  dcRefdr[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c/mOxTot
                - mOx[i]*bulkSystem[i].gk_c*dmOxTotdr[j]/(mOxTot*mOxTot);
                dcRefdr[j]   += (iOxAl2O3 != -1) ? ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
                - 2.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot) : 0.0;
  d2cdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt/mOxTot
                - mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotdr[j]/(mOxTot*mOxTot);
                d3vdrdp2[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
                d4vdrdp3[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
                d5vdrdp4[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
            }
            if (iCmpAl2O3 != -1) dcRefdr[iCmpAl2O3] += ((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*mOx[i]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot);

            for (j=0; j<NR; j++) {
  for (k=j; k<NR; k++) {
    d2cRefdr2[j][k] += -((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c*dmOxTotdr[k]/(mOxTot*mOxTot)
                - ((liquid[k+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c*dmOxTotdr[j]/(mOxTot*mOxTot)
            + 2.0*mOx[i]*bulkSystem[i].gk_c*dmOxTotdr[j]*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot);
    d2cRefdr2[j][k] += (iOxAl2O3 != -1) ?
                - 2.0*((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot)
                - 2.0*((liquid[k+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot)
            + 6.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot*mOxTot) : 0.0;

    d3cdr2dt[j][k]  += -((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt*dmOxTotdr[k]/(mOxTot*mOxTot)
                - ((liquid[k+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt*dmOxTotdr[j]/(mOxTot*mOxTot)
            + 2.0*mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotdr[j]*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot);
  }
  if (iCmpAl2O3 != -1) d2cRefdr2[(j<iCmpAl2O3) ? j : iCmpAl2O3][(j>iCmpAl2O3) ? j : iCmpAl2O3] +=
    ((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
    - 2.0*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*mOx[i]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot);
            }
            if (iCmpAl2O3 != -1) d2cRefdr2[iCmpAl2O3][iCmpAl2O3] +=
  ((liquid[iCmpAl2O3+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
  - 2.0*mOx[i]*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*bulkSystem[i].gk_cXal2o3*dmOxTotdr[iCmpAl2O3]/(mOxTot*mOxTot*mOxTot);

            for (j=0; j<NR; j++) {
  for (k=j; k<NR; k++) {
    for (l=k; l<NR; l++) {
      d3cRefdr3[j][k][l] += 2.0*((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c*dmOxTotdr[k]*dmOxTotdr[l]/(mOxTot*mOxTot*mOxTot)
                    + 2.0*((liquid[k+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c*dmOxTotdr[j]*dmOxTotdr[l]/(mOxTot*mOxTot*mOxTot)
                + 2.0*((liquid[l+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c*dmOxTotdr[j]*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot)
                - 6.0*mOx[i]*bulkSystem[i].gk_c*dmOxTotdr[j]*dmOxTotdr[k]*dmOxTotdr[l]/(mOxTot*mOxTot*mOxTot*mOxTot);
      d3cRefdr3[j][k][l] += (iOxAl2O3 != -1) ?
                    + 6.0*((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[k]*dmOxTotdr[l]/(mOxTot*mOxTot*mOxTot*mOxTot)
                    + 6.0*((liquid[k+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]*dmOxTotdr[l]/(mOxTot*mOxTot*mOxTot*mOxTot)
                + 6.0*((liquid[l+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot*mOxTot)
                - 24.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]*dmOxTotdr[k]*dmOxTotdr[l]/(mOxTot*mOxTot*mOxTot*mOxTot*mOxTot) : 0.0;

      d4cdr3dt[j][k][l] += 2.0*((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt*dmOxTotdr[k]*dmOxTotdr[l]/(mOxTot*mOxTot*mOxTot)
                        + 2.0*((liquid[k+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt*dmOxTotdr[j]*dmOxTotdr[l]/(mOxTot*mOxTot*mOxTot)
                    + 2.0*((liquid[l+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt*dmOxTotdr[j]*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot)
                    - 6.0*mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotdr[j]*dmOxTotdr[k]*dmOxTotdr[l]/(mOxTot*mOxTot*mOxTot*mOxTot);
    }
    /* case for j || k || l == iCmpAl2O3 */
    if (iCmpAl2O3 != -1) d3cRefdr3[(j<iCmpAl2O3) ? j : iCmpAl2O3][(j>iCmpAl2O3) ? j : ((k<iCmpAl2O3) ? k : iCmpAl2O3)][(k>iCmpAl2O3) ? k : iCmpAl2O3] +=
      - 2.0*((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*bulkSystem[i].gk_cXal2o3*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot)
      - 2.0*((liquid[k+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot)
      + 6.0*mOx[i]*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot*mOxTot);
  }
  /* case for j,k || j,l || k,l == iCmpAl2O3 */
  if (iCmpAl2O3 != -1) d3cRefdr3[(j<iCmpAl2O3) ? j : iCmpAl2O3][iCmpAl2O3][(j>iCmpAl2O3) ? j : iCmpAl2O3] +=
    - 2.0*((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*bulkSystem[i].gk_cXal2o3*dmOxTotdr[iCmpAl2O3]/(mOxTot*mOxTot*mOxTot)
    - 2.0*((liquid[iCmpAl2O3+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot)
    + 6.0*mOx[i]*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]*dmOxTotdr[iCmpAl2O3]/(mOxTot*mOxTot*mOxTot*mOxTot);
            }
            /* Case for j,k,l == iCmpAl2O3*/
            if (iCmpAl2O3 != -1) d3cRefdr3[iCmpAl2O3][iCmpAl2O3][iCmpAl2O3] +=
  - 2.0*((liquid[iCmpAl2O3+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*bulkSystem[i].gk_cXal2o3*dmOxTotdr[iCmpAl2O3]/(mOxTot*mOxTot*mOxTot)
  - 2.0*((liquid[iCmpAl2O3+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*bulkSystem[i].gk_cXal2o3*dmOxTotdr[iCmpAl2O3]/(mOxTot*mOxTot*mOxTot)
  + 6.0*mOx[i]*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*bulkSystem[i].gk_cXal2o3*dmOxTotdr[iCmpAl2O3]*dmOxTotdr[iCmpAl2O3]/(mOxTot*mOxTot*mOxTot*mOxTot);
        }
        if (v == 0.0) return;

        alpha   = dvdt/v;
        v0      = v*exp(alpha*(t-tr));
        v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
        c       = cRef + (t-tr)*dcdt;
        v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
        v2      = d2vdp2;
        denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
        a        = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
        b        = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;
        sum     = a*a - 4.0*b;

        for (i=0; i<NR; i++) {
            double dalphadri   = d2vdrdt[i]/v - dvdt*dvdr[i]/(v*v);
            double dv0dri      = (dvdr[i] + v*dalphadri*(t-tr))*exp(alpha*(t-tr));
            double dv1Refdri   = -2.0*v*dvdr[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                           -v*v*(- 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
                        + 2.0*tr*alpha*dalphadri/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp));
            double dcdri       = dcRefdr[i] + (t-tr)*d2cdrdt[i];
            double dv1dri      = -2.0*v0*dv0dri*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                           -v0*v0*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)/(mw*mw*c*c*c)
                            + 2.0*t*alpha*dalphadri/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp));
            double ddenomdri   = 2.0*dv1Refdri*d3vdp3+2.0*v1Ref*d4vdrdp3[i]-6.0*d2vdp2*d3vdrdp2[i];
            double dadri       = (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdri*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)/denom
                         - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdri;
            double dbdri       = (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])/denom
                         - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdri;

            for (j=i; j<NR; j++) {
                double dalphadrj   = d2vdrdt[j]/v - dvdt*dvdr[j]/(v*v);
                double dv0drj       = (dvdr[j] + v*dalphadrj*(t-tr))*exp(alpha*(t-tr));
                double dv1Refdrj   = -2.0*v*dvdr[j]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                 -v*v*(- 1000.0*(cRef*dmwdr[j]+2.0*mw*dcRefdr[j])/(mw*mw*cRef*cRef*cRef)
                            + 2.0*tr*alpha*dalphadrj/(cp) - tr*alpha*alpha*dcpdr[j]/(cp*cp));
                double dcdrj       = dcRefdr[j] + (t-tr)*d2cdrdt[j];
                double dv1drj       = -2.0*v0*dv0drj*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                 -v0*v0*(- 1000.0*(c*dmwdr[j]+2.0*mw*dcdrj)/(mw*mw*c*c*c)
                                + 2.0*t*alpha*dalphadrj/(cp) - t*alpha*alpha*dcpdr[j]/(cp*cp));
                double ddenomdrj   = 2.0*dv1Refdrj*d3vdp3+2.0*v1Ref*d4vdrdp3[j]-6.0*d2vdp2*d3vdrdp2[j];
                double dadrj       = (d3vdrdp2[j]*d3vdp3+d2vdp2*d4vdrdp3[j]-dv1Refdrj*d4vdp4/2.0-v1Ref*d5vdrdp4[j]/2.0)/denom
               - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdrj;
                double dbdrj       = (d3vdrdp2[j]*d4vdp4/4.0+d2vdp2*d5vdrdp4[j]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[j])/denom
                         - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdrj;

                double d2alphadridrj = -d2vdrdt[i]*dvdr[j]/(v*v)-d2vdrdt[j]*dvdr[i]/(v*v)+2.0*dvdt*dvdr[i]*dvdr[j]/(v*v*v);
                double d2v0dridrj    = dvdr[i]*dalphadrj*(t-tr)*exp(alpha*(t-tr))
                + dvdr[j]*dalphadri*(t-tr)*exp(alpha*(t-tr))
                + v*d2alphadridrj*(t-tr)*exp(alpha*(t-tr))
                + v*dalphadri*pow(t-tr,2.0)*dalphadrj*exp(alpha*(t-tr));
                double d2v1Refdridrj = -2.0*dvdr[j]*dvdr[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                                            - 2.0*v*dvdr[i]*(- 1000.0*(cRef*dmwdr[j]+2.0*mw*dcRefdr[j])/(mw*mw*cRef*cRef*cRef)
                            + 2.0*tr*alpha*dalphadrj/(cp) - tr*alpha*alpha*dcpdr[j]/(cp*cp))
                 - 2.0*v*dvdr[j]*(- 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
                            + 2.0*tr*alpha*dalphadri/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp))
                - v*v*(- 1000.0*(dcRefdr[j]*dmwdr[i]+2.0*dmwdr[j]*dcRefdr[i]+2.0*mw*d2cRefdr2[i][j])/(mw*mw*cRef*cRef*cRef)
                            + 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])*(2.0*dmwdr[j]*cRef+3.0*mw*dcRefdr[j])/(mw*mw*mw*cRef*cRef*cRef*cRef)
                            + 2.0*tr*alpha*d2alphadridrj/(cp)
           + 2.0*tr*dalphadrj*dalphadri/(cp) - 2.0*tr*alpha*dalphadri*dcpdr[j]/(cp*cp)
           - 2.0*tr*alpha*dalphadrj*dcpdr[i]/(cp*cp) + 2.0*tr*alpha*alpha*dcpdr[i]*dcpdr[j]/(cp*cp*cp));
                double d2cdridrj       = d2cRefdr2[i][j] + (t-tr)*d3cdr2dt[i][j];
                double d2v1dridrj       = -2.0*dv0drj*dv0dri*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                      - 2.0*v0*d2v0dridrj*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                      - 2.0*v0*dv0dri*(- 1000.0*(c*dmwdr[j]+2.0*mw*dcdrj)/(mw*mw*c*c*c)
                                                + 2.0*t*alpha*dalphadrj/(cp) - t*alpha*alpha*dcpdr[j]/(cp*cp))
                            - 2.0*v0*dv0drj*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)/(mw*mw*c*c*c)
                                                + 2.0*t*alpha*dalphadri/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp))
                -v0*v0*(- 1000.0*(dcdrj*dmwdr[i]+2.0*dmwdr[j]*dcdri+2.0*mw*d2cdridrj)/(mw*mw*c*c*c)
                                + 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)*(2.0*dmwdr[j]*c+3.0*mw*dcdrj)/(mw*mw*mw*c*c*c*c)
             + 2.0*t*alpha*d2alphadridrj/(cp)
                                + 2.0*t*dalphadrj*dalphadri/(cp) - 2.0*t*alpha*dalphadri*dcpdr[j]/(cp*cp)
             - 2.0*t*alpha*dalphadrj*dcpdr[i]/(cp*cp) + 2.0*t*alpha*alpha*dcpdr[i]*dcpdr[j]/(cp*cp*cp));
                double d2denomdridrj   = 2.0*d2v1Refdridrj*d3vdp3 + 2.0*dv1Refdri*d4vdrdp3[j] + 2.0*dv1Refdrj*d4vdrdp3[i] - 6.0*d3vdrdp2[j]*d3vdrdp2[i];
                double d2adridrj       = (d3vdrdp2[i]*d4vdrdp3[j]+d3vdrdp2[j]*d4vdrdp3[i]-d2v1Refdridrj*d4vdp4/2.0-dv1Refdri*d5vdrdp4[j]/2.0-dv1Refdrj*d5vdrdp4[i]/2.0)/denom
                                        - (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdri*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)*ddenomdrj/(denom*denom)
               - (d3vdrdp2[j]*d3vdp3+d2vdp2*d4vdrdp3[j]-dv1Refdrj*d4vdp4/2.0-v1Ref*d5vdrdp4[j]/2.0)*ddenomdri/(denom*denom)
                           + 2.0*(d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*ddenomdri*ddenomdrj/(denom*denom*denom)
                           - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*d2denomdridrj/(denom*denom);
                double d2bdridrj       = (d3vdrdp2[i]*d5vdrdp4[j]/4.0+d3vdrdp2[j]*d5vdrdp4[i]/4.0-2.0/3.0*d4vdrdp3[j]*d4vdrdp3[i])/denom
                           - (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])*ddenomdrj/(denom*denom)
               - (d3vdrdp2[j]*d4vdp4/4.0+d2vdp2*d5vdrdp4[j]/4.0-2.0*d4vdrdp3[j]*d3vdp3/3.0)*ddenomdri/(denom*denom)
            + 2.0*(d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)*ddenomdri*ddenomdrj/(denom*denom*denom)
            - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)*d2denomdridrj/(denom*denom);

                for (k=j; k<NR; k++) {
                    double dalphadrk   = d2vdrdt[k]/v - dvdt*dvdr[k]/(v*v);
                    double dv0drk      = (dvdr[k] + v*dalphadrk*(t-tr))*exp(alpha*(t-tr));
                    double dv1Refdrk   = -2.0*v*dvdr[k]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                                    -v*v*(- 1000.0*(cRef*dmwdr[k]+2.0*mw*dcRefdr[k])/(mw*mw*cRef*cRef*cRef)
                    + 2.0*tr*alpha*dalphadrk/(cp) - tr*alpha*alpha*dcpdr[k]/(cp*cp));
                    double dcdrk       = dcRefdr[k] + (t-tr)*d2cdrdt[k];
                    double dv1drk      = -2.0*v0*dv0drk*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                                    -v0*v0*(- 1000.0*(c*dmwdr[k]+2.0*mw*dcdrk)/(mw*mw*c*c*c)
                      + 2.0*t*alpha*dalphadrk/(cp) - t*alpha*alpha*dcpdr[k]/(cp*cp));
                    double ddenomdrk   = 2.0*dv1Refdrk*d3vdp3+2.0*v1Ref*d4vdrdp3[k]-6.0*d2vdp2*d3vdrdp2[k];
                    double dadrk       = (d3vdrdp2[k]*d3vdp3+d2vdp2*d4vdrdp3[k]-dv1Refdrk*d4vdp4/2.0-v1Ref*d5vdrdp4[k]/2.0)/denom
                                - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdrk;
                    double dbdrk       = (d3vdrdp2[k]*d4vdp4/4.0+d2vdp2*d5vdrdp4[k]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[k])/denom
                             - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdrk;

                    double d2alphadridrk = -d2vdrdt[i]*dvdr[k]/(v*v)-d2vdrdt[k]*dvdr[i]/(v*v)+2.0*dvdt*dvdr[i]*dvdr[k]/(v*v*v);
                    double d2v0dridrk    = dvdr[i]*dalphadrk*(t-tr)*exp(alpha*(t-tr))
              + dvdr[k]*dalphadri*(t-tr)*exp(alpha*(t-tr))
              + v*d2alphadridrk*(t-tr)*exp(alpha*(t-tr))
              + v*dalphadri*pow(t-tr,2.0)*dalphadrk*exp(alpha*(t-tr));
                    double d2v1Refdridrk = -2.0*dvdr[k]*dvdr[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
              - 2.0*v*dvdr[i]*(- 1000.0*(cRef*dmwdr[k]+2.0*mw*dcRefdr[k])/(mw*mw*cRef*cRef*cRef)
             + 2.0*tr*alpha*dalphadrk/(cp) - tr*alpha*alpha*dcpdr[k]/(cp*cp))
                                    - 2.0*v*dvdr[k]*(- 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
                            + 2.0*tr*alpha*dalphadri/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp))
              - v*v*(- 1000.0*(dcRefdr[k]*dmwdr[i]+2.0*dmwdr[k]*dcRefdr[i]+2.0*mw*d2cRefdr2[i][k])/(mw*mw*cRef*cRef*cRef)
               + 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])*(2.0*dmwdr[k]*cRef+3.0*mw*dcRefdr[k])/(mw*mw*mw*cRef*cRef*cRef*cRef)
               + 2.0*tr*alpha*d2alphadridrk/(cp)
               + 2.0*tr*dalphadrk*dalphadri/(cp) - 2.0*tr*alpha*dalphadri*dcpdr[k]/(cp*cp)
               - 2.0*tr*alpha*dalphadrk*dcpdr[i]/(cp*cp) + 2.0*tr*alpha*alpha*dcpdr[i]*dcpdr[k]/(cp*cp*cp));
                    double d2cdridrk          = d2cRefdr2[i][k] + (t-tr)*d3cdr2dt[i][k];
                    double d2v1dridrk          = -2.0*dv0drk*dv0dri*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                    - 2.0*v0*d2v0dridrk*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                    - 2.0*v0*dv0dri*(- 1000.0*(c*dmwdr[k]+2.0*mw*dcdrk)/(mw*mw*c*c*c)
                  + 2.0*t*alpha*dalphadrk/(cp) - t*alpha*alpha*dcpdr[k]/(cp*cp))
                   - 2.0*v0*dv0drk*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)/(mw*mw*c*c*c)
                                        + 2.0*t*alpha*dalphadri/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp))
              -v0*v0*(- 1000.0*(dcdrk*dmwdr[i]+2.0*dmwdr[k]*dcdri+2.0*mw*d2cdridrk)/(mw*mw*c*c*c)
                        + 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)*(2.0*dmwdr[k]*c+3.0*mw*dcdrk)/(mw*mw*mw*c*c*c*c)
                        + 2.0*t*alpha*d2alphadridrk/(cp)
                        + 2.0*t*dalphadrk*dalphadri/(cp) - 2.0*t*alpha*dalphadri*dcpdr[k]/(cp*cp)
                        - 2.0*t*alpha*dalphadrk*dcpdr[i]/(cp*cp) + 2.0*t*alpha*alpha*dcpdr[i]*dcpdr[k]/(cp*cp*cp));
                    double d2denomdridrk = 2.0*d2v1Refdridrk*d3vdp3 + 2.0*dv1Refdri*d4vdrdp3[k] + 2.0*dv1Refdrk*d4vdrdp3[i] - 6.0*d3vdrdp2[k]*d3vdrdp2[i];
                    double d2adridrk     = (d3vdrdp2[i]*d4vdrdp3[k]+d3vdrdp2[k]*d4vdrdp3[i]-d2v1Refdridrk*d4vdp4/2.0-dv1Refdri*d5vdrdp4[k]/2.0-dv1Refdrk*d5vdrdp4[i]/2.0)/denom
              - (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdri*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)*ddenomdrk/(denom*denom)
                                    - (d3vdrdp2[k]*d3vdp3+d2vdp2*d4vdrdp3[k]-dv1Refdrk*d4vdp4/2.0-v1Ref*d5vdrdp4[k]/2.0)*ddenomdri/(denom*denom)
                                    + 2.0*(d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*ddenomdri*ddenomdrk/(denom*denom*denom)
                                    - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*d2denomdridrk/(denom*denom);
                    double d2bdridrk     = (d3vdrdp2[i]*d5vdrdp4[k]/4.0+d3vdrdp2[k]*d5vdrdp4[i]/4.0-2.0/3.0*d4vdrdp3[k]*d4vdrdp3[i])/denom
                                    - (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])*ddenomdrk/(denom*denom)
                                    - (d3vdrdp2[k]*d4vdp4/4.0+d2vdp2*d5vdrdp4[k]/4.0-2.0*d4vdrdp3[k]*d3vdp3/3.0)*ddenomdri/(denom*denom)
              + 2.0*(d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)*ddenomdri*ddenomdrk/(denom*denom*denom)
              - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)*d2denomdridrk/(denom*denom);

                    double d2alphadrjdrk = -d2vdrdt[j]*dvdr[k]/(v*v)-d2vdrdt[k]*dvdr[j]/(v*v)+2.0*dvdt*dvdr[j]*dvdr[k]/(v*v*v);
                    double d2v0drjdrk    = dvdr[j]*dalphadrk*(t-tr)*exp(alpha*(t-tr))
              + dvdr[k]*dalphadrj*(t-tr)*exp(alpha*(t-tr))
              + v*d2alphadrjdrk*(t-tr)*exp(alpha*(t-tr))
              + v*dalphadrj*pow(t-tr,2.0)*dalphadrk*exp(alpha*(t-tr));
                    double d2v1Refdrjdrk = -2.0*dvdr[k]*dvdr[j]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
              - 2.0*v*dvdr[j]*(- 1000.0*(cRef*dmwdr[k]+2.0*mw*dcRefdr[k])/(mw*mw*cRef*cRef*cRef)
             + 2.0*tr*alpha*dalphadrk/(cp) - tr*alpha*alpha*dcpdr[k]/(cp*cp))
                                    - 2.0*v*dvdr[k]*(- 1000.0*(cRef*dmwdr[j]+2.0*mw*dcRefdr[j])/(mw*mw*cRef*cRef*cRef)
                            + 2.0*tr*alpha*dalphadrj/(cp) - tr*alpha*alpha*dcpdr[j]/(cp*cp))
              - v*v*(- 1000.0*(dcRefdr[k]*dmwdr[j]+2.0*dmwdr[k]*dcRefdr[j]+2.0*mw*d2cRefdr2[j][k])/(mw*mw*cRef*cRef*cRef)
               + 1000.0*(cRef*dmwdr[j]+2.0*mw*dcRefdr[j])*(2.0*dmwdr[k]*cRef+3.0*mw*dcRefdr[k])/(mw*mw*mw*cRef*cRef*cRef*cRef)
               + 2.0*tr*alpha*d2alphadrjdrk/(cp)
               + 2.0*tr*dalphadrk*dalphadrj/(cp) - 2.0*tr*alpha*dalphadrj*dcpdr[k]/(cp*cp)
               - 2.0*tr*alpha*dalphadrk*dcpdr[j]/(cp*cp) + 2.0*tr*alpha*alpha*dcpdr[j]*dcpdr[k]/(cp*cp*cp));
                    double d2cdrjdrk     = d2cRefdr2[j][k] + (t-tr)*d3cdr2dt[j][k];
                    double d2v1drjdrk    = -2.0*dv0drk*dv0drj*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
              - 2.0*v0*d2v0drjdrk*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
              - 2.0*v0*dv0drj*(- 1000.0*(c*dmwdr[k]+2.0*mw*dcdrk)/(mw*mw*c*c*c)
                    + 2.0*t*alpha*dalphadrk/(cp) - t*alpha*alpha*dcpdr[k]/(cp*cp))
                                    - 2.0*v0*dv0drk*(- 1000.0*(c*dmwdr[j]+2.0*mw*dcdrj)/(mw*mw*c*c*c)
                         + 2.0*t*alpha*dalphadrj/(cp) - t*alpha*alpha*dcpdr[j]/(cp*cp))
              -v0*v0*(- 1000.0*(dcdrk*dmwdr[j]+2.0*dmwdr[k]*dcdrj+2.0*mw*d2cdrjdrk)/(mw*mw*c*c*c)
                        + 1000.0*(c*dmwdr[j]+2.0*mw*dcdrj)*(2.0*dmwdr[k]*c+3.0*mw*dcdrk)/(mw*mw*mw*c*c*c*c)
                        + 2.0*t*alpha*d2alphadrjdrk/(cp)
                        + 2.0*t*dalphadrk*dalphadrj/(cp) - 2.0*t*alpha*dalphadrj*dcpdr[k]/(cp*cp)
                        - 2.0*t*alpha*dalphadrk*dcpdr[j]/(cp*cp) + 2.0*t*alpha*alpha*dcpdr[j]*dcpdr[k]/(cp*cp*cp));
                    double d2denomdrjdrk = 2.0*d2v1Refdrjdrk*d3vdp3 + 2.0*dv1Refdrj*d4vdrdp3[k] + 2.0*dv1Refdrk*d4vdrdp3[j] - 6.0*d3vdrdp2[k]*d3vdrdp2[j];
                    double d2adrjdrk     = (d3vdrdp2[j]*d4vdrdp3[k]+d3vdrdp2[k]*d4vdrdp3[j]-d2v1Refdrjdrk*d4vdp4/2.0-dv1Refdrj*d5vdrdp4[k]/2.0-dv1Refdrk*d5vdrdp4[j]/2.0)/denom
              - (d3vdrdp2[j]*d3vdp3+d2vdp2*d4vdrdp3[j]-dv1Refdrj*d4vdp4/2.0-v1Ref*d5vdrdp4[j]/2.0)*ddenomdrk/(denom*denom)
                                    - (d3vdrdp2[k]*d3vdp3+d2vdp2*d4vdrdp3[k]-dv1Refdrk*d4vdp4/2.0-v1Ref*d5vdrdp4[k]/2.0)*ddenomdrj/(denom*denom)
                                    + 2.0*(d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*ddenomdrj*ddenomdrk/(denom*denom*denom)
                                    - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*d2denomdrjdrk/(denom*denom);
                    double d2bdrjdrk     = (d3vdrdp2[j]*d5vdrdp4[k]/4.0+d3vdrdp2[k]*d5vdrdp4[j]/4.0-2.0/3.0*d4vdrdp3[k]*d4vdrdp3[j])/denom
                                    - (d3vdrdp2[j]*d4vdp4/4.0+d2vdp2*d5vdrdp4[j]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[j])*ddenomdrk/(denom*denom)
                                    - (d3vdrdp2[k]*d4vdp4/4.0+d2vdp2*d5vdrdp4[k]/4.0-2.0*d4vdrdp3[k]*d3vdp3/3.0)*ddenomdrj/(denom*denom)
              + 2.0*(d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)*ddenomdrj*ddenomdrk/(denom*denom*denom)
              - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)*d2denomdrjdrk/(denom*denom);


                    double d3alphadr3 = 2.0*d2vdrdt[i]/(v*v*v)*dvdr[j]*dvdr[k]
                            + 2.0*d2vdrdt[j]/(v*v*v)*dvdr[i]*dvdr[k]
                            + 2.0*d2vdrdt[k]/(v*v*v)*dvdr[i]*dvdr[j]
                            - 6.0*dvdt/(v*v*v*v)*dvdr[i]*dvdr[j]*dvdr[k];
                    double d3v0dr3 = dvdr[i]*d2alphadrjdrk*(t-tr)*exp(alpha*(t-tr))
             + dvdr[i]*dalphadrj*pow(t-tr,2.0)*dalphadrk*exp(alpha*(t-tr))
             + dvdr[j]*d2alphadridrk*(t-tr)*exp(alpha*(t-tr))
             + dvdr[j]*dalphadri*pow(t-tr,2.0)*dalphadrk*exp(alpha*(t-tr))
             + dvdr[k]*d2alphadridrj*(t-tr)*exp(alpha*(t-tr))
             + v*d3alphadr3*(t-tr)*exp(alpha*(t-tr))
             + v*d2alphadridrj*pow(t-tr,2.0)*dalphadrk*exp(alpha*(t-tr))
             + dvdr[k]*dalphadri*pow(t-tr,2.0)*dalphadrj*exp(alpha*(t-tr))
             + v*d2alphadridrk*pow(t-tr,2.0)*dalphadrj*exp(alpha*(t-tr))
             + v*dalphadri*pow(t-tr,2.0)*d2alphadrjdrk*exp(alpha*(t-tr))
             + v*dalphadri*pow(t-tr,3.0)*dalphadrj*dalphadrk*exp(alpha*(t-tr));
    double d3v1Refdr3 = -2.0*dvdr[j]*dvdr[i]*(-1000.0*dmwdr[k]/(mw*mw*cRef*cRef) - 2000.0*dcRefdr[k]/(mw*cRef*cRef*cRef) + 2.0*tr*alpha*dalphadrk/(cp))
                                            - 2.0*dvdr[k]*dvdr[i]*(- 1000.0*(cRef*dmwdr[j]+2.0*mw*dcRefdr[j])/(mw*mw*cRef*cRef*cRef)
                                  + 2.0*tr*alpha*dalphadrj/(cp) - tr*alpha*alpha*dcpdr[j]/(cp*cp))
                - 2.0*v*dvdr[i]*(- 1000.0*(dcRefdr[k]*dmwdr[j]+2.0*dmwdr[k]*dcRefdr[j]+2.0*mw*d2cRefdr2[j][k])/(mw*mw*cRef*cRef*cRef)
                            + 2000.0*(cRef*dmwdr[j]+2.0*mw*dcRefdr[j])*dmwdr[k]/(mw*mw*mw*cRef*cRef*cRef)
                            + 3000.0*(cRef*dmwdr[j]+2.0*mw*dcRefdr[j])*dcRefdr[k]/(mw*mw*cRef*cRef*cRef*cRef)
                            + 2.0*tr*dalphadrk*dalphadrj/(cp)
                + 2.0*tr*alpha*d2alphadrjdrk/(cp)
                - 2.0*tr*alpha*dalphadrj*dcpdr[k]/(cp*cp)
                - tr*2.0*alpha*dalphadrk*dcpdr[j]/(cp*cp)
                + 2.0*tr*alpha*alpha*dcpdr[j]*dcpdr[k]/(cp*cp*cp)
                    )
                 - 2.0*dvdr[k]*dvdr[j]*(- 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
                                  + 2.0*tr*alpha*dalphadri/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp))
                - 2.0*v*dvdr[j]*(- 1000.0*(dcRefdr[k]*dmwdr[i]+2.0*dmwdr[k]*dcRefdr[i]+2.0*mw*d2cRefdr2[i][k])/(mw*mw*cRef*cRef*cRef)
                            + 2000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])*dmwdr[k]/(mw*mw*mw*cRef*cRef*cRef)
                            + 3000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])*dcRefdr[k]/(mw*mw*cRef*cRef*cRef*cRef)
                            + 2.0*tr*dalphadrk*dalphadri/(cp)
                + 2.0*tr*alpha*d2alphadridrk/(cp)
                - 2.0*tr*alpha*dalphadri*dcpdr[k]/(cp*cp)
                - tr*2.0*alpha*dalphadrk*dcpdr[i]/(cp*cp)
                + 2.0*tr*alpha*alpha*dcpdr[i]*dcpdr[k]/(cp*cp*cp)
                    )
                - 2.0*v*dvdr[k]*(- 1000.0*(dcRefdr[j]*dmwdr[i]+2.0*dmwdr[j]*dcRefdr[i]+2.0*mw*d2cRefdr2[i][j])/(mw*mw*cRef*cRef*cRef)
                                            + 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])*(2.0*dmwdr[j]*cRef+3.0*mw*dcRefdr[j])/(mw*mw*mw*cRef*cRef*cRef*cRef)
                                            + 2.0*tr*alpha*d2alphadridrj/(cp)
                    + 2.0*tr*dalphadrj*dalphadri/(cp)
                    - 2.0*tr*alpha*dalphadri*dcpdr[j]/(cp*cp)
                    - 2.0*tr*alpha*dalphadrj*dcpdr[i]/(cp*cp)
                    + 2.0*tr*alpha*alpha*dcpdr[i]*dcpdr[j]/(cp*cp*cp)
                            )
            - v*v*(- 1000.0*(d2cRefdr2[j][k]*dmwdr[i]+2.0*dmwdr[j]*d2cRefdr2[i][k]+2.0*dmwdr[k]*d2cRefdr2[i][j]+2.0*mw*d3cRefdr3[i][j][k])/(mw*mw*cRef*cRef*cRef)
                                + 2000.0*(dcRefdr[j]*dmwdr[i]+2.0*dmwdr[j]*dcRefdr[i]+2.0*mw*d2cRefdr2[i][j])*dmwdr[k]/(mw*mw*mw*cRef*cRef*cRef)
                                + 3000.0*(dcRefdr[j]*dmwdr[i]+2.0*dmwdr[j]*dcRefdr[i]+2.0*mw*d2cRefdr2[i][j])*dcRefdr[k]/(mw*mw*cRef*cRef*cRef*cRef)
                                + 1000.0*(dcRefdr[k]*dmwdr[i]+2.0*dmwdr[k]*dcRefdr[i]+2.0*mw*d2cRefdr2[i][k])*(2.0*dmwdr[j]*cRef+3.0*mw*dcRefdr[j])/(mw*mw*mw*cRef*cRef*cRef*cRef)
             + 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])*(2.0*dmwdr[j]*dcRefdr[k]+3.0*dmwdr[k]*dcRefdr[j]+3.0*mw*d2cRefdr2[j][k])/(mw*mw*mw*cRef*cRef*cRef*cRef)
             - 3000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])*(2.0*dmwdr[j]*cRef+3.0*mw*dcRefdr[j])*dmwdr[k]/(mw*mw*mw*mw*cRef*cRef*cRef*cRef)
             - 4000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])*(2.0*dmwdr[j]*cRef+3.0*mw*dcRefdr[j])*dcRefdr[k]/(mw*mw*mw*cRef*cRef*cRef*cRef*cRef)
                                + 2.0*tr*dalphadrk*d2alphadridrj/(cp) + 2.0*tr*alpha*d3alphadr3/(cp) - 2.0*tr*alpha*d2alphadridrj*dcpdr[k]/(cp*cp)
             + 2.0*tr*d2alphadrjdrk*dalphadri/(cp) + 2.0*tr*dalphadrj*d2alphadridrk/(cp) - 2.0*tr*dalphadrj*dalphadri*dcpdr[k]/(cp*cp)
             - 2.0*tr*dalphadrk*dalphadri*dcpdr[j]/(cp*cp) - 2.0*tr*alpha*d2alphadridrk*dcpdr[j]/(cp*cp) + 4.0*tr*alpha*dalphadri*dcpdr[j]*dcpdr[k]/(cp*cp*cp)
             - 2.0*tr*dalphadrk*dalphadrj*dcpdr[i]/(cp*cp) - 2.0*tr*alpha*d2alphadrjdrk*dcpdr[i]/(cp*cp) + 4.0*tr*alpha*dalphadrj*dcpdr[i]*dcpdr[k]/(cp*cp*cp)
             + 4.0*tr*alpha*dalphadrk*dcpdr[i]*dcpdr[j]/(cp*cp*cp) - 6.0*tr*alpha*alpha*dcpdr[i]*dcpdr[j]*dcpdr[k]/(cp*cp*cp*cp)
                );
    double d3cdr3     = d3cRefdr3[i][j][k] + (t-tr)*d4cdr3dt[i][j][k];
                    double d3v1dr3    = - 2.0*(d2v0drjdrk*dv0dri + dv0drj*d2v0dridrk)*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                        - 2.0*dv0drj*dv0dri*(- 1000.0*dmwdr[k]/(mw*mw*c*c) - 2000.0*dcdrk/(mw*c*c*c) + t*2.0*alpha*dalphadrk/(cp) - t*alpha*alpha*dcpdr[k]/(cp*cp))
                        - 2.0*(dv0drk*d2v0dridrj + v0*d3v0dr3)*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
            - 2.0*v0*d2v0dridrj*(- 1000.0*dmwdr[k]/(mw*mw*c*c) - 2000.0*dcdrk/(mw*c*c*c) + 2.0*t*alpha*dalphadrk/(cp) - t*alpha*alpha*dcpdr[k]/(cp*cp))
                        - 2.0*(dv0drk*dv0dri + v0*d2v0dridrk)*(- 1000.0*(c*dmwdr[j]+2.0*mw*dcdrj)/(mw*mw*c*c*c)
                                                                                                + 2.0*t*alpha*dalphadrj/(cp)
                                                                    - t*alpha*alpha*dcpdr[j]/(cp*cp)
                                      )
            - 2.0*v0*dv0dri*(- 1000.0*(dcdrk*dmwdr[j]+2.0*dmwdr[k]*dcdrj+2.0*mw*d2cdrjdrk)/(mw*mw*c*c*c)
                                                    + 2000.0*(c*dmwdr[j]+2.0*mw*dcdrj)*dmwdr[k]/(mw*mw*mw*c*c*c)
                        + 3000.0*(c*dmwdr[j]+2.0*mw*dcdrj)*dcdrk/(mw*mw*c*c*c*c)
                                                    + 2.0*t*(dalphadrk*dalphadrj + alpha*d2alphadrjdrk)/(cp) - 2.0*t*alpha*dalphadrj*dcpdr[k]/(cp*cp)
                        - 2.0*t*alpha*dalphadrk*dcpdr[j]/(cp*cp) + 2.0*t*alpha*alpha*dcpdr[j]*dcpdr[k]/(cp*cp*cp)
                )
                                - 2.0*(dv0drk*dv0drj + v0*d2v0drjdrk)*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)/(mw*mw*c*c*c)
                                                                                                + 2.0*t*alpha*dalphadri/(cp)
                                                                    - t*alpha*alpha*dcpdr[i]/(cp*cp)
                                      )
            - 2.0*v0*dv0drj*(- 1000.0*(dcdrk*dmwdr[i]+2.0*dmwdr[k]*dcdri+2.0*mw*d2cdridrk)/(mw*mw*c*c*c)
                                                    + 2000.0*(c*dmwdr[i]+2.0*mw*dcdri)*dmwdr[k]/(mw*mw*mw*c*c*c)
                        + 3000.0*(c*dmwdr[i]+2.0*mw*dcdri)*dcdrk/(mw*mw*c*c*c*c)
                                                    + 2.0*t*(dalphadrk*dalphadri + alpha*d2alphadridrk)/(cp) - 2.0*t*alpha*dalphadri*dcpdr[k]/(cp*cp)
                        - 2.0*t*alpha*dalphadrk*dcpdr[i]/(cp*cp) + 2.0*t*alpha*alpha*dcpdr[i]*dcpdr[k]/(cp*cp*cp)
                )
            - 2.0*v0*dv0drk*(- 1000.0*(dcdrj*dmwdr[i]+2.0*dmwdr[j]*dcdri+2.0*mw*d2cdridrj)/(mw*mw*c*c*c)
                       + 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)*(2.0*dmwdr[j]*c+3.0*mw*dcdrj)/(mw*mw*mw*c*c*c*c)
                       + 2.0*t*alpha*d2alphadridrj/(cp)
                       + 2.0*t*dalphadrj*dalphadri/(cp) - 2.0*t*alpha*dalphadri*dcpdr[j]/(cp*cp)
                       - 2.0*t*alpha*dalphadrj*dcpdr[i]/(cp*cp) + 2.0*t*alpha*alpha*dcpdr[i]*dcpdr[j]/(cp*cp*cp)
                            )
            - v0*v0*(- 1000.0*(d2cdrjdrk*dmwdr[i]+2.0*dmwdr[j]*d2cdridrk+2.0*dmwdr[k]*d2cdridrj+2.0*mw*d3cdr3)/(mw*mw*c*c*c)
                                    + 2000.0*(dcdrj*dmwdr[i]+2.0*dmwdr[j]*dcdri+2.0*mw*d2cdridrj)*dmwdr[k]/(mw*mw*mw*c*c*c)
               + 3000.0*(dcdrj*dmwdr[i]+2.0*dmwdr[j]*dcdri+2.0*mw*d2cdridrj)*dcdrk/(mw*mw*c*c*c*c)
                                    + 1000.0*(dcdrk*dmwdr[i]+2.0*dmwdr[k]*dcdri+2.0*mw*d2cdridrk)*(2.0*dmwdr[j]*c+3.0*mw*dcdrj)/(mw*mw*mw*c*c*c*c)
               + 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)*(2.0*dmwdr[j]*dcdrk+3.0*dmwdr[k]*dcdrj+3.0*mw*d2cdrjdrk)/(mw*mw*mw*c*c*c*c)
               - 3000.0*(c*dmwdr[i]+2.0*mw*dcdri)*(2.0*dmwdr[j]*c+3.0*mw*dcdrj)*dmwdr[k]/(mw*mw*mw*mw*c*c*c*c)
               - 4000.0*(c*dmwdr[i]+2.0*mw*dcdri)*(2.0*dmwdr[j]*c+3.0*mw*dcdrj)*dcdrk/(mw*mw*mw*c*c*c*c*c)
                                    + 2.0*t*(dalphadrk*d2alphadridrj + alpha*d3alphadr3)/(cp) - 2.0*t*alpha*d2alphadridrj*dcpdr[k]/(cp*cp)
                                    + 2.0*t*(d2alphadrjdrk*dalphadri + dalphadrj*d2alphadridrk)/(cp) - 2.0*t*dalphadrj*dalphadri*dcpdr[k]/(cp*cp)
               - 2.0*t*(dalphadrk*dalphadri + alpha*d2alphadridrk)*dcpdr[j]/(cp*cp) + 4.0*t*alpha*dalphadri*dcpdr[j]*dcpdr[k]/(cp*cp*cp)
                                    - 2.0*t*(dalphadrk*dalphadrj + alpha*d2alphadrjdrk)*dcpdr[i]/(cp*cp) + 4.0*t*alpha*dalphadrj*dcpdr[i]*dcpdr[k]/(cp*cp*cp)
               + 4.0*t*alpha*dalphadrk*dcpdr[i]*dcpdr[j]/(cp*cp*cp) - 6.0*t*alpha*alpha*dcpdr[i]*dcpdr[j]*dcpdr[k]/(cp*cp*cp*cp)
                    );
    double d3denomdr3 = 2.0*d3v1Refdr3*d3vdp3 + 2.0*d2v1Refdridrj*d4vdrdp3[k] + 2.0*d2v1Refdridrk*d4vdrdp3[j] + 2.0*d2v1Refdrjdrk*d4vdrdp3[i];
    double d3adr3 = (-d3v1Refdr3*d4vdp4/2.0-d2v1Refdridrj*d5vdrdp4[k]/2.0-d2v1Refdridrk*d5vdrdp4[j]/2.0-d2v1Refdrjdrk*d5vdrdp4[i]/2.0)/denom
      - (d3vdrdp2[i]*d4vdrdp3[j]+d3vdrdp2[j]*d4vdrdp3[i]-d2v1Refdridrj*d4vdp4/2.0-dv1Refdri*d5vdrdp4[j]/2.0-dv1Refdrj*d5vdrdp4[i]/2.0)*ddenomdrk/(denom*denom)
                  - (d3vdrdp2[i]*d4vdrdp3[k]+d3vdrdp2[k]*d4vdrdp3[i]-d2v1Refdridrk*d4vdp4/2.0-dv1Refdri*d5vdrdp4[k]/2.0-dv1Refdrk*d5vdrdp4[i]/2.0)*ddenomdrj/(denom*denom)
      - (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdri*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)*d2denomdrjdrk/(denom*denom)
      + 2.0*(d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdri*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)*ddenomdrj*ddenomdrk/(denom*denom*denom)
                    - (d3vdrdp2[j]*d4vdrdp3[k]+d3vdrdp2[k]*d4vdrdp3[j]-d2v1Refdrjdrk*d4vdp4/2.0-dv1Refdrj*d5vdrdp4[k]/2.0-dv1Refdrk*d5vdrdp4[j]/2.0)*ddenomdri/(denom*denom)
      - (d3vdrdp2[j]*d3vdp3+d2vdp2*d4vdrdp3[j]-dv1Refdrj*d4vdp4/2.0-v1Ref*d5vdrdp4[j]/2.0)*d2denomdridrk/(denom*denom)
      + 2.0*(d3vdrdp2[j]*d3vdp3+d2vdp2*d4vdrdp3[j]-dv1Refdrj*d4vdp4/2.0-v1Ref*d5vdrdp4[j]/2.0)*ddenomdri*ddenomdrk/(denom*denom*denom)
                                                + 2.0*(d3vdrdp2[k]*d3vdp3+d2vdp2*d4vdrdp3[k]-dv1Refdrk*d4vdp4/2.0-v1Ref*d5vdrdp4[k]/2.0)*ddenomdri*ddenomdrj/(denom*denom*denom)
      + 2.0*(d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*(d2denomdridrk*ddenomdrj+ddenomdri*d2denomdrjdrk)/(denom*denom*denom)
      - 6.0*(d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*ddenomdri*ddenomdrj*ddenomdrk/(denom*denom*denom*denom)
                                                - (d3vdrdp2[k]*d3vdp3+d2vdp2*d4vdrdp3[k]-dv1Refdrk*d4vdp4/2.0-v1Ref*d5vdrdp4[k]/2.0)*d2denomdridrj/(denom*denom)
      - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*d3denomdr3/(denom*denom)
      + 2.0*(d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*d2denomdridrj*ddenomdrk/(denom*denom*denom);
                    double d3bdr3 = -(d3vdrdp2[i]*d5vdrdp4[j]/4.0+d3vdrdp2[j]*d5vdrdp4[i]/4.0-2.0/3.0*d4vdrdp3[j]*d4vdrdp3[i])/(denom*denom)*ddenomdrk
                        -(d3vdrdp2[i]*d5vdrdp4[k]/4.0+d3vdrdp2[k]*d5vdrdp4[i]/4.0-2.0/3.0*d4vdrdp3[k]*d4vdrdp3[i])/(denom*denom)*ddenomdrj
                    + 2.0*(d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])/(denom*denom*denom)*ddenomdrj*ddenomdrk
                    - (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])/(denom*denom)*d2denomdrjdrk
                    - (d3vdrdp2[j]*d5vdrdp4[k]/4.0+d3vdrdp2[k]*d5vdrdp4[j]/4.0-2.0/3.0*d4vdrdp3[k]*d4vdrdp3[j])/(denom*denom)*ddenomdri
                    + 2.0*(d3vdrdp2[j]*d4vdp4/4.0+d2vdp2*d5vdrdp4[j]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[j])/(denom*denom*denom)*ddenomdri*ddenomdrk
                    - (d3vdrdp2[j]*d4vdp4/4.0+d2vdp2*d5vdrdp4[j]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[j])/(denom*denom)*d2denomdridrk
                    + 2.0*(d3vdrdp2[k]*d4vdp4/4.0+d2vdp2*d5vdrdp4[k]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[k])/(denom*denom*denom)*ddenomdri*ddenomdrj
                    - 6.0*(d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom*denom*denom)*ddenomdri*ddenomdrj*ddenomdrk
                    + 2.0*(d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom*denom)*d2denomdridrk*ddenomdrj
                    + 2.0*(d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom*denom)*ddenomdri*d2denomdrjdrk
                    - (d3vdrdp2[k]*d4vdp4/4.0+d2vdp2*d5vdrdp4[k]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[k])/(denom*denom)*d2denomdridrj
                    + 2.0*(d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom*denom)*d2denomdridrj*ddenomdrk;

                    double d3gIntdr3 = 0.0;

                    if ((a == 0.0) && (b == 0.0)) {
                        d3gIntdr3  = d3v0dr3*(p-pr) + d3v1dr3*(p-pr)*(p-pr)/2.0;

                    } else if ((a != 0.0) && (b == 0.0)) {
                        printf("*-->Exception in fillD3GDR3 (liquid.c). a is not equal to zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
                        d3gIntdr3  = 0.0;

                    } else if ((a == 0.0) && (b != 0.0)) {
                        printf("*-->Exception in fillD3GDR3 (liquid.c). a is zero, b is not equal to zero.\n"); liqERRstate = ERR_A_ZERO;
                        d3gIntdr3  = 0.0;

                    } else if (sum > 0.0) {
                        d3gIntdr3  = d3gdr3GMAP(p/10000.0, pr/10000.0,
                              v0,                 v1*10000.0,      d2vdp2*10000.0*10000.0,         a*10000.0,         b*10000.0*10000.0,
                              dv0dri,         dv1dri*10000.0, d3vdrdp2[i]*10000.0*10000.0,     dadri*10000.0,     dbdri*10000.0*10000.0,
                              dv0drj,         dv1drj*10000.0, d3vdrdp2[j]*10000.0*10000.0,     dadrj*10000.0,     dbdrj*10000.0*10000.0,
                dv0drk,         dv1drk*10000.0, d3vdrdp2[k]*10000.0*10000.0,     dadrk*10000.0,     dbdrk*10000.0*10000.0,
                                                                        d2v0dridrj, d2v1dridrj*10000.0,                  0.0, d2adridrj*10000.0, d2bdridrj*10000.0*10000.0,
                d2v0dridrk, d2v1dridrk*10000.0,                  0.0, d2adridrk*10000.0, d2bdridrk*10000.0*10000.0,
                d2v0drjdrk, d2v1drjdrk*10000.0,                  0.0, d2adrjdrk*10000.0, d2bdrjdrk*10000.0*10000.0,
                                                                        d3v0dr3,       d3v1dr3*10000.0,                  0.0,    d3adr3*10000.0,    d3bdr3*10000.0*10000.0)*10000.0;

                    } else if (sum == 0.0) {
                        printf("*-->Exception in fillD3GDR3 (liquid.c). a*a-4*b is equal to zero.\n"); liqERRstate = ERR_SUM_ZERO;
                        d3gIntdr3  = 0.0;

                    } else if(sum < 0.0) {
                        d3gIntdr3  = d3gdr3LMAP(p/10000.0, pr/10000.0,
                              v0,                 v1*10000.0,      d2vdp2*10000.0*10000.0,         a*10000.0,         b*10000.0*10000.0,
                              dv0dri,         dv1dri*10000.0, d3vdrdp2[i]*10000.0*10000.0,     dadri*10000.0,     dbdri*10000.0*10000.0,
                              dv0drj,         dv1drj*10000.0, d3vdrdp2[j]*10000.0*10000.0,     dadrj*10000.0,     dbdrj*10000.0*10000.0,
                dv0drk,         dv1drk*10000.0, d3vdrdp2[k]*10000.0*10000.0,     dadrk*10000.0,     dbdrk*10000.0*10000.0,
                                                                        d2v0dridrj, d2v1dridrj*10000.0,                  0.0, d2adridrj*10000.0, d2bdridrj*10000.0*10000.0,
                d2v0dridrk, d2v1dridrk*10000.0,                  0.0, d2adridrk*10000.0, d2bdridrk*10000.0*10000.0,
                d2v0drjdrk, d2v1drjdrk*10000.0,                  0.0, d2adrjdrk*10000.0, d2bdrjdrk*10000.0*10000.0,
                                                                        d3v0dr3,       d3v1dr3*10000.0,                  0.0,    d3adr3*10000.0,    d3bdr3*10000.0*10000.0)*10000.0;

                    }

                    result[i][j][k] += coeffCN*d3gIntdr3;

                    if      (i == j && j != k) { result[i][k][j] += coeffCN*d3gIntdr3; result[k][i][j] += coeffCN*d3gIntdr3; }
    else if (i != j && j == k) { result[j][i][k] += coeffCN*d3gIntdr3; result[j][k][i] += coeffCN*d3gIntdr3; }
    else if (i != j && j != k) {
      result[i][k][j] += coeffCN*d3gIntdr3;
      result[j][i][k] += coeffCN*d3gIntdr3;
      result[j][k][i] += coeffCN*d3gIntdr3;
      result[k][i][j] += coeffCN*d3gIntdr3;
      result[k][j][i] += coeffCN*d3gIntdr3;
    }
  }
            }
        }
    }

#endif /* USE_GHIORSO_KRESS_MODEL */

}

static void fillD3GDR2DS (double r[NR], double s[NT], double t, double p, double result[NR][NR][NT]) {
    int i, j, k, l;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) for (j=0; j<NR; j++) memset(result[i][j], '\0', (size_t) NT*sizeof(double));

    /* Configurational entropy terms */
    for (j=0; j<NR; j++) {
        for (k=j; k<NR; k++) {
            for (l=0; l<NS; l++) {
                double config = 0.0;
                for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
                    config += nSpecies*(
                        (d2xSpeciesdrds[i][j][l]*dxSpeciesdr[i][k]+dxSpeciesdr[i][j]*d2xSpeciesdrds[i][k][l])/xSpecies[i]
                        - dxSpeciesdr[i][j]*dxSpeciesdr[i][k]*dxSpeciesds[i][l]/(xSpecies[i]*xSpecies[i]) )
                        + dnSpeciesds[l]*dxSpeciesdr[i][j]*dxSpeciesdr[i][k]/xSpecies[i];
                if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
                    config += nSpecies*( (d2xSpeciesdrds[nH2O][j][l]*dxSpeciesdr[nH2O][k] + dxSpeciesdr[nH2O][j]*d2xSpeciesdrds[nH2O][k][l])
                            *(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]))
                        - dxSpeciesdr[nH2O][j]*dxSpeciesdr[nH2O][k]*dxSpeciesds[nH2O][l]
                            *(1.0/(xSpecies[nH2O]*xSpecies[nH2O]) - 1.0/((1.0-xSpecies[nH2O])*(1.0-xSpecies[nH2O]))) )
                        + dnSpeciesds[l]*dxSpeciesdr[nH2O][j]*dxSpeciesdr[nH2O][k]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]));

                result[j][k][l] += R*t*config;
                result[k][j][l]  = result[j][k][l];
            }
        }
    }

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double d2gIntIVdr2[NR][NR];

        d2IntegralV_GKdr2(r, s, t, p, d2gIntIVdr2);
        for (i=0; i<NR; i++) {
            for (j=i; j<NR; j++) {
                for (k=0; k<NY; k++) {
    result[i][j][NS+k] += (fCN[k]-1.0)*d2gIntIVdr2[i][j];
    result[j][i][NS+k]  = result[i][j][NS+k];
  }
            }
        }
    }
#endif /* USE_GHIORSO_KRESS_MODEL */

}

static void fillD3GDR2DT (double r[NR], double s[NT], double t, double p, double result[NR][NR]) {
    int i, j, k;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) for (j=i; j<NR; j++) {
        result[i][j] = (i == j) ? -2.0*srr[i][i] : -srr[i][j];
        result[j][i] = result[i][j];
    }

    /* Configurational entropy terms */
    for (j=0; j<NR; j++) {
        for (k=j; k<NR; k++) {
            double config = 0.0;
            for (i=0; i<NE; i++) if (xSpecies[i] > 0.0 && dxSpeciesdr[i][j] != 0.0 && dxSpeciesdr[i][k] != 0.0)
                config += dxSpeciesdr[i][j]*dxSpeciesdr[i][k]/xSpecies[i];
            if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0 && dxSpeciesdr[nH2O][j] != 0.0 && dxSpeciesdr[nH2O][k] != 0.0)
                config += dxSpeciesdr[nH2O][j]*dxSpeciesdr[nH2O][k]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]));
            result[j][k] += R*nSpecies*config;
            result[k][j]  = result[j][k];
        }
    }

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        const double pr         = 1.0;
        const double tr         = 1673.15;
        double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum, coeffCN;
        double dvdr[NR], d2vdrdt[NR], dcRefdr[NR], d2cdrdt[NR], dmwdr[NR], dcpdr[NR], d3vdrdp2[NR], d4vdrdp3[NR], d5vdrdp4[NR], dmOxTotdr[NR], denom;
        double d2cRefdr2[NR][NR], d3cdr2dt[NR][NR];
        double dv0dt, dv1dt;

        if (fabs(p-pr) < 10.0*DBL_EPSILON) return;

        for (i=0, coeffCN=1.0; i<NY; i++) coeffCN += (fCN[i]-1.0)*s[NS+i];

        /* Convert input composition (r) to liquid moles (m)  */
        for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

        /* Compute moles and total moles of oxides */
        for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
        if (mOxTot == 0.0) return;

        /* Deal with the special case of FeO1.3 */
        mOx[NA] = 0.0;
        if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
            const double y = 0.3;
            mOx[iOxFeO1_3] = 0.0;
            if (iCmpFe2SiO4_6 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
            }
            if (iCmpFe2AlO4_1 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
            }
        }

        for (i=0; i<NR; i++) {
            dvdr[i]    = 0.0; d2vdrdt[i]  = 0.0; dcRefdr[i]  = 0.0; d2cdrdt[i]  = 0.0; dmwdr[i] = 0.0;
            dcpdr[i]   = 0.0; d3vdrdp2[i] = 0.0; d4vdrdp3[i] = 0.0; d5vdrdp4[i] = 0.0;
            for (j=0, dmOxTotdr[i]=0.0; j<nc; j++) dmOxTotdr[i] += (liquid[i+1].liqToOx)[j] - (liquid[0].liqToOx)[j];
            for (j=i; j<NR; j++) { d2cRefdr2[i][j] = 0.0; d3cdr2dt[i][j] = 0.0; }
        }

        for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
            v       += mOx[i]*bulkSystem[i].gk_v;
            dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
            cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
            dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
            cp      += mOx[i]*bulkSystem[i].gk_cp;
            d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
            mw      += mOx[i]*bulkSystem[i].mw;

            for (j=0; j<NR; j++) {
                dvdr[j]      += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_v;
                d2vdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dvdt;
                dmwdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].mw;
                dcpdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_cp;
  dcRefdr[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c/mOxTot
                - mOx[i]*bulkSystem[i].gk_c*dmOxTotdr[j]/(mOxTot*mOxTot);
  dcRefdr[j]   += (iOxAl2O3 != -1) ? ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
                - 2.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot) : 0.0;
  d2cdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt/mOxTot
                - mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotdr[j]/(mOxTot*mOxTot);
                d3vdrdp2[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
                d4vdrdp3[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
                d5vdrdp4[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
            }
            if (iCmpAl2O3 != -1) dcRefdr[iCmpAl2O3] += ((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*mOx[i]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot);

            for (j=0; j<NR; j++) {
  for (k=j; k<NR; k++) {
    d2cRefdr2[j][k] += -((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c*dmOxTotdr[k]/(mOxTot*mOxTot)
                - ((liquid[k+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c*dmOxTotdr[j]/(mOxTot*mOxTot)
            + 2.0*mOx[i]*bulkSystem[i].gk_c*dmOxTotdr[j]*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot);
    d2cRefdr2[j][k] += (iOxAl2O3 != -1) ?
                - 2.0*((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot)
                - 2.0*((liquid[k+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot)
            + 6.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot*mOxTot)
            : 0.0;

    d3cdr2dt[j][k]  += -((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt*dmOxTotdr[k]/(mOxTot*mOxTot)
                - ((liquid[k+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt*dmOxTotdr[j]/(mOxTot*mOxTot)
            + 2.0*mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotdr[j]*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot);
  }
  if (iCmpAl2O3 != -1) d2cRefdr2[(j<iCmpAl2O3) ? j : iCmpAl2O3][(j>iCmpAl2O3) ? j : iCmpAl2O3] +=
    ((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
    - 2.0*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*mOx[i]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot);
            }
            if (iCmpAl2O3 != -1) d2cRefdr2[iCmpAl2O3][iCmpAl2O3] +=
  ((liquid[iCmpAl2O3+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
  - 2.0*mOx[i]*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*bulkSystem[i].gk_cXal2o3*dmOxTotdr[iCmpAl2O3]/(mOxTot*mOxTot*mOxTot);

        }
        if (v == 0.0) return;

        alpha   = dvdt/v;
        v0      = v*exp(alpha*(t-tr));
        v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
        c       = cRef + (t-tr)*dcdt;
        v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
        v2      = d2vdp2;
        denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
        a        = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
        b        = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;
        sum     = a*a - 4.0*b;
        dv0dt   = alpha*v0;
        dv1dt   = - 2.0*v0*v0*alpha*(1000.0/(mw*c*c) + t*alpha*alpha/(cp)) - v0*v0*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp));

        for (i=0; i<NR; i++) {
            double dalphadri   = d2vdrdt[i]/v - dvdt*dvdr[i]/(v*v);
            double dv0dri      = (dvdr[i] + v*dalphadri*(t-tr))*exp(alpha*(t-tr));
            double dv1Refdri   = -2.0*v*dvdr[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                           -v*v*(- 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
                        + 2.0*tr*alpha*dalphadri/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp));
            double dcdri       = dcRefdr[i] + (t-tr)*d2cdrdt[i];
            double dv1dri      = -2.0*v0*dv0dri*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                           -v0*v0*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)/(mw*mw*c*c*c)
                            + 2.0*t*alpha*dalphadri/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp));
            double ddenomdri   = 2.0*dv1Refdri*d3vdp3+2.0*v1Ref*d4vdrdp3[i]-6.0*d2vdp2*d3vdrdp2[i];
            double dadri       = (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdri*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)/denom
                         - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdri;
            double dbdri       = (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])/denom
                         - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdri;
            double d2v0dridt   = dvdr[i]*alpha*exp(alpha*(t-tr))+v*(d2vdrdt[i]/v-dvdt/(v*v)*dvdr[i])*exp(alpha*(t-tr))
                         + v*(d2vdrdt[i]/v-dvdt/(v*v)*dvdr[i])*(t-tr)*alpha*exp(alpha*(t-tr));
            double d2v1dridt   = -2.0*(dv0dt*dv0dri + v0*d2v0dridt)*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                         - 2.0*v0*dv0dri*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp))
                         - 2.0*v0*dv0dt*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)/(mw*mw*c*c*c)
                                        + 2.0*t*alpha*dalphadri/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp))
             - v0*v0*(- 1000.0*(dcdt*dmwdr[i]+2.0*mw*d2cdrdt[i])/(mw*mw*c*c*c) + 3000.0*(c*dmwdr[i]+2.0*mw*dcdri)*dcdt/(mw*mw*c*c*c*c)
                                        + 2.0*alpha*dalphadri/(cp) - alpha*alpha*dcpdr[i]/(cp*cp));

            for (j=i; j<NR; j++) {
                double dalphadrj   = d2vdrdt[j]/v - dvdt*dvdr[j]/(v*v);
                double dv0drj       = (dvdr[j] + v*dalphadrj*(t-tr))*exp(alpha*(t-tr));
                double dv1Refdrj   = -2.0*v*dvdr[j]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                 -v*v*(- 1000.0*(cRef*dmwdr[j]+2.0*mw*dcRefdr[j])/(mw*mw*cRef*cRef*cRef)
                            + 2.0*tr*alpha*dalphadrj/(cp) - tr*alpha*alpha*dcpdr[j]/(cp*cp));
                double dcdrj       = dcRefdr[j] + (t-tr)*d2cdrdt[j];
                double dv1drj       = -2.0*v0*dv0drj*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                 -v0*v0*(- 1000.0*(c*dmwdr[j]+2.0*mw*dcdrj)/(mw*mw*c*c*c)
                                + 2.0*t*alpha*dalphadrj/(cp) - t*alpha*alpha*dcpdr[j]/(cp*cp));
                double ddenomdrj   = 2.0*dv1Refdrj*d3vdp3+2.0*v1Ref*d4vdrdp3[j]-6.0*d2vdp2*d3vdrdp2[j];
                double dadrj       = (d3vdrdp2[j]*d3vdp3+d2vdp2*d4vdrdp3[j]-dv1Refdrj*d4vdp4/2.0-v1Ref*d5vdrdp4[j]/2.0)/denom
               - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdrj;
                double dbdrj       = (d3vdrdp2[j]*d4vdp4/4.0+d2vdp2*d5vdrdp4[j]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[j])/denom
                         - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdrj;
                double d2v0drjdt   = dvdr[j]*alpha*exp(alpha*(t-tr))+v*(d2vdrdt[j]/v-dvdt/(v*v)*dvdr[j])*exp(alpha*(t-tr))
               + v*(d2vdrdt[j]/v-dvdt/(v*v)*dvdr[j])*(t-tr)*alpha*exp(alpha*(t-tr));
                double d2v1drjdt   = -2.0*(dv0dt*dv0drj + v0*d2v0drjdt)*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                           - 2.0*v0*dv0drj*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp))
                           - 2.0*v0*dv0dt*(- 1000.0*(c*dmwdr[j]+2.0*mw*dcdrj)/(mw*mw*c*c*c)
                                            + 2.0*t*alpha*dalphadrj/(cp) - t*alpha*alpha*dcpdr[j]/(cp*cp))
               - v0*v0*(- 1000.0*(dcdt*dmwdr[j]+2.0*mw*d2cdrdt[j])/(mw*mw*c*c*c) + 3000.0*(c*dmwdr[j]+2.0*mw*dcdrj)*dcdt/(mw*mw*c*c*c*c)
                                            + 2.0*alpha*dalphadrj/(cp) - alpha*alpha*dcpdr[j]/(cp*cp));

                double d2alphadr2  = -d2vdrdt[i]*dvdr[j]/(v*v)-d2vdrdt[j]*dvdr[i]/(v*v)+2.0*dvdt*dvdr[i]*dvdr[j]/(v*v*v);
                double d2v0dr2     = dvdr[i]*dalphadrj*(t-tr)*exp(alpha*(t-tr))
            + dvdr[j]*dalphadri*(t-tr)*exp(alpha*(t-tr))
            + v*d2alphadr2*(t-tr)*exp(alpha*(t-tr))
            + v*dalphadri*pow(t-tr,2.0)*dalphadrj*exp(alpha*(t-tr));
                double d2v1Refdr2  = -2.0*dvdr[j]*dvdr[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                      - 2.0*v*dvdr[i]*(- 1000.0*(cRef*dmwdr[j]+2.0*mw*dcRefdr[j])/(mw*mw*cRef*cRef*cRef)
                                                + 2.0*tr*alpha*dalphadrj/(cp) - tr*alpha*alpha*dcpdr[j]/(cp*cp))
                            - 2.0*v*dvdr[j]*(- 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
                                                + 2.0*tr*alpha*dalphadri/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp))
                -v*v*(- 1000.0*(dcRefdr[j]*dmwdr[i]+2.0*dmwdr[j]*dcRefdr[i]+2.0*mw*d2cRefdr2[i][j])/(mw*mw*cRef*cRef*cRef)
                            + 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])*(2.0*dmwdr[j]*cRef+3.0*mw*dcRefdr[j])/(mw*mw*mw*cRef*cRef*cRef*cRef)
                            + 2.0*tr*alpha*d2alphadr2/(cp)
           + 2.0*tr*dalphadrj*dalphadri/(cp) - 2.0*tr*alpha*dalphadri*dcpdr[j]/(cp*cp)
           - 2.0*tr*alpha*dalphadrj*dcpdr[i]/(cp*cp) + 2.0*tr*alpha*alpha*dcpdr[i]*dcpdr[j]/(cp*cp*cp));
                double d2cdr2       = d2cRefdr2[i][j] + (t-tr)*d3cdr2dt[i][j];
                double d2v1dr2       = -2.0*dv0drj*dv0dri*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                      - 2.0*v0*d2v0dr2*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                      - 2.0*v0*dv0dri*(- 1000.0*(c*dmwdr[j]+2.0*mw*dcdrj)/(mw*mw*c*c*c)
                                                + 2.0*t*alpha*dalphadrj/(cp) - t*alpha*alpha*dcpdr[j]/(cp*cp))
                            - 2.0*v0*dv0drj*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)/(mw*mw*c*c*c)
                                                + 2.0*t*alpha*dalphadri/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp))
                -v0*v0*(- 1000.0*(dcdrj*dmwdr[i]+2.0*dmwdr[j]*dcdri+2.0*mw*d2cdr2)/(mw*mw*c*c*c)
                                + 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)*(2.0*dmwdr[j]*c+3.0*mw*dcdrj)/(mw*mw*mw*c*c*c*c)
             + 2.0*t*alpha*d2alphadr2/(cp)
                                + 2.0*t*dalphadrj*dalphadri/(cp) - 2.0*t*alpha*dalphadri*dcpdr[j]/(cp*cp)
             - 2.0*t*alpha*dalphadrj*dcpdr[i]/(cp*cp) + 2.0*t*alpha*alpha*dcpdr[i]*dcpdr[j]/(cp*cp*cp));
                double d2denomdr2  = 2.0*d2v1Refdr2*d3vdp3 + 2.0*dv1Refdri*d4vdrdp3[j] + 2.0*dv1Refdrj*d4vdrdp3[i] - 6.0*d3vdrdp2[j]*d3vdrdp2[i];
                double d2adr2       = (d3vdrdp2[i]*d4vdrdp3[j]+d3vdrdp2[j]*d4vdrdp3[i]-d2v1Refdr2*d4vdp4/2.0-dv1Refdri*d5vdrdp4[j]/2.0-dv1Refdrj*d5vdrdp4[i]/2.0)/denom
                                        - (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdri*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)*ddenomdrj/(denom*denom)
               - (d3vdrdp2[j]*d3vdp3+d2vdp2*d4vdrdp3[j]-dv1Refdrj*d4vdp4/2.0-v1Ref*d5vdrdp4[j]/2.0)*ddenomdri/(denom*denom)
                           + 2.0*(d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*ddenomdri*ddenomdrj/(denom*denom*denom)
                           - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*d2denomdr2/(denom*denom);
                double d2bdr2       = (d3vdrdp2[i]*d5vdrdp4[j]/4.0+d3vdrdp2[j]*d5vdrdp4[i]/4.0-2.0/3.0*d4vdrdp3[j]*d4vdrdp3[i])/denom
                           - (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])*ddenomdrj/(denom*denom)
               - (d3vdrdp2[j]*d4vdp4/4.0+d2vdp2*d5vdrdp4[j]/4.0-2.0*d4vdrdp3[j]*d3vdp3/3.0)*ddenomdri/(denom*denom)
            + 2.0*(d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)*ddenomdri*ddenomdrj/(denom*denom*denom)
            - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)*d2denomdr2/(denom*denom);
                double d3v0dr2dt   = dvdr[i]*dalphadrj*exp(alpha*(t-tr))+dvdr[i]*dalphadrj*(t-tr)*alpha*exp(alpha*(t-tr))
                           + dvdr[j]*(d2vdrdt[i]/v-dvdt/(v*v)*dvdr[i])*exp(alpha*(t-tr))
            + dvdr[j]*(d2vdrdt[i]/v-dvdt/(v*v)*dvdr[i])*(t-tr)*alpha*exp(alpha*(t-tr))
                                        + v*d2alphadr2*exp(alpha*(t-tr))+v*d2alphadr2*(t-tr)*alpha*exp(alpha*(t-tr))
                                        + 2.0*v*(d2vdrdt[i]/v-dvdt/(v*v)*dvdr[i])*(t-tr)*dalphadrj*exp(alpha*(t-tr))
                                        + v*(d2vdrdt[i]/v-dvdt/(v*v)*dvdr[i])*pow(t-tr,2.0)*dalphadrj*alpha*exp(alpha*(t-tr));
                double d3v1dr2dt   = -2.0*(d2v0drjdt*dv0dri+dv0drj*d2v0dridt)*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                                        - 2.0*dv0drj*dv0dri*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp))
            - 2.0*(dv0dt*d2v0dr2+v0*d3v0dr2dt)*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
            - 2.0*v0*d2v0dr2*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp))
                                        - 2.0*(dv0dt*dv0dri+v0*d2v0dridt)*(- 1000.0*(c*dmwdr[j]+2.0*mw*dcdrj)/(mw*mw*c*c*c)
                        + 2.0*t*alpha*dalphadrj/(cp) - t*alpha*alpha*dcpdr[j]/(cp*cp))
            - 2.0*v0*dv0dri*(- 1000.0*(dcdt*dmwdr[j]+2.0*mw*d2cdrdt[j])/(mw*mw*c*c*c) + 3000.0*dcdt*(c*dmwdr[j]+2.0*mw*dcdrj)/(mw*mw*c*c*c*c)
                        + 2.0*alpha*dalphadrj/(cp) - alpha*alpha*dcpdr[j]/(cp*cp))
               - 2.0*(dv0dt*dv0drj+v0*d2v0drjdt)*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)/(mw*mw*c*c*c)
                          + 2.0*t*alpha*dalphadri/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp))
            - 2.0*v0*dv0drj*(- 1000.0*(dcdt*dmwdr[i]+2.0*mw*d2cdrdt[i])/(mw*mw*c*c*c) + 3000.0*dcdt*(c*dmwdr[i]+2.0*mw*dcdri)/(mw*mw*c*c*c*c)
                          + 2.0*alpha*dalphadri/(cp) - alpha*alpha*dcpdr[i]/(cp*cp))
            - 2.0*v0*dv0dt*(- 1000.0*(dcdrj*dmwdr[i]+2.0*dmwdr[j]*dcdri+2.0*mw*d2cdr2)/(mw*mw*c*c*c)
                            + 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)*(2.0*dmwdr[j]*c+3.0*mw*dcdrj)/(mw*mw*mw*c*c*c*c)
                + 2.0*t*alpha*d2alphadr2/(cp)
                            + 2.0*t*dalphadrj*dalphadri/(cp) - 2.0*t*alpha*dalphadri*dcpdr[j]/(cp*cp)
                                    - 2.0*t*alpha*dalphadrj*dcpdr[i]/(cp*cp) + 2.0*t*alpha*alpha*dcpdr[i]*dcpdr[j]/(cp*cp*cp)
                            )
            - v0*v0*(- 1000.0*(d2cdrdt[j]*dmwdr[i]+2.0*dmwdr[j]*d2cdrdt[i]+2.0*mw*d3cdr2dt[i][j])/(mw*mw*c*c*c)
                  + 3000.0*(dcdrj*dmwdr[i]+2.0*dmwdr[j]*dcdri+2.0*mw*d2cdr2)*dcdt/(mw*mw*c*c*c*c)
               + 1000.0*(dcdt*dmwdr[i]+2.0*mw*d2cdrdt[i])*(2.0*dmwdr[j]*c+3.0*mw*dcdrj)/(mw*mw*mw*c*c*c*c)
                + 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)*(2.0*dmwdr[j]*dcdt+3.0*mw*d2cdrdt[j])/(mw*mw*mw*c*c*c*c)
                - 4000.0*(c*dmwdr[i]+2.0*mw*dcdri)*(2.0*dmwdr[j]*c+3.0*mw*dcdrj)*dcdt/(mw*mw*mw*c*c*c*c*c)
                + 2.0*alpha*d2alphadr2/(cp)
               + 2.0*dalphadrj*dalphadri/(cp) - 2.0*alpha*dalphadri*dcpdr[j]/(cp*cp)
                - 2.0*alpha*dalphadrj*dcpdr[i]/(cp*cp) + 2.0*alpha*alpha*dcpdr[i]*dcpdr[j]/(cp*cp*cp)
           );
                double d3gIntdr2dt = 0.0;

                if ((a == 0.0) && (b == 0.0)) {
                    d3gIntdr2dt  = d3v0dr2dt*(p-pr) + d3v1dr2dt*(p-pr)*(p-pr)/2.0;

                } else if ((a != 0.0) && (b == 0.0)) {
                    printf("*-->Exception in fillD3GDR2DT (liquid.c). a is not equal to zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
                    d3gIntdr2dt  = 0.0;

                } else if ((a == 0.0) && (b != 0.0)) {
                    printf("*-->Exception in fillD3GDR2DT (liquid.c). a is zero, b is not equal to zero.\n"); liqERRstate = ERR_A_ZERO;
                    d3gIntdr2dt  = 0.0;

                } else if (sum > 0.0) {
                    d3gIntdr2dt  = d3gdr2dtGMAP(p/10000.0, pr/10000.0,
                                v0,               v1*10000.0,      d2vdp2*10000.0*10000.0,      a*10000.0,      b*10000.0*10000.0,
                    dv0dri,       dv1dri*10000.0, d3vdrdp2[i]*10000.0*10000.0,  dadri*10000.0,  dbdri*10000.0*10000.0,
                    dv0drj,       dv1drj*10000.0, d3vdrdp2[j]*10000.0*10000.0,  dadrj*10000.0,  dbdrj*10000.0*10000.0,
                    d2v0dr2,     d2v1dr2*10000.0,                         0.0, d2adr2*10000.0, d2bdr2*10000.0*10000.0,
                    dv0dt,         dv1dt*10000.0,
                    d2v0dridt, d2v1dridt*10000.0,
                    d2v0drjdt, d2v1drjdt*10000.0,
                    d3v0dr2dt, d3v1dr2dt*10000.0)*10000.0;

                } else if (sum == 0.0) {
                    printf("*-->Exception in fillD3GDR2DT (liquid.c). a*a-4*b is equal to zero.\n"); liqERRstate = ERR_SUM_ZERO;
                    d3gIntdr2dt  = 0.0;

                } else if(sum < 0.0) {
                    d3gIntdr2dt  = d3gdr2dtLMAP(p/10000.0, pr/10000.0,
                                v0,               v1*10000.0,      d2vdp2*10000.0*10000.0,      a*10000.0,      b*10000.0*10000.0,
                    dv0dri,       dv1dri*10000.0, d3vdrdp2[i]*10000.0*10000.0,  dadri*10000.0,  dbdri*10000.0*10000.0,
                    dv0drj,       dv1drj*10000.0, d3vdrdp2[j]*10000.0*10000.0,  dadrj*10000.0,  dbdrj*10000.0*10000.0,
                    d2v0dr2,     d2v1dr2*10000.0,                         0.0, d2adr2*10000.0, d2bdr2*10000.0*10000.0,
                    dv0dt,         dv1dt*10000.0,
                    d2v0dridt, d2v1dridt*10000.0,
                    d2v0drjdt, d2v1drjdt*10000.0,
                    d3v0dr2dt, d3v1dr2dt*10000.0)*10000.0;

                }

                result[i][j] += coeffCN*d3gIntdr2dt;
  if (i != j) result[j][i] += coeffCN*d3gIntdr2dt;
            }
        }
    }

#endif /* USE_GHIORSO_KRESS_MODEL */

}

static void fillD3GDR2DP (double r[NR], double s[NT], double t, double p, double result[NR][NR]) {
    int i, j;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) for (j=i; j<NR; j++) {
        result[i][j] = (i == j) ? 2.0*vrr[i][i] : vrr[i][j];
        result[j][i] = result[i][j];
    }

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        const double pr         = 1.0;
        const double tr         = 1673.15;
        double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum, coeffCN;
        double dvdr[NR], d2vdrdt[NR], dcRefdr[NR], d2cdrdt[NR], dmwdr[NR], dcpdr[NR], d3vdrdp2[NR], d4vdrdp3[NR], d5vdrdp4[NR], dmOxTotdr[NR], denom;
        double d2cRefdr2[NR][NR], d3cdr2dt[NR][NR];
        int k;

        for (i=0, coeffCN=1.0; i<NY; i++) coeffCN += (fCN[i]-1.0)*s[NS+i];

        /* Convert input composition (r) to liquid moles (m)  */
        for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

        /* Compute moles and total moles of oxides */
        for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
        if (mOxTot == 0.0) return;

        /* Deal with the special case of FeO1.3 */
        mOx[NA] = 0.0;
        if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
            const double y = 0.3;
            mOx[iOxFeO1_3] = 0.0;
            if (iCmpFe2SiO4_6 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
            }
            if (iCmpFe2AlO4_1 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
            }
        }

        for (i=0; i<NR; i++) {
            dvdr[i]    = 0.0; d2vdrdt[i]  = 0.0; dcRefdr[i]  = 0.0; d2cdrdt[i]  = 0.0; dmwdr[i] = 0.0;
            dcpdr[i]   = 0.0; d3vdrdp2[i] = 0.0; d4vdrdp3[i] = 0.0; d5vdrdp4[i] = 0.0;
            for (j=0, dmOxTotdr[i]=0.0; j<nc; j++) dmOxTotdr[i] += (liquid[i+1].liqToOx)[j] - (liquid[0].liqToOx)[j];
            for (j=i; j<NR; j++) { d2cRefdr2[i][j] = 0.0; d3cdr2dt[i][j] = 0.0; }
        }

        for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
            v       += mOx[i]*bulkSystem[i].gk_v;
            dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
            cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
            dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
            cp      += mOx[i]*bulkSystem[i].gk_cp;
            d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
            mw      += mOx[i]*bulkSystem[i].mw;

            for (j=0; j<NR; j++) {
                dvdr[j]      += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_v;
                d2vdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dvdt;
                dmwdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].mw;
                dcpdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_cp;
  dcRefdr[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c/mOxTot
                - mOx[i]*bulkSystem[i].gk_c*dmOxTotdr[j]/(mOxTot*mOxTot);
  dcRefdr[j]   += (iOxAl2O3 != -1) ? ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
                - 2.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot) : 0.0;
  d2cdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt/mOxTot
                - mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotdr[j]/(mOxTot*mOxTot);
                d3vdrdp2[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
                d4vdrdp3[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
                d5vdrdp4[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
            }
            if (iCmpAl2O3 != -1) dcRefdr[iCmpAl2O3] += ((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*mOx[i]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot);

            for (j=0; j<NR; j++) {
  for (k=j; k<NR; k++) {
    d2cRefdr2[j][k] += -((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c*dmOxTotdr[k]/(mOxTot*mOxTot)
                -  ((liquid[k+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c*dmOxTotdr[j]/(mOxTot*mOxTot)
            + 2.0*mOx[i]*bulkSystem[i].gk_c*dmOxTotdr[j]*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot);
    d2cRefdr2[j][k] += (iOxAl2O3 != -1) ?
                - 2.0*((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot)
                - 2.0*((liquid[k+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot)
            + 6.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot*mOxTot)
            : 0.0;

    d3cdr2dt[j][k]  += -((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt*dmOxTotdr[k]/(mOxTot*mOxTot)
                - ((liquid[k+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt*dmOxTotdr[j]/(mOxTot*mOxTot)
            + 2.0*mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotdr[j]*dmOxTotdr[k]/(mOxTot*mOxTot*mOxTot);
  }
  if (iCmpAl2O3 != -1) d2cRefdr2[(j<iCmpAl2O3) ? j : iCmpAl2O3][(j>iCmpAl2O3) ? j : iCmpAl2O3] +=
    ((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
    - 2.0*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*mOx[i]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot);
            }
            if (iCmpAl2O3 != -1) d2cRefdr2[iCmpAl2O3][iCmpAl2O3] +=
  ((liquid[iCmpAl2O3+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
  - 2.0*mOx[i]*((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*bulkSystem[i].gk_cXal2o3*dmOxTotdr[iCmpAl2O3]/(mOxTot*mOxTot*mOxTot);

        }
        if (v == 0.0) return;

        alpha   = dvdt/v;
        v0      = v*exp(alpha*(t-tr));
        v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
        c       = cRef + (t-tr)*dcdt;
        v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
        v2      = d2vdp2;
        denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
        a        = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
        b        = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;
        sum     = a*a - 4.0*b;

        for (i=0; i<NR; i++) {
            double dalphadri   = d2vdrdt[i]/v - dvdt*dvdr[i]/(v*v);
            double dv0dri      = (dvdr[i] + v*dalphadri*(t-tr))*exp(alpha*(t-tr));
            double dv1Refdri   = -2.0*v*dvdr[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                           -v*v*(- 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
                        + 2.0*tr*alpha*dalphadri/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp));
            double dcdri       = dcRefdr[i] + (t-tr)*d2cdrdt[i];
            double dv1dri      = -2.0*v0*dv0dri*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                           -v0*v0*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)/(mw*mw*c*c*c)
                            + 2.0*t*alpha*dalphadri/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp));
            double ddenomdri   = 2.0*dv1Refdri*d3vdp3+2.0*v1Ref*d4vdrdp3[i]-6.0*d2vdp2*d3vdrdp2[i];
            double dadri       = (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdri*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)/denom
                         - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdri;
            double dbdri       = (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])/denom
                         - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdri;

            for (j=i; j<NR; j++) {
                double dalphadrj   = d2vdrdt[j]/v - dvdt*dvdr[j]/(v*v);
                double dv0drj       = (dvdr[j] + v*dalphadrj*(t-tr))*exp(alpha*(t-tr));
                double dv1Refdrj   = -2.0*v*dvdr[j]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                 -v*v*(- 1000.0*(cRef*dmwdr[j]+2.0*mw*dcRefdr[j])/(mw*mw*cRef*cRef*cRef)
                            + 2.0*tr*alpha*dalphadrj/(cp) - tr*alpha*alpha*dcpdr[j]/(cp*cp));
                double dcdrj       = dcRefdr[j] + (t-tr)*d2cdrdt[j];
                double dv1drj       = -2.0*v0*dv0drj*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                 -v0*v0*(- 1000.0*(c*dmwdr[j]+2.0*mw*dcdrj)/(mw*mw*c*c*c)
                                + 2.0*t*alpha*dalphadrj/(cp) - t*alpha*alpha*dcpdr[j]/(cp*cp));
                double ddenomdrj   = 2.0*dv1Refdrj*d3vdp3+2.0*v1Ref*d4vdrdp3[j]-6.0*d2vdp2*d3vdrdp2[j];
                double dadrj       = (d3vdrdp2[j]*d3vdp3+d2vdp2*d4vdrdp3[j]-dv1Refdrj*d4vdp4/2.0-v1Ref*d5vdrdp4[j]/2.0)/denom
               - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdrj;
                double dbdrj       = (d3vdrdp2[j]*d4vdp4/4.0+d2vdp2*d5vdrdp4[j]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[j])/denom
                         - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdrj;

                double d2alphadr2   = -d2vdrdt[i]*dvdr[j]/(v*v)-d2vdrdt[j]*dvdr[i]/(v*v)+2.0*dvdt*dvdr[i]*dvdr[j]/(v*v*v);
                double d2v0dr2      = dvdr[i]*dalphadrj*(t-tr)*exp(alpha*(t-tr))
          + dvdr[j]*dalphadri*(t-tr)*exp(alpha*(t-tr))
          + v*d2alphadr2*(t-tr)*exp(alpha*(t-tr))
          + v*dalphadri*pow(t-tr,2.0)*dalphadrj*exp(alpha*(t-tr));
                double d2v1Refdr2   = -2.0*dvdr[j]*dvdr[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                      - 2.0*v*dvdr[i]*(- 1000.0*(cRef*dmwdr[j]+2.0*mw*dcRefdr[j])/(mw*mw*cRef*cRef*cRef)
                                                + 2.0*tr*alpha*dalphadrj/(cp) - tr*alpha*alpha*dcpdr[j]/(cp*cp))
                            - 2.0*v*dvdr[j]*(- 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
                                                + 2.0*tr*alpha*dalphadri/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp))
                -v*v*(- 1000.0*(dcRefdr[j]*dmwdr[i]+2.0*dmwdr[j]*dcRefdr[i]+2.0*mw*d2cRefdr2[i][j])/(mw*mw*cRef*cRef*cRef)
                            + 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])*(2.0*dmwdr[j]*cRef+3.0*mw*dcRefdr[j])/(mw*mw*mw*cRef*cRef*cRef*cRef)
                            + 2.0*tr*alpha*d2alphadr2/(cp)
           + 2.0*tr*dalphadrj*dalphadri/(cp) - 2.0*tr*alpha*dalphadri*dcpdr[j]/(cp*cp)
           - 2.0*tr*alpha*dalphadrj*dcpdr[i]/(cp*cp) + 2.0*tr*alpha*alpha*dcpdr[i]*dcpdr[j]/(cp*cp*cp));
                double d2cdr2       = d2cRefdr2[i][j] + (t-tr)*d3cdr2dt[i][j];
                double d2v1dr2       = -2.0*dv0drj*dv0dri*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                      - 2.0*v0*d2v0dr2*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                      - 2.0*v0*dv0dri*(- 1000.0*(c*dmwdr[j]+2.0*mw*dcdrj)/(mw*mw*c*c*c)
                                                + 2.0*t*alpha*dalphadrj/(cp) - t*alpha*alpha*dcpdr[j]/(cp*cp))
                            - 2.0*v0*dv0drj*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)/(mw*mw*c*c*c)
                                                + 2.0*t*alpha*dalphadri/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp))
                -v0*v0*(- 1000.0*(dcdrj*dmwdr[i]+2.0*dmwdr[j]*dcdri+2.0*mw*d2cdr2)/(mw*mw*c*c*c)
                                + 1000.0*(c*dmwdr[i]+2.0*mw*dcdri)*(2.0*dmwdr[j]*c+3.0*mw*dcdrj)/(mw*mw*mw*c*c*c*c)
             + 2.0*t*alpha*d2alphadr2/(cp)
                                + 2.0*t*dalphadrj*dalphadri/(cp) - 2.0*t*alpha*dalphadri*dcpdr[j]/(cp*cp)
             - 2.0*t*alpha*dalphadrj*dcpdr[i]/(cp*cp) + 2.0*t*alpha*alpha*dcpdr[i]*dcpdr[j]/(cp*cp*cp));
                double d2denomdr2   = 2.0*d2v1Refdr2*d3vdp3 + 2.0*dv1Refdri*d4vdrdp3[j] + 2.0*dv1Refdrj*d4vdrdp3[i] - 6.0*d3vdrdp2[j]*d3vdrdp2[i];
                double d2adr2       = (d3vdrdp2[i]*d4vdrdp3[j]+d3vdrdp2[j]*d4vdrdp3[i]-d2v1Refdr2*d4vdp4/2.0-dv1Refdri*d5vdrdp4[j]/2.0-dv1Refdrj*d5vdrdp4[i]/2.0)/denom
                                        - (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdri*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)*ddenomdrj/(denom*denom)
               - (d3vdrdp2[j]*d3vdp3+d2vdp2*d4vdrdp3[j]-dv1Refdrj*d4vdp4/2.0-v1Ref*d5vdrdp4[j]/2.0)*ddenomdri/(denom*denom)
                           + 2.0*(d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*ddenomdri*ddenomdrj/(denom*denom*denom)
                           - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)*d2denomdr2/(denom*denom);
                double d2bdr2       = (d3vdrdp2[i]*d5vdrdp4[j]/4.0+d3vdrdp2[j]*d5vdrdp4[i]/4.0-2.0/3.0*d4vdrdp3[j]*d4vdrdp3[i])/denom
                           - (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])*ddenomdrj/(denom*denom)
               - (d3vdrdp2[j]*d4vdp4/4.0+d2vdp2*d5vdrdp4[j]/4.0-2.0*d4vdrdp3[j]*d3vdp3/3.0)*ddenomdri/(denom*denom)
            + 2.0*(d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)*ddenomdri*ddenomdrj/(denom*denom*denom)
            - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)*d2denomdr2/(denom*denom);
                double d3gIntdr2dp;

                if ((a == 0.0) && (b == 0.0)) {
                    d3gIntdr2dp  = d2v0dr2 + d2v1dr2*(p-pr);

                } else if ((a != 0.0) && (b == 0.0)) {
                    printf("*-->Exception in fillD3GDR2DP (liquid.c). a is not equal to zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
                    d3gIntdr2dp  = 0.0;

                } else if ((a == 0.0) && (b != 0.0)) {
                    printf("*-->Exception in fillD3GDR2DP (liquid.c). a is zero, b is not equal to zero.\n"); liqERRstate = ERR_A_ZERO;
                    d3gIntdr2dp  = 0.0;

                } else {
                    d3gIntdr2dp  = d3gdr2dpMAP(p, pr, v0, v1, d2vdp2, a, b, dv0dri, dv1dri, d3vdrdp2[i], dadri, dbdri, dv0drj, dv1drj, d3vdrdp2[j], dadrj, dbdrj,
                                                            d2v0dr2, d2v1dr2, 0.0, d2adr2, d2bdr2);

                }

                result[i][j] += coeffCN*d3gIntdr2dp;
  if (i != j) result[j][i] += coeffCN*d3gIntdr2dp;
            }
        }
    }

#endif /* USE_GHIORSO_KRESS_MODEL */

}

static void fillD3GDRDS2 (double r[NR], double s[NT], double t, double p, double result[NR][NT][NT]) {
    int i, j, k, l;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) for (j=0; j<NT; j++) memset(result[i][j], '\0', (size_t) NT*sizeof(double));

    /* Configurational entropy terms */
    for (j=0; j<NR; j++) {
        for (k=0; k<NS; k++) {
            for (l=k; l<NS; l++) {
                double config = 0.0;
                for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
                    config += nSpecies*(
                        (d2xSpeciesdrds[i][j][k]*dxSpeciesds[i][l] + dxSpeciesds[i][k]*d2xSpeciesdrds[i][j][l])/xSpecies[i]
                        - dxSpeciesdr[i][j]*dxSpeciesds[i][k]*dxSpeciesds[i][l]/(xSpecies[i]*xSpecies[i]) )
                        + dnSpeciesds[l]*(d2xSpeciesdrds[i][j][k]*log(xSpecies[i]) + dxSpeciesdr[i][j]*dxSpeciesds[i][k]/xSpecies[i])
                        + dnSpeciesds[k]*(d2xSpeciesdrds[i][j][l]*log(xSpecies[i]) + dxSpeciesdr[i][j]*dxSpeciesds[i][l]/xSpecies[i])
                        + d2nSpeciesds2[k][l]*(dxSpeciesdr[i][j]*(1.0 + log(xSpecies[i])));
                if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
                    config += d2nSpeciesds2[k][l]*(dxSpeciesdr[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O])))
                        + dnSpeciesds[k]*(d2xSpeciesdrds[nH2O][j][l]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                        + dxSpeciesdr[nH2O][j]*dxSpeciesds[nH2O][l]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O])))
                        + dnSpeciesds[l]*(d2xSpeciesdrds[nH2O][j][k]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                        + dxSpeciesds[nH2O][k]*dxSpeciesdr[nH2O][j]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O])))
                        + nSpecies*(
                            d2xSpeciesdrds[nH2O][j][k]*dxSpeciesds[nH2O][l]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]))
                            + dxSpeciesds[nH2O][k]*d2xSpeciesdrds[nH2O][j][l]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]))
                            - dxSpeciesds[nH2O][k]*dxSpeciesdr[nH2O][j]*dxSpeciesds[nH2O][l]
                                *(1.0/(xSpecies[nH2O]*xSpecies[nH2O]) - 1.0/((1.0-xSpecies[nH2O])*(1.0-xSpecies[nH2O]))) );
                result[j][k][l] += R*t*config;
                result[j][l][k]  = result[j][k][l];
            }
        }
    }
}

static void fillD3GDRDSDT (double r[NR], double s[NT], double t, double p, double result[NR][NT]) {
    int i, j, k;

    for (i=0; i<NR; i++) memset(result[i], '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) for (j=0; j<NS; j++) result[i][j] = -srs[i][j];

    /* Configurational entropy terms */
    for (j=0; j<NR; j++) {
        for (k=0; k<NS; k++) {
            double config = 0.0;
            for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
                config += nSpecies*(d2xSpeciesdrds[i][j][k]*log(xSpecies[i]) + dxSpeciesdr[i][j]*dxSpeciesds[i][k]/xSpecies[i])
                    + dnSpeciesds[k]*(dxSpeciesdr[i][j]*(1.0 + log(xSpecies[i])));
            if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
                config += dnSpeciesds[k]*dxSpeciesdr[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                    + nSpecies*(d2xSpeciesdrds[nH2O][j][k]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                    + dxSpeciesds[nH2O][k]*dxSpeciesdr[nH2O][j]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O])));
            result[j][k] += R*config;
        }
    }

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double d2gIntIVdrdt[NR];

        d2IntegralV_GKdrdt(r, s, t, p, d2gIntIVdrdt);
        for (i=0; i<NR; i++) for (j=0; j<NY; j++) result[i][NS+j] += (-ssCnModelParameters[j*nc+i+1].entropy + ssCnModelParameters[j*nc  +0].entropy) + (fCN[j]-1.0)*d2gIntIVdrdt[i];
    }
#endif /* USE_GHIORSO_KRESS_MODEL */
}

static void fillD3GDRDSDP (double r[NR], double s[NT], double t, double p, double result[NR][NT]) {
    int i, j;

    for (i=0; i<NR; i++) memset(result[i], '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) for (j=0; j<NS; j++) result[i][j] = vrs[i][j];

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double d2gIntIVdrdp[NR];

        d2IntegralV_GKdrdp(r, s, t, p, d2gIntIVdrdp);
        for (i=0; i<NR; i++) for (j=0; j<NY; j++) result[i][NS+j] += (fCN[j]-1.0)*d2gIntIVdrdp[i];
    }
#endif /* USE_GHIORSO_KRESS_MODEL */
}

static void fillD3GDS3 (double r[NR], double s[NT], double t, double p, double result[NT][NT][NT]) {
    int i, j, k, l;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NT; i++) for (j=0; j<NT; j++) for (k=0; k<NT; k++) result[i][j][k] = 0.0;

    /* Configurational entropy terms */
    for (j=0; j<NS; j++) {
        for (k=j; k<NS; k++) {
            for (l=k; l<NS; l++) {
                double config = 0.0;
                for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
                    config += -nSpecies*dxSpeciesds[i][j]*dxSpeciesds[i][k]*dxSpeciesds[i][l]/(xSpecies[i]*xSpecies[i])
                        + dnSpeciesds[l]*dxSpeciesds[i][j]*dxSpeciesds[i][k]/xSpecies[i]
                        + dnSpeciesds[k]*dxSpeciesds[i][j]*dxSpeciesds[i][l]/xSpecies[i]
                        + d2nSpeciesds2[k][l]*dxSpeciesds[i][j]*(1.0 + log(xSpecies[i]))
                        + dnSpeciesds[j]*dxSpeciesds[i][k]*dxSpeciesds[i][l]/xSpecies[i]
                        + d2nSpeciesds2[j][l]*dxSpeciesds[i][k]*(1.0 + log(xSpecies[i]))
                        + d2nSpeciesds2[j][k]*dxSpeciesds[i][l]*(1.0 + log(xSpecies[i]))
                        + d3nSpeciesds3[j][k][l]*xSpecies[i]*log(xSpecies[i]);
                if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
                    config += d3nSpeciesds3[j][k][l]*(xSpecies[nH2O]*log(xSpecies[nH2O]) + (1.0-xSpecies[nH2O])*log(1.0-xSpecies[nH2O]))
                        + d2nSpeciesds2[j][k]*dxSpeciesds[nH2O][l]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                        + d2nSpeciesds2[j][l]*dxSpeciesds[nH2O][k]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                        + dnSpeciesds[j]*dxSpeciesds[nH2O][k]*dxSpeciesds[nH2O][l]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]))
                        + d2nSpeciesds2[k][l]*dxSpeciesds[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                        + dnSpeciesds[k]*dxSpeciesds[nH2O][j]*dxSpeciesds[nH2O][l]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]))
                        + dnSpeciesds[l]*dxSpeciesds[nH2O][j]*dxSpeciesds[nH2O][k]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]))
                        - nSpecies*dxSpeciesds[nH2O][j]*dxSpeciesds[nH2O][k]*dxSpeciesds[nH2O][l]
                        *(1.0/(xSpecies[nH2O]*xSpecies[nH2O]) - 1.0/((1.0-xSpecies[nH2O])*(1.0-xSpecies[nH2O])));
                result[j][k][l] += R*t*config;
                result[k][j][l]  = result[j][k][l];
                result[l][j][k]  = result[j][k][l];
                result[l][k][j]  = result[j][k][l];
                result[j][l][k]  = result[j][k][l];
                result[k][l][j]  = result[j][k][l];
            }
        }
    }

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double yIV = 1.0;

        for (i=0; i<NY; i++) yIV -= (s[NS+i] > 0.0) ? s[NS+i] : DBL_EPSILON;

        for (i=0; i<NY; i++) {
            double y = (s[NS+i] > 0.0) ? s[NS+i] : DBL_EPSILON;
            for (j=i; j<NY; j++) {
                for (k=j; k<NY; k++) {
                    if ((i == j) && (j == k)) result[NS+i][NS+i][NS+i] = R*t*(-1.0/(y*y) + 1.0/(yIV*yIV));
    else                      result[NS+i][NS+j][NS+k] = R*t/(yIV*yIV);
    result[NS+j][NS+i][NS+k] = result[NS+i][NS+j][NS+k];
    result[NS+k][NS+i][NS+j] = result[NS+i][NS+j][NS+k];
    result[NS+k][NS+j][NS+i] = result[NS+i][NS+j][NS+k];
    result[NS+i][NS+k][NS+j] = result[NS+i][NS+j][NS+k];
    result[NS+j][NS+k][NS+i] = result[NS+i][NS+j][NS+k];

  }
            }
        }

    }
#endif /* USE_GHIORSO_KRESS_MODEL */
}

static void fillD3GDS2DT (double r[NR], double s[NT], double t, double p, double result[NT][NT]) {
    int i, j, k;

    for (i=0; i<NS; i++) memset(result[i], '\0', (size_t) NS*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) for (j=i; j<NS; j++) {
            result[i][j] = (i == j) ? -2.0*sss[i][i] : -sss[i][j];
            result[j][i] = result[i][j];
    }

    /* Configurational entropy terms */
    for (j=0; j<NS; j++) {
        for (k=j; k<NS; k++) {
            double config = 0.0;
            for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
                config += nSpecies*dxSpeciesds[i][j]*dxSpeciesds[i][k]/xSpecies[i]
                    + dnSpeciesds[k]*dxSpeciesds[i][j]*(1.0 + log(xSpecies[i]))
                    + dnSpeciesds[j]*dxSpeciesds[i][k]*(1.0 + log(xSpecies[i]))
                    + d2nSpeciesds2[j][k]*xSpecies[i]*log(xSpecies[i]);
            if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
                config += d2nSpeciesds2[j][k]*(xSpecies[nH2O]*log(xSpecies[nH2O]) + (1.0-xSpecies[nH2O])*log(1.0-xSpecies[nH2O]))
                    + dnSpeciesds[j]*dxSpeciesds[nH2O][k]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                    + dnSpeciesds[k]*dxSpeciesds[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                    + nSpecies*dxSpeciesds[nH2O][j]*dxSpeciesds[nH2O][k]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]));
            result[j][k] += R*config;
            result[k][j]  = result[j][k];
        }
    }

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double yIV = 1.0;

        for (i=0; i<NY; i++) yIV -= (s[NS+i] > 0.0) ? s[NS+i] : DBL_EPSILON;

        for (i=0; i<NY; i++) {
            double y = (s[NS+i] > 0.0) ? s[NS+i] : DBL_EPSILON;
            result[NS+i][NS+i] = R*(1.0/y + 1.0/yIV);
            for (j=i+1; j<NY; j++) { result[NS+i][NS+j] = R/yIV; result[NS+j][NS+i] = result[NS+i][NS+j]; }
        }

    }
#endif /* USE_GHIORSO_KRESS_MODEL */

}

static void fillD3GDS2DP (double r[NR], double s[NT], double t, double p, double result[NT][NT]) {
    int i, j;

    for (i=0; i<NS; i++) memset(result[i], '\0', (size_t) NS*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) for (j=i; j<NS; j++) {
            result[i][j] = (i == j) ? 2.0*vss[i][i] : vss[i][j];
            result[j][i] = result[i][j];
    }
}

static void fillD3GDSDT2 (double r[NR], double s[NT], double t, double p, double *result) {
    int i;

    memset(result, '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) result[i] = -cps[i]/t;

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double d2gIntIVdt2 = d2IntegralV_GKdT2(r, s, t, p);
        for (i=0; i<NY; i++) result[NS+i] = (fCN[i]-1.0)*d2gIntIVdt2;
    }
#endif /* USE_GHIORSO_KRESS_MODEL */

}

static void fillD3GDSDTDP (double r[NR], double s[NT], double t, double p, double *result) {
    int i;

    memset(result, '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) result[i] = dvdts[i];

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double d2gIntIVdtdp = d2IntegralV_GKdTdP(r, s, t, p);
        for (i=0; i<NY; i++) result[NS+i] = (fCN[i]-1.0)*d2gIntIVdtdp;
    }
#endif /* USE_GHIORSO_KRESS_MODEL */

}

static void fillD3GDSDP2 (double r[NR], double s[NT], double t, double p, double *result) {
    int i;

    memset(result, '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) result[i] = dvdps[i];

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        double d2gIntIVdp2 = d2IntegralV_GKdP2(r, s, t, p);
        for (i=0; i<NY; i++) result[NS+i] = (fCN[i]-1.0)*d2gIntIVdp2;
    }
#endif /* USE_GHIORSO_KRESS_MODEL */

}

static void fillD3GDRDT2 (double r[NR], double s[NT], double t, double p, double *result) {
    int i;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) result[i] = -cpr[i]/t;

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        const double pr        = 1.0;
        const double tr        = 1673.15;
        double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum, coeffCN;
        double dvdr[NR], d2vdrdt[NR], dcRefdr[NR], d2cdrdt[NR], dmwdr[NR], dcpdr[NR], d3vdrdp2[NR], d4vdrdp3[NR], d5vdrdp4[NR], dmOxTotdr[NR], denom;
        int j;

        if (fabs(p-pr) < 10.0*DBL_EPSILON) return;

        for (i=0, coeffCN=1.0; i<NY; i++) coeffCN += (fCN[i]-1.0)*s[NS+i];

        /* Convert input composition (r) to liquid moles (m)  */
        for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

        /* Compute moles and total moles of oxides */
        for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
        if (mOxTot == 0.0) return;

        /* Deal with the special case of FeO1.3 */
        mOx[NA] = 0.0;
        if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
            const double y = 0.3;
            mOx[iOxFeO1_3] = 0.0;
            if (iCmpFe2SiO4_6 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
            }
            if (iCmpFe2AlO4_1 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
            }
        }

        for (i=0; i<NR; i++) {
            dvdr[i]    = 0.0; d2vdrdt[i]  = 0.0; dcRefdr[i]  = 0.0; d2cdrdt[i]  = 0.0; dmwdr[i] = 0.0;
            dcpdr[i]   = 0.0; d3vdrdp2[i] = 0.0; d4vdrdp3[i] = 0.0; d5vdrdp4[i] = 0.0;
            for (j=0, dmOxTotdr[i]=0.0; j<nc; j++) dmOxTotdr[i] += (liquid[i+1].liqToOx)[j] - (liquid[0].liqToOx)[j];
        }

        for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
            v       += mOx[i]*bulkSystem[i].gk_v;
            dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
            cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
            dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
            cp      += mOx[i]*bulkSystem[i].gk_cp;
            d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
            mw      += mOx[i]*bulkSystem[i].mw;
            for (j=0; j<NR; j++) {
                dvdr[j]      += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_v;
                d2vdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dvdt;
                dmwdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].mw;
                dcpdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_cp;
  dcRefdr[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c/mOxTot
                - mOx[i]*bulkSystem[i].gk_c*dmOxTotdr[j]/(mOxTot*mOxTot);
                dcRefdr[j]   += (iOxAl2O3 != -1) ? ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
                - 2.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot) : 0.0;
  d2cdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt/mOxTot
                - mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotdr[j]/(mOxTot*mOxTot);
                d3vdrdp2[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
                d4vdrdp3[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
                d5vdrdp4[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
            }
            if (iCmpAl2O3 != -1) dcRefdr[iCmpAl2O3] += ((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*mOx[i]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot);
        }
        if (v == 0.0) return;

        alpha   = dvdt/v;
        v0      = v*exp(alpha*(t-tr));
        v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
        c       = cRef + (t-tr)*dcdt;
        v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
        v2      = d2vdp2;
        denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
        a        = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
        b        = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;
        sum     = a*a - 4.0*b;

        for (i=0; i<NR; i++) {
            double dalphadr = d2vdrdt[i]/v - dvdt*dvdr[i]/(v*v);
            double dv0dr    = (dvdr[i] + v*dalphadr*(t-tr))*exp(alpha*(t-tr));
            double dv1Refdr = -2.0*v*dvdr[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                                                -v*v*(- 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
                                + 2.0*tr*alpha*dalphadr/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp));
            double dcdr     = dcRefdr[i] + (t-tr)*d2cdrdt[i];
            double dv1dr    = -2.0*v0*dv0dr*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                                                -v0*v0*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdr)/(mw*mw*c*c*c)
              + 2.0*t*alpha*dalphadr/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp));
            double dv2dr    = d3vdrdp2[i];
            double ddenomdr = 2.0*dv1Refdr*d3vdp3+2.0*v1Ref*d4vdrdp3[i]-6.0*d2vdp2*d3vdrdp2[i];
            double dadr     = (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdr*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)/denom
                                            - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdr;
            double dbdr     = (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])/denom
                                            - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdr;
            double dv0dt    = v*alpha*exp(alpha*(t-tr));
            double dv1dt    = - 2.0*v0*v0*alpha*(1000.0/(mw*c*c) + t*alpha*alpha/(cp)) - v0*v0*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp));
            double d2v0drdt = v*dalphadr*exp(alpha*(t-tr)) + (dvdr[i] + v*dalphadr*(t-tr))*alpha*exp(alpha*(t-tr));

            double d2v1drdt = -2.0*(dv0dt*dv0dr + v0*d2v0drdt)*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                                            - 2.0*v0*dv0dr*(- 2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp))
                                            - 2.0*v0*dv0dt*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdr)/(mw*mw*c*c*c)
                    + 2.0*t*alpha*dalphadr/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp))
                - v0*v0*(- 1000.0*(dcdt*dmwdr[i]+2.0*mw*d2cdrdt[i])/(mw*mw*c*c*c) + 3000.0*(c*dmwdr[i]+2.0*mw*dcdr)*dcdt/(mw*mw*c*c*c*c)
                   + 2.0*alpha*dalphadr/(cp) - alpha*alpha*dcpdr[i]/(cp*cp));
            double d2v0dt2  = v*alpha*alpha*exp(alpha*(t-tr));
            double d2v1dt2  = - 2.0*(dv0dt*dv0dt+v0*d2v0dt2)*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                                            - 4.0*v0*dv0dt*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp))
                - v0*v0*(6000.0*dcdt*dcdt/(mw*c*c*c*c));
            double d3v0drdt2  = 2.0*v*dalphadr*alpha*exp(alpha*(t-tr)) + (dvdr[i] + v*dalphadr*(t-tr))*alpha*alpha*exp(alpha*(t-tr));

            double d3v1drdt2  = - 2.0*(d2v0dt2*dv0dr + dv0dt*d2v0drdt + dv0dt*d2v0drdt + v0*d3v0drdt2)*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                                                    - 2.0*(dv0dt*dv0dr + v0*d2v0drdt)*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp))
                                                    - 2.0*(dv0dt*dv0dr + v0*d2v0drdt)*(- 2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp)) - 2.0*v0*dv0dr*6000.0*dcdt*dcdt/(mw*c*c*c*c)
                                                    - 2.0*(dv0dt*dv0dt + v0*d2v0dt2)*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdr)/(mw*mw*c*c*c) + 2.0*t*alpha*dalphadr/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp))
                        - 2.0*v0*dv0dt*(- 1000.0*(dcdt*dmwdr[i]+2.0*mw*d2cdrdt[i])/(mw*mw*c*c*c) + 3000.0*(c*dmwdr[i]+2.0*mw*dcdr)*dcdt/(mw*mw*c*c*c*c)
                                                        + 2.0*alpha*dalphadr/(cp) - alpha*alpha*dcpdr[i]/(cp*cp))
                        - 2.0*v0*dv0dt*(- 1000.0*(dcdt*dmwdr[i]+2.0*mw*d2cdrdt[i])/(mw*mw*c*c*c) + 3000.0*(c*dmwdr[i]+2.0*mw*dcdr)*dcdt/(mw*mw*c*c*c*c)
                                                        + 2.0*alpha*dalphadr/(cp) - alpha*alpha*dcpdr[i]/(cp*cp))
        - v0*v0*(6000.0*(dcdt*dmwdr[i]+2.0*mw*d2cdrdt[i])*dcdt/(mw*mw*c*c*c*c) - 12000.0*(c*dmwdr[i]+2.0*mw*dcdr)*dcdt*dcdt/(mw*mw*c*c*c*c*c));

            double d3gIntdrdt2 = 0.0;

            if ((a == 0.0) && (b == 0.0)) {
                d3gIntdrdt2 = d3v0drdt2*(p-pr) + d3v1drdt2*(p-pr)*(p-pr)/2.0;

            } else if ((a != 0.0) && (b == 0.0)) {
  printf("*-->Exception in fillD3GDRDT2 (liquid.c). a is greater than zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
                d3gIntdrdt2 = 0.0;

            } else if ((a == 0.0) && (b != 0.0)) {
  printf("*-->Exception in fillD3GDRDT2 (liquid.c). a is zero, b is greater than zero.\n"); liqERRstate = ERR_A_ZERO;
                d3gIntdrdt2 = 0.0;

            } else if (sum > 0.0) {
                d3gIntdrdt2 = d3gdrdt2GMAP(p/10000.0, pr/10000.0,
                                                        v0,               v1*10000.0,    v2*10000.0*10000.0,    a*10000.0,    b*10000.0*10000.0,
           dv0dr,         dv1dr*10000.0, dv2dr*10000.0*10000.0, dadr*10000.0, dbdr*10000.0*10000.0,
           dv0dt,         dv1dt*10000.0,
                                                        d2v0drdt,   d2v1drdt*10000.0,
           d2v0dt2,     d2v1dt2*10000.0,
           d3v0drdt2, d3v1drdt2*10000.0)*10000.0;

            } else if (sum == 0.0) {
  printf("*-->Exception in fillD3GDRDT2 (liquid.c). a*a-4*b is equal to zero.\n"); liqERRstate = ERR_SUM_ZERO;
                d3gIntdrdt2 = 0.0;

            } else if(sum < 0.0) {
                d3gIntdrdt2 = d3gdrdt2LMAP(p/10000.0, pr/10000.0,
                                                        v0,               v1*10000.0,    v2*10000.0*10000.0,    a*10000.0,    b*10000.0*10000.0,
           dv0dr,         dv1dr*10000.0, dv2dr*10000.0*10000.0, dadr*10000.0, dbdr*10000.0*10000.0,
           dv0dt,         dv1dt*10000.0,
                                                        d2v0drdt,   d2v1drdt*10000.0,
           d2v0dt2,     d2v1dt2*10000.0,
           d3v0drdt2, d3v1drdt2*10000.0)*10000.0;

            }

            result[i] += coeffCN*d3gIntdrdt2;
        }
    }

#endif /* USE_GHIORSO_KRESS_MODEL */

}

static void fillD3GDRDTDP (double r[NR], double s[NT], double t, double p, double *result) {
    int i;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) result[i] = dvdtr[i];

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        const double pr        = 1.0;
        const double tr        = 1673.15;
        double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum, coeffCN;
        double dvdr[NR], d2vdrdt[NR], dcRefdr[NR], d2cdrdt[NR], dmwdr[NR], dcpdr[NR], d3vdrdp2[NR], d4vdrdp3[NR], d5vdrdp4[NR], dmOxTotdr[NR], denom;
        int j;

        for (i=0, coeffCN=1.0; i<NY; i++) coeffCN += (fCN[i]-1.0)*s[NS+i];

        /* Convert input composition (r) to liquid moles (m)  */
        for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

        /* Compute moles and total moles of oxides */
        for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
        if (mOxTot == 0.0) return;

        /* Deal with the special case of FeO1.3 */
        mOx[NA] = 0.0;
        if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
            const double y = 0.3;
            mOx[iOxFeO1_3] = 0.0;
            if (iCmpFe2SiO4_6 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
            }
            if (iCmpFe2AlO4_1 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
            }
        }

        for (i=0; i<NR; i++) {
            dvdr[i]    = 0.0; d2vdrdt[i]  = 0.0; dcRefdr[i]  = 0.0; d2cdrdt[i]  = 0.0; dmwdr[i] = 0.0;
            dcpdr[i]   = 0.0; d3vdrdp2[i] = 0.0; d4vdrdp3[i] = 0.0; d5vdrdp4[i] = 0.0;
            for (j=0, dmOxTotdr[i]=0.0; j<nc; j++) dmOxTotdr[i] += (liquid[i+1].liqToOx)[j] - (liquid[0].liqToOx)[j];
        }

        for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
            v       += mOx[i]*bulkSystem[i].gk_v;
            dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
            cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
            dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
            cp      += mOx[i]*bulkSystem[i].gk_cp;
            d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
            mw      += mOx[i]*bulkSystem[i].mw;
            for (j=0; j<NR; j++) {
                dvdr[j]      += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_v;
                d2vdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dvdt;
                dmwdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].mw;
                dcpdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_cp;
  dcRefdr[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c/mOxTot
                - mOx[i]*bulkSystem[i].gk_c*dmOxTotdr[j]/(mOxTot*mOxTot);
  dcRefdr[j]   += (iOxAl2O3 != -1) ? ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
                - 2.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot) : 0.0;
  d2cdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt/mOxTot
                - mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotdr[j]/(mOxTot*mOxTot);
                d3vdrdp2[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
                d4vdrdp3[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
                d5vdrdp4[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
            }
            if (iCmpAl2O3 != -1) dcRefdr[iCmpAl2O3] += ((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*mOx[i]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot);
        }
        if (v == 0.0) return;

        alpha   = dvdt/v;
        v0      = v*exp(alpha*(t-tr));
        v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
        c       = cRef + (t-tr)*dcdt;
        v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
        v2      = d2vdp2;
        denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
        a        = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
        b        = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;
        sum     = a*a - 4.0*b;

        for (i=0; i<NR; i++) {
            double dalphadr = d2vdrdt[i]/v - dvdt*dvdr[i]/(v*v);
            double dv0dr    = (dvdr[i] + v*dalphadr*(t-tr))*exp(alpha*(t-tr));
            double dv1Refdr = -2.0*v*dvdr[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                                                -v*v*(- 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
                                + 2.0*tr*alpha*dalphadr/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp));
            double dcdr     = dcRefdr[i] + (t-tr)*d2cdrdt[i];
            double dv1dr    = -2.0*v0*dv0dr*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                                                -v0*v0*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdr)/(mw*mw*c*c*c)
              + 2.0*t*alpha*dalphadr/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp));
            double dv2dr    = d3vdrdp2[i];
            double ddenomdr = 2.0*dv1Refdr*d3vdp3+2.0*v1Ref*d4vdrdp3[i]-6.0*d2vdp2*d3vdrdp2[i];
            double dadr     = (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdr*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)/denom
                                            - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdr;
            double dbdr     = (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])/denom
                                            - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdr;
            double dv0dt    = v*alpha*exp(alpha*(t-tr));
            double dv1dt    = - 2.0*v0*v0*alpha*(1000.0/(mw*c*c) + t*alpha*alpha/(cp)) - v0*v0*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp));

            double d2v0drdt = v*dalphadr*exp(alpha*(t-tr)) + (dvdr[i] + v*dalphadr*(t-tr))*alpha*exp(alpha*(t-tr));
            double d2v1drdt = -2.0*(dv0dt*dv0dr + v0*d2v0drdt)*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                                            - 2.0*v0*dv0dr*(- 2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp))
                                            - 2.0*v0*dv0dt*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdr)/(mw*mw*c*c*c)
                    + 2.0*t*alpha*dalphadr/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp))
                - v0*v0*(- 1000.0*(dcdt*dmwdr[i]+2.0*mw*d2cdrdt[i])/(mw*mw*c*c*c) + 3000.0*(c*dmwdr[i]+2.0*mw*dcdr)*dcdt/(mw*mw*c*c*c*c)
                   + 2.0*alpha*dalphadr/(cp) - alpha*alpha*dcpdr[i]/(cp*cp));
            double d3gdrdtdp;

            if ((a == 0.0) && (b == 0.0)) {
                d3gdrdtdp = d2v0drdt + d2v1drdt*(p-pr);

            } else if ((a != 0.0) && (b == 0.0)) {
  printf("*-->Exception in fillD2GDRDT (liquid.c). a is greater than zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
                d3gdrdtdp = 0.0;

            } else if ((a == 0.0) && (b != 0.0)) {
  printf("*-->Exception in fillD2GDRDT (liquid.c). a is zero, b is greater than zero.\n"); liqERRstate = ERR_A_ZERO;
                d3gdrdtdp = 0.0;

            } else {
                d3gdrdtdp = d3gdrdtdpMAP(p, pr, v0, v1, v2, a, b, dv0dr, dv1dr, dv2dr, dadr, dbdr, dv0dt, dv1dt, d2v0drdt, d2v1drdt);

            }

            result[i] += coeffCN*d3gdrdtdp;
        }
    }

#endif /* USE_GHIORSO_KRESS_MODEL */

}

static void fillD3GDRDP2 (double r[NR], double s[NT], double t, double p, double *result) {
    int i;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) result[i] = dvdpr[i];

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        const double pr        = 1.0;
        const double tr        = 1673.15;
        double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum, coeffCN;
        double dvdr[NR], d2vdrdt[NR], dcRefdr[NR], d2cdrdt[NR], dmwdr[NR], dcpdr[NR], d3vdrdp2[NR], d4vdrdp3[NR], d5vdrdp4[NR], dmOxTotdr[NR], denom;
        int j;

        for (i=0, coeffCN=1.0; i<NY; i++) coeffCN += (fCN[i]-1.0)*s[NS+i];

        /* Convert input composition (r) to liquid moles (m)  */
        for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

        /* Compute moles and total moles of oxides */
        for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
        if (mOxTot == 0.0) return;

        /* Deal with the special case of FeO1.3 */
        mOx[NA] = 0.0;
        if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
            const double y = 0.3;
            mOx[iOxFeO1_3] = 0.0;
            if (iCmpFe2SiO4_6 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
            }
            if (iCmpFe2AlO4_1 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
            }
        }

        for (i=0; i<NR; i++) {
            dvdr[i]    = 0.0; d2vdrdt[i]  = 0.0; dcRefdr[i]  = 0.0; d2cdrdt[i]  = 0.0; dmwdr[i] = 0.0;
            dcpdr[i]   = 0.0; d3vdrdp2[i] = 0.0; d4vdrdp3[i] = 0.0; d5vdrdp4[i] = 0.0;
            for (j=0, dmOxTotdr[i]=0.0; j<nc; j++) dmOxTotdr[i] += (liquid[i+1].liqToOx)[j] - (liquid[0].liqToOx)[j];
        }

        for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
            v       += mOx[i]*bulkSystem[i].gk_v;
            dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
            cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
            dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
            cp      += mOx[i]*bulkSystem[i].gk_cp;
            d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
            mw      += mOx[i]*bulkSystem[i].mw;
            for (j=0; j<NR; j++) {
                dvdr[j]      += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_v;
                d2vdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dvdt;
                dmwdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].mw;
                dcpdr[j]     += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_cp;
  dcRefdr[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_c/mOxTot
                - mOx[i]*bulkSystem[i].gk_c*dmOxTotdr[j]/(mOxTot*mOxTot);
  dcRefdr[j]   += (iOxAl2O3 != -1) ? ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot)
                - 2.0*mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3*dmOxTotdr[j]/(mOxTot*mOxTot*mOxTot) : 0.0;
  d2cdrdt[j]   += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*bulkSystem[i].gk_dcdt/mOxTot
                - mOx[i]*bulkSystem[i].gk_dcdt*dmOxTotdr[j]/(mOxTot*mOxTot);
                d3vdrdp2[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
                d4vdrdp3[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
                d5vdrdp4[j]  += ((liquid[j+1].liqToOx)[i]-(liquid[0].liqToOx)[i])*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
            }
            if (iCmpAl2O3 != -1) dcRefdr[iCmpAl2O3] += ((liquid[iCmpAl2O3+1].liqToOx)[iOxAl2O3]-(liquid[0].liqToOx)[iOxAl2O3])*mOx[i]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot);
        }
        if (v == 0.0) return;

        alpha   = dvdt/v;
        v0      = v*exp(alpha*(t-tr));
        v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
        c       = cRef + (t-tr)*dcdt;
        v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
        v2      = d2vdp2;
        denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
        a        = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
        b        = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;
        sum     = a*a - 4.0*b;

        for (i=0; i<NR; i++) {
            double dalphadr   = d2vdrdt[i]/v - dvdt*dvdr[i]/(v*v);
            double dv0dr      = (dvdr[i] + v*dalphadr*(t-tr))*exp(alpha*(t-tr));
            double dv1Refdr   = -2.0*v*dvdr[i]*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp))
                                                    -v*v*(- 1000.0*(cRef*dmwdr[i]+2.0*mw*dcRefdr[i])/(mw*mw*cRef*cRef*cRef)
              + 2.0*tr*alpha*dalphadr/(cp) - tr*alpha*alpha*dcpdr[i]/(cp*cp));
            double dcdr       = dcRefdr[i] + (t-tr)*d2cdrdt[i];
            double dv1dr      = -2.0*v0*dv0dr*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                                                    -v0*v0*(- 1000.0*(c*dmwdr[i]+2.0*mw*dcdr)/(mw*mw*c*c*c)
                + 2.0*t*alpha*dalphadr/(cp) - t*alpha*alpha*dcpdr[i]/(cp*cp));
            double dv2dr      = d3vdrdp2[i];
            double ddenomdr   = 2.0*dv1Refdr*d3vdp3+2.0*v1Ref*d4vdrdp3[i]-6.0*d2vdp2*d3vdrdp2[i];
            double dadr       = (d3vdrdp2[i]*d3vdp3+d2vdp2*d4vdrdp3[i]-dv1Refdr*d4vdp4/2.0-v1Ref*d5vdrdp4[i]/2.0)/denom
                                                - (d2vdp2*d3vdp3-v1Ref*d4vdp4/2.0)/(denom*denom)*ddenomdr;
            double dbdr       = (d3vdrdp2[i]*d4vdp4/4.0+d2vdp2*d5vdrdp4[i]/4.0-2.0/3.0*d3vdp3*d4vdrdp3[i])/denom
                                                - (d2vdp2*d4vdp4/4.0-d3vdp3*d3vdp3/3.0)/(denom*denom)*ddenomdr;
            double d3gdrdp2;

            if ((a == 0.0) && (b == 0.0)) {
                d3gdrdp2 = dv2dr;

            } else if ((a != 0.0) && (b == 0.0)) {
  printf("*-->Exception in fillD3GDRDP2 (liquid.c). a is greater than zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
                d3gdrdp2 = 0.0;

            } else if ((a == 0.0) && (b != 0.0)) {
  printf("*-->Exception in fillD3GDRDP2 (liquid.c). a is zero, b is greater than zero.\n"); liqERRstate = ERR_A_ZERO;
                d3gdrdp2 = 0.0;

            } else {
                d3gdrdp2 = d3gdrdp2MAP(p, pr, v0, v1, v2, a, b, dv0dr, dv1dr, dv2dr, dadr, dbdr);

            }

            result[i] += coeffCN*d3gdrdp2;
        }
    }

#endif /* USE_GHIORSO_KRESS_MODEL */

}

static double fillD3GDT3 (double r[NR], double s[NT], double t, double p) {
    double result;
    int i;

    /* Taylor expansion and standard state terms */
    result = CPconst/(t*t) - DCPDTconst/t;
    for (i=0; i<NR; i++) result += cpr[i]*r[i]/(t*t) - dcpdtr[i]*r[i]/t;
    for (i=0; i<NS; i++) result += cps[i]*s[i]/(t*t) - dcpdts[i]*s[i]/t;

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        int j;
        const double pr         = 1.0;
        const double tr         = 1673.15;
        double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum, denom, coeffCN;
        double dv0dt, dv1dt, d2v0dt2, d2v1dt2, d3v0dt3, d3v1dt3, d3gIntdt3 = 0.0;

        if (fabs(p-pr) < 10.0*DBL_EPSILON) return 0.0;

        for (i=0, coeffCN=1.0; i<NY; i++) coeffCN += (fCN[i]-1.0)*s[NS+i];

        /* Convert input composition (r) to liquid moles (m)  */
        for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

        /* Compute moles and total moles of oxides */
        for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
        if (mOxTot == 0.0) return result;

        /* Deal with the special case of FeO1.3 */
        mOx[NA] = 0.0;
        if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
            const double y = 0.3;
            mOx[iOxFeO1_3] = 0.0;
            if (iCmpFe2SiO4_6 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
            }
            if (iCmpFe2AlO4_1 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
            }
        }

        for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
            v       += mOx[i]*bulkSystem[i].gk_v;
            dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
            cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
            dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
            cp      += mOx[i]*bulkSystem[i].gk_cp;
            d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
            mw      += mOx[i]*bulkSystem[i].mw;
        }
        if (v == 0.0) return result;

        alpha   = dvdt/v;
        v0      = v*exp(alpha*(t-tr));
        v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
        c       = cRef + (t-tr)*dcdt;
        v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
        v2      = d2vdp2;
        denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
        a        = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
        b        = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;
        sum     = a*a - 4.0*b;

        dv0dt   = alpha*v*exp(alpha*(t-tr));
        dv1dt   = - 2.0*v0*dv0dt*(1000.0/(mw*c*c) + t*alpha*alpha/(cp)) - v0*v0*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp));
        d2v0dt2 = alpha*alpha*v*exp(alpha*(t-tr));
        d2v1dt2 = - 2.0*(dv0dt*dv0dt+v0*d2v0dt2)*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                            - 4.0*v0*dv0dt*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp))
        - v0*v0*(6000.0*dcdt*dcdt/(mw*c*c*c*c));
        d3v0dt3   = alpha*alpha*alpha*v*exp(alpha*(t-tr));
        d3v1dt3   = - 2.0*(3.0*dv0dt*d2v0dt2 + v0*d3v0dt3)*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                                - 2.0*(dv0dt*dv0dt+v0*d2v0dt2)*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp))
                                - 4.0*(dv0dt*dv0dt + v0*d2v0dt2)*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp))
          - 4.0*v0*dv0dt*6000.0*dcdt*dcdt/(mw*c*c*c*c)
          - 2.0*v0*dv0dt*(6000.0*dcdt*dcdt/(mw*c*c*c*c)) - v0*v0*(-24000.0*dcdt*dcdt*dcdt/(mw*c*c*c*c*c));

        if ((a == 0.0) && (b == 0.0)) {
            d3gIntdt3    = d3v0dt3*(p-pr) + d3v1dt3*(p-pr)*(p-pr)/2.0;

        } else if ((a != 0.0) && (b == 0.0)) {
            printf("*-->Exception in fillD3GDT3 (liquid.c). a is greater than zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
            d3gIntdt3       = 0.0;

        } else if ((a == 0.0) && (b != 0.0)) {
            printf("*-->Exception in fillD3GDT3 (liquid.c). a is zero, b is greater than zero.\n"); liqERRstate = ERR_A_ZERO;
            d3gIntdt3       = 0.0;

        } else if (sum > 0.0) {
            d3gIntdt3 = d3gdt3GMAP(p/10000.0, pr/10000.0,
                             v0,           v1*10000.0, v2*10000.0*10000.0, a*10000.0, b*10000.0*10000.0,
                dv0dt,     dv1dt*10000.0,
                d2v0dt2, d2v1dt2*10000.0,
                d3v0dt3, d3v1dt3*10000.0)*10000.0;

        } else if (sum == 0.0) {
            printf("*-->Exception in fillD3GDT3 (liquid.c). a*a-4*b is equal to zero.\n"); liqERRstate = ERR_SUM_ZERO;
            d3gIntdt3       = 0.0;

        } else if(sum < 0.0) {
            d3gIntdt3 = d3gdt3LMAP(p/10000.0, pr/10000.0,
                             v0,           v1*10000.0, v2*10000.0*10000.0, a*10000.0, b*10000.0*10000.0,
                dv0dt,     dv1dt*10000.0,
                d2v0dt2, d2v1dt2*10000.0,
                d3v0dt3, d3v1dt3*10000.0)*10000.0;

        }

        result += coeffCN*d3gIntdt3;
    }

#endif /* USE_GHIORSO_KRESS_MODEL */

    return result;
}

static double fillD3GDT2DP (double r[NR], double s[NT], double t, double p) {
    double result;
    int i;

    /* Taylor expansion and standard state terms */
    result = D2VDT2const;
    for (i=0; i<NR; i++) result += d2vdt2r[i]*r[i];
    for (i=0; i<NS; i++) result += d2vdt2s[i]*s[i];

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        int j;
        const double pr         = 1.0;
        const double tr         = 1673.15;
        double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum, denom, coeffCN;
        double dv0dt, dv1dt, d2v0dt2, d2v1dt2, d3gdt2dp;

        for (i=0, coeffCN=1.0; i<NY; i++) coeffCN += (fCN[i]-1.0)*s[NS+i];

        /* Convert input composition (r) to liquid moles (m)  */
        for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

        /* Compute moles and total moles of oxides */
        for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
        if (mOxTot == 0.0) return 0.0;

        /* Deal with the special case of FeO1.3 */
        mOx[NA] = 0.0;
        if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
            const double y = 0.3;
            mOx[iOxFeO1_3] = 0.0;
            if (iCmpFe2SiO4_6 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
            }
            if (iCmpFe2AlO4_1 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
            }
        }

        for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
            v       += mOx[i]*bulkSystem[i].gk_v;
            dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
            cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
            dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
            cp      += mOx[i]*bulkSystem[i].gk_cp;
            d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
            mw      += mOx[i]*bulkSystem[i].mw;
        }
        if (v == 0.0) return result;

        alpha   = dvdt/v;
        v0      = v*exp(alpha*(t-tr));
        v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
        c       = cRef + (t-tr)*dcdt;
        v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
        v2      = d2vdp2;
        denom   = 2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2;
        a        = (denom != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /denom : 0.0;
        b        = (denom != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/denom : 0.0;
        sum     = a*a - 4.0*b;

        dv0dt   = alpha*v*exp(alpha*(t-tr));
        dv1dt   = - 2.0*v0*dv0dt*(1000.0/(mw*c*c) + t*alpha*alpha/(cp)) - v0*v0*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp));
        d2v0dt2 = alpha*alpha*v*exp(alpha*(t-tr));
        d2v1dt2 = - 2.0*(dv0dt*dv0dt+v0*d2v0dt2)*(1000.0/(mw*c*c) + t*alpha*alpha/(cp))
                            - 4.0*v0*dv0dt*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp))
        - v0*v0*(6000.0*dcdt*dcdt/(mw*c*c*c*c));
        d3gdt2dp = 0.0;

        if ((a == 0.0) && (b == 0.0)) {
            d3gdt2dp    += d2v0dt2 + d2v1dt2*(p-pr);

        } else if ((a != 0.0) && (b == 0.0)) {
            printf("*-->Exception in fillD3GDT2DP (liquid.c). a is greater than zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
            d3gdt2dp       = 0.0;

        } else if ((a == 0.0) && (b != 0.0)) {
            printf("*-->Exception in fillD3GDT2DP (liquid.c). a is zero, b is greater than zero.\n"); liqERRstate = ERR_A_ZERO;
            d3gdt2dp       = 0.0;

        } else {
            d3gdt2dp = d3gdt2dpMAP(p, pr, v0, v1, v2, a, b, dv0dt, dv1dt, d2v0dt2, d2v1dt2);

        }

        result += coeffCN*d3gdt2dp;
    }

#endif /* USE_GHIORSO_KRESS_MODEL */

    return result;
}

static double fillD3GDTDP2 (double r[NR], double s[NT], double t, double p) {
    double result;
    int i;

    /* Taylor expansion and standard state terms */
    result = D2VDTDPconst;
    for (i=0; i<NR; i++) result += d2vdtdpr[i]*r[i];
    for (i=0; i<NS; i++) result += d2vdtdps[i]*s[i];

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        int j;
        const double pr         = 1.0;
        const double tr         = 1673.15;
        double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum, coeffCN;
        double dv0dt, dv1dt, d3gdtdp2;

        for (i=0, coeffCN=1.0; i<NY; i++) coeffCN += (fCN[i]-1.0)*s[NS+i];

        /* Convert input composition (r) to liquid moles (m)  */
        for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

        /* Compute moles and total moles of oxides */
        for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
        if (mOxTot == 0.0) return result;

        /* Deal with the special case of FeO1.3 */
        mOx[NA] = 0.0;
        if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
            const double y = 0.3;
            mOx[iOxFeO1_3] = 0.0;
            if (iCmpFe2SiO4_6 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
            }
            if (iCmpFe2AlO4_1 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
            }
        }

        for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
            v       += mOx[i]*bulkSystem[i].gk_v;
            dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
            cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
            dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
            cp      += mOx[i]*bulkSystem[i].gk_cp;
            d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
            mw      += mOx[i]*bulkSystem[i].mw;
        }
        if (v == 0.0) return result;

        alpha   = dvdt/v;
        v0      = v*exp(alpha*(t-tr));
        v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
        c       = cRef + (t-tr)*dcdt;
        v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
        v2      = d2vdp2;
        a        = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;
        b        = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;
        sum     = a*a - 4.0*b;
        dv0dt   = alpha*v0;
        dv1dt   = - 2.0*v0*v0*alpha*(1000.0/(mw*c*c) + t*alpha*alpha/(cp)) - v0*v0*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(cp));

        if ((a == 0.0) && (b == 0.0)) {
            d3gdtdp2 = dv1dt;

        } else if ((a != 0.0) && (b == 0.0)) {
            printf("*-->Exception in fillD3GDTDP2 (liquid.c). a is greater than zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
            d3gdtdp2 = 0.0;

        } else if ((a == 0.0) && (b != 0.0)) {
            printf("*-->Exception in fillD3GDTDP2 (liquid.c). a is zero, b is greater than zero.\n"); liqERRstate = ERR_A_ZERO;
            d3gdtdp2 = 0.0;

        } else {
            d3gdtdp2 = d3gdtdp2MAP(p, pr, v0, v1, v2, a, b, dv0dt, dv1dt);

        }

        result += coeffCN*d3gdtdp2;
    }

#endif /* USE_GHIORSO_KRESS_MODEL */

    return result;
}

static double fillD3GDP3 (double r[NR], double s[NT], double t, double p) {
    double result;
    int i;

    /* Taylor expansion and standard state terms */
    result = D2VDP2const;
    for (i=0; i<NR; i++) result += d2vdp2r[i]*r[i];
    for (i=0; i<NS; i++) result += d2vdp2s[i]*s[i];

#ifdef USE_GHIORSO_KRESS_MODEL
    {
        int j;
        const double pr         = 1.0;
        const double tr         = 1673.15;
        double m[NA], mOx[NA+1], mOxTot, v, dvdt, c, cRef, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum, coeffCN;
        double d3gdp3;

        for (i=0, coeffCN=1.0; i<NY; i++) coeffCN += (fCN[i]-1.0)*s[NS+i];

        /* Convert input composition (r) to liquid moles (m)  */
        for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

        /* Compute moles and total moles of oxides */
        for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
        if (mOxTot == 0.0) return result;

        /* Deal with the special case of FeO1.3 */
        mOx[NA] = 0.0;
        if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
            const double y = 0.3;
            mOx[iOxFeO1_3] = 0.0;
            if (iCmpFe2SiO4_6 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
            }
            if (iCmpFe2AlO4_1 != -1) {
        mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
        mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
            }
        }

        for (i=0, v=0.0, dvdt=0.0, cRef=0.0, dcdt=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
            v       += mOx[i]*bulkSystem[i].gk_v;
            dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
            cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
            dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
            cp      += mOx[i]*bulkSystem[i].gk_cp;
            d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
            d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
            d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
            mw      += mOx[i]*bulkSystem[i].mw;
        }
        if (v == 0.0) return result;

        alpha   = dvdt/v;
        v0      = v*exp(alpha*(t-tr));
        v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
        c       = cRef + (t-tr)*dcdt;
        v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
        v2      = d2vdp2;
        a        = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;
        b        = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;
        sum     = a*a - 4.0*b;

        if ((a == 0.0) && (b == 0.0)) {
            d3gdp3 = v2;

        } else if ((a != 0.0) && (b == 0.0)) {
            printf("*-->Exception in fillD3GDP3 (liquid.c). a is greater than zero, b is zero.\n"); liqERRstate = ERR_B_ZERO;
            d3gdp3 = 0.0;

        } else if ((a == 0.0) && (b != 0.0)) {
            printf("*-->Exception in fillD3GDP3 (liquid.c). a is zero, b is greater than zero.\n"); liqERRstate = ERR_A_ZERO;
            d3gdp3 = 0.0;

        } else {
            d3gdp3 = d3gdp3MAP(p, pr, v0, v1, v2, a, b);

        }

        result += coeffCN*d3gdp3;
    }

#endif /* USE_GHIORSO_KRESS_MODEL */

    return result;
}

static void fillD3GDS2DW (double r[NR], double s[NT], double t, double p, double result[NT][NT][3*NP]) {
    int i, j, k, l, m, n, ii, iii;

    /*******************************
    * Parameters: NW WH(), NE H() *
    *             NW WS(), NE S() *
    *             NW WV(), NE V() *
    *******************************/
    for (k=0; k<NT; k++) for (j=0; j<NT; j++) memset(result[k][j], '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
    * NW W parameters solution are first *
    **************************************/
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        for (ii=0; ii<NS; ii++) {
            for (iii=ii; iii<NS; iii++) {
                m = ii*NS+(iii+1)-(ii+1)*(ii+2)/2+(ii+1)-1+NR*(NR-1)/2+NR+NR*NS;
                if (taylorCoeff[n+NE][1+NR+NS+m] != 0.0) {
                    result[ii][iii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m];
                    result[ii][iii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m];
#ifndef USE_GHIORSO_KRESS_MODEL
                    result[ii][iii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m];
#endif /* USE_GHIORSO_KRESS_MODEL */
                    result[iii][ii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m];
                    result[iii][ii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m];
#ifndef USE_GHIORSO_KRESS_MODEL
                    result[iii][ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m];
#endif /* USE_GHIORSO_KRESS_MODEL */
                }
            }
        }
    }
}

static void fillD3GDSDTDW (double r[NR], double s[NT], double t, double p, double result[NT][3*NP]) {
    int i, j, k, l, m, n, ii;

    /*******************************
    * Parameters: NW WH(), NE H() *
    *             NW WS(), NE S() *
    *             NW WV(), NE V() *
    *******************************/
    for (ii=0; ii<NT; ii++) memset(result[ii], '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
    * NW W parameters solution are first *
    **************************************/
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        for (ii=0; ii<NS; ii++) {
            for (j=0, m=0; j<NR; j++) {
                m += NR - j;
                m += ii;
                result[ii][  NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*r[j];
                m += NS - ii;
            }
            result[ii][  NP+n] += -taylorCoeff[n+NE][1+NR+ii];
            for (k=ii; k<NS; k++) {
                m = ii*NS+(k+1)-(ii+1)*(ii+2)/2+(ii+1)-1+NR*(NR-1)/2+NR+NR*NS;
                result[ii][  NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*s[k];
            }
            for (j=0; j<=ii; j++) {
                m = j*NS+(ii+1)-(j+1)*(j+2)/2+(j+1)-1+NR*(NR-1)/2+NR+NR*NS;
                result[ii][  NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*s[j];
            }
        }
    }

    /**************************************
    * NE standard state terms are second *
    **************************************/
    for (i=0; i<NE; i++, n++) {
        for (ii=0; ii<NS; ii++) result[ii][  NP+n] += -taylorCoeff[i][1+NR+ii];
    }
}

static void fillD3GDRDTDW (double r[NR], double s[NT], double t, double p, double result[NR][3*NP]) {
    int i, j, k, l, m, n, ii;

    /*******************************
    * Parameters: NW WH(), NE H() *
    *             NW WS(), NE S() *
    *             NW WV(), NE V() *
    *******************************/
    for (ii=0; ii<NR; ii++) memset(result[ii], '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
    * NW W parameters solution are first *
    **************************************/
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        for (ii=0; ii<NR; ii++) {
            result[ii][  NP+n] += -taylorCoeff[n+NE][1+ii];
            m = 0;
            for (j=0, m=0; j<NR; j++) {
                for (k=j; k<NR; k++, m++) {
                    if (j == ii) result[ii][  NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*r[k];
                    if (k == ii) result[ii][  NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*r[j];
                }
                for (k=0; k<NS; k++, m++) {
                    if (j == ii) result[ii][  NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*s[k];
                }
            }
        }
    }

    /**************************************
    * NE standard state terms are second *
    **************************************/
    for (i=0; i<NE; i++, n++) {
        for (ii=0; ii<NR; ii++) result[ii][  NP+n] += -taylorCoeff[i][1+ii];
    }
}

static void fillD3GDRDSDW (double r[NR], double s[NT], double t, double p, double result[NR][NT][3*NP]) {
    int i, l, m, n, ii, iii;

    /*******************************
    * Parameters: NW WH(), NE H() *
    *             NW WS(), NE S() *
    *             NW WV(), NE V() *
    *******************************/
    for (iii=0; iii<NR; iii++) for (ii=0; ii<NT; ii++) memset(result[iii][ii], '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
    * NW W parameters solution are first *
    **************************************/
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        for (ii=0; ii<NR; ii++) {
            for (iii=0; iii<NS; iii++) {
                m = (ii+1)*NR-ii*(ii+1)/2+iii*(ii+1)+(NS-iii)*ii;
                result[ii][iii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m];
                result[ii][iii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m];
#ifndef USE_GHIORSO_KRESS_MODEL
                result[ii][iii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m];
#endif /* USE_GHIORSO_KRESS_MODEL */
            }
        }
    }
}

static void initialGuessOrdering(double r[NR], double s[NT]) {
    int i;
    static double *sCorr;
    double factor = 1.0;
#ifndef TESTDYNAMICLIB
    double rSum, s0;
#endif

    /* Update to match Thermoengine */
    double denom = 1.0 + (((r[9] > 0.0) && (r[13] > 0.0)) ? 1.0 : 0.0); /* [ 0] SiO2 [10] CaO [14] CO2 */
    for (i=0; i<NS; i++)  s[i] = r[13]/denom; // if CaO  is present 1/denom of all CO2 is CaCO3

    if (NT == 0) return;

#ifdef DEBUG
    printf("Call to initialGuessOrdering in liquid.c\n");
#endif
    if (sCorr == NULL) sCorr = (double *) malloc((size_t) NS*sizeof(double));

    /*for (i=0; i<NS; i++)  s[i] = 0.0;
    for (i=NS; i<NT; i++) s[i] = 1.0/(((double) NY)+1.0);*/

#ifndef TESTDYNAMICLIB
    /* A tweak for alkali-rich liquids where SiO2 tends to have negative mole fraction
        in the initial guess. Used for dynamically loaded library to help avoid doAbort. */
    for (i=0, rSum=0.0; i<NR; i++) rSum += r[i];
    /*coeff = 1.0;
    xSpecies[ 0] = 1.0 - rSum*coeff;*/                   /* SiO2  */
    if (rSum > 1.0) s[0] = rSum - 1.0 + sqrt(DBL_EPSILON);
#endif

    if (!rANDsTOx (r, s)) {
        /* [ 0] SiO2 [10] CaO [14] CO2 */
        static double tolerance = 0.0;

        if (tolerance == 0.0) tolerance = pow(DBL_EPSILON, (double) (2.0/3.0));

        /***********************************/
        /* Find minimum length Soln vector */
        /***********************************/

        /* objective fucntion row for finding minimum length solution vector that satisfies constraints */

        /***********************************/
        /* Find maximum length Soln vector */
        /***********************************/

        /* objective fucntion row for finding maximum length solution vector that satisfies constraints */

        /***************************************************************/
        /* Return average of two solutions - should always be feasible */
        /***************************************************************/

        for (i=0; i<NS; i++) {
            //s[i] = (sMin[i] + sMax[i])/2.0;
            if (fabs(s[i]) < tolerance) s[i] = 0.0;
        }

        if(!rANDsTOx (r, s)) {
#ifdef DEBUG
            fprintf(stderr, "Failed to find feasible solution in initialGuessOrdering.\n");
#endif
        }

        //for (i=NS; i<NT; i++) s[i] = 1.0/(((double) NY)+1.0);

#ifdef DEBUG
        printf("Results of call to initialGuessOrdering:\n");
        printf("   %20.20s %13.13s %13.13s %13.13s\n", "Species", "Mole frac", "r", "s");
        printf("   %20.20s %13.6g\n", liquid[0].label, xSpecies[0]);
        for (i=0;  i<NR; i++) printf("   %20.20s %13.6g %13.6g\n", liquid[i+1].label, xSpecies[i+1], r[i]);
        for (i=0;  i<NS; i++) printf("   %20.20s %13.6g %13.13s %13.6g\n", liquid[i+NA].label, xSpecies[i+NA], "", s[i]);
        for (i=NS; i<NT; i++) printf("   %20.20s %13.13s %13.13s %13.6g\n", "order CN[*]", "", "", s[i]);
#endif

        //return;
    } /* end block on simplex method */

    for (i=0; i<NS; i++) {
#ifndef TESTDYNAMICLIB
        s0 = MAX(rSum - 1.0, 0.0);
        sCorr[i] = s0; s[i] = sCorr[i] + sqrt(DBL_EPSILON);
        if (!rANDsTOx (r, s)) s[i] = sCorr[i];
        else {
            sCorr[i] = 0.5; s[i] = 1.0;
            while (sCorr[i] > sqrt(DBL_EPSILON)) {
                if(!rANDsTOx (r, s)) s[i] -= sCorr[i]; else s[i] += sCorr[i];
                sCorr[i] /= 2.0;
            }
        }
        sCorr[i] = (s[i] + s0)/2.0;
        s[i]     = 0.0;
    }

    for (i=0; i<NS; i++) s[i] = sCorr[i];
    for (i=0; i<NS; i++) sCorr[i] -= s0;
    while (factor > sqrt(DBL_EPSILON) && !rANDsTOx (r, s)) {
        factor /= 2.0;
        for (i=0; i<NS; i++) s[i] = s0 + sCorr[i]*factor;
    }

#else
        sCorr[i] = 0.0; s[i] = sqrt(DBL_EPSILON);
        if (!rANDsTOx (r, s)) s[i] = 0.0;
        else {
            sCorr[i] = 0.5; s[i] = 1.0;
            while (sCorr[i] > sqrt(DBL_EPSILON)) {
                if(!rANDsTOx (r, s)) s[i] -= sCorr[i]; else s[i] += sCorr[i];
                sCorr[i] /= 2.0;
            }
        }
        sCorr[i] = s[i]/2.0;
        s[i]     = 0.0;
    }

    for (i=0; i<NS; i++) s[i] = sCorr[i];
    while (factor > sqrt(DBL_EPSILON) && !rANDsTOx (r, s)) {
        factor /= 2.0;
        for (i=0; i<NS; i++) s[i] = sCorr[i]*factor;
    }

#endif
    //for (i=NS; i<NT; i++) s[i] = 1.0/(((double) NY)+1.0);

#ifdef DEBUG
    printf("Results of call to initialGuessOrdering:\n");
    printf("   %20.20s %13.13s %13.13s %13.13s\n", "Species", "Mole frac", "r", "s");
    printf("   %20.20s %13.6g\n", liquid[0].label, xSpecies[0]);
    for (i=0;  i<NR; i++) printf("   %20.20s %13.6g %13.6g\n", liquid[i+1].label, xSpecies[i+1], r[i]);
    for (i=0;  i<NS; i++) printf("   %20.20s %13.6g %13.13s %13.6g\n", liquid[i+NA].label, xSpecies[i+NA], "", s[i]);
    for (i=NS; i<NT; i++) printf("   %20.20s %13.13s %13.13s %13.6g\n", "order CN[*]", "", "", s[i]);
#endif
}

static void
order(int mask, double t, double p, double r[NR],
            double s[NT],            /* s[NT]                  BINARY MASK: 000000000001 */
            double dr[NT][NR] ,      /* ds[NT]/dr[NR]          BINARY MASK: 000000000010 */
            double dt[NT],           /* ds[NT]/dt              BINARY MASK: 000000000100 */
            double dp[NT],           /* ds[NT]/dp              BINARY MASK: 000000001000 */
            double dr2[NT][NR][NR],  /* d2s[NT]/dr[NR]dr[NR]   BINARY MASK: 000000010000 */
            double drt[NT][NR],      /* d2s[NT]/dr[NR]dt       BINARY MASK: 000000100000 */
            double drp[NT][NR],      /* d2s[NT]/dr[NR]dp       BINARY MASK: 000001000000 */
            double dt2[NT],          /* d2s[NT]/dt2            BINARY MASK: 000010000000 */
            double dtp[NT],          /* d2s[NT]/dtp            BINARY MASK: 000100000000 */
            double dp2[NT],          /* d2s[NT]/dp2            BINARY MASK: 001000000000 */
            double dw[NT][3*NP],     /* ds[NT]/dw[3*NP]        BINARY MASK: 010000000000 */
            double dtw[NT][3*NP]     /* ds[NT]/dtdw[3*NP]      BINARY MASK: 100000000000 */
     )
{
    double tOld         = getTOld();
    double pOld         = getPOld();
    double *rOld        = getROld();
    double *sOld        = getSOld();
    double **d2gds2     = getD2gds2();
    gsl_matrix       *ptToD2gds2 = getPtToD2gds2();
    gsl_matrix       *ptToV = getPtToV();
    gsl_vector       *ptToS = getPtToS();
    gsl_permutation  *indexLU = getIndexD2gds2();
    static ModelParameters *shadowParameters = NULL;
    static EosModelParameters *eosShadowParameters;
    int i, j, iter = 0, doAbort=FALSE, update=FALSE, loop, dLU;

    if (shadowParameters == NULL) {
        shadowParameters = (ModelParameters *) calloc((size_t) NP, sizeof(ModelParameters));
        eosShadowParameters = (EosModelParameters *) calloc((size_t) nc, sizeof(EosModelParameters));
    }

    /* default value for update is FALSE */
    for (i=0; i<NP; i++) {                             /* check if modelparameters are changing (independent of t, p) */
        if (modelParameters[i].activeH && (modelParameters[i].enthalpy != shadowParameters[i].enthalpy)) { update |= TRUE;
            shadowParameters[i].enthalpy = modelParameters[i].enthalpy; }
        if (modelParameters[i].activeS && (modelParameters[i].entropy  != shadowParameters[i].entropy )) { update |= TRUE;
            shadowParameters[i].entropy  = modelParameters[i].entropy;  }
        if (modelParameters[i].activeV && (modelParameters[i].volume   != shadowParameters[i].volume  )) { update |= TRUE;
            shadowParameters[i].volume   = modelParameters[i].volume ;  }
    }

    for (i=0; i<nc; i++) {                             /* check if modelparameters are changing (independent of t, p) */
        if (eosModelParameters[i].activeKp   && (eosModelParameters[i].Kp   != eosShadowParameters[i].Kp  )) { update |= TRUE; eosShadowParameters[i].Kp   = eosModelParameters[i].Kp  ; }
        if (eosModelParameters[i].activeKpp  && (eosModelParameters[i].Kpp  != eosShadowParameters[i].Kpp )) { update |= TRUE; eosShadowParameters[i].Kpp  = eosModelParameters[i].Kpp ; }
        if (eosModelParameters[i].activeKppp && (eosModelParameters[i].Kppp != eosShadowParameters[i].Kppp)) { update |= TRUE; eosShadowParameters[i].Kppp = eosModelParameters[i].Kppp; }
    }

    update |= (t != tOld);
    update |= (p != pOld);
    if (update) loadTaylorCoefficients(t, p);         /* if T, P, or meltsAndCO2_H2OModelParameters change                    */
    for (i=0; i<NR; i++) update |= (r[i] != rOld[i]);

    /* look-up or compute the current ordering state */
    if (update) {                                     /* if T, P, meltsAndCO2_H2OModelParameters or liquid composition change */
        double dgds[NT], sNew[NT], dgdsNORM=0.0, residual[NT];
        gsl_matrix *ptToD2gds2Copy = gsl_matrix_alloc((size_t) NT, (size_t) NT);
        gsl_vector_view vvToDgds = gsl_vector_view_array(dgds, (size_t) NT);
        for (i=0; i<NT; i++) { sNew[i] = sOld[i]; sOld[i] = 2.0; }
        convergedInOrder = TRUE;

        initialGuessOrdering(r, sNew);
        if (!rANDsTOx (r, sNew)) {
            printf("Initial guess to ordering state iteration is infeasible\n");
            doAbort = TRUE;
        }
        if (doAbort) {
            printf("Results of call to ordering with bad initial guess:\n");
            printf("   %20.20s %13.13s %13.13s %13.13s\n", "Species", "Mole frac", "r", "s");
            printf("   %20.20s %13.6g\n", liquid[0].label, xSpecies[0]);
            for (i=0;  i<NR; i++) printf("   %20.20s %13.6g %13.6g\n",         liquid[i+1].label,  xSpecies[i+1],  r[i]);
            for (i=0;  i<NS; i++) printf("   %20.20s %13.6g %13.13s %13.6g\n", liquid[i+NA].label, xSpecies[i+NA], "", sNew[i]);
            for (i=NS; i<NT; i++) printf("   %20.20s %13.13s %13.13s %13.6g\n", "order CN[*]", "", "", sNew[i]);
#ifdef USESEH
            doInterrupt = TRUE;
            raise_sigabrt(EXCEPTION_FLT_INVALID_OPERATION);
#else
            (void) raise(SIGABRT);
#endif
        }

        loop = TRUE;
        while (loop) {
            double s[NT], deltaS[NT], lambda=1.0, dgdsCopy[NT]; // dgdsCopy used as work space
            gsl_vector_view vvToDeltaS = gsl_vector_view_array(deltaS, (size_t) NT);
            gsl_vector_view vvToDgdsCopy = gsl_vector_view_array(dgdsCopy, (size_t) NT);
            int cycle = TRUE /*, errGaussj */;

            for (i=0; i<NT; i++) s[i] = sNew[i];

            fillDGDS   (r, s, t, p, dgds);
            fillD2GDS2New (r, s, t, p, d2gds2);

            //dgdsNORM=0.0
            for (i=0; i<NT; i++) {
                sOld[i] = s[i];
                for (j=0; j<NT; j++) if (s[i] == 0.0 || s[j] == 0.0) d2gds2[i][j] = 0.0;
                if (s[i] == 0.0) {
                    d2gds2[i][i] = 1.0;
                    dgds[i]      = 0.0;
                }
                //dgdsNORM += pow(dgds[i]/MAX(1.0, fabs(eosIntDGDS[i])), (double) 2.0);
            }
            //dgdsNORM = sqrt(dgdsNORM);

            gsl_matrix_scale(ptToD2gds2, -1.0);
            gsl_matrix_memcpy(ptToD2gds2Copy, ptToD2gds2);

            if (USE_SVD) {
                gsl_linalg_SV_decomp(ptToD2gds2, ptToV, ptToS, &vvToDgdsCopy.vector);
                gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToDgds.vector, &vvToDeltaS.vector);
            }
            else {
                melts_LU_decomp(ptToD2gds2, indexLU, &dLU);
                melts_LU_solve(ptToD2gds2, indexLU, &vvToDgds.vector, &vvToDeltaS.vector);
                melts_LU_refine(ptToD2gds2Copy, ptToD2gds2, indexLU,  &vvToDgds.vector, &vvToDeltaS.vector, &vvToDgdsCopy.vector);
            }
            for (i=0; i<NS; i++) s[i] += deltaS[i];


#ifdef DEBUG
            printf("--->dgds:   ");
            for (i=0; i<NT; i++) if (fabs(s[i]) > 10.0*DBL_EPSILON) printf("%20.13g", dgds[i]);
            printf("\n");
            for (j=0; j<NT; j++) if (fabs(s[j]) > 10.0*DBL_EPSILON) {
                printf("--->d2gds2[%d][]: ", j);
                for (i=0; i<NT; i++) if (fabs(s[i]) > 10.0*DBL_EPSILON) printf("%20.13g", d2gds2Copy[j][i]);
                printf("\n");
            }
            printf("--->s:      ");
            for (i=0; i<NT; i++) if (fabs(s[i]) > 10.0*DBL_EPSILON) printf("%20.13g", s[i]);
            printf("\n");
            printf("--->dels:   ");
            for (i=0; i<NT; i++) if (fabs(s[i]) > 10.0*DBL_EPSILON) printf("%20.13g", deltaS[i]);
            printf("\n");
#endif

            while (cycle && !rANDsTOx (r, s)) {
                lambda /= 2.0;
                for (j=0; j<NT; j++) s[j] = sOld[j] + lambda*deltaS[j];
                /* if (lambda < DBL_EPSILON) { */
                if (lambda < DBL_MIN) {
                    cycle = FALSE;
                    s[0] = (double) iter;
                    iter = MAX_ITER - 1;
                    fprintf(stderr, "\n*****lambda -> zero in ORDER. Terminating search loop.\n");
                  }
            }
#ifdef DEBUG
            printf("steplength correction:  = %20.13g\n", lambda);
            printf("--->s(adj): ");
            for (i=0; i<NT; i++) if (fabs(s[i]) > 10.0*DBL_EPSILON) printf("%20.13g", s[i]);
            printf("\n----- end -----\n");
#endif

            for (i=0; i<NT; i++) sNew[i] = s[i];
            iter++;
            loop = FALSE;
            if (iter < MAX_ITER) for (i=0; i<NT; i++) loop |= (fabs(sNew[i]-sOld[i]) > 10.0*DBL_EPSILON);
        }
        tOld = t;
        pOld = p;
        for (i=0; i<NR; i++) rOld[i] = r[i];

        (void) rANDsTOx (rOld, sOld);

        if (iter == MAX_ITER) {
            double sNorm;
            for (i=0, sNorm=0.0; i<NT; i++) sNorm += pow(sNew[i]-sOld[i], (double) 2.0);
            sNorm = sqrt(sNorm);
            if (sNorm > sqrt(DBL_EPSILON)) {
                /* convergedInOrder = FALSE; */
                fprintf(stderr, "ERROR in LIQUID.C (function ORDER). Failed to converge!\n");
                if (iter >= MAX_ITER) fprintf(stderr, " Iteration limit (%4d) exceeded.\n", iter);
                fprintf(stderr, "   T (C) = %8.2f, P (GPa) = %10.4f\n", t-273.15, p/10000.0);
                fprintf(stderr, "   %20.20s %13.13s %13.13s %13.13s %13.13s %13.13s\n", "Species", "Mole frac", "r", "s", "dgds", "deltaS");
                fprintf(stderr, "   %20.20s %13.6g\n", liquid[0].label, xSpecies[0]);
                for (i=0;  i<NR; i++) fprintf(stderr, "   %20.20s %13.6g %13.6g\n",         liquid[i+1].label,  xSpecies[i+1],  r[i]);
                for (i=0;  i<NS; i++) fprintf(stderr, "   %20.20s %13.6g %13.13s %13.6g %13.6g %13.6g\n", liquid[i+NA].label, xSpecies[i+NA], "", sOld[i], dgds[i], sNew[i]-sOld[i]);
                for (i=NS; i<NT; i++) fprintf(stderr, "   %20.20s %13.13s %13.13s %13.6g %13.6g %13.6g\n", "order CN[*]", "", "", sOld[i], dgds[i], sNew[i]-sOld[i]);
                fprintf(stderr, " sNorm             = %20.13g\n", sNorm);
                fprintf(stderr, " dgdsNorm          = %20.13g\n", dgdsNORM);
                fprintf(stderr, " 10*DBL_EPSILON    = %20.13g\n", 10.0*DBL_EPSILON);
                fprintf(stderr, " DBL_EPSILON^(2/3) = %20.13g\n", pow(DBL_EPSILON, 2.0/3.0));
                fprintf(stderr, " DBL_EPSILON^(1/2) = %20.13g\n", sqrt(DBL_EPSILON));
                //fprintf(stderr, " eosIntegralBranch = %s\n", (eosIntegralBranch == GMAPeosBRANCH) ? "GMAP" : "LMAP");
            } else if (sNorm > pow(DBL_EPSILON, 2.0/3.0)) {
                fprintf(stderr, "WARNING in LIQUID.C (function ORDER). sNorm = %g, dgdsNorm = %g [eps = %g, sqrt(eps) = %g]\n", sNorm, dgdsNORM, DBL_EPSILON, sqrt(DBL_EPSILON));
            }
        }

#ifdef DEBUG
        printf("Results of ordering state calculation:\n");
        printf("   T (C) = %8.2f, P (GPa) = %10.4f\n", t-273.15, p/10000.0);
        printf("   %20.20s %13.13s %13.13s %13.13s %13.13s %13.13s %13.13s\n", "Species", "Mole frac", "r", "s", "dgds", "deltaS", "eosIntDGDS");
        printf("   %20.20s %13.6g\n", liquid[0].label, xSpecies[0]);
        for (i=0;  i<NR; i++) printf("   %20.20s %13.6g %13.6g\n", liquid[i+1].label, xSpecies[i+1], r[i]);
        for (i=0;  i<NS; i++) printf("   %20.20s %13.6g %13.13s %13.6g %13.6g %13.6g %13.6g\n", liquid[i+NA].label, xSpecies[i+NA], "", sOld[i], dgds[i], sNew[i]-sOld[i], eosIntDGDS[i]);
        for (i=NS; i<NT; i++) printf("   %20.20s %13.13s %13.13s %13.6g %13.6g %13.6g\n", "order CN[*]", "", "", sOld[i], dgds[i], sNew[i]-sOld[i]);
        printf(" 10*DBL_EPSILON    = %20.13g\n", 10.0*DBL_EPSILON);
        printf(" DBL_EPSILON^(2/3) = %20.13g\n", pow(DBL_EPSILON, 2.0/3.0));
        printf(" DBL_EPSILON^(1/2) = %20.13g\n", sqrt(DBL_EPSILON));
        printf(" eosIntegralBranch = %s\n", (eosIntegralBranch == GMAPeosBRANCH) ? "GMAP" : "LMAP");
        for (i=0; i<NS; i++) {
            double s[NT];
            for (j=0; j<NT; j++) s[j] = sOld[j];
            s[i] += 30.0*DBL_EPSILON;
            fillDGDS(r, s, t, p, dgds);
            printf(" s[%d]+30eps gives dgds = %20.13g\n", i, dgds[i]);
            s[i] -= 10.0*DBL_EPSILON;
            fillDGDS(r, s, t, p, dgds);
            printf(" s[%d]+20eps gives dgds = %20.13g\n", i, dgds[i]);
            s[i] -= 10.0*DBL_EPSILON;
            fillDGDS(r, s, t, p, dgds);
            printf(" s[%d]+10eps gives dgds = %20.13g\n", i, dgds[i]);
            printf(" ----------------------\n");
            s[i] = sOld[i];
            fillDGDS(r, s, t, p, dgds);
            printf(" s[%d]       gives dgds = %20.13g\n", i, dgds[i]);
            printf(" ----------------------\n");
            s[i] -= 10.0*DBL_EPSILON;
            fillDGDS(r, s, t, p, dgds);
            printf(" s[%d]-10eps gives dgds = %20.13g\n", i, dgds[i]);
            s[i] -= 10.0*DBL_EPSILON;
            fillDGDS(r, s, t, p, dgds);
            printf(" s[%d]-20eps gives dgds = %20.13g\n", i, dgds[i]);
            s[i] -= 10.0*DBL_EPSILON;
            fillDGDS(r, s, t, p, dgds);
            printf(" s[%d]-30eps gives dgds = %20.13g\n", i, dgds[i]);
            printf(" ++++++++++++++++++++++\n");
        }
#endif

        gsl_matrix_free(ptToD2gds2Copy);

        setTOld(tOld);
        setPOld(pOld);
        /* arrays (rOld, sOld, d2gds2, etc) should be preserved automatically */

    }

    if (mask & FIRST  ) {   /* return s        */
        for (i=0; i<NT; i++) s[i] = sOld[i];
    }
    if (mask & SECOND ) {   /* compute ds/dr:  */
        double *s = sOld;
        double d2gdrds[NR][NT];
        gsl_matrix_view mvToDr = gsl_matrix_view_array((double *) dr, (size_t) NT, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NT);

        fillD2GDRDS (r, s, t, p,  d2gdrds);

        for (j=0; j<NR; j++) {
            gsl_vector_view vvToDr = gsl_matrix_column(&mvToDr.matrix, j),
                vvToD2gdrds = gsl_matrix_row(&mvToD2gdrds.matrix, j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdrds.vector, &vvToDr.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdrds.vector, &vvToDr.vector);
            for (i=0; i<NT; i++) dr[i][j] = (s[i] > 0.0) ? dr[i][j] : 0.0;
        }
    }

    if (mask & THIRD  ) {   /* compute ds/dt:  */
        double *s = sOld;
        double d2gdsdt[NT];
        gsl_vector_view vvToDt = gsl_vector_view_array(dt, (size_t) NT),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NT);

        fillD2GDSDT (r, s, t, p, d2gdsdt);

        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdt.vector, &vvToDt.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdt.vector, &vvToDt.vector);
        for (i=0; i<NT; i++) dt[i] = (s[i] > 0.0) ? dt[i] : 0.0;

    }
    if (mask & FOURTH ) {   /* compute ds/dp:  */
        double *s = sOld;
        double d2gdsdp[NT];
        gsl_vector_view vvToDp = gsl_vector_view_array(dp, (size_t) NT),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NT);

        fillD2GDSDP (r, s, t, p, d2gdsdp);

        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdp.vector, &vvToDp.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdp.vector, &vvToDp.vector);
        for (i=0; i<NT; i++) dp[i] = (s[i] > 0.0) ? dp[i] : 0.0;
    }
    if (mask & FIFTH  ) {   /* compute d2s/dr2 */
        double *s = sOld;
        double d2gdrds[NR][NT];
        double d3gdr2ds[NR][NR][NT];
        double d3gdrds2[NR][NT][NT];
        double d3gds3[NT][NT][NT];
        double dsdr[NT][NR], temp[NT], temp2[NT];
        gsl_matrix_view mvToDsdr = gsl_matrix_view_array((double *) dsdr, (size_t) NT, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NT);
        gsl_vector_view vvToTemp = gsl_vector_view_array(temp, (size_t) NT), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NT);
        int k, l, m, n;

        fillD2GDRDS  (r, s, t, p,  d2gdrds);
        fillD3GDR2DS (r, s, t, p,  d3gdr2ds);
        fillD3GDRDS2 (r, s, t, p,  d3gdrds2);
        fillD3GDS3   (r, s, t, p,  d3gds3);

        /* compute dsdr matrix */
        for (j=0; j<NR; j++) {
            gsl_vector_view vvToDsdr = gsl_matrix_column(&mvToDsdr.matrix, j),
              vvToD2gdrds = gsl_matrix_row(&mvToD2gdrds.matrix, j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdrds.vector, &vvToDsdr.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdrds.vector, &vvToDsdr.vector);
            for (i=0; i<NT; i++) dsdr[i][j] = (s[i] > 0.0) ? dsdr[i][j] : 0.0;
        }

        /* compute dsdr2 cube */
        for (j=0; j<NR; j++) {
            for (k=0; k<NR; k++) {
                for (l=0; l<NT; l++) {
                    temp[l] = d3gdr2ds[j][k][l];
                    for (m=0; m<NT; m++) {
                        temp[l] += d3gdrds2[j][l][m]*dsdr[m][k]
                            + d3gdrds2[k][l][m]*dsdr[m][j];
                        for (n=0; n<NT; n++)
                            temp[l] += d3gds3[l][m][n]*dsdr[m][j]*dsdr[n][k];
                    }
                }
                if (USE_SVD) {
                    gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
                    for (i=0; i<NT; i++) dr2[i][j][k] = (s[i] > 0.0) ? temp2[i] : 0.0;
                }
                else {
                    melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
                    for (i=0; i<NT; i++) dr2[i][j][k] = (s[i] > 0.0) ? temp[i] : 0.0;
                }
            }
        }
    }

    if (mask & SIXTH  ) {   /* compute d2s/drt */
        double *s = sOld;
        double d2gdrds[NR][NT];
        double d2gdsdt[NT];
        double d3gdrds2[NR][NT][NT];
        double d3gdrdsdt[NR][NT];
        double d3gds2dt[NT][NT];
        double d3gds3[NT][NT][NT];
        double dsdr[NT][NR], dsdt[NT], temp[NT], temp2[NT];
        gsl_matrix_view mvToDsdr = gsl_matrix_view_array((double *) dsdr, (size_t) NT, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NT);
        gsl_vector_view vvToDsdt = gsl_vector_view_array(dsdt, (size_t) NT),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NT),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NT), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NT);
        int k, l, m;

        fillD2GDRDS   (r, s, t, p,  d2gdrds);
        fillD2GDSDT   (r, s, t, p,  d2gdsdt);
        fillD3GDRDS2  (r, s, t, p,  d3gdrds2);
        fillD3GDRDSDT (r, s, t, p,  d3gdrdsdt);
        fillD3GDS3    (r, s, t, p,  d3gds3);
        fillD3GDS2DT  (r, s, t, p,  d3gds2dt);

        /* compute dsdr matrix */
        for (j=0; j<NR; j++) {
            gsl_vector_view vvToDsdr = gsl_matrix_column(&mvToDsdr.matrix, j),
              vvToD2gdrds = gsl_matrix_row(&mvToD2gdrds.matrix, j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdrds.vector, &vvToDsdr.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdrds.vector, &vvToDsdr.vector);
            for (i=0; i<NT; i++) dsdr[i][j] = (s[i] > 0.0) ? dsdr[i][j] : 0.0;
        }

        /* compute dsdt vector */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        for (i=0; i<NT; i++) dsdt[i] = (s[i] > 0.0) ? dsdt[i] : 0.0;

        /* compute dsdrdt matrix */
        for (j=0; j<NR; j++) {
            for (k=0; k<NT; k++) {
                temp[k] = d3gdrdsdt[j][k];
                for (l=0; l<NT; l++) {
                    temp[k] += d3gdrds2[j][k][l]*dsdt[l] + d3gds2dt[k][l]*dsdr[l][j];
                    for (m=0; m<NT; m++) temp[k] += d3gds3[k][l][m]*dsdr[l][j]*dsdt[m];
                }
            }
            if (USE_SVD) {
                gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
                for (i=0; i<NT; i++) drt[i][j] = (s[i] > 0.0) ? temp2[i] : 0.0;
            }
            else {
                melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
                for (i=0; i<NT; i++) drt[i][j] = (s[i] > 0.0) ? temp[i] : 0.0;
            }
        }

    }
    if (mask & SEVENTH) {   /* compute d2s/drp */
        double *s = sOld;
        double d2gdrds[NR][NT];
        double d2gdsdp[NT];
        double d3gdrds2[NR][NT][NT];
        double d3gdrdsdp[NR][NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dp[NT][NT];
        double dsdr[NT][NR], dsdp[NT], temp[NT], temp2[NT];
        gsl_matrix_view mvToDsdr = gsl_matrix_view_array((double *) dsdr, (size_t) NT, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NT);
        gsl_vector_view vvToDsdp = gsl_vector_view_array(dsdp, (size_t) NT),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NT),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NT), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NT);

        int k, l, m;

        fillD2GDRDS   (r, s, t, p,  d2gdrds);
        fillD2GDSDP   (r, s, t, p,  d2gdsdp);
        fillD3GDRDS2  (r, s, t, p,  d3gdrds2);
        fillD3GDRDSDP (r, s, t, p,  d3gdrdsdp);
        fillD3GDS3    (r, s, t, p,  d3gds3);
        fillD3GDS2DP  (r, s, t, p,  d3gds2dp);

        /* compute dsdr matrix */
        for (j=0; j<NR; j++) {
            gsl_vector_view vvToDsdr = gsl_matrix_column(&mvToDsdr.matrix, j),
                vvToD2gdrds = gsl_matrix_row(&mvToD2gdrds.matrix, j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdrds.vector, &vvToDsdr.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdrds.vector, &vvToDsdr.vector);
            for (i=0; i<NT; i++) dsdr[i][j] = (s[i] > 0.0) ? dsdr[i][j] : 0.0;
        }

        /* compute dsdp vector */
        if (USE_SVD)  gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdp.vector, &vvToDsdp.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdp.vector, &vvToDsdp.vector);
        for (i=0; i<NT; i++) dsdp[i] = (s[i] > 0.0) ? dsdp[i] : 0.0;

        /* compute dsdrdp matrix */
        for (j=0; j<NR; j++) {
            for (k=0; k<NT; k++) {
                temp[k] = d3gdrdsdp[j][k];
                for (l=0; l<NT; l++) {
                        temp[k] += d3gdrds2[j][k][l]*dsdp[l] + d3gds2dp[k][l]*dsdr[l][j];
                        for (m=0; m<NT; m++) temp[k] += d3gds3[k][l][m]*dsdr[l][j]*dsdp[m];
                }
            }
            if (USE_SVD) {
                gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
                for (i=0; i<NT; i++) drp[i][j] = (s[i] > 0.0) ? temp2[i] : 0.0;
            }
            else {
                melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
                for (i=0; i<NT; i++) drp[i][j] = (s[i] > 0.0) ? temp[i] : 0.0;
            }
        }
    }

    if (mask & EIGHTH ) {   /* compute d2s/dt2 */
        double *s = sOld;
        double d2gdsdt[NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dt[NT][NT];
        double d3gdsdt2[NT];
        double dsdt[NT], temp[NT], temp2[NT];
        gsl_vector_view vvToDsdt = gsl_vector_view_array(dsdt, (size_t) NT),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NT),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NT), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NT);
        int k, l;

        fillD2GDSDT  (r, s, t, p,  d2gdsdt);
        fillD3GDS3   (r, s, t, p,  d3gds3);
        fillD3GDS2DT (r, s, t, p,  d3gds2dt);
        fillD3GDSDT2 (r, s, t, p,  d3gdsdt2);

        /* compute dsdt vector */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        for (i=0; i<NT; i++) dsdt[i] = (s[i] > 0.0) ? dsdt[i] : 0.0;

        /* compute dsdt2 vector */
        for (j=0; j<NT; j++) {
            temp[j] = d3gdsdt2[j];
            for (k=0; k<NT; k++) {
                temp[j] +=  2.0*d3gds2dt[j][k]*dsdt[k];
                for (l=0; l<NT; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdt[l];
            }
        }
        if (USE_SVD) {
            gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
            for (i=0; i<NT; i++) dt2[i] = (s[i] > 0.0) ? temp2[i] : 0.0;
        }
        else {
            melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
            for (i=0; i<NT; i++) dt2[i] = (s[i] > 0.0) ? temp[i] : 0.0;
        }
    }

    if (mask & NINTH  ) {   /* compute d2s/dtp */
        double *s = sOld;
        double d2gdsdt[NT];
        double d2gdsdp[NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dt[NT][NT];
        double d3gds2dp[NT][NT];
        double d3gdsdtdp[NT];
        double dsdt[NT], dsdp[NT], temp[NT], temp2[NT];
        gsl_vector_view vvToDsdt = gsl_vector_view_array(dsdt, (size_t) NT),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NT),
            vvToDsdp = gsl_vector_view_array(dsdp, (size_t) NT),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NT),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NT), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NT);
        int k, l;

        fillD2GDSDT   (r, s, t, p,  d2gdsdt);
        fillD2GDSDP   (r, s, t, p,  d2gdsdp);
        fillD3GDS3    (r, s, t, p,  d3gds3);
        fillD3GDS2DT  (r, s, t, p,  d3gds2dt);
        fillD3GDS2DP  (r, s, t, p,  d3gds2dp);
        fillD3GDSDTDP (r, s, t, p,  d3gdsdtdp);

        /* compute dsdt vector */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        for (i=0; i<NT; i++) dsdt[i] = (s[i] > 0.0) ? dsdt[i] : 0.0;

        /* compute dsdp vector */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdp.vector, &vvToDsdp.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdp.vector, &vvToDsdp.vector);
        for (i=0; i<NT; i++) dsdp[i] = (s[i] > 0.0) ? dsdp[i] : 0.0;

        /* compute dsdtp vector */
        for (j=0; j<NT; j++) {
            temp[j] = d3gdsdtdp[j];
            for (k=0; k<NT; k++) {
                temp[j] += d3gds2dt[j][k]*dsdp[k] + d3gds2dp[j][k]*dsdt[k];
                for (l=0; l<NT; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdp[l];
            }
        }
        if (USE_SVD) {
            gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
            for (i=0; i<NT; i++) dtp[i] = (s[i] > 0.0) ? temp2[i] : 0.0;
        }
        else {
            melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
            for (i=0; i<NT; i++) dtp[i] = (s[i] > 0.0) ? temp[i] : 0.0;
        }
    }

    if (mask & TENTH  ) {   /* compute d2s/dp2 */
        double *s = sOld;
        double d2gdsdp[NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dp[NT][NT];
        double d3gdsdp2[NT];
        double dsdp[NT], temp[NT], temp2[NT];
        gsl_vector_view vvToDsdp = gsl_vector_view_array(dsdp, (size_t) NT),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NT),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NT), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NT);
        int k, l;

        fillD2GDSDP  (r, s, t, p,  d2gdsdp);
        fillD3GDS3   (r, s, t, p,  d3gds3);
        fillD3GDS2DP (r, s, t, p,  d3gds2dp);
        fillD3GDSDP2 (r, s, t, p,  d3gdsdp2);

        /* compute dsdp vector */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdp.vector, &vvToDsdp.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdp.vector, &vvToDsdp.vector);
        for (i=0; i<NT; i++) dsdp[i] = (s[i] > 0.0) ? dsdp[i] : 0.0;

        /* compute dsdp2 vector */
        for (j=0; j<NT; j++) {
            temp[j] = d3gdsdp2[j];
            for (k=0; k<NT; k++) {
                temp[j] +=  2.0*d3gds2dp[j][k]*dsdp[k];
                for (l=0; l<NT; l++) temp[j] += d3gds3[j][k][l]*dsdp[k]*dsdp[l];
            }
        }
        if (USE_SVD) {
            gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
            for (i=0; i<NT; i++) dp2[i] = (s[i] > 0.0) ? temp2[i] : 0.0;
        }
        else {
            melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
            for (i=0; i<NT; i++) dp2[i] = (s[i] > 0.0) ? temp[i] : 0.0;
        }
    }

    if (mask & ELEVENTH ) {   /* compute ds/dw:  */
        double *s = sOld;
        double d2gdsdw[NT][3*NP];
        gsl_matrix_view mvToDw = gsl_matrix_view_array((double *) dw, (size_t) NT, (size_t) 3*NP),
            mvToD2gdsdw = gsl_matrix_view_array((double *) d2gdsdw, (size_t) NT, (size_t) 3*NP);

        fillD2GDSDW (r, s, t, p,  d2gdsdw);
        for (j=0; j<NP; j++) {
            gsl_vector_view vvToDw, vvToD2gdsdw;
            vvToDw = gsl_matrix_column(&mvToDw.matrix,      j);
            vvToD2gdsdw = gsl_matrix_column(&mvToD2gdsdw.matrix,      j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdw.vector, &vvToDw.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdw.vector, &vvToDw.vector);
            for (i=0; i<NT; i++) dw[i][     j] = (s[i] > 0.0) ? dw[i][     j] : 0.0;
            vvToDw = gsl_matrix_column(&mvToDw.matrix,   NP+j);
            vvToD2gdsdw = gsl_matrix_column(&mvToD2gdsdw.matrix,   NP+j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdw.vector, &vvToDw.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdw.vector, &vvToDw.vector);
            for (i=0; i<NT; i++) dw[i][  NP+j] = (s[i] > 0.0) ? dw[i][  NP+j] : 0.0;
            vvToDw = gsl_matrix_column(&mvToDw.matrix, 2*NP+j);
            vvToD2gdsdw = gsl_matrix_column(&mvToD2gdsdw.matrix, 2*NP+j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdw.vector, &vvToDw.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdw.vector, &vvToDw.vector);
            for (i=0; i<NT; i++) dw[i][2*NP+j] = (s[i] > 0.0) ? dw[i][2*NP+j] : 0.0;
        }
    }

    if (mask & TWELFTH ) {   /* compute ds/dtdw:  */
        double *s = sOld;
        double d2gdsdt[NT];
        double d2gdsdw[NT][3*NP];
        double d3gds3[NT][NT][NT];
        double d3gds2dt[NT][NT];
        double d3gds2dw[NT][NT][3*NP];
        double d3gdsdtdw[NT][3*NP];
        double dsdt[NT], dsdw[NT][3*NP], temp[NT], temp2[NT];
        gsl_matrix_view mvToDsdw = gsl_matrix_view_array((double *) dsdw, (size_t) NT, (size_t) 3*NP),
            mvToD2gdsdw = gsl_matrix_view_array((double *) d2gdsdw, (size_t) NT, (size_t) 3*NP);
        gsl_vector_view vvToDsdt = gsl_vector_view_array(dsdt, (size_t) NT),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NT),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NT), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NT);

        fillD2GDSDW (r, s, t, p,  d2gdsdw);
        for (j=0; j<NP; j++) {
            gsl_vector_view vvToDsdw, vvToD2gdsdw;
            vvToDsdw = gsl_matrix_column(&mvToDsdw.matrix,      j);
            vvToD2gdsdw = gsl_matrix_column(&mvToD2gdsdw.matrix,      j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdw.vector, &vvToDsdw.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdw.vector, &vvToDsdw.vector);
            for (i=0; i<NT; i++) dsdw[i][     j] = (s[i] > 0.0) ? dsdw[i][     j] : 0.0;
            vvToDsdw = gsl_matrix_column(&mvToDsdw.matrix,   NP+j);
            vvToD2gdsdw = gsl_matrix_column(&mvToD2gdsdw.matrix,   NP+j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdw.vector, &vvToDsdw.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdw.vector, &vvToDsdw.vector);
            for (i=0; i<NT; i++) dsdw[i][  NP+j] = (s[i] > 0.0) ? dsdw[i][  NP+j] : 0.0;
            vvToDsdw = gsl_matrix_column(&mvToDsdw.matrix, 2*NP+j);
            vvToD2gdsdw = gsl_matrix_column(&mvToD2gdsdw.matrix, 2*NP+j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdw.vector, &vvToDsdw.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdw.vector, &vvToDsdw.vector);
            for (i=0; i<NT; i++) dsdw[i][2*NP+j] = (s[i] > 0.0) ? dsdw[i][2*NP+j] : 0.0;
        }

        fillD2GDSDT (r, s, t, p, d2gdsdt);
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        for (i=0; i<NT; i++) dsdt[i] = (s[i] > 0.0) ? dsdt[i] : 0.0;

        fillD3GDSDTDW (r, s, t, p,  d3gdsdtdw);
        fillD3GDS2DW  (r, s, t, p,  d3gds2dw);
        fillD3GDS2DT  (r, s, t, p,  d3gds2dt);
        fillD3GDS3    (r, s, t, p,  d3gds3);
        for (j=0; j<NP; j++) {
            for (i=0; i<NT; i++) {
                int k, l;
                temp[i] = d3gdsdtdw[i][j];
                for (k=0; k<NT; k++) {
                    if (dsdw[k][j] != 0.0) temp[i] += d3gds2dt[i][k]*dsdw[k][j];
                    if (dsdt[k] != 0.0) {
                        temp[i] += d3gds2dw[i][k][j]*dsdt[k];
                        for (l=0; l<NT; l++) if (dsdw[l][j] != 0.0) temp[i] += d3gds3[i][k][l]*dsdt[k]*dsdw[l][j];
                    }
                }
            }
            if (USE_SVD) {
                gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
                for (i=0; i<NT; i++) dtw[i][j] = (s[i] > 0.0) ? temp2[i] : 0.0;
            }
            else {
                melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
                for (i=0; i<NT; i++) dtw[i][j] = (s[i] > 0.0) ? temp[i] : 0.0;
            }
            for (i=0; i<NT; i++) {
                int k, l;
                temp[i] = d3gdsdtdw[i][j+NP];
                for (k=0; k<NT; k++) {
                    if (dsdw[k][j+NP] != 0.0) temp[i] += d3gds2dt[i][k]*dsdw[k][j+NP];
                    if (dsdt[k] != 0.0) {
                        temp[i] += d3gds2dw[i][k][j+NP]*dsdt[k];
                        for (l=0; l<NT; l++) if (dsdw[l][j+NP] != 0.0) temp[i] += d3gds3[i][k][l]*dsdt[k]*dsdw[l][j+NP];
                    }
                }
            }
            if (USE_SVD) {
                gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
                for (i=0; i<NT; i++) dtw[i][j+NP] = (s[i] > 0.0) ? temp2[i] : 0.0;
            }
            else {
                melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
                for (i=0; i<NT; i++) dtw[i][j+NP] = (s[i] > 0.0) ? temp[i] : 0.0;
            }
            for (i=0; i<NT; i++) {
                int k, l;
                temp[i] = d3gdsdtdw[i][j+2*NP];
                for (k=0; k<NT; k++) {
                    if (dsdw[k][j+2*NP] != 0.0) temp[i] += d3gds2dt[i][k]*dsdw[k][j+2*NP];
                    if (dsdt[k] != 0.0) {
                        temp[i] += d3gds2dw[i][k][j+2*NP]*dsdt[k];
                        for (l=0; l<NT; l++) if (dsdw[l][j+2*NP] != 0.0) temp[i] += d3gds3[i][k][l]*dsdt[k]*dsdw[l][j+2*NP];
                    }
                }
            }
            if (USE_SVD) {
                gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
                for (i=0; i<NT; i++) dtw[i][j+2*NP] = (s[i] > 0.0) ? temp2[i] : 0.0;
            }
            else {
                melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
                for (i=0; i<NT; i++) dtw[i][j+2*NP] = (s[i] > 0.0) ? temp[i] : 0.0;
            }
        }

    }

}

/*
 *=============================================================================
 * Public functions:
 *    t       -  Temperature (K)
 *    p       -  Pressure (bars)
 *    *r      -  (pointer to x[]) Array of independent compositional variables
 */

double /* returns 1 + a (p-pr) + b (p-pr)^2 */
retLiqEosParam (double r[NR], double t, double p,
    double *v1Ret, /* dv/dp   at reference temperature and 1 bar */
    double *v2Ret, /* d2v/dp2 at reference temperature and 1 bar */
    double *v3Ret, /* d3v/dp3 at reference temperature and 1 bar */
    double *v4Ret, /* d4v/dp4 at reference temperature and 1 bar */
    double *aRet,  /* a = f(v1,v2,v3,v4)                         */
    double *bRet)  /* b = f(v1,v2,v3,v4)                         */
{
    const double pr      = 1.0;
    const double tr      = 1673.15;
    double m[NA], mOx[NA+1], mOxTot, v, dvdt, cRef, mw, cp, d2vdp2, d3vdp3, d4vdp4, v1Ref, alpha, a, b;
    double s[NT];
    int i, j;
    double result = 0.0;

    liqERRstate = ERR_NONE;
    MTHREAD_ONCE(&initThreadBlock, threadInit);

    /* Convert input composition (r) to liquid moles (m)  */
    for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* Compute moles and total moles of oxides */
    for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
    if (mOxTot == 0.0) return result;

    order(FIRST, t, p, r, s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    /* Deal with the special case of FeO1.3 */
    /* Deal with the special case of FeO1.3 */
    mOx[NA] = 0.0;
    if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
        const double y = 0.3;
        mOx[iOxFeO1_3] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
        }
    }

    for (i=0, v=0.0, dvdt=0.0, cRef=0.0, mw=0.0, cp=0.0, d2vdp2=0.0, d3vdp3=0.0, d4vdp4=0.0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
        v        += mOx[i]*bulkSystem[i].gk_v;
        dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
        cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
        cp      += mOx[i]*bulkSystem[i].gk_cp;
        d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
        d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
        d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        mw      += mOx[i]*bulkSystem[i].mw;
    }
    if (v == 0.0) return result;

    alpha   = dvdt/v;
    v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
    a      = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;
    b      = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;

    *v1Ret = v1Ref;
    *v2Ret = d2vdp2;
    *v3Ret = d3vdp3;
    *v4Ret = d4vdp4;
    *aRet  = a;
    *bRet  = b;
    return (1.0 + a*(p-pr) + b*(p-pr)*(p-pr));
}

typedef struct _eosRefParameters {
    double v, dvdt, cRef, dcdt, cp, d2vdp2, d3vdp3, d4vdp4, mw;
    double alpha, v0, v1, v2, v3, v4, a, b;
    double K, Kp, Kpp, Kppp;
} EosRefParameters;

EosRefParameters *getEosRefParameters(double *r) {
    static EosRefParameters eos;
    double s[NT];
    const double tr = 1673.15;
    const double pr = 1.0;
    double m[NA], mOx[NA+1], mOxTot, denom;
    int i, j;

    liqERRstate = ERR_NONE;
    MTHREAD_ONCE(&initThreadBlock, threadInit);

    eos.v      = 0.0; eos.dvdt   = 0.0; eos.cRef   = 0.0; eos.dcdt = 0.0; eos.cp = 0.0;
    eos.d2vdp2 = 0.0; eos.d3vdp3 = 0.0; eos.d4vdp4 = 0.0; eos.mw   = 0.0;
    eos.alpha  = 0.0; eos.v0     = 0.0; eos.v1     = 0.0; eos.v2   = 0.0; eos.a  = 0.0; eos.b = 0.0;
    eos.K      = 0.0; eos.Kp     = 0.0; eos.Kpp    = 0.0; eos.Kppp = 0.0;

    /* Convert input composition (r) to liquid moles (m)  */
    for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* Compute moles and total moles of oxides */
    for (i=0, mOxTot=0.0; i<nc; i++) { for (j=0, mOx[i]=0.0; j<NA; j++) mOx[i] += m[j]*(liquid[j].liqToOx)[i]; mOxTot += mOx[i]; }
    if (mOxTot < 100.0*DBL_EPSILON) return &eos;

    order(FIRST, tr, pr, r, s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    /* Deal with the special case of FeO1.3 */
    mOx[NA] = 0.0;
    if ((iOxFe2O3 != -1) && (iOxFeO != -1) && (iOxFeO1_3 != -1)) {
        const double y = 0.3;
        mOx[iOxFeO1_3] = 0.0;
        if (iCmpFe2SiO4_6 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2SiO4_6]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2SiO4_6]*nSpecies;
        }
        if (iCmpFe2AlO4_1 != -1) {
            mOx[iOxFeO1_3] += 2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFe2O3]  -= y*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOx[iOxFeO]    -= (1.0-2.0*y)*2.0*s[iCmpFe2AlO4_1]*nSpecies;
            mOxTot         += 2.0*y*s[iCmpFe2AlO4_1]*nSpecies;
        }
    }

    for (i=0; i<((iOxFeO1_3 != -1) ? nc+1 : nc); i++) {
        eos.v    += mOx[i]*bulkSystem[i].gk_v;
        eos.dvdt    += mOx[i]*bulkSystem[i].gk_dvdt;
        eos.cRef    += mOx[i]*bulkSystem[i].gk_c/mOxTot + ((iOxAl2O3 != -1) ? mOx[i]*mOx[iOxAl2O3]*bulkSystem[i].gk_cXal2o3/(mOxTot*mOxTot) : 0.0);
        eos.dcdt    += mOx[i]*bulkSystem[i].gk_dcdt/mOxTot;
        eos.cp      += mOx[i]*bulkSystem[i].gk_cp;
        eos.d2vdp2  += mOx[i]*(bulkSystem[i].gk_d2vdp2 + eosModelParameters[i].v2);
        eos.d3vdp3  += mOx[i]*(bulkSystem[i].gk_d3vdp3 + eosModelParameters[i].v3);
        eos.d4vdp4  += mOx[i]*(bulkSystem[i].gk_d4vdp4 + eosModelParameters[i].v4);
        eos.mw      += mOx[i]*bulkSystem[i].mw;
    }
    if (eos.v < 100.0*DBL_EPSILON) return &eos;

    eos.alpha = eos.dvdt/eos.v;
    eos.v0    = eos.v;
    eos.v1    = -eos.v*eos.v*(1000.0/(eos.mw*eos.cRef*eos.cRef) + tr*eos.alpha*eos.alpha/eos.cp);
    eos.v2    = eos.d2vdp2;
    eos.v3    = eos.d3vdp3;
    eos.v4    = eos.d4vdp4;
    denom     = 2.0*eos.v1*eos.d3vdp3-3.0*eos.d2vdp2*eos.d2vdp2;
    eos.a        = (denom) ? (eos.d2vdp2*eos.d3vdp3 - eos.v1*eos.d4vdp4/2.0)        /denom : 0.0;
    eos.b        = (denom) ? (eos.d2vdp2*eos.d4vdp4/4.0 - eos.d3vdp3*eos.d3vdp3/3.0)/denom : 0.0;

    eos.K    = -eos.v0/eos.v1;
    eos.Kp   = eos.v2*eos.K*eos.K/eos.v0 - 1.0;
    eos.Kpp  = (eos.v3*eos.K*eos.K*eos.K/eos.v0 + (2.0*eos.Kp+1.0)*(eos.Kp+1.0))/eos.K;
    eos.Kppp = (eos.v4*eos.K*eos.K*eos.K*eos.K/eos.v0 + eos.Kpp*eos.K*(4.0+6.0*eos.Kp)
                            - (3.0*eos.Kp+1.0)*(2.0*eos.Kp+1.0)*(eos.Kp+1.0))
            /(eos.K*eos.K);

    return &eos;
}

/* function utilized by conLiq and test_eos for high-P ferric/ferrous calculation */

double integralV_GKsp(int index, double t, double p) {
    const double pr = 1.0;
    const double tr = 1673.15;
    double v, dvdt, cRef, c, dcdt, mw, cp, d2vdp2, d3vdp3, d4vdp4, v0, v1, v1Ref, v2, alpha, a, b, sum;
    double gInt;

    if (fabs(p-pr) < 10.0*DBL_EPSILON) return (double) 0.0;

    v      = bulkSystem[index].gk_v;
    dvdt    = bulkSystem[index].gk_dvdt;
    cRef    = bulkSystem[index].gk_c;
    dcdt    = bulkSystem[index].gk_dcdt;
    cp      = bulkSystem[index].gk_cp;
    d2vdp2  = bulkSystem[index].gk_d2vdp2 + eosModelParameters[index].v2;
    d3vdp3  = bulkSystem[index].gk_d3vdp3 + eosModelParameters[index].v3;
    d4vdp4  = bulkSystem[index].gk_d4vdp4 + eosModelParameters[index].v4;
    mw      = bulkSystem[index].mw;

    if (v == 0.0) return 0.0;

    alpha   = dvdt/v;
    v0      = v*exp(alpha*(t-tr));
    v1Ref   = -v*v*(1000.0/(mw*cRef*cRef) + tr*alpha*alpha/(cp));
    c      = cRef + (t-tr)*dcdt;
    v1      = -v0*v0*(1000.0/(mw*c*c) + t*alpha*alpha/(cp));
    v2      = d2vdp2;
    a      = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d3vdp3 - v1Ref*d4vdp4/2.0)     /(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;
    b      = (2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2 != 0.0) ? (d2vdp2*d4vdp4/4.0 - d3vdp3*d3vdp3/3.0)/(2.0*v1Ref*d3vdp3-3.0*d2vdp2*d2vdp2) : 0.0;
    sum      = a*a - 4.0*b;
    gInt    = 0.0;

    if ((a == 0.0) && (b == 0.0)) {
        gInt      += v0*(p-pr) + v1*(p-pr)*(p-pr)/2.0 + v2*(p-pr)*(p-pr)*(p-pr)/6.0;

    } else if ((a != 0.0) && (b == 0.0)) {
        gInt      += (v0 - v2/(2.0*a*a))*(p-pr) + (v1 + v2/(2.0*a))*(p-pr)*(p-pr)/2.0 + v2*log(1.0+a*(p-pr))/(2.0*a*a*a);

    } else if ((a == 0.0) && (b != 0.0)) {
        gInt      += (v0 + v2/(2.0*b))*(p-pr) + v1*log(1.0 + b*(p-pr)*(p-pr))/(2.0*b);
        gInt      += (b > 0.0) ? -v2*atan(sqrt(b)*(p-pr))/(2.0*b*sqrt(b)) : -v2*log((1.0+sqrt(-b)*(p-pr))/(1.0-sqrt(-b)*(p-pr)))/(4.0*b*sqrt(-b));

    } else if (sum > 0.0) {
        double x = sqrt(sum);
        double y = (a + 2.0*b*(p-pr))/x;
        double z = a/x;
        double PcA = (2.0*b*pr - a + x)/(2.0*b);
        double PcB = (2.0*b*pr - a - x)/(2.0*b);
        double arg = (1.0 + y)*(1.0 - z)/((1.0 - y)*(1.0 + z));

        if (((pr < PcA) && (PcA < p)) || ((pr < PcB) && (PcB < p))) fprintf(stderr, "index %d, v %g,alpha %g, v1Ref %g, v1 %g\nc %g, v2 %g, a %g, b %g\n", index, v, alpha, v1Ref, v1, c, v2, a, b);
        if (((pr < PcA) && (PcA < p)) || ((pr < PcB) && (PcB < p))) printERR("integralV_GKsp", ERR_SUM_GT_ZERO, "1.0+a*(p-pr)+b*(p-pr)*(p-pr)", 1.0+a*(p-pr)+b*(p-pr)*(p-pr))
        if (arg <= 0.0)                        printERR("integralV_GKsp", ERR_SUM_GT_ZERO, "arg", arg)

        gInt      += (v0 + a*v1/b + v2/(2.0*b))*(p-pr);
        gInt      += (v1*(1.0-a*a/b) - v2*a/(2.0*b))*log(1.0+a*(p-pr)+b*(p-pr)*(p-pr))/(2.0*b);
        gInt      += (v1*a*(3.0-a*a/b) + v2*(1.0 - a*a/(2.0*b)))*log(arg)/(2.0*b*x);

    } else if (sum == 0.0) {
        gInt      += (v0 + 4.0*v1/a + 2.0*v2/(a*a))*(p-pr);
        gInt      += -8.0*(v2/a + v1)/(a*a*(2.0+a*(p-pr))) + 4.0*(v2/a + v1)/(a*a);
        gInt      += -4.0*(3.0*v1 + 2.0*v2/a)*log(1.0 + a*(p-pr)/2.0)/(a*a);

    } else if(sum < 0.0) {
        double x = sqrt(-sum);
        double y = (a + 2.0*b*(p-pr))/x;
        double z = a/x;

        gInt      += (v0 + a*v1/b + v2/(2.0*b))*(p-pr);
        gInt      += (v1*(1.0-a*a/b) - v2*a/(2.0*b))*log(1.0+a*(p-pr)+b*(p-pr)*(p-pr))/(2.0*b);
        gInt      += (v1*a*(3.0-a*a/b) + v2*(1.0 - a*a/(2.0*b)))*atan((z-y)/(1.0+z*y))/(b*x);

    }

    return gInt;
}

/*
 *=============================================================================
 * Public functions:
 *    inpMask -  bitwise mask for specifying input parameters
 *    outMask -  bitwise mask for selecting output
 *    mask    -  bitwise mask for selecting output
 *    t       -  Temperature (K)
 *    p       -  Pressure (bars)
 *    *r      -  (pointer to x[]) Array of independent compositional variables
 */

void
conLiq(int inpMask, int outMask, double t, double p,
    double *o,      /* comp of liquid in moles of oxides                        */
    double *m,      /* comp of liquid in moles of endmember components          */
    double *r,      /* comp of liquid in terms of the independent comp var      */
    double *x,      /* comp of liquid in mole fractions of endmember components */
    double **dm,    /* Jacobian matrix: dm[i][j] = dr[i]/dm[j]                  */
    double ***d2m,  /* vector of matrices: d2m[i][j][k] = d2r[i]/dm[j]dm[k]     */
    double *logfo2) /* base 10 logarithm of the oxygen fugacity                 */
{
    /*---------------------------------------------------------------------------
    Not all combinations of inpMask and outMask are feasible. Valid
        combinations are:

       inpMask          outMask
    (1)  FIRST | SEVENTH  FIRST
    (2)  FIRST            SEVENTH
    (3)  SECOND           THIRD | FOURTH | FIFTH | SIXTH
    (4)  THIRD            FOURTH
    (5)  THIRD            FOURTH | EIGHTH
    (6)  THIRD            FOURTH | NINTH

    (1) converts a vector of moles of oxides into a vector of moles of oxides
            with the correct redox state for the given t, p, and logfo2. Note that
            the original vector is used as output.
    (2) calculates from a vector of moles of oxides and the given t and p, the
            appropriate logfo2
    (3) calculates from a vector of moles of endmember components, one or
            all of: r[], x[], dr[]/dm[], or d2r[]/dm[]dm[]
    (4) calculates from a vector of independent compositional variables
            mole fractions of endmember components
    (5) calculates from a vector of independent compositional variables
            mole fractions of endmember species
    (6) calculates from a vector of independent compositional variables
            fraction of system in CN state (vector 0 ... NT-1)
    ----------------------------------------------------------------------------*/

    int i, j, k;

    liqERRstate = ERR_NONE;

#if defined(PHMELTS_ADJUSTMENTS) && ! defined(TRUE_xMELTS)
    if ((silminState->fo2Alt) && (calculationMode != MODE_xMELTS)) {
        //      conLiq_Alt(inpMask, outMask, t, p, o, m, r, x, dm, d2m, logfo2);
            return;
    } else
#endif
    if ((calculationMode == MODE__MELTS) || (calculationMode == MODE_pMELTS) ) {
            conLiq_v34(inpMask, outMask, t, p, o, m, r, x, dm, d2m, logfo2);
            return;
    } else if (calculationMode == MODE__MELTSandCO2) {
            conLiq_CO2(inpMask, outMask, t, p, o, m, r, x, dm, d2m, logfo2);
            return;
    } else if (calculationMode == MODE__MELTSandCO2_H2O) {
            conLiq_CO2_H2O(inpMask, outMask, t, p, o, m, r, x, dm, d2m, logfo2);
            return;
    }

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    if (inpMask & FIRST) {
        /*-------------------------------------------------------------------------
            Oxide+logfo2 -> oxide or oxide -> logfo2. The algorithm used is that
            given by:
            Kress, VC and Carmichael, ISE (1991) The compressibility of silicate
                liquids containing Fe2O3 and the effect of composition, temperature,
                oxygen fugacity and pressure on their redox states.
                Contributions to Mineralogy and Petrology (in press)
            Coefficients for the oxides are initialized in LIQ_STRUCT_DATA.H
        --------------------------------------------------------------------------*/
        static int indexFeO = 0, indexFe2O3 = 0;
#ifdef USE_KRESS_CARMICHAEL_FO2
        static const double t0 = 1673.15,                      /* K       */
                         a =    0.196,
                         b =    1.1492e4,                  /* K       */
                         c =   -6.675,
                         e =   -3.364,
                         f =   -7.01e-7  * 1.0e5,          /* K/bar   */
                         g =   -1.54e-10 * 1.0e5,          /* 1/bar   */
                         h =    3.85e-17 * 1.0e5 * 1.0e5;  /* K/bar^2 */
#endif /* USE_KRESS_CARMICHAEL_FO2 */
        double sum = 0.0, temp;

        if (indexFeO == 0 || indexFe2O3 == 0) {
            for (i=0; i<nc; i++) {
                if (bulkSystem[i].type == FEO)   indexFeO   = i;
                if (bulkSystem[i].type == FE2O3) indexFe2O3 = i;
            }
            if (indexFeO == 0 || indexFe2O3 == 0) {
         printf("Fatal error in conLiq (LIQUID.C)\n");
         printf("The oxides FeO and Fe2O3 cannot be identified.\n");
         return;
            }
        }

        if (inpMask == (FIRST | SEVENTH)  && outMask == FIRST) {
            /*----------------------------------------------------------------------
                Converts a vector of moles of oxides (as defined in LIQ_STRUCT_DATA.H
                for the structure bulkSystem) into a vector of moles of oxides with
                the correct redox state (ferric/ferrous ratio) for the given bulk
                composition, t and p.
            ------------------------------------------------------------------------*/

            o[indexFeO]   += 2.0*o[indexFe2O3];
            o[indexFe2O3]  = 0.0;
            if (o[indexFeO] == 0.0) return;

            for (i=0; i<nc; i++) sum += o[i];
            if (sum == 0.0) return;

#ifdef USE_KRESS_CARMICHAEL_FO2
            if (p < 50000.0) temp = a*log(10.0)*(*logfo2) + b/t + c + e*(1.0 - t0/t - log(t/t0))
                                                        + f*p/t + g*(t-t0)*p/t + h*SQUARE(p)/t;
            else             temp = a*log(10.0)*(*logfo2) + b/t + c + e*(1.0 - t0/t - log(t/t0))
                                                        + f*50000.0/t + g*(t-t0)*50000.0/t + h*SQUARE(50000.0)/t
                - a*log(10.0)*(608.966*p/10000.0-608.966*5.0)/t;

            for (i=0; i<nc; i++) temp += bulkSystem[i].coeff*o[i]/sum;
            temp = exp(temp);

            o[indexFe2O3]  = temp*o[indexFeO]/(1.0 + 2.0*temp);
            o[indexFeO]   -= 2.0*o[indexFe2O3];
#else
            {
                double y = 0.3;
                double intVFe2O3  = integralV_GKsp(iOxFe2O3,  t, p);
                double intVFeO    = integralV_GKsp(iOxFeO,    t, p);
                double intVFeO1_3 = integralV_GKsp(iOxFeO1_3, t, p);
                double deltaG = -106200.0 - t*(-55.1) + 31.86*(t - 1673.15 - t*log(t/1673.15)) + intVFe2O3/2.0 - intVFeO;
                double KD1 = exp(-deltaG/(R*t)
                            -(39860.0*o[iOxAl2O3] - 62520.0*o[iOxCaO] - 102000.0*o[iOxNa2O] - 119000.0*o[iOxK2O])/(R*t*sum));
                double K2 = 0.4*exp(-(intVFeO1_3 - (1.0-2.0*y)*intVFeO - y*intVFe2O3)/(R*t));
                double fo2 = exp((*logfo2)*log(10.0));

                temp =  (KD1*pow(fo2, (double) 0.25) + 2.0*y*K2*pow(KD1, 2.0*y)*pow(fo2, y/2.0))
                /(1.0 + (1.0-2.0*y)*K2*pow(KD1, 2.0*y)*pow(fo2, y/2.0));

                o[indexFe2O3]  = o[indexFeO]*(1.0-1.0/(1.0+temp))/2.0;
                o[indexFeO]   -= 2.0*o[indexFe2O3];
            }
#endif /* USE_KRESS_CARMICHAEL_FO2 */

        } else if (inpMask == FIRST && outMask == SEVENTH) {
            /*----------------------------------------------------------------------
                Calculates from the given t and p and a vector of moles of oxides
                (as defined in LIQ_STRUCT_DATA.H for the structure bulkSystem) the
                appropriate log10fo2 for the given t and p.
            ------------------------------------------------------------------------*/

            if (o[indexFeO] == 0.0 || o[indexFe2O3] == 0.0) { *logfo2 = 0.0; return; }
            for (i=0; i<nc; i++) sum += o[i];
            sum += o[indexFe2O3];
            if (sum == 0.0) { *logfo2 = 0.0; return; }

#ifdef USE_KRESS_CARMICHAEL_FO2
            if (p< 50000.0) temp = b/t + c + e*(1.0 - t0/t - log(t/t0)) + f*p/t + g*(t-t0)*p/t + h*SQUARE(p)/t;
            else            temp = b/t + c + e*(1.0 - t0/t - log(t/t0)) + f*50000.0/t + g*(t-t0)*50000.0/t + h*SQUARE(50000.0)/t
                             - a*log(10.0)*(608.966*p/10000.0-608.966*5.0)/t;
            for (i=0; i<nc; i++) temp += bulkSystem[i].coeff*o[i]/sum;
            temp += 2.0*bulkSystem[indexFeO].coeff*o[indexFe2O3]/sum
                                - bulkSystem[indexFe2O3].coeff*o[indexFe2O3]/sum;
            *logfo2 = (log(o[indexFe2O3]/o[indexFeO]) - temp)/(a*log(10.0));
#else
            {
                double y = 0.3;
                double intVFe2O3  = integralV_GKsp(iOxFe2O3,  t, p);
                double intVFeO    = integralV_GKsp(iOxFeO,    t, p);
                double intVFeO1_3 = integralV_GKsp(iOxFeO1_3, t, p);
                double deltaG = -106200.0 - t*(-55.1) + 31.86*(t - 1673.15 - t*log(t/1673.15)) + intVFe2O3/2.0 - intVFeO;
                double KD1 = exp(-deltaG/(R*t)
                            -(39860.0*o[iOxAl2O3] - 62520.0*o[iOxCaO] - 102000.0*o[iOxNa2O] - 119000.0*o[iOxK2O])/(R*t*sum));
                double K2 = 0.4*exp(-(intVFeO1_3 - (1.0-2.0*y)*intVFeO - y*intVFe2O3)/(R*t));
                int converged = FALSE;
                int iter = 0;

                *logfo2 = -10.0;

                while (!converged && (iter < 200)) {
                double fo2    = exp((*logfo2)*log(10.0));
                double dfo2   = exp((*logfo2)*log(10.0))*log(10.0);
                double numer  = KD1*pow(fo2, (double) 0.25) + 2.0*y*K2*pow(KD1, 2.0*y)*pow(fo2, y/2.0);
                double dnumer = 0.25*KD1*dfo2/pow(fo2, (double) 0.75) + 2.0*y*K2*pow(KD1, 2.0*y)*(y/2.0)*dfo2/pow(fo2, 1.0-y/2.0);
                double denom  = 1.0 + (1.0-2.0*y)*K2*pow(KD1, 2.0*y)*pow(fo2, y/2.0);
                double ddenom = (1.0-2.0*y)*K2*pow(KD1, 2.0*y)*(y/2.0)*dfo2/pow(fo2, 1.0-y/2.0);

                double f  =  numer/denom - 2.0*o[indexFe2O3]/o[indexFeO];
                double df = dnumer/denom - numer*ddenom/(denom*denom);
                double corr = -f/df;

                if (fabs(corr) > sqrt(DBL_EPSILON)) *logfo2 += corr; else converged = TRUE;
                if (*logfo2 >  10.0) *logfo2 =  10.0;
                if (*logfo2 < -50.0) *logfo2 = -50.0;
                iter++;
            }

    if (!converged) printf("Convergence failure in conLiq_CO2_H2O(FIRST,SEVENTH)\n");
            }
#endif /* USE_KRESS_CARMICHAEL_FO2 */

        } else
            printf("Illegal call to conLiq with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

    } else if (inpMask == SECOND) {
        double sum;

        if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH))
            printf("Illegal call to conLiq with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        for (i=0, sum=0.0; i<NA; i++) sum += m[i];

        if (outMask & THIRD) {
            /* Converts a vector of moles of end-member components (m) into a vector
         of independent compositional variables (r) required as input for the
         remaining public functions.
         The dependent variable is taken to be SiO2 (1st component), as this
         component will never have a mole fraction of zero.                   */

            for (i=0; i<NR; i++) r[i] = (sum != 0.0) ? m[i+1]/sum : 0.0;
        }

        if (outMask & FOURTH) {
            /* Converts a vector of moles of end-member components (m) into a vector
         of mole fractions of endmember components                            */

            for (i=0; i<NA; i++) x[i] = (sum != 0.0) ? m[i]/sum : 0.0;
        }

        if (outMask & FIFTH) {
            /* Calculates the matrix dr[i]/dm[j] using m[] as input                 */

            if (sum == 0.0) {
                for (i=0; i<NR; i++) for (j=0; j<NA; j++) dm[i][j] = 0.0;
            } else {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++) {
                        dm[i][j] = (i+1 == j) ? (1.0-m[i+1]/sum)/sum : - m[i+1]/SQUARE(sum);
                    }
                }
            }
        }

        if (outMask & SIXTH) {
            /* Calculates the matrix d2r[i]/dm[j]dm[k] using m[] as input           */

            if (sum == 0.0) {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++)  {
                        for (k=0; k<NA; k++) d2m[i][j][k] = 0.0;
                    }
                }
            } else {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++) {
                        for (k=0; k<NA; k++) {
                            d2m[i][j][k]  = 2.0*m[i+1]/CUBE(sum);
                            d2m[i][j][k] -= (i+1 == j) ? 1.0/SQUARE(sum) : 0.0;
                            d2m[i][j][k] -= (i+1 == k) ? 1.0/SQUARE(sum) : 0.0;
                        }
                    }
                }
            }

        }

    } else if (inpMask == THIRD && outMask == FOURTH) {
        /* Converts a vector of independent compositional variables (r)
            into a vector of mole fractions of end-member components (x)            */

        for (i=0, x[0] = 1.0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }

    } else if (inpMask == THIRD && outMask == (FOURTH | EIGHTH)) {
        /* Converts a vector of independent compositional variables (r)
            into a vector of mole fractions of endmember species (x)                */

        MTHREAD_MUTEX_LOCK(&global_data_mutex);
        order(0, t, p, r,
                NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
        for (i=0; i<NE; i++) x[i] = xSpecies[i];
        MTHREAD_MUTEX_UNLOCK(&global_data_mutex);

    } else if (inpMask == THIRD && outMask == (FOURTH | NINTH)) {
        double s[NT];
        /* Converts a vector of independent compositional variables (r)
            into a vector of fractions of system in CN state * (x)                */

        MTHREAD_MUTEX_LOCK(&global_data_mutex);
        order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
        for (i=0; i<NY; i++) x[i] = s[NS+i];
        MTHREAD_MUTEX_UNLOCK(&global_data_mutex);

    } else {
        printf("Illegal call to conLiq with inpMask = %o and outMask = %o\n",
            inpMask, outMask);
    }

}

int
testLiq(int mask, double t, double p,
    int na,          /* Expected number of endmember components                 */
    int nr,          /* Expected number of independent compositional variables  */
    char **names,    /* array of strings of names of endmember oxides           */
    char **formulas, /* array of strings of formulas of endmember components    */
    double *r,       /* array of indepependent compos variables, check bounds   */
    double *m)       /* array of moles of endmember components, check bounds    */
{
    const char *phase = "liquid.c";
    int result = TRUE, i;

    if ((calculationMode == MODE__MELTS) || (calculationMode == MODE_pMELTS)) return testLiq_v34(mask, t, p, na, nr, names, formulas, r, m);
        else if (calculationMode == MODE__MELTSandCO2) return testLiq_CO2(mask, t, p, na, nr, names, formulas, r, m);
        else if (calculationMode == MODE__MELTSandCO2_H2O) return testLiq_CO2_H2O(mask, t, p, na, nr, names, formulas, r, m);

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    if (mask & FIRST) {
        result = result && (na == NA);
        if (!result) printf("<<%s>> Wrong number of components!\n", phase);
    }
    if (mask & SECOND) {
        result = result && (nr == NR);
        if (!result) printf("<<%s>> Wrong number of indep variables!\n", phase);
    }
    if (mask & THIRD) {
        for (i=0; i<NA; i++) {
            result = result && (strcmp(names[i],bulkSystem[i].label) == 0);
            if (!result)
                printf("<<%s>> Oxide[%d] should be %s not %s.\n",
                    phase, i, bulkSystem[i].label, names[i]);
        }
    }
    if (mask & FOURTH) {
        for (i=0; i<NA; i++) {
            result = result && (strcmp(formulas[i],liquid[i].label) == 0);
            if (!result)
                printf("<<%s>> Component[%d] should have formula %s not %s.\n",
                    phase, i, liquid[i].label, formulas[i]);
        }
    }
    /* Check bounds on the independent compositional variables */
    if (mask & FIFTH) {
        double s[NT];
        MTHREAD_MUTEX_LOCK(&global_data_mutex);
        initialGuessOrdering(r, s);
        result = rANDsTOx (r, s);
        MTHREAD_MUTEX_UNLOCK(&global_data_mutex);
    }
    /* Check bounds on moles of endmember components */
    if (mask & SIXTH) {
        double rTemp[NR], s[NT], sum;
        for (i=0, sum=0.0; i<NA; i++) sum += m[i];
        for (i=0; i<NR; i++) rTemp[i] = (sum != 0.0) ? m[i+1]/sum : 0.0;
        if (sum > 0.0) {
            MTHREAD_MUTEX_LOCK(&global_data_mutex);
            initialGuessOrdering(rTemp, s);
            result = rANDsTOx (rTemp, s);
            MTHREAD_MUTEX_UNLOCK(&global_data_mutex);
        } else {
            result = FALSE;
        }
    }
    /* Check if ordering state calculation has converged */
    if (mask & SEVENTH) { /* This call is NOT thread safe, but it should never be called from a threaded app */
        result = convergedInOrder;
    }
    /* Check if EOS calculation is valid */
    if (mask & EIGHTH) { /* This call is NOT thread safe, but it should never be called from a threaded app */
        result = ((liqERRstate == ERR_NONE) ? TRUE : FALSE);
    }
    /* Check if number of CN states is corrcet */
    if (mask & NINTH) {
        result = result && (NY == nCN);
        if (!result) printf("<<%s>> Wrong number of CN states!\n", phase);
    }

    return result;
}

void
dispLiq(int mask, double t, double p, double *x,
    char **formula            /* Mineral formula for interface display MASK: 1 */
    )
{
    double *r = x;

    if ((calculationMode == MODE__MELTS) || (calculationMode == MODE_pMELTS)) {
            dispLiq_v34(mask, t, p, x, formula);
            return;
    } else if (calculationMode == MODE__MELTSandCO2) {
            dispLiq_CO2(mask, t, p, x, formula);
            return;
    } else if (calculationMode == MODE__MELTSandCO2_H2O) {
            dispLiq_CO2_H2O(mask, t, p, x, formula);
            return;
    }

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    if (mask & FIRST) {    /* assume maximum string length is 5 */
        char *string = (char *) malloc((unsigned) (7+NA*12+1)*sizeof(char));;
        double m[NA], oxVal[NA], oxSum;
        int i, j, n;

        (void) snprintf(string, 8, "wt%% ox:");

        for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }
        for (i=0, oxSum=0.0; i<NA; i++) {
            for (j=0, oxVal[i]=0.0; j<NA; j++) oxVal[i] += m[j]*(liquid[j].liqToOx)[i];
            oxVal[i] *= bulkSystem[i].mw;
            oxSum    += oxVal[i];
        }

        if (oxSum != 0.0) for (i=0, n=7; i<NA; i++)
            if (oxVal[i] != 0.0) {
                double w = 100.0*oxVal[i]/oxSum;
                int nn = snprintf(&string[n], 13, " %s %.2f", bulkSystem[i].label, w);
                n += (nn < 13) ? nn : 12;
            }

        *formula = string;
    }
}

static int returnMixingProperties = TRUE;
void setModeToMixingLiq(int flag) { returnMixingProperties = flag; }

void
actLiq(int mask, double t, double p, double *x,
    double *a,   /* (pointer to a[]) activities              BINARY MASK: 000001 */
    double *mu,  /* (pointer to mu[]) chemical potentials    BINARY MASK: 000010 */
    double **dx, /* (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 000100 */
    double **dw  /* (pointer to dw[][]) d(mu[])/d(w[])       BINARY MASK: 001000 */
    )
{
    double *r = x;
    double s[NT], g, dgdr[NR];
    double fr[NA][NR];
    int i, j;

    liqERRstate = ERR_NONE;

    if ((calculationMode == MODE__MELTS) || (calculationMode == MODE_pMELTS)) {
            actLiq_v34(mask, t, p, x, a, mu, dx);
            return;
    } else if (calculationMode == MODE__MELTSandCO2) {
            actLiq_CO2(mask, t, p, x, a, mu, dx);
            return;
    } else if (calculationMode == MODE__MELTSandCO2_H2O) {
            actLiq_CO2_H2O(mask, t, p, x, a, mu, dx);
            return;
    }

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    MTHREAD_MUTEX_LOCK(&global_data_mutex);
    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    for(i=0; i<NA; i++) for (j=0; j<NR; j++) fr[i][j] = rsEndmembers[i][j] - r[j];

    g = fillG(r, s, t, p);
    fillDGDR (r, s, t, p, dgdr);

    if (!mask && a != NULL) {
        for (i=0; i<NT; i++) a[i] = s[i];
    }

    if (mask & FIRST) {
        for(i=0; i<NA; i++) {
       for (a[i]=g, j=0; j<NR; j++) a[i] += fr[i][j]*dgdr[j];
       if (returnMixingProperties) a[i] -= G(i);
       a[i] = exp(a[i]/(R*t));
        }
    }

    if (mask & SECOND) {
        for(i=0; i<NA; i++) {
            for (mu[i]=g, j=0; j<NR; j++) mu[i] += fr[i][j]*dgdr[j];
            if (returnMixingProperties) mu[i] -= G(i);
        }

        g = gT;
        for (j=0; j<NR; j++) g += fr[nCO2][j]*dgdrT[j];
        // printf("X H2O, XCO2 = %g %g, muTernary CO2 %g %g\n", r[13], r[17], g, mu[nCO2]);
    }

    if (mask & THIRD) {
        double d2gdr2[NR][NR];
        double d2gdrds[NR][NT];
        double d2gds2[NT][NT];

        double dsdr[NT][NR], dfrdr[NA][NR], gs[NA][NT], aref[NA], sum;
        int k, l;

        fillD2GDR2  (r, s, t, p,  d2gdr2);
        fillD2GDRDS (r, s, t, p,  d2gdrds);
        fillD2GDS2  (r, s, t, p,  d2gds2);

        /* fill Darken structures */
        for(i=0; i<NA; i++) {
            for (j=0;  j<NR; j++) dfrdr[i][j] = -1.0;
            for (j=0;  j<NS; j++) gs[i][j] = rsEndmembers[i][NR+j] - s[j];
            for (j=NS; j<NT; j++) gs[i][j] = 0.0;
        }

        order(SECOND, t, p, r,
                    NULL, dsdr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        /* get reference activities */
        if (!(mask & FIRST)) {
            for(i=0; i<NA; i++) {
                for (aref[i]=g, j=0; j<NR; j++) aref[i] += fr[i][j]*dgdr[j];
                if (returnMixingProperties) aref[i] -= G(i);
                aref[i] = exp(aref[i]/(R*t));
            }
        } else for (i=0; i<NA; i++) aref[i] = a[i];

        /* Compute derivatives of the chemical potentials */
        for (i=0; i<NA; i++) {
            for (k=0; k<NR; k++) {
                sum = (1.0+dfrdr[i][k])*dgdr[k];
                for (j=0; j<NR; j++) {
                    sum += fr[i][j]*d2gdr2[j][k];
                    for (l=0; l<NT; l++) sum += fr[i][j]*d2gdrds[j][l]*dsdr[l][k];
                }
                for (j=0; j<NT; j++) {
                    sum += gs[i][j]*d2gdrds[k][j];
                    for (l=0; l<NT; l++) sum += gs[i][j]*d2gds2[j][l]*dsdr[l][k];
                }
                dx[i][k] = sum; /* This is d mu/d r */
            }
        }

        /* convert result to d a/d r */
        for (i=0; i<NA; i++) for (j=0; j<NR; j++) dx[i][j] *= aref[i]/(R*t);

    }

    if (mask & FOURTH) {
        double dgdw[3*NP];
        double d2gdrds[NR][NT];
        double d2gdrdw[NR][3*NP];
        double dsdw[NT][3*NP], sum;
        int k, l;

        fillDGDW    (r, s, t, p,  dgdw);
        fillD2GDRDS (r, s, t, p,  d2gdrds);
        fillD2GDRDW (r, s, t, p,  d2gdrdw);

        order(ELEVENTH, t, p, r,
                    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, dsdw, NULL);

        for (i=0; i<NA; i++) {
            for (j=0; j<NP; j++) {
                if (modelParameters[j].activeH) {
                    for (k=0, dw[i][j]=dgdw[j]; k<NR; k++) {
                        for (l=0, sum=d2gdrdw[k][j]; l<NT; l++) sum += d2gdrds[k][l]*dsdw[l][j];
                        dw[i][j] += fr[i][k]*sum;
                    }
                } else dw[i][j] = 0.0;

                if (modelParameters[j].activeS) {
                    for (k=0, dw[i][NP+j]=dgdw[NP+j]; k<NR; k++) {
                        for (l=0, sum=d2gdrdw[k][NP+j]; l<NT; l++) sum += d2gdrds[k][l]*dsdw[l][NP+j];
                        dw[i][NP+j] += fr[i][k]*sum;
                    }
                } else dw[i][NP+j] = 0.0;

                if (modelParameters[j].activeV) {
                    for (k=0, dw[i][2*NP+j]=dgdw[2*NP+j]; k<NR; k++) {
                        for (l=0, sum=d2gdrdw[k][2*NP+j]; l<NT; l++) sum += d2gdrds[k][l]*dsdw[l][2*NP+j];
                        dw[i][2*NP+j] += fr[i][k]*sum;
                    }
                } else dw[i][2*NP+j] = 0.0;
            }

            if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
                if (modelParameters[NW+i].activeH) dw[i][     NW+i] -= 1.0;
                if (modelParameters[NW+i].activeS) dw[i][  NP+NW+i] -= -t;
#ifndef USE_GHIORSO_KRESS_MODEL
                if (modelParameters[NW+i].activeV) dw[i][2*NP+NW+i] -= (p-1.0);
#endif /* USE_GHIORSO_KRESS_MODEL */
            }
        }
    }

    /************************************************************************
   * This return is not-public and is used for testing of the derivatives *
   * of G with respect to w[] internal to test_liquid.                    *
   ************************************************************************/
    if (mask & FIFTH) {
        fillDGDW (r, s, t, p, mu);

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            mu[     NW+0] -= 1.0;
            mu[  NP+NW+0] -= -t;
#ifndef USE_GHIORSO_KRESS_MODEL
            mu[2*NP+NW+0] -= (p-1.0);
#endif /* USE_GHIORSO_KRESS_MODEL */
            for (i=0; i<NR; i++) {
                mu[     NW  +0] +=         r[i];
                mu[  NP+NW  +0] +=      -t*r[i];
#ifndef USE_GHIORSO_KRESS_MODEL
                mu[2*NP+NW  +0] += (p-1.0)*r[i];
#endif /* USE_GHIORSO_KRESS_MODEL */
                mu[     NW+i+1] -=         r[i];
                mu[  NP+NW+i+1] -=      -t*r[i];
#ifndef USE_GHIORSO_KRESS_MODEL
                mu[2*NP+NW+i+1] -= (p-1.0)*r[i];
#endif /* USE_GHIORSO_KRESS_MODEL */
            }
        }

        for (i=0; i<(3*NP); i++) if (fabs(mu[i]) < 10.0*DBL_EPSILON) mu[i] = 0.0;
    }

    /*****************************************************************************
   * This return is used by the preclb_support.c functions to return the       *
   * configurational part of the activity.  It is now obsolete and is included *
   * for backwards compatibility.                                              *
   * ***************************************************************************/
    if (mask & SIXTH) {
        for (i=0; i<NA; i++) {
            a[i] = x[i];
            a[i] = (i != NA-1)   ? (1.0 - x[NA-1])*a[i] : x[NA-1]*a[NA-1];
        }
    }

    /*****************************************************************************
   * This return is not-public and is used within preclb_slave for             *
   * computation of the parameter derivatives of dgdr which are returned in dw.*
   *****************************************************************************/
    if (mask & SEVENTH) {
        double d2gdrds[NR][NT];
        double d2gdrdw[NR][3*NP];
        double dsdw[NT][3*NP];
        int k;

        fillD2GDRDS (r, s, t, p,  d2gdrds);
        fillD2GDRDW (r, s, t, p,  d2gdrdw);

        order(ELEVENTH, t, p, r,
                    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, dsdw, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NP; j++) {
                if (modelParameters[j].activeH) {
                    for (k=0, dw[i][j]=d2gdrdw[i][j]; k<NT; k++) dw[i][j] += d2gdrds[i][k]*dsdw[k][j];
                } else dw[i][j] = 0.0;
                if (modelParameters[j].activeS) {
                    for (k=0, dw[i][NP+j]=d2gdrdw[i][NP+j]; k<NT; k++) dw[i][NP+j] += d2gdrds[i][k]*dsdw[k][NP+j];
                } else dw[i][NP+j] = 0.0;
                if (modelParameters[j].activeV) {
                    for (k=0, dw[i][2*NP+j]=d2gdrdw[i][2*NP+j]; k<NT; k++) dw[i][2*NP+j] += d2gdrds[i][k]*dsdw[k][2*NP+j];
                } else dw[i][2*NP+j] = 0.0;
            }

            if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties dr[i] += (G(0)-G(i+1)); */
                if (modelParameters[NW  +0].activeH) dw[i][     NW  +0] += 1.0;
                if (modelParameters[NW  +0].activeS) dw[i][  NP+NW  +0] += -t;
#ifndef USE_GHIORSO_KRESS_MODEL
                if (modelParameters[NW  +0].activeV) dw[i][2*NP+NW  +0] += (p-1.0);
#endif /* USE_GHIORSO_KRESS_MODEL */
                if (modelParameters[NW+i+1].activeH) dw[i][     NW+i+1] -= 1.0;
                if (modelParameters[NW+i+1].activeS) dw[i][  NP+NW+i+1] -= -t;
#ifndef USE_GHIORSO_KRESS_MODEL
                if (modelParameters[NW+i+1].activeV) dw[i][2*NP+NW+i+1] -= (p-1.0);
#endif /* USE_GHIORSO_KRESS_MODEL */
            }
        }
    }

    /*******************************************************************************
    * This return is non-public and is designed for xMELTS calibration purposes.   *
    * The negative of the partial molar entropy (i.e. d mu dT) is returned in mu[] *
    ********************************************************************************/
    if (mask & EIGHTH) {
        double d2gdrds[NR][NT];
        double d2gdrdt[NR];
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double dsdt[NT], gs[NA][NT], sum, dgdt;
        int k;

        dgdt = -fillS(r, s, t, p);
        fillD2GDRDS (r, s, t, p,  d2gdrds);
        fillD2GDRDT (r, s, t, p,  d2gdrdt);
        fillD2GDS2  (r, s, t, p,  d2gds2);
        fillD2GDSDT (r, s, t, p,  d2gdsdt);

        /* fill Darken structures */
        for(i=0; i<NA; i++) {
            for (j=0;  j<NS; j++) gs[i][j] = rsEndmembers[i][NR+j] - s[j];
            for (j=NS; j<NT; j++) gs[i][j] = 0.0;
        }

        order(THIRD, t, p, r,
                    NULL, NULL, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        /* Compute derivatives of the chemical potentials */
        for (i=0; i<NA; i++) {
            sum = dgdt;
            for (j=0; j<NR; j++) {
                sum += fr[i][j]*d2gdrdt[j];
                for (k=0; k<NT; k++) sum += fr[i][j]*d2gdrds[j][k]*dsdt[k];
            }
            for (j=0; j<NT; j++) {
                sum += gs[i][j]*d2gdsdt[j];
                for (k=0; k<NT; k++) sum += gs[i][j]*d2gds2[j][k]*dsdt[k];
            }
            mu[i] = sum; /* This is d mu/d t */
            if (returnMixingProperties) mu[i] -= -S(i);
        }

    }

    /*******************************************************************************
    * This return is non-public and is designed for xMELTS calibration purposes.   *
    * The model parameter derivative of the negative of the partial molar entropy  *
    * (i.e. d mu / dT dW) is returned in dw[][]                                    *
    ********************************************************************************/
    if (mask & NINTH) {
        double d2gdsdw[NT][3*NP];
        double d2gdtdw[3*NP];
        double d2gdrds[NR][NT];
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double d3gdrds2[NR][NT][NT];
        double d3gdrdsdt[NR][NT];
        double d3gdrdsdw[NR][NT][3*NP];
        double d3gdrdtdw[NR][3*NP];
        double d3gds3[NT][NT][NT];
        double d3gds2dt[NT][NT];
        double d3gds2dw[NT][NT][3*NP];
        double d3gdsdtdw[NT][3*NP];
        double dsdt[NT], dsdw[NT][3*NP], d2sdtdw[NT][3*NP], gs[NA][NT];
        int k, l, m;

        fillD2GDSDW   (r, s, t, p,  d2gdsdw);
        fillD2GDTDW   (r, s, t, p,  d2gdtdw);
        fillD2GDRDS   (r, s, t, p,  d2gdrds);
        fillD2GDS2    (r, s, t, p,  d2gds2);
        fillD2GDSDT   (r, s, t, p,  d2gdsdt);
        fillD3GDRDS2  (r, s, t, p,  d3gdrds2);
        fillD3GDRDSDT (r, s, t, p,  d3gdrdsdt);
        fillD3GDRDSDW (r, s, t, p,  d3gdrdsdw);
        fillD3GDRDTDW (r, s, t, p,  d3gdrdtdw);
        fillD3GDS3    (r, s, t, p,  d3gds3);
        fillD3GDS2DT  (r, s, t, p,  d3gds2dt);
        fillD3GDS2DW  (r, s, t, p,  d3gds2dw);
        fillD3GDSDTDW (r, s, t, p,  d3gdsdtdw);

        /* fill Darken structures */
        for(i=0; i<NA; i++) {
            for (j=0;  j<NS; j++) gs[i][j] = rsEndmembers[i][NR+j] - s[j];
            for (j=NS; j<NT; j++) gs[i][j] = 0.0;
        }

        order(THIRD | ELEVENTH | TWELFTH, t, p, r,
                    NULL, NULL, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL, dsdw, d2sdtdw);

        for (i=0; i<NA; i++) {
            for (j=0; j<NP; j++) {
                if (modelParameters[j].activeH) {
                    dw[i][j] = d2gdtdw[j];
                    for (l=0; l<NT; l++) dw[i][j] += d2gdsdt[l]*dsdw[l][j];
                    for (k=0; k<NR; k++) {
                        dw[i][j] += fr[i][k]*d3gdrdtdw[k][j];
                        for (l=0; l<NT; l++) dw[i][j] += fr[i][k]*d3gdrdsdt[k][l]*dsdw[l][j];
                        for (m=0; m<NT; m++) {
                            dw[i][j] += fr[i][k]*(d3gdrdsdw[k][m][j]*dsdt[m]+d2gdrds[k][m]*d2sdtdw[m][j]);
                            for (l=0; l<NT; l++) if (dsdw[l][j] != 0.0) dw[i][j] += fr[i][k]*d3gdrds2[k][m][l]*dsdw[l][j]*dsdt[m];
                        }
                    }
                    for (k=0; k<NT; k++) {
                        dw[i][j] += gs[i][k]*d3gdsdtdw[k][j];
                        for (l=0; l<NT; l++) dw[i][j] += gs[i][k]*d3gds2dt[k][l]*dsdw[l][j];
                        for (m=0; m<NT; m++) {
                            dw[i][j] += gs[i][k]*(d3gds2dw[k][m][j]*dsdt[m]+d2gds2[k][m]*d2sdtdw[m][j]);
                            for (l=0; l<NT; l++) if (dsdw[l][j] != 0.0) dw[i][j] += gs[i][k]*d3gds3[k][m][l]*dsdw[l][j]*dsdt[m];
                        }
                    }
                } else dw[i][j] = 0.0;

                if (modelParameters[j].activeS) {
                    dw[i][NP+j] = d2gdtdw[NP+j];
                    for (l=0; l<NT; l++) dw[i][NP+j] += d2gdsdt[l]*dsdw[l][NP+j];
                    for (k=0; k<NR; k++) {
                        dw[i][NP+j] += fr[i][k]*d3gdrdtdw[k][NP+j];
                        for (l=0; l<NT; l++) dw[i][NP+j] += fr[i][k]*d3gdrdsdt[k][l]*dsdw[l][NP+j];
                        for (m=0; m<NT; m++) {
                            dw[i][NP+j] += fr[i][k]*(d3gdrdsdw[k][m][NP+j]*dsdt[m]+d2gdrds[k][m]*d2sdtdw[m][NP+j]);
                            for (l=0; l<NT; l++) if (dsdw[l][NP+j] != 0.0) dw[i][NP+j] += fr[i][k]*d3gdrds2[k][m][l]*dsdw[l][NP+j]*dsdt[m];
                        }
                    }
                    for (k=0; k<NT; k++) {
                        dw[i][NP+j] += gs[i][k]*d3gdsdtdw[k][NP+j];
                        for (l=0; l<NT; l++) dw[i][NP+j] += gs[i][k]*d3gds2dt[k][l]*dsdw[l][NP+j];
                        for (m=0; m<NT; m++) {
                            dw[i][NP+j] += gs[i][k]*(d3gds2dw[k][m][NP+j]*dsdt[m]+d2gds2[k][m]*d2sdtdw[m][NP+j]);
                            for (l=0; l<NT; l++) if (dsdw[l][NP+j] != 0.0) dw[i][NP+j] += gs[i][k]*d3gds3[k][m][l]*dsdw[l][NP+j]*dsdt[m];
                        }
                    }
                } else dw[i][NP+j] = 0.0;

                if (modelParameters[j].activeV) {
                    dw[i][2*NP+j] = d2gdtdw[2*NP+j];
                    for (l=0; l<NT; l++) dw[i][2*NP+j] += d2gdsdt[l]*dsdw[l][2*NP+j];
                    for (k=0; k<NR; k++) {
                        dw[i][2*NP+j] += fr[i][k]*d3gdrdtdw[k][2*NP+j];
                        for (l=0; l<NT; l++) dw[i][2*NP+j] += fr[i][k]*d3gdrdsdt[k][l]*dsdw[l][2*NP+j];
                        for (m=0; m<NT; m++) {
                            dw[i][2*NP+j] += fr[i][k]*(d3gdrdsdw[k][m][2*NP+j]*dsdt[m]+d2gdrds[k][m]*d2sdtdw[m][2*NP+j]);
                            for (l=0; l<NT; l++) if (dsdw[l][2*NP+j] != 0.0) dw[i][2*NP+j] += fr[i][k]*d3gdrds2[k][m][l]*dsdw[l][2*NP+j]*dsdt[m];
                        }
                    }
                    for (k=0; k<NT; k++) {
                        dw[i][2*NP+j] += gs[i][k]*d3gdsdtdw[k][2*NP+j];
                        for (l=0; l<NT; l++) dw[i][2*NP+j] += gs[i][k]*d3gds2dt[k][l]*dsdw[l][2*NP+j];
                        for (m=0; m<NT; m++) {
                            dw[i][2*NP+j] += gs[i][k]*(d3gds2dw[k][m][2*NP+j]*dsdt[m]+d2gds2[k][m]*d2sdtdw[m][2*NP+j]);
                            for (l=0; l<NT; l++) if (dsdw[l][2*NP+j] != 0.0) dw[i][2*NP+j] += gs[i][k]*d3gds3[k][m][l]*dsdw[l][2*NP+j]*dsdt[m];
                        }
                    }
                } else dw[i][2*NP+j] = 0.0;
            }

            if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
                if (modelParameters[NW+i].activeS) dw[i][  NP+NW+i] -= -1.0;
            }

        }
    }

    MTHREAD_MUTEX_UNLOCK(&global_data_mutex);
}

void
gmixLiq(int mask, double t, double p, double *x,
    double *gmix, /* Gibbs energy of mixing             BINARY MASK: 0001 */
    double *dx,   /* (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010 */
    double **dx2  /* (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100 */
    )
{
    double *r = x;
    double s[NT];
    int i;

    liqERRstate = ERR_NONE;

    if ((calculationMode == MODE__MELTS) || (calculationMode == MODE_pMELTS)) {
            gmixLiq_v34(mask, t, p, x, gmix, dx, dx2);
            return;
    } else if (calculationMode == MODE__MELTSandCO2) {
            gmixLiq_CO2(mask, t, p, x, gmix, dx, dx2);
            return;
    } else if (calculationMode == MODE__MELTSandCO2_H2O) {
            gmixLiq_CO2_H2O(mask, t, p, x, gmix, dx, dx2);
            return;
    }

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    MTHREAD_MUTEX_LOCK(&global_data_mutex);
    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    if (mask & FIRST) {
        *gmix  = fillG (r, s, t, p);
        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *gmix -= G(0);
            for (i=0; i<NR; i++) *gmix += r[i]*(G(0)-G(i+1));
        }
    }

    if(mask & SECOND) {
        fillDGDR (r, s, t, p, dx);
        if (returnMixingProperties) for (i=0; i<NR; i++) dx[i] += (G(0)-G(i+1));
    }

    if(mask & THIRD) {
        double d2gdr2[NR][NR];
        double d2gdrds[NR][NT];
        double d2gds2[NT][NT];
        double dsdr[NT][NR];
        int j, k, l;

        fillD2GDR2  (r, s, t, p,  d2gdr2);
        fillD2GDRDS (r, s, t, p,  d2gdrds);
        fillD2GDS2  (r, s, t, p,  d2gds2);

        order(SECOND, t, p, r,
                    NULL, dsdr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dx2[i][j] = d2gdr2[i][j];
                for (k=0; k<NT; k++) {
                    dx2[i][j] += d2gdrds[i][k]*dsdr[k][j] + d2gdrds[j][k]*dsdr[k][i];
                    for (l=0; l<NT; l++) dx2[i][j] += d2gds2[k][l]*dsdr[k][i]*dsdr[l][j];
                }
            }
        }
    }

    if(mask & FOURTH) {
        double d3gdr3[NR][NR][NR];
        double d3gdr2ds[NR][NR][NT];
        double d3gdrds2[NR][NT][NT];
        double d3gds3[NT][NT][NT];
        double dx3[NR][NR][NR]; /* This should be passed to the function */
        double dsdr[NT][NR], d2sdr2[NT][NR][NR];
        int i, j, k, l, m, n;

        fillD3GDR3   (r, s, t, p,  d3gdr3);
        fillD3GDR2DS (r, s, t, p,  d3gdr2ds);
        fillD3GDRDS2 (r, s, t, p,  d3gdrds2);
        fillD3GDS3   (r, s, t, p,  d3gds3);

        order(SECOND, t, p, r,
                    NULL, dsdr, NULL, NULL, d2sdr2, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                for (k=0; k<NR; k++) {
                    dx3[i][j][k] = d3gdr3[i][j][k];
                    for (l=0; l<NT; l++) {
                        dx3[i][j][k] += d3gdr2ds[i][j][l]*dsdr[l][k] +
                            d3gdr2ds[j][k][l]*dsdr[l][i] + d3gdr2ds[k][i][l]*dsdr[l][j];
                        for (m=0; m<NT; m++) {
                            dx3[i][j][k] +=
                                d3gdrds2[i][l][m]*dsdr[l][j]*dsdr[m][k] +
                                d3gdrds2[j][l][m]*dsdr[l][k]*dsdr[m][i] +
                                d3gdrds2[k][l][m]*dsdr[l][i]*dsdr[m][j];
                            for (n=0; n<NT; n++)
                                dx3[i][j][k] +=
                                    d3gds3[l][m][n]*dsdr[l][i]*dsdr[m][j]*dsdr[n][k];
                        }
                    }
                }
            }
        }
    }

    MTHREAD_MUTEX_UNLOCK(&global_data_mutex);
}

void
hmixLiq(int mask, double t, double p, double *x,
    double *hmix, /* Enthalpy of mixing BINARY MASK:                 01 */
    double *dw    /* (pointer to dw[]) d(H)/d(x[])      BINARY MASK: 10 */
    )
{
    double *r = x;
    double s[NT], dsdt[NT];
    int i;

    liqERRstate = ERR_NONE;

    if ((calculationMode == MODE__MELTS) || (calculationMode == MODE_pMELTS)) {
            hmixLiq_v34 (mask, t, p, x, hmix);
            return;
    } else if (calculationMode == MODE__MELTSandCO2) {
            hmixLiq_CO2(mask, t, p, x, hmix);
            return;
    } else if (calculationMode == MODE__MELTSandCO2_H2O) {
            hmixLiq_CO2_H2O(mask, t, p, x, hmix);
            return;
    }

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    MTHREAD_MUTEX_LOCK(&global_data_mutex);
    order(FIRST | THIRD, t, p, r,
                s, NULL, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    if (mask & FIRST) {
        *hmix = fillG (r, s, t, p) + t*fillS (r, s, t, p);                            /* was: *hmix = fillH (r, s, t, p);  */

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *hmix -= (G(0) + t*S(0));                                                   /* was: *hmix -= H(0);               */
            for (i=0; i<NR; i++) *hmix += r[i]*((G(0) + t*S(0))-(G(i+1) + t*S(i+1)));   /* was: *hmix += r[i]*(H(0)-H(i+1)); */
        }
    }

    if (mask & SECOND) {
        double dgdw[3*NP];
        double d2gdsdw[NT][3*NP];
        double d2gdtdw[3*NP];
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double dsdw[NT][3*NP], dsdt[NT];
        int i, k, l;

        fillDGDW    (r, s, t, p, dgdw);
        fillD2GDSDW (r, s, t, p, d2gdsdw);
        fillD2GDTDW (r, s, t, p, d2gdtdw);
        fillD2GDS2  (r, s, t, p, d2gds2);
        fillD2GDSDT (r, s, t, p, d2gdsdt);

        order(THIRD | ELEVENTH, t, p, r,
                    NULL, NULL, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL, dsdw, NULL);

        for (i=0; i<NP; i++) {
            if (modelParameters[i].activeH) {
                dw[i] = d2gdtdw[i];
                for (k=0; k<NT; k++) {
                    dw[i] += d2gdsdw[k][i]*dsdt[k] + d2gdsdt[k]*dsdw[k][i];
                    for (l=0; l<NT; l++) dw[i] += d2gds2[k][l]*dsdt[k]*dsdw[l][i] ;
                }
                dw[i] = dgdw[i] - t*dw[i];
            } else dw[i] = 0.0;

            if (modelParameters[i].activeS) {
                dw[NP+i] = d2gdtdw[NP+i];
                for (k=0; k<NT; k++) {
                    dw[NP+i] += d2gdsdw[k][NP+i]*dsdt[k] + d2gdsdt[k]*dsdw[k][NP+i];
                    for (l=0; l<NT; l++) dw[NP+i] += d2gds2[k][l]*dsdt[k]*dsdw[l][NP+i] ;
                }
                dw[NP+i] = dgdw[NP+i] - t*dw[NP+i];
            } else dw[NP+i] = 0.0;

            if (modelParameters[i].activeV) {
                dw[2*NP+i] = d2gdtdw[2*NP+i];
                for (k=0; k<NT; k++) {
                    dw[2*NP+i] += d2gdsdw[k][2*NP+i]*dsdt[k] + d2gdsdt[k]*dsdw[k][2*NP+i];
                    for (l=0; l<NT; l++) dw[2*NP+i] += d2gds2[k][l]*dsdt[k]*dsdw[l][2*NP+i] ;
                }
                dw[2*NP+i] = dgdw[2*NP+i] - t*dw[2*NP+i];
            } else dw[2*NP+i] = 0.0;
        }

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            if (modelParameters[NW+0].activeH) dw[     NW+0] -= 1.0;
#ifndef USE_GHIORSO_KRESS_MODEL
            if (modelParameters[NW+0].activeV) dw[2*NP+NW+0] -= (p-1.0);
#endif /* USE_GHIORSO_KRESS_MODEL */
            for (i=0; i<NR; i++) {
                if (modelParameters[NW  +0].activeH) dw[     NW  +0] +=         r[i];
#ifndef USE_GHIORSO_KRESS_MODEL
                if (modelParameters[NW  +0].activeV) dw[2*NP+NW  +0] += (p-1.0)*r[i];
#endif /* USE_GHIORSO_KRESS_MODEL */
                if (modelParameters[NW+i+1].activeH) dw[     NW+i+1] -=         r[i];
#ifndef USE_GHIORSO_KRESS_MODEL
                if (modelParameters[NW+i+1].activeV) dw[2*NP+NW+i+1] -= (p-1.0)*r[i];
#endif /* USE_GHIORSO_KRESS_MODEL */
            }
        }

        for (i=0; i<(3*NP); i++) if (fabs(dw[i]) < 10.0*DBL_EPSILON) dw[i] = 0.0;
    }

    MTHREAD_MUTEX_UNLOCK(&global_data_mutex);
}

void
smixLiq(int mask, double t, double p, double *x,
    double *smix, /* Entropy of mixing                  BINARY MASK: 0001 */
    double *dx,   /* (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 0010 */
    double **dx2, /* (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 0100 */
    double *dw    /* (pointer to dw[]) d(s)/d(x[])      BINARY MASK: 1000 */
    )
{
    double *r = x;
    double s[NT];

    liqERRstate = ERR_NONE;

    if ((calculationMode == MODE__MELTS) || (calculationMode == MODE_pMELTS)) {
            smixLiq_v34(mask, t, p, x, smix, dx, dx2);
            return;
    } else if (calculationMode == MODE__MELTSandCO2) {
            smixLiq_CO2(mask, t, p, x, smix, dx, dx2);
            return;
    } else if (calculationMode == MODE__MELTSandCO2_H2O) {
            smixLiq_CO2_H2O(mask, t, p, x, smix, dx, dx2);
            return;
    }

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    MTHREAD_MUTEX_LOCK(&global_data_mutex);
    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    if (mask & FIRST) {
        int i;

        *smix = fillS (r, s, t, p);

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *smix -= S(0);
            for (i=0; i<NR; i++) *smix += r[i]*(S(0)-S(i+1));
        }
    }

    if(mask & SECOND) {
        double d2gdrds[NR][NT];
        double d2gdrdt[NR];
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double dsdr[NT][NR], dsdt[NT];
        int i, k, l;

        fillD2GDRDS (r, s, t, p, d2gdrds);
        fillD2GDRDT (r, s, t, p, d2gdrdt);
        fillD2GDS2  (r, s, t, p, d2gds2);
        fillD2GDSDT (r, s, t, p, d2gdsdt);

        order(SECOND | THIRD, t, p, r,
                    NULL, dsdr, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            dx[i] = d2gdrdt[i];
            for (k=0; k<NT; k++) {
                dx[i] += d2gdrds[i][k]*dsdt[k] + d2gdsdt[k]*dsdr[k][i];
                for (l=0; l<NT; l++) dx[i] += d2gds2[k][l]*dsdt[k]*dsdr[l][i] ;
            }
            dx[i] *= -1.0;
        }

        /* Convert Solution Properties -> Mixing Properties */
        if (returnMixingProperties) for (i=0; i<NR; i++) dx[i] += (S(0)-S(i+1));
    }

    if(mask & THIRD) {
        double d2gdrds[NR][NT];
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double d3gdr2ds[NR][NR][NT];
        double d3gdr2dt[NR][NR];
        double d3gdrds2[NR][NT][NT];
        double d3gdrdsdt[NR][NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dt[NT][NT];
        double dsdr[NT][NR], dsdt[NT], d2sdr2[NT][NR][NR], d2sdrdt[NT][NR];
        int i, j, k, l, m;

        fillD2GDRDS   (r, s, t, p, d2gdrds);
        fillD2GDS2    (r, s, t, p, d2gds2);
        fillD2GDSDT   (r, s, t, p, d2gdsdt);
        fillD3GDR2DS  (r, s, t, p, d3gdr2ds);
        fillD3GDR2DT  (r, s, t, p, d3gdr2dt);
        fillD3GDRDS2  (r, s, t, p, d3gdrds2);
        fillD3GDRDSDT (r, s, t, p, d3gdrdsdt);
        fillD3GDS3    (r, s, t, p, d3gds3);
        fillD3GDS2DT  (r, s, t, p, d3gds2dt);

        order(SECOND | THIRD | FIFTH | SIXTH, t, p, r,
                    NULL, dsdr, dsdt, NULL, d2sdr2, d2sdrdt, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dx2[i][j] = d3gdr2dt[i][j];
                for (k=0; k<NT; k++) {
                    dx2[i][j] += d3gdr2ds[i][j][k]*dsdt[k]
                     + d3gdrdsdt[i][k]*dsdr[k][j]
                     + d3gdrdsdt[j][k]*dsdr[k][i]
                     + d2gdsdt[k]*d2sdr2[k][i][j]
                     + d2gdrds[i][k]*d2sdrdt[k][j]
                     + d2gdrds[j][k]*d2sdrdt[k][i];
                    for (l=0; l<NT; l++) {
                        dx2[i][j] += d3gdrds2[i][k][l]*dsdr[k][j]*dsdt[l]
                       + d3gdrds2[j][k][l]*dsdr[k][i]*dsdt[l]
                       + d2gds2[k][l]*d2sdr2[k][i][j]*dsdt[l]
                       + d3gds2dt[k][l]*dsdr[k][i]*dsdr[l][j]
                       + d2gds2[k][l]*dsdr[k][i]*d2sdrdt[l][j]
                       + d2gds2[k][l]*dsdr[k][j]*d2sdrdt[l][i];
                        for (m=0; m<NT; m++)
                            dx2[i][j] += d3gds3[k][l][m]*dsdr[k][i]*dsdr[l][j]*dsdt[m];
                    }
                }
                dx2[i][j] *= -1.0;
            }
        }
    }

    if(mask & FOURTH) {
        double d2gdsdw[NT][3*NP];
        double d2gdtdw[3*NP];
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double dsdw[NT][3*NP], dsdt[NT];
        int i, k, l;

        fillD2GDSDW (r, s, t, p, d2gdsdw);
        fillD2GDTDW (r, s, t, p, d2gdtdw);
        fillD2GDS2  (r, s, t, p, d2gds2);
        fillD2GDSDT (r, s, t, p, d2gdsdt);

        order(THIRD | ELEVENTH, t, p, r,
                    NULL, NULL, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL, dsdw, NULL);

        for (i=0; i<NP; i++) {
            if (modelParameters[i].activeH) {
                dw[i] = d2gdtdw[i];
                for (k=0; k<NT; k++) {
                    dw[i] += d2gdsdw[k][i]*dsdt[k] + d2gdsdt[k]*dsdw[k][i];
                    for (l=0; l<NT; l++) dw[i] += d2gds2[k][l]*dsdt[k]*dsdw[l][i] ;
                }
                dw[i] *= -1.0;
            } else dw[i] = 0.0;

            if (modelParameters[i].activeS) {
                dw[NP+i] = d2gdtdw[NP+i];
                for (k=0; k<NT; k++) {
                    dw[NP+i] += d2gdsdw[k][NP+i]*dsdt[k] + d2gdsdt[k]*dsdw[k][NP+i];
                    for (l=0; l<NT; l++) dw[NP+i] += d2gds2[k][l]*dsdt[k]*dsdw[l][NP+i] ;
                }
                dw[NP+i] *= -1.0;
            } else dw[NP+i] = 0.0;

            if (modelParameters[i].activeV) {
                dw[2*NP+i] = d2gdtdw[2*NP+i];
                for (k=0; k<NT; k++) {
                    dw[2*NP+i] += d2gdsdw[k][2*NP+i]*dsdt[k] + d2gdsdt[k]*dsdw[k][2*NP+i];
                    for (l=0; l<NT; l++) dw[2*NP+i] += d2gds2[k][l]*dsdt[k]*dsdw[l][2*NP+i] ;
                }
                dw[2*NP+i] *= -1.0;
            } else dw[2*NP+i] = 0.0;
        }

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            if (modelParameters[NW+0].activeS) dw[NP+NW+0] -= 1.0;
            for (i=0; i<NR; i++) {
                if (modelParameters[NW  +0].activeS) dw[NP+NW  +0] += r[i];
                if (modelParameters[NW+i+1].activeS) dw[NP+NW+i+1] -= r[i];
            }
        }

        for (i=0; i<(3*NP); i++) if (fabs(dw[i]) < 10.0*DBL_EPSILON) dw[i] = 0.0;
    }

    MTHREAD_MUTEX_UNLOCK(&global_data_mutex);
}

void
cpmixLiq(int mask, double t, double p, double *x,
    double *cpmix, /* Heat capacity of mixing               BINARY MASK: 001 */
    double *dt,    /* d(cp)/d(t)                            BINARY MASK: 010 */
    double *dx     /* d(cp)/d(x[])d(t)                      BINARY MASK: 100 */
    )
{
    double *r = x;
    double d2gdsdt[NT];
    double d2gds2[NT][NT];
    double s[NT], dsdt[NT];
    double d2gdt2;
    int i, j;

    liqERRstate = ERR_NONE;

    if ((calculationMode == MODE__MELTS) || (calculationMode == MODE_pMELTS)) {
            cpmixLiq_v34(mask, t, p, x, cpmix, dt, dx);
            return;
    } else if (calculationMode == MODE__MELTSandCO2) {
            cpmixLiq_CO2(mask, t, p, x, cpmix, dt, dx);
            return;
    } else if (calculationMode == MODE__MELTSandCO2_H2O) {
            cpmixLiq_CO2_H2O(mask, t, p, x, cpmix, dt, dx);
            return;
    }

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    MTHREAD_MUTEX_LOCK(&global_data_mutex);
    order(FIRST | THIRD, t, p, r,
                s,    NULL, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    fillD2GDSDT (r, s, t, p, d2gdsdt);
    fillD2GDS2  (r, s, t, p, d2gds2);
    d2gdt2 = fillD2GDT2  (r, s, t, p);

    if (mask & FIRST) {
        *cpmix = d2gdt2;
        for (i=0; i<NT; i++) {
            *cpmix += 2.0*d2gdsdt[i]*dsdt[i];
            for (j=0; j<NT; j++) *cpmix += d2gds2[i][j]*dsdt[i]*dsdt[j];
        }
        *cpmix *= -t;

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *cpmix -= CP(0);
            for (i=0; i<NR; i++) *cpmix += r[i]*(CP(0)-CP(i+1));
        }
    }

    if(mask & SECOND) {
        double d3gds3[NT][NT][NT];
        double d3gds2dt[NT][NT];
        double d3gdsdt2[NT];
        double d3gdt3 = fillD3GDT3   (r, s, t, p);
        double d2sdt2[NT], temp;
        int k;

        fillD3GDS3   (r, s, t, p, d3gds3);
        fillD3GDS2DT (r, s, t, p, d3gds2dt);
        fillD3GDSDT2 (r, s, t, p, d3gdsdt2);

        order(EIGHTH, t, p, r,
                    NULL, NULL, NULL, NULL, NULL, NULL, NULL, d2sdt2, NULL, NULL, NULL, NULL);

        /* compute d2gdt2 */
        temp = d2gdt2;
        for (i=0; i<NT; i++) {
            temp += 2.0*d2gdsdt[i]*dsdt[i];
            for (j=0; j<NT; j++) temp += d2gds2[i][j]*dsdt[i]*dsdt[j];
        }

        *dt = d3gdt3;
        for (i=0; i<NT; i++) {
            *dt += 3.0*d3gdsdt2[i]*dsdt[i] + 3.0*d2gdsdt[i]*d2sdt2[i];
            for (j=0; j<NT; j++) {
                *dt += 3.0*d2gds2[i][j]*dsdt[i]*d2sdt2[j]
             + 3.0*d3gds2dt[i][j]*dsdt[i]*dsdt[j];
                for (k=0; k<NT; k++) *dt += d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdt[k];
            }
        }
        *dt = -t*(*dt) - temp;

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *dt -= DCPDT(0);
            for (i=0; i<NR; i++) *dt += r[i]*(DCPDT(0)-DCPDT(i+1));
        }
    }

    if(mask & THIRD) {
        double d3gds3[NT][NT][NT];
        double d3gdrds2[NR][NT][NT];
        double d3gdrdsdt[NR][NT];
        double d3gds2dt[NT][NT];
        double d2gdrds[NR][NT];
        double d3gdrdt2[NR];
        double d3gdsdt2[NT];
        double dsdr[NT][NR], d2sdrdt[NT][NR], d2sdt2[NT];
        int k, l;

        fillD3GDS3    (r, s, t, p, d3gds3);
        fillD3GDRDS2  (r, s, t, p, d3gdrds2);
        fillD3GDRDSDT (r, s, t, p, d3gdrdsdt);
        fillD3GDS2DT  (r, s, t, p, d3gds2dt);
        fillD2GDRDS   (r, s, t, p, d2gdrds);
        fillD3GDRDT2  (r, s, t, p, d3gdrdt2);
        fillD3GDSDT2  (r, s, t, p, d3gdsdt2);

        order(SECOND | SIXTH | EIGHTH, t, p, r,
                    NULL, dsdr, NULL, NULL, NULL, d2sdrdt, NULL, d2sdt2, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0,dx[i]=d3gdrdt2[i]; j<NT; j++) {
                dx[i] += d3gdsdt2[j]*dsdr[j][i] + 2.0*d2gdsdt[j]*d2sdrdt[j][i] +
                 2.0*d3gdrdsdt[i][j]*dsdt[j] + d2gdrds[i][j]*d2sdt2[j];
                for (k=0; k<NT; k++) {
                    dx[i] += d3gdrds2[i][j][k]*dsdt[j]*dsdt[k] +
                   2.0*d2gds2[j][k]*dsdt[j]*d2sdrdt[k][i] +
                   2.0*d3gds2dt[j][k]*dsdr[j][i]*dsdt[k] +
                   d2gds2[j][k]*dsdr[j][i]*d2sdt2[k];
                    for (l=0; l<NT; l++)
                        dx[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdt[k]*dsdt[l];
                }
            }
            dx[i] *= -t;
        }

        if (returnMixingProperties) for (i=0; i<NR; i++) dx[i] += (CP(0)-CP(i+1));
    }

    MTHREAD_MUTEX_UNLOCK(&global_data_mutex);
}

void
vmixLiq(int mask, double t, double p, double *x,
    double *vmix, /* Volume of mixing                BINARY MASK: 00000000001 */
    double *dx,   /* (pointer to dx[]) d(v)/d(x[])   BINARY MASK: 00000000010 */
    double **dx2, /* (point to dx2[][]) d(v)/d(x[])2 BINARY MASK: 00000000100 */
    double *dt,   /* d(v)/d(t)                       BINARY MASK: 00000001000 */
    double *dp,   /* d(v)/d(p)                       BINARY MASK: 00000010000 */
    double *dt2,  /* d2(v)/d(t)2                     BINARY MASK: 00000100000 */
    double *dtdp, /* d2(v)/d(t)d(p)                  BINARY MASK: 00001000000 */
    double *dp2,  /* d2(v)/d(p)2                     BINARY MASK: 00010000000 */
    double *dxdt, /* d2(v)/d(x[])d(t)                BINARY MASK: 00100000000 */
    double *dxdp, /* d2(v)/d(x[])d(p)                BINARY MASK: 01000000000 */
    double *dw    /* (pointer to dw[]) d(v)/d(x[])   BINARY MASK: 10000000000 */
    )
{
    double *r = x;
    double s[NT];

    liqERRstate = ERR_NONE;

    if ((calculationMode == MODE__MELTS) || (calculationMode == MODE_pMELTS)) {
        vmixLiq_v34 (mask, t, p, x, vmix, dx, dx2, dt, dp, dt2, dtdp, dp2, dxdt, dxdp);
        return;
    } else if (calculationMode == MODE__MELTSandCO2) {
            vmixLiq_CO2 (mask, t, p, x, vmix, dx, dx2, dt, dp, dt2, dtdp, dp2, dxdt, dxdp);
            return;
    } else if (calculationMode == MODE__MELTSandCO2_H2O) {
            vmixLiq_CO2_H2O (mask, t, p, x, vmix, dx, dx2, dt, dp, dt2, dtdp, dp2, dxdt, dxdp);
            return;
    }

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    MTHREAD_MUTEX_LOCK(&global_data_mutex);
    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

#ifdef USE_GHIORSO_KRESS_MODEL
#undef  V
#undef  DVDT
#undef  DVDP
#undef  D2VDT2
#undef  D2VDTDP
#undef  D2VDP2
#define V(i)       ((liquid[i].cur).v + modelParameters[NW+i].volume)
#define DVDT(i)    ((liquid[i].cur).dvdt)
#define DVDP(i)    ((liquid[i].cur).dvdp)
#define D2VDT2(i)  ((liquid[i].cur).d2vdt2)
#define D2VDTDP(i) ((liquid[i].cur).d2vdtdp)
#define D2VDP2(i)  ((liquid[i].cur).d2vdp2)
#endif /* USE_GHIORSO_KRESS_MODEL */

    if (mask & FIRST) {
        int i;
        *vmix = fillV (r, s, t, p);
        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *vmix -= V(0);
            for (i=0; i<NR; i++) *vmix += r[i]*(V(0)-V(i+1));
        }
    }

    if(mask & SECOND) {
        double d2gdrds[NR][NT];
        double d2gdrdp[NR];
        double d2gds2[NT][NT];
        double d2gdsdp[NT];
        double dsdr[NT][NR], dsdp[NT];
        int i, j, k;

        fillD2GDRDS (r, s, t, p, d2gdrds);
        fillD2GDRDP (r, s, t, p, d2gdrdp);
        fillD2GDS2  (r, s, t, p, d2gds2);
        fillD2GDSDP (r, s, t, p, d2gdsdp);

        order(SECOND | FOURTH, t, p, r,
                    NULL, dsdr, NULL, dsdp, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            dx[i] = d2gdrdp[i];
            for (j=0; j<NT; j++) {
                dx[i] += d2gdrds[i][j]*dsdp[j] + d2gdsdp[j]*dsdr[j][i];
                for (k=0; k<NT; k++) dx[i] += d2gds2[j][k]*dsdp[j]*dsdr[k][i];
            }
        }

        if (returnMixingProperties) for (i=0; i<NR; i++) dx[i] += (V(0)-V(i+1));
    }

    if(mask & THIRD) {
        double d2gdrds[NR][NT];
        double d2gds2[NT][NT];
        double d2gdsdp[NT];
        double d3gdr2ds[NR][NR][NT];
        double d3gdr2dp[NR][NR];
        double d3gdrds2[NR][NT][NT];
        double d3gdrdsdp[NR][NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dp[NT][NT];
        double dsdr[NT][NR], dsdp[NT], d2sdr2[NT][NR][NR], d2sdrdp[NT][NR];
        int i, j, k, l, m;

        fillD2GDRDS   (r, s, t, p, d2gdrds);
        fillD2GDS2    (r, s, t, p, d2gds2);
        fillD2GDSDP   (r, s, t, p, d2gdsdp);
        fillD3GDR2DS  (r, s, t, p, d3gdr2ds);
        fillD3GDR2DP  (r, s, t, p, d3gdr2dp);
        fillD3GDRDS2  (r, s, t, p, d3gdrds2);
        fillD3GDRDSDP (r, s, t, p, d3gdrdsdp);
        fillD3GDS3    (r, s, t, p, d3gds3);
        fillD3GDS2DP  (r, s, t, p, d3gds2dp);

        order(SECOND | FOURTH | FIFTH | SEVENTH, t, p, r,
                    NULL, dsdr, NULL, dsdp, d2sdr2, NULL, d2sdrdp,  NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dx2[i][j] = d3gdr2dp[i][j];
                for (k=0; k<NT; k++) {
                    dx2[i][j] += d3gdr2ds[i][j][k]*dsdp[k]
                     + d3gdrdsdp[i][k]*dsdr[k][j]
                     + d3gdrdsdp[j][k]*dsdr[k][i]
                     + d2gdsdp[k]*d2sdr2[k][i][j]
                     + d2gdrds[i][k]*d2sdrdp[k][j]
                     + d2gdrds[j][k]*d2sdrdp[k][i];
                    for (l=0; l<NT; l++) {
                        dx2[i][j] += d3gdrds2[i][k][l]*dsdr[k][j]*dsdp[l]
                       + d3gdrds2[j][k][l]*dsdr[k][i]*dsdp[l]
                       + d2gds2[k][l]*d2sdr2[k][i][j]*dsdp[l]
                       + d3gds2dp[k][l]*dsdr[k][i]*dsdr[l][j]
                       + d2gds2[k][l]*dsdr[k][i]*d2sdrdp[l][j]
                       + d2gds2[k][l]*dsdr[k][j]*d2sdrdp[l][i];
                        for (m=0; m<NT; m++)
                            dx2[i][j] += d3gds3[k][l][m]*dsdr[k][i]*dsdr[l][j]*dsdp[m];
                    }
                }
            }
        }

    }

    if(mask & FOURTH) {
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double d2gdsdp[NT];
        double d2gdtdp = fillD2GDTDP (r, s, t, p);
        double dsdt[NT], dsdp[NT];
        int i, j;

        fillD2GDS2  (r, s, t, p, d2gds2);
        fillD2GDSDT (r, s, t, p, d2gdsdt);
        fillD2GDSDP (r, s, t, p, d2gdsdp);

        order(THIRD | FOURTH, t, p, r,
                    NULL, NULL, dsdt, dsdp, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        *dt = d2gdtdp;
        for (i=0; i<NT; i++) {
            *dt += d2gdsdt[i]*dsdp[i] + d2gdsdp[i]*dsdt[i];
            for (j=0; j<NT; j++) *dt += d2gds2[i][j]*dsdt[i]*dsdp[j];
        }

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *dt -= DVDT(0);
            for (i=0; i<NR; i++) *dt += r[i]*(DVDT(0)-DVDT(i+1));
        }
    }

    if(mask & FIFTH) {
        double d2gds2[NT][NT];
        double d2gdsdp[NT];
        double d2gdp2  = fillD2GDP2  (r, s, t, p);
        double dsdp[NT];
        int i,j;

        fillD2GDS2  (r, s, t, p, d2gds2);
        fillD2GDSDP (r, s, t, p, d2gdsdp);

        order(FOURTH, t, p, r,
                    NULL, NULL, NULL, dsdp, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        *dp = d2gdp2;
        for (i=0; i<NT; i++) {
            *dp += 2.0*d2gdsdp[i]*dsdp[i];
            for (j=0; j<NT; j++) *dp += d2gds2[i][j]*dsdp[i]*dsdp[j];
        }

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *dp -= DVDP(0);
            for (i=0; i<NR; i++) *dp += r[i]*(DVDP(0)-DVDP(i+1));
        }
    }

    if(mask & SIXTH) {
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double d2gdsdp[NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dp[NT][NT];
        double d3gds2dt[NT][NT];
        double d3gdsdtdp[NT];
        double d3gdsdt2[NT];
        double d3gdt2dp  = fillD3GDT2DP  (r, s, t, p);
        double dsdt[NT], dsdp[NT], d2sdt2[NT], d2sdtdp[NT];
        int i, j, k;

        fillD2GDS2    (r, s, t, p, d2gds2);
        fillD2GDSDT   (r, s, t, p, d2gdsdt);
        fillD2GDSDP   (r, s, t, p, d2gdsdp);
        fillD3GDS3    (r, s, t, p, d3gds3);
        fillD3GDS2DP  (r, s, t, p, d3gds2dp);
        fillD3GDS2DT  (r, s, t, p, d3gds2dt);
        fillD3GDSDTDP (r, s, t, p, d3gdsdtdp);
        fillD3GDSDT2  (r, s, t, p, d3gdsdt2);

        order(THIRD | FOURTH | EIGHTH | NINTH, t, p, r,
                    NULL, NULL, dsdt, dsdp, NULL, NULL, NULL, d2sdt2, d2sdtdp, NULL, NULL, NULL);

        *dt2 = d3gdt2dp;
        for (i=0; i<NT; i++) {
            *dt2 += d3gdsdt2[i]*dsdp[i] + 2.0*d2gdsdt[i]*d2sdtdp[i]
                        + d2gdsdp[i]*d2sdt2[i] + 2.0*d3gdsdtdp[i]*dsdt[i];
            for (j=0; j<NT; j++) {
                *dt2 += 2.0*d3gds2dt[i][j]*dsdt[i]*dsdp[j]
                            + d2gds2[i][j]*d2sdt2[i]*dsdp[j]
                            + 2.0*d2gds2[i][j]*dsdt[i]*d2sdtdp[j]
                            + d3gds2dp[i][j]*dsdt[i]*dsdt[j];
                for (k=0; k<NT; k++) *dt2 += d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdp[k];
            }
        }

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *dt2 -= D2VDT2(0);
            for (i=0; i<NR; i++) *dt2 += r[i]*(D2VDT2(0)-D2VDT2(i+1));
        }
    }

    if(mask & SEVENTH) {
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double d2gdsdp[NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dt[NT][NT];
        double d3gds2dp[NT][NT];
        double d3gdsdtdp[NT];
        double d3gdsdp2[NT];
        double    d3gdtdp2  = fillD3GDTDP2  (r, s, t, p);
        double dsdt[NT], dsdp[NT], d2sdtdp[NT], d2sdp2[NT];
        int i, j, k;

        fillD2GDS2    (r, s, t, p, d2gds2);
        fillD2GDSDT   (r, s, t, p, d2gdsdt);
        fillD2GDSDP   (r, s, t, p, d2gdsdp);
        fillD3GDS3    (r, s, t, p, d3gds3);
        fillD3GDS2DT  (r, s, t, p, d3gds2dt);
        fillD3GDS2DP  (r, s, t, p, d3gds2dp);
        fillD3GDSDTDP (r, s, t, p, d3gdsdtdp);
        fillD3GDSDP2  (r, s, t, p, d3gdsdp2);

        order(THIRD | FOURTH | NINTH | TENTH, t, p, r,
                    NULL, NULL, dsdt, dsdp, NULL, NULL, NULL, NULL, d2sdtdp, d2sdp2, NULL, NULL);

        *dtdp = d3gdtdp2;
        for (i=0; i<NT; i++) {
            *dtdp += 2.0*d3gdsdtdp[i]*dsdp[i] + d2gdsdt[i]*d2sdp2[i]
             + 2.0*d2gdsdp[i]*d2sdtdp[i] + d3gdsdp2[i]*dsdt[i];
            for (j=0; j<NT; j++) {
                *dtdp += 2.0*d3gds2dp[i][j]*dsdt[i]*dsdp[j]
               + d2gds2[i][j]*dsdt[i]*d2sdp2[j]
               + 2.0*d2gds2[i][j]*d2sdtdp[i]*dsdp[j]
               + d3gds2dt[i][j]*dsdp[i]*dsdp[j];
                for (k=0; k<NT; k++) *dtdp += d3gds3[i][j][k]*dsdt[i]*dsdp[j]*dsdp[k];
            }
        }

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *dtdp -= D2VDTDP(0);
            for (i=0; i<NR; i++) *dtdp += r[i]*(D2VDTDP(0)-D2VDTDP(i+1));
        }
    }

    if(mask & EIGHTH) {
        double d2gds2[NT][NT];
        double d2gdsdp[NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dp[NT][NT];
        double d3gdsdp2[NT];
        double    d3gdp3   = fillD3GDP3   (r, s, t, p);
        double dsdp[NT], d2sdp2[NT];
        int i, j, k;

        fillD2GDS2   (r, s, t, p, d2gds2);
        fillD2GDSDP  (r, s, t, p, d2gdsdp);
        fillD3GDS3   (r, s, t, p, d3gds3);
        fillD3GDS2DP (r, s, t, p, d3gds2dp);
        fillD3GDSDP2 (r, s, t, p, d3gdsdp2);

        order(FOURTH | TENTH, t, p, r,
                    NULL, NULL, NULL, dsdp, NULL, NULL, NULL, NULL, NULL, d2sdp2, NULL, NULL);

        *dp2 = d3gdp3;
        for (i=0; i<NT; i++) {
            *dp2 += 3.0*d3gdsdp2[i]*dsdp[i] + 3.0*d2gdsdp[i]*d2sdp2[i];
            for (j=0; j<NT; j++) {
                *dp2 += 3.0*d2gds2[i][j]*dsdp[i]*d2sdp2[j]
                            + 3.0*d3gds2dp[i][j]*dsdp[i]*dsdp[j];
                for (k=0; k<NT; k++) *dp2 += d3gds3[i][j][k]*dsdp[i]*dsdp[j]*dsdp[k];
            }
        }

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *dp2 -= D2VDP2(0);
            for (i=0; i<NR; i++) *dp2 += r[i]*(D2VDP2(0)-D2VDP2(i+1));
        }
    }

    if(mask & NINTH) {
        double d3gds3[NT][NT][NT];
        double d3gdrds2[NR][NT][NT];
        double d3gdrdsdt[NR][NT];
        double d3gds2dp[NT][NT];
        double d2gdrds[NR][NT];
        double d3gdrdtdp[NR];
        double d3gdsdtdp[NT];
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double d3gdrdsdp[NR][NT];
        double d2gdsdp[NT];
        double d3gds2dt[NT][NT];
        double dsdt[NT], dsdp[NT], dsdr[NT][NR], d2sdrdt[NT][NR], d2sdrdp[NT][NR],
            d2sdtdp[NT];
        int i, j, k, l;

        fillD3GDS3    (r, s, t, p, d3gds3);
        fillD3GDRDS2  (r, s, t, p, d3gdrds2);
        fillD3GDRDSDT (r, s, t, p, d3gdrdsdt);
        fillD3GDS2DP  (r, s, t, p, d3gds2dp);
        fillD2GDRDS   (r, s, t, p, d2gdrds);
        fillD3GDRDTDP (r, s, t, p, d3gdrdtdp);
        fillD3GDSDTDP (r, s, t, p, d3gdsdtdp);
        fillD2GDS2    (r, s, t, p, d2gds2);
        fillD2GDSDT   (r, s, t, p, d2gdsdt);
        fillD3GDRDSDP (r, s, t, p, d3gdrdsdp);
        fillD2GDSDP   (r, s, t, p, d2gdsdp);
        fillD3GDS2DT  (r, s, t, p, d3gds2dt);

        order(SECOND | THIRD | FOURTH | SIXTH | SEVENTH | NINTH, t, p, r,
                    NULL, dsdr, dsdt, dsdp, NULL, d2sdrdt, d2sdrdp, NULL, d2sdtdp, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0,dxdt[i]=d3gdrdtdp[i]; j<NT; j++) {
                dxdt[i] += d3gdsdtdp[j]*dsdr[j][i] + d2gdsdt[j]*d2sdrdp[j][i] +
                   d3gdrdsdt[i][j]*dsdp[j] + d2gdrds[i][j]*d2sdtdp[j] +
                   d3gdrdsdp[i][j]*dsdt[j] + d2gdsdp[j]*d2sdrdt[j][i];
                for (k=0; k<NT; k++) {
                    dxdt[i] += d3gdrds2[i][j][k]*dsdt[j]*dsdp[k] +
                     d2gds2[j][k]*dsdt[j]*d2sdrdp[k][i] +
                     d2gds2[j][k]*dsdp[j]*d2sdrdt[k][i] +
                     d3gds2dt[j][k]*dsdr[j][i]*dsdp[k] +
                     d3gds2dp[j][k]*dsdr[j][i]*dsdt[k] +
                     d2gds2[j][k]*dsdr[j][i]*d2sdtdp[k];
                    for (l=0; l<NT; l++)
                        dxdt[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdt[k]*dsdp[l];
                }
            }
        }

        if (returnMixingProperties) for (i=0; i<NR; i++) dxdt[i] += (DVDT(0)-DVDT(i+1));
    }

    if(mask & TENTH) {
        double d3gds3[NT][NT][NT];
        double d3gdrds2[NR][NT][NT];
        double d3gds2dp[NT][NT];
        double d2gdrds[NR][NT];
        double d3gdrdsdp[NR][NT];
        double d3gdrdp2[NR];
        double d3gdsdp2[NT];
        double d2gdsdp[NT];
        double d2gds2[NT][NT];
        double dsdr[NT][NR], dsdp[NT], d2sdrdp[NT][NR], d2sdp2[NT];
        int i, j, k, l;

        fillD3GDS3    (r, s, t, p, d3gds3);
        fillD3GDRDS2  (r, s, t, p, d3gdrds2);
        fillD3GDS2DP  (r, s, t, p, d3gds2dp);
        fillD2GDRDS   (r, s, t, p, d2gdrds);
        fillD3GDRDSDP (r, s, t, p, d3gdrdsdp);
        fillD3GDRDP2  (r, s, t, p, d3gdrdp2);
        fillD3GDSDP2  (r, s, t, p, d3gdsdp2);
        fillD2GDSDP   (r, s, t, p, d2gdsdp);
        fillD2GDS2    (r, s, t, p, d2gds2);

        order(SECOND | FOURTH | SEVENTH | TENTH, t, p, r,
                    NULL, dsdr, NULL, dsdp, NULL, NULL, d2sdrdp, NULL, NULL, d2sdp2, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0,dxdp[i]=d3gdrdp2[i]; j<NT; j++) {
                dxdp[i] += d3gdsdp2[j]*dsdr[j][i] + d2gdsdp[j]*d2sdrdp[j][i] +
                   2.0*d3gdrdsdp[i][j]*dsdp[j] + d2gdrds[i][j]*d2sdp2[j] +
                   d2gdsdp[j]*d2sdrdp[j][i];
                for (k=0; k<NT; k++) {
                    dxdp[i] += d3gdrds2[i][j][k]*dsdp[j]*dsdp[k] +
                     2.0*d2gds2[j][k]*dsdp[j]*d2sdrdp[k][i] +
                     2.0*d3gds2dp[j][k]*dsdr[j][i]*dsdp[k] +
                     d2gds2[j][k]*dsdr[j][i]*d2sdp2[k];
                    for (l=0; l<NT; l++)
                        dxdp[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdp[k]*dsdp[l];
                }
            }
        }

        if (returnMixingProperties) for (i=0; i<NR; i++) dxdp[i] += (DVDP(0)-DVDP(i+1));
    }

    if(mask & ELEVENTH) {
        double d2gdsdw[NT][3*NP];
        double d2gdpdw[3*NP];
        double d2gds2[NT][NT];
        double d2gdsdp[NT];
        double dsdw[NT][3*NP], dsdp[NT];
        int i, k, l;

        fillD2GDSDW (r, s, t, p, d2gdsdw);
        fillD2GDPDW (r, s, t, p, d2gdpdw);
        fillD2GDS2  (r, s, t, p, d2gds2);
        fillD2GDSDP (r, s, t, p, d2gdsdp);

        order(FOURTH | ELEVENTH, t, p, r,
                    NULL, NULL, NULL, dsdp, NULL, NULL, NULL, NULL, NULL, NULL, dsdw, NULL);

        for (i=0; i<NP; i++) {
            if (modelParameters[i].activeH) {
                dw[i] = d2gdpdw[i];
                for (k=0; k<NT; k++) {
                    dw[i] += d2gdsdw[k][i]*dsdp[k] + d2gdsdp[k]*dsdw[k][i];
                    for (l=0; l<NT; l++) dw[i] += d2gds2[k][l]*dsdp[k]*dsdw[l][i] ;
                }
                dw[i] *= -1.0;
            } else dw[i] = 0.0;

            if (modelParameters[i].activeS) {
                dw[NP+i] = d2gdpdw[NP+i];
                for (k=0; k<NT; k++) {
                    dw[NP+i] += d2gdsdw[k][NP+i]*dsdp[k] + d2gdsdp[k]*dsdw[k][NP+i];
                    for (l=0; l<NT; l++) dw[NP+i] += d2gds2[k][l]*dsdp[k]*dsdw[l][NP+i] ;
                }
                dw[NP+i] *= -1.0;
            } else dw[NP+i] = 0.0;

            if (modelParameters[i].activeV) {
                dw[2*NP+i] = d2gdpdw[2*NP+i];
                for (k=0; k<NT; k++) {
                    dw[2*NP+i] += d2gdsdw[k][2*NP+i]*dsdp[k] + d2gdsdp[k]*dsdw[k][2*NP+i];
                    for (l=0; l<NT; l++) dw[2*NP+i] += d2gds2[k][l]*dsdp[k]*dsdw[l][2*NP+i] ;
                }
                dw[2*NP+i] *= -1.0;
            } else dw[2*NP+i] = 0.0;
        }

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            if (modelParameters[NW+0].activeS) dw[NP+NW+0] -= 1.0;
            for (i=0; i<NR; i++) {
                if (modelParameters[NW  +0].activeS) dw[NP+NW  +0] += r[i];
                if (modelParameters[NW+i+1].activeS) dw[NP+NW+i+1] -= r[i];
            }
        }

        for (i=0; i<(3*NP); i++) if (fabs(dw[i]) < 10.0*DBL_EPSILON) dw[i] = 0.0;
    }

#ifdef USE_GHIORSO_KRESS_MODEL
#undef  V
#undef  DVDT
#undef  DVDP
#undef  D2VDT2
#undef  D2VDTDP
#undef  D2VDP2
#define V(i)       0.0
#define DVDT(i)    0.0
#define DVDP(i)    0.0
#define D2VDT2(i)  0.0
#define D2VDTDP(i) 0.0
#define D2VDP2(i)  0.0
#endif  /* USE_GHIORSO_KRESS_MODEL */

    MTHREAD_MUTEX_UNLOCK(&global_data_mutex);
}

/* ============================================================================
   In the following public routine:
   m = m[i] = moles of the ith component in the liquid, and
   mu*O2    = mu O2 - mu0 O2, defined from the vector o[]
   This routine should be deprecated as it is inconsistent with conLiq
   ==========================================================================*/

void
muO2Liq(int mask, double t, double p, double *m,
    double *muO2, /* muO2      = mu*O2                 BINARY MASK: 0000000001 */
    double *dm,   /* dm[i]     = d mu*O2/dm[i]         BINARY MASK: 0000000010 */
    double *dt,   /* dt        = d mu*O2/d T           BINARY MASK: 0000000100 */
    double *dp,   /* dp        = d mu*O2/d P           BINARY MASK: 0000001000 */
    double **d2m, /* d2m[i][j] = d mu*O2/dm[i][j]      BINARY MASK: 0000010000 */
    double *d2mt, /* d2mt[i]   = d mu*O2/dm[i]dt       BINARY MASK: 0000100000 */
    double *d2mp, /* d2mp[i]   = d mu*O2/dm[i]dp       BINARY MASK: 0001000000 */
    double *d2t2, /* d2t2      = d mu*O2/dt2           BINARY MASK: 0010000000 */
    double *d2tp, /* d2tp      = d mu*O2/dtdp          BINARY MASK: 0100000000 */
    double *d2p2) /* d2p2      = d mu*O2/dp2           BINARY MASK: 1000000000 */
{
    /*-------------------------------------------------------------------------
        The algorithm used is that given by:
        Kress, VC and Carmichael, ISE (1991) The compressibility of silicate
            liquids containing Fe2O3 and the effect of composition, temperature,
            oxygen fugacity and pressure on their redox states.
            Contributions to Mineralogy and Petrology (in press)
        Coefficients for the oxides are initialized in LIQ_STRUCT_DATA.H
    --------------------------------------------------------------------------*/
    static int indexFeO = -1, indexFe2O3 = -1;
    static const double t0 = 1673.15,                      /* K       */
                       a =    0.196,
                       b =    1.1492e4,                  /* K       */
                       c =   -6.675,
                       e =   -3.364,
                       f =   -7.01e-7  * 1.0e5,          /* K/bar   */
                       g =   -1.54e-10 * 1.0e5,          /* 1/bar   */
                       h =    3.85e-17 * 1.0e5 * 1.0e5;  /* K/bar^2 */
    double mOx[NA], vTemp[NA], mTemp[NA][NA];
    int i, j;
    double sum;

    liqERRstate = ERR_NONE;

#ifdef PHMELTS_ADJUSTMENTS
    if ((silminState->fo2Alt) && (calculationMode != MODE_xMELTS)) {
        //    muO2Liq_Alt(mask, t, p, m, muO2, dm, dt, dp, d2m, d2mt, d2mp, d2t2, d2tp, d2p2);
        return;
    } else
#endif
    if ((calculationMode == MODE__MELTS) || (calculationMode == MODE_pMELTS)) {
        muO2Liq_v34(mask, t, p, m, muO2, dm, dt, dp, d2m, d2mt, d2mp, d2t2, d2tp, d2p2);
        return;
    } else if (calculationMode == MODE__MELTSandCO2) {
            muO2Liq_CO2(mask, t, p, m, muO2, dm, dt, dp, d2m, d2mt, d2mp, d2t2, d2tp, d2p2);
            return;
    } else if (calculationMode == MODE__MELTSandCO2_H2O) {
            muO2Liq_CO2_H2O(mask, t, p, m, muO2, dm, dt, dp, d2m, d2mt, d2mp, d2t2, d2tp, d2p2);
            return;
    }

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    MTHREAD_MUTEX_LOCK(&global_data_mutex);
    if ((indexFeO == -1) || (indexFe2O3 == -1)) {
        for (i=0; i<nc; i++) {
            if (bulkSystem[i].type == FEO)   indexFeO   = i;
            if (bulkSystem[i].type == FE2O3) indexFe2O3 = i;
        }
        if (indexFeO == 0 || indexFe2O3 == 0) {
       printf("Fatal error in muO2Liq (LIQUID.C)\n");
       printf("The oxides FeO and Fe2O3 cannot be identified.\n");
       return;
        }
    }
    MTHREAD_MUTEX_UNLOCK(&global_data_mutex);

    for (i=0; i<nc; i++)
        for (j=0, mOx[i]=0.0; j<nlc; j++) mOx[i] += (liquid[j].liqToOx)[i]*m[j];

    for (i=0, sum=0.0; i<nc; i++) { sum += mOx[i]; } sum += mOx[indexFe2O3];
    if (sum == 0.0 || mOx[indexFeO] == 0.0 || mOx[indexFe2O3] == 0.0) {
        if (mask & FIRST)   *muO2 = 0.0;
        if (mask & SECOND)  for (i=0; i<nlc; i++) dm[i] = 0.0;
        if (mask & THIRD)   *dt = 0.0;
        if (mask & FOURTH)  *dp = 0.0;
        if (mask & FIFTH)   {
            for (i=0; i<nlc; i++) for (j=0; j<nlc; j++) d2m[i][j] = 0.0;
        }
        if (mask & SIXTH)   for (i=0; i<nlc; i++) d2mt[i] = 0.0;
        if (mask & SEVENTH) for (i=0; i<nlc; i++) d2mp[i] = 0.0;
        if (mask & EIGHTH)  *d2t2 = 0.0;
        if (mask & NINTH)   *d2tp = 0.0;
        if (mask & TENTH)   *d2p2 = 0.0;
        return;
    }

    /*-------------------------------------------------------------------------*/

    if (mask & FIRST) {
        double temp;
        if (p < 50000.0) temp = b/t + c + e*(1.0 - t0/t - log(t/t0)) + f*p/t + g*(t-t0)*p/t + h*SQUARE(p)/t;
        else             temp = b/t + c + e*(1.0 - t0/t - log(t/t0)) + f*50000.0/t + g*(t-t0)*50000.0/t + h*SQUARE(50000.0)/t
                                                    - a*log(10.0)*(608.966*p/10000.0-608.966*5.0)/t;
        for (i=0; i<nc; i++) temp += bulkSystem[i].coeff*mOx[i]/sum;
        temp += 2.0*bulkSystem[indexFeO].coeff*mOx[indexFe2O3]/sum - bulkSystem[indexFe2O3].coeff*mOx[indexFe2O3]/sum;
        *muO2 = R*t*(log(mOx[indexFe2O3]/mOx[indexFeO]) - temp)/a;
    }

    if (mask & SECOND) {
        for (j=0; j<nc; j++) {
            double factor = (j == indexFe2O3) ? 2.0 : 1.0;
            for (i=0, vTemp[j]=0.0; i<nc; i++)
                vTemp[j] -= (i == j) ?
                    bulkSystem[i].coeff*(1.0-factor*mOx[i]/sum)/sum :
                    - bulkSystem[i].coeff*factor*mOx[i]/SQUARE(sum);
            vTemp[j] += - (factor*mOx[indexFe2O3]/SQUARE(sum))
                *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
            if      (j == indexFeO) vTemp[j] += -1.0/mOx[indexFeO];
            else if (j == indexFe2O3) {
                vTemp[j] += 1.0/mOx[indexFe2O3];
                vTemp[j] += (1.0/sum)
                    *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
            }
            vTemp[j] *= R*t/a;
        }
        for (i=0; i<nlc; i++)
            for (j=0, dm[i]=0.0; j<nc; j++) dm[i] += vTemp[j]*(liquid[i].liqToOx)[j];
    }

    if (mask & THIRD) {
        double temp;
        if (p < 50000.0) temp = b/t + c + e*(1.0 - t0/t - log(t/t0)) + f*p/t + g*(t-t0)*p/t + h*SQUARE(p)/t;
        else             temp = b/t + c + e*(1.0 - t0/t - log(t/t0)) + f*50000.0/t + g*(t-t0)*50000.0/t + h*SQUARE(50000.0)/t
                                                    - a*log(10.0)*(608.966*p/10000.0-608.966*5.0)/t;
        for (i=0; i<nc; i++) temp += bulkSystem[i].coeff*mOx[i]/sum;
        temp += 2.0*bulkSystem[indexFeO].coeff*mOx[indexFe2O3]/sum - bulkSystem[indexFe2O3].coeff*mOx[indexFe2O3]/sum;
        *dt = R*(log(mOx[indexFe2O3]/mOx[indexFeO]) - temp)/a + R*t*(b/SQUARE(t) - e*(t0/t-1.0)*(1.0/t))/a;
        if (p < 50000.0) *dt += R*t*(f*p/SQUARE(t) - g*(t0/t)*(p/t) + h*SQUARE(p/t))/a;
        else             *dt += R*t*(f*50000.0/SQUARE(t) - g*(t0/t)*(50000.0/t) + h*SQUARE(50000.0/t)
                                 - a*log(10.0)*(608.966*p/10000.0-608.966*5.0)/SQUARE(t))/a;
    }

    if (mask & FOURTH) {
        if (p < 50000.0) *dp = R*t*(-f/t - g*(t-t0)/t - 2.0*h*p/t)/a;
        else             *dp = R*t*(a*log(10.0)*608.966/10000.0/t)/a;
    }

    if (mask & FIFTH) {
        int k, l;
        for (k=0; k<nc; k++) {
            double factorK = (k == indexFe2O3) ? 2.0 : 1.0;
            for (j=0; j<nc; j++) {
                double factorJ = (j == indexFe2O3) ? 2.0 : 1.0;
                for (i=0, mTemp[k][j]=0.0; i<nc; i++) {
                    mTemp[k][j] -=
                        2.0*factorJ*factorK*bulkSystem[i].coeff*mOx[i]/CUBE(sum);
                    if (i == j) mTemp[k][j] -= - factorK*bulkSystem[i].coeff/SQUARE(sum);
                    if (i == k) mTemp[k][j] -= - factorJ*bulkSystem[i].coeff/SQUARE(sum);
                }
                mTemp[k][j] += 2.0*(factorJ*factorK*mOx[indexFe2O3]/CUBE(sum))
                    *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
                if      (j == indexFeO && k == indexFeO)
                    mTemp[k][j] += 1.0/SQUARE(mOx[indexFeO]);
                else if (j == indexFe2O3 && k == indexFe2O3) {
                    mTemp[k][j] += -1.0/SQUARE(mOx[indexFe2O3]);
                    mTemp[k][j] += -((factorJ+factorK)/SQUARE(sum))
                        *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
                }
                else if (j == indexFe2O3) mTemp[k][j] += -(factorK/SQUARE(sum))
                        *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
                else if (k == indexFe2O3) mTemp[k][j] += -(factorJ/SQUARE(sum))
                        *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
                mTemp[k][j] *= R*t/a;
            }
        }
        for (i=0; i<nlc; i++) for (j=0; j<nlc; j++)
            for (k=0, d2m[i][j]=0.0; k<nc; k++) for (l=0; l<nc; l++)
                d2m[i][j] += mTemp[k][l]*(liquid[i].liqToOx)[k]*(liquid[j].liqToOx)[l];
    }

    if (mask & SIXTH) {
        for (j=0; j<nc; j++) {
            double factor = (j == indexFe2O3) ? 2.0 : 1.0;
            for (i=0, vTemp[j]=0.0; i<nc; i++)
                vTemp[j] -= (i == j) ?
                    bulkSystem[i].coeff*(1.0-factor*mOx[i]/sum)/sum :
                    - bulkSystem[i].coeff*factor*mOx[i]/SQUARE(sum);
            vTemp[j] += - (factor*mOx[indexFe2O3]/SQUARE(sum))
                *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
            if      (j == indexFeO) vTemp[j] += -1.0/mOx[indexFeO];
            else if (j == indexFe2O3) {
                vTemp[j] += 1.0/mOx[indexFe2O3];
                vTemp[j] += (1.0/sum)
                    *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
            }
            vTemp[j] *= R/a;
        }
        for (i=0; i<nlc; i++) for (j=0, d2mt[i]=0.0; j<nc; j++)
            d2mt[i] += vTemp[j]*(liquid[i].liqToOx)[j];
    }

    if (mask & SEVENTH) {
        for (i=0; i<nlc; i++) d2mp[i] = 0.0;
    }

    if (mask & EIGHTH) {

        *d2t2 = 2.0*R*(b/SQUARE(t) - e*(t0/t-1.0)*(1.0/t))/a
                    + R*t*(-2.0*b/CUBE(t) - e*(1.0 - 2.0*t0/t)/SQUARE(t))/a;
    }

    if (mask & NINTH) {
        if (p < 50000.0) *d2tp = -R*(f/t + g*(t-t0)/t + 2.0*h*p/t)/a + R*t*(f/SQUARE(t) - g*(t0/t)/t + 2.0*h*p/SQUARE(t))/a;
        else             *d2tp = -R*(-a*log(10.0)*608.966/10000.0/t)/a + R*t*(-a*log(10.0)*608.966/10000.0/SQUARE(t))/a;
    }

    if (mask & TENTH) {
        if (p < 50000.0) *d2p2 = - (R*t*2.0*h/t)/a;
        else             *d2p2 = 0.0;
    }

}

void
visLiq(int mask, double t, double p, double *r,
    double *viscosity  /* log(10) viscosity            BINARY MASK: 00000001 */
    )
{
    double coeff[NA], factor[NA], m[NA], x[NA], sum;
    int nSiO2 = -1, i, j;

    liqERRstate = ERR_NONE;

    if ((calculationMode == MODE__MELTS) || (calculationMode == MODE_pMELTS)) {
            visLiq_v34(mask, t, p, r, viscosity);
            return;
    } else if (calculationMode == MODE__MELTSandCO2) {
            visLiq_CO2(mask, t, p, r, viscosity);
            return;
    } else if (calculationMode == MODE__MELTSandCO2_H2O) {
            visLiq_CO2_H2O(mask, t, p, r, viscosity);
            return;
    }

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    struct _shawModel {
        char   *oxide;
        double coeff;
        double factor;
    } shawModel[] = {
        { "TiO2",    4.5, 1.0 }, { "Al2O3",  6.7, 2.0 },
        { "Fe2O3",  3.4, 2.0 }, { "FeO",    3.4, 1.0 },
        { "MgO",    3.4, 1.0 }, { "CaO",    4.5, 1.0 },
        { "Na2O",    2.8, 1.0 }, { "K2O",    2.8, 1.0 },
        { "H2O",    2.0, 1.0 }
    };
    const int nShaw = (sizeof shawModel / sizeof(struct _shawModel));

    for (j=0; j<NA; j++) { coeff[j] = 0.0; factor[j] = 0.0; }
    for (i=0; i<nShaw; i++) {
        for (j=0; j<NA; j++) if (strcmp(shawModel[i].oxide, bulkSystem[j].label) == 0) {
            coeff[j]  = shawModel[i].coeff;
            factor[j] = shawModel[i].factor;
            break;
        }
    }
    for (i=0; i<NA; i++) if (strcmp("SiO2", bulkSystem[i].label) == 0) { nSiO2 = i; break; }

    if (nSiO2 == -1) { *viscosity = 0.0; return; }

    /* m[0] --> m[NA-1] is an array of mole fractions of liquid components      */
    for (m[0]=1.0, i=0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* convert m[] -> x[] : mole fractions of liquid comp -> moles of oxides    */
    for (i=0; i<NA; i++)
        for (j=0, x[i]=0.0; j<NA; j++) x[i] += (liquid[j].liqToOx)[i]*m[j];

    /* Convert to the Shaw mole fractions                                       */
    for (i=0, sum=0.0; i<NA; i++) {
        if (factor[i] > 0.0) x[i] *= factor[i];
        sum += x[i];
    }
    for (i=0; i<NA; i++) x[i] /= (sum != 0.0) ? sum : 1.0;

    if (mask & FIRST) {
        for (i=0, *viscosity=0.0; i<NA; i++) *viscosity += coeff[i]*x[nSiO2]*x[i];
        *viscosity /= (x[nSiO2] < 1.0) ? 1.0 - x[nSiO2] : 1.0;
        *viscosity  = (*viscosity)*(10000.0/t - 1.50)  - 6.40;
        *viscosity /= log(10.0);
    }
}

/* end of file LIQUID.C */
