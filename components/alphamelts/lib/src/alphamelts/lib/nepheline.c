const char *nepheline_ver(void) { return "$Id: nepheline.c,v 1.4 2007/03/12 20:06:36 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: nepheline.c,v $
MELTS Source Code: RCS Revision 1.4  2007/03/12 20:06:36  ghiorso
MELTS Source Code: RCS Changed fO2 inclusion criteria and adjusted some mineral endmember inclusion
MELTS Source Code: RCS tolerances.
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
MELTS Source Code: RCS Revision 1.1.1.1  2004/01/02 19:21:49  cvsaccount
MELTS Source Code: RCS CTserver University of Chicago
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:58:12  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
 * Revision 1.6  1997/06/21  22:49:37  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 1.5  1997/05/03  20:23:13  ghiorso
 * *** empty log message ***
 *
 * Revision 1.4  1997/03/27  17:03:21  ghiorso
 * *** empty log message ***
 *
 * Revision 1.3  1996/09/24  20:33:28  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 1.2  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 1.1  1995/11/01  22:40:27  ghiorso
 * Initial revision
 *
 * Revision 1.1  1995/11/01  22:40:27  ghiorso
 * Initial revision
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Routines to compute nepheline solution properties
**      (file: NEPHELINE.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso October 14, 1995 Original Version - started
**              Modified from SPINEL.C
**--
*/

#include "melts_gsl.h"
#include "silmin.h"  /* Structure definitions foor SILMIN package */

#ifdef DEBUG
#undef DEBUG
#endif

#define SQUARE(x)  ((x)*(x))
#define CUBE(x)    ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

/*
 *=============================================================================
 * Nepheline solution parameters:
 * Sack, R.O., Ghiorso, M.S. (1995)
 * Thermodynamics of Nepheline Solutions
 */

                                        /* kals - neph */
#define DH1                  32145.67 /* Na4Al4Si4O16 joules */
#define DS1                     10.46 /* joules/T */
#define DV1                       0.0 /* joules/bar */
#define DH2                   -5648.4 /* K4Al4Si4O16 joules */
#define DS2                       0.0 /* joules/T */
#define DV2                     -0.04 /* joules/bar */
#define DH3                   14644.0 /* []Na3Al3Si5O16 joules */
#define DS3                     -18.7 /* joules/T */
#define DV3                       0.0 /* joules/bar */
#define DH4                   35184.0 /* []CaNa2Al2Si2O16 joules */
#define DS4                     -7.17 /* joules/T */
#define DV4                       0.0 /* joules/bar */

#define HEX                 -31744.63 /* joules */
#define SEX                    -20.92 /* joules/T */
#define VEX                     -1.05 /* joules/bar */
#define HX                  -13893.18 /* joules */
#define SX                      12.55 /* joules/T */
#define VX                        0.0 /* joules/bar */

#define H23                   73520.0 /* joules */
#define H24                   42560.0 /* joules */
#define S23                       0.0 /* joules/K */
#define S24                       0.0 /* joules/K */
                                        /* nepheline structure */
#define WHNAKLS               6861.76 /* W(LS)Na-K joules */
#define WVNAKLS                  0.33 /* joules/bar */
#define WHNAKSS              51002.96 /* W(SS)Na-K joules */
#define WVNAKSS                  0.54 /* joules/bar */
#define WVN                   14644.0 /* W[]Si-NaAl joules */
#define WVK                    8368.0 /* W[]Si-KAl joules */
#define WCANA                     0.0 /* W[]Ca-Na2 joules */
#define WCAK                      0.0 /* W[]Ca-K2 joules */
#define WPLAG                     0.0 /* WCaAl-NaSi joules */
                                        /* kalsilite structure */
#define WKALS                 29288.0 /* WNa4-K4 joules */
#define WKVKALS               30334.0 /* WK4-[]Na3 joules */
#define WNVKALS               14646.0 /* WNa4-[]Na3 joules */
#define WNCKALS                   0.0 /* WNa4-[]CaNa2 joules */
#define WKCKALS               14644.0 /* WK4-[]CaNa2 joules */
#define WVVKALS                   0.0 /* W[]Na3-[]CaNa2 joules */

#define DWKNALS                   0.0 /* joules */
#define DWKNASS                   0.0 /* joules */
#define DWVNASS                   0.0 /* joules */
#define DWVNALS                   0.0 /* joules */
#define DWVKLS                    0.0 /* joules */
#define AWLS                      0.0 /* joules */
#define AWSS                      0.0 /* joules */
#define AWVNA                     0.0 /* joules */
#define AWVK                      0.0 /* joules */

/*
 * Global (to this file): variables
 */

#define R  8.3143
#define NR         3    /* Three independent composition variables */
#define NS         1    /* One ordering parameter                  */
#define NA         4    /* Four endmember compositions             */
                   /* Penalty constant for approaching the zero
                                            concentration of r[1] - see #define G        */
#define PENALTY    sqrt(DBL_EPSILON)

                                /* correction to zero point entropy of []CaNa2Al4Si4O16 */
#define S4      -15.8765
#define G4      -t*(S4)

#define GEX     (HEX)-t*(SEX)+(p-1.0)*(VEX)
#define GX      (HX)-t*(SX)+(p-1.0)*(VX)
#define WNAKLS  (WHNAKLS)+(p-1.0)*(WVNAKLS)
#define WNAKSS  (WHNAKSS)+(p-1.0)*(WVNAKSS)
#define G23     (H23)-t*(S23)
#define G24     (H24)-t*(S24)

#define MAX_ITER 200    /* Maximum number of iterations allowed in order */

/*************************************/
/* Statics for Ordering Calculations */
/*************************************/

static MTHREAD_ONCE_T initThreadOBlock = MTHREAD_ONCE_INIT;

static MTHREAD_KEY_T tOldKey;
static MTHREAD_KEY_T pOldKey;
static MTHREAD_KEY_T rOldKey;
static MTHREAD_KEY_T sOldKey;
static MTHREAD_KEY_T ptToD2gds2Key;
static MTHREAD_KEY_T d2gds2Key;
static MTHREAD_KEY_T indexD2gds2Key;

static void freeNSarray(void *NSarray) {
    gsl_vector_free((gsl_vector *) NSarray);
}

static void freePtToD2gds2(void *ptToD2gds2) {
    gsl_matrix_free((gsl_matrix *) ptToD2gds2);
}

static void freeIndexD2gds2(void *indexD2gds2) {
    gsl_permutation_free((gsl_permutation *) indexD2gds2);
}

static void threadOInit(void) {
    MTHREAD_KEY_CREATE(&tOldKey,       free);
    MTHREAD_KEY_CREATE(&pOldKey,       free);
    MTHREAD_KEY_CREATE(&rOldKey,       freeNSarray);
    MTHREAD_KEY_CREATE(&sOldKey,       freeNSarray);
    MTHREAD_KEY_CREATE(&ptToD2gds2Key, freePtToD2gds2);
    MTHREAD_KEY_CREATE(&d2gds2Key,     free);
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
        sOldPt = gsl_vector_alloc((size_t) NS);
        gsl_vector_set_all(sOldPt, 2.0);
        MTHREAD_SETSPECIFIC(sOldKey, (void *) sOldPt);
    }
    return sOldPt->data;
}

static gsl_matrix *getPtToD2gds2() {
    gsl_matrix *ptToD2gds2Pt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    ptToD2gds2Pt = (gsl_matrix *) MTHREAD_GETSPECIFIC(ptToD2gds2Key);
    if (ptToD2gds2Pt == NULL) {
        ptToD2gds2Pt  = gsl_matrix_alloc((size_t) NS, (size_t) NS);
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
        d2gds2Pt  = (double **) malloc((size_t) NS*sizeof(double *));
        for (i=0; i<NS; i++) d2gds2Pt[i] = &(ptToD2gds2Pt->data[i*NS]);
        MTHREAD_SETSPECIFIC(d2gds2Key, (void *) d2gds2Pt);
    }
    return d2gds2Pt;
}

static gsl_permutation *getIndexD2gds2() {
    gsl_permutation *indexD2gds2Pt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    indexD2gds2Pt = (gsl_permutation *) MTHREAD_GETSPECIFIC(indexD2gds2Key);
    if (indexD2gds2Pt == NULL) {
        indexD2gds2Pt = gsl_permutation_alloc((size_t) NS);
        gsl_permutation_init(indexD2gds2Pt);
        MTHREAD_SETSPECIFIC(indexD2gds2Key, (void *) indexD2gds2Pt);
    }
    return indexD2gds2Pt;
}

/***********************************/
/* Statics for Site Mole Fractions */
/***********************************/

#define DECLARE_SITE_FRACTIONS \
    double xkls, xvcls, xnals, xkss, xnass, xcass;

#define XKLS  0
#define XVCLS 1
#define XNALS 2
#define XKSS  3
#define XNASS 4
#define XCASS 5

#define NX    6

static MTHREAD_ONCE_T initThreadXBlock = MTHREAD_ONCE_INIT;

static MTHREAD_KEY_T xKey[NX];

static void threadXInit(void) {
    int i;
    for (i=0; i<NX; i++) MTHREAD_KEY_CREATE(&xKey[i], free);
}

static double getX(int n) {
    double *xPt;
    MTHREAD_ONCE(&initThreadXBlock, threadXInit);

    xPt = (double *) MTHREAD_GETSPECIFIC(xKey[n]);
    if (xPt == NULL) {
        xPt = (double *) malloc(sizeof(double)); *xPt = 0.0;
        MTHREAD_SETSPECIFIC(xKey[n], (void *) xPt);
    }
    return *xPt;
}

static void setX(int n, double x) {
    double *xPt;
    MTHREAD_ONCE(&initThreadXBlock, threadXInit);

    xPt = (double *) MTHREAD_GETSPECIFIC(xKey[n]);
    if (xPt == NULL) {
        xPt = (double *) malloc(sizeof(double));
        MTHREAD_SETSPECIFIC(xKey[n], (void *) xPt);
    }
    *xPt = x;
}

#define GET_SITE_FRACTIONS \
    xkls  = getX(XKLS);  \
    xvcls = getX(XVCLS); \
    xnals = getX(XNALS); \
    xkss  = getX(XKSS);  \
    xnass = getX(XNASS); \
    xcass = getX(XCASS);

#define SET_SITE_FRACTIONS \
    setX(XKLS,  xkls ); \
    setX(XVCLS, xvcls); \
    setX(XNALS, xnals); \
    setX(XKSS,  xkss ); \
    setX(XNASS, xnass); \
    setX(XCASS, xcass);

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conSpn defines the conversion from m[i], to r[j]
 */
                   /* Order: X2, X3, X4 */
#define FR2(i)     (i == 1) ? 1.0 - r[0] : - r[0]
#define FR3(i)     (i == 2) ? 1.0 - r[1] : - r[1]
#define FR4(i)     (i == 3) ? 1.0 - r[2] : - r[2]

                   /* Order: S1 */
#define GS1(i)     - s[0]

#define DFR2DR2(i) - 1.0
#define DFR3DR3(i) - 1.0
#define DFR4DR4(i) - 1.0

#define DGS1DS1(i) - 1.0

/*
 * Global (to this file): derivative definitions
 */

#define S (S4)*r[2] - R*(xkls*log(xkls) + xvcls*log(xvcls) + \
                    xnals*log(xnals) - (1.0-3.0*xcass)*log(1.0-3.0*xcass) + \
                    3.0*xkss*log(xkss) + 3.0*xcass*log(xcass) + 3.0*xnass*log(xnass)) \
                    + 0.25*(2.0*(SEX)+(SX))*s[0] + (SX)*r[0]*(1.0-r[0]) \
                    + (3.0/16.0)*(SX)*s[0]*s[0] - 0.5*(SX)*r[0]*s[0] \
                    + 0.5*(2.0*(S23)+(SEX)-(SX))*r[0]*r[1] + 0.125*((SX)-(SEX)-2.0*(S23))*r[1]*s[0] \
                    + (6.0*(S23)+3.0*(SEX)-15.0*(SX))*r[0]*r[2]/18.0 \
                    + (-9.0*(SEX)-3.0*(SX)+6.0*(S23)-36.0*(S24))*r[2]*s[0]/24.0
#define H 0.25*(2.0*(HEX)+(HX)+3.0*(WHNAKLS)-(WHNAKSS))*s[0] \
                    + ((HX)+(WHNAKLS)+(WHNAKSS))*r[0]*(1.0-r[0]) \
                    + (1.0/16.0)*(3.0*(HX)-9.0*(WHNAKLS)-(WHNAKSS))*s[0]*s[0] \
                    - 0.5*((HX)+3.0*(WHNAKLS)-(WHNAKSS))*r[0]*s[0] + (WVN)*r[1]*(1.0-r[1]) \
                    + 0.5*(2.0*(H23)+(HEX)-(HX)-2.0*(WHNAKLS)-2.0*(WVN)+2.0*(WVK))*r[0]*r[1] \
                    + 0.125*((HX)-(HEX)-2.0*(H23)-6.0*(WHNAKLS)-6.0*(WVN)+6.0*(WVK))*r[1]*s[0] \
                    + (DWKNALS)*(r[0]*r[0]+0.5*r[0]*s[0]-3.0*s[0]*s[0]/16.0-r[0]*r[0]*r[0] \
                       -5.0*r[0]*r[0]*s[0]/4.0-3.0*r[0]*s[0]*s[0]/16.0 \
                       +9.0*s[0]*s[0]*s[0]/64.0-r[0]*r[1]/3.0-r[1]*s[0]/4.0 \
                       -2.0*r[0]*r[0]*r[1]/3.0+r[0]*r[1]*r[1]/3.0+r[1]*r[1]*s[0]/4.0 \
                       +3.0*r[1]*s[0]*s[0]/8.0) \
                    + (DWKNASS)*(r[0]*r[0]+0.5*r[0]*s[0]-3.0*s[0]*s[0]/16.0-r[0]*r[0]*r[0] \
                       -r[0]*r[0]*s[0]/4.0+5.0*r[0]*s[0]*s[0]/16.0 \
                       -3.0*s[0]*s[0]*s[0]/64.0) \
                    + (DWVNASS)*(4.0*r[0]*r[1]/3.0-4.0*r[0]*r[0]*r[1]/3.0-r[0]*r[1]*r[1]/3.0 \
                       -r[1]*r[1]*s[0]/4.0-r[1]*s[0]*s[0]/4.0) \
                    + (DWVNALS)*(2.0*r[0]*r[1]/3.0-r[1]*s[0]/2.0-2.0*r[0]*r[0]*r[1]/3.0 \
                       -2.0*r[0]*r[1]*r[1]/3.0+0.5*r[1]*r[1]*s[0] \
                       +3.0*r[1]*s[0]*s[0]/8.0) \
                    + (DWVKLS)*(r[0]*r[1]/3.0+r[1]*s[0]/4.0+2.0*r[0]*r[0]*r[1]/3.0 \
                                            -r[0]*r[1]*r[1]/3.0-r[1]*r[1]*s[0]/4.0 \
                                            -3.0*r[1]*s[0]*s[0]/8.0) \
                    + (AWLS)*(r[0]+3.0*s[0]/4.0-3.0*r[0]*r[0]-9.0*r[0]*s[0]/2.0 \
                                        -27.0*s[0]*s[0]/16.0+2.0*r[0]*r[0]*r[0]+9.0*r[0]*r[0]*s[0]/2.0 \
                                        +27.0*r[0]*s[0]*s[0]/8.0+27.0*s[0]*s[0]*s[0]/32.0 + r[0]*r[1] \
                                        +3.0*r[1]*s[0]/4.0-2.0*r[0]*r[1]*r[1]-3.0*r[1]*r[1]*s[0]/2.0) \
                    + (AWSS)*(r[0]-s[0]/4.0-3.0*r[0]*r[0]+3.0*r[0]*s[0]/2.0-3.0*s[0]*s[0]/16.0 \
                                        +2.0*r[0]*r[0]*r[0]-3.0*r[0]*r[0]*s[0]/2.0+3.0*r[0]*s[0]*s[0]/8.0 \
                                        -s[0]*s[0]*s[0]/32.0) \
                    + (AWVNA)*(r[1]-r[0]*r[1]-3.0*r[1]*s[0]/4.0-3.0*r[1]*r[1]+2.0*r[0]*r[1]*r[1] \
                     + 3.0*r[1]*r[1]*s[0]/2.0+2.0*r[1]*r[1]*r[1]) \
                    + (AWVK)*(r[0]*r[1]+3.0*r[1]*s[0]/4.0-2.0*r[0]*r[1]*r[1] \
                                        -3.0*r[1]*r[1]*s[0]/2.0) \
                    + (WCANA)*r[2]*(1.0-r[2]) \
                    + (6.0*(H23)+3.0*(HEX)-15.0*(HX)-18.0*(WHNAKLS)-4.0*(WHNAKSS)+36.0*(WCAK) \
             -36.0*(WCANA)+18.0*(WVN)-18.0*(WVK))*r[0]*r[2]/18.0 \
                    + ((WPLAG)-(WCANA)-(WVN))*r[1]*r[2] \
                    + (18.0*(WVN)-18.0*(WVK)+36.0*(WCAK)-36.0*(WCANA)-18.0*(WHNAKLS) \
             +4.0*(WHNAKSS)-9.0*(HEX)-3.0*(HX)+6.0*(H23)-36.0*(H24))*r[2]*s[0]/24.0
#define V 0.25*(2.0*(VEX)+(VX)+3.0*(WVNAKLS)-(WVNAKSS))*s[0] \
                    + ((VX)+(WVNAKLS)+(WVNAKSS))*r[0]*(1.0-r[0]) \
                    + (1.0/16.0)*(3.0*(VX)-9.0*(WVNAKLS)-(WVNAKSS))*s[0]*s[0] \
                    - 0.5*((VX)+3.0*(WVNAKLS)-(WVNAKSS))*r[0]*s[0] \
                    + 0.5*((VEX)-(VX)-2.0*(WVNAKLS))*r[0]*r[1] \
                    + 0.125*((VX)-(VEX)-6.0*(WVNAKLS))*r[1]*s[0] \
                    + (3.0*(VEX)-15.0*(VX)-18.0*(WVNAKLS)-4.0*(WVNAKSS))*r[0]*r[2]/18.0 \
                    + (-18.0*(WVNAKLS)+4.0*(WVNAKSS)-9.0*(VEX)-3.0*(VX))*r[2]*s[0]/24.0
#define G (H) - t*(S) + (p-1.0)*(V) + (PENALTY)/r[1]

/*----------------------------------------------------------------------------*/

#define DGDR0 R*t*(log(xkls) - log(xnals) + 3.0*log(xkss) - 3.0*log(xnass)) \
                            + ((GX)+(WNAKLS)+(WNAKSS))*(1.0-2.0*r[0]) \
                            - 0.5*((GX)+3.0*(WNAKLS)-(WNAKSS))*s[0] \
                            + 0.5*(2.0*(G23)+(GEX)-(GX)-2.0*(WNAKLS)-2.0*(WVN)+2.0*(WVK))*r[1] \
                            + (DWKNALS)*(2.0*r[0]+0.5*s[0]-3.0*r[0]*r[0]-5.0*r[0]*s[0]/2.0 \
                           -3.0*s[0]*s[0]/16.0-r[1]/3.0-4.0*r[0]*r[1]/3.0 \
                           +r[1]*r[1]/3.0+3.0*r[1]*s[0]*s[0]/8.0) \
                            + (DWKNASS)*(2.0*r[0]+0.5*s[0]-3.0*r[0]*r[0]-r[0]*s[0]/2.0 \
                           +5.0*s[0]*s[0]/16.0) \
                            + (DWVNASS)*(4.0*r[1]/3.0-8.0*r[0]*r[1]/3.0-r[1]*r[1]/3.0) \
                            + (DWVNALS)*(2.0*r[1]/3.0-4.0*r[0]*r[1]/3.0-2.0*r[1]*r[1]/3.0) \
                            + (DWVKLS)*(r[1]/3.0+4.0*r[0]*r[1]/3.0-r[1]*r[1]/3.0) \
                            + (AWLS)*(1.0-6.0*r[0]-9.0*s[0]/2.0+6.0*r[0]*r[0]+9.0*r[0]*s[0] \
                                                +27.0*s[0]*s[0]/8.0+r[1]-2.0*r[1]*r[1]) \
                            + (AWSS)*(1.0-6.0*r[0]+3.0*s[0]/2.0+6.0*r[0]*r[0] \
                                                -6.0*r[0]*s[0]/2.0+3.0*s[0]*s[0]/8.0) \
                            + (AWVNA)*(-r[1]+2.0*r[1]*r[1]) + (AWVK)*(r[1]-2.0*r[1]*r[1]) \
                            + (6.0*(G23)+3.0*(GEX)-15.0*(GX)-18.0*(WNAKLS)-4.0*(WNAKSS)+36.0*(WCAK) \
                 -36.0*(WCANA)+18.0*(WVN)-18.0*(WVK))*r[2]/18.0
#define DGDR1 R*t*(log(xvcls) - log(xnals)) + (WVN)*(1.0-2.0*r[1]) \
                            + 0.5*(2.0*(G23)+(GEX)-(GX)-2.0*(WNAKLS)-2.0*(WVN)+2.0*(WVK))*r[0] \
                            + 0.125*((GX)-(GEX)-2.0*(G23)-6.0*(WNAKLS)-6.0*(WVN)+6.0*(WVK))*s[0] \
                            + (DWKNALS)*(-r[0]/3.0-s[0]/4.0-2.0*r[0]*r[0]/3.0+2.0*r[0]*r[1]/3.0 \
                           +r[1]*s[0]/2.0+3.0*s[0]*s[0]/8.0) \
                            + (DWVNASS)*(4.0*r[0]/3.0-4.0*r[0]*r[0]/3.0-2.0*r[0]*r[1]/3.0 \
                           -r[1]*s[0]/2.0-s[0]*s[0]/4.0) \
                            + (DWVNALS)*(2.0*r[0]/3.0-s[0]/2.0-2.0*r[0]*r[0]/3.0 \
                           -4.0*r[0]*r[1]/3.0+r[1]*s[0]+3.0*s[0]*s[0]/8.0) \
                            + (DWVKLS)*(r[0]/3.0+s[0]/4.0+2.0*r[0]*r[0]/3.0 \
                                            -2.0*r[0]*r[1]/3.0-r[1]*s[0]/2.0-3.0*s[0]*s[0]/8.0) \
                            + (AWLS)*(r[0]+3.0*s[0]/4.0-4.0*r[0]*r[1]-6.0*r[1]*s[0]/2.0) \
                            + (AWVNA)*(1.0-r[0]-3.0*s[0]/4.0-6.0*r[1]+4.0*r[0]*r[1] \
                         + 6.0*r[1]*s[0]/2.0+6.0*r[1]*r[1]) \
                            + (AWVK)*(r[0]+3.0*s[0]/4.0-4.0*r[0]*r[1]-6.0*r[1]*s[0]/2.0) \
                            + ((WPLAG)-(WCANA)-(WVN))*r[2] \
                            - (PENALTY)/(r[1]*r[1])
#define DGDR2 (G4) + R*t*(log(xvcls) - log(xnals) + log(1.0-3.0*xcass) \
                       + log(xcass) - log(xnass) + 1.0) \
                            + (WCANA)*(1.0-2.0*r[2]) \
                            + (6.0*(G23)+3.0*(GEX)-15.0*(GX)-18.0*(WNAKLS)-4.0*(WNAKSS)+36.0*(WCAK) \
                 -36.0*(WCANA)+18.0*(WVN)-18.0*(WVK))*r[0]/18.0 \
                            + ((WPLAG)-(WCANA)-(WVN))*r[1] \
                            + (18.0*(WVN)-18.0*(WVK)+36.0*(WCAK)-36.0*(WCANA)-18.0*(WNAKLS) \
                 +4.0*(WNAKSS)-9.0*(GEX)-3.0*(GX)+6.0*(G23)-36.0*(G24))*s[0]/24.0
#define DGDS0 R*t*(0.75*log(xkls) - 0.75*log(xnals) \
                            - 3.0*0.25*log(xkss) + 3.0*0.25*log(xnass)) \
                            + 0.25*(2.0*(GEX)+(GX)+3.0*(WNAKLS)-(WNAKSS)) \
                            + (1.0/8.0)*(3.0*(GX)-9.0*(WNAKLS)-(WNAKSS))*s[0] \
                            - 0.5*((GX)+3.0*(WNAKLS)-(WNAKSS))*r[0] \
                            + 0.125*((GX)-(GEX)-2.0*(G23)-6.0*(WNAKLS)-6.0*(WVN)+6.0*(WVK))*r[1] \
                            + (DWKNALS)*(0.5*r[0]-3.0*s[0]/8.0-5.0*r[0]*r[0]/4.0 \
                           -3.0*r[0]*s[0]/8.0+27.0*s[0]*s[0]/64.0-r[1]/4.0 \
                           +r[1]*r[1]/4.0+3.0*r[1]*s[0]/4.0) \
                            + (DWKNASS)*(0.5*r[0]-3.0*s[0]/8.0-r[0]*r[0]/4.0+5.0*r[0]*s[0]/8.0 \
                           -9.0*s[0]*s[0]/64.0) \
                            + (DWVNASS)*(-r[1]*r[1]/4.0-r[1]*s[0]/2.0) \
                            + (DWVNALS)*(-r[1]/2.0+0.5*r[1]*r[1]+3.0*r[1]*s[0]/4.0) \
                            + (DWVKLS)*(r[1]/4.0-r[1]*r[1]/4.0-3.0*r[1]*s[0]/4.0) \
                            + (AWLS)*(3.0/4.0-9.0*r[0]/2.0-27.0*s[0]/8.0+9.0*r[0]*r[0]/2.0 \
                                                +27.0*r[0]*s[0]/4.0+81.0*s[0]*s[0]/32.0+3.0*r[1]/4.0 \
                                                -3.0*r[1]*r[1]/2.0) \
                            + (AWSS)*(-1.0/4.0+3.0*r[0]/2.0-3.0*s[0]/8.0-3.0*r[0]*r[0]/2.0 \
                                                +3.0*r[0]*s[0]/4.0-3.0*s[0]*s[0]/32.0) \
                            + (AWVNA)*(-3.0*r[1]/4.0+3.0*r[1]*r[1]/2.0) \
                            + (AWVK)*(3.0*r[1]/4.0-3.0*r[1]*r[1]/2.0) \
                            + (18.0*(WVN)-18.0*(WVK)+36.0*(WCAK)-36.0*(WCANA)-18.0*(WNAKLS) \
                 +4.0*(WNAKSS)-9.0*(GEX)-3.0*(GX)+6.0*(G23)-36.0*(G24))*r[2]/24.0
#define DGDT  - (S)
#define DGDP  (V)

/*----------------------------------------------------------------------------*/

#define D2GDR0R0 R*t*(1.0/xkls + 1.0/xnals + 3.0/xkss + 3.0/xnass) \
                 - 2.0*((GX)+(WNAKLS)+(WNAKSS)) \
                 + (DWKNALS)*(2.0-6.0*r[0]-5.0*s[0]/2.0-4.0*r[1]/3.0) \
                 + (DWKNASS)*(2.0-6.0*r[0]-s[0]/2.0) + (DWVNASS)*(-8.0*r[1]/3.0) \
                 - (DWVNALS)*(4.0*r[1]/3.0) + (DWVKLS)*(4.0*r[1]/3.0) \
                 + (AWLS)*(-6.0+12.0*r[0]+9.0*s[0]) \
                 + (AWSS)*(-6.0+12.0*r[0]-6.0*s[0]/2.0)
#define D2GDR0R1 R*t/xnals \
                 + 0.5*(2.0*(G23)+(GEX)-(GX)-2.0*(WNAKLS)-2.0*(WVN)+2.0*(WVK)) \
                 + (DWKNALS)*(-1.0/3.0-4.0*r[0]/3.0+2.0*r[1]/3.0+3.0*s[0]*s[0]/8.0) \
                 + (DWVNASS)*(4.0/3.0-8.0*r[0]/3.0-2.0*r[1]/3.0) \
                 + (DWVNALS)*(2.0/3.0-4.0*r[0]/3.0-4.0*r[1]/3.0) \
                 + (DWVKLS)*(1.0/3.0+4.0*r[0]/3.0-2.0*r[1]/3.0) \
                 + (AWLS)*(1.0-4.0*r[1]) + (AWVNA)*(-1.0+4.0*r[1]) \
                 + (AWVK)*(1.0-4.0*r[1])
#define D2GDR0R2 R*t*(1.0/xnals + 1.0/xnass) \
                 + (6.0*(G23)+3.0*(GEX)-15.0*(GX)-18.0*(WNAKLS)-4.0*(WNAKSS)+36.0*(WCAK) \
                                        -36.0*(WCANA)+18.0*(WVN)-18.0*(WVK))/18.0
#define D2GDR0S0 R*t*(0.75/xkls + 0.75/xnals - 3.0*0.25/xkss \
                                                                - 3.0*0.25/xnass) \
                 - 0.5*((GX)+3.0*(WNAKLS)-(WNAKSS)) \
                 + (DWKNALS)*(0.5-5.0*r[0]/2.0-6.0*s[0]/16.0+6.0*r[1]*s[0]/8.0) \
                 + (DWKNASS)*(0.5-r[0]/2.0+10.0*s[0]/16.0) \
                 + (AWLS)*(-9.0/2.0+9.0*r[0]+27.0*s[0]/4.0+r[1]) \
                 + (AWSS)*(3.0/2.0-6.0*r[0]/2.0+3.0*s[0]/4.0)
#define D2GDR0DT R*(log(xkls) - log(xnals) + 3.0*log(xkss) \
                                                    - 3.0*log(xnass)) - (SX)*(1.0-2.0*r[0]) \
                 + 0.5*(SX)*s[0] - 0.5*(2.0*(S23)+(SEX)-(SX))*r[1] \
                 + (-6.0*(S23)-3.0*(SEX)+15.0*(SX))*r[2]/18.0
#define D2GDR0DP ((VX)+(WVNAKLS)+(WVNAKSS))*(1.0-2.0*r[0]) \
                 - 0.5*((VX)+3.0*(WVNAKLS)-(WVNAKSS))*s[0] \
                 + 0.5*((VEX)-(VX)-2.0*(WVNAKLS))*r[1] \
                 + (3.0*(VEX)-15.0*(VX)-18.0*(WVNAKLS)-4.0*(WVNAKSS))*r[2]/18.0

#define D2GDR1R1 R*t*(1.0/xvcls + 1.0/xnals) - 2.0*(WVN) \
                 + (DWKNALS)*(2.0*r[0]/3.0+s[0]/2.0) \
                 + (DWVNASS)*(-2.0*r[0]/3.0-s[0]/2.0) \
                 + (DWVNALS)*(-4.0*r[0]/3.0+s[0]) \
                 + (DWVKLS)*(-2.0*r[0]/3.0-s[0]/2.0) \
                 + (AWLS)*(-4.0*r[0]-6.0*s[0]/2.0) \
                 + (AWVNA)*(-6.0+4.0*r[0]+6.0*s[0]/2.0+12.0*r[1]) \
                 + (AWVK)*(-4.0*r[0]-6.0*s[0]/2.0) \
                 + 2.0*(PENALTY)/(r[1]*r[1]*r[1])
#define D2GDR1R2 R*t*(1.0/xvcls + 1.0/xnals) + ((WPLAG)-(WCANA)-(WVN))
#define D2GDR1S0 R*t*0.75/xnals \
                 + 0.125*((GX)-(GEX)-2.0*(G23)-6.0*(WNAKLS)-6.0*(WVN)+6.0*(WVK)) \
                 + (DWKNALS)*(-1.0/4.0+r[1]/2.0+3.0*s[0]/4.0) \
                 + (DWVNASS)*(-r[1]/2.0-s[0]/2.0) \
                 + (DWVNALS)*(-1.0/2.0+r[1]+3.0*s[0]/4.0) \
                 + (DWVKLS)*(1.0/4.0-r[1]/2.0-3.0*s[0]/4.0) \
                 + (AWLS)*(3.0/4.0-6.0*r[1]/2.0) \
                 + (AWVNA)*(-3.0/4.0+6.0*r[1]/2.0) \
                 + (AWVK)*(3.0/4.0-6.0*r[1]/2.0)
#define D2GDR1DT R*(log(xvcls) - log(xnals)) \
                 - 0.5*(2.0*(S23)+(SEX)-(SX))*r[0] - 0.125*((SX)-(SEX)-2.0*(S23))*s[0]
#define D2GDR1DP 0.5*((VEX)-(VX)-2.0*(WVNAKLS))*r[0] \
                 + 0.125*((VX)-(VEX)-6.0*(WVNAKLS))*s[0]

#define D2GDR2R2 R*t*(1.0/xvcls + 1.0/xnals - 1.0/(1.0-3.0*xcass) \
                       + 1.0/(3.0*xcass) + 1.0/(3.0*xnass)) - 2.0*(WCANA)
#define D2GDR2S0 R*t*(0.75/xnals - 0.25/xnass) \
                 + (18.0*(WVN)-18.0*(WVK)+36.0*(WCAK)-36.0*(WCANA)-18.0*(WNAKLS) \
                                        +4.0*(WNAKSS)-9.0*(GEX)-3.0*(GX)+6.0*(G23)-36.0*(G24))/24.0
#define D2GDR2DT - (S4) + R*(log(xvcls) - log(xnals) + log(1.0-3.0*xcass) \
                             + log(xcass) - log(xnass) + 1.0) \
                 + (-6.0*(S23)-3.0*(SEX)+15.0*(SX))*r[0]/18.0 \
                 + (9.0*(SEX)+3.0*(SX)-6.0*(S23)+36.0*(S24))*s[0]/24.0
#define D2GDR2DP (3.0*(VEX)-15.0*(VX)-18.0*(WVNAKLS)-4.0*(WVNAKSS))*r[0]/18.0 \
                 + (-18.0*(WVNAKLS)+4.0*(WVNAKSS)-9.0*(VEX)-3.0*(VX))*s[0]/24.0

#define D2GDS0S0 R*t*(0.75*0.75/xkls + 0.75*0.75/xnals \
                 + 3.0*0.25*0.25/xkss + 3.0*0.25*0.25/xnass) \
                 + (1.0/8.0)*(3.0*(GX)-9.0*(WNAKLS)-(WNAKSS)) \
                 + (DWKNALS)*(-3.0/8.0-3.0*r[0]/8.0+27.0*s[0]/32.0 \
                                                            +3.0*r[1]/4.0) \
                 + (DWKNASS)*(-3.0/8.0+5.0*r[0]/8.0-9.0*s[0]/32.0) \
                 + (DWVNASS)*(-r[1]/2.0) + (DWVNALS)*(3.0*r[1]/4.0) \
                 + (DWVKLS)*(-3.0*r[1]/4.0) \
                 + (AWLS)*(-27.0/8.0+27.0*r[0]/4.0+81.0*s[0]/16.0) \
                 + (AWSS)*(-3.0/8.0+3.0*r[0]/4.0-3.0*s[0]/16.0)
#define D2GDS0DT R*(0.75*log(xkls) - 0.75*log(xnals) \
                 - 3.0*0.25*log(xkss) + 3.0*0.25*log(xnass)) \
                 - 0.25*(2.0*(SEX)+(SX)) - (1.0/8.0)*(3.0*(SX))*s[0] \
                 + 0.5*((SX))*r[0] - 0.125*((SX)-(SEX)-2.0*(S23))*r[1] \
                 + (9.0*(SEX)+3.0*(SX)-6.0*(S23)+36.0*(S24))*r[2]/24.0
#define D2GDS0DP 0.25*(2.0*(VEX)+(VX)+3.0*(WVNAKLS)-(WVNAKSS)) \
                 + (1.0/8.0)*(3.0*(VX)-9.0*(WVNAKLS)-(WVNAKSS))*s[0] \
                 - 0.5*((VX)+3.0*(WVNAKLS)-(WVNAKSS))*r[0] \
                 + 0.125*((VX)-(VEX)-6.0*(WVNAKLS))*r[1] \
                 + (-18.0*(WVNAKLS)+4.0*(WVNAKSS)-9.0*(VEX)-3.0*(VX))*r[2]/24.0

#define D2GDT2   0.0
#define D2GDTDP  0.0
#define D2GDP2   0.0

/*----------------------------------------------------------------------------*/

#define D3GDR0R0R0 - R*t*(1.0/(xkls*xkls) - 1.0/(xnals*xnals) \
                   + 3.0/(xkss*xkss) - 3.0/(xnass*xnass)) - 6.0*(DWKNALS) \
                   - 6.0*(DWKNASS) + (AWLS)*(12.0) + (AWSS)*(12.0)
#define D3GDR0R0R1 R*t/(xnals*xnals) \
                   + (DWKNALS)*(-4.0/3.0) + (DWVNASS)*(-8.0/3.0) \
                   - (DWVNALS)*(4.0/3.0) + (DWVKLS)*(4.0/3.0)

#define D3GDR0R0R2 R*t*(1.0/(xnals*xnals) + 1.0/(xnass*xnass))
#define D3GDR0R0S0 - R*t*(0.75/(xkls*xkls) - 0.75/(xnals*xnals) \
                   - 3.0*0.25/(xkss*xkss) + 3.0*0.25/(xnass*xnass)) \
                   + (DWKNALS)*(-5.0/2.0) + (DWKNASS)*(-1.0/2.0) \
                   + (AWLS)*(9.0) + (AWSS)*(-6.0/2.0)
#define D3GDR0R0DT R*(1.0/xkls + 1.0/xnals + 3.0/xkss + 3.0/xnass) \
                   + 2.0*(SX)
#define D3GDR0R0DP - 2.0*((VX)+(WVNAKLS)+(WVNAKSS))

#define D3GDR0R1R1 R*t/(xnals*xnals) \
                   + (DWKNALS)*(2.0/3.0) + (DWVNASS)*(-2.0/3.0) \
                   + (DWVNALS)*(-4.0/3.0) + (DWVKLS)*(-2.0/3.0) \
                   + (AWLS)*(-4.0) + (AWVNA)*(4.0) + (AWVK)*(-4.0)
#define D3GDR0R1R2 R*t/(xnals*xnals)
#define D3GDR0R1S0 R*t*0.75/(xnals*xnals) \
                   + (DWKNALS)*(3.0*s[0]/4.0)
#define D3GDR0R1DT R/xnals - 0.5*(2.0*(S23)+(SEX)-(SX))
#define D3GDR0R1DP 0.5*((VEX)-(VX)-2.0*(WVNAKLS))

#define D3GDR0R2R2 R*t*(1.0/(xnals*xnals) + 1.0/(3.0*xnass*xnass))
#define D3GDR0R2S0 R*t*(0.75/(xnals*xnals) - 0.25/(xnass*xnass))
#define D3GDR0R2DT R*(1.0/xnals + 1.0/xnass) \
                                                                            + (-6.0*(S23)-3.0*(SEX)+15.0*(SX))/18.0
#define D3GDR0R2DP (3.0*(VEX)-15.0*(VX)-18.0*(WVNAKLS)-4.0*(WVNAKSS))/18.0

#define D3GDR0S0S0 - R*t*(0.75*0.75/(xkls*xkls) - 0.75*0.75/(xnals*xnals) \
                   + 3.0*0.25*0.25/(xkss*xkss) - 3.0*0.25*0.25/(xnass*xnass)) \
                   + (DWKNALS)*(-6.0/16.0+6.0*r[1]/8.0) \
                   + (DWKNASS)*(10.0/16.0) + (AWLS)*(27.0/4.0) + (AWSS)*(3.0/4.0)
#define D3GDR0S0DT R*(0.75/xkls + 0.75/xnals - 3.0*0.25/xkss \
                   - 3.0*0.25/xnass) + 0.5*(SX)
#define D3GDR0S0DP - 0.5*((VX)+3.0*(WVNAKLS)-(WVNAKSS))

#define D3GDR1R1R1 - R*t*(1.0/(xvcls*xvcls) - 1.0/(xnals*xnals)) \
                   + (AWVNA)*(12.0) \
                   - 6.0*(PENALTY)/(r[1]*r[1]*r[1]*r[1])
#define D3GDR1R1R2 R*t*(-1.0/(xvcls*xvcls) + 1.0/(xnals*xnals))
#define D3GDR1R1S0 R*t*0.75/(xnals*xnals) \
                   + (DWKNALS)*(1.0/2.0) + (DWVNASS)*(-1.0/2.0) \
                   + (DWVNALS) + (DWVKLS)*(-1.0/2.0) + (AWLS)*(-6.0/2.0) \
                   + (AWVNA)*(6.0/2.0)+ (AWVK)*(-6.0/2.0)
#define D3GDR1R1DT R*(1.0/xvcls + 1.0/xnals)
#define D3GDR1R1DP 0.0

#define D3GDR1R2R2 R*t*(-1.0/(xvcls*xvcls) + 1.0/(xnals*xnals))
#define D3GDR1R2S0 R*t*0.75/(xnals*xnals)
#define D3GDR1R2DT  R*(1.0/xvcls + 1.0/xnals)
#define D3GDR1R2DP 0.0

#define D3GDR1S0S0 R*t*0.75*0.75/(xnals*xnals) \
                   + (DWKNALS)*(3.0/4.0) + (DWVNASS)*(-1.0/2.0) \
                   + (DWVNALS)*(3.0/4.0) + (DWVKLS)*(-3.0/4.0)
#define D3GDR1S0DT R*0.75/xnals - 0.125*(-2.0*(S23)+(SX)-(SEX))
#define D3GDR1S0DP 0.125*((VX)-(VEX)-6.0*(WVNAKLS))

#define D3GDR2R2R2 R*t*(-1.0/(xvcls*xvcls) + 1.0/(xnals*xnals) \
                                                        - 1.0/SQUARE(1.0-3.0*xcass) \
                                                        - 1.0/SQUARE(3.0*xcass) + 1.0/SQUARE(3.0*xnass))
#define D3GDR2R2S0 R*t*(0.75/(xnals*xnals) - 0.25/(3.0*xnass*xnass))
#define D3GDR2R2DT R*(1.0/xvcls + 1.0/xnals - 1.0/(1.0-3.0*xcass) \
                               + 1.0/(3.0*xcass) + 1.0/(3.0*xnass))
#define D3GDR2R2DP 0.0

#define D3GDR2S0S0 R*t*(0.75*0.75/(xnals*xnals) + 0.25*0.25/(xnass*xnass))
#define D3GDR2S0DT R*(0.75/xnals - 0.25/xnass) \
                                                                    + (9.0*(SEX)+3.0*(SX)-6.0*(S23)+36.0*(S24))/24.0
#define D3GDR2S0DP (-18.0*(WVNAKLS)+4.0*(WVNAKSS)-9.0*(VEX)-3.0*(VX))/24.0

#define D3GDS0S0S0 - R*t*(0.75*0.75*0.75/(xkls*xkls) \
                   - 0.75*0.75*0.75/(xnals*xnals) - 3.0*0.25*0.25*0.25/(xkss*xkss) \
                   + 3.0*0.25*0.25*0.25/(xnass*xnass)) + (DWKNALS)*(27.0/32.0) \
                   + (DWKNASS)*(-9.0/32.0) + (AWLS)*(81.0/16.0) \
                   + (AWSS)*(-3.0/16.0)
#define D3GDS0S0DT R*(0.75*0.75/xkls + 0.75*0.75/xnals \
                   + 3.0*0.25*0.25/xkss + 3.0*0.25*0.25/xnass) \
                   - (1.0/8.0)*(3.0*(SX))
#define D3GDS0S0DP  (1.0/8.0)*(3.0*(VX)-9.0*(WVNAKLS)-(WVNAKSS))
#define D3GDS0DT2   0.0
#define D3GDS0DTDP  0.0
#define D3GDS0DP2   0.0

#define D3GDT3      0.0
#define D3GDT2DP    0.0
#define D3GDTDP2    0.0
#define D3GDP3      0.0

#define D3GDR0DT2   0.0
#define D3GDR0DTDP  0.0
#define D3GDR0DP2   0.0
#define D3GDR1DT2   0.0
#define D3GDR1DTDP  0.0
#define D3GDR1DP2   0.0
#define D3GDR2DT2   0.0
#define D3GDR2DTDP  0.0
#define D3GDR2DP2   0.0

/*
 *=============================================================================
 * Macros for automatic array initialization
 */
#define fillD2GDR2 \
 d2gdr2[0][0] = D2GDR0R0;     d2gdr2[0][1] = D2GDR0R1;     d2gdr2[0][2] = D2GDR0R2; \
 d2gdr2[1][0] = d2gdr2[0][1]; d2gdr2[1][1] = D2GDR1R1;     d2gdr2[1][2] = D2GDR1R2; \
 d2gdr2[2][0] = d2gdr2[0][2]; d2gdr2[2][1] = d2gdr2[1][2]; d2gdr2[2][2] = D2GDR2R2;

#define fillD2GDRDS d2gdrds[0][0] = D2GDR0S0; d2gdrds[1][0] = D2GDR1S0; d2gdrds[2][0] = D2GDR2S0;
#define fillD2GDRDT d2gdrdt[0] = D2GDR0DT;    d2gdrdt[1] = D2GDR1DT;    d2gdrdt[2] = D2GDR2DT;
#define fillD2GDRDP d2gdrdp[0] = D2GDR0DP;    d2gdrdp[1] = D2GDR1DP;    d2gdrdp[2] = D2GDR2DP;
#define fillD2GDS2  d2gds2[0][0] = D2GDS0S0;
#define fillD2GDSDT d2gdsdt[0] = D2GDS0DT;
#define fillD2GDSDP d2gdsdp[0] = D2GDS0DP;

#define fillD3GDR3 \
 d3gdr3[0][0][0] = D3GDR0R0R0;		d3gdr3[0][0][1] = D3GDR0R0R1;      d3gdr3[0][0][2] = D3GDR0R0R2; \
 d3gdr3[0][1][0] = d3gdr3[0][0][1];	d3gdr3[0][1][1] = D3GDR0R1R1;      d3gdr3[0][1][2] = D3GDR0R1R2; \
 d3gdr3[0][2][0] = d3gdr3[0][0][2];     d3gdr3[0][2][1] = d3gdr3[0][1][2]; d3gdr3[0][2][2] = D3GDR0R2R2; \
 d3gdr3[1][0][0] = d3gdr3[0][0][1];	d3gdr3[1][0][1] = d3gdr3[0][1][1]; d3gdr3[1][0][2] = d3gdr3[0][1][2]; \
 d3gdr3[1][1][0] = d3gdr3[0][1][1];	d3gdr3[1][1][1] = D3GDR1R1R1;      d3gdr3[1][1][2] = D3GDR1R1R2; \
 d3gdr3[1][2][0] = d3gdr3[0][1][2];	d3gdr3[1][2][1] = d3gdr3[1][1][2]; d3gdr3[1][2][2] = D3GDR1R2R2; \
 d3gdr3[2][0][0] = d3gdr3[0][0][2];	d3gdr3[2][0][1] = d3gdr3[0][1][2]; d3gdr3[2][0][2] = d3gdr3[0][2][2]; \
 d3gdr3[2][1][0] = d3gdr3[0][1][2];	d3gdr3[2][1][1] = d3gdr3[1][1][2]; d3gdr3[2][1][2] = d3gdr3[1][2][2]; \
 d3gdr3[2][2][0] = d3gdr3[0][2][2];	d3gdr3[2][2][1] = d3gdr3[1][2][2]; d3gdr3[2][2][2] = D3GDR2R2R2; \

#define fillD3GDR2DS \
 d3gdr2ds[0][0][0] = D3GDR0R0S0;        d3gdr2ds[0][1][0] = D3GDR0R1S0;        d3gdr2ds[0][2][0] = D3GDR0R2S0; \
 d3gdr2ds[1][0][0] = d3gdr2ds[0][1][0]; d3gdr2ds[1][1][0] = D3GDR1R1S0;        d3gdr2ds[1][2][0] = D3GDR1R2S0; \
 d3gdr2ds[2][0][0] = d3gdr2ds[0][2][0]; d3gdr2ds[2][1][0] = d3gdr2ds[1][2][0]; d3gdr2ds[2][2][0] = D3GDR2R2S0;

#define fillD3GDR2DT \
 d3gdr2dt[0][0] = D3GDR0R0DT;     d3gdr2dt[0][1] = D3GDR0R1DT;     d3gdr2dt[0][2] = D3GDR0R2DT; \
 d3gdr2dt[1][0] = d3gdr2dt[0][1]; d3gdr2dt[1][1] = D3GDR1R1DT;     d3gdr2dt[1][2] = D3GDR1R2DT; \
 d3gdr2dt[2][0] = d3gdr2dt[0][2]; d3gdr2dt[2][1] = d3gdr2dt[1][2]; d3gdr2dt[2][2] = D3GDR2R2DT;

#define fillD3GDR2DP \
 d3gdr2dp[0][0] = D3GDR0R0DP;     d3gdr2dp[0][1] = D3GDR0R1DP;     d3gdr2dp[0][2] = D3GDR0R2DP; \
 d3gdr2dp[1][0] = d3gdr2dp[0][1]; d3gdr2dp[1][1] = D3GDR1R1DP;     d3gdr2dp[1][2] = D3GDR1R2DP; \
 d3gdr2dp[2][0] = d3gdr2dp[0][2]; d3gdr2dp[2][1] = d3gdr2dp[1][2]; d3gdr2dp[2][2] = D3GDR2R2DP;

#define fillD3GDRDS2 \
 d3gdrds2[0][0][0] = D3GDR0S0S0;  d3gdrds2[1][0][0] = D3GDR1S0S0;  d3gdrds2[2][0][0] = D3GDR2S0S0;

#define fillD3GDRDSDT \
 d3gdrdsdt[0][0] = D3GDR0S0DT;    d3gdrdsdt[1][0] = D3GDR1S0DT;    d3gdrdsdt[2][0] = D3GDR2S0DT;

#define fillD3GDRDSDP \
 d3gdrdsdp[0][0] = D3GDR0S0DP;    d3gdrdsdp[1][0] = D3GDR1S0DP;    d3gdrdsdp[2][0] = D3GDR2S0DP;

#define fillD3GDS3    d3gds3[0][0][0] = D3GDS0S0S0;
#define fillD3GDS2DT  d3gds2dt[0][0] = D3GDS0S0DT;
#define fillD3GDS2DP  d3gds2dp[0][0] = D3GDS0S0DP;
#define fillD3GDSDT2  d3gdsdt2[0] = D3GDS0DT2;
#define fillD3GDSDTDP d3gdsdtdp[0] = D3GDS0DTDP;
#define fillD3GDSDP2  d3gdsdp2[0] = D3GDS0DP2;
#define fillD3GDRDT2  d3gdrdt2[0] = D3GDR0DT2;   d3gdrdt2[1] = D3GDR1DT2;   d3gdrdt2[2] = D3GDR2DT2;
#define fillD3GDRDTDP d3gdrdtdp[0] = D3GDR0DTDP; d3gdrdtdp[1] = D3GDR1DTDP; d3gdrdtdp[2] = D3GDR2DTDP;
#define fillD3GDRDP2  d3gdrdp2[0] = D3GDR0DP2;   d3gdrdp2[1] = D3GDR1DP2;   d3gdrdp2[2] = D3GDR2DP2;

/*
 *=============================================================================
 * Local function to compute ordering state and associated derivatives
 */

static void
order(int mask, double t, double p, double r[NR],
            double s[NS],           /* s[NS]                BINARY MASK: 0000000001 */
            double dr[NS][NR] ,     /* ds[NS]/dr[NR]        BINARY MASK: 0000000010 */
            double dt[NS],          /* ds[NS]/dt            BINARY MASK: 0000000100 */
            double dp[NS],          /* ds[NS]/dp            BINARY MASK: 0000001000 */
            double dr2[NS][NR][NR], /* d2s[NS]/dr[NR]dr[NR] BINARY MASK: 0000010000 */
            double drt[NS][NR],     /* d2s[NS]/dr[NR]dt     BINARY MASK: 0000100000 */
            double drp[NS][NR],     /* d2s[NS]/dr[NR]dp     BINARY MASK: 0001000000 */
            double dt2[NS],         /* d2s[NS]/dt2          BINARY MASK: 0010000000 */
            double dtp[NS],         /* d2s[NS]/dtp          BINARY MASK: 0100000000 */
            double dp2[NS]          /* d2s[NS]/dp2          BINARY MASK: 1000000000 */
            )
{
    DECLARE_SITE_FRACTIONS
    double tOld         = getTOld();
    double pOld         = getPOld();
    double *rOld        = getROld();
    double *sOld        = getSOld();
    double **d2gds2     = getD2gds2();
    gsl_matrix      *ptToD2gds2  = getPtToD2gds2();
    gsl_permutation *indexD2gds2 = getIndexD2gds2();
    int i, j, iter = 0, signum;

    GET_SITE_FRACTIONS

    /* look-up or compute the current ordering state */
    if ( (t != tOld)       || (p != pOld) ||
       (r[0] != rOld[0]) || (r[1] != rOld[1]) || (r[2] != rOld[2]) ) {
        double dgds[NS], sNew[NS];
        int skipCheck = FALSE;
        double xk  = r[0];
        double xvc = (r[1]+r[2])/4.0;
        double xna = 1.0 - r[0] - r[1]/4.0 - r[2]/2.0;
        gsl_vector_view vvToDgds = gsl_vector_view_array(dgds, (size_t) NS);

        for (i=0; i<NS; i++) { sOld[i] = 2.0; sNew[i] = 0.0; dgds[i] = 0.0; }

        /* Initial guess assumes random distribution */
        sNew[0] = -4.0*r[0]*(r[1]+2.0*r[2]/3.0)/(4.0-r[1]-2.0*r[2]);

        xkls  = r[0] + 3.0*sNew[0]/4.0;
        xvcls = r[1] + r[2];
        xnals = 1.0 - r[0] - r[1] - r[2] - 3.0*sNew[0]/4.0;
        xkss  = r[0] - sNew[0]/4.0;
        xcass = r[2]/3.0;
        xnass = 1.0 - r[0] - r[2]/3.0 + sNew[0]/4.0;

        if (xkls  <= DBL_EPSILON) xkls  = DBL_EPSILON;
        if (xvcls <= DBL_EPSILON) xvcls = DBL_EPSILON;
        if (xnals <= DBL_EPSILON) xnals = DBL_EPSILON;
        if (xkss  <= DBL_EPSILON) xkss  = DBL_EPSILON;
        if (xcass <= DBL_EPSILON) xcass = DBL_EPSILON;
        if (xnass <= DBL_EPSILON) xnass = DBL_EPSILON;

        if (xkls  >= 1.0     - DBL_EPSILON) xkls  = 1.0     - DBL_EPSILON;
        if (xvcls >= 1.0     - DBL_EPSILON) xvcls = 1.0     - DBL_EPSILON;
        if (xnals >= 1.0     - DBL_EPSILON) xnals = 1.0     - DBL_EPSILON;
        if (xkss  >= 1.0     - DBL_EPSILON) xkss  = 1.0     - DBL_EPSILON;
        if (xcass >= 1.0/3.0 - DBL_EPSILON) xcass = 1.0/3.0 - DBL_EPSILON;
        if (xnass >= 1.0     - DBL_EPSILON) xnass = 1.0     - DBL_EPSILON;

        sNew[0] = xkls - xkss;

        /* Check for irrelevant compositions */
        if (xk  < sqrt(DBL_EPSILON) || xna < sqrt(DBL_EPSILON) ||
                xvc > ((double) (1.0/4.0) - sqrt(DBL_EPSILON)) ) {
            for (i=0; i<NS; i++) sOld[i] = sNew[i];
            skipCheck = TRUE;
        }

        while ( (ABS(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON) && (iter < MAX_ITER)) {
            double s[NS], sCorr[NS], lambda = 1.0;
            gsl_vector_view vvToSCorr = gsl_vector_view_array(sCorr, (size_t) NS);

            for (i=0; i<NS; i++) s[i] = sNew[i];

            dgds[0] = DGDS0;
            d2gds2[0][0] = D2GDS0S0;

            for (i=0; i<NS; i++) sOld[i] = s[i];

            /* original: gaussj(ptToD2gds2, NS, (double **) NULL, 0);
              	for (i=0; i<NS; i++) {
                    for(j=0, sCorr[i]=0.0; j<NS; j++) sCorr[i] += - d2gds2[i][j]*dgds[j];
                    s[i] = sOld[i] + lambda*sCorr[i];
        	  }*/

            gsl_matrix_scale(ptToD2gds2, -1.0);
            melts_LU_decomp(ptToD2gds2, indexD2gds2, &signum);

            melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToDgds.vector, &vvToSCorr.vector);
            for (i=0; i<NS; i++) s[i] += lambda*sCorr[i];

            xkls  = r[0] + 3.0*s[0]/4.0;
            xvcls = r[1] + r[2];
            xnals = 1.0 - r[0] - r[1] - r[2] - 3.0*s[0]/4.0;
            xkss  = r[0] - s[0]/4.0;
            xcass = r[2]/3.0;
            xnass = 1.0 - r[0] - r[2]/3.0 + s[0]/4.0;

            while ((   xkls  < 0.0 || xkls  > 1.0     || xvcls < 0.0 || xvcls > 1.0
                            || xnals < 0.0 || xnals > 1.0     || xkss  < 0.0 || xkss  > 1.0
                            || xcass < 0.0 || xcass > 1.0/3.0 || xnass < 0.0 || xnass > 1.0)
             && lambda > DBL_EPSILON ) {
                lambda /= 2.0;
                for (i=0; i<NS; i++) s[i] = sOld[i] + lambda*sCorr[i];
                xkls  = r[0] + 3.0*s[0]/4.0;
                xvcls = r[1] + r[2];
                xnals = 1.0 - r[0] - r[1] - r[2] - 3.0*s[0]/4.0;
                xkss  = r[0] - s[0]/4.0;
                xcass = r[2]/3.0;
                xnass = 1.0 - r[0] - r[2]/3.0 + s[0]/4.0;
            }

            if (xkls  <= DBL_EPSILON) xkls  = DBL_EPSILON;
            if (xvcls <= DBL_EPSILON) xvcls = DBL_EPSILON;
            if (xnals <= DBL_EPSILON) xnals = DBL_EPSILON;
            if (xkss  <= DBL_EPSILON) xkss  = DBL_EPSILON;
            if (xcass <= DBL_EPSILON) xcass = DBL_EPSILON;
            if (xnass <= DBL_EPSILON) xnass = DBL_EPSILON;

            if (xkls  >= 1.0     - DBL_EPSILON) xkls  = 1.0     - DBL_EPSILON;
            if (xvcls >= 1.0     - DBL_EPSILON) xvcls = 1.0     - DBL_EPSILON;
            if (xnals >= 1.0     - DBL_EPSILON) xnals = 1.0     - DBL_EPSILON;
            if (xkss  >= 1.0     - DBL_EPSILON) xkss  = 1.0     - DBL_EPSILON;
            if (xcass >= 1.0/3.0 - DBL_EPSILON) xcass = 1.0/3.0 - DBL_EPSILON;
            if (xnass >= 1.0     - DBL_EPSILON) xnass = 1.0     - DBL_EPSILON;

            s[0] = xkls - xkss;

            for (i=0; i<NS; i++) sNew[i] = s[i];
            iter++;
        }
        tOld = t;
        pOld = p;
        for (i=0; i<NR; i++) rOld[i] = r[i];

#ifdef DEBUG
        for (i=0; i<NS; i++) {
            if (!skipCheck && ABS(dgds[i]) > sqrt(DBL_EPSILON)
                    && (ABS(          xk *dgds[i]) > sqrt(DBL_EPSILON) &&
                            ABS(          xna*dgds[i]) > sqrt(DBL_EPSILON) &&
                            ABS((1.0/4.0-xvc)*dgds[i]) > sqrt(DBL_EPSILON) )
                    && ABS(sOld[i]) > DBL_EPSILON) {
                printf("ERROR in NEPHELINE.C (function ORDER). Failed to converge!\n");
                if (iter >= MAX_ITER)
                    printf("  Iteration limit (%4d) exceeded.\n", iter);
                printf("  X2    = %13.6g, X3    = %13.6g\n", r[0], r[1]);
                printf("  X4    = %13.6g, s1    = %13.6g\n", r[2], sOld[0]);
                printf("  dgds1 = %13.6g\n", dgds[0]);
                printf("  X K  ls: %13.6g  X K  ss: %13.6g\n", xkls,  xkss);
                printf("  X Na ls: %13.6g  X Na ss: %13.6g\n", xnals, xnass);
                printf("  X Vc ls: %13.6g  X Ca ss: %13.6g\n", xvcls, xcass);
                break;
            }
        }
#endif

        setTOld(tOld);
        setPOld(pOld);
        /* arrays (rOld, sOld, d2gds2, indexD2gds2) should be preserved automatically */

        SET_SITE_FRACTIONS
    }

    if (mask & FIRST  ) {   /* return s        */
        for (i=0; i<NS; i++) s[i] = sOld[i];
    }

    if (mask & SECOND) {   /* compute ds/dr:  */
        double *s = sOld;
        double d2gdrds[NR][NS];
        gsl_matrix_view mvToDr = gsl_matrix_view_array((double *) dr, (size_t) NS, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NS);

        fillD2GDRDS

        /* original: dr[i][j] += - d2gds2[i][k]*d2gdrds[j][k]; */
        for (j=0; j<NR; j++) {
            gsl_vector_view vvToDr = gsl_matrix_column(&mvToDr.matrix, j),
            	vvToD2gdrds = gsl_matrix_row(&mvToD2gdrds.matrix, j);
            melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdrds.vector, &vvToDr.vector);
        }
    }

    if (mask & THIRD) {   /* compute ds/dt:  */
        double *s = sOld;
        double d2gdsdt[NS];
        gsl_vector_view vvToDt = gsl_vector_view_array(dt, (size_t) NS),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NS);

        fillD2GDSDT

        /* original: dt[i] += - d2gds2[i][j]*d2gdsdt[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdt.vector, &vvToDt.vector);

    }

    if (mask & FOURTH) {   /* compute ds/dp:  */
        double *s = sOld;
        double d2gdsdp[NS];
        gsl_vector_view vvToDp = gsl_vector_view_array(dp, (size_t) NS),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NS);

        fillD2GDSDP

        /* original: dp[i] += - d2gds2[i][j]*d2gdsdp[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdp.vector, &vvToDp.vector);

    }

    if (mask & FIFTH) {   /* compute d2s/dr2 */
        double *s = sOld;
        double d2gdrds[NR][NS], d3gdr2ds[NR][NR][NS], d3gdrds2[NR][NS][NS],
            d3gds3[NS][NS][NS], dsdr[NS][NR], temp[NS];
        gsl_matrix_view mvToDsdr = gsl_matrix_view_array((double *) dsdr, (size_t) NS, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NS);
        gsl_vector_view vvToTemp = gsl_vector_view_array(temp, (size_t) NS);
        int k, l, m, n;

        fillD2GDRDS
        fillD3GDR2DS
        fillD3GDRDS2
        fillD3GDS3

        /* compute dsdr matrix */
        /* original: dsdr[i][j] += - d2gds2[i][k]*d2gdrds[j][k]; */
        for (j=0; j<NR; j++) {
            gsl_vector_view vvToDsdr = gsl_matrix_column(&mvToDsdr.matrix, j),
            	vvToD2gdrds = gsl_matrix_row(&mvToD2gdrds.matrix, j);
            melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdrds.vector, &vvToDsdr.vector);
        }

        /* compute dsdr2 cube */
        for (j=0; j<NR; j++) {
            for (k=0; k<NR; k++) {
                for (l=0; l<NS; l++) {
                    temp[l] = d3gdr2ds[j][k][l];
                    for (m=0; m<NS; m++) {
                        temp[l] += d3gdrds2[j][l][m]*dsdr[m][k]
                            + d3gdrds2[k][l][m]*dsdr[m][j];
                        for (n=0; n<NS; n++)
                            temp[l] += d3gds3[l][m][n]*dsdr[m][j]*dsdr[n][k];
                    }
                    }
                /* original: dr2[i][j][k] += - d2gds2[i][l]*temp[l]; */
                melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
                for (l=0; l<NS; l++) dr2[l][j][k] = temp[l];
            }
        }

    }

    if (mask & SIXTH) {   /* compute d2s/drt */
        double *s = sOld;
        double d2gdrds[NR][NS], d2gdsdt[NS], d3gdrds2[NR][NS][NS],
            d3gdrdsdt[NR][NS], d3gds3[NS][NS][NS], d3gds2dt[NS][NS], dsdr[NS][NR],
            dsdt[NS], temp[NS];
        gsl_matrix_view mvToDsdr = gsl_matrix_view_array((double *) dsdr, (size_t) NS, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NS);
        gsl_vector_view vvToDsdt = gsl_vector_view_array(dsdt, (size_t) NS),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NS),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NS);
        int k, l, m;

        fillD2GDRDS
        fillD2GDSDT
        fillD3GDRDS2
        fillD3GDRDSDT
        fillD3GDS3
        fillD3GDS2DT

        /* compute dsdr matrix */
        /* original: dsdr[i][j] += - d2gds2[i][k]*d2gdrds[j][k]; */
        for (j=0; j<NR; j++) {
            gsl_vector_view vvToDsdr = gsl_matrix_column(&mvToDsdr.matrix, j),
            	vvToD2gdrds = gsl_matrix_row(&mvToD2gdrds.matrix, j);
            melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdrds.vector, &vvToDsdr.vector);
        }

        /* compute dsdt vector */
        /* original: dsdt[i] += - d2gds2[i][j]*d2gdsdt[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdt.vector, &vvToDsdt.vector);

        /* compute dsdrdt matrix */
        for (j=0; j<NR; j++) {
            for (k=0; k<NS; k++) {
                temp[k] = d3gdrdsdt[j][k];
                for (l=0; l<NS; l++) {
                    temp[k] += d3gdrds2[j][k][l]*dsdt[l] + d3gds2dt[k][l]*dsdr[l][j];
                    for (m=0; m<NS; m++) temp[k] += d3gds3[k][l][m]*dsdr[l][j]*dsdt[m];
                }
            }
            /* original: drt[i][j] += - d2gds2[i][k]*temp[k]; */
            melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
            for (k=0; k<NS; k++) drt[k][j] = temp[k];
        }

    }
    if (mask & SEVENTH) {   /* compute d2s/drp */
        double *s = sOld;
        double d2gdrds[NR][NS], d2gdsdp[NS], d3gdrds2[NR][NS][NS],
            d3gdrdsdp[NR][NS], d3gds3[NS][NS][NS], d3gds2dp[NS][NS], dsdr[NS][NR],
            dsdp[NS], temp[NS];
        gsl_matrix_view mvToDsdr = gsl_matrix_view_array((double *) dsdr, (size_t) NS, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NS);
        gsl_vector_view vvToDsdp = gsl_vector_view_array(dsdp, (size_t) NS),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NS),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NS);
        int k, l, m;

        fillD2GDRDS
        fillD2GDSDP
        fillD3GDRDS2
        fillD3GDRDSDP
        fillD3GDS3
        fillD3GDS2DP

        /* compute dsdr matrix */
        /* original: dsdr[i][j] += - d2gds2[i][k]*d2gdrds[j][k]; */
        for (j=0; j<NR; j++) {
            gsl_vector_view vvToDsdr = gsl_matrix_column(&mvToDsdr.matrix, j),
            	vvToD2gdrds = gsl_matrix_row(&mvToD2gdrds.matrix, j);
            melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdrds.vector, &vvToDsdr.vector);
        }

        /* compute dsdp vector */
        /* original: dsdp[i] += - d2gds2[i][j]*d2gdsdp[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdp.vector, &vvToDsdp.vector);

        /* compute dsdrdp matrix */
        for (j=0; j<NR; j++) {
            for (k=0; k<NS; k++) {
                temp[k] = d3gdrdsdp[j][k];
                for (l=0; l<NS; l++) {
                    temp[k] += d3gdrds2[j][k][l]*dsdp[l] + d3gds2dp[k][l]*dsdr[l][j];
                    for (m=0; m<NS; m++) temp[k] += d3gds3[k][l][m]*dsdr[l][j]*dsdp[m];
                }
            }
            /* original: drp[i][j] += - d2gds2[i][k]*temp[k]; */
            melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
            for (k=0; k<NS; k++) drp[k][j] = temp[k];
        }

    }
    if (mask & EIGHTH) {   /* compute d2s/dt2 */
        double *s = sOld;
        double d2gdsdt[NS], d3gds3[NS][NS][NS], d3gds2dt[NS][NS], d3gdsdt2[NS],
            dsdt[NS], temp[NS];
        gsl_vector_view vvToDsdt = gsl_vector_view_array(dsdt, (size_t) NS),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NS),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NS);
        int k, l;

        fillD2GDSDT
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDSDT2

        /* compute dsdt vector */
        /* original: dsdt[i] += - d2gds2[i][j]*d2gdsdt[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdt.vector, &vvToDsdt.vector);

        /* compute dsdt2 vector */
        for (j=0; j<NS; j++) {
            temp[j] = d3gdsdt2[j];
            for (k=0; k<NS; k++) {
                temp[j] +=  2.0*d3gds2dt[j][k]*dsdt[k];
                for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdt[l];
            }
        }
        /* original: dt2[i] += - d2gds2[i][j]*temp[j]; */
        melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
        for (j=0; j<NS; j++) dt2[j] = temp[j];

    }
    if (mask & NINTH) {   /* compute d2s/dtp */
        double *s = sOld;
        double d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS][NS][NS], d3gds2dt[NS][NS],
            d3gds2dp[NS][NS], d3gdsdtdp[NS], dsdt[NS], dsdp[NS], temp[NS];
        gsl_vector_view vvToDsdt = gsl_vector_view_array(dsdt, (size_t) NS),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NS),
            vvToDsdp = gsl_vector_view_array(dsdp, (size_t) NS),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NS),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NS);
        int k, l;

        fillD2GDSDT
        fillD2GDSDP
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDS2DP
        fillD3GDSDTDP

        /* compute dsdt vector */
        /* original: dsdt[i] += - d2gds2[i][j]*d2gdsdt[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdt.vector, &vvToDsdt.vector);

        /* compute dsdp vector */
        /* original: dsdp[i] += - d2gds2[i][j]*d2gdsdp[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdp.vector, &vvToDsdp.vector);

        /* compute dsdtp vector */
        for (j=0; j<NS; j++) {
            temp[j] = d3gdsdtdp[j];
            for (k=0; k<NS; k++) {
                temp[j] += d3gds2dt[j][k]*dsdp[k] + d3gds2dp[j][k]*dsdt[k];
                for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdp[l];
            }
        }
        /* original: dtp[i] += - d2gds2[i][j]*temp[j]; */
        melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
        for (j=0; j<NS; j++) dtp[j] = temp[j];

    }
    if (mask & TENTH) {   /* compute d2s/dp2 */
        double *s = sOld;
        double d2gdsdp[NS], d3gds3[NS][NS][NS], d3gds2dp[NS][NS], d3gdsdp2[NS],
            dsdp[NS], temp[NS];
        gsl_vector_view vvToDsdp = gsl_vector_view_array(dsdp, (size_t) NS),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NS),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NS);
        int k, l;

        fillD2GDSDP
        fillD3GDS3
        fillD3GDS2DP
        fillD3GDSDP2

        /* compute dsdp vector */
        /* original: dsdp[i] += - d2gds2[i][j]*d2gdsdp[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdp.vector, &vvToDsdp.vector);

        /* compute dsdp2 vector */
        for (j=0; j<NS; j++) {
            temp[j] = d3gdsdp2[j];
            for (k=0; k<NS; k++) {
                temp[j] +=  2.0*d3gds2dp[j][k]*dsdp[k];
                for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdp[k]*dsdp[l];
            }
        }
        /* original: dp2[i] += - d2gds2[i][j]*temp[j]; */
        melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
        for (j=0; j<NS; j++) dp2[j] = temp[j];

    }

}

/*
 *=============================================================================
 * Public functions:
 *    mask  -  bitwise mask for selecting output
 *    t     -  Temperature (K)
 *    p     -  Pressure (bars)
 *    *x    -  (pointer to x[]) Array of independent compositional variables
 */
int
testNph(int mask, double t, double p,
    int na,          /* Expected number of endmember components                 */
    int nr,          /* Expected number of independent compositional variables  */
    char **names,    /* array of strings of names of endmember components       */
    char **formulas, /* array of strings of formulas of endmember components    */
    double *r,       /* array of indepependent compos variables, check bounds   */
    double *m)       /* array of moles of endmember components, check bounds    */
{
    const char *phase = "nepheline.c";
    const char *NAMES[NA]    = { "na-nepheline", "k-nepheline", "vc-nepheline", "ca-nepheline" };
    const char *FORMULAS[NA] = { "Na4Al4Si4O16", "K4Al4Si4O16", "Na3Al3Si5O16", "CaNa2Al4Si4O16" };
    int result = TRUE, i;
    double sum;

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
            result = result && (strcmp(names[i],NAMES[i]) == 0);
            if (!result)
                printf("<<%s>> Component[%d] should be %s not %s.\n",
                    phase, i, NAMES[i], names[i]);
        }
    }
    if (mask & FOURTH) {
        for (i=0; i<NA; i++) {
            result = result && (strcmp(formulas[i],FORMULAS[i]) == 0);
            if (!result)
                printf("<<%s>> Component[%d] should have formula %s not %s.\n",
                    phase, i, FORMULAS[i], formulas[i]);
        }
    }
    /* Check bounds on the independent compositional variables */
    if (mask & FIFTH) {
        result = result && (r[0] >= 0.0) && (r[0] <= 1.0-r[1]/4.0-r[2]/2.0);
        result = result && (r[1] >= 0.0) && (r[1] <= 1.0);
        result = result && (r[2] >= 0.0) && (r[2] <= 1.0);

        result = result && (4.0*r[0] >= 0.0) && (4.0*r[0] <= 4.0+sqrt(DBL_EPSILON));               /* tot K  */
        result = result && (4.0-4.0*r[0]-r[1]-2.0*r[2] >= 0.0)                                     /* tot Na */
                                        && (4.0-4.0*r[0]-r[1]-2.0*r[2] <= 4.0+sqrt(DBL_EPSILON));
        result = result && (r[2] >= 0.0) && (r[2] <= 1.0+sqrt(DBL_EPSILON));                       /* tot Ca */
        result = result && (r[1]+r[2] >= 0.0) && (r[1]+r[2] <= 1.0+sqrt(DBL_EPSILON));             /* tot Vc */
        result = result && (4.0-r[1] >= 3.0-sqrt(DBL_EPSILON)) && (4.0-r[1] <= 5.0+sqrt(DBL_EPSILON));   /* tot Al */
        result = result && (4.0+r[1] >= 3.0-sqrt(DBL_EPSILON)) && (4.0+r[1] <= 5.0+sqrt(DBL_EPSILON));   /* tot Si */
    }
    /* Check bounds on moles of endmember components */
    if (mask & SIXTH) {
        int pResult;
        for (i=0, sum=0.0; i<NA; i++) sum += m[i];
        result = result && (sum >= 0.0);
        if (sum > 0.0) {
            pResult = (m[0] >= -3.0*m[2]/4.0-m[3]/2.0-DBL_EPSILON) && (m[0] <= sum+DBL_EPSILON);
#ifdef DEBUG
            if (!pResult) printf("Bound check m[0] : %g <= %g <= %g\n", -3.0*m[2]/4.0-m[3]/2.0, m[0], sum);
#endif
            result = result && pResult;

            pResult = (m[1] >= 0.0) && (m[1] <= sum-m[2]/4.0-m[3]/2.0+DBL_EPSILON);
#ifdef DEBUG
            if (!pResult) printf("Bound check m[1] : %g <= %g <= %g\n", 0.0, m[1], sum-m[2]/4.0-m[3]/2.0);
#endif
            result = result && pResult;

            pResult = (m[2] >= 0.0) && (m[2] <= sum+DBL_EPSILON);
#ifdef DEBUG
            if (!pResult) printf("Bound check m[2] : %g <= %g <= %g\n", 0.0, m[2], sum);
#endif
            result = result && pResult;

            pResult = (m[3] >= 0.0) && (m[3] <= sum+DBL_EPSILON);
#ifdef DEBUG
            if (!pResult) printf("Bound check m[3] : %g <= %g <= %g\n", 0.0, m[3], sum);
#endif
            result = result && pResult;

            pResult = (4.0*m[1] >= 0.0) && (4.0*m[1] <= 4.0*sum+sqrt(DBL_EPSILON));       /* tot K  */
#ifdef DEBUG
            if (!pResult) printf("Bound check tot K : %g <= %g <= %g\n", 0.0, 4.0*m[1], 4.0*sum);
#endif
            result = result && pResult;

            pResult = (4.0*m[0]+3.0*m[2]+2.0*m[3] >= 0.0)                                 /* tot Na */
             && (4.0*m[0]+3.0*m[2]+2.0*m[3] <= 4.0*sum+sqrt(DBL_EPSILON));
#ifdef DEBUG
            if (!pResult) printf("Bound check tot Na: %g <= %g <= %g\n", 0.0, 4.0*m[0]+3.0*m[2]+2.0*m[3], 4.0*sum);
#endif
            result = result && pResult;

            pResult = (m[3] >= 0.0) && (m[3] <= sum+sqrt(DBL_EPSILON));                   /* tot Ca */
#ifdef DEBUG
            if (!pResult) printf("Bound check tot Ca: %g <= %g <= %g\n", 0.0, m[3], sum);
#endif
            result = result && pResult;

            pResult = (m[2]+m[3] >= 0.0) && (m[2]+m[3] <= sum+sqrt(DBL_EPSILON));         /* tot Vc */
#ifdef DEBUG
            if (!pResult) printf("Bound check tot Vc: %g <= %g <= %g\n", 0.0, m[2]+m[3], sum);
#endif
            result = result && pResult;

            pResult = (4.0*m[0]+4.0*m[1]+3.0*m[2]+4.0*m[3] >= 3.0*sum-sqrt(DBL_EPSILON))  /* tot Al */
             && (4.0*m[0]+4.0*m[1]+3.0*m[2]+4.0*m[3] <= 5.0*sum+sqrt(DBL_EPSILON));
#ifdef DEBUG
            if (!pResult) printf("Bound check tot Al: %g <= %g <= %g\n", 3.0*sum,
                4.0*m[0]+4.0*m[1]+3.0*m[2]+4.0*m[3], 5.0*sum);
#endif
            result = result && pResult;

            pResult = (4.0*m[0]+4.0*m[1]+5.0*m[2]+4.0*m[3] >= 3.0*sum-sqrt(DBL_EPSILON))  /* tot Si */
             && (4.0*m[0]+4.0*m[1]+5.0*m[2]+4.0*m[3] <= 5.0*sum+sqrt(DBL_EPSILON));
#ifdef DEBUG
            if (!pResult) printf("Bound check tot Si: %g <= %g <= %g\n", 3.0*sum,
                4.0*m[0]+4.0*m[1]+5.0*m[2]+4.0*m[3], 5.0*sum);
#endif
            result = result && pResult;
        }
    }

    return result;
}

void
conNph(int inpMask, int outMask, double t, double p,
    double *e,      /* comp of nepheline in moles of elements                      */
    double *m,      /* comp of nepheline in moles of endmember components          */
    double *r,      /* comp of nepheline in terms of the independent comp var      */
    double *x,      /* comp of nepheline in mole fractions of endmember comp       */
    double **dm,    /* Jacobian matrix: dm[i][j] = dr[i]/dm[j]                  */
    double ***d2m,  /* vector of matrices: d2m[i][j][k] = d2r[i]/dm[j]dm[k]     */
    double **dr,    /* Jacobian matrix: dr[i][j] = dx[i]/dr[j]                  */
    double ****d3m) /* 3rd deriv matrix: d3m[i][j][k][l]=d3r[i]/dm[j]dm[k]dm[l] */
{
    /*---------------------------------------------------------------------------
    Not all combinations of inpMask and outMask are feasible. Valid
        combinations are:

       inpMask          outMask
    (1)  FIRST            SECOND
    (2)  SECOND           THIRD  | FOURTH  | FIFTH | SIXTH | EIGHTH
    (3)  THIRD            FOURTH | SEVENTH

    (1) converts a vector of moles of elements into a vector of moles of
            endmember nepheline components.
    (2) calculates from a vector of moles of endmember components, one or
            all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
    (3) calculates from a vector of independent compositional variables
            mole fractions of endmember components and/or the Jacobian matrix
            dx[]/dr[]

    In this routine it is assumed that the elements are in the order of atomic
    numbers and that the order of nepheline components has been verified as:
            m[0] = na-nepheline   (Na4Al4Si4O16) ,
            m[1] = k-nepheline    (K4Al4Si4O16),
            m[2] = vc-nepheline   (Na3Al3Si5O16),
            m[3] = ca-nepheline   (CaNa2Al4Si4O16),

    ----------------------------------------------------------------------------*/

    int i, j, k;

    if (inpMask == FIRST && outMask == SECOND) {
        /* Converts a vector of moles of elements into a vector of moles of
       end-member components.                                                 */
        static const int Na = 11;
        static const int Mg = 12;
        static const int Al = 13;
        static const int Si = 14;
        static const int K  = 19;
        static const int Ca = 20;
        static const int Fe = 26;
        double divalent = e[Ca] + e[Mg];
        double vacancy  = (e[Al]+e[Fe]+e[Si])/2.0 - e[Na] - e[K] - divalent; /* assume: Fe is trivalent */

        m[0] = (e[Na] - 3.0*vacancy + divalent)/4.0;  /* Moles of Na4Al4Si4O16 */
        m[1] = e[K]/4.0;                              /* Moles of K4Al4Si4O16  */
        m[2] = vacancy-divalent;                      /* Moles of Na3Al3Si5O16 */
        m[3] = divalent;                              /* Moles of CaNa2Al4Si4O16 */

    } else if (inpMask == SECOND) {
        double sum;

        if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
            printf("Illegal call to conNph with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        for (i=0, sum=0.0; i<NA; i++) sum += m[i];

        if (outMask & THIRD) {
            /* Converts a vector of moles of end-member components (m) into a vector
         of independent compositional variables (r) required as input for the
         remaining public functions.                                          */
            r[0] = (sum != 0.0) ? m[1]/sum : 0.0;  /* X2 = X K4Al4Si4O16    */
            r[1] = (sum != 0.0) ? m[2]/sum : 0.0;  /* X3 = X Na3Al3Si5O16   */
            r[2] = (sum != 0.0) ? m[3]/sum : 0.0;  /* X4 = X CaNa2Al4Si4O16 */
        }

        if (outMask & FOURTH) {
            /* Converts a vector of moles of end-member components (m) into a vector
         of mole fractions of endmember components                            */
            for (i=0; i<NA; i++) x[i] = (sum != 0.0) ? m[i]/sum : 0.0;
        }

        if (outMask & FIFTH) {
            /* Calculates the matrix dr[i]/dm[j] using m[] as input                 */

            if (sum == 0.0) {
                for (i=0; i<NR; i++) { for (j=0; j<NA; j++) dm[i][j] = 0.0; }
            } else {
                for (j=0; j<NA; j++) {
                    dm[0][j] = (j == 1) ? (1.0 - m[1]/sum)/sum : -m[1]/SQUARE(sum);
                    dm[1][j] = (j == 2) ? (1.0 - m[2]/sum)/sum : -m[2]/SQUARE(sum);
                    dm[2][j] = (j == 3) ? (1.0 - m[3]/sum)/sum : -m[3]/SQUARE(sum);
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
                for (j=0; j<NA; j++) {
                    for (k=0; k<NA; k++) {
                        d2m[0][j][k]  = 2.0*m[1]/CUBE(sum);
                        d2m[0][j][k] -= (j == 1) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[0][j][k] -= (k == 1) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[1][j][k]  = 2.0*m[2]/CUBE(sum);
                        d2m[1][j][k] -= (j == 2) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[1][j][k] -= (k == 2) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[2][j][k]  = 2.0*m[3]/CUBE(sum);
                        d2m[2][j][k] -= (j == 3) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[2][j][k] -= (k == 3) ? 1.0/SQUARE(sum) : 0.0;
                    }
                }
            }

        }

        if (outMask & EIGHTH) {
            /* calculates the matrix d3r[i]/dm[j]dm[k]dm[l] using m[] as input        */
            int l;

            if (sum == 0.0) {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++)  {
                        for (k=0; k<NA; k++)  {
                            for (l=0; l<NA; l++) d3m[i][j][k][l] = 0.0;
                        }
                    }
                }
            } else {
                for (j=0; j<NA; j++)  {
                    for (k=0; k<NA; k++)  {
                        for (l=0; l<NA; l++)  {
                            d3m[0][j][k][l]  = -6.0*m[1]/QUARTIC(sum);
                            d3m[0][j][k][l] += (j == 1) ? 2.0/CUBE(sum) : 0.0;
                            d3m[0][j][k][l] += (k == 1) ? 2.0/CUBE(sum) : 0.0;
                            d3m[0][j][k][l] += (l == 1) ? 2.0/CUBE(sum) : 0.0;
                            d3m[1][j][k][l]  = -6.0*m[2]/QUARTIC(sum);
                            d3m[1][j][k][l] += (j == 2) ? 2.0/CUBE(sum) : 0.0;
                            d3m[1][j][k][l] += (k == 2) ? 2.0/CUBE(sum) : 0.0;
                            d3m[1][j][k][l] += (l == 2) ? 2.0/CUBE(sum) : 0.0;
                            d3m[2][j][k][l]  = -6.0*m[3]/QUARTIC(sum);
                            d3m[2][j][k][l] += (j == 3) ? 2.0/CUBE(sum) : 0.0;
                            d3m[2][j][k][l] += (k == 3) ? 2.0/CUBE(sum) : 0.0;
                            d3m[2][j][k][l] += (l == 3) ? 2.0/CUBE(sum) : 0.0;
                        }
                    }
                }
            }
        }

    } else if (inpMask == THIRD) {

        if (outMask & ~(FOURTH | SEVENTH))
            printf("Illegal call to conNph with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        if (outMask & FOURTH) {
            /* Converts a vector of independent compositional variables (r) into a
         vector of mole fractions of endmember components (x).                */
            x[0] = 1.0 - r[0] - r[1] - r[2];
            x[1] = r[0];
            x[2] = r[1];
            x[3] = r[2];
        }

        if (outMask & SEVENTH) {
            /* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
            for (i=0; i<NA; i++) for (j=0; j<NR; j++) dr[i][j] = 0.0;
            dr[0][0] = -1.0; dr[0][1] = -1.0; dr[0][2] = -1.0;
            dr[1][0] =  1.0;
            dr[2][1] =  1.0;
            dr[3][2] =  1.0;
        }

    } else  {
        printf("Illegal call to conNph with inpMask = %o and outMask = %o\n",
            inpMask, outMask);
    }

}

void
dispNph(int mask, double t, double p, double *x,
    char **formula            /* Mineral formula for interface display MASK: 1 */
    )
{
    double *r = x;
    static char masterString[] = {
/*             1111111111222222222233333333334444444444555555555566666666667
     01234567890123456789012345678901234567890123456789012345678901234567890 */
        "neph Na_.__K_.__Ca_.__[]_.__Al_.__Si_.__O16" };

    if (mask & FIRST) {
        char *string = strcpy((char *) malloc((size_t) (strlen(masterString)+1)*sizeof(char)), masterString);
        double totNa, totK, totCa, totVc, totAl, totSi;
        char n[5];
        int i;

        totNa = 4.0*(1.0-r[0]-r[1]-r[2]) + 3.0*r[1] + 2.0*r[2];
        totK  = 4.0*r[0];
        totCa = r[2];
        totVc = r[1] + r[2];
        totAl = 4.0*(1.0-r[0]-r[1]-r[2]) + 4.0*r[0] + 3.0*r[1] + 4.0*r[2];
        totSi = 4.0*(1.0-r[0]-r[1]-r[2]) + 4.0*r[0] + 5.0*r[1] + 4.0*r[2];

        (void) snprintf(n, 5, "%4.2f", totNa); for (i=0; i<4; i++) string[ 7+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totK);  for (i=0; i<4; i++) string[12+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totCa); for (i=0; i<4; i++) string[18+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totVc); for (i=0; i<4; i++) string[24+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totAl); for (i=0; i<4; i++) string[30+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totSi); for (i=0; i<4; i++) string[36+i] = n[i];

        *formula = string;
    }
}

void
actNph(int mask, double t, double p, double *x,
    double *a,  /* (pointer to a[]) activities              BINARY MASK: 0001 */
    double *mu, /* (pointer to mu[]) chemical potentials    BINARY MASK: 0010 */
    double **dx /* (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100 */
    )           /* exclusion criteria applied to results if BINARY MASK: 1000 */
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS], g, dgdr[NR];
    double fr[NA][NR];
    int i, j;

    for(i=0; i<NA; i++) {
     fr[i][0] = FR2(i); /* X2 */
     fr[i][1] = FR3(i); /* X3 */
     fr[i][2] = FR4(i); /* X4 */
    }

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    g       = G;
    dgdr[0] = DGDR0;
    dgdr[1] = DGDR1;
    dgdr[2] = DGDR2;

    /* activities for library */
    if (!mask && a != NULL) {
        for(i=0; i<NA; i++) {
       for (a[i]=g, j=0; j<NR; j++) a[i] += fr[i][j]*dgdr[j];
       a[i] = exp(a[i]/(R*t));
        }
    }

    if (mask & FIRST) {
        for(i=0; i<NA; i++) {
       for (a[i]=g, j=0; j<NR; j++) a[i] += fr[i][j]*dgdr[j];
       a[i] = exp(a[i]/(R*t));
        }
    }

    if (mask & SECOND) {
        for(i=0; i<NA; i++)
       for (mu[i]=g, j=0; j<NR; j++) mu[i] += fr[i][j]*dgdr[j];
    }

    if (mask & THIRD) {
        double d2gdr2[NR][NR], d2gdrds[NR][NS], d2gds2[NS][NS], dsdr[NS][NR],
            dfrdr[NA][NR], gs[NA][NS], dgsds[NA][NS], sum;
        int k, l;

        fillD2GDR2
        fillD2GDRDS
        fillD2GDS2

        for(i=0; i<NA; i++) {
       gs[i][0] = GS1(i); /* s1 */
       dfrdr[i][0] = DFR2DR2(i); /* X2 */
       dfrdr[i][1] = DFR3DR3(i); /* X3 */
       dfrdr[i][2] = DFR4DR4(i); /* X4 */
       dgsds[i][0] = DGS1DS1(i); /* s1 */
        }

        order(SECOND, t, p, r,
                    NULL, dsdr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NA; i++) {
            for (k=0; k<NR; k++) {
                /* compute activity of the i-th component */
                for (dx[i][k]=g, j=0; j<NR; j++) dx[i][k] += fr[i][j]*dgdr[j];
                dx[i][k] = exp(dx[i][k]/(R*t));

                /* compute derivative of i-th activity with respect to r(k) */
                sum = (1.0+dfrdr[i][k])*dgdr[k];
                for (j=0; j<NR; j++) {
                    sum += fr[i][j]*d2gdr2[j][k];
                    for (l=0; l<NS; l++) sum += fr[i][j]*d2gdrds[j][l]*dsdr[l][k];
                }
                for (j=0; j<NS; j++) {
                    sum += gs[i][j]*d2gdrds[k][j];
                    for (l=0; l<NS; l++) sum += gs[i][j]*d2gds2[j][l]*dsdr[l][k];
                }
                dx[i][k] *= sum/(R*t);
            }
        }
    }

    if (mask & FOURTH) {
        /* implement exclusion criteria on quantities for preclb routines         */
        static const double exclusion[NA] = {
       0.05,  /* exclusion criteria on the mole fraction of Na+  */
       0.05,  /* exclusion criteria on the mole fraction of K+   */
       0.05,  /* exclusion criteria on the mole fraction of Vc   */
       0.05,  /* exclusion criteria on the mole fraction of Ca++ */
        };
        double x[NA], totNa, totK, totCa, totVc;

        totNa = 4.0*(1.0-r[0]-r[1]-r[2]) + 3.0*r[1] + 2.0*r[2];
        totK  = 4.0*r[0];
        totCa = r[2];
        totVc = r[1] + r[2];

        x[0] = totNa/4.0;   /* Na+ averaged on both sites */
        x[1] = totK/4.0;    /* K+ averaged on both sites  */
        x[2] = totVc-totCa; /* Vc on large site           */
        x[3] = totCa;       /* Ca on small site           */

        for (i=0; i<NA; i++) {
            if (x[i] < exclusion[i]) {
                if (mask & FIRST)  a[i]  = 0.0;
                if (mask & SECOND) mu[i] = 0.0;
                if (mask & THIRD)  for (j=0; j<NR; j++) dx[i][j] = 0.0;
            }
        }
    }

}

void
gmixNph(int mask, double t, double p, double *x,
    double *gmix, /* Gibbs energy of mixing             BINARY MASK: 0001 */
    double *dx,   /* (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010 */
    double **dx2, /* (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100 */
    double ***dx3 /* (pointer to dx3[][][]) d3(g)/d(x[])3 NARY MASK: 1000 */
    )
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS];

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    if (mask & FIRST) {
        *gmix = G;
    }

    if(mask & SECOND) {
        dx[0] = DGDR0;
        dx[1] = DGDR1;
        dx[2] = DGDR2;
    }

    if(mask & THIRD) {
        double d2gdr2[NR][NR], d2gdrds[NR][NS], d2gds2[NS][NS], dsdr[NS][NR];
        int i, j, k, l;

        fillD2GDR2
        fillD2GDRDS
        fillD2GDS2

        order(SECOND, t, p, r,
                    NULL, dsdr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dx2[i][j] = d2gdr2[i][j];
                for (k=0; k<NS; k++) {
                    dx2[i][j] += d2gdrds[i][k]*dsdr[k][j] + d2gdrds[j][k]*dsdr[k][i];
                    for (l=0; l<NS; l++) dx2[i][j] += d2gds2[k][l]*dsdr[k][i]*dsdr[l][j];
                }
            }
        }

    }

    if(mask & FOURTH) {
        double d3gdr3[NR][NR][NR], d3gdr2ds[NR][NR][NS], d3gdrds2[NR][NS][NS];
        double d3gds3[NS][NS][NS], dsdr[NS][NR];
        int i, j, k, l, m, n;

        fillD3GDR3
        fillD3GDR2DS
        fillD3GDRDS2
        fillD3GDS3

        order(SECOND, t, p, r,
                    NULL, dsdr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                for (k=0; k<NR; k++) {
                    dx3[i][j][k] = d3gdr3[i][j][k];
                    for (l=0; l<NS; l++) {
                        dx3[i][j][k] += d3gdr2ds[i][j][l]*dsdr[l][k] +
                            d3gdr2ds[j][k][l]*dsdr[l][i] + d3gdr2ds[k][i][l]*dsdr[l][j];
                        for (m=0; m<NS; m++) {
                            dx3[i][j][k] +=
                                d3gdrds2[i][l][m]*dsdr[l][j]*dsdr[m][k] +
                                d3gdrds2[j][l][m]*dsdr[l][k]*dsdr[m][i] +
                                d3gdrds2[k][l][m]*dsdr[l][i]*dsdr[m][j];
                            for (n=0; n<NS; n++)
                                dx3[i][j][k] +=
                                    d3gds3[l][m][n]*dsdr[l][i]*dsdr[m][j]*dsdr[n][k];
                        }
                    }
                }
            }
        }

    }

}

void
hmixNph(int mask, double t, double p, double *x,
    double *hmix /* Enthalpy of mixing BINARY MASK: 1 */
    )
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS];

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    *hmix = (G) + t*(S);
}

void
smixNph(int mask, double t, double p, double *x,
    double *smix, /* Entropy of mixing                  BINARY MASK: 001 */
    double *dx,   /* (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010 */
    double **dx2  /* (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100 */
    )
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS];

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    if (mask & FIRST) {
        *smix = S;
    }

    if(mask & SECOND) {
        double d2gdrds[NR][NS], d2gdrdt[NR], d2gds2[NS][NS], d2gdsdt[NS],
            dsdr[NS][NR], dsdt[NS];
        int i, k, l;

        fillD2GDRDS
        fillD2GDRDT
        fillD2GDS2
        fillD2GDSDT

        order(SECOND | THIRD, t, p, r,
                    NULL, dsdr, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            dx[i] = d2gdrdt[i];
            for (k=0; k<NS; k++) {
                dx[i] += d2gdrds[i][k]*dsdt[k] + d2gdsdt[k]*dsdr[k][i];
                for (l=0; l<NS; l++) dx[i] += d2gds2[k][l]*dsdt[k]*dsdr[l][i] ;
            }
            dx[i] *= -1.0;
        }
    }

    if(mask & THIRD) {
        double d2gdrds[NR][NS], d2gds2[NS][NS], d2gdsdt[NS], d3gdr2ds[NR][NR][NS],
            d3gdr2dt[NR][NR], d3gdrds2[NR][NS][NS], d3gdrdsdt[NR][NS],
            d3gds3[NS][NS][NS], d3gds2dt[NS][NS], dsdr[NS][NR], dsdt[NS],
            d2sdr2[NS][NR][NR], d2sdrdt[NS][NR];
        int i, j, k, l, m;

        fillD2GDRDS
        fillD2GDS2
        fillD2GDSDT
        fillD3GDR2DS
        fillD3GDR2DT
        fillD3GDRDS2
        fillD3GDRDSDT
        fillD3GDS3
        fillD3GDS2DT

        order(SECOND | THIRD | FIFTH | SIXTH, t, p, r,
                    NULL, dsdr, dsdt, NULL, d2sdr2, d2sdrdt, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dx2[i][j] = d3gdr2dt[i][j];
                for (k=0; k<NS; k++) {
                    dx2[i][j] += d3gdr2ds[i][j][k]*dsdt[k]
                     + d3gdrdsdt[i][k]*dsdr[k][j]
                     + d3gdrdsdt[j][k]*dsdr[k][i]
                     + d2gdsdt[k]*d2sdr2[k][i][j]
                     + d2gdrds[i][k]*d2sdrdt[k][j]
                     + d2gdrds[j][k]*d2sdrdt[k][i];
                    for (l=0; l<NS; l++) {
                        dx2[i][j] += d3gdrds2[i][k][l]*dsdr[k][j]*dsdt[l]
                       + d3gdrds2[j][k][l]*dsdr[k][i]*dsdt[l]
                       + d2gds2[k][l]*d2sdr2[k][i][j]*dsdt[l]
                       + d3gds2dt[k][l]*dsdr[k][i]*dsdr[l][j]
                       + d2gds2[k][l]*dsdr[k][i]*d2sdrdt[l][j]
                       + d2gds2[k][l]*dsdr[k][j]*d2sdrdt[l][i];
                        for (m=0; m<NS; m++)
                            dx2[i][j] += d3gds3[k][l][m]*dsdr[k][i]*dsdr[l][j]*dsdt[m];
                    }
                }
                dx2[i][j] *= -1.0;
            }
        }

    }

}

void
cpmixNph(int mask, double t, double p, double *x,
    double *cpmix, /* Heat capacity of mixing               BINARY MASK: 001 */
    double *dt,    /* d(cp)/d(t)                            BINARY MASK: 010 */
    double *dx     /* d(cp)/d(x[])d(t)                      BINARY MASK: 100 */
    )
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS], dsdt[NS], d2gdsdt[NS], d2gds2[NS][NS], d2gdt2;
    int i, j;

    order(FIRST | THIRD, t, p, r,
                s, NULL, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    fillD2GDS2
    fillD2GDSDT
    d2gdt2  = D2GDT2;

    if (mask & FIRST) {

        *cpmix = d2gdt2;
        for (i=0; i<NS; i++) {
            *cpmix += 2.0*d2gdsdt[i]*dsdt[i];
            for (j=0; j<NS; j++) *cpmix += d2gds2[i][j]*dsdt[i]*dsdt[j];
        }
        *cpmix *= -t;
    }

    if(mask & SECOND) {
        double d3gds3[NS][NS][NS], d3gds2dt[NS][NS], d3gdsdt2[NS], d2sdt2[NS],
            temp;
        double d3gdt3 = D3GDT3;
        int k;

        fillD3GDS3
        fillD3GDS2DT
        fillD3GDSDT2

        order(EIGHTH, t, p, r,
                    NULL, NULL, NULL, NULL, NULL, NULL, NULL, d2sdt2, NULL, NULL);

        /* compute d2gdt2 */
        temp = d2gdt2;
        for (i=0; i<NS; i++) {
            temp += 2.0*d2gdsdt[i]*dsdt[i];
            for (j=0; j<NS; j++) temp += d2gds2[i][j]*dsdt[i]*dsdt[j];
        }

        *dt = d3gdt3;
        for (i=0; i<NS; i++) {
            *dt += 3.0*d3gdsdt2[i]*dsdt[i] + 3.0*d2gdsdt[i]*d2sdt2[i];
            for (j=0; j<NS; j++) {
                *dt += 3.0*d2gds2[i][j]*dsdt[i]*d2sdt2[j]
             + 3.0*d3gds2dt[i][j]*dsdt[i]*dsdt[j];
                for (k=0; k<NS; k++) *dt += d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdt[k];
            }
        }
        *dt = -t*(*dt) - temp;
    }

    if(mask & THIRD) {
        double d3gds3[NS][NS][NS], d3gdrds2[NR][NS][NS], d3gdrdsdt[NR][NS],
            d3gds2dt[NS][NS], d2gdrds[NR][NS], d3gdrdt2[NR], d3gdsdt2[NS],
            dsdr[NS][NR], d2sdrdt[NS][NR], d2sdt2[NS];
        int k, l;

        fillD2GDRDS
        fillD3GDRDS2
        fillD3GDRDSDT
        fillD3GDRDT2
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDSDT2

        order(SECOND | SIXTH | EIGHTH, t, p, r,
                    NULL, dsdr, NULL, NULL, NULL, d2sdrdt, NULL, d2sdt2, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0,dx[i]=d3gdrdt2[i]; j<NS; j++) {
                dx[i] += d3gdsdt2[j]*dsdr[j][i] + 2.0*d2gdsdt[j]*d2sdrdt[j][i] +
                 2.0*d3gdrdsdt[i][j]*dsdt[j] + d2gdrds[i][j]*d2sdt2[j];
                for (k=0; k<NS; k++) {
                    dx[i] += d3gdrds2[i][j][k]*dsdt[j]*dsdt[k] +
                   2.0*d2gds2[j][k]*dsdt[j]*d2sdrdt[k][i] +
                   2.0*d3gds2dt[j][k]*dsdr[j][i]*dsdt[k] +
                   d2gds2[j][k]*dsdr[j][i]*d2sdt2[k];
                    for (l=0; l<NS; l++)
                        dx[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdt[k]*dsdt[l];
                }
            }
            dx[i] *= -t;
        }
    }

}

void
vmixNph(int mask, double t, double p, double *x,
    double *vmix, /* Volume of mixing                BINARY MASK: 0000000001 */
    double *dx,   /* (pointer to dx[]) d(v)/d(x[])   BINARY MASK: 0000000010 */
    double **dx2, /* (point to dx2[][]) d(v)/d(x[])2 BINARY MASK: 0000000100 */
    double *dt,   /* d(v)/d(t)                       BINARY MASK: 0000001000 */
    double *dp,   /* d(v)/d(p)                       BINARY MASK: 0000010000 */
    double *dt2,  /* d2(v)/d(t)2                     BINARY MASK: 0000100000 */
    double *dtdp, /* d2(v)/d(t)d(p)                  BINARY MASK: 0001000000 */
    double *dp2,  /* d2(v)/d(p)2                     BINARY MASK: 0010000000 */
    double *dxdt, /* d2(v)/d(x[])d(t)                BINARY MASK: 0100000000 */
    double *dxdp  /* d2(v)/d(x[])d(p)                BINARY MASK: 1000000000 */
    )
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS];

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    if (mask & FIRST) {
        *vmix = DGDP;
    }

    if(mask & SECOND) {
        double d2gdrds[NR][NS], d2gdrdp[NR], d2gds2[NS][NS], d2gdsdp[NS],
            dsdr[NS][NR], dsdp[NS];
        int i, j, k;

        fillD2GDRDS
        fillD2GDRDP
        fillD2GDS2
        fillD2GDSDP

        order(SECOND | FOURTH, t, p, r,
                    NULL, dsdr, NULL, dsdp, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            dx[i] = d2gdrdp[i];
            for (j=0; j<NS; j++) {
                dx[i] += d2gdrds[i][j]*dsdp[j] + d2gdsdp[j]*dsdr[j][i];
                for (k=0; k<NS; k++) dx[i] += d2gds2[j][k]*dsdp[j]*dsdr[k][i];
            }
        }
    }

    if(mask & THIRD) {
        double d2gdrds[NR][NS], d2gds2[NS][NS], d2gdsdp[NS], d3gdr2ds[NR][NR][NS],
            d3gdr2dp[NR][NR], d3gdrds2[NR][NS][NS], d3gdrdsdp[NR][NS],
            d3gds3[NS][NS][NS], d3gds2dp[NS][NS], dsdr[NS][NR], dsdp[NS],
            d2sdr2[NS][NR][NR], d2sdrdp[NS][NR];
        int i, j, k, l, m;

        fillD2GDRDS
        fillD2GDS2
        fillD2GDSDP
        fillD3GDR2DS
        fillD3GDR2DP
        fillD3GDRDS2
        fillD3GDRDSDP
        fillD3GDS3
        fillD3GDS2DP

        order(SECOND | FOURTH | FIFTH | SEVENTH, t, p, r,
                    NULL, dsdr, NULL, dsdp, d2sdr2, NULL, d2sdrdp,  NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dx2[i][j] = d3gdr2dp[i][j];
                for (k=0; k<NS; k++) {
                    dx2[i][j] += d3gdr2ds[i][j][k]*dsdp[k]
                     + d3gdrdsdp[i][k]*dsdr[k][j]
                     + d3gdrdsdp[j][k]*dsdr[k][i]
                     + d2gdsdp[k]*d2sdr2[k][i][j]
                     + d2gdrds[i][k]*d2sdrdp[k][j]
                     + d2gdrds[j][k]*d2sdrdp[k][i];
                    for (l=0; l<NS; l++) {
                        dx2[i][j] += d3gdrds2[i][k][l]*dsdr[k][j]*dsdp[l]
                       + d3gdrds2[j][k][l]*dsdr[k][i]*dsdp[l]
                       + d2gds2[k][l]*d2sdr2[k][i][j]*dsdp[l]
                       + d3gds2dp[k][l]*dsdr[k][i]*dsdr[l][j]
                       + d2gds2[k][l]*dsdr[k][i]*d2sdrdp[l][j]
                       + d2gds2[k][l]*dsdr[k][j]*d2sdrdp[l][i];
                        for (m=0; m<NS; m++)
                            dx2[i][j] += d3gds3[k][l][m]*dsdr[k][i]*dsdr[l][j]*dsdp[m];
                    }
                }
            }
        }

    }

    if(mask & FOURTH) {
        double d2gdtdp = D2GDTDP;
        double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], dsdt[NS], dsdp[NS];
        int i, j;

        fillD2GDS2
        fillD2GDSDT
        fillD2GDSDP

        order(THIRD | FOURTH, t, p, r,
                    NULL, NULL, dsdt, dsdp, NULL, NULL, NULL, NULL, NULL, NULL);

        *dt = d2gdtdp;
        for (i=0; i<NS; i++) {
            *dt += d2gdsdt[i]*dsdp[i] + d2gdsdp[i]*dsdt[i];
            for (j=0; j<NS; j++) *dt += d2gds2[i][j]*dsdt[i]*dsdp[j];
        }
    }

    if(mask & FIFTH) {
        double d2gdp2 = D2GDP2;
        double d2gds2[NS][NS], d2gdsdp[NS], dsdp[NS];
        int i,j;

        fillD2GDS2
        fillD2GDSDP

        order(FOURTH, t, p, r,
                    NULL, NULL, NULL, dsdp, NULL, NULL, NULL, NULL, NULL, NULL);

        *dp = d2gdp2;
        for (i=0; i<NS; i++) {
            *dp += 2.0*d2gdsdp[i]*dsdp[i];
            for (j=0; j<NS; j++) *dp += d2gds2[i][j]*dsdp[i]*dsdp[j];
        }
    }

    if(mask & SIXTH) {
        double d3gdt2dp = D3GDT2DP;
        double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS][NS][NS],
            d3gds2dp[NS][NS], d3gds2dt[NS][NS], d3gdsdtdp[NS], d3gdsdt2[NS],
            dsdt[NS], dsdp[NS], d2sdt2[NS], d2sdtdp[NS];
        int i, j, k;

        fillD2GDS2
        fillD2GDSDT
        fillD2GDSDP
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDS2DP
        fillD3GDSDT2
        fillD3GDSDTDP

        order(THIRD | FOURTH | EIGHTH | NINTH, t, p, r,
                    NULL, NULL, dsdt, dsdp, NULL, NULL, NULL, d2sdt2, d2sdtdp, NULL);

        *dt2 = d3gdt2dp;
        for (i=0; i<NS; i++) {
            *dt2 += d3gdsdt2[i]*dsdp[i] + 2.0*d2gdsdt[i]*d2sdtdp[i]
                        + d2gdsdp[i]*d2sdt2[i] + 2.0*d3gdsdtdp[i]*dsdt[i];
            for (j=0; j<NS; j++) {
                *dt2 += 2.0*d3gds2dt[i][j]*dsdt[i]*dsdp[j]
                            + d2gds2[i][j]*d2sdt2[i]*dsdp[j]
                            + 2.0*d2gds2[i][j]*dsdt[i]*d2sdtdp[j]
                            + d3gds2dp[i][j]*dsdt[i]*dsdt[j];
                for (k=0; k<NS; k++) *dt2 += d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdp[k];
            }
        }
    }

    if(mask & SEVENTH) {
        double d3gdtdp2 = D3GDTDP2;
        double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS][NS][NS],
            d3gds2dt[NS][NS], d3gds2dp[NS][NS], d3gdsdtdp[NS], d3gdsdp2[NS],
            dsdt[NS], dsdp[NS], d2sdtdp[NS], d2sdp2[NS];
        int i, j, k;

        fillD2GDS2
        fillD2GDSDT
        fillD2GDSDP
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDS2DP
        fillD3GDSDTDP
        fillD3GDSDP2

        order(THIRD | FOURTH | NINTH | TENTH, t, p, r,
                    NULL, NULL, dsdt, dsdp, NULL, NULL, NULL, NULL, d2sdtdp, d2sdp2);

        *dtdp = d3gdtdp2;
        for (i=0; i<NS; i++) {
            *dtdp += 2.0*d3gdsdtdp[i]*dsdp[i] + d2gdsdt[i]*d2sdp2[i]
             + 2.0*d2gdsdp[i]*d2sdtdp[i] + d3gdsdp2[i]*dsdt[i];
            for (j=0; j<NS; j++) {
                *dtdp += 2.0*d3gds2dp[i][j]*dsdt[i]*dsdp[j]
               + d2gds2[i][j]*dsdt[i]*d2sdp2[j]
               + 2.0*d2gds2[i][j]*d2sdtdp[i]*dsdp[j]
               + d3gds2dt[i][j]*dsdp[i]*dsdp[j];
                for (k=0; k<NS; k++) *dtdp += d3gds3[i][j][k]*dsdt[i]*dsdp[j]*dsdp[k];
            }
        }
    }

    if(mask & EIGHTH) {
        double d3gdp3 = D3GDP3;
        double d2gds2[NS][NS], d2gdsdp[NS], d3gds3[NS][NS][NS], d3gds2dp[NS][NS],
            d3gdsdp2[NS], dsdp[NS], d2sdp2[NS];
        int i, j, k;

        fillD2GDS2
        fillD2GDSDP
        fillD3GDS3
        fillD3GDS2DP
        fillD3GDSDP2

        order(FOURTH | TENTH, t, p, r,
                    NULL, NULL, NULL, dsdp, NULL, NULL, NULL, NULL, NULL, d2sdp2);

        *dp2 = d3gdp3;
        for (i=0; i<NS; i++) {
            *dp2 += 3.0*d3gdsdp2[i]*dsdp[i] + 3.0*d2gdsdp[i]*d2sdp2[i];
            for (j=0; j<NS; j++) {
                *dp2 += 3.0*d2gds2[i][j]*dsdp[i]*d2sdp2[j]
                            + 3.0*d3gds2dp[i][j]*dsdp[i]*dsdp[j];
                for (k=0; k<NS; k++) *dp2 += d3gds3[i][j][k]*dsdp[i]*dsdp[j]*dsdp[k];
            }
        }
    }

    if(mask & NINTH) {
        double d3gds3[NS][NS][NS], d3gdrds2[NR][NS][NS], d3gdrdsdt[NR][NS],
            d3gds2dp[NS][NS], d2gdrds[NR][NS], d3gdrdtdp[NR], d3gdsdtdp[NS],
            dsdt[NS], dsdp[NS], dsdr[NS][NR], d2sdrdt[NS][NR], d2sdrdp[NS][NR],
            d2gds2[NS][NS], d2gdsdt[NS], d3gdrdsdp[NR][NS], d2gdsdp[NS],
            d2sdtdp[NS], d3gds2dt[NS][NS];
        int i, j, k, l;

        fillD2GDRDS
        fillD2GDS2
        fillD2GDSDT
        fillD2GDSDP
        fillD3GDRDS2
        fillD3GDRDSDT
        fillD3GDRDSDP
        fillD3GDRDTDP
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDSDTDP
        fillD3GDS2DP

        order(SECOND | THIRD | FOURTH | SIXTH | SEVENTH | NINTH, t, p, r,
            NULL, dsdr, dsdt, dsdp, NULL, d2sdrdt, d2sdrdp, NULL, d2sdtdp, NULL);

        for (i=0; i<NR; i++) {
            for (j=0,dxdt[i]=d3gdrdtdp[i]; j<NS; j++) {
                dxdt[i] += d3gdsdtdp[j]*dsdr[j][i] + d2gdsdt[j]*d2sdrdp[j][i] +
                   d3gdrdsdt[i][j]*dsdp[j] + d2gdrds[i][j]*d2sdtdp[j] +
                   d3gdrdsdp[i][j]*dsdt[j] + d2gdsdp[j]*d2sdrdt[j][i];
                for (k=0; k<NS; k++) {
                    dxdt[i] += d3gdrds2[i][j][k]*dsdt[j]*dsdp[k] +
                     d2gds2[j][k]*dsdt[j]*d2sdrdp[k][i] +
                     d2gds2[j][k]*dsdp[j]*d2sdrdt[k][i] +
                     d3gds2dt[j][k]*dsdr[j][i]*dsdp[k] +
                     d3gds2dp[j][k]*dsdr[j][i]*dsdt[k] +
                     d2gds2[j][k]*dsdr[j][i]*d2sdtdp[k];
                    for (l=0; l<NS; l++)
                        dxdt[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdt[k]*dsdp[l];
                }
            }
        }
    }

    if(mask & TENTH) {
        double d3gds3[NS][NS][NS], d3gdrds2[NR][NS][NS], d3gds2dp[NS][NS],
            d2gdrds[NR][NS], dsdp[NS], dsdr[NS][NR], d2sdrdp[NS][NR], d2gds2[NS][NS],
            d3gdrdsdp[NR][NS], d3gdrdp2[NR], d3gdsdp2[NS], d2gdsdp[NS], d2sdp2[NS];
        int i, j, k, l;

        fillD2GDRDS
        fillD2GDS2
        fillD2GDSDP
        fillD3GDRDS2
        fillD3GDRDSDP
        fillD3GDRDP2
        fillD3GDS3
        fillD3GDS2DP
        fillD3GDSDP2

        order(SECOND | FOURTH | SEVENTH | TENTH, t, p, r,
            NULL, dsdr, NULL, dsdp, NULL, NULL, d2sdrdp, NULL, NULL, d2sdp2);

        for (i=0; i<NR; i++) {
            for (j=0,dxdp[i]=d3gdrdp2[i]; j<NS; j++) {
                dxdp[i] += d3gdsdp2[j]*dsdr[j][i] + d2gdsdp[j]*d2sdrdp[j][i] +
                   2.0*d3gdrdsdp[i][j]*dsdp[j] + d2gdrds[i][j]*d2sdp2[j] +
                   d2gdsdp[j]*d2sdrdp[j][i];
                for (k=0; k<NS; k++) {
                    dxdp[i] += d3gdrds2[i][j][k]*dsdp[j]*dsdp[k] +
                     2.0*d2gds2[j][k]*dsdp[j]*d2sdrdp[k][i] +
                     2.0*d3gds2dp[j][k]*dsdr[j][i]*dsdp[k] +
                     d2gds2[j][k]*dsdr[j][i]*d2sdp2[k];
                    for (l=0; l<NS; l++)
                        dxdp[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdp[k]*dsdp[l];
                }
            }
        }
    }

}

/* end of file NEPHELINE.C */
