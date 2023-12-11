const char *cummingtonite_ver(void) { return "$Id: cummingtonite.c,v 1.3 2006/10/20 00:59:22 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: cummingtonite.c,v $
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
MELTS Source Code: RCS Revision 1.1.1.1  2006/08/15 16:57:35  ghiorso
MELTS Source Code: RCS xMELTS gcc 3.x sources
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2004/01/02 19:21:49  cvsaccount
MELTS Source Code: RCS CTserver University of Chicago
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2003/05/01 17:33:54  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:58:12  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
 * Revision 3.7  1997/06/21  22:49:59  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.6  1997/05/03  20:23:38  ghiorso
 * *** empty log message ***
 *
 * Revision 3.5  1997/03/27  17:03:41  ghiorso
 * *** empty log message ***
 *
 * Revision 3.4  1996/09/24  20:33:44  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 3.3  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 3.2  1995/09/09  23:37:35  ghiorso
 * Modifications by Asimow to include new Derivatives of gmix and convert
 * options for liquid-absent fO2 buffering. Also new olivine components
 * implemented.
 *
 * Revision 3.2  1995/09/09  23:37:35  ghiorso
 * Modifications by Asimow to include new Derivatives of gmix and convert
 * options for liquid-absent fO2 buffering. Also new olivine components
 * implemented.
 *
 * Revision 3.1  1995/08/18  17:28:52  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Routines to compute cummingtonite solution properties
**      (file: CUMMINGTONITE.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  April 16, 1991   Original Version
**      V1.0-2  Mark S. Ghiorso  April 17, 1992
**              (1) Modified bound constraint algorithm in the calculation
**                  of the ordering parameter
**      V1.1-1  Mark S. Ghiorso  July 13, 1992
**              Added (*display) function
**      V2.0-1  Mark S. Ghiorso  May 9, 1994
**              Updated three site formulation and added revised solution
**              model
**      V3.0-1  Mark S. Ghiorso  May 10, 1994
**              (1) Began modifications for new isenthalpic, isentropic,
**                  isochoric derivatives
**      V3.0-2  Mark S. Ghiorso  May 17, 1994
**              (1) Changed definition of the enthalpy of mixing
**      V3.1-1  Mark S. Ghiorso  February 28, 1995
**              (1) Updated mixing parameters to conform to Ghiorso et al.
**                  (1995) (Note: also updated sol_struct_data.h)
**	V4.0-1  Paul D. Asimow  July 27, 1995
**		(1) added d3rdm3 to (*convert) and d3gdr3 to (*gmix)
**--
*/

#include "melts_gsl.h"
#include "silmin.h"  /* Structure definitions foor SILMIN package */

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))

#ifdef PRINT_NONCONVERGENCE_IN_ORDER
#undef PRINT_NONCONVERGENCE_IN_ORDER
#endif

#ifdef DEBUG
#undef DEBUG
#endif

/*
 *=============================================================================
 * Cummingtonite solution parameters:
 * Ghiorso, M.S., Evans, B.W., Hirschmann, M, Yang H. (1994)
 *   Thermodynamics of the Amphiboles: I. (Fe2+,Mg) Cummingtonite Solid
 *   Solutions. American Mineralogist (submitted)
 */

#define HEX     (-25856.0  ) /* J - fit to the 3-site model */
#define HDEX    (  5905.0  )
#define WH123   (  5115.64 )
#define WH13    (  4592.11 )
#define WH2     (  4324.28 )
#define WH4     (  3559.52 )
#define HX123_4 ( 12102.8  )
#define HX134_2 ( 15113.3  )
#define HX13_24 (  3296.10 )

#define VEX     0.0      /* J/bar */
#define VDEX    0.0
#define WV123   0.0
#define WV13    0.0
#define WV2     0.0
#define WV4     0.0
#define VX123_4 0.0
#define VX134_2 0.0
#define VX13_24 0.0

#define GEX     (HEX     + (p-1.0)*VEX)
#define GDEX    (HDEX    + (p-1.0)*VDEX)
#define WG123   (WH123   + (p-1.0)*WV123)
#define WG13    (WH13    + (p-1.0)*WV13)
#define WG2     (WH2     + (p-1.0)*WV2)
#define WG4     (WH4     + (p-1.0)*WV4)
#define GX123_4 (HX123_4 + (p-1.0)*VX123_4)
#define GX134_2 (HX134_2 + (p-1.0)*VX134_2)
#define GX13_24 (HX13_24 + (p-1.0)*VX13_24)


/*
 * Definitions of Taylor expansion coefficients in terms of solution
 * parameters.
 */

#define H0      (5.0*WH123/4.0+WH4/2.0+HX123_4/4.0)
#define HR      0.0
#define HS13    (3.0*HEX/5.0+HDEX)
#define HS2     (2.0*HEX/5.0-HDEX)
#define HRR     (-5.0*WH123/4.0-WH4/2.0-HX123_4/4.0)
#define HRS13   (-15.0*WH123/7.0+3.0*WH13-6.0*WH4/7.0-3.0*HX123_4/7.0 \
                 +HX13_24/2.0)
#define HRS2    (-10.0*WH123/7.0+2.0*WH2-4.0*WH4/7.0-2.0*HX123_4/7.0 \
                 +HX134_2/2.0)
#define HS13S13 (-45.0*WH123/49.0-3.0*WH13/7.0-18.0*WH4/49.0-9.0*HX123_4/49.0 \
                 +3.0*HX13_24/7.0)
#define HS13S2  (-60.0*WH123/49.0+12.0*WH13/7.0+12.0*WH2/7.0-24.0*WH4/49.0 \
                 +25.0*HX123_4/98.0-3.0*HX13_24/14.0-HX134_2/14.0)
#define HS2S2   (-20.0*WH123/49.0-6.0*WH2/7.0-8.0*WH4/49.0-4.0*HX123_4/49.0 \
                 +2.0*HX134_2/7.0)

#define V0      (5.0*WV123/4.0+WV4/2.0+VX123_4/4.0)
#define VR      0.0
#define VS13    (3.0*VEX/5.0+VDEX)
#define VS2     (2.0*VEX/5.0-VDEX)
#define VRR     (-5.0*WV123/4.0-WV4/2.0-VX123_4/4.0)
#define VRS13   (-15.0*WV123/7.0+3.0*WV13-6.0*WV4/7.0-3.0*VX123_4/7.0 \
                 +VX13_24/2.0)
#define VRS2    (-10.0*WV123/7.0+2.0*WV2-4.0*WV4/7.0-2.0*VX123_4/7.0 \
                 +VX134_2/2.0)
#define VS13S13 (-45.0*WV123/49.0-3.0*WV13/7.0-18.0*WV4/49.0-9.0*VX123_4/49.0 \
                 +3.0*VX13_24/7.0)
#define VS13S2  (-60.0*WV123/49.0+12.0*WV13/7.0+12.0*WV2/7.0-24.0*WV4/49.0 \
                 +25.0*VX123_4/98.0-3.0*VX13_24/14.0-VX134_2/14.0)
#define VS2S2   (-20.0*WV123/49.0-6.0*WV2/7.0-8.0*WV4/49.0-4.0*VX123_4/49.0 \
                 +2.0*VX134_2/7.0)

#define G0      (5.0*WG123/4.0+WG4/2.0+GX123_4/4.0)
#define GR      0.0
#define GS13    (3.0*GEX/5.0+GDEX)
#define GS2     (2.0*GEX/5.0-GDEX)
#define GRR     (-5.0*WG123/4.0-WG4/2.0-GX123_4/4.0)
#define GRS13   (-15.0*WG123/7.0+3.0*WG13-6.0*WG4/7.0-3.0*GX123_4/7.0 \
                 +GX13_24/2.0)
#define GRS2    (-10.0*WG123/7.0+2.0*WG2-4.0*WG4/7.0-2.0*GX123_4/7.0 \
                 +GX134_2/2.0)
#define GS13S13 (-45.0*WG123/49.0-3.0*WG13/7.0-18.0*WG4/49.0-9.0*GX123_4/49.0 \
                 +3.0*GX13_24/7.0)
#define GS13S2  (-60.0*WG123/49.0+12.0*WG13/7.0+12.0*WG2/7.0-24.0*WG4/49.0 \
                 +25.0*GX123_4/98.0-3.0*GX13_24/14.0-GX134_2/14.0)
#define GS2S2   (-20.0*WG123/49.0-6.0*WG2/7.0-8.0*WG4/49.0-4.0*GX123_4/49.0 \
                 +2.0*GX134_2/7.0)

/*
 * Global (to this file): variables
 */

#define R   8.3143
#define NR  1           /* One independent composition variable */
#define NS  2           /* Two ordering parameters              */
#define NA  2           /* Two endmember compositions           */

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
    double xfe2m13, xmg2m13, xfe2m2, xmg2m2, xfe2m4, xmg2m4;

#define XFE2M13 0
#define XMG2M13 1
#define XFE2M2  2
#define XMG2M2  3
#define XFE2M4  4
#define XMG2M4  5

#define NX      6

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
    xfe2m13 = getX(XFE2M13); \
    xmg2m13 = getX(XMG2M13); \
    xfe2m2  = getX(XFE2M2); \
    xmg2m2  = getX(XMG2M2); \
    xfe2m4  = getX(XFE2M4); \
    xmg2m4  = getX(XMG2M4);

#define SET_SITE_FRACTIONS \
    setX(XFE2M13, xfe2m13); \
    setX(XMG2M13, xmg2m13); \
    setX(XFE2M2,  xfe2m2); \
    setX(XMG2M2,  xmg2m2); \
    setX(XFE2M4,  xfe2m4); \
    setX(XMG2M4,  xmg2m4);

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conCum defines the conversion from m[i], to r[j]
 */

#define FR0(i)        (i == 0) ? -(1.0+r[0]) : (1.0-r[0])

#define GS0(i)       - s[0]
#define GS1(i)       - s[1]

#define DFR0DR0(i)   - 1.0

#define DGS0DS0(i)   - 1.0
#define DGS1DS1(i)   - 1.0

/*
 * Global (to this file): derivative definitions
 */

#define S  -R*(2.0*xfe2m4 *log(xfe2m4)  + 2.0*xmg2m4 *log(xmg2m4)  \
             + 3.0*xfe2m13*log(xfe2m13) + 3.0*xmg2m13*log(xmg2m13) \
             + 2.0*xfe2m2 *log(xfe2m2)  + 2.0*xmg2m2 *log(xmg2m2)  )
#define H  H0 + HR*r[0] + HS13*s[0] + HS2*s[1] + HRR*r[0]*r[0] \
           + HRS13*r[0]*s[0] + HRS2*r[0]*s[1] + HS13S13*s[0]*s[0] \
           + HS13S2*s[0]*s[1] + HS2S2*s[1]*s[1]
#define V  V0 + VR*r[0] + VS13*s[0] + VS2*s[1] + VRR*r[0]*r[0] \
           + VRS13*r[0]*s[0] + VRS2*r[0]*s[1] + VS13S13*s[0]*s[0] \
           + VS13S2*s[0]*s[1] + VS2S2*s[1]*s[1]
#define G  H - t*(S) + (p-1.0)*(V)

/*----------------------------------------------------------------------------*/

#define DGDR0  R*t*(      log(xfe2m4)      -     log(xmg2m4) \
                                        + 3.0*log(xfe2m13)/2.0 - 3.0*log(xmg2m13)/2.0 \
                                        +     log(xfe2m2)      -     log(xmg2m2) ) \
               + GR + 2.0*GRR*r[0] + GRS13*s[0] + GRS2*s[1]
#define DGDS0  R*t*(   6.0*log(xfe2m4)/7.0  -  6.0*log(xmg2m4)/7.0  \
                                        - 12.0*log(xfe2m13)/7.0 + 12.0*log(xmg2m13)/7.0 \
                                        +  6.0*log(xfe2m2)/7.0  -  6.0*log(xmg2m2)/7.0 ) \
               + GS13 + GRS13*r[0] + 2.0*GS13S13*s[0] + GS13S2*s[1]
#define DGDS1  R*t*(   4.0*log(xfe2m4)/7.0  -  4.0*log(xmg2m4)/7.0  \
                                        +  6.0*log(xfe2m13)/7.0 -  6.0*log(xmg2m13)/7.0 \
                                        - 10.0*log(xfe2m2)/7.0  + 10.0*log(xmg2m2)/7.0 ) \
               + GS2 + GRS2*r[0] + GS13S2*s[0] + 2.0*GS2S2*s[1]
#define DGDT   -(S)
#define DGDP   (V)

/*----------------------------------------------------------------------------*/

#define D2GDR0R0  R*t*(0.5/xfe2m4 + 0.5/xmg2m4 + 0.75/xfe2m13 + 0.75/xmg2m13 \
                       + 0.5/xfe2m2 + 0.5/xmg2m2) + 2.0*GRR
#define D2GDR0S0  R*t*(  (3.0/7.0)/xfe2m4  + (3.0/7.0)/xmg2m4  \
                       - (6.0/7.0)/xfe2m13 - (6.0/7.0)/xmg2m13 \
                       + (3.0/7.0)/xfe2m2  + (3.0/7.0)/xmg2m2 ) + GRS13
#define D2GDR0S1  R*t*(  (2.0/7.0)/xfe2m4  + (2.0/7.0)/xmg2m4  \
                       + (3.0/7.0)/xfe2m13 + (3.0/7.0)/xmg2m13 \
                       - (5.0/7.0)/xfe2m2  - (5.0/7.0)/xmg2m2 ) + GRS2
#define D2GDR0DT  R*(      log(xfe2m4)      -     log(xmg2m4) \
                     + 3.0*log(xfe2m13)/2.0 - 3.0*log(xmg2m13)/2.0 \
                     +     log(xfe2m2)      -     log(xmg2m2) )
#define D2GDR0DP  VR + 2.0*VRR*r[0] + VRS13*s[0] + VRS2*s[1]

#define D2GDS0S0  R*t*(  (18.0/49.0)/xfe2m4  + (18.0/49.0)/xmg2m4  \
                       + (48.0/49.0)/xfe2m13 + (48.0/49.0)/xmg2m13 \
                       + (18.0/49.0)/xfe2m2  + (18.0/49.0)/xmg2m2 ) \
                                    + 2.0*GS13S13
#define D2GDS0S1  R*t*(  (12.0/49.0)/xfe2m4  + (12.0/49.0)/xmg2m4  \
                       - (24.0/49.0)/xfe2m13 - (24.0/49.0)/xmg2m13 \
                       - (30.0/49.0)/xfe2m2  - (30.0/49.0)/xmg2m2 ) \
                                    + GS13S2
#define D2GDS0DT  R*(   6.0*log(xfe2m4)/7.0  -  6.0*log(xmg2m4)/7.0  \
                     - 12.0*log(xfe2m13)/7.0 + 12.0*log(xmg2m13)/7.0 \
                     +  6.0*log(xfe2m2)/7.0  -  6.0*log(xmg2m2)/7.0 )
#define D2GDS0DP  VS13 + VRS13*r[0] + 2.0*VS13S13*s[0] + VS13S2*s[1]

#define D2GDS1S1  R*t*(  (8.0/49.0)/xfe2m4   + (8.0/49.0)/xmg2m4   \
                       + (12.0/49.0)/xfe2m13 + (12.0/49.0)/xmg2m13 \
                       + (50.0/49.0)/xfe2m2  + (50.0/49.0)/xmg2m2 ) \
                                    + 2.0*GS2S2
#define D2GDS1DT  R*(   4.0*log(xfe2m4)/7.0  -  4.0*log(xmg2m4)/7.0  \
                     +  6.0*log(xfe2m13)/7.0 -  6.0*log(xmg2m13)/7.0 \
                     - 10.0*log(xfe2m2)/7.0  + 10.0*log(xmg2m2)/7.0 )
#define D2GDS1DP  VS2 + VRS2*r[0] + VS13S2*s[0] + 2.0*VS2S2*s[1]

#define D2GDT2    0.0
#define D2GDTDP   0.0
#define D2GDP2    0.0

/*----------------------------------------------------------------------------*/

#define D3GDR0R0R0  R*t*(-0.25/SQUARE(xfe2m2)+0.25/SQUARE(xmg2m2)- \
            0.375/SQUARE(xfe2m13)+0.375/SQUARE(xmg2m13)- \
            0.25/SQUARE(xfe2m4)+0.25/SQUARE(xmg2m4))
#define D3GDR0R0S0  -R*t*((3.0/14.0)/SQUARE(xfe2m4)-(3.0/14.0)/SQUARE(xmg2m4) \
                                        - (12.0/28.0)/SQUARE(xfe2m13)+(12.0/28.0)/SQUARE(xmg2m13) \
                                        + (3.0/14.0)/SQUARE(xfe2m2) - (3.0/14.0)/SQUARE(xmg2m2) )
#define D3GDR0R0S1  -R*t*((1.0/7.0)/SQUARE(xfe2m4)-(1.0/7.0)/SQUARE(xmg2m4) \
                                        + (6.0/28.0)/SQUARE(xfe2m13)-(6.0/28.0)/SQUARE(xmg2m13) \
                                        - (5.0/14.0)/SQUARE(xfe2m2) + (5.0/14.0)/SQUARE(xmg2m2) )
#define D3GDR0R0DT  R*(0.5/xfe2m4 + 0.5/xmg2m4 + 0.75/xfe2m13 + 0.75/xmg2m13 \
                       + 0.5/xfe2m2 + 0.5/xmg2m2)
#define D3GDR0R0DP  2.0*VRR

#define D3GDR0S0S0  -R*t*((9.0/49.0)/SQUARE(xfe2m4)-(9.0/49.0)/SQUARE(xmg2m4) \
                                        + (24.0/49.0)/SQUARE(xfe2m13)-(24.0/49.0)/SQUARE(xmg2m13) \
                                        + (9.0/49.0)/SQUARE(xfe2m2)-(9.0/49.0)/SQUARE(xmg2m2) )
#define D3GDR0S0S1  -R*t*((6.0/49.0)/SQUARE(xfe2m4)-(6.0/49.0)/SQUARE(xmg2m4) \
                                        - (12.0/49.0)/SQUARE(xfe2m13)+(12.0/49.0)/SQUARE(xmg2m13) \
                                        - (15.0/49.0)/SQUARE(xfe2m2)+(15.0/49.0)/SQUARE(xmg2m2) )
#define D3GDR0S0DT  R*(  (3.0/7.0)/xfe2m4  + (3.0/7.0)/xmg2m4  \
                       - (6.0/7.0)/xfe2m13 - (6.0/7.0)/xmg2m13 \
                       + (3.0/7.0)/xfe2m2  + (3.0/7.0)/xmg2m2 )
#define D3GDR0S0DP  VRS13

#define D3GDR0S1S1  -R*t*((4.0/49.0)/SQUARE(xfe2m4)-(4.0/49.0)/SQUARE(xmg2m4) \
                                        + (6.0/49.0)/SQUARE(xfe2m13)-(6.0/49.0)/SQUARE(xmg2m13) \
                                        + (25.0/49.0)/SQUARE(xfe2m2)-(25.0/49.0)/SQUARE(xmg2m2) )
#define D3GDR0S1DT  R*(  (2.0/7.0)/xfe2m4  + (2.0/7.0)/xmg2m4  \
                       + (3.0/7.0)/xfe2m13 + (3.0/7.0)/xmg2m13 \
                       - (5.0/7.0)/xfe2m2  - (5.0/7.0)/xmg2m2 )
#define D3GDR0S1DP  VRS2

#define D3GDS0S0S0  -R*t*(  (54.0/343.0)/SQUARE(xfe2m4)   \
                                                    - (54.0/343.0)/SQUARE(xmg2m4)   \
                                                    - (192.0/343.0)/SQUARE(xfe2m13) \
                                                    + (192.0/343.0)/SQUARE(xmg2m13) \
                                                    + (54.0/343.0)/SQUARE(xfe2m2)   \
                                                    - (54.0/343.0)/SQUARE(xmg2m2) )
#define D3GDS0S0S1  -R*t*(  (36.0/343.0)/SQUARE(xfe2m4)   \
                                                    - (36.0/343.0)/SQUARE(xmg2m4)   \
                                                    + (96.0/343.0)/SQUARE(xfe2m13)  \
                                                    - (96.0/343.0)/SQUARE(xmg2m13)  \
                                                    - (90.0/343.0)/SQUARE(xfe2m2)   \
                                                    + (90.0/343.0)/SQUARE(xmg2m2) )
#define D3GDS0S0DT  R*(  (18.0/49.0)/xfe2m4  + (18.0/49.0)/xmg2m4  \
                       + (48.0/49.0)/xfe2m13 + (48.0/49.0)/xmg2m13 \
                       + (18.0/49.0)/xfe2m2  + (18.0/49.0)/xmg2m2 )
#define D3GDS0S0DP  2.0*VS13S13
#define D3GDS0S1S1  -R*t*(  (24.0/343.0)/SQUARE(xfe2m4)  \
                                                    - (24.0/343.0)/SQUARE(xmg2m4)  \
                                                    - (48.0/343.0)/SQUARE(xfe2m13) \
                                                    + (48.0/343.0)/SQUARE(xmg2m13) \
                                                    + (150.0/343.0)/SQUARE(xfe2m2) \
                                                    - (150.0/343.0)/SQUARE(xmg2m2) )
#define D3GDS0S1DT  R*(  (12.0/49.0)/xfe2m4  + (12.0/49.0)/xmg2m4  \
                       - (24.0/49.0)/xfe2m13 - (24.0/49.0)/xmg2m13 \
                       - (30.0/49.0)/xfe2m2  - (30.0/49.0)/xmg2m2 )
#define D3GDS0S1DP  VS13S2
#define D3GDS0DT2   0.0
#define D3GDS0DTDP  0.0
#define D3GDS0DP2   0.0

#define D3GDS1S1S1  -R*t*(  (16.0/343.0)/SQUARE(xfe2m4)  \
                                                    - (16.0/343.0)/SQUARE(xmg2m4)  \
                                                    + (24.0/343.0)/SQUARE(xfe2m13) \
                                                    - (24.0/343.0)/SQUARE(xmg2m13) \
                                                    - (250.0/343.0)/SQUARE(xfe2m2) \
                                                    + (250.0/343.0)/SQUARE(xmg2m2) )
#define D3GDS1S1DT  R*(  (8.0/49.0)/xfe2m4   + (8.0/49.0)/xmg2m4   \
                       + (12.0/49.0)/xfe2m13 + (12.0/49.0)/xmg2m13 \
                       + (50.0/49.0)/xfe2m2  + (50.0/49.0)/xmg2m2 )
#define D3GDS1S1DP  2.0*VS2S2
#define D3GDS1DT2   0.0
#define D3GDS1DTDP  0.0
#define D3GDS1DP2   0.0

#define D3GDT3      0.0
#define D3GDT2DP    0.0
#define D3GDTDP2    0.0
#define D3GDP3      0.0

#define D3GDR0DT2   0.0
#define D3GDR0DTDP  0.0
#define D3GDR0DP2   0.0

/*
 *=============================================================================
 * Macros for automatic array initialization
 */
#define fillD2GDR2  d2gdr2[0][0]  = D2GDR0R0;

#define fillD2GDRDS d2gdrds[0][0] = D2GDR0S0;     d2gdrds[0][1] = D2GDR0S1;

#define fillD2GDRDT d2gdrdt[0]    = D2GDR0DT;

#define fillD2GDRDP d2gdrdp[0]    = D2GDR0DP;

#define fillD2GDS2  d2gds2[0][0]  = D2GDS0S0;     d2gds2[0][1] = D2GDS0S1; \
                                        d2gds2[1][0]  = d2gds2[0][1]; d2gds2[1][1] = D2GDS1S1;

#define fillD2GDSDT d2gdsdt[0]    = D2GDS0DT;     d2gdsdt[1] = D2GDS1DT;

#define fillD2GDSDP d2gdsdp[0]    = D2GDS0DP;     d2gdsdp[1] = D2GDS1DP;

#define fillD3GDR3 d3gdr3[0][0][0] = D3GDR0R0R0;

#define fillD3GDR2DS \
 d3gdr2ds[0][0][0] = D3GDR0R0S0; d3gdr2ds[0][0][1] = D3GDR0R0S1; \

#define fillD3GDR2DT d3gdr2dt[0][0] = D3GDR0R0DT;

#define fillD3GDR2DP d3gdr2dp[0][0] = D3GDR0R0DP;

#define fillD3GDRDS2 \
 d3gdrds2[0][0][0] = D3GDR0S0S0;        d3gdrds2[0][0][1] = D3GDR0S0S1; \
 d3gdrds2[0][1][0] = d3gdrds2[0][0][1]; d3gdrds2[0][1][1] = D3GDR0S1S1;

#define fillD3GDRDSDT \
 d3gdrdsdt[0][0] = D3GDR0S0DT; d3gdrdsdt[0][1] = D3GDR0S1DT;

#define fillD3GDRDSDP \
 d3gdrdsdp[0][0] = D3GDR0S0DP; d3gdrdsdp[0][1] = D3GDR0S1DP; \

#define fillD3GDS3 \
 d3gds3[0][0][0] = D3GDS0S0S0;      d3gds3[0][0][1] = D3GDS0S0S1; \
 d3gds3[0][1][0] = d3gds3[0][0][1]; d3gds3[0][1][1] = D3GDS0S1S1; \
 d3gds3[1][0][0] = d3gds3[0][0][1]; d3gds3[1][0][1] = d3gds3[0][1][1]; \
 d3gds3[1][1][0] = d3gds3[0][1][1]; d3gds3[1][1][1] = D3GDS1S1S1;

#define fillD3GDS2DT \
 d3gds2dt[0][0] = D3GDS0S0DT;     d3gds2dt[0][1] = D3GDS0S1DT; \
 d3gds2dt[1][0] = d3gds2dt[0][1]; d3gds2dt[1][1] = D3GDS1S1DT;

#define fillD3GDS2DP \
 d3gds2dp[0][0] = D3GDS0S0DP;     d3gds2dp[0][1] = D3GDS0S1DP; \
 d3gds2dp[1][0] = d3gds2dp[0][1]; d3gds2dp[1][1] = D3GDS1S1DP;

#define fillD3GDSDT2 \
 d3gdsdt2[0] = D3GDS0DT2; d3gdsdt2[1] = D3GDS1DT2;

#define fillD3GDSDTDP \
 d3gdsdtdp[0] = D3GDS0DTDP; d3gdsdtdp[1] = D3GDS1DTDP;

#define fillD3GDSDP2 \
 d3gdsdp2[0] = D3GDS0DP2; d3gdsdp2[1] = D3GDS1DP2;

#define fillD3GDRDT2 \
 d3gdrdt2[0] = D3GDR0DT2;

#define fillD3GDRDTDP \
 d3gdrdtdp[0] = D3GDR0DTDP;

#define fillD3GDRDP2 \
 d3gdrdp2[0] = D3GDR0DP2;

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
    int i, j, iter=0, signum;

    GET_SITE_FRACTIONS

    /* look-up or compute the current ordering state */
    if ( (t != tOld)       || (p != pOld) || (r[0] != rOld[0]) ) {
        double dgds[NS], sNew[NS];
        gsl_vector_view vvToDgds = gsl_vector_view_array(dgds, (size_t) NS);

        for (i=0; i<NS; i++) sOld[i] = 2.0;

        sNew[0] = 0.0;
        sNew[1] = 0.0;

        while ( ((ABS(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON) ||
             (ABS(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON)    ) &&
             (iter < MAX_ITER)) {
            double s[NS], deltaS[NS], lambda;
            gsl_vector_view vvToDeltaS = gsl_vector_view_array(deltaS, (size_t) NS);

            for (i=0; i<NS; i++) s[i] = sNew[i];

            xfe2m4  = (1.0+r[0]+6.0*s[0]/7.0+ 4.0*s[1]/7.0)/2.0;
            xfe2m13 = (1.0+r[0]-8.0*s[0]/7.0+ 4.0*s[1]/7.0)/2.0;
            xfe2m2  = (1.0+r[0]+6.0*s[0]/7.0-10.0*s[1]/7.0)/2.0;
            xmg2m4  = 1.0 - xfe2m4;
            xmg2m13 = 1.0 - xfe2m13;
            xmg2m2  = 1.0 - xfe2m2;

            if (xfe2m4  <= 0.0) xfe2m4  = DBL_EPSILON;
            if (xmg2m4  <= 0.0) xmg2m4  = DBL_EPSILON;
            if (xfe2m13 <= 0.0) xfe2m13 = DBL_EPSILON;
            if (xmg2m13 <= 0.0) xmg2m13 = DBL_EPSILON;
            if (xfe2m2  <= 0.0) xfe2m2  = DBL_EPSILON;
            if (xmg2m2  <= 0.0) xmg2m2  = DBL_EPSILON;

            if (xfe2m4  >= 1.0) xfe2m4  = 1.0 - DBL_EPSILON;
            if (xmg2m4  >= 1.0) xmg2m4  = 1.0 - DBL_EPSILON;
            if (xfe2m13 >= 1.0) xfe2m13 = 1.0 - DBL_EPSILON;
            if (xmg2m13 >= 1.0) xmg2m13 = 1.0 - DBL_EPSILON;
            if (xfe2m2  >= 1.0) xfe2m2  = 1.0 - DBL_EPSILON;
            if (xmg2m2  >= 1.0) xmg2m2  = 1.0 - DBL_EPSILON;

            dgds[0] = DGDS0;
            dgds[1] = DGDS1;

            d2gds2[0][0] = D2GDS0S0;
            d2gds2[0][1] = D2GDS0S1;
            d2gds2[1][0] = d2gds2[0][1];
            d2gds2[1][1] = D2GDS1S1;

            for (i=0; i<NS; i++) sOld[i] = s[i];

            /* original: gaussj(ptToD2gds2, NS, (double **) NULL, 0);
            for (i=0; i<NS; i++) {
            for(j=0; j<NS; j++) s[i] += - d2gds2[i][j]*dgds[j];
            deltaS[i] = s[i] - sOld[i];
            } */

            gsl_matrix_scale(ptToD2gds2, -1.0);
            melts_LU_decomp(ptToD2gds2, indexD2gds2, &signum);

            melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToDgds.vector, &vvToDeltaS.vector);
            for (i=0; i<NS; i++) s[i] += deltaS[i];


            /* Default steplength along search direction */
            lambda = 1.0;
            /* Test Projected value of X Fe2 M4 */
            if      (xfe2m4+lambda*(2.0*deltaS[1]/7.0+3.0*deltaS[0]/7.0) < 0.0)
                lambda = -xfe2m4/(2.0*deltaS[1]/7.0+3.0*deltaS[0]/7.0);
            else if (xfe2m4+lambda*(2.0*deltaS[1]/7.0+3.0*deltaS[0]/7.0) > 1.0)
                lambda = (1.0-xfe2m4)/(2.0*deltaS[1]/7.0+3.0*deltaS[0]/7.0);
            /* Test Projected value of X Fe2 M13 */
            if      (xfe2m13+lambda*(2.0*deltaS[1]/7.0-4.0*deltaS[0]/7.0) < 0.0)
                lambda = -xfe2m13/(2.0*deltaS[1]/7.0-4.0*deltaS[0]/7.0);
            else if (xfe2m13+lambda*(2.0*deltaS[1]/7.0-4.0*deltaS[0]/7.0) > 1.0)
                lambda = (1.0-xfe2m13)/(2.0*deltaS[1]/7.0-4.0*deltaS[0]/7.0);
            /* Test Projected value of X Fe2 M2 */
            if      (xfe2m2+lambda*(-5.0*deltaS[1]/7.0+3.0*deltaS[0]/7.0) < 0.0)
                lambda = -xfe2m2/(-5.0*deltaS[1]/7.0+3.0*deltaS[0]/7.0);
            else if (xfe2m2+lambda*(-5.0*deltaS[1]/7.0+3.0*deltaS[0]/7.0) > 1.0)
                lambda = (1.0-xfe2m2)/(-5.0*deltaS[1]/7.0+3.0*deltaS[0]/7.0);
            /* Modify steplength if required to maintain feasibility */
            if (lambda < 1.0) for (i=0; i<NS; i++) s[i] = sOld[i] + lambda*deltaS[i];

            for (i=0; i<NS; i++) sNew[i] = s[i];
            iter++;
        }
        tOld = t;
        pOld = p;
        for (i=0; i<NR; i++) rOld[i] = r[i];

#ifdef DEBUG
        for (i=0; i<NS; i++) {
            if (dgds[i] > sqrt(DBL_EPSILON) && ABS(sOld[i]) > DBL_EPSILON) {
                printf("ERROR in CUMMINGTONITE.C (func ORDER). Failed to converge!\n");
                printf("  r      = %13.6g\n", r[0]);
                printf("  s13    = %13.6g, s2    = %13.6g\n", sOld[0], sOld[1]);
                printf("  dgds13 = %13.6g, dgds2 = %13.6g\n", dgds[0], dgds[1]);
                printf("  X Fe2+ M4 : %13.6g  X Mg   M4 : %13.6g\n", xfe2m4 , xmg2m4 );
                printf("  X Fe2+ M13: %13.6g  X Mg   M13: %13.6g\n", xfe2m13, xmg2m13);
                printf("  X Fe2+ M2 : %13.6g  X Mg   M2 : %13.6g\n", xfe2m2 , xmg2m2 );
                break;
            }
        }
#endif /* PRINT_NONCONVERGENCE_IN_ORDER */

        setTOld(tOld);
        setPOld(pOld);
        /* arrays (rOld, sOld, d2gds2, indexD2gds2) should be preserved automatically */

        SET_SITE_FRACTIONS
    }

    if (mask & FIRST  ) {   /* return s        */
        for (i=0; i<NS; i++) s[i] = sOld[i];
    }
    if (mask & SECOND ) {   /* compute ds/dr:  */
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
    if (mask & THIRD  ) {   /* compute ds/dt:  */
        double d2gdsdt[NS];
        gsl_vector_view vvToDt = gsl_vector_view_array(dt, (size_t) NS),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NS);

        fillD2GDSDT

        /* original: dt[i] += - d2gds2[i][j]*d2gdsdt[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdt.vector, &vvToDt.vector);
    }
    if (mask & FOURTH ) {   /* compute ds/dp:  */
        double *s = sOld;
        double d2gdsdp[NS];
        gsl_vector_view vvToDp = gsl_vector_view_array(dp, (size_t) NS),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NS);

        fillD2GDSDP

        /* original: dp[i] += - d2gds2[i][j]*d2gdsdp[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdp.vector, &vvToDp.vector);
    }
    if (mask & FIFTH  ) {   /* compute d2s/dr2 */
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
    if (mask & SIXTH  ) {   /* compute d2s/drt */
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
    if (mask & EIGHTH ) {   /* compute d2s/dt2 */
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
    if (mask & NINTH  ) {   /* compute d2s/dtp */
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
    if (mask & TENTH  ) {   /* compute d2s/dp2 */
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
testCum(int mask, double t, double p,
    int na,          /* Expected number of endmember components                 */
    int nr,          /* Expected number of independent compositional variables  */
    char **names,    /* array of strings of names of endmember components       */
    char **formulas, /* array of strings of formulas of endmember components    */
    double *r,       /* array of indepependent compos variables, check bounds   */
    double *m)       /* array of moles of endmember components, check bounds    */
{
    const char *phase = "cummingtonite.c";
    const char *NAMES[NA]    = { "cummingtonite", "grunerite" };
    const char *FORMULAS[NA] = { "Mg7Si8O22(OH)2", "Fe7Si8O22(OH)2" };
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
        result = result && (r[0] >= -1.0) && (r[0] <= 1.0);
    }
    /* Check bounds on moles of endmember components */
    if (mask & SIXTH) {
        for (i=0, sum=0.0; i<NA; i++) sum += m[i];
        result = result && (sum >= 0.0);
        if (sum > 0.0) {
            result = result && (m[0]/sum >= 0.0) && (m[0]/sum <= 1.0);
            result = result && (m[1]/sum >= 0.0) && (m[1]/sum <= 1.0);
        }
    }

    return result;
}

void
conCum(int inpMask, int outMask, double t, double p,
    double *e,      /* comp of spinel in moles of elements                      */
    double *m,      /* comp of spinel in moles of endmember components          */
    double *r,      /* comp of spinel in terms of the independent comp var      */
    double *x,      /* comp of spinel in mole fractions of endmember comp       */
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
            endmember spinel components.
    (2) calculates from a vector of moles of endmember components, one or
            all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
    (3) calculates from a vector of independent compositional variables
            mole fractions of endmember components and/or the Jacobian matrix
            dx[]/dr[]

    In this routine it is assumed that the elements are in the order of atomic
    numbers and that the order of cummingtonite components has been verified as:
            m[0] = cummingtonite (Mg7Si8O22(OH)2),
            m[1] = grunerite     (Fe7Si8O22(OH)2)

    ----------------------------------------------------------------------------*/

    int i, j, k;

    if (inpMask == FIRST && outMask == SECOND) {
        /* Converts a vector of moles of elements into a vector of moles of
       end-member components.                                                 */
        static const int Mg = 12;
        static const int Fe = 26;

        m[0] = e[Mg]/7.0;
        m[1] = e[Fe]/7.0;

    } else if (inpMask == SECOND) {
        double sum;

        if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
            printf("Illegal call to conCum with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        for (i=0, sum=0.0; i<NA; i++) sum += m[i];

        if (outMask & THIRD) {
            /* Converts a vector of moles of end-member components (m) into a vector
         of independent compositional variables (r) required as input for the
         remaining public functions.                                          */
            r[0] = (sum != 0.0) ? 2.0*m[1]/sum - 1.0 : 0.0;
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
                    dm[0][j] = (j == 1) ? 2.0*(1.0-m[1]/sum)/sum : -2.0*m[1]/SQUARE(sum);
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
                        d2m[0][j][k]  = 4.0*m[1]/CUBE(sum);
                        d2m[0][j][k] -= (j == 1) ? 2.0/SQUARE(sum) : 0.0;
                        d2m[0][j][k] -= (k == 1) ? 2.0/SQUARE(sum) : 0.0;
                    }
                }
            }

        }

        if (outMask & EIGHTH) {
            /* Calculates the matrix d3r[i]/dm[j]dm[k]dm[l] using m[] as input         */
            int l;

            if (sum == 0.0) {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++)  {
                        for (k=0; k<NA; k++)   {
	      for (l=0; l<NA; l++) d3m[i][j][k][l] = 0.0;
	    }
                    }
                }
            } else {
                for (j=0; j<NA; j++) {
                    for (k=0; k<NA; k++) {
	    for (l=0; l<NA; l++) {
                            d3m[0][j][k][l]  = -12.0*m[1]/SQUARE(SQUARE(sum));
                            d3m[0][j][k][l] += (j == 1) ? 4.0/CUBE(sum) : 0.0;
                            d3m[0][j][k][l] += (k == 1) ? 4.0/CUBE(sum) : 0.0;
                            d3m[0][j][k][l] += (l == 1) ? 4.0/CUBE(sum) : 0.0;
                        }
	  }
                }
            }

        }
    } else if (inpMask == THIRD) {

        if (outMask & ~(FOURTH | SEVENTH))
            printf("Illegal call to conCum with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        if (outMask & FOURTH) {
            /* Converts a vector of independent compositional variables (r) into a
         vector of mole fractions of endmember components (x).                */
            x[0] = (1.0-r[0])/2.0;
            x[1] = (1.0+r[0])/2.0;
        }

        if (outMask & SEVENTH) {
            /* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
            for (i=0; i<NA; i++) for (j=0; j<NR; j++) dr[i][j] = 0.0;
            dr[0][0] = -1.0/2.0;
            dr[1][0] =  1.0/2.0;
        }

    } else  {
        printf("Illegal call to conCum with inpMask = %o and outMask = %o\n",
            inpMask, outMask);
    }

}

void
dispCum(int mask, double t, double p, double *x,
    char **formula            /* Mineral formula for interface display MASK: 1 */
    )
{
    double *r = x;
    static char masterString[] = {
/*             1111111111222222222233333333334444444444555555555566666666667
     01234567890123456789012345678901234567890123456789012345678901234567890 */
        "(Fe''_.__Mg_.__)7Si8O22(OH)2" };

    if (mask & FIRST) {
        char *string = strcpy((char *) malloc((size_t) (strlen(masterString)+1)*sizeof(char)), masterString);
        double totFe2, totMg;
        char n[5];
        int i;

        totFe2 = (1.0+r[0])/2.0;
        totMg  = (1.0-r[0])/2.0;

        (void) snprintf(n, 5, "%4.2f", totFe2); for (i=0; i<4; i++) string[ 5+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totMg);  for (i=0; i<4; i++) string[11+i] = n[i];

        *formula = string;
    }
}

void
actCum(int mask, double t, double p, double *x,
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
     fr[i][0] = FR0(i); /* r */
    }

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    g       = G;
    dgdr[0] = DGDR0;

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

        for(i=0; i<NA; i++) {
       for (mu[i]=g, j=0; j<NR; j++) mu[i] += fr[i][j]*dgdr[j];
        }

    }

    if (mask & THIRD) {
        double d2gdr2[NR][NR], d2gdrds[NR][NS], d2gds2[NS][NS], dsdr[NS][NR],
            dfrdr[NA][NR], gs[NA][NS], dgsds[NA][NS], sum;
        int k, l;

        fillD2GDR2
        fillD2GDRDS
        fillD2GDS2

        for(i=0; i<NA; i++) {
       gs[i][0]    = GS0(i);     /* s13 */
       gs[i][1]    = GS1(i);     /* s2  */
       dfrdr[i][0] = DFR0DR0(i); /* r   */
       dgsds[i][0] = DGS0DS0(i); /* s13 */
       dgsds[i][1] = DGS1DS1(i); /* s2  */
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
       0.05,  /* exclusion criteria on the mole fraction of cummingtonite     */
       0.05   /* exclusion criteria on the mole fraction of grunerite         */
        };
        double x[NA];

        x[0] = (1.0-r[0])/2.0;
        x[1] = (1.0+r[0])/2.0;

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
gmixCum(int mask, double t, double p, double *x,
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
hmixCum(int mask, double t, double p, double *x,
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
smixCum(int mask, double t, double p, double *x,
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
cpmixCum(int mask, double t, double p, double *x,
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
vmixCum(int mask, double t, double p, double *x,
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

/* end of file CUMMINGTONITE.C */
