const char *albite_ver(void) { return "$Id: albite.c,v 1.4 2007/11/22 04:08:13 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: albite.c,v $
MELTS Source Code: RCS Revision 1.4  2007/11/22 04:08:13  ghiorso
MELTS Source Code: RCS Corrected infinite loop error in order() in albite.c
MELTS Source Code: RCS Removed arbitrary volume corrections in sol_struct_data.h
MELTS Source Code: RCS Turned on non-quadrilateral cpx endmembers for regression.
MELTS Source Code: RCS Added MgSiO3 species to liquid model.
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
MELTS Source Code: RCS Revision 1.1.1.1  2006/08/15 16:57:35  ghiorso
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
 * Revision 3.6  1997/06/21  22:50:11  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.5  1997/05/03  20:23:45  ghiorso
 * *** empty log message ***
 *
 * Revision 3.4  1997/03/27  17:03:53  ghiorso
 * *** empty log message ***
 *
 * Revision 3.3  1996/09/24  20:31:03  ghiorso
 * *** empty log message ***
 *
 * Revision 3.2  1995/11/23  22:37:42  ghiorso
 * Final implementation of subsolidus fO2 buffering.
 *
 * Revision 3.2  1995/11/23  22:37:42  ghiorso
 * Final implementation of subsolidus fO2 buffering.
 *
 * Revision 3.1  1995/08/18  17:12:33  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**  Routine to compute thermodynamic properties of albite ordering
**  relative to low albite standard state function of Salje (1985)
**  uses monalbite std state (file: ALBITE.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  February 7, 1995
**              (1) Original Version
**--
*/

#include "melts_gsl.h"
#include "silmin.h"

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))

#define MAX_ITER 200    /* Maximum number of iterations allowed in order */
#define NS       2      /* Two ordering parameters                       */

static const double a0   =     5.479;
static const double b	 =  6854.0;
static const double aod0 =    41.620;
static const double bod  = -9301.0;
static const double cod  = 43600.0;
static const double d0   =    -2.171;
static const double d1   =    -3.043;
static const double d2   =    -0.001569;
static const double d3   =     0.000002109;
static const double tc   =  1251.0;
static const double tod  =   824.1;

/*************************************/
/* Statics for Ordering Calculations */
/*************************************/

static MTHREAD_ONCE_T initThreadOBlock = MTHREAD_ONCE_INIT;

static MTHREAD_KEY_T tOldKey;
static MTHREAD_KEY_T pOldKey;
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

/* q = s[0], qod = s[1] */

#define S   -(0.5*a0*s[0]*s[0] + 0.5*aod0*s[1]*s[1] \
                            + (d1+2.0*d2*t+3.0*d3*t*t)*s[0]*s[1])
#define H   - 0.5*a0*tc*s[0]*s[0] + 0.25 *b*s[0]*s[0]*s[0]*s[0] \
                        - 0.5*aod0*tod*s[1]*s[1] + 0.25*bod*s[1]*s[1]*s[1]*s[1] \
                        + (cod/6.0)*s[1]*s[1]*s[1]*s[1]*s[1]*s[1] \
                        + (d0 - d2*t*t - 2.0*d3*t*t*t)*s[0]*s[1]
#define V   (H)/335282.925
#define G   (H) - t*(S) + (p-1.0)*(V)

/*----------------------------------------------------------------------------*/

#define DGDS0 (-a0*tc*s[0] + b*s[0]*s[0]*s[0] \
               + (d0-d2*t*t-2.0*d3*t*t*t)*s[1])*(1.0+(p-1.0)/335282.925) \
               + t*(a0*s[0] + (d1+2.0*d2*t+3.0*d3*t*t)*s[1])
#define DGDS1 (-aod0*tod*s[1] + bod*s[1]*s[1]*s[1] \
               + cod*s[1]*s[1]*s[1]*s[1]*s[1] \
               + (d0-d2*t*t-2.0*d3*t*t*t)*s[0])*(1.0+(p-1.0)/335282.925) \
               + t*(aod0*s[1] + (d1+2.0*d2*t+3.0*d3*t*t)*s[0])
#define DGDT  -(S) + (-2.0*d2*t-6.0*d3*t*t)*s[0]*s[1]*(p-1.0)/335282.925
#define DGDP  (V)

/*----------------------------------------------------------------------------*/

#define D2GDS0S0 (-a0*tc + 3.0*b*s[0]*s[0])*(1.0+(p-1.0)/335282.925) + t*a0
#define D2GDS0S1 (d0-d2*t*t-2.0*d3*t*t*t)*(1.0+(p-1.0)/335282.925) \
                 + t*(d1+2.0*d2*t+3.0*d3*t*t)
#define D2GDS0DT -2.0*(d2*t+3.0*d3*t*t)*s[1]*(p-1.0)/335282.925 \
                 + a0*s[0] + (d1+2.0*d2*t+3.0*d3*t*t)*s[1]
#define D2GDS0DP (-a0*tc*s[0] + b*s[0]*s[0]*s[0] \
                                    + (d0-d2*t*t-2.0*d3*t*t*t)*s[1])/335282.925

#define D2GDS1S1 (-aod0*tod + 3.0*bod*s[1]*s[1] + 5.0*cod*s[1]*s[1]*s[1]*s[1]) \
                                    *(1.0+(p-1.0)/335282.925) + t*aod0
#define D2GDS1DT -2.0*(d2*t+3.0*d3*t*t)*s[0]*(p-1.0)/335282.925 \
                 + aod0*s[1] + (d1+2.0*d2*t+3.0*d3*t*t)*s[0]
#define D2GDS1DP (-aod0*tod*s[1] + bod*s[1]*s[1]*s[1] \
                                    + cod*s[1]*s[1]*s[1]*s[1]*s[1] \
                                    + (d0-d2*t*t-2.0*d3*t*t*t)*s[0])/335282.925

#define D2GDT2   -2.0*(d2+6.0*d3*t)*s[0]*s[1]*(p-1.0)/335282.925 \
                 + (2.0*d2+6.0*d3*t)*s[0]*s[1]
#define D2GDTDP  -2.0*(d2*t+3.0*d3*t*t)*s[0]*s[1]/335282.925
#define D2GDP2   0.0

/*----------------------------------------------------------------------------*/

#define D3GDS0S0S0 6.0*b*s[0]*(1.0+(p-1.0)/335282.925)
#define D3GDS0S0S1 0.0
#define D3GDS0S0DT a0
#define D3GDS0S0DP (-a0*tc + 3.0*b*s[0]*s[0])/335282.925
#define D3GDS0S1S1 0.0
#define D3GDS0S1DT (-2.0*d2*t-6.0*d3*t*t)*(1.0+(p-1.0)/335282.925) \
                   + d1 + 2.0*d2*t + 3.0*d3*t*t + t*(2.0*d2+6.0*d3*t)
#define D3GDS0S1DP (d0-d2*t*t-2.0*d3*t*t*t)/335282.925
#define D3GDS0DT2  -2.0*(d2+6.0*d3*t)*s[1]*(p-1.0)/335282.925 \
                   + (2.0*d2+6.0*d3*t)*s[1]
#define D3GDS0DTDP -2.0*(d2*t+3.0*d3*t*t)*s[1]/335282.925
#define D3GDS0DP2  0.0

#define D3GDS1S1S1 (6.0*bod*s[1] + 20.0*cod*s[1]*s[1]*s[1]) \
                                        *(1.0+(p-1.0)/335282.925)
#define D3GDS1S1DT aod0
#define D3GDS1S1DP (-aod0*tod + 3.0*bod*s[1]*s[1] \
                                        + 5.0*cod*s[1]*s[1]*s[1]*s[1])/335282.925
#define D3GDS1DT2  -2.0*(d2+6.0*d3*t)*s[0]*(p-1.0)/335282.925 \
                   + (2.0*d2+6.0*d3*t)*s[0]
#define D3GDS1DTDP -2.0*(d2*t+3.0*d3*t*t)*s[0]/335282.925
#define D3GDS1DP2  0.0

#define D3GDT3     -12.0*d3*s[0]*s[1]*(p-1.0)/335282.925 + 6.0*d3*s[0]*s[1]
#define D3GDT2DP   -2.0*(d2+6.0*d3*t)*s[0]*s[1]/335282.925
#define D3GDTDP2   0.0
#define D3GDP3     0.0

#define fillD2GDS2 \
 d2gds2[0][0] = D2GDS0S0;     d2gds2[0][1] = D2GDS0S1; \
 d2gds2[1][0] = d2gds2[0][1]; d2gds2[1][1] = D2GDS1S1;

#define fillD2GDSDT \
 d2gdsdt[0] = D2GDS0DT;  d2gdsdt[1] = D2GDS1DT;

#define fillD2GDSDP \
 d2gdsdp[0] = D2GDS0DP;  d2gdsdp[1] = D2GDS1DP;

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

static void
order(double t, double p,
            double s[NS],           /* s[NS]       */
            double dt[NS],          /* ds[NS]/dt   */
            double dp[NS],          /* ds[NS]/dp   */
            double dt2[NS],         /* d2s[NS]/dt2 */
            double dtp[NS],         /* d2s[NS]/dtp */
            double dp2[NS]          /* d2s[NS]/dp2 */
            )
{
    double tOld         = getTOld();
    double pOld         = getPOld();
    double *sOld        = getSOld();
    double **d2gds2     = getD2gds2();
    gsl_matrix      *ptToD2gds2  = getPtToD2gds2();
    gsl_permutation *indexD2gds2 = getIndexD2gds2();

    int i, j, k, l, iter=0, signum;
    double d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS][NS][NS], d3gds2dt[NS][NS],
        d3gdsdt2[NS], temp[NS], d3gds2dp[NS][NS], d3gdsdtdp[NS], d3gdsdp2[NS];

    gsl_vector_view vvToDt = gsl_vector_view_array(dt, (size_t) NS),
        vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NS),
        vvToDp = gsl_vector_view_array(dp, (size_t) NS),
        vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NS),
        vvToTemp = gsl_vector_view_array(temp, (size_t) NS);

    /* look-up or compute the current ordering state */
    if ( (t != tOld) || (p != pOld) ) {
        double dgds[NS], sNew[NS];
        gsl_vector_view vvToDgds = gsl_vector_view_array(dgds, (size_t) NS);

        for (i=0; i<NS; i++) sOld[i] = 2.0;

        sNew[0] = 0.6;
        sNew[1] = 0.9;

        while (   ((ABS(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON) || (ABS(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON))
           && (iter < MAX_ITER)) {
            double s[NS], deltaS[NS];
            gsl_vector_view vvToDeltaS = gsl_vector_view_array(deltaS, (size_t) NS);

            for (i=0; i<NS; i++) s[i] = sNew[i];

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

            for (i=0; i<NS; i++) sNew[i] = s[i];
            iter++;
        }
        tOld = t;
        pOld = p;
        setTOld(tOld);
        setPOld(pOld);
    }

    /* s */
    for (i=0; i<NS; i++) s[i] = sOld[i];

    /* dsdt */
    fillD2GDSDT

    /* original: dt[i] += - d2gds2[i][j]*d2gdsdt[j]; */
    melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdt.vector, &vvToDt.vector);

    /* dsdp */
    fillD2GDSDP

    /* original: dp[i] += - d2gds2[i][j]*d2gdsdp[j]; */
    melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdp.vector, &vvToDp.vector);

    /* d2sdt2 */
    fillD2GDSDT
    fillD3GDS3
    fillD3GDS2DT
    fillD3GDSDT2

    for (j=0; j<NS; j++) {
        temp[j] = d3gdsdt2[j];
        for (k=0; k<NS; k++) {
            temp[j] +=  2.0*d3gds2dt[j][k]*dt[k];
            for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dt[k]*dt[l];
        }
    }
    /* original: dt2[i] += - d2gds2[i][j]*temp[j]; */
    melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
    for (j=0; j<NS; j++) dt2[j] = temp[j];

    /* d2sdtdp */
    fillD2GDSDT
    fillD2GDSDP
    fillD3GDS3
    fillD3GDS2DT
    fillD3GDS2DP
    fillD3GDSDTDP

    for (j=0; j<NS; j++) {
        temp[j] = d3gdsdtdp[j];
        for (k=0; k<NS; k++) {
            temp[j] += d3gds2dt[j][k]*dp[k] + d3gds2dp[j][k]*dt[k];
            for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dt[k]*dp[l];
        }
    }
    /* original: dtp[i] += - d2gds2[i][j]*temp[j]; */
    melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
    for (j=0; j<NS; j++) dtp[j] = temp[j];

    /* d2sdp2 */
    fillD2GDSDP
    fillD3GDS3
    fillD3GDS2DP
    fillD3GDSDP2

    for (j=0; j<NS; j++) {
        temp[j] = d3gdsdp2[j];
        for (k=0; k<NS; k++) {
            temp[j] +=  2.0*d3gds2dp[j][k]*dp[k];
            for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dp[k]*dp[l];
        }
    }
    /* original: dp2[i] += - d2gds2[i][j]*temp[j]; */
    melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
    for (j=0; j<NS; j++) dp2[j] = temp[j];

}

void albite(double p, double t, double *gDis, double *hDis, double *sDis,
   double *cpDis, double *dcpdt, double *vDis, double *dvdt, double *dvdp,
   double *d2vdt2, double *d2vdtdp, double *d2vdp2)
{
   double s[NS], dsdt[NS], dsdp[NS], d2sdt2[NS], d2sdtdp[NS], d2sdp2[NS];
   int skip;
   *gDis    = 0.0;
   *hDis    = 0.0;
   *sDis    = 0.0;
   *cpDis   = 0.0;
   *dcpdt   = 0.0;
   *vDis    = 0.0;
   *dvdt    = 0.0;
   *dvdp    = 0.0;
   *d2vdt2  = 0.0;
   *d2vdtdp = 0.0;
   *d2vdp2  = 0.0;

   order(t, p, s, dsdt, dsdp, d2sdt2, d2sdtdp, d2sdp2);
   skip = (fabs(s[0]) < sqrt(DBL_EPSILON)) && (fabs(s[1]) < sqrt(DBL_EPSILON));

/* if (t > 1290.0) { This was Berman's test - change on 8/19/09 */
   if (skip) {
            return;
   } else {
            int i, j, k;
            double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], d2gdt2, d2gdtdp, d2gdp2;
            double d3gds3[NS][NS][NS], d3gds2dt[NS][NS], d3gds2dp[NS][NS],
                d3gdsdt2[NS], d3gdsdtdp[NS], d3gdsdp2[NS], d3gdt3, d3gdt2dp, d3gdtdp2,
                d3gdp3;
            double temp;

/* ---------- */
            *gDis = (G);
/* ---------- */
            *hDis = *gDis - t*(DGDT);
/* ---------- */
            *sDis = -(DGDT);
/* ---------- */
            fillD2GDS2
            fillD2GDSDT
            d2gdt2  = D2GDT2;

            *cpDis = d2gdt2;
            for (i=0; i<NS; i++) {
                *cpDis += 2.0*d2gdsdt[i]*dsdt[i];
                for (j=0; j<NS; j++) *cpDis += d2gds2[i][j]*dsdt[i]*dsdt[j];
            }
            temp = *cpDis;
            *cpDis *= -t;
/* ---------- */
            fillD3GDS3
            fillD3GDS2DT
            fillD3GDSDT2
            d3gdt3 = D3GDT3;

            *dcpdt = d3gdt3;
            for (i=0; i<NS; i++) {
                *dcpdt += 3.0*d3gdsdt2[i]*dsdt[i] + 3.0*d2gdsdt[i]*d2sdt2[i];
                for (j=0; j<NS; j++) {
                    *dcpdt += 3.0*d2gds2[i][j]*dsdt[i]*d2sdt2[j]
                 + 3.0*d3gds2dt[i][j]*dsdt[i]*dsdt[j];
                    for (k=0; k<NS; k++) *dcpdt +=
                        d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdt[k];
                }
            }
            *dcpdt = -t*(*dcpdt) - temp;
/* ---------- */
            *vDis = (V);
/* ---------- */
            fillD2GDSDP
            d2gdtdp = D2GDTDP;

            *dvdt = d2gdtdp;
            for (i=0; i<NS; i++) {
                *dvdt += d2gdsdt[i]*dsdp[i] + d2gdsdp[i]*dsdt[i];
                for (j=0; j<NS; j++) *dvdt += d2gds2[i][j]*dsdt[i]*dsdp[j];
            }
/* ---------- */
            d2gdp2 = D2GDP2;

            *dvdp = d2gdp2;
            for (i=0; i<NS; i++) {
                *dvdp += 2.0*d2gdsdp[i]*dsdp[i];
                for (j=0; j<NS; j++) *dvdp += d2gds2[i][j]*dsdp[i]*dsdp[j];
            }
/* ---------- */
            fillD3GDS2DP
            fillD3GDSDTDP
            d3gdt2dp = D3GDT2DP;

            *d2vdt2 = d3gdt2dp;
            for (i=0; i<NS; i++) {
                *d2vdt2 += d3gdsdt2[i]*dsdp[i] + 2.0*d2gdsdt[i]*d2sdtdp[i]
                                + d2gdsdp[i]*d2sdt2[i] + 2.0*d3gdsdtdp[i]*dsdt[i];
                for (j=0; j<NS; j++) {
                    *d2vdt2 += 2.0*d3gds2dt[i][j]*dsdt[i]*dsdp[j]
                                + d2gds2[i][j]*d2sdt2[i]*dsdp[j]
                                + 2.0*d2gds2[i][j]*dsdt[i]*d2sdtdp[j]
                                + d3gds2dp[i][j]*dsdt[i]*dsdt[j];
                    for (k=0; k<NS; k++) *d2vdt2
                        += d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdp[k];
                }
            }
/* ---------- */
            fillD3GDSDP2
            d3gdtdp2 = D3GDTDP2;

            *d2vdtdp = d3gdtdp2;
            for (i=0; i<NS; i++) {
                *d2vdtdp += 2.0*d3gdsdtdp[i]*dsdp[i] + d2gdsdt[i]*d2sdp2[i]
                 + 2.0*d2gdsdp[i]*d2sdtdp[i] + d3gdsdp2[i]*dsdt[i];
                for (j=0; j<NS; j++) {
                    *d2vdtdp += 2.0*d3gds2dp[i][j]*dsdt[i]*dsdp[j]
                   + d2gds2[i][j]*dsdt[i]*d2sdp2[j]
                   + 2.0*d2gds2[i][j]*d2sdtdp[i]*dsdp[j]
                   + d3gds2dt[i][j]*dsdp[i]*dsdp[j];
                    for (k=0; k<NS; k++) *d2vdtdp
                        += d3gds3[i][j][k]*dsdt[i]*dsdp[j]*dsdp[k];
                }
            }
/* ---------- */
            d3gdp3 = D3GDP3;

            *d2vdp2 = d3gdp3;
            for (i=0; i<NS; i++) {
                *d2vdp2 += 3.0*d3gdsdp2[i]*dsdp[i] + 3.0*d2gdsdp[i]*d2sdp2[i];
                for (j=0; j<NS; j++) {
                    *d2vdp2 += 3.0*d2gds2[i][j]*dsdp[i]*d2sdp2[j]
                                    + 3.0*d3gds2dp[i][j]*dsdp[i]*dsdp[j];
                    for (k=0; k<NS; k++) *d2vdp2
                        += d3gds3[i][j][k]*dsdp[i]*dsdp[j]*dsdp[k];
                }
            }
   }
}

/* end of file ALBITE.C */
