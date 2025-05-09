const char *check_coexisting_liquids_ver(void) { return "$Id: check_coexisting_liquids.c,v 1.2 2006/08/17 16:47:18 ghiorso Exp $"; }
/*
 MELTS Source Code: RCS $Log: check_coexisting_liquids.c,v $
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
 * Revision 1.5  1997/06/21  22:50:08  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 1.4  1997/05/03  20:23:43  ghiorso
 * *** empty log message ***
 *
 * Revision 1.3  1997/03/27  17:03:51  ghiorso
 * *** empty log message ***
 *
 * Revision 1.2  1996/09/24  20:33:52  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 1.1  1995/12/17  00:19:08  ghiorso
 * Initial revision
 *
 */

/*
 **++
 **  FACILITY:  Silicate Melts Regression/Crystallization Package
 **
 **  MODULE DESCRIPTION:
 **
 **      Routines to determine if immiscible solid phases are present in the
 **      system.
 **      (file: CHECK_COEXISTING_LIQUIDS.C)
 **
 **  MODIFICATION HISTORY:
 **
 **      V1.0-1  Mark S. Ghiorso  December 16, 1995  Original Version
 **              Canabalized from check_coexisting_solids.c.
 **--
 */

#include "melts_gsl.h"

#include "silmin.h"
#include "nash.h"
#include "lawson_hanson.h"

#define SQUARE(x) ((x)*(x))
#define REALLOC(x, y) (((x) == NULL) ? malloc(y) : realloc((x), (y)))

#ifdef DEBUG
#undef DEBUG
#endif

/*****************************************************************************
 * Globally static variables
 *****************************************************************************/

static double  *bRef;     /* reference initial composition                    */
static double  *dg;       /* temporary storage for gradient                   */
static double  *dgRef;    /* reference gradient for input composition         */
static double **dxdr;     /* temporary storage for d(mole frac)/d(indep)      */
static double   gRef;     /* reference function value for input composition   */
static double **gHess;    /* Hessian of g at input composition for spinode    */
static int      hasNull;  /* TRUE if some endmembers have zero concen         */
static double  *hVec;     /* pivot elements for H12 decomp of eq constr       */
static int      na;       /* number of endmembers in liquid phase             */
static int      ne;       /* number of equality constraints                   */
static int      nr;       /* number of independent variables in liquid phase  */
static int      nz;       /* number of non-zero endmember concentrations      */
static int     *nullComp; /* vector of TRUE/FALSE flags for zero concen       */
static double   p;        /* pressure (bars)                                  */
static double   t;        /* temperature (K)                                  */
static double  *tVec;     /* temporary soln for fminfn and fmingr routines    */
static double  *yVec;     /* range space soln if prob has eq constr           */

/*****************************************************************************
 * Private function definitions:
 *****************************************************************************/

static int noFminfnCalls;
static int noFmingrCalls;

static double fminfn(int n, double *bVec, int *notcomp);
static void   fmingr(int n, double *bVec, double *g);

/*****************************************************************************
 * Global function declarations:
 *****************************************************************************/

#define SUCCESS TRUE
#define FAILURE FALSE

static int checkpd(gsl_matrix *a, int n) {
    gsl_matrix_view p = gsl_matrix_submatrix(a, (size_t) 0, (size_t) 0, (size_t) n, (size_t) n);
    int status;

    (void) gsl_set_error_handler_off();
    status = gsl_linalg_cholesky_decomp1(&p.matrix);
    (void) gsl_set_error_handler(NULL);

    return (status == GSL_SUCCESS);
}

int checkForCoexistingLiquids(  /* returns a MODE flag for success or failure */
                              void)
{
    static gsl_matrix *gHessMat;
    static double **bVec, *Fmin, **mVec, *rTr;
    double reltest;
    int i, j, k, l, mode, np, nl, nHess, result;

#ifdef DEBUG
    printf("Call to checkForCoexistingLiquids\n");
#endif

    if (silminState->liquidMass == 0.0) return FAILURE;

    if (dgRef == NULL) {
        na        = nlc;
        nr        = nlc - 1;
        bRef     = (double *)  malloc((size_t) nr*sizeof(double));
        bVec     = (double **) malloc((size_t) na*sizeof(double *));
        for (i=0; i<na; i++) bVec[i]  = (double *)  malloc((size_t) nr*sizeof(double));
        dg       = (double *)  malloc((size_t) nr*sizeof(double));
        dgRef    = (double *)  malloc((size_t) nr*sizeof(double));
        gHess    = (double **) malloc((size_t) nr*sizeof(double *));
        for (i=0; i<nr; i++) gHess[i] = (double *)  malloc((size_t) nr*sizeof(double));
        gHessMat = gsl_matrix_alloc((size_t) nr, (size_t) nr);
        dxdr     = (double **) malloc((size_t) na*sizeof(double *));
        for (i=0; i<na; i++) dxdr[i]  = (double *)  malloc((size_t) nr*sizeof(double));
        Fmin     = (double *)  malloc((size_t) na*sizeof(double));
        hVec     = (double *)  malloc((size_t) nr*sizeof(double));
        mVec     = (double **) malloc((size_t) na*sizeof(double *));
        for (i=0; i<na; i++) mVec[i]  = (double *)  malloc((size_t) na*sizeof(double));
        nullComp = (int *)	malloc((size_t) na*sizeof(int));
        rTr      = (double *)  malloc((size_t) na*sizeof(double));
        tVec     = (double *)  malloc((size_t) nr*sizeof (double));
        yVec     = (double *)  malloc((size_t) nr*sizeof(double));
    }

    t       = silminState->T; if (t <= 0.0) return FAILURE;
    p       = silminState->P; if (p <  0.0) return FAILURE;
    reltest = sqrt(DBL_EPSILON);
    result  = FAILURE;

    for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
        for (i=0; i<na; i++) mVec[0][i] = (silminState->liquidComp)[nl][i];

        conLiq(SECOND, THIRD, t, p, NULL, mVec[0], bRef, NULL, NULL, NULL, NULL);
        gmixLiq(FIRST | SECOND | THIRD, t, p, bRef, &gRef, dgRef, gHess);

        /* Determine if the Hessian is positive definite */
        for (i=0, k=0; i<nr; i++) if (bRef[i] != 0.0) {
            for (j=0, l=0; j<nr; j++) if (bRef[j] != 0.0) {
                gsl_matrix_set(gHessMat, k, l, gHess[i][j]);
                gsl_matrix_set(gHessMat, l, k, gHess[i][j]);
                l++;
            }
            k++;
        }
        nHess = k;

        if (!checkpd(gHessMat, nHess)) {
#ifdef DEBUG
            printf("Check for coexisting liquids called inside spinodal!\n");
#endif
            /* return result; */
        }

        /* Determine if any endmembers have zero concentration */
        for (i=0, hasNull=FALSE, nz=0; i<na; i++) if ((silminState->liquidComp)[nl][i] == 0.0) { hasNull = TRUE; nullComp[i] = TRUE;
        } else { nullComp[i] = FALSE; nz++; }
        nz += nr - na;

        /* Form the orthogonal projection operator for the equality constraints. The static ne is initialized here. */
        if (hasNull) {

            /* we need a general proceedure for dxdr */
            for (i=0; i<na; i++) for (j=0; j<nr; j++) dxdr[i][j] = 0.0;
            for (j=0; j<nr; j++)  { dxdr[0][j] = -1.0; dxdr[j+1][j] = 1.0; }
            /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

            for (i=0, ne=0; i<na; i++) if(nullComp[i]) { for (j=0; j<nr; j++) dxdr[ne][j] = dxdr[i][j]; ne++; }
            for (i=0; i<ne; i++) householderRowRow(HOUSEHOLDER_CALC_MODE_H1, i, i+1, nr-1, dxdr, i, &hVec[i], dxdr, i+1, ne-1);
        } else ne = 0;

        for (i=0, np=0; i<na; i++) if (!nullComp[i]) {

            for (j=0; j<na; j++) mVec[np][j] = (nullComp[j]) ? 0.0 : 1.0;
            mVec[np][i] *= 10.0*(na-ne);
            conLiq(SECOND, THIRD, t, p, NULL, mVec[np], bVec[np], NULL, NULL, NULL, NULL);

            if (hasNull) {
                for (j=0; j<ne; j++) householderRowRow(HOUSEHOLDER_CALC_MODE_H2, j, j+1, nr-1, dxdr, j, &hVec[j], bVec, np, np);
                for (j=0; j<ne; j++) yVec[j] = (fabs(bVec[np][j]) < 10.0*DBL_EPSILON) ? 0.0 : bVec[np][j];
            }

#ifdef DEBUG
            printf("Making call to vmmin in checkForCoexistingLiquids().\n");
#endif
            noFminfnCalls = 0; noFmingrCalls = 0;
            mode = vmmin(nz, &bVec[np][ne], &Fmin[np], reltest, fminfn, fmingr);
#ifdef DEBUG
            printf("Return from checkForCoexistingLiquids(). Calls to fminfn = %d, calls to fmingr = %d.\n", noFminfnCalls, noFmingrCalls);
#endif

            if (hasNull) {
                for (j=0; j<ne; j++) bVec[np][j] = yVec[j];
                for (j=(ne-1); j>=0; j--) householderRowRow(HOUSEHOLDER_CALC_MODE_H2, j, j+1, nr-1, dxdr, j, &hVec[j], bVec, np, np);
                for (j=0; j<nr; j++) if (fabs(bVec[np][j]) < 10.0*DBL_EPSILON) bVec[np][j] = 0.0;
            }

            if (mode == VMMIN_SUCCESS && Fmin[np] < -reltest*fabs(gRef)) {
                rTr[np]=0.0;
                for (j=0; j<nr; j++) rTr[np] += SQUARE(bVec[np][j]-bRef[j]);
                conLiq(THIRD, FOURTH, t, p, NULL, NULL, bVec[np], mVec[np], NULL, NULL, NULL);
                np++;
            } else if (mode != VMMIN_SUCCESS) {
#ifdef DEBUG
                printf("  FAILURE in vmmin for liquid phase from direction %s\n",
                       liquid[i].label);
                printf("  Initial gRef = %g, composition:\n", gRef);
                for (j=0; j<nr; j++) printf("  r[%1d] = %g",  j, bRef[j]);  printf("\n");
                for (j=0; j<nr; j++) printf("  dg[%1d] = %g", j, dgRef[j]); printf("\n");
                printf("  Fmin = %g\n", Fmin[np]);
                for (j=0; j<nr; j++) printf("  r[%1d] = %g", j, bVec[np][j]); printf("\n");
#endif
            }
        }

        if (np > 0) {
            double inmass;
#ifdef DEBUG
            {
                double *oxVal = (double *) malloc((size_t) nc*sizeof(double));
                double sum;
                int l;

                printf("  Coexisting liquids found.\n");
                printf("    %13.13s %13.13s\n", "Fmin", "rTr");
                for (j=0; j<nc; j++) printf(" %5.5s", bulkSystem[j].label); printf("\n");

                printf("    %13.6g %13.13s\n", gRef, "");
                for (j=0, sum=0.0; j<nc; j++) {
                    for (k=0, oxVal[j]=0.0; k<na; k++) oxVal[j] += bulkSystem[j].mw*(liquid[k].liqToOx)[j]*(silminState->liquidComp)[nl][k];
                    sum += oxVal[j];
                }
                if (sum != 0.0) for (j=0; j<nc; j++) oxVal[j] *= 100.0/sum;
                for (j=0; j<nc; j++) printf(" %5.2f", oxVal[j]); printf("\n");

                for (j=0; j<np; j++) {
                    printf("    %13.6g %13.6g\n", Fmin[j], rTr[j]);
                    for (k=0, sum=0.0; k<nc; k++) {
                        for (l=0, oxVal[k]=0.0; l<na; l++) oxVal[k] += bulkSystem[k].mw*(liquid[l].liqToOx)[k]*mVec[j][l];
                        sum += oxVal[k];
                    }
                    if (sum != 0.0) for (k=0; k<nc; k++) oxVal[k] *= 100.0/sum;
                    for (k=0; k<nc; k++) printf(" %5.2f", oxVal[k]); printf("\n");
                }
                free(oxVal);
            }
#endif
            /* Stop the search through the immiscible solid phases. This insures
             that we only add one new phase at a time which prevents duplicate
             entries. */
            nl = silminState->nLiquidCoexist;
            /* Allocate space to store the new compositional data */

            silminState->liquidComp  = (double **) REALLOC(silminState->liquidComp,  (size_t) (nl+1)*sizeof(double *));
            silminState->liquidDelta = (double **) REALLOC(silminState->liquidDelta, (size_t) (nl+1)*sizeof(double *));
            (silminState->liquidComp)[nl]  = (double *) malloc((size_t) na*sizeof(double));
            (silminState->liquidDelta)[nl] = (double *) malloc((size_t) na*sizeof(double));

            /* Find the composition most distant from the initial composition */
            for (i=0, j=0; i<np; i++) if (rTr[i] > rTr[j]) j = i;
            np = j;

            /* Add the new phase to the system */
            inmass = MASSIN;
            silminState->nLiquidCoexist++;
            for (i=0; i<na; i++) {
                (silminState->liquidComp)[nl][i] = mVec[np][i]*inmass;
#ifdef DEBUG
                printf("%20.20s Before %13.6g After %13.6g New %13.6g\n", liquid[i].label,
                       (silminState->liquidComp)[ 0][i], (silminState->liquidComp)[ 0][i]-(silminState->liquidComp)[nl][i],
                       (silminState->liquidComp)[nl][i]);
#endif
                (silminState->liquidComp)[ 0][i] -= (silminState->liquidComp)[nl][i];
            }
#ifdef DEBUG
            if (!testLiq(SIXTH, t, p, 0, 0, NULL, NULL, NULL, (silminState->liquidComp)[ 0])) printf("... primary liquid infeasible in check_coexisting_liquid!\n");
            if (!testLiq(SIXTH, t, p, 0, 0, NULL, NULL, NULL, (silminState->liquidComp)[nl])) printf("... derived liquid infeasible in check_coexisting_liquid!\n");
#endif
            result = SUCCESS;
        }
    } /* loop on nl (number of coexisting liquids             */

#ifdef DEBUG
    printf("Normal exit from checkForCoexistingLiquids\n");
#endif

    return result;
}

/*****************************************************************************
 * Private function declarations:
 *****************************************************************************/

static double fminfn(        /* returned function value, if notcomp == FALSE */
                     int n,         /* number of independent variables                          */
                     double *bVec,  /* vector of independent variables, length n                */
                     int *notcomp)  /* returned flag, TRUE is current parameters are infeasible */
{
    double result;
    int i;

    noFminfnCalls++;

    for (i=0; i<nz; i++) tVec[i+ne] = bVec[i];
    if (hasNull) {
        for (i=0; i<ne; i++) tVec[i] = yVec[i];
        for (i=(ne-1); i>=0; i--) householderRowRow(HOUSEHOLDER_CALC_MODE_H2, i, i+1, nr-1, dxdr, i, &hVec[i], &tVec, 0, 0);
        for (i=0; i<nr; i++) if (fabs(tVec[i]) < 10.0*DBL_EPSILON) tVec[i] = 0.0;
    }

    if (silminState->fo2Path != FO2_NONE) {
        int j;
        double *mOxd = (double *) malloc((size_t)  nc*sizeof(double));
        double *mCmp = (double *) malloc((size_t) nlc*sizeof(double));

        conLiq(THIRD, FOURTH, t, p, NULL, NULL, tVec, mCmp, NULL, NULL, NULL);
        for (i=0; i<nc; i++) for (j=0, mOxd[i]=0.0; j<nlc; j++) mOxd[i] += mCmp[j]*(liquid[j].liqToOx)[i];
        conLiq(FIRST | SEVENTH, FIRST, t, p, mOxd, NULL, NULL, NULL, NULL, NULL, &(silminState->fo2));
        for (i=0; i<nlc; i++) for (j=0, mCmp[i]=0.0; j<nc; j++) mCmp[i] += mOxd[j]*(bulkSystem[j].oxToLiq)[i];
        free(mOxd);

        if (!testLiq(SIXTH, t, p, 0, 0, NULL, NULL, NULL, mCmp)) {
#ifdef DEBUG
            printf("... Infeasible fO2 corrected composition in fminfn (check_coexisting_liquids.c) at iteration %d\n", noFminfnCalls);
#endif
            *notcomp = TRUE;
            free(mCmp);
            return 0.0;
        }
        conLiq(SECOND, THIRD, t, p, NULL, mCmp, tVec, NULL, NULL, NULL, NULL);
        free(mCmp);
    } else if (!testLiq(FIFTH, t, p, 0, 0, NULL, NULL, tVec, NULL)) {
#ifdef DEBUG
        printf("... Infeasible composition in fminfn (check_coexisting_liquids.c) at iteration %d\n", noFminfnCalls);
#endif
        *notcomp = TRUE;
        return 0.0;
    }
    *notcomp = FALSE;

    gmixLiq(FIRST, t, p, tVec, &result, NULL, NULL);

    result -= gRef;
    for (i=0; i<nr; i++) result -= dgRef[i]*(tVec[i]-bRef[i]);

    return result;
}

/*****************************************************************************/

static void fmingr(
                   int n,              /* number of independent variables                     */
                   double *bVec,       /* vector of independent variables, length n           */
                   double *g)          /* returned gradient vector                            */
{
    int i;

    noFmingrCalls++;

    gmixLiq(SECOND, t, p, tVec, NULL, dg, NULL);
    for (i=0; i<nr; i++) dg[i] -= dgRef[i];

    if (hasNull) for (i=0; i<ne; i++) householderRowRow(HOUSEHOLDER_CALC_MODE_H2, i, i+1, nr-1, dxdr, i, &hVec[i], &dg, 0, 0);

    for (i=0; i<nz; i++) g[i] = dg[i+ne];
}
#undef REALLOC
#undef SQUARE

/* end of file CHECK_COEXISTING_LIQUIDS.C */
