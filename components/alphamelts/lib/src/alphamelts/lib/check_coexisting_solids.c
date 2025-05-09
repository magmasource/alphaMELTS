const char *check_coexisting_solids_ver(void) { return "$Id: check_coexisting_solids.c,v 1.2 2006/08/17 16:47:18 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: check_coexisting_solids.c,v $
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
 * Revision 3.7  1997/06/21  22:50:08  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.6  1997/05/03  20:23:43  ghiorso
 * *** empty log message ***
 *
 * Revision 3.5  1997/03/27  17:03:46  ghiorso
 * *** empty log message ***
 *
 * Revision 3.4  1996/09/24  20:33:52  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 3.3  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 3.2  1995/11/01  22:40:27  ghiorso
 * Implementation of subsolidus options after Asimow.
 * Additional implementation of nepheline solid solutions.
 *
 * Revision 3.2  1995/11/01  22:40:27  ghiorso
 * Implementation of subsolidus options after Asimow.
 * Additional implementation of nepheline solid solutions.
 *
 * Revision 3.1  1995/08/18  17:20:27  ghiorso
 * MELTS Version 3 - Initial Entry
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
**      (file: CHECK_COEXISTING_SOLIDS.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  October 12, 1991  Original Version
**      V1.0-2  Mark S. Ghiorso  October 14, 1991
**              (1) Added equality constraints for absent endmember components
**      V1.0-3  Mark S. Ghiorso  October 17, 1991
**              (1) First complete version of logic for detecting and
**                  adding an immiscible phase
**      V1.0-4  Mark S. Ghiorso  November 4, 1991
**              Added test on backprojected solution vector from the null
**              space constraints to evaluate whether values are sufficiently
**              close to zero (+/- 10 machine precision).
**      V1.0-5  Mark S. Ghiorso  January 28, 1992
**              changed index to Index to avoid conflict with UNIX string
**              function index()
**      V1.0-6  Mark S. Ghiorso  March 23, 1992
**              Replaced call to pow() function with macro SQUARE
**      V1.1-1  Mark S. Ghiorso  June 5, 1992
**              Corrected initialization of range-space vector for zero
**                concentrations of endmembers. Error was causing erroneous
**                compositions for centro-symmetric solid solutions.
**      V1.1-2  Mark S. Ghiorso  September 29, 1993
**              Modified call to realloc to catch zero pointer (SPARC port)
**      V2.1-0  Paul D. Asimow   August 3, 1994
**              Begin optimizing for speed in realistic mantle calculations:
**              Turn off everything but pyroxene
**      V2.2-0  Paul D. Asimow  August 12, 1994
**              Restrict pyroxene search directions to cpx/opx to save time
**      V3.0-1  Paul D. Asimow  March 22, 1995
**              Add solids-only capability
**      V3.1-1  Paul D. Asimow  August 7, 1995
**              Add subsolidus fO2 buffering (compatibility changes only)
**      V3.2-1  Paul D. Asimow  September 7, 1995
**              reject 2.1-0
**--
*/
//#include "melts_gsl.h"

#include "silmin.h"
#include "nash.h"
#include "lawson_hanson.h"

#define SQUARE(x) ((x)*(x))
#define REALLOC(x, y) (((x) == NULL) ? malloc(y) : realloc((x), (y)))

/* If change this would need to set up gsl matrix like in check coexisting liquids */
#define _NO_SPINODE_TEST

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
static int      Index;    /* solid phase index                                */
static int      na;       /* number of endmembers in phase index              */
static int      ne;       /* number of equality constraints                   */
static int      nr;       /* number of independent variables in phase index   */
static int      nz;       /* number of non-zero endmember concentrations      */
static int     *nullComp; /* vector of TRUE/FALSE flags for zero concen       */
static double   p;        /* pressure (bars)                                  */
static double   t;        /* temperature (K)                                  */
static double  *tVec;     /* temporary soln for fminfn and fmingr routines    */
static double  *yVec;     /* range space soln if prob has eq constr           */

/*****************************************************************************
 * Private function definitions:
 *****************************************************************************/

static double fminfn(int n, double *bVec, int *notcomp);
static void   fmingr(int n, double *bVec, double *g);
static void   initializeGlobalStatics(void);

/*****************************************************************************
 * Global function declarations:
 *****************************************************************************/

#define SUCCESS TRUE
#define FAILURE FALSE

#ifdef SPINODE_TEST
static int checkpd(gsl_matrix *a, int n) {
    gsl_matrix_view p = gsl_matrix_submatrix(a, (size_t) 0, (size_t) 0, (size_t) n, (size_t) n);
    int status;

    (void) gsl_set_error_handler_off();
    status = gsl_linalg_cholesky_decomp1(&p.matrix);
    (void) gsl_set_error_handler(NULL);

    return (status == GSL_SUCCESS);
}
}
#endif

int spinodeTest(void)
{
  int result, Index, ns, i;

  if (dgRef == NULL) initializeGlobalStatics();

  t       = silminState->T; if (t <= 0.0) return FAILURE;
  p       = silminState->P; if (p <  0.0) return FAILURE;
  result  = SUCCESS;

  for (Index=0; Index<npc; Index++) if (solids[Index].na > 1) {
    for (ns=0; ns<(silminState->nSolidCoexist)[Index]; ns++) {
      double *mVec;
      na = solids[Index].na;
      nr = solids[Index].nr;

      mVec  = (double *) malloc((unsigned) na*sizeof(double));
      for (i=0; i<na; i++) mVec[i] = (silminState->solidComp)[Index+1+i][ns];
      (*solids[Index].convert)(SECOND, THIRD, t, p, NULL, mVec, bRef, NULL, NULL, NULL, NULL, NULL);
      free(mVec);

      (*solids[Index].gmix)(FIRST | SECOND | THIRD, t, p, bRef, &gRef, dgRef, gHess, NULL);

#ifdef SPINODE_TEST
      /* Determine if the Hessian is positive definite */
      if (!checkpd(gHess, nr)) {
#ifdef DEBUG
        printf("Check for coexisting solids called inside spinodal of %s!\n", solids[Index].label);
#endif
        result = FAILURE;
      }
#endif

    }
  }
  return result;
}


int checkForCoexistingSolids(  /* returns a MODE flag for success or failure */
  void)
{
  static double **bVec, *Fmin, **mVec, *rTr;
  double reltest, inmass;
  int i, j, mode, np, ns, result, acceptable;
  int hasLiquid = (silminState->liquidMass != 0.0);

#ifdef DEBUG
  printf("Call to checkForCoexistingSolids\n");
#endif

  /* Global statics */
  if (dgRef == NULL) initializeGlobalStatics();
  /* Local statics */
  if (bVec  == NULL) {
    bVec      = (double **) malloc((unsigned) na*sizeof(double *));
    for (i=0; i<na; i++)
      bVec[i] = (double *)  malloc((unsigned) nr*sizeof(double));
    Fmin      = (double *)  malloc((unsigned) na*sizeof(double));
    mVec      = (double **) malloc((unsigned) na*sizeof(double *));
    for (i=0; i<na; i++)
      mVec[i] = (double *)  malloc((unsigned) na*sizeof(double));
    rTr	      = (double *)  malloc((unsigned) na*sizeof(double));
  }

  t       = silminState->T; if (t <= 0.0) return FAILURE;
  p       = silminState->P; if (p <  0.0) return FAILURE;
  reltest = sqrt(DBL_EPSILON);
  result  = FAILURE;

  for (Index=0; Index<npc; Index++) if (solids[Index].na > 1) {
#ifdef PHMELTS_ADJUSTMENTS
  /* solids[Index].type == PHASE if na > 1 */
  for (i=0, j=0; i<npc; i++) if (solids[i].type == PHASE) { if (i == Index) break; else j++; }
#endif
#ifdef RHYOLITE_ADJUSTMENTS
    if (!strcmp(solids[Index].label, "orthopyroxene")) {
#ifdef DEBUG
    printf("Check for coexisting solids is blocked for phase %s\n", solids[Index].label);
#endif
#ifndef TESTDYNAMICLIB
    } else if (!strcmp(solids[Index].label, "alkali-feldspar") && ((silminState->nSolidCoexist)[Index] >= 2)) {
#ifdef DEBUG
    printf("Check for coexisting solids is blocked for phase %s which already has %d entries.\n", solids[Index].label, (silminState->nSolidCoexist)[Index]);
#endif
    } else if (!strcmp(solids[Index].label, "plagioclase") && ((silminState->nSolidCoexist)[Index] >= 2)) {
#ifdef DEBUG
    printf("Check for coexisting solids is blocked for phase %s which already has %d entries.\n", solids[Index].label, (silminState->nSolidCoexist)[Index]);
#endif
#else
    } else if (!strcmp(solids[Index].label, "feldspar") && ((silminState->nSolidCoexist)[Index] >= 2)) {
#ifdef DEBUG
    printf("Check for coexisting solids is blocked for phase %s which already has %d entries.\n", solids[Index].label, (silminState->nSolidCoexist)[Index]);
#endif
#endif /* TESTDYNAMICLIB */
    } else if (!strcmp(solids[Index].label, "clinopyroxene") && ((silminState->nSolidCoexist)[Index] >= 3)) {
#ifdef DEBUG
    printf("Check for coexisting solids is blocked for phase %s which already has %d entries.\n", solids[Index].label, (silminState->nSolidCoexist)[Index]);
#endif
#ifdef PHMELTS_ADJUSTMENTS
    } else if ((silminState->nSolidCoexist)[Index] >= (silminState->incSolids)[j]) {
#ifdef DEBUG
    printf("Check for coexisting solids is blocked for phase %s which already has %d entries.\n", solids[Index].label, (silminState->nSolidCoexist)[Index]);
#endif
#endif /* PHMELTS_ADJUSTMENTS */
    } else {
#endif
#ifdef DEBUG
    printf("Check for coexisting solids called for phase %s\n", solids[Index].label);
#endif
    for (ns=0; ns<(silminState->nSolidCoexist)[Index]; ns++) {
      na = solids[Index].na;
      nr = solids[Index].nr;

      for (i=0; i<na; i++) mVec[0][i] = (silminState->solidComp)[Index+1+i][ns];
      (*solids[Index].convert)(SECOND, THIRD, t, p, NULL, mVec[0], bRef, NULL, NULL, NULL, NULL, NULL);
      (*solids[Index].gmix)(FIRST | SECOND | THIRD, t, p, bRef, &gRef, dgRef, gHess, NULL);

#ifdef SPINODE_TEST
      /* Determine if the Hessian is positive definite */
      if (!choldc(gHess, nr)) {
#ifdef DEBUG
        printf("Check for coexisting solids called inside spinodal!\n");
#endif
        return result;
      }
#endif

      /* Determine if any enemembers have zero concentration */
      for (i=0, hasNull=FALSE, nz=0; i<na; i++)
        if ((silminState->solidComp)[Index+1+i][ns] == 0.0) { hasNull = TRUE; nullComp[i] = TRUE; } else { nullComp[i] = FALSE; nz++; }
      nz += nr - na;

      /* Form the orthogonal projection operator for the equality constraints.
         the static ne is initialized here.                                   */
      if (hasNull) {
        (*solids[Index].convert)(THIRD, SEVENTH, t, p, NULL, NULL, bRef, NULL, NULL, NULL, dxdr, NULL);
        for (i=0, ne=0; i<na; i++) if(nullComp[i]) { for (j=0; j<nr; j++) dxdr[ne][j] = dxdr[i][j]; ne++; }
        for (i=0; i<ne; i++) householderRowRow(HOUSEHOLDER_CALC_MODE_H1, i, i+1, nr-1, dxdr, i, &hVec[i], dxdr, i+1, ne-1);
      } else ne = 0;

      for (i=0, np=0; i<na; i++) if (!nullComp[i]) {

        for (j=0; j<na; j++) mVec[np][j] = (nullComp[j]) ? 0.0 : 1.0;
        mVec[np][i] *= 10.0*(na-ne);
        (*solids[Index].convert)(SECOND, THIRD, t, p, NULL, mVec[np], bVec[np], NULL, NULL, NULL, NULL, NULL);

        if (hasNull) {
          for (j=0; j<ne; j++) householderRowRow(HOUSEHOLDER_CALC_MODE_H2, j, j+1, nr-1, dxdr, j, &hVec[j], bVec, np, np);
          for (j=0; j<ne; j++) yVec[j] = (fabs(bVec[np][j]) < 10.0*DBL_EPSILON) ? 0.0 : bVec[np][j];
        }

        mode = vmmin(nz, &bVec[np][ne], &Fmin[np], reltest, fminfn, fmingr);
        if (hasNull) {
          for (j=0; j<ne; j++) bVec[np][j] = yVec[j];
          for (j=(ne-1); j>=0; j--) householderRowRow(HOUSEHOLDER_CALC_MODE_H2, j, j+1, nr-1, dxdr, j, &hVec[j], bVec, np, np);
          for (j=0; j<nr; j++) if (fabs(bVec[np][j]) < 10.0*DBL_EPSILON) bVec[np][j] = 0.0;
        }

        if (mode == VMMIN_SUCCESS && Fmin[np] < -reltest*fabs(gRef)) {
          rTr[np]=0.0;
          for (j=0; j<nr; j++) rTr[np] += SQUARE(bVec[np][j]-bRef[j]);
          (*solids[Index].convert)(THIRD, FOURTH, t, p, NULL, NULL, bVec[np], mVec[np], NULL, NULL, NULL, NULL);
          np++;
        } else if (mode != VMMIN_SUCCESS) {
#ifdef DEBUG
  printf("  FAILURE in vmmin for phase %s from direction %s\n", solids[Index].label, solids[Index+1+i].label);
  printf("  Initial gRef = %g, composition:\n", gRef);
  for (j=0; j<nr; j++) printf("  r[%1d] = %g",  j, bRef[j]);  printf("\n");
  for (j=0; j<nr; j++) printf("  dg[%1d] = %g", j, dgRef[j]); printf("\n");
  printf("  Fmin = %g\n", Fmin[np]);
  for (j=0; j<nr; j++) printf("  r[%1d] = %g", j, bVec[np][j]); printf("\n");
#endif
        }
      }

      if (np > 0) {
#ifdef DEBUG
  printf("  Coexisting solids found for phase %s.\n", solids[Index].label);
  printf("    %13.13s %13.13s", "Fmin", "rTr");
  for (j=0; j<na; j++) printf(" %9.9s", solids[Index+1+j].label); printf("\n");
  printf("    %13.6g %13.13s", gRef, "");
  for (j=0; j<na; j++) printf("%10.6f", (silminState->solidComp)[Index+1+j][ns]/(silminState->solidComp)[Index][ns]);
  printf("\n");
  for (j=0; j<np; j++) {
    int k;
    printf("    %13.6g %13.6g", Fmin[j], rTr[j]);
    for (k=0; k<na; k++) printf("%10.6f", mVec[j][k]); printf("\n");
  }
#endif
        /* Stop the search through the immiscible solid phases. This insures
           that we only add one new phase at a time which prevents duplicate
           entries. */
        ns = (silminState->nSolidCoexist)[Index];
        /* Allocate space to store the new compositional data */
        for (i=0; i<=na; i++) {
          (silminState->solidComp)[Index+i]  = (double *) REALLOC((silminState->solidComp)[Index+i],  (size_t) (ns+1)*sizeof(double));
          (silminState->solidDelta)[Index+i] = (double *) REALLOC((silminState->solidDelta)[Index+i], (size_t) (ns+1)*sizeof(double));
        }
        /* Find the composition most distant from the initial composition */
        for (i=0, j=0; i<np; i++) if (rTr[i] > rTr[j]) j = i;
        np = j;
        /* Add the new phase to the system */
        inmass = MASSIN; acceptable = FALSE;
        if (hasLiquid) {  /* test to see whether withdrawing from liquid will cause crash (use first liquid) */
          while (!acceptable) {
            double *dummyComp = (double *) malloc((size_t) nlc*sizeof(double));
            for (i=0; i<nlc; i++) dummyComp[i] = silminState->liquidComp[0][i];
            acceptable = TRUE;
            for (i=0; i<na; i++) for (j=0; j<nlc; j++) dummyComp[j] -= (solids[Index+1+i].solToLiq)[j]*mVec[np][i]*inmass;
	    acceptable = testLiq(SIXTH, t, p, 0, 0, NULL, NULL, NULL, dummyComp);
            if (acceptable == FALSE) inmass *= 0.5;
            free(dummyComp);
          }
        }
        (silminState->nSolidCoexist)[Index]++;
        (silminState->solidComp)[Index][ns] = inmass;
        if (!hasLiquid) (silminState->solidComp)[Index][ns-1] -= inmass;
        for (i=0; i<na; i++) {
          (silminState->solidComp)[Index+1+i][ns] = mVec[np][i]*inmass;
          if (hasLiquid) for (j=0; j<nlc; j++) (silminState->liquidComp)[0][j] -= (solids[Index+1+i].solToLiq)[j]*(silminState->solidComp)[Index+1+i][ns];
          else (silminState->solidComp)[Index+1+i][ns-1] -= mVec[np][i]*inmass;
        }
        if (hasLiquid) {
          for (i=0, silminState->liquidMass=0.0; i<nlc; i++) for (j=0; j<nc; j++) silminState->liquidMass +=
            silminState->liquidComp[0][i]*(liquid[i].liqToOx)[j]*bulkSystem[j].mw;
        }
        result = SUCCESS;
      }

    } /* loop on ns (number of coexisting solids             */
#ifdef RHYOLITE_ADJUSTMENTS
    } /* end if block on coexisting phase exceptions */
#endif
  }   /* loop on Index (solids with more than one endmember) */

  return result;
}

/*****************************************************************************
 * Private function declarations:
 *****************************************************************************/

static void initializeGlobalStatics(void) {
  int i;

  for (i=0, nr=0, na=0; i<npc; i++) {
    if (solids[i].nr > nr) nr = solids[i].nr;
    if (solids[i].na > na) na = solids[i].na;
  }

  bRef      = (double *)  malloc((unsigned) nr*sizeof(double));
  dg	    = (double *)  malloc((unsigned) nr*sizeof(double));
  dgRef     = (double *)  malloc((unsigned) nr*sizeof(double));
  gHess     = (double **) malloc((unsigned) nr*sizeof(double *));
  for (i=0; i<nr; i++)
    gHess[i] = (double *) malloc((unsigned) nr*sizeof(double));
  dxdr      = (double **) malloc((unsigned) na*sizeof(double *));
  for (i=0; i<na; i++)
    dxdr[i] = (double *)  malloc((unsigned) nr*sizeof(double));
  hVec      = (double *)  malloc((unsigned) nr*sizeof(double));
  nullComp  = (int *)	  malloc((unsigned) na*sizeof(int));
  tVec      = (double *)  malloc((unsigned) nr*sizeof (double));
  yVec      = (double *)  malloc((unsigned) nr*sizeof(double));
}

static double fminfn(        /* returned function value, if notcomp == FALSE */
  int n,         /* number of independent variables                          */
  double *bVec,  /* vector of independent variables, length n                */
  int *notcomp)  /* returned flag, TRUE is current parameters are infeasible */
{
  double result;
  int i;

  for (i=0; i<nz; i++) tVec[i+ne] = bVec[i];
  if (hasNull) {
    for (i=0; i<ne; i++) tVec[i] = yVec[i];
    for (i=(ne-1); i>=0; i--) householderRowRow(HOUSEHOLDER_CALC_MODE_H2, i, i+1, nr-1, dxdr, i, &hVec[i], &tVec, 0, 0);
    for (i=0; i<nr; i++) if (fabs(tVec[i]) < 10.0*DBL_EPSILON) tVec[i] = 0.0;
  }

  if (!(*solids[Index].test)(FIFTH, t, p, 0, 0, NULL, NULL, tVec, NULL)) { *notcomp = TRUE; return 0.0; }
  *notcomp = FALSE;

  (*solids[Index].gmix)(FIRST, t, p, tVec, &result, NULL, NULL, NULL);

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

  (*solids[Index].gmix)(SECOND, t, p, tVec, (double *) NULL, dg, NULL, NULL);
  for (i=0; i<nr; i++) dg[i] -= dgRef[i];

  if (hasNull) for (i=0; i<ne; i++) householderRowRow(HOUSEHOLDER_CALC_MODE_H2, i, i+1, nr-1, dxdr, i, &hVec[i], &dg, 0, 0);

  for (i=0; i<nz; i++) g[i] = dg[i+ne];
}
#undef REALLOC
#undef SQUARE

/* end of file CHECK_COEXISTING_SOLIDS.C */
