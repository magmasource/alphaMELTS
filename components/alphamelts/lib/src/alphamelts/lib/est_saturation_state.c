const char *est_saturation_state_ver(void) { return "$Id: est_saturation_state.c,v 1.4 2009/04/16 16:35:23 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: est_saturation_state.c,v $
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
MELTS Source Code: RCS Revision 1.1.1.1  2004/01/02 19:21:49  cvsaccount
MELTS Source Code: RCS CTserver University of Chicago
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:58:12  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
 * Revision 3.7  1997/06/21  22:49:58  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.6  1997/05/03  20:23:36  ghiorso
 * *** empty log message ***
 *
 * Revision 3.5  1997/03/27  17:03:40  ghiorso
 * *** empty log message ***
 *
 * Revision 3.4  1996/09/24  20:33:42  ghiorso
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
 * Revision 3.1  1995/08/18  17:30:23  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Routines to estimate the composition and chemical affinity for
**      saturation with a solid solution at a particular T and P.
**      (file: EST_SATURATION_STATE.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  August 14, 1991  Original Version
**                               code completed August 16, 1991
**      V1.0-2  Mark S. Ghiorso  August 21, 1991
**                               final debugging  of working version
**      V1.0-3  Mark S. Ghiorso  September 14, 1991
**              Changed call to (*solids[].convert) to reflect new
**              parameter list
**      V2.0-1  Mark S. Ghiorso  September 24, 1991
**              Reorganized and expanded code to allow for zero concentrations
**              of endmembers. Revised argument list for (*solids.convert)
**              routines
**      V2.0-2  Mark S. Ghiorso  September 26, 1991
**              Revised initial guess algorithm to avoid zero endmembers
**      V3.0-1  Mark S. Ghiorso  September 28, 1991
**              Revised code to detect and correct non-zero sum-of-squares
**              conditions by providing alternative initial guesses
**      V3.0-2  Mark S. Ghiorso  October 1, 1991
**              Implemented conjugate gradient search and minimization
**              algorithm to get initial guess out of the spinoidal region
**      V3.0-3  Mark S. Ghiorso  October 7, 1991
**              Revised algorithm for guess in spinoidal regions for
**              searching in two directions, if the first attempt fails
**      V4.0-1  Mark S. Ghiorso  October 8, 1991
**              Revised modmrt code for equality constraints and reorganized
**              structure to accomodate these modifications
**      V4.0-2  Mark S. Ghiorso  November 4, 1991
**              (1) Upon backprojection of the null space solution vector into
**                  the range space, elements are checked to so that non-zero
**                  quantities on the order of the machine precision are zeroed
**                  explicitly.
**              (2) in computeSolnWithinSolvus, changed double loop to nz
**                  in returned hessian to double loop to (nz-1)
**      V4.0-3  Mark S. Ghiorso March 4, 1992
**              Added code to scale the hessian in computeSolnWithinSolvus
**              in order to prevent overflow in the call to rqmcg()
**      V4.1-1  Mark S. Ghiorso  June 5, 1992
**              Corrected initialization of range-space vector for zero
**                concentrations of endmembers. Error was causing erroneous
**                compositions for centro-symmetric solid solutions.
**      V4.1-2  Mark S. Ghiorso  July 22, 1992
**              (1) Added warning message for return FAILURE of external
**                  function
**              (2) Added compilation flag (SPECIAL) to allow a successful
**                  return if the calculated affinity is negative but
**                  the algorithm indicates non-convergence
**      V4.2-1  Mark S. Ghiorso  March 23, 1995
**              (1) Added failure return to getAffinityAndComposition to deal
**                  with the case that nz == na
**      V5.0-1  Paul D. Asimow  March 23, 1995
**              (1) fix for nz==na suggested by MSG
**              (2) enable liquid saturation estimate for subsolidus assemblage
**      V5.0-2  Paul D. Asimow August 7, 1995
**              (1) Compatibility changes for subsolidus buffering
**--
*/

#include "silmin.h"
#include "nash.h"
#include "lawson_hanson.h"

#ifdef DEBUG
#undef DEBUG
#endif

#ifdef SPECIAL  /* If defined, return SUCCESS if affinity is negative     */
#undef SPECIAL  /* even if convergence criteria (Fmin < sqrt(DBL_EPSILON) */
#endif          /* is not satisfied                                       */

/*****************************************************************************
 * External variables whose scope is confined to this file:
 *****************************************************************************/

static double *activity;     /* temporary storage for endmember activities    */
static double **dadr;        /* temporary storage for d(activity)/d(indep)    */
static double **dxdr;        /* temporary storage for d(mole frac)/d(indep)   */
static double **d2gdr2;      /* matrix of second derivatives of G mixing      */
static double *finalR;       /* r vector which min Gmix in solvus search      */
static int    hasNull;       /* TRUE if some endmembers have zero concen      */
static double *hVec;         /* pivot elements for H12 decomp of eq constr    */
static double **identity;    /* the identity matrix                           */
static double *moleFrac;     /* temporray storage for endmember mole fractions*/
static double *mu;           /* vector of chemical potentials of end-members  */
static int     na;           /* number of end-members in solid solution       */
static int     ne;           /* number of equality constraints                */
static int     nr;           /* number of indepen parameters in solid solution*/
static int     nz;           /* number of non-zero endmember concentrations   */
static int    *nullComp;     /* vector of TRUE/FALSE flags for zero concen    */
static int    *nullList;     /* list of indices for non-nulls of length nz    */
static double *refR;         /* refernece r vector for solvus search          */
static double *searchR;      /* search direction for solvus search            */
static int     solidID = -1; /* index number of solid phase in solids[0:npc]  */
static double  t;            /* temperature (K)                               */
static double *tVec;         /* temporary soln for func and grad routines     */
static double  p;            /* pressure (bars)                               */
static double *yVec;         /* range space soln if prob has eq constr        */
static int     liquidMode;   /* TRUE if we are looking for solidus            */

/*****************************************************************************
 * Private function definitions:
 *****************************************************************************/

static double func(int i, int n, double *bVec, int *notcomp);
static void   grad(int i, int n, double *bVec, double *X);
static double gmix(double lambda, int *notcomp);
static void   initialGuess(double *bVec);
static int    computeSolnWithinSolvus(double *bVec, double *Fmin);

/*****************************************************************************
 * Global function declarations:
 *****************************************************************************/

#define SCALE    1000.0  /* Scaling factor for chemical affinities           */

#define SUCCESS TRUE
#define FAILURE FALSE

int getAffinityAndComposition( /* Returns a MODE flag for success or failure */
    double temp,        /* temperature (K)                                     */
    double pres,        /* pressure (bars)                                     */
    int index,          /* index of solid phase in the solids[] structure      */
                                            /* -1 indicates liquid                                 */
    int    *zeroX,      /* TRUE if endmember component has zero mole fraction  */
    double *muMinusMu0, /* vector of end-member mu - mu0 i.e. A + RTln(a)      */
    double *affinity,   /* returned value, chemical affinity (J)               */
    double *indepVar)   /* returned vector, composition of phase (length nr)   */
{
    int i, j, maxIter, mode;
    static double *bVec = NULL;
    double Fmin, reltest;

    /* Check static vector size allocation */
    if (dadr == NULL) {
        activity  = (double *)  malloc((unsigned) nlc*sizeof (double));
        bVec      = (double *)  malloc((unsigned) nlc*sizeof (double));
        dadr      = (double **) malloc((unsigned) nlc*sizeof (double *));
        dxdr      = (double **) malloc((unsigned) nlc*sizeof (double *));
        for (i=0; i<nlc; i++) {
            dadr[i] = (double *)  malloc((unsigned) (nlc-1)*sizeof (double));
            dxdr[i] = (double *)  malloc((unsigned) (nlc-1)*sizeof (double));
        }
        d2gdr2    = (double **) malloc((unsigned) (nlc-1)*sizeof (double *));
        finalR    = (double *)  malloc((unsigned) (nlc-1)*sizeof (double));
        hVec      = (double *)  malloc((unsigned) (nlc-1)*sizeof (double));
        identity  = (double **) malloc((unsigned) (nlc-1)*sizeof (double *));
        for (i=0; i<(nlc-1); i++) {
            d2gdr2[i]   = (double *)  malloc((unsigned) (nlc-1)*sizeof (double));
            identity[i] = (double *)  malloc((unsigned) (nlc-1)*sizeof (double));
            for (j=0; j<(nlc-1); j++) identity[i][j] = 0.0;
            identity[i][i] = 1.0;
        }
        moleFrac  = (double *)  malloc((unsigned) nlc*sizeof (double));
        mu        = (double *)  malloc((unsigned) nlc*sizeof (double));
        nullComp  = (int *)     malloc((unsigned) nlc*sizeof (int));
        nullList  = (int *)     malloc((unsigned) nlc*sizeof (int));
        refR      = (double *)  malloc((unsigned) (nlc-1)*sizeof (double));
        searchR   = (double *)  malloc((unsigned) (nlc-1)*sizeof (double));
        tVec      = (double *)  malloc((unsigned) (nlc-1)*sizeof (double));
        yVec      = (double *)  malloc((unsigned) (nlc-1)*sizeof (double));
    }

    /* Test input parameters */
    t       = temp;  if (t <= 0.0)                   return FAILURE;
    p       = pres;  if (p <  0.0)                   return FAILURE;
                   if (index < -1 || index >= npc) return FAILURE;

    /* Check parameters for algorithmic assumptions */
    liquidMode = (index == -1);
    if (!liquidMode) {
        na = solids[index].na; nr = solids[index].nr;
        if (na != (nr+1)) return FAILURE;
        solidID = index;
    } else {
        na = nlc; nr = nlc-1;
    }

    /* remove excluded components for computation of mole fraction estimates.
     the static nz is initialized here.                                    */
    for (i=0, hasNull=FALSE, nz=0; i<na; i++) {
        nullComp[i] = zeroX[i]; hasNull |= zeroX[i];
        if (!nullComp[i]) { nullList[nz] = i; mu[nz++] = muMinusMu0[i]; }
    }
    if (nz == 0) return FAILURE;

    /* obtain an initial guess for the solution */
    initialGuess(bVec);

    /* Form the orthogonal projection operator for the equality constraints.
     the static ne is initialized here.                                    */
    if (hasNull) {
        if (!liquidMode)
            (*solids[solidID].convert)(THIRD, SEVENTH, t, p, (double *) NULL,
                (double *) NULL, bVec, (double *) NULL, (double **) NULL,
                (double ***) NULL, dxdr, (double ****) NULL);
        else
            for (i=0;i<na;i++) for (j=0;j<nr;j++) dxdr[i][j] = (double) (i == j+1);
        for (i=0, ne=0; i<na; i++) if(nullComp[i])
            { for (j=0; j<nr; j++) dxdr[ne][j] = dxdr[i][j]; ne++; }
        for (i=0; i<ne; i++) householderRowRow(HOUSEHOLDER_CALC_MODE_H1, i, i+1,
            nr-1, dxdr, i, &hVec[i], dxdr, i+1, ne-1);
    } else ne = 0;

    reltest = 1.0e-6;
    maxIter = 250;

    /* Project the initial guess into the Null space of constraints */
    if (hasNull) {
        for (i=0; i<ne; i++) householderRowRow(HOUSEHOLDER_CALC_MODE_H2, i, i+1,
            nr-1, dxdr, i, &hVec[i], &bVec, 0, 0);
        for (i=0; i<ne; i++) yVec[i] = (fabs(bVec[i]) < 10.0*DBL_EPSILON) ? 0.0 :
            bVec[i];
    }

    if((mode = modmrt(nz, nz, &bVec[ne], &Fmin, reltest, maxIter, func, grad))
        == MODMRT_SUCCESS) {
        /* Reassemble solution if equality constraints are present */
        if (hasNull) {
            for (i=0; i<ne; i++) bVec[i] = yVec[i];
            for (i=(ne-1); i>=0; i--) householderRowRow(HOUSEHOLDER_CALC_MODE_H2,
                i, i+1, nr-1, dxdr, i, &hVec[i], &bVec, 0, 0);
            for (i=0; i<nr; i++) if (fabs(bVec[i]) < 10.0*DBL_EPSILON) bVec[i] = 0.0;
        }
        /* Test if solution is on the wrong side of a solvus */
        if (Fmin > sqrt(DBL_EPSILON)) {
#ifdef DEBUG
            if (!liquidMode)
                printf("FAILURE for %s in getAffinityAndComposition. Refining guess.\n",
                    solids[index].label);
            else printf("FAILURE for liquid in getAffinityAndComposition. Refining guess.\n");
            printf("  T (K) = %8.2f, P (bars) = %8.1f, Fmin = %13.6g\n", t, p, Fmin);
            for (i=0; i<na; i++) printf("  bVec[%d] = %f, mu-mu0[%d] = %f\n",
                i, bVec[i], i, muMinusMu0[i]);
            if (!liquidMode) {
                char *formula;
                (*solids[solidID].display)(FIRST, t, p, bVec, &formula);
                printf("  Formula: %s\n", formula);
	free(formula);
            } else
                for (i=0;i<nlc;i++) printf("%s %f\n", liquid[i].label, bVec[i]);
#endif
            if (computeSolnWithinSolvus(bVec, &Fmin) == FAILURE) {
#ifdef DEBUG
                printf("  SECOND FAILURE after attempt with refined initial guess.\n");
                printf("  T (K) = %8.2f, P (bars) = %8.1f, Fmin = %13.6g\n", t, p,
                    Fmin);
                for (i=0; i<na; i++) printf("  bVec[%d] = %f, mu-mu0[%d] = %f\n",
                    i, bVec[i], i, muMinusMu0[i]);
                if (!liquidMode) {
                    char *formula;
                    (*solids[solidID].display)(FIRST, t, p, bVec, &formula);
                    printf("return FAILURE (type 1) for %s in getAffinityAndComposition.\n",
                        solids[index].label);
                    printf(" Affinity = %6.0f,", bVec[nr]*SCALE);
                    printf(" T (C) = %8.2f, P (bars) = %8.1f, Fmin = %6.0f\n", t-273.15,
                        p, Fmin);
                    printf(" Formula: %s\n", formula);
	  free(formula);
                } else {
                    printf("return FAILURE (type 1) for liquid in getAffinityAndComposition.\n");
                    printf(" Affinity = %6.0f,", bVec[nr]*SCALE);
                    printf(" T (C) = %8.2f, P (bars) = %8.1f, Fmin = %6.0f\n", t-273.15,
                        p, Fmin);
                    for (i=0;i<nlc;i++) printf("%s %f\n", liquid[i].label, bVec[i]);
                }
#endif

#ifdef SPECIAL
                if (bVec[nr]*SCALE < 0.0) {
                    for (i=0; i<nr; i++) indepVar[i] = bVec[i];
                    *affinity = bVec[nr]*SCALE;
                    return SUCCESS;
                } else return FAILURE;
#else
                return FAILURE;
#endif
            }
        }
        for (i=0; i<nr; i++) indepVar[i] = bVec[i];
        *affinity = bVec[nr]*SCALE;
        return SUCCESS;
    } else {
        /* Reassemble solution if equality constraints are present */
        if (hasNull) {
            for (i=0; i<ne; i++) bVec[i] = yVec[i];
            for (i=(ne-1); i>=0; i--) householderRowRow(HOUSEHOLDER_CALC_MODE_H2,
                i, i+1, nr-1, dxdr, i, &hVec[i], &bVec, 0, 0);
            for (i=0; i<nr; i++) if (fabs(bVec[i]) < 10.0*DBL_EPSILON) bVec[i] = 0.0;
        }
#ifdef DEBUG
        if (mode == MODMRT_BAD_INITIAL) {
            printf("Bad initial guess in call to MODMRT.\n");
        } else if (mode == MODMRT_ITERS_EXCEEDED) {
            printf("Function call limit exceeded in call to MODMRT.\n");
        }
        printf("  T (K) = %8.2f, P (bars) = %8.1f, Fmin = %13.6g\n", t, p, Fmin);
        for (i=0; i<na; i++) printf("  bVec[%d] = %f, mu-mu0[%d] = %f\n",
            i, bVec[i], i, muMinusMu0[i]);
#endif
        if (!liquidMode) {
#ifdef DEBUG
            char *formula;
            (*solids[solidID].display)(FIRST, t, p, bVec, &formula);
            printf("return FAILURE (type 2) for %s in getAffinityAndComposition.\n",
                solids[index].label);
            printf(" Affinity = %6.0f,", bVec[nr]*SCALE);
            printf(" T (C) = %8.2f, P (bars) = %8.1f, Fmin = %6.0f\n", t-273.15,
                p, Fmin);
            printf(" Formula: %s\n", formula);
            free(formula);
#endif
            return FAILURE;
#ifdef DEBUG
        } else {
            printf("return FAILURE (type 2) for liquid in getAffinityAndComposition.\n");
            printf(" Affinity = %6.0f,", bVec[nr]*SCALE);
            printf(" T (C) = %8.2f, P (bars) = %8.1f, Fmin = %6.0f\n", t-273.15,
                p, Fmin);
            for (i=0;i<nlc;i++) printf("%s %f\n", liquid[i].label, bVec[i]);
            return FAILURE;
#endif
        }
    }
    return FAILURE;
}

/*****************************************************************************
 * Private function declarations:
 *****************************************************************************/

static void initialGuess(
    double *bVec)       /* returned vector, composition (length nr) + affinity */
{
    /***************************************************************************
   * Consider an ideal solution of n+1 endmembers:
   * - A + RT ln x[1]                  = - mu[1],
   * - A + RT ln x[2]                  = - mu[2],
   *   ...
   * - A + RT ln x[n]                  = - mu[n],
   * - A + RT ln(1-x[1]-x[2]-...-x[n]) = - mu[n+1].
   *
   * Setting: f[1] = exp(-(mu[1]-mu[n])/(R*t)),
   *          f[2] = exp(-(mu[2]-mu[n])/(R*t)),
   *          ...
   *          f[n] = exp(-(mu[n]-mu[n+1])/(R*t))
   *
   * A solution is given by:
   * x[n]   = f[n]/(1+f[n]*(1+f[1]+f[2]+...+f[n-1])
   * x[n-1] = f[n-1]*x[n]
   * x[n-2] = f[n-2]*x[n]
   * ....
   * x[2]   = f[2]*x[n]
   * x[1]   = f[1]*x[n]
   *
   * Any of the original equations may be used to compute A
   ***************************************************************************/
    int i, j;
    double sum;

    /* take appropriate action for the number of non-zero components */
       if (nz == 0)                      bVec[na-1] = 0.0;
    else if (nz == 1) { moleFrac[0] = 1.0; bVec[na-1] = mu[0]/SCALE; }
    else {

        /* Compute the f[n] terms and store them temporarily in moleFrac[] */
        sum = 1.0;
        if (nz > 2) for (i=0; i<(nz-2); i++) {
            moleFrac[i] = exp(-(mu[i]-mu[nz-2])/(R*t)); sum += moleFrac[i];
        }
        moleFrac[nz-2] = exp(-(mu[nz-2]-mu[nz-1])/(R*t));

        /* Solve for the composition variables (mole fractions) */
        moleFrac[nz-2] /= 1.0 + moleFrac[nz-2]*sum;
        moleFrac[nz-1] = 1.0 - moleFrac[nz-2];
        if (nz > 2) for (i=0; i<(nz-2); i++) {
            moleFrac[i] *= moleFrac[nz-2]; moleFrac[nz-1] -= moleFrac[i];
        }

        /* compute the chemical affinity (choice of mu[] is arbitrary) */
        bVec[na-1] = (mu[0] + R*t*log(moleFrac[0]))/SCALE;

    }

    /* Reassemble the mole fraction and chemical potential vectors with zeros
     for the absent endmembers                                              */
    if (hasNull) for (i=na-1, j=nz; i>=0; i--) {
        if(!nullComp[i]) moleFrac[i] = moleFrac[--j];
        else             moleFrac[i] = 0.0;
    }

    /* convert mole fractions of endmembers into independent compos var */
    if (!liquidMode)
        (*solids[solidID].convert)(SECOND, THIRD, t, p, (double *) NULL,
            moleFrac, bVec, (double *) NULL, (double **) NULL, (double ***) NULL,
            (double **) NULL, (double ****) NULL);
    else
        conLiq(SECOND, THIRD, t, p, (double *) NULL, moleFrac, bVec,
            (double *) NULL, (double **) NULL, (double ***) NULL, (double *) NULL);

}

/*****************************************************************************/

static int computeSolnWithinSolvus(            /* Returns SUCCESS or FAILURE */
    double *bVec,           /* i/o: initial guess on input, solution on output */
    double *Fmin)           /* minimum sum-of-squares value from modmrt        */
{
    int i, j, maxIter, mode;
    double lambda, stepSize;
    double small = sqrt(DBL_EPSILON);
    double scale = DBL_EPSILON;

    /* Redetermine the initial guess based upon ideal solution behavior */
    initialGuess(bVec);
    for (i=0; i<nr; i++) refR[i] = bVec[i];
#ifdef DEBUG
    for (i=0; i<na; i++) printf("  Initial guess: bVec[%d] = %f\n", i, bVec[i]);
    if (!liquidMode) {
        char *formula;
        (*solids[solidID].display)(FIRST, t, p, bVec, &formula);
        printf("  Formula: %s\n", formula);
        free(formula);
    } else for (i=0;i<nlc;i++) printf("%s %f\n", liquid[i].label, bVec[i]);
#endif

    /* Obtain the Hessian of G */
    if (!liquidMode)
        (*solids[solidID].gmix)(THIRD, t, p, bVec, (double *) NULL, (double *) NULL,
            d2gdr2, (double ***) NULL);
    else
     gmixLiq(THIRD, t, p, bVec, (double *) NULL, (double *) NULL, d2gdr2);

    /* Project the Hessian into the Null space of the equality constraints */
    if(hasNull) {
        for (i=0; i<ne; i++) householderRowCol(HOUSEHOLDER_CALC_MODE_H2, i, i+1,
            nr-1, dxdr, i, &hVec[i], d2gdr2, 0, nr-1);
        for (i=0; i<ne; i++) householderRowRow(HOUSEHOLDER_CALC_MODE_H2, i, i+1,
            nr-1, dxdr, i, &hVec[i], d2gdr2, ne, nr-1);
        for (i=0; i<(nz-1); i++)
            for (j=0; j<(nz-1); j++) d2gdr2[i][j] = d2gdr2[ne+i][ne+j];
    }

    /* Project the initial guess into the Null space of constraints */
    if (hasNull) for (i=0; i<ne; i++) householderRowRow(HOUSEHOLDER_CALC_MODE_H2,
        i, i+1, nr-1, dxdr, i, &hVec[i], &bVec, 0, 0);

    /* Scale the hessian */
    for (i=0; i<(nr-ne); i++) for (j=0; j<(nr-ne); j++)
        if (scale < fabs(d2gdr2[i][j])) scale = d2gdr2[i][j];
    for (i=0; i<(nr-ne); i++) for (j=0; j<(nr-ne); j++) d2gdr2[i][j] /= scale;

    /* Find the direction of most negative curvature on the Gibbs surface */
    maxIter = 100;
    if(rqmcg(nr-ne, d2gdr2, identity, &bVec[ne], &maxIter, Fmin)
        == RQMCG_SUCCESS) {
        /* Reassemble the solution if equality constraints are present */
        if (hasNull) {
            for (i=0; i<ne; i++) bVec[i] = yVec[i];
            for (i=(ne-1); i>=0; i--) householderRowRow(HOUSEHOLDER_CALC_MODE_H2,
                i, i+1, nr-1, dxdr, i, &hVec[i], &bVec, 0, 0);
            for (i=0; i<nr; i++) if (fabs(bVec[i]) < 10.0*DBL_EPSILON) bVec[i] = 0.0;
        }
        for (i=0; i<nr; i++) searchR[i] = bVec[i];
#ifdef DEBUG
        for (i=0; i<nr; i++) printf("  Search direction: bVec[%d] = %f\n", i,
            searchR[i]);
#endif

        /* find the lowest Gibbs energy along the direction of negative curvature */
        lambda   = 1.0;
        stepSize = 0.10;
        do {
            maxIter = 100;
            mode = min1d(&lambda, &stepSize, small, &maxIter, Fmin, gmix);
            if (mode == MIN1D_BAD_INITIAL) {
                lambda *= 0.9;
                if (fabs(lambda) < DBL_EPSILON) mode = MIN1D_SUCCESS;
            } else mode = MIN1D_SUCCESS;
        } while (mode != MIN1D_SUCCESS);
        for (i=0; i<nr; i++) bVec[i] = finalR[i];
#ifdef DEBUG
        printf("  G at minimum: %f (lambda: %f)\n", *Fmin, lambda);
        for (i=0; i<nr; i++) printf("  Gmin guess: bVec[%d] = %f\n", i, finalR[i]);
        if (!liquidMode) {
            char *formula;
            (*solids[solidID].display)(FIRST, t, p, bVec, &formula);
            printf("  Formula: %s\n", formula);
            free(formula);
        } else for (i=0;i<nlc;i++) printf("%s %f\n", liquid[i].label, bVec[i]);
#endif

        /* armed with the new initial guess, attempt the solution again */
        maxIter = 250;

        /* Project the new guess into the Null space of constraints */
        if (hasNull) for (i=0; i<ne; i++)
            householderRowRow(HOUSEHOLDER_CALC_MODE_H2, i, i+1, nr-1, dxdr, i,
                &hVec[i], &bVec, 0, 0);

        if (modmrt(nz, nz, &bVec[ne], Fmin, small, maxIter, func, grad) !=
                MODMRT_SUCCESS || *Fmin > small) {
            /* Reassemble the solution if equality constraints are present */
            if (hasNull) {
                for (i=0; i<ne; i++) bVec[i] = yVec[i];
                for (i=(ne-1); i>=0; i--) householderRowRow(HOUSEHOLDER_CALC_MODE_H2,
                    i, i+1, nr-1, dxdr, i, &hVec[i], &bVec, 0, 0);
                for (i=0; i<nr; i++) if(fabs(bVec[i]) < 10.0*DBL_EPSILON) bVec[i] = 0.0;
            }
#ifdef DEBUG
            printf("  FAILURE after attempt with refined initial guess.\n");
            printf("  T (K) = %8.2f, P (bars) = %8.1f, Fmin = %13.6g\n", t, p, *Fmin);
            for (i=0, j=0; i<na; i++) printf("  bVec[%d] = %f, mu-mu0[%d] = %f\n",
                i, bVec[i], i, (nullComp[i]) ? 0.0 : mu[j++]);
            if (!liquidMode) {
                char *formula;
                (*solids[solidID].display)(FIRST, t, p, bVec, &formula);
                printf("  Formula: %s\n", formula);
	free(formula);
            } else for (i=0;i<nlc;i++) printf("%s %f\n", liquid[i].label, bVec[i]);
#endif
            /* Try one last time by looking in the opposite direction */
            lambda   = -1.0;
            stepSize = -0.10;
            do {
                maxIter = 100;
                mode = min1d(&lambda, &stepSize, small, &maxIter, Fmin, gmix);
                if (mode == MIN1D_BAD_INITIAL) {
                    lambda *= 0.9;
                    if (fabs(lambda) < DBL_EPSILON) mode = MIN1D_SUCCESS;
                } else mode = MIN1D_SUCCESS;
            } while (mode != MIN1D_SUCCESS);
            for (i=0; i<nr; i++) bVec[i] = finalR[i];
#ifdef DEBUG
            printf("  G at minimum: %f (lambda: %f)\n", *Fmin, lambda);
            for (i=0; i<nr; i++) printf("  Gmin guess: bVec[%d] = %f\n", i,
                finalR[i]);
            if (!liquidMode) {
                char *formula;
                (*solids[solidID].display)(FIRST, t, p, bVec, &formula);
                printf("  Formula: %s\n", formula);
	free(formula);
            } else for (i=0;i<nlc;i++) printf("%s %f\n", liquid[i].label, bVec[i]);
#endif
            /* armed with the new initial guess, attempt the solution again */
            maxIter = 250;

            /* Project the new guess into the Null space of constraints */
            if (hasNull) for (i=0; i<ne; i++)
                householderRowRow(HOUSEHOLDER_CALC_MODE_H2, i, i+1, nr-1, dxdr, i,
                    &hVec[i], &bVec, 0, 0);

            if (modmrt(nz, nz, &bVec[ne], Fmin, small, maxIter, func, grad) !=
                    MODMRT_SUCCESS || *Fmin > small) {
                /* Reassemble the solution if equality constraints are present */
                if (hasNull) {
                    for (i=0; i<ne; i++) bVec[i] = yVec[i];
                    for (i=(ne-1); i>=0; i--) householderRowRow(HOUSEHOLDER_CALC_MODE_H2,
                        i, i+1, nr-1, dxdr, i, &hVec[i], &bVec, 0, 0);
                    for (i=0; i<nr; i++)
                        if (fabs(bVec[i]) < 10.0*DBL_EPSILON) bVec[i] = 0.0;
                }
                return FAILURE;
            }
        }

    } else {
#ifdef DEBUG
        printf("  FAILED to find a negative curavture search direction\n");
#endif
        return FAILURE;
    }

    /* Reassemble the solution if equality constraints are present */
    if (hasNull) {
        for (i=0; i<ne; i++) bVec[i] = yVec[i];
        for (i=(ne-1); i>=0; i--) householderRowRow(HOUSEHOLDER_CALC_MODE_H2,
            i, i+1, nr-1, dxdr, i, &hVec[i], &bVec, 0, 0);
        for (i=0; i<nr; i++) if (fabs(bVec[i]) < 10.0*DBL_EPSILON) bVec[i] = 0.0;
    }
#ifdef DEBUG
    printf("  SUCCESS after attempt with refined initial guess.\n");
    printf("  T (K) = %8.2f, P (bars) = %8.1f, Fmin = %13.6g\n", t, p, *Fmin);
    for (i=0, j=0; i<na; i++) printf("  bVec[%d] = %f, mu-mu0[%d] = %f\n",
        i, bVec[i], i, (nullComp[i]) ? 0.0 : mu[j++]);
    if (!liquidMode) {
        char *formula;
        (*solids[solidID].display)(FIRST, t, p, bVec, &formula);
        printf("  Formula: %s\n", formula);
        free(formula);
    } else for (i=0;i<nlc;i++) printf("%s %f\n", liquid[i].label, bVec[i]);
#endif
    return SUCCESS;
}

/*****************************************************************************/

static double func(      /* returned value, ith residual if notcomp == FALSE */
    int i,         /* 0 <= i < m, index of residual to be computed             */
    int n,         /* number of independent variables                          */
    double *bVec,  /* vector of independent variables, length n                */
    int *notcomp)  /* returned flag, TRUE is current parameters are infeasible */
{
    double result;
    int j;

    /***************************************************************************
   * When modmrt computes a new estimate to the solution, it calls this
   * routine na times. On the first call we compute the activities of the
   * na endmembers in the solid solution after first testing for violation
   * of bound constraints
   ***************************************************************************/
    *notcomp = FALSE;
    if (i == 0) {
        for (j=0; j<(nz-1); j++) tVec[j+ne] = bVec[j];
        /* Reassemble the solution if equality constraints are present */
        if (hasNull) {
            for (j=0; j<ne; j++) tVec[j] = yVec[j];
            for (j=(ne-1); j>=0; j--) householderRowRow(HOUSEHOLDER_CALC_MODE_H2,
                j, j+1, nr-1, dxdr, j, &hVec[j], &tVec, 0, 0);
            for (j=0; j<nr; j++) if (fabs(tVec[j]) < 10.0*DBL_EPSILON) tVec[j] = 0.0;
        }
        if (!liquidMode) {
            if(!(*solids[solidID].test)(FIFTH, t, p, 0, 0,
                        (char **) NULL, (char **) NULL, tVec, (double *) NULL)) {
                *notcomp = TRUE; return 0.0;
            }
        (*solids[solidID].activity)(FIRST, t, p, tVec, activity, (double *) NULL,
            (double **) NULL);
        } else {
            if (!testLiq(FIFTH, t, p, 0, 0,
                    (char **) NULL, (char **) NULL, tVec, (double *) NULL))
                { *notcomp = TRUE; return 0.0; }
            actLiq(FIRST, t, p, tVec, activity, NULL, NULL, NULL);
        }
    }
    result = (- R*t*log(activity[nullList[i]]) + bVec[nz-1]*SCALE - mu[i])/SCALE;
    return result;
}

/*****************************************************************************/

static void grad(
    int i,              /* 0 <= i < n, index of row of Jacobian to be computed */
    int n,              /* number of independent variables                     */
    double *bVec,       /* vector of independent variables, length n           */
    double *X)          /* returned vector, ith row of the Jacobian            */
{
    int j;

    /***************************************************************************
   * When modmrt computes a new estimate to the Jacobian, it calls this
   * routine na times. On the first call we compute the derivatives of all
   * of the activities with respect to the independent compositional
   * variables. The parameter vector is guaranteed to be feasible on all
   * calls and consequently the value of the activity cannot be zero.
   ***************************************************************************/
    if (i == 0) {
        if (!liquidMode)
            (*solids[solidID].activity)(THIRD, t, p, tVec, (double *) NULL,
                (double *) NULL, dadr);
        else
            actLiq(THIRD, t, p, tVec, NULL, NULL, dadr, NULL);

        /* project the gradient into the Null space of the equality constraints */
        if (hasNull) for (j=0; j<ne; j++)
            householderRowRow(HOUSEHOLDER_CALC_MODE_H2, j, j+1, nr-1, dxdr, j,
                &hVec[j], dadr, 0, na-1);
    }

    for (j=0; j<(nz-1); j++)
        X[j] = - R*t*dadr[nullList[i]][j+ne]/(activity[nullList[i]]*SCALE);
    X[nz-1] = 1.0;
}

/*****************************************************************************/

static double gmix(double lambda, int *notcomp)
{
    int i;
    double gmix;

    for (i=0; i<nr; i++) finalR[i] = refR[i] + lambda*searchR[i];
    if (!liquidMode) {
        if (!(*solids[solidID].test)(FIFTH, t, p, 0, 0,
         (char **) NULL, (char **) NULL, finalR, (double *) NULL))
            { *notcomp = TRUE; return 0.0; }
    } else {
        if (!testLiq(FIFTH, t, p, 0, 0,
         (char **) NULL, (char **) NULL, finalR, (double *) NULL))
            { *notcomp = TRUE; return 0.0; }
    }
    *notcomp = FALSE;
    if (!liquidMode)
        (*solids[solidID].gmix)(FIRST, t, p, finalR, &gmix, (double *) NULL,
            (double **) NULL, (double ***) NULL);
    else gmixLiq(FIRST, t, p, finalR, &gmix, (double *) NULL, (double **) NULL);
    return gmix;
}

/* end of file EST_SATURATION_STATE.C */
