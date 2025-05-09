const char *equality_constraints_ver(void) { return "$Id: equality_constraints.c,v 1.5 2007/06/08 17:25:42 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: equality_constraints.c,v $
MELTS Source Code: RCS Revision 1.5  2007/06/08 17:25:42  ghiorso
MELTS Source Code: RCS Added code to allow regression of Ghiorso EOS parameters
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.4  2006/10/20 00:59:22  ghiorso
MELTS Source Code: RCS (1) Made initial modifications for thread safe code.
MELTS Source Code: RCS (2) Added support for XML I/O in batch mode
MELTS Source Code: RCS (3) Added support for Melts-batch listener for eventual integration into VIGMCS
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.3  2006/08/17 20:47:54  ghiorso
MELTS Source Code: RCS Clarified variable initialization issues in routines.  Problems discovered
MELTS Source Code: RCS when compiler optimization is turned on.
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
MELTS Source Code: RCS Revision 1.2  2005/01/08 22:21:02  cvsaccount
MELTS Source Code: RCS
MELTS Source Code: RCS Set tolerance in silmin (before HFTI call) to 10*DBL_EPSILON to insure
MELTS Source Code: RCS catching phase rule violations in simple system crystallization.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2004/01/02 19:21:49  cvsaccount
MELTS Source Code: RCS CTserver University of Chicago
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.4  2003/04/28 20:44:46  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.3  2002/02/25 17:45:46  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2002/02/18 19:37:55  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:58:12  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
 * Revision 3.8  1997/06/21  22:49:58  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.7  1997/05/03  20:23:37  ghiorso
 * *** empty log message ***
 *
 * Revision 3.6  1997/03/27  17:03:40  ghiorso
 * *** empty log message ***
 *
 * Revision 3.5  1996/09/24  20:33:43  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 3.4  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 3.3  1995/11/23  22:37:42  ghiorso
 * Final implementation of subsolidus fO2 buffering.
 *
 * Revision 3.2  1995/11/01  22:40:27  ghiorso
 * Implementation of subsolidus options after Asimow.
 * Additional implementation of nepheline solid solutions.
 *
 * Revision 3.1  1995/08/18  17:29:58  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Function to construct and project the bulk composition equality
**      constraint matrix
**      (file: EQUALITY_CONSTRAINTS.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  September 18, 1991
**              Extracted from file SILMIN_SUPPORT.C
**      V1.0-2  Mark S. Ghiorso  September 28, 1991
**              Corrected error in column counting for a single component
**              solid
**      V1.0-3  Mark S. Ghiorso  October 15, 1991
**              (1) Modified silminState->solidComp component to prepare for
**                  immiscible solid code
**              (2) Implemented code for immiscible solid phases
**      V2.0-1  Mark S. Ghiorso  November 14, 1991
**              Conversion to OSF Motif V1.1.1/X11 Release 4
**      V3.0-1  Mark S. Ghiorso  April 28, 1992
**              (1) Implementation of mu O2 constraints
**      V3.0-2  Mark S. Ghiorso  September 29, 1993
**              Modified call to realloc to catch zero pointer (SPARC port)
**      V4.0-1  Mark S. Ghiorso  May 11, 1994
**              (1) Added test for isentropic constraints
**                               May 17, 1994
**              (2) Added new code for is***ic constraints (coded on the
**                  assumption that isenthalpic, isentropic, and isochoric
**                  are mutually exclusive) (done May 18, 1994)
**      V4.0-2  Mark S. Ghiorso  May 26, 1994
**              (1) Fixed some minor bugs
**              (2) implemented local isenthalpic, isentropic, isochoric
**                  variables that test for initial state (i.e. whether
**                  refEnthalpy, etc. have been set)
**      V4.0-3  Mark S. Ghiorso  June 14, 1994
**              (1) Corrected test for solid solution inclusion, i.e.
**                  mSol != 0.0
**              (2) Corrected T,P entry for multicomponent solids
**      V4.0-4  Mark S. Ghiorso  June 16, 1994
**              (1) Corrected construction of oxygen potential derivatives
**                  for isenthalpic and isochoric cases
**                               June 17, 1994
**              (2) Fixed V4.0-4.1 for reference state T,P dependence
**                               July 2, 1994
**              (3) Revoked 4.0-4.2 (TEMPORARY)
**                               July 4, 1994
**              (4) Back to 4.0-4.2
**      V4.0-5  Mark S. Ghiorso  March 24, 1995
**              (1) Corrected construction of oxygen potential derivatives
**                  for isenthalpic and isochoric cases in reference state T,P
**                  dependence
**      V5.0-1  Paul D. Asimow  April 26, 1995
**              (1) Enabled subsolidus operation
**      V5.1-1  Paul D. Asimow  July 31, 1995
**              (1) Enabled subsolidus fO2 buffered operation
**--
*/

#include <stdlib.h>
#include <stdio.h>

#ifndef BATCH_VERSION
#include <Xm/Xm.h>
#include "interface.h"            /*Specific external declarations          */
#endif

#include "silmin.h"               /*SILMIN structures include file          */
#include "lawson_hanson.h"        /*Lawson and Hanson, i.e. householder     */

#ifdef DEBUG
#undef DEBUG
#endif

#define REALLOC(x, y) (((x) == NULL) ? malloc(y) : realloc((x), (y)))

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

void getEqualityConstraints(int *conRows, int *conCols, double ***cMatrixPt,
    double **hVectorPt, double **dVectorPt, double **yVectorPt)
{
    static int indFe2O3 = -1, indFeO = -1;
    int i, j, k, l, m, n, nl, ns, rows, cols, kMin = -1, indexT = -1, indexP = -1;
#ifdef PHMELTS_ADJUSTMENTS
    int indH2O = -1, indexH2O = -1;
#endif
    double **cMatrix = *cMatrixPt,
            *hVector = *hVectorPt,
            *dVector = *dVectorPt,
            *yVector = *yVectorPt;
    int isenthalpic = (silminState->refEnthalpy != 0.0) && silminState->isenthalpic;
    int isentropic  = (silminState->refEntropy  != 0.0) && silminState->isentropic;
    int isochoric   = (silminState->refVolume   != 0.0) && silminState->isochoric;
    int hasLiquid   = (silminState->liquidMass  != 0.0);
#ifdef PHMELTS_ADJUSTMENTS
    int H2Obuffer = (silminState->aH2O != 0.0) && silminState->H2Obuffer && hasLiquid; /* don't try to buffer H2O w/o liquid */
#endif

#ifdef PHMELTS_ADJUSTMENTS
    if ((isenthalpic || isochoric || isentropic  || H2Obuffer || (silminState->fo2Path != FO2_NONE)) && constraints == NULL) {
#else
    if ((isenthalpic || isochoric || isentropic  || (silminState->fo2Path != FO2_NONE)) && constraints == NULL) {
#endif

        int na = 0;
        constraints = (Constraints *) malloc((unsigned) sizeof(Constraints));
        constraints->lambda = (double *) malloc((size_t) (2*nc+4)*sizeof(double)); /* nc components, + nc fO2 constraints + 1 + V + S + H  (there are at most nc
                                                                                                                                                                    coexisting liquids, each with an fO2 constraint plus - perhaps - and
                        additional fO2 constraint for the solid phase)    */
        constraints->liquidDelta = (double **) malloc((size_t) nlc*sizeof(double *));                       /* nlc coexisting liquids   */
        for (i=0; i<nlc; i++) constraints->liquidDelta[i] = (double *) malloc((size_t) nlc*sizeof(double)); /* each with nlc components */
        constraints->solidDelta = (double **) malloc((size_t) npc*sizeof(double *));                        /* npc solid phases         */
        for (i=0; i<npc; i++) {
            if (solids[i].type == PHASE) na = MAX(solids[i].na, 1);
            constraints->solidDelta[i] = (double *) malloc((size_t) na*sizeof(double));                       /* na coexisting solids     */
        }
        constraints->lambdaO2 = (double *) malloc((size_t)   (nc+1)*sizeof(double));                        /* There are at most nc (liquid) + 1 (solid) fO2 cons  */
        for (i=0; i<nc; i++) {
            if      (bulkSystem[i].type == FEO)   indFeO   = i;
            else if (bulkSystem[i].type == FE2O3) indFe2O3 = i;
        }
    }

#ifdef PHMELTS_ADJUSTMENTS
    for (i=0; i<nc; i++) if (!strcmp(bulkSystem[i].label, "H2O")) indH2O = i;
#endif

    if (silminState->fo2Path != FO2_NONE && (indFe2O3 == -1 || indFeO == -1)) {
     printf("FATAL ERROR getEqualityConstraints() (equality_constraints.c).\n");
     printf("  Cannot determine array index for Fe2O3 or FeO.\n");
     printf("  Oxygen buffering constraints cannot be imposed.\n");
     exit(1);
    }
#ifdef PHMELTS_ADJUSTMENTS
    if (H2Obuffer && indH2O == -1) {
        printf("FATAL ERROR getEqualityConstraints() (equality_constraints.c).\n");
        printf("  Cannot determine array index for H2O component.\n");
        printf("  Water buffering constraints cannot be imposed.\n");
        exit(1);
    }
    if (H2Obuffer && (isenthalpic || isentropic || isochoric)) {
        printf("FATAL ERROR getEqualityConstraints() (equality_constraints.c).\n");
        printf("  Water buffering constraints cannot be imposed on this path.\n");
        exit(1);
    }
    if (H2Obuffer && silminState->incLiquids > 1) {
        /* Know that H2OBuffer works for one liquid so safest to stick to that */
        printf("FATAL ERROR getEqualityConstraints() (equality_constraints.c).\n");
        printf("  Water buffering constraints not enabled for multiple liquids.\n");
        exit(1);
    }
    /*  if ((silminState->fo2Path != FO2_NONE) && silminState->fo2Alt && (silminState->incLiquids > 1)) {*/
        /* Know that alternative fO2 buffer works for one liquid so safest to stick to that */
        /*printf("FATAL ERROR getEqualityConstraints() (equality_constraints.c).\n");
        printf("  Alternative fO2 buffering constraints not enabled for multiple liquids.\n");
        exit(1);
    }*/
#endif

    /* count rows (number of constraints) and columns (number of variables) */
    for (i=0, rows=0; i<nc; i++) if ((silminState->bulkComp)[i] != 0.0) rows++; /* was > */
    for (i=0, cols=(hasLiquid ? rows*silminState->nLiquidCoexist : 0); i<npc; i++) if ((ns = (silminState->nSolidCoexist)[i]) > 0) {
        if (solids[i].na == 1) cols++;
        else for (j=i+1; j<=i+solids[i].na; j++) for (k=0; k<ns; k++) if ((silminState->solidComp)[j][k] != 0.0) cols++;
    }
        /* f O2 constraint takes a column and row from FeO+Fe2O3 combination; extra liquids add a row */
        /* aH2O constraint replaces xH2O constraint in same row and column */
    if ((silminState->fo2Path != FO2_NONE) && (silminState->nLiquidCoexist > 1)) rows += silminState->nLiquidCoexist-1;
        /* add an extra row and column for temperature */
    if (isenthalpic || isentropic) { indexT = cols; cols++; rows++; }
        /* add an extra row and column for pressure */
    if (isochoric) { indexP = cols; cols++; rows++; }

    /* (re)allocate space for the constraint matrix and solution vectors    */
    /* for (i=rows; i<(*conRows); i++) free(cMatrix[i]);        SPECIAL FIX */

    /* SPECIAL FIX */
    for (i=0; i<(*conRows); i++) free(cMatrix[i]);
    free(cMatrix);
    cMatrix = NULL;
    *conRows = 0;
    /* END SPECIAL FIX */

    cMatrix = (double **) REALLOC(cMatrix, (unsigned) rows*sizeof(double *));
    for (i=0; i<(*conRows); i++)    cMatrix[i] = (double *) REALLOC(cMatrix[i], (unsigned) cols*sizeof(double));
    for (i=(*conRows); i<rows; i++) cMatrix[i] = (double *) malloc((unsigned) cols*sizeof(double));
    hVector = (double *) REALLOC(hVector, (unsigned) rows*sizeof(double));
    dVector = (double *) REALLOC(dVector, (unsigned) rows*sizeof(double));
    yVector = (double *) REALLOC(yVector, (unsigned) rows*sizeof(double));

    *conRows   = rows;
    *conCols   = cols;
    *cMatrixPt = cMatrix;
    *hVectorPt = hVector;
    *dVectorPt = dVector;
    *yVectorPt = yVector;

    /* construct cMatrix and dVector                                        */

    for (i=0, k=0; i<nc; i++) {
        if ((silminState->bulkComp)[i] != 0.0) { /* was > */
#ifdef VERBOSE_DEBUG
            printf(" %15.15s %15.15s %20.13g\n", bulkSystem[i].label, "", (silminState->bulkComp)[i]);
#endif
            l = 0;
            /* liquid first */
            if (hasLiquid) {
                for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                    for (j=0; j<nlc; j++) {
                        if ((silminState->liquidComp)[nl][j] != 0.0) {
#ifdef VERBOSE_DEBUG
            printf(" %15.15s %15.15s %20.13g\n", "", liquid[j].label, (silminState->liquidComp)[nl][j]);
#endif
                            cMatrix[k][l] = (liquid[j].liqToOx)[i];
                            l++;
                        }
                    }
            	}
            }
            /* now solids */
            for (j=0; j<npc; j++) {
                if ((ns = (silminState->nSolidCoexist)[j]) > 0) {
                    if (solids[j].na == 1) {
                        cMatrix[k][l] = (solids[j].solToOx)[i];
                        l++;
                    } else {
                        for (n=0; n<ns; n++) for (m=j+1; m<=j+solids[j].na; m++) {
                            if((silminState->solidComp)[m][n] != 0.0) {
                                cMatrix[k][l] = (solids[m].solToOx)[i];
                                l++;
                            }
                        }
                        j += solids[j].na;
                    }
                }
            }
            /* zero additional columns for T and P constraints */
            while (l<cols) cMatrix[k][l++] = 0.0;
            /* total moles for the ith oxide in system */
#ifdef PHMELTS_ADJUSTMENTS
            dVector[k] = (isenthalpic || isentropic  || isochoric   || H2Obuffer || (silminState->fo2Path != FO2_NONE)) ? 0.0 : (silminState->bulkComp)[i];
#else
            dVector[k] = (isenthalpic || isentropic  || isochoric   || (silminState->fo2Path != FO2_NONE)) ? 0.0 : (silminState->bulkComp)[i];
#endif
            /* If the reaction path is mu O2 constrained, combine FeO and Fe2O3 */
            if (silminState->fo2Path != FO2_NONE) {
                if      (i == indFe2O3) for (j=0; j<cols; j++) cMatrix[k][j] *= 2.0;
                if      (i == MIN(indFe2O3, indFeO)) kMin = k;
                else if (i == MAX(indFe2O3, indFeO)) {
                    for (j=0; j<cols; j++) cMatrix[kMin][j] += cMatrix[k][j];
                    k--;
                }
            }
#ifdef PHMELTS_ADJUSTMENTS
            /* record line holding H2O constraint */
            if (i == indH2O) indexH2O = k;
#endif
            k++;
        }
    }

    /* If non-linear constraints exist, add the appropriate gradient row    */
    if (silminState->fo2Path != FO2_NONE) {
        int blockIndex = 0;
        for (nl=0; nl<((hasLiquid) ? silminState->nLiquidCoexist : 1); nl++) {
            if (hasLiquid) {
                double *gradO2 = (double *) malloc((unsigned) nlc*sizeof(double));
                muO2Liq(SECOND, silminState->T, silminState->P, silminState->liquidComp[nl], NULL, gradO2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
                for (i=0; i<cols; i++) cMatrix[k][i] = 0.0;
                for (i=0, j=blockIndex; i<nlc; i++) if ((silminState->liquidComp)[nl][i] != 0.0) cMatrix[k][j++] = gradO2[i];
                blockIndex = j;
                free(gradO2);
            } else {
                double *gradO2 = (double *) malloc((unsigned) cols*sizeof(double));
                subsolidusmuO2(SECOND, NULL, gradO2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
                for (j=0;j<(cols-((isentropic||isenthalpic)?1:0)-(isochoric?1:0));j++) cMatrix[k][j] = gradO2[j];
                for (i=j; i<cols; i++) cMatrix[k][i] = 0.0;
                free(gradO2);
            }
            if (isenthalpic || isentropic) {
                double dmuO2dt;
                if (hasLiquid) muO2Liq(THIRD, silminState->T, silminState->P, silminState->liquidComp[nl], NULL, NULL, &dmuO2dt, NULL, NULL, NULL, NULL,  NULL, NULL, NULL);
                else subsolidusmuO2(THIRD, NULL, NULL, &dmuO2dt, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
                cMatrix[k][indexT] = dmuO2dt - R*log((double) 10.0)*silminState->fo2
                    - R*silminState->T*log((double) 10.0)*getdlog10fo2dt(silminState->T, silminState->P, silminState->fo2Path);
            }
            if (isochoric) {
                double dmuO2dp;
                if (hasLiquid) muO2Liq(FOURTH, silminState->T, silminState->P, silminState->liquidComp[nl], NULL, NULL, NULL, &dmuO2dp, NULL, NULL, NULL, NULL, NULL, NULL);
                else subsolidusmuO2(FOURTH, NULL, NULL, NULL, &dmuO2dp, NULL, NULL, NULL, NULL, NULL, NULL);
                cMatrix[k][indexP] = dmuO2dp - R*silminState->T*log((double) 10.0)*getdlog10fo2dp(silminState->T, silminState->P, silminState->fo2Path);
            }
            dVector[k++] = 0.0;
        }
    }

#ifdef PHMELTS_ADJUSTMENTS
    /* replace xH2O constraint with aH2O constraint if so desired; for now only implement isothermal, isobaric case */
    if (H2Obuffer) {
        double *gradH2O = vector_alloc(nlc);
        for (n=0, j=0; n<silminState->nLiquidCoexist; n++) {
            muH2OLiq(SECOND, silminState->T, silminState->P, silminState->liquidComp[n], NULL, gradH2O, NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL);
            for (i=0; i<nlc; i++) if ((silminState->liquidComp)[n][i] > 0.0) cMatrix[indexH2O][j++] = gradH2O[i];
        }
        for (i=j; i<cols; i++) cMatrix[indexH2O][i] = 0.0;
        vector_free(gradH2O, nlc);
    }
#endif

    if (isenthalpic || isentropic || isochoric) {
        double *rLiq, **drdmLiq, ***d2rdm2Liq, mTotal, pMixLiq, *dpMixLiq, *dpLiq,
            dpMixLiqdT, dpLiqdT = 0.0, dpMixLiqdP, dpLiqdP = 0.0, *mSol, *rSol, **drdmSol,
            ***d2rdm2Sol, pMixSol, *dpMixSol, *dpSol, dpMixSoldT, dpSoldT = 0.0, dpMixSoldP,
            dpSoldP = 0.0;

        if (isenthalpic || isentropic) cMatrix[k][indexT] = 0.0;
        if (isochoric)                 cMatrix[k][indexP] = 0.0;
        l = 0;

        if (hasLiquid) {
            rLiq      = vector_alloc(nlc-1);
            drdmLiq   = matrix_alloc(nlc-1, nlc);
            d2rdm2Liq = (double ***) malloc((unsigned) (nlc-1)*sizeof(double **));
            for (i=0; i<(nlc-1); i++) d2rdm2Liq[i] = matrix_alloc(nlc, nlc);
            dpMixLiq  = vector_alloc(nlc-1);
            dpLiq     = vector_alloc(nlc);

            for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                conLiq(SECOND, THIRD | FIFTH | SIXTH, silminState->T, silminState->P, NULL, silminState->liquidComp[nl], rLiq, NULL, drdmLiq, d2rdm2Liq, NULL);
                for (i=0, mTotal=0.0; i<nlc; i++) mTotal += (silminState->liquidComp)[nl][i];

                if (isenthalpic) {
                    gmixLiq(SECOND, silminState->T, silminState->P, rLiq, NULL, dpMixLiq, NULL);
                    smixLiq(SECOND, silminState->T, silminState->P, rLiq, NULL, dpLiq, NULL, NULL);
                    for (i=0; i<(nlc-1); i++) dpMixLiq[i] += (silminState->T)*dpLiq[i];
                    hmixLiq(FIRST,  silminState->T, silminState->P, rLiq, &pMixLiq, NULL);
                    cpmixLiq(FIRST, silminState->T, silminState->P, rLiq, &dpMixLiqdT, NULL, NULL);
                    intenToExtenGradient(pMixLiq, dpMixLiq, nlc-1, dpLiq, nlc, mTotal, drdmLiq);
                    for (i=0, dpLiqdT=mTotal*dpMixLiqdT; i<nlc; i++) {
                        dpLiq[i] += (liquid[i].cur).h;
                        dpLiqdT  += (silminState->liquidComp)[nl][i]*(liquid[i].cur).cp;
                    }
                } else if (isentropic) {
                    smixLiq(FIRST | SECOND, silminState->T, silminState->P, rLiq, &pMixLiq, dpMixLiq, NULL, NULL);
                    cpmixLiq(FIRST, silminState->T, silminState->P, rLiq, &dpMixLiqdT, NULL, NULL);
                    intenToExtenGradient(pMixLiq, dpMixLiq, nlc-1, dpLiq, nlc, mTotal, drdmLiq);
                    for (i=0, dpLiqdT=mTotal*dpMixLiqdT/silminState->T; i<nlc; i++) {
                        dpLiq[i] += (liquid[i].cur).s;
                        dpLiqdT  += (silminState->liquidComp)[nl][i]*(liquid[i].cur).cp/silminState->T;
                    }
                } else if (isochoric) {
                    vmixLiq(FIRST | SECOND | FIFTH, silminState->T, silminState->P, rLiq, &pMixLiq, dpMixLiq, NULL, NULL, &dpMixLiqdP, NULL, NULL, NULL, NULL, NULL, NULL);
                    intenToExtenGradient(pMixLiq, dpMixLiq, nlc-1, dpLiq, nlc, mTotal, drdmLiq);
                    for (i=0, dpLiqdP=mTotal*dpMixLiqdP; i<nlc; i++) {
                        dpLiq[i] += (liquid[i].cur).v;
                        dpLiqdP  += (silminState->liquidComp)[nl][i]*(liquid[i].cur).dvdp;
                    }
                }

                for (i=0; i<nlc; i++) if ((silminState->liquidComp)[nl][i] != 0.0) cMatrix[k][l++] = dpLiq[i];
                if (isenthalpic || isentropic) cMatrix[k][indexT] = dpLiqdT;
                if (isochoric) cMatrix[k][indexP] = dpLiqdP;
            }

            vector_free(dpMixLiq, nlc-1);
            vector_free(dpLiq, nlc);
            vector_free(rLiq, nlc-1);
            matrix_free(drdmLiq, nlc-1, nlc);
            for (i=0; i<(nlc-1); i++) matrix_free(d2rdm2Liq[i], nlc, nlc);
            free(d2rdm2Liq);
        }

        for (j=0; j<npc; j++) {
            if ((ns = (silminState->nSolidCoexist)[j]) > 0) {
                if (solids[j].na == 1) {
                    mTotal = (silminState->solidComp)[j][0];
                    if (isenthalpic) {
                        cMatrix[k][l++]     = (solids[j].cur).h;
                        cMatrix[k][indexT] += mTotal*(solids[j].cur).cp;
                    } else if (isentropic) {
                        cMatrix[k][l++]     = (solids[j].cur).s;
                        cMatrix[k][indexT] += mTotal*(solids[j].cur).cp/silminState->T;
                    } else if (isochoric) {
                        cMatrix[k][l++]     = (solids[j].cur).v;
                        cMatrix[k][indexP] += mTotal*(solids[j].cur).dvdp;
                    }
                } else {
                    int nr = solids[j].nr;
                    int na = solids[j].na;

                    rSol      = vector_alloc(nr);
                    mSol      = vector_alloc(na);
                    drdmSol   = matrix_alloc(nr, na);
                    d2rdm2Sol = (double ***) malloc((unsigned) nr*sizeof(double **));
                    for (i=0; i<nr; i++) d2rdm2Sol[i] = matrix_alloc(na, na);
                    dpMixSol  = vector_alloc(nr);
                    dpSol     = vector_alloc(na);

                    for (n=0; n<ns; n++) {
                        mTotal = (silminState->solidComp)[j][n];
                        for (i=0; i<na; i++) mSol[i] = (silminState->solidComp)[j+1+i][n];

                        (*solids[j].convert)(SECOND, THIRD | FIFTH | SIXTH, silminState->T, silminState->P, NULL, mSol, rSol, NULL, drdmSol, d2rdm2Sol, NULL, NULL);

                        if (isenthalpic) {
                            (*solids[j].gmix)(SECOND, silminState->T, silminState->P, rSol, NULL, dpMixSol, NULL, NULL);
                            (*solids[j].smix)(SECOND, silminState->T, silminState->P, rSol, NULL, dpSol, NULL);
                            for (i=0; i<nr; i++) dpMixSol[i] += (silminState->T)*dpSol[i];
                            (*solids[j].hmix)(FIRST,  silminState->T, silminState->P, rSol, &pMixSol);
                            (*solids[j].cpmix)(FIRST, silminState->T, silminState->P, rSol, &dpMixSoldT, NULL, NULL);
                            intenToExtenGradient(pMixSol, dpMixSol, nr, dpSol, na, mTotal, drdmSol);
                            for (i=0, dpSoldT=mTotal*dpMixSoldT; i<na; i++) {
                                dpSol[i] += (solids[j+1+i].cur).h;
                                dpSoldT  += mSol[i]*(solids[j+1+i].cur).cp;
                            }
                        } else if (isentropic) {
                            (*solids[j].smix)(FIRST | SECOND, silminState->T, silminState->P, rSol, &pMixSol, dpMixSol, NULL);
                            (*solids[j].cpmix)(FIRST, silminState->T, silminState->P, rSol, &dpMixSoldT, NULL, NULL);
                            intenToExtenGradient(pMixSol, dpMixSol, nr, dpSol, na, mTotal, drdmSol);
                            for (i=0, dpSoldT=mTotal*dpMixSoldT/silminState->T; i<na; i++) {
                                dpSol[i] += (solids[j+1+i].cur).s;
                                dpSoldT  += mSol[i]*(solids[j+1+i].cur).cp/silminState->T;
                            }
                        } else if (isochoric) {
                            (*solids[j].vmix)(FIRST | SECOND | FIFTH, silminState->T, silminState->P, rSol, &pMixSol, dpMixSol, NULL, NULL, &dpMixSoldP, NULL, NULL, NULL, NULL, NULL);
                            intenToExtenGradient(pMixSol, dpMixSol, nr, dpSol, na, mTotal, drdmSol);
                            for (i=0, dpSoldP=mTotal*dpMixSoldP; i<na; i++) {
                                dpSol[i] += (solids[j+1+i].cur).v;
                                dpSoldP  += mSol[i]*(solids[j+1+i].cur).dvdp;
                            }
                        }

                        for (i=0; i<na; i++) if (mSol[i] != 0.0) cMatrix[k][l++] = dpSol[i];
                        if (isenthalpic || isentropic) cMatrix[k][indexT] += dpSoldT;
                        if (isochoric) cMatrix[k][indexP] += dpSoldP;

                    }
                    vector_free(rSol, nr);
                    vector_free(mSol, na);
                    matrix_free(drdmSol, nr, na);
                    for (i=0; i<nr; i++) matrix_free(d2rdm2Sol[i], na, na);
                    free(d2rdm2Sol);
                    vector_free(dpMixSol, nr);
                    vector_free(dpSol, na);
                }
            }
        }
        dVector[k++] = 0.0;
    }

#ifdef DEBUG
    printf("d:c from getEqualityConstraints (rows = %d, cols = %d)\n", rows, cols);
    for (i=0; i<rows; i++) {
        printf(" %9.2g  ", dVector[i]);
        for (j=0; j<cols; j++) printf(" %9.2g", cMatrix[i][j]);
        printf("\n");
    }
#endif

    /* Householder decompose the constraint matrix                          */

    for (i=0; i<rows; i++) householderRowRow(HOUSEHOLDER_CALC_MODE_H1, i, i+1, cols-1, cMatrix, i, &hVector[i], cMatrix, i+1, rows-1);

    /* back-solve the lower-triangular system and store the soln in yVector */

    for (i=0; i<rows; i++) {
        yVector[i] = dVector[i];
        for (j=0; j<i; j++) yVector[i] -= cMatrix[i][j]*yVector[j];
        yVector[i] /= cMatrix[i][i];
    }

#ifdef DEBUG
    printf("y:c' from getEqualityConstraints (rows = %d, cols = %d)\n", rows, cols);
    for (i=0; i<rows; i++) {
        printf(" %9.2g  ", yVector[i]);
        for (j=0; j<cols; j++) printf(" %9.2g", cMatrix[i][j]);
        printf("\n");
    }
#endif

}

/* end of file EQUALITY_CONSTRAINTS.C */
