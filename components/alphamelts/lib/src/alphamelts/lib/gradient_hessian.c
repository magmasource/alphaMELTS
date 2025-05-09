const char *gradient_hessian_ver(void) { return "$Id: gradient_hessian.c,v 1.5 2007/06/08 17:25:42 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: gradient_hessian.c,v $
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
MELTS Source Code: RCS Revision 1.3  2003/05/01 17:33:54  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2003/04/28 20:44:46  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:58:12  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
 * Revision 3.8  1997/06/21  22:49:55  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.7  1997/05/03  20:23:33  ghiorso
 * *** empty log message ***
 *
 * Revision 3.6  1997/03/27  17:03:36  ghiorso
 * *** empty log message ***
 *
 * Revision 3.5  1996/09/24  20:33:40  ghiorso
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
 * Revision 3.1  1995/08/18  17:33:25  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      intenToExtenGradient()
**      intenToExtenHessian()
**
**      These two functions convert partial derivatives of the intensive
**      (molar) thermodynamic potential with respect to a set of independent
**      compositional variables (r), into partial derivatives of the extensive
**      (total) thermodynamic potential with respect to moles of a set of
**      thermodynamic components. pMix, dpMix and d2pMix are the zeroth,
**      first and second derivatives of the intensive potential with respect
**      to r[i], where 0 <= i < nr dp and d2p are the first and second
**      derivatives of the extensive potential with respect to m[i], where
**      0 <= i < na. Note that generally nr+1 == na. mTotal is the total
**      number of moles, i.e. sum of all m[i], and drdm and d2rdm2 are the
**      jacobian and its derivative of r[i] with respect to m[j], and r[i]
**      with respect to m[j] and m[k], respectively.
**
**      getProjGradientAndHessian()
**
**      Function to create the projected gradient and hessian for the system
**      Storage is (re)allocated and contents for eMatrix and bMatrix are
**      set by this routine. The integers conRows and conCols as well as the
**      arrays cMatrix, hVector, dVector and yVector must have been set in a
**      previous call to getEqualityConstraints()
**
**      (file: GRADIENT_HESSIAN.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  September 18, 1991
**              Extracted from file SILMIN_SUPPORT.C
**      V1.0-2  Mark S. Ghiorso  September 24, 1991
**              Altered parameter list to (*solids[].convert)
**      V1.0-3  Mark S. Ghiorso  October 9, 1991
**              Corrected error in gradient construction for pure-component
**              solid phases
**      V1.0-4  Mark S. Ghiorso  October 15, 1991
**              (1) Modified references to silminState->solidComp to reflect
**                  changes to allow for immiscible solid phases
**              (2) Implemented multiple solid phases of the same type
**      V2.0-1  Mark S. Ghiorso  November 14, 1991
**              Conversion to OSF Motif V1.1.1/X11 Release 4
**      V2.0-2  Mark S. Ghiorso  February 13, 1992
**              Corrected reallocatio algorithm for eMatrix and bMatrix to
**              avoid access violation when conCols is less than colRows
**              i.e. when a phase is dropped
**      V3.0-1  Mark S. Ghiorso  April 29, 1992
**              (1) Added support for mu O2 path constraints
**      V3.0-2  Mark S. Ghiorso  May 1, 1992
**              (1) Removed reference to silminState->oxygen
**              (2) Modified algorithm to utilize current - reference oxygen
**                  content when constructing the Korzhinskii potential
**              (3) 3.0-2.2 revoked
**      V3.0-3  Mark S. Ghiorso  September 29, 1993
**              Modified call to realloc to catch zero pointer (SPARC port)
**      V4.0-1  Mark S. Ghiorso  May 11, 1994
**              (1) Added reference to isentropic constraints
**      V4.0-2  Mark S. Ghiorso  May 20, 1994
**              (1) Began modifications for is*ic constraints
**                               May 26, 1994
**              (2) Continue. Also implemented local isenthalpic, isentropic,
**                  isochoric flags to test for initialization of refEnthalpy,
**                  refEntropy, and refVolume
**                               May 27, 1994
**              (3) Finished.
**      V4.0-3  Mark S. Ghiorso  June 17, 1994
**              (1) Corrected error in constructing TT and PP hessian terms
**                  for oxygen constraints
**                               July 2, 1994
**              (2) Revoked 4.0-3.1 (TEMPORARY)
**                               July 4, 1994
**              (3) Back to 4.0-3.1
**      V4.0-4  Mark S. Ghiorso  July 5, 1994
**              (1) Corrected error in oxygen hessians for solid phases
**      V4.0-5  Mark S. Ghiorso  August 3, 1994
**              (1) Changed some debug printouts
**      V4.0-6  Mark S. Ghiorso  March 24, 1995
**              (1) Corrected construction of oxygen potential derivatives
**                  for isenthalpic and isochoric cases in reference state T,P
**                  dependence
**      V5.0-1  Paul D. Asimow  April 26, 1995
**              (1) Enable subsolidus calculation
**      V5.1-1  Paul D. Asimow  August 1, 1995
**              (1) Enable subsolidus fO2 buffering
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

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))

#ifdef DEBUG
#undef DEBUG
#endif

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

void intenToExtenGradient(double pMix, double *dpMix, int nr,  double *dp,
    int na, double mTotal, double **drdm)
{
    int i, j;

    for (j=0; j<na; j++) {
        for (i=0, dp[j] = 0.0; i<nr; i++) dp[j] += dpMix[i]*drdm[i][j];
        dp[j] = mTotal*dp[j] + pMix;
    }
}

void intenToExtenHessian(double pMix, double *dpMix, double **d2pMix,
    int nr, double **d2p, int na, double mTotal, double **drdm,
    double ***d2rdm2)
{
    int i, j, k, l;
    double temp;

    for (j=0; j<na; j++) {
        for (l=0; l<na; l++) {
            for (i=0, d2p[j][l] = 0.0; i<nr; i++) {
                for (k=0, temp=0.0; k<nr; k++) temp += d2pMix[i][k]*drdm[k][l];
                d2p[j][l] += dpMix[i]*(drdm[i][j] + drdm[i][l] + mTotal*d2rdm2[i][j][l])
                   + mTotal*drdm[i][j]*temp;
            }
        }
    }
}

#define REALLOC(x, y) (((x) == NULL) ? malloc(y) : realloc((x), (y)))

int getProjGradientAndHessian(int conRows, int conCols, double ***eMatrixPt,
       double ***bMatrixPt, double **cMatrix, double *hVector, double *dVector,
       double *yVector)
{
    static int colRow = 0;
    double **eMatrix = *eMatrixPt,
         **bMatrix = *bMatrixPt;
    double pMixLiq, *dpMixLiq, **d2pMixLiq, pMixSol, *dpMixSol, **d2pMixSol,
         pMixTmp, *dpMixTmp, **d2pMixTmp, pMixCon, *dpMixCon, d2pMixCon;
    double *dpLiq, **d2pLiq, *dpSol, **d2pSol, *dpTmp, **d2pTmp, pCon = 0.0, *dpCon, d2pCon = 0.0;
    double *rLiq, **drdmLiq, ***d2rdm2Liq;
    double *rSol, *mSol, **drdmSol, ***d2rdm2Sol;
    double *muO2L = NULL, *dmuO2Ldt = NULL, *dmuO2Ldp = NULL, *d2muO2Ldt2 = NULL, *d2muO2Ldp2 = NULL,
         **dmuO2Ldm = NULL, ***d2muO2Ldm2 = NULL, **d2muO2Ldmdt = NULL, **d2muO2Ldmdp = NULL,
    *mO2L = NULL, dmuO2LdpAlt = 0.0;
    double  muO2S = 0.0,  dmuO2Sdt,  dmuO2Sdp,  d2muO2Sdt2,  d2muO2Sdp2,  *dmuO2Sdm = NULL,  **d2muO2Sdm2 = NULL,
         *d2muO2Sdmdt = NULL,  *d2muO2Sdmdp = NULL,  mO2S = 0.0;
    double  muO2E = 0.0,  dmuO2Edt = 0.0,  dmuO2Edp = 0.0,  d2muO2Edt2 = 0.0,  d2muO2Edp2 = 0.0, mO2T;
#ifdef PHMELTS_ADJUSTMENTS
    double muH2O, **dmuH2Odm, mH2O, ***d2muH2Odm2;
    double ***eMatrixfO2 = NULL, **eMatrixT = NULL, **eMatrixP = NULL, **eMatrixH2O = NULL;
    double mTotal;
    int i, j, k, l, m, n, nl, ns, hasNlCon, indexT = -1, indexP = -1, indexH2O = -1, iLiqH2O, iBulkH2O, iPhaseH2O;
#else
    double ***eMatrixfO2 = NULL, **eMatrixT = NULL, **eMatrixP = NULL;
    double mTotal;
    int i, j, k, l, m, n, nl, ns, hasNlCon, indexT = -1, indexP = -1;
#endif
    int isenthalpic = (silminState->refEnthalpy != 0.0) && silminState->isenthalpic;
    int isentropic  = (silminState->refEntropy  != 0.0) && silminState->isentropic;
    int isochoric   = (silminState->refVolume   != 0.0) && silminState->isochoric;
    int hasLiquid   = (silminState->liquidMass  != 0.0);
#ifdef PHMELTS_ADJUSTMENTS
    int H2Obuffer = (silminState->aH2O != 0.0) && silminState->H2Obuffer && hasLiquid; /* don't try to buffer H2O w/o liquid */
#endif

    int numLiquids  = 0;
    int numEMsolids = 0;
    int numSSsolids = 0;

    int hessianType = HESSIAN_TYPE_NORMAL;

    /* (re)allocate storage for eMatrix and bMatrix                             */
    for (i=conCols; i<colRow; i++) { free(eMatrix[i]); free(bMatrix[i]); }
    eMatrix = (double **) REALLOC(eMatrix, (size_t) conCols*sizeof(double *));
    bMatrix = (double **) REALLOC(bMatrix, (size_t) conCols*sizeof(double *));
    for (i=0; i<MIN(colRow, conCols); i++) {
        eMatrix[i] = (double *) REALLOC(eMatrix[i], (size_t) conCols*sizeof(double));
        bMatrix[i] = (double *) REALLOC(bMatrix[i], (size_t) sizeof(double));
    }
    for (i=colRow; i<conCols; i++) {
        eMatrix[i] = (double *) malloc((size_t) conCols*sizeof(double));
        bMatrix[i] = (double *) malloc((size_t) sizeof(double));
    }
    *eMatrixPt = eMatrix;
    *bMatrixPt = bMatrix;

    for (i=0; i<conCols; i++) {
        bMatrix[i][0] = 0.0;
        for (j=0; j<conCols; j++) eMatrix[i][j] = 0.0;
    }

    /****************************************************************************
     Allocate temporary storage for nonlinear constraints
   ***************************************************************************/

#ifdef PHMELTS_ADJUSTMENTS
    hasNlCon = (silminState->fo2Path != FO2_NONE) || H2Obuffer || isenthalpic || isentropic || isochoric;
#else
    hasNlCon = (silminState->fo2Path != FO2_NONE) || isenthalpic || isentropic || isochoric;
#endif
    if (isentropic || isenthalpic) {
        indexT   = conCols - 1;
        eMatrixT = matrix_alloc(conCols, conCols);
        for (i=0; i<conCols; i++) for (j=0; j<conCols; j++) eMatrixT[i][j] = 0.0;
    }
    if (isochoric) {
        indexP   = conCols - 1;
        eMatrixP = matrix_alloc(conCols, conCols);
        for (i=0; i<conCols; i++) for (j=0; j<conCols; j++) eMatrixP[i][j] = 0.0;
    }
#ifdef PHMELTS_ADJUSTMENTS
    if (H2Obuffer) {
        // Make sure only get here for pHMELTS (or Rhyolite-MELTS 1.0.2)
        for (i=0; i<nc; i++) if (!strcmp(bulkSystem[i].label, "H2O") && (silminState->bulkComp[i] != 0.0)) indexH2O = i;
        for (iLiqH2O=0; iLiqH2O<nlc; iLiqH2O++) if (!strcmp(liquid[iLiqH2O].label, "H2O")) break;
        for (iBulkH2O=0; iBulkH2O<nc; iBulkH2O++)if (!strcmp(bulkSystem[iBulkH2O].label, "H2O")) break;
        for (iPhaseH2O=0; iPhaseH2O<npc; iPhaseH2O++) if (!strcmp(solids[iPhaseH2O].label, "fluid")) break;
    }
#endif

    /****************************************************************************
     We now obtain first and second compositional derivatives of the
     thermodynamic potential we seek to minimize
   ****************************************************************************/
    if (hasLiquid) {
        rLiq      = vector_alloc(nlc-1);
        drdmLiq   = matrix_alloc(nlc-1, nlc);
        d2rdm2Liq = (double ***) malloc((size_t) (nlc-1)*sizeof(double **));
        for (i=0; i<(nlc-1); i++) d2rdm2Liq[i] = matrix_alloc(nlc, nlc);

        dpMixLiq  = vector_alloc(nlc-1);
        dpLiq     = vector_alloc(nlc);
        d2pMixLiq = matrix_alloc(nlc-1, nlc-1);
        d2pLiq    = matrix_alloc(nlc, nlc);

        if (isenthalpic || isentropic || isochoric) {
            dpMixTmp  = vector_alloc(nlc-1);
            dpMixCon  = vector_alloc(nlc-1);
            d2pMixTmp = matrix_alloc(nlc-1, nlc-1);
            dpTmp     = vector_alloc(nlc);
            dpCon     = vector_alloc(nlc);
            d2pTmp    = matrix_alloc(nlc, nlc);
        } else {
            dpMixTmp  = NULL;
            dpMixCon  = NULL;
            d2pMixTmp = NULL;
            dpTmp     = NULL;
            dpCon     = NULL;
            d2pTmp    = NULL;
        }

        /***************************************************************************
           nl indexes the loop over all liquids in the assemblage.
           k indexes rows of the Hessian.
           colRow keeps track of the last index for each liquid block in the hessian
        ****************************************************************************/
        for (nl=0, k=0, colRow=0; nl<silminState->nLiquidCoexist; nl++) {
            numLiquids++;
            /* Obtain r (indep var) for liquid, as well as dr/dm and d2r/dm2  	  */
            conLiq(SECOND, THIRD | FIFTH | SIXTH, silminState->T, silminState->P, NULL, silminState->liquidComp[nl], rLiq, NULL, drdmLiq, d2rdm2Liq, NULL);

            /* Obtain the molar potential and its associated derivatives for the liquid */
            if (isenthalpic) {
    	    /**************************************************************************
                *Liq quantities contain negative entropy and it's derivatives
                *Tmp quantities contain enthalpy and it's derivatives
                *Con quantities contain cpMix, dcpMix/dn, dcpMix/dt
            **************************************************************************/
                smixLiq(FIRST | SECOND | THIRD, silminState->T, silminState->P, rLiq, &pMixLiq, dpMixLiq, d2pMixLiq, NULL);
                gmixLiq(FIRST | SECOND | THIRD, silminState->T, silminState->P, rLiq, &pMixTmp, dpMixTmp, d2pMixTmp);
                pMixTmp += silminState->T*pMixLiq;
                for (i=0; i<(nlc-1); i++) {
                    dpMixTmp[i] += silminState->T*dpMixLiq[i];
                    for(j=0; j<(nlc-1); j++) d2pMixTmp[i][j] += silminState->T*d2pMixLiq[i][j];
                }
                pMixLiq *= -1.0;
                for (i=0; i<(nlc-1); i++) {
                    dpMixLiq[i] *= -1.0;
                    for (j=0; j<(nlc-1); j++) d2pMixLiq[i][j] *= -1.0;
                }
                cpmixLiq(FIRST | SECOND | THIRD, silminState->T, silminState->P, rLiq, &pMixCon, &d2pMixCon, dpMixCon);
            } else if (isentropic) {
            /**************************************************************************
                *Liq quantities contain enthalpy and it's derivatives
                *Tmp quantities contain entropy  and it's derivatives
                *Con quantities contain cpMix, dcpMix/dn, dcpMix/dt
            **************************************************************************/
                gmixLiq(FIRST | SECOND | THIRD, silminState->T, silminState->P, rLiq, &pMixLiq, dpMixLiq, d2pMixLiq);
                smixLiq(FIRST | SECOND | THIRD, silminState->T, silminState->P, rLiq, &pMixTmp, dpMixTmp, d2pMixTmp, NULL);
                pMixLiq += silminState->T*pMixTmp;
                for (i=0; i<(nlc-1); i++) {
                    dpMixLiq[i] += silminState->T*dpMixTmp[i];
                    for(j=0; j<(nlc-1); j++) d2pMixLiq[i][j] += silminState->T*d2pMixTmp[i][j];
                }
                cpmixLiq(FIRST | SECOND | THIRD, silminState->T, silminState->P, rLiq, &pMixCon, &d2pMixCon, dpMixCon);
            } else if (isochoric) {
        	/**************************************************************************
                *Liq quantities contain Helmholtz free energy and it's derivatives
                *Tmp quantities contain volume and it's derivatives
                *Con quantities contain dvMix/dp, d2vMix/dpdn, d2vMix/dp2
            **************************************************************************/
                gmixLiq(FIRST | SECOND | THIRD, silminState->T, silminState->P, rLiq, &pMixLiq, dpMixLiq, d2pMixLiq);
                vmixLiq(FIRST | SECOND | THIRD | FIFTH | EIGHTH | TENTH, silminState->T, silminState->P, rLiq, &pMixTmp, dpMixTmp, d2pMixTmp,
            	    NULL, &pMixCon, NULL, NULL, &d2pMixCon, NULL, dpMixCon, NULL);
    	        pMixLiq -= silminState->P*pMixTmp;
                for (i=0; i<(nlc-1); i++) {
                    dpMixLiq[i] -= silminState->P*dpMixTmp[i];
                    for(j=0; j<(nlc-1); j++) d2pMixLiq[i][j] -= silminState->P*d2pMixTmp[i][j];
                }
            } else {
            /**************************************************************************
                *Liq quantities contain Gibbs free energy and it's derivatives
                *Tmp quantities are undefined
                *Con quantities are undefined
            **************************************************************************/
            	gmixLiq(FIRST | SECOND | THIRD, silminState->T, silminState->P, rLiq, &pMixLiq, dpMixLiq, d2pMixLiq);
            }

            /* Obtain the extensive potential gradient				  */
            for (i=0, mTotal=0.0; i<nlc; i++) mTotal += (silminState->liquidComp)[nl][i];
            intenToExtenGradient(pMixLiq, dpMixLiq, nlc-1, dpLiq, nlc, mTotal, drdmLiq);

            if (dpTmp != NULL) intenToExtenGradient(pMixTmp, dpMixTmp, nlc-1, dpTmp, nlc, mTotal, drdmLiq);
            if (dpCon != NULL) intenToExtenGradient(pMixCon, dpMixCon, nlc-1, dpCon, nlc, mTotal, drdmLiq);

            /* Add in the standard state contribution 				  */
            if (isenthalpic) {
                for (i=0, pCon=mTotal*pMixCon, d2pCon=mTotal*d2pMixCon; i<nlc; i++) {
                    dpLiq[i] -= (liquid[i].cur).s;
                    dpTmp[i] += (liquid[i].cur).h;
                    pCon     += (silminState->liquidComp)[nl][i]*(liquid[i].cur).cp;
                    dpCon[i] += (liquid[i].cur).cp;
                    d2pCon   += (silminState->liquidComp)[nl][i]*(liquid[i].cur).dcpdt;
                }
            } else if (isentropic) {
                for (i=0, pCon=mTotal*pMixCon, d2pCon=mTotal*d2pMixCon; i<nlc; i++) {
                    dpLiq[i] += (liquid[i].cur).h;
                    dpTmp[i] += (liquid[i].cur).s;
                    pCon     += (silminState->liquidComp)[nl][i]*(liquid[i].cur).cp;
                    dpCon[i] += (liquid[i].cur).cp;
                    d2pCon   += (silminState->liquidComp)[nl][i]*(liquid[i].cur).dcpdt;
                }
            } else if (isochoric) {
                for (i=0, pCon=mTotal*pMixCon, d2pCon=mTotal*d2pMixCon; i<nlc; i++) {
                    dpLiq[i] += (liquid[i].cur).g - silminState->P*(liquid[i].cur).v;
                    dpTmp[i] += (liquid[i].cur).v;
                    pCon     += (silminState->liquidComp)[nl][i]*(liquid[i].cur).dvdp;
                    dpCon[i] += (liquid[i].cur).dvdp;
                    d2pCon   += (silminState->liquidComp)[nl][i]*(liquid[i].cur).d2vdp2;
                }
            	dmuO2LdpAlt += 2.0*(dpTmp[0]+dpTmp[3]-dpTmp[5]);
            } else {
            	for (i=0; i<nlc; i++) dpLiq[i] += (liquid[i].cur).g;
            }

            /* Obtain the extensive potential Hessian 				  */
            intenToExtenHessian(pMixLiq, dpMixLiq, d2pMixLiq, nlc-1, d2pLiq, nlc, mTotal, drdmLiq, d2rdm2Liq);
            if (d2pTmp != NULL) intenToExtenHessian(pMixTmp, dpMixTmp, d2pMixTmp, nlc-1, d2pTmp, nlc, mTotal, drdmLiq, d2rdm2Liq);

            /****************************************************************************
                Construct and store the gradient and Hessian contribution to the total
                potential for the liquid components.
                eMatrix[][]  contains the Hessian and
                bMatrix[][0] contains the - gradient[] + Hessian * liquid composition or
                only the - gradient[] if non-linear constraints are active
               ***************************************************************************/

            for (i=0; i<nlc; i++) {
            	if ((silminState->liquidComp)[nl][i] != 0.0) {
                    bMatrix[k][0] += -dpLiq[i];
                    for (j=0, l=colRow; j<nlc; j++) {
                        if ((silminState->liquidComp)[nl][j] != 0.0) {
                            eMatrix[k][l] += d2pLiq[i][j];
                            bMatrix[k][0] += (hasNlCon) ? 0.0 : (silminState->liquidComp)[nl][j]*d2pLiq[i][j];
                            if (isenthalpic || isentropic) eMatrixT[k][l] -= d2pTmp[i][j];
                            if (isochoric)		     eMatrixP[k][l] -= d2pTmp[i][j];
                            l++;
                        }
                    }
                    if (isenthalpic) {
                        eMatrix[k][indexT]  += -dpCon[i]/silminState->T;
                        eMatrix[indexT][k]  += -dpCon[i]/silminState->T;
                        eMatrixT[k][indexT] += -dpCon[i];
                        eMatrixT[indexT][k] += -dpCon[i];
                    } else if (isentropic) {
                        eMatrix[k][indexT]  +=  dpCon[i];
                        eMatrix[indexT][k]  +=  dpCon[i];
                        eMatrixT[k][indexT] += -dpCon[i]/silminState->T;
                        eMatrixT[indexT][k] += -dpCon[i]/silminState->T;
                    } else if (isochoric) {
                        eMatrix[k][indexP]  += -silminState->P*dpCon[i];
                        eMatrix[indexP][k]  += -silminState->P*dpCon[i];
                        eMatrixP[k][indexP] += -dpCon[i];
                        eMatrixP[indexP][k] += -dpCon[i];
                    }
                    k++;
                }
            }
            colRow = k; /* save current column/row counter */

            if (isenthalpic) {
                bMatrix[indexT][0]	 += pCon/silminState->T;
                eMatrix[indexT][indexT]  += pCon/SQUARE(silminState->T) - d2pCon/silminState->T;
                eMatrixT[indexT][indexT] += -d2pCon;
            } else if (isentropic) {
                bMatrix[indexT][0]	 += -pCon;
                eMatrix[indexT][indexT]  += d2pCon;
                eMatrixT[indexT][indexT] += pCon/SQUARE(silminState->T) - d2pCon/silminState->T;
            } else if (isochoric) {
                bMatrix[indexP][0]	 += silminState->P*pCon;
                eMatrix[indexP][indexP]  += -pCon - silminState->P*d2pCon;
                eMatrixP[indexP][indexP] += -d2pCon;
            }
        } /* End loop on number of liquids */

        /* Delete storage for the liquid                                            */
        matrix_free(d2pLiq, nlc, nlc);
        vector_free(dpLiq, nlc);
        matrix_free(d2pMixLiq, nlc-1, nlc-1);
        vector_free(dpMixLiq, nlc-1);
        for (i=0; i<(nlc-1); i++) matrix_free(d2rdm2Liq[i], nlc, nlc);
        free(d2rdm2Liq);
        matrix_free(drdmLiq, nlc-1, nlc);
        vector_free(rLiq, nlc-2);
        if (dpMixTmp  != (double *)  NULL) vector_free(dpMixTmp,  nlc-1);
        if (dpMixCon  != (double *)  NULL) vector_free(dpMixCon,  nlc-1);
        if (d2pMixTmp != (double **) NULL) matrix_free(d2pMixTmp, nlc-1, nlc-1);
        if (dpTmp     != (double *)  NULL) vector_free(dpTmp,     nlc);
        if (dpCon     != (double *)  NULL) vector_free(dpCon,     nlc);
        if (d2pTmp    != (double **) NULL) matrix_free(d2pTmp,    nlc, nlc);
    } else colRow = 0;

    /* mu O2 path constraints - revised formulation 2/16/00 */
    if (silminState->fo2Path != FO2_NONE) {
        int muO2Flags = FIRST | SECOND | FIFTH;

        /**************************************************************************
         Compute the total amount of oxygen in the system. Note that the total
            amount of oxygen includes contributions from both the liquid and solid.
        *************************************************************************/
        mO2T = 0.0;
        if (hasLiquid) {
            mO2L = vector_alloc(silminState->nLiquidCoexist);
            for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                for (i=0, mO2L[nl]=0.0; i<nlc; i++) mO2L[nl] += (oxygen.liqToOx)[i]*(silminState->liquidComp)[nl][i];
            	mO2T += mO2L[nl];
            }
        }
        for (i=0, mO2S=0.0; i<npc; i++) for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
            if (solids[i].na == 1)              mO2S += (oxygen.solToOx)[i]*    (silminState->solidComp)[    i][ns];
            else for (j=0; j<solids[i].na; j++) mO2S += (oxygen.solToOx)[i+1+j]*(silminState->solidComp)[i+1+j][ns];
        }
        mO2T += mO2S;

        /**************************************************************************
         Compute the appropriate thermodynamic potential of oxygen, its gradient
            and hessian. Note that the potential and its derivatives are only
            functions of the number of moles of liquid components. Also note that
            the oxygen potential and its T and P derivatives (muO2, dmuO2dt, etc) are
            referenced outside this if block.
        *************************************************************************/
        if (hasLiquid) {
            muO2L      = vector_alloc(silminState->nLiquidCoexist);
            dmuO2Ldm   = matrix_alloc(silminState->nLiquidCoexist, nlc);
            d2muO2Ldm2 = (double ***) malloc((size_t) silminState->nLiquidCoexist*sizeof(double **));
            for (nl=0; nl<silminState->nLiquidCoexist; nl++) d2muO2Ldm2[nl] = matrix_alloc(nlc, nlc);
        } else {
            dmuO2Sdm   = vector_alloc(2*npc);
            d2muO2Sdm2 = matrix_alloc(2*npc, 2*npc);
        }
        if (isenthalpic || isentropic) {
            muO2Flags |= THIRD | SIXTH | EIGHTH;
            if (hasLiquid) {
                d2muO2Ldmdt = matrix_alloc(silminState->nLiquidCoexist, nlc);
                dmuO2Ldt    = vector_alloc(silminState->nLiquidCoexist);
                d2muO2Ldt2  = vector_alloc(silminState->nLiquidCoexist);
            } else {
                d2muO2Sdmdt = vector_alloc(2*npc);
            }
        } else if (isochoric) {
            muO2Flags |= FOURTH | SEVENTH | TENTH;
            if (hasLiquid) {
                d2muO2Ldmdp = matrix_alloc(silminState->nLiquidCoexist, nlc);
                dmuO2Ldp    = vector_alloc(silminState->nLiquidCoexist);
                d2muO2Ldp2  = vector_alloc(silminState->nLiquidCoexist);
            } else {
                d2muO2Sdmdp = vector_alloc(2*npc);
            }
        }
        if (hasLiquid) for (nl=0; nl<silminState->nLiquidCoexist; nl++) muO2Liq(muO2Flags, silminState->T, silminState->P, silminState->liquidComp[nl],
            &muO2L[nl], dmuO2Ldm[nl], (muO2Flags & THIRD) ? &dmuO2Ldt[nl] : NULL, (muO2Flags & FOURTH) ? &dmuO2Ldp[nl] : NULL, d2muO2Ldm2[nl],
            (muO2Flags & SIXTH) ? d2muO2Ldmdt[nl] : NULL, (muO2Flags & SEVENTH) ? d2muO2Ldmdp[nl] : NULL, (muO2Flags & EIGHTH) ? &d2muO2Ldt2[nl] : NULL, NULL,
            (muO2Flags & TENTH) ? &d2muO2Ldp2[nl] : NULL);
        else subsolidusmuO2(muO2Flags, &muO2S, dmuO2Sdm, &dmuO2Sdt, &dmuO2Sdp, d2muO2Sdm2, d2muO2Sdmdt, d2muO2Sdmdp, &d2muO2Sdt2, NULL, &d2muO2Sdp2);

#ifdef DEBUG
        printf("Call to %s from getProjGradientAndHessian with %o\n", (hasLiquid ? "muO2Liq" : "subsolidusmuO2"), muO2Flags);
        if (muO2Flags & FIRST)  {
            if (hasLiquid) {
                printf("  muO2     ="); for (nl=0; nl<silminState->nLiquidCoexist; nl++) printf(" %13.6g", muO2L[nl]);
            	printf(", (E) = %13.6g\n", R*silminState->T*log((double) 10.0)*getlog10fo2(silminState->T, silminState->P, silminState->fo2Path));
            } else printf("  muO2     = %13.6g\n",  muO2S);
        }
        if (muO2Flags & SECOND)  {
            if (hasLiquid) {
                for (i=0; i<nlc; i++) if ((silminState->liquidComp)[0][i] != 0.0) {
                    printf("dmuO2dm[%2.2d]=", i); for (nl=0; nl<silminState->nLiquidCoexist; nl++) printf(" %13.6g", dmuO2Ldm[nl][i]); printf("\n");
            	}
            }
        }
        if (muO2Flags & THIRD)  {
            if (hasLiquid) {
                printf("  dmuO2dt  ="); for (nl=0; nl<silminState->nLiquidCoexist; nl++) printf(" %13.6g", dmuO2Ldt[nl]);   printf("\n");
            } else printf("  dmuO2dt  = %13.6g\n", dmuO2Sdt);
        }
        if (muO2Flags & FOURTH) {
            if (hasLiquid) {
                printf("  dmuO2dp  ="); for (nl=0; nl<silminState->nLiquidCoexist; nl++) printf(" %13.6g", dmuO2Ldp[nl]);
            	printf(", (Alt) = %13.6g\n", dmuO2LdpAlt);
            } else printf("  dmuO2dp  = %13.6g\n", dmuO2Sdp);
        }
        if (muO2Flags & EIGHTH) {
            if (hasLiquid) {
                printf(" d2muO2dt2 ="); for (nl=0; nl<silminState->nLiquidCoexist; nl++) printf(" %13.6g", d2muO2Ldt2[nl]); printf("\n");
            } else printf(" d2muO2dt2 = %13.6g\n", d2muO2Sdt2);
        }
        if (muO2Flags & TENTH)  {
            if (hasLiquid) {
                printf(" d2muO2dp2 ="); for (nl=0; nl<silminState->nLiquidCoexist; nl++) printf(" %13.6g", d2muO2Ldp2[nl]); printf("\n");
            } else printf(" d2muO2dp2 = %13.6g\n", d2muO2Sdp2);
        }
#endif

        if (hasLiquid) {
            for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
            	muO2L[nl] += (oxygen.cur).g;
                if (isenthalpic || isentropic) {
                    dmuO2Ldt[nl]   += -(oxygen.cur).s;
                    d2muO2Ldt2[nl] += -(oxygen.cur).cp/silminState->T;
                }
            	if (isochoric) {
#ifdef PHMELTS_ADJUSTMENTS
                    if (silminState->fo2Alt) {
                        /* no V-fO2 fix */
                        dmuO2Ldp[nl]   += (oxygen.cur).v;
                        d2muO2Ldp2[nl] += (oxygen.cur).dvdp;
                    } else {
#endif
                    dmuO2Ldp[nl]   += 0.0; /* V-fO2 fix: Standard state properties are independent of pressure */
                    d2muO2Ldp2[nl] += 0.0; /* V-fO2 fix: Standard state properties are independent of pressure */
#ifdef PHMELTS_ADJUSTMENTS
                    }
#endif
            	}
            }
        } else {
            muO2S += (oxygen.cur).g;
            if (isenthalpic || isentropic) {
                dmuO2Sdt   += -(oxygen.cur).s;
                d2muO2Sdt2 += -(oxygen.cur).cp/silminState->T;
            }
            if (isochoric) {
                dmuO2Sdp   += (oxygen.cur).v;
                d2muO2Sdp2 += (oxygen.cur).dvdp;
            }
        }

        /* These are the properties of the externally imposed buffer */
        muO2E =  (oxygen.cur).g + R*silminState->T*log((double) 10.0)*getlog10fo2(silminState->T, silminState->P, silminState->fo2Path);
        if (isenthalpic || isentropic) {
            dmuO2Edt   = -(oxygen.cur).s
                 + R*log((double) 10.0)*getlog10fo2(silminState->T, silminState->P, silminState->fo2Path)
                 + R*silminState->T*log((double) 10.0)*getdlog10fo2dt(silminState->T, silminState->P, silminState->fo2Path);
            d2muO2Edt2 = -(oxygen.cur).cp/silminState->T
                 + 2.0*R*log((double) 10.0)*getdlog10fo2dt(silminState->T, silminState->P, silminState->fo2Path)
                 + R*silminState->T*log((double) 10.0)*getd2log10fo2dt2(silminState->T, silminState->P, silminState->fo2Path);
        }
        if (isochoric) {
#ifdef PHMELTS_ADJUSTMENTS
            if (silminState->fo2Alt) {
                /*  no V-fO2 fix */
                dmuO2Edp   = (oxygen.cur).v
                    + R*silminState->T*log((double) 10.0)*getdlog10fo2dp(silminState->T, silminState->P, silminState->fo2Path);
                d2muO2Edp2 = (oxygen.cur).dvdp
                    + R*silminState->T*log((double) 10.0)*getd2log10fo2dp2(silminState->T, silminState->P, silminState->fo2Path);
            } else {
#endif
            /* V-fO2 fix: Standard state properties are independent of pressure */
            dmuO2Edp   = R*silminState->T*log((double) 10.0)*getdlog10fo2dp(silminState->T, silminState->P, silminState->fo2Path);
            d2muO2Edp2 = R*silminState->T*log((double) 10.0)*getd2log10fo2dp2(silminState->T, silminState->P, silminState->fo2Path);
#ifdef PHMELTS_ADJUSTMENTS
            }
#endif
        }

#ifdef DEBUG
        printf("After Correction:\n");
        if (muO2Flags & FIRST)  {
            if (hasLiquid) {
                printf("  muO2     ="); for (nl=0; nl<silminState->nLiquidCoexist; nl++) printf(" %13.6g", muO2L[nl]);      printf(", (E) = %13.6g\n", muO2E);
            } else printf("  muO2     = %13.6g\n",  muO2S);
        }
        if (muO2Flags & THIRD)  {
            if (hasLiquid) {
                printf("  dmuO2dt  ="); for (nl=0; nl<silminState->nLiquidCoexist; nl++) printf(" %13.6g", dmuO2Ldt[nl]);   printf(", (E) = %13.6g\n", dmuO2Edt);
            } else printf("  dmuO2dt  = %13.6g\n", dmuO2Sdt);
        }
        if (muO2Flags & FOURTH) {
            if (hasLiquid) {
                printf("  dmuO2dp  ="); for (nl=0; nl<silminState->nLiquidCoexist; nl++) printf(" %13.6g", dmuO2Ldp[nl]);   printf(", (E) = %13.6g\n", dmuO2Edp);
            } else printf("  dmuO2dp  = %13.6g\n", dmuO2Sdp);
        }
        if (muO2Flags & EIGHTH) {
            if (hasLiquid) {
                printf(" d2muO2dt2 ="); for (nl=0; nl<silminState->nLiquidCoexist; nl++) printf(" %13.6g", d2muO2Ldt2[nl]); printf(", (E) = %13.6g\n", d2muO2Edt2);
            } else printf(" d2muO2dt2 = %13.6g\n", d2muO2Sdt2);
        }
        if (muO2Flags & TENTH)  {
            if (hasLiquid) {
                printf(" d2muO2dp2 ="); for (nl=0; nl<silminState->nLiquidCoexist; nl++) printf(" %13.6g", d2muO2Ldp2[nl]); printf(", (E) = %13.6g\n", d2muO2Edp2);
            } else printf(" d2muO2dp2 = %13.6g\n", d2muO2Sdp2);
        }
#endif

        eMatrixfO2 = (double ***) malloc((size_t) ((hasLiquid) ? silminState->nLiquidCoexist : 1)*sizeof(double **));
        for (nl=0; nl<((hasLiquid) ? silminState->nLiquidCoexist : 1); nl++) {
            eMatrixfO2[nl] = matrix_alloc(conCols, conCols);
            for (i=0; i<conCols; i++) for (j=0; j<conCols; j++) eMatrixfO2[nl][i][j] = 0.0;
        }

        if (hasLiquid) {
            int topOfBlock;
            /**************************************************************************
             Correct the Liquid gradient vector (bMatrix) and liquid Hessian (eMatrix)
            for the oxygen contribution and store extra oxygen Hessian which will be
            used after the Lagrange multipliers are estimated.
            *************************************************************************/
            for (nl=0, k=0, topOfBlock=0; nl<silminState->nLiquidCoexist; nl++) {
                for (i=0; i<nlc; i++) if ((silminState->liquidComp)[nl][i] != 0.0) {
                    if (isenthalpic) {
                        bMatrix[k][0]	      += (oxygen.liqToOx)[i]*muO2E/silminState->T;
                        eMatrix[indexT][k]        += (oxygen.liqToOx)[i]*muO2E/SQUARE(silminState->T) - (oxygen.liqToOx)[i]*dmuO2Edt/silminState->T;
                        eMatrix[k][indexT]        += (oxygen.liqToOx)[i]*muO2E/SQUARE(silminState->T) - (oxygen.liqToOx)[i]*dmuO2Edt/silminState->T;
                        eMatrixfO2[nl][indexT][k] += -d2muO2Ldmdt[nl][i];
                        eMatrixfO2[nl][k][indexT] += -d2muO2Ldmdt[nl][i];
                    } else if (isentropic) {
                        bMatrix[k][0]	      +=  (oxygen.liqToOx)[i]*muO2E;
                        eMatrix[indexT][k]        += -(oxygen.liqToOx)[i]*dmuO2Edt;
                        eMatrix[k][indexT]        += -(oxygen.liqToOx)[i]*dmuO2Edt;
                        eMatrixfO2[nl][indexT][k] += -d2muO2Ldmdt[nl][i];
                        eMatrixfO2[nl][k][indexT] += -d2muO2Ldmdt[nl][i];
                    } else if (isochoric) {
#ifdef PHMELTS_ADJUSTMENTS
                        if (silminState->fo2Alt) {
                            bMatrix[k][0]           +=  mO2T*dmuO2Ldm[nl][i]    + (oxygen.liqToOx)[i]*muO2E;
                            eMatrix[indexP][k]        += -mO2T*d2muO2Ldmdp[nl][i] - (oxygen.liqToOx)[i]*dmuO2Edp;
                            eMatrix[k][indexP]        += -mO2T*d2muO2Ldmdp[nl][i] - (oxygen.liqToOx)[i]*dmuO2Edp;
                            eMatrixfO2[nl][indexP][k] += -d2muO2Ldmdp[nl][i];
                            eMatrixfO2[nl][k][indexP] += -d2muO2Ldmdp[nl][i];
                        } else {
#endif
                        bMatrix[k][0]	      +=  mO2T*dmuO2Ldm[nl][i]    + (oxygen.liqToOx)[i]*muO2L[nl];    /* V-fO2 fix was:  (oxygen.liqToOx)[i]*muO2E;    */
                        eMatrix[indexP][k]        += -mO2T*d2muO2Ldmdp[nl][i] - (oxygen.liqToOx)[i]*dmuO2Ldp[nl]; /* V-fO2 fix was: -(oxygen.liqToOx)[i]*dmuO2Edp; */
                        eMatrix[k][indexP]        += -mO2T*d2muO2Ldmdp[nl][i] - (oxygen.liqToOx)[i]*dmuO2Ldp[nl]; /* V-fO2 fix was: -(oxygen.liqToOx)[i]*dmuO2Edp; */
                        eMatrixfO2[nl][indexP][k] += -d2muO2Ldmdp[nl][i];                                         /* V-fO2 fix was: -d2muO2Ldmdp[nl][i];  	       */
                        eMatrixfO2[nl][k][indexP] += -d2muO2Ldmdp[nl][i];                                         /* V-fO2 fix was: -d2muO2Ldmdp[nl][i];  	       */
#ifdef PHMELTS_ADJUSTMENTS
                        }
#endif
                    } else {
#ifdef PHMELTS_ADJUSTMENTS
                        if (silminState->fo2Alt) {
                            bMatrix[k][0]       +=  mO2T*dmuO2Ldm[nl][i]    + (oxygen.liqToOx)[i]*muO2E;
                        } else {
#endif
                        bMatrix[k][0]       +=  mO2T*dmuO2Ldm[nl][i]    + (oxygen.liqToOx)[i]*muO2L[nl];    /* V-fO2 fix was:  (oxygen.liqToOx)[i]*muO2E;    */
#ifdef PHMELTS_ADJUSTMENTS
                        }
#endif
                    }

#ifdef PHMELTS_ADJUSTMENTS
                    if (silminState->fo2Alt) {
                        for (j=0, l=topOfBlock; j<nlc; j++) if ((silminState->liquidComp)[nl][j] != 0.0) {
                                eMatrixfO2[nl][k][l] += -d2muO2Ldm2[nl][i][j];
                                l++;
                        }
                    } else {
#endif
                    for (j=0, l=topOfBlock; j<nlc; j++) if ((silminState->liquidComp)[nl][j] != 0.0) {
                	    eMatrix[k][l]        += -(oxygen.liqToOx)[i]*dmuO2Ldm[nl][j] - (oxygen.liqToOx)[j]*dmuO2Ldm[nl][i] - mO2T*d2muO2Ldm2[nl][i][j]; /* V-fO2 fix new */
                        eMatrixfO2[nl][k][l] += -d2muO2Ldm2[nl][i][j];
                        l++;
                    }

                    /* These entries in E pertain to solids and will be multiplied below by appropriate stoichiometric coefficients  */
                    for (j=colRow; j<((isenthalpic || isentropic || isochoric) ? conCols-1 : conCols); j++) { /* V-fO2 fix new */
                        eMatrix[k][j] += -dmuO2Ldm[nl][i];                                                      /* V-fO2 fix new */
                        eMatrix[j][k] += -dmuO2Ldm[nl][i];                                                      /* V-fO2 fix new */
                    }                                                                                         /* V-fO2 fix new */
#ifdef PHMELTS_ADJUSTMENTS
                    }
#endif

                    k++;
                }
                topOfBlock = k;
            } /* end loop on number of liquids */
        } /* end hasLiquid block */

        /* Note that the terms in this block pertain to BOTH the liquid and solid phases */
        if (isenthalpic) {
            if (hasLiquid) {
                for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                    bMatrix[indexT][0]	         += -mO2T*muO2E/SQUARE(silminState->T) + mO2T*dmuO2Edt/silminState->T;
                    eMatrix[indexT][indexT]        += -2.0*mO2T*muO2E/CUBE(silminState->T) + 2.0*mO2T*dmuO2Edt/SQUARE(silminState->T) - mO2T*d2muO2Edt2/silminState->T;
                    eMatrixfO2[nl][indexT][indexT] += -d2muO2Ldt2[nl] + d2muO2Edt2;
            	}
            } else {
                bMatrix[indexT][0]	         += -mO2S*muO2S/SQUARE(silminState->T) + mO2S*dmuO2Sdt/silminState->T;
                eMatrix[indexT][indexT]          += -2.0*mO2S*muO2S/CUBE(silminState->T) + 2.0*mO2S*dmuO2Sdt/SQUARE(silminState->T) - mO2S*d2muO2Sdt2/silminState->T;
                eMatrixfO2[0][indexT][indexT]    += -d2muO2Sdt2 + d2muO2Edt2;
            }
        } else if (isentropic) {
            if (hasLiquid) {
                for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                    bMatrix[indexT][0]	         +=  mO2T*dmuO2Edt;
                    eMatrix[indexT][indexT]        += -mO2T*d2muO2Edt2;
                    eMatrixfO2[nl][indexT][indexT] += -d2muO2Ldt2[nl] + d2muO2Edt2;
            	}
            } else {
                bMatrix[indexT][0]	         +=  mO2S*dmuO2Sdt;
                eMatrix[indexT][indexT]          += -mO2S*d2muO2Sdt2;
                eMatrixfO2[0][indexT][indexT]    += -d2muO2Sdt2 + d2muO2Edt2;
            }
        } else if (isochoric) {
            if (hasLiquid) {
                for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
#ifdef PHMELTS_ADJUSTMENTS
                    if (silminState->fo2Alt) {
                        bMatrix[indexP][0]           +=  mO2T*dmuO2Edp;
                        eMatrix[indexP][indexP]      += -mO2T*d2muO2Edp2;
                    } else {
#endif
                    bMatrix[indexP][0]	         +=  mO2T*dmuO2Ldp[nl];            /* V-fO2 fix was:  mO2T*dmuO2Edp;   */
                    eMatrix[indexP][indexP]        += -mO2T*d2muO2Ldp2[nl];          /* V-fO2 fix was: -mO2T*d2muO2Edp2; */
#ifdef PHMELTS_ADJUSTMENTS
                    }
#endif
    	  eMatrixfO2[nl][indexP][indexP] += -d2muO2Ldp2[nl] + d2muO2Edp2;
                }
            } else {
                bMatrix[indexP][0]	         +=  mO2S*dmuO2Sdp;
                eMatrix[indexP][indexP]          += -mO2S*d2muO2Sdp2;
                eMatrixfO2[0][indexP][indexP]    += -d2muO2Sdp2 + d2muO2Edp2;
            }
        }
    }

#ifdef PHMELTS_ADJUSTMENTS
    /* mu H2O path constraints                                                  */
    if (H2Obuffer) {
        int muH2OFlags = FIRST | SECOND | FIFTH;
        int colRowH2O = 0;
        /**************************************************************************
         Compute the total amount of water in the system. Note that the total
            amount of water includes contributions from both the liquid and solid.
        *************************************************************************/
        for (nl=0, mH2O=0.0; nl<(silminState->nLiquidCoexist); nl++) mH2O += (silminState->liquidComp)[nl][iLiqH2O];
        for (i=0; i<npc; i++)
            for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
                if (solids[i].na == 1) mH2O += (solids[i].solToOx)[iBulkH2O]*(silminState->solidComp)[i][ns];
                else for (j=0; j<solids[i].na; j++) mH2O += (solids[i+1+j].solToOx)[iBulkH2O]*(silminState->solidComp)[i+1+j][ns];
        }
        /**************************************************************************
         Compute the appropriate thermodynamic potential of water, its gradient
            and hessian.
        *************************************************************************/
        dmuH2Odm   = matrix_alloc(silminState->nLiquidCoexist, nlc);
        d2muH2Odm2 = (double ***) malloc((size_t) silminState->nLiquidCoexist*sizeof(double **));
        for (nl=0; nl<silminState->nLiquidCoexist; nl++) d2muH2Odm2[nl] = matrix_alloc(nlc, nlc);
        eMatrixH2O = matrix_alloc(conCols, conCols);
        for (i=0; i<conCols; i++) for (j=0; j<conCols; j++) eMatrixH2O[i][j] = 0.0;

        for (nl=0, k=0; nl< silminState->nLiquidCoexist; nl++) {
            muH2OLiq(muH2OFlags, silminState->T, silminState->P, silminState->liquidComp[nl], &muH2O, dmuH2Odm[nl], NULL, NULL, d2muH2Odm2[nl], NULL, NULL, NULL, NULL, NULL);
            muH2O += (liquid[iLiqH2O].cur).g;

            /**************************************************************************
             Correct the Liquid gradient vector (bMatrix) and liquid Hessian (eMatrix)
                for the water contribution and store extra water Hessian which will be
                used after the Lagrange multipliers are estimated.
            *************************************************************************/
            for (i=0; i<nlc; i++) if ((silminState->liquidComp)[nl][i] > 0.0) {
                bMatrix[k][0] += (i==iLiqH2O ? muH2O : 0) + mH2O*dmuH2Odm[nl][i];

                for (j=0, l=0; j<nlc; j++) if ((silminState->liquidComp)[nl][j] > 0.0) {
                    eMatrix[k][colRowH2O+l] += -(i==iLiqH2O ? 1 : 0)*dmuH2Odm[nl][j] - (j==iLiqH2O?1:0)*dmuH2Odm[nl][i] - mH2O*d2muH2Odm2[nl][i][j];
                    eMatrixH2O[k][colRowH2O+l] += -d2muH2Odm2[nl][i][j];
                    l++;
                }
                k++;
            }
            colRowH2O=k;
        }
    }
#endif

    /* Loop over all solids in the system, adding those as required             */
    for (i=0; i<npc; i++) {
        if ((ns = (silminState->nSolidCoexist)[i]) > 0) {
            if (solids[i].na == 1) {
                numEMsolids++;
                /* add a simple one-component solid to the gradient and Hessian       */
                if (isenthalpic) {
                    mTotal                    = (silminState->solidComp)[i][0];
                    bMatrix[colRow][0]       += (solids[i].cur).s;
                    bMatrix[indexT][0]       += mTotal*(solids[i].cur).cp/silminState->T;
                    eMatrix[colRow][indexT]  += -(solids[i].cur).cp/silminState->T;
                    eMatrix[indexT][colRow]  += -(solids[i].cur).cp/silminState->T;
                    eMatrixT[colRow][indexT] += -(solids[i].cur).cp;
                    eMatrixT[indexT][colRow] += -(solids[i].cur).cp;
                    eMatrix[indexT][indexT]  += mTotal*( - (solids[i].cur).dcpdt/silminState->T + (solids[i].cur).cp/SQUARE(silminState->T) );
                    eMatrixT[indexT][indexT] += -mTotal*(solids[i].cur).dcpdt;
                } else if (isentropic) {
                    mTotal                    = (silminState->solidComp)[i][0];
                    bMatrix[colRow][0]       += -(solids[i].cur).h;
                    bMatrix[indexT][0]       += -mTotal*(solids[i].cur).cp;
                    eMatrix[colRow][indexT]  += (solids[i].cur).cp;
                    eMatrix[indexT][colRow]  += (solids[i].cur).cp;
                    eMatrixT[colRow][indexT] += -(solids[i].cur).cp/silminState->T;
                    eMatrixT[indexT][colRow] += -(solids[i].cur).cp/silminState->T;
                    eMatrix[indexT][indexT]  += mTotal*(solids[i].cur).dcpdt;
                    eMatrixT[indexT][indexT] += mTotal*( - (solids[i].cur).dcpdt/silminState->T + (solids[i].cur).cp/SQUARE(silminState->T) );
                } else if (isochoric) {
                    mTotal                    = (silminState->solidComp)[i][0];
                    bMatrix[colRow][0]       += - ((solids[i].cur).g - silminState->P*(solids[i].cur).v);
                    bMatrix[indexP][0]       += mTotal*silminState->P*(solids[i].cur).dvdp;
                    eMatrix[colRow][indexP]  += -silminState->P*(solids[i].cur).dvdp;
                    eMatrix[indexP][colRow]  += -silminState->P*(solids[i].cur).dvdp;
                    eMatrixP[colRow][indexP] += -(solids[i].cur).dvdp;
                    eMatrixP[indexP][colRow] += -(solids[i].cur).dvdp;
                    eMatrix[indexP][indexP]  += -mTotal*((solids[i].cur).dvdp + silminState->P*(solids[i].cur).d2vdp2);
                    eMatrixP[indexP][indexP] += -mTotal*(solids[i].cur).d2vdp2;
                } else bMatrix[colRow][0]  += - (solids[i].cur).g;

                /**********************************************************************************************************************/
                /* Solid contribution: mu O2 path constraints:                                                                        */
                /*   If a liquid is present, the mu f O2 is imposed on the solid and is given by the "external" properties stored in  */
                /*     muO2E, dmuO2Edt, dmuO2Edp, d2muO2Edt2, d2muO2Edp2.  There are no solid-composition derivatives.                */
                /*   If a liquid is not present, the mu f O2 is a function of solid composition and appropriate derivatives are given */
                /*     in subSolidusfO2(...).                                                                                         */
                /**********************************************************************************************************************/
                if (silminState->fo2Path != FO2_NONE) {
                    if (isenthalpic)       {
                        bMatrix[colRow][0]      += hasLiquid ? (oxygen.solToOx)[i]*muO2E/silminState->T
                                                                                    : (oxygen.solToOx)[i]*muO2S/silminState->T + mO2S*dmuO2Sdm[colRow]/silminState->T;
                        if (!hasLiquid) {
                            for (k=0; k<conCols-1; k++) {
                                eMatrix[colRow][k] += -(oxygen.solToOx)[i]*dmuO2Sdm[k]/silminState->T;
                                eMatrix[k][colRow] += -(oxygen.solToOx)[i]*dmuO2Sdm[k]/silminState->T;
                            }
                            eMatrixfO2[0][indexT][colRow] += -d2muO2Sdmdt[colRow];
                            eMatrixfO2[0][colRow][indexT] += -d2muO2Sdmdt[colRow];
                        }
                        eMatrix[colRow][indexT] += hasLiquid ? (oxygen.solToOx)[i]*(muO2E/SQUARE(silminState->T) - dmuO2Edt/silminState->T)
                                                 : (oxygen.solToOx)[i]*(muO2S/SQUARE(silminState->T) - dmuO2Sdt/silminState->T)
						 + mO2S*(dmuO2Sdm[colRow]/SQUARE(silminState->T) - d2muO2Sdmdt[colRow]/silminState->T);
                        eMatrix[indexT][colRow] += hasLiquid ? (oxygen.solToOx)[i]*(muO2E/SQUARE(silminState->T) - dmuO2Edt/silminState->T)
                                                                                    : (oxygen.solToOx)[i]*(muO2S/SQUARE(silminState->T) - dmuO2Sdt/silminState->T)
						 + mO2S*(dmuO2Sdm[colRow]/SQUARE(silminState->T) - d2muO2Sdmdt[colRow]/silminState->T);
                    } else if (isentropic) {
                        bMatrix[colRow][0]      += hasLiquid ? (oxygen.solToOx)[i]*muO2E : (oxygen.solToOx)[i]*muO2S + mO2S*dmuO2Sdm[colRow];
                        if (!hasLiquid) {
                            for (k=0; k<conCols-1; k++) {
                                eMatrix[colRow][k] += -(oxygen.solToOx)[i]*dmuO2Sdm[k];
                                eMatrix[k][colRow] += -(oxygen.solToOx)[i]*dmuO2Sdm[k];
                            }
                            eMatrixfO2[0][indexT][colRow] += -d2muO2Sdmdt[colRow];
                            eMatrixfO2[0][colRow][indexT] += -d2muO2Sdmdt[colRow];
                        }
                        eMatrix[colRow][indexT] += hasLiquid ? -(oxygen.solToOx)[i]*dmuO2Edt : -(oxygen.solToOx)[i]*dmuO2Sdt - mO2S*d2muO2Sdmdt[colRow];
                        eMatrix[indexT][colRow] += hasLiquid ? -(oxygen.solToOx)[i]*dmuO2Edt : -(oxygen.solToOx)[i]*dmuO2Sdt - mO2S*d2muO2Sdmdt[colRow];
                    } else if (isochoric)  {
#ifdef PHMELTS_ADJUSTMENTS
                        if (silminState->fo2Alt) {
                            bMatrix[colRow][0]      += hasLiquid ? (oxygen.solToOx)[i]*muO2E
                                                        : (oxygen.solToOx)[i]*muO2S + mO2S*dmuO2Sdm[colRow];
                        } else {
#endif
                        bMatrix[colRow][0]      += hasLiquid ? (oxygen.solToOx)[i]*muO2L[0]                       /* V-fO2 fix was: (oxygen.solToOx)[i]*muO2E */
                                                    : (oxygen.solToOx)[i]*muO2S + mO2S*dmuO2Sdm[colRow];
#ifdef PHMELTS_ADJUSTMENTS
                        }
#endif
                        if (!hasLiquid) {
                            for (k=0; k<conCols-1; k++) {
                                eMatrix[colRow][k] += -(oxygen.solToOx)[i]*dmuO2Sdm[k];
                                eMatrix[k][colRow] += -(oxygen.solToOx)[i]*dmuO2Sdm[k];
                            }
                            eMatrixfO2[0][indexP][colRow] += -d2muO2Sdmdp[colRow];
                            eMatrixfO2[0][colRow][indexP] += -d2muO2Sdmdp[colRow];
                        }
#ifdef PHMELTS_ADJUSTMENTS
                    if (silminState->fo2Alt) {
                        eMatrix[colRow][indexP] += hasLiquid ? -(oxygen.solToOx)[i]*dmuO2Edp
                                                 : -(oxygen.solToOx)[i]*dmuO2Sdp - mO2S*d2muO2Sdmdp[colRow];
                        eMatrix[indexP][colRow] += hasLiquid ? -(oxygen.solToOx)[i]*dmuO2Edp
                                                 : -(oxygen.solToOx)[i]*dmuO2Sdp - mO2S*d2muO2Sdmdp[colRow];
                    } else {
#endif
                        eMatrix[colRow][indexP] += hasLiquid ? -(oxygen.solToOx)[i]*dmuO2Ldp[0]                          /* V-fO2 fix was: -(oxygen.solToOx)[i]*dmuO2Edp */
                                                        : -(oxygen.solToOx)[i]*dmuO2Sdp - mO2S*d2muO2Sdmdp[colRow];
                        eMatrix[indexP][colRow] += hasLiquid ? -(oxygen.solToOx)[i]*dmuO2Ldp[0]                          /* V-fO2 fix was: -(oxygen.solToOx)[i]*dmuO2Edp */
                                                        : -(oxygen.solToOx)[i]*dmuO2Sdp - mO2S*d2muO2Sdmdp[colRow];
#ifdef PHMELTS_ADJUSTMENTS
                    }
#endif
                    } else {
#ifdef PHMELTS_ADJUSTMENTS
                        if (silminState->fo2Alt) {
                            bMatrix[colRow][0] += hasLiquid ? (oxygen.solToOx)[i]*muO2E
                                                        : (oxygen.solToOx)[i]*muO2S + mO2S*dmuO2Sdm[colRow];
                        } else {
#endif
                        bMatrix[colRow][0] += hasLiquid ? (oxygen.solToOx)[i]*muO2L[0]                                   /* V-fO2 fix was: (oxygen.solToOx)[i]*muO2E */
	                                    : (oxygen.solToOx)[i]*muO2S + mO2S*dmuO2Sdm[colRow];
#ifdef PHMELTS_ADJUSTMENTS
                        }
#endif
                        if (!hasLiquid) {
                            for (k=0; k<conCols-1; k++) {
                                eMatrix[colRow][k] += -(oxygen.solToOx)[i]*dmuO2Sdm[k];
                                eMatrix[k][colRow] += -(oxygen.solToOx)[i]*dmuO2Sdm[k];
                            }
                        }
                    }

                    if (hasLiquid) {                                           /* V-fO2 fix new */
#ifdef PHMELTS_ADJUSTMENTS
                        if (!silminState->fo2Alt) {
                            for (nl=0, k=0; nl<silminState->nLiquidCoexist; nl++) { /* V-fO2 fix new */
                                for (j=0; j<nlc; j++) {                                /* V-fO2 fix new */
                                    if ((silminState->liquidComp)[nl][j] != 0.0) {       /* V-fO2 fix new */
                                        eMatrix[colRow][k] *= (oxygen.solToOx)[i];         /* V-fO2 fix new */
                                        eMatrix[k][colRow] *= (oxygen.solToOx)[i];         /* V-fO2 fix new */
                                        k++;                                               /* V-fO2 fix new */
                                    }                                                    /* V-fO2 fix new */
                                }                                                      /* V-fO2 fix new */
                            }                                                        /* V-fO2 fix new */
                        }
#endif
                    /***********************************************************************************************
                     Asimow executed this for all cases. Ghiorso changed it for the liquid absent case. His l -> j
                    ***********************************************************************************************/
                    } else /* (!hasLiquid) */ {                                                                       /* V-fO2 fix was: if (!hasLiquid) { */
                        for (j=0; j<conCols-(isenthalpic||isentropic||isochoric ? 1 : 0); j++) {
                            if (isenthalpic) {
                                eMatrix[colRow][j]       -= mO2S*d2muO2Sdm2[colRow][j]/silminState->T;
                                eMatrixfO2[0][colRow][j] += -d2muO2Sdm2[colRow][j];
                            } else {
                                eMatrix[colRow][j]       -= mO2S*d2muO2Sdm2[colRow][j];
                                eMatrixfO2[0][colRow][j] += -d2muO2Sdm2[colRow][j];
                            }
                        }
                    }
                    /* END */
                }
#ifdef PHMELTS_ADJUSTMENTS
                if (H2Obuffer) {
                    bMatrix[colRow][0] += (solids[i].solToOx)[iBulkH2O]*muH2O;
                    for (nl=0,l=0; nl< silminState->nLiquidCoexist ;nl++) {
                        for (k=0; k<nlc; k++) if ((silminState->liquidComp)[nl][k] > 0.0) {
                            eMatrix[colRow][l] += -(solids[i].solToOx)[iBulkH2O]*dmuH2Odm[nl][k];
                            eMatrix[l][colRow] += -(solids[i].solToOx)[iBulkH2O]*dmuH2Odm[nl][k];
                            l++;
                        }
                    }
                }
#endif

                colRow++;

            } else {
                int nr = solids[i].nr;
                int na = solids[i].na;

                /* Allocate storage for this solid */
                rSol      = vector_alloc(nr);
                mSol      = vector_alloc(na);
                drdmSol   = matrix_alloc(nr, na);
                d2rdm2Sol = (double ***) malloc((size_t) nr*sizeof(double **));
                for (j=0; j<nr; j++) d2rdm2Sol[j] = matrix_alloc(na, na);
                dpMixSol  = vector_alloc(nr);
                d2pMixSol = matrix_alloc(nr, nr);
                dpSol  = vector_alloc(na);
                d2pSol = matrix_alloc(na, na);
                if (isenthalpic || isentropic || isochoric) {
                    dpMixTmp  = vector_alloc(nr);
                    dpMixCon  = vector_alloc(nr);
                    d2pMixTmp = matrix_alloc(nr, nr);
                    dpTmp     = vector_alloc(na);
                    dpCon     = vector_alloc(na);
                    d2pTmp    = matrix_alloc(na, na);
                } else {
                    dpMixTmp  = (double *)  NULL;
                    dpMixCon  = (double *)  NULL;
                    d2pMixTmp = (double **) NULL;
                    dpTmp     = (double *)  NULL;
                    dpCon     = (double *)  NULL;
                    d2pTmp    = (double **) NULL;
                }

                /* loop over all coexisting phases of this type */
                for (n=0; n<ns; n++) {
                    numSSsolids++;
                    /* Obtain r (indep var) for solid, as well as dr/dm and d2r/dm2 */
                    for (j=0; j<na; j++) mSol[j] = (silminState->solidComp)[i+1+j][n];
                    (*solids[i].convert)(SECOND, THIRD | FIFTH | SIXTH, silminState->T, silminState->P, NULL, mSol, rSol, NULL, drdmSol, d2rdm2Sol, NULL, NULL);

                    /* Obtain the molar potential and its associated derivatives for the solid */
                    if (isenthalpic) {
                        /******************************************************************
                         *Sol quantities contain negative entropy and it's derivatives
                         *Tmp quantities contain enthalpy and it's derivatives
                         *Con quantities contain cpMix, dcpMix/dn, dcpMix/dt
                        ******************************************************************/
                        (*solids[i].smix)(FIRST | SECOND | THIRD, silminState->T, silminState->P, rSol, &pMixSol, dpMixSol, d2pMixSol);
                        (*solids[i].gmix)(FIRST | SECOND | THIRD, silminState->T, silminState->P, rSol, &pMixTmp, dpMixTmp, d2pMixTmp, NULL);
                        pMixTmp += silminState->T*pMixSol;
                        for (j=0; j<nr; j++) {
                            dpMixTmp[j] += silminState->T*dpMixSol[j];
                            for(k=0; k<nr; k++) d2pMixTmp[j][k] += silminState->T*d2pMixSol[j][k];
                        }
                        pMixSol *= -1.0;
                        for (j=0; j<nr; j++) {
                            dpMixSol[j] *= -1.0;
                            for (k=0; k<nr; k++) d2pMixSol[j][k] *= -1.0;
                        }
                        (*solids[i].cpmix)(FIRST | SECOND | THIRD, silminState->T, silminState->P, rSol, &pMixCon, &d2pMixCon, dpMixCon);
                    } else if (isentropic) {
                        /******************************************************************
                         *Sol quantities contain enthalpy and it's derivatives
                         *Tmp quantities contain entropy  and it's derivatives
                         *Con quantities contain cpMix, dcpMix/dn, dcpMix/dt
                        ******************************************************************/
                        (*solids[i].gmix)(FIRST | SECOND | THIRD, silminState->T, silminState->P, rSol, &pMixSol, dpMixSol, d2pMixSol, NULL);
                        (*solids[i].smix)(FIRST | SECOND | THIRD, silminState->T, silminState->P, rSol, &pMixTmp, dpMixTmp, d2pMixTmp);
                        pMixSol += silminState->T*pMixTmp;
                        for (j=0; j<nr; j++) {
                            dpMixSol[j] += silminState->T*dpMixTmp[j];
                            for(k=0; k<nr; k++) d2pMixSol[j][k] += silminState->T*d2pMixTmp[j][k];
                        }
                        (*solids[i].cpmix)(FIRST | SECOND | THIRD, silminState->T, silminState->P, rSol, &pMixCon, &d2pMixCon, dpMixCon);
                    } else if (isochoric) {
                        /******************************************************************
                         *Sol quantities contain Helmholtz free energy and it's derivatives
                         *Tmp quantities contain volume and it's derivatives
                         *Con quantities contain dvMix/dp, d2vMix/dpdn, d2vMix/dp2
                        ******************************************************************/
                        (*solids[i].gmix)(FIRST | SECOND | THIRD, silminState->T, silminState->P, rSol, &pMixSol, dpMixSol, d2pMixSol, NULL);
                        (*solids[i].vmix)(FIRST | SECOND | THIRD | FIFTH | EIGHTH | TENTH,
                            silminState->T, silminState->P, rSol, &pMixTmp, dpMixTmp, d2pMixTmp, NULL, &pMixCon, NULL, NULL, &d2pMixCon, NULL, dpMixCon);
                        pMixSol -= silminState->P*pMixTmp;
                        for (j=0; j<nr; j++) {
                            dpMixSol[j] -= silminState->P*dpMixTmp[j];
                            for(k=0; k<nr; k++) d2pMixSol[j][k] -= silminState->P*d2pMixTmp[j][k];
                        }
                    } else {
                        /******************************************************************
                         *Sol quantities contain Gibbs free energy and it's derivatives
                         *Tmp quantities are undefined
                         *Con quantities are undefined
                        ******************************************************************/
                        (*solids[i].gmix)(FIRST | SECOND | THIRD, silminState->T, silminState->P, rSol, &pMixSol, dpMixSol, d2pMixSol, NULL);
                    }

                    /* Obtain the extensive potential gradient                          */
                    mTotal = (silminState->solidComp)[i][n];
                    intenToExtenGradient(pMixSol, dpMixSol, nr, dpSol, na, mTotal, drdmSol);

                    if (dpTmp != (double *) NULL) intenToExtenGradient(pMixTmp, dpMixTmp, nr, dpTmp, na, mTotal, drdmSol);
                    if (dpCon != (double *) NULL) intenToExtenGradient(pMixCon, dpMixCon, nr, dpCon, na, mTotal, drdmSol);

                    /* Add in the standard state contribution                           */
                    if (isenthalpic) {
                        for (j=0, pCon=mTotal*pMixCon, d2pCon=mTotal*d2pMixCon; j<na; j++) {
                            dpSol[j] -= (solids[i+1+j].cur).s;
                            dpTmp[j] += (solids[i+1+j].cur).h;
                            pCon     += mSol[j]*(solids[i+1+j].cur).cp;
                            dpCon[j] += (solids[i+1+j].cur).cp;
                            d2pCon   += mSol[j]*(solids[i+1+j].cur).dcpdt;
                        }
                    } else if (isentropic) {
                        for (j=0, pCon=mTotal*pMixCon, d2pCon=mTotal*d2pMixCon; j<na; j++) {
                            dpSol[j] += (solids[i+1+j].cur).h;
                            dpTmp[j] += (solids[i+1+j].cur).s;
                            pCon     += mSol[j]*(solids[i+1+j].cur).cp;
                            dpCon[j] += (solids[i+1+j].cur).cp;
                            d2pCon   += mSol[j]*(solids[i+1+j].cur).dcpdt;
                        }
                    } else if (isochoric) {
                        for (j=0, pCon=mTotal*pMixCon, d2pCon=mTotal*d2pMixCon; j<na; j++) {
                            dpSol[j] += (solids[i+1+j].cur).g
                                                - silminState->P*(solids[i+1+j].cur).v;
                            dpTmp[j] += (solids[i+1+j].cur).v;
                            pCon     += mSol[j]*(solids[i+1+j].cur).dvdp;
                            dpCon[j] += (solids[i+1+j].cur).dvdp;
                            d2pCon   += mSol[j]*(solids[i+1+j].cur).d2vdp2;
                        }
                    } else {
                        for (j=0; j<na; j++) dpSol[j] += (solids[i+1+j].cur).g;
                    }

                    /* Obtain the extensive potential Hessian                           */
                    intenToExtenHessian(pMixSol, dpMixSol, d2pMixSol, nr, d2pSol, na, mTotal, drdmSol, d2rdm2Sol);

                    if (d2pTmp != (double **) NULL) intenToExtenHessian(pMixTmp, dpMixTmp, d2pMixTmp, nr, d2pTmp, na, mTotal, drdmSol, d2rdm2Sol);

                    /********************************************************************
                     Construct and store the gradient and Hessian contribution to the
                     total potential for the solid components.
                     eMatrix[][]  contains the Hessian and
                     bMatrix[][0] contains the - gradient[] + Hessian * solid composition
                            or only the - gradient[] if non-linear constraints are active
                    *******************************************************************/

                    for (j=0, m=0; j<na; j++) {
                        if ((silminState->solidComp)[i+1+j][n] != 0.0) {

                            bMatrix[colRow+m][0] += -dpSol[j];

                            /* mu O2 constraint contribution to the solid phase             */
                            if (silminState->fo2Path != FO2_NONE) {
                                if (isenthalpic)       {
                                    bMatrix[colRow+m][0] += hasLiquid ? (oxygen.solToOx)[i+1+j]*muO2E/silminState->T
                                                                            : (oxygen.solToOx)[i+1+j]*muO2S/silminState->T + mO2S*dmuO2Sdm[colRow+m]/silminState->T;
                                    if (!hasLiquid) {
                                        for (k=0; k<conCols-1; k++) {
                                            eMatrix[colRow+m][k] += -(oxygen.solToOx)[i+1+j]*dmuO2Sdm[k]/silminState->T;
                                            eMatrix[k][colRow]   += -(oxygen.solToOx)[i+1+j]*dmuO2Sdm[k]/silminState->T;
                                        }
                                        eMatrixfO2[0][indexT][colRow+m] += -d2muO2Sdmdt[colRow+m];
                                        eMatrixfO2[0][colRow+m][indexT] += -d2muO2Sdmdt[colRow+m];
                                    }
                                    eMatrix[colRow+m][indexT] += hasLiquid ? (oxygen.solToOx)[i+1+j]*(muO2E/SQUARE(silminState->T) - dmuO2Edt/silminState->T)
		                                         : (oxygen.solToOx)[i+1+j]*(muO2S/SQUARE(silminState->T) - dmuO2Sdt/silminState->T)
		                                         + mO2S*(dmuO2Sdm[colRow+m]/SQUARE(silminState->T) - d2muO2Sdmdt[colRow+m]/silminState->T);
                                    eMatrix[indexT][colRow+m] += hasLiquid ? (oxygen.solToOx)[i+1+j]*(muO2E/SQUARE(silminState->T) - dmuO2Edt/silminState->T)
		                                         : (oxygen.solToOx)[i+1+j]*(muO2S/SQUARE(silminState->T) - dmuO2Sdt/silminState->T)
		                                         + mO2S*(dmuO2Sdm[colRow+m]/SQUARE(silminState->T) - d2muO2Sdmdt[colRow+m]/silminState->T);
                                } else if (isentropic) {
                                    bMatrix[colRow+m][0] += hasLiquid ? (oxygen.solToOx)[i+1+j]*muO2E : (oxygen.solToOx)[i+1+j]*muO2S + mO2S*dmuO2Sdm[colRow+m];
                                    if (!hasLiquid) {
                                        for (k=0; k<conCols-1; k++) {
                                            eMatrix[colRow+m][k] += -(oxygen.solToOx)[i+1+j]*dmuO2Sdm[k];
                                            eMatrix[k][colRow]   += -(oxygen.solToOx)[i+1+j]*dmuO2Sdm[k];
                                        }
                                        eMatrixfO2[0][indexT][colRow+m] += -d2muO2Sdmdt[colRow+m];
                                        eMatrixfO2[0][colRow+m][indexT] += -d2muO2Sdmdt[colRow+m];
                                    }
                                    eMatrix[colRow+m][indexT] += hasLiquid ? -(oxygen.solToOx)[i+1+j]*dmuO2Edt : -(oxygen.solToOx)[i+1+j]*dmuO2Sdt - mO2S*d2muO2Sdmdt[colRow+m];
                                    eMatrix[indexT][colRow+m] += hasLiquid ? -(oxygen.solToOx)[i+1+j]*dmuO2Edt : -(oxygen.solToOx)[i+1+j]*dmuO2Sdt - mO2S*d2muO2Sdmdt[colRow+m];
                                } else if (isochoric)  {
#ifdef PHMELTS_ADJUSTMENTS
                                    if (silminState->fo2Alt) {
                                        bMatrix[colRow+m][0] += hasLiquid ? (oxygen.solToOx)[i+1+j]*muO2E
                                                                    : (oxygen.solToOx)[i+1+j]*muO2S + mO2S*dmuO2Sdm[colRow+m];
                                    } else {
#endif
                                    bMatrix[colRow+m][0] += hasLiquid ? (oxygen.solToOx)[i+1+j]*muO2L[0]                         /* V-fO2 fix was: (oxygen.solToOx)[i+1+j]*muO2E */
                                                                    : (oxygen.solToOx)[i+1+j]*muO2S + mO2S*dmuO2Sdm[colRow+m];
#ifdef PHMELTS_ADJUSTMENTS
                                    }
#endif
                                    if (!hasLiquid) {
                                        for (k=0; k<conCols-1; k++) {
                                            eMatrix[colRow+m][k] += -(oxygen.solToOx)[i+1+j]*dmuO2Sdm[k];
                                            eMatrix[k][colRow]   += -(oxygen.solToOx)[i+1+j]*dmuO2Sdm[k];
                                        }
                                        eMatrixfO2[0][indexP][colRow+m] += -d2muO2Sdmdp[colRow+m];
                                        eMatrixfO2[0][colRow+m][indexP] += -d2muO2Sdmdp[colRow+m];
                                    }
#ifdef PHMELTS_ADJUSTMENTS
                                    if (silminState->fo2Alt) {
                                        eMatrix[colRow+m][indexP] += hasLiquid ? -(oxygen.solToOx)[i+1+j]*dmuO2Edp
                                                            : -(oxygen.solToOx)[i+1+j]*dmuO2Sdp - mO2S*d2muO2Sdmdp[colRow+m];
                                        eMatrix[indexP][colRow+m] += hasLiquid ? -(oxygen.solToOx)[i+1+j]*dmuO2Edp
                                                            : -(oxygen.solToOx)[i+1+j]*dmuO2Sdp - mO2S*d2muO2Sdmdp[colRow+m];
                                    } else {
#endif
                                    eMatrix[colRow+m][indexP] += hasLiquid ? -(oxygen.solToOx)[i+1+j]*dmuO2Ldp[0]                            /* V-fO2 fix was: -(oxygen.solToOx)[i+1+j]*dmuO2Edp */
		                                         : -(oxygen.solToOx)[i+1+j]*dmuO2Sdp - mO2S*d2muO2Sdmdp[colRow+m];
                                    eMatrix[indexP][colRow+m] += hasLiquid ? -(oxygen.solToOx)[i+1+j]*dmuO2Ldp[0]                            /* V-fO2 fix was: -(oxygen.solToOx)[i+1+j]*dmuO2Edp */
		                                         : -(oxygen.solToOx)[i+1+j]*dmuO2Sdp - mO2S*d2muO2Sdmdp[colRow+m];
#ifdef PHMELTS_ADJUSTMENTS
                                    }
#endif
                                } else {
#ifdef PHMELTS_ADJUSTMENTS
                                    if (silminState->fo2Alt) {
                                        bMatrix[colRow+m][0] += hasLiquid ? (oxygen.solToOx)[i+1+j]*muO2E
                                                                                                            : (oxygen.solToOx)[i+1+j]*muO2S + mO2S*dmuO2Sdm[colRow+m];
                                    } else {
#endif
                                    bMatrix[colRow+m][0] += hasLiquid ? (oxygen.solToOx)[i+1+j]*muO2L[0]                                     /* V-fO2 fix was: (oxygen.solToOx)[i+1+j]*muO2E*/
                                                                            : (oxygen.solToOx)[i+1+j]*muO2S + mO2S*dmuO2Sdm[colRow+m];
#ifdef PHMELTS_ADJUSTMENTS
                                    }
#endif
                                    if (!hasLiquid) {
                                        for (k=0;k<conCols-1; k++) {
                                            eMatrix[colRow+m][k] += -(oxygen.solToOx)[i+1+j]*dmuO2Sdm[k];
                                            eMatrix[k][colRow]   += -(oxygen.solToOx)[i+1+j]*dmuO2Sdm[k];
                                        }
                                    }
                                }

                                if (hasLiquid) {                                          /* V-fO2 fix new */
#ifdef PHMELTS_ADJUSTMENTS
                                if (!silminState->fo2Alt) {
                                    int ii = 0;                                             /* V-fO2 fix new */
                                    for (nl=0; nl<silminState->nLiquidCoexist; nl++) {     /* V-fO2 fix new */
                                        for (k=0; k<nlc; k++) {                               /* V-fO2 fix new */
                                            if ((silminState->liquidComp)[nl][k] != 0.0) {      /* V-fO2 fix new */
                                                eMatrix[colRow+m][ii] *= (oxygen.solToOx)[i+1+j]; /* V-fO2 fix new */
                                                eMatrix[ii][colRow+m] *= (oxygen.solToOx)[i+1+j]; /* V-fO2 fix new */
                                                ii++;                                             /* V-fO2 fix new */
                                            }                                                   /* V-fO2 fix new */
                                        }                                                     /* V-fO2 fix new */
                                    }                                                       /* V-fO2 fix new */
                                }
#endif
                                /*********************************************************************************************************
                                 Asimow executed this for all cases. Ghiorso changed it for the liquid absent case. His j -> k; His l -> k
                                *********************************************************************************************************/
                                } else /* (!hasLiquid) */ {                                                                                 /* V-fO2 fix was: if (!hasLiquid) */
                                    for(k=0; k<conCols-(isenthalpic||isentropic||isochoric ? 1 : 0); k++){
                                        if (isenthalpic) {
                                            eMatrix[colRow+m][k]       -= mO2S*d2muO2Sdm2[colRow+m][k]/silminState->T;
                                            eMatrixfO2[0][colRow+m][k] += -d2muO2Sdm2[colRow+m][k];
                                        } else {
                                            eMatrix[colRow+m][k]       -= mO2S*d2muO2Sdm2[colRow+m][k];
                                            eMatrixfO2[0][colRow+m][k] += -d2muO2Sdm2[colRow+m][k];
                                        }
                                    }
                                }
                                /* END */
                            }

#ifdef PHMELTS_ADJUSTMENTS
                            /* muH2O contribution to the solid phase */
                            if (H2Obuffer) {
                                bMatrix[colRow+m][0] += (solids[i+1+j].solToOx)[iBulkH2O]*muH2O;
                                for (nl=0, l=0; nl< silminState->nLiquidCoexist ; nl++) {
                                    for (k=0; k<nlc; k++) {
                                        if ((silminState->liquidComp)[nl][k] > 0.0) {
                                            eMatrix[colRow+m][l] +=
                                                        -(solids[i+1+j].solToOx)[iBulkH2O]*dmuH2Odm[nl][k];
                                            eMatrix[l][colRow+m] +=
                                                        -(solids[i+1+j].solToOx)[iBulkH2O]*dmuH2Odm[nl][k];
                                            l++;
                                        }
                                    }
                                }
                            }
#endif

                            for (k=0, l=0; k<na; k++) {
                                if ((silminState->solidComp)[i+1+k][n] != 0.0) {
                                    eMatrix[colRow+m][colRow+l] += d2pSol[j][k];
                                    bMatrix[colRow+m][0]        += (hasNlCon) ? 0.0 : (silminState->solidComp)[i+1+k][n]*d2pSol[j][k];
                                    if (isenthalpic || isentropic) eMatrixT[colRow+m][colRow+l] -= d2pTmp[j][k];
                                    if (isochoric)                 eMatrixP[colRow+m][colRow+l] -= d2pTmp[j][k];
                                    l++;
                                }
                            }

                            if (isenthalpic) {
                                eMatrix[colRow+m][indexT]  += -dpCon[j]/silminState->T;
                                eMatrix[indexT][colRow+m]  += -dpCon[j]/silminState->T;
                                eMatrixT[colRow+m][indexT] += -dpCon[j];
                                eMatrixT[indexT][colRow+m] += -dpCon[j];
                            } else if (isentropic) {
                                eMatrix[colRow+m][indexT]  +=  dpCon[j];
                                eMatrix[indexT][colRow+m]  +=  dpCon[j];
                                eMatrixT[colRow+m][indexT] += -dpCon[j]/silminState->T;
                                eMatrixT[indexT][colRow+m] += -dpCon[j]/silminState->T;
                            } else if (isochoric) {
                                eMatrix[colRow+m][indexP]  += -silminState->P*dpCon[j];
                                eMatrix[indexP][colRow+m]  += -silminState->P*dpCon[j];
                                eMatrixP[colRow+m][indexP] += -dpCon[j];
                                eMatrixP[indexP][colRow+m] += -dpCon[j];
                            }

                            m++;
                        }
                    }
                    colRow += m;

                    if (isenthalpic) {
                        bMatrix[indexT][0]       += pCon/silminState->T;
                        eMatrix[indexT][indexT]  += pCon/SQUARE(silminState->T) - d2pCon/silminState->T;
                        eMatrixT[indexT][indexT] += -d2pCon;
                    } else if (isentropic) {
                        bMatrix[indexT][0]       += -pCon;
                        eMatrix[indexT][indexT]  += d2pCon;
                        eMatrixT[indexT][indexT] += pCon/SQUARE(silminState->T) - d2pCon/silminState->T;
                    } else if (isochoric) {
                        bMatrix[indexP][0]       += silminState->P*pCon;
                        eMatrix[indexP][indexP]  += -pCon - silminState->P*d2pCon;
                        eMatrixP[indexP][indexP] += -d2pCon;
                    }

                } /* end loop on coexisting phases */

                /* Delete storage for the solid                                       */
                matrix_free(d2pSol, na, na);
                vector_free(dpSol, na);
                matrix_free(d2pMixSol, nr, nr);
                vector_free(dpMixSol, nr);
                for (j=0; j<nr; j++) matrix_free(d2rdm2Sol[j], na, na);
                free(d2rdm2Sol);
                matrix_free(drdmSol, nr, na);
                vector_free(rSol, nr);
                vector_free(mSol, na);
                if (dpMixTmp  != (double *)  NULL) vector_free(dpMixTmp,  nr);
                if (dpMixCon  != (double *)  NULL) vector_free(dpMixCon,  nr);
                if (d2pMixTmp != (double **) NULL) matrix_free(d2pMixTmp, nr, nr);
                if (dpTmp     != (double *)  NULL) vector_free(dpTmp,     na);
                if (dpCon     != (double *)  NULL) vector_free(dpCon,     na);
                if (d2pTmp    != (double **) NULL) matrix_free(d2pTmp,    na, na);

            }
        }
    }
    if (silminState->fo2Path != FO2_NONE) {
        /**************************************************************************
         Free oxygen constraint specific storage
         *************************************************************************/
        if (hasLiquid) {
            vector_free(mO2L,     silminState->nLiquidCoexist);
            vector_free(muO2L,    silminState->nLiquidCoexist);
            matrix_free(dmuO2Ldm, silminState->nLiquidCoexist, nlc);
            for (nl=0; nl<silminState->nLiquidCoexist; nl++) matrix_free(d2muO2Ldm2[nl], nlc, nlc);
            free(d2muO2Ldm2);
        } else {
            vector_free(dmuO2Sdm,   2*npc);
            matrix_free(d2muO2Sdm2, 2*npc, 2*npc);
        }
        if (isenthalpic || isentropic) {
            if (hasLiquid) {
                matrix_free(d2muO2Ldmdt, silminState->nLiquidCoexist, nlc);
                vector_free(dmuO2Ldt,    silminState->nLiquidCoexist);
                vector_free(d2muO2Ldt2,  silminState->nLiquidCoexist);
            } else {
                vector_free(d2muO2Sdmdt, 2*npc);
            }
        } else if (isochoric) {
            if (hasLiquid) {
                matrix_free(d2muO2Ldmdp, silminState->nLiquidCoexist, nlc);
                vector_free(dmuO2Ldp,    silminState->nLiquidCoexist);
                vector_free(d2muO2Ldp2,  silminState->nLiquidCoexist);
            } else {
                vector_free(d2muO2Sdmdp, 2*npc);
            }
        }
    }
#ifdef PHMELTS_ADJUSTMENTS
    if (H2Obuffer) {
        /* Free H2O storage */
        for (nl=0; nl<silminState->nLiquidCoexist; nl++) matrix_free(d2muH2Odm2[nl], nlc, nlc);
        free(d2muH2Odm2);
    }
#endif

    /*****************************************************************************************************
     Compute the Lagrange multipliers by backsubstitution. g = (C^^T) lambda -> (K^^T) g = (R^^T) lambda
   *****************************************************************************************************/

    if (hasNlCon) {
        double **gMatrix = matrix_alloc(conCols, 1);
        for (i=0; i<conCols; i++) gMatrix[i][0] = - bMatrix[i][0];

        /* Form (K^^T) g */
        for (i=0; i<conRows; i++) householderRowCol(HOUSEHOLDER_CALC_MODE_H2, i, i+1, conCols-1, cMatrix, i, &hVector[i], gMatrix, 0, 0);

        /* Backsolve the Upper triangular system */
        for (i=(conRows-1); i>=0; i--) {
            constraints->lambda[i] = gMatrix[i][0];
            for (j=(conRows-1); j>i; j--) constraints->lambda[i] -= cMatrix[i][j]*constraints->lambda[j];
            constraints->lambda[i] /= cMatrix[i][i];
        }

        i = conRows - 1;
        if (isochoric)   constraints->lambdaV  = constraints->lambda[i--];
        if (isenthalpic) constraints->lambdaH  = constraints->lambda[i--];
        if (isentropic)  constraints->lambdaS  = constraints->lambda[i--];
        if (silminState->fo2Path != FO2_NONE) {
            if (hasLiquid) for (nl=(silminState->nLiquidCoexist-1); nl>=0; nl--) constraints->lambdaO2[nl] = constraints->lambda[i--];
            else                                                                 constraints->lambdaO2[ 0] = constraints->lambda[i--];
        }
#ifdef PHMELTS_ADJUSTMENTS
        if (H2Obuffer) constraints->lambdaH2O = constraints->lambda[indexH2O];
#endif
        matrix_free(gMatrix, conCols, 1);
    }

    /****************************************************************************
     DEBUG mode printouts
   ***************************************************************************/

#ifdef DEBUG
    printf("hasNlCon = %s, numLiquids = %d, numEMsolids = %d, numSSsolids = %d\n", (hasNlCon) ? "T" : "F", numLiquids, numEMsolids, numSSsolids);
    printf("conRows = %3d, conCols = %3d, colRow = %3d\n", conRows, conCols, colRow);
    for (i=0; i<conCols; i++) {
        printf("  ");
        for (j=0; j<conCols; j++) printf("%s", (eMatrix[i][j] != 0.0) ? "x" : ".");
        printf("\n");
    }
    for (i=0; i<conCols; i++) {
        printf("e[%3d][%3d] = %13.6g\n", i, i, eMatrix[i][i]);
        for (j=0; j<i; j++) printf("e[%3d][%3d] = %13.6g, e[%3d][%3d] = %13.6g, del = %13.6g\n",
            i, j, eMatrix[i][j], j, i, eMatrix[j][i], eMatrix[i][j]-eMatrix[j][i]);
    }
    if (hasNlCon) {
        for (i=0; i<conRows; i++) printf("lambda[%3d] = %13.6g\n", i, (constraints->lambda)[i]);
        if (silminState->fo2Path != FO2_NONE)
            for (nl=0; nl<(hasLiquid ? silminState->nLiquidCoexist : 1); nl++) printf("lambda[O2][%d] = %13.6g\n", nl, constraints->lambdaO2[nl]);
        if (isenthalpic) printf("lambda[H] = %13.6g\n", constraints->lambdaH);
        if (isentropic)  printf("lambda[S] = %13.6g\n", constraints->lambdaS);
        if (isochoric)   printf("lambda[V] = %13.6g\n", constraints->lambdaV);
#ifdef PHMELTS_ADJUSTMENTS
        if (H2Obuffer)   printf("lambda[H2O] = %13.6g\n", constraints->lambdaH2O);
#endif
    }
#endif

    /****************************************************************************
     Modify the system Hessian for the Lagrange multiplier estimates and free
     the constraint Hessians, if they were allocated.
   ****************************************************************************/

#ifdef PHMELTS_ADJUSTMENTS
    if (H2Obuffer) {
        matrix_free(dmuH2Odm, silminState->nLiquidCoexist, nlc);
        for (i=0; i<conCols; i++) for (j=0; j<conCols; j++) eMatrix[i][j] += constraints->lambdaH2O*eMatrixH2O[i][j];
        matrix_free(eMatrixH2O, conCols-1, conCols-1);
    }
#endif
    if (silminState->fo2Path != FO2_NONE) {
        for (nl=0; nl<(hasLiquid ? silminState->nLiquidCoexist : 1); nl++) {
            for (i=0; i<conCols; i++) for (j=0; j<conCols; j++) eMatrix[i][j] += constraints->lambdaO2[nl]*eMatrixfO2[nl][i][j];

#ifdef DEBUG
    for (i=0; i<conCols; i++) {
        printf("  ");
        for (j=0; j<conCols; j++) printf("%s", (eMatrixfO2[nl][i][j] != 0.0) ? "x" : ".");
        printf("\n");
    }
    for (i=0; i<conCols; i++) {
        printf("e[%3d][%3d] = %13.6g\n", i, i, eMatrixfO2[nl][i][i]);
        for (j=0; j<i; j++) printf("e[%3d][%3d] = %13.6g, e[%3d][%3d] = %13.6g, del = %13.6g\n",
            i, j, eMatrixfO2[nl][i][j], j, i, eMatrixfO2[nl][j][i], eMatrixfO2[nl][i][j]-eMatrixfO2[nl][j][i]);
    }
#endif

            matrix_free(eMatrixfO2[nl], conCols-1, conCols-1);
        }
        free(eMatrixfO2);
    }

    if (isenthalpic) {
        colRow++;
        for (i=0; i<conCols; i++) for (j=0; j<conCols; j++) eMatrix[i][j] += constraints->lambdaH*eMatrixT[i][j];

#ifdef DEBUG
    for (i=0; i<conCols; i++) {
        printf("  ");
        for (j=0; j<conCols; j++) printf("%s", (eMatrixT[i][j] != 0.0) ? "x" : ".");
        printf("\n");
    }
    for (i=0; i<conCols; i++) {
        printf("e[%3d][%3d] = %13.6g\n", i, i, eMatrixT[i][i]);
        for (j=0; j<i; j++) printf("e[%3d][%3d] = %13.6g, e[%3d][%3d] = %13.6g, del = %13.6g\n",
            i, j, eMatrixT[i][j], j, i, eMatrixT[j][i], eMatrixT[i][j]-eMatrixT[j][i]);
    }
#endif

        matrix_free(eMatrixT, conCols-1, conCols-1);
    } else if (isentropic) {
        colRow++;
        for (i=0; i<conCols; i++) for (j=0; j<conCols; j++) eMatrix[i][j] += constraints->lambdaS*eMatrixT[i][j];

#ifdef DEBUG
    for (i=0; i<conCols; i++) {
        printf("  ");
        for (j=0; j<conCols; j++) printf("%s", (eMatrixT[i][j] != 0.0) ? "x" : ".");
        printf("\n");
    }
    for (i=0; i<conCols; i++) {
        printf("e[%3d][%3d] = %13.6g\n", i, i, eMatrixT[i][i]);
        for (j=0; j<i; j++) printf("e[%3d][%3d] = %13.6g, e[%3d][%3d] = %13.6g, del = %13.6g\n",
            i, j, eMatrixT[i][j], j, i, eMatrixT[j][i], eMatrixT[i][j]-eMatrixT[j][i]);
    }
#endif

        matrix_free(eMatrixT, conCols-1, conCols-1);
    } else if (isochoric) {
        colRow++;
        for (i=0; i<conCols; i++) for (j=0; j<conCols; j++) eMatrix[i][j] += constraints->lambdaV*eMatrixP[i][j];

#ifdef DEBUG
    for (i=0; i<conCols; i++) {
        printf("  ");
        for (j=0; j<conCols; j++) printf("%s", (eMatrixP[i][j] != 0.0) ? "x" : ".");
        printf("\n");
    }
    for (i=0; i<conCols; i++) {
        printf("e[%3d][%3d] = %13.6g\n", i, i, eMatrixP[i][i]);
        for (j=0; j<i; j++) printf("e[%3d][%3d] = %13.6g, e[%3d][%3d] = %13.6g, del = %13.6g\n",
            i, j, eMatrixP[i][j], j, i, eMatrixP[j][i], eMatrixP[i][j]-eMatrixP[j][i]);
    }
#endif

        matrix_free(eMatrixP, conCols-1, conCols-1);
    }

    /*****************************************************************************************************
     Compute e~ = (K^^T) E
    *****************************************************************************************************/

    for (i=0; i<conRows; i++) householderRowCol(HOUSEHOLDER_CALC_MODE_H2, i, i+1, conCols-1, cMatrix, i, &hVector[i], eMatrix, 0, conCols-1);

    /******************************************************************************************************
     Form the last conCols - conRows of e^ = e~ K, i.e. compute e21^ = K2^^T E K1  and  e22^ = K2^^T E K2
    ******************************************************************************************************/

    for (i=0; i<conRows; i++) householderRowRow(HOUSEHOLDER_CALC_MODE_H2, i, i+1, conCols-1, cMatrix, i, &hVector[i], eMatrix, conRows, conCols-1);

    /******************************************************************************************************
     Compute b ~ = (K^^T) B
    ******************************************************************************************************/

    for (i=0; i<conRows; i++) householderRowCol(HOUSEHOLDER_CALC_MODE_H2, i, i+1, conCols-1, cMatrix, i, &hVector[i], bMatrix, 0, 0);

    /******************************************************************************************************
     Compute b2^ = - b2~ - e21^ Y1^
    ******************************************************************************************************/

    for (i=conRows; i<conCols; i++) for (j=0; j<conRows; j++) bMatrix[i][0] -= eMatrix[i][j]*yVector[j];

    /******************************************************************************************************
     The least squares problem is now stored in the matrix and vector
     e[conRows:conCols-1][conRows:conCols-1] and b[conRows:conCols-1][0]
   ******************************************************************************************************/

#ifdef DEBUG
    for (i=0; i<conRows; i++) printf("bMatrix[%2d][0] = %13.6g  eMatrix[%2d][%2d] = %13.6g\n", i, bMatrix[i][0], i, i, eMatrix[i][i]);
    for (i=conRows; i<conCols; i++) {
        printf("bMatrix[%2d][0] = %13.6g\neMatrix[%2d][%2d:%2d] =", i, bMatrix[i][0], i, conRows, conCols-1);
        for (j=conRows; j<conCols; j++) printf(" %13.6g", eMatrix[i][j]);
        printf("\n");
    }
#endif

    if (!hasNlCon && (numLiquids == 1) && (numSSsolids == 0)) hessianType = HESSIAN_TYPE_ONE;
    return hessianType;
}

/* end of file GRADIENT_HESSIAN.C */

