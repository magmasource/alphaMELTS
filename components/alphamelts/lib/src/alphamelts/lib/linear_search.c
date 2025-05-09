const char *linear_search_ver(void) { return "$Id: linear_search.c,v 1.6 2007/11/29 05:32:12 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: linear_search.c,v $
MELTS Source Code: RCS Revision 1.6  2007/11/29 05:32:12  ghiorso
MELTS Source Code: RCS Majorite testMaj corrections.  Fixed "O2" in sol_struct_data.h to "o2".
MELTS Source Code: RCS
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
MELTS Source Code: RCS Revision 1.2  2003/04/28 20:44:46  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:58:12  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
 * Revision 3.10  1997/06/21  22:49:49  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.9  1997/05/03  20:23:28  ghiorso
 * *** empty log message ***
 *
 * Revision 3.8  1997/03/27  17:03:31  ghiorso
 * *** empty log message ***
 *
 * Revision 3.7  1996/09/24  20:33:36  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 3.6  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 3.5  1995/11/23  22:37:42  ghiorso
 * Final implementation of subsolidus fO2 buffering.
 *
 * Revision 3.4  1995/11/22  02:26:17  ghiorso
 * minor correction
 *
 * Revision 3.3  1995/11/19  18:49:02  ghiorso
 * Corrections to Asimow code for subsolidus calculations
 *
 * Revision 3.2  1995/11/01  22:40:27  ghiorso
 * Implementation of subsolidus options after Asimow.
 * Additional implementation of nepheline solid solutions.
 *
 * Revision 3.1  1995/08/18  17:46:04  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Function called by min1d to evaluate the total thermodynamic potential
**      of the system at the "point"
**        silminState->liquidComp + lambda * silminState->liquidDelta    and
**        silminState->solidComp  + lambda * silminState->solidDelta
**      where *Delta are correction vectors determined by the latest quadratic
**      minimization attempt.
**      notcomp is set to TRUE if the value of lambda causes the resulting
**      system composition vector to be infeasible.
**      (file: LINEAR_SEARCH.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  September 18, 1991
**              Extracted from file SILMIN_SUPPORT.C
**      V1.0-2  Mark S. Ghiorso  September 24, 1991
**              Altered parameter list to (*solids[].convert)
**      V1.0-3  Mark S. Ghiorso  September 28, 1991
**              Implemented accurate bound constraints for solid solutions:
**              i.e., calls to (*solids[].test)
**      V1.0-4  Mark S. Ghiorso  October 10, 1991
**              Revised output format in debug mode
**      V1.0-5  Mark S. Ghiorso  October 15, 1991
**              (1) Modified silminState->solidComp and ->solidDelta arrays
**                  to allow for immiscible solids
**              (2) Reorganization for immiscible solid code
**      V2.0-1  Mark S. Ghiorso  November 14, 1991
**              Conversion to OSF Motif V1.1.1/X11 Release 4
**      V2.0-2  Implemented bound constraint on magnitude of search
**              increment (empirical value of 2.0 adopted)
**      V3.0-1  Mark S. Ghiorso  April 29, 1992
**              (1) Added support for mu O2 constraint paths
**      V3.0-2  Mark S. Ghiorso  May 1, 1992
**              (1) Modified computation of Korzhinskii potential to account
**                  for excess oxygen (mO2 - silminState->oxygen) only
**              (2) 3.0-2.1 revoked
**      V4.0-1  Mark S. Ghiorso  May 27, 1994
**              (1) Began modifications for isenthalpic, isentropic and
**                  isochoric constraints
**              (2) Done
**      V4.0-2  Mark S. Ghiorso  June 16, 1994
**              (1) Corrected error in termination of enthalpy, entropy,
**                  and volume correction algorithms: negative residuals
**                  were causing abnormal termination
**      V4.0-3  Mark S. Ghiorso  July 2, 1994
**              (1) Alterted convergence criteria for isenthalpic, isentropic
**                  and isochoric nonlinear corrections
**      V4.1-0  Mark S. Ghiorso  February 2, 1995
**              (1) Experiments with convergence criteria in entropy matching
**      V4.1-1  Mark S. Ghiorso  February 27, 1995
**              (1) Alterted convergence criteria for isenthalpic, isentropic
**                  and isochoric nonlinear corrections
**      V5.0-1  Paul D. Asimow  April 26, 1995
**              (1) Enable subsolidus operation
**      V5.1-1  Paul D. Asimow  August 1, 1995
**              (1) Enable subsolidus fo2 buffering
**
**--
*/

#include <stdlib.h>
#include <stdio.h>

#ifndef BATCH_VERSION
#include <Xm/Xm.h>
#include <Xm/ToggleBG.h>
#include "interface.h"            /*Specific external declarations          */
#endif

#include "silmin.h"               /*SILMIN structures include file          */

#ifdef DEBUG
#undef DEBUG
#endif
#ifdef DETAIL_DEBUG
#undef DETAIL_DEBUG
#endif

#define REALLOC(x, y) (((x) == NULL) ? malloc(y) : realloc((x), (y)))

double linearSearch(double lambda, int *notcomp)
{
  static double *mSol, **mLiq, **mLiqRef, *rLiq, *rSol, *mOx, *mO2L, *muO2L;
  static SilminState *hypoState, *saveSilminState, *saveHypoState;
  static int maxNliq = 1;
  double pTotal, pTemp, mTotal, mO2S, muO2S, mO2T = 0.0, muO2E, currentT, currentP;
  int i, j = -1, k, nl, ns;
  int isenthalpic = (silminState->refEnthalpy != 0.0) && silminState->isenthalpic;
  int isentropic  = (silminState->refEntropy  != 0.0) && silminState->isentropic;
  int isochoric   = (silminState->refVolume   != 0.0) && silminState->isochoric;
  int hasLiquid   = (silminState->liquidMass  != 0.0);
#ifdef PHMELTS_ADJUSTMENTS
  int iLiqH2O, iBulkH2O, iPhaseH2O;
  int H2Obuffer   = silminState->H2Obuffer && (silminState->aH2O != 0.0) && hasLiquid;
  double mH2O, muH2O;
#endif

  if (mSol == NULL) {
    for (i=0, j=1, k=1; i<npc; i++) if (solids[i].type == PHASE) { j = MAX(j, solids[i].nr); k = MAX(k, solids[i].na); }
    mSol    = (double  *) malloc((size_t)       k*sizeof(double));
    mLiq    = (double **) malloc((size_t) maxNliq*sizeof(double *));
    mLiqRef = (double **) malloc((size_t) maxNliq*sizeof(double *));
    rLiq    = (double  *) malloc((size_t) (nlc-1)*sizeof(double));
    rSol    = (double  *) malloc((size_t)       j*sizeof(double));
    mOx     = (double  *) malloc((size_t)      nc*sizeof(double));
    mO2L    = (double  *) malloc((size_t) maxNliq*sizeof(double));
    muO2L   = (double  *) malloc((size_t) maxNliq*sizeof(double));
    for (i=0; i<maxNliq; i++) {
      mLiq[i]    = (double *) malloc((size_t) nlc*sizeof(double));
      mLiqRef[i] = (double *) malloc((size_t) nlc*sizeof(double));
    }
  }

  if (maxNliq < silminState->nLiquidCoexist) {
    mLiq    = (double **) realloc(mLiq,    (size_t) silminState->nLiquidCoexist*sizeof(double *));
    mLiqRef = (double **) realloc(mLiqRef, (size_t) silminState->nLiquidCoexist*sizeof(double *));
    for (i=maxNliq; i<silminState->nLiquidCoexist; i++) {
      mLiq[i]    = (double *) malloc((size_t) nlc*sizeof(double));
      mLiqRef[i] = (double *) malloc((size_t) nlc*sizeof(double));
    }
    mO2L  = (double *) realloc(mO2L,  (size_t) silminState->nLiquidCoexist*sizeof(double));
    muO2L = (double *) realloc(muO2L, (size_t) silminState->nLiquidCoexist*sizeof(double));
    maxNliq = silminState->nLiquidCoexist;
  }

#ifdef PHMELTS_ADJUSTMENTS
  // only for pMELTS and rhyolite-MELTS 1.0.2
  for (iLiqH2O=0; iLiqH2O<nlc; iLiqH2O++) if (!strcmp(liquid[iLiqH2O].label, "H2O")) break;
  for (iBulkH2O=0; iBulkH2O<nc; iBulkH2O++)if (!strcmp(bulkSystem[iBulkH2O].label, "H2O")) break;
  for (iPhaseH2O=0; iPhaseH2O<npc; iPhaseH2O++) if (!strcmp(solids[iPhaseH2O].label, "fluid")) break;
#endif

  *notcomp = FALSE;

  /* Reject solution if steplength exceeds an arbitrary upper bound */
  if (lambda < -2.0 || lambda > 2.0) {
#ifdef DEBUG
    printf("LinearSearch: Exiting with lambda = %g out-of-bounds.\n", lambda);
#endif
    *notcomp = TRUE;
    return 0.0;
  }

  /* compute linearSearch approximation to the solution vector */
  if (hasLiquid) {
    for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
      for (i=0; i<nlc; i++) mLiq[nl][i] = (silminState->liquidComp)[nl][i] + lambda*(silminState->liquidDelta)[nl][i];
      if(!testLiq(SIXTH, silminState->T, silminState->P, 0, 0, NULL, NULL, NULL, mLiq[nl])) {
#ifdef DEBUG
        printf("LinearSearch: Bad liquid composition (No: %d) at lambda = %g.\n", nl, lambda);
#endif
        *notcomp = TRUE;
        return 0.0;
      }
      for (i=0; i<nlc; i++) mLiqRef[nl][i] = mLiq[nl][i];
    }
  }

  if (!hasLiquid) {
    hypoState       = copySilminStateStructure(silminState, hypoState);
    saveSilminState = copySilminStateStructure(silminState, saveSilminState);
  }

  for (i=0; i<npc; i++) for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
    mTotal = (silminState->solidComp)[i][ns] + lambda*(silminState->solidDelta)[i][ns];
    if (!hasLiquid) hypoState->solidComp[i][ns] = mTotal;
    if (mTotal < 0.0) {
#ifdef DEBUG
      printf("LinearSearch: Negative Total moles of %s at lambda = %g.\n", solids[i].label, lambda);
#endif
      *notcomp = TRUE;
      return 0.0;
    } else if (mTotal > 0.0 && solids[i].na > 1) {
      for (j=0; j<solids[i].na; j++) {
        mSol[j] = (silminState->solidComp)[i+1+j][ns] + lambda*(silminState->solidDelta)[i+1+j][ns];
        if (!hasLiquid) hypoState->solidComp[i+1+j][ns] = mSol[j];
      }
      if(!(*solids[i].test)(SIXTH, silminState->T, silminState->P, 0, 0, NULL, NULL, NULL, mSol)) {
#ifdef DEBUG
        printf("LinearSearch: Bad composition for %s at lambda = %g.\n", solids[i].label, lambda);
#endif
        *notcomp = TRUE;
        return 0.0;
      }
    }
  }

  if ( (isenthalpic || isentropic) &&
       ( (silminState->T+lambda*silminState->tDelta <  298.15) ||
         (silminState->T+lambda*silminState->tDelta > 3000.00) ) )
    { *notcomp = TRUE; return 0.0;}

  if ( isochoric &&
       ( (silminState->P+lambda*silminState->pDelta <     1.013) ||
         (silminState->P+lambda*silminState->pDelta > 50000.000) ) )
    { *notcomp = TRUE; return 0.0;}

  /****************************************************************************
   Compute correction vector if constraints are non-linear
   ****************************************************************************/

#ifdef PHMELTS_ADJUSTMENTS
  if (H2Obuffer) {
    double oldH2O, mubuff, *gradH2O;
    mubuff = (solids[iPhaseH2O].cur).g + R*(silminState->T)*log(silminState->aH2O) - (liquid[iLiqH2O].cur).g;
    gradH2O = (double *) malloc(nlc*sizeof(double));
    for (ns=0, mH2O = 0.0; ns<silminState->nLiquidCoexist; ns++) {
      oldH2O = mLiq[ns][iLiqH2O];
      /* COMMENTED OUT UNTIL FULLY IMPLEMENTED IN LIQUID_H2O.C */
      /* PRODUCES COMPILER WARNING (GOOD!) */
      muH2OLiq(FIRST | SECOND, silminState->T, silminState->P, mLiq[ns], &muH2O, gradH2O, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
      while (fabs(muH2O - mubuff) > sqrt(DBL_EPSILON)) {
        mLiq[ns][iLiqH2O] -= (muH2O-mubuff)/gradH2O[iLiqH2O];
        muH2OLiq(FIRST | SECOND, silminState->T, silminState->P, mLiq[ns], &muH2O, gradH2O, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
      }
      constraints->liquidDelta[ns][iLiqH2O] = mLiq[ns][iLiqH2O] - oldH2O;
      mH2O += mLiq[ns][iLiqH2O];
    }
    muH2O += (liquid[iLiqH2O].cur).g;
    free(gradH2O);
  }
#endif

  if ((silminState->fo2Path != FO2_NONE) || isenthalpic || isentropic || isochoric) {

    if (silminState->fo2Path != FO2_NONE && !isenthalpic && !isentropic && !isochoric) {
      if (hasLiquid) {
        for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
          for (i=0; i<nc; i++) for (j=0, mOx[i]=0.0; j<nlc; j++) mOx[i] += (liquid[j].liqToOx)[i]*mLiq[nl][j];
          conLiq(FIRST | SEVENTH, FIRST, silminState->T, silminState->P, mOx, NULL, NULL, NULL, NULL, NULL, &(silminState->fo2));
          for (i=0, mO2L[nl]=0.0; i<nlc; i++) {
            for (j=0, (constraints->liquidDelta)[nl][i] = -mLiq[nl][i]; j<nc; j++) (constraints->liquidDelta)[nl][i] += (bulkSystem[j].oxToLiq)[i]*mOx[j];
            mLiq[nl][i] += (constraints->liquidDelta)[nl][i];
            mO2L[nl]    += (oxygen.liqToOx)[i]*mLiq[nl][i];
          }
          if (!testLiq(SIXTH, silminState->T, silminState->P, 0, 0, NULL, NULL, NULL, mLiq[nl])) { *notcomp = TRUE; return 0.0; }
	}
      } else {
        muO2S = silminState->fo2*(R*silminState->T*log(10.0));
        saveHypoState = copySilminStateStructure(hypoState, saveHypoState);
        silminState   = copySilminStateStructure(hypoState, silminState);
        if (!subsolidusmuO2(0, &muO2S, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)) {
	  printf("FAILURE in subsolidus fO2 buffering routine.  Removing buffer constraint.\n");
#ifdef PHMELTS_ADJUSTMENTS
	  silminState->fo2Liq = silminState->fo2Path;
#endif
	  silminState->fo2Path = FO2_NONE;
#ifndef BATCH_VERSION
          XmToggleButtonGadgetSetState(tg_path_none, True, True);
#endif
	}
        hypoState   = copySilminStateStructure(silminState,     hypoState);
        silminState = copySilminStateStructure(saveSilminState, silminState);
        for (i=0; i<npc; i++) {
          if (hypoState->nSolidCoexist[i] != 0) for (j=0;j<=solids[i].na;j++) (constraints->solidDelta)[i+j] =
              REALLOC((constraints->solidDelta)[i+j], hypoState->nSolidCoexist[i]*sizeof(double));
          for (j=0; j<hypoState->nSolidCoexist[i]; j++) {
            (constraints->solidDelta)[i][j] = hypoState->solidComp[i][j] - saveHypoState->solidComp[i][j];
            if (solids[i].na > 1) for (k=0; k<solids[i].na; k++) {
              (constraints->solidDelta)[i+1+k][j] = hypoState->solidComp[i+1+k][j] - saveHypoState->solidComp[i+1+k][j];
            }
          }
        }
      }
    }

    if (isenthalpic) {
      double hTotal = 0.0, cpTotal = 0.0, residual = DBL_MAX;
      int iter = 0;
      constraints->T = silminState->T + lambda*silminState->tDelta;

      while ((fabs(residual) > 10.0*DBL_EPSILON*fabs(silminState->refEnthalpy)) && (iter < 50)) {

        if (silminState->fo2Path != FO2_NONE && hasLiquid) {
          constraints->fo2 = getlog10fo2(constraints->T, silminState->P, silminState->fo2Path);
	  for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
            for (i=0; i<nc; i++) for (j=0, mOx[i]=0.0; j<nlc; j++) mOx[i] += (liquid[j].liqToOx)[i]*mLiq[nl][j];
            conLiq(FIRST | SEVENTH, FIRST, constraints->T, silminState->P, mOx, NULL, NULL, NULL, NULL, NULL, &(constraints->fo2));
            for (i=0, mO2L[nl]=0.0; i<nlc; i++) {
              for (j=0, mLiq[nl][i]=0.0; j<nc; j++) mLiq[nl][i] += (bulkSystem[j].oxToLiq)[i]*mOx[j];
	      (constraints->liquidDelta)[nl][i] = mLiq[nl][i] - mLiqRef[nl][i];
              mO2L[nl] += (oxygen.liqToOx)[i]*mLiq[nl][i];
            }
	  }
        } else if (silminState->fo2Path != FO2_NONE && !hasLiquid) {
          constraints->fo2 = getlog10fo2(constraints->T, silminState->P, silminState->fo2Path);
          muO2S = constraints->fo2*(R*silminState->T*log(10.0));
          hypoState->T  = constraints->T;
          saveHypoState = copySilminStateStructure(hypoState, saveHypoState);
          silminState   = copySilminStateStructure(hypoState, silminState);
          if (!subsolidusmuO2(0, &muO2S, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)) {
	    printf("FAILURE in subsolidus fO2 buffering routine.  Removing buffer constraint.\n");
#ifdef PHMELTS_ADJUSTMENTS
	    silminState->fo2Liq = silminState->fo2Path;
#endif
	    silminState->fo2Path = FO2_NONE;
#ifndef BATCH_VERSION
            XmToggleButtonGadgetSetState(tg_path_none, True, True);
#endif
	  }
          hypoState   = copySilminStateStructure(silminState,     hypoState);
          silminState = copySilminStateStructure(saveSilminState, silminState);
          for (i=0; i<npc; i++) {
	    if (hypoState->nSolidCoexist[i] != 0) for (j=0;j<=solids[i].na;j++) (constraints->solidDelta)[i+j] =
              REALLOC((constraints->solidDelta)[i+j], hypoState->nSolidCoexist[i]*sizeof(double));
            for (j=0; j<hypoState->nSolidCoexist[i]; j++) {
              (constraints->solidDelta)[i][j] = hypoState->solidComp[i][j] - saveHypoState->solidComp[i][j];
              if (solids[i].na > 1) for (k=0; k<solids[i].na; k++) {
                (constraints->solidDelta)[i+1+k][j] = hypoState->solidComp[i+1+k][j] - saveHypoState->solidComp[i+1+k][j];
              }
            }
          }
        }

        hTotal  = 0.0;
        cpTotal = 0.0;
        if (hasLiquid) {
	  for (i=0; i<nlc; i++) if (mLiq[0][i] != 0.0) gibbs(constraints->T, silminState->P, (char *) liquid[i].label,
	    &(liquid[i].ref), &(liquid[i].liq), &(liquid[i].fus), &(liquid[i].cur));
          for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
            for (i=0, mTotal= 0.0; i<nlc; i++) {
              mTotal  += mLiq[nl][i];
              hTotal  += mLiq[nl][i]*(liquid[i].cur).h;
              cpTotal += mLiq[nl][i]*(liquid[i].cur).cp;
            }
            conLiq(SECOND, THIRD, constraints->T, silminState->P, NULL, mLiq[nl], rLiq, NULL, NULL, NULL, NULL);
            hmixLiq(FIRST,  constraints->T, silminState->P, rLiq, &pTemp, NULL);
            hTotal += mTotal*pTemp;
            cpmixLiq(FIRST, constraints->T, silminState->P, rLiq, &pTemp, NULL, NULL);
            cpTotal += mTotal*pTemp;
	  }
        }
#ifdef DETAIL_DEBUG
  printf("LinearSearch: H Liq = %20.13g, Cp Liq = %20.13g\n", hTotal, cpTotal);
#endif

        for (i=0; i<npc; i++)
          for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
            mTotal = (silminState->solidComp)[i][ns]  + lambda*(silminState->solidDelta)[i][ns]
	           + (silminState->fo2Path != FO2_NONE && !hasLiquid ? constraints->solidDelta[i+1+j][ns] : 0.0);
            if (solids[i].na == 1) {
              gibbs(constraints->T, silminState->P, (char *) solids[i].label, &(solids[i].ref), NULL, NULL, &(solids[i].cur));
              hTotal  += mTotal*(solids[i].cur).h;
              cpTotal += mTotal*(solids[i].cur).cp;
#ifdef DETAIL_DEBUG
  printf("LinearSearch: H+ %s = %20.13g, Cp+ = %20.13g\n", solids[i].label, hTotal, cpTotal);
#endif
            } else {
              for (j=0; j<solids[i].na; j++) {
                mSol[j] = (silminState->solidComp)[i+1+j][ns] + lambda*(silminState->solidDelta)[i+1+j][ns]
                        + (silminState->fo2Path != FO2_NONE && !hasLiquid ? constraints->solidDelta[i+1+j][ns] : 0.0);
                gibbs(constraints->T, silminState->P, (char *) solids[i+1+j].label, &(solids[i+1+j].ref), NULL, NULL, &(solids[i+1+j].cur));
                hTotal  += mSol[j]*(solids[i+1+j].cur).h;
                cpTotal += mSol[j]*(solids[i+1+j].cur).cp;
              }
              (*solids[i].convert)(SECOND,THIRD,constraints->T,silminState->P, NULL, mSol, rSol, NULL, NULL, NULL, NULL, NULL);
              (*solids[i].hmix)(FIRST, constraints->T, silminState->P, rSol, &pTemp);
              hTotal += mTotal*pTemp;
              (*solids[i].cpmix)(FIRST, constraints->T, silminState->P, rSol, &pTemp, (double *) NULL, (double *) NULL);
              cpTotal += mTotal*pTemp;
#ifdef DETAIL_DEBUG
  printf("LinearSearch: H+ %s = %20.13g, Cp+ = %20.13g\n", solids[i].label, hTotal, cpTotal);
#endif
            }
          }

        residual = hTotal - silminState->refEnthalpy;
        constraints->T -= residual/cpTotal;
        iter++;

#ifdef DEBUG
  printf("LinearSearch: Iter = %d, T(C) = %.6f, H-Href = %g, Cp = %g\n", iter, constraints->T-273.15, residual, cpTotal);
#endif

      }

#ifdef DEBUG
  printf("LinearSearch: Iter = %d, T(C) = %.6f, H-Href = %g, Cp = %g\n", iter, constraints->T-273.15, residual, cpTotal);
#endif

    } else if (isentropic) {
      double sTotal = 0.0, cpTotal = 0.0, residual = DBL_MAX;
      int iter = 0;
      constraints->T = silminState->T + lambda*silminState->tDelta;

      while ((fabs(residual) > 10.0*DBL_EPSILON*fabs(silminState->refEntropy)) && (iter < 50)) {

#ifdef DEBUG
  printf("LinearSearch: Iter = %d, T(C) = %.2f (top of loop)\n", iter, constraints->T-273.15);
#endif

        if (silminState->fo2Path != FO2_NONE && hasLiquid) {
          constraints->fo2 = getlog10fo2(constraints->T, silminState->P, silminState->fo2Path);
	  for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
            for (i=0; i<nc; i++) for (j=0, mOx[i]=0.0; j<nlc; j++) mOx[i] += (liquid[j].liqToOx)[i]*mLiq[nl][j];
            conLiq(FIRST | SEVENTH, FIRST, constraints->T, silminState->P, mOx, NULL, NULL, NULL, NULL, NULL, &(constraints->fo2));
            for (i=0, mO2L[nl]=0.0; i<nlc; i++) {
              for (j=0, mLiq[nl][i]=0.0; j<nc; j++) mLiq[nl][i] += (bulkSystem[j].oxToLiq)[i]*mOx[j];
	      (constraints->liquidDelta)[nl][i] = mLiq[nl][i] - mLiqRef[nl][i];
              mO2L[nl] += (oxygen.liqToOx)[i]*mLiq[nl][i];
            }
	  }
        } else if (silminState->fo2Path != FO2_NONE && !hasLiquid) {
          constraints->fo2 = getlog10fo2(constraints->T, silminState->P, silminState->fo2Path);
          muO2S = constraints->fo2*(R*silminState->T*log(10.0));
          hypoState->T  = constraints->T;
          saveHypoState = copySilminStateStructure(hypoState, saveHypoState);
          silminState   = copySilminStateStructure(hypoState, silminState);
          if (!subsolidusmuO2(0, &muO2S, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)) {
	    printf("FAILURE in subsolidus fO2 buffering routine.  Removing buffer constraint.\n");
#ifdef PHMELTS_ADJUSTMENTS
      	    silminState->fo2Liq = silminState->fo2Path;
#endif
	    silminState->fo2Path = FO2_NONE;
#ifndef BATCH_VERSION
            XmToggleButtonGadgetSetState(tg_path_none, True, True);
#endif
	  }
          hypoState   = copySilminStateStructure(silminState,     hypoState);
          silminState = copySilminStateStructure(saveSilminState, silminState);
          for (i=0; i<npc; i++) {
	    if (hypoState->nSolidCoexist[i] != 0) for (j=0;j<=solids[i].na;j++) (constraints->solidDelta)[i+j] =
              REALLOC((constraints->solidDelta)[i+j], hypoState->nSolidCoexist[i]*sizeof(double));
            for (j=0; j<hypoState->nSolidCoexist[i]; j++) {
              (constraints->solidDelta)[i][j] = hypoState->solidComp[i][j] - saveHypoState->solidComp[i][j];
              if (solids[i].na > 1) for (k=0; k<solids[i].na; k++) {
                (constraints->solidDelta)[i+1+k][j] = hypoState->solidComp[i+1+k][j] - saveHypoState->solidComp[i+1+k][j];
              }
            }
          }
        }

#ifdef DEBUG
  printf("LinearSearch: Iter = %d, T(C) = %.2f (middle of loop)\n", iter, constraints->T-273.15);
#endif

        sTotal  = 0.0;
        cpTotal = 0.0;
        if (hasLiquid) {
	  for (i=0; i<nlc; i++) if (mLiq[0][i] != 0.0) gibbs(constraints->T, silminState->P, (char *) liquid[i].label,
	     &(liquid[i].ref), &(liquid[i].liq), &(liquid[i].fus), &(liquid[i].cur));
          for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
	    for (i=0, mTotal= 0.0; i<nlc; i++) {
              mTotal  += mLiq[nl][i];
              sTotal  += mLiq[nl][i]*(liquid[i].cur).s;
              cpTotal += mLiq[nl][i]*(liquid[i].cur).cp;
            }
            conLiq(SECOND, THIRD, constraints->T, silminState->P, NULL, mLiq[nl], rLiq, NULL, NULL, NULL, NULL);
            smixLiq(FIRST,  constraints->T, silminState->P, rLiq, &pTemp, NULL, NULL, NULL);
            sTotal += mTotal*pTemp;
            cpmixLiq(FIRST, constraints->T, silminState->P, rLiq, &pTemp, NULL, NULL);
            cpTotal += mTotal*pTemp;
	  }
        }

        for (i=0; i<npc; i++)
          for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
            mTotal = (silminState->solidComp)[i][ns] + lambda*(silminState->solidDelta)[i][ns]
                   + (silminState->fo2Path != FO2_NONE && !hasLiquid ? constraints->solidDelta[i+1+j][ns] : 0.0);
            if (solids[i].na == 1) {
              gibbs(constraints->T, silminState->P, (char *) solids[i].label, &(solids[i].ref), NULL, NULL, &(solids[i].cur));
              sTotal  += mTotal*(solids[i].cur).s;
              cpTotal += mTotal*(solids[i].cur).cp;
            } else {
              for (j=0; j<solids[i].na; j++) {
                mSol[j] = (silminState->solidComp)[i+1+j][ns] + lambda*(silminState->solidDelta)[i+1+j][ns]
                        + (silminState->fo2Path != FO2_NONE && !hasLiquid ? constraints->solidDelta[i+1+j][ns] : 0.0);
                gibbs(constraints->T, silminState->P, (char *) solids[i+1+j].label, &(solids[i+1+j].ref), NULL, NULL, &(solids[i+1+j].cur));
                sTotal  += mSol[j]*(solids[i+1+j].cur).s;
                cpTotal += mSol[j]*(solids[i+1+j].cur).cp;
              }
              (*solids[i].convert)(SECOND, THIRD, constraints->T, silminState->P, NULL, mSol, rSol, NULL, NULL, NULL, NULL, NULL);
              (*solids[i].smix)(FIRST, constraints->T, silminState->P, rSol, &pTemp, NULL, NULL);
              sTotal += mTotal*pTemp;
              (*solids[i].cpmix)(FIRST, constraints->T, silminState->P, rSol, &pTemp, NULL, NULL);
              cpTotal += mTotal*pTemp;
            }
        }

        residual = sTotal - silminState->refEntropy;
        constraints->T -= residual*constraints->T/cpTotal;
        iter++;

#ifdef DEBUG
  printf("LinearSearch: Iter = %d, T(C) = %.2f, S-Sref = %g, Cp = %g\n", iter, constraints->T-273.15, residual, cpTotal);
#endif

      }

#ifdef DEBUG
  printf("LinearSearch: Iter = %d, T(C) = %.2f, S-Sref = %g, Cp = %g\n", iter, constraints->T-273.15, residual, cpTotal);
#endif

    } else if (isochoric) {
      double vTotal = 0.0, dvdpTotal = 0.0, dpTemp, residual = DBL_MAX;
      int iter = 0;
#ifdef DEBUG
      double vOnInput = 0.0, vOnOutput = 0.0;
#endif
      constraints->P = silminState->P + lambda*silminState->pDelta;

#ifdef DEBUG
  printf("LinearSearch: Iter = %d, P = %.2f, Vref = %e (100*eps*Vref = %g)\n", iter, constraints->P, silminState->refVolume,
         100.0*DBL_EPSILON*fabs(silminState->refVolume));
#endif

      while ((fabs(residual) > 100.0*DBL_EPSILON*fabs(silminState->refVolume)) && (iter < 50)) {

        if (silminState->fo2Path != FO2_NONE && hasLiquid) {
          constraints->fo2 = getlog10fo2(silminState->T, constraints->P, silminState->fo2Path);
	  for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
#ifdef DEBUG
	    for (i=0, vOnInput=0.0; i<nlc; i++) if (mLiq[nl][i] != 0.0) {
	      gibbs(silminState->T, constraints->P, (char *) liquid[i].label, &(liquid[i].ref), &(liquid[i].liq), &(liquid[i].fus), &(liquid[i].cur));
	      vOnInput += mLiq[nl][i]*(liquid[i].cur).v;
	    }
#endif
            for (i=0; i<nc; i++) for (j=0, mOx[i]=0.0; j<nlc; j++) mOx[i] += (liquid[j].liqToOx)[i]*mLiq[nl][j];
            conLiq(FIRST | SEVENTH, FIRST, silminState->T, constraints->P, mOx, NULL, NULL, NULL, NULL, NULL, &(constraints->fo2));
            for (i=0, mO2L[nl]=0.0; i<nlc; i++) {
              for (j=0, mLiq[nl][i]=0.0; j<nc; j++) mLiq[nl][i] += (bulkSystem[j].oxToLiq)[i]*mOx[j];
	      (constraints->liquidDelta)[nl][i] = mLiq[nl][i] - mLiqRef[nl][i];
              mO2L[nl] += (oxygen.liqToOx)[i]*mLiq[nl][i];
            }
#ifdef DEBUG
	    for (i=0, vOnOutput=0.0; i<nlc; i++) if (mLiq[nl][i] != 0.0) {
	      gibbs(silminState->T, constraints->P, (char *) liquid[i].label, &(liquid[i].ref), &(liquid[i].liq), &(liquid[i].fus), &(liquid[i].cur));
	      vOnOutput += mLiq[nl][i]*(liquid[i].cur).v;
	    }
#endif
	  }
        } else if (silminState->fo2Path != FO2_NONE && !hasLiquid) {
          constraints->fo2 = getlog10fo2(constraints->T, silminState->P, silminState->fo2Path);
          muO2S = constraints->fo2*(R*silminState->T*log(10.0));
          hypoState->P  = constraints->P;
          saveHypoState = copySilminStateStructure(hypoState, saveHypoState);
          silminState   = copySilminStateStructure(hypoState, silminState);
          if (!subsolidusmuO2(0, &muO2S, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)) {
	    printf("FAILURE in subsolidus fO2 buffering routine.  Removing buffer constraint.\n");
#ifdef PHMELTS_ADJUSTMENTS
	    silminState->fo2Liq = silminState->fo2Path;
#endif
	    silminState->fo2Path = FO2_NONE;
#ifndef BATCH_VERSION
            XmToggleButtonGadgetSetState(tg_path_none, True, True);
#endif
	  }
          hypoState   = copySilminStateStructure(silminState,     hypoState);
          silminState = copySilminStateStructure(saveSilminState, silminState);
          for (i=0; i<npc; i++) {
	    if (hypoState->nSolidCoexist[i] != 0) for (j=0;j<=solids[i].na;j++) (constraints->solidDelta)[i+j] =
              REALLOC((constraints->solidDelta)[i+j], hypoState->nSolidCoexist[i]*sizeof(double));
            for (j=0; j<hypoState->nSolidCoexist[i]; j++) {
              (constraints->solidDelta)[i][j] = hypoState->solidComp[i][j] - saveHypoState->solidComp[i][j];
              if (solids[i].na > 1) for (k=0; k<solids[i].na; k++) {
                (constraints->solidDelta)[i+1+k][j] = hypoState->solidComp[i+1+k][j] - saveHypoState->solidComp[i+1+k][j];
              }
            }
          }
        }

        vTotal  = 0.0;
        dvdpTotal = 0.0;
        if (hasLiquid) {
	  for (i=0; i<nlc; i++) if (mLiq[0][i] != 0.0) gibbs(silminState->T, constraints->P, (char *) liquid[i].label,
	     &(liquid[i].ref), &(liquid[i].liq), &(liquid[i].fus), &(liquid[i].cur));
          for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
	    for (i=0, mTotal= 0.0; i<nlc; i++) {
              mTotal    += mLiq[nl][i];
              vTotal    += mLiq[nl][i]*(liquid[i].cur).v;
              dvdpTotal += mLiq[nl][i]*(liquid[i].cur).dvdp;
            }
            conLiq(SECOND, THIRD, silminState->T, constraints->P, NULL, mLiq[nl], rLiq, NULL, NULL, NULL, NULL);
            vmixLiq(FIRST | FIFTH,  silminState->T, constraints->P, rLiq, &pTemp, NULL, NULL, NULL, &dpTemp, NULL, NULL, NULL, NULL, NULL, NULL);
            vTotal    += mTotal*pTemp;
            dvdpTotal += mTotal*dpTemp;
	  }
        }

        for (i=0; i<npc; i++)
          for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
            mTotal = (silminState->solidComp)[i][ns] + lambda*(silminState->solidDelta)[i][ns]
                   + (silminState->fo2Path != FO2_NONE && !hasLiquid ? constraints->solidDelta[i+1+j][ns] : 0.0);
            if (solids[i].na == 1) {
              gibbs(silminState->T, constraints->P, (char *) solids[i].label, &(solids[i].ref), NULL, NULL, &(solids[i].cur));
              vTotal    += mTotal*(solids[i].cur).v;
              dvdpTotal += mTotal*(solids[i].cur).dvdp;
            } else {
              for (j=0; j<solids[i].na; j++) {
                mSol[j] = (silminState->solidComp)[i+1+j][ns] + lambda*(silminState->solidDelta)[i+1+j][ns]
                        + (silminState->fo2Path != FO2_NONE && !hasLiquid ? constraints->solidDelta[i+1+j][ns] : 0.0);
                gibbs(silminState->T, constraints->P, (char *) solids[i+1+j].label, &(solids[i+1+j].ref), NULL, NULL, &(solids[i+1+j].cur));
                vTotal    += mSol[j]*(solids[i+1+j].cur).v;
                dvdpTotal += mSol[j]*(solids[i+1+j].cur).dvdp;
              }
              (*solids[i].convert)(SECOND, THIRD, silminState->T, constraints->P, NULL, mSol, rSol, NULL, NULL, NULL, NULL, NULL);
              (*solids[i].vmix)(FIRST | FIFTH, silminState->T, constraints->P, rSol, &pTemp, NULL, NULL, NULL, &dpTemp, NULL, NULL, NULL, NULL, NULL);
              vTotal    += mTotal*pTemp;
              dvdpTotal += mTotal*dpTemp;
            }
          }

        residual = vTotal - silminState->refVolume;
        constraints->P -= residual/dvdpTotal;
        iter++;

#ifdef DEBUG
  printf("LinearSearch: Iter = %d, P = %.2f, V-Vref = %13.6e, dVdP = %13.6e", iter, constraints->P, residual, dvdpTotal);
  if (silminState->fo2Path != FO2_NONE && hasLiquid) printf(", V(post-fO2)-V(pre-fO2) = %13.6e\n", vOnOutput-vOnInput); else printf("\n");
#endif

      }

#ifdef DEBUG
  printf("LinearSearch: Iter = %d, P = %.2f, V-Vref = %13.6e, dVdP = %13.6e", iter, constraints->P, residual, dvdpTotal);
  if (silminState->fo2Path != FO2_NONE && hasLiquid) printf(", V(post-fO2)-V(pre-fO2) = %13.6e\n", vOnOutput-vOnInput); else printf("\n");
#endif

    }

    if (silminState->fo2Path != FO2_NONE && (isenthalpic || isentropic))
      gibbs(constraints->T, silminState->P, "o2", &(oxygen.ref), NULL, NULL, &(oxygen.cur));
    else if (silminState->fo2Path != FO2_NONE && isochoric)
      gibbs(silminState->T, constraints->P, "o2", &(oxygen.ref), NULL, NULL, &(oxygen.cur));
  }

  /****************************************************************************
   Compute potential for minimization
   ****************************************************************************/

  currentT = (isenthalpic || isentropic) ? constraints->T : silminState->T;
  currentP = (isochoric)                 ? constraints->P : silminState->P;

  /* Compute the chemical potential if path is O2 buffered                    */
  if (silminState->fo2Path != FO2_NONE) {
    if (hasLiquid) {
      for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
        muO2Liq(FIRST, currentT, currentP, mLiq[nl], &muO2L[nl], NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
	muO2L[nl] += (oxygen.cur).g;
      }
    }
    muO2S = (isentropic || isenthalpic || isochoric ? constraints->fo2 : silminState->fo2)*R*currentT*log(10.0) + (oxygen.cur).g;
  }
  /* These are the properties of the externally imposed buffer */
  muO2E =  (oxygen.cur).g + R*currentT*log((double) 10.0)*getlog10fo2(currentT, currentP, silminState->fo2Path);

#ifdef DEBUG
  printf("muO2L = %g, muO2E = %g, diff = %g\n", muO2L[0], muO2E, muO2L[0]-muO2E);
#endif

  pTotal = 0.0;

  /* compute total potential for the liquid -> first standard state */
  if (hasLiquid) {
    for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
      for (i=0, mTotal = 0.0; i<nlc; i++) {
        mTotal += mLiq[nl][i];
        if      (isenthalpic) pTotal -= mLiq[nl][i]*(liquid[i].cur).s;
        else if (isentropic)  pTotal += mLiq[nl][i]*(liquid[i].cur).h;
        else if (isochoric)   pTotal += mLiq[nl][i]*((liquid[i].cur).g - currentP*(liquid[i].cur).v);
        else                  pTotal += mLiq[nl][i]*(liquid[i].cur).g;
      }

      /*                                        -> mixing properties    */
      conLiq(SECOND, THIRD, currentT, currentP, NULL, mLiq[nl], rLiq, NULL, NULL, NULL, NULL);

      if (isenthalpic) {
        smixLiq(FIRST, currentT, currentP, rLiq, &pTemp, NULL, NULL, NULL);
        pTotal -= mTotal*pTemp;
      } else if (isentropic) {
        hmixLiq(FIRST, currentT, currentP, rLiq, &pTemp, NULL);
        pTotal += mTotal*pTemp;
      } else if (isochoric) {
        gmixLiq(FIRST, currentT, currentP, rLiq, &pTemp, NULL, NULL);
        pTotal += mTotal*pTemp;
        vmixLiq(FIRST, currentT, currentP, rLiq, &pTemp, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
        pTotal -= mTotal*currentP*pTemp;
      } else {
        gmixLiq(FIRST, currentT, currentP, rLiq, &pTemp, NULL, NULL);
        pTotal += mTotal*pTemp;
      }
    }
  }

  /* compute total potential for the solids */
  for (i=0, mO2S=0.0; i<npc; i++) for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
    mTotal = (silminState->solidComp)[i][ns] + lambda*(silminState->solidDelta)[i][ns]
           + (silminState->fo2Path != FO2_NONE && !hasLiquid ? constraints->solidDelta[i][ns] : 0.0);
    if (solids[i].na == 1) {
      if (isenthalpic)     pTotal -= mTotal*(solids[i].cur).s;
      else if (isentropic) pTotal += mTotal*(solids[i].cur).h;
      else if (isochoric)  pTotal += mTotal*((solids[i].cur).g - currentP*(solids[i].cur).v);
      else                 pTotal += mTotal*(solids[i].cur).g;
      if (silminState->fo2Path != FO2_NONE) mO2S += mTotal*(oxygen.solToOx)[i];
#ifdef PHMELTS_ADJUSTMENTS
      if (H2Obuffer) mH2O += mTotal*(solids[i].solToOx)[iBulkH2O];
#endif
    } else {
      for (j=0; j<solids[i].na; j++) {
        mSol[j] = (silminState->solidComp)[i+1+j][ns] + lambda*(silminState->solidDelta)[i+1+j][ns]
                + (silminState->fo2Path != FO2_NONE && !hasLiquid ? constraints->solidDelta[i+1+j][ns] : 0.0);
        if (isenthalpic)     pTotal -= mSol[j]*(solids[i+1+j].cur).s;
        else if (isentropic) pTotal += mSol[j]*(solids[i+1+j].cur).h;
        else if (isochoric)  pTotal += mSol[j]*((solids[i+1+j].cur).g - currentP*(solids[i+1+j].cur).v);
        else                 pTotal += mSol[j]*(solids[i+1+j].cur).g;
        if (silminState->fo2Path != FO2_NONE) mO2S += mSol[j]*(oxygen.solToOx)[i+1+j];
#ifdef PHMELTS_ADJUSTMENTS
        if (H2Obuffer) mH2O += mSol[j]*(solids[i+1+j].solToOx)[iBulkH2O];
#endif
      }
      (*solids[i].convert)(SECOND, THIRD, currentT, currentP, NULL, mSol, rSol, NULL, NULL, NULL, NULL, NULL);
      if (isenthalpic) {
        (*solids[i].smix)(FIRST, currentT, currentP, rSol, &pTemp, NULL, NULL);
        pTotal -= mTotal*pTemp;
      } else if (isentropic) {
        (*solids[i].hmix)(FIRST, currentT, currentP, rSol, &pTemp);
        pTotal += mTotal*pTemp;
      } else if (isochoric) {
        (*solids[i].gmix)(FIRST, currentT, currentP, rSol, &pTemp, NULL, NULL, NULL);
        pTotal += mTotal*pTemp;
        (*solids[i].vmix)(FIRST, currentT, currentP, rSol, &pTemp, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
        pTotal -= mTotal*currentP*pTemp;
      } else {
        (*solids[i].gmix)(FIRST, currentT, currentP, rSol, &pTemp, NULL, NULL, NULL);
        pTotal += mTotal*pTemp;
      }
    }
  }

  /* Add term for oxygen if path is buffered */
  if (silminState->fo2Path != FO2_NONE) {
    mO2T = mO2S;
    if (hasLiquid) for (nl=0; nl<silminState->nLiquidCoexist; nl++) mO2T += mO2L[nl];
    if (isenthalpic) pTotal -= mO2T*muO2E/currentT;
#ifdef PHMELTS_ADJUSTMENTS
    else if (silminState->fo2Alt)
                     pTotal -= mO2T*muO2E;
#endif
    else             pTotal -= mO2T*muO2L[0];                                          /* V-fO2 fix was: mO2T*muO2E */
  }
#ifdef PHMELTS_ADJUSTMENTS
  /* Add term for H2O if path is buffered */
  if (H2Obuffer) pTotal -= mH2O*muH2O;
#endif

#ifdef DEBUG
  printf("linearSearch: pTotal = %.*g, lambda = %.*g\n", DBL_DIG, pTotal, DBL_DIG, lambda);
  if (silminState->fo2Path != FO2_NONE) {
    printf("  moles of O2 in system(-) = %.*g (mu = %.*g)\n", DBL_DIG, mO2T, DBL_DIG, muO2E);
    printf("  moles of O2 in solids(-) = %.*g (mu = %.*g)\n", DBL_DIG, mO2S, DBL_DIG, muO2S);
    for (nl=0; nl<silminState->nLiquidCoexist; nl++)
    printf("  moles of O2 in liquid(%d) = %.*g (mu = %.*g)\n", nl, DBL_DIG, mO2L[nl], DBL_DIG, muO2L[nl]);
    printf("  delta moles of O2 = %.*g\n", DBL_DIG, mO2T-silminState->oxygen);
    if (isenthalpic) printf("  delta H of O2 = %.*g\n", DBL_DIG, (mO2T-silminState->oxygen)*(oxygen.cur).h);
    if (isentropic)  printf("  delta S of O2 = %.*g\n", DBL_DIG, (mO2T-silminState->oxygen)*(oxygen.cur).s);
    if (isochoric)   printf("  delta V of O2 = %.*g\n", DBL_DIG, (mO2T-silminState->oxygen)*(oxygen.cur).v);
  }
#endif

  return pTotal;
}

/* end of file LINEAR_SEARCH.C */
