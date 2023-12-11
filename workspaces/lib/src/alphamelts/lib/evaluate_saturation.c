const char *evaluate_saturation_ver(void) { return "$Id: evaluate_saturation.c,v 1.4 2009/04/16 16:35:23 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: evaluate_saturation.c,v $
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
MELTS Source Code: RCS Revision 1.2  2002/10/12 16:41:29  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:58:12  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
 * Revision 3.8  1997/06/21  22:49:57  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.7  1997/05/03  20:23:36  ghiorso
 * *** empty log message ***
 *
 * Revision 3.6  1997/03/27  17:03:39  ghiorso
 * *** empty log message ***
 *
 * Revision 3.5  1996/09/24  20:33:42  ghiorso
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
 * Revision 3.1  1995/08/18  17:30:53  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Function to evaluate the saturation state with respect to solid phases
**      returns TRUE  if a supersaturated solid has been detected
**           FALSE if no solids are supersaturated
**      In either event the first npc elements of the double array rSol
**      are returned. The elements correspond to:
**        rSol[i]   = chemical affinity for the ith solid if
**                    solids[i].type == PHASE, else
**                  = estimated comosition for the ith solid in terms of
**                    the independent compositional variables if
**                    solids[i].type == COMPONENT. (Note only the first
**                    solids[i].nr elements are used to store the composition,
**                    the solids[i].na location is undefined).
**      (file: EVALUATE_SATURATION.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  September 18, 1991
**              Extracted from file SILMIN_SUPPORT.C
**      V1.1-1  Mark S. Ghiorso  September 23, 1991
**              Changed logic to allow for estimation of saturation state
**              when a solid component is absent from the liquid phase
**              Altered call to getAffinityAndComposition
**      V1.1-2  Mark S. Ghiorso  September 28, 1991
**              Added summary printout prior to return (debugging purposes)
**      V1.1-3  Mark S. Ghiorso  October 15, 1991
**              (1) Modiified references to silminState->solidComp array
**                  to reflect immiscible solids structures
**      V1.1-4  Mark S. Ghiorso  Noovember 2, 1991
**              Added a temorpary code fragment to prevent precipitation of
**              non-quadrilateral pyroxene components
**      V2.0-1  Mark S. Ghiorso  November 14, 1991
**              Conversion to OSF Motif V1.1.1/X11 Release 4
**      V2.0-2  Mark S. Ghiorso  May 4, 1992
**              Added cycling flag to computation to prevent infinite
**              loops in constrained path crystallization
**      V2.0-3  Mark S. Ghiorso  July 1, 1992
**              Removed non-quad pyroxene provision
**      V2.0-4  Mark S. Ghiorso  July 3, 1992
**              Removed Na from pyroxene
**      V2.0-5  Mark S. Ghiorso  July 20, 1992
**              Revoked 2.0-4
**      V2.1-1  Mark S. Ghiorso  June 15, 1993
**              Added special case for failure of pyroxene. Second try attempts
**              to detect saturation for a purely quadrilateral pyroxene
**      V3.0-1  Paul D. Asimow  March 17, 1995
**              If liquid is absent, construct fictive liquid component
**              chemical potentials by assuming solid assemblage is in equilibrium
**      V3.1-1  Paul D. Asimow March 22, 1995
**              If liquid is absent, attempt to identify solidus by saturation test
**      V3.1-2  Paul D. Asimow March 29, 1995
**              silminState->incSolids[npc] is now used to store whether liquid is allowed
**      V3.1-3  Paul D. Asimow August 7, 1995
**              compatibility changes for subsolidus buffering
**--
*/

#include <stdlib.h>
#include <stdio.h>

#include "melts_gsl.h"            /* Needs to go before silmin.h */
#include "silmin.h"               /*SILMIN structures include file          */

int getAffinityAndCompositionGeneric(double t, double p, int index,
    int *zeroX, double *muMinusMu0, double *affinity, double *indepVar);

int getAffinityAndCompositionPyroxene(double t, double p, int index,
    int *zeroX, double *muMinusMu0, double *affinity, double *indepVar);

#ifdef DEBUG
#undef DEBUG
#endif

#define DO_PYROXENE_COMPROMISE

#define CYL 2                 /* Limit on the number of phase in/out cycles */

int evaluateSaturationState(double *rSol, double *rLiq)
{
    static int *zeroX;
    static double *muSol, *xSol, *muLiq, *muTemp, *liquidComp;
    int i, j, k, l = 0, hasSupersat;
    double t, p;
#ifdef PHMELTS_ADJUSTMENTS
    int iLiqH2O, iBulkH2O, iPhaseH2O;
#endif

    if (muLiq == NULL) {
        muLiq      = (double *) malloc((unsigned) nlc*sizeof(double));
        muSol      = (double *) malloc((unsigned) nlc*sizeof(double));
        xSol       = (double *) malloc((unsigned) nlc*sizeof(double));
        zeroX      = (int *)    malloc((unsigned) nlc*sizeof(int));
        muTemp     = (double *) calloc((unsigned) nlc, sizeof(double));
        liquidComp = (double *) calloc((unsigned) nlc, sizeof(double));
    }

#ifdef PHMELTS_ADJUSTMENTS
    for (iLiqH2O=0; iLiqH2O<nlc; iLiqH2O++) if (!strcmp(liquid[iLiqH2O].label, "H2O")) break;
    for (iBulkH2O=0; iBulkH2O<nc; iBulkH2O++) if (!strcmp(bulkSystem[iBulkH2O].label, "H2O")) break;
    if (calculationMode == MODE__MELTSandCO2_H2O || calculationMode == MODE__MELTSandCO2) {
        for (iPhaseH2O=0; iPhaseH2O<npc; iPhaseH2O++) if (!strcmp(solids[iPhaseH2O].label, "h2oduan")) break;
    }
    else {
        for (iPhaseH2O=0; iPhaseH2O<npc; iPhaseH2O++) if (!strcmp(solids[iPhaseH2O].label, "fluid")) break;
    }
#endif

    t = silminState->T;
    p = silminState->P;

    if (silminState->liquidMass != 0.0) {

        /* obtain liquid chemical potentials - use the first liquid */
        conLiq(SECOND, THIRD, t, p, (double *) NULL, silminState->liquidComp[0], rLiq,
            (double *) NULL, (double **) NULL, (double ***) NULL,  &(silminState->fo2));
        actLiq(SECOND, t, p, rLiq, NULL, muLiq, NULL, NULL);
        for (i=0; i<nlc; i++) if ((silminState->liquidComp)[0][i] != 0.0) muLiq[i] += (liquid[i].cur).g;

        // Correction to Fe2O3 chemical potential due to the fact that when the fO2 buffer path is constrained,
        // homogeneous redox equilibrium in the melt is not consistent with Kress and Carmichael.
        // This results in a saturation state convergence to an incorrect composition/affinity for phases that
        // contain ferric iron if the Fe2O3 chemical potential is not adjusted for internal consistency with K&C.

        if ( (silminState->fo2Path != FO2_NONE) && (
#ifndef TRUE_xMELTS
                    (calculationMode == MODE_xMELTS) ||
#endif
                    (calculationMode == MODE__MELTS) ||
                    (calculationMode == MODE__MELTSandCO2) ||
                    (calculationMode == MODE__MELTSandCO2_H2O)
                )) {
            // This is the chemical potential of oxygen calculated from the imposed buffer
            double muO2 = R*(silminState->T)*(silminState->fo2)*log(10.0) + oxygen.cur.g;
            // This is the free energy change of the reaction: Fe2SiO4 + 1/2 O2 = Fe2O3 + SiO2 (MELTS component choice)
            // If Kress & Carmichael were consistent with homogeneous equilibrium in MELTS, then this number would be zero
            double deltaG = muLiq[0] + muLiq[3] - muLiq[5] - muO2/2.0;
#ifdef DEBUG
            printf("The chemical potential of oxygen is %g, deltaG is %g\n", muO2, deltaG);
#endif
            // adjust the chemical potential of Fe2O3 so that the fO2 constraint will propagate
            // into the correct selection of saturation phases
            muLiq[3] -= deltaG;

        } else if ( (silminState->fo2Path != FO2_NONE) && (calculationMode == MODE_pMELTS) ) {
            // This is the chemical potential of oxygen calculated from the imposed buffer
            double muO2 = R*(silminState->T)*(silminState->fo2)*log(10.0) + oxygen.cur.g;
            // This is the free energy change of the reaction: Fe2SiO4 + 1/2 O2 = Fe2O3 + 1/4 Si4O8 (pMELTS component choice)
            // If Kress & Carmichael were consistent with homogeneous equilibrium in MELTS, then this number would be zero
            double deltaG = muLiq[0]/4.0 + muLiq[3] - muLiq[5] - muO2/2.0;
#ifdef DEBUG
            printf("The chemical potential of oxygen is %g, deltaG is %g\n", muO2, deltaG);
#endif
            // adjust the chemical potential of Fe2O3 so that the fO2 constraint will propagate
            // into the correct selection of saturation phases
            muLiq[3] -= deltaG;
        }

#ifdef PHMELTS_ADJUSTMENTS
        /* this is a convenient place to set aH2O if not buffered */
        if (silminState->H2Obuffer == FALSE) {
            if (silminState->liquidComp[0][iLiqH2O] != 0.0) silminState->aH2O = exp((muLiq[iLiqH2O] - (solids[iPhaseH2O].cur).g)/(R*silminState->T));
            else silminState->aH2O = 0.0;
        }
#endif

        for (i=0;i<nlc;i++) liquidComp[i] = silminState->liquidComp[0][i];

    } else {    /* liquid is absent */

        gsl_vector *muAllSol, *muAllLiq, *S;
        gsl_matrix *stoichMatrix, *V;
        gsl_matrix_view A;
        gsl_vector_view b, x;
        int nonZero, m, n;

        muAllSol = gsl_vector_alloc((size_t) npc);
        muAllLiq = gsl_vector_alloc((size_t) nlc+1);
        stoichMatrix = gsl_matrix_alloc((size_t) npc, (size_t) nlc+1);

        /* obtain solid chemical potentials */
        for (i=0,j=0;i<npc;i++) {
            if (solids[i].type == PHASE) {
                if (silminState->nSolidCoexist[i]) {
                    if (solids[i].na == 1) {
                        gsl_vector_set(muAllSol, j, (solids[i].cur).g);
                        j++;
                    }
                    else {
                        for (k=0;k<solids[i].na;k++) xSol[k] = silminState->solidComp[i+1+k][0];
                        (*solids[i].convert)(SECOND, THIRD, t, p, NULL, xSol, rSol, NULL, NULL, NULL, NULL, NULL);
                        (*solids[i].activity)(SECOND, t, p, rSol, NULL, muTemp, NULL);
                        for (k=0;k<solids[i].na;k++) {
                            if (silminState->solidComp[i+1+k][0] != 0.0) {
                                gsl_vector_set(muAllSol, j, (solids[i+1+k].cur).g + muTemp[k]);
                                j++;
                            }
                        }
                    }
                }
            }
        }
        /* obtain stoichiometry matrix */
        for (i=0,j=0;i<npc;i++) {
            if (solids[i].type == PHASE) {
                if (silminState->nSolidCoexist[i]) {
                    if (solids[i].na == 1) {
                        for (k=0,l=0;k<nlc;k++) {
                            for (n=0,nonZero=1;n<nc;n++) if (silminState->bulkComp[n] == 0.0 && liquid[k].liqToOx[n] != 0) nonZero = 0;
                            if (nonZero) {
                                gsl_matrix_set(stoichMatrix, j, l,  solids[i].solToLiq[k]);
                                l++;
                            }
                        }
                        j++;
                    } else {
                        for (m=0;m<solids[i].na;m++) {
                            if (silminState->solidComp[i+1+m][0] != 0.0) {
                                for (k=0,l=0;k<nlc;k++) {
                                    for (n=0,nonZero=1;n<nc;n++) if (silminState->bulkComp[n] == 0.0 && liquid[k].liqToOx[n] != 0) nonZero = 0;
                                    if (nonZero) {
                                        gsl_matrix_set(stoichMatrix, j, l, solids[i+1+m].solToLiq[k]);
                                        l++;
                                    }
                                }
                                j++;
                            }
                        }
                    }
                }
            }
        }
        /* l now holds number of columns, j number of rows */
        /* obtain liquid chemical potentials by least squares; if solids are in
           equilibrium the system should be overdetermined but exactly consistent */
        S = gsl_vector_alloc((size_t) l);
        V = gsl_matrix_alloc((size_t) l, (size_t) l);
        A = gsl_matrix_submatrix(stoichMatrix, (size_t) 0, (size_t) 0, (size_t) j, (size_t) l);
        b = gsl_vector_subvector(muAllSol, (size_t) 0, (size_t) j);
        x = gsl_vector_subvector(muAllLiq, (size_t) 0, (size_t) l);

        gsl_linalg_SV_decomp(&A.matrix, V, S, &x.vector); // muAllLiq used as work space
        for (i=0;i<l;i++) if (gsl_vector_get(S, i) < 1.0e-08) gsl_vector_set(S, i, 0.0);
        gsl_linalg_SV_solve(&A.matrix, V, S, &b.vector, &x.vector);

        for (i=0,l=0;i<nlc;i++) {
            for (n=0,nonZero=1;n<nc;n++) if (silminState->bulkComp[n] == 0.0 && liquid[i].liqToOx[n] != 0) nonZero = 0;
            if (!nonZero) {
                muLiq[i] = 0.0;
                liquidComp[i] = 0.0;
            } else {
                muLiq[i] = gsl_vector_get(&x.vector, l);
                liquidComp[i] = 1.0;
                l++;
            }
        }

        gsl_matrix_free(V); gsl_vector_free(S);
        gsl_matrix_free(stoichMatrix);
        gsl_vector_free(muAllLiq); gsl_vector_free(muAllSol);

#ifdef PHMELTS_ADJUSTMENTS
        /* this is a convenient place to set aH2O if not buffered */
        if (silminState->H2Obuffer == FALSE) { // only get here if pMELTS or rhyolite-MELTS 1.0.2
            if (silminState->bulkComp[iBulkH2O] != 0.0) silminState->aH2O = exp((muLiq[iLiqH2O] - (solids[iPhaseH2O].cur).g)/(R*silminState->T));
            else silminState->aH2O = 0.0;
        }
        /* allow saturation with hydrous phases other than H2O even if H2O is absent from bulk composition, if aH2O buffered */
        if (silminState->H2Obuffer && (silminState->aH2O != 0.0)) {
            muLiq[iLiqH2O] = (solids[iPhaseH2O].cur).g + R*silminState->T*log(silminState->aH2O);
            liquidComp[iLiqH2O] = 1.0;
        }
#endif
    }

    hasSupersat = FALSE;

    /* obtain solid chemical potentials, affinities and composition estimates  */
    for (i=0, j=0; i<npc; i++) {
        rSol[i] = 0.0;
        if (solids[i].type == PHASE) {
            if ((silminState->incSolids)[j] && ((silminState->cylSolids)[i] < CYL) &&
                    (silminState->solidComp)[i][0] == 0.0) {
                if (solids[i].na == 1) {

                    muSol[0] = (solids[i].cur).g;
                    for (k=0; k<nlc; k++) {
                        if( (solids[i].solToLiq)[k] != 0.0) {
                            if (liquidComp[k] != 0.0)
                                muSol[0] -= (solids[i].solToLiq)[k] * muLiq[k];
                            else {
                                muSol[0] = 0.0;
                                break;
                            }
                        }
                    }
                    rSol[i] = muSol[0]; /* Affinity is < 0 if phase is supersaturated */
                    hasSupersat |= (rSol[i] < 0.0);

                } else if (solids[i].na > 1) {

                    for (k=0; k<solids[i].na; k++) {
                        muSol[k] = (solids[i+1+k].cur).g;
                        zeroX[k] = FALSE;
                        for (l=0; l<nlc; l++) {
                            if( (solids[i+1+k].solToLiq)[l] != 0.0) {
                                if (liquidComp[l] != 0.0)
                                    muSol[k] -= (solids[i+1+k].solToLiq)[l] * muLiq[l];
                                else {
                                    muSol[k] = 0.0;
                                    zeroX[k] = TRUE;
                                    break;
                                }
                            }
                        }
                    }

                    /* Affinity is returned in rSol[i]. It is < 0 if phase is
                    supersaturated. The composition in terms of independent
                    variables (ie. solids[i].nr of them) is returned in
                    rSol[i+1] to rSol[i+1+solids[i].nr]. It may be converted
                    subsequently to moles of endmembers if necessary                 */
#ifdef RHYOLITE_ADJUSTMENTS
#ifndef TESTDYNAMICLIB
                	if (!strcmp(solids[i].label, "alkali-feldspar")) {
#ifdef DEBUG
                        printf("Using generic method for alkali-feldspar.\n");
#endif
                        getAffinityAndCompositionGeneric(t, p, i, zeroX, muSol, &rSol[i], &rSol[i+1]);
                    } else if (!strcmp(solids[i].label, "plagioclase")) {
#ifdef DEBUG
                        printf("Using generic method for plagioclase.\n");
#endif
                	    getAffinityAndCompositionGeneric(t, p, i, zeroX, muSol, &rSol[i], &rSol[i+1]);
                    }
#else
                	if (!strcmp(solids[i].label, "feldspar")) {
#ifdef DEBUG
                        printf("Using generic method for feldspar.\n");
#endif
	                    getAffinityAndCompositionGeneric(t, p, i, zeroX, muSol, &rSol[i], &rSol[i+1]);
                    }
#endif /* TESTDYNAMICLIB */
        	        else if (!strcmp(solids[i].label, "orthopyroxene")) {
#ifdef DEBUG
                        printf("Using speciation method for orthopyroxene.\n");
#endif
                        getAffinityAndCompositionPyroxene(t, p, i, zeroX, muSol, &rSol[i], &rSol[i+1]);

                    } else if (!strcmp(solids[i].label, "clinopyroxene")) {
#ifdef DEBUG
                        printf("Using speciation method for clinopyroxene.\n");
#endif
                	    getAffinityAndCompositionPyroxene(t, p, i, zeroX, muSol, &rSol[i], &rSol[i+1]);
/*
            	    } else if (!strcmp(solids[i].label, "spinel")) {
#ifdef DEBUG
                        printf("Using speciation method for spinel.\n");
#endif
	                    getAffinityAndCompositionSpinel(t, p, i, zeroX, muSol, &rSol[i], &rSol[i+1]);
*/
                    } else
#endif /* RHYOLITE_ADJUSTMENTS */
        	        if (!getAffinityAndComposition(t, p, i, zeroX, muSol, &rSol[i], &rSol[i+1])) {
                        if (!strcmp(solids[i].label, "clinopyroxene") ||
                                !strcmp(solids[i].label, "orthopyroxene")    ) {
#ifdef DO_PYROXENE_COMPROMISE
                	        if (muSol[2] != 0.0) {
                                int tempZeroX[8];
                                for (k=0; k<solids[i].na; k++) tempZeroX[k] = zeroX[k];
                                muSol[3] = 0.0; tempZeroX[3] = TRUE;  /* Ca(Ti,Mg)(Al,Si)2O6   */
                                muSol[4] = 0.0; tempZeroX[4] = TRUE;  /* Ca(Ti,Mg)(Fe3+,Si)2O6 */
                                muSol[5] = 0.0; tempZeroX[5] = TRUE;  /* CaFeAlSiO6            */
                                muSol[6] = 0.0; tempZeroX[6] = TRUE;  /* NaAlSi2O6             */
#ifdef DEBUG
                                printf("Making the PYROXENE compromise in evaluateSaturationState!\n");
#endif
                                if (getAffinityAndComposition(t, p, i, tempZeroX, muSol, &rSol[i],
                                        &rSol[i+1])) {
                                    if (tempZeroX[3] && !zeroX[3]) rSol[i+2] = 0.0001;
                                    if (tempZeroX[4] && !zeroX[4]) rSol[i+3] = 0.0001;
                                    if (tempZeroX[5] && !zeroX[5]) rSol[i+4] = 0.0001;
                                    if (tempZeroX[6] && !zeroX[6]) rSol[i+5] = 0.0001;
                                } else for (k=0; k<=solids[i].na; k++) rSol[i+k] = 0.0;
                            } else {
#ifdef DEBUG
                                printf("Making the Mg-binary PYROXENE compromise in evaluateSaturationState!\n");
#endif
                                rSol[i]   = -100000.0;
                                rSol[i+1] = 0.0;
                                rSol[i+2] = 0.0;
                                rSol[i+3] = 0.0;
                                rSol[i+4] = 0.0;
                                rSol[i+5] = 0.0;
                                rSol[i+6] = (!strcmp(solids[i].label, "clinopyroxene")) ? 0.5 : 0.5;
	      }
#endif /* DO_PYROXENE_COMPROMISE */
                        } else if (!strcmp(solids[i].label, "nepheline")) {
                            int tempZeroX[4];
                            for (k=0; k<solids[i].na; k++) tempZeroX[k] = zeroX[k];
                            /* muSol[2] = 0.0; tempZeroX[2] = TRUE; */ /* vc-nepheline */
                            /* muSol[3] = 0.0; tempZeroX[3] = TRUE; */ /* ca-nepheline */
#ifdef DEBUG
                            printf("Making the NEPHELINE compromise in evaluateSaturationState!\n");
#endif
                            if (getAffinityAndComposition(t, p, i, tempZeroX, muSol, &rSol[i],
                                    &rSol[i+1])) {
                                /* if (tempZeroX[2] && !zeroX[2]) rSol[i+2] = 0.0001; */
                                /* if (tempZeroX[3] && !zeroX[3]) rSol[i+3] = 0.0001; */
                            } else for (k=0; k<=solids[i].na; k++) rSol[i+k] = 0.0;

                        } else if (!strcmp(solids[i].label, "kalsilite")) {
                            int tempZeroX[4];
                            for (k=0; k<solids[i].na; k++) tempZeroX[k] = zeroX[k];
                            /* muSol[2] = 0.0; tempZeroX[2] = TRUE; */ /* vc-nepheline */
                            /* muSol[3] = 0.0; tempZeroX[3] = TRUE; */ /* ca-nepheline */
#ifdef DEBUG
                            printf("Making the KALSILITE compromise in evaluateSaturationState!\n");
#endif
                            if (getAffinityAndComposition(t, p, i, tempZeroX, muSol, &rSol[i],
                   &rSol[i+1])) {
                                /* if (tempZeroX[2] && !zeroX[2]) rSol[i+2] = 0.0001; */
                                /* if (tempZeroX[3] && !zeroX[3]) rSol[i+3] = 0.0001; */
                            } else for (k=0; k<=solids[i].na; k++) rSol[i+k] = 0.0;

                        } else for (k=0; k<=solids[i].na; k++) rSol[i+k] = 0.0;
                    }  /* end of if-block for getAffinityAndComposition() */

                    hasSupersat |= (rSol[i] < 0.0);
                    i += solids[i].na;
                }
            }
            j++;
        }
    }

#ifdef DEBUG
    printf("Results of call to evaluateSaturationState:\n");
#endif
    /* obtain liquid chemical potentials, affinity and composition estimate
     silminState->incSolids[npc] is now used to store whether liquid is allowed  */
    rLiq[nlc-1] = 0.0;
#ifdef PHMELTS_ADJUSTMENTS
    if (silminState->liquidMass == 0.0 && silminState->incLiquids) {
#else
    if (silminState->liquidMass == 0.0 && silminState->incSolids[npc]) {
#endif
        for (k=0; k<nlc; k++) {
            muSol[k] = (liquid[k].cur).g;
            zeroX[k] = FALSE;
            if (liquidComp[k] != 0.0) {
                muSol[k] -= muLiq[k];
            } else {
                muSol[k] = 0.0;
                zeroX[k] = TRUE;
            }
        }

        /* Affinity is returned in rLiq[nlc-1]. It is < 0 if phase is
       supersaturated. The composition in terms of independent
       variables (ie. NR of them) is returned in
       rLiq[0] to rLiq[NR-1]. It may be converted
       subsequently to moles of endmembers if necessary                 */

        if (!getAffinityAndComposition(t, p, -1, zeroX, muSol, &rLiq[nlc-1], rLiq))
            for (k=0; k<(nlc-1); k++) rLiq[k] = 0.0;
        hasSupersat |= (rLiq[nlc-1] < 0.0);

#ifdef DEBUG
        printf("%-10.10s A = %13.6g\n", "melt", rLiq[nlc-1]);
        for (j=0; j<(nlc-1); j++) {
            printf("  r[%2d] = %7.4f", j, rLiq[j]);
            if (j == 3 || j == 7 || j == 11 || j == 15) printf("\n");
        }
        printf("\n");
#endif
    }

#ifdef DEBUG
    for (i=0; i<npc; i++) {
        if(solids[i].type == PHASE && rSol[i] != 0.0) {
            printf("%-10.10s A = %13.6g", solids[i].label, rSol[i]);
            if(solids[i].na > 1) for (j=0; j<solids[i].nr; j++)
                printf("  r[%1d] = %7.4f", j, rSol[i+1+j]);
            printf("\n");
        }
    }
#endif

    return hasSupersat;
}

/* end of file EVALUATE_SATURATION.C */
