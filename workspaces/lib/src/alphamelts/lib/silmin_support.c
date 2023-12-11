const char *silmin_support_ver(void) { return "$Id: silmin_support.c,v 1.8 2009/04/24 20:51:09 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: silmin_support.c,v $
MELTS Source Code: RCS Revision 1.6  2008/03/06 17:51:23  ghiorso
MELTS Source Code: RCS New fluid fractionation mode and other enhancements.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.5  2007/08/23 16:09:40  ghiorso
MELTS Source Code: RCS Database updates.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.4  2007/06/08 17:25:43  ghiorso
MELTS Source Code: RCS Added code to allow regression of Ghiorso EOS parameters
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
MELTS Source Code: RCS Revision 1.4  2003/05/03 18:43:56  ghiorso
MELTS Source Code: RCS *** empty log message ***
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
 * Revision 3.12  1997/06/21  22:49:26  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.11  1997/05/03  20:23:06  ghiorso
 * *** empty log message ***
 *
 * Revision 3.10  1997/03/27  17:03:10  ghiorso
 * *** empty log message ***
 *
 * Revision 3.9  1996/09/24  20:33:20  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 3.8  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 3.7  1995/11/23  22:37:42  ghiorso
 * Final implementation of subsolidus fO2 buffering.
 *
 * Revision 3.6  1995/11/01  22:40:27  ghiorso
 * Implementation of subsolidus options after Asimow.
 * Additional implementation of nepheline solid solutions.
 *
 * Revision 3.5  1995/09/04  23:30:20  ghiorso
 * Added new function copySilminStateStructure, for copy the silminState
 * structure. Called in create_xy_plot.padb.c for allocating history records.
 *
 * Revision 3.4  1995/09/04  20:01:28  ghiorso
 * Update to allow display of bulk composition (in grams) in the text entry
 * fields of the main silmin display. Liquid composition is no longer
 * display here, and is available only through the popup selection.
 *
 * Revision 3.3  1995/09/01  23:53:03  ghiorso
 * Modifications made to update interface for V3.x and consolidate
 * Graph Widgets
 *
 * Revision 3.2  1995/08/31  00:25:55  ghiorso
 * Removed References to multiple user plots
 *
 * Revision 3.1  1995/08/18  19:14:22  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Collection of support functions for crystallization calculations
**
**      allocSilminStatePointer()
**      Function to allocate space for the SilminState structure
**
**      checkStateAgainstInterface()
**      Function to evaluate state of interface display contained in the
**      INTERFACE.H structures:
**        (1) state of toggle button widgets:
**            tg_path_none, tg_path_hm, tg_path_nno, tg_path_fmq, tg_path_iw,
**            tg_isenthalpic, tg_isochoric, tg_fractionation,
**        (2) compositionValues
**        (3) debugEntries
**        (4) tpValues
**        (5) assimilantValues and assimilantUnits
**        (6) magmaValues
**        (7) includedSolids
**      and compare it to the state of the system contained in the SILMIN.H
**      structures:
**        (1) silminState.
**      returning a flag (see global macros (defined below) indicating what has
**      been done.
**
**      update*()
**      Functions to update the various display elements, including the
**      text widgets, label widgets and graphs. The "member" argument
**      corresponds to a partricular entry to update. These are defined
**      by macros in INTERFACE.H. Passing a negative argument will
**      update all entries.
**
**      (file: SILMIN_SUPPORT.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  September 5, 1991
**              function: checkStateAgainstInterface()
**                preliminary version finished September 10, 1991
**              function: evaluateSaturationState()
**                preliminary version finished September 11, 1991
**              function getEqualityConstraints()
**                preliminary version finished September 12, 1991
**              function getProjGradientAndHessian()
**                altered call to liqCon on September 13, 1991 as
**                  part of code development for the function
**                create code to transform intensive to extensive
**                  derivatives of the thermodynamic potential
**                altered call to (*solids[].convert) on September 14, 1991
**                  as part of code development for the function
**                preliminary version finished September 17, 1991
**              function linearSearch()
**                preliminary version finished September 18, 1991
**      V1.0-2  Mark S. Ghiorso  September 18, 1991
**              Removed major numerical routines, left those pertaining
**                to interface update and display
**      V1.0-3  Mark S. Ghiorso  September 19, 1991
**              Began work on update*() functions for interface display
**                               September 20, 1991
**              continues ...
**      V1.0-4  Mark S. Ghiorso  September 24, 1991
**              Altered parameter list to (*solids[].convert)
**      V1.0-5  Mark S. Ghiorso  October 15, 1991
**              (1) Modifications to references for arrays
**                  silminState->solidComp and silminState->solidDelta
**                  to allow for immiscible solid phases
**              (2) Initialized silminState->nSolidCoexist to count number
**                  of solid phases
**              (3) The module updateSolidGW needs to be rewritten in light
**                  of these and other functionality considerations
**              (4) Complete reorganization of updateSolidGW logic (includes
**                  extensions for immiscible solid phases)
**      V1.0-6  Mark S. Ghiorso  October 18, 1991
**              (1) Added graph legend code
**      V2.0-1  Mark S. Ghiorso  November 14, 1991
**              Conversion to OSF Motif V1.1.1/X11 Release 4
**      V2.0-2  Mark S. Ghiorso  November 23, 1991
**              (1) Removed references to debug entries on interface
**                  display and menu bar
**      V2.0-3  Mark S. Ghiorso  November 27, 1991
**              (1) Complete rewrite of updateSolidADB to utilize the
**                  Vlist widget structure and output affinities
**      V2.0-4  Mark S. Ghiorso  December 10, 1991
**              (1) Reorganized assimilantValues references and removed
**                  magmaValues and magma_padb references
**      V2.0-5  Mark S. Ghiorso  December 14, 1991
**              (1) Revised assimilation algorithm and eliminated magma
**                  mixing code
**      V2.0-6  Mark S. Ghiorso  December 21, 1991
**              (1) Inserted code from silmin (SILMIN.C) to update
**                  silminState->liquidComp when display composition of
**                  liquid is altered
**      V2.0-7  Mark S. Ghiorso  January 2, 1992
**              (1) Implemented changes for fractionation mode
**      V2.0-8  Mark S. Ghiorso  January 4, 1992
**              (1) Altered resource declarations for graph widget from
**                  XtN -> GwN
**              (2) Altered calls to polyline() and polymarker() to
**                  GwPolyline() and GwPolymarker()
**              (3) Added calls to GwStartBatchUpdate() and GwEndBatchUpdate()
**                  surrounding the liquid and solid graph poly* calls
**              (4) Added support for display on user_graph (AFM only)
**      V2.0-9  Mark S. Ghiorso  January 6, 1992
**              (1) Corrected bound constraint algorithm in userGraph by
**                  eliminating multiple "identical" points connected by
**                  polylines.
**      V2.0-10 Mark S. Ghiorso  January 10, 1992
**              (1) Added wt % assimilation support
**      V2.0-11 Mark S. Ghiorso  January 15, 1992
**              (1) removed definition of updateAssimilantPADB()
**                  (now defined in create_assimilant_padb.c - we needed
**                  locally static widget names)
**      V2.0-12 Mark S. Ghiorso  February 19, 1992
**              (1) Corrected minor casting violations for ANSI C compliance
**              (2) Removed global dependence on arg_set
**      V3.0-1  Mark S. Ghiorso  April 27, 1992
**              (1) Begin nodifications for f O2 buffering of reaction path
**      V3.0-2  Mark S. Ghiorso  May 1, 1992
**              (1) Added correction to silminState->oxygen value when
**                  user modifies bulk composition
**      V3.0-3  Mark S. Ghiorso  May 4, 1992
**              (1) Added initialization of space for cylSolids in
**                  allocSilminStatePointer()
**      V3.1-1  Mark S. Ghiorso  July 13, 1992
**              (1) Added a nrw column (mineral formulas) to the phases
**                  vList widget
**      V3.1-2  Mark S. Ghiorso  September 29, 1992
**              Converted TextField to Text widgets as a bug workaround
**              for DECWindows Motif V 1.1
**      V3.1-3  Mark S. Ghiorso  June 14, 1993
**              Added code for new fo2 path constraints
**      V3.1-4  Mark S. Ghiorso  September 21, 1993
**              XtFree -> XmStringFree
**      V3.1-5  Mark S. Ghiorso  September 29, 1993
**              Modified call to realloc to catch zero pointer (SPARC port)
**      V3.1-6  Mark S. Ghiorso  April 5, 1994
**              Added #ifdef __osf__ to correct 64 bit errors
**      V4.0-1  Mark S. Ghiorso  May 11, 1994
**              (1) Modified calls to *vmix to reflect additional derivatives
**              (2) Removed reference to grove, walker plots
**              (3) Added reference to isentropic constraints
**              (4) Removed warning dialogs if isenthalpic, isentropic or
**                  isochoric constraints are switched on
**                               June 9, 1994
**              (5) Added section to compute thermodynamic properties of
**                  assimilant and store them in silminState->assimTD
**                               June 11, 1994
**              (6) Fixed logic error in computing assimTD structure
**                               June 13, 1994
**              (7) Added three new functions to adjust T and P for "update"
**                  changes in enthalpy, entropy or volume
**                  void correctTforChangeInEnthalpy(void)
**                  void correctTforChangeInEntropy(void)
**                  void correctPforChangeInVolume(void)
**                               June 15, 1994
**              (8) Added interface checks for enthalpy, entropy and
**                  volume terms in tp_values[]
**                               June 16, 1994
**              (9) Corrected error in 4.0-1.7 routines on abnormal
**                  termination for negative residuals
**                               July 2, 1994
**              (10) Altered convergence criteria for 4.0-1.7 routines
**                   (scaled test on residuals)
**      V4.1-1  Mark S. Ghiorso
**              (1) Experiments with convergence in correctTforChangeInEntropy
**      V4.1-2  Mark S. Ghiorso  February 27, 1995
**              (1) Altered convergence criteria in correctTforChangeInEnthalpy,
**                  correctTforChangeInEntropy, correctPforChangeInVolume
**      V4.1-3  Mark S. Ghiorso  March 25, 1995
**              (1) Disallowed constraints on both fo2 and enthalpy or
**                  entropy or volume
**      V5.1-1  Paul D. Asimow  April 26, 1995
**              (1) Increase storage size of incSolids to allow for liquid
**      V5.2-1  Paul D. Asimow  August 2, 1995
**              (1) Enable subsolidus buffering -- changes to correctTforChange
**              InEnthalpy, correctTforChangeInEntropy, and correctPforChangeIn
**              Volume.
**              (2) New functions copyStateInfo and copyThermoData
**              (3) New function subsolidusmuO2 for subsolidus buffering
**              (4) New function addOrDropLiquid supports liquid-absent
**                  phase drops and adds
**--
*/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "melts_gsl.h"
#include "silmin.h"                 /*SILMIN structures include file        */

#ifndef BATCH_VERSION
#include <Xm/PushBG.h>
#include <Xm/RowColumn.h>
#include <Xm/Text.h>
#include <Xm/ToggleBG.h>

#include "interface.h"              /*Specific external declarations        */
#include "vframe.h"                 /*vframe widget include file            */
#include "vheader.h"                /*vheader widget include file           */
#include "vlist.h"                  /*vlist widget include file             */
#else
#define True '\001'
#endif

#ifdef DEBUG
#undef DEBUG
#endif

#define REALLOC(x, y) (((x) == NULL) ? malloc(y) : realloc((x), (y)))

static void multiplyThermoData(ThermoData *target, double factor)
{
    target->g *= factor;
    target->h *= factor;
    target->s *= factor;
    target->v *= factor;
    target->cp *= factor;
    target->dcpdt *= factor;
    target->dvdt *= factor;
    target->dvdp *= factor;
    target->d2vdt2 *= factor;
    target->d2vdtdp *= factor;
    target->d2vdp2 *= factor;
}

Constraints *allocConstraintsPointer(void)
{
#ifdef PHMELTS_ADJUSTMENTS
    int i, na = 0;
#endif
    Constraints *p;

    p = (Constraints *) malloc((unsigned) sizeof(Constraints));
#ifdef PHMELTS_ADJUSTMENTS
    p = (Constraints *) malloc((unsigned) sizeof(Constraints));
    p->lambda = (double *) malloc((size_t) (2*nc+4)*sizeof(double));
    p->lambdaO2 = (double *) malloc((size_t)   (nc+1)*sizeof(double));
    p->liquidDelta = (double **) malloc((size_t) nlc*sizeof(double *));
    for (i=0; i<nlc; i++) p->liquidDelta[i] = (double *) malloc((size_t) nlc*sizeof(double));
    p->solidDelta = (double **) malloc((size_t) npc*sizeof(double *));
    for (i=0; i<npc; i++) {
        if (solids[i].type == PHASE) na = MAX(solids[i].na, 1);
        p->solidDelta[i] = (double *) malloc((size_t) na*sizeof(double));
    }
#endif
    return p;
}

void destroyConstraintsStructure(void *pt)
{
    Constraints *p = pt;
#ifdef PHMELTS_ADJUSTMENTS
    int i;

    for (i=0; i<nlc; i++) free((p->liquidDelta)[i]);
    for (i=0; i<npc; i++) free((p->solidDelta)[i]);

    free(p->lambda);
    free(p->lambdaO2);
    free(p->liquidDelta);
    free(p->solidDelta);
#endif
    free(p);
}

SilminState *allocSilminStatePointer(void)
{
    int i;
    SilminState *p;

    /* assume initialization with zero bytes yields default entries */
    p = (SilminState *) calloc((size_t) 1, sizeof(SilminState));
#ifdef PHMELTS_ADJUSTMENTS
    p->nLiquidCoexist = 1;
    p->liquidMass     = 0.0;
#endif

    /* allocate space for mandatory entries with zero bytes         */
    p->bulkComp       = (double *)     calloc((size_t)  nc,     sizeof(double));
    p->dspBulkComp    = (double *)     calloc((size_t)  nc,     sizeof(double));
    p->liquidComp     = (double **)    calloc((size_t)   1,     sizeof(double *));
    p->liquidComp[0]  = (double *)     calloc((size_t) nlc,     sizeof(double));
    p->liquidDelta    = (double **)    calloc((size_t)   1,     sizeof(double *));
    p->liquidDelta[0] = (double *)     calloc((size_t) nlc,     sizeof(double));
    p->solidComp      = (double **)    calloc((size_t) npc,     sizeof(double *));
    p->nSolidCoexist  = (int *)        calloc((size_t) npc,     sizeof(int));
    p->solidDelta     = (double **)    calloc((size_t) npc,     sizeof(double *));
#ifdef PHMELTS_ADJUSTMENTS
    p->incSolids      = (int *)        calloc((size_t) npc,     sizeof(int));
#else
    p->incSolids      = (int *)        calloc((size_t) (npc+1), sizeof(int));
#endif
    p->cylSolids      = (int *)        calloc((size_t) npc,     sizeof(int));

    for (i=0; i<npc; i++) {
        (p->solidComp)[i]  = (double *) calloc((size_t) 1, sizeof(double));
        (p->solidDelta)[i] = (double *) calloc((size_t) 1, sizeof(double));
    }

    return p;
}

void destroySilminStateStructure(void *pt)
{
    SilminState *p = pt;
    int i;

    for (i=0; i<MAX(1,p->nLiquidCoexist); i++) {
        free((p->liquidComp)[i]);
        free((p->liquidDelta)[i]);
    }

    for (i=0; i<npc; i++) {
        free((p->solidComp)[i]);
        free((p->solidDelta)[i]);
    }

    free(p->bulkComp);
    free(p->dspBulkComp);
    free(p->liquidComp);
    free(p->liquidDelta);
    free(p->solidComp);
    free(p->nSolidCoexist);
    free(p->solidDelta);
    free(p->incSolids);
    free(p->cylSolids);
#ifdef PHMELTS_ADJUSTMENTS
    if (p->fracSComp != NULL) {
        for (i=0; i<npc; i++) if ((p->fracSComp)[i] != NULL) free((p->fracSComp)[i]);
        free(p->fracSComp);
    }
    if (p->nFracCoexist != NULL) free(p->nFracCoexist);

    if (p->fracLComp != NULL) free(p->fracLComp);

    if (p->dspAssimComp != NULL) {
        for (i=0; i<(npc+nc); i++)  if ((p->dspAssimComp)[i] != NULL) free((p->dspAssimComp)[i]);
        free(p->dspAssimComp);
    }
    if (p->assimComp    != NULL) {
        for (i=0; i<(npc+nlc); i++) if ((p->assimComp)[i]    != NULL) free((p->assimComp)[i]);
        free(p->assimComp);
    }
    if (p->nDspAssimComp != NULL) free(p->nDspAssimComp);
    if (p->nAssimComp    != NULL) free(p->nAssimComp);

    if ((p->ySol) != NULL) free(p->ySol);
    if ((p->yLiq) != NULL) free(p->yLiq);
#endif
    free(p);
}

SilminState *copySilminStateStructure(SilminState *pOld, SilminState *pNew)
{
    int i, j, k, ns;

    if (pNew == (SilminState *) NULL) {
        /* allocate space for the SilminState structure and fill with defaults   */
        pNew = (SilminState *) calloc((size_t) 1, sizeof(SilminState));

        /* copy across non-pointer specific entries                              */
        memcpy(pNew, pOld, sizeof(SilminState));

        /* annihilate copies pointers of mandatory arrays                        */
        pNew->bulkComp       = (double *)     calloc((size_t)  nc,     sizeof(double));
        pNew->dspBulkComp    = (double *)     calloc((size_t)  nc,     sizeof(double));
        pNew->liquidComp     = (double **)    calloc((size_t)   1,     sizeof(double *));
        pNew->liquidComp[0]  = (double *)     calloc((size_t) nlc,     sizeof(double));
        pNew->liquidDelta    = (double **)    calloc((size_t)   1,     sizeof(double *));
        pNew->liquidDelta[0] = (double *)     calloc((size_t) nlc,     sizeof(double));
        pNew->solidComp      = (double **)    calloc((size_t) npc,     sizeof(double *));
        pNew->nSolidCoexist  = (int *)        calloc((size_t) npc,     sizeof(int));
        pNew->solidDelta     = (double **)    calloc((size_t) npc,     sizeof(double *));
#ifdef PHMELTS_ADJUSTMENTS
        pNew->incSolids      = (int *)        calloc((size_t) npc,     sizeof(int));
#else
        pNew->incSolids      = (int *)        calloc((size_t) (npc+1), sizeof(int));
#endif
        pNew->cylSolids      = (int *)        calloc((size_t) npc,     sizeof(int));

        /* copy across mandatory entry arrays (space has already been allocated) */
        for (i=0; i<nc; i++) {
            (pNew->bulkComp)[i]      = (pOld->bulkComp)[i];
            (pNew->dspBulkComp)[i]   = (pOld->dspBulkComp)[i];
        }
        for (i=0; i<nlc; i++) {
            (pNew->liquidComp)[0][i]    = (pOld->liquidComp)[0][i];
            (pNew->liquidDelta)[0][i]   = (pOld->liquidDelta)[0][i];
        }
        for (i=0; i<npc; i++) {
            (pNew->nSolidCoexist)[i] = (pOld->nSolidCoexist)[i];
            (pNew->incSolids)[i]     = (pOld->incSolids)[i];
            (pNew->cylSolids)[i]     = (pOld->cylSolids)[i];
        }
#ifndef PHMELTS_ADJUSTMENTS
        (pNew->incSolids)[npc]     = (pOld->incSolids)[npc];
#endif

        /* Now copy the problem specific entries                                 */
        if (pNew->nLiquidCoexist > 1) {
            pNew->liquidComp  = (double **) realloc(pNew->liquidComp,  (size_t) (pNew->nLiquidCoexist)*sizeof(double *));
            pNew->liquidDelta = (double **) realloc(pNew->liquidDelta, (size_t) (pNew->nLiquidCoexist)*sizeof(double *));
            for (i=1; i<pNew->nLiquidCoexist; i++) {
                (pNew->liquidComp)[i]  = (double *) calloc((size_t) nlc, sizeof(double));
                (pNew->liquidDelta)[i] = (double *) calloc((size_t) nlc, sizeof(double));
                for (j=0; j<nlc; j++) {
                    (pNew->liquidComp)[i][j]  = (pOld->liquidComp)[i][j];
                    (pNew->liquidDelta)[i][j] = (pOld->liquidDelta)[i][j];
                }
            }
        }

        for (i=0; i<npc; i++) {
            if ((ns = (pOld->nSolidCoexist)[i]) > 0) {
                if (solids[i].type == PHASE) {
                    (pNew->solidComp)[i]  = (double *) malloc((unsigned) ns*sizeof(double));
                    (pNew->solidDelta)[i] = (double *) malloc((unsigned) ns*sizeof(double));
                    for (j=0; j<ns; j++) {
                        (pNew->solidComp)[i][j]  = (pOld->solidComp)[i][j];
                        (pNew->solidDelta)[i][j] = (pOld->solidDelta)[i][j];
                    }
                    if (solids[i].na > 1) {
      	    for (k=0; k<solids[i].na; k++) {
                            (pNew->solidComp)[i+1+k]  = (double *) malloc((unsigned) ns*sizeof(double));
                            (pNew->solidDelta)[i+1+k] = (double *) malloc((unsigned) ns*sizeof(double));
                            for (j=0; j<ns; j++) {
                                (pNew->solidComp)[i+1+k][j]  = (pOld->solidComp)[i+1+k][j];
                                (pNew->solidDelta)[i+1+k][j] = (pOld->solidDelta)[i+1+k][j];
                            }
                        }
                        i += solids[i].na;
	        }
                }
            } else {
                (pNew->solidComp)[i]  = (double *) calloc((unsigned) 1, sizeof(double));
                (pNew->solidDelta)[i] = (double *) calloc((unsigned) 1, sizeof(double));
            }
        }

        if ((pOld->fractionateSol || pOld->fractionateFlu) && (pOld->nFracCoexist != NULL) && (pOld->fracSComp != NULL)) {
            pNew->fracSComp    = (double **) calloc((unsigned) npc, sizeof(double *));
            pNew->nFracCoexist = (int *)     calloc((unsigned) npc, sizeof(int));
            for (i=0; i<npc; i++) {
                (pNew->nFracCoexist)[i] = (pOld->nFracCoexist)[i];
                if ((ns = (pOld->nFracCoexist)[i]) > 0) {
                    (pNew->fracSComp)[i]  = (double *) malloc((unsigned) ns*sizeof(double));
                    for (j=0; j<ns; j++) (pNew->fracSComp)[i][j] = (pOld->fracSComp)[i][j];
	  if (solids[i].na > 1) for (k=0; k<solids[i].na; k++) {
                            (pNew->fracSComp)[i+1+k]  = (double *) malloc((unsigned) ns*sizeof(double));
                            for (j=0; j<ns; j++) (pNew->fracSComp)[i+1+k][j] = (pOld->fracSComp)[i+1+k][j];
	  }
                }
            }
        }

        if (pOld->fractionateLiq && (pOld->fracLComp != NULL)) {
            pNew->fracLComp = (double *) calloc((unsigned) nlc, sizeof(double));
            for (i=0; i<nlc; i++) (pNew->fracLComp)[i] = (pOld->fracLComp)[i];
        }

        if (pOld->assimilate) {
            pNew->dspAssimComp
                = (double **) calloc((unsigned) (npc+nc), sizeof(double *));
            pNew->assimComp
                = (double **) calloc((unsigned) (npc+nlc), sizeof(double *));
            pNew->nDspAssimComp
                = (int *)     calloc((unsigned) (npc+nc), sizeof(int));
            pNew->nAssimComp
                = (int *)     calloc((unsigned) (npc+nlc), sizeof(int));
            for (i=0; i<(npc+nc); i++) {
                (pNew->nDspAssimComp)[i] = (pOld->nDspAssimComp)[i];
                if ((ns = (pOld->nDspAssimComp)[i]) > 0) {
                    (pNew->dspAssimComp)[i]
                        = (double *) malloc((unsigned) ns*sizeof(double));
                    for (j=0; j<ns; j++)
                        (pNew->dspAssimComp)[i][j] = (pOld->dspAssimComp)[i][j];
                }
            }
            for (i=0; i<(npc+nlc); i++) {
                (pNew->nAssimComp)[i] = (pOld->nAssimComp)[i];
                if ((ns = (pOld->nAssimComp)[i]) > 0) {
                    (pNew->assimComp)[i] = (double *) malloc((unsigned) ns*sizeof(double));
                    for (j=0; j<ns; j++) (pNew->assimComp)[i][j] = (pOld->assimComp)[i][j];
                }
            }
        }
    } else {
        /* copy across mandatory entries (space has already been allocated)      */
        for (i=0; i<nc; i++) (pNew->bulkComp)[i]    = (pOld->bulkComp)[i];
        for (i=0; i<nc; i++) (pNew->dspBulkComp)[i] = (pOld->dspBulkComp)[i];
        memcpy(&(pNew->bulkTD), &(pOld->bulkTD), sizeof(ThermoData));
        for (i=0; i<nlc; i++) (pNew->liquidComp)[0][i]  = (pOld->liquidComp)[0][i];
        for (i=0; i<nlc; i++) (pNew->liquidDelta)[0][i] = (pOld->liquidDelta)[0][i];
        pNew->liquidMass = pOld->liquidMass;
        memcpy(&(pNew->liquidTD), &(pOld->liquidTD), sizeof(ThermoData));

        if (pNew->nLiquidCoexist != pOld->nLiquidCoexist) {
            if (pNew->nLiquidCoexist > pOld->nLiquidCoexist) for (i=pOld->nLiquidCoexist; i<pNew->nLiquidCoexist; i++) {
                free((pNew->liquidComp)[i]);
                free((pNew->liquidDelta)[i]);
            }
            pNew->liquidComp  = (double **) realloc(pNew->liquidComp,  (size_t) (pOld->nLiquidCoexist)*sizeof(double *));
            pNew->liquidDelta = (double **) realloc(pNew->liquidDelta, (size_t) (pOld->nLiquidCoexist)*sizeof(double *));
            if (pNew->nLiquidCoexist < pOld->nLiquidCoexist) for (i=pNew->nLiquidCoexist; i<pOld->nLiquidCoexist; i++) {
                pNew->liquidComp[i]  = (double *) calloc((size_t) nlc, sizeof(double));
                pNew->liquidDelta[i] = (double *) calloc((size_t) nlc, sizeof(double));
            }
            for (i=1; i<pOld->nLiquidCoexist; i++) for (j=0; j<nlc; j++) {
                pNew->liquidComp[i][j]  = pOld->liquidComp[i][j];
                pNew->liquidDelta[i][j] = pOld->liquidDelta[i][j];
            }
            pNew->nLiquidCoexist = pOld->nLiquidCoexist;
        }

        for (i=0; i<npc; i++) {
            if ((solids[i].type == PHASE) && ((pNew->nSolidCoexist)[i] != (pOld->nSolidCoexist)[i])) {
                int newSize = MAX((pOld->nSolidCoexist)[i], 1);
                pNew->solidComp[i]  = (double *) realloc(pNew->solidComp[i], (unsigned) newSize*sizeof(double));
                pNew->solidDelta[i] = (double *) realloc(pNew->solidDelta[i],(unsigned) newSize*sizeof(double));
	if (solids[i].na > 1) {
	  for (j=0; j<solids[i].na; j++) {
                        pNew->solidComp[i+1+j]  = (double *) realloc(pNew->solidComp[i+1+j], (unsigned) newSize*sizeof(double));
                        pNew->solidDelta[i+1+j] = (double *) realloc(pNew->solidDelta[i+1+j],(unsigned) newSize*sizeof(double));
	  }
	  i += solids[i].na;
	}
            }
        }

        for (i=0; i<npc; i++) {
            for (j=0; j<(pOld->nSolidCoexist)[i]; j++) {
                pNew->solidComp[i][j]  = pOld->solidComp[i][j];
                pNew->solidDelta[i][j] = pOld->solidDelta[i][j];
                if (solids[i].type == PHASE && solids[i].na > 1)
                    for (k=0; k<solids[i].na; k++) {
                        (pNew->solidComp)[i+1+k][j]  = (pOld->solidComp)[i+1+k][j];
                        (pNew->solidDelta)[i+1+k][j] = (pOld->solidDelta)[i+1+k][j];
                    }
            }
            (pNew->nSolidCoexist)[i] = (pOld->nSolidCoexist)[i];
        }

        pNew->solidMass = pOld->solidMass;

        for (i=0; i<npc; i++) (pNew->incSolids)[i] = (pOld->incSolids)[i];
        for (i=0; i<npc; i++) (pNew->cylSolids)[i] = (pOld->cylSolids)[i];
#ifdef PHMELTS_ADJUSTMENTS
        pNew->incLiquids           = pOld->incLiquids;
        pNew->cylLiquids           = pOld->cylLiquids;
#else
        (pNew->incSolids)[npc]     = (pOld->incSolids)[npc];
        (pNew->cylSolids)[npc]     = (pOld->cylSolids)[npc];
#endif

        memcpy(&(pNew->solidTD), &(pOld->solidTD), sizeof(ThermoData));
        pNew->dspTstart = pOld->dspTstart;
        pNew->dspTstop  = pOld->dspTstop;
        pNew->dspTinc   = pOld->dspTinc;
        pNew->dspHstop  = pOld->dspHstop;
        pNew->dspHinc   = pOld->dspHinc;
        pNew->dspSstop  = pOld->dspSstop;
        pNew->dspSinc   = pOld->dspSinc;
        pNew->T         = pOld->T;
        pNew->dspPstart = pOld->dspPstart;
        pNew->dspPstop  = pOld->dspPstop;
        pNew->dspPinc   = pOld->dspPinc;
        pNew->dspVstop  = pOld->dspVstop;
        pNew->dspVinc   = pOld->dspVinc;
        pNew->dspDPDt   = pOld->dspDPDt;
        pNew->dspDPDH   = pOld->dspDPDH;
        pNew->dspDPDS   = pOld->dspDPDS;
#ifdef PHMELTS_ADJUSTMENTS
        pNew->dspDTDV   = pOld->dspDTDV;
#else
        pNew->dspDVDt   = pOld->dspDVDt;
#endif
        pNew->P         = pOld->P;
        pNew->fo2       = pOld->fo2;
        pNew->fo2Path   = pOld->fo2Path;
        pNew->fo2Delta  = pOld->fo2Delta;
        pNew->oxygen    = pOld->oxygen;

#ifdef PHMELTS_ADJUSTMENTS
        pNew->fo2Alt    = pOld->fo2Alt;
        pNew->fo2Liq    = pOld->fo2Liq;
        pNew->fo2Sol    = pOld->fo2Sol;
        pNew->fo2Iter   = pOld->fo2Iter;
        pNew->H2Obuffer = pOld->H2Obuffer;
        pNew->aH2O      = pOld->aH2O;

#endif

        if ((pOld->fractionateSol || pOld->fractionateFlu) && (pOld->nFracCoexist != NULL) && (pOld->fracSComp != NULL)) {
            if (pNew->nFracCoexist == NULL) {
                pNew->fracSComp    = (double **) calloc((unsigned) npc, sizeof(double *));
                pNew->nFracCoexist = (int *)     calloc((unsigned) npc, sizeof(int));
            }
            for (i=0; i<npc; i++) {
                if ((pNew->nFracCoexist)[i] != (pOld->nFracCoexist)[i]) {
                    int newSize = MAX((pOld->nFracCoexist)[i], 1);
                    pNew->fracSComp[i]  = (double *) REALLOC(pNew->fracSComp[i],(unsigned) newSize*sizeof(double));
	        if (solids[i].na > 1) for (k=0; k<solids[i].na; k++) {
	          pNew->fracSComp[i+1+k]  = (double *) REALLOC(pNew->fracSComp[i+1+k],(unsigned) newSize*sizeof(double));
      	  }
                }
                if ((pOld->nFracCoexist)[i] > 0) for (j=0; j<(pOld->nFracCoexist)[i]; j++) {
	        pNew->fracSComp[i][j] = pOld->fracSComp[i][j];
                    if (solids[i].na > 1) for (k=0; k<solids[i].na; k++) pNew->fracSComp[i+1+k][j] = pOld->fracSComp[i+1+k][j];
                }
                (pNew->nFracCoexist)[i] = (pOld->nFracCoexist)[i];
            }
        } else if ((pNew->fractionateSol || pNew->fractionateFlu) && (pNew->nFracCoexist != NULL) && (pNew->fracSComp != NULL)) {
            for (i=0; i<npc; i++) if (pNew->fracSComp[i] != NULL) free (pNew->fracSComp[i]);
            free(pNew->nFracCoexist); pNew->nFracCoexist = NULL;
            free(pNew->fracSComp);    pNew->fracSComp    = NULL;
        }
        pNew->fractionateSol = pOld->fractionateSol;
        pNew->fractionateFlu = pOld->fractionateFlu;

        if (pOld->fractionateLiq && (pOld->fracLComp != NULL)) {
            if (pNew->fracLComp == NULL) pNew->fracLComp = (double *) calloc((unsigned) nlc, sizeof(double));
            for (i=0; i<nlc; i++) pNew->fracLComp[i] = pOld->fracLComp[i];
        } else if (pNew->fractionateLiq && (pNew->fracLComp != NULL)) {
            free(pNew->fracLComp);
            pNew->fracLComp = NULL;
        }
        pNew->fractionateLiq = pOld->fractionateLiq;

        pNew->fracMass       = pOld->fracMass;
#ifdef PHMELTS_ADJUSTMENTS
        pNew->incLiquids     = pOld->incLiquids;
#else
        pNew->multipleLiqs   = pOld->multipleLiqs;
#endif
        pNew->isenthalpic    = pOld->isenthalpic;
        pNew->refEnthalpy    = pOld->refEnthalpy;
        pNew->isentropic     = pOld->isentropic;
        pNew->refEntropy     = pOld->refEntropy;
        pNew->tDelta         = pOld->tDelta;

        pNew->isochoric      = pOld->isochoric;
        pNew->refVolume      = pOld->refVolume;
        pNew->pDelta         = pOld->pDelta;

        if (pOld->assimilate) {
            pNew->dspAssimComp
                = (double **) calloc((unsigned) (npc+nc), sizeof(double *));
            pNew->assimComp
                = (double **) calloc((unsigned) (npc+nlc), sizeof(double *));
            pNew->nDspAssimComp
                = (int *)     calloc((unsigned) (npc+nc), sizeof(int));
            pNew->nAssimComp
                = (int *)     calloc((unsigned) (npc+nlc), sizeof(int));
            for (i=0; i<(nlc+npc); i++) {
                if ((pNew->nAssimComp)[i] != (pOld->nAssimComp)[i]) {
                    int newSize = MAX((pOld->nAssimComp)[i], 1);
                    (pNew->assimComp)[i]   = (double *) REALLOC((pNew->assimComp)[i],
                        (unsigned) newSize*sizeof(double));
                }
            }
            for (i=0; i<(nc+npc); i++) {
                if ((pNew->nDspAssimComp)[i] != (pOld->nDspAssimComp)[i]) {
                    int newSize = MAX((pOld->nDspAssimComp)[i], 1);
                    (pNew->dspAssimComp)[i]   = (double *) REALLOC((pNew->dspAssimComp)[i],
                        (unsigned) newSize*sizeof(double));
                }
            }
            for (i=0; i<(nc+npc); i++) {
                for (j=0; j<(pOld->nAssimComp)[i]; j++)
                    (pNew->assimComp)[i][j]    = (pOld->assimComp)[i][j];
                for (j=0; j<(pOld->nDspAssimComp)[i]; j++)
                    (pNew->dspAssimComp)[i][j] = (pOld->dspAssimComp)[i][j];
                (pNew->nAssimComp)[i]        = (pOld->nAssimComp)[i];
                (pNew->nDspAssimComp)[i]     = (pOld->nDspAssimComp)[i];
            }
        } else if (pNew->assimilate) {
            for (i=0; i<(nlc+npc); i++)
                if (pNew->assimComp[i] != NULL)    free (pNew->assimComp[i]);
            for (i=0; i<(nc+npc); i++)
                if (pNew->dspAssimComp[i] != NULL) free (pNew->dspAssimComp[i]);
            free(pNew->nAssimComp);    pNew->nAssimComp    = NULL;
            free(pNew->assimComp);     pNew->assimComp     = NULL;
            free(pNew->nDspAssimComp); pNew->nDspAssimComp = NULL;
            free(pNew->dspAssimComp);  pNew->dspAssimComp  = NULL;
        }
        pNew->assimilate    = pOld->assimilate;
        pNew->dspAssimUnits = pOld->dspAssimUnits;
        pNew->dspAssimMass  = pOld->dspAssimMass;
        pNew->dspAssimT     = pOld->dspAssimT;
        pNew->dspAssimInc   = pOld->dspAssimInc;
        pNew->dspAssimLiqM  = pOld->dspAssimLiqM;
        pNew->assimMass     = pOld->assimMass;
        pNew->assimT        = pOld->assimT;
        pNew->assimInc      = pOld->assimInc;
        memcpy(&(pNew->assimTD), &(pOld->assimTD), sizeof(ThermoData));

#ifdef PHMELTS_ADJUSTMENTS
        pNew->txtOutput     = pOld->txtOutput;
#else
        pNew->plotState     = pOld->plotState;
#endif

    }

    return pNew;
}

#ifndef BATCH_VERSION

/*
#define ERROR(string) \
{ \
    XmString csString = XmStringCreateLtoR(string, "ISO8859-1"); \
    XtVaSetValues(message, XmNmessageString, csString, NULL); \
    XtManageChild(message); \
    XmStringFree(csString); \
    free(diff); \
    return mask | SILMIN_STATE_CHANGE_FATAL_ERROR; \
}

int checkStateAgainstInterface(void) {}

#undef ERROR
*/

#endif /* BATCH_VERSION */

/******************************************************************************
 * The screen update functions:
 ******************************************************************************/

void updateBulkADB(void)
{
#ifndef BATCH_VERSION
    static char compositionEntry[9];
#endif
    double sum, *temporary;
    int i, j, nl;
    int hasLiquid = (silminState->liquidMass != 0.0);

    temporary = (double *) calloc((size_t) nc, sizeof(double));

    /* convert moles of liquid components into grams of oxides */
    if (hasLiquid) for (i=0; i<nc; i++) for (nl=0, temporary[i]=0.0; nl<silminState->nLiquidCoexist; nl++) for (j=0; j<nlc; j++)
        temporary[i] += (silminState->liquidComp)[nl][j]*(liquid[j].liqToOx)[i];

    if (silminState->fo2Path == FO2_NONE && hasLiquid) { /* This returns the AVERAGE log fO2 */
        double *oxides = (double *) malloc((size_t) nc*sizeof(double));
        double logfO2 = 0.0;
        for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
            for (i=0; i<nc; i++) for (j=0, oxides[i]=0.0; j<nlc; j++) oxides[i] += (silminState->liquidComp)[nl][j]*(liquid[j].liqToOx)[i];
            conLiq(FIRST, SEVENTH, silminState->T, silminState->P, oxides, NULL, NULL, NULL, NULL, NULL, &(silminState->fo2));
            logfO2 += silminState->fo2;
        }
        silminState->fo2 = logfO2/((double) nl);
        free(oxides);
    } else if (silminState->fo2Path == FO2_NONE && !hasLiquid) {
        double muO2;
        subsolidusmuO2(FIRST, &muO2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
        silminState->fo2 = muO2/(R*silminState->T*log(10.0));
    }
    for (i=0, sum=0.0; i<nc; i++) {
        temporary[i] *= bulkSystem[i].mw;     /* grams of oxide components */
        sum          += temporary[i];         /* total grams of oxides     */
    }
    silminState->liquidMass = sum;

    /* display grams oxides and preserve an exact reference copy */
    for (i=0; i<nc; i++) {
        (silminState->dspBulkComp)[i] = (silminState->bulkComp)[i]*bulkSystem[i].mw;
    }
    free(temporary);
}

#ifndef BATCH_VERSION

/*
typedef struct _clsStable {
    Widget   name;
    XmString label;
} ClsStable;

static ClsStable *updateCompADBaddButton(char *label)
static void updateCompADBoxides(double *oxVal)
void updateCompADB(void)
void updateSolidADB(double *ySol, double *yLiq)
void updateTpPADB(int member)
*/

#endif /* BATCH_VERSION */

void correctTforChangeInEnthalpy(void)
{
    static double *mSol, *rLiq, *rSol, *mOx;
    double pTemp, mTotal, hTotal, cpTotal, residual = DBL_MAX;
    int i, j, k, nl, ns, iter = 0;
    int hasLiquid = (silminState->liquidMass != 0.0);

    if (mSol == NULL) {
        for (i=0, j=1, k=1; i<npc; i++) if (solids[i].type == PHASE)
            { j = MAX(j, solids[i].nr); k = MAX(k, solids[i].na); }

        mSol = (double *) malloc((size_t)       k*sizeof(double));
        rLiq = (double *) malloc((size_t) (nlc-1)*sizeof(double));
        rSol = (double *) malloc((size_t)       j*sizeof(double));
        mOx  = (double *) malloc((size_t)      nc*sizeof(double));
    }

    while ((fabs(residual) > 10.0*DBL_EPSILON*fabs(silminState->refEnthalpy)) && (iter < 50)) {

        if (silminState->fo2Path != FO2_NONE && hasLiquid) {
            silminState->fo2 = getlog10fo2(silminState->T, silminState->P, silminState->fo2Path);
            for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                for (i=0; i<nc; i++) for (j=0, mOx[i]=0.0; j<nlc; j++) mOx[i] += (liquid[j].liqToOx)[i]*(silminState->liquidComp)[nl][j];
                conLiq(FIRST | SEVENTH, FIRST, silminState->T, silminState->P, mOx, NULL, NULL, NULL, NULL, NULL, &(silminState->fo2));
                for (i=0; i<nlc; i++) for (j=0, (silminState->liquidComp)[nl][i]=0.0; j<nc; j++) (silminState->liquidComp)[nl][i] += (bulkSystem[j].oxToLiq)[i]*mOx[j];
            }
        } else if (silminState->fo2Path != FO2_NONE && !hasLiquid) {
            double muO2;
            silminState->fo2 = getlog10fo2(silminState->T, silminState->P, silminState->fo2Path);
            muO2 = silminState->fo2*(R*silminState->T*log(10.0));
            if (!subsolidusmuO2(0,  &muO2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)) {
                printf("Failure to impose fo2 buffer in subsolidus.  Releasing buffer constraint from system.\n");
#ifdef PHMELTS_ADJUSTMENTS
	silminState->fo2Liq = silminState->fo2Path;
#endif
	silminState->fo2Path = FO2_NONE;
#ifndef BATCH_VERSION
                XmToggleButtonGadgetSetState(tg_path_none, True, True);
#endif
            }
        }

        hTotal  = 0.0;
        cpTotal = 0.0;
        if (hasLiquid) {
            for (i=0; i<nlc; i++) if ((silminState->liquidComp)[0][i] != 0.0)
                gibbs(silminState->T, silminState->P, (char *) liquid[i].label, &(liquid[i].ref), &(liquid[i].liq), &(liquid[i].fus), &(liquid[i].cur));
            for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                for (i=0, mTotal= 0.0; i<nlc; i++) {
                    mTotal  += (silminState->liquidComp)[nl][i];
                    hTotal  += (silminState->liquidComp)[nl][i]*(liquid[i].cur).h;
                    cpTotal += (silminState->liquidComp)[nl][i]*(liquid[i].cur).cp;
                }
                conLiq(SECOND, THIRD, silminState->T, silminState->P, NULL, (silminState->liquidComp)[nl], rLiq, NULL, NULL, NULL, NULL);
                hmixLiq(FIRST,  silminState->T, silminState->P, rLiq, &pTemp, NULL);
                hTotal += mTotal*pTemp;
                cpmixLiq(FIRST, silminState->T, silminState->P, rLiq, &pTemp, NULL, NULL);
                cpTotal += mTotal*pTemp;
            }
        }

        for (i=0; i<npc; i++)
            for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
                mTotal = (silminState->solidComp)[i][ns];
                if (solids[i].na == 1) {
                    gibbs(silminState->T, silminState->P, (char *) solids[i].label, &(solids[i].ref), NULL, NULL, &(solids[i].cur));
                    hTotal  += mTotal*(solids[i].cur).h;
                    cpTotal += mTotal*(solids[i].cur).cp;
                } else {
                    for (j=0; j<solids[i].na; j++) {
                        mSol[j] = (silminState->solidComp)[i+1+j][ns];
                        gibbs(silminState->T, silminState->P, (char *) solids[i+1+j].label, &(solids[i+1+j].ref), NULL, NULL, &(solids[i+1+j].cur));
                        hTotal  += mSol[j]*(solids[i+1+j].cur).h;
                        cpTotal += mSol[j]*(solids[i+1+j].cur).cp;
                    }
                    (*solids[i].convert)(SECOND, THIRD, silminState->T, silminState->P, NULL, mSol, rSol, NULL, NULL, NULL, NULL, NULL);
                    (*solids[i].hmix)(FIRST, silminState->T, silminState->P, rSol, &pTemp);
                    hTotal += mTotal*pTemp;
                    (*solids[i].cpmix)(FIRST, silminState->T, silminState->P, rSol, &pTemp, NULL, NULL);
                    cpTotal += mTotal*pTemp;
                }
            }

        residual = hTotal - silminState->refEnthalpy;
        silminState->T -= residual/cpTotal;
        iter++;

#ifdef DEBUG
    printf("Enthalpy corr: Iter = %d, T(C) = %.2f, H-Href = %g, Cp = %g\n", iter, silminState->T-273.15, residual, cpTotal);
#endif
    }
#ifdef DEBUG
    printf("Enthalpy corr: Iter = %d, T(C) = %.2f, H-Href = %g, Cp = %g\n", iter, silminState->T-273.15, residual, cpTotal);
#endif
}

void correctTforChangeInEntropy(void)
{
    static double *mSol, *rLiq, *rSol, *mOx;
    double pTemp, mTotal, sTotal, cpTotal, residual = DBL_MAX;
    int i, j, k, nl, ns, iter = 0;
    int hasLiquid = (silminState->liquidMass != 0.0);

    if (mSol == NULL) {
        for (i=0, j=1, k=1; i<npc; i++) if (solids[i].type == PHASE)
            { j = MAX(j, solids[i].nr); k = MAX(k, solids[i].na); }

        mSol = (double *) malloc((size_t)       k*sizeof(double));
        rLiq = (double *) malloc((size_t) (nlc-1)*sizeof(double));
        rSol = (double *) malloc((size_t)       j*sizeof(double));
        mOx  = (double *) malloc((size_t)      nc*sizeof(double));
    }

    while ((fabs(residual) > 10.0*DBL_EPSILON*fabs(silminState->refEntropy)) && (iter < 50)) {

        if (silminState->fo2Path != FO2_NONE && hasLiquid) {
            silminState->fo2 = getlog10fo2(silminState->T, silminState->P, silminState->fo2Path);
            for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                for (i=0; i<nc; i++) for (j=0, mOx[i]=0.0; j<nlc; j++) mOx[i] += (liquid[j].liqToOx)[i]*(silminState->liquidComp)[nl][j];
                conLiq(FIRST | SEVENTH, FIRST, silminState->T, silminState->P, mOx, NULL, NULL, NULL, NULL, NULL, &(silminState->fo2));
                for (i=0; i<nlc; i++) for (j=0, (silminState->liquidComp)[nl][i]=0.0; j<nc; j++) (silminState->liquidComp)[nl][i] += (bulkSystem[j].oxToLiq)[i]*mOx[j];
            }
        } else if (silminState->fo2Path != FO2_NONE && !hasLiquid) {
            double muO2;
            silminState->fo2 = getlog10fo2(silminState->T, silminState->P, silminState->fo2Path);
            muO2 = silminState->fo2*(R*silminState->T*log(10.0));
            if (!subsolidusmuO2(0, &muO2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)) {
                printf("Failure to impose fo2 buffer in subsolidus.  Releasing buffer constraint on the system.\n");
#ifdef PHMELTS_ADJUSTMENTS
	silminState->fo2Liq = silminState->fo2Path;
#endif
	silminState->fo2Path = FO2_NONE;
#ifndef BATCH_VERSION
                XmToggleButtonGadgetSetState(tg_path_none, True, True);
#endif
            }
        }

        sTotal  = 0.0;
        cpTotal = 0.0;
        if (hasLiquid) {
            for (i=0; i<nlc; i++) if ((silminState->liquidComp)[0][i] != 0.0)
                gibbs(silminState->T, silminState->P, (char *) liquid[i].label, &(liquid[i].ref), &(liquid[i].liq), &(liquid[i].fus), &(liquid[i].cur));
            for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                for (i=0, mTotal= 0.0; i<nlc; i++) {
                    mTotal  += (silminState->liquidComp)[nl][i];
                    sTotal  += (silminState->liquidComp)[nl][i]*(liquid[i].cur).s;
                    cpTotal += (silminState->liquidComp)[nl][i]*(liquid[i].cur).cp;
                }
                conLiq(SECOND, THIRD, silminState->T, silminState->P, NULL, (silminState->liquidComp)[nl], rLiq, NULL, NULL, NULL, NULL);
                smixLiq(FIRST,  silminState->T, silminState->P, rLiq, &pTemp, NULL, NULL, NULL);
                sTotal += mTotal*pTemp;
                cpmixLiq(FIRST, silminState->T, silminState->P, rLiq, &pTemp, NULL, NULL);
                cpTotal += mTotal*pTemp;
            }
        }

        for (i=0; i<npc; i++)
            for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
                mTotal = (silminState->solidComp)[i][ns];
                if (solids[i].na == 1) {
                    gibbs(silminState->T, silminState->P, (char *) solids[i].label, &(solids[i].ref), NULL, NULL, &(solids[i].cur));
                    sTotal  += mTotal*(solids[i].cur).s;
                    cpTotal += mTotal*(solids[i].cur).cp;
                } else {
                    for (j=0; j<solids[i].na; j++) {
                        mSol[j] = (silminState->solidComp)[i+1+j][ns];
                        gibbs(silminState->T, silminState->P, (char *) solids[i+1+j].label, &(solids[i+1+j].ref), NULL, NULL, &(solids[i+1+j].cur));
                        sTotal  += mSol[j]*(solids[i+1+j].cur).s;
                        cpTotal += mSol[j]*(solids[i+1+j].cur).cp;
                    }
                    (*solids[i].convert)(SECOND, THIRD, silminState->T, silminState->P, NULL, mSol, rSol, NULL, NULL, NULL, NULL, NULL);
                    (*solids[i].smix)(FIRST, silminState->T, silminState->P, rSol, &pTemp, NULL, NULL);
                    sTotal += mTotal*pTemp;
                    (*solids[i].cpmix)(FIRST, silminState->T, silminState->P, rSol, &pTemp, NULL, NULL);
                    cpTotal += mTotal*pTemp;
                }
            }

        residual = sTotal - silminState->refEntropy;
        silminState->T -= residual*silminState->T/cpTotal;
        iter++;

#ifdef DEBUG
    printf("Entropy corr: Iter = %d, T(C) = %.2f, S-Sref = %g, Cp = %g\n", iter, silminState->T-273.15, residual, cpTotal);
#endif

    }
#ifdef DEBUG
    printf("Entropy corr: Iter = %d, T(C) = %.2f, S-Sref = %g, Cp = %g\n", iter, silminState->T-273.15, residual, cpTotal);
#endif
}

void correctPforChangeInVolume(void)
{
    static double *mSol, *rLiq, *rSol, *mOx;
    double pTemp, dpTemp, mTotal, vTotal, dvdpTotal, residual = DBL_MAX;
    int i, j, k, nl, ns, iter = 0;
    int hasLiquid = (silminState->liquidMass != 0.0);

    if (mSol == NULL) {
        for (i=0, j=1, k=1; i<npc; i++) if (solids[i].type == PHASE) { j = MAX(j, solids[i].nr); k = MAX(k, solids[i].na); }

        mSol = (double *) malloc((size_t)       k*sizeof(double));
        rLiq = (double *) malloc((size_t) (nlc-1)*sizeof(double));
        rSol = (double *) malloc((size_t)       j*sizeof(double));
        mOx  = (double *) malloc((size_t)      nc*sizeof(double));
    }

    while ((fabs(residual) > 10.0*DBL_EPSILON*fabs(silminState->refVolume)) && (iter < 50)) {

        if (silminState->fo2Path != FO2_NONE && hasLiquid) {
            silminState->fo2 = getlog10fo2(silminState->T, silminState->P, silminState->fo2Path);
            for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                for (i=0; i<nc; i++) for (j=0, mOx[i]=0.0; j<nlc; j++) mOx[i] += (liquid[j].liqToOx)[i]*(silminState->liquidComp)[nl][j];
                conLiq(FIRST | SEVENTH, FIRST, silminState->T, silminState->P, mOx, NULL, NULL, NULL, NULL, NULL, &(silminState->fo2));
                for (i=0; i<nlc; i++) for (j=0, (silminState->liquidComp)[nl][i]=0.0; j<nc; j++) (silminState->liquidComp)[nl][i] += (bulkSystem[j].oxToLiq)[i]*mOx[j];
            }
        } else if (silminState->fo2Path != FO2_NONE && !hasLiquid) {
            double muO2;
            silminState->fo2 = getlog10fo2(silminState->T, silminState->P, silminState->fo2Path);
            muO2 = silminState->fo2*(R*silminState->T*log(10.0));
            if (!subsolidusmuO2(0, &muO2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)) {
                printf("Failure to impose fo2 buffer in subsolidus.  Releasing buffer constraint from the system.\n");
#ifdef PHMELTS_ADJUSTMENTS
                silminState->fo2Liq = silminState->fo2Path;
#endif
	silminState->fo2Path = FO2_NONE;
#ifndef BATCH_VERSION
                XmToggleButtonGadgetSetState(tg_path_none, True, True);
#endif
            }
        }

        vTotal  = 0.0;
        dvdpTotal = 0.0;
        if (hasLiquid) {
            for (i=0, mTotal= 0.0; i<nlc; i++) if ((silminState->liquidComp)[0][i] != 0.0)
                gibbs(silminState->T, silminState->P, (char *) liquid[i].label, &(liquid[i].ref), &(liquid[i].liq), &(liquid[i].fus), &(liquid[i].cur));
            for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                for (i=0, mTotal= 0.0; i<nlc; i++) {
                    mTotal    += (silminState->liquidComp)[nl][i];
                    vTotal    += (silminState->liquidComp)[nl][i]*(liquid[i].cur).v;
                    dvdpTotal += (silminState->liquidComp)[nl][i]*(liquid[i].cur).dvdp;
                }
                conLiq(SECOND, THIRD, silminState->T, silminState->P, NULL, (silminState->liquidComp)[nl], rLiq, NULL, NULL, NULL, NULL);
                vmixLiq(FIRST | FIFTH,  silminState->T, silminState->P, rLiq, &pTemp, NULL, NULL, NULL, &dpTemp, NULL, NULL, NULL, NULL, NULL, NULL);
                vTotal    += mTotal*pTemp;
                dvdpTotal += mTotal*dpTemp;
            }
        }

        for (i=0; i<npc; i++)
            for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
                mTotal = (silminState->solidComp)[i][ns];
                if (solids[i].na == 1) {
                    gibbs(silminState->T, silminState->P, (char *) solids[i].label, &(solids[i].ref), NULL, NULL, &(solids[i].cur));
                    vTotal    += mTotal*(solids[i].cur).v;
                    dvdpTotal += mTotal*(solids[i].cur).dvdp;
                } else {
                    for (j=0; j<solids[i].na; j++) {
                        mSol[j] = (silminState->solidComp)[i+1+j][ns];
                        gibbs(silminState->T, silminState->P, (char *) solids[i+1+j].label, &(solids[i+1+j].ref), NULL, NULL, &(solids[i+1+j].cur));
                        vTotal    += mSol[j]*(solids[i+1+j].cur).v;
                        dvdpTotal += mSol[j]*(solids[i+1+j].cur).dvdp;
                    }
                    (*solids[i].convert)(SECOND, THIRD, silminState->T, silminState->P, NULL, mSol, rSol, NULL, NULL, NULL, NULL, NULL);
                    (*solids[i].vmix)(FIRST | FIFTH, silminState->T, silminState->P, rSol, &pTemp, NULL, NULL, NULL, &dpTemp, NULL, NULL, NULL, NULL, NULL);
                    vTotal    += mTotal*pTemp;
                    dvdpTotal += mTotal*dpTemp;
                }
            }

        residual = vTotal - silminState->refVolume;
        silminState->P -= residual/dvdpTotal;
        iter++;

#ifdef DEBUG
    printf("Volume corr: Iter = %d, P = %.2f, V-Vref = %g, dVdP = %g\n", iter, silminState->P, residual, dvdpTotal);
#endif

    }
#ifdef DEBUG
    printf("Volume corr: Iter = %d, P = %.2f, V-Vref = %g, dVdP = %g\n", iter, silminState->P, residual, dvdpTotal);
#endif
}

#ifdef PHMELTS_ADJUSTMENTS
void correctTforChangeInVolume(void)
{
    static double *mSol, *rLiq, *rSol, *mOx;
    double pTemp, dpTemp, mTotal, vTotal, dvdtTotal, residual = DBL_MAX;
    int i, j, k, nl, ns, iter = 0;
    int hasLiquid = (silminState->liquidMass != 0.0);

    if (mSol == NULL) {
        for (i=0, j=1, k=1; i<npc; i++) if (solids[i].type == PHASE) { j = MAX(j, solids[i].nr); k = MAX(k, solids[i].na); }

        mSol = (double *) malloc((size_t)       k*sizeof(double));
        rLiq = (double *) malloc((size_t) (nlc-1)*sizeof(double));
        rSol = (double *) malloc((size_t)       j*sizeof(double));
        mOx  = (double *) malloc((size_t)      nc*sizeof(double));
    }

    while ((fabs(residual) > 10.0*DBL_EPSILON*fabs(silminState->refVolume)) && (iter < 50)) {

        if (silminState->fo2Path != FO2_NONE && hasLiquid) {
            silminState->fo2 = getlog10fo2(silminState->T, silminState->P, silminState->fo2Path);
            for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                for (i=0; i<nc; i++) for (j=0, mOx[i]=0.0; j<nlc; j++) mOx[i] += (liquid[j].liqToOx)[i]*(silminState->liquidComp)[nl][j];
                conLiq(FIRST | SEVENTH, FIRST, silminState->T, silminState->P, mOx, NULL, NULL, NULL, NULL, NULL, &(silminState->fo2));
                for (i=0; i<nlc; i++) for (j=0, (silminState->liquidComp)[nl][i]=0.0; j<nc; j++) (silminState->liquidComp)[nl][i] += (bulkSystem[j].oxToLiq)[i]*mOx[j];
            }
        } else if (silminState->fo2Path != FO2_NONE && !hasLiquid) {
            double muO2;
            silminState->fo2 = getlog10fo2(silminState->T, silminState->P, silminState->fo2Path);
            muO2 = silminState->fo2*(R*silminState->T*log(10.0));
            if (!subsolidusmuO2(0, &muO2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)) {
                printf("Failure to impose fo2 buffer in subsolidus.  Releasing buffer constraint from the system.\n");
#ifdef PHMELTS_ADJUSTMENTS
                silminState->fo2Liq = silminState->fo2Path;
#endif
	silminState->fo2Path = FO2_NONE;
#ifndef BATCH_VERSION
                XmToggleButtonGadgetSetState(tg_path_none, True, True);
#endif
            }
        }

        vTotal  = 0.0;
        dvdtTotal = 0.0;
        if (hasLiquid) {
            for (i=0, mTotal= 0.0; i<nlc; i++) if ((silminState->liquidComp)[0][i] != 0.0)
                gibbs(silminState->T, silminState->P, (char *) liquid[i].label, &(liquid[i].ref), &(liquid[i].liq), &(liquid[i].fus), &(liquid[i].cur));
            for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                for (i=0, mTotal= 0.0; i<nlc; i++) {
                    mTotal    += (silminState->liquidComp)[nl][i];
                    vTotal    += (silminState->liquidComp)[nl][i]*(liquid[i].cur).v;
                    dvdtTotal += (silminState->liquidComp)[nl][i]*(liquid[i].cur).dvdt;
                }
                conLiq(SECOND, THIRD, silminState->T, silminState->P, NULL, (silminState->liquidComp)[nl], rLiq, NULL, NULL, NULL, NULL);
                vmixLiq(FIRST | FOURTH,  silminState->T, silminState->P, rLiq, &pTemp, NULL, NULL, &dpTemp, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
                vTotal    += mTotal*pTemp;
                dvdtTotal += mTotal*dpTemp;
            }
        }

        for (i=0; i<npc; i++)
            for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
                mTotal = (silminState->solidComp)[i][ns];
                if (solids[i].na == 1) {
                    gibbs(silminState->T, silminState->P, (char *) solids[i].label, &(solids[i].ref), NULL, NULL, &(solids[i].cur));
                    vTotal    += mTotal*(solids[i].cur).v;
                    dvdtTotal += mTotal*(solids[i].cur).dvdt;
                } else {
                    for (j=0; j<solids[i].na; j++) {
                        mSol[j] = (silminState->solidComp)[i+1+j][ns];
                        gibbs(silminState->T, silminState->P, (char *) solids[i+1+j].label, &(solids[i+1+j].ref), NULL, NULL, &(solids[i+1+j].cur));
                        vTotal    += mSol[j]*(solids[i+1+j].cur).v;
                        dvdtTotal += mSol[j]*(solids[i+1+j].cur).dvdt;
                    }
                    (*solids[i].convert)(SECOND, THIRD, silminState->T, silminState->P, NULL, mSol, rSol, NULL, NULL, NULL, NULL, NULL);
                    (*solids[i].vmix)(FIRST | FOURTH, silminState->T, silminState->P, rSol, &pTemp, NULL, NULL, &dpTemp, NULL, NULL, NULL, NULL, NULL, NULL);
                    vTotal    += mTotal*pTemp;
                    dvdtTotal += mTotal*dpTemp;
                }
            }

        residual = vTotal - silminState->refVolume;
        silminState->T -= residual/dvdtTotal;
        iter++;

#ifdef DEBUG
    printf("Volume corr: Iter = %d, T = %.2f, V-Vref = %g, dVdT = %g\n", iter, silminState->T, residual, dvdtTotal);
#endif

    }
#ifdef DEBUG
    printf("Volume corr: Iter = %d, T = %.2f, V-Vref = %g, dVdT = %g\n", iter, silminState->T, residual, dvdtTotal);
#endif
}
#endif

/*************************************************************************
   This (misnamed) function finds a way to adjust the solid assemblage
   so as to change the bulk solid composition by *deltaBulkComp (a vector
   of moles of oxides).  This is used to add or drop solids when liquid is
   absent or to add or drop liquid, without kludging the bulk composition.

   Additional twist, added 7/24/95 -- each column of SolToOx is weighted
   by the silminState->solidComp[], and the solution vector is then un-
   weighted after the SVD.  This forces a solution which changes most those
   phase-components which are most abundant, hence minimizing the reulting
   departure from equilibrium.

   8/10/95 -- add test of each solid, return failure if infeasible
*************************************************************************/

int addOrDropLiquid(double *deltaBulkComp) {
    int i, j, k, l, n, result;
    gsl_matrix *reducedSolToOx, *V;
    gsl_vector *deltaSol, *deltaLiq, *S;
    double *weights;

    result = TRUE;
    /* count solid phase-components in assemblage */
    for (i=0, n=0; i<npc; i++) if (silminState->nSolidCoexist[i]) n+= solids[i].na;
    if (n == 0) return FALSE;

    /* allocate memory */
    if (nc >= n) {
        reducedSolToOx = gsl_matrix_alloc((size_t) nc, (size_t) n);
        deltaLiq = gsl_vector_alloc((size_t) nc);
    }
    else {
        reducedSolToOx = gsl_matrix_alloc((size_t) n, (size_t) n);
        deltaLiq = gsl_vector_alloc((size_t) n);
        gsl_matrix_set_zero(reducedSolToOx); // pad with zeros
        gsl_vector_set_zero(deltaLiq);
    }
    deltaSol = gsl_vector_alloc((size_t) n);
    weights =  (double *) malloc((size_t) n*sizeof(double));

    /* construct reduced SolToLiq matrix */
    for (i=0, k=0; i<npc; i++) {
        if (silminState->nSolidCoexist[i]) {
            if (solids[i].na == 1) {
                weights[k] = silminState->solidComp[i][0];
                if (nc >= n) for (j=0; j<nc; j++) gsl_matrix_set(reducedSolToOx, j, k, (solids[i].solToOx)[j]*weights[k]);
                else for (j=0; j<nc; j++) gsl_matrix_set(reducedSolToOx, k, j, (solids[i].solToOx)[j]*weights[k]);
                k++;
            } else for (l=0; l<solids[i].na; l++) {
                weights[k] = silminState->solidComp[i+1+l][0];
                if (nc >= n) for (j=0; j<nc; j++) gsl_matrix_set(reducedSolToOx, j, k, (solids[i+1+l].solToOx)[j]*weights[k]);
                else for (j=0; j<nc; j++) gsl_matrix_set(reducedSolToOx, k, j, (solids[i+1+l].solToOx)[j]*weights[k]);
                k++;
            }
        }
    }

    V = gsl_matrix_alloc((size_t) n, (size_t) n);
    S = gsl_vector_alloc((size_t) n);

    /* SVD it and backsub deltaBulkComp to get correction vector */
    gsl_linalg_SV_decomp(reducedSolToOx, V, S, deltaSol); // deltaSol used as work space
    for (i=0;i<n;i++) if (gsl_vector_get(S, i) < 1.0e-08) gsl_vector_set(S, i, 0.0);

    for (i=0; i<nc; i++) gsl_vector_set(deltaLiq, i, deltaBulkComp[i]);

    // SVD not implemented for M < N in GSL
    if (nc >= n) gsl_linalg_SV_solve(reducedSolToOx, V, S, deltaLiq, deltaSol);
    else gsl_linalg_SV_solve(V, reducedSolToOx, S, deltaLiq, deltaSol);

    /* apply the correction */
    for (i=0, k=0; i<npc; i++) {
        if (silminState->nSolidCoexist[i]) {
            if (solids[i].na == 1) {
                silminState->solidComp[i][0] -= weights[k]*gsl_vector_get(deltaSol, k);
                k++;
                if (silminState->solidComp[i][0] <= 0.0) result = FALSE;
            } else {
                double *mSol = (double *) malloc((size_t) solids[i].na*sizeof(double));
                for (j=0;j<solids[i].na;j++) {
                    mSol[j] = silminState->solidComp[i+1+j][0];
                    silminState->solidComp[i][0] -= weights[k]*gsl_vector_get(deltaSol, k);
                    silminState->solidComp[i+1+j][0] -= weights[k]*gsl_vector_get(deltaSol, k);
                    k++;
                }
                if (silminState->solidComp[i][0] <= 0.0) result = FALSE;
                else result = result & (*solids[i].test)(SIXTH, silminState->T, silminState->P, 0, 0, NULL, NULL, NULL, mSol);
                free(mSol);
            }
        }
    }

    /* undo the correction on failure */
    if (result == FALSE) {
        for (i=0, k=0; i<npc; i++) {
            if (silminState->nSolidCoexist[i]) {
                if (solids[i].na == 1) {
                    silminState->solidComp[i][0] += weights[k]*gsl_vector_get(deltaSol, k);
                    k++;
               } else {
                    for (j=0; j<solids[i].na; j++) {
                        silminState->solidComp[i][0] += weights[k]*gsl_vector_get(deltaSol, k);
                        silminState->solidComp[i+1+j][0] += weights[k]*gsl_vector_get(deltaSol, k);
                        k++;
                    }
                }
            }
        }
    }
    gsl_matrix_free(reducedSolToOx);
    gsl_matrix_free(V); gsl_vector_free(S);
    gsl_vector_free(deltaLiq); gsl_vector_free(deltaSol);
    free(weights);

    return result;
}

#ifdef PHMELTS_ADJUSTMENTS
/***********************************************************************/

/* silmin assumes that the initial guess satisfies the bulk composition
   constraints.  Hence for routines that change the bulk composition or
   monkey with the phases, we need a routine to check that the masses add
   up before we call silmin().  This function considers the bulk comp
   recorded in silminState->bulkComp[] to be gospel and does whatever surgery
   is needed to make the phases add up to it */
void correctXforChangeInBulkComp() {
    //static int beenHere = 0;
    //static double **liqKernel;

    int i, j, k, l, ns, hasLiquid = (silminState->liquidMass != 0.0);
    double *deltaBulkComp = (double *) malloc((size_t) nc*sizeof(double));
    int enoughLiquid;
    int iBulkH2O, iPhaseH2O;

    /* Previously: was obtaining LiqToOx, SVDing it, then backsubstituting with unit
        RHS vectors to get liquid kernels; stored in rows of **liqKernel         */
    /* This gives the same values as bulkSystem[].oxToLiq[] (maybe that was added later?) */

    /* Check no fossil components in solids that have been zeroed in bulk */
    for (i=0; i<npc; i++) {
        for (ns=0; ns<silminState->nSolidCoexist[i]; ns++) {
            if (solids[i].na == 1) {
                for (j=0; j<nc; j++)
	                if (solids[i].solToOx[j] != 0.0 && silminState->bulkComp[j] == 0.0 && silminState->solidComp[i][ns] != 0.0)
	                    silminState->solidComp[i][ns] = 0.0;
            } else {
                for (k=0; k<solids[i].na; k++) {
                    for (j=0; j<nc; j++) {
                        if (solids[i+1+k].solToOx[j] != 0.0 && silminState->bulkComp[j] == 0.0 && silminState->solidComp[i+1+k][ns] != 0.0) {
                            silminState->solidComp[i][ns] -= silminState->solidComp[i+1+k][ns];
                            silminState->solidComp[i+1+k][ns] = 0.0;
                        }
                    }
                }
            }
        }
    }

    /* check no fossil components in liquid */
    for (i=0; i<nlc; i++) {
        for (j=0; j<nc; j++) {
            if (liquid[i].liqToOx[j] != 0.0 && silminState->bulkComp[j] == 0.0) {
                silminState->liquidComp[0][i] = 0.0;
                for (k=1; k<silminState->nLiquidCoexist; k++) {
                    silminState->liquidComp[k][i] = 0.0;
                }
            }
        }
    }
    //if (silminState->nLiquidCoexist == 0) {
    if (!hasLiquid) {
        for (i=0; i<nlc; i++) silminState->liquidComp[0][i] = 0.0;
    }


    /* compute deltaBulkComp, residual between bulk and sum of phases */
    for (i=0; i<nc; i++) {
        deltaBulkComp[i] = -silminState->bulkComp[i];
        if (hasLiquid) {
            for (k=0; k<silminState->nLiquidCoexist; k++)
                for (j=0; j<nlc; j++) deltaBulkComp[i] += silminState->liquidComp[k][j] * (liquid[j].liqToOx)[i];
        }
        for (j=0; j<npc; j++) {
            for (k=0; k<silminState->nSolidCoexist[j]; k++) {
                if (solids[j].na == 1) deltaBulkComp[i] += silminState->solidComp[j][k] * (solids[j].solToOx)[i];
                else for (l=0;l<solids[j].na;l++) deltaBulkComp[i] += silminState->solidComp[j+1+l][k] * (solids[j+1+l].solToOx)[i];
            }
        }
    }

    /* check that H2O component has a phase to go into */
    for (iBulkH2O=0; iBulkH2O<nc; iBulkH2O++) if (!strcmp(bulkSystem[iBulkH2O].label, "H2O")) break;
    for (iPhaseH2O=0; iPhaseH2O<npc; iPhaseH2O++) if (!strcmp(solids[iPhaseH2O].formula, "H2O")) break;

    if(silminState->bulkComp[iBulkH2O] > 0.0) {
        if (fabs(deltaBulkComp[iBulkH2O] + silminState->bulkComp[iBulkH2O]) < 100.0*DBL_EPSILON) {
            silminState->nSolidCoexist[iPhaseH2O] = 1;
            silminState->solidComp[iPhaseH2O][0] = silminState->bulkComp[iBulkH2O];
            deltaBulkComp[iBulkH2O] = 0.0; /* fix double */
        }
    }

    /* we use the local function if first liquid is sufficient else use addOrDropLiquid().
     First criterion is that the change should not reduce any liquid component by
     more than 95%. */
    enoughLiquid = FALSE;
    //if(silminState->nLiquidCoexist) {
    if (hasLiquid) {
        for (i=0,enoughLiquid=TRUE; i<nlc; i++) {
            double dumb;
            for (j=0,dumb=0.0; j<nc; j++) dumb += bulkSystem[j].oxToLiq[i] * deltaBulkComp[j];
            if (dumb > 0.0 && dumb > 0.95*silminState->liquidComp[0][i]) enoughLiquid = FALSE;
        }
    }
    if (enoughLiquid) {
        /* If deltaBulkComp[i] != 0.0, use liqKernel to adjust liquidComp */
        for (i=0; i<nc; i++) {
            if (deltaBulkComp[i] != 0.0) {
                for (j=0; j<nlc; j++) silminState->liquidComp[0][j] -= deltaBulkComp[i] * bulkSystem[i].oxToLiq[j];
            }
        }
        for (i=0; i<nlc; i++) if (silminState->liquidComp[0][i] < 100.0*DBL_EPSILON)
            silminState->liquidComp[0][i] = 0.0;
    } else if (!addOrDropLiquid(deltaBulkComp)) {
        /* failed once -- try again with half of the residual assigned to each method */
        printf("Correct X for Bulk Composition FAILED!\n");
        for (i=0; i<nc; i++) deltaBulkComp[i] /= 2.0;
        enoughLiquid = FALSE;
        //if(silminState->nLiquidCoexist) {
        if (hasLiquid) {
            for (i=0,enoughLiquid=TRUE; i<nlc; i++) {
                double dumb;
                for (j=0,dumb=0.0; j<nc; j++) dumb += bulkSystem[j].oxToLiq[i] * deltaBulkComp[j];
                if (dumb > 0.0 && dumb > 0.95*silminState->liquidComp[0][i]) enoughLiquid = FALSE;
            }
        }
        if (enoughLiquid && addOrDropLiquid(deltaBulkComp)) {
            /* If deltaBulkComp[i] != 0.0, use liqKernel to adjust liquidComp */
            for (i=0; i<nc; i++) {
                if (deltaBulkComp[i] != 0.0) {
                    for (j=0; j<nlc; j++) silminState->liquidComp[0][j] -= deltaBulkComp[i] * bulkSystem[i].oxToLiq[j];
                }
            }
            for (i=0; i<nlc; i++) if (silminState->liquidComp[0][i] < DBL_EPSILON) silminState->liquidComp[0][i] = 0.0;
        } else {
            /* failed again -- dump liquid, run pdaNorm for solid initial guess */

            /* check that H2O component has a phase to go into */
            if(silminState->bulkComp[iBulkH2O] > 0.0) {
                /* may not be exact if within parallel_adiabat */
                if (fabs(silminState->solidComp[iPhaseH2O][0] - silminState->bulkComp[iBulkH2O]) < 100.0*DBL_EPSILON) {
                    silminState->nSolidCoexist[iPhaseH2O] = 0;
                    silminState->solidComp[iPhaseH2O][0] = 0.0;
                    deltaBulkComp[iBulkH2O] = -silminState->bulkComp[iBulkH2O];
                }
            }
            printf("Correct X for Bulk Composition FAILED again - running norm!\n");
            if (!pdaNorm()) {
                /* failed again -- melt system for sure-fire liquid initial guess */
                printf("Correct X for Bulk Composition FAILED three times - melting system!\n");
                for (i=0; i<nlc; i++)
                    for (j=0, silminState->liquidComp[0][i]=0.0; j<nc; j++)
                            silminState->liquidComp[0][i] += bulkSystem[j].oxToLiq[i] * silminState->bulkComp[j];
                for (i=0,silminState->liquidMass=0.0; i<nc; i++)
                    silminState->liquidMass += silminState->bulkComp[i]*bulkSystem[i].mw;
                for (ns=1;ns<silminState->nLiquidCoexist;ns++)
                    for (i=0; i<nlc; i++) silminState->liquidComp[ns][i] = 0.0;
                silminState->nLiquidCoexist = 1;
                for (i=0; i<npc; i++) {
                    for (ns=0;ns<silminState->nSolidCoexist[i];ns++) {
                        silminState->solidComp[i][ns] = 0.0;
                        if (solids[i].na > 1) for (j=0; j<solids[i].na; j++)
                            silminState->solidComp[i+1+j][ns] = 0.0;
                    }
                    silminState->nSolidCoexist[i] = 0;
                }
                silminState->solidMass = 0.0;
                silminState->solidTD.s = 0.0;
            } else {
                if (hasLiquid) {
                    for (j=1; j<silminState->nLiquidCoexist; j++) {
                            free((silminState->liquidComp )[j]); (silminState->liquidComp )[j] = NULL;
                            free((silminState->liquidDelta)[j]); (silminState->liquidDelta)[j] = NULL;
                    }
                    for (i=0;i<nlc;i++) silminState->liquidComp[0][i] = 0.0;
                    //for (i=0;i<nc;i++)  silminState->dspLiquidComp[i] = 0.0;

                    //for (ns=0; ns<silminState->nLiquidCoexist;ns++)
                   // for (i=0; i<nlc; i++) silminState->liquidComp[ns][i] = 0.0;
                }
                //silminState->nLiquidCoexist = 0;
	            silminState->nLiquidCoexist = 1; // else memory leak
                silminState->liquidMass = 0.0;
                multiplyThermoData(&(silminState->liquidTD), 0.0);

            }
        }
    }

    free(deltaBulkComp);

}

/***********************************************************************/
#endif /* PHMELTS_ADJUSTMENTS */

/* end of file SILMIN_SUPPORT.C */
