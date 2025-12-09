const char *silmin_ver(void) { return "$Id: silmin.c,v 1.12 2009/04/24 20:51:09 ghiorso Exp $"; }
/*
 MELTS Source Code: RCS $Log: silmin.c,v $
 MELTS Source Code: RCS Revision 1.9  2008/03/06 17:51:23  ghiorso
 MELTS Source Code: RCS New fluid fractionation mode and other enhancements.
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.8  2007/12/22 22:43:30  ghiorso
 MELTS Source Code: RCS Fixed error in BM integration in Gibbs.c
 MELTS Source Code: RCS Updated param_struct_data.h file for AGU 2007 xMELTS parameters
 MELTS Source Code: RCS Added support for status file production in MELTS-batch (XML)
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.7  2007/11/29 05:32:14  ghiorso
 MELTS Source Code: RCS Majorite testMaj corrections.  Fixed "O2" in sol_struct_data.h to "o2".
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.6  2007/08/23 16:09:39  ghiorso
 MELTS Source Code: RCS Database updates.
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.5  2007/06/08 17:25:43  ghiorso
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
 MELTS Source Code: RCS Revision 1.3  2005/01/08 22:21:02  cvsaccount
 MELTS Source Code: RCS
 MELTS Source Code: RCS Set tolerance in silmin (before HFTI call) to 10*DBL_EPSILON to insure
 MELTS Source Code: RCS catching phase rule violations in simple system crystallization.
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.2  2005/01/08 03:14:02  cvsaccount
 MELTS Source Code: RCS *** empty log message ***
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
 * Revision 3.10  1997/06/21  22:49:27  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.9  1997/05/03  20:23:07  ghiorso
 * *** empty log message ***
 *
 * Revision 3.8  1997/03/27  17:03:11  ghiorso
 * *** empty log message ***
 *
 * Revision 3.7  1996/09/24  20:33:21  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 3.6  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 3.5  1995/11/23  22:37:42  ghiorso
 * Final implementation of subsolidus fO2 buffering.
 *
 * Revision 3.4  1995/11/01  22:40:27  ghiorso
 * Implementation of subsolidus options after Asimow.
 * Additional implementation of nepheline solid solutions.
 *
 * Revision 3.3  1995/09/04  20:01:28  ghiorso
 * Update to allow display of bulk composition (in grams) in the text entry
 * fields of the main silmin display. Liquid composition is no longer
 * display here, and is available only through the popup selection.
 *
 * Revision 3.2  1995/09/01  23:53:03  ghiorso
 * Modifications made to update interface for V3.x and consolidate
 * Graph Widgets
 *
 * Revision 3.1  1995/08/18  19:13:18  ghiorso
 * MELTS Version 3 - Initial Entry
 *
 */

/*
 **++
 **  FACILITY:  Silicate Melts Regression/Crystallization Package
 **
 **  MODULE DESCRIPTION:
 **
 **      Toolkit work proceedure to drive a SILMIN iteration (file: SILMIN.C)
 **
 **  MODIFICATION HISTORY:
 **
 **      V1.0-1  Mark S. Ghiorso  August 24, 1990   Original Version - test case
 **      V1.1-1  Mark S. Ghiorso  April 27, 1991
 **              New solid and liquid include file dependencies added
 **      V1.1-2  Mark S. Ghiorso  August 10, 1991
 **              Initial structure defined.
 **      V1.1-3  Mark S. Ghiorso  September 4, 1991
 **              Added initialization of global structures
 **      V1.1-4  Mark S. Ghiorso  September 6, 1991
 **              Added WorkProcData structure and set active to status on return
 **              (check is only preliminary)
 **      V1.1-5  Mark S. Ghiorso  September 10, 1991
 **              (1) Began implementation of CHANGE_COMPOSITION code and
 **                  CHANGE_TP code
 **              Mark S. Ghiorso  September 12, 1991
 **              (1) Implemented evaluateSaturationState and
 **                  getEqualityConstraints
 **      V1.1-6  Mark S. Ghiorso  September 13, 1991
 **              Altered parameter list to conLiq. Part of development
 **              process for code in getProjGradientAndHessian()
 **      V1.1-7  Mark S. Ghiorso  September 14, 1991
 **              Altered call to (*solids[].convert) to reflect change
 **              in parameter list
 **      V1.1-8  Mark S. Ghiorso  September 18, 1991
 **              (1) Finished preliminary implementation of linearSearch (min1d)
 **                  and getProjGradientAndHessian. Implemented new saturation
 **                  check.
 **              (2) Create update* routines and began to construct output
 **                  to interface routines
 **              (3) Removed definition of wkarea and wksize
 **      V1.1-9  Mark S. Ghiorso  September 24, 1991
 **              Altered parameter list to (*solids[].convert)
 **      V1.1-10 Mark S. Ghiorso  September 28, 1991
 **              Changed treatment of stepsize for linear search procedure
 **              and provided debug check on pseudoRank determination in HFTI
 **      V1.1-11 Mark S. Ghiorso  October 1, 1991
 **              Corrected error in modification of length of search direction
 **              if initial guess to min1d fails
 **      V1.1-12 Mark S. Ghiorso  October 9, 1991
 **              Corrected error in addition of 1-component solid
 **      V1.1-13 Mark S. Ghiorso  October 12, 1991
 **              Added call to checkForCoexistingSolids
 **      V1.1-14 Mark S. Ghiorso  October 15, 1991
 **              (1) Modified references to silminState->solidComp and
 **                  silminState->solidDelta arrays to allow for precipitation
 **                  of immiscible solid phases
 **              (2) Conversion of code to allow for multiple instances of
 **                  solids (silminState->nSolidCoexist)[]
 **      V2.0-1  Mark S. Ghiorso  November 14, 1991
 **              Conversion to OSF Motif V1.1.1/X11 Release 4
 **      V2.0-2  Mark S. Ghiorso  November 23, 1991
 **              (1) Converted method of display for the
 **                  statusEntries[STATUS_ADB_INDEX_STATUS].name widget to
 **                  reflect its change to a scrolled text window
 **      V2.0-3  Mark S. Ghiorso  November 29, 1991
 **              (1) Added numerous informational messages via wprintf
 **              (2) Added additional calls to updateSolidADB to display
 **                  chemical affinities
 **      V2.0-4  Mark S. Ghiorso  December 10, 1991
 **              (1) Corrected indexing problem for loops on npc in linear
 **                  search and phase dropping algorithm
 **              (2) Corrected evaluation of state change fatal error and
 **                  forced a Return of TRUE when constraints are invalid
 **              (3) Removed all references to second magma and magma
 **                  mixing code. (These are now all dealt with by the
 **                  generic assimilation routines)
 **      V2.0-5  Mark S. Ghiorso  December 21, 1991
 **              (1) Removed composition update to silminState->liquidComp,
 **                  transfering the code to checkStateAgainstInterface
 **                  (SILMIN_SUPPORT.C)
 **              (2) Added error handler
 **      V2.1-1  Mark S. Ghiorso  January 2, 1992
 **              (1) Added support for crystal fractionation
 **      V2.1-2  Mark S. Ghiorso  January 6, 1992
 **              (1) Added call to putOutputDataToFile()
 **      V2.1-3  Mark S. Ghiorso  January 10, 1992
 **              (1) Added code for variation in pressure
 **              (2) Added support for assimilation
 **      V2.1-4  Mark S. Ghiorso  February 13, 1992
 **              (1) Added scaling to call to HFTI to prevent overflow
 **      V2.1-5  Mark S. Ghiorso  March 235, 1992
 **              (1) Replaced illegal calls to pow() with macro SQUARE
 **              (2) Removed redundant static qualifier in enum statement
 **      V2.1-6  Mark S. Ghiorso  April 1, 1992
 **              (1) Modified quadratic convergence algorithm to accept
 **                  convergence at sqrt(tau) (non-optimal, but acceptable)
 **                  and fail with a dialog box if itermx is exceeded under
 **                  other conditions
 **      V3.0-1  Mark S. Ghiorso  April 27, 1992
 **              (1) Begin modifications for f O2 buffering
 **                               April 30, 1992 - done
 **      V3.0-2  Mark S. Ghiorso  May 1, 1992
 **              (1) Added computation of reference O2 content, including
 **                  subtraction of O2 during fractionation and addition of
 **                  O2 during assimilation
 **              (2) Cosmetic changes to wprintf 's
 **      V3.0-3  Mark S. Ghiorso  May 2, 1992
 **              (1) Corrected bulkComp vector on exit from linear search
 **                  proceedure to reflect change in Fe2O3/FeO ratio due to
 **                  formation of solid phases
 **      V3.0-4  Mark S. Ghiorso  May 4, 1992
 **              (1) Added a solid cycling algorthm, such that if a phase is
 **                  added, then dropped, it is not considered again until
 **                  the system is updated (i.e. change of T, P, composition).
 **                  This should eliminate infinite looping during constrained
 **                  "path" calculations.
 **      V3.0-5  Mark S. Ghiorso  September 29, 1993
 **              Modified call to realloc to catch zero pointer (SPARC port)
 **      V3.1-1  Mark S. Ghiorso  October 5, 1993
 **              (1) Added code to check for dissapearance of a component from
 **                  the system due to fractionation of a highly incompatible
 **                  element (i.e. Cr from liquid -> spinel)
 **      V4.0-1  Mark S. Ghiorso  May 11, 1994
 **              (1) Added constraint tests for isentropic constraints
 **                               May 30, 1994
 **              (2) Changed silminState->isenthalpic, etc tests to account for
 **                  first pass (refEnthalpy, etc at zero)
 **                               June 6, 1994
 **              (3) Added T and P to reassemble solution dection for
 **                  isenthalpic, isentropic and isochoric constraints
 **                               June 9, 1994
 **              (4) Added T and P updates to linear search section
 **                               June 10, 1994
 **              (5) Modified UPDATE_SYSTEM case for isenthalpic, isentropic
 **                  and isochoric constraints
 **                               June 16, 1994
 **              (6) Corrected error in assimilation update
 **      V5.0-1  Paul D. Asimow  April 14, 1995
 **              Enable liquid-absent operation -- algorithmic changes only;
 **              no subsolidus oxygen buffering
 **      V5.1-1  Paul D. Asimow  August 1, 1995
 **              (1) Enable subsolidus fO2 buffering
 **              (2) New liquid-absent phase add/drop procedure
 **              (3) Insert check that adding phase from liquid will not
 **                  push any liquid component less than zero
 **              (4) Retain best quadratic iteration in case last iteration
 **                  fails convergence test but earlier one passes
 **--
 */

#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

#include "lawson_hanson.h"        /*func decl for Lawson and Hanson routines*/
#include "silmin.h"               /*SILMIN structures include file          */
#include "nash.h"                 /*func decl for Nash routines             */

#ifndef BATCH_VERSION
#include <Xm/Xm.h>
#include <Xm/ToggleBG.h>
#include "interface.h"            /*Specific external declarations          */
#else
#include <gsl/gsl_errno.h>
#include "status.h"               /*Status of calculation in batch mode     */
#endif


#define SQUARE(x) ((x)*(x))
#define REALLOC(x, y) (((x) == NULL) ? malloc(y) : realloc((x), (y)))

#ifdef DEBUG
#undef DEBUG
#endif

#ifdef PRINT_ENERGY_AT_EACH_QUAD_ITERATION
#undef PRINT_ENERGY_AT_EACH_QUAD_ITERATION
#endif

/*
 *=============================================================================
 * Global variables declared extern in SILMIN.H
 */

SilminState     *silminState;
SilminInputData silminInputData = {NULL, NULL};
SilminHistory   *silminHistory;
Constraints     *constraints;

/*
 *=============================================================================
 * Global variables that were initialized in interface.c
 */

#ifdef BATCH_VERSION
#include "status.h"
MeltsStatus meltsStatus;
#endif

int calculationMode = MODE_DEFAULT;
int quad_tol_modifier = 1;

void (*additionalOutput) (char *filename); // not used?
char *addOutputFileName;

/*
 *=============================================================================
 * Global variables needed in rMELTSframework to change calculation modes in
 * webservivces
 */

extern SilminState *bestState;
SilminState *bestState = NULL;

/*
 *=============================================================================
 * Global variables needed in Magma Chamber Simulator / standalone Melts-batch
 * to mimic webservices output
 */

/*
 *=============================================================================
 * Executable code
 */

#ifndef SIG_ERR                /* For VAX C implementations (define BADSIG) */
#define SIG_ERR (int (*)())-1
#endif

#define UPDATE_ALL_ENTRIES -1

#ifndef BATCH_VERSION
#define RELAY_ERROR_COND(string) \
{ \
XmString csString = XmStringCreateLtoR(string, "ISO8859-1"); \
XtVaSetValues(message, XmNmessageString, csString, NULL); \
XtManageChild(message); \
XmStringFree(csString); \
}

/*
 * Static variables and functions designed to capture arithmetic errors and
 *   return control as gracefully as possible.
 */

static void newErrorHandler(int sig);  /* new error handler function */
static void (*oldErrorHandler)();      /* old error handler function */
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

static int *permutation_alloc(int n) {
    int *p = (int *) malloc((size_t) n*sizeof(int));
    return p;
}

static void permutation_free(int *p, int n) {
    free(p);
}

#ifndef BATCH_VERSION
Boolean silmin(XtPointer client_data)
#else
int silmin(void)
#endif /* BATCH_VERSION */
{
    enum steps {
        CHANGE_COMPOSITION,      CHANGE_TP,               CHECK_SATURATION,
        ADD_PHASE,               PROJECT_CONSTRAINTS,     PRE_QUADRATIC,
        CONSTRUCT_QUADRATIC,     SOLVE_QUADRATIC,         REASSEMBLE_SOLUTION,
        LINEAR_SEARCH,           DROP_PHASE,              CONVERGENCE_TEST,
        VERIFY_SATURATION,       OUTPUT_RESULTS,          UPDATE_SYSTEM
    };
    static int curStep = 0;
    enum stages {
        PRE_STAGE_ZERO,
#ifdef PHMELTS_ADJUSTMENTS
        L_X_STAGE_ONE, L_X_STAGE_TWO, L_X_STAGE_THREE,
#endif
        L_H_STAGE_ONE, L_H_STAGE_TWO, L_H_STAGE_THREE,
        L_S_STAGE_ONE, L_S_STAGE_TWO, L_S_STAGE_THREE,
        L_V_STAGE_ONE, L_V_STAGE_TWO, L_V_STAGE_THREE
    };
    static int curStage = 0;
    static int hasSupersaturation, conCols, conRows, iterQuad;
    static double **cMatrix, *hVector, *dVector, *yVector;
    static double **eMatrix, **bMatrix, rNorm, sNorm;
    double mTotal;
    int i, j, k, nl, ns, stateChange, hasNlCon;
    int hasLiquid = ((silminState != NULL) && (silminState->liquidMass != 0.0));
    static double bestrNorm;
    static int acceptable = FALSE, bestIter, hessianType = HESSIAN_TYPE_NORMAL;
/* Does this get initialised outside the rMELTS framework?
    static SilminState *bestState = NULL; */
#ifdef PHMELTS_ADJUSTMENTS
    int heteroFlag, iLiqH2O = 0, iBulkH2O = 0, iPhaseH2O = 0;
    int H2Obuffer = silminState->H2Obuffer;
    double *saveT, *saveP, *saveLiqComp;

    /* use bestState as thread-safe storage */
    /* works because allocSilminStatePointer() uses calloc, not malloc */
    if (bestState == (SilminState *) NULL) {
        bestState = allocSilminStatePointer();
    }
    saveT = &(bestState->T);
    saveP = &(bestState->P);
    saveLiqComp = bestState->liquidComp[0];

    // pMELTS and rholite-MELTS 1.0.2 only
    for (iLiqH2O=0; iLiqH2O<nlc; iLiqH2O++) if (!strcmp(liquid[iLiqH2O].label, "H2O")) break;
    for (iBulkH2O=0; iBulkH2O<nc; iBulkH2O++) if (!strcmp(bulkSystem[iBulkH2O].label, "H2O")) break;
    for (iPhaseH2O=0; iPhaseH2O<npc; iPhaseH2O++) if (!strcmp(solids[iPhaseH2O].label, "fluid")) break;
#endif

#ifndef BATCH_VERSION
    /*WorkProcData *workProcData = (WorkProcData *) client_data;*/

    /******************************************************************************
   * On entry, check status of calculation and check workProcData.mode for tag:
   *   FALSE   return call without interface modification
   *   TRUE    initial call to invoke a crystallization run
   ******************************************************************************/

    /*if (workProcData->mode) {
        } */

    /******************************************************************************
   * Step to current phase of the calculation
   ******************************************************************************/

    /*updateStatusADB(STATUS_ADB_INDEX_PHASE, &curStep);*/
#else /* BATCH_VERSION */
    (void) gsl_set_error_handler_off();
    if (curStep == 0) curStep = CHANGE_COMPOSITION;
#endif /* BATCH_VERSION */

    if (curStage == 0) curStage = PRE_STAGE_ZERO;

    switch(curStep) {
            /* ======================================================================== */
        case CHANGE_COMPOSITION:

            updateBulkADB();
#ifndef BATCH_VERSION
            updateStatusADB(STATUS_ADB_INDEX_MASS_LIQUID, &(silminState->liquidMass));
            if(silminState->fo2Path == FO2_NONE) updateStatusADB(STATUS_ADB_INDEX_LOGFO2, &(silminState->fo2));

            workProcData->active = TRUE;
#endif

            curStep++;
            return FALSE;
            /* ------------------------------------------------------------------------ */
        case CHANGE_TP:

            /* -> Calculate liquid end-member properties                                  */
            for (i=0; i<nlc; i++)
                gibbs(silminState->T, silminState->P, (char *) liquid[i].label, &(liquid[i].ref), &(liquid[i].liq), &(liquid[i].fus), &(liquid[i].cur));

            /* -> Calculate solid  end-member properties                                  */
            for (i=0, j=0; i<npc; i++) {
                if (solids[i].type == PHASE) {
                    if ((silminState->incSolids)[j]) {
                        if(solids[i].na == 1)
                            gibbs(silminState->T, silminState->P, (char *) solids[i].label, &(solids[i].ref), NULL, NULL, &(solids[i].cur));
                        else {
                            for (k=0; k<solids[i].na; k++) {
                                gibbs(silminState->T, silminState->P, (char *) solids[i+1+k].label, &(solids[i+1+k].ref), NULL, NULL, &(solids[i+1+k].cur));
                            }
                            i += solids[i].na;
                        }
                    }
                    j++;
                }
            }

            /* -> Calculate O2 end-member properties if path is buffered                  */
            if (silminState->fo2Path != FO2_NONE) gibbs(silminState->T, silminState->P, "o2", &(oxygen.ref), NULL, NULL, &(oxygen.cur));


            /* -> Redistribute Fe2O3 and FeO in liquid phase to establish buffer at this T and P                                                            */
            if (silminState->fo2Path != FO2_NONE && hasLiquid) {
                double *moles = (double *) malloc((size_t) nc*sizeof(double));
                silminState->fo2 = getlog10fo2(silminState->T, silminState->P, silminState->fo2Path);
                for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                    for (i=0; i<nc; i++) {
                        for (j=0, moles[i]=0.0; j<nlc; j++) moles[i] += (silminState->liquidComp)[nl][j]*(liquid[j].liqToOx)[i];
                        (silminState->bulkComp)[i] -= moles[i];
                    }
                    conLiq(FIRST | SEVENTH, FIRST, silminState->T, silminState->P, moles, NULL, NULL, NULL, NULL, NULL, &(silminState->fo2));
                    for (i=0; i<nc; i++) (silminState->bulkComp)[i] += moles[i];
                    for (i=0; i<nlc; i++) for (j=0, (silminState->liquidComp)[nl][i]=0.0; j<nc; j++) (silminState->liquidComp)[nl][i] += moles[j]*(bulkSystem[j].oxToLiq)[i];
                }
                free(moles);

                if (silminState->oxygen == 0.0) {
                    /* This block should really be before the buffer is imposed, but makes no difference for alphaMELTS or Melts-batch software */
                    for (i=0, silminState->oxygen=0.0; i<nlc; i++) for (nl=0; nl<silminState->nLiquidCoexist; nl++)
                        silminState->oxygen += (oxygen.liqToOx)[i]*(silminState->liquidComp)[nl][i];
                    for (i=0; i<npc; i++) for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
                        if (solids[i].na == 1) silminState->oxygen += (oxygen.solToOx)[i]*(silminState->solidComp)[i][ns];
                        else {
                            for (j=0; j<solids[i].na; j++) silminState->oxygen += (oxygen.solToOx)[i+1+j]*(silminState->solidComp)[i+1+j][ns];
                        }
                    }
                }

                updateBulkADB();
#ifndef BATCH_VERSION
                updateStatusADB(STATUS_ADB_INDEX_MASS_LIQUID, &(silminState->liquidMass));
                updateStatusADB(STATUS_ADB_INDEX_LOGFO2, &(silminState->fo2));
#endif

                /* -> If no liquid, find a subsolidus buffer reaction, run it to establish buffer at this T and P                                       */
            } else if (silminState->fo2Path != FO2_NONE && !hasLiquid) {
                double muO2;
                silminState->fo2 = getlog10fo2(silminState->T, silminState->P, silminState->fo2Path);
                muO2 = silminState->fo2*(R*silminState->T*log(10.0));
                if (silminState->oxygen == 0.0) /* Atomic Number of O is 8; we want O2 */
                    for (i=0, silminState->oxygen=0.0; i<nc; i++) silminState->oxygen += (silminState->bulkComp)[i]*(bulkSystem[i].oxToElm)[8]/2.0;

#ifdef DEBUG
                printf("\nBefore entry to subsolidusmuO2:\n");
                for (i=0; i<npc; i++) {
                    for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
                        printf("soln %-15.15s\n", solids[i].label);
                        if (solids[i].na == 1) printf("soln %-15.15s = %13.6g\n", "Total moles", (silminState->solidComp)[i][ns]);
                        else {
                            for (k=0, mTotal=0.0; k<solids[i].na; k++) {
                                if ((silminState->solidComp)[i+1+k][ns] != 0.0) {
                                    printf("soln %-15.15s = %13.6g\n", solids[i+1+k].label, (silminState->solidComp)[i+1+k][ns]);
                                }
                            }
                            printf("soln %-15.15s = %13.6g\n", "Total moles", (silminState->solidComp)[i][ns]);
                        }
                    }
                }
                printf("oxygen content = %13.6g\n", silminState->oxygen);
#endif /* DEBUG */

                if (!subsolidusmuO2(0, &muO2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)) {
                    printf("Failure to impose fO2 buffer in subsolidus. Releasing buffer constraint from the system.\n");
#ifdef PHMELTS_ADJUSTMENTS
            silminState->fo2Liq = silminState->fo2Path;
#endif
                    silminState->fo2Path = FO2_NONE;
#ifndef BATCH_VERSION
                    XmToggleButtonGadgetSetState(tg_path_none, True, True);
#endif
                }

#ifdef DEBUG
                printf("\nAfter exit from subsolidusmuO2:\n");
                for (i=0; i<npc; i++) {
                    for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
                        printf("soln %-15.15s\n", solids[i].label);
                        if (solids[i].na == 1) printf("soln %-15.15s = %13.6g\n", "Total moles", (silminState->solidComp)[i][ns]);
                        else {
                            for (k=0, mTotal=0.0; k<solids[i].na; k++) {
                                if ((silminState->solidComp)[i+1+k][ns] != 0.0) {
                                    printf("soln %-15.15s = %13.6g\n", solids[i+1+k].label, (silminState->solidComp)[i+1+k][ns]);
                                }
                            }
                            printf("soln %-15.15s = %13.6g\n", "Total moles", (silminState->solidComp)[i][ns]);
                        }
                    }
                }
                printf("oxygen content = %13.6g\n", silminState->oxygen);
#endif /* DEBUG */

                updateBulkADB();
#ifndef BATCH_VERSION
                updateStatusADB(STATUS_ADB_INDEX_LOGFO2, &(silminState->fo2));
#endif
            }

#ifndef BATCH_VERSION
            wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...Checking saturation state of potential solids.\n");
            workProcData->active = TRUE;
#elif defined(ALPHAMELTS_UPDATE_SYSTEM)
            printf("...Checking saturation state of potential solids.\n");
#else
            fprintf(stderr, "...Checking saturation state of potential solids.\n");
#endif

            curStep++;
            return FALSE;
            /* ------------------------------------------------------------------------ */
        case CHECK_SATURATION:

            /*  Note: This routine uses the array (silminState->ySol)[]. It returns a result in
                the the first npc elements. For liquids it uses (silminState->yLiq)[] and returns
                nlc elements; (silminState->yLiq)[nlc-1] is the liquid affinity.

                Note: The storage mode of the array silminState->solidComp obeys the
                following convention:

                if solids[i].type == PHASE then
                (silminState->solidComp)[i][*] = total moles of the phase in the system

                if solids[i].type == PHASE and solids[i].na > 1,
                i.e. an solids[i].na-component solid solution, then
                (silminState->solidComp)[i+1][*] through
                (silminState->solidComp)[i+solids[i].na][*] contain the composition of
                the solid solution in terms of moles of the endmember components

                For all solids if solids[i].type == PHASE, then a non-zero value in
                (silminState->solidComp)[i][*] means the solid phase is present. This
                test is used in evaluateSaturationState()                            */

            if ((silminState->ySol) == NULL) {
                (silminState->ySol) = (double *) malloc((size_t) npc*sizeof(double));
                (silminState->yLiq) = (double *) malloc((size_t) nlc*sizeof(double));
            }
            for (i=0; i<npc; i++) (silminState->cylSolids)[i] = 0;

#ifdef PHMELTS_ADJUSTMENTS
            heteroFlag = hasLiquid;

            if (silminState->H2Obuffer && (silminState->aH2O != 0.0) && hasLiquid) {
                heteroFlag = heteroFlag & (silminState->T == *saveT);
                heteroFlag = heteroFlag & (silminState->P == *saveP);
                for (i=0; i<nlc ;i++) heteroFlag = heteroFlag & (silminState->liquidComp[0][i] == saveLiqComp[i]);
                for (i=0; i<nlc; i++) saveLiqComp[i] = 0.0; /* in case we exit early */
            }

            /* Only call at this stage if we are starting from (at least metastable) heterogeneous equilibrium */
            if (heteroFlag) hasSupersaturation = evaluateSaturationState((silminState->ySol), (silminState->yLiq));
#else
            /* Only call at this stage if we are starting from liquid */
            if (hasLiquid) hasSupersaturation = evaluateSaturationState((silminState->ySol), (silminState->yLiq));
#endif
            else           hasSupersaturation = FALSE;

#ifndef BATCH_VERSION
            updateSolidADB((silminState->ySol), (silminState->yLiq));
            workProcData->active = TRUE;
#endif

            curStep++;
            return FALSE;
            /* ------------------------------------------------------------------------ */
        case ADD_PHASE:

            /* Note: The storage (silminState->ySol)[0:npc-1] must contain results set by the
            previous call to evaluateSaturationState. The space may be
            reused upon exit from this case.  If liquid is absent, (silminState->yLiq)[0:nlc-1]
            should contain results for liquid.  SilminState->incSolids[npc] is
            used to store inclusion or suppression of liquid phase.          */

            if (hasSupersaturation) {
                double minAffinity = 0.0;
                int    index = -9999;
                for (i=0; i<npc; i++) {
                    if (solids[i].type == PHASE && (silminState->ySol)[i] < 0.0) {
                        if ((silminState->ySol)[i] < minAffinity) { minAffinity = (silminState->ySol)[i]; index = i; }
                    }
                }
#ifdef PHMELTS_ADJUSTMENTS
                if ((silminState->yLiq)[nlc-1] < minAffinity && silminState->incLiquids)
#else
                if ((silminState->yLiq)[nlc-1] < minAffinity && silminState->incSolids[npc])
#endif
                { minAffinity = (silminState->yLiq)[nlc-1]; index = -1; }
                if (index >= 0) {
#ifdef PHMELTS_ADJUSTMENTS
                    //		    double inmass = MASSIN*silminState->refMass/100.0;
                    double inmass = MASSIN;
#else
                    double inmass = MASSIN;
#endif
                    int acceptable = FALSE;

#ifndef BATCH_VERSION
                    wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...Adding the solid phase %s to the assemblage.\n", solids[index].label);
#elif defined(ALPHAMELTS_UPDATE_SYSTEM)
                    printf("...Adding the solid phase %s to the assemblage.\n", solids[index].label);
#else
                    fprintf(stderr, "...Adding the solid phase %s to the assemblage.\n", solids[index].label);
#endif

                    /* test to see how much we can safely add - use the first liquid as a test case */
                    if (hasLiquid) {
                        while (!acceptable) {
                            double *dummyComp = (double *) calloc((size_t) nlc, sizeof(double));
                            for (i=0;i<nlc;i++) dummyComp[i] = silminState->liquidComp[0][i];
                            acceptable = TRUE;
                            if (solids[index].na == 1) for (j=0; j<nlc; j++) dummyComp[j]-=(solids[index].solToLiq)[j]*inmass;
                            else {
                                int na = solids[index].na;
                                double *mSol = (double *) malloc((size_t) na*sizeof(double));
                                (*solids[index].convert)(THIRD, FOURTH, silminState->T, silminState->P, NULL, NULL,&(silminState->ySol)[index+1],mSol, NULL, NULL,  NULL,  NULL);
                                for (i=0; i<na; i++) for (j=0; j<nlc; j++) dummyComp[j] -= (solids[index+1+i].solToLiq)[j]*mSol[i]*inmass;
                                free(mSol);
                            }
                            // line removed and replaced with the one below it to deal with new liquid speciation model
                            // for (j=0;j<nlc;j++) if (dummyComp[j] < 0.0) { acceptable = FALSE; printf("...addPhase: inmass reduced, %s < 0\n", liquid[j].label); }
                            acceptable = testLiq(SIXTH, silminState->T, silminState->P, 0, 0, NULL, NULL, NULL, dummyComp);
                            if (acceptable == FALSE) {
                                inmass *= 0.5;
#ifdef BATCH_VERSION
                                if (inmass < 10.0*DBL_EPSILON) {
                                    meltsStatus.status = SILMIN_ADD_LIQUID_1;
                                    curStage = 0;
                                    curStep = 0;
                                    hasSupersaturation = 0;
                                    return TRUE;
                                }
#endif
                            }
                            free(dummyComp);
                        }
                    }
#ifdef DEBUG
                    printf("...addPhase: inmass(a) = %g\n", inmass);
#endif
                    acceptable = FALSE;
                    while (!acceptable) {
                        double *deltaBulkComp = (double *) calloc((size_t) nc,sizeof(double));
                        acceptable = TRUE;
                        (silminState->solidComp)[index][0] = inmass;
                        if (solids[index].na == 1) {
                            if (hasLiquid) for (j=0; j<nlc; j++) (silminState->liquidComp)[0][j] -= (solids[index].solToLiq)[j] * (silminState->solidComp)[index][0];
                            else {
                                for (j=0; j<nc; j++) deltaBulkComp[j] = (solids[index].solToOx)[j] * (silminState->solidComp)[index][0];
#ifndef PHMELTS_ADJUSTMENTS
                                addOrDropLiquid(deltaBulkComp);
#endif
                            }
                        } else {
                            int na = solids[index].na;
                            double *mSol = (double *) malloc((size_t) na*sizeof(double));
                            (*solids[index].convert)(THIRD, FOURTH, silminState->T, silminState->P, NULL, NULL, &(silminState->ySol)[index+1], mSol, NULL, NULL,  NULL, NULL);
                            for (i=0; i<na; i++) {
                                (silminState->solidComp)[index+1+i][0] = mSol[i]*(silminState->solidComp)[index][0];
                                if (hasLiquid)
                                    for (j=0; j<nlc; j++) (silminState->liquidComp)[0][j] -= (solids[index+1+i].solToLiq)[j] * (silminState->solidComp)[index+1+i][0];
                                else
                                    for (j=0; j<nc; j++) deltaBulkComp[j] += (solids[index+1+i].solToOx)[j] * (silminState->solidComp)[index+1+i][0];
                            }
#ifndef PHMELTS_ADJUSTMENTS
                            if (!hasLiquid) acceptable = addOrDropLiquid(deltaBulkComp);
#endif
                            free (mSol);
                        }
#ifdef PHMELTS_ADJUSTMENTS
                        /* modification 4/11/02 PDA: allow saturation with hydrous phases other than H2O
                        even if H2O is absent from bulk composition if aH2O buffered */
                        if (!hasLiquid) {
                            if (silminState->H2Obuffer && (silminState->aH2O != 0.0) &&
                            (silminState->bulkComp[iBulkH2O] == 0.0) && (deltaBulkComp[iBulkH2O] != 0.0)) {
                                double saveDeltaH2O = deltaBulkComp[iBulkH2O];

                                deltaBulkComp[iBulkH2O] = 0.0;
                                acceptable = addOrDropLiquid(deltaBulkComp);

                                if(acceptable) {

                            silminState->bulkComp[iBulkH2O] = saveDeltaH2O;

                            /* set water in NAM to have same specific entropy as vapour */
                            if     (silminState->isentropic  && (silminState->refEntropy  != 0.0))
                                silminState->refEntropy  += silminState->bulkComp[iBulkH2O]*(solids[iPhaseH2O].cur).s;
                            else if(silminState->isenthalpic && (silminState->refEnthalpy != 0.0))
                                silminState->refEnthalpy += silminState->bulkComp[iBulkH2O]*(solids[iPhaseH2O].cur).h;
                            else if(silminState->isochoric   && (silminState->refVolume   != 0.0))
                                silminState->refVolume   += silminState->bulkComp[iBulkH2O]*(solids[iPhaseH2O].cur).v;

                            if(silminState->fo2Path != FO2_NONE && silminState->oxygen != 0.0)
                                silminState->oxygen += silminState->bulkComp[iBulkH2O]*(oxygen.liqToOx)[iLiqH2O];
                                }

                            } else {
                            acceptable = addOrDropLiquid(deltaBulkComp);
                            }
                        }
#else
                        if (!hasLiquid) acceptable = acceptable & spinodeTest();
#endif
                        if (acceptable == FALSE) {
                            inmass /= 2.0; /* note this only happens without liquid */
#ifdef BATCH_VERSION
                            if (inmass < 10.0*DBL_EPSILON) {
                                meltsStatus.status = SILMIN_ADD_LIQUID_2;
                                curStage = 0;
                                curStep = 0;
                                hasSupersaturation = 0;
                                return TRUE;
                            }
#endif
                        }
                        free(deltaBulkComp);
                    }
#ifdef DEBUG
                    printf("...addPhase: inmass(b) = %g\n", inmass);
#endif
                    (silminState->nSolidCoexist)[index] = 1;

                } else {  /* adding liquid to subsolidus assemblage */
#ifdef PHMELTS_ADJUSTMENTS
                    //		    double inmass = MASSIN*silminState->refMass/100.0;
                    double inmass = MASSIN;
#else
                    double inmass = MASSIN;
#endif
                    int acceptable = FALSE;
                    double *mLiq = (double *) malloc((size_t) nlc*sizeof(double));
#ifndef BATCH_VERSION
                    wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...Adding liquid to the assemblage.\n");
#elif defined(ALPHAMELTS_UPDATE_SYSTEM)
                    printf("...Adding liquid to the assemblage.\n");
#else
                    fprintf(stderr, "...Adding liquid to the assemblage.\n");
#endif
                    conLiq(THIRD, FOURTH, silminState->T, silminState->P, NULL, NULL, (silminState->yLiq), mLiq, NULL, NULL, NULL);
                    acceptable = FALSE;
                    while (!acceptable) {
                        double *deltaBulkComp = (double *) calloc((size_t) nc,sizeof(double));
                        acceptable = TRUE;
                        hasLiquid = TRUE; silminState->liquidMass = inmass;
                        for (i=0;i<nlc;i++) {
                            (silminState->liquidComp)[0][i] = mLiq[i]*inmass;
                            for (j=0;j<nc;j++) deltaBulkComp[j] += (liquid[i].liqToOx)[j] * (silminState->liquidComp)[0][i];
                        }
#ifdef PHMELTS_ADJUSTMENTS
                        /* modification 4/11/02 PDA: allow saturation with hydrous phases other than H2O
                        even if H2O is absent from bulk composition if aH2O buffered */
                        if (silminState->H2Obuffer && (silminState->aH2O != 0.0) &&
                            (silminState->bulkComp[iBulkH2O] == 0.0) && (deltaBulkComp[iBulkH2O] != 0.0)) {
                            double saveDeltaH2O = deltaBulkComp[iBulkH2O];

                            deltaBulkComp[iBulkH2O] = 0.0;
                            acceptable = addOrDropLiquid(deltaBulkComp);

                            if(acceptable) {

                                silminState->bulkComp[iBulkH2O] = saveDeltaH2O;

                                /* set water in NAM to have same specific entropy as vapour */
                                if     (silminState->isentropic  && (silminState->refEntropy  != 0.0))
                                    silminState->refEntropy  += silminState->bulkComp[iBulkH2O]*(solids[iPhaseH2O].cur).s;
                                else if(silminState->isenthalpic && (silminState->refEnthalpy != 0.0))
                                    silminState->refEnthalpy += silminState->bulkComp[iBulkH2O]*(solids[iPhaseH2O].cur).h;
                                else if(silminState->isochoric   && (silminState->refVolume   != 0.0))
                                    silminState->refVolume   += silminState->bulkComp[iBulkH2O]*(solids[iPhaseH2O].cur).v;

                                if(silminState->fo2Path != FO2_NONE && silminState->oxygen != 0.0)
                                    silminState->oxygen += silminState->bulkComp[iBulkH2O]*(oxygen.liqToOx)[iLiqH2O];
                            }

                        }
                        else {
                            acceptable = addOrDropLiquid(deltaBulkComp);
                        }
#else
                        acceptable = addOrDropLiquid(deltaBulkComp);
#endif
                        free(deltaBulkComp);
                        if (acceptable == FALSE) {
                            inmass /= 2.0;
#ifdef BATCH_VERSION
                            if (inmass < 10.0*DBL_EPSILON) {
                                meltsStatus.status = SILMIN_ADD_LIQUID_3;
                                curStage = 0;
                                curStep = 0;
                                hasSupersaturation = 0;
                                return TRUE;
                            }
#endif
                        }
                    }
#ifdef PHMELTS_ADJUSTMENTS
                    if (silminState->fo2Liq != FO2_NONE) silminState->fo2Path = silminState->fo2Liq;

                    /* Note: if initially above anhydrous solidus in pHMELTS then H2Obuffer won't be on yet. */
                    /* Note: on josquin was this just
                       if (silminState->H2Obuffer && silminState->isentropic) silminState->H2Obuffer = FALSE; */

                    /* Experimental - turn isentropic etc. on and off in the same way as with fO2 buffer,
                       unless using Kress & Carmichael 1991. */
                    if (silminState->fo2Path != FO2_NONE && !silminState->fo2Alt) {
                        if (silminState->H2Obuffer && (silminState->aH2O != 0.0) && hasLiquid &&
                            ((silminState->isenthalpic && silminState->refEnthalpy != 0.0) ||
                            (silminState->isentropic  && silminState->refEntropy  != 0.0) ||
                            (silminState->isochoric   && silminState->refVolume   != 0.0))) silminState->H2Obuffer = FALSE;
                    }
                    H2Obuffer = silminState->H2Obuffer;
#endif
                    silminState->nLiquidCoexist = 1;
                    free (mLiq);
                }

                if (silminState->fo2Path != FO2_NONE && hasLiquid) {
                    double *moles = (double *) malloc((size_t) nc*sizeof(double));
                    for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                        for (i=0; i<nc; i++) {
                            for (j=0, moles[i]=0.0; j<nlc; j++) moles[i] += (silminState->liquidComp)[nl][j]*(liquid[j].liqToOx)[i];
                            (silminState->bulkComp)[i] -= moles[i];
                        }
                        conLiq(FIRST | SEVENTH, FIRST, silminState->T, silminState->P, moles, NULL, NULL, NULL, NULL, NULL, &(silminState->fo2));
                        for (i=0; i<nc; i++) (silminState->bulkComp)[i] += moles[i];
                        for (i=0; i<nlc; i++) for (j=0, (silminState->liquidComp)[nl][i]=0.0; j<nc; j++) (silminState->liquidComp)[nl][i] += moles[j]*(bulkSystem[j].oxToLiq)[i];
                    }
                    free(moles);
                } else if (silminState->fo2Path != FO2_NONE && !hasLiquid) {
                    double muO2 = silminState->fo2*(R*silminState->T*log(10.0));
                    if (!subsolidusmuO2(0, &muO2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)) {
                        printf("Failure to impose fO2 buffer in subsolidus.  Releasing buffer constraint from the system.\n");
#ifdef PHMELTS_ADJUSTMENTS
                        silminState->fo2Liq = silminState->fo2Path;
#endif
                        silminState->fo2Path = FO2_NONE;
#ifndef BATCH_VERSION
                        XmToggleButtonGadgetSetState(tg_path_none, True, True);
#endif
                    }
                }
            }

            iterQuad = 0; bestrNorm = 1.0; acceptable = FALSE;

#ifndef BATCH_VERSION
            updateStatusADB(STATUS_ADB_INDEX_QUADRATIC, &iterQuad);
            workProcData->active = TRUE;
#endif

            curStep++;
            return FALSE;
            /* ------------------------------------------------------------------------ */
        case PROJECT_CONSTRAINTS:

            /*  Note: Storage for cMatrix, hVector, dVector and yVector is allocated
                within getEqualityConstraints. Hence, they are passed as pointers.
                The sizes of the resulting arrays are returned in conRows and conCols  */

#ifdef PHMELTS_ADJUSTMENTS
            /*  Buffering fO2 with Kress & Carmichael 1991 and nonlinear constraints does not work, even with the V-fO2 fix
                (see linear_search_H2O.c and gradient_hessian_H2O.c), so the constraints are turned on and off, here and below.
                Buffering fO2 with subsolidusmuo2 (including fo2Alt) and nonlinear constraints is correct if there is no V-fO2 fix.
                Buffering H2O with non-linear constraints is not fully implemented so we 'borrow' the fO2 switch needed for K&C. */

            if      (   (silminState->isenthalpic && (silminState->refEnthalpy != 0.0)  ) && (H2Obuffer && (silminState->aH2O != 0.0) && hasLiquid)
                && (curStage == PRE_STAGE_ZERO) ) { curStage = L_H_STAGE_ONE; silminState->isenthalpic = FALSE; }
            else if (     (silminState->isentropic  && (silminState->refEntropy  != 0.0)) && (H2Obuffer && (silminState->aH2O != 0.0) && hasLiquid)
                && (curStage == PRE_STAGE_ZERO) ) { curStage = L_S_STAGE_ONE; silminState->isentropic  = FALSE; }
            else if (   (silminState->isochoric   && (silminState->refVolume   != 0.0)  ) && (H2Obuffer && (silminState->aH2O != 0.0) && hasLiquid)
                && (curStage == PRE_STAGE_ZERO) ) { curStage = L_V_STAGE_ONE; silminState->isochoric   = FALSE; }
            else if (silminState->fo2Iter && !silminState->fo2Alt && hasLiquid) {
                if      (   (silminState->isenthalpic && (silminState->refEnthalpy != 0.0)  ) && (silminState->fo2Path != FO2_NONE)
                            && (curStage == PRE_STAGE_ZERO) ) { curStage = L_H_STAGE_ONE; }
                else if (     (silminState->isentropic  && (silminState->refEntropy  != 0.0)) && (silminState->fo2Path != FO2_NONE)
                            && (curStage == PRE_STAGE_ZERO) ) { curStage = L_S_STAGE_ONE; }
                else if (   (silminState->isochoric   && (silminState->refVolume   != 0.0)  ) && (silminState->fo2Path != FO2_NONE)
                            && (curStage == PRE_STAGE_ZERO) ) { curStage = L_V_STAGE_ONE; }
                else if (   (silminState->fo2Path != FO2_NONE)
                            && (curStage == PRE_STAGE_ZERO) ) { curStage = L_X_STAGE_ONE; }
            } else
#endif
            if      (   (silminState->isenthalpic && (silminState->refEnthalpy != 0.0)  ) && (silminState->fo2Path != FO2_NONE)
               && (curStage == PRE_STAGE_ZERO) ) { curStage = L_H_STAGE_ONE; silminState->isenthalpic = FALSE; }
            else if (     (silminState->isentropic  && (silminState->refEntropy  != 0.0)) && (silminState->fo2Path != FO2_NONE)
               && (curStage == PRE_STAGE_ZERO) ) { curStage = L_S_STAGE_ONE; silminState->isentropic  = FALSE; }
            else if (   (silminState->isochoric   && (silminState->refVolume   != 0.0)  ) && (silminState->fo2Path != FO2_NONE)
               && (curStage == PRE_STAGE_ZERO) ) { curStage = L_V_STAGE_ONE; silminState->isochoric   = FALSE; }

#ifndef BATCH_VERSION
            if (iterQuad == 0) wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...Projecting equality constraints.\n");
#elif defined(ALPHAMELTS_UPDATE_SYSTEM)
            if (iterQuad == 0) printf("...Projecting equality constraints.\n");
#else
            if (iterQuad == 0) fprintf(stderr, "...Projecting equality constraints.\n");
#endif
#ifdef DEBUG
            printf("\nMaking call to getEqualityConstraints(...) with curStage = ");
            if      (curStage == PRE_STAGE_ZERO ) printf("PRE_STAGE_ZERO");
#ifdef PHMELTS_ADJUSTMENTS
            else if (curStage == L_X_STAGE_ONE  ) printf("L_X_STAGE_ONE");
            else if (curStage == L_X_STAGE_TWO  ) printf("L_X_STAGE_TWO");
            else if (curStage == L_X_STAGE_THREE) printf("L_X_STAGE_THREE");
#endif
            else if (curStage == L_H_STAGE_ONE  ) printf("L_H_STAGE_ONE");
            else if (curStage == L_H_STAGE_TWO  ) printf("L_H_STAGE_TWO");
            else if (curStage == L_H_STAGE_THREE) printf("L_H_STAGE_THREE");
            else if (curStage == L_S_STAGE_ONE  ) printf("L_S_STAGE_ONE");
            else if (curStage == L_S_STAGE_TWO  ) printf("L_S_STAGE_TWO");
            else if (curStage == L_S_STAGE_THREE) printf("L_S_STAGE_THREE");
            else if (curStage == L_V_STAGE_ONE  ) printf("L_V_STAGE_ONE");
            else if (curStage == L_V_STAGE_TWO  ) printf("L_V_STAGE_TWO");
            else if (curStage == L_V_STAGE_THREE) printf("L_V_STAGE_THREE");
            printf(". T = %g, P = %g\n", silminState->T-273.15, silminState->P);
#endif
            getEqualityConstraints(&conRows, &conCols, &cMatrix, &hVector, &dVector, &yVector);

#ifndef BATCH_VERSION
            workProcData->active = TRUE;
#endif

            curStep++;
            return FALSE;
            /* ------------------------------------------------------------------------ */
        case PRE_QUADRATIC:

#ifndef BATCH_VERSION
            if (iterQuad == 0) wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...Minimizing the thermodynamic potential.\n");
            workProcData->active = TRUE;
#elif defined(ALPHAMELTS_UPDATE_SYSTEM)
            if (iterQuad == 0) printf("...Minimizing the thermodynamic potential.\n");
#else
            if (iterQuad == 0) fprintf(stderr, "...Minimizing the thermodynamic potential.\n");
#endif

            curStep++;
            return FALSE;
            /* ------------------------------------------------------------------------ */
        case CONSTRUCT_QUADRATIC:

            /* Note: Storage for eMatrix, and bMatrix is allocated within getProjGradientAndHessian. Hence, they are passed as pointers. */

            iterQuad++;
#ifndef BATCH_VERSION
            updateStatusADB(STATUS_ADB_INDEX_QUADRATIC, &iterQuad);
#endif

#ifdef DEBUG
            printf("\nMaking call to getProjGradientAndHessian(...) with conRows = %d and conCols = %d\n", conRows, conCols);
#endif
            hessianType = getProjGradientAndHessian(conRows, conCols, &eMatrix, &bMatrix, cMatrix, hVector, dVector, yVector);

#ifndef BATCH_VERSION
            wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...-->Solving quadratic minimization Iter: %d.\n", iterQuad);
            workProcData->active = TRUE;
#else
            fprintf(stderr, "...-->Solving quadratic minimization Iter: %d.\n", iterQuad);
#endif
            curStep++;
            return FALSE;
            /* ------------------------------------------------------------------------ */
        case SOLVE_QUADRATIC:
            if (conRows < conCols) {
                double **aMatrix, *hVector, *gVector;
                int    *pVector, pseudoRank, nLiqs, nSols, nCmps;
                double tolerance = 10.0*DBL_EPSILON /**DBL_EPSILON */, scale = DBL_MIN, fNORM;

                aMatrix =  matrix_alloc(conCols-conRows, conCols-conRows);
                hVector =  vector_alloc(conCols-conRows);
                gVector =  vector_alloc(conCols-conRows);
                pVector =  permutation_alloc(conCols-conRows);

#ifdef PHMELTS_ADJUSTMENTS
                nLiqs = (silminState->incLiquids > 1) ? silminState->nLiquidCoexist : 1;
#else
                nLiqs = (silminState->multipleLiqs) ? silminState->nLiquidCoexist : 1;
#endif
                for (i=0, nCmps=0; i<nc;  i++) if ((silminState->bulkComp)[i] != 0.0) nCmps++;
                for (i=0, nSols=0; i<npc; i++) if (solids[i].type == PHASE) nSols += (silminState->nSolidCoexist)[i];
#ifdef DEBUG
                printf("\nIteration %d. Entering HFTI with %d component(s), %d liquid(s) and %d solid(s)\n", iterQuad, nCmps, nLiqs, nSols);
                if (nCmps < (nLiqs+nSols))
                    printf("*****Assemblage input to HFTI violates phase rule for divariant assemblage.\n");
#endif

                for (i=0; i<(conCols-conRows); i++) {
                    if (fabs(bMatrix[conRows+i][0]) > scale) scale = fabs(bMatrix[conRows+i][0]);
                    for (j=0; j<(conCols-conRows); j++) if (fabs(eMatrix[i+conRows][j+conRows]) > scale) scale = fabs(eMatrix[i+conRows][j+conRows]);
                }
                for (i=0; i<(conCols-conRows); i++) {
                    bMatrix[conRows+i][0] /= scale;
                    for (j=0; j<(conCols-conRows); j++) aMatrix[i][j] = eMatrix[i+conRows][j+conRows]/scale;
                }

                for (i=0, fNORM=0.0; i<(conCols-conRows); i++) for (j=0; j<(conCols-conRows); j++) fNORM += aMatrix[i][j]*aMatrix[i][j];
                fNORM = sqrt(fNORM);
#ifdef DEBUG
                printf("Frobenious norm of HFTI matrix = %g\n", fNORM);
#endif
                hfti(aMatrix, conCols-conRows, conCols-conRows, &bMatrix[conRows], 1, tolerance, &pseudoRank, &rNorm, hVector, gVector, pVector);

                if (pseudoRank < conCols-conRows) {
#ifndef BATCH_VERSION
                    wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...-->Rank deficiency detected by HFTI, rank = %d\n", pseudoRank);
#else
                    fprintf(stderr, "...-->Rank deficiency detected by HFTI, rank = %d\n", pseudoRank);
                    meltsStatus.status = SILMIN_RANK;
#endif
                }

#ifdef DEBUG
                printf("HFTI: scale factor: %20.13g\n", scale);
                /*
         for (i=0; i<(conCols-conRows); i++) printf("HFTI Soln: bMatrix[%d] = %20.13g\n", conRows+i, bMatrix[conRows+i][0]);
         */
#endif

                permutation_free(pVector, conCols-conRows);
                vector_free (gVector, conCols-conRows);
                vector_free (hVector, conCols-conRows);
                matrix_free (aMatrix, conCols-conRows, conCols-conRows);
            } else if (silminState->isenthalpic && (silminState->refEnthalpy != 0.0)) {
                correctTforChangeInEnthalpy();
#ifndef BATCH_VERSION
                updateStatusADB(STATUS_ADB_INDEX_T, &(silminState->T));
                tpValues[TP_PADB_INDEX_T_INITIAL].value = silminState->T - 273.15;
#endif
#ifdef ALPHAMELTS_UPDATE_SYSTEM
                // CHECK WHAT SHOULD HAPPEN FOR MATLAB/Python
                //silminState->dspTstart = silminState->T - 273.15; Fix: MSG 2/11/15
#else
                silminState->dspTstart = silminState->T;
#endif
            } else if (silminState->isentropic  && (silminState->refEntropy  != 0.0)) { correctTforChangeInEntropy();
#ifndef BATCH_VERSION
                updateStatusADB(STATUS_ADB_INDEX_T, &(silminState->T));
                tpValues[TP_PADB_INDEX_T_INITIAL].value = silminState->T - 273.15;
#endif
#ifdef ALPHAMELTS_UPDATE_SYSTEM
                //silminState->dspTstart = silminState->T - 273.15; Fix: MSG 2/11/15
#else
                silminState->dspTstart = silminState->T;
#endif
            } else if (silminState->isochoric   && (silminState->refVolume   != 0.0)) {
                correctPforChangeInVolume();
#ifndef BATCH_VERSION
                updateStatusADB(STATUS_ADB_INDEX_P, &(silminState->P));
                tpValues[TP_PADB_INDEX_P_INITIAL].value = silminState->P;
#endif
#ifndef ALPHAMELTS_UPDATE_SYSTEM
                silminState->dspPstart = silminState->P;
#endif
            }

#ifndef BATCH_VERSION
            workProcData->active = TRUE;
#endif

            curStep++;
            return FALSE;
            /* ------------------------------------------------------------------------ */
        case REASSEMBLE_SOLUTION:
        {

            hasNlCon =  (silminState->isenthalpic && (silminState->refEnthalpy != 0.0))
            || (silminState->isentropic  && (silminState->refEntropy  != 0.0))
            || (silminState->isochoric   && (silminState->refVolume   != 0.0))
#ifdef PHMELTS_ADJUSTMENTS
            || (silminState->H2Obuffer && (silminState->aH2O != 0.0) && hasLiquid)
#endif
            || (silminState->fo2Path != FO2_NONE);
            for (i=0; i<conRows; i++) bMatrix[i][0] = yVector[i];
            for (i=(conRows-1); i>=0; i--) householderRowCol(HOUSEHOLDER_CALC_MODE_H2, i, i+1, conCols-1, cMatrix, i, &hVector[i], bMatrix, 0, 0);
#ifdef DEBUG
            /*
       for (i=0; i<conCols; i++) printf("unProj HFTI Soln: bMatrix[%d] = %20.13g\n", i, bMatrix[i][0]);
       */
#endif

            j=0; rNorm=0.0; sNorm = 0.0;
            if (hasLiquid) {
                for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                    for (i=0; i<nlc; i++) {
                        if((silminState->liquidComp)[nl][i] != 0.0) {
                            (silminState->liquidDelta)[nl][i] = (hasNlCon) ? bMatrix[j][0] : bMatrix[j][0] - (silminState->liquidComp)[nl][i];
                            sNorm += (hasNlCon) ? SQUARE(bMatrix[j][0]+(silminState->liquidComp)[nl][i]) : SQUARE(bMatrix[j][0]);
                            rNorm += SQUARE((silminState->liquidDelta)[nl][i]);
#ifdef DEBUG
                            printf("soln %-15.15s = %13.6g  ref = %13.6g  delta = %13.6g\n", liquid[i].label,
                   (hasNlCon) ? bMatrix[j][0]+(silminState->liquidComp)[nl][i] : bMatrix[j][0], (silminState->liquidComp)[nl][i], (silminState->liquidDelta)[nl][i]);
#endif
                            j++;
                        } else (silminState->liquidDelta)[nl][i] = 0.0;
                    }
                    if (hasNlCon) for (i=0; i<nlc; i++) (constraints->liquidDelta)[nl][i] = 0.0;
                }
            }

#ifdef PRINT_ENERGY_AT_EACH_QUAD_ITERATION
            if (iterQuad == 1) printf("T = %g, P = %g\n", silminState->T - 273.15, silminState->P);
#endif
            for (i=0; i<npc; i++) {
                for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
                    (silminState->solidDelta)[i][ns] = 0.0;
#ifdef DEBUG
                    printf("soln %-15.15s\n", solids[i].label);
#endif
#ifdef PRINT_ENERGY_AT_EACH_QUAD_ITERATION
                    if (iterQuad == 1) printf("%-15.15s ", solids[i].label);
#endif
                    if (solids[i].na == 1) {
                        (silminState->solidDelta)[i][ns] = (hasNlCon) ? bMatrix[j][0] : bMatrix[j][0] - (silminState->solidComp)[i][ns];
                        sNorm  += (hasNlCon) ? SQUARE(bMatrix[j][0]+(silminState->solidComp)[i][ns]) : SQUARE(bMatrix[j][0]);
                        rNorm  += SQUARE((silminState->solidDelta)[i][ns]);
#ifdef DEBUG
                        printf("soln %-15.15s = %13.6g  ref = %13.6g  delta = %13.6g\n", "Total moles",
                         (hasNlCon) ? bMatrix[j][0]+(silminState->solidComp)[i][ns] : bMatrix[j][0], (silminState->solidComp)[i][ns], (silminState->solidDelta)[i][ns]);
#endif
                        j++;
                    } else {
                        for (k=0, mTotal=0.0; k<solids[i].na; k++) {
                            if ((silminState->solidComp)[i+1+k][ns] != 0.0) {
                                (silminState->solidDelta)[i+1+k][ns] = (hasNlCon) ? bMatrix[j][0] : bMatrix[j][0] - (silminState->solidComp)[i+1+k][ns];
                                mTotal += (hasNlCon) ? bMatrix[j][0]+(silminState->solidComp)[i+1+k][ns] : bMatrix[j][0];
                                sNorm  += (hasNlCon) ? SQUARE(bMatrix[j][0] +(silminState->solidComp)[i+1+k][ns]) : SQUARE(bMatrix[j][0]);
                                rNorm  += SQUARE((silminState->solidDelta)[i+1+k][ns]);
#ifdef DEBUG
                                printf("soln %-15.15s = %13.6g  ref = %13.6g  delta = %13.6g\n", solids[i+1+k].label,
                                 (hasNlCon) ? bMatrix[j][0]+(silminState->solidComp)[i+1+k][ns] : bMatrix[j][0], (silminState->solidComp)[i+1+k][ns], (silminState->solidDelta)[i+1+k][ns]);
#endif
                                j++;
                            } else (silminState->solidDelta)[i+1+k][ns] = 0.0;
                        }
                        (silminState->solidDelta)[i][ns] = mTotal - (silminState->solidComp)[i][ns];
#ifdef DEBUG
                        printf("soln %-15.15s = %13.6g  ref = %13.6g  delta = %13.6g\n", "Total moles", mTotal, (silminState->solidComp)[i][ns], (silminState->solidDelta)[i][ns]);
#endif
                    }
                }
            }
#ifdef PRINT_ENERGY_AT_EACH_QUAD_ITERATION
            if (iterQuad == 1) printf("\n");
#endif

            if ((silminState->isenthalpic && (silminState->refEnthalpy != 0.0)) || (silminState->isentropic  && (silminState->refEntropy  != 0.0))) {
                silminState->tDelta = bMatrix[j][0];
                constraints->T = 0.0;
                sNorm += SQUARE((bMatrix[j][0]+silminState->T)/SCALET);
                rNorm += SQUARE(bMatrix[j][0]/SCALET);
#ifdef DEBUG
                printf("soln %-15.15s = %13.6g, delta = %13.6g\n", "T (C)", bMatrix[j][0]+silminState->T-273.15, bMatrix[j][0]);
#endif
                j++;
            } else if (silminState->isochoric && (silminState->refVolume != 0.0)) {
                silminState->pDelta = bMatrix[j][0];
                constraints->P = 0.0;
                sNorm += SQUARE((bMatrix[j][0]+silminState->P)/SCALEP);
                rNorm += SQUARE(bMatrix[j][0]/SCALEP);
#ifdef DEBUG
                printf("soln %-15.15s = %13.6g, delta = %13.6g\n", "P (bars)", bMatrix[j][0]+silminState->P, bMatrix[j][0]);
#endif
                j++;
            }

            rNorm = sqrt(rNorm);  sNorm = sqrt(sNorm);

            /* New code to save best iteration that meets non-optimal criteria */
            if (rNorm < sqrt(DBL_EPSILON)*sNorm && rNorm < bestrNorm) {
                int fractionateSol = silminState->fractionateSol, fractionateFlu = silminState->fractionateFlu,
                    fractionateLiq = silminState->fractionateLiq;
                silminState->fractionateSol = FALSE;
                silminState->fractionateFlu = FALSE;
                silminState->fractionateLiq = FALSE;

                bestrNorm = rNorm; bestIter = iterQuad;
                acceptable = TRUE;
                bestState = copySilminStateStructure(silminState, bestState);

                silminState->fractionateSol = fractionateSol;
                silminState->fractionateFlu = fractionateFlu;
                silminState->fractionateLiq = fractionateLiq;

            }

#ifndef BATCH_VERSION
            wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...-->rNorm = %13.6g, sNorm = %13.6g\n", rNorm, sNorm);
#else
            fprintf(stderr, "...-->rNorm = %13.6g, sNorm = %13.6g\n", rNorm, sNorm);
#endif

#ifdef PRINT_ENERGY_AT_EACH_QUAD_ITERATION
            printf("iter = %3.3d rNorm = %13.6e, sNorm = %13.6e", iterQuad, rNorm, sNorm);
#endif

            if (rNorm < pow(DBL_EPSILON, (double) 0.75)*sNorm*((double) quad_tol_modifier)) curStep = CONVERGENCE_TEST;
            else if (iterQuad > ITERMX) {
                if (rNorm < sqrt(DBL_EPSILON)*sNorm*((double) quad_tol_modifier)) {
                    curStep = CONVERGENCE_TEST;
#ifndef BATCH_VERSION
                    wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...Quadratic convergence accepted, but non-optimal.\n");
#elif defined(ALPHAMELTS_UPDATE_SYSTEM)
                    printf("...Quadratic convergence accepted, but non-optimal.\n");
#else
                    fprintf(stderr, "...Quadratic convergence accepted, but non-optimal.\n");
#endif
                } else if (acceptable) {
                    int fractionateSol = silminState->fractionateSol, fractionateFlu = silminState->fractionateFlu,
                        fractionateLiq = silminState->fractionateLiq;
                    silminState->fractionateSol = FALSE;
                    silminState->fractionateFlu = FALSE;
                    silminState->fractionateLiq = FALSE;

                    curStep = CONVERGENCE_TEST;
#ifndef BATCH_VERSION
                    wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...Quadratic convergence was acceptable at iterQuad = %d (rNorm = %g).\n", bestIter, bestrNorm);
#elif defined(ALPHAMELTS_UPDATE_SYSTEM)
                    printf("...Quadratic convergence was acceptable at iterQuad = %d (rNorm = %g).\n", bestIter, bestrNorm);
#else
                    fprintf(stderr, "...Quadratic convergence was acceptable at iterQuad = %d (rNorm = %g).\n", bestIter, bestrNorm);
#endif
                    silminState = copySilminStateStructure(bestState, silminState);

                    silminState->fractionateSol = fractionateSol;
                    silminState->fractionateFlu = fractionateFlu;
                    silminState->fractionateLiq = fractionateLiq;

                } else {
#ifndef BATCH_VERSION
                    XmString csString1, csString2, csString3;

                    csString1 = XmStringCreateLtoR("THE QUADRATIC MINIMIZATION ALGORITHM HAS FAILED TO CONVERGE.\n", "ISO8859-1");
                    csString2 = XmStringCreateLtoR("Try restarting the calculation at the current point with\n", "ISO8859-1");
                    csString3 = XmStringConcat(csString1, csString2);
                    XmStringFree(csString1); XmStringFree(csString2);

                    csString1 = XmStringCreateLtoR("slightly modified system parameters. You may wish to save the\n", "ISO8859-1");
                    csString2 = XmStringConcat(csString3, csString1);
                    XmStringFree(csString3); XmStringFree(csString1);

                    csString1 = XmStringCreateLtoR("current system state before proceeding.", "ISO8859-1");
                    csString3 = XmStringConcat(csString2, csString1);
                    XmStringFree(csString1); XmStringFree(csString2);

                    XtVaSetValues(error_message, XmNmessageString, csString3, NULL);
                    XtManageChild(error_message);
                    XmStringFree(csString3);

                    wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...Quadratic convergence failure. Aborting.\n");
                    workProcData->active = FALSE;
#elif defined(ALPHAMELTS_UPDATE_SYSTEM)
                    printf("...Quadratic convergence failure. Aborting.\n");
                    meltsStatus.status = SILMIN_QUAD_MAX;
#else
                    fprintf(stderr, "...Quadratic convergence failure. Aborting.\n");
                    meltsStatus.status = SILMIN_QUAD_MAX;
#endif
                    iterQuad = 0; curStep = 0; return TRUE;
                }
            } else curStep++;
        }
        jumpFromLinSearch:
            if ((curStage != PRE_STAGE_ZERO) && (curStep == CONVERGENCE_TEST)) {
                static int fo2PathOld = FO2_NONE;
#ifdef PHMELTS_ADJUSTMENTS
                if (H2Obuffer && (silminState->aH2O != 0.0) && hasLiquid) {
                    if      (curStage == L_H_STAGE_ONE)   { silminState->H2Obuffer = FALSE; silminState->isenthalpic = TRUE;  curStage = L_H_STAGE_TWO;   }
                    else if (curStage == L_H_STAGE_TWO)   { silminState->H2Obuffer = TRUE;  silminState->isenthalpic = FALSE; curStage = L_H_STAGE_THREE; }
                    else if (curStage == L_H_STAGE_THREE) { silminState->H2Obuffer = TRUE;  silminState->isenthalpic = TRUE;  curStage = PRE_STAGE_ZERO ; }
                    else if (curStage == L_H_STAGE_ONE)   { silminState->H2Obuffer = FALSE; silminState->isentropic  = TRUE;  curStage = L_H_STAGE_TWO;   }
                    else if (curStage == L_H_STAGE_TWO)   { silminState->H2Obuffer = TRUE;  silminState->isentropic  = FALSE; curStage = L_H_STAGE_THREE; }
                    else if (curStage == L_H_STAGE_THREE) { silminState->H2Obuffer = TRUE;  silminState->isentropic  = TRUE;  curStage = PRE_STAGE_ZERO ; }
                    else if (curStage == L_H_STAGE_ONE)   { silminState->H2Obuffer = FALSE; silminState->isochoric   = TRUE;  curStage = L_H_STAGE_TWO;   }
                    else if (curStage == L_H_STAGE_TWO)   { silminState->H2Obuffer = TRUE;  silminState->isochoric   = FALSE; curStage = L_H_STAGE_THREE; }
                    else if (curStage == L_H_STAGE_THREE) { silminState->H2Obuffer = TRUE;  silminState->isochoric   = TRUE;  curStage = PRE_STAGE_ZERO ; }

                    /* CHECK THE MUH2O BUFFER */

                } else if (silminState->fo2Iter && !silminState->fo2Alt && hasLiquid) {
                    if      (curStage == L_X_STAGE_ONE  ) { fo2PathOld = silminState->fo2Path;                                   curStage = L_X_STAGE_TWO;   }
                    else if (curStage == L_X_STAGE_TWO  ) { silminState->fo2Path = FO2_NONE;                                     curStage = L_X_STAGE_THREE; }
                    else if (curStage == L_X_STAGE_THREE) { silminState->fo2Path = fo2PathOld; 				       curStage = PRE_STAGE_ZERO ; }
                    else if (curStage == L_H_STAGE_ONE  ) { fo2PathOld = silminState->fo2Path; silminState->isenthalpic = FALSE; curStage = L_H_STAGE_TWO;   }
                    else if (curStage == L_H_STAGE_TWO  ) { silminState->fo2Path = FO2_NONE;   silminState->isenthalpic = TRUE;  curStage = L_H_STAGE_THREE; }
                    else if (curStage == L_H_STAGE_THREE) { silminState->fo2Path = fo2PathOld; silminState->isenthalpic = TRUE;  curStage = PRE_STAGE_ZERO ; }
                    else if (curStage == L_S_STAGE_ONE  ) { fo2PathOld = silminState->fo2Path; silminState->isentropic  = FALSE; curStage = L_S_STAGE_TWO;   }
                    else if (curStage == L_S_STAGE_TWO  ) { silminState->fo2Path = FO2_NONE;   silminState->isentropic  = TRUE;  curStage = L_S_STAGE_THREE; }
                    else if (curStage == L_S_STAGE_THREE) { silminState->fo2Path = fo2PathOld; silminState->isentropic  = TRUE;  curStage = PRE_STAGE_ZERO ; }
                    else if (curStage == L_V_STAGE_ONE  ) { fo2PathOld = silminState->fo2Path; silminState->isochoric   = FALSE; curStage = L_V_STAGE_TWO;   }
                    else if (curStage == L_V_STAGE_TWO  ) { silminState->fo2Path = FO2_NONE;   silminState->isochoric   = TRUE;  curStage = L_V_STAGE_THREE; }
                    else if (curStage == L_V_STAGE_THREE) { silminState->fo2Path = fo2PathOld; silminState->isochoric   = TRUE;  curStage = PRE_STAGE_ZERO ; }
                } else
#endif
                if      (curStage == L_H_STAGE_ONE  ) { fo2PathOld = silminState->fo2Path; silminState->fo2Path = FO2_NONE; silminState->isenthalpic = TRUE;  curStage = L_H_STAGE_TWO;   }
                else if (curStage == L_H_STAGE_TWO  ) { silminState->fo2Path = fo2PathOld;                                  silminState->isenthalpic = FALSE; curStage = L_H_STAGE_THREE; }
                else if (curStage == L_H_STAGE_THREE) { silminState->fo2Path = fo2PathOld;                                  silminState->isenthalpic = TRUE;  curStage = PRE_STAGE_ZERO ; }
                else if (curStage == L_S_STAGE_ONE  ) { fo2PathOld = silminState->fo2Path; silminState->fo2Path = FO2_NONE; silminState->isentropic  = TRUE;  curStage = L_S_STAGE_TWO;	}
                else if (curStage == L_S_STAGE_TWO  ) { silminState->fo2Path = fo2PathOld;				  silminState->isentropic  = FALSE; curStage = L_S_STAGE_THREE; }
                else if (curStage == L_S_STAGE_THREE) { silminState->fo2Path = fo2PathOld;				  silminState->isentropic  = TRUE;  curStage = PRE_STAGE_ZERO ; }
                else if (curStage == L_V_STAGE_ONE  ) { fo2PathOld = silminState->fo2Path; silminState->fo2Path = FO2_NONE; silminState->isochoric   = TRUE;  curStage = L_V_STAGE_TWO;	}
                else if (curStage == L_V_STAGE_TWO  ) { silminState->fo2Path = fo2PathOld;				  silminState->isochoric   = FALSE; curStage = L_V_STAGE_THREE; }
                else if (curStage == L_V_STAGE_THREE) { silminState->fo2Path = fo2PathOld;				  silminState->isochoric   = TRUE;  curStage = PRE_STAGE_ZERO ; }

                if (curStage != PRE_STAGE_ZERO) { curStep = PROJECT_CONSTRAINTS; iterQuad = 0; }
#ifdef PHMELTS_ADJUSTMENTS
                if (!silminState->fo2Iter || silminState->fo2Alt || !hasLiquid) {
#endif
                    if ( (curStage == L_H_STAGE_THREE) || (curStage == L_S_STAGE_THREE) || (curStage == L_V_STAGE_THREE) ) {
                        silminState->fo2 = getlog10fo2(silminState->T, silminState->P, silminState->fo2Path);
                        gibbs(silminState->T, silminState->P, "o2", &(oxygen.ref), NULL, NULL, &(oxygen.cur));
                    }
#ifdef PHMELTS_ADJUSTMENTS
                }
#endif
            }

#ifndef BATCH_VERSION
            workProcData->active = TRUE;
#endif
            return FALSE;
            /* ------------------------------------------------------------------------ */
        case LINEAR_SEARCH:
        {
            static double pTotalLast = 0.0, pTotalAveLast = 0.0, pTotalHistory[ITERMX+1];
            double lambda = 1.0, stepSize = 0.10, pTotal = 0.0, pTotalAverage = 0.0;
            int iter, status;

            do {
                double reltest = MAX(10.0*DBL_EPSILON, DBL_EPSILON*sNorm);
                iter   = ITERMX;
                status = min1d(&lambda, &stepSize, reltest, &iter, &pTotal, linearSearch);
                if        (status == MIN1D_BAD_INITIAL) {
                    lambda /= 2.0; /* was 10.0 */
                    if (lambda < DBL_EPSILON) {
                        status = MIN1D_SUCCESS;
#ifdef BATCH_VERSION
                        meltsStatus.status = SILMIN_LIN_ZERO;
#endif
                    }
                } else if (status == MIN1D_ITERS_EXCEEDED) {
#ifndef BATCH_VERSION
                    wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...-->Linear search: Iteration exceeded - aborting.\n");
#else
                    fprintf(stderr, "...-->Linear search: Iteration exceeded - aborting.\n");
                    meltsStatus.status = SILMIN_LIN_MAX;
#endif
                    status = MIN1D_SUCCESS;
                }
            } while (status != MIN1D_SUCCESS);
#ifdef DEBUG
            printf("On exit from Linear Search routine, iter = %d, lambda = %g", iter, lambda);
            printf(" pTotal = %20.13e (%20.13e, %21.1f)\n", pTotal, (pTotal-pTotalLast)/pTotal, (pTotal-pTotalLast)/(pTotal*DBL_EPSILON));
#endif
#ifdef PRINT_ENERGY_AT_EACH_QUAD_ITERATION
            printf(" lambda = %13.6e pTotal = %20.13e (%20.13e, %21.1f)\n", lambda, pTotal, (pTotal-pTotalLast)/pTotal, (pTotal-pTotalLast)/(pTotal*DBL_EPSILON));
#endif

#ifndef BATCH_VERSION
            wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...-->Linear search: Min = %13.6g, step = %13.6g\n", pTotal, lambda);
            updateStatusADB(STATUS_ADB_INDEX_LINEAR, &iter);
#else
            fprintf(stderr, "...-->Linear search: Min = %13.6g, step = %13.6g\n", pTotal, lambda);
#endif

            if (hasLiquid) {
                silminState->liquidMass = 0.0;
                for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                    for (i=0; i<nlc; i++) {
                        (silminState->liquidComp)[nl][i] += lambda*(silminState->liquidDelta)[nl][i];
#ifdef PHMELTS_ADJUSTMENTS
                        if ((silminState->fo2Path != FO2_NONE) || (silminState->H2Obuffer && (silminState->aH2O != 0.0))) (silminState->liquidComp)[nl][i] += (constraints->liquidDelta)[nl][i];
#else
                        if (silminState->fo2Path != FO2_NONE) (silminState->liquidComp)[nl][i] += (constraints->liquidDelta)[nl][i];
#endif
                        for (j=0; j<nc; j++) silminState->liquidMass += silminState->liquidComp[nl][i]*(liquid[i].liqToOx)[j]*bulkSystem[j].mw;
                    }
                }
            }
            for (i=0; i<npc; i++)
                for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
                    (silminState->solidComp)[i][ns] += lambda*(silminState->solidDelta)[i][ns];
                    if (silminState->fo2Path != FO2_NONE && !hasLiquid) (silminState->solidComp)[i][ns] += (constraints->solidDelta)[i][ns];
                    if (solids[i].na > 1) for (j=0; j<solids[i].na; j++) {
                        (silminState->solidComp)[i+1+j][ns] += lambda*(silminState->solidDelta)[i+1+j][ns];
                        if (silminState->fo2Path != FO2_NONE && !hasLiquid) (silminState->solidComp)[i+1+j][ns] += (constraints->solidDelta)[i+1+j][ns];
                    }
                }
#ifndef PHMELTS_ADJUSTMENTS
            if (silminState->fo2Path != FO2_NONE) {
                for (i=0; i<nc; i++) {
                    (silminState->bulkComp)[i] = 0.0;
                    for (nl=0; nl<silminState->nLiquidCoexist; nl++) for (j=0; j<nlc; j++) (silminState->bulkComp)[i] += (silminState->liquidComp)[nl][j]*(liquid[j].liqToOx)[i];
                    for (j=0; j<npc; j++) for (ns=0; ns<(silminState->nSolidCoexist)[j]; ns++) {
                        if (solids[j].na == 1) (silminState->bulkComp)[i] += (silminState->solidComp)[j][ns]*(solids[j].solToOx)[i];
                        else for (k=0; k<solids[j].na; k++) (silminState->bulkComp)[i] += (silminState->solidComp)[j+1+k][ns]*(solids[j+1+k].solToOx)[i];
                    }
                }
            }
#else
            if ((silminState->fo2Path != FO2_NONE) || (silminState->H2Obuffer && (silminState->aH2O != 0.0) && hasLiquid)) {
                for (i=0; i<nc; i++) {
                    (silminState->bulkComp)[i] = 0.0;
                    for (nl=0; nl<silminState->nLiquidCoexist; nl++) for (j=0; j<nlc; j++) (silminState->bulkComp)[i] += (silminState->liquidComp)[nl][j]*(liquid[j].liqToOx)[i];
                    for (j=0; j<npc; j++) for (ns=0; ns<(silminState->nSolidCoexist)[j]; ns++) {
                        if (solids[j].na == 1) (silminState->bulkComp)[i] += (silminState->solidComp)[j][ns]*(solids[j].solToOx)[i];
                        else for (k=0; k<solids[j].na; k++) (silminState->bulkComp)[i] += (silminState->solidComp)[j+1+k][ns]*(solids[j+1+k].solToOx)[i];
                    }
                }
            }

#ifdef DEBUG
            { /* check totals */
                double sum;
                for (j=0, sum=0.0; j<npc; j++) {
                    for (ns=0; ns<(silminState->nSolidCoexist)[j]; ns++) {
                    if (solids[j].na == 1) {
                        sum += (silminState->solidComp)[j][ns]*solids[j].mw;

                    } else {
                    for (k=0; k<solids[j].na; k++) {
                        sum += (silminState->solidComp)[j+1+k][ns]*solids[j+1+k].mw;
                    }
                    }
                    }
                }
                silminState->solidMass = sum;
            }
#endif
#endif

            /* -> Change system temperature                                         */
            if ((silminState->isenthalpic && (silminState->refEnthalpy != 0.0)) || (silminState->isentropic  && (silminState->refEntropy  != 0.0))) {
#ifdef DEBUG
                printf("Temperature reset from %g -> %g\n", silminState->T, constraints->T);
#endif
                silminState->T = constraints->T;
#ifndef BATCH_VERSION
                updateStatusADB(STATUS_ADB_INDEX_T, &(silminState->T));
                tpValues[TP_PADB_INDEX_T_INITIAL].value = silminState->T - 273.15;
#endif
#ifdef ALPHAMELTS_UPDATE_SYSTEM
                //silminState->dspTstart = silminState->T - 273.15; Fix: MSG 2/11/15
#else
                silminState->dspTstart = silminState->T;
#endif

                if (silminState->fo2Path != FO2_NONE) {
#ifdef DEBUG
                    printf("log10 fO2 reset from %g -> %g\n", silminState->fo2, constraints->fo2);
#endif
                    silminState->fo2 = constraints->fo2;
#ifndef BATCH_VERSION
                    updateStatusADB(STATUS_ADB_INDEX_LOGFO2, &(silminState->fo2));
#endif
                }
            }

            /* -> Change system pressure                                            */
            if (silminState->isochoric   && (silminState->refVolume   != 0.0)) {
#ifdef DEBUG
                printf("Pressure reset from %g -> %g\n", silminState->P, constraints->P);
#endif
                silminState->P = constraints->P;
#ifndef BATCH_VERSION
                updateStatusADB(STATUS_ADB_INDEX_P, &(silminState->P));
                tpValues[TP_PADB_INDEX_P_INITIAL].value = silminState->P;
#endif
#ifndef ALPHAMELTS_UPDATE_SYSTEM
                silminState->dspPstart = silminState->P;
#endif

                if (silminState->fo2Path != FO2_NONE) {
#ifdef DEBUG
                    printf("log10 fO2 reset from %g -> %g\n", silminState->fo2, constraints->fo2);
#endif
                    silminState->fo2 = constraints->fo2;
#ifndef BATCH_VERSION
                    updateStatusADB(STATUS_ADB_INDEX_LOGFO2, &(silminState->fo2));
#endif
                }
            }

            /* -> Calculate liquid and solid end-member properties                  */
            if ((silminState->isenthalpic && (silminState->refEnthalpy != 0.0))
                || (silminState->isentropic  && (silminState->refEntropy  != 0.0))
                || (silminState->isochoric   && (silminState->refVolume   != 0.0))) {
                for (i=0; i<nlc; i++) gibbs(silminState->T, silminState->P, (char *) liquid[i].label, &(liquid[i].ref), &(liquid[i].liq), &(liquid[i].fus), &(liquid[i].cur));
                for (i=0, j=0; i<npc; i++) {
                    if (solids[i].type == PHASE) {
                        if ((silminState->incSolids)[j]) {
                            if(solids[i].na == 1) gibbs(silminState->T, silminState->P, (char *) solids[i].label,  &(solids[i].ref), NULL, NULL, &(solids[i].cur));
                            else {
                                for (k=0; k<solids[i].na; k++) gibbs(silminState->T, silminState->P, (char *) solids[i+1+k].label, &(solids[i+1+k].ref), NULL, NULL, &(solids[i+1+k].cur));
                                i += solids[i].na;
                            }
                        }
                        j++;
                    }
                }
                if (silminState->fo2Path != FO2_NONE) gibbs(silminState->T, silminState->P, "o2", &(oxygen.ref), NULL, NULL, &(oxygen.cur));
            }

#ifdef DEBUG
            if ((silminState->fo2Path != FO2_NONE) && hasLiquid) {
                double muO2E =  (oxygen.cur).g + R*silminState->T*log((double) 10.0)*getlog10fo2(silminState->T, silminState->P, silminState->fo2Path);
                for (i=0; i<silminState->nLiquidCoexist; i++) {
                    double muO2L;
                    muO2Liq(FIRST, silminState->T, silminState->P, (silminState->liquidComp)[i], &muO2L, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
                    muO2L += (oxygen.cur).g;
                    printf("muO2L = %g, muO2E = %g, diff = %g\n", muO2L, muO2E, muO2L-muO2E);
                }
            }
#endif

            pTotalHistory[iterQuad] = pTotal;
            if (iterQuad > 5) {
                pTotalAverage = (pTotalHistory[iterQuad-5]+pTotalHistory[iterQuad-4]+pTotalHistory[iterQuad-3]+pTotalHistory[iterQuad-2]+pTotalHistory[iterQuad-1])/5.0;
                if (fabs((pTotalAverage-pTotalAveLast)/pTotalAverage) < 10.0*DBL_EPSILON) {
                    curStep = CONVERGENCE_TEST;
#ifndef BATCH_VERSION
                    wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...-->Linear search: Convergence detected.\n");
#else
                    fprintf(stderr, "...-->Linear search: Convergence detected.\n");
#endif
#ifdef DEBUG
                    printf("<><><> Optimal energy exit. <><><>\n");
#endif
#ifdef PRINT_ENERGY_AT_EACH_QUAD_ITERATION
                    printf("<><><> Optimal energy exit. <><><>\n");
#endif
                    goto jumpFromLinSearch;
                }
            } else pTotalAverage = 0.0;

            pTotalLast = pTotal;
            pTotalAveLast = pTotalAverage;

        }

#ifndef BATCH_VERSION
            workProcData->active = TRUE;
#endif

            curStep++;
            return FALSE;
            /* ------------------------------------------------------------------------ */
        case DROP_PHASE:
    {
#ifdef PHMELTS_ADJUSTMENTS
        //           double outmass = MASSOUT*silminState->refMass/100.0;
     double outmass = MASSOUT;
#else
     double outmass = MASSOUT;
#endif

            curStep = ((silminState->isenthalpic && (silminState->refEnthalpy != 0.0))
             || (silminState->isentropic  && (silminState->refEntropy  != 0.0))
             || (silminState->isochoric   && (silminState->refVolume   != 0.0))
#ifdef PHMELTS_ADJUSTMENTS
         || (silminState->H2Obuffer && (silminState->aH2O != 0.0) && hasLiquid)
#endif
             || (silminState->fo2Path != FO2_NONE)
             ) ? PROJECT_CONSTRAINTS : CONSTRUCT_QUADRATIC;

            /* check if a solid phase is trying to disappear - add the solid components to the 1st liquid */
            for (i=0; i<npc; i++) {
                int nDrop = 0;

                for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
                    if ( (silminState->solidComp)[i][ns] < outmass) {
#ifndef BATCH_VERSION
                        wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...Dropping phase %s from the assemblage.\n", solids[i].label);
#elif defined(ALPHAMELTS_UPDATE_SYSTEM)
                        printf("...Dropping phase %s from the assemblage.\n", solids[i].label);
#else
                        fprintf(stderr, "...Dropping phase %s from the assemblage.\n", solids[i].label);
#endif
                        curStep = PROJECT_CONSTRAINTS;
                        if (solids[i].na == 1) {
                            if (hasLiquid) {
                                for (j=0; j<nlc; j++) (silminState->liquidComp)[0][j] += (solids[i].solToLiq)[j] * (silminState->solidComp)[i][ns];
                                (silminState->nSolidCoexist)[i]--;
                            } else {
                                double *deltaBulkComp = (double *) calloc((size_t) nc,sizeof(double));
                                for (j=0; j<nc; j++) deltaBulkComp[j] -= (solids[i].solToOx)[j] * (silminState->solidComp)[i][ns];
#ifdef PHMELTS_ADJUSTMENTS
                                if (fabs(deltaBulkComp[iBulkH2O] + silminState->bulkComp[iBulkH2O]) < DBL_EPSILON) {

                                    /* set water in NAM to have same specific entropy as vapour */
                                    if (silminState->isentropic && (silminState->refEntropy != 0.0))
                                    silminState->refEntropy  -= silminState->bulkComp[iBulkH2O]*(solids[iPhaseH2O].cur).s;
                                    else if(silminState->isenthalpic && (silminState->refEnthalpy != 0.0))
                                    silminState->refEnthalpy -= silminState->bulkComp[iBulkH2O]*(solids[iPhaseH2O].cur).h;
                                    else if(silminState->isochoric && (silminState->refVolume != 0.0))
                                    silminState->refEnthalpy -= silminState->bulkComp[iBulkH2O]*(solids[iPhaseH2O].cur).v;

                                    if(silminState->fo2Path != FO2_NONE && silminState->oxygen != 0.0)
                                    silminState->oxygen -= silminState->bulkComp[iBulkH2O]*(oxygen.liqToOx)[iLiqH2O];

                                    deltaBulkComp[iBulkH2O] = 0.0;
                                    silminState->bulkComp[iBulkH2O] = 0.0;
                                }
#endif
                                (silminState->nSolidCoexist)[i]--;
                                addOrDropLiquid(deltaBulkComp);
                                free(deltaBulkComp);
                            }
                            (silminState->solidComp)[i][ns] = 0.0;
                        } else {
                            double *deltaBulkComp = (double *) calloc((size_t) nc,sizeof(double));
                            (silminState->solidComp)[i][ns] = 0.0;
                            for (j=0; j<solids[i].na; j++) {
                                if (hasLiquid) for (k=0; k<nlc; k++) (silminState->liquidComp)[0][k] += (solids[i+1+j].solToLiq)[k]*(silminState->solidComp)[i+1+j][ns];
                                else for (k=0; k<nc; k++) deltaBulkComp[k] -= (solids[i+1+j].solToOx)[k]*(silminState->solidComp)[i+1+j][ns];
                                (silminState->solidComp)[i+1+j][ns] = 0.0;
                            }
#ifdef PHMELTS_ADJUSTMENTS
                if (fabs(deltaBulkComp[iBulkH2O] + silminState->bulkComp[iBulkH2O]) < DBL_EPSILON) {

                    /* set water in NAM to have same specific entropy as vapour */
                    if (silminState->isentropic && (silminState->refEntropy != 0.0))
                    silminState->refEntropy -= silminState->bulkComp[iBulkH2O]*(solids[iPhaseH2O].cur).s;
                    else if(silminState->isenthalpic && (silminState->refEnthalpy != 0.0))
                    silminState->refEnthalpy -= silminState->bulkComp[iBulkH2O]*(solids[iPhaseH2O].cur).h;
                    else if(silminState->isochoric && (silminState->refVolume != 0.0))
                    silminState->refEnthalpy -= silminState->bulkComp[iBulkH2O]*(solids[iPhaseH2O].cur).v;

                    if(silminState->fo2Path != FO2_NONE && silminState->oxygen != 0.0)
                    silminState->oxygen -= silminState->bulkComp[iBulkH2O]*(oxygen.liqToOx)[iLiqH2O];

                    deltaBulkComp[iBulkH2O] = 0.0;
                    silminState->bulkComp[iBulkH2O] = 0.0;
                }
#endif
                            silminState->nSolidCoexist[i]--;
                            if (!hasLiquid) addOrDropLiquid(deltaBulkComp);
                            silminState->nSolidCoexist[i]++;
                            free(deltaBulkComp);
                            nDrop++;
                        }
                        if (silminState->fo2Path != FO2_NONE && hasLiquid) {
                            double *moles = (double *) malloc((size_t) nc*sizeof(double));
                            /* Just the first liquid, because only its composition was changed */
                            for (j=0; j<nc; j++) {
                                for (k=0, moles[j]=0.0; k<nlc; k++) moles[j] += (silminState->liquidComp)[0][k]*(liquid[k].liqToOx)[j];
                                (silminState->bulkComp)[j] -= moles[j];
                            }
                            conLiq(FIRST | SEVENTH, FIRST, silminState->T, silminState->P, moles, NULL, NULL, NULL, NULL, NULL, &(silminState->fo2));
                            for (j=0; j<nc; j++) (silminState->bulkComp)[j] += moles[j];
                            for (j=0; j<nlc; j++) for (k=0, (silminState->liquidComp)[0][j]=0.0; k<nc; k++) (silminState->liquidComp)[0][j] += moles[k]*(bulkSystem[k].oxToLiq)[j];
                            free(moles);
                        } else if (silminState->fo2Path != FO2_NONE && !hasLiquid) {
                            double muO2;
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
                        (silminState->cylSolids)[i]++;
                        iterQuad = 0; bestrNorm = 1.0; acceptable = FALSE;
                    }
                }

                if (nDrop > 0) {
                    for (ns=0, k=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
                        if ((silminState->solidComp)[i][ns] != 0.0) {
                            (silminState->solidComp)[i][k] = (silminState->solidComp)[i][ns];
                            if (solids[i].na > 1) for (j=0; j<solids[i].na; j++) (silminState->solidComp)[i+1+j][k] = (silminState->solidComp)[i+1+j][ns];
                            k++;
#ifdef PHMELTS_ADJUSTMENTS
                            /* allow selective fractionation of one instance of a phase (e.g. cpx)? */
                            /* move fractionated phases if needed, as well as current composition */
                            /* CODE IN HERE */
#endif
                        }
                    }
                    (silminState->nSolidCoexist)[i] -= nDrop;
                }

            }

            /* check if a liquid phase is trying to disappear - add the liquid components to the 1st liquid */
#ifdef PHMELTS_ADJUSTMENTS
            if (hasLiquid && (silminState->incLiquids > 1) && (silminState->nLiquidCoexist > 1)) {
#else
            if (hasLiquid && silminState->multipleLiqs && (silminState->nLiquidCoexist > 1)) {
#endif
                int nDrop = 0;

                for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                    double mass = 0.0;
                    for (i=0; i<nlc; i++) for (j=0; j<nc; j++) mass += silminState->liquidComp[nl][i]*(liquid[i].liqToOx)[j]*bulkSystem[j].mw;

                    if (mass < outmass) {
#ifndef BATCH_VERSION
                        wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...Dropping liquid %d from the assemblage.\n", nl);
#elif defined(ALPHAMELTS_UPDATE_SYSTEM)
                        printf("...Dropping multiple liquid %d (of %d) from the assemblage.\n", nl, silminState->nLiquidCoexist);
#else
                        fprintf(stderr, "...Dropping multiple liquid %d (of %d) from the assemblage.\n", nl, silminState->nLiquidCoexist);
#endif

//#ifdef DEBUG
//                        printf("...Dropping multiple liquid %d (of %d) from the assemblage.\n", nl, silminState->nLiquidCoexist);
//#endif

                        curStep = PROJECT_CONSTRAINTS;
                        nDrop++;

                        if (nl > 0) {
#ifdef PHMELTS_ADJUSTMENTS
                            double *deltaBulkComp = (double *) calloc((size_t) nc,sizeof(double));

                            /* for multiple liquid case, continue to deposit mass from dropped liquid in solid assemblage
                                rather than in other liquids, since the latter operation may cause cycling across the solvus */
                            for (i=0; i<nlc; i++) {
                                for (k=0; k<nc; k++) deltaBulkComp[k] -= (liquid[i].liqToOx)[k]*(silminState->liquidComp)[nl][i];
                                silminState->liquidComp[nl][i] = 0.0;
                            }
                            addOrDropLiquid(deltaBulkComp);
                            free(deltaBulkComp);

                            silminState->cylLiquids++;
#else
                            for (i=0; i<nlc; i++) {
                                (silminState->liquidComp)[0][i] += (silminState->liquidComp)[nl][i];
                                (silminState->liquidComp)[nl][i] = 0.0;
                            }
#endif
                        }
                        for (i=nl+1; i<silminState->nLiquidCoexist; i++) for (j=0; j<nlc; j++) (silminState->liquidComp)[i-1][j] += (silminState->liquidComp)[i][j];
                        free((silminState->liquidComp )[silminState->nLiquidCoexist - 1]); (silminState->liquidComp )[silminState->nLiquidCoexist - 1] = NULL;
                        free((silminState->liquidDelta)[silminState->nLiquidCoexist - 1]); (silminState->liquidDelta)[silminState->nLiquidCoexist - 1] = NULL;
                        silminState->nLiquidCoexist--;
                        iterQuad = 0; bestrNorm = 1.0; acceptable = FALSE;
                    }
                }

                if ((nDrop > 0) && (silminState->fo2Path != FO2_NONE)) {
                    double *moles = (double *) malloc((size_t) nc*sizeof(double));
                    /* Just the first liquid, because only its composition was changed */
                    for (j=0; j<nc; j++) {
                        for (k=0, moles[j]=0.0; k<nlc; k++) moles[j] += (silminState->liquidComp)[0][k]*(liquid[k].liqToOx)[j];
                        (silminState->bulkComp)[j] -= moles[j];
                    }
                    conLiq(FIRST | SEVENTH, FIRST, silminState->T, silminState->P, moles, NULL, NULL, NULL, NULL, NULL, &(silminState->fo2));
                    for (j=0; j<nc; j++) (silminState->bulkComp)[j] += moles[j];
                    for (j=0; j<nlc; j++) for (k=0, (silminState->liquidComp)[0][j]=0.0; k<nc; k++) (silminState->liquidComp)[0][j] += moles[k]*(bulkSystem[k].oxToLiq)[j];
                    free(moles);
                }
            }

            /* check if liquid is trying to disappear - if the liquid mass is this low, there will only be one liquid */
            if (hasLiquid) {
                if (silminState->liquidMass < outmass) {
                    int success;
                    double *deltaBulkComp = (double *) calloc((size_t) nc,sizeof(double));
#ifndef BATCH_VERSION
                    wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...Dropping liquid from the assemblage.\n");
#elif defined(ALPHAMELTS_UPDATE_SYSTEM)
                    printf("...Dropping liquid from the assemblage.\n");
#else
                    fprintf(stderr, "...Dropping liquid from the assemblage.\n");
#endif

//#ifdef DEBUG
//                    printf("...Dropping liquid from the assemblage.\n");
//#endif

                    curStep = PROJECT_CONSTRAINTS;
                    silminState->liquidMass = 0.0; hasLiquid = FALSE;
#ifdef PHMELTS_ADJUSTMENTS
                    /* important to eliminate water from system when using aH2Obuffer */
                    if (fabs(silminState->bulkComp[iBulkH2O] - silminState->liquidComp[0][iLiqH2O]) < DBL_EPSILON) {

                        /* set water in NAM to have same specific entropy as vapour */
                        if (silminState->isentropic && (silminState->refEntropy != 0.0))
                        silminState->refEntropy  -= silminState->bulkComp[iBulkH2O]*(solids[iPhaseH2O].cur).s;
                        else if(silminState->isenthalpic && (silminState->refEnthalpy != 0.0))
                        silminState->refEnthalpy -= silminState->bulkComp[iBulkH2O]*(solids[iPhaseH2O].cur).h;
                        else if(silminState->isochoric && (silminState->refVolume != 0.0))
                        silminState->refEnthalpy -= silminState->bulkComp[iBulkH2O]*(solids[iPhaseH2O].cur).v;

                        if(silminState->fo2Path != FO2_NONE && silminState->oxygen != 0.0)
                        silminState->oxygen -= silminState->bulkComp[iBulkH2O]*(oxygen.liqToOx)[iLiqH2O];

                        deltaBulkComp[iBulkH2O] = 0.0;
                        silminState->bulkComp[iBulkH2O] = 0.0;
                    }
#endif
                    for (i=0; i<nlc; i++) {
                        for (k=0; k<nc; k++) deltaBulkComp[k] -= (liquid[i].liqToOx)[k]*(silminState->liquidComp)[0][i];
                        silminState->liquidComp[0][i] = 0.0;
                    }
                    success = addOrDropLiquid(deltaBulkComp);
                    iterQuad = 0; bestrNorm = 1.0; acceptable = FALSE;
                    free(deltaBulkComp);

#ifdef PHMELTS_ADJUSTMENTS
#ifdef DEBUG
                    { /* check totals */
                        double sum;
                        for (j=0, sum=0.0; j<npc; j++) {
                            for (ns=0; ns<(silminState->nSolidCoexist)[j]; ns++) {
                            if (solids[j].na == 1) {
                                sum += (silminState->solidComp)[j][ns]*solids[j].mw;

                            } else {
                            for (k=0; k<solids[j].na; k++) {
                                sum += (silminState->solidComp)[j+1+k][ns]*solids[j+1+k].mw;
                            }
                            }
                            }
                        }
                        silminState->solidMass = sum;
                    }
#endif

                    //silminState->nLiquidCoexist = 0;
                    if (silminState->fo2Liq != FO2_NONE) {
                        silminState->fo2Liq = silminState->fo2Path;
                        silminState->fo2Path = FO2_NONE;
                    }

                    /* reset aH2O before turning on buffer */
                    /*	  if ((getenv("ALPHAMELTS_DO_TRACE_H2O") != NULL) && (silminState->bulkComp[iBulkH2O] == 0.0)) {
                        (void) evaluateSaturationState(ySol, yLiq);
                        if(silminState->aH2O != 0.0) silminState->H2Obuffer = TRUE;
                        }*/
#endif

                }
            }

    }
#ifndef BATCH_VERSION
            workProcData->active = TRUE;
#endif

            return FALSE;
            /* ------------------------------------------------------------------------ */
        case CONVERGENCE_TEST:

            iterQuad = 0;

            updateBulkADB();
#ifndef BATCH_VERSION
            updateSolidADB((double *) NULL, (double *) NULL);
            updateCompADB();
            updateStatusADB(STATUS_ADB_INDEX_MASS_LIQUID, &(silminState->liquidMass));
            updateStatusADB(STATUS_ADB_INDEX_MASS_SOLID, &(silminState->solidMass));
            if(silminState->fo2Path == FO2_NONE) updateStatusADB(STATUS_ADB_INDEX_LOGFO2, &(silminState->fo2));

            wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...Checking saturation state of potential solids.\n");
            workProcData->active = TRUE;
#elif defined(PHMELTS_ADJUSTMENTS)
            printf("...Checking saturation state of potential solids.\n");
            /* PUT CALL TO SUBTHERMO_CALC IN HERE */
            /* NO NEED TO UPDATE fO2 */
            (void) getSystemProperties(silminState, TRUE);
#else
            fprintf(stderr, "...Checking saturation state of potential solids.\n");
#endif

            curStep++;
            return FALSE;
            /* ------------------------------------------------------------------------ */
        case VERIFY_SATURATION:

            if ((hasSupersaturation = evaluateSaturationState((silminState->ySol), (silminState->yLiq)))) curStep = ADD_PHASE;
            else if ((hasSupersaturation = checkForCoexistingSolids())) {
                curStep = PROJECT_CONSTRAINTS;
#ifndef BATCH_VERSION
                wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...One of the solid phases has undergone phase separation.\n");
#elif defined(ALPHAMELTS_UPDATE_SYSTEM)
                printf("...One of the solid phases has undergone phase separation.\n");
#else
                fprintf(stderr, "...One of the solid phases has undergone phase separation.\n");
#endif
                if (silminState->fo2Path != FO2_NONE && hasLiquid) {
                    double *moles = (double *) malloc((size_t) nc*sizeof(double));
                    for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                        for (i=0; i<nc; i++) {
                            for (j=0, moles[i]=0.0; j<nlc; j++) moles[i] += (silminState->liquidComp)[nl][j]*(liquid[j].liqToOx)[i];
                            (silminState->bulkComp)[i] -= moles[i];
                        }
                        conLiq(FIRST | SEVENTH, FIRST, silminState->T, silminState->P, moles, NULL, NULL, NULL, NULL, NULL, &(silminState->fo2));
                        for (i=0; i<nc; i++) (silminState->bulkComp)[i] += moles[i];
                        for (i=0; i<nlc; i++) for (j=0, (silminState->liquidComp)[nl][i]=0.0; j<nc; j++) (silminState->liquidComp)[nl][i] += moles[j]*(bulkSystem[j].oxToLiq)[i];
                    }
                    free(moles);
                } else if (silminState->fo2Path != FO2_NONE && !hasLiquid) {
                    double muO2 = silminState->fo2*(R*silminState->T*log(10.0));
                    if (!subsolidusmuO2(0, &muO2, NULL,  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)) {
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
                iterQuad = 0; acceptable = FALSE; bestrNorm = 1.0;
#ifdef PHMELTS_ADJUSTMENTS
            } else if (hasLiquid && (silminState->nLiquidCoexist < silminState->incLiquids) && (hasSupersaturation = checkForCoexistingLiquids())) {
#else
            } else if (hasLiquid && silminState->multipleLiqs && (hasSupersaturation = checkForCoexistingLiquids())) {
#endif
                curStep = PROJECT_CONSTRAINTS;
#ifndef BATCH_VERSION
                wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...A liquid has undergone phase separation.\n");
#elif defined(ALPHAMELTS_UPDATE_SYSTEM)
                printf("...A liquid has undergone phase separation.\n");
#else
                fprintf(stderr, "...A liquid has undergone phase separation.\n");
#endif
                if (silminState->fo2Path != FO2_NONE) {
                    double *moles = (double *) malloc((size_t) nc*sizeof(double));
                    for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                        for (i=0; i<nc; i++) {
                            for (j=0, moles[i]=0.0; j<nlc; j++) moles[i] += (silminState->liquidComp)[nl][j]*(liquid[j].liqToOx)[i];
                            (silminState->bulkComp)[i] -= moles[i];
                        }
                        conLiq(FIRST | SEVENTH, FIRST, silminState->T, silminState->P, moles, NULL, NULL, NULL, NULL, NULL, &(silminState->fo2));
                        for (i=0; i<nc; i++) (silminState->bulkComp)[i] += moles[i];
                        for (i=0; i<nlc; i++) for (j=0, (silminState->liquidComp)[nl][i]=0.0; j<nc; j++) (silminState->liquidComp)[nl][i] += moles[j]*(bulkSystem[j].oxToLiq)[i];
                    }
                    free(moles);
                }
                iterQuad = 0; acceptable = FALSE; bestrNorm = 1.0;
            } else {
#ifndef BATCH_VERSION
                wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "<> Stable liquid solid assemblage achieved.\n");
#elif defined(ALPHAMELTS_UPDATE_SYSTEM)
                printf("<> Stable liquid solid assemblage achieved.\n");
#else
                fprintf(stderr, "<> Stable liquid solid assemblage achieved.\n");
#endif
                curStep++;
            }

#ifndef BATCH_VERSION
            updateSolidADB((silminState->ySol), (silminState->yLiq));
            workProcData->active = TRUE;
#endif

            return FALSE;
            /* ------------------------------------------------------------------------ */
        case OUTPUT_RESULTS:

#ifndef BATCH_VERSION
            updateUserGraphGW();
#endif

#ifdef PHMELTS_ADJUSTMENTS

            /* Reset fracSComp and fracLComp (previousState is used for accumulation instead) */
            if (silminState->fractionateSol || silminState->fractionateFlu || silminState->fractionateLiq) (silminState->fracMass) = 0.0;
            if ((silminState->fractionateSol || silminState->fractionateFlu) && silminState->fracSComp != NULL) {
                for (i=0; i<npc; i++) if (solids[i].type == PHASE) {
                    if ((silminState->nFracCoexist)[i] > 0) {
                        int nf = (silminState->nFracCoexist)[i];
                        for (j=0; j<nf; j++) (silminState->fracSComp)[i][j] = 0.0;
                        if (solids[i].na > 1) for (j=0; j<solids[i].na; j++) {
                            for (k=0; k<nf; k++) (silminState->fracSComp)[i+1+j][k] = 0.0;
                        }
                    }
                }
            }
            if (silminState->fractionateLiq && silminState->fracLComp != NULL) {
                /* SEPARATE OUT MULTIPLE LIQUIDS (JUST LIKE SOLIDS) */
                for (i=0; i<nlc; i++) {
                    (silminState->fracLComp)[i] = 0.0;
                }
            }

            if (silminState->txtOutput != TEXT_NONE && silminState->txtOutput < TEXT_ALPHA_0) {
                if ((silminState->fractionateSol || silminState->fractionateFlu) && silminState->fracSComp == NULL) {
                    silminState->fracSComp    = (double **) calloc((unsigned) npc, sizeof(double *));
                    silminState->nFracCoexist = (int *) calloc((unsigned) npc, sizeof(int));
                }
                if (silminState->fractionateLiq && silminState->fracLComp == NULL) {
                    silminState->fracLComp = (double *) calloc((unsigned) nlc, sizeof(double));
                }
                /* NOT NEEDED? */
                /*if  ((silminState->fractionateSol || silminState->fractionateFlu || silminState->fractionateLiq) && previousSilminState == NULL) {
                    previousSilminState = copySilminStateStructure(silminState, previousSilminState);
                }*/



                (void) putOutputDataToFile((char *) NULL);
                if (additionalOutput != NULL) (*additionalOutput)(addOutputFileName);
            }

    /*******************/
    //if (strstr(silminInputData.name, ".xml")   != NULL) putSequenceDataToXmlFile(TRUE);
    // NEED TO CLOSE OUT FILE
    //(void) putSequenceDataToXmlFile(TRUE);

#else
#ifndef DO_NOT_PRODUCE_OUTPUT_FILES
            (void) putOutputDataToFile((char *) NULL);
            if (additionalOutput != NULL) (*additionalOutput)(addOutputFileName);
#endif
#endif

#ifndef BATCH_VERSION
            workProcData->active = TRUE;
#endif

            curStep++;
            return FALSE;
            /* ------------------------------------------------------------------------ */
        case UPDATE_SYSTEM:

#ifndef BATCH_VERSION
            if (oldErrorHandler != NULL && signal(SIGFPE, oldErrorHandler) == SIG_ERR)
                wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...Error in installing old signal handler.\n");
            oldErrorHandler = NULL;
#else
            (void) gsl_set_error_handler(NULL);
#endif

#ifndef PHMELTS_ADJUSTMENTS
            /* Solid Phase Fractionation */
#ifndef BATCH_VERSION
            if ((silminState->fractionateSol || silminState->fractionateFlu) && !hasLiquid) wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name,
                                                                                                    "...Cannot do solid/fluid fractionation without a liquid phase.\n");
#else
            if ((silminState->fractionateSol || silminState->fractionateFlu) && !hasLiquid) fprintf(stderr, "...Cannot do solid/fluid fractionation without a liquid phase.\n");
#endif

            if ((silminState->fractionateSol || silminState->fractionateFlu) && hasLiquid) {
                double *m = (double *) malloc((size_t) nlc*sizeof(double));
                double *r = (double *) malloc((size_t) nlc*sizeof(double));
                for (i=0; i<npc; i++) if (solids[i].type == PHASE) {
                    if ((silminState->nSolidCoexist)[i] > (silminState->nFracCoexist)[i]) {
                        int ns = (silminState->nSolidCoexist)[i];
                        int nf = (silminState->nFracCoexist)[i];
                        (silminState->nFracCoexist)[i] = ns;
                        if (nf == 0) {
                            (silminState->fracSComp)[i] = (double *) calloc((size_t) ns, sizeof(double));
                            if (solids[i].na > 1) for (j=0; j<solids[i].na; j++) (silminState->fracSComp)[i+1+j] = (double *) calloc((size_t) ns, sizeof(double));
                        } else {
                            (silminState->fracSComp)[i] = (double *) REALLOC((silminState->fracSComp)[i], (size_t) ns*sizeof(double));
                            for (j=nf; j<ns; j++) (silminState->fracSComp)[i][j] = 0.0;
                            if (solids[i].na > 1) for (j=0; j<solids[i].na; j++) {
                                (silminState->fracSComp)[i+1+j] = (double *) REALLOC((silminState->fracSComp)[i+1+j], (size_t) ns*sizeof(double));
                                for (k=nf; k<ns; k++) (silminState->fracSComp)[i+1+j][k] = 0.0;
                            }
                        }
                    }
                }
                int haveWater = ((calculationMode == MODE__MELTS) || (calculationMode == MODE_pMELTS));
                for (i=0; i<npc; i++) {
                    if ( haveWater &&  silminState->fractionateSol && !silminState->fractionateFlu && !strcmp((char *) solids[i].label, "water")) continue;
                    if ( haveWater && !silminState->fractionateSol &&  silminState->fractionateFlu &&  strcmp((char *) solids[i].label, "water")) continue;
                    if (!haveWater &&  silminState->fractionateSol && !silminState->fractionateFlu && !strcmp((char *) solids[i].label, "fluid")) continue;
                    if (!haveWater && !silminState->fractionateSol &&  silminState->fractionateFlu &&  strcmp((char *) solids[i].label, "fluid")) continue;
                    for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
                        if (solids[i].na == 1) {
                            (silminState->fracSComp)[i][ns] += (silminState->solidComp)[i][ns]-MASSIN;
                            if (silminState->fo2Path != FO2_NONE) silminState->oxygen -= (oxygen.solToOx)[i]*((silminState->solidComp)[i][ns]-MASSIN);
                            silminState->fracMass += ((silminState->solidComp)[i][ns]-MASSIN)*solids[i].mw;
                            for (j=0; j<nc; j++) (silminState->bulkComp)[j] -= (solids[i].solToOx)[j]*((silminState->solidComp)[i][ns]-MASSIN);

                            /* Subtract off H, S or V if appropriate                          */
                            if (silminState->isenthalpic && (silminState->refEnthalpy != 0.0))
                                silminState->refEnthalpy -= ((silminState->solidComp)[i][ns]-MASSIN)*(solids[i].cur).h;
                            if (silminState->isentropic && (silminState->refEntropy != 0.0))
                                silminState->refEntropy -= ((silminState->solidComp)[i][ns]-MASSIN)*(solids[i].cur).s;
                            if (silminState->isochoric && (silminState->refVolume != 0.0))
                                silminState->refVolume -= ((silminState->solidComp)[i][ns]-MASSIN)*(solids[i].cur).v;

                            (silminState->solidComp)[i][ns] = MASSIN;
                        } else {
                            double moleF, totalMoles=0.0;
                            (silminState->fracSComp)[i][ns] += (silminState->solidComp)[i][ns] - MASSIN;
                            for (j=0; j<solids[i].na; j++) {
                                moleF = (silminState->solidComp)[i+1+j][ns]/(silminState->solidComp)[i][ns];
                                m[j] = (silminState->solidComp)[i+1+j][ns] - MASSIN*moleF;
                                totalMoles += m[j];
                                (silminState->fracSComp)[i+1+j][ns] += m[j];
                                if (silminState->fo2Path != FO2_NONE) silminState->oxygen -= (oxygen.solToOx)[i+1+j]*m[j];
                                silminState->fracMass += m[j]*solids[i+1+j].mw;
                                for (k=0; k<nc; k++) (silminState->bulkComp)[k] -= (solids[i+1+j].solToOx)[k]*m[j];
                                (silminState->solidComp)[i+1+j][ns] = MASSIN*moleF;

                                /* Subtract off H, S or V if appropriate                        */
                                if (silminState->isenthalpic && (silminState->refEnthalpy != 0.0)) silminState->refEnthalpy -= m[j]*(solids[i+1+j].cur).h;
                                if (silminState->isentropic && (silminState->refEntropy != 0.0))   silminState->refEntropy  -= m[j]*(solids[i+1+j].cur).s;
                                if (silminState->isochoric && (silminState->refVolume != 0.0))     silminState->refVolume   -= m[j]*(solids[i+1+j].cur).v;
                            }
                            (silminState->solidComp)[i][ns] = MASSIN;

                            /* Subtract off H, S or V if appropriate                          */
                            if (silminState->isenthalpic && (silminState->refEnthalpy != 0.0)) {
                                double enthalpy;
                                (*solids[i].convert)(SECOND, THIRD, silminState->T,silminState->P, NULL, m, r, NULL,  NULL, NULL, NULL, NULL);
                                (*solids[i].hmix)(FIRST, silminState->T, silminState->P, r, &enthalpy);
                                silminState->refEnthalpy -= totalMoles*enthalpy;
                            }
                            if (silminState->isentropic && (silminState->refEntropy != 0.0)) {
                                double entropy;
                                (*solids[i].convert)(SECOND, THIRD,silminState->T,silminState->P, NULL, m, r, NULL, NULL, NULL, NULL, NULL);
                                (*solids[i].smix)(FIRST, silminState->T, silminState->P, r, &entropy, (double *) NULL, (double **) NULL);
                                silminState->refEntropy  -= totalMoles*entropy;
                            }
                            if (silminState->isochoric && (silminState->refVolume != 0.0)) {
                                double volume;
                                (*solids[i].convert)(SECOND, THIRD, silminState->T,silminState->P, NULL, m, r, NULL, NULL, NULL, NULL, NULL);
                                (*solids[i].vmix)(FIRST, silminState->T, silminState->P, r, &volume, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
                                silminState->refVolume   -= totalMoles*volume;
                            }

                        }
                    }
                }

#ifndef BATCH_VERSION
                updateStatusADB(STATUS_ADB_INDEX_MASS_FRAC, &(silminState->fracMass));
#endif

                for (i=0; i<nc; i++) {
                    if ((silminState->bulkComp)[i] != 0.0 && (silminState->bulkComp)[i] <  MASSOUT && bulkSystem[i].type != FE2O3) {
#ifndef BATCH_VERSION
                        wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "  Moles of %5.5s in system (%g) < %g\n.", bulkSystem[i].label, (silminState->bulkComp)[i], MASSOUT);
#else
                        fprintf(stderr, "  Moles of %5.5s in system (%g) < %g\n.", bulkSystem[i].label, (silminState->bulkComp)[i], MASSOUT);
#endif
                        (silminState->bulkComp)[i] = 0.0;
                        for (j=0; j<nlc; j++) if ((liquid[j].liqToOx)[i] != 0.0) {
                            for (nl=0; nl<silminState->nLiquidCoexist; nl++) (silminState->liquidComp)[nl][j] = 0.0;
#ifndef BATCH_VERSION
                            wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "    Moles of %s in liquid(s) set to zero.\n", liquid[j].label);
#else
                            fprintf(stderr, "    Moles of %s in liquid(s) set to zero.\n", liquid[j].label);
#endif
                        }
                        for (j=0; j<npc; j++) {
                            for (ns=0; ns<(silminState->nSolidCoexist)[j]; ns++) {
                                if (solids[j].na == 1) {
                                    if ((solids[j].solToOx)[i] != 0.0) {
                                        (silminState->solidComp)[j][ns] = 0.0;
#ifndef BATCH_VERSION
                                        wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "    Moles of %s in solid set to zero.\n", solids[j].label);
#else
                                        fprintf(stderr, "    Moles of %s in solid set to zero.\n", solids[j].label);
#endif
                                    }
                                } else {
                                    for (k=0; k<solids[j].na; k++) {
                                        if ((solids[j+1+k].solToOx)[i] != 0.0) {
                                            (silminState->solidComp)[j][ns] -= (silminState->solidComp)[j+1+k][ns];
                                            (silminState->solidComp)[j+1+k][ns] = 0.0;
#ifndef BATCH_VERSION
                                            wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "    Moles of %s in %s solid set to zero.\n", solids[j+1+k].label, solids[j].label);
#else
                                            fprintf(stderr, "    Moles of %s in %s solid set to zero.\n", solids[j+1+k].label, solids[j].label);
#endif
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                free(m);
                free(r);
            }

            /* Liquid Phase Fractionation */
#ifndef BATCH_VERSION
            if (silminState->fractionateLiq && !hasLiquid) wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name,
                                   "...Cannot do liquid fractionation without a liquid phase.\n");
#else
            if (silminState->fractionateLiq && !hasLiquid) fprintf(stderr, "...Cannot do liquid fractionation without a liquid phase.\n");
#endif

            if (silminState->fractionateLiq && hasLiquid) {
                double *m = (double *) malloc((size_t) nlc*sizeof(double));
                double *r = (double *) malloc((size_t) nlc*sizeof(double));
                for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
                    double refMoles, totalMoles;
                    for (i=0, refMoles=0.0; i<nlc; i++) refMoles += (silminState->liquidComp)[nl][i];

                    for (i=0, totalMoles=0.0; i<nlc; i++) {
                        if (((silminState->liquidComp)[nl][i] != 0.0) && (refMoles != 0.0)) {
                            double mw;
                            double moleF = (silminState->liquidComp)[nl][i]/refMoles;

                            for (j=0, mw = 0.0; j<nc; j++) mw += (liquid[i].liqToOx)[j]*bulkSystem[j].mw;
                            m[i] = (silminState->liquidComp)[nl][i] - MASSIN*moleF;
                            totalMoles += m[i];
                            (silminState->fracLComp)[i] += m[i];
                            if (silminState->fo2Path != FO2_NONE) silminState->oxygen -= (oxygen.liqToOx)[i]*m[i];
                            silminState->fracMass += m[i]*mw;
                            for (j=0; j<nc; j++) (silminState->bulkComp)[j] -= (liquid[i].liqToOx)[j]*m[i];
                            (silminState->liquidComp)[nl][i] = MASSIN*moleF;

                            /* Subtract off H, S or V if appropriate			    */
                            if (silminState->isenthalpic && (silminState->refEnthalpy != 0.0)) silminState->refEnthalpy -= m[i]*(liquid[i].cur).h;
                            if (silminState->isentropic  && (silminState->refEntropy  != 0.0)) silminState->refEntropy  -= m[i]*(liquid[i].cur).s;
                            if (silminState->isochoric   && (silminState->refVolume   != 0.0)) silminState->refVolume	-= m[i]*(liquid[i].cur).v;
                        } else m[i] = 0.0;
                    }

                    /* Subtract off H, S or V if appropriate			  */
                    if (silminState->isenthalpic && (silminState->refEnthalpy != 0.0)) {
                        double enthalpy;
                        conLiq (SECOND, THIRD, silminState->T,silminState->P, NULL, m, r, NULL,  NULL, NULL, NULL);
                        hmixLiq(FIRST, silminState->T, silminState->P, r, &enthalpy, NULL);
                        silminState->refEnthalpy -= totalMoles*enthalpy;
                    }
                    if (silminState->isentropic && (silminState->refEntropy != 0.0)) {
                        double entropy;
                        conLiq (SECOND, THIRD,silminState->T,silminState->P, NULL, m, r, NULL, NULL, NULL, NULL);
                        smixLiq(FIRST, silminState->T, silminState->P, r, &entropy, NULL, NULL, NULL);
                        silminState->refEntropy  -= totalMoles*entropy;
                    }
                    if (silminState->isochoric && (silminState->refVolume != 0.0)) {
                        double volume;
                        conLiq (SECOND, THIRD, silminState->T,silminState->P, NULL, m, r, NULL, NULL, NULL, NULL);
                        vmixLiq(FIRST, silminState->T, silminState->P, r, &volume, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
                        silminState->refVolume   -= totalMoles*volume;
                    }

                }

#ifndef BATCH_VERSION
                updateStatusADB(STATUS_ADB_INDEX_MASS_FRAC, &(silminState->fracMass));
#endif

                for (i=0; i<nc; i++) {
                    if ((silminState->bulkComp)[i] != 0.0 && (silminState->bulkComp)[i] <  MASSOUT && bulkSystem[i].type != FE2O3) {
#ifndef BATCH_VERSION
                        wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "  Moles of %5.5s in system (%g) < %g\n.", bulkSystem[i].label, (silminState->bulkComp)[i], MASSOUT);
#else
                        fprintf(stderr, "  Moles of %5.5s in system (%g) < %g\n.", bulkSystem[i].label, (silminState->bulkComp)[i], MASSOUT);
#endif
                        (silminState->bulkComp)[i] = 0.0;
                        for (j=0; j<nlc; j++) if ((liquid[j].liqToOx)[i] != 0.0) {
                            for (nl=0; nl<silminState->nLiquidCoexist; nl++) (silminState->liquidComp)[nl][j] = 0.0;
#ifndef BATCH_VERSION
                            wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "    Moles of %s in liquid(s) set to zero.\n", liquid[j].label);
#else
                            fprintf(stderr, "    Moles of %s in liquid(s) set to zero.\n", liquid[j].label);
#endif
                        }
                        for (j=0; j<npc; j++) {
                            for (ns=0; ns<(silminState->nSolidCoexist)[j]; ns++) {
                                if (solids[j].na == 1) {
                                    if ((solids[j].solToOx)[i] != 0.0) {
                                        (silminState->solidComp)[j][ns] = 0.0;
#ifndef BATCH_VERSION
                                        wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "    Moles of %s in solid set to zero.\n", solids[j].label);
#else
                                        fprintf(stderr, "    Moles of %s in solid set to zero.\n", solids[j].label);
#endif
                                    }
                                } else {
                                    for (k=0; k<solids[j].na; k++) {
                                        if ((solids[j+1+k].solToOx)[i] != 0.0) {
                                            (silminState->solidComp)[j][ns] -= (silminState->solidComp)[j+1+k][ns];
                                            (silminState->solidComp)[j+1+k][ns] = 0.0;
#ifndef BATCH_VERSION
                                            wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "    Moles of %s in %s solid set to zero.\n", solids[j+1+k].label, solids[j].label);
#else
                                            fprintf(stderr, "    Moles of %s in %s solid set to zero.\n", solids[j+1+k].label, solids[j].label);
#endif
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                free(m);
                free(r);
            }

            stateChange = FALSE;

            /* Changing T ? */
            if (fabs(silminState->dspTstart - silminState->dspTstop) >=  (silminState->dspTinc != 0.0 ? silminState->dspTinc : 0.001)
                && !(silminState->isenthalpic && (silminState->refEnthalpy != 0.0))
                && !(silminState->isentropic  && (silminState->refEntropy  != 0.0))) {
                stateChange = TRUE;
#ifndef BATCH_VERSION
                workProcData->mode = TRUE;
                if (silminState->dspTstart - silminState->dspTstop < 0.0) tpValues[TP_PADB_INDEX_T_INITIAL].value += silminState->dspTinc;
                else                                                      tpValues[TP_PADB_INDEX_T_INITIAL].value -= silminState->dspTinc;
#else
                if (silminState->dspTstart - silminState->dspTstop < 0.0) {
                    silminState->T += silminState->dspTinc; silminState->dspTstart += silminState->dspTinc;
                } else {
                    silminState->T -= silminState->dspTinc; silminState->dspTstart -= silminState->dspTinc;
                }
#endif
                /* Changing H ? */
            } else if ((fabs(silminState->dspHstop) > 0.0) &&
             (fabs(silminState->dspHstop - silminState->refEnthalpy) >= (silminState->dspHinc != 0.0 ? fabs(silminState->dspHinc) : 0.001)) &&
             (silminState->isenthalpic && (silminState->refEnthalpy != 0.0))) {
                stateChange = TRUE;
                silminState->refEnthalpy += silminState->dspHinc;
                correctTforChangeInEnthalpy();
#ifndef BATCH_VERSION
                workProcData->mode = TRUE;
                tpValues[TP_PADB_INDEX_T_INITIAL].value = silminState->T - 273.15;
#endif
            } else if ((silminState->dspHstop == 0.0) &&
             (fabs(silminState->dspPstart - silminState->dspPstop) >=  (silminState->dspPinc != 0.0 ? silminState->dspPinc : 0.001)) &&
             (silminState->isenthalpic && (silminState->refEnthalpy != 0.0))) {
                stateChange = TRUE;
                silminState->refEnthalpy += silminState->dspHinc;
                correctTforChangeInEnthalpy();
#ifndef BATCH_VERSION
                workProcData->mode = TRUE;
                tpValues[TP_PADB_INDEX_T_INITIAL].value = silminState->T - 273.15;
#endif
        /* Could get into loop if dH/dP = 0.0 (e.g. isenthalpic degassing) or for isenthalpic assimilation once all material assimilated */
        /* (fabs(silminState->dspTstart - silminState->dspTstop) >=  (silminState->dspTinc != 0.0 ? silminState->dspTinc : 0.001)) */
            } else if ((silminState->dspHstop == 0.0) && (silminState->dspHinc != 0.0) &&
             (silminState->isenthalpic && (silminState->refEnthalpy != 0.0))) {
            int changeH = (silminState->dspHinc < 0.0) ? (silminState->dspTstart - silminState->dspTstop < 0.0) : (silminState->dspTstart - silminState->dspTstop > 0.0);
            if (changeH) {
                stateChange = TRUE;
                silminState->refEnthalpy += silminState->dspHinc;
                correctTforChangeInEnthalpy();
#ifndef BATCH_VERSION
                workProcData->mode = TRUE;
                tpValues[TP_PADB_INDEX_T_INITIAL].value = silminState->T - 273.15;
#endif
            }
                /* Changing S ? */
            } else if ((fabs(silminState->dspSstop) > 0.0) &&
             (fabs(silminState->dspSstop - silminState->refEntropy) >= (silminState->dspSinc != 0.0 ? fabs(silminState->dspSinc) : 0.001)) &&
             (silminState->isentropic && (silminState->refEntropy != 0.0))) {
                stateChange = TRUE;
                silminState->refEntropy += silminState->dspSinc;
                correctTforChangeInEntropy();
#ifndef BATCH_VERSION
                workProcData->mode = TRUE;
                tpValues[TP_PADB_INDEX_T_INITIAL].value = silminState->T - 273.15;
#endif
            } else if ((silminState->dspSstop == 0.0) &&
             (fabs(silminState->dspPstart - silminState->dspPstop) >=  (silminState->dspPinc != 0.0 ? silminState->dspPinc : 0.001)) &&
             (silminState->isentropic  && (silminState->refEntropy  != 0.0))) {
                stateChange = TRUE;
                silminState->refEntropy += silminState->dspSinc;
                correctTforChangeInEntropy();
#ifndef BATCH_VERSION
                workProcData->mode = TRUE;
                tpValues[TP_PADB_INDEX_T_INITIAL].value = silminState->T - 273.15;
#endif
        /* (fabs(silminState->dspTstart - silminState->dspTstop) >=  (silminState->dspTinc != 0.0 ? silminState->dspTinc : 0.001)) */
            } else if ((silminState->dspSstop == 0.0) && (silminState->dspSinc != 0.0) &&
             (silminState->isentropic  && (silminState->refEntropy  != 0.0))) {
            int changeS = (silminState->dspSinc > 0.0) ? (silminState->dspTstart - silminState->dspTstop < 0.0) : (silminState->dspTstart - silminState->dspTstop > 0.0);
            if (changeS) {
                stateChange = TRUE;
                silminState->refEntropy += silminState->dspSinc;
                correctTforChangeInEntropy();
#ifndef BATCH_VERSION
                workProcData->mode = TRUE;
                tpValues[TP_PADB_INDEX_T_INITIAL].value = silminState->T - 273.15;
#endif
            }
            }
            /* Changing P ? */
            if (fabs(silminState->dspPstart - silminState->dspPstop) >=  (silminState->dspPinc != 0.0 ? silminState->dspPinc : 0.001)
                && !(silminState->isochoric && (silminState->refVolume != 0.0))) {
                stateChange = TRUE;
#ifndef BATCH_VERSION
                workProcData->mode = TRUE;
                if (silminState->dspPstart - silminState->dspPstop < 0.0) tpValues[TP_PADB_INDEX_P_INITIAL].value += silminState->dspPinc;
                else                                                      tpValues[TP_PADB_INDEX_P_INITIAL].value -= silminState->dspPinc;
#else
                if (silminState->dspPstart - silminState->dspPstop < 0.0) {
                    silminState->P += silminState->dspPinc; silminState->dspPstart += silminState->dspPinc;
                } else {
                    silminState->P -= silminState->dspPinc; silminState->dspPstart -= silminState->dspPinc;
                }
#endif
                /* Changing V ? */
            } else if ((fabs(silminState->dspVstop) > 0.0) &&
             (fabs(silminState->dspVstop - silminState->refVolume) >= (silminState->dspVinc != 0.0 ? fabs(silminState->dspVinc) : 0.0001)) &&
             (silminState->isochoric && (silminState->refVolume != 0.0))) {
                stateChange = TRUE;
                silminState->refVolume += silminState->dspVinc;
                correctPforChangeInVolume();
#ifndef BATCH_VERSION
                workProcData->mode = TRUE;
                tpValues[TP_PADB_INDEX_P_INITIAL].value = silminState->P;
#endif
                /*    } else if (fabs(silminState->dspPstart - silminState->dspPstop) >=  (silminState->dspPinc != 0.0 ? silminState->dspPinc : 0.001) FIX-fO2 */
            } else if ((silminState->dspVstop == 0.0) &&
             (fabs(silminState->dspTstart - silminState->dspTstop) >=  (silminState->dspTinc != 0.0 ? silminState->dspTinc : 0.001)) &&
             (silminState->isochoric && (silminState->refVolume != 0.0))) {
                stateChange = TRUE;
                silminState->refVolume += silminState->dspVinc;
                correctPforChangeInVolume();
#ifndef BATCH_VERSION
                workProcData->mode = TRUE;
                tpValues[TP_PADB_INDEX_P_INITIAL].value = silminState->P;
#endif
        /* (fabs(silminState->dspPstart - silminState->dspPstop) >=  (silminState->dspPinc != 0.0 ? silminState->dspPinc : 0.001) */
            } else if ((silminState->dspVstop == 0.0) && (silminState->dspVinc != 0.0) &&
             (silminState->isochoric && (silminState->refVolume != 0.0))) {
            int changeV = (silminState->dspVinc < 0.0) ? (silminState->dspPstart - silminState->dspPstop < 0.0) : (silminState->dspPstart - silminState->dspPstop > 0.0);
            if (changeV) {
                stateChange = TRUE;
                silminState->refVolume += silminState->dspVinc;
                correctPforChangeInVolume();
#ifndef BATCH_VERSION
                workProcData->mode = TRUE;
                tpValues[TP_PADB_INDEX_P_INITIAL].value = silminState->P;
#endif
            }
            }
            /* Assimilation ? - add assimilant to the first liquid, if a liquid is present */
            if (silminState->assimilate && silminState->assimMass < silminState->dspAssimMass) {
                double fraction = 1.0/silminState->dspAssimInc;

                stateChange = TRUE;
                silminState->assimMass += fraction*silminState->dspAssimMass;
                if (silminState->isenthalpic && (silminState->refEnthalpy != 0.0)) silminState->refEnthalpy += fraction*(silminState->assimTD).h;
                if (silminState->isentropic && (silminState->refEntropy != 0.0))   silminState->refEntropy += fraction*(silminState->assimTD).s;
                if (silminState->isochoric && (silminState->refVolume != 0.0))     silminState->refVolume += fraction*(silminState->assimTD).v;
                silminState->assimInc++;
#ifndef BATCH_VERSION
                workProcData->mode = TRUE;
                updateStatusADB(STATUS_ADB_INDEX_MASS_ASSIM, &(silminState->assimMass));
#endif

                for (i=0; i<npc; i++) for (ns=0; ns<(silminState->nAssimComp)[i]; ns++) {
                    if (solids[i].type == PHASE && solids[i].na == 1) {
                        for (k=0; k<nc; k++) (silminState->bulkComp)[k] += fraction*(silminState->assimComp)[i][ns]*(solids[i].solToOx)[k];
                        if (hasLiquid) for (k=0; k<nlc; k++) (silminState->liquidComp)[0][k] += fraction*(silminState->assimComp)[i][ns]*(solids[i].solToLiq)[k];
                        if (!hasLiquid) {  /* if liquid is absent, assimilate as solid */
                            if (!(silminState->nSolidCoexist)[i]) (silminState->nSolidCoexist)[i] = 1;
                            (silminState->solidComp)[i][0] += fraction*(silminState->assimComp)[i][ns];
                        }
                        if (silminState->fo2Path != FO2_NONE) silminState->oxygen += (oxygen.solToOx)[i]*fraction*(silminState->assimComp)[i][ns];
                    } else if (solids[i].type == PHASE && solids[i].na > 1) {
                        for (j=0; j<solids[i].na; j++) {
                            for (k=0; k<nc; k++) (silminState->bulkComp)[k] += fraction*(silminState->assimComp)[i+1+j][ns]*(solids[i+1+j].solToOx)[k];
                            if (hasLiquid) for (k=0; k<nlc; k++) (silminState->liquidComp)[0][k] += fraction*(silminState->assimComp)[i+1+j][ns]*(solids[i+1+j].solToLiq)[k];
                            if (!hasLiquid) {  /* if liquid is absent, assimilate as solid -- risky ? */
                                if (!(silminState->nSolidCoexist)[i]) (silminState->nSolidCoexist)[i] = 1;
                                (silminState->solidComp)[i+1+j][0] += fraction*(silminState->assimComp)[i+1+j][ns];
                            }
                            if (silminState->fo2Path != FO2_NONE) silminState->oxygen += (oxygen.solToOx)[i+1+j]*fraction*(silminState->assimComp)[i+1+j][ns];
                        }
                    }
                }
                for (ns=0; ns<(silminState->nAssimComp)[npc+1]; ns++) {
                    if (!hasLiquid) hasLiquid = TRUE;
                    for (i=0; i<nlc; i++) {
                        (silminState->liquidComp)[0][i] += fraction*(silminState->assimComp)[npc+i][ns];
                        for (j=0; j<nc; j++) (silminState->bulkComp)[j] += fraction*(silminState->assimComp)[npc+i][ns]*(liquid[i].liqToOx)[j];
                        if (silminState->fo2Path != FO2_NONE) silminState->oxygen += (oxygen.liqToOx)[i]*fraction*(silminState->assimComp)[npc+i][ns];
                    }
                    for (i=0,silminState->liquidMass=0.0; i<nc; i++) {
                        for (j=0; j<nlc; j++) silminState->liquidMass += silminState->liquidComp[0][j]*(liquid[j].liqToOx)[i]*bulkSystem[i].mw;
                    }
                }
                if (silminState->isenthalpic && (silminState->refEnthalpy != 0.0)) {
                    correctTforChangeInEnthalpy();
#ifndef BATCH_VERSION
                    tpValues[TP_PADB_INDEX_T_INITIAL].value = silminState->T - 273.15;
#endif
                }
                if (silminState->isentropic && (silminState->refEntropy != 0.0)) {
                    correctTforChangeInEntropy();
#ifndef BATCH_VERSION
                    tpValues[TP_PADB_INDEX_T_INITIAL].value = silminState->T - 273.15;
#endif
                }
                if (silminState->isochoric && (silminState->refVolume != 0.0)) {
                    correctPforChangeInVolume();
#ifndef BATCH_VERSION
                    tpValues[TP_PADB_INDEX_P_INITIAL].value = silminState->P;
#endif
                }
            }

#else
            stateChange = FALSE;
#endif /* PHMELTS_ADJUSTMENTS */

#ifdef PHMELTS_ADJUSTMENTS
        *saveT = silminState->T; *saveP = silminState->P;
        if(hasLiquid) for (i=0; i<nlc; i++) saveLiqComp[i] = silminState->liquidComp[0][i];
        else for (i=0; i<nlc; i++) saveLiqComp[i] = 0.0;
#endif

#ifndef BATCH_VERSION
            workProcData->active = stateChange;
#else
            meltsStatus.status = SILMIN_SUCCESS;
#endif

            curStep = 0;
            return (!stateChange);
            /* ======================================================================== */
    } /* end switch */

#ifdef PHMELTS_ADJUSTMENTS
    *saveT = silminState->T; *saveP = silminState->P;
    if(hasLiquid) for (i=0; i<nlc; i++) saveLiqComp[i] = silminState->liquidComp[0][i];
    else for (i=0; i<nlc; i++) saveLiqComp[i] = 0.0;
#endif

#ifdef BATCH_VERSION
    meltsStatus.status = SILMIN_SUCCESS;
#endif
    return TRUE;
}

#ifndef BATCH_VERSION
static void newErrorHandler(int sig)
{
    XmString csString1, csString2, csString3;
    XEvent event;
    XtAppContext app;

    csString1 = XmStringCreateLtoR("A FATAL ARITHMETIC ERROR HAS OCCURRED.\n",          "ISO8859-1");
    csString2 = XmStringCreateLtoR("The program has been left in an unstable state.\n", "ISO8859-1");
    csString3 = XmStringConcat(csString1, csString2);
    XmStringFree(csString1); XmStringFree(csString2);

    csString1 = XmStringCreateLtoR("You must try to save your work, and restart the calculation.\n", "ISO8859-1");
    csString2 = XmStringConcat(csString3, csString1);
    XmStringFree(csString3); XmStringFree(csString1);

    csString1 = XmStringCreateLtoR("Please report the sequence of events which led to this disaster.", "ISO8859-1");
    csString3 = XmStringConcat(csString2, csString1);
    XmStringFree(csString1); XmStringFree(csString2);

    XtVaSetValues(error_message, XmNmessageString, csString3, NULL);
    XtManageChild(error_message);
    XmStringFree(csString3);

    updateBulkADB();
    updateSolidADB((double *) NULL, (double *) NULL);
    updateCompADB();
    updateStatusADB(STATUS_ADB_INDEX_MASS_LIQUID, &(silminState->liquidMass));
    updateStatusADB(STATUS_ADB_INDEX_MASS_SOLID, &(silminState->solidMass));
    if(silminState->fo2Path == FO2_NONE) updateStatusADB(STATUS_ADB_INDEX_LOGFO2, &(silminState->fo2));
    wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "...FATAL arithmetic error detected. Aborting.\n");

    app = XtDisplayToApplicationContext(XtDisplay(topLevel));
    for(;;) {
        XtAppNextEvent(app, &event);
        XtDispatchEvent(&event);
    }

}
#endif /* BATCH_VERSION */

/* end of file SILMIN.C */
