const char *liquidus_ver(void) { return "$Id: liquidus.c,v 1.4 2007/12/22 22:43:30 ghiorso Exp $"; }
/*
    MELTS Source Code: RCS $Log: liquidus.c,v $
    MELTS Source Code: RCS Revision 1.4  2007/12/22 22:43:30  ghiorso
    MELTS Source Code: RCS Fixed error in BM integration in Gibbs.c
    MELTS Source Code: RCS Updated param_struct_data.h file for AGU 2007 xMELTS parameters
    MELTS Source Code: RCS Added support for status file production in MELTS-batch (XML)
    MELTS Source Code: RCS
    MELTS Source Code: RCS Revision 1.3  2007/11/29 05:32:13  ghiorso
    MELTS Source Code: RCS Majorite testMaj corrections.  Fixed "O2" in sol_struct_data.h to "o2".
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
    * Revision 3.9  1997/06/21  22:49:46  ghiorso
    * June 1997 MELTS 3.0.x release
    * (prior to new entropy and regression model being introduced)
    *
    * Revision 3.8  1997/05/03  20:23:25  ghiorso
    * *** empty log message ***
    *
    * Revision 3.7  1997/03/27  17:03:29  ghiorso
    * *** empty log message ***
    *
    * Revision 3.6  1996/09/24  20:33:35  ghiorso
    * Version modified for OSF/1 4.0
    *
    * Revision 3.5  1995/12/09  19:26:38  ghiorso
    * Interface revisions for status box and graphics display
    *
    * Revision 3.4  1995/11/23  22:37:42  ghiorso
    * Final implementation of subsolidus fO2 buffering.
    *
    * Revision 3.3  1995/11/01  22:40:27  ghiorso
    * Implementation of subsolidus options after Asimow.
    * Additional implementation of nepheline solid solutions.
    *
    * Revision 3.2  1995/09/04  20:01:28  ghiorso
    * Update to allow display of bulk composition (in grams) in the text entry
    * fields of the main silmin display. Liquid composition is no longer
    * display here, and is available only through the popup selection.
    *
    * Revision 3.1  1995/08/18  17:47:56  ghiorso
    * MELTS Version 3 - Initial Entry
    *
    */

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Toolkit work proceedure to determine the liquidis temperature of
**      a given bulk composition (file: LIQUIDUS.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  September 7, 1991
**              Original (dummy version)
**      V2.0-1  Mark S. Ghiorso  November 14, 1991
**              Conversion to OSF Motif V1.1.1/X11 Release 4
**      V2.0-2  Mark S. Ghiorso  January 16, 1992
**              First working version
**      V2.0-3  Mark S. Ghiorso  March 23, 1992
**              Removed redundant static qualifier for enum statement
**      V3.0-1  Mark S. Ghiorso  May 11, 1994
**              (1) Added constraint checking on isentropic constraint
**      V4.0-1  Paul D. Asimow  April 26, 1995
**              (1) Changed call to evaluateSaturationState to reflect
**              new arguments
*--
*/

#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef BATCH_VERSION
#include <Xm/Xm.h>
#include "interface.h"            /*Specific external declarations          */
#endif

#include "silmin.h"               /*SILMIN structures include file          */

#ifdef BATCH_VERSION
#include "status.h"               /*Status of calculation                   */
#endif

/*
 *=============================================================================
 * Executable code
 */

#ifndef BATCH_VERSION
Boolean liquidus(XtPointer client_data)
#else
int liquidus(void)
#endif
{
    enum steps {
        CHECK_CONDITION, CHANGE_TP, CHECK_SATURATION, UPDATE_SYSTEM
    };
    static int curStep = 0;
    static int hasSupersaturation, tState = -1;
    static double tInterval;
    int i, j, k, stateChange;

#ifndef BATCH_VERSION
    /*WorkProcData *workProcData = (WorkProcData *) client_data;*/

    /******************************************************************************
   * On entry, check status of calculation and check workProcData.mode for tag:
   *   FALSE   return call without interface modification
   *   TRUE    initial call to invoke a liquidus determination
   ******************************************************************************/

    /*if (workProcData->mode) {
        } */
#else /* BATCH_VERSION */
    if (curStep == 0) curStep = CHECK_CONDITION;
#endif /* BATCH_VERSION */

    /******************************************************************************
   * Step to current phase of the calculation
   ******************************************************************************/

    switch(curStep) {
        /* ======================================================================== */
    case CHECK_CONDITION:
        if ((silminState->solidMass > 0.0) || (silminState->nLiquidCoexist > 1)) {
#ifndef BATCH_VERSION
            XmString csString = XmStringCreateLtoR("Cannot execute Find Liquidus with solids/multiple liquids present.", "ISO8859-1");
            XtVaSetValues(message, XmNmessageString, csString, NULL);
            XtManageChild(message);
            XmStringFree(csString);
            workProcData->active = FALSE;
#else
            meltsStatus.status = LIQUIDUS_MULTIPLE;
            printf("Cannot execute Find Liquidus with solids/multiple liquids present.\n");
#endif
            curStep = 0;
            return TRUE;
        }

#ifndef BATCH_VERSION
        workProcData->active = TRUE;
#endif

        curStep++;
        return FALSE;
        /* ------------------------------------------------------------------------ */
    case CHANGE_TP:

        /* -> Calculate liquid end-member properties                                  */
        for (i=0; i<nlc; i++) gibbs(silminState->T, silminState->P, (char *) liquid[i].label, &(liquid[i].ref),
                                                                &(liquid[i].liq), &(liquid[i].fus), &(liquid[i].cur));

        /* -> Calculate solid  end-member properties                                  */
        for (i=0, j=0; i<npc; i++) {
            if (solids[i].type == PHASE) {
                if ((silminState->incSolids)[j]) {
                    if(solids[i].na == 1) gibbs(silminState->T, silminState->P, (char *) solids[i].label, &(solids[i].ref), NULL, NULL, &(solids[i].cur));
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


        /* -> Redistribute Fe2O3 and FeO in liquid phase to establish buffer at this T and P */
        if (silminState->fo2Path != FO2_NONE) {
            double *moles = (double *) malloc((unsigned) nc*sizeof(double));
            gibbs(silminState->T, silminState->P, "o2", &(oxygen.ref), NULL, NULL, &(oxygen.cur));
            silminState->fo2 = getlog10fo2(silminState->T, silminState->P, silminState->fo2Path);
            for (i=0; i<nc; i++) {
                for (j=0, moles[i]=0.0; j<nlc; j++) moles[i] += (silminState->liquidComp)[0][j]*(liquid[j].liqToOx)[i];
                (silminState->bulkComp)[i] -= moles[i];
            }
            conLiq(FIRST | SEVENTH, FIRST, silminState->T, silminState->P, moles, NULL, NULL, NULL, NULL, NULL, &(silminState->fo2));
            for (i=0; i<nc; i++) (silminState->bulkComp)[i] += moles[i];
            for (i=0; i<nlc; i++) for (j=0, (silminState->liquidComp)[0][i]=0.0; j<nc; j++)
                                                            (silminState->liquidComp)[0][i] += moles[j]*(bulkSystem[j].oxToLiq)[i];
            free(moles);
#ifndef BATCH_VERSION
            updateBulkADB();
            updateStatusADB(STATUS_ADB_INDEX_MASS_LIQUID, &(silminState->liquidMass));
            updateStatusADB(STATUS_ADB_INDEX_LOGFO2, &(silminState->fo2));
        } else {
            updateBulkADB();
            if(silminState->fo2Path == FO2_NONE) updateStatusADB(STATUS_ADB_INDEX_LOGFO2, &(silminState->fo2));
#endif
        }

#ifndef BATCH_VERSION
        workProcData->active = TRUE;
#endif
        curStep++;
        return FALSE;
        /* ------------------------------------------------------------------------ */
    case CHECK_SATURATION:

        if ((silminState->ySol) == NULL) {
            (silminState->ySol) = (double *) malloc((unsigned) npc*sizeof(double));
            (silminState->yLiq) = (double *) malloc((unsigned) nlc*sizeof(double));
        }

        hasSupersaturation = evaluateSaturationState((silminState->ySol), (silminState->yLiq));

#ifndef BATCH_VERSION
        updateSolidADB((silminState->ySol), (silminState->yLiq));
        workProcData->active = TRUE;
#endif
        curStep++;
        return FALSE;
        /* ------------------------------------------------------------------------ */
    case UPDATE_SYSTEM:

        stateChange = FALSE;

        if (tState < 0) {
            stateChange = TRUE;
            tState    = hasSupersaturation;
            tInterval = (hasSupersaturation) ? 50.0 : -50.0;
#ifndef BATCH_VERSION
            workProcData->mode = TRUE;
            tpValues[TP_PADB_INDEX_T_INITIAL].value += tInterval;
#else
            silminState->T += tInterval;
            if (silminState->T > 2500.0) { meltsStatus.status = LIQUIDUS_MAX_T; return TRUE; }
            if (silminState->T <  500.0) { meltsStatus.status = LIQUIDUS_MIN_T; return TRUE; }
#endif
        } else {
            if (tState != hasSupersaturation) {
                tState = hasSupersaturation; tInterval *= -0.5;
            }
            if (fabs(tInterval) > 0.1) {
                stateChange = TRUE;
#ifndef BATCH_VERSION
                workProcData->mode = TRUE;
                tpValues[TP_PADB_INDEX_T_INITIAL].value += tInterval;
#else
                silminState->T += tInterval;
#endif
            } else {
#ifndef BATCH_VERSION
                wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "<> Found the liquidus at T = %.2f (C).\n", silminState->T-273.15);
#else
                printf("<> Found the liquidus at T = %.2f (C).\n", silminState->T-273.15);
                meltsStatus.status = LIQUIDUS_SUCCESS;
#endif
                tState = -1;
            }
        }

#ifndef BATCH_VERSION
        workProcData->active = stateChange;
#endif
        curStep = 0;
        return (!stateChange);
        /* ======================================================================== */
    } /* end switch */
    return FALSE;
}

#ifdef PHMELTS_ADJUSTMENTS
/* Adapted from rMELTSframework.m */
int findWetLiquidus(void) {
    int *oldIncSolids = (int *) malloc((size_t) npc*sizeof(int));
    double oldT = 0.0;
    int i, j, iter = 0, failure = TRUE;

    if ((silminState->solidMass > 0.0) || (silminState->nLiquidCoexist > 1)) {
        meltsStatus.status = LIQUIDUS_MULTIPLE;
        printf("<><><><> !!! <><><><> Unable to locate (Wet) liquidus with solids/multiple liquids present.\n");
    }
    else {
        do {
            oldT = silminState->T;
            for (i=0, j=0; i<npc; i++) if ((solids[i].type == PHASE) && (solids[i].nr == 0 || (solids[i].nr > 0 && solids[i].convert != NULL))) {
                    oldIncSolids[j] = (silminState->incSolids)[j];
                    (silminState->incSolids)[j] = FALSE;
                    if (!strcmp("fluid", solids[i].label)) (silminState->incSolids)[j] = TRUE;
                    j++;
                }
            while(!silmin());
            fprintf(stderr, "<><><><> !!! <><><><> T after equilibrate: %lf\n", silminState->T);
            if (meltsStatus.status != SILMIN_SUCCESS) {
                meltsStatus.status = LIQUIDUS_SILMIN_ERROR;
                printf("<><><><> !!! <><><><>  Failure in silmin() stage of findWetLiquidus. Aborting ...\n");
                break;
            }

            for (i=0, j=0; i<npc; i++) if ((solids[i].type == PHASE) && (solids[i].nr == 0 || (solids[i].nr > 0 && solids[i].convert != NULL))) {
                    (silminState->incSolids)[j] = oldIncSolids[j];
                    j++;
                }
            silminState->solidMass = 0.0; // added to make sure liquidus() doesn't fail
            while(!liquidus());
            fprintf(stderr, "<><><><> !!! <><><><> T after liquidus: %lf, iter: %d\n", silminState->T, iter);
            if (meltsStatus.status != LIQUIDUS_SUCCESS) {
                printf("<><><><> !!! <><><><>  Failure in liquidus() stage of findWetLiquidus. Aborting ...\n");
                break;
            }

#ifndef ALPHAMELTS_UPDATE_SYSTEM
            /* check what should happen if coming from MATLAB/Python */
            silminState->dspTstart = silminState->T;
            silminState->dspTstop  = silminState->T;
#endif
            iter++;
        } while ((fabs(oldT - silminState->T) > 0.5) && (iter < 50));

        if (iter == 50) {
            meltsStatus.status = GENERIC_INTERNAL_ERROR;
            printf("<><><><> !!! <><><><> Unable to locate (Wet) liquidus.\n");
        } else {
            printf("<><><><> !!! <><><><> Found (Wet) liquidus at %.2f (C).\n", silminState->T-273.15);
            failure = FALSE;
        }
    }
    free(oldIncSolids);

    return failure;

}
#endif

/* end of file LIQUIDUS.C */
