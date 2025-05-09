const char *interface_ver(void) { return "$Id: interface.c,v 1.13 2009/05/14 04:23:59 ghiorso Exp $"; }
/*
    MELTS Source Code: RCS $Log: interface.c,v $
    MELTS Source Code: RCS Revision 1.9  2008/03/06 17:51:23  ghiorso
    MELTS Source Code: RCS New fluid fractionation mode and other enhancements.
    MELTS Source Code: RCS
    MELTS Source Code: RCS Revision 1.8  2007/12/22 22:43:30  ghiorso
    MELTS Source Code: RCS Fixed error in BM integration in Gibbs.c
    MELTS Source Code: RCS Updated param_struct_data.h file for AGU 2007 xMELTS parameters
    MELTS Source Code: RCS Added support for status file production in MELTS-batch (XML)
    MELTS Source Code: RCS
    MELTS Source Code: RCS Revision 1.7  2007/10/29 19:58:02  ghiorso
    MELTS Source Code: RCS Updated ss volume terms for new regression.
    MELTS Source Code: RCS Miscellaneous fixes.
    MELTS Source Code: RCS
    MELTS Source Code: RCS Revision 1.6  2007/08/23 16:09:36  ghiorso
    MELTS Source Code: RCS Database updates.
    MELTS Source Code: RCS
    MELTS Source Code: RCS Revision 1.5  2007/06/14 16:41:58  ghiorso
    MELTS Source Code: RCS Corrected error in Batch MELTS (xml, Tinc)
    MELTS Source Code: RCS
    MELTS Source Code: RCS Revision 1.4  2007/06/08 17:25:42  ghiorso
    MELTS Source Code: RCS Added code to allow regression of Ghiorso EOS parameters
    MELTS Source Code: RCS
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
    MELTS Source Code: RCS Revision 1.3  2005/01/24 03:38:04  cvsaccount
    MELTS Source Code: RCS
    MELTS Source Code: RCS Added new files and modifications to perform builds for MgO-SiO2 system
    MELTS Source Code: RCS
    MELTS Source Code: RCS Revision 1.2  2004/12/04 19:10:36  cvsaccount
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
    * Revision 3.6  1997/06/21  22:49:51  ghiorso
    * June 1997 MELTS 3.0.x release
    * (prior to new entropy and regression model being introduced)
    *
    * Revision 3.5  1997/05/03  20:23:30  ghiorso
    * *** empty log message ***
    *
    * Revision 3.4  1997/03/27  17:03:33  ghiorso
    * *** empty log message ***
    *
    * Revision 3.3  1996/09/24  20:33:38  ghiorso
    * Version modified for OSF/1 4.0
    *
    * Revision 3.2  1995/12/09  19:26:38  ghiorso
    * Interface revisions for status box and graphics display
    *
    * Revision 3.1  1995/08/18  17:41:13  ghiorso
    * MELTS Version 3 - Initial Entry
    *
    * Revision 3.1  1995/08/18  17:41:13  ghiorso
    * MELTS Version 3 - Initial Entry
    *
    */

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      User Interface main program (file: INTERFACE.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  August 20, 1990   Original Version
**              Canabalized from TEST.C (graph widget test facility)
**              Testing of basic Widget structures
**      V1.0-2  Mark S. Ghiorso  August 23, 1990
**              First modularization
**      V1.0-3  Mark S. Ghiorso  August 24, 1990
**              Second modularization
**      V1.0-4  Mark S. Ghiorso  September 8, 1990
**              Stabilitized logic for initial conditions and looping
**      V1.0-5  Mark S. Ghiorso  August 30, 1991
**              Added global definitions of extern widgets
**      V1.0-6  Mark S. Ghiorso  October 18, 1991
**              Added management of silmin_padb widget
**      V2.0-1  Mark S. Ghiorso  November 14, 1991
**              Conversion to OSF Motif V1.1.1/X11 Release 4
**      V2.0-2  Mark S. Ghiorso  November 23, 1991
**              Removed silmin_padb widget
**      V2.0-3  Mark S. Ghiorso  November 29, 1991
**              (1) Added icon pixmap, icon name and toLevel shell
**                  title. The latter defines the application release
**      V2.0-4  Mark S. Ghiorso  February 18, 1992
**              (1) Minor changes for ANSI compliance
**      V2.0-5  Mark S. Ghiorso  March 17, 1992
**              Version 1.0 - alpha release
**      V2.0-6  Mark S. Ghiorso  June 14, 1993
**              Version 1.1.1 - beta release
**      V2.0-7  Mark S. Ghiorso  September 21, 1993
**              Corrected fifth argument to XtVaAppInitialize
**              Version 1.1.2 - beta release
**      V2.0-8  Mark S. Ghiorso  October 5, 1993
**              Version 1.1.3 - beta release
**      V2.0-9  Mark S. Ghiorso  March 16, 1993
**              Added cast to (Cardinal *) for R4 OpenVMS alpha
**              platform
**      V3.0-1  Mark S. Ghiorso  May 11, 1994
**              (1) Added banner for Version 2.0.1 alpha release
**      V3.0-2  Mark S. Ghiorso  February, 18, 1995
**              (1) Added banner for Version 2.0.1 beta release
**              (2) Added demo banner
**--
*/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "interface.h"

/*
 *=============================================================================
 * Initialize global SILMIN structures
 */

#include "silmin.h"

/* Defined in silmin.c instead
   #ifdef BATCH_VERSION
   #include "status.h"
   MeltsStatus meltsStatus;
   #endif

   #include "liq_struct_data.h"
   #include "sol_struct_data.h"

   int calculationMode = MODE_DEFAULT;
   int quad_tol_modifier = 1;

   void (*additionalOutput) (char *filename) = NULL;
   char *addOutputFileName = NULL;
*/

// NEED TO FIX UP BELOW
#ifdef ALPHAMELTS_UPDATE_SYSTEM
#define READ_ERROR printf("Error during input file read, offending record: %s", line); return GET_INPUT_ERROR_BAD_READ;
#elif defined(PHMELTS_ADJUSTMENTS)
#define READ_ERROR printf("Error during settings read, offending record: %s", line); return GET_INPUT_ERROR_BAD_READ;
#else
#define READ_ERROR printf("Error during input file read, offending record: %s", line); return GET_INPUT_ERROR_BAD_READ;
#endif

#ifndef BATCH_VERSION
/*
    static void getCommandLineAndEnvironment (int argc, char *argv[])
*/
#else
/* BATCH_VERSION */

extern SilminState *previousSilminState;

#define REALLOC(x, y) (((x) == NULL) ? malloc(y) : realloc((x), (y)))
#define REC 134

/* returns TRUE if successful */
/*         FALSE if not       */

#ifdef PHMELTS_ADJUSTMENTS
int getInputDataFromFile(char *fileName)
#else
    int batchInputDataFromFile(char *fileName)
#endif
{
    static FILE *input = NULL;
    static char line[REC];
    static double *oxideWt = NULL;
#ifdef PHMELTS_ADJUSTMENTS
    static int npOld = -1, nsOld = 0; /* liquid or bulk */
#endif
    int i, j, np, ns;
    double sum;

    /* -> Check validity of file name */

    if ((input = fopen (fileName, "r")) == NULL) {
        printf("Error in SILMIN file input procedure. Cannot open file: %s\n", fileName);
        return GET_INPUT_ERROR_BAD_FILE;
    }

    /* -> Initial global and static storage (label moved to getInputDataFromLine, silminInputData added from read_write.c) */

    if (silminInputData.name == NULL)  silminInputData.name  = (char *) malloc((unsigned) (REC+1)*sizeof(char));
    if (silminInputData.title == NULL) silminInputData.title = (char *) malloc((unsigned) (REC+1)*sizeof(char));

    if (oxideWt == NULL) oxideWt = (double *) malloc((unsigned) nc*sizeof(double));

    /* -> Set default values */

    /* from read_write.c (should strn)*/
    (void) strcpy(silminInputData.name, fileName);
    (void) strcpy(silminInputData.title, "");

    sum = 0.0;
    for (i=0; i<nc; i++) oxideWt[i] = 0.0;

#ifdef PHMELTS_ADJUSTMENTS
    for (i=0, np=0; i<npc; i++) if (solids[i].type == PHASE) {
        if (solids[i].inStdSet)
            (silminState->incSolids)[np] = solids[i].na;
        else
            (silminState->incSolids)[np] = FALSE;
        np++;
    }
    //(silminState->incSolids)[npc] = nlc;
    silminState->incLiquids = 1;
    silminState->fracOut = 1.0;
    silminState->maxF = 1.0;

    np = -1; /* liquid or bulk */
    ns = 0; /* number of assimilant */
#else
    for (i=0, np=0; i<npc; i++) if (solids[i].type == PHASE) { (silminState->incSolids)[np] = TRUE; np++; }
    (silminState->incSolids)[npc] = TRUE;
#endif

    silminState->nLiquidCoexist  = 1;
    silminState->fo2Path  = FO2_NONE;

    /* -> Read input data file */

    for (;;) {
        /* -> title record */
        if (fgets(line, REC, input) == NULL) break;
        if (getInputDataFromLine(line, &np, &ns, oxideWt) == GET_INPUT_ERROR_BAD_READ) {
            fclose(input);
            return GET_INPUT_ERROR_BAD_READ;
        }
        else if (np != npOld || (npOld >= 0 && ns > (nsOld+1))) {
            if (npOld < 0) {
                for (i=0; i<nc; i++) sum += oxideWt[i];
                if (sum != 0.0) {

                    for (i=0, silminState->liquidMass=0.0; i<nc; i++) {
                        silminState->liquidMass    += oxideWt[i];
                        (silminState->bulkComp)[i]  = oxideWt[i]/bulkSystem[i].mw;
                    }

                    for (i=0; i<nlc; i++) {
                        for ((silminState->liquidComp)[0][i]=0.0, silminState->oxygen=0.0, j=0; j<nc; j++) {
                            (silminState->liquidComp)[0][i] += (silminState->bulkComp)[j]*(bulkSystem[j].oxToLiq)[i];
                            silminState->oxygen += (silminState->bulkComp)[j]*(bulkSystem[j].oxToLiq)[i]*(oxygen.liqToOx)[i];
                        }
                    }
                }
            }
            else {
                /* should already be (re)allocated */
                if ((silminState->nDspAssimComp)[npOld] == 0) (silminState->dspAssimComp)[npOld] = (double *) malloc((unsigned) sizeof(double));
                else  (silminState->dspAssimComp)[npOld] = (double *) REALLOC((silminState->dspAssimComp)[npOld], (unsigned) (nsOld+1)*sizeof(double));

                if (solids[npOld].na != 1) { /* ignore composition if pure phase */
                    double *e = calloc((size_t) 106, sizeof(double)), *m = calloc((size_t) solids[npOld].na, sizeof(double));
                    for (i=0; i<nc; i++) {
                        double mOx = oxideWt[i]/bulkSystem[i].mw;
                        for (j=0; j<106; j++) e[j] += mOx*(bulkSystem[i].oxToElm)[j];
                    }

                    if (!strncmp(solids[npOld].label,"clinopyroxene", MIN((int) strlen(solids[npOld].label), 13)) ||
                            !strncmp(solids[npOld].label,"orthopyroxene", MIN((int) strlen(solids[npOld].label), 13))) {
                        // use input Fe2O3 for Opx/Cpx as entered (see marc.c)
                        (*solids[npOld].convert)(FIRST, FIRST, silminState->assimT+273.15, silminState->P, e, m, NULL, NULL, NULL, NULL, NULL, NULL);
                    }
                    else {
                        (*solids[npOld].convert)(FIRST, SECOND, silminState->assimT+273.15, silminState->P, e, m, NULL, NULL, NULL, NULL, NULL, NULL);
                    }

                    for (j=0; j<solids[npOld].na; j++) {
                        if ((silminState->nDspAssimComp)[npOld+1+j] == 0) (silminState->dspAssimComp)[npOld+1+j] = (double *) malloc((unsigned) sizeof(double));
                        else  (silminState->dspAssimComp)[npOld+1+j] = (double *) REALLOC((silminState->dspAssimComp)[npOld+1+j], (unsigned) (nsOld+1)*sizeof(double));
                        (silminState->dspAssimComp)[npOld+1+j][nsOld] *= 100.0*m[j];
                        (silminState->nDspAssimComp)[npOld+1+j] = nsOld+1;
                    }
                    free(m);
                    free(e);

                }
                nsOld = ns;
            }
            npOld = np;
            sum = 0.0;
            for (i=0; i<nc; i++) oxideWt[i] = 0.0;
        }

    }

    /* -> Close and discard file and return */
    fclose(input);

#ifdef PHMELTS_ADJUSTMENTS
    for (i=0; i<nc-3; i++) {  /* SO3, Cl2O-1, F2O-1 */
        if (calculationMode != MODE_xMELTS) { /* i.e. calculationMode has been set */
            if ((calculationMode == MODE_pMELTS) &&  /* not calibrated */
	            (!strncmp(bulkSystem[i].label, "MnO", 3) || !strncmp(bulkSystem[i].label, "NiO", 3) ||
                !strncmp(bulkSystem[i].label, "CoO", 3) || !strncmp(bulkSystem[i].label, "CO2", 3))) {
                    if (oxideWt[i] != 0.0) {
                        printf("Oxide is not calibrated for the pMELTS model, zeroing composition: %s\n", bulkSystem[i].label); \
                        fprintf(stderr, "Ignored line... Initial Composition: %s %f\n", bulkSystem[i].label, oxideWt[i]);
                        oxideWt[i] = 0.0;
                    }
                }
            else if ((calculationMode == MODE__MELTS) && !strncmp(bulkSystem[i].label, "CO2", 3)) {
                printf("Oxide is not calibrated for the Rhyolite-MELTS 1.0.2 model, zeroing composition: %s\n", bulkSystem[i].label); \
                fprintf(stderr, "Ignored line... Initial Composition: %s %f\n", bulkSystem[i].label, oxideWt[i]);
                oxideWt[i] = 0.0;
            }
        }
    }
    for (i=nc-3; i<nc; i++) {  /* SO3, Cl2O-1, F2O-1 */
        if (oxideWt[i] != 0.0) {
            printf("Oxide is not calibrated for the chosen MELTS model, zeroing composition: %s\n", bulkSystem[i].label); \
            fprintf(stderr, "Ignored line... Initial Composition: %s %f\n", bulkSystem[i].label, oxideWt[i]);
            oxideWt[i] = 0.0;
        }
    }
#endif

    if (npOld < 0) {
        for (i=0; i<nc; i++) sum += oxideWt[i];
        if (sum != 0.0) { /* sum will be zero if bulkComp has already been set */

            for (i=0, silminState->liquidMass=0.0; i<nc; i++) {
                silminState->liquidMass    += oxideWt[i];
                (silminState->bulkComp)[i]  = oxideWt[i]/bulkSystem[i].mw;
            }

            for (i=0, silminState->oxygen=0.0; i<nlc; i++) {
                for ((silminState->liquidComp)[0][i]=0.0, j=0; j<nc; j++) {
                    (silminState->liquidComp)[0][i] += (silminState->bulkComp)[j]*(bulkSystem[j].oxToLiq)[i];
                    silminState->oxygen += (silminState->bulkComp)[j]*(bulkSystem[j].oxToLiq)[i]*(oxygen.liqToOx)[i];
                }
            }
        }
    } else {
        /* should already be (re)allocated */
        if ((silminState->nDspAssimComp)[npOld] == 0) (silminState->dspAssimComp)[npOld] = (double *) malloc((unsigned) sizeof(double));
        else  (silminState->dspAssimComp)[npOld] = (double *) REALLOC((silminState->dspAssimComp)[npOld], (unsigned) (nsOld+1)*sizeof(double));

        if (solids[npOld].na != 1) { /* ignore composition if pure phase */
            double *e = calloc((size_t) 106, sizeof(double)), *m = calloc((size_t) solids[npOld].na, sizeof(double));
            for (i=0; i<nc; i++) {
                double mOx = oxideWt[i]/bulkSystem[i].mw;
                for (j=0; j<106; j++) e[j] += mOx*(bulkSystem[i].oxToElm)[j];
            }

            if (!strncmp(solids[npOld].label,"clinopyroxene", MIN((int) strlen(solids[npOld].label), 13)) ||
                    !strncmp(solids[npOld].label,"orthopyroxene", MIN((int) strlen(solids[npOld].label), 13))) {
                // use input Fe2O3 for Opx/Cpx as entered (see marc.c)
                (*solids[npOld].convert)(FIRST, FIRST, silminState->assimT+273.15, silminState->P, e, m, NULL, NULL, NULL, NULL, NULL, NULL);
            }
            else {
                (*solids[npOld].convert)(FIRST, SECOND, silminState->assimT+273.15, silminState->P, e, m, NULL, NULL, NULL, NULL, NULL, NULL);
            }

            for (j=0; j<solids[npOld].na; j++) {
                if ((silminState->nDspAssimComp)[npOld+1+j] == 0) (silminState->dspAssimComp)[npOld+1+j] = (double *) malloc((unsigned) sizeof(double));
                else  (silminState->dspAssimComp)[npOld+1+j] = (double *) REALLOC((silminState->dspAssimComp)[npOld+1+j], (unsigned) (nsOld+1)*sizeof(double));
                (silminState->dspAssimComp)[npOld+1+j][nsOld] *= 100.0*m[j];
                (silminState->nDspAssimComp)[npOld+1+j] = nsOld+1;
            }
            free(m);
            free(e);

        }

    }

    if (silminState->assimilate) {
        if (silminState->nAssimComp == NULL) silminState->nAssimComp = (int *)     calloc((unsigned) (npc+nc), sizeof(int));
        if (silminState->assimComp  == NULL) silminState->assimComp  = (double **) calloc((unsigned) (npc+nc), sizeof(double *));

        { /* adapted from silmin_support - this should be moved to its own function */

            if      (silminState->dspAssimUnits == ASSIM_PADB_UNITS_WEIGHT) {
                double *m = (double *) malloc((unsigned) nc*sizeof(double));
                for (i=0; i<npc; i++) if (solids[i].type == PHASE) {
                        for (ns=0; ns<silminState->nDspAssimComp[i]; ns++) {
                            if ((silminState->nAssimComp)[i] == 0) (silminState->assimComp)[i] = (double *) malloc((unsigned) sizeof(double));
                            else  (silminState->assimComp)[i] = (double *) REALLOC((silminState->assimComp)[i], (unsigned) (ns+1)*sizeof(double));
                            if (solids[i].na == 1) {
                                (silminState->assimComp)[i][ns] = (silminState->dspAssimComp)[i][ns]*(silminState->dspAssimMass)/(100.0*solids[i].mw);
                            } else {
                                double mw = 0.0;
                                for (j=0; j<solids[i].na; j++) {
                                    if ((silminState->nAssimComp)[i] == 0) (silminState->assimComp)[i+1+j] = (double *) malloc((unsigned) sizeof(double));
                                    else  (silminState->assimComp)[i+1+j] = (double *) REALLOC((silminState->assimComp)[i+1+j], (unsigned) (ns+1)*sizeof(double));
                                }
                                for (j=0; j<solids[i].na; j++) mw += (silminState->dspAssimComp)[i+1+j][ns]*solids[i+1+j].mw/100.0;
                                (silminState->assimComp)[i][ns] = (mw == 0.0) ? 0.0 : (silminState->dspAssimComp)[i][ns]*(silminState->dspAssimMass)/(100.0*mw);
                                for (j=0; j<solids[i].na; j++) (silminState->assimComp)[i+1+j][ns] = (silminState->dspAssimComp)[i+1+j][ns]*(silminState->assimComp)[i][ns]/100.0;
                            }
                            silminState->nAssimComp[i] = ns+1;
                        }
                    }
                for (ns=0; ns<silminState->nDspAssimComp[1+npc]; ns++) {
                    if ((silminState->nAssimComp)[i] == 0) (silminState->assimComp)[i] = (double *) malloc((unsigned) sizeof(double));
                    else  (silminState->assimComp)[i] = (double *) REALLOC((silminState->assimComp)[i], (unsigned) (ns+1)*sizeof(double));
                    for (i=0; i<nc; i++) m[i] = (silminState->dspAssimComp)[i+npc][ns]*(silminState->dspAssimLiqM)*(silminState->dspAssimMass)/(100.0*100.0*bulkSystem[i].mw);
                    for (i=0; i<nlc; i++)
                        for (j=0, (silminState->assimComp)[i+npc][ns]=0.0; j<nc; j++) (silminState->assimComp)[i+npc][ns] += m[j]*(bulkSystem[j].oxToLiq)[i];
                    silminState->nAssimComp[1+npc] = ns+1;
                }
                free(m);
            } else if (silminState->dspAssimUnits == ASSIM_PADB_UNITS_VOLUME) {
                printf("Code not yet implemented to convert a vol %% assimilant !");
            }

        }


        { /* from silmin_support - this should be moved to its own function */
            double *m = (double *) malloc((unsigned)     nlc*sizeof(double));
            double *r = (double *) malloc((unsigned) (nlc-1)*sizeof(double));
            double enthalpy, entropy, volume, totalEnthalpy, totalEntropy,
                totalVolume;

            totalEnthalpy = 0.0;
            totalEntropy  = 0.0;
            totalVolume   = 0.0;

            for (j=0; j<npc; j++) {
                if (solids[j].type == PHASE) {
                    for (ns=0; ns<(silminState->nAssimComp)[j]; ns++) {
                        if (solids[j].na == 1) {
                            gibbs(silminState->assimT+273.15, silminState->P, (char *) solids[j].label, &(solids[j].ref), NULL, NULL, &(silminState->assimTD));
                            enthalpy = (silminState->assimComp)[j][ns]*(silminState->assimTD).h;
                            entropy  = (silminState->assimComp)[j][ns]*(silminState->assimTD).s;
                            volume   = (silminState->assimComp)[j][ns]*(silminState->assimTD).v;
                        } else {
                            for (i=0; i<solids[j].na; i++) m[i] = (silminState->assimComp)[j+1+i][ns];
                            (*solids[j].convert)(SECOND, THIRD, silminState->assimT+273.15, silminState->P, NULL, m, r, NULL, NULL, NULL, NULL, NULL);
                            (*solids[j].hmix)(FIRST, silminState->assimT+273.15, silminState->P, r, &enthalpy);
                            (*solids[j].smix)(FIRST, silminState->assimT+273.15, silminState->P, r, &entropy, NULL, NULL);
                            (*solids[j].vmix)(FIRST, silminState->assimT+273.15, silminState->P, r, &volume, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
                            enthalpy     *= (silminState->assimComp)[j][ns];
                            entropy      *= (silminState->assimComp)[j][ns];
                            volume       *= (silminState->assimComp)[j][ns];
                            for (i=0; i<solids[j].na; i++) {
                                gibbs(silminState->assimT+273.15, silminState->P, (char *) solids[j+1+i].label, &(solids[j+1+i].ref), NULL, NULL, &(silminState->assimTD));
                                enthalpy += m[i]*(silminState->assimTD).h;
                                entropy  += m[i]*(silminState->assimTD).s;
                                volume   += m[i]*(silminState->assimTD).v;
                            }
                        }
                        totalEnthalpy += enthalpy;
                        totalEntropy  += entropy;
                        totalVolume   += volume;
                    }
                }
            }

            for (ns=0; ns<(silminState->nAssimComp)[npc+1]; ns++) {
                double totalMoles;
                for (i=0, totalMoles=0.0; i<nlc; i++) {
                    m[i] = (silminState->assimComp)[npc+i][ns];
                    totalMoles += m[i];
                }
                conLiq(SECOND, THIRD, silminState->assimT+273.15, silminState->P, NULL, m, r, NULL, NULL, NULL, NULL);
                hmixLiq(FIRST, silminState->assimT+273.15, silminState->P, r, &enthalpy, NULL);
                smixLiq(FIRST, silminState->assimT+273.15, silminState->P, r, &entropy,  NULL, NULL, NULL);
                vmixLiq(FIRST, silminState->assimT+273.15, silminState->P, r, &volume,   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
                enthalpy *= totalMoles;
                entropy  *= totalMoles;
                volume   *= totalMoles;
                for (i=0; i<nlc; i++) {
                    gibbs(silminState->assimT+273.15, silminState->P, (char *) liquid[i].label, &(liquid[i].ref), &(liquid[i].liq), &(liquid[i].fus), &(silminState->assimTD));
                    enthalpy += m[i]*(silminState->assimTD).h;
                    entropy  += m[i]*(silminState->assimTD).s;
                    volume   += m[i]*(silminState->assimTD).v;
                }
                totalEnthalpy += enthalpy;
                totalEntropy  += entropy;
                totalVolume   += volume;
            }

            (silminState->assimTD).g       = 0.0;
            (silminState->assimTD).h       = totalEnthalpy;
            (silminState->assimTD).s       = totalEntropy;
            (silminState->assimTD).v       = totalVolume;
            (silminState->assimTD).cp      = 0.0;
            (silminState->assimTD).dcpdt   = 0.0;
            (silminState->assimTD).dvdt    = 0.0;
            (silminState->assimTD).dvdp    = 0.0;
            (silminState->assimTD).d2vdt2  = 0.0;
            (silminState->assimTD).d2vdp2  = 0.0;
            (silminState->assimTD).d2vdtdp = 0.0;

            free(m);
            free(r);

        }

    }

    fflush(stderr);
    fflush(stdout);
    printf("Input file read. Waiting for command or user input.\n");
    return GET_INPUT_SUCCESS;
}

int getInputDataFromLine(char *line, int *npOld, int *nsOld, double *oxideWt) {

    static char *label;

    size_t len;
    int i, j, np, ns;
    float temporary;

    /* -> Initial global and static storage */

    if (label == NULL) {
        for (i=0, len=0; i<nc; i++) len = MAX(len, (int) strlen(bulkSystem[i].label));
        label = (char *) malloc((unsigned) (len+1)*sizeof(char));
    }

    len = strlen(line); for (i=0; i<len; i++) line[i] = tolower(line[i]);

    if        (!strncmp(line, "title: ",                 MIN(len, 7))) {

#ifdef ALPHAMELTS_UPDATE_SYSTEM
            if (sscanf(&line[6], "%s", silminInputData.title) == EOF) { READ_ERROR }
#endif
        /* -> initial composition record */
    } else if (!strncmp(line, "initial composition: ",   MIN(len,21))) {
#ifdef PHMELTS_ADJUSTMENTS
        if (npOld == NULL || nsOld == NULL || oxideWt == NULL) {READ_ERROR}
        *npOld = -1;
#endif
        for (i=0; i<nc; i++) {
            for (j=0; j < ((int)strlen(bulkSystem[i].label)); j++) label[j] = tolower((bulkSystem[i].label)[j]);
            label[j] = '\0';
            if (!strncmp(&line[21], label, MIN((len-21), (int) strlen(label)))) {
                if (sscanf(&line[21 + (int) strlen(label)], "%f", &temporary) == EOF) { READ_ERROR }
                oxideWt[i]  = (double) temporary;
                //                    sum        += (double) temporary;
                break;
            }
        }
        if (i == nc) { READ_ERROR }

        /* -> initial temperature record */
    } else if (!strncmp(line, "initial temperature: ",   MIN(len,21))) {
        if (sscanf(&line[21], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->T         = (double) temporary + 273.15;
#ifdef ALPHAMELTS_UPDATE_SYSTEM
        silminState->dspTstart = (double) temporary; /* could be changed */
#elif !defined(PHMELTS_ADJUSTMENTS)
        silminState->dspTstart = (double) temporary + 273.15;
#endif
        /* -> final temperature record */
    } else if (!strncmp(line, "final temperature: ",     MIN(len,19))) {
        if (sscanf(&line[19], "%f", &temporary) == EOF) { READ_ERROR }
#ifdef ALPHAMELTS_UPDATE_SYSTEM
        silminState->dspTstop = (double) temporary; /* could be changed */
#elif !defined(PHMELTS_ADJUSTMENTS)
        silminState->dspTstop = (double) temporary + 273.15;
#endif
        /* -> increment temperature record */
    } else if (!strncmp(line, "increment temperature: ", MIN(len,23))) {
        if (sscanf(&line[23], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->dspTinc = (double) temporary;

        /* -> initial pressure record */
    } else if (!strncmp(line, "initial pressure: ",      MIN(len,18))) {
        if (sscanf(&line[18], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->P         = (double) temporary;
        silminState->dspPstart = (double) temporary;

        /* -> final pressure record */
    } else if (!strncmp(line, "final pressure: ",        MIN(len,16))) {
        if (sscanf(&line[16], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->dspPstop = (double) temporary;

#ifdef PHMELTS_ADJUSTMENTS
        /* -> min f record for silicate liquid(s) only */
    } else if (!strncmp(line, "set min f: ",        MIN(len,11))) {
        if (sscanf(&line[11], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->minF = (double) temporary;
        /* -> min phi record for fluid or liquid */
    } else if (!strncmp(line, "set min phi: fluid",        MIN(len,13))) {
        if (sscanf(&line[19], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->minFluPhi = (double) temporary;
    } else if (!strncmp(line, "set min phi: water",        MIN(len,13))) {
        if (sscanf(&line[19], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->minFluPhi = (double) temporary;
    } else if (!strncmp(line, "set min phi: liquid",        MIN(len,13))) {
        if (sscanf(&line[20], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->minLiqPhi = (double) temporary;
    } else if (!strncmp(line, "set min phi: ",        MIN(len,13))) {
        if (sscanf(&line[13], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->minLiqPhi = (double) temporary;
        /* -> continuous ratio record for all phases */
    } else if (!strncmp(line, "set fractionation coefficient: ",        MIN(len,15))) {
        if (sscanf(&line[22], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->fracOut = (double) temporary;
        /* -> min f record for silicate liquid(s) only */
    } else if (!strncmp(line, "set max f: ",        MIN(len,11))) {
        if (sscanf(&line[11], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->maxF  = (double) temporary;

        /* -> initial mass scale record */
    } else if        (!strncmp(line, "initial mass: ",          MIN(len, 14))) {
        if (sscanf(&line[14], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->refMass   = (double) temporary;

#endif


        /* -> increment pressure record */
    } else if (!strncmp(line, "increment pressure: ",    MIN(len,20))) {
        if (sscanf(&line[20], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->dspPinc = (double) temporary;

        /* -> dp/dt record */
    } else if (!strncmp(line, "dp/dt: ",                 MIN(len, 7))) {
        if (sscanf(&line[7], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->dspDPDt = (double) temporary;

#ifdef ALPHAMELTS_UPDATE_SYSTEM
        /* -> initial enthalpy record */
    } else if (!strncmp(line, "initial enthalpy: ",      MIN(len,18))) {
        if (sscanf(&line[18], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->refEnthalpy = (double) temporary;

        /* -> final enthalpy record */
    } else if (!strncmp(line, "final enthalpy: ",        MIN(len,16))) {
        if (sscanf(&line[16], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->dspHstop = (double) temporary;
#endif
        /* -> increment enthalpy record */
    } else if (!strncmp(line, "increment enthalpy: ", MIN(len,20))) {
        if (sscanf(&line[20], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->dspHinc = (double) temporary;

#ifdef ALPHAMELTS_UPDATE_SYSTEM
        /* -> initial entropy record */
    } else if (!strncmp(line, "initial entropy: ",      MIN(len,17))) {
        if (sscanf(&line[17], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->refEntropy = (double) temporary;

        /* -> final entropy record */
    } else if (!strncmp(line, "final entropy: ",        MIN(len,15))) {
        if (sscanf(&line[15], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->dspSstop = (double) temporary;
#endif
        /* -> increment entropy record */
    } else if (!strncmp(line, "increment entropy: ", MIN(len,19))) {
        if (sscanf(&line[19], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->dspSinc = (double) temporary;

#ifdef ALPHAMELTS_UPDATE_SYSTEM
        /* -> initial volume record */
    } else if (!strncmp(line, "initial volume: ",      MIN(len,16))) {
        if (sscanf(&line[16], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->refVolume = (double) temporary;

        /* -> final volume record */
    } else if (!strncmp(line, "final volume: ",        MIN(len,14))) {
        if (sscanf(&line[14], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->dspVstop = (double) temporary;
#endif
        /* -> increment volume record */
    } else if (!strncmp(line, "increment volume: ", MIN(len,18))) {
        if (sscanf(&line[18], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->dspVinc = (double) temporary;

        /* -> dp/dH record */
    } else if (!strncmp(line, "dp/dh: ",                 MIN(len, 7))) {
        if (sscanf(&line[7], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->dspDPDH = (double) temporary;

        /* -> dp/dS record */
    } else if (!strncmp(line, "dp/ds: ",                 MIN(len, 7))) {
        if (sscanf(&line[7], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->dspDPDS = (double) temporary;

        /* -> dt/dV record */
    } else if (!strncmp(line, "dt/dv: ",                 MIN(len, 7))) {
        if (sscanf(&line[7], "%f", &temporary) == EOF) { READ_ERROR }
#ifdef PHMELTS_ADJUSTMENTS
        silminState->dspDTDV = (double) temporary;
#else
        silminState->dspDVDt = (temporary != 0.0) ? (double) 1.0/temporary : 0.0;
#endif
        /* -> log fo2 path record */
    } else if (!strncmp(line, "log fo2 path: ",          MIN(len,14))) {
        if      (!strncmp(&line[14], "none",  MIN((len-14), 4))) { silminState->fo2Path  = FO2_NONE; silminState->fo2Delta =  0.0; }
        else if (!strncmp(&line[14], "fmq",   MIN((len-14), 3))) { silminState->fo2Path  = FO2_QFM; silminState->fo2Delta =  0.0; }
        else if (!strncmp(&line[14], "coh",   MIN((len-14), 3))) { silminState->fo2Path  = FO2_COH; silminState->fo2Delta =  0.0; }
        else if (!strncmp(&line[14], "nno",   MIN((len-14), 3))) { silminState->fo2Path  = FO2_NNO; silminState->fo2Delta =  0.0; }
        else if (!strncmp(&line[14], "iw",    MIN((len-14), 2))) { silminState->fo2Path  = FO2_IW; silminState->fo2Delta =  0.0; }
        else if (!strncmp(&line[14], "hm",    MIN((len-14), 2))) { silminState->fo2Path  = FO2_HM; silminState->fo2Delta =  0.0; }
        else if (!strncmp(&line[14], "+3fmq", MIN((len-14), 5))) { silminState->fo2Path  = FO2_QFM; silminState->fo2Delta =  3.0; }
        else if (!strncmp(&line[14], "+2fmq", MIN((len-14), 5))) { silminState->fo2Path  = FO2_QFM; silminState->fo2Delta =  2.0; }
        else if (!strncmp(&line[14], "+1fmq", MIN((len-14), 5))) { silminState->fo2Path  = FO2_QFM; silminState->fo2Delta =  1.0; }
        else if (!strncmp(&line[14], "-1fmq", MIN((len-14), 5))) { silminState->fo2Path  = FO2_QFM; silminState->fo2Delta = -1.0; }
        else if (!strncmp(&line[14], "-2fmq", MIN((len-14), 5))) { silminState->fo2Path  = FO2_QFM; silminState->fo2Delta = -2.0; }
        else if (!strncmp(&line[14], "-3fmq", MIN((len-14), 5))) { silminState->fo2Path  = FO2_QFM; silminState->fo2Delta = -3.0; }
        else if (!strncmp(&line[14], "-4fmq", MIN((len-14), 5))) { silminState->fo2Path  = FO2_QFM; silminState->fo2Delta = -4.0; }
        else if (!strncmp(&line[14], "-5fmq", MIN((len-14), 5))) { silminState->fo2Path  = FO2_QFM; silminState->fo2Delta = -5.0; }
        else if (!strncmp(&line[14], "-6fmq", MIN((len-14), 5))) { silminState->fo2Path  = FO2_QFM; silminState->fo2Delta = -6.0; }
        else if (!strncmp(&line[14], "-7fmq", MIN((len-14), 5))) { silminState->fo2Path  = FO2_QFM; silminState->fo2Delta = -7.0; }
        else if (!strncmp(&line[14], "-8fmq", MIN((len-14), 5))) { silminState->fo2Path  = FO2_QFM; silminState->fo2Delta = -8.0; }
        else if (!strncmp(&line[14], "-9fmq", MIN((len-14), 5))) { silminState->fo2Path  = FO2_QFM; silminState->fo2Delta = -9.0; }
        else if (!strncmp(&line[14], "+0.5fmq", MIN((len-14), 7))) { silminState->fo2Path  = FO2_QFM; silminState->fo2Delta = +0.5; }
        else if (!strncmp(&line[14], "+1.5fmq", MIN((len-14), 7))) { silminState->fo2Path  = FO2_QFM; silminState->fo2Delta = +1.5; }
#ifdef PHMELTS_ADJUSTMENTS
        else if (!strncmp(&line[14], "abs", MIN((len-14), 3))) { silminState->fo2Path  = FO2_ABS; }
#endif
        else { READ_ERROR }



        /* DELTA RECORD */
    } else if (!strncmp(line, "log fo2 offset:", MIN(len, 15))) {
        if (sscanf(&line[15], "%f", &temporary) == EOF) { READ_ERROR }
        silminState->fo2Delta = (double) temporary;


        // ADD LIQUID AND FIX RHM OXIDE ETC

        /* -> suppress a solid phase record */
    } else if (!strncmp(line, "suppress: ",              MIN(len,10))) {
        if (!strncmp(&line[10], "none", MIN(len-10, 4))) {
            for (i=0, np=0; i<npc; i++) if (solids[i].type == PHASE) { (silminState->incSolids)[np] = TRUE; np++; }
#ifdef PHMELTS_ADJUSTMENTS
            silminState->incLiquids = TRUE;
#else
            (silminState->incSolids)[npc] = TRUE;
#endif
        }
        else if (!strncmp(&line[10], "liquid", MIN(len-10, 6))) {
#ifdef PHMELTS_ADJUSTMENTS
            silminState->incLiquids = FALSE;
#else
            (silminState->incSolids)[npc] = FALSE;
#endif
        }
        else {
            int j;
            for (i=0, j=0; i<npc; i++) {
                if (solids[i].type == PHASE) {
                    int phaseStrLen = (int) strlen(solids[i].label);
                    // Was == but this causes problems on Windows.
                    if (((len-10-phaseStrLen-1)  >= 0) && !strncmp(&line[10], solids[i].label, phaseStrLen)) {
                        if ( solids[i].nr == 0 || (solids[i].nr > 0 && solids[i].convert != NULL)) {
                            silminState->incSolids[j] = FALSE;
                        }
                        break;
                    }
                    j++;
                }
            }
            if (i == npc) { READ_ERROR }
        }

#ifdef PHMELTS_ADJUSTMENTS
        /* -> suppress exsolution for a phase record */
    } else if (!strncmp(line, "limit number: ", MIN(len,14))) {
        if (!strncmp(&line[14], "liquid", MIN(len-14, 6))) {
            if (sscanf(&line[20], "%f", &temporary) != EOF) {
                silminState->incLiquids = (int) temporary;
                if (silminState->incLiquids <= 0) silminState->incLiquids = FALSE;
            }
        }
        else {
            int j;
            for (i=0, j=0; i<npc; i++) {
                if (solids[i].type == PHASE) {
                    int phaseStrLen = (int) strlen(solids[i].label);
                    if (!strncmp(&line[14], solids[i].label, MIN(len-14, phaseStrLen)) && (line[14+phaseStrLen] = ' ')) {
                        break;
                    }
                    j++;
                }
            }
            if (i == npc) { READ_ERROR }
            else if (sscanf(&line[14 + (int) strlen(solids[i].label)], "%f", &temporary) != EOF) {
                silminState->incSolids[j] = (int) temporary;
                if (silminState->incSolids[j] <= 0) silminState->incSolids[j] = FALSE;
            }
        }
#endif

        /* -> mode record */
    } else if (!strncmp(line, "mode: ",                   MIN(len, 6))) {
        int i;
        if      (!strncmp(&line[6],  "fractionate solids",  MIN((len-6), 18))) {
            silminState->fractionateSol = TRUE;
#ifdef PHMELTS_ADJUSTMENTS
            if (silminState->fracSolids == NULL) silminState->fracSolids = (double *) calloc((unsigned) npc, sizeof(double));
            for (i=0; i<npc; i++) silminState->fracSolids[i] = silminState->fracOut;
#endif
        }
        else if (!strncmp(&line[6],  "fractionate liquids", MIN((len-6), 19))) {
            silminState->fractionateLiq = TRUE;
#ifdef PHMELTS_ADJUSTMENTS
            if (silminState->fracLiquids == NULL) silminState->fracLiquids = (double *) calloc((unsigned) nlc, sizeof(double));
            for (i=0; i<nlc; i++) silminState->fracLiquids[i] = silminState->fracOut;
#endif
        }
        else if (!strncmp(&line[6],  "fractionate fluids",  MIN((len-6), 18))) {
            silminState->fractionateFlu = TRUE;
#ifdef PHMELTS_ADJUSTMENTS
            if (silminState->fracFluids == NULL) silminState->fracFluids = (double *) calloc((unsigned) 2, sizeof(double));
            for (i=0; i<2; i++) silminState->fracFluids[i] = silminState->fracOut;
#endif
        }
        else if (!strncmp(&line[6],  "fractionate none",  MIN((len-6), 16))) {
            silminState->fractionateSol = FALSE; silminState->fractionateLiq = FALSE; silminState->fractionateFlu = FALSE;
#ifdef PHMELTS_ADJUSTMENTS
            if (silminState->fracLiquids != NULL) for (i=0; i<nlc; i++) silminState->fracLiquids[i] = 0.0;
            if (silminState->fracFluids != NULL) for (i=0; i<2; i++) silminState->fracFluids[i] = 0.0;
            if (silminState->fracSolids != NULL) for (i=0; i<npc; i++) silminState->fracSolids[i] = 0.0;
#endif
        }
#ifdef PHMELTS_ADJUSTMENTS
        /* To catch typos (especially if multiple liquids and/or mixed fluid is off) and to choose optionally to choose liquid/fluid */
        else if (!strncmp(&line[6],  "fractionate liquid", MIN((len-6), 18))) {
            silminState->fractionateLiq = TRUE;
            if (silminState->fracLiquids == NULL) silminState->fracLiquids = (double *) calloc((unsigned) nlc, sizeof(double));
            if (sscanf(&line[24], "%d", &i) == EOF) i = 0; i = MIN(MAX(i, 1), nlc);
            silminState->fracLiquids[i-1] = silminState->fracOut;
        }
        else if (!strncmp(&line[6],  "fractionate fluid",  MIN((len-6), 17))) {
            silminState->fractionateFlu = TRUE;
            if (silminState->fracFluids == NULL) silminState->fracFluids = (double *) calloc((unsigned) 2, sizeof(double));
            if (sscanf(&line[23], "%d", &i) == EOF) i = 0; i = MIN(MAX(i, 1), 2);
            silminState->fracFluids[i-1] = silminState->fracOut;
        }
        /* -> fractionate a phase record */
        else if (!strncmp(&line[6], "fractionate ", MIN((len-6),12))) {
            // EVENTUALLY PACK THIS THE SAME AS INCSOLIDS? */
            if (silminState->fracSolids == NULL) silminState->fracSolids = (double *) calloc((unsigned) npc, sizeof(double));
            if (silminState->fracFluids == NULL) silminState->fracFluids = (double *) calloc((unsigned) 2, sizeof(double));
            silminState->fractionateSol = TRUE;

            for (i=0; i<npc; i++) {
                if (solids[i].type == PHASE) {
                    int phaseStrLen = (int) strlen(solids[i].label);
                    if (((len-18-phaseStrLen-1)  >= 0) && !strncmp(&line[18], solids[i].label, phaseStrLen)) {
                        if ( solids[i].nr == 0 || (solids[i].nr > 0 && solids[i].convert != NULL)) {
                            (silminState->fracSolids)[i] = silminState->fracOut;
                        }
                        break;
                    }
                }
            }
#endif
        }

#ifdef PHMELTS_ADJUSTMENTS
        else if (!strncmp(&line[6],  "multiple liquids",    MIN((len-6), 16))) silminState->incLiquids     = nlc;
#else
        else if (!strncmp(&line[6],  "multiple liquids",    MIN((len-6), 16))) silminState->multipleLiqs   = TRUE;
#endif
#ifdef ALPHAMELTS_UPDATE_SYSTEM
        else if (!strncmp(&line[6],  "isenthalpic",   MIN((len-6), 11))) { silminState->isenthalpic  = TRUE;
            silminState->isentropic  = FALSE; silminState->isochoric  = FALSE; }
        else if (!strncmp(&line[6],  "isentropic",    MIN((len-6), 10))) { silminState->isentropic   = TRUE;
            silminState->isenthalpic = FALSE; silminState->isochoric  = FALSE; }
        else if (!strncmp(&line[6],  "isochoric",     MIN((len-6),  9))) { silminState->isochoric    = TRUE;
            silminState->isenthalpic = FALSE; silminState->isentropic = FALSE; }
#elif !defined(PHMELTS_ADJUSTMENTS)
        /* If coming from MATLAB/Python then these constraints are handled differently */
        else if (!strncmp(&line[6],  "isenthalpic",   MIN((len-6), 11))) silminState->isenthalpic    = TRUE;
        else if (!strncmp(&line[6],  "isentropic",    MIN((len-6), 10))) silminState->isentropic     = TRUE;
        else if (!strncmp(&line[6],  "isochoric",     MIN((len-6),  9))) silminState->isochoric      = TRUE;
#endif
        else { READ_ERROR }
#ifdef PHMELTS_ADJUSTMENTS
        /* -> output record */
    } else if (!strncmp(line, "output: ",                   MIN(len, 8)) ||
                            !strncmp(line, "set output: ",                   MIN(len, 12))) {

        if      (!strncmp(&line[8],  "none",   MIN((len-8), 4))) { silminState->txtOutput = TEXT_NONE; }
        else if (!strncmp(&line[8],  "tbl",    MIN((len-8), 3))) { silminState->txtOutput = TEXT_TABLE; }
        else if (!strncmp(&line[8],  "both_0", MIN((len-8), 6))) { silminState->txtOutput = TEXT_BOTH_0; }
        else if (!strncmp(&line[8],  "both",   MIN((len-8), 4))) { silminState->txtOutput = TEXT_BOTH_1; }
        else if (!strncmp(&line[8],  "txt_0",  MIN((len-8), 5))) { silminState->txtOutput = TEXT_ALPHA_0; }
        else if (!strncmp(&line[8],  "txt",    MIN((len-8), 3))) { silminState->txtOutput = TEXT_ALPHA_1; }
        else { READ_ERROR }
#endif
        /* -> assimilate a solid phase record */
    } else if (!strncmp(line, "assimilant: ",            MIN(len,12))) {

        silminState->assimilate = TRUE;
        if (silminState->nDspAssimComp == NULL) silminState->nDspAssimComp = (int *)     calloc((unsigned) (npc+nc), sizeof(int));
        if (silminState->dspAssimComp  == NULL) silminState->dspAssimComp  = (double **) calloc((unsigned) (npc+nc), sizeof(double *));

        if        (!strncmp(&line[12],  "units ",       MIN((len-12),  6))) {
            if      (!strncmp(&line[18],  "vol %%", MIN((len-18), 5))) silminState->dspAssimUnits = ASSIM_PADB_UNITS_VOLUME;
            else if (!strncmp(&line[18],  "wt %%",  MIN((len-18), 4))) silminState->dspAssimUnits = ASSIM_PADB_UNITS_WEIGHT;
            else { READ_ERROR }
        } else if (!strncmp(&line[12],  "temperature ", MIN((len-12), 12))) {
            if (sscanf(&line[24], "%f", &temporary) == EOF) { READ_ERROR }
            silminState->dspAssimT = (double) temporary;
            silminState->assimT = (double) temporary; /* for now keep in C */
        } else if (!strncmp(&line[12],  "mass ",        MIN((len-12),  5))) {
            if (sscanf(&line[17], "%f", &temporary) == EOF) { READ_ERROR }
            silminState->dspAssimMass = (double) temporary;
        } else if (!strncmp(&line[12],  "increments ",  MIN((len-12), 11))) {
            if (sscanf(&line[23], "%f", &temporary) == EOF) { READ_ERROR }
            silminState->dspAssimInc = (int) round(temporary);
        } else if (!strncmp(&line[12],  "liquid mass ",  MIN((len-12), 12))) {
            if (sscanf(&line[24], "%f", &temporary) == EOF) { READ_ERROR }
#ifdef PHMELTS_ADJUSTMENTS
            if (npOld == NULL || nsOld == NULL || oxideWt == NULL) {READ_ERROR}
            *npOld = -1;
#endif
            silminState->dspAssimLiqM = (double) temporary;
        } else {
#ifdef PHMELTS_ADJUSTMENTS
            if (npOld == NULL || nsOld == NULL || oxideWt == NULL) {READ_ERROR}
            for (i=0; i<npc; i++) {
                if (solids[i].type == PHASE) np = i;
                if (!strncmp(&line[12], solids[i].label, MIN((len-12), (int) strlen(solids[i].label)))) {
                    int l = (int) strlen(solids[i].label);
                    if (!strncmp(&line[12+l],  " mass ",  MIN((len-12-l), 6))) {
                        *npOld = i; *nsOld = (silminState->nDspAssimComp)[i]+1;
                        if (sscanf(&line[12 + (int) strlen(solids[i].label) + 6], "%f", &temporary) == EOF) { READ_ERROR }
                    }
                    else {
                        /* original code would not distinguish between cpx and opx endmembers? */
                        if ((solids[i].type != PHASE) && ((silminState->nDspAssimComp)[i] == (silminState->nDspAssimComp)[np])) continue;
                        if (sscanf(&line[12 + (int) strlen(solids[i].label)], "%f", &temporary) == EOF) { READ_ERROR }
                    }
#else
            for (i=0; i<npc; i++) {
                if (!strncmp(&line[12], solids[i].label, MIN((len-12), (int) strlen(solids[i].label)))) {
                    if (sscanf(&line[12 + (int) strlen(solids[i].label)], "%f", &temporary) == EOF) { READ_ERROR }
#endif
                    if ((ns = (silminState->nDspAssimComp)[i]) == 0) (silminState->dspAssimComp)[i] = (double *) malloc((unsigned) sizeof(double));
                    else (silminState->dspAssimComp)[i] = (double *) REALLOC((silminState->dspAssimComp)[i], (unsigned) (ns+1)*sizeof(double));
                    (silminState->nDspAssimComp)[i]++;
                    (silminState->dspAssimComp)[i][ns] = (double) temporary;
                    break;
                }
            }
            if (i == npc) {
                for (i=0; i<nc; i++) {
                    for (j=0; j < ((int)strlen(bulkSystem[i].label)); j++) label[j] = tolower((bulkSystem[i].label)[j]);
                    label[j] = '\0';
                    if (!strncmp(&line[12], label, MIN((len-12), (int) strlen(label)))) {
                        if (sscanf(&line[12 + (int) strlen(label)], "%f", &temporary) == EOF) { READ_ERROR }
#ifdef PHMELTS_ADJUSTMENTS
                        if (*npOld < 0) { /* liquid or bulk */
#endif
                            if ((ns = (silminState->nDspAssimComp)[npc+i]) == 0) (silminState->dspAssimComp)[npc+i] = (double *) malloc((unsigned) sizeof(double));
                            else  (silminState->dspAssimComp)[npc+i] = (double *) REALLOC((silminState->dspAssimComp)[npc+i], (unsigned)(ns+1)*sizeof(double));
                            (silminState->nDspAssimComp)[npc+i]++;
                            (silminState->dspAssimComp)[npc+i][ns] = (double) temporary;
                            break;
#ifdef PHMELTS_ADJUSTMENTS
                        } else {
                            oxideWt[i]  = (double) temporary; /* solid assimilant expressed in oxides */
                            break;
                        }
#endif
                    }
                }
                if (i == nc) { READ_ERROR }
            }
        }
    } else if (len > 0 && (line[0] == ' ' || line[0] == '\r' || line[0] == '\n')) {
        /* -> illegal record */
    } else {
        { READ_ERROR }
    }
#ifdef PHMELTS_ADJUSTMENTS
    if (len > 0 && (line[0] == ' ' || line[0] == '\r' || line[0] == '\n')) {
        printf("Blank line(?) encountered during input file read will be ignored, offending record: %s\n", line); \
        fprintf(stderr, "Ignored line... %s", line);
    }
    else {
        fprintf(stderr, "Processed line... %s", line);
    }
#endif
    return GET_INPUT_SUCCESS;

}

/*
    static int batchInputDataFromXmlFile(char *fileName)
    static void putOutputDataToXmlFile(char *outputFile)
    static void putStatusDataToXmlFile(char *statusFile)
*/

/* based on doBatchFractionation(void) */
//void doBatchFractionation(void) {

#ifdef PHMELTS_ADJUSTMENTS
#define MASSINS ((fracIn) ? (silminState->solidComp)[i][ns]*fracIn : MASSIN)
#define MASSINL ((fracIn) ? (silminState->liquidComp)[nl][i]*fracIn : MASSIN)
#else
#define MASSINS MASSIN
#define MASSINL MASSIN
#endif

void doBatchFractionation(double fracOut) {
    int i, j, k, ns, nl;
    int hasLiquid = ((silminState != NULL) && (silminState->liquidMass != 0.0));
    double fracIn = 0.0;

#ifdef PHMELTS_ADJUSTMENTS
    /* Use fracMass, fracSComp and fracLComp only for most recent fractionation (zeroed in silmin.c)*/
    /* Frac options may have changed at menu (unlike Melts-batch) */
    if (previousSilminState == NULL) {
        previousSilminState = copySilminStateStructure(silminState, previousSilminState);
        previousSilminState->fracMass = 0.0;/* should already be zero but just in case */
    }
    else {
        previousSilminState->fractionateSol = silminState->fractionateSol;
        previousSilminState->fractionateFlu = silminState->fractionateFlu;
        previousSilminState->fractionateLiq = silminState->fractionateLiq;
    }
    if ((previousSilminState->fractionateSol || previousSilminState->fractionateFlu) && previousSilminState->fracSComp == NULL) {
        previousSilminState->fracSComp    = (double **) calloc((unsigned) npc, sizeof(double *));
        previousSilminState->nFracCoexist = (int *) calloc((unsigned) npc, sizeof(int));
    }
    if (previousSilminState->fractionateLiq && previousSilminState->fracLComp == NULL) {
        previousSilminState->fracLComp = (double *) calloc((unsigned) nlc, sizeof(double));
    }

    /* Allow fractionation of fluid if no liquid */
    if (silminState->fractionateSol && !silminState->fractionateFlu && !hasLiquid) printf("...Cannot do solid fractionation without a liquid phase.\n");

    /* Solid Phase Fractionation */
    if (silminState->fractionateFlu || (silminState->fractionateSol && hasLiquid)) {
#else
    if ((silminState->fractionateSol || silminState->fractionateFlu) && !hasLiquid) fprintf(stderr, "...Cannot do solid/fluid fractionation without a liquid phase.\n");

    /* Solid Phase Fractionation */
    if ((silminState->fractionateSol || silminState->fractionateFlu) && hasLiquid) {
#endif
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
#ifdef PHMELTS_ADJUSTMENTS
                    (previousSilminState->nFracCoexist)[i] = ns;
                    if (nf == 0) {
                        (previousSilminState->fracSComp)[i] = (double *) calloc((size_t) ns, sizeof(double));
                        if (solids[i].na > 1) for (j=0; j<solids[i].na; j++) (previousSilminState->fracSComp)[i+1+j] = (double *) calloc((size_t) ns, sizeof(double));
                    } else {
                        (previousSilminState->fracSComp)[i] = (double *) REALLOC((previousSilminState->fracSComp)[i], (size_t) ns*sizeof(double));
                        for (j=nf; j<ns; j++) (previousSilminState->fracSComp)[i][j] = 0.0;
                        if (solids[i].na > 1) for (j=0; j<solids[i].na; j++) {
                                (previousSilminState->fracSComp)[i+1+j] = (double *) REALLOC((previousSilminState->fracSComp)[i+1+j], (size_t) ns*sizeof(double));
                            for (k=nf; k<ns; k++) (previousSilminState->fracSComp)[i+1+j][k] = 0.0;
                        }
                    }
#endif
                }
        }
        for (i=0; i<npc; i++) {
            if (!silminState->fractionateFlu && !strcmp((char *) solids[i].label, "fluid")) continue;
            if (!silminState->fractionateSol &&  strcmp((char *) solids[i].label, "fluid")) continue;

            for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {

#ifdef PHMELTS_ADJUSTMENTS
                if (silminState->fractionateFlu && !silminState->fracFluids[ns] && !strcmp((char *) solids[i].label, "fluid")) continue;

                if (fracOut >= 1.0) { /* individual fractionation coefficients or default MASSIN */
                    if (!strcmp((char *) solids[i].label, "fluid")) fracIn = MAX(1.0 - silminState->fracFluids[ns], 0.0);
                    else    fracIn = MAX(1.0 - silminState->fracSolids[i], 0.0); /* may evenutally pack like incSolids */
                }
                else { /* fractionation controlled by mass or vol frac, 'coefficients' act as flags */
                         if (!strcmp((char *) solids[i].label, "fluid") && (silminState->fracFluids[ns] == 0.0)) fracIn = 1.0;
                    else if (silminState->fracSolids[i] == 0.0) fracIn = 1.0; /* may evenutally pack like incSolids */
                    else    fracIn = MIN(1.0 - fracOut, 1.0);
                }
#endif
                if (!fracIn && (silminState->solidComp)[i][ns] < MASSIN) continue;       // from MCS

                if (solids[i].na == 1) {
                    (silminState->fracSComp)[i][ns] += (silminState->solidComp)[i][ns]-MASSINS;
#ifdef PHMELTS_ADJUSTMENTS
                    (previousSilminState->fracSComp)[i][ns] += (silminState->solidComp)[i][ns]-MASSINS;
#endif
                    if (silminState->fo2Path != FO2_NONE) silminState->oxygen -= (oxygen.solToOx)[i]*((silminState->solidComp)[i][ns]-MASSINS);

                    silminState->fracMass += ((silminState->solidComp)[i][ns]-MASSINS)*solids[i].mw;
#ifdef PHMELTS_ADJUSTMENTS
                    previousSilminState->fracMass += ((silminState->solidComp)[i][ns]-MASSINS)*solids[i].mw;
#endif
                    for (j=0; j<nc; j++) (silminState->bulkComp)[j] -= (solids[i].solToOx)[j]*((silminState->solidComp)[i][ns]-MASSINS);

                    /* Subtract off H, S or V if appropriate                          */
                    if (silminState->isenthalpic && (silminState->refEnthalpy != 0.0))
                        silminState->refEnthalpy -= ((silminState->solidComp)[i][ns]-MASSINS)*(solids[i].cur).h;
                    if (silminState->isentropic && (silminState->refEntropy != 0.0))
                        silminState->refEntropy -= ((silminState->solidComp)[i][ns]-MASSINS)*(solids[i].cur).s;
                    if (silminState->isochoric && (silminState->refVolume != 0.0))
                        silminState->refVolume -= ((silminState->solidComp)[i][ns]-MASSINS)*(solids[i].cur).v;

                    (silminState->solidComp)[i][ns] = MASSINS;
                } else {
                    double moleF, totalMoles=0.0;
                    (silminState->fracSComp)[i][ns] += (silminState->solidComp)[i][ns] - MASSINS;
#ifdef PHMELTS_ADJUSTMENTS
                    (previousSilminState->fracSComp)[i][ns] += (silminState->solidComp)[i][ns] - MASSINS;
#endif
                    for (j=0; j<solids[i].na; j++) {
                        moleF = (silminState->solidComp)[i+1+j][ns]/(silminState->solidComp)[i][ns];
                        m[j] = (silminState->solidComp)[i+1+j][ns] - MASSINS*moleF;
                        totalMoles += m[j];
                        (silminState->fracSComp)[i+1+j][ns] += m[j];
#ifdef PHMELTS_ADJUSTMENTS
                        (previousSilminState->fracSComp)[i+1+j][ns] += m[j];
#endif
                        if (silminState->fo2Path != FO2_NONE) silminState->oxygen -= (oxygen.solToOx)[i+1+j]*m[j];
                        silminState->fracMass += m[j]*solids[i+1+j].mw;
#ifdef PHMELTS_ADJUSTMENTS
                        previousSilminState->fracMass += m[j]*solids[i+1+j].mw;
#endif
                        for (k=0; k<nc; k++) (silminState->bulkComp)[k] -= (solids[i+1+j].solToOx)[k]*m[j];
                        (silminState->solidComp)[i+1+j][ns] = MASSINS*moleF;

                        /* Subtract off H, S or V if appropriate                        */
                        if (silminState->isenthalpic && (silminState->refEnthalpy != 0.0)) silminState->refEnthalpy -= m[j]*(solids[i+1+j].cur).h;
                        if (silminState->isentropic && (silminState->refEntropy != 0.0))   silminState->refEntropy  -= m[j]*(solids[i+1+j].cur).s;
                        if (silminState->isochoric && (silminState->refVolume != 0.0))     silminState->refVolume   -= m[j]*(solids[i+1+j].cur).v;
                    }
                    (silminState->solidComp)[i][ns] = MASSINS;

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

        for (i=0; i<nc; i++) {
            if ((silminState->bulkComp)[i] != 0.0 && (silminState->bulkComp)[i] <  MASSOUT && bulkSystem[i].type != FE2O3) {
#ifdef ALPHAMELTS_UPDATE_SYSTEM
                printf("  Moles of %5.5s in system (%g) < %g\n.", bulkSystem[i].label, (silminState->bulkComp)[i], MASSOUT);
#else
                fprintf(stderr, "  Moles of %5.5s in system (%g) < %g\n.", bulkSystem[i].label, (silminState->bulkComp)[i], MASSOUT);
#endif
                (silminState->bulkComp)[i] = 0.0;
                for (j=0; j<nlc; j++) if ((liquid[j].liqToOx)[i] != 0.0) {
                        for (nl=0; nl<silminState->nLiquidCoexist; nl++) (silminState->liquidComp)[nl][j] = 0.0;
                        fprintf(stderr, "    Moles of %s in liquid(s) set to zero.\n", liquid[j].label);
                    }
                for (j=0; j<npc; j++) {
                    for (ns=0; ns<(silminState->nSolidCoexist)[j]; ns++) {
                        if (solids[j].na == 1) {
                            if ((solids[j].solToOx)[i] != 0.0) {
                                (silminState->solidComp)[j][ns] = 0.0;
                                fprintf(stderr, "    Moles of %s in solid set to zero.\n", solids[j].label);
                            }
                        } else {
                            for (k=0; k<solids[j].na; k++) {
                                if ((solids[j+1+k].solToOx)[i] != 0.0) {
                                    (silminState->solidComp)[j][ns] -= (silminState->solidComp)[j+1+k][ns];
                                    (silminState->solidComp)[j+1+k][ns] = 0.0;
                                    fprintf(stderr, "    Moles of %s in %s solid set to zero.\n", solids[j+1+k].label, solids[j].label);
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
#ifdef ALPHAMELTS_UPDATE_SYSTEM
    if (silminState->fractionateLiq && !hasLiquid) printf("...Cannot do liquid fractionation without a liquid phase.\n");
#else
    if (silminState->fractionateLiq && !hasLiquid) fprintf(stderr, "...Cannot do liquid fractionation without a liquid phase.\n");
#endif

    if (silminState->fractionateLiq && hasLiquid) {
        double *m = (double *) malloc((size_t) nlc*sizeof(double));
        double *r = (double *) malloc((size_t) nlc*sizeof(double));
        for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
            double refMoles, totalMoles;
            for (i=0, refMoles=0.0; i<nlc; i++) refMoles += (silminState->liquidComp)[nl][i];

#ifdef PHMELTS_ADJUSTMENTS
            if (fracOut >= 1.0) {
                fracIn = MAX(1.0 - silminState->fracLiquids[nl], 0.0);
            }
            else {
                if (silminState->fracLiquids[nl] == 0.0) fracIn = 1.0;
                else fracIn = MIN(1.0 - fracOut, 1.0);
            }
#endif
            if (!fracIn && refMoles < MASSIN) continue;       // from MCS

            for (i=0, totalMoles=0.0; i<nlc; i++) {
                if (((silminState->liquidComp)[nl][i] != 0.0) && (refMoles != 0.0)) {
                    double mw;
                    double moleF = (silminState->liquidComp)[nl][i]/refMoles;

                    for (j=0, mw = 0.0; j<nc; j++) mw += (liquid[i].liqToOx)[j]*bulkSystem[j].mw;
                    m[i] = (silminState->liquidComp)[nl][i] - MASSINL*moleF;
                    totalMoles += m[i];
#ifdef PHMELTS_ADJUSTMENTS
                    /* Use fracLComp only for most recent fractionation (like 'extractions' before) */
                    (silminState->fracLComp)[i] = m[i];
                    (previousSilminState->fracLComp)[i] += m[i];
#else
                    (silminState->fracLComp)[i] += m[i];
#endif
                    if (silminState->fo2Path != FO2_NONE) silminState->oxygen -= (oxygen.liqToOx)[i]*m[i];
                    silminState->fracMass += m[i]*mw;
#ifdef PHMELTS_ADJUSTMENTS
                    previousSilminState->fracMass += m[i]*mw;
#endif
                    for (j=0; j<nc; j++) (silminState->bulkComp)[j] -= (liquid[i].liqToOx)[j]*m[i];
                    (silminState->liquidComp)[nl][i] = MASSINL*moleF;

                    /* Subtract off H, S or V if appropriate                            */
                    if (silminState->isenthalpic && (silminState->refEnthalpy != 0.0)) silminState->refEnthalpy -= m[i]*(liquid[i].cur).h;
                    if (silminState->isentropic  && (silminState->refEntropy  != 0.0)) silminState->refEntropy  -= m[i]*(liquid[i].cur).s;
                    if (silminState->isochoric   && (silminState->refVolume   != 0.0)) silminState->refVolume       -= m[i]*(liquid[i].cur).v;
                } else m[i] = 0.0;
            }

            /* Subtract off H, S or V if appropriate                      */
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

        for (i=0; i<nc; i++) {
            if ((silminState->bulkComp)[i] != 0.0 && (silminState->bulkComp)[i] <  MASSOUT && bulkSystem[i].type != FE2O3) {
#ifdef ALPHAMELTS_UPDATE_SYSTEM
                printf("  Moles of %5.5s in system (%g) < %g\n.", bulkSystem[i].label, (silminState->bulkComp)[i], MASSOUT);
#else
                fprintf(stderr, "  Moles of %5.5s in system (%g) < %g\n.", bulkSystem[i].label, (silminState->bulkComp)[i], MASSOUT);
#endif

                (silminState->bulkComp)[i] = 0.0;
                for (j=0; j<nlc; j++) if ((liquid[j].liqToOx)[i] != 0.0) {
                        for (nl=0; nl<silminState->nLiquidCoexist; nl++) (silminState->liquidComp)[nl][j] = 0.0;
#ifdef ALPHAMELTS_UPDATE_SYSTEM
                        printf("    Moles of %s in liquid(s) set to zero.\n", liquid[j].label);
#else
                        fprintf(stderr, "    Moles of %s in liquid(s) set to zero.\n", liquid[j].label);
#endif
                    }
                for (j=0; j<npc; j++) {
                    for (ns=0; ns<(silminState->nSolidCoexist)[j]; ns++) {
                        if (solids[j].na == 1) {
                            if ((solids[j].solToOx)[i] != 0.0) {
                                (silminState->solidComp)[j][ns] = 0.0;
#ifdef ALPHAMELTS_UPDATE_SYSTEM
                                printf("    Moles of %s in solid set to zero.\n", solids[j].label);
#else
                                fprintf(stderr, "    Moles of %s in solid set to zero.\n", solids[j].label);
#endif
                            }
                        } else {
                            for (k=0; k<solids[j].na; k++) {
                                if ((solids[j+1+k].solToOx)[i] != 0.0) {
                                    (silminState->solidComp)[j][ns] -= (silminState->solidComp)[j+1+k][ns];
                                    (silminState->solidComp)[j+1+k][ns] = 0.0;
#ifdef ALPHAMELTS_UPDATE_SYSTEM
                                    printf("    Moles of %s in %s solid set to zero.\n", solids[j+1+k].label, solids[j].label);
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
}

int doBatchStateChange(void) {

    int stateChange = FALSE;

    /* Changing T ? */
    if (fabs(silminState->dspTstart - silminState->dspTstop) >=  (silminState->dspTinc != 0.0 ? silminState->dspTinc : 0.001)
            && !(silminState->isenthalpic && (silminState->refEnthalpy != 0.0))
            && !(silminState->isentropic  && (silminState->refEntropy  != 0.0))) {
        stateChange = TRUE;

        if (silminState->dspTstart - silminState->dspTstop < 0.0) { // alphaMELTS has signed increments
            silminState->T += silminState->dspTinc; silminState->dspTstart += fabs(silminState->dspTinc);
        } else {
            silminState->T -= silminState->dspTinc; silminState->dspTstart -= fabs(silminState->dspTinc);
        }

        /* Changing H ? */
    } else if ((fabs(silminState->dspHstop) > 0.0) &&
             (fabs(silminState->dspHstop - silminState->refEnthalpy) >= (silminState->dspHinc != 0.0 ? fabs(silminState->dspHinc) : 0.001)) &&
             (silminState->isenthalpic && (silminState->refEnthalpy != 0.0))) {
        stateChange = TRUE;
        silminState->refEnthalpy += silminState->dspHinc;
        correctTforChangeInEnthalpy();
    } else if ((silminState->dspHstop == 0.0) &&
             (fabs(silminState->dspPstart - silminState->dspPstop) >=  (silminState->dspPinc != 0.0 ? silminState->dspPinc : 0.001)) &&
             (silminState->isenthalpic && (silminState->refEnthalpy != 0.0))) {
        stateChange = TRUE;
        silminState->refEnthalpy += silminState->dspHinc;
        correctTforChangeInEnthalpy();

        /* Could get into loop if dH/dP = 0.0 (e.g. isenthalpic degassing) or for isenthalpic assimilation once all material assimilated */
        /* (fabs(silminState->dspTstart - silminState->dspTstop) >=  (silminState->dspTinc != 0.0 ? silminState->dspTinc : 0.001)) */
    } else if ((silminState->dspHstop == 0.0) && (silminState->dspHinc != 0.0) &&
             (silminState->isenthalpic && (silminState->refEnthalpy != 0.0))) {
        int changeH = (silminState->dspHinc < 0.0) ? (silminState->dspTstart - silminState->dspTstop < 0.0) : (silminState->dspTstart - silminState->dspTstop > 0.0);
        if (changeH) {
            stateChange = TRUE;
            silminState->refEnthalpy += silminState->dspHinc;
            correctTforChangeInEnthalpy();
        }

        /* Changing S ? */
    } else if ((fabs(silminState->dspSstop) > 0.0) &&
             (fabs(silminState->dspSstop - silminState->refEntropy) >= (silminState->dspSinc != 0.0 ? fabs(silminState->dspSinc) : 0.001)) &&
             (silminState->isentropic && (silminState->refEntropy != 0.0))) {
        stateChange = TRUE;
        silminState->refEntropy += silminState->dspSinc;
        correctTforChangeInEntropy();

    } else if ((silminState->dspSstop == 0.0) &&
             (fabs(silminState->dspPstart - silminState->dspPstop) >=  (silminState->dspPinc != 0.0 ? silminState->dspPinc : 0.001)) &&
             (silminState->isentropic  && (silminState->refEntropy  != 0.0))) {
        stateChange = TRUE;
        silminState->refEntropy += silminState->dspSinc;
        correctTforChangeInEntropy();

        /* (fabs(silminState->dspTstart - silminState->dspTstop) >=  (silminState->dspTinc != 0.0 ? silminState->dspTinc : 0.001)) */
    } else if ((silminState->dspSstop == 0.0) && (silminState->dspSinc != 0.0) &&
             (silminState->isentropic  && (silminState->refEntropy  != 0.0))) {
        int changeS = (silminState->dspSinc > 0.0) ? (silminState->dspTstart - silminState->dspTstop < 0.0) : (silminState->dspTstart - silminState->dspTstop > 0.0);
        if (changeS) {
            stateChange = TRUE;
            silminState->refEntropy += silminState->dspSinc;
            correctTforChangeInEntropy();
        }
    }

    /* Changing P ? */
    if (fabs(silminState->dspPstart - silminState->dspPstop) >=  (silminState->dspPinc != 0.0 ? silminState->dspPinc : 0.001)
            && !(silminState->isochoric && (silminState->refVolume != 0.0))) {
        stateChange = TRUE;
        if (silminState->dspPstart - silminState->dspPstop < 0.0) { // alphaMELTS has signed increments
            silminState->P += silminState->dspPinc; silminState->dspPstart += fabs(silminState->dspPinc);
        } else {
            silminState->P -= silminState->dspPinc; silminState->dspPstart -= fabs(silminState->dspPinc);
        }

        /* Changing V ? */
    } else if ((fabs(silminState->dspVstop) > 0.0) &&
             (fabs(silminState->dspVstop - silminState->refVolume) >= (silminState->dspVinc != 0.0 ? fabs(silminState->dspVinc) : 0.0001)) &&
             (silminState->isochoric && (silminState->refVolume != 0.0))) {
        stateChange = TRUE;
        silminState->refVolume += silminState->dspVinc;
        correctPforChangeInVolume();

        /*    } else if (fabs(silminState->dspPstart - silminState->dspPstop) >=  (silminState->dspPinc != 0.0 ? silminState->dspPinc : 0.001) FIX-fO2 */
    } else if ((silminState->dspVstop == 0.0) &&
             (fabs(silminState->dspTstart - silminState->dspTstop) >=  (silminState->dspTinc != 0.0 ? silminState->dspTinc : 0.001)) &&
             (silminState->isochoric && (silminState->refVolume != 0.0))) {
        stateChange = TRUE;
        silminState->refVolume += silminState->dspVinc;
        correctPforChangeInVolume();

        /* (fabs(silminState->dspPstart - silminState->dspPstop) >=  (silminState->dspPinc != 0.0 ? silminState->dspPinc : 0.001) */
    } else if ((silminState->dspVstop == 0.0) && (silminState->dspVinc != 0.0) &&
             (silminState->isochoric && (silminState->refVolume != 0.0))) {
        int changeV = (silminState->dspVinc < 0.0) ? (silminState->dspPstart - silminState->dspPstop < 0.0) : (silminState->dspPstart - silminState->dspPstop > 0.0);
        if (changeV) {
            stateChange = TRUE;
            silminState->refVolume += silminState->dspVinc;
            correctPforChangeInVolume();
        }
    }

    return stateChange;

}

void doBatchAssimilation(void) {

    int i, j, k, ns;
    int hasLiquid = ((silminState != NULL) && (silminState->liquidMass != 0.0));

    /* Assimilation ? - add assimilant to the first liquid, if a liquid is present */
#ifndef PHMELTS_ADJUSTMENTS
    if (silminState->assimilate && silminState->assimMass < silminState->dspAssimMass) {
        double fraction = 1.0/silminState->dspAssimInc;
#else
    if (silminState->assimilate && ((silminState->dspAssimInc == 0 && silminState->dspAssimMass > 0.0) ||
                                                                    (silminState->assimMass < silminState->dspAssimMass))) {
        /* if the number of increments is zero then add mass each time */
        double fraction = (silminState->dspAssimInc) ? 1.0/silminState->dspAssimInc : 1.0;
#endif
        //stateChange = TRUE;
        silminState->assimMass += fraction*silminState->dspAssimMass;
        if (silminState->isenthalpic && (silminState->refEnthalpy != 0.0)) silminState->refEnthalpy += fraction*(silminState->assimTD).h;
        if (silminState->isentropic && (silminState->refEntropy != 0.0))   silminState->refEntropy += fraction*(silminState->assimTD).s;
        if (silminState->isochoric && (silminState->refVolume != 0.0))     silminState->refVolume += fraction*(silminState->assimTD).v;
        silminState->assimInc++;

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
        }
        if (silminState->isentropic && (silminState->refEntropy != 0.0)) {
            correctTforChangeInEntropy();
        }
        if (silminState->isochoric && (silminState->refVolume != 0.0)) {
            correctPforChangeInVolume();
        }
#ifndef PHMELTS_ADJUSTMENTS
    }
#else
    }
#endif
}

#endif /* BATCH_VERSION */

    /*****************/
    /* MAIN FUNCTION */
    /****************/

    /* int main (int argc, char *argv[]) */

    /* end of file INTERFACE.C */
