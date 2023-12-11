/*
**++
**  FACILITY:  alphamelts: Silicate Melts batch path finder
**
**  MODULE DESCRIPTION:
**
**      Calculates isentropic, isothermal, or isobaric
**        equilibrium paths for according to SILMIN,
**        for fixed bulk composition.
**
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "silmin.h"                /* SILMIN structures include file        */
#include "interface.h"

#include "adiabat.h"
#include "alphamelts.h"
#include "phmelts.h"

#define FREE(x) free(x); x = NULL

static int iAmInitialized = FALSE;

static void initializeMenu(void) {

    /* Test for correct allocation in here */
    if (silminState == (SilminState *) NULL) silminState = createSilminState();
    if (startState  == (SilminState *) NULL) startState  = allocSilminStatePointer();
    if (oldState    == (SilminState *) NULL) oldState    = allocSilminStatePointer();
    if (states      == (SilminState *) NULL) states      = allocSilminStatePointer();
	iAmInitialized = TRUE;

}





void menu_option0(void) {

    if (states      != (SilminState *) NULL) {
        adiabat_1ph(0, 0, 0, -1); //clean up states - does not work if last iteration failed???
        destroySilminStateStructure(states);
    }
    if (startState  != (SilminState *) NULL) destroySilminStateStructure(startState);
    if (oldState    != (SilminState *) NULL) destroySilminStateStructure(oldState);
    if (silminState != (SilminState *) NULL) destroySilminStateStructure(silminState);

}

int menu_option1(char *filename) {

    if (silminState != (SilminState *) NULL) {
        destroySilminStateStructure(silminState);
        silminState = createSilminState();
    }

    if (!iAmInitialized) initializeMenu();
    guessFlag = FALSE;

    if(getInputDataFromFile(filename) == GET_INPUT_ERROR_BAD_FILE) return FALSE;
    else assignProblemStatics();

    return TRUE;

}

int menu_option2(double dspTstart, double dspTstop, double dspTinc,
    double dspPstart, double dspPstop, double dspPinc) {

    if(!iAmInitialized) return FALSE;

    silminState->T = dspTstart + 273.15;
    silminState->dspTstart = dspTstart;
    silminState->dspTstop = dspTstop;
    silminState->dspTinc = dspTinc;

    silminState->P = dspPstart;
    silminState->dspPstart = dspPstart;
    silminState->dspPstop = dspPstop;
    silminState->dspPinc = dspPinc;

    assignProblemStatics();

    return TRUE;

}

int menu_option3(int equilibriumGuess) {

    if(!iAmInitialized) return FALSE;

	(void) adiabat_0ph(equilibriumGuess);

    return TRUE;


}


int menu_option4(int equilibriumGuess, int saveAll, int iterMax) {

    if(!iAmInitialized) return FALSE;

    int success = FALSE;
     /* Execute (follow path, mineral isograd or melt contour) */


    success = adiabat_1ph(equilibriumGuess, saveAll, iterMax, FALSE);

	/* If have assimilation do single calc to get assimilant in (liquidus will be before though) */

    /*use allocProblemStatics after first calc (once 'Hstart' has been set) */
    /* turn off ref quantities afterwards ????*/



            switch (success) {
                case 0: printf("Not all calculations performed!\n"); break;
                default: printf("Successful return from alphamelts...\n");
            }


    return success;

}



int menu_option5(int fo2Path, double fo2Delta) {

    if(!iAmInitialized) return FALSE;

    // Need some sanity checks for NONE and ABS etc.
    if(fo2Path < -1 || fo2Path > 5) {
        printf("WARNING: Invalid choice for fO2 Path, setting to NONE.\n");
        fo2Delta = 0.0;
    }
    if(fo2Path == FO2_NONE && fo2Delta != 0.0) {
        printf("WARNING: fO2 Path set to NONE, so setting fO2 Offset to 0.0.\n");
        fo2Delta = 0.0;
    }

    silminState->fo2Path = fo2Path;
    silminState->fo2Delta = fo2Delta;

    // Absolute in here too
    return TRUE;

}

int menu_option6() {

    if(!iAmInitialized) return FALSE;


            /*else {
	printf("Activity of Water to impose (0 < aH2O <= 1) ");
	printf("or 0 for unbuffered: ");
	scanf("%lf", &silminState->aH2O);
	if (silminState->aH2O <= 0.0) {
	  silminState->aH2O = 0.0; silminState->H2Obuffer = FALSE;
	} else {
	  if (silminState->aH2O > 1.0) silminState->aH2O = 1.0;
	  silminState->H2Obuffer = TRUE;
	  }*/

    return TRUE;

}

int menu_option7(int constraint, double newRef, double newStop, double newInc) {

    /* Impose initial entropy, enthalpy or volume */
    int isenthalpic = FALSE, isentropic = FALSE, isochoric = FALSE;

    if(!iAmInitialized) return FALSE;

    switch (constraint) {
    case 0:
        isenthalpic = FALSE; isentropic = FALSE; isochoric = FALSE;
        break;
    case 1:
        isenthalpic = TRUE; isentropic = FALSE; isochoric = FALSE;
        break;
    case 2:
        isenthalpic = FALSE; isentropic = TRUE; isochoric = FALSE;
        break;
    case 3:
        isenthalpic = FALSE; isentropic = FALSE; isochoric = TRUE;
        break;
    }

    silminState->isenthalpic = isenthalpic;
    silminState->isentropic  = isentropic;
    silminState->isochoric   = isochoric;

    // NEEDS CHECKING FOR INTERPLAY BETWEEN DIFFERENT MENU OPTIONS
    if (isenthalpic) {
        silminState->refEnthalpy = newRef;
        silminState->dspHstop = newStop;
        silminState->dspHinc = newInc;
    }
    else if (isentropic) {
        silminState->refEntropy = newRef;
        silminState->dspSstop = newStop;
        silminState->dspSinc = newInc;
    }
    else if (isochoric) {
        silminState->refVolume = newRef;
        silminState->dspVstop = newStop;
        silminState->dspVinc = newInc;
    }

    return TRUE;

}


int menu_option8(int index, int jphase, int incSolids, int fracSolids, int minType, double newMin) {

    int i, j, k, np;

    if(!iAmInitialized) return FALSE;

    if (silminState->fracSolids == NULL) silminState->fracSolids = (double *) calloc((unsigned) npc, sizeof(double));

    /* all */
    if (index < 0) {
        if (incSolids > 0) {
            for (i=0, np=0; i<npc; i++) if (solids[i].type == PHASE) {
                int len = strlen(solids[i].label);
                if (strncmp(solids[i].label, "water", MIN(len, 5)) && strncmp(solids[i].label, "fluid", MIN(len, 5))) {
                    if (solids[i].inStdSet)
                        (silminState->incSolids)[np] = MIN(incSolids, solids[i].na);
                    else
                        (silminState->incSolids)[np] = FALSE;
                }
                np++;
            }
        }
        else if (incSolids == -1) { /* unlimited */
            for (i=0, np=0; i<npc; i++) if (solids[i].type == PHASE) {
                int len = strlen(solids[i].label);
                if (strncmp(solids[i].label, "water", MIN(len, 5)) && strncmp(solids[i].label, "fluid", MIN(len, 5))) {
                    if (solids[i].inStdSet)
                        (silminState->incSolids)[np] = solids[i].na;
                    else
                        (silminState->incSolids)[np] = FALSE;
                }
                np++;
            }
        }
        else if (!incSolids) {
            for (i=0, np=0; i<npc; i++) if (solids[i].type == PHASE) {
                int len = strlen(solids[i].label);
                if (strncmp(solids[i].label, "water", MIN(len, 5)) && strncmp(solids[i].label, "fluid", MIN(len, 5))) {
                    (silminState->incSolids)[np] = FALSE; np++;
                }
            }
        }

        silminState->fractionateSol = fracSolids;

        switch (minType) {
            case 1: /* by mass */
                silminState->maxF = 1.0 - newMin;
                break;
            case 2: /* global ratio */
                if (newMin > 0.0) silminState->fracOut = newMin;
                for (i = 0; i<npc; i++) if (solids[i].type == PHASE) {
                    int len = strlen(solids[i].label);
                    if (strncmp(solids[i].label, "water", MIN(len, 5)) && strncmp(solids[i].label, "fluid", MIN(len, 5))) {
                        silminState->fracSolids[i] = newMin;
                    }
                }
                break;
            case 3: /* individual ratio */
                if (newMin > 0.0) silminState->fracOut = newMin;
                break;
            default:
                silminState->maxF = 1.0;
                silminState->fracOut = 1.0;
                for (i = 0; i<npc; i++) if (solids[i].type == PHASE) {
                    int len = strlen(solids[i].label);
                    if (strncmp(solids[i].label, "water", MIN(len, 5)) && strncmp(solids[i].label, "fluid", MIN(len, 5))) {
                        silminState->fracSolids[i] = 1.0;
                    }
                }
                break;
        }

    }
    /* specific phase */
    else {

        if (incSolids >= 0) {
            silminState->incSolids[jphase] = MIN(incSolids, solids[index].na);
        }
        else if (incSolids == -1) { /* unlimited */
            silminState->incSolids[jphase] = solids[index].na;
        }

        if(!silminState->incSolids[jphase]) {
	    for (j=0; j<(silminState->nSolidCoexist)[index]; j++) {
	      (silminState->nSolidCoexist)[index] = 0;
	      (silminState->solidComp)[index][j] = 0.0;
	      if (solids[index].na != 1) for (k=0; k<solids[index].na; k++) (silminState->solidComp)[index+1+k][j] = 0.0;
  	  }
	    correctXforChangeInBulkComp();
        }

        /* fluid moved into the liquid menu */
        if (!silminState->fractionateSol) silminState->fractionateSol = fracSolids;

        switch (minType) {
            case 1: /* by mass */
                silminState->maxF = 1.0 - newMin;
                break;
            case 2: /* global ratio */
                if (newMin > 0.0) silminState->fracOut = newMin;
                silminState->fracSolids[index] = newMin;
                break;
            case 3: /* individual ratio */
                if (newMin > 0.0) silminState->fracOut = newMin;
                break;
            default:
                silminState->maxF = 1.0;
                silminState->fracOut = 1.0;
                silminState->fracSolids[index] = 1.0;
                break;
        }

    }
    return TRUE;

}

int menu_option9(int index, int jphase, int incLiquids, int fracLiquids, int minType, double newMin) {

    int i, j, k;

    if(!iAmInitialized) return FALSE;

    if (index < 0) {

        if (silminState->fracLiquids == NULL) silminState->fracLiquids = (double *) calloc((unsigned) nlc, sizeof(double));

        if (jphase) {/* a particular liquid (jphase is number of liquids) */
            if (incLiquids >= 0) {
                silminState->incLiquids = MIN(incLiquids, nlc);
            }
            else if (incLiquids == -1) { /* unlimited */
                silminState->incLiquids = nlc;
            }

            if (!silminState->incLiquids) {
                for (j=1; j<silminState->nLiquidCoexist; j++) {
                        free((silminState->liquidComp )[j]); (silminState->liquidComp )[j] = NULL;
                        free((silminState->liquidDelta)[j]); (silminState->liquidDelta)[j] = NULL;
                }
                for (i=0;i<nlc;i++) silminState->liquidComp[0][i] = 0.0;
                //for (i=0;i<nc;i++)  silminState->dspLiquidComp[i] = 0.0;
                silminState->nLiquidCoexist = 1; // else memory leak
                silminState->liquidMass = 0.0;
                correctXforChangeInBulkComp();
            }

            if (!silminState->fractionateLiq) silminState->fractionateLiq = fracLiquids;

            switch (minType) {
                case 1: /* by vol */
                    silminState->minLiqPhi = newMin;
                    break;
                case 2: /* by mass */
                    silminState->minF = newMin;
                    break;
                case 3: /* global ratio */
                    silminState->fracOut = newMin;
                    silminState->fracLiquids[jphase-1] = newMin*fracLiquids;
                    break;
                case 4: /* individual ratio */
                    silminState->fracLiquids[jphase-1] = newMin*fracLiquids;
                    break;
                default:
                    silminState->minLiqPhi = 0.0; silminState->minF = 0.0;
                    silminState->fracLiquids[jphase-1] = 1.0;
                    break;
            }

        }
        else { /* all liquid */
            if (incLiquids >= 0) {
                silminState->incLiquids = MIN(incLiquids, nlc);
            }
            else if (incLiquids == -1) {
                silminState->incLiquids = nlc;
            }

            if (!silminState->incLiquids) {
                for (j=1; j<silminState->nLiquidCoexist; j++) {
                    free((silminState->liquidComp )[j]); (silminState->liquidComp )[j] = NULL;
                    free((silminState->liquidDelta)[j]); (silminState->liquidDelta)[j] = NULL;
                }
                for (i=0;i<nlc;i++) silminState->liquidComp[0][i] = 0.0;
                //for (i=0;i<nc;i++)  silminState->dspLiquidComp[i] = 0.0;
                silminState->nLiquidCoexist = 1; // else memory leak
                silminState->liquidMass = 0.0;
                correctXforChangeInBulkComp();
            }

            silminState->fractionateLiq = fracLiquids;

            switch (minType) {
                case 1: /* by vol */
                    silminState->minLiqPhi = newMin;
                    break;
                case 2: /* by mass */
                    silminState->minF = newMin;
                    break;
                case 3: /* global ratio */
                    if (newMin > 0.0) silminState->fracOut = newMin;
                    for (j = 0; j<nlc; j++) silminState->fracLiquids[j] = newMin;
                    break;
                case 4: /* individual ratio */
                    if (newMin > 0.0) silminState->fracOut = newMin;
                    break;
                default:
                    silminState->minLiqPhi = 0.0; silminState->minF = 0.0;
                    silminState->fracOut = 1.0;
                    for (j = 0; j<nlc; j++) silminState->fracLiquids[j] = 1.0;
                    break;
            }

        }

    }
    else { /* jphase is used to locate the correct fluid phase */

        if (silminState->fracFluids == NULL) silminState->fracFluids = (double *) calloc((unsigned) 2, sizeof(double));

        if (incLiquids >= 0) silminState->incSolids[jphase] = MIN(incLiquids, solids[index].na);
        else if (incLiquids == -1) silminState->incSolids[jphase] = solids[index].na;

        if(!silminState->incSolids[jphase]) {
	    for (j=0; j<(silminState->nSolidCoexist)[index]; j++) {
	      (silminState->nSolidCoexist)[index] = 0;
	      (silminState->solidComp)[index][j] = 0.0;
	      if (solids[index].na != 1) for (k=0; k<solids[index].na; k++) (silminState->solidComp)[index+1+k][j] = 0.0;
  	  }
	    correctXforChangeInBulkComp();
        }

        /* fluid moved into the liquid menu */
        silminState->fractionateFlu = fracLiquids;

        switch (minType) {
            case 1: /* by vol */
                silminState->minFluPhi = newMin;
                break;
            case 2: /* by mass */
                silminState->maxF = 1.0 - newMin;
                break;
            case 3: /* global ratio */
                if (newMin > 0.0) silminState->fracOut = newMin;
                for (j = 0; j<2; j++) silminState->fracFluids[j] = newMin;
                break;
            case 4: /* individual ratio */
                if (newMin > 0.0) silminState->fracOut = newMin;
                /* will need to change this if more than one fluid possible */
                for (j = 0; j<2; j++) silminState->fracFluids[j] = newMin;
                break;
            default:
                silminState->minFluPhi = 0.0; silminState->maxF = 1.0;
                silminState->fracOut = 1.0;
                for (j = 0; j<2; j++) silminState->fracFluids[j] = 1.0;
                break;
        }

    }
    return TRUE;

}



int menu_option13(char *filename) {

    /* At the moment can only read in .melts file */
    /* eventually can read in trace data etc. */

    if (!iAmInitialized) initializeMenu();

    if(getInputDataFromFile(filename) == GET_INPUT_ERROR_BAD_FILE) return FALSE;
    else assignProblemStatics();

    return TRUE;

}
