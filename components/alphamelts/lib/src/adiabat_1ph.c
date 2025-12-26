#include <math.h>
#include <string.h>
#include <stdio.h>

#include "silmin.h"    /* SILMIN structures include file    */

#include "adiabat.h"
#include "phmelts.h"
#include "interface.h"

#ifdef DEBUG
#undef DEBUG
#endif

#define REALLOC(x, y) (((x) == NULL) ? malloc(y) : realloc((x), (y)))

// NEED TO DO SOMETHING ABOUT ALLOWING P < 1
// Double check Final T for isenthalpic etc.

#define DELTAP ((modeFlag != ISOBARIC && modeFlag != ISOCHORIC) ? (silminState->dspPinc) : 0.0 )
#define DELTAT ((modeFlag != ISENTHALPIC && modeFlag != ISENTROPIC) ? (silminState->dspTinc) : 0.0 )
#define MAX_P  ((silminState->dspPstop > 0 && silminState->dspPinc > 0.0) ? ((silminState->dspPstop < Pmax) ? silminState->dspPstop : Pmax) : Pmax)
#define MIN_P  ((silminState->dspPstop > 0 && silminState->dspPinc < 0.0) ? \
    ((silminState->dspPstop > Pmin) ? silminState->dspPstop : ((silminState->incLiquids) ? Pmin : 1)) : \
    ((silminState->incLiquids) ? Pmin : 1))
#define MAX_T  ((silminState->dspTstop > 0 && silminState->dspTinc > 0.0) ? ((silminState->dspTstop + 273.15 < Tmax) ? silminState->dspTstop + 273.15 : Tmax) : Tmax)
#define MIN_T  ((silminState->dspTstop > 0 && silminState->dspTinc < 0.0) ? ((silminState->dspTstop + 273.15 > Tmin) ? silminState->dspTstop + 273.15 : Tmin) : Tmin)

#define REC 134
#define APPEND 0

static int ista = 0;
int guessFlag = FALSE, modeFlag, ptpathFlag;
double Pmax, Pmin, Tmax, Tmin;
SilminState *startState = NULL, *oldState = NULL, *states = NULL;

int adiabat_0ph(int equilibriumGuess)
{
    int i, j, failed = FALSE;
    int isentropic = silminState->isentropic, isenthalpic = silminState->isenthalpic,
        isochoric = silminState->isochoric, txtOutput = silminState->txtOutput;

    /* first, need initial guess assemblage.  Three methods are available:
   all liquid starting guess, as in MELTS; and a customized norm calculation,
   pdaNorm, in which case a list of phases must be specified. */

    if(!guessFlag && !equilibriumGuess) {
    startingSolution();
	}
    else if  (!guessFlag && equilibriumGuess > 2) {
        findWetLiquidus();

        silminState->dspTstart = silminState->T - 273.15;
        silminState->dspPstart = silminState->P;

    }
    else if (!guessFlag && equilibriumGuess > 1) {
        while(!liquidus());

    silminState->dspTstart = silminState->T - 273.15;
  	silminState->dspPstart = silminState->P;

    }
    else if(!equilibriumGuess){ /* set liquid components to zero */
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

    silminState->isentropic = FALSE;
    silminState->isenthalpic = FALSE;
    silminState->isochoric = FALSE;
    silminState->txtOutput = TEXT_NONE;

    if(!adiabatFunc()) {
        printf("Initial calculation failed (%f bars, %f C)!\n", silminState->dspPstart, silminState->dspTstart);
        failed = TRUE;
    }

    if (failed) {
        copyStateInfo(silminState, oldState);
        return 0;
    }
    guessFlag = TRUE;

    printf("Initial alphaMELTS calculation at: P %f (bars), T %f (C)\n", (silminState->dspPstart),
        (silminState->dspTstart));
    printPhases();

    copyStateInfo(oldState, silminState);

    silminState->isentropic  = isentropic;
    silminState->isenthalpic = isenthalpic;
    silminState->isochoric   = isochoric;
    silminState->txtOutput = txtOutput;

    if (silminState->P > MAX_P) printf("WARNING: Maximum Pressure exceeded\n");
    if (silminState->P < MIN_P) printf("WARNING: Minimum Pressure exceeded\n");
    if (silminState->T > MAX_T) printf("WARNING: Maximum Temperature exceeded\n");
    if (silminState->T < MIN_T) printf("WARNING: Minimum Temperature exceeded\n");

    return 1;

}

int adiabat_1ph(int equilibriumGuess, int saveAll, int iterMax, int noprint)
{
    char line[REC], *pt_filename;
    static int lineno = 0;
    int i, failed = FALSE, iter = 0, ifail = 0, max_iter = 0, max_fails = 0;
    int isentropic = silminState->isentropic, isenthalpic = silminState->isenthalpic, isochoric = silminState->isochoric;
    double refPress, refTemp;
    FILE *fp;

#ifdef ALPHAMELTS_UPDATE_SYSTEM
    pt_filename = silminState->PTfile;
#endif
    max_iter = (iterMax < 0) ? 0 : iterMax;
    max_fails = (iterMax < 0) ? abs(iterMax) : 0;

    /* first, need initial guess assemblage.  Two methods are available:
   all liquid starting guess, as in MELTS; and a customized norm calculation,
   pdaNorm, in which case a list of phases must be specified.  The list of
   phases must include two pyroxenes. */

    if(!saveAll && (ista != 0)) {
        /* normally, unless coming from Amoeba noprint is FALSE */
        states = reAllocSilminStatePointer(states, (size_t) (ista + 1), (size_t) MAX(noprint, 1));
        ista = MAX(noprint-1, 0);
    }
    /* set noprint < 0 to clean up before exit */
    if (noprint < 0) return 0;

    if(modeFlag == PTFILE) {

        if((fp = fopen(pt_filename,"r")) == NULL) {
            printf("Could not open PTpath file %s.\n", pt_filename);
            return 0;
        }

        if(lineno < 0) lineno = 0;
        for(i = 0; i<=lineno; i++) {
            if (fgetstring(line, REC, fp) == NULL) {
                lineno = -1;
                printf("End of PTfile reached.\n");
                fclose(fp);
                return 0;
            }
        }

        lineno++;
        for (i=0; i<(int) strlen(line); i++)
            if (line[i] == '\t') line[i] = ' ';
        sscanf(line, "%lf %lf", &silminState->dspPstart, &silminState->dspTstart);
        silminState->T += silminState->dspTstart + 273.15;

    }

    if(guessFlag) {

        if((isentropic && (silminState->refEntropy != 0.0)) || (isenthalpic && (silminState->refEnthalpy != 0.0)) ||
            (isochoric && (silminState->refVolume != 0.0))) {

            if (isentropic && (silminState->refEntropy != 0.0)) {
                silminState->isentropic = TRUE;
                correctTforChangeInEntropy();
            }
            else if (isenthalpic && (silminState->refEnthalpy != 0.0)) {
                silminState->isenthalpic = TRUE;
                correctTforChangeInEnthalpy();
            }
            else if (isochoric && (silminState->refVolume != 0.0)) {
                silminState->isochoric = TRUE;
                correctPforChangeInVolume();
            }


        /**********/
        if (silminState->assimilate) doBatchAssimilation();


        /* If fO2 is not 'alternative fO2' will need to turn it off to use the H2O buffer with
        constraints like isentropic. */

    }


    }
    else {

        /* First calculation gives reference entropy unless user assigned.
     Then we step in pressure and/or temperature along the desired path... */

        adiabat_0ph(equilibriumGuess);

    }

    if ((silminState->fractionateSol || silminState->fractionateFlu) && silminState->fracSComp == NULL) {
        silminState->fracSComp  = (double **) calloc((unsigned) npc, sizeof(double *));
        silminState->nFracCoexist = (int *) calloc((unsigned) npc, sizeof(int));
    }
    if (silminState->fractionateLiq && silminState->fracLComp == NULL) {
        silminState->fracLComp = (double *) calloc((unsigned) nlc, sizeof(double));
    }


    /*if((Tmax > 0.0) && getenv("ALPHAMELTS_FRACTIONATE_TARGET") != NULL)
        fractionatingSolids = -1;*/




    if(!adiabatFunc()) {
        printf("Initial calculation failed (%f bars, %f C)!\n", silminState->dspPstart, silminState->dspTstart);
        failed = TRUE;
        ifail++;
    }
    else iter++;

    /* The behaviour that we return if no calculations are successful, even with SKIP_FAILURE,
   predates being able to make multipe option 4 calls. Trying to fix it. */
    if (failed) {
        if (ifail <= max_fails) {
            copyStateInfo(silminState, oldState);
            printf("WARNING: Failed once but continuing anyway...\n");
        }
        else {
            if(modeFlag == PTFILE) {
                lineno--;
                fclose(fp);
            }
            return 0;
        }
    }

    refTemp = silminState->T;
    refPress = silminState->P;

    copyStateInfo(oldState, silminState);

    /* Put this before any fO2 buffer is turned off so that it will be recorded appropriately */
    copyStateInfo(&states[ista++], silminState);
    states = reAllocSilminStatePointer(states,(size_t) ista, (size_t) (ista+1));

    /* May need WRITE, but actually should be OK becuase will WRITE if the static storage is NULL */

    if(!noprint && (silminState->txtOutput > TEXT_TABLE))
        putMultipleDataToFile("majors_tbl.txt", "traces_tbl.txt", states, ista, !saveAll);

    if(!noprint) {
        printf("Initial alphaMELTS calculation at: P %f (bars), T %f (C)\n", (silminState->dspPstart),
        (silminState->dspTstart));
        printPhases();
    }


    if(isentropic || isenthalpic || isochoric) {

        if (isentropic && (silminState->refEntropy == 0.0)) {
            printf("Setting reference Entropy to current entropy = %g\n", silminState->bulkTD.s);
            silminState->refEntropy = silminState->bulkTD.s;
        }
        else if (isenthalpic && (silminState->refEnthalpy == 0.0)) {
            //if(DELTAP != 0.0)	printf("WARNING: isenthalpic path chosen for non-zero DELTAP!\n");
            printf("Setting reference Enthalpy to current enthalpy = %g\n", silminState->bulkTD.h);
            silminState->refEnthalpy = silminState->bulkTD.h;
        }
        else if (isochoric && (silminState->refVolume == 0.0)) {
            printf("Setting reference Volume to current volume = %g\n", silminState->bulkTD.v);
            silminState->refVolume = silminState->bulkTD.v;
        }

                /* If fO2 is not 'alternative fO2' will need to turn it off to use the H2O buffer with
        constraints like isentropic. */

    }

    silminState->isentropic  = isentropic;
    silminState->isenthalpic = isenthalpic;
    silminState->isochoric   = isochoric;




#ifdef NEVER_DEFINED
 if((MAX_T > 0.0) && getenv("ALPHAMELTS_FRACTIONATE_TARGET") != NULL) {

        if (targetNo == 0.0) {
            if (silminState->dspLiquidComp[0][mgoindex] == targetWt) {
	if(!noprint) printf("Target MgO wt%% reached: %f\n", silminState->dspLiquidComp[0][mgoindex]);
	return ista;
            }
            else if (silminState->dspLiquidComp[0][mgoindex] > targetWt) {
	if(!noprint) printf("Forward fractionation...MgO = %f\n",
                silminState->dspLiquidComp[0][mgoindex]);
	fractionatingSolids = 1;
            } else {
	if(!noprint) printf("Back fractionation...MgO = %f\n",
                silminState->dspLiquidComp[0][mgoindex]);
	fractionatingSolids = -1;
            }
        } else {
            if (mgnumber == targetNo) {
	if(!noprint) printf("Target Mg# reached: %f\n", mgnumber);
	return ista;
            }
            else if (mgnumber > targetNo) {
	if(!noprint) printf("Forward fractionation...Mg# = %f\n",mgnumber);
	fractionatingSolids = 1;
            } else {
	if(!noprint) printf("Back fractionation...Mg# = %f\n",mgnumber);
	fractionatingSolids = -1;
            }
        }
    }
 #endif







    if (silminState->P > MAX_P) printf("WARNING: Maximum Pressure exceeded\n");
    if (silminState->P < MIN_P) printf("WARNING: Minimum Pressure exceeded\n");
    if (silminState->T > MAX_T) printf("WARNING: Maximum Temperature exceeded\n");
    if (silminState->T < MIN_T) printf("WARNING: Minimum Temperature exceeded\n");

    while (silminState->P <= MAX_P && silminState->P >= MIN_P && silminState->T <= MAX_T && silminState->T >= MIN_T
        && (max_iter == 0 || iter < max_iter)) {

        /********** POST-ALPHAMELTSFUNC *************/

        /* option to do fractional or continuous melting at evenly spaced increments */
        if (silminState->fractionateLiq || silminState->fractionateFlu) {
            (void) extractMelt();
        }
        /* option to fractionate solids */
        if (silminState->fractionateSol || silminState->fractionateFlu) {
            (void) fractionateSolids();
        }

        /********** TEST WHETHER SHOULD CONTINUE ***************/

        /*if ((DELTAT == 0.0 || modeFlag == ISOTHERMAL ||
    modeFlag == ISENTROPIC ||  modeFlag == ISENTHALPIC) &&
	(DELTAP == 0.0 || modeFlag == ISOBARIC || modeFlag == ISOCHORIC)) {
            if (ista == max_iter) {
	printf("Maximum number of iterations reached.\n");
	break;
            }
            }*/

        if (modeFlag == PTFILE) {
            if (fgetstring(line, REC, fp) == NULL) {
    	lineno = -1;
                printf("End of PTfile reached.\n");
                fclose(fp);
                break;

            } else {

    	lineno++;
                for (i=0; i<(int) strlen(line); i++)
                    if (line[i] == '\t') line[i] = ' ';
                sscanf(line, "%lf %lf", &silminState->P, &silminState->T);
                silminState->T += 273.15;

                if (silminState->P > MAX_P) printf("WARNING: Maximum Pressure exceeded\n");
                if (silminState->P < MIN_P) printf("WARNING: Minimum Pressure exceeded\n");
                if (silminState->T > MAX_T) printf("WARNING: Maximum Temperature exceeded\n");
                if (silminState->T < MIN_T) printf("WARNING: Minimum Temperature exceeded\n");

            }
        }
        else {

            double saveP, saveT;

            saveP = silminState->P; saveT = silminState->T;

        // THIS NEEDS FIXING FOR TPGRID EQUIVALENT
            if (ptpathFlag == PTGRID && ((silminState->T + DELTAT) > MAX_T || (silminState->T + DELTAT) < MIN_T)) silminState->P += DELTAP;
            else if (ptpathFlag != PTGRID && modeFlag != ISOBARIC && modeFlag != ISOCHORIC) silminState->P += DELTAP;
            if (fabs(silminState->P) < 1.0) silminState->P = 1.0;


// ISOCHORIC etc if Vinc etc non-zero...
// Hstop  etc.


            /*if (getenv("ALPHAMELTS_FRACTIONATE_TARGET") != NULL) silminState->T -= fabs(DELTAT);
	else*/ if (ptpathFlag == PTGRID && ((silminState->T + DELTAT) > MAX_T || (silminState->T + DELTAT) < MIN_T))
	    silminState->T = refTemp;
            else if (modeFlag == PTPATH || ptpathFlag == PTGRID || modeFlag == ISOBARIC || ptpathFlag == GEOTHERMAL ||
                modeFlag == ISOCHORIC) silminState->T += DELTAT;





#ifdef NEVER_DEFINED
            if (getenv("ALPHAMELTS_FRACTIONATE_TARGET") != NULL) silminState->T -= fabs(DELTAT);
            else if ((modeFlag == PTGRID) && ptlimits) silminState->T = refTemp;
            else if ((modeFlag == TPGRID) && ptlimits) {
                copyStateInfo(silminState, startState);
	silminState->T = saveT + DELTAT;
            }
            else if (modeFlag != TPGRID && modeFlag != ISENTROPIC && modeFlag != ISOTHERMAL &&
                modeFlag != ISENTHALPIC) silminState->T += DELTAT;

            if (getenv("ALPHAMELTS_FRACTIONATE_TARGET") != NULL) {

	if(fractionatingSolids > 0) {
	  if((targetNo != 0.0) && (mgnumber < targetNo)) {
	  if(!noprint) printf("Target Mg# reached: %f\n", mgnumber);
	  saveP = silminState->P; saveT = silminState->T;
	  break;
	  }
	  else if((targetWt != 0.0) && (silminState->dspLiquidComp[0][mgoindex] < targetWt)) {
	  if(!noprint) printf("Target MgO wt%% reached: %f\n", silminState->dspLiquidComp[0][mgoindex]);
	  saveP = silminState->P; saveT = silminState->T;
	  break;
	  }
	}
#endif











            if (silminState->P > MAX_P || silminState->P < MIN_P || silminState->T > MAX_T || silminState->T < MIN_T)  {

            if (silminState->P > MAX_P) printf("Maximum Pressure reached\n");
            if (silminState->P < MIN_P) printf("Minimum Pressure reached\n");
            if (silminState->T > MAX_T) printf("Maximum Temperature reached\n");
            if (silminState->T < MIN_T) printf("Minimum Temperature reached\n");

            silminState->P = saveP; silminState->T = saveT;
            break;

            }

        }

        /************ PRE-ALPHAMELTSFUNC **************/

        failed = FALSE;
        copyStateInfo(oldState, silminState);


        if (silminState->assimilate) doBatchAssimilation();

        if (isentropic) correctTforChangeInEntropy();
        else if (isenthalpic) correctTforChangeInEnthalpy();
        else if (isochoric) correctPforChangeInVolume();


        /************* ALPHAMELTSFUNC ***************/

        if (!failed) {
            /* option to distribute trace elements; loop to converge water kludge */
            if (!adiabatFunc()) {
                printf("Failure in silmin (%f bars, %f K)\n", silminState->P, silminState->T);
                failed = TRUE;
                ifail++;
            }
            else {
                ifail = 0;
                iter++;
            }
        }

        if (failed) {


            if (ifail <= max_fails) {
                copyStateInfo(silminState, oldState);
                printf("WARNING: Failed %s but continuing anyway...\n", (ifail == 1) ? "once" : "again");

                if (silminState->P >= MAX_P) printf("WARNING: Maximum Pressure reached\n");
                if (silminState->P <= MIN_P) printf("WARNING: Minimum Pressure reached\n");
                if (silminState->T >= MAX_T) printf("WARNING: Maximum Temperature reached\n");
                if (silminState->T <= MIN_T) printf("WARNING: Minimum Temperature reached\n");


                // Put failed state in output file - should be obvious as will be identical
                if (saveAll) {

                    double saveP, saveT;

                    saveP = silminState->P; saveT = silminState->T;


                    silminState->P = oldState->P; silminState->T = oldState->T;
                    copyStateInfo(&states[ista++], silminState);
                    silminState->P = saveP; silminState->T = saveT;

                    states = reAllocSilminStatePointer(states,(size_t) ista, (size_t) (ista+1));

                    if(!noprint && (silminState->txtOutput > TEXT_TABLE))
                        putMultipleDataToFile("majors_tbl.txt", "traces_tbl.txt", states, ista, APPEND);
                }
                continue;
            }
            else {
                if(modeFlag == PTFILE) {
                    lineno--;
                    fclose(fp);
                }
                return 0;
            }
  	}

        /********** POST-ALPHAMELTSFUNC *************/

        /* save history of path in array states[] -- note order of operations changed to record
                state before fractionations */
        copyStateInfo(&states[ista++], silminState);
        states = reAllocSilminStatePointer(states,(size_t) ista, (size_t) (ista+1));

        if(!noprint && (silminState->txtOutput > TEXT_TABLE))
            putMultipleDataToFile("majors_tbl.txt", "traces_tbl.txt", states, ista, APPEND);

        if(!noprint) {
            printf("alphaMELTS at: P %f (bars), T %.2f (C)\n", (silminState->P),(silminState->dspTstart));
            printPhases();
        }

        /* Special case: break before fractionations i.e. just after it finds liquidus again. */
        /* This ensures that we have the correct T, bulk composition etc. for any subsequent
     calculations, regardless of which direction things are going in. */

    }

    /* may be kicked out because reached max_iter */
    if ((MAX_T > 0.0) && (modeFlag == PTFILE) &&
            (DELTAT == 0.0) && (DELTAP == 0.0)) {
        if (fgetstring(line, REC, fp) == NULL) {
            lineno = -1;
            printf("End of PTfile reached.\n");
            fclose(fp);
        }
    }
    return ista;
}

#ifdef NEVER_DEFINED

/* function to construct phase diagram boundaries by following
   appearance/disappearance of phases.  Can do minerals or liquid,
   including exsolution boundaries as well as saturation boundaries.
   For liquid, can follow constant melt fraction contours too.
   New twist 1/19/04 -- deal with one-sided phases that have trouble saturating */
int isograd(void)
{
    int i,j;  /* counting indices */
    char phasename[40];
    int iter, phage, tSide, nPhage, upswitch, searchstage; /* metaswitch,*/
    int isentropic = (modeFlag == ISENTROPIC), isenthalpic = (modeFlag == ISENTHALPIC),
        isochoric = (modeFlag == ISOCHORIC);
    int equilibriumGuess, max_iter = 500;
    double bracket, inT, outT, satT, eRef, *silminT, tol = 0.01;
    double F, fTarget = 0.0, fType, affinity, aTarget = 0.0, aType = 0;

    while(TRUE) {
        printf("Phase to track boundary of (by name, lower case) or 'x' to return to menu: ");
        scanf("%s", phasename);

        if (!strcmp(phasename, "x")) {
	return ista;
        }
        else {
            guessFlag = 1;
        }
        if (!strcmp(phasename, "liquid")) {
            phage = npc;
            printf("Type of melt contour to track:\n0. Phi (melt fraction by volume)\n");
            printf("1. F (melt fraction by mass)\n2. aH2O (activity of water in the melt)\n");
            printf("Choose: ");
            scanf("%lf", &fType);

            if(fType == 2) {
	aType = 1;
	fType = 1;
            }
            fType = !fType;
            printf("Type the %s value to set (or < 0.0 for default): ", (aType ? "aH2O" : (fType ? "Phi" : "F")));
            scanf("%lf", &fTarget);

            if (fTarget > 1.0) fTarget = 1.0;
            if (fTarget < 0.0) {
                if(aType) fTarget = 0.0;
                else if(getenv("ALPHAMELTS_CONTINUOUS_MELTING") == NULL) fTarget = 0.0;
                else if(fType) fTarget = MINphi;
                else fTarget = MINf;
            }
            aTarget = fTarget;
            tSide = 1;

            break;

        }
        else {
            //if (!strcmp(phasename, "rhm")) strcpy(phasename, "rhm oxide");
            for(i=0;i<npc;i++) {
	if (!strcmp(phasename, solids[i].label)) {
	  if (!strcmp(phasename,"clinopyroxene")) {
	  printf("How many coexisting phases to track? ");
	  scanf("%d", &nPhage);
	  } else nPhage = 1;
	  phage = i;
	  if((DELTAP == 0.0) && (DELTAT != 0.0))
	  printf("Phase is in on high-pressure (1) or low-pressure (0) side? ");
	  else printf("Phase is in on high-temperature (1) or low temperature (0) side? ");
	  scanf("%d", &tSide);
	  break;
	}
            }
            if (i==npc) printf("Phase not recognized!\n");
            else break;
        }
    }
    if (tSide != 1) tSide = -1;
    printf("Use special monotonic search for troublesome phases (1) or quick search (0)? ");
    scanf("%d", &upswitch);
    if (upswitch != 1) upswitch = 0;

    if(!saveAll) {
        if (!quickWrite) states = reAllocSilminStatePointer(states,(size_t) (ista + 1), (size_t) 1);
        ista=0;
    }

    /* these get switched on again on execution */
    silminState->isentropic = FALSE;
    silminState->isenthalpic = FALSE;
    silminState->isochoric = FALSE;
    if(isentropic) {
        eRef = silminState->refEntropy;
        silminState->refEntropy = 0.0;
    }
    if(isenthalpic) {
        eRef = silminState->refEnthalpy;
        silminState->refEnthalpy = 0.0;
    }
    if(isochoric) {
        eRef = silminState->refVolume;
        silminState->refVolume = 0.0;
    }

    /* optional loop to converge water kludge */
    if (!adiabatFunc()) {
        printf("Initial calculation failed (%f bars, %f K)!\n", silminState->P, silminState->T);
        return 0;
    }

    printf("Initial Guess (not an isograd solution): P %f, T %f\n", (silminState->P),
    (silminState->dspTstart));
    printPhases();

    if(fType) F = silminState->meltFraction;
    else F = silminState->liquidMass/(silminState->liquidMass + silminState->solidMass);

    if(aType) {
        if((F < 0.00001) && (silminState->aH2O < aTarget))
            fTarget = 0.0;
        else if ((F > 0.99999) && (silminState->aH2O > aTarget))
            fTarget = 1.0;
        else {
            F = -silminState->aH2O;
            fTarget = -aTarget;
        }
    }

    if((DELTAP == 0.0) && (DELTAT != 0.0) && (phage != npc)) silminT = &(silminState->P);
    else silminT = &(silminState->T);

    if ((DELTAP == 0.0) && (DELTAT != 0.0) && (phage != npc)) {
        tSide *= 100;
        tol *= 100.0;
    }

    bracket = 2.0;
    while (TRUE) {
        inT = 0.0; outT = 1.0; satT = 0.0;
        iter = 0; searchstage = 0;
        while(TRUE) {

            /* for upswitch mode, find satT -- first get phase in, then get it out */
            if (searchstage == 0 && upswitch == 1 &&
	  ((phage != npc && silminState->nSolidCoexist[phage] < nPhage) || (phage == npc && F < fTarget))) {

	if(phage == npc) *silminT += 5.0*tSide*bracket;
	else *silminT += tSide*bracket;

            } else if (searchstage == 0 && upswitch == 1) {

                satT = *silminT;
                inT = *silminT;
                searchstage = 1;
                *silminT -= tSide*bracket;

            } else if (searchstage == 2 && upswitch == 1) {
	/* double check on satT */
	if((phage != npc && silminState->nSolidCoexist[phage] < nPhage) || (phage == npc && F < fTarget)) {
	  searchstage = 0;

	  if(phage == npc) *silminT += 5.0*tSide*bracket;
	  else *silminT += tSide*bracket;
	  inT = 0.0; satT = 0.0;
	}
	else {

	  searchstage = 1;
	  *silminT = (inT + outT)/2.0;

	}
            } else if ((phage != npc && silminState->nSolidCoexist[phage] < nPhage) || (phage == npc && F < fTarget)) {

                outT = *silminT;
                if (upswitch == 0 && inT != 0.0) *silminT = (inT + outT)/2.0;
                else if (upswitch == 0) {
	  if(phage == npc) *silminT += 5.0*tSide*bracket;
	  else *silminT += tSide*bracket;
    }
                else {
	  searchstage = 2;
                    *silminT = satT;
                }

   } else {

                inT = *silminT;
                if (outT != 1.0) *silminT = (inT + outT)/2.0;
                else *silminT -= tSide*bracket;

            }

            /* optional loop to converge water kludge */
            if (!adiabatFunc()) {
	printf("Failure in silmin (%f bars, %f K)\n", silminState->P, silminState->T);
	return 0;
            }
            affinity = isoFunc1(phage);

            iter++;

            if(fType) F = silminState->meltFraction;
            else F = silminState->liquidMass/(silminState->liquidMass + silminState->solidMass);

            if(aType) {
	if((F < 0.00001) && (silminState->aH2O < aTarget))
	  fTarget = 0.0;
	else if ((F > 0.99999) && (silminState->aH2O > aTarget))
	  fTarget = 1.0;
	else {
	  F = -silminState->aH2O;
	  fTarget = -aTarget;
	}
            }

            /* convergence criteria */
            if(iter == max_iter) {
	  printf("Failed to find isograd solution!\n");
	  if ((DELTAP == 0.0) && (DELTAT != 0.0) && (phage != npc))
	  printf("Please try twiddling starting pressure.\n");
	  else
	  printf("Please try twiddling starting temperature.\n");
	  return 0;
            }
            if (phage==npc) {
                if (fTarget == 0.0) {
	  if((affinity != 0.0) && (affinity < 0.0001) && (silminState->liquidMass == 0.0)) break;
	  if((affinity == 0.0) && (fabs(inT-outT) < tol) && (silminState->liquidMass == 0.0)) break;
	  if ((fabs(F - fTarget) < 0.00001) && (fabs(inT-outT) < tol)) break;
	}
	/* add similar convergence criteria for liquidus to that for solidus */
	else if (fTarget == 1.0) {
	  double minAffinity = isoFunc1(-1);
	  double water = silminState->solidComp[iPhaseH2O][0]*solids[iPhaseH2O].mw;
	  if((minAffinity != 0.0) && (minAffinity < 0.0001) &&
            (silminState->solidMass == 0.0)) break;
	  if((minAffinity == 0.0) && (fabs(inT-outT) < tol) &&
            (silminState->solidMass == 0.0)) break;
	  if ((fabs(F - fTarget) < 0.00001) && (fabs(inT-outT) < tol)) break;
	  /* this could fail if searching in P instead of T? */
	  if((water != 0.0) && (silminState->solidMass - water < MASSIN) && (silminState->T > MAX_T)) {
	  printf("Failed to find isograd solution, due to water saturation!\n");
	  return 0;
	  }
	}
                else if (fabs(F - fTarget) < 0.00001) break;
	else if (fabs(inT-outT) < tol) { /* not totally out to lunch? */
	  if(fabs(F - fTarget) < (aType ? 0.0001 : 0.001)) break;
	  else { /* found a wrong solution, possibly due to hysteresis */
	  outT -= tSide*bracket;
	  inT += 5.0*tSide*bracket;
	  }
	}
            } else { /* looking for a solid */
	if (nPhage == 1 && upswitch == 0) {
	  if ((affinity != 0.0) && (affinity < 0.0001) && silminState->nSolidCoexist[phage] == 0) break;
	  if ((affinity == 0.0) && (fabs(inT-outT) < tol) && silminState->nSolidCoexist[phage] == 0) break;
	}
	else if (fabs(inT-outT) < tol) break;
            }
        }
        if(outT == 1.0) outT = *silminT; /* just in case e.g. for second calling */
        if(inT == 0.0) inT = *silminT;

        /* has now converged so (used to) get rid of phase again and / or make sure not at satT */
        if(phage == npc) {
            /* liquid would be low P side but this should still be O.K. as tSide will be -ve */
            if (fTarget == 1.0) *silminT = outT; //*silminT = inT;
            else if (fTarget == 0.0) *silminT = outT;
            /* for all other F probably best to keep current value? */
        }
        else {
            //*silminT = outT;
            *silminT = inT;
        }

        /* optional loop to converge water kludge */
        if (!adiabatFunc()) {
            printf("Failure in silmin (%f bars, %f K)\n", silminState->P, silminState->T);
            return 0;
        }

        if (doingTraces)
            (void) doTraceElements(SECOND, silminState);


        if((MAX_T > 0.0) || saveAll) {
            putMultipleDataToFile(FALSE, silminState, &states, &ista, APPEND); /* reallocates states, if appropriate */
        }


        printf("Isograd 1 solution at: P %f, T %f\n", (silminState->P), (silminState->dspTstart));
        printPhases();



        if(MAX_T < 0.0) {
            if (silminState->P > MAX_P) printf("WARNING: Maximum Pressure exceeded\n");
            if (silminState->P < MIN_P) printf("WARNING: Minimum Pressure exceeded\n");
            if (silminState->T > -MAX_T) printf("WARNING: Maximum Temperature exceeded\n");
            if (silminState->T < MIN_T) printf("WARNING: Minimum Temperature exceeded\n");
            break;
        }
        else {

            double saveP = silminState->P, saveT = silminState->T;

            silminState->P += DELTAP;
            if (fabs(silminState->P) < 1.0) silminState->P = 1.0;
            if((DELTAP == 0.0) && (DELTAT != 0.0) && (phage != npc)) silminState->T += DELTAT;
            iter = 0;
            if (silminState->P > MAX_P || silminState->P < MIN_P || silminState->T > MAX_T ||
	  silminState->T < MIN_T)  {
                if (silminState->P > MAX_P) printf("Maximum Pressure reached\n");
                if (silminState->P < MIN_P) printf("Minimum Pressure reached\n");
                if (silminState->T > MAX_T) printf("Maximum Temperature reached\n");
                if (silminState->T < MIN_T) printf("Minimum Temperature reached\n");

	silminState->P = saveP; silminState->T = saveT;
                break;
            }
        }

        /* silminState->refEntropy = 0.0; */
        /* if (getenv("ALPHAMELTS_EMERGENCY_FO2") != NULL) fo2Loop(); */
        /* optional loop to converge water kludge */
        if (!adiabatFunc()) {
            printf("Failure in silmin (%f bars, %f K)\n", silminState->P, silminState->T);
            return 0;
        }

        if(fType) F = silminState->meltFraction;
        else F = silminState->liquidMass/(silminState->liquidMass + silminState->solidMass);
    }

    if(isentropic) {
        silminState->isentropic = TRUE;
        silminState->refEntropy = eRef;
    }
    if(isenthalpic) {
        silminState->isenthalpic = TRUE;
        silminState->refEnthalpy = eRef;
    }
    if(isochoric) {
        silminState->isochoric = TRUE;
        silminState->refVolume = eRef;
    }

    return ista;
}

/* New twist 8/28/02 -- just follow affinity=0 contour, suppress
   all phases but liquid */
int isograd2(SilminState *silminState)
{
    int i, j, phage, equilibriumGuess = TRUE;
    char phasename[40];
    int tSide;
    double bracket, inT, outT, affinity, *silminT, tol = 0.01;
    double saveP, saveT;

    while(TRUE) {
        printf("Phase to track boundary of (by name, lower case) or 'x' to return to menu: ");
        scanf("%s", phasename);

        if (!strcmp(phasename, "x")) return ista;
        if (!strncmp(phasename, "liquid", 6)) {
            phage = npc;
            break;
        }
        else {
            if (!strcmp(phasename, "rhm")) strcpy(phasename, "rhm oxide");
            for(i=0;i<npc;i++) {
	if (!strcmp(phasename, solids[i].label)) {
	  phage = i;
	  if((DELTAP == 0.0) && (DELTAT != 0.0))
	  printf("Phase is in on high-pressure (1) or low-pressure (0) side? ");
	  else printf("Phase is in on high-temperature (1) or low temperature (0) side? ");
	  scanf("%d", &tSide);
	  break;
	}
            }
            if (i==npc) printf("Phase not recognized!\n");
            else break;
        }
    }

    if (phage != npc) {

        if (tSide != 1) tSide = -1;
        if((DELTAP == 0.0) && (DELTAT != 0.0)) silminT = &(silminState->P);
        else silminT = &(silminState->T);
        if ((DELTAP == 0.0) && (DELTAT != 0.0)) {
            tSide *= 100;
            tol *= 100.0;
        }
        bracket = 2.0;

    }

    if(!saveAll) {
        if (!quickWrite) states = reAllocSilminStatePointer(states,(size_t) (ista + 1), (size_t) 1);
        ista=0;
    }

    while (TRUE) {

        if (phage == npc) {
            int temp = fractionatingSolids;

            fractionatingSolids = -1;
            equilibriumGuess = silminState->incLiquids;
            silminState->incLiquids = TRUE; // Can't find liquidus without liquid!

            /* failure */
            saveP = silminState->P, saveT = silminState->T;
            if(!adiabatFunc()) {
	printf("WARNING: Failed to find isograd solution!\n");
	silminState->P = saveP; silminState->P = saveT;
            }
            fractionatingSolids = temp;

        }
        else {

            inT = 0.0; outT = 1.0;
            while (inT == 0.0 || outT == 1.0) {
	affinity = isoFunc2(silminState, phage, inT, outT, 0.0);
	if ((phage != npc && affinity*tSide > 0) || (phage == npc && affinity*tSide > 0)) {
	  outT = *silminT;
	  *silminT -= tSide*bracket;
	} else {
	  inT = *silminT;
	  *silminT += tSide*bracket;
	}
            }

            /*was zbrent(satfunc, inT, outT, 0.01);*/

            /* failure */
            saveP = silminState->P, saveT = silminState->T;
            if(isoFunc2(silminState, phage, inT, outT, tol) > 99999.9) {
	printf("WARNING: Failed to find isograd solution!\n");
	silminState->P = saveP; silminState->P = saveT;
            }

        }

        putMultipleDataToFile(FALSE, silminState, &states, &ista, APPEND); /* reallocates states, if appropriate */


        if((MAX_T > 0.0) || saveAll) {
            putMultipleDataToFile(FALSE, silminState, &states, &ista, APPEND); /* reallocates states, if appropriate */
        }


        printf("Isograd 2 solution at: P %f, T %f\n", (silminState->P), (silminState->dspTstart));
        if (phage == npc) printPhases();


        if(MAX_T < 0.0) {
	if (silminState->P > MAX_P) printf("WARNING: Maximum Pressure exceeded\n");
	if (silminState->P < MIN_P) printf("WARNING: Minimum Pressure exceeded\n");
	if (silminState->T > -MAX_T) printf("WARNING: Maximum Temperature exceeded\n");
	if (silminState->T < MIN_T) printf("WARNING: Minimum Temperature exceeded\n");
	break;
        }
        else {

            double saveP = silminState->P, saveT = silminState->T;

            silminState->P += DELTAP;
            if (fabs(silminState->P) < 1.0) silminState->P = 1.0;
            if((DELTAP == 0.0) && (DELTAT != 0.0) && (phage != npc)) silminState->T += DELTAT;
            //iter = 0;
            if (silminState->P > MAX_P || silminState->P < MIN_P || silminState->T > MAX_T ||
	  silminState->T < MIN_T)  {
                if (silminState->P > MAX_P) printf("Maximum Pressure reached\n");
                if (silminState->P < MIN_P) printf("Minimum Pressure reached\n");
                if (silminState->T > MAX_T) printf("Maximum Temperature reached\n");
                if (silminState->T < MIN_T) printf("Minimum Temperature reached\n");

	silminState->P = saveP; silminState->T = saveT;
                break;
            }
        }

    }

    if(!equilibriumGuess){ /* set liquid components to zero */
        silminState->incLiquids = FALSE;
        for (j=0; j<silminState->nLiquidCoexist; j++) {
            for (i=0;i<nlc;i++) silminState->liquidComp[j][i] = 0.0;
            for (i=0;i<nc;i++)  silminState->dspLiquidComp[j][i] = 0.0;
        }
        silminState->nLiquidCoexist = 0;
        silminState->liquidMass = 0.0;
        correctXforChangeInBulkComp();
    }

    return ista;
}



SilminState *reMix(SilminState *silminState, int z) {
    static char filename[300], phasename[40] = "startvalue";
    static int nmelts = 1; /* the number of melts files */
    int i, j, index, ns, flag = (getenv("ALPHAMELTS_ASSIMILATE") != NULL);
    int isentropic = (modeFlag == ISENTROPIC), isenthalpic = (modeFlag == ISENTHALPIC),
        isochoric = (modeFlag == ISOCHORIC);
    double ttemp, ptemp, inmass;
    static double *eTrace = NULL, *eMajor = NULL, Xe, eRef, eO2;
    TraceElements *traceElements;
    SilminState *assimState;

    traceElements = silminState->traceElements;

    /* first time with this bunch of melts files */
    if ((eTrace == (double *) NULL) && (eMajor == (double *) NULL)) {

        if (getenv("ALPHAMELTS_ASSIMILATE") != NULL) {
            /*previously asked if text or binary for (isentropic || isenthalpic || isochoric) */
            printf("Number of assimilant MELTS files (or '0' for binary input file\n");
            printf("or negative value to use separate MELTS files for each phase): ");
            scanf("%d", &nmelts);
            if((nmelts > 1) && (nmelts + z < max_iter)) {
	printf("WARNING: not enough MELTS files (%d) for the requested number of iterations (%d).",
    nmelts, max_iter - z);
	max_iter = nmelts + z;
            }
        }

        /* Note change in order of operations! */
        /* printf("Trace and major (1) or trace only (0): ") */
        if (flag) {
            printf("Mass of enriching agent to be added in grams per cycle");
            if (nmelts > 0) printf("\n(or negative value to use 'Initial Mass:' from file)");
            printf (": ");
            scanf("%lf", &Xe);
        }
        else {
            printf("Mass proportion of enriching agent to be added\n");
            if (nmelts > 0) printf("(or negative value to use 'Initial Mass:' from file)");
            printf (": ");
            scanf("%lf", &Xe);
        }

        if (getenv("ALPHAMELTS_ASSIMILATE") != NULL) {
            if (nmelts > 1) printf("Assimilant filenames (use '?' as wildcard): ");
            else if (nmelts < 0) {
	while(TRUE) {
	  //printf("Note: if assimilating liquid choose this before solid phases.\n");
	  if (doingTraces) printf("Please include bulk trace elements in first assimilant file.\n");
	  printf("First phase to assimilate (by name, lower case; or 'x' to return): ");
	  scanf("%s", phasename);
	  for(i=0;i<npc;i++) {
	  /* should really be strncmp (and same in original) */
	  if (!strcmp(phasename, solids[i].label)) break;
	  }
	  if (!strcmp(phasename,"x")) {
	  freeSilminStatePointer(silminState);
	  return  (SilminState *) GET_INPUT_ERROR_BAD_FILE;
	  }
	  if (i==npc && strcmp(phasename, "liquid")) printf("Phase not recognized!\n");
	  else break;
	}
	index = i;
	printf("First assimilant filename: ");
            }
            else printf("Assimilant filename: ");
        }
        else if ((getenv("ALPHAMELTS_FLUX_MELTING") != NULL) && doingTraces) {
            printf("Enriched composition filename for flux melting: ");
        }
        else if (getenv("ALPHAMELTS_FLUX_MELTING") != NULL) {
            printf("Cannot do flux melting if ALPHAMELTS_DO_TRACE is not set!");
            freeSilminStatePointer(silminState);
            return (SilminState *) NULL;
        }

        scanfilename(filename);

        if (nmelts) { /* MELTS file */

            assimState = allocSilminStatePointer();
            copyStateInfo(assimState, silminState);
            if (nmelts > 1)
	assimState = getMultipleDataFromFile(filename, assimState, z);
            else
	assimState = getInputDataFromFile(filename, assimState);

        }
        else { /* binary file */

            assimState = readSilminStatesFromFile(filename, &nmelts);
            /* reset refMass so assimilant is scaled properly */
            assimState->refMass = assimState->liquidMass + assimState->solidMass;
            nmelts = 0;

        }

        /* Note change in order of operations! */
        /* printf("Trace and major (1) or trace only (0): ") */
        if (assimState == (SilminState *) GET_INPUT_ERROR_BAD_FILE) {
            if (nmelts) printf("Could not open enrichment file.");
            else printf("Could not open restart file.");
            freeSilminStatePointer(silminState);
            return  (SilminState *) GET_INPUT_ERROR_BAD_FILE;
        }
        else if (flag) {
            if (Xe < 0.0) Xe = (nmelts > 0) ? assimState->refMass : 0.0;
        }
        else {
            if (Xe < 0.0) Xe = (nmelts > 0) ? assimState->refMass /
                (assimState->refMass + silminState->liquidMass+silminState->solidMass) : 0.0;
        }

        if((nmelts < 0) && (index != npc)) {

            inmass = assimState->refMass;
            assimState->liquidMass = 0.0;
            assimState->solidMass = inmass;
            assimState->fo2Path = TRUE; // don't try to calculate fO2 for mineral files

            /* from silmin() add phase */
            ns = assimState->nSolidCoexist[index];
            /* ignore given composition */
            if (solids[index].na == 1) {
	(assimState->solidComp)[index][0] = inmass /  solids[index].mw;
	for (j=0; j<nlc; j++) (assimState->liquidComp)[0][j] -=
        (solids[index].solToLiq)[j] * (assimState->solidComp)[index][ns];
	assimState->nSolidCoexist[index] = 1;
            }
            /* onus on user to give something stoichiometric */
            else {
	int     na = solids[index].na;
	double *wt = (double *) calloc((size_t) nc, sizeof(double));
	double molesElmSol[107];
	double *mSol = (double *) calloc((size_t) solids[index].na, sizeof(double));

	/* from preclb.c */
	for (i=0; i<nc; i++) wt[i] = assimState->bulkComp[i]; /* /= bulkSystem[i].mw; */
	/* Convert moles of oxides to moles of elements */
	for (i=0; i<107; i++) {
	  molesElmSol[i] = 0.0;
	  for (j=0; j<nc; j++)
	  if ((bulkSystem[j].oxToElm)[i] != 0)
	    molesElmSol[i] += ((double) (bulkSystem[j].oxToElm)[i])*wt[j];
	}

	if (!strcmp(solids[index].label,"clinopyroxene") || !strcmp(solids[index].label,"orthopyroxene")) {
	  // use input Fe2O3 for Opx and Cpx as entered
	  double *e = molesElmSol;
	  double *m = mSol;
	  double sumcat, sumchg, fe2, fe3;
	  static const int Na = 11;
	  static const int Mg = 12;
	  static const int Al = 13;
	  static const int Si = 14;
	  static const int Ca = 20;
	  static const int Ti = 22;
	  static const int Cr = 24;
	  static const int Mn = 25;
	  static const int Fe = 26;

	  /* Sum the cations and correct the analysis for silica deficiency */
	  sumcat  = e[Na] +   e[Mg] +   e[Al] +   e[Si] +   e[Ca] +   e[Ti] +   e[Cr] +   e[Mn] + e[Fe];
	  sumchg  = e[Na] + 2.0*e[Mg] + 3.0*e[Al] + 4.0*e[Si] + 2.0*e[Ca] + 4.0*e[Ti] + 3.0*e[Cr] + 2.0*e[Mn];

	  /* Compute the ferric/ferrous ratio */
	  fe3 = 3.0*sumcat - sumchg - 2.0*e[Fe];
	  fe2 = e[Fe] - fe3;
	  if (fe3 < 0.0) { fe3 = 0.0; fe2 = e[Fe]; }
	  if (fe2 < 0.0) { fe2 = 0.0; fe3 = e[Fe]; }

	  /* Assign moles of endmembers */
	  m[0] = -fe3/2.0 - fe2 - e[Mn] - e[Al]/2.0 - e[Cr]/2.0 + e[Ca] + e[Na]/2.0 - e[Ti];
	  m[1] =  fe3/4.0 + fe2/2.0 + e[Mn]/2.0 + e[Al]/4.0 + e[Cr]/4.0 - e[Ca]/2.0 + e[Mg]/2.0 - e[Na]/4.0;
	  m[2] =  fe2 + e[Mn];
	  m[3] = -fe3/2.0 + e[Al]/2.0 + e[Cr]/2.0 - e[Na]/2.0 + e[Ti];
	  m[4] =  fe3/2.0 - e[Al]/2.0 - e[Cr]/2.0 + e[Na]/2.0 + e[Ti];
	  m[5] =  fe3/2.0 + e[Al]/2.0 + e[Cr]/2.0 - e[Na]/2.0 - e[Ti];
	  m[6] =  e[Na];

	} else {
	  (*solids[index].convert)(FIRST, SECOND, assimState->T, assimState->P,
				   molesElmSol, mSol, (double *) NULL, (double *) NULL, (double **) NULL,
				   (double ***) NULL, (double **) NULL, (double ****) NULL);
	}

	(assimState->solidComp)[index][ns] = 0.0; /* moles rather than mole frac */
	for (i=0; i<na; i++) {
	  (assimState->solidComp)[index][ns] += mSol[i];
	  (assimState->solidComp)[index+1+i][ns] = mSol[i]; /* *(assimState->solidComp)[index][ns];*/
	  for (j=0; j<nlc; j++) (assimState->liquidComp)[0][j] -= (solids[index+1+i].solToLiq)[j] * mSol[i];
	}
	free(wt);
	free(mSol);
	/* solution phases (except liquid) are allowed multiple files */
	assimState->nSolidCoexist[index] += 1;

	/* ensures bulk composition etc. are updated */
	(void) extractMelt(assimState, 1.0); /* pass non-zero liquidout */

            }

        }

    }

    /* first time through or multiple melts files */
    if (((eTrace == (double *) NULL) && (eMajor == (double *) NULL)) || (nmelts > 1)) {

        if ((eTrace == (double *) NULL) && (eMajor == (double *) NULL)) {

            /* first time we have already read a melts file */
            if(nTraceElements) eTrace = (double *) dvector(0, nTraceElements);
            if(flag) eMajor = (double *) dvector(0, nc);

            if (nmelts < 0) {

	/* first phase is in assimState, now read the rest */
	while(strcmp(phasename,"x")) {
	  printf("Next phase to assimilate (by name, lower case; 'x' when done): ");
	  scanf("%s", phasename);
	  for(i=0;i<npc;i++) {
	  /* should really be strncmp (and same in original) */
	  if (!strcmp(phasename, solids[i].label)) break;
	  }
	  if ((i==npc) && strcmp(phasename, "liquid") && strcmp(phasename, "x")) {
	  printf("Phase not recognized!\n");
	  }
	  else if(strcmp(phasename, "x")) {

	  SilminState *silminTemp;

	  index = i;

	  printf("Assimilant filename: ");
	  scanfilename(filename);

	  silminTemp = allocSilminStatePointer();
	  copyStateInfo(silminTemp, silminState);
	  silminTemp = getInputDataFromFile(filename, silminTemp);

	  if (silminTemp == (SilminState *) GET_INPUT_ERROR_BAD_FILE) {
	    printf("Could not open enrichment file.");
	    freeSilminStatePointer(silminState);
	    return  (SilminState *) GET_INPUT_ERROR_BAD_FILE;
	  }

	  //for (j=0; j<nc; j++) (assimState->bulkComp)[j] += (silminTemp->bulkComp)[j];
	  //for (j=0; j<nlc; j++) (assimState->liquidComp)[0][j] += (silminTemp->liquidComp)[0][j];

	  /* should correct for preclb bit... */
	  /* may be slightly off if not stoichiometric */
	  if(silminState->fo2Path != FO2_NONE) assimState->oxygen += silminTemp->oxygen;

	  inmass = silminTemp->refMass;
	  assimState->refMass += inmass;

	  if (index != npc) {

	    assimState->fo2Path = TRUE; // don't try to calculate fO2 for mineral files
	    assimState->solidMass += inmass;

	    /* from silmin() add phase */
	    ns = assimState->nSolidCoexist[index];
	    /* ignore given composition */
	    if (solids[index].na == 1) {
    /* If more than one file for pure phases, just sum them */
    (silminTemp->solidComp)[index][0] = inmass /  solids[index].mw;
    for (j=0; j<nlc; j++) (silminTemp->liquidComp)[0][j] -=
					(solids[index].solToLiq)[j] * (silminTemp->solidComp)[index][0];
    if (ns == 0) (assimState->solidComp)[index][0] = 0.0;
    (assimState->solidComp)[index][0] += (silminTemp->solidComp)[index][0];
    assimState->nSolidCoexist[index] = 1;
	    }
	    /* onus on user to give something stoichiometric */
	    else {
    int     na = solids[index].na;
    double *wt = (double *) calloc((size_t) nc, sizeof(double));
    double molesElmSol[107];
    double *mSol = (double *) calloc((size_t) solids[index].na, sizeof(double));

    /* from check_coexisting_solids.c */
    /* Allocate space to store the new compositional data */
    for (i=0; i<=na; i++) {
        (assimState->solidComp)[index+i]
            = (double *) REALLOC((assimState->solidComp)[index+i], (size_t) (ns+1)*sizeof(double));
        (assimState->solidDelta)[index+i]
            = (double *) REALLOC((assimState->solidDelta)[index+i], (size_t) (ns+1)*sizeof(double));
    }

    /* from preclb.c */
    for (i=0; i<nc; i++) wt[i] = silminTemp->bulkComp[i]; /*/= bulkSystem[i].mw;*/
    /* Convert moles of oxides to moles of elements */
    for (i=0; i<107; i++) {
        molesElmSol[i] = 0.0;
        for (j=0; j<nc; j++)
            if ((bulkSystem[j].oxToElm)[i] != 0)
                molesElmSol[i] += ((double) (bulkSystem[j].oxToElm)[i])*wt[j];
    }

    if (!strcmp(solids[index].label,"orthopyroxene")) { // use input Fe2O3 for Opx as entered

        double *e = molesElmSol;
        double *m = mSol;
        double sumcat, sumchg, fe2, fe3;
        static const int Na = 11;
        static const int Mg = 12;
        static const int Al = 13;
        static const int Si = 14;
        static const int Ca = 20;
        static const int Ti = 22;
        static const int Cr = 24;
        static const int Mn = 25;
        static const int Fe = 26;

        /* Sum the cations and correct the analysis for silica deficiency */
        sumcat  = e[Na] +   e[Mg] +   e[Al] +   e[Si] +   e[Ca] +   e[Ti] +   e[Cr] +   e[Mn] + e[Fe];
        sumchg  = e[Na] + 2.0*e[Mg] + 3.0*e[Al] + 4.0*e[Si] + 2.0*e[Ca] + 4.0*e[Ti] + 3.0*e[Cr] + 2.0*e[Mn];

        /* Compute the ferric/ferrous ratio */
        fe3 = 3.0*sumcat - sumchg - 2.0*e[Fe];
        fe2 = e[Fe] - fe3;
        if (fe3 < 0.0) { fe3 = 0.0; fe2 = e[Fe]; }
        if (fe2 < 0.0) { fe2 = 0.0; fe3 = e[Fe]; }

        /* Assign moles of endmembers */
        m[0] = -fe3/2.0 - fe2 - e[Mn] - e[Al]/2.0 - e[Cr]/2.0 + e[Ca] + e[Na]/2.0 - e[Ti];
        m[1] =  fe3/4.0 + fe2/2.0 + e[Mn]/2.0 + e[Al]/4.0 + e[Cr]/4.0 - e[Ca]/2.0 + e[Mg]/2.0 - e[Na]/4.0;
        m[2] =  fe2 + e[Mn];
        m[3] = -fe3/2.0 + e[Al]/2.0 + e[Cr]/2.0 - e[Na]/2.0 + e[Ti];
        m[4] =  fe3/2.0 - e[Al]/2.0 - e[Cr]/2.0 + e[Na]/2.0 + e[Ti];
        m[5] =  fe3/2.0 + e[Al]/2.0 + e[Cr]/2.0 - e[Na]/2.0 - e[Ti];
        m[6] =  e[Na];

    } else {

        (*solids[index].convert)(FIRST, SECOND, silminTemp->T, silminTemp->P,
                molesElmSol, mSol, (double *) NULL, (double *) NULL, (double **) NULL,
                (double ***) NULL, (double **) NULL, (double ****) NULL);

    }

    (silminTemp->solidComp)[index][0] = 0.0; /* moles rather than mole frac */
    for (i=0; i<na; i++) {
        (silminTemp->solidComp)[index][0] += mSol[i];
        (silminTemp->solidComp)[index+1+i][0] = mSol[i]; /* *(silminTemp->solidComp)[index][ns];*/
        for (j=0; j<nlc; j++) (silminTemp->liquidComp)[0][j] -= (solids[index+1+i].solToLiq)[j] * mSol[i];
    }
    for (i=0; i<=na; i++) (assimState->solidComp)[index+i][ns] = (silminTemp->solidComp)[index+i][0];
    free(wt);
    free(mSol);
    assimState->nSolidCoexist[index] += 1;
	    }

	    /* ensures bulk composition etc. are updated */
	    (void) extractMelt(silminTemp, 1.0); /* pass non-zero liquidout */

	  }
	  else {
	    /* If more than one file for liquid, just sum them */
	    assimState->liquidMass += inmass;
	    if (assimState->nLiquidCoexist == 0)
    for (j=0; j<nlc; j++) (assimState->liquidComp)[0][j] = 0.0;
	    for (j=0; j<nlc; j++) (assimState->liquidComp)[0][j] += (silminTemp->liquidComp)[0][j];
	    assimState->nLiquidCoexist = 1;
	  }

	  for (j=0; j<nc; j++) (assimState->bulkComp)[j] += (silminTemp->bulkComp)[j];
	  freeSilminStatePointer(silminTemp);
	  }
	}
            }
        }
        else if (nmelts > 1) {
            /* read a new melts file */
            assimState = allocSilminStatePointer();
            copyStateInfo(assimState, silminState);
            assimState = getMultipleDataFromFile(filename, assimState, z);

            if (assimState == (SilminState *) GET_INPUT_ERROR_BAD_FILE) {
	printf("Could not open enrichment file.");
	freeSilminStatePointer(silminState);
	return  (SilminState *) GET_INPUT_ERROR_BAD_FILE;
            }

        }

        if(nTraceElements) for (i=0;i<nTraceElements; i++)
        eTrace[i] = assimState->traceElements[i].bulkComp;

        if(flag) {

            for (i=0;i<nc;i++)
	eMajor[i] = (Xe / assimState->refMass) * assimState->bulkComp[i];

            if ((getenv("ALPHAMELTS_ASSIMILATE") != NULL) &&
	  (isentropic || isenthalpic || isochoric)) {

	if ((isentropic && (silminState->refEntropy == 0.0)) ||
	  (isenthalpic && (silminState->refEnthalpy == 0.0)) ||
	  (isochoric && (silminState->refVolume == 0.0))) {

	  if (!guessFlag) silminState->incLiquids = TRUE; /* superliquidus start */
	  if(!adiabatFunc(silminState)) {
	  if(isentropic) printf("Initial entropy calculation failed (%f bars, %f K)!\n", silminState->P, silminState->T);
	  else if (isenthalpic) printf("Initial enthalpy calculation failed (%f bars, %f K)!\n", silminState->P, silminState->T);
	  else if (isochoric) printf("Initial volume calculation failed (%f bars, %f K)!\n", silminState->P, silminState->T);
	  printf ("Try using a subsolidus start (option 3) and\nimposing the reference ");
	  if(isentropic) printf("entropy");
	  else if (isenthalpic) printf("enthalpy");
	  else if (isochoric) printf("volume");
	  printf(" (option 7) before execution.\n");
	  freeSilminStatePointer(silminState);
	  return (SilminState *) NULL;
	  }
	  else {
	  if (isentropic) {
	    printf("Setting reference Entropy to current entropy = %g\n", silminState->bulkTD.s);
	    silminState->refEntropy = silminState->bulkTD.s;
	    silminState->isentropic = TRUE;
	  }
	  else if (isenthalpic) {
	    if(DELTAP != 0.0)	printf("WARNING: isenthalpic path chosen for non-zero DELTAP!\n");
	    printf("Setting reference Enthalpy to current enthalpy = %g\n", silminState->bulkTD.h);
	    silminState->refEnthalpy = silminState->bulkTD.h;
	    silminState->isenthalpic = TRUE;
	  }
	  else if (isochoric) {
	    printf("Setting reference Volume to current volume = %g\n", silminState->bulkTD.v);
	    silminState->refVolume = silminState->bulkTD.v;
	    silminState->isochoric = TRUE;
	  }

	  }
	}

	ttemp = assimState->T;
	ptemp = assimState->P;
	if(nmelts > 0) { /* i.e. we don't have a starting solution */

	  assimState->refEntropy = 0.0;
	  assimState->refEnthalpy = 0.0;
	  assimState->refVolume = 0.0;

	  assimState->isenthalpic = FALSE;
	  assimState->isenthalpic = FALSE;
	  assimState->isochoric = FALSE;

	  assimState->incLiquids = TRUE; /* superliquidus start */
	  assimState->fo2Path = silminState->fo2Path;
	  assimState->fo2Delta = silminState->fo2Delta;

	}
	if (nmelts) {
	  /* if no starting solution then use current PT for estimate */
	  assimState->P = silminState->P;
	  assimState->T = silminState->T;
	}

	if (isentropic || isenthalpic) {
	  if (silminState->P != assimState->P) {
	  printf("Assimilant P must be the same as current P!\n");
	  return 0;
	  }
	}
	else if (isochoric) {
	  if (silminState->T != assimState->T) {
	  printf("Assimilant T must be the same as current T!\n");
	  return 0;
	  }
	}

	if((nmelts > 0)  && !adiabatFunc(assimState)) {
	  if(isentropic) printf("Assimilant entropy calculation failed (%f bars, %f K)!\n", silminState->P, silminState->T);
	  else if (isenthalpic) printf("Assimilant enthalpy calculation failed (%f bars, %f K)!\n", silminState->P, silminState->T);
	  else if (isochoric) printf("Assimilant volume calculation failed (%f bars, %f K)!\n", silminState->P, silminState->T);
	  freeSilminStatePointer(silminState);
	  return (SilminState *) NULL;
	}

	/* this may be incorrect sometimes(?) but we don't output it anyway */
	if(silminState->fo2Path != FO2_NONE)
	  eO2 = assimState->oxygen *(Xe / assimState->refMass);

	/* we have a starting successful solution */
	/* recalculate thermodynamic properties for new (probably lower) T */
	assimState->P = ptemp;
	assimState->T = ttemp;

	updateThermoProp(assimState);
	getSystemProperties(assimState, TRUE); /* check O.K. if no fO2 path */

	if (nmelts < 0)
	  /* reset refMass in case some material 'lost' by non-stoichiometry */
	  assimState->refMass = assimState->liquidMass + assimState->solidMass;

	if(isentropic) eRef = assimState->bulkTD.s;
	else if(isenthalpic) eRef = assimState->bulkTD.h;
	else if(isochoric) eRef = assimState->bulkTD.v;
	eRef *= (Xe / assimState->refMass);

            }

        }

        if(Xe == 0.0) {
            if(getenv("ALPHAMELTS_ASSIMILATE") != NULL) printf("WARNING: Assimilant mass is 0.0!\n");
            else printf("WARNING: Fluid proportion is 0.0!\n");
        }

        freeSilminStatePointer(assimState);
    }

    if (eMajor != (double *) NULL) {

        for (i=0;i<nTraceElements; i++)
            traceElements[i].bulkComp =
	((silminState->solidMass+silminState->liquidMass)*traceElements[i].bulkComp + Xe*eTrace[i])/
	(silminState->liquidMass+silminState->solidMass+Xe);

        for (i=0;i<nc;i++) {
            silminState->bulkComp[i] += eMajor[i];
            if (silminState->liquidMass != 0.0) {
	for (j=0;j<nlc;j++) silminState->liquidComp[0][j] += eMajor[i]*(bulkSystem[i].oxToLiq)[j];
	silminState->liquidMass += eMajor[i]*bulkSystem[i].mw;
            } else correctXforChangeInBulkComp(silminState);
        }
        silminState->refMass += Xe;
        if(silminState->fo2Path != FO2_NONE && silminState->oxygen != 0.0)
            silminState->oxygen += eO2; /* if zero it will be included later */

        if ((getenv("ALPHAMELTS_ASSIMILATE") != NULL) &&
	(isentropic || isenthalpic || isochoric)) {

            if(isentropic) {
	silminState->refEntropy += eRef;
	correctTforChangeInEntropy(silminState);
            }
            else if(isenthalpic) {
	silminState->refEnthalpy += eRef;
	correctTforChangeInEnthalpy(silminState);
            }
            else if(isochoric) {
	silminState->refVolume += eRef;
	correctPforChangeInVolume(silminState);
            }

            updateThermoProp(silminState);
            (void) getSystemProperties(silminState, TRUE);

        }
    }
    else if (eTrace != (double *) NULL) {
        /* Use Source Mixer equation of traces only */
        for (i=0;i<nTraceElements; i++)
            traceElements[i].bulkComp = (1.0-Xe)*traceElements[i].bulkComp + Xe*eTrace[i];

    }

    if ((nmelts > 1) && (z == nmelts)) { /* last melts file */
        if (eMajor != ( double * ) NULL) free(eMajor);
        if (eTrace != ( double * ) NULL) free(eTrace);
        eTrace = NULL;
        eMajor = NULL;
    }

    return silminState;

}

SilminState *sourceMix(SilminState *silminState) {

    char filename[300];
    int i, j, flag = 0, ns;
    double Xe;
    TraceElements *traceElements;
    SilminState *silminTemp;

    if (silminState == (SilminState *) NULL) {

        printf("MELTS filename: ");
        scanfilename(filename);

        guessFlag = 0;
        silminState = getInputDataFromFile(filename, (SilminState *) NULL);

        if(silminState != (SilminState *) GET_INPUT_ERROR_BAD_FILE) {

            assignProblemStatics(silminState);
            if((doingTraces) && (getenv("ALPHAMELTS_TRACE_INPUT_FILE") != NULL))
	readinTraceInputFile(getenv("ALPHAMELTS_TRACE_INPUT_FILE"),
                silminState->nTraces, silminState->traceElements);
            if(getenv("ALPHAMELTS_SKIP_FAILURE") != NULL) copyStateInfo(oldState, silminState);

        }

    }
    else {

        traceElements = silminState->traceElements;

        printf("MELTS filename for source mixer: ");
        scanfilename(filename);

        silminTemp = allocSilminStatePointer();
        copyStateInfo(silminTemp, silminState);
        silminTemp = getInputDataFromFile(filename, silminTemp);

        if (silminTemp == (SilminState *) GET_INPUT_ERROR_BAD_FILE) {
            printf("Could not open enrichment file.");
            freeSilminStatePointer(silminState);
            return  (SilminState *) GET_INPUT_ERROR_BAD_FILE;
        }

        printf("Mix traces and majors (1) or traces only (0)? ");
        scanf("%d", &flag);

        printf("Mass proportion of enriching agent to be added\n");
        printf("(or negative value to use 'Initial Mass:' from file): ");
        scanf("%lf", &Xe);

        if(flag) {
            if(Xe < 0.0) {

	Xe = (silminTemp->refMass + silminState->liquidMass+silminState->solidMass) /
	  (silminState->liquidMass+silminState->solidMass);

	for (i=0;i<nc;i++)
	  silminState->bulkComp[i] *= Xe;
	silminState->liquidMass *= Xe;
	silminState->solidMass *= Xe;
	silminState->refMass *= Xe;

	Xe = silminTemp->refMass /
	  (silminState->liquidMass+silminState->solidMass);

            }
            /* scale to current mass */
            for (i=0;i<nc;i++)
	silminTemp->bulkComp[i] *= (silminState->liquidMass+silminState->solidMass)
	  / silminTemp->refMass;
        }
        else if (Xe < 0.0) {
            Xe = silminTemp->refMass /
	(silminTemp->refMass + silminState->liquidMass+silminState->solidMass);
        }

        for (i=0;i<nTraceElements; i++)
            traceElements[i].bulkComp = (1.0-Xe)*traceElements[i].bulkComp
	+ Xe*silminTemp->traceElements[i].bulkComp;

        if (flag) { /* majors were involved */
            for (i=0;i<nc;i++)
	silminState->bulkComp[i] =
	  (1.0-Xe)*silminState->bulkComp[i] + Xe*silminTemp->bulkComp[i];
            for (j=0;j<nlc;j++) {
	silminState->liquidComp[0][j] = 0.0;
	for (i=0;i<nc;i++) silminState->liquidComp[0][j] +=
                silminState->bulkComp[i]*(bulkSystem[i].oxToLiq)[j];
            }
            for (i=0;i<npc;i++) {
	for (ns=0;ns<silminState->nSolidCoexist[i];ns++) {
	  silminState->solidComp[i][ns] = 0.0;
	  if (solids[i].na > 1) for (j=0;j<solids[i].na;j++)
	  silminState->solidComp[i+1+j][ns] = 0.0;
	}
	silminState->nSolidCoexist[i] = 0;
            }

            /* need to reset in case any calculations have already been made */
            silminState->nLiquidCoexist = 1;
            silminState->liquidMass = silminState->liquidMass+silminState->solidMass;
            silminState->solidMass = 0.0;
            /* reset refMass */
            silminState->refMass = silminState->liquidMass;

            guessFlag = 0;

            if (modeFlag == ISENTROPIC) {
	if(silminState->refEntropy != 0.0)
	  printf ("Old reference Entropy (S0): %g\n", silminState->refEntropy);
	else
	  printf ("Old reference Entropy (S0): not set yet\n");
	if(silminState->bulkTD.s != 0.0)
	  printf ("Old total Entropy: %g\n", silminState->bulkTD.s);
	else
	  printf ("Old total Entropy: not calculated yet\n");
	if(silminState->refEntropy != 0.0)
	  printf ("WARNING: Unsetting reference Entropy as bulk composition has changed in Source Mixer!\n");
	silminState->refEntropy = 0.0;
	silminState->bulkTD.s = 0.0;
            }
            else if (modeFlag == ISENTHALPIC){
	if(silminState->refEnthalpy != 0.0)
	  printf ("Old reference Enthalpy (H0): %g\n", silminState->refEnthalpy);
	else
	  printf ("Old reference Enthalpy (H0): not set yet\n");
	if(silminState->bulkTD.h != 0.0)
	  printf ("Old total Enthalpy: %g\n", silminState->bulkTD.h);
	else
	  printf ("Old total Enthalpy: not calculated yet\n");
	if(silminState->refEnthalpy != 0.0)
	  printf ("WARNING: Unsetting reference Enthalpy as bulk composition has changed in Source Mixer!\n");
	silminState->refEnthalpy = 0.0;
	silminState->bulkTD.h = 0.0;
            }
            else if (modeFlag == ISOCHORIC) {
	if(silminState->refVolume != 0.0)
	  printf ("Old reference Volume (V0): %g\n", silminState->refVolume);
	else
	  printf ("Old reference Volume (V0): not set yet\n");
	if(silminState->bulkTD.v != 0.0)
	  printf ("Old total Volume: %g\n", silminState->bulkTD.v);
	else
	  printf ("Old total Volume: not calculated yet\n");
	if(silminState->refVolume != 0.0)
	  printf ("WARNING: Unsetting reference Volume as bulk composition has changed in Source Mixer!\n");
	silminState->refVolume = 0.0;
	silminState->bulkTD.v = 0.0;

            }
        } /* majors were involved */

        freeSilminStatePointer(silminTemp);

    }
    return silminState;

}













SilminState *sourceAdd(SilminState *silminState) {

    char filename[300];
    int i, j, flag = 0, ns, nsmax;
    double Xe, reference;
    TraceElements *traceElements;
    SilminState *silminTemp;

    if (silminState != (SilminState *) NULL) {

        printf("Source Adder option:\n0. Keep old state (but scale mass)\n1. Add new state to old state\n");
        printf("2. Replace old state with new state ");
        scanf("%d", &flag);
        if (flag < 0) flag = 0;
        else if (flag > 2) flag = 2;

    }

    if (flag) {

        traceElements = silminState->traceElements;

        printf("Binary input filename for source adder: ");
        scanfilename(filename);

        i = 0;
        silminTemp = readSilminStatesFromFile(filename, &i);

        if (silminTemp == (SilminState *) GET_INPUT_ERROR_BAD_FILE) {
            printf("Could not open restart file.");
            freeSilminStatePointer(silminState);
            return (SilminState *) GET_INPUT_ERROR_BAD_FILE;
        }

        if (flag == 2) printf ("Current reference mass: %g\n", silminTemp->refMass);
        printf ("Current total mass (liquid + solid): %g\n", silminTemp->liquidMass + silminTemp->solidMass);
        if (flag == 1) printf("Type new mass (or negative to use current total): ");
        else printf("Type new mass (or 0 to keep old values, or negative value to use current total): ");
        scanf("%lf", &reference);

        if (reference > 0.0) silminTemp->refMass = reference;
        else if (reference < 0.0) silminTemp->refMass = silminTemp->liquidMass + silminTemp->solidMass;
        else if (flag == 1) silminTemp->refMass = reference; /* can add zero amount */

        Xe = silminTemp->refMass / (silminTemp->liquidMass + silminTemp->solidMass);

        printf("Add majors and traces (1) or majors only (0)? ");
        scanf("%d", &flag);

        if(Xe == 0.0) {
            printf("WARNING: Scaling factor is 0.0! No phases added.\n");
        }
        else {

            if (flag) {
	double Xe2 = Xe;
	Xe2 *= silminTemp->liquidMass + silminTemp->solidMass;
	for (i=0;i<nTraceElements; i++)
	  traceElements[i].bulkComp =
	  ((silminState->solidMass+silminState->liquidMass)*traceElements[i].bulkComp +
            Xe2*silminTemp->traceElements[i].bulkComp)/
	  (silminState->liquidMass+silminState->solidMass+Xe2);
            }

            for (i=0;i<nc;i++)
	silminState->bulkComp[i] += Xe*silminTemp->bulkComp[i];

            nsmax = MAX(silminState->nLiquidCoexist, silminTemp->nLiquidCoexist);
            for (j=0;j<nlc;j++) {
	for (ns=0;ns<nsmax;ns++) {
	  silminState->liquidComp[ns][j] += Xe*silminTemp->liquidComp[ns][j];
	}
            }
            silminState->liquidMass += Xe*silminTemp->liquidMass;
            silminState->nLiquidCoexist = nsmax;

            for (i=0;i<npc;i++) {
	nsmax = MAX(silminState->nSolidCoexist[i], silminTemp->nSolidCoexist[i]);
	for (ns=0;ns<nsmax;ns++) {
	  silminState->solidComp[i][ns] += Xe*silminTemp->solidComp[i][ns];
	  if (solids[i].na > 1) for (j=0;j<solids[i].na;j++)
            silminState->solidComp[i+1+j][ns] += Xe*silminTemp->solidComp[i+1+j][ns];
	}
	silminState->nSolidCoexist[i] = nsmax;
            }
            silminState->solidMass += Xe*silminTemp->solidMass;
            silminState->refMass += Xe*silminTemp->refMass;

            if (doingTraces)
	(void) doTraceElements(SECOND, silminState);

            guessFlag = 1;

            if (modeFlag == ISENTROPIC) {
	if(silminState->refEntropy != 0.0)
	  printf ("Old reference Entropy (S0): %g\n", silminState->refEntropy);
	else
	  printf ("Old reference Entropy (S0): not set yet\n");
	if(silminState->bulkTD.s != 0.0)
	  printf ("Old total Entropy: %g\n", silminState->bulkTD.s);
	else
	  printf ("Old total Entropy: not calculated yet\n");
	if(silminState->refEntropy != 0.0)
	  printf ("WARNING: Unsetting reference Entropy as bulk composition has changed in Source Adder!\n");
	silminState->refEntropy = 0.0;
	silminState->bulkTD.s = 0.0;
            }
            else if (modeFlag == ISENTHALPIC){
	if(silminState->refEnthalpy != 0.0)
	  printf ("Old reference Enthalpy (H0): %g\n", silminState->refEnthalpy);
	else
	  printf ("Old reference Enthalpy (H0): not set yet\n");
	if(silminState->bulkTD.h != 0.0)
	  printf ("Old total Enthalpy: %g\n", silminState->bulkTD.h);
	else
	  printf ("Old total Enthalpy: not calculated yet\n");
	if(silminState->refEnthalpy != 0.0)
	  printf ("WARNING: Unsetting reference Enthalpy as bulk composition has changed in Source Adder!\n");
	silminState->refEnthalpy = 0.0;
	silminState->bulkTD.h = 0.0;
            }
            else if (modeFlag == ISOCHORIC) {
	if(silminState->refVolume != 0.0)
	  printf ("Old reference Volume (V0): %g\n", silminState->refVolume);
	else
	  printf ("Old reference Volume (V0): not set yet\n");
	if(silminState->bulkTD.v != 0.0)
	  printf ("Old total Volume: %g\n", silminState->bulkTD.v);
	else
	  printf ("Old total Volume: not calculated yet\n");
	if(silminState->refVolume != 0.0)
	  printf ("WARNING: Unsetting reference Volume as bulk composition has changed in Source Adder!\n");
	silminState->refVolume = 0.0;
	silminState->bulkTD.v = 0.0;

            }
        }

        freeSilminStatePointer(silminTemp);

    }
    else { /* either picked keep old state or there is no old state */

        if (silminState == (SilminState *) NULL) {

            printf("Restart filename: ");
            scanfilename(filename);

            i = 0;
            guessFlag = 0;
            silminState = readSilminStatesFromFile(filename, &i);

            if(silminState != (SilminState *) GET_INPUT_ERROR_BAD_FILE) {

	assignProblemStatics(silminState);
	guessFlag = 1;
	if((doingTraces) && (getenv("ALPHAMELTS_TRACE_INPUT_FILE") != NULL))
	  readinTraceInputFile(getenv("ALPHAMELTS_TRACE_INPUT_FILE"),
                    silminState->nTraces, silminState->traceElements);
	if(getenv("ALPHAMELTS_SKIP_FAILURE") != NULL) copyStateInfo(oldState, silminState);
            }

        }

        printf ("Current reference mass: %g\n", silminState->refMass);
        printf ("Current total mass (liquid + solid): %g\n", silminState->liquidMass + silminState->solidMass);
        printf("Type new mass (or 0 to keep old values, or negative value to use current total): ");

        scanf("%lf", &reference);
        if (reference > 0.0) silminState->refMass = reference;
        else if (reference < 0.0) silminState->refMass = silminState->liquidMass + silminState->solidMass;

        if (reference != 0.0) {

            Xe = silminState->refMass / (silminState->liquidMass + silminState->solidMass);

            for (i=0;i<nc;i++) silminState->bulkComp[i] *= Xe;

            for (j=0;j<nlc;j++)
	for (ns=0;ns<silminState->nLiquidCoexist;ns++)
	  silminState->liquidComp[ns][j] *= Xe;
            silminState->liquidMass *= Xe;

            for (i=0;i<npc;i++) {
	for (ns=0;ns<silminState->nSolidCoexist[i];ns++) {
	  silminState->solidComp[i][ns] *= Xe;
	  if (solids[i].na > 1) for (j=0;j<solids[i].na;j++) silminState->solidComp[i+1+j][ns] *= Xe;
	}
            }
            silminState->solidMass *= Xe;

        }

    }

    /* not strictly necessary but shouldn't do any harm */
    updateThermoProp(silminState);
    (void) getSystemProperties(silminState, TRUE);

    return silminState;

}











#endif


#undef REALLOC
