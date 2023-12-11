#include <math.h>
#include <string.h>
#include <stdio.h>

#include "silmin.h"                /* SILMIN structures include file        */
#include "alphamelts.h"

#include "adiabat.h"
#include "phmelts.h"

#ifdef DEBUG
#undef DEBUG
#endif

#define REALLOC(x, y) (((x) == NULL) ? malloc(y) : realloc((x), (y)))

int modeFlag, ptpathFlag;
double Pmax, Pmin, Tmax, Tmin;

#define DATESTAMP " (" __DATE__ " " __TIME__ ")"
double thisversion = 2.0201; char *filestatus = DATESTAMP;

void splashScreen(double newestversion) {

    int i, j, k;
    char vstring[10], nstring[10];

    if(newestversion != 0.0) newestversion += 1.0e-8;
    thisversion += 1.0e-8; /* nudge to avoid rounding error on certain doubles */

    i = ( int ) newestversion;
    j = ( int ) floor(100.0*newestversion) - 100*i;
    k = ( int ) floor(10000.0*newestversion) - 10000*i - 100*j;

    ( void ) sprintf(nstring, "%i.%i.%i", i, j, k);

    i = ( int ) thisversion;
    j = ( int ) floor(100.0*thisversion) - 100*i;
    k = ( int ) floor(10000.0*thisversion) - 10000*i - 100*j;
    ( void ) sprintf(vstring, "%i.%i.%i", i, j, k);

    printf("\n*** alphaMELTS %s%s -- ", vstring, filestatus);
    //if(getenv("ALPHAMELTS_DO_TRACE_H2O") != NULL) printf("pHMELTS ");
    if (calculationMode == MODE_pMELTS) printf("pMELTS ");
    else printf("rhyolite-MELTS ");

       if (ptpathFlag == PSEUDOSECTION) printf("pseudosection w/ or w/o liquid ***\n");
    else if (ptpathFlag == GEOTHERMAL)    printf("geotherm w/ or w/o liquid ***\n");
    else if (ptpathFlag == PTGRID) {
         if (modeFlag == ISENTROPIC)  printf("P-S grid w/ or w/o liquid ***\n");
        else if (modeFlag == ISENTHALPIC) printf("P-H grid w/ or w/o liquid ***\n");
        else if (modeFlag == ISOCHORIC)   printf("V-T grid w/ or w/o liquid ***\n");
        else                              printf("P-T grid w/ or w/o liquid ***\n");
    }
    else if (modeFlag == ISENTROPIC)  printf("P-S path w/ or w/o liquid ***\n");
    else if (modeFlag == ISOBARIC)    printf("isobar w/ or w/o liquid ***\n");
    else if (modeFlag == ISOTHERMAL)  printf("isotherm w/ or w/o liquid ***\n");
    else if (modeFlag == ISENTHALPIC) printf("P-H path w/ or w/o liquid ***\n");
    else if (modeFlag == ISOCHORIC)   printf("V-T path w/ or w/o liquid ***\n");
    else                              printf("P-T path w/ or w/o liquid ***\n");

    printf("\nThis front end is the work of Paul Asimow and Paula Antoshechkina (nee Smith) and it\n");
    printf("uses the rhyolite-MELTS and pMELTS algorithms developed by Mark Ghiorso & co-workers.\n");
    printf("You are welcome to use and distribute this program, under the condition that you\n");
    printf("acknowledge all the contributors by citing the appropriate references with any results.\n");
    printf("See Smith & Asimow (2005), documentation and the forum for details.\n\n");

    if (newestversion > thisversion)
        printf("This is version %s; updated version %s is available on GitLab and described at\n", vstring, nstring);
    else if (newestversion == 0.0)
        printf("Unable to automatically check for updates; check for updates at\n");
    else
        printf("Links to binaries, documentation and example files are available at\n");

    printf("http://magmasource.caltech.edu/alphamelts/version2.php and\n");
    printf("http://magmasource.caltech.edu/forum/index.php/board,31.0.html\n");

}

int assignGlobalStatics(int mode) {

    /* If successful, this runs InitComputeDataStruct() etc. */
    if (setCalculationMode(mode)) {
        if      (mode == MODE__MELTS)           printf("---> Calculation mode is rhyolite-MELTS (public release v 1.0.2).\n");
        else if (mode == MODE__MELTSandCO2)     printf("---> Calculation mode is rhyolite-MELTS (public release v 1.1.0).\n");
        else if (mode == MODE__MELTSandCO2_H2O) printf("---> Calculation mode is rhyolite-MELTS (public release v 1.2.0).\n");
        else if (mode == MODE_pMELTS)           printf("---> Calculation mode is pMELTS (public release v 5.6.1).\n");
        //if (outsw == (int *) NULL) setoutput();
    }
    else {
        printf("Unable to set calculation mode!\n");
        return FALSE;
    }
    return TRUE;

}

int getEnvironmentSettings(void) {

    modeFlag = 0;

    /* get parameters from environment variables, if set */
    /* need to check for buffer overflow */
    /* need to check the values are sensible too */

    if (getenv("ALPHAMELTS_MAXP") != NULL) Pmax = atof(getenv("ALPHAMELTS_MAXP"));
    else if (calculationMode == MODE_pMELTS) Pmax = 40000.0;
    else Pmax = 30000;
    if(Pmax < 1.0) Pmax = 1.0;

    if (getenv("ALPHAMELTS_MAXT") != NULL) Tmax = atof(getenv("ALPHAMELTS_MAXT"));
    else if (calculationMode == MODE_pMELTS) Tmax = 2500.0;
    else Tmax = 2000.0;
    if (Tmax < 0.0) Tmax = 0.0;
    Tmax += 273.15;

    if (getenv("ALPHAMELTS_MINP") != NULL) Pmin = atof(getenv("ALPHAMELTS_MINP"));
    //else if (calculationMode == MODE_pMELTS) Pmin = 10000.0;
    else Pmin = 1.0;
    if(Pmin < 1.0) Pmin = 1.0;

    if (getenv("ALPHAMELTS_MINT") != NULL) Tmin = atof(getenv("ALPHAMELTS_MINT"));
    else if (calculationMode == MODE_pMELTS) Tmin = 1000.0;
    else Tmin = 500.0;
    if (Tmin < 0.0) Tmin = 0.0;
    Tmin += 273.15;

    modeFlag = PTPATH; /* dspPinc (or dspVinc) and dspTinc (or dspHinc, dspSinc) are respected */
    ptpathFlag = PTPATH;

    /* For display purposes */
    if (getenv("ALPHAMELTS_RUN_MODE") != NULL) {
        char line[20];
        int i, len = MIN(strlen(getenv("ALPHAMELTS_RUN_MODE")), strlen("isenthalpic"));

        ( void ) strncpy(line, getenv("ALPHAMELTS_RUN_MODE"), len);
        for (i=0; i<len; i++) line[i] = tolower(line[i]);

         if (!strncmp(line, "isenthalpic", 11)) modeFlag = ISENTHALPIC;
        else if (!strncmp(line, "isentropic", 10))  modeFlag = ISENTROPIC;
        else if (!strncmp(line, "isochoric", 9))   modeFlag = ISOCHORIC;
        else if (!strncmp(line, "isothermal", 10)) modeFlag = ISOTHERMAL; // Tinc ignored
        else if (!strncmp(line, "isobaric", 8))   modeFlag = ISOBARIC; // Pinc ignored
        else if (!strncmp(line, "ptpath", 6)) modeFlag = PTPATH; // Pinc and Tinc observed
    }
    if (getenv("ALPHAMELTS_PTPATH_MODE") != NULL) {
        char line[20];
        int i, len = MIN(strlen(getenv("ALPHAMELTS_PTPATH_MODE")), strlen("pseudosection"));

        ( void ) strncpy(line, getenv("ALPHAMELTS_PTPATH_MODE"), len);
        for (i=0; i<len; i++) line[i] = tolower(line[i]);

         if (!strncmp(line, "pseudosection", 13)) ptpathFlag = PSEUDOSECTION;
        else if (!strncmp(line, "geothermal", 10))    ptpathFlag = GEOTHERMAL; // Use gradient instead of inc
        else if (!strncmp(line, "grid", 4))           ptpathFlag = PTGRID;
        else if (!strncmp(line, "file", 4))           ptpathFlag = PTFILE;
        else if (!strncmp(line, "path", 4))           ptpathFlag = PTPATH;
    }

    // return value not actually used
    return modeFlag;
}

void assignProblemStatics(void) {

    /* for display and backwards compatability */
    if (modeFlag == ISENTHALPIC)     { silminState->isenthalpic = TRUE; modeFlag = PTPATH; }
    else if (modeFlag == ISENTROPIC) { silminState->isentropic  = TRUE; modeFlag = PTPATH; }
    else if (modeFlag == ISOCHORIC)  { silminState->isochoric   = TRUE; modeFlag = PTPATH; }

    silminState->T = silminState->dspTstart; /* output in C */
    if (silminState->T != 0.0) silminState->T += 273.15;
    if (silminState->dspPstart == 0.0) silminState->P = 1.0;
    else silminState->P = silminState->dspPstart;

    /* alphaMELTS uses signed increments */
    if ((silminState->T != 0.0) && (silminState->dspTstop != 0))
        silminState->dspTinc = (silminState->dspTstart < silminState->dspTstop) ? fabs(silminState->dspTinc) : -fabs(silminState->dspTinc);
    if ((silminState->P != 0.0) && (silminState->dspPstop != 0))
        silminState->dspPinc = (silminState->dspPstart < silminState->dspPstop) ? fabs(silminState->dspPinc) : -fabs(silminState->dspPinc);

    if (silminState->isenthalpic && (silminState->refEnthalpy != 0.0) && (silminState->dspHstop != 0.0))
        silminState->dspHinc = (silminState->refEnthalpy < silminState->dspHstop) ? fabs(silminState->dspHinc) : -fabs(silminState->dspHinc);
    else if (silminState->isentropic && (silminState->refEntropy != 0.0) && (silminState->dspSstop != 0.0))
        silminState->dspSinc = (silminState->refEntropy < silminState->dspSstop) ? fabs(silminState->dspSinc) : -fabs(silminState->dspSinc);
    else if (silminState->isochoric && (silminState->refVolume != 0.0) && (silminState->dspVstop != 0.0))
        silminState->dspVinc = (silminState->refVolume < silminState->dspVstop) ? fabs(silminState->dspVinc) : -fabs(silminState->dspVinc);


    /* THIS NEEDS FIXING ***********/
    if (ptpathFlag == GEOTHERMAL) {
        /* see create_tp_padb.c */
        if (silminState->isochoric) {
            double vInc = silminState->dspVinc, dtdv = silminState->dspDTDV;
            if (dtdv != 0.0) silminState->dspVinc = dtdv;
            else silminState->dspTinc = vInc*dtdv;
            printf("NOTE: resetting Tinc = %f C\n", silminState->dspTinc);
        }
        else {
            if (silminState->isenthalpic) {
  	    double hInc = silminState->dspHinc, dpdh = silminState->dspDPDH;
                if (dpdh == 0.0) silminState->dspHinc = dpdh;
    	  else silminState->dspPinc = hInc*dpdh;
            }
            else if (silminState->isentropic) {
	      double sInc = silminState->dspSinc, dpds = silminState->dspDPDS;
      	if (dpds == 0.0) silminState->dspSinc = dpds;
      	else silminState->dspPinc = sInc*dpds;
            }
            else {
                double tInc = silminState->dspTinc, dpdt = silminState->dspDPDt;
                if (dpdt == 0.0) silminState->dspTinc = dpdt;
                else silminState->dspPinc = tInc*dpdt;
            }
            printf("NOTE: resetting Pinc = %f bars\n", silminState->dspPinc);
        }

    }

}

void startingSolution() {

    int i,j,k, hasLiquid = (silminState->liquidMass != 0.0);    /* counting indices */
    //char phasename[40] = "startvalue";

    if (!pdaNorm()) {
        printf("Error in initial guess routine -- using superliquidus start\n");
        /*equilibriumGuess = TRUE;*/
        silminState->incLiquids = 1;
        for (i=0;i<npc;i++) {
            for (j=0;j<silminState->nSolidCoexist[i];j++) {
                silminState->solidComp[i][j] = 0.0;
                if (solids[i].na > 1) {
                    for (k=0;k<solids[i].na;k++) silminState->solidComp[i+1+k][j] = 0.0;
                }
            }
            silminState->nSolidCoexist[i] = 0;
        }

        printf("Unbuffered initial guess at: P %f, T %f\n", (silminState->P), silminState->dspTstart);
        printPhases();
        return;

    }

    /* set liquid components to zero */
    if (hasLiquid) {
        for (j=1; j<silminState->nLiquidCoexist; j++) {
            free((silminState->liquidComp )[j]); (silminState->liquidComp )[j] = NULL;
            free((silminState->liquidDelta)[j]); (silminState->liquidDelta)[j] = NULL;
        }
        for (i=0;i<nlc;i++) silminState->liquidComp[0][i] = 0.0;
        //for (i=0;i<nc;i++)  silminState->dspLiquidComp[i] = 0.0;
        silminState->nLiquidCoexist = 1; // else memory leak
        silminState->liquidMass = 0.0;
        hasLiquid = FALSE;
    }

    printf("Norm calculated initial guess at: P %f, T %f\n", (silminState->P), silminState->dspTstart);
    printPhases();

    if((silminState->fo2Path != FO2_NONE) && !silminState->fo2Liq) {

        int fo2path = silminState->fo2Path, isentropic = silminState->isentropic,
            isenthalpic = silminState->isenthalpic, isochoric = silminState->isochoric;
        double fo2delta = silminState->fo2Delta;

        silminState->fo2Path = FO2_NONE;
        silminState->fo2Delta = 0.0;
        if(isentropic) silminState->isentropic = FALSE; /* this was in case ALPHAMELTS_IMPOSE_FO2 */
        else if (isenthalpic) silminState->isenthalpic = FALSE;
        else if (isochoric) silminState->isochoric = FALSE;

        if(!silmin()) {
            printf("Error in initial guess routine -- using superliquidus start\n");

            silminState->incLiquids = 1;
            for (i=0;i<npc;i++) {
                for (j=0;j<silminState->nSolidCoexist[i];j++) {
                    silminState->solidComp[i][j] = 0.0;
                    if (solids[i].na > 1) {
                        for (k=0;k<solids[i].na;k++) silminState->solidComp[i+1+k][j] = 0.0;
                    }
                }
                silminState->nSolidCoexist[i] = 0;
            }
                    silminState->solidMass = 0;
            for (i=0; i<nlc; i++) {
                (silminState->liquidComp)[0][i] = 0.0;
                for (j=0; j<nc; j++) (silminState->liquidComp)[0][i] +=
                               (silminState->bulkComp)[j]*(bulkSystem[j].oxToLiq)[i];
            }
            /*for (j=0; j<nc; j++)
                (silminState->dspLiquidComp)[j] = (silminState->bulkComp)[j]*bulkSystem[j].mw;*/
            silminState->nLiquidCoexist = 1;
            silminState->liquidMass = silminState->refMass;
            silminState->liquidMass = silminState->solidMass;
            silminState->solidMass = 0;
        }
        silminState->fo2Path = fo2path;
        silminState->fo2Delta = fo2delta;
        if(isentropic) silminState->isentropic = TRUE;
        else if (isenthalpic) silminState->isenthalpic = TRUE;
        else if (isochoric) silminState->isochoric = TRUE;

        printf("Unbuffered initial guess at: P %f, T %f\n", (silminState->P), silminState->dspTstart);
        printPhases();
    }

}

/* print out phase compositions */
void printPhases(void) {
    int i, j, k, width, hasLiquid = (silminState->liquidMass != 0.0);
    double *mm, *r, mass, *dspLiquidComp;
    char *formula;

    if (outsw == (int *) NULL) setoutput();

    if (hasLiquid) {
        for (k=0;k<silminState->nLiquidCoexist;k++) {
            dspLiquidComp = (double *) calloc(nlc, sizeof(double));
            //    for (i=0,mass=0.0;i<nlc;i++) mass += silminState->liquidComp[k][i]*liquid[i].mw;
            for (i=0,mass=0.0;i<nc;i++) {
                for (j=0, dspLiquidComp[i]=0.0; j<nlc; j++) {
                    dspLiquidComp[i] += silminState->liquidComp[k][j]*(liquid[j].liqToOx)[i]*bulkSystem[i].mw;
                }
                mass += dspLiquidComp[i];
            }
            for (i=0; i<nc; i++) if (outsw[i]) dspLiquidComp[i] *= ((mass != 0.0) ? 100.0/mass : 0);

            /* Just print once for multiple liquids */
            printf("liquid:   ");
            for (i=0; i<nc; i++) if (outsw[i]) {
                width = MAX(outsw[i], ((dspLiquidComp[i] >= 10.0) ?  floor(log10(dspLiquidComp[i])) + 4 : 4));
                printf("%*.*s ", width, width, bulkSystem[i].label); /* SO3, Cl2O-1, F2O-1 */
            }
            width = 5 - ((mass >= 10.0) ? floor(log10(mass)) : 0);
            printf("\n%7.*f g ", width, mass);
            //for (i=0;i<nc;i++) if (outsw[i]) printf("%.2f ", silminState->dspLiquidComp[k][i]);
            for (i=0; i<nc; i++) if (outsw[i]) {
                width = MAX(outsw[i], ((dspLiquidComp[i] >= 10.0) ?  floor(log10(dspLiquidComp[i])) + 4 : 4));
                printf("%*.2f ", width, dspLiquidComp[i]); /* SO3, Cl2O-1, F2O-1 */
            }
            printf("\n");
            printf("Activity of H2O = %g  Melt fraction = %g\n", silminState->aH2O,
            silminState->liquidMass/(silminState->liquidMass + silminState->solidMass));
            free(dspLiquidComp);
        }
    }

    for (i=0;i<npc;i++) {
        for (j=0;j<(silminState->nSolidCoexist)[i];j++) {
            if (solids[i].na == 1)
                printf("%s: %f g, composition %s\n", solids[i].label,
                    (silminState->solidComp)[i][j]*solids[i].mw, solids[i].formula);
            else {
                mm = (double *) calloc(solids[i].na, sizeof(double));
                r = (double *) calloc(solids[i].na, sizeof(double));
                for (k=0, mass=0.0; k<solids[i].na; k++) {
                    mm[k] = (silminState->solidComp)[i+1+k][j];
                    mass += (silminState->solidComp)[i+1+k][j]*solids[i+1+k].mw;
                }
                (*solids[i].convert)(SECOND, THIRD, silminState->T, silminState->P,
                    (double *) NULL, mm, r, (double *) NULL, (double **) NULL,
                    (double ***) NULL, (double **) NULL, (double ****) NULL);
                (*solids[i].display)(FIRST, silminState->T, silminState->P,
                    r, &formula);
                for (k=0; strlen(formula); k++) {
                    if (formula[k] != ' ') break;
                }
                printf("%s: %f g, composition %s\n", solids[i].label, mass, &formula[k]);
                free(formula);
                if (!(*solids[i].test)(SIXTH, silminState->T, silminState->P, 0, 0, (char **) NULL,
                        (char **) NULL, (double *) NULL, mm))
                        printf("Phase %s is infeasible\n", solids[i].label);
                free(mm); free(r);
            }
        }
    }

}




#undef REALLOC
