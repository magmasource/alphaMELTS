#include <math.h>
#include <string.h>
#include <stdio.h>

#include "silmin.h"                /* SILMIN structures include file        */
#include "alphamelts.h"

//#include "adiabat.h"
#include "phmelts.h"

#ifdef DEBUG
#undef DEBUG
#endif

#define REALLOC(x, y) (((x) == NULL) ? malloc(y) : realloc((x), (y)))

#define DATESTAMP " (" __DATE__ " " __TIME__ ")"
/* double thisversion = 2.0300; */
double thisversion = THIS_VERSION; char *filestatus = DATESTAMP;

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

    printf("\nThis front end is the work of Paula Antoshechkina (nee Smith) and Paul Asimow and it\n");
    printf("uses the rhyolite-MELTS and pMELTS algorithms developed by Mark Ghiorso & co-workers.\n");
    printf("You are welcome to use and distribute this program, under the condition that you\n");
    printf("acknowledge all the contributors by citing the appropriate references with any results.\n");
    printf("See Smith & Asimow (2005), documentation and the GitHub Wiki for details.\n\n");

    if (newestversion > thisversion)
        printf("This is version %s; updated version %s is available at\n", vstring, nstring);
    else if (newestversion == 0.0)
        printf("Unable to automatically check for updates; check for updates at\n");
    else
        printf("Links to binaries, documentation and example files are available at\n");

    printf("https://github.com/magmasource/alphaMELTS\n");

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

#undef REALLOC
