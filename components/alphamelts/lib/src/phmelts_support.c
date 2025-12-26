/*
**        Support routines for adiabatic path finders
**
**        Version 5 for compatibility with Melts 3.0.x  May 1996
**        Version 6 for multiple liquids  July 1998
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "silmin.h"
#include "adiabat.h"
#include "phmelts.h"

#include "liq_struct_data.h"
#include "sol_struct_data.h"

#define REALLOC(x, y) (((x) == NULL) ? malloc(y) : realloc((x), (y)))
#define SQUARE(x) ((x)*(x))

#define REC   134

void initializeLibrary(void) {

    if (silminInputData.name == NULL)  silminInputData.name  = (char *) calloc((unsigned) (REC+1), sizeof(char));
    if (silminInputData.title == NULL) silminInputData.title = (char *) calloc((unsigned) (REC+1), sizeof(char));

    if ((calculationMode == MODE__MELTS) || (calculationMode > MODE__MELTSandCO2_H2O)) {
        liquid = meltsLiquid;
        solids = meltsSolids;
        nlc    = meltsNlc;
        nls    = meltsNls;
        npc    = meltsNpc;
    } else if ((calculationMode == MODE__MELTSandCO2) || (calculationMode == MODE__MELTSandCO2_H2O)) {
        liquid = meltsFluidLiquid;
        solids = meltsFluidSolids;
        nlc = meltsFluidNlc;
        nls = meltsFluidNls;
        npc = meltsFluidNpc;
    } else if (calculationMode == MODE_pMELTS) {
        liquid = pMeltsLiquid;
        solids = pMeltsSolids;
        nlc    = pMeltsNlc;
        nls    = pMeltsNls;
        npc    = pMeltsNpc;
    } else {
        // It should not be possible to get here through normal interfaces;
        printf("ALPHAMELTS_CALC_MODE not set in initializeLibrary!\n\n");
        printf("Please report the sequence of events which led to this disaster.\n\n");
        calculationMode = MODE_xMELTS;
    }
    InitComputeDataStruct();
}

SilminState *reAllocSilminStatePointer(SilminState *q, size_t zOld, size_t z)
{
    int i, j;
    SilminState *p = NULL;

    if(z > zOld) {

        p = (SilminState *) realloc(q, z*sizeof(SilminState));
        for(i= (int) zOld; i< (int) z; i++) {
            /* assume initialization with zero bytes yields default entries */

            p[i].nLiquidCoexist = 1;
            p[i].liquidMass     = 0.0;

            /* allocate space for mandatory entries with zero bytes         */
            p[i].bulkComp      = (double *)     calloc((size_t)  nc, sizeof(double));
            p[i].dspBulkComp   = (double *)     calloc((size_t)  nc, sizeof(double *));
            p[i].liquidComp    = (double **)    calloc((size_t)  1, sizeof(double *));
            p[i].liquidComp[0]  = (double *)     calloc((size_t) nlc,     sizeof(double));
            p[i].liquidDelta   = (double **)    calloc((size_t)  1, sizeof(double *));
            p[i].liquidDelta[0] = (double *)     calloc((size_t) nlc,     sizeof(double));
            p[i].solidComp     = (double **)    calloc((size_t) npc, sizeof(double *));
            p[i].nSolidCoexist = (int *)        calloc((size_t) npc, sizeof(int));
            p[i].solidDelta    = (double **)    calloc((size_t) npc, sizeof(double *));
            p[i].incSolids     = (int *)        calloc((size_t) npc, sizeof(int));
            p[i].cylSolids     = (int *)        calloc((size_t) npc, sizeof(int));

            for (j=0; j<npc; j++) {
                (p[i].solidComp)[j]  = (double *) calloc((size_t)   1, sizeof(double));
      	(p[i].solidDelta)[j] = (double *) calloc((size_t)   1, sizeof(double));
            }
            p[i].fractionateSol = 0;
            p[i].fractionateLiq = 0;
            p[i].fractionateFlu = 0;
            p[i].assimilate     = 0;

            p[i].nFracCoexist   = NULL;
            p[i].fracSComp      = NULL;
            p[i].fracLComp      = NULL;

            p[i].dspAssimComp = NULL;
            p[i].assimComp    = NULL;
            p[i].nDspAssimComp = NULL;
            p[i].nAssimComp    = NULL;

            p[i].ySol = NULL;
            p[i].yLiq = NULL;

        }
    }
    else if(z < zOld) {

        for(i = (int) z; i < (int) zOld; i++) {

            for (j=0; j<MAX(1,q[i].nLiquidCoexist);j++) {
                free((q[i].liquidComp)[j]);
                free((q[i].liquidDelta)[j]);
            }

            for (j=0;j<npc;j++) {
                free((q[i].solidDelta)[j]);
                free((q[i].solidComp)[j]);
            }

            free(q[i].incSolids);
            free(q[i].cylSolids);
            free(q[i].solidDelta);
            free(q[i].nSolidCoexist);
            free(q[i].solidComp);
            free(q[i].liquidDelta);
            free(q[i].liquidComp);
            free(q[i].dspBulkComp);
            free(q[i].bulkComp);

            /* NEED TO FREE FRAC AND ASSIM COMPS IN HERE */
#ifdef PHMELTS_ADJUSTMENTS
    if (q[i].fracSComp != NULL) {
        for (j=0; j<npc; j++) if ((q[i].fracSComp)[j] != NULL) free((q[i].fracSComp)[j]);
        free(q[i].fracSComp);
    }
    if (q[i].nFracCoexist != NULL) free(q[i].nFracCoexist);

    if (q[i].fracLComp != NULL) free(q[i].fracLComp);

    if (q[i].dspAssimComp != NULL) {
        for (j=0; j<(npc+nc); j++)  if ((q[i].dspAssimComp)[j] != NULL) free((q[i].dspAssimComp)[j]);
        free(q[i].dspAssimComp);
    }
    if (q[i].assimComp != NULL) {
        for (j=0; j<(npc+nlc); j++) if ((q[i].assimComp)[j]    != NULL) free((q[i].assimComp)[j]);
        free(q[i].assimComp);
    }
    if (q[i].nDspAssimComp != NULL) free(q[i].nDspAssimComp);
    if (q[i].nAssimComp    != NULL) free(q[i].nAssimComp);

    if ((q[i].ySol) != NULL) free(q[i].ySol);
    if ((q[i].yLiq) != NULL) free(q[i].yLiq);
#endif

        }
        p = (SilminState *) realloc(q, z*sizeof(SilminState));

    }
    else p = q;
    return p;
}

void copyStateInfo(SilminState *target, SilminState *source)
{
   target = copySilminStateStructure(source, target);
}

/*void freeSilminStatePointer(SilminState *p) {

    destroySilminStateStructure((void *) p);

    }*/

#undef REALLOC


