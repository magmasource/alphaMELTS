/*
**        Support routines for adiabatic path finders
**
**        Version 5 for compatibility with Melts 3.0.x  May 1996
**        Version 6 for multiple liquids  July 1998
*/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "interface.h"
#include "silmin.h"                /* SILMIN structures include file        */

#include "status.h"

#ifdef DEBUG
#undef DEBUG
#endif

double extractMelt()       /* returns new reference entropy */
{

    int fractionateLiq = silminState->fractionateLiq, fractionateSol = silminState->fractionateSol,
        fractionateFlu = silminState->fractionateFlu;
    double liquidout, fluidout;

    if (silminState->fractionateLiq && silminState->fracLComp == NULL) {
        silminState->fracLComp = (double *) calloc((unsigned) nlc, sizeof(double));
    }
    if (silminState->fractionateFlu && silminState->fracSComp == NULL) {
        silminState->fracSComp    = (double **) calloc((unsigned) npc, sizeof(double *));
        silminState->nFracCoexist = (int *) calloc((unsigned) npc, sizeof(int));
    }

    if (silminState->fractionateLiq && (silminState->liquidMass != 0.0)) {

        //(void) doTraceElements(THIRD, silminState);

        if (silminState->minLiqPhi > 0.0) {
            double meltFraction = (silminState->liquidTD.v)/(silminState->bulkTD.v);
            /*liquidout = 1.0 - (MINphi / state->meltFraction);*/
            liquidout = 1.0 - (silminState->minLiqPhi / meltFraction) * (1.0 - meltFraction) / (1.0 - silminState->minLiqPhi);
        }
        else if (silminState->minF > 0.0) {
            double F = silminState->liquidMass / (silminState->liquidMass + silminState->solidMass);
            double minF = silminState->minF;

            /*liquidout = 1.0 - (MINf / F);*/
            liquidout = 1.0 - (minF/ F) * (1.0 - F) / (1.0 - minF);
        }
        else liquidout = 1.0; /* individual fractionation coeffs or default to old definition of MASSIN */
    }
    else liquidout = 0.0;
    if (liquidout < 0.0) liquidout = 0.0;

    if (silminState->fractionateFlu) {

        int i, j, haveWater = FALSE; //((calculationMode == MODE__MELTS) || (calculationMode == MODE_pMELTS));

        if (haveWater) {
            for (j=0; j<npc; j++) {
                if (solids[j].type == PHASE) {
                    // Check min phase strlen and false positives
                    if (!strncmp("water", solids[j].label, 5)) {
                        break;
                    }
                }
            }
        }
        else {
            for (j=0; j<npc; j++) {
                if (solids[j].type == PHASE) {
                    if (!strncmp("fluid", solids[j].label, 5)) {
                        break;
                    }
                }
            }
        }
        /* check j==npc */

        if (silminState->minFluPhi != 0.0) {
            double V, fluidFraction;

            if (haveWater) {
                V = (silminState->solidComp)[j][0]*(solids[j].cur).v;
            }
            else {
                /* Should really loop in here */
                double *m = (double *) malloc((size_t)      2*sizeof(double));
                double *r = (double *) malloc((size_t)      1*sizeof(double));

                for (i=0; i<solids[j].na; i++) m[i] = (silminState->solidComp)[j+1+i][0];

                (*solids[j].convert)(SECOND, THIRD, silminState->T, silminState->P, NULL, m, r, NULL, NULL, NULL, NULL, NULL);
                (*solids[j].vmix) (FIRST,
                    silminState->T, silminState->P, r, &V, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

                for (i=0; i<solids[j].na; i++) {
                    V       += m[i]*(solids[j+1+i].cur).v;
                }
                free(m);
                free(r);
            }
            fluidFraction = V/silminState->bulkTD.v;
            fluidout = 1.0 - (silminState->minFluPhi / fluidFraction) * (1.0 - fluidFraction) / (1.0 - silminState->minFluPhi);
        }
        else fluidout = 0.0; /* will be fractionated in fracSolids instead */
    }
    else fluidout = 0.0;
    if (fluidout < 0.0) fluidout = 0.0;

    if (silminState->fractionateLiq && liquidout != 0.0) {
        silminState->fractionateSol = FALSE; silminState->fractionateFlu = FALSE;
        //  doBatchFractionation((1.0-liquidout)*silminState->liquidMass);
        doBatchFractionation(liquidout);
        silminState->fractionateSol = fractionateSol; silminState->fractionateFlu = fractionateFlu;
    }

    /* here is the place to save information for the integrator */

    if (silminState->fractionateFlu && fluidout != 0.0) {
        silminState->fractionateLiq = FALSE; silminState->fractionateSol = FALSE;
        doBatchFractionation(fluidout);
        silminState->fractionateLiq = fractionateLiq; silminState->fractionateSol = fractionateSol;
    }

    return silminState->bulkTD.s;

}

double fractionateSolids() {

    int fractionateLiq = silminState->fractionateLiq, fractionateSol = silminState->fractionateSol,
        fractionateFlu = silminState->fractionateFlu;
    double solidout;

    if ((silminState->fractionateSol || silminState->fractionateFlu) && silminState->fracSComp == NULL) {
        silminState->fracSComp    = (double **) calloc((unsigned) npc, sizeof(double *));
        silminState->nFracCoexist = (int *) calloc((unsigned) npc, sizeof(int));
    }
#ifdef PHMELTS_ADJUSTMENTS
    if ((silminState->fractionateSol || silminState->fractionateFlu) && silminState->fracSolids == NULL) {
        //int i, haveWater = ((calculationMode == MODE__MELTS) || (calculationMode == MODE_pMELTS));

        /* so far only have 2 component fluid; assume at most one fluid for fractionation */
        silminState->fracLiquids = (double *) calloc((unsigned) nlc, sizeof(double));
        silminState->fracFluids = (double *) calloc((unsigned) 2, sizeof(double));

        // EVENTUALLY PACK THIS THE SAME AS INCSOLIDS? */
        silminState->fracSolids = (double *) calloc((unsigned) npc, sizeof(double));
        //if (silminState->fractionateSol) for (i=0; i<npc; i++) silminState->fracSolids[i] = 1.0;
    }
#endif

    //if (silminState->solidMass != 0.0) (void) doTraceElements(FOURTH, silminState);

    /* NEED TO CHANGE MASSIN IN SILMIN TOO */

    /* will be fractionated in extractMelt() instead */
    if (silminState->fractionateFlu && silminState->minFluPhi != 0.0) silminState->fractionateFlu = FALSE;

    if (silminState->fractionateSol || silminState->fractionateFlu) {
        if (silminState->maxF < 1.0) {
            double F = silminState->solidMass / (silminState->liquidMass + silminState->solidMass);
            double minF = 1.0 - silminState->maxF;
            solidout = 1.0 - (minF/ F) * (1.0 - F) / (1.0 - minF);
        }
        else solidout = 1.0; /* individual fractionation coeffs or default to old definition of MASSIN */
    }
    else solidout = 0.0;

    silminState->fractionateLiq = FALSE;
    if (solidout != 0.0) doBatchFractionation(solidout);
    silminState->fractionateLiq = fractionateLiq; silminState->fractionateFlu = fractionateFlu;

    return silminState->bulkTD.s;

}

/* Phase diagram and liquidus stuff in here */


/* pHMELTS part will go in here */
int adiabatFunc() {
    int success;

    /* ALPHAMELTS_IMPOSE_FO2 only needed for aH2O + Kress & Carmichael */
    while(!silmin());

    silminState->dspTstart = silminState->T - 273.15;
    silminState->dspPstart = silminState->P;

    success = (meltsStatus.status == 5);
    return success;

}


/* This is getLiquidProperties in create_managed.c, modified to be passed a state */
/* Note: returned volume is in cc */
double getDspLiquidProperties(SilminState *state, int nl, double *value)
{
    double *m, mass, moles, *r;
    int i, j, na, nr;

    na = nlc;
    nr = nlc - 1;

    m = (double *) malloc((size_t) na*sizeof(double));
    r = (double *) malloc((size_t) nr*sizeof(double));
    //mass =  state->liquidMass; could be more than one liquid
    for (i=0, mass=0.0, moles=0.0; i<na; i++) {
        m[i]  = ( state->liquidComp)[nl][i]; moles += m[i];
        for (j=0; j<nc; j++) mass += m[i]*(liquid[i].liqToOx)[j]*bulkSystem[j].mw;
    }
    (*conLiq)(SECOND, THIRD,  state->T,  state->P, NULL, m, r, NULL, NULL, NULL, NULL);

    /* The Gibbs energy of the phase (J)    */
    (*gmixLiq)(FIRST,  state->T,  state->P, r, &value[0], NULL, NULL);
    for (i=0, value[0] *= moles; i<na; i++) value[0] += m[i]*liquid[i].cur.g;

    /* The Enthalpy of the phase (J)        */
    (*hmixLiq)(FIRST,  state->T,  state->P, r, &value[1], NULL);
    for (i=0, value[1] *= moles; i<na; i++) value[1] += m[i]*liquid[i].cur.h;

    /* The Entropy of the phase (J/K)       */
    (*smixLiq)(FIRST,  state->T,  state->P, r, &value[2], NULL, NULL, NULL);
    for (i=0, value[2] *= moles; i<na; i++) value[2] += m[i]*liquid[i].cur.s;

    /* The Volume of the phase (cc)         */
    (*vmixLiq)(FIRST,  state->T,  state->P, r, &value[3], NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    for (i=0, value[3] *= moles; i<na; i++) value[3] += m[i]*liquid[i].cur.v;
    value[3] *= 10.0; /* joules/bar -> cc */

    /* The Heat Capacity of the phase (J/K) */
    (*cpmixLiq)(FIRST,  state->T,  state->P, r, &value[4], NULL, NULL);
    for (i=0, value[4] *= moles; i<na; i++) value[4] += m[i]*liquid[i].cur.cp;

    /* The Density of the phase (gm/cc)     */
    value[5] = (value[3] == 0.0) ? 0.0 : mass/value[3];

    /* The Viscosity of the phase (log poise)   */
    (*visLiq)(FIRST,  state->T,  state->P, r, &value[6]);

    free(m); if (na > 1) free(r);
    return mass;
}

/* This is getSolidProperties in the original create_managed.c, modified to be passed a state */
/* Note: returned volume is in cc */
double getDspSolidProperties(SilminState *state, int index, int ns, double *value)
{
    double *m, mass, moles, *r = NULL;
    int i, na, nr;

    na = solids[index].na;
    nr = solids[index].nr;

    m = (double *) malloc((size_t) na*sizeof(double));
    if (na == 1) {
        m[0]  = ( state->solidComp)[index][ns];
        if (( state->fractionateSol ||  state->fractionateFlu) && (( state->fracSComp)[index] != NULL))
            m[0] += ( state->fracSComp)[index][ns];
        mass  = m[0]*solids[index].mw;
        moles = m[0];
    } else {
        for (i=0, mass=0.0, moles=0.0; i<na; i++) {
            m[i]   = ( state->solidComp)[index+1+i][ns];
            if (( state->fractionateSol ||  state->fractionateFlu) && (( state->fracSComp)[index+1+i] != NULL))
                m[i] += ( state->fracSComp)[index+1+i][ns];
            mass  += m[i]*solids[index+1+i].mw;
            moles += m[i];
        }
        r = (double *) malloc((size_t) nr*sizeof(double));
        (*solids[index].convert)(SECOND, THIRD,  state->T,  state->P, NULL, m, r, NULL, NULL, NULL, NULL, NULL);
    }

    /* The Gibbs energy of the phase (J)    */
    if (na == 1) value[0] = m[0]*solids[index].cur.g;
    else {
        (*solids[index].gmix)(FIRST,  state->T,  state->P, r, &value[0], NULL, NULL, NULL);
        value[0] *= moles;
        for (i=0; i<na; i++) value[0] += m[i]*solids[index+1+i].cur.g;
    }

    /* The Enthalpy of the phase (J)        */
    if (na == 1) value[1] = m[0]*solids[index].cur.h;
    else {
        (*solids[index].hmix)(FIRST,  state->T,  state->P, r, &value[1]);
        value[1] *= moles;
        for (i=0; i<na; i++) value[1] += m[i]*solids[index+1+i].cur.h;
    }

    /* The Entropy of the phase (J/K)       */
    if (na == 1) value[2] = m[0]*solids[index].cur.s;
    else {
        (*solids[index].smix)(FIRST,  state->T,  state->P, r, &value[2], NULL, NULL);
        value[2] *= moles;
        for (i=0; i<na; i++) value[2] += m[i]*solids[index+1+i].cur.s;
    }

    /* The Volume of the phase (cc)         */
    if (na == 1) value[3] = m[0]*solids[index].cur.v;
    else {
        (*solids[index].vmix)(FIRST,  state->T,  state->P, r, &value[3], NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL);
        value[3] *= moles;
        for (i=0; i<na; i++) value[3] += m[i]*solids[index+1+i].cur.v;
    }
    value[3] *= 10.0; /* joules/bar -> cc */

    /* The Heat Capacity of the phase (J/K) */
    if (na == 1) value[4] = m[0]*solids[index].cur.cp;
    else {
        (*solids[index].cpmix)(FIRST,  state->T,  state->P, r, &value[4], NULL, NULL);
        value[4] *= moles;
        for (i=0; i<na; i++) value[4] += m[i]*solids[index+1+i].cur.cp;
    }

    /* The Density of the phase (gm/cc)     */
    value[5] = (value[3] == 0.0) ? 0.0 : mass/value[3];

    /* The Viscosity of the phase (poise)   */
    value[6] = 0.0;

    free(m); if (na > 1) free(r);
    return mass;
}

/* this is based on updateBulkADB from silmin_support.c */
void updateDspBulk(void)
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
