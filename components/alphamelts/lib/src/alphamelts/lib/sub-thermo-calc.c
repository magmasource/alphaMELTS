/*
**++
**  FACILITY:  subthermo_calc.c
**
**  MODULE DESCRIPTION:
**
**      routines to return thermodynamic properties (descended from
**        create_managed.c in Melts package; now updated to use silmin.c code)
**
**  MODIFICATION HISTORY:
**
**        V4.0-1        Paul D. Asimow  July 24, 1994
**                deleted all X11 dependencies, create_managed(), etc.
**                keep only useful init routines for non-X11 subroutine version
**        V4.0-2  Paul D. Asimow  July 27, 1994
**                add getSystemProperties(), getBulkSolidProperties(),
**                and addThermoData().  Modify all to pass ThermoData structures.
**        V5.0-1  Paul D. Asimow  November 11, 1994
**                Modify calls to *vmix and *cpmix for isentropic derivatives.
**                Gibbs now returns S(P,T)
**        V6.0-1  Paul D. Asimow  June 3, 1996
**                Compatible with Melts 3
**        V7.0-1  Paul D. Asimow  July 29, 1998
**        Multiple liquids.
**--
*/

#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "interface.h"
#include "silmin.h"
#include "adiabat.h"

void multiplyThermoData(ThermoData *target, double factor)
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

void addThermoData(ThermoData *value, ThermoData temp)
{
    value->g += temp.g;
    value->h += temp.h;
    value->s += temp.s;
    value->v += temp.v;
    value->cp += temp.cp;
    value->dcpdt += temp.dcpdt;
    value->dvdt += temp.dvdt;
    value->dvdp += temp.dvdp;
    value->d2vdt2 += temp.d2vdt2;
    value->d2vdp2 += temp.d2vdp2;
    value->d2vdtdp += temp.d2vdtdp;
}

/* From pHMELTS 3 but volume in J/bar */
double getLiquidProperties(SilminState *state, int ns, ThermoData *value, double *viscosity, int derivatives)
{
    double *m, mass, moles, *r;
    int i, j, na, nr;

    na = nlc;
    nr = nlc - 1;
    // mass = silminState->liquidMass;
    m = (double *) malloc((size_t) na*sizeof(double));
    r = (double *) malloc((size_t) nr*sizeof(double));

    for (i=0, mass=0.0, moles=0.0; i<na; i++) {
        m[i]  = (state->liquidComp)[ns][i]; moles += m[i];
        for (j=0; j<nc; j++) mass += m[i]*(liquid[i].liqToOx)[j]*bulkSystem[j].mw;
    }

    if(!testLiq(SIXTH, state->T, state->P, 0, 0, NULL, NULL, NULL, m)) {
        free(m);
        return -1.0;
    }

    for (i=0; i<nlc; i++) {
        gibbs(state->T, state->P, (char *) liquid[i].label, &(liquid[i].ref), &(liquid[i].liq), &(liquid[i].fus), &(liquid[i].cur));
    }

    (*conLiq)(SECOND, THIRD, state->T, state->P, NULL, m, r, NULL, NULL, NULL, NULL);

    /* The Gibbs energy of the phase (J)    */
    (*gmixLiq)(FIRST, state->T, state->P, r, &(value->g), NULL, NULL);
    for (i=0, value->g *= moles; i<na; i++) value->g += m[i]*(liquid[i].cur).g;

    /* The Enthalpy of the phase (J)        */
    (*hmixLiq)(FIRST, state->T, state->P, r, &(value->h), NULL);
    for (i=0, value->h *= moles; i<na; i++) value->h += m[i]*(liquid[i].cur).h;

    /* The Entropy of the phase (J/K) at P,T    */
    (*smixLiq)(FIRST, state->T, state->P, r, &(value->s), NULL, NULL, NULL);
    for (i=0, value->s *= moles; i<na; i++) value->s += m[i]*(liquid[i].cur).s;

    /* The Volume of the phase (cc)         */
    (*vmixLiq)(FIRST, state->T, state->P, r, &(value->v), NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    for (i=0, value->v *= moles; i<na; i++) value->v += m[i]*(liquid[i].cur).v;
    /*value->v *= 10.0;*/ /* joules/bar -> cc */

    /* The Heat Capacity of the phase (J/K) */
    (*cpmixLiq)(FIRST, state->T, state->P, r, &(value->cp),
        (double *) NULL, (double *) NULL);
    for (i=0, value->cp *= moles; i<na; i++) value->cp += m[i]*(liquid[i].cur).cp;

    /* Various Derivatives */
    if (derivatives) {
        value->dcpdt = 0.0;
        for (i=0;i<na;i++) value->dcpdt += m[i]*(liquid[i].cur).dcpdt; /* dcpdt */
        value->dvdt = 0.0;
        for (i=0;i<na;i++) value->dvdt += m[i]*(liquid[i].cur).dvdt;   /* dvdt  */
        value->dvdp = 0.0;
        for (i=0;i<na;i++) value->dvdp += m[i]*(liquid[i].cur).dvdp;   /* dvdp  */
        value->d2vdt2 = 0.0;
        for (i=0;i<na;i++) value->d2vdt2 += m[i]*(liquid[i].cur).d2vdt2;   /* d2vdt2  */
        value->d2vdp2 = 0.0;
        for (i=0;i<na;i++) value->d2vdp2 += m[i]*(liquid[i].cur).d2vdp2;   /* d2vdp2  */
        value->d2vdtdp = 0.0;
        for (i=0;i<na;i++) value->d2vdtdp += m[i]*(liquid[i].cur).d2vdtdp; /* d2vdtdp */
    }

    /* The Viscosity of the phase (log poise)   */
    (*visLiq)(FIRST, state->T, state->P, r, viscosity);

    free(m); if (na > 1) free(r);
    return mass;
}

double getSolidProperties(SilminState *state, int index, int ns, ThermoData *value, int derivatives)
{
    double *m, mass, moles, *r;
    int i, na, nr;

    na = solids[index].na;
    nr = solids[index].nr;

    m = (double *) malloc((size_t) na*sizeof(double));
    if (na == 1) {
        m[0]  = (state->solidComp)[index][ns];
        mass  = m[0]*solids[index].mw;
        moles = m[0];
    } else {
        for (i=0, mass=0.0, moles=0.0; i<na; i++) {
            m[i]   = (state->solidComp)[index+1+i][ns];
            mass  += m[i]*solids[index+1+i].mw;
            moles += m[i];
        }
        if(!(*solids[index].test)(SIXTH, state->T, state->P,
			      0, 0, NULL, NULL, NULL, m)) {
            free(m);
            return -1.0;
        }
        r = (double *) malloc((size_t) nr*sizeof(double));
        (*solids[index].convert)(SECOND, THIRD, state->T, state->P, NULL, m, r, NULL, NULL, NULL, NULL, NULL);
    }

    if (na == 1) {
        gibbs(state->T, state->P, (char *) solids[index].label, &(solids[index].ref), NULL, NULL, &(solids[index].cur));
    } else {
        for (i=0; i<na; i++) {
            gibbs(state->T, state->P, (char *) solids[index+1+i].label, &(solids[index+1+i].ref), NULL, NULL, &(solids[index+1+i].cur));
        }
    }

    /* The Gibbs energy of the phase (J)    */
    if (na == 1) value->g = m[0]*(solids[index].cur).g;
    else {
        (*solids[index].gmix)(FIRST, state->T, state->P, r, &(value->g), NULL, NULL, NULL);
        value->g *= moles;
        for (i=0; i<na; i++) value->g += m[i]*(solids[index+1+i].cur).g;
    }

    /* The Enthalpy of the phase (J)        */
    if (na == 1) value->h = m[0]*(solids[index].cur).h;
    else {
        (*solids[index].hmix)(FIRST, state->T, state->P, r, &(value->h));
        value->h *= moles;
        for (i=0; i<na; i++) value->h += m[i]*(solids[index+1+i].cur).h;
    }

    /* The Entropy of the phase (J/K) at P,T     */
    if (na == 1) value->s = m[0]*(solids[index].cur).s;
    else {
        (*solids[index].smix)(FIRST, state->T, state->P, r, &(value->s), NULL, NULL);
        value->s *= moles;
        for (i=0; i<na; i++) value->s += m[i]*(solids[index+1+i].cur).s;
    }

    /* The Volume of the phase (cc)         */
    if (na == 1) value->v = m[0]*(solids[index].cur).v;
    else {
        (*solids[index].vmix)(FIRST, state->T, state->P, r, &(value->v), NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL);
        value->v *= moles;
        for (i=0; i<na; i++) value->v += m[i]*(solids[index+1+i].cur).v;
    }
    /*value->v *= 10.0; *//* joules/bar -> cc */

    /* The Heat Capacity of the phase (J/K) */
    if (na == 1) value->cp = m[0]*(solids[index].cur).cp;
    else {
        (*solids[index].cpmix)(FIRST, state->T, state->P, r, &(value->cp), NULL, NULL);
        value->cp *= moles;
        for (i=0; i<na; i++) value->cp += m[i]*(solids[index+1+i].cur).cp;
    }

    /* Various Derivatives */
    if (derivatives) {
        if (na == 1) {
                value->dcpdt = m[0]*(solids[index].cur).dcpdt;
                value->dvdt = m[0]*(solids[index].cur).dvdt;
                value->dvdp = m[0]*(solids[index].cur).dvdp;
                value->d2vdt2 = m[0]*(solids[index].cur).d2vdt2;
                value->d2vdp2 = m[0]*(solids[index].cur).d2vdp2;
                value->d2vdtdp = m[0]*(solids[index].cur).d2vdtdp;
        } else {
            (*solids[index].cpmix)(SECOND, state->T, state->P, r, NULL, &(value->dcpdt), NULL);
            value->dcpdt *= moles;
            for (i=0; i<na; i++) value->dcpdt += m[i]*(solids[index+1+i].cur).dcpdt;
            (*solids[index].vmix)(FOURTH | FIFTH | SIXTH | SEVENTH | EIGHTH, state->T, state->P, r, NULL, NULL, NULL,
                &(value->dvdt), &(value->dvdp), &(value->d2vdt2), &(value->d2vdtdp), &(value->d2vdp2), NULL, NULL);
            value->dvdt *= moles;
            for (i=0; i<na; i++) value->dvdt += m[i]*(solids[index+1+i].cur).dvdt;
            value->dvdp *= moles;
            for (i=0; i<na; i++) value->dvdp += m[i]*(solids[index+1+i].cur).dvdp;
            value->d2vdt2 *= moles;
            for (i=0; i<na; i++) value->d2vdt2 += m[i]*(solids[index+1+i].cur).d2vdt2;
            value->d2vdp2 *= moles;
            for (i=0; i<na; i++) value->d2vdp2 += m[i]*(solids[index+1+i].cur).d2vdp2;
            value->d2vdtdp *= moles;
            for (i=0; i<na; i++) value->d2vdtdp += m[i]*(solids[index+1+i].cur).d2vdtdp;
        }
    }

    free(m); if (na > 1) free(r);
    return mass;
}

double getSystemProperties(SilminState *state, int derivatives)
{
    double mass, totalMass, volSol; // Lviscosity;
    int i, j;
    /*double *moles;*/

    /* these pointers used to be passed instead (plus viscosity and phi) */
    ThermoData *value = &(state->bulkTD);

    totalMass = 0.0; multiplyThermoData(value, 0.0);

    if (state->liquidMass != 0.0) {
        mass = getBulkLiquidProperties(state, derivatives);
        if(mass < 0.0) return -1.0;
        else state->liquidMass = mass;
        addThermoData(value, state->liquidTD);
    }

    mass = getBulkSolidProperties(state, derivatives);
    if(mass < 0.0) return -1.0;
    else state->solidMass = mass;
    totalMass = state->liquidMass + state->solidMass;
    if (mass > 0.0) addThermoData(value, state->solidTD);

    volSol = state->solidTD.v;
    //Lviscosity = state->liqViscosity;

    /* The Viscosity of the system (log poise)   A == 0.5; not meaningful for two liquids */
    //state->viscosity = (volSol < 0.5*value->v) ?
    //               Lviscosity - 2.0*log10(1.0 - 2.0*volSol/value->v) : 0.0;
    /* Melt Fraction */
    //state->meltFraction = (value->v == 0.0) ? 0.0 : (1.0 - volSol/value->v);

    /* fO2 */
    if (state->fo2Path == FO2_NONE) {
        double muO2;
        double *moles = (double *) malloc(nc*sizeof(double));
        if (state->liquidMass != 0.0) {// && getenv("ALPHAELTS_ALTERNATIVE_FO2") == NULL) {
            for (i=0;i<nc;i++) {
                for (j=0, moles[i]=0.0; j<nlc; j++) moles[i] += (state->liquidComp)[0][j]*(liquid[j].liqToOx)[i];
            }
            conLiq(FIRST, SEVENTH, state->T, state->P, moles, NULL, NULL, NULL, NULL, NULL, &(state->fo2));
        }
        else { //if(getenv("ALPHAMELTS_ALTERNATIVE_FO2") != NULL || getenv("ALPHAMELTS_LIQUID_FO2") == NULL) {
            subsolidusmuO2(FIRST, &muO2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
            state->fo2 = muO2/(R*state->T*log(10.0));
        }
        free(moles);
    }
    /*free(moles);*/
    return totalMass;
}

double getBulkSolidProperties(SilminState *state, int derivatives)
{
    double mass, totalMass;
    int i, j;
    ThermoData temp;

    /* these pointers used to be passed instead */
    ThermoData *value = &(state->solidTD);

    totalMass = 0.0; multiplyThermoData(value, 0.0);

    for (i=0; i<npc; i++) {
        if (solids[i].type == PHASE) {
            for (j=0; j<(state->nSolidCoexist)[i]; j++) {
                mass = getSolidProperties(state, i, j, &temp, derivatives);
                if(mass < 0.0) return -1.0;
                else totalMass += mass;
                addThermoData(value, temp);
            }
        }
    }

    state->solidMass = totalMass;
    return totalMass;
}

double getBulkLiquidProperties(SilminState *state, int derivatives)
{
    double mass, totalMass, viscosity;
    int j, hasLiquid = (state->liquidMass != 0.0);
    ThermoData temp;

    /* these pointers used to be passed instead */
    ThermoData *value = &(state->liquidTD);
    //double *Lviscosity = &(state->liqViscosity);

    totalMass = 0.0; multiplyThermoData(value, 0.0);

    if (hasLiquid) {
        for (j=0; j<state->nLiquidCoexist; j++) {
            mass = getLiquidProperties(state,  j, &temp, &viscosity, derivatives);
            if(mass < 0.0) return -1.0;
            else totalMass += mass;
            //if (j==0) *Lviscosity = viscosity;
            addThermoData(value, temp);
        }
    }

    state->liquidMass = totalMass;
    return totalMass;
}