/* Module initial_guess.c contains pdaNorm, a method for obtaining an
   initial guess in the subsolidus
**        V1.0-1        Paul D. Asimow  March 22, 1995
**                extracted from solid_support.c
*/

/*#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "silmin.h"
#include "interface.h"

int pdaNorm() {
    double SiO2, TiO2, Al2O3, Fe2O3, Cr2O3, FeO, MgO, CaO, Na2O, MnO, CoO, NiO, K2O, P2O5, H2O, CO2;
    int garnet, spinel, plag, olivine, quartz;
    int spindex, gtindex, plindex, olindex, qindex, cpxindex, opxindex, whitindex, windex;
    int i, j, k, l, success = FALSE, *incSolids, hasLiquid = (silminState->liquidMass != 0.0);

    incSolids = calloc((size_t) (npc+1), sizeof(int));

    /* find the phases and store their indices.  This avoid the assumption that the list of phases hasn't changed. */
    for (i=0, j=0; i<npc; i++) {
        if (solids[i].type == PHASE) {

            if (silminState->nSolidCoexist[i]) for (k = 0; k<silminState->nSolidCoexist[i]; k++)
                for (l = 0; l<=solids[i].na; l++) silminState->solidComp[i+l][k] = 0.0;
            if (silminState->nSolidCoexist[i] > 1) silminState->nSolidCoexist[i] = 1;

            if      (!strcmp(solids[i].label, "olivine"))       olindex = i;
            else if (!strcmp(solids[i].label, "garnet"))        gtindex = i;
            else if (!strcmp(solids[i].label, "clinopyroxene")) cpxindex = i;
            else if (!strcmp(solids[i].label, "orthopyroxene")) opxindex = i;
            else if (!strcmp(solids[i].label, "spinel"))        spindex = i;
            else if (!strcmp(solids[i].label, "plagioclase"))   plindex = i;
            else if (!strcmp(solids[i].label, "quartz"))        qindex = i;
            else if (!strcmp(solids[i].label, "whitlockite"))   whitindex = i;
            else if (!strcmp(solids[i].label, "fluid"))         windex = i;
            /* need to zero all other phases in case a hydrous phase is already present */
            else if (silminState->nSolidCoexist[i]) silminState->nSolidCoexist[i] = 0;

            incSolids[i] = silminState->incSolids[j++];
        }
    }

    spinel  = silminState->nSolidCoexist[spindex];
    plag    = silminState->nSolidCoexist[plindex];
    garnet  = silminState->nSolidCoexist[gtindex];
    olivine = silminState->nSolidCoexist[olindex];
    quartz  = silminState->nSolidCoexist[qindex];

    while (!success) {

        /* find the oxides and store their indices; this avoids the assumption that their order hasn't changed. */
        for (i=0; i<nc; i++) {
            if      (!strcmp(bulkSystem[i].label, "SiO2"))  SiO2  = silminState->bulkComp[i];
            else if (!strcmp(bulkSystem[i].label, "TiO2"))  TiO2  = silminState->bulkComp[i];
            else if (!strcmp(bulkSystem[i].label, "Al2O3")) Al2O3 = silminState->bulkComp[i];
            else if (!strcmp(bulkSystem[i].label, "Fe2O3")) Fe2O3 = silminState->bulkComp[i];
            else if (!strcmp(bulkSystem[i].label, "Cr2O3")) Cr2O3 = silminState->bulkComp[i];
            else if (!strcmp(bulkSystem[i].label, "FeO"))   FeO   = silminState->bulkComp[i];
            else if (!strcmp(bulkSystem[i].label, "MgO"))   MgO   = silminState->bulkComp[i];
            else if (!strcmp(bulkSystem[i].label, "NiO"))   NiO   = silminState->bulkComp[i];
            else if (!strcmp(bulkSystem[i].label, "MnO"))   MnO   = silminState->bulkComp[i];
            else if (!strcmp(bulkSystem[i].label, "CoO"))   CoO   = silminState->bulkComp[i];
            else if (!strcmp(bulkSystem[i].label, "CaO"))   CaO   = silminState->bulkComp[i];
            else if (!strcmp(bulkSystem[i].label, "Na2O"))  Na2O  = silminState->bulkComp[i];
            else if (!strcmp(bulkSystem[i].label, "K2O"))   K2O   = silminState->bulkComp[i];
            else if (!strcmp(bulkSystem[i].label, "P2O5"))  P2O5  = silminState->bulkComp[i];
            else if (!strcmp(bulkSystem[i].label, "H2O"))   H2O   = silminState->bulkComp[i];
            else if (!strcmp(bulkSystem[i].label, "CO2"))   CO2   = silminState->bulkComp[i];
        }

        /* olivine or quartz ? How many pyroxenes? */
        //if (!quartz && !olivine) printf("WARNING: Not sure how to solve w/o quartz or olivine!\n");

        /* Rule 0 choose aluminous phase(s) */
        if (silminState->P < 10000) {
                if (!incSolids[plindex]) {
                    printf("Need plagioclase but it's suppressed; trying spinel instead!\n");
                    if (!incSolids[spindex]) {
                        printf("Need spinel but it's suppressed!\n");
                        break;
                    }
                    if(!spinel) printf("Adding spinel to assemblage.\n");
                    spinel = silminState->nSolidCoexist[spindex] = 1;
                }
                else {
                    if(!plag) printf("Adding plagioclase to assemblage.\n");
                    plag = silminState->nSolidCoexist[plindex] = 1;
                }
        }
        else if (silminState->P > 20000) {
                if (!incSolids[gtindex]) {
                    printf("Need garnet but it's suppressed; trying spinel instead!\n");
                    if (!incSolids[spindex]) {
                        printf("Need spinel but it's suppressed!\n");
                        break;
                    }
                    if(!spinel) printf("Adding spinel to assemblage.\n");
                    spinel = silminState->nSolidCoexist[spindex] = 1;
                }
                else {
                    if(!garnet) printf("Adding garnet to assemblage.\n");
                    garnet = silminState->nSolidCoexist[gtindex] = 1;
                }
        }
        if (silminState->P >  5000 && silminState->P < 30000) {
                if (!incSolids[spindex]) {
                    printf("Need spinel but it's suppressed!\n");
                }
                if(!spinel) printf("Adding spinel to assemblage.\n");
                spinel = silminState->nSolidCoexist[spindex] = 1;
        }
        if (!plag && !spinel && !garnet) break;

        if (!silminState->nSolidCoexist[cpxindex] || !silminState->nSolidCoexist[opxindex]) {

            //printf("WARNING: Not sure how to solve without opx and cpx!\n");

            if(!incSolids[cpxindex]) {
                printf("Need clinopyroxene but it's suppressed!\n");
                break;
            }
            if(!silminState->nSolidCoexist[cpxindex]) printf("Adding clinopyroxene to assemblage.\n");
            silminState->nSolidCoexist[cpxindex] = 1;

            if(!incSolids[opxindex]) {
                printf("Need orthopyroxene but it's suppressed!\n");
                break;
            }
            if(!silminState->nSolidCoexist[opxindex]) printf("Adding orthopyroxene to assemblage.\n");
            silminState->nSolidCoexist[opxindex] = 1;

        }

        /* Rule 1 Cr2O3 and P2O5 and transition metals */
        if (Cr2O3 != 0.0) {
            if (!incSolids[spindex]) {
                printf("Need spinel but it's suppressed!\n");
                break;
            }
            if(!spinel) printf("Adding spinel to handle Cr2O3.\n");
            spinel = silminState->nSolidCoexist[spindex] = 1;

            for (i=1; i<=solids[spindex].na; i++) if (!strcmp(solids[spindex+i].label, "chromite")) break;
            silminState->solidComp[spindex+i][0] = Cr2O3;      /* chromite */
            Cr2O3 -= silminState->solidComp[spindex+i][0];
            FeO -= silminState->solidComp[spindex+i][0];
        }

        if (P2O5 != 0.0) {  /* whitlockite */
            if (!incSolids[whitindex]) {
                printf("Need whitlockite but it's suppressed!\n");
                break;
            }
            if(!silminState->nSolidCoexist[whitindex]) printf("Adding whitlockite to handle P2O5.\n");
            silminState->nSolidCoexist[whitindex] = 1;
            silminState->solidComp[whitindex][0] = P2O5;
            P2O5 -= silminState->solidComp[whitindex][0];
            CaO -= 3.0*silminState->solidComp[whitindex][0];
        }

        if (MnO != 0.0 || NiO != 0.0 || CoO != 0.0) {
            if (!incSolids[olindex]) {
                printf("Need olivine but it's suppressed!\n");
                break;
            }
            if(!olivine) printf("Adding olivine to handle transition metals.\n");
            olivine = silminState->nSolidCoexist[olindex] = 1;
        }

        /* Rule 2 Na2O, K2O */
        if (K2O != 0.0) {
            if (!incSolids[plindex]) {
                printf("Need feldspar but it's suppressed!\n");
                break;
            }
            if (!plag) printf("Adding feldspar to handle K2O.\n");
            plag = silminState->nSolidCoexist[plindex] = 1;
        }

        if (plag) {
            for (i=1; i<=solids[plindex].na; i++) if (!strcmp(solids[plindex+i].label, "albite")) break;
            silminState->solidComp[plindex+i][0] = 1.8*Na2O;   /* albite */
            Na2O  -= 0.5*silminState->solidComp[plindex+i][0];
            Al2O3 -= 0.5*silminState->solidComp[plindex+i][0];
            SiO2  -= 3.0*silminState->solidComp[plindex+i][0];

            for (i=1; i<=solids[plindex].na; i++) if (!strcmp(solids[plindex+i].label, "highsanidine")) break;
            silminState->solidComp[plindex+i][0] = K2O;/* sanidine */
            Na2O  -= 0.5*silminState->solidComp[plindex+i][0];
            Al2O3 -= 0.5*silminState->solidComp[plindex+i][0];
            SiO2  -= 3.0*silminState->solidComp[plindex+i][0];
        }

        for (i=1; i<=solids[cpxindex].na; i++) if (!strcmp(solids[cpxindex+i].label, "jadeite")) break;
        silminState->solidComp[cpxindex+i][0] = 1.32*Na2O;    /* cpx jadeite */
        Al2O3 -= 0.5*silminState->solidComp[cpxindex+i][0];
        Na2O  -= 0.5*silminState->solidComp[cpxindex+i][0];
        SiO2  -= 2.0*silminState->solidComp[cpxindex+i][0];
        silminState->solidComp[opxindex+i][0] = 2.0*Na2O;             /* opx jadeite */
        Na2O  -= 0.5*silminState->solidComp[opxindex+i][0];
        Al2O3 -= 0.5*silminState->solidComp[opxindex+i][0];
        SiO2  -= 2.0*silminState->solidComp[opxindex+i][0];

        /* Rule 3 TiO2 and Fe2O3 in spinel */
        if (spinel) {
            for (i=1;i<=solids[spindex].na;i++) if (!strcmp(solids[spindex+i].label, "ulvospinel")) break;
            silminState->solidComp[spindex+i][0] = 0.25*TiO2;  /* ulvospinel */
            TiO2 -= silminState->solidComp[spindex+i][0];
            FeO  -= 2.0*silminState->solidComp[spindex+i][0];

            for (i=1; i<=solids[spindex].na; i++) if (!strcmp(solids[spindex+i].label, "magnetite")) break;
            silminState->solidComp[spindex+i][0] = 0.4*Fe2O3; /* magnetite */
            Fe2O3 -= silminState->solidComp[spindex+i][0];
            FeO  -= silminState->solidComp[spindex+i][0];
        }

        /* Rule 4 Al2O3 */
        if (spinel) {
            int spinelcomp, hercynite, magnetite, chromite, ulvospinel;

            for (i=1;i<=solids[spindex].na;i++){
                if (!strcmp(solids[spindex+i].label, "spinel")) spinelcomp = i;
                else if (!strcmp(solids[spindex+i].label, "hercynite")) hercynite = i;
                else if (!strcmp(solids[spindex+i].label, "magnetite")) magnetite = i;
                else if (!strcmp(solids[spindex+i].label, "chromite")) chromite = i;
                else if (!strcmp(solids[spindex+i].label, "ulvospinel")) ulvospinel = i;
            }
            silminState->solidComp[spindex+spinelcomp][0] =  0.5*Al2O3/(spinel + plag + garnet); /* spinel */
            silminState->solidComp[spindex+hercynite][0]  = -0.2*Al2O3/(spinel + plag + garnet); /* hercynite */

            if (silminState->solidComp[spindex+hercynite][0] + silminState->solidComp[spindex+magnetite][0] +
      	  silminState->solidComp[spindex+chromite][0] + silminState->solidComp[spindex+ulvospinel][0]*2.0 < 0.0)
                silminState->solidComp[spindex+hercynite][0] = 0.001 - silminState->solidComp[spindex+magnetite][0] +
	        silminState->solidComp[spindex+chromite][0] + silminState->solidComp[spindex+ulvospinel][0]*2.0; /* enforce composition limit */
            MgO   -= silminState->solidComp[spindex+spinelcomp][0];
            FeO   -= silminState->solidComp[spindex+hercynite][0];
            Al2O3 -= (silminState->solidComp[spindex+hercynite][0] + silminState->solidComp[spindex+spinelcomp][0]);
        }

        if (plag) {
            for (i=1;i<=solids[plindex].na;i++) if (!strcmp(solids[plindex+i].label, "anorthite")) break;
            silminState->solidComp[plindex+i][0] = 0.64*Al2O3; /* anorthite */
            CaO -= silminState->solidComp[plindex+i][0];
            Al2O3 -= silminState->solidComp[plindex+i][0];
            SiO2 -= 2.0*silminState->solidComp[plindex+i][0];
        }

        if (garnet) {
            double denominator = 2.0*FeO + 3.0*MgO + 0.5*CaO;

            for (i=1;i<=solids[gtindex].na;i++) if (!strcmp(solids[gtindex+i].label, "almandine")) break;
            silminState->solidComp[gtindex+i][0] = 0.6*Al2O3*(2.0*FeO/denominator);  /* almandine */
            FeO -= 3.0*silminState->solidComp[gtindex+i][0];

            for (i=1;i<=solids[gtindex].na;i++) if (!strcmp(solids[gtindex+i].label, "grossular")) break;
            silminState->solidComp[gtindex+i][0] = 0.6*Al2O3*(0.5*CaO/denominator);  /* grossular */
            CaO -= 3.0*silminState->solidComp[gtindex+i][0];

            for (i=1;i<=solids[gtindex].na;i++) if (!strcmp(solids[gtindex+i].label, "pyrope")) break;
            silminState->solidComp[gtindex+i][0] = 0.6*Al2O3*(3.0*MgO/denominator);       /* pyrope */
            MgO   -= 3.0*silminState->solidComp[gtindex+i][0];
            SiO2  -= 3.0*0.6*Al2O3;
            Al2O3 -= 0.6*Al2O3;
        }

        /* Rule 5 TiO2, Fe2O3, Al2O3 */
        for (i=1;i<=solids[cpxindex].na;i++) if (!strcmp(solids[cpxindex+i].label, "alumino-buffonite")) break;
        silminState->solidComp[opxindex+i][0] = 0.75*(Al2O3 + TiO2 - Fe2O3);  /* opx aluminobuffonite */
        silminState->solidComp[cpxindex+i][0] = 0.25*(Al2O3 + TiO2 - Fe2O3);  /* cpx aluminobuffonite */

        for (i=1;i<=solids[cpxindex].na;i++) if (!strcmp(solids[cpxindex+i].label, "buffonite")) break;
        silminState->solidComp[opxindex+i][0] = 0.75*(-Al2O3 + TiO2 + Fe2O3); /* opx buffonite */
        silminState->solidComp[cpxindex+i][0] = 0.25*(-Al2O3 + TiO2 + Fe2O3); /* cpx buffonite */

        for (i=1;i<=solids[cpxindex].na;i++) if (!strcmp(solids[cpxindex+i].label, "essenite")) break;
        silminState->solidComp[opxindex+i][0] = 0.75*(Al2O3 - TiO2 + Fe2O3);  /* opx essenite */
        silminState->solidComp[cpxindex+i][0] = 0.25*(Al2O3 - TiO2 + Fe2O3);  /* cpx essenite */

        CaO   -= (Al2O3 + TiO2 + Fe2O3);
        MgO   -= TiO2;
        SiO2  -= (Al2O3 + TiO2 + Fe2O3);
        Al2O3 -= Al2O3;
        TiO2  -= TiO2;
        Fe2O3 -= Fe2O3;

        /* Rule 6 CaO */
        {
            int hd, di;

            for (i=1; i<=solids[cpxindex].na; i++) {
                if (!strcmp(solids[cpxindex+i].label, "hedenbergite")) hd = i;
                else if (!strcmp(solids[cpxindex+i].label, "diopside")) di = i;
            }
            silminState->solidComp[cpxindex+hd][0] = 0.15*CaO;                /* cpx hedenbergite */
            FeO -= 0.15*CaO; SiO2 -= 0.30*CaO;

            silminState->solidComp[opxindex+hd][0] = (olivine ? CaO : FeO);                        /* opx hedenbergite */
            silminState->solidComp[opxindex+di][0] = -1.05*silminState->solidComp[opxindex+hd][0]; /* opx diopside */
            CaO -= silminState->solidComp[opxindex+hd][0] + silminState->solidComp[opxindex+di][0] + silminState->solidComp[cpxindex+hd][0];
            MgO -= silminState->solidComp[opxindex+di][0];
            SiO2 -= 2.0*(silminState->solidComp[opxindex+hd][0] + silminState->solidComp[opxindex+di][0]);
            FeO -= silminState->solidComp[opxindex+hd][0];

            silminState->solidComp[cpxindex+di][0] = (olivine ? 0.833*CaO : CaO);    /* cpx diopside */
            SiO2 -= 2.0*silminState->solidComp[cpxindex+di][0];
            MgO -= silminState->solidComp[cpxindex+di][0];
            CaO -= silminState->solidComp[cpxindex+di][0];

            for (i=1;i<=solids[olindex].na;i++) if (!strcmp(solids[olindex+i].label, "monticellite")) break;
            silminState->solidComp[olindex+i][0] = ((CaO > 0.0 && olivine)?CaO:0.0); /* monticellite */
            CaO -= silminState->solidComp[olindex+i][0];
            MgO -= silminState->solidComp[olindex+i][0];
            SiO2 -= silminState->solidComp[olindex+i][0];
        }

        /* Rule 7 MnO, NiO, CoO */
        if (olivine) {
            for (i=1; i<=solids[olindex].na; i++) if (!strcmp(solids[olindex+i].label, "tephroite")) break;
            silminState->solidComp[olindex+i][0] = MnO/2.0;
            SiO2 -= MnO/2.0; MnO -= MnO;

            for (i=1; i<=solids[olindex].na; i++) if (!strcmp(solids[olindex+i].label, "co-olivine")) break;
            silminState->solidComp[olindex+i][0] = CoO/2.0;
            SiO2 -= CoO/2.0; CoO -= CoO;

            for (i=1; i<=solids[olindex].na; i++) if (!strcmp(solids[olindex+i].label, "ni-olivine")) break;
            silminState->solidComp[olindex+i][0] = NiO/2.0;
            SiO2 -= NiO/2.0; NiO -= NiO;
        }

        /* Rule 8 FeO, MgO, SiO2 */
        for (i=1; i<=solids[olindex].na; i++) if (!strcmp(solids[olindex+i].label, "fayalite")) break;
        silminState->solidComp[olindex+i][0] = (olivine ? FeO/2.0 : 0.0); /* fayalite */
        FeO  -= 2.0*silminState->solidComp[olindex+i][0];
        SiO2 -= silminState->solidComp[olindex+i][0];

        for (i=1; i<=solids[olindex].na; i++) if (!strcmp(solids[olindex+i].label, "forsterite")) break;
        silminState->solidComp[olindex+i][0] = (((MgO-SiO2)>0 && olivine) ? (MgO-SiO2) : 0.0);    /* forsterite */
        MgO  -= 2.0*silminState->solidComp[olindex+i][0];
        SiO2 -= silminState->solidComp[olindex+i][0];

        for (i=1; i<=solids[cpxindex].na; i++) if (!strcmp(solids[cpxindex+i].label, "clinoenstatite")) break;
        silminState->solidComp[opxindex+i][0] = 0.5*0.95*MgO; /* opx enstatite */
        silminState->solidComp[cpxindex+i][0] = 0.5*0.05*MgO; /* cpx enstatite */
        SiO2 -= MgO;
        MgO  -= MgO;

        if (SiO2 > 0.0) {
            if (!incSolids[qindex]) {
                printf("Need quartz but it's suppressed!\n");
                break;
            }
            if(!quartz) printf("Adding quartz to handle excess SiO2.\n");
            silminState->solidComp[qindex][0] = SiO2;
            silminState->nSolidCoexist[qindex] = 1;
            SiO2 -= SiO2;
        }

        /* Rule 9 H2O - adding CO2 in here too*/
        if (H2O != 0.0) {
            if (!incSolids[windex]) {
                printf("Need fluid but it's suppressed!\n");
                break;
            }
            if(!silminState->nSolidCoexist[windex]) printf("Adding fluid to initial guess assemblage\n");
            incSolids[windex] = 1;
            silminState->nSolidCoexist[windex] = 1;
            silminState->solidComp[windex][0] = H2O;
            H2O -= H2O;
        }

        if (CO2 != 0.0) {
            if (!incSolids[windex]) {
                printf("Need fluid but it's suppressed!\n");
                break;
            }
            if(!silminState->nSolidCoexist[windex]) {
                printf("Adding fluid to initial guess assemblage\n");
            silminState->solidComp[windex][0] = 0.0;
            }
            incSolids[windex] = 1;
            silminState->nSolidCoexist[windex] = 1;
            silminState->solidComp[windex][1] = CO2;
            CO2 -= CO2;
        }

        /* Total minerals */
        for (i=0, silminState->solidComp[olindex][0]=0.0; i<solids[olindex].na; i++)
            silminState->solidComp[olindex][0] += silminState->solidComp[olindex+i+1][0];
        for (i=0, silminState->solidComp[gtindex][0]=0.0; i<solids[gtindex].na; i++)
            silminState->solidComp[gtindex][0] += silminState->solidComp[gtindex+i+1][0];
        for (i=0, silminState->solidComp[opxindex][0]=0.0, silminState->solidComp[cpxindex][0]=0.0; i<solids[cpxindex].na; i++) {
            silminState->solidComp[opxindex][0] += silminState->solidComp[opxindex+i+1][0];
            silminState->solidComp[cpxindex][0] += silminState->solidComp[cpxindex+i+1][0];
        }
        for (i=0, silminState->solidComp[plindex][0]=0.0; i<solids[plindex].na; i++)
            silminState->solidComp[plindex][0] += silminState->solidComp[plindex+i+1][0];
        for (i=0, silminState->solidComp[spindex][0]=0.0; i<solids[spindex].na; i++)
            silminState->solidComp[spindex][0] += silminState->solidComp[spindex+i+1][0];

        /* Check sums */
        if (SiO2 != 0.0) {
            printf("Silica error of %f\n", SiO2);
        }
        if (MgO != 0.0) {
            printf("Magnesia error of %f\n", MgO);
        }

        if(!olivine && (SiO2 != 0.0 || MgO != 0.0)) {
            printf("Adding olivine to the assemblage.\n");
            olivine = silminState->nSolidCoexist[olindex] = 1;
            continue;
        }

        /* check minerals */
        for (i=0; i<npc; i++) {
            for (j=0; j<silminState->nSolidCoexist[i]; j++) {
                if (solids[i].na == 1) {
                    if (silminState->solidComp[i][j] < 0.0) break;
                } else {
                    double *mSol = calloc((size_t) (solids[i].na), sizeof(double));
                    for (k=0;k<solids[i].na;k++) mSol[k] = silminState->solidComp[i+1+k][j];
                    if (silminState->solidComp[i][j] < 0.0) {
                        free(mSol);
                        break;
                    }
                    if (!(*solids[i].test)(SIXTH, 1000.0, 1.0, 0, 0, NULL, NULL, NULL, mSol)) {
                        free(mSol);
                        break;
                    }
                    free(mSol);
                }
            }
        }

        /* no liquid allowed */
        if (hasLiquid) {
            for (j=1; j<silminState->nLiquidCoexist; j++) {
                free((silminState->liquidComp )[j]); (silminState->liquidComp )[j] = NULL;
                free((silminState->liquidDelta)[j]); (silminState->liquidDelta)[j] = NULL;
            }
            for (i=0;i<nlc;i++) silminState->liquidComp[0][i] = 0.0;
            //for (i=0;i<nc;i++)  silminState->dspLiquidComp[i] = 0.0;
        }
        silminState->nLiquidCoexist = 1; // else memory leak
        silminState->liquidMass = 0.0;

        success = TRUE;

    }

    free(incSolids);
    return success;

}
