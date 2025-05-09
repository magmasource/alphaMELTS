const char *liquid_v34_ver(void) { return "$Id: liquid_v34.c,v 1.3 2006/10/20 00:59:22 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: liquid_v34.c,v $
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
MELTS Source Code: RCS Revision 1.1.1.1  2004/01/02 19:21:49  cvsaccount
MELTS Source Code: RCS CTserver University of Chicago
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2002/04/06 00:51:39  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
*/

/*
**++
**  FACILITY:  Thread Safe Silicate Melts Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Routines to compute liquid solution properties
**      (file: LIQUID_V34.C)
**--
*/

#include "silmin.h"  /* Structure definitions for SILMIN package */
#ifdef USESEH
#include <windows.h>
void raise_sigabrt(DWORD dwType);
extern int doInterrupt;
#else
#include <signal.h>
#endif

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))

#include "param_struct_data_v34.h"

#define NR   18             /* Number of independent mole fraction variables */
#define NA   19             /* Number of liquid components                   */
#define SMX  SHRT_MAX

/*
 * Array to convert W(i,j) indexes to entries in the array of structures
 * modelParameters[k], where k is:
 */

static const short int Index[NA][NA] = {
 { SMX,   0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,  16,  17 },
 {   0, SMX,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,  32,  33,  34 },
 {   1,  18, SMX,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50 },
 {   2,  19,  35, SMX,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,  65 },
 {   3,  20,  36,  51, SMX,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79 },
 {   4,  21,  37,  52,  66, SMX,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92 },
 {   5,  22,  38,  53,  67,  80, SMX,  93,  94,  95,  96,  97,  98,  99, 100, 101, 102, 103, 104 },
 {   6,  23,  39,  54,  68,  81,  93, SMX, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115 },
 {   7,  24,  40,  55,  69,  82,  94, 105, SMX, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125 },
 {   8,  25,  41,  56,  70,  83,  95, 106, 116, SMX, 126, 127, 128, 129, 130, 131, 132, 133, 134 },
 {   9,  26,  42,  57,  71,  84,  96, 107, 117, 126, SMX, 135, 136, 137, 138, 139, 140, 141, 142 },
 {  10,  27,  43,  58,  72,  85,  97, 108, 118, 127, 135, SMX, 143, 144, 145, 146, 147, 148, 149 },
 {  11,  28,  44,  59,  73,  86,  98, 109, 119, 128, 136, 143, SMX, 150, 151, 152, 153, 154, 155 },
 {  12,  29,  45,  60,  74,  87,  99, 110, 120, 129, 137, 144, 150, SMX, 156, 157, 158, 159, 160 },
 {  13,  30,  46,  61,  75,  88, 100, 111, 121, 130, 138, 145, 151, 156, SMX, 161, 162, 163, 164 },
 {  14,  31,  47,  62,  76,  89, 101, 112, 122, 131, 139, 146, 152, 157, 161, SMX, 165, 166, 167 },
 {  15,  32,  48,  63,  77,  90, 102, 113, 123, 132, 140, 147, 153, 158, 162, 165, SMX, 168, 169 },
 {  16,  33,  49,  64,  78,  91, 103, 114, 124, 133, 141, 148, 154, 159, 163, 166, 168, SMX, 170 },
 {  17,  34,  50,  65,  79,  92, 104, 115, 125, 134, 142, 149, 155, 160, 164, 167, 169, 170, SMX }
};

#undef SMX
#define WH(i,j) ((calculationMode == MODE__MELTS) ? meltsModelParameters[Index[i][j]].enthalpy : pMeltsModelParameters[Index[i][j]].enthalpy )
#define WS(i,j) ((calculationMode == MODE__MELTS) ? meltsModelParameters[Index[i][j]].entropy  : pMeltsModelParameters[Index[i][j]].entropy  )
#define WV(i,j) ((calculationMode == MODE__MELTS) ? meltsModelParameters[Index[i][j]].volume   : pMeltsModelParameters[Index[i][j]].volume   )

/* possible types */

#define ENTHALPY 0
#define ENTROPY  1
#define VOLUME   2

double outputParamLiq(int type, int index1, int index2) {
    double value = 0.0;

    // Will have full test on type and set here (or may just use type and index)...
    if      (type == ENTHALPY) value = WH(index1,index2);
    else if (type == ENTROPY)  value = WS(index1,index2);
    else if (type == VOLUME)   value = WV(index1,index2);

    return value;
}

int changeParamLiq(int type, int index1, int index2, double value) {

    if      (type == ENTHALPY) meltsModelParameters[Index[index1][index2]].enthalpy = value;
    else if (type == ENTROPY)  meltsModelParameters[Index[index1][index2]].entropy = value;
    else if (type == VOLUME)   meltsModelParameters[Index[index1][index2]].volume = value;

    return TRUE;
}

/*
 *=============================================================================
 * Public functions:
 *    inpMask -  bitwise mask for specifying input parameters
 *    outMask -  bitwise mask for selecting output
 *    mask    -  bitwise mask for selecting output
 *    t       -  Temperature (K)
 *    p       -  Pressure (bars)
 *    *r      -  (pointer to x[]) Array of independent compositional variables
 */

void
conLiq_v34(int inpMask, int outMask, double t, double p,
    double *o,      /* comp of liquid in moles of oxides                        */
    double *m,      /* comp of liquid in moles of endmember components          */
    double *r,      /* comp of liquid in terms of the independent comp var      */
    double *x,      /* comp of liquid in mole fractions of endmember components */
    double **dm,    /* Jacobian matrix: dm[i][j] = dr[i]/dm[j]                  */
    double ***d2m,  /* vector of matrices: d2m[i][j][k] = d2r[i]/dm[j]dm[k]     */
    double *logfo2) /* base 10 logarithm of the oxygen fugacity                 */
{
    /*---------------------------------------------------------------------------
    Not all combinations of inpMask and outMask are feasible. Valid
        combinations are:

       inpMask          outMask
    (1)  FIRST | SEVENTH  FIRST
    (2)  FIRST            SEVENTH
    (3)  SECOND           THIRD | FOURTH | FIFTH | SIXTH
    (4)  THIRD            FOURTH

    (1) converts a vector of moles of oxides into a vector of moles of oxides
            with the correct redox state for the given t, p, and logfo2. Note that
            the original vector is used as output.
    (2) calculates from a vector of moles of oxides and the given t and p, the
            appropriate logfo2
    (3) calculates from a vector of moles of endmember components, one or
            all of: r[], x[], dr[]/dm[], or d2r[]/dm[]dm[]
    (4) calculates from a vector of independent compositional variables
            mole fractions of endmember components
    ----------------------------------------------------------------------------*/

    int i, j, k;

    if (inpMask & FIRST) {
        /*-------------------------------------------------------------------------
            Oxide+logfo2 -> oxide or oxide -> logfo2. The algorithm used is that
            given by:
            Kress, VC and Carmichael, ISE (1991) The compressibility of silicate
                liquids containing Fe2O3 and the effect of composition, temperature,
                oxygen fugacity and pressure on their redox states.
                Contributions to Mineralogy and Petrology (in press)
            Coefficients for the oxides are initialized in LIQ_STRUCT_DATA.H
        --------------------------------------------------------------------------*/
        static const double t0 = 1673.15,                      /* K       */
                         a =    0.196,
                         b =    1.1492e4,                  /* K       */
                         c =   -6.675,
                         e =   -3.364,
                         f =   -7.01e-7  * 1.0e5,          /* K/bar   */
                         g =   -1.54e-10 * 1.0e5,          /* 1/bar   */
                         h =    3.85e-17 * 1.0e5 * 1.0e5;  /* K/bar^2 */
        double sum = 0.0, temp;
        int indexFeO = -1, indexFe2O3 = -1;

        for (i=0; i<NA; i++) {
            if (bulkSystem[i].type == FEO)   indexFeO   = i;
            if (bulkSystem[i].type == FE2O3) indexFe2O3 = i;
        }
        if (indexFeO == -1 || indexFe2O3 == -1) {
       printf("Fatal error in conLiq_v34 (LIQUID_V34.C)\n");
       printf("The oxides FeO and Fe2O3 cannot be identified.\n");
       return;
        }

        if (inpMask == (FIRST | SEVENTH)  && outMask == FIRST) {
            /*----------------------------------------------------------------------
                Converts a vector of moles of oxides (as defined in LIQ_STRUCT_DATA.H
                for the structure bulkSystem) into a vector of moles of oxides with
                the correct redox state (ferric/ferrous ratio) for the given bulk
                composition, t and p.
            ------------------------------------------------------------------------*/

            o[indexFeO]   += 2.0*o[indexFe2O3];
            o[indexFe2O3]  = 0.0;
            if (o[indexFeO] == 0.0) return;

            for (i=0; i<NA; i++) sum += o[i];
            if (sum == 0.0) return;

            temp = a*log(10.0)*(*logfo2) + b/t + c + e*(1.0-t0/t - log(t/t0))
           + f*p/t + g*(t-t0)*p/t + h*SQUARE(p)/t;
            for (i=0; i<NA; i++) temp += bulkSystem[i].coeff*o[i]/sum;
            temp = exp(temp);

            o[indexFe2O3]  = temp*o[indexFeO]/(1.0 + 2.0*temp);
            o[indexFeO]   -= 2.0*o[indexFe2O3];

        } else if (inpMask == FIRST && outMask == SEVENTH) {
            /*----------------------------------------------------------------------
                Calculates from the given t and p and a vector of moles of oxides
                (as defined in LIQ_STRUCT_DATA.H for the structure bulkSystem) the
                appropriate log10fo2 for the given t and p.
            ------------------------------------------------------------------------*/

            if (o[indexFeO] == 0.0 || o[indexFe2O3] == 0.0) { *logfo2 = 0.0; return; }
            for (i=0; i<NA; i++) sum += o[i];
            sum += o[indexFe2O3];
            if (sum == 0.0) { *logfo2 = 0.0; return; }
            temp = b/t + c + e*(1.0-t0/t - log(t/t0)) + f*p/t + g*(t-t0)*p/t
           + h*SQUARE(p)/t;
            for (i=0; i<NA; i++) temp += bulkSystem[i].coeff*o[i]/sum;
            temp += 2.0*bulkSystem[indexFeO].coeff*o[indexFe2O3]/sum
                                - bulkSystem[indexFe2O3].coeff*o[indexFe2O3]/sum;
            *logfo2 = (log(o[indexFe2O3]/o[indexFeO]) - temp)/(a*log(10.0));

        } else
            printf("Illegal call to conLiq_v34 with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

    } else if (inpMask == SECOND) {
        double sum;

        if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH))
            printf("Illegal call to conLiq_v34 with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        for (i=0, sum=0.0; i<NA; i++) sum += m[i];

        if (outMask & THIRD) {
            /* Converts a vector of moles of end-member components (m) into a vector
         of independent compositional variables (r) required as input for the
         remaining public functions.
         The dependent variable is taken to be SiO2 (1st component), as this
         component will never have a mole fraction of zero.                   */

            for (i=0; i<NR; i++) r[i] = (sum != 0.0) ? m[i+1]/sum : 0.0;
        }

        if (outMask & FOURTH) {
            /* Converts a vector of moles of end-member components (m) into a vector
         of mole fractions of endmember components                            */

            for (i=0; i<NA; i++) x[i] = (sum != 0.0) ? m[i]/sum : 0.0;
        }

        if (outMask & FIFTH) {
            /* Calculates the matrix dr[i]/dm[j] using m[] as input                 */

            if (sum == 0.0) {
                for (i=0; i<NR; i++) { for (j=0; j<NA; j++) dm[i][j] = 0.0; }
            } else {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++)
                        dm[i][j] = (i+1 == j) ? (1.0-m[i+1]/sum)/sum : - m[i+1]/SQUARE(sum);
                }
            }
        }

        if (outMask & SIXTH) {
            /* Calculates the matrix d2r[i]/dm[j]dm[k] using m[] as input           */

            if (sum == 0.0) {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++)  {
                        for (k=0; k<NA; k++) d2m[i][j][k] = 0.0;
                    }
                }
            } else {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++)  {
                        for (k=0; k<NA; k++) {
                            d2m[i][j][k]  = 2.0*m[i+1]/CUBE(sum);
                            d2m[i][j][k] -= (i+1 == j) ? 1.0/SQUARE(sum) : 0.0;
                            d2m[i][j][k] -= (i+1 == k) ? 1.0/SQUARE(sum) : 0.0;
                        }
                    }
                }
            }

        }

    } else if (inpMask == THIRD && outMask == FOURTH) {
   /* Converts a vector of independent compositional variables (r)
            into a vector of mole fractions of end-member components (x)            */

        for (i=0, x[0] = 1.0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }

    } else {
        printf("Illegal call to conLiq_v34 with inpMask = %o and outMask = %o\n",
            inpMask, outMask);
    }

}

int
testLiq_v34(int mask, double t, double p,
    int na,          /* Expected number of endmember components                 */
    int nr,          /* Expected number of independent compositional variables  */
    char **names,    /* array of strings of names of endmember oxides           */
    char **formulas, /* array of strings of formulas of endmember components    */
    double *r,       /* array of indepependent compos variables, check bounds   */
    double *m)       /* array of moles of endmember components, check bounds    */
{
    const char *phase = "liquid.c";
    const char *NAMES[NA]    = { "SiO2"  , "TiO2"  , "Al2O3" , "Fe2O3" , "Cr2O3" ,
                               "FeO"   , "MnO"   , "MgO"   , "NiO"   , "CoO"   ,
                               "CaO"   , "Na2O"  , "K2O"   , "P2O5"  , "H2O"   ,
                               "CO2"   , "SO3"   , "Cl2O-1", "F2O-1" };
    const char *FORMULAS[NA] = { "SiO2"     , "TiO2"     , "Al2O3"    , "Fe2O3"    ,
                               "MgCr2O4"  , "Fe2SiO4"  , "MnSi0.5O2", "Mg2SiO4"  ,
                               "NiSi0.5O2", "CoSi0.5O2", "CaSiO3"   , "Na2SiO3"  ,
                               "KAlSiO4"  , "Ca3(PO4)2", "CO2"      , "SO3"      ,
                               "Cl2O-1"   , "F2O-1"    , "H2O"      };
    int result = TRUE, i;
    double sum;

    if (mask & FIRST) {
        result = result && (na == NA);
        if (!result) printf("<<%s>> Wrong number of components!\n", phase);
    }
    if (mask & SECOND) {
        result = result && (nr == NR);
        if (!result) printf("<<%s>> Wrong number of indep variables!\n", phase);
    }
    if (mask & THIRD) {
        for (i=0; i<NA; i++) {
            result = result && (strcmp(names[i],NAMES[i]) == 0);
            if (!result)
                printf("<<%s>> Oxide[%d] should be %s not %s.\n",
                    phase, i, NAMES[i], names[i]);
        }
    }
    if (mask & FOURTH) {
        for (i=0; i<NA; i++) {
            result = result && (strcmp(formulas[i],FORMULAS[i]) == 0);
            if (!result)
                printf("<<%s>> Component[%d] should have formula %s not %s.\n",
                    phase, i, FORMULAS[i], formulas[i]);
        }
    }
    /* Check bounds on the independent compositional variables */
    if (mask & FIFTH) {
        for (i=0, sum=0.0; i<NR; i++) {
            result = result && (r[i] >= 0.0) && (r[i] <= 1.0);
            sum += r[i];
        }
        result = result && (sum <= 1.0);
    }
    /* Check bounds on moles of endmember components */
    if (mask & SIXTH) {
        for (i=0; i<NA; i++) result = result && (m[i] >= 0.0);
    }

    return result;
}

void
dispLiq_v34(int mask, double t, double p, double *x,
    char **formula            /* Mineral formula for interface display MASK: 1 */
    )
{
    double *r = x;
    char *string = (char *) malloc((unsigned) (7+NA*12+1)*sizeof(char));
    (void) snprintf(string, 8, "wt%% ox:");

    if (mask & FIRST) {    /* assume maximum string length is 5 */
        const char *NAMES[NA]    = { "SiO2"  , "TiO2"  , "Al2O3" , "Fe2O3" , "Cr2O3" ,
                                 "FeO"   , "MnO"   , "MgO"   , "NiO"   , "CoO"   ,
                                 "CaO"   , "Na2O"  , "K2O"   , "P2O5"  , "H2O"   ,
                                 "CO2"   , "S"     , "Cl"    , "F"     };
        double m[NA], oxVal[NA], oxSum;
        int i, j, n;

        for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }
        for (i=0, oxSum=0.0; i<NA; i++) {
            for (j=0, oxVal[i]=0.0; j<NA; j++) oxVal[i] += m[j]*(liquid[j].liqToOx)[i];
            oxVal[i] *= bulkSystem[i].mw;
            oxSum    += oxVal[i];
        }

        if (oxSum != 0.0) for (i=0, n=7; i<NA; i++)
            if (oxVal[i] != 0.0) {
                double w = 100.0*oxVal[i]/oxSum;
                int nn = snprintf(&string[n], 13, " %s %.2f", NAMES[i], w);
                n += (nn < 13) ? nn : 12;
            }

        *formula = string;
    }
}

#ifdef PHMELTS_ADJUSTMENTS
#define W(i,j) ((WH(i,j)-t*WS(i,j)+(p-1.0)*WV(i,j)))
/* all in terms of mu - mu0 */
void muH2OLiq(int mask, double t, double p, double *m,
    double *muH2O,/* muH2O     = mu*H2O                 BINARY MASK: 0000000001 */
    double *dm,   /* dm[i]     = d mu*H2O/dm[i]         BINARY MASK: 0000000010 */
    double *dt,   /* dt        = d mu*H2O/d T           BINARY MASK: 0000000100 */
    double *dp,   /* dp        = d mu*H2O/d P           BINARY MASK: 0000001000 */
    double **d2m, /* d2m[i][j] = d mu*H2O/dm[i][j]      BINARY MASK: 0000010000 */
    double *d2mt, /* d2mt[i]   = d mu*H2O/dm[i]dt       BINARY MASK: 0000100000 */
    double *d2mp, /* d2mp[i]   = d mu*H2O/dm[i]dp       BINARY MASK: 0001000000 */
    double *d2t2, /* d2t2      = d mu*H2O/dt2           BINARY MASK: 0010000000 */
    double *d2tp, /* d2tp      = d mu*H2O/dtdp          BINARY MASK: 0100000000 */
    double *d2p2) /* d2p2      = d mu*H2O/dp2           BINARY MASK: 1000000000 */
{
    double *a = (double *) malloc((size_t) NA*sizeof(double));
    double *mu = (double *) malloc((size_t) NA*sizeof(double));
    double **dadr = (double **) malloc((size_t) NA*sizeof (double *));
    double **drdm = (double **) malloc((size_t) NA*sizeof (double *));
    double dadm, sum;
    int indH2O, i, j, k, l;
    double *r = (double *) malloc(NR*sizeof(double));
    int iw = NA-1; /* Only valid for Rhyolite-MELTS 1.0.2 and pMELTS */

    for (i=0; i<NA; i++) {
        dadr[i] = (double *)  malloc((size_t) NR*sizeof (double));
        drdm[i] = (double *)  malloc((size_t) NA*sizeof (double));
    }

    for (indH2O=0; indH2O<nlc; indH2O++) if (!strcmp(liquid[indH2O].label, "H2O")) break;
    conLiq(SECOND, THIRD | FIFTH, t, p, (double *) NULL, m, r, (double *) NULL, drdm, (double ***) NULL, (double *) NULL);
    actLiq(FIRST | SECOND | THIRD, t, p, r, a, mu, dadr, NULL);

    for (i=0,sum=0.0; i<NA; i++) sum += m[i];
    if (mask & FIRST) {
        *muH2O = mu[indH2O];
    }
    if (mask & SECOND) {
        for (i=0; i<NA; i++) {
            for (j=0, dadm=0.0; j<NR; j++) dadm += dadr[indH2O][j] * drdm[j][i];
            dm[i] = R*t*dadm/a[indH2O];
        }

    }
    if (mask & FIFTH) {
        double general;

        for (i=0,sum=0.0; i<NA; i++) sum += m[i];
        for (i=0,general=0.0; i<NA; i++)
            for (j=i+1; j<NA; j++) general -= 6.0*W(i,j)*m[i]*m[j];
        for (k=0; k<NA; k++) {
            for (l=0; l<NA; l++) {
                d2m[k][l] = general - sum*sum*((k==iw?0.0:W(k,iw)) + (l==iw?0.0:W(l,iw)) + (k==l?0.0:W(k,l)));
                d2m[k][l] += sum*sum*2.0*R*t*(1.0-((k==iw && l==iw)?sum*sum/(m[iw]*m[iw]):0.0));
                for (i=0;i<NA;i++) d2m[k][l]+=2.0*sum*m[i]*((i==iw?0.0:W(i,iw)) + (i==k?0.0:W(i,k)) + (i==l?0.0:W(i,l)));
                d2m[k][l] /= sum*sum*sum*sum;
            }
        }
    }
    free(a); free(mu); free(r);
    for (i=0; i<NA; i++) {  free(dadr[i]); free(drdm[i]);  }
    free(dadr); free(drdm);
}
#undef W
#endif

void
actLiq_v34(int mask, double t, double p, double *r,
    double *a,  /* (pointer to a[]) activities           BINARY MASK: 001 */
    double *mu, /* (pointer to mu[]) chemical potentials BINARY MASK: 010 */
    double **dr /* (pointer to dr[][]) d(a[])/d(x[])     BINARY MASK: 100 */
    )
{
    double x[NA], gex;
    int i, j;

    /* x[0] --> x[NA-1] is an array of mole fractions of liquid components      */
    for (x[0] = 1.0, i=0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }
    for (i=0; i<NA; i++) if (x[i] < 0.0)
        printf("CAUTION: Liquid component %s has negative mole fraction.\n",
            liquid[i].label);
    for (gex=0.0, i=0; i<NA; i++) {
        for (j=i+1; j<NA; j++) gex += x[i]*x[j]*(WH(i,j)-t*WS(i,j)+(p-1.0)*WV(i,j));
    }

    if (mask & FIRST) {
        for (i=0; i<NA; i++) {
            a[i] = - gex;
            for (j=0;   j<i;  j++) a[i] += (WH(i,j)-t*WS(i,j)+(p-1.0)*WV(i,j))*x[j];
            for (j=i+1; j<NA; j++) a[i] += (WH(i,j)-t*WS(i,j)+(p-1.0)*WV(i,j))*x[j];

            a[i] = (x[i] != 0.0) ? x[i]*exp(a[i]/(R*t)) : 0.0;
            a[i] = (i != NA-1)   ? (1.0 - x[NA-1])*a[i] : x[NA-1]*a[NA-1];
        }
    }

    if (mask & SECOND) {
        for (i=0; i<NA; i++) {
            mu[i] = - gex;
            for (j=0;   j<i;  j++) mu[i] += (WH(i,j)-t*WS(i,j)+(p-1.0)*WV(i,j))*x[j];
            for (j=i+1; j<NA; j++) mu[i] += (WH(i,j)-t*WS(i,j)+(p-1.0)*WV(i,j))*x[j];
#ifdef USESEH
            if (x[i] < 0.0) { doInterrupt = TRUE; raise_sigabrt(EXCEPTION_FLT_INVALID_OPERATION); }
#else
            if (x[i] < 0.0) (void) raise(SIGABRT);
#endif
            mu[i] = (x[i] != 0.0) ? mu[i] + R*t*log(x[i]) : 0.0;
            if (i != NA-1)           mu[i]    += R*t*log(1.0-x[NA-1]);
            else if (x[NA-1] != 0.0) mu[NA-1] += R*t*log(x[NA-1]);
        }
    }

    if (mask & THIRD) {
        double a[NA], dgexdr[NR];

        for (i=0; i<NA; i++) {
            a[i] = - gex;
            for (j=0;   j<i;  j++) a[i] += (WH(i,j)-t*WS(i,j)+(p-1.0)*WV(i,j))*x[j];
            for (j=i+1; j<NA; j++) a[i] += (WH(i,j)-t*WS(i,j)+(p-1.0)*WV(i,j))*x[j];

            a[i] = (x[i] != 0.0) ? x[i]*exp(a[i]/(R*t)) : 0.0;
            a[i] = (i != NA-1)   ? (1.0 - x[NA-1])*a[i] : x[NA-1]*a[NA-1];
        }

        for (i=0; i<NR; i++) {
            for (dgexdr[i] = x[0]*(WH(0,i+1)-t*WS(0,i+1)+(p-1.0)*WV(0,i+1)), j=0;
                j<NR; j++) dgexdr[i] += (i != j) ? r[j]*((WH(i+1,j+1)-t*WS(i+1,j+1)
                +(p-1.0)*WV(i+1,j+1)) - (WH(0,j+1)-t*WS(0,j+1)+(p-1.0)*WV(0,j+1)))
                : - r[j]*(WH(0,j+1)-t*WS(0,j+1)+(p-1.0)*WV(0,j+1));
        }

        /* Special case for component 0 (SiO2)                                    */
        for (j=0; j<NR; j++) {
            dr[0][j] = - R*t/x[0] + (WH(0,j+1)-t*WS(0,j+1)+(p-1.0)*WV(0,j+1)) - dgexdr[j];
            if (j == NR-1) dr[0][j] += - R*t/(1.0-x[NA-1]);

            dr[0][j] *= a[0]/(R*t);
        }

        /* All other cases                                                        */
        for (i=1; i<NA; i++) {
            for (j=0; j<NR; j++) {
                dr[i][j] = (i == (j+1)) ? ((r[j] != 0.0) ? R*t/r[j] : 0.0) :
                    WH(i,j+1)-t*WS(i,j+1)+(p-1.0)*WV(i,j+1);
                dr[i][j] += - (WH(0,i)-t*WS(0,i)+(p-1.0)*WV(0,i)) - dgexdr[j];
                if (i != NA-1 && j == NR-1) dr[i][NR-1] += - R*t/(1.0-x[NA-1]);
                else if (i == NA-1 && j == NR-1 && x[NA-1] != 0.0)
                    dr[NA-1][NR-1] += R*t/x[NA-1];

                dr[i][j] *= a[i]/(R*t);
            }
        }

    }

    if (mask & FOURTH) {
        for (i=0; i<NA; i++) {
            a[i] = x[i];
            a[i] = (i != NA-1)   ? (1.0 - x[NA-1])*a[i] : x[NA-1]*a[NA-1];
        }
    }

}

void
gmixLiq_v34(int mask, double t, double p, double *r,
    double *gmix, /* Gibbs energy of mixing             BINARY MASK: 001 */
    double *dr,   /* (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 010 */
    double **dr2  /* (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 100 */
    )
{
    double x[NA];
    int i, j;

    /* x[0] --> x[NA-1] is an array of mole fractions of liquid components      */
    for (x[0] = 1.0, i=0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }
    for (i=0; i<NA; i++) if (x[i] < 0.0)
        printf("CAUTION: Liquid component %s has negative mole fraction.\n",
            liquid[i].label);

    if (mask & FIRST) {
        for (*gmix = 0.0, i=0; i<NA; i++) {
            for (j=i+1; j<NA; j++) *gmix += x[i]*x[j]*(WH(i,j) - t*WS(i,j)
                                                                        + (p-1.0)*WV(i,j));
#ifdef USESEH
            if (x[i] < 0.0) { doInterrupt = TRUE; raise_sigabrt(EXCEPTION_FLT_INVALID_OPERATION); }
#else
            if (x[i] < 0.0) (void) raise(SIGABRT);
#endif
            *gmix += (x[i] != 0.0) ? R*t*x[i]*log(x[i]) : 0.0;
        }
        *gmix += (x[NA-1] != 0.0) ?
            R*t*(x[NA-1]*log(x[NA-1]) + (1.0-x[NA-1])*log(1.0-x[NA-1])) : 0.0;
    }

    if(mask & SECOND) {
        for (i=0; i<NR; i++) {
#ifdef USESEH
            if ((x[0] < 0.0) || (r[i] < 0.0)) { doInterrupt = TRUE; raise_sigabrt(EXCEPTION_FLT_INVALID_OPERATION); }
#else
            if ((x[0] < 0.0) || (r[i] < 0.0)) (void) raise(SIGABRT);
#endif
            dr[i] = (r[i] != 0.0) ? R*t*(log(r[i]) - log(x[0])) : 0.0;
            dr[i] += x[0]*(WH(0,i+1)-t*WS(0,i+1)+(p-1.0)*WV(0,i+1));
            for (j=0; j<NR; j++)
                dr[i] += (i != j) ? r[j]*((WH(i+1,j+1) - t*WS(i+1,j+1)
               + (p-1.0)*WV(i+1,j+1))
               - (WH(0,j+1)-t*WS(0,j+1)+(p-1.0)*WV(0,j+1)))
               : - r[j]*(WH(0,j+1)-t*WS(0,j+1)+(p-1.0)*WV(0,j+1));
        }
        dr[NR-1] += (x[NA-1] != 0.0) ? R*t*(log(x[NA-1])-log(1.0-x[NA-1])) : 0.0;

    }

    if(mask & THIRD) {
        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dr2[i][j] = R*t/x[0] - (WH(0,i+1)-t*WS(0,i+1)+(p-1.0)*WV(0,i+1))
                    - (WH(0,j+1)-t*WS(0,j+1)+(p-1.0)*WV(0,j+1));
                dr2[i][j] += (i != j) ? (WH(i+1,j+1)-t*WS(i+1,j+1)+(p-1.0)*WV(i+1,j+1))
                   : 0.0;
            }
            dr2[i][i] += (r[i] != 0.0) ? R*t/r[i] : 0.0;
        }
        dr2[NR-1][NR-1] +=
            (x[NA-1] != 0.0) ? R*t*(1.0/x[NA-1] + 1.0/(1.0-x[NA-1])) : 0.0;

    }
}

void
hmixLiq_v34(int mask, double t, double p, double *r,
    double *hmix /* Enthalpy of mixing BINARY MASK: 1 */
    )
{
    double x[NA];
    int i, j;

    /* x[0] --> x[NA-1] is an array of mole fractions of liquid components      */
    for (x[0] = 1.0, i=0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }
    for (i=0; i<NA; i++) if (x[i] < 0.0)
        printf("CAUTION: Liquid component %s has negative mole fraction.\n",
            liquid[i].label);

    for (*hmix = 0.0, i=0; i<NA; i++) {
        for (j=i+1; j<NA; j++) *hmix += x[i]*x[j]*(WH(i,j)+(p-1.0)*WV(i,j));
    }

}

void
smixLiq_v34(int mask, double t, double p, double *r,
    double *smix, /* Entropy of mixing                  BINARY MASK: 001 */
    double *dr,   /* (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010 */
    double **dr2  /* (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100 */
    )
{
    double x[NA];
    int i, j;

    /* x[0] --> x[NA-1] is an array of mole fractions of liquid components      */
    for (x[0] = 1.0, i=0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }
    for (i=0; i<NA; i++) if (x[i] < 0.0)
        printf("CAUTION: Liquid component %s has negative mole fraction.\n",
            liquid[i].label);

    if (mask & FIRST) {
        for (*smix = 0.0, i=0; i<NA; i++) {
            for (j=i+1; j<NA; j++) *smix += x[i]*x[j]*WS(i,j);
#ifdef USESEH
            if (x[i] < 0.0) { doInterrupt = TRUE; raise_sigabrt(EXCEPTION_FLT_INVALID_OPERATION); }
#else
            if (x[i] < 0.0) (void) raise(SIGABRT);
#endif
            *smix += (x[i] != 0.0) ? - R*x[i]*log(x[i]) : 0.0;
        }
        *smix += (x[NA-1] != 0.0) ?
            -R*(x[NA-1]*log(x[NA-1]) + (1.0-x[NA-1])*log(1.0-x[NA-1])) : 0.0;
    }

    if(mask & SECOND) {
        for (i=0; i<NR; i++) {
#ifdef USESEH
            if ((x[0] < 0.0) || (r[i] < 0.0)) { doInterrupt = TRUE; raise_sigabrt(EXCEPTION_FLT_INVALID_OPERATION); }
#else
            if ((x[0] < 0.0) || (r[i] < 0.0)) (void) raise(SIGABRT);
#endif
            dr[i] = (r[i] != 0.0) ? R*(log(x[0]) - log(r[i])) : 0.0;
            dr[i] += x[0]*WS(0,i+1);
            for (j=0; j<NR; j++)
                dr[i] += (i != j) ? r[j]*(WS(i+1,j+1) - WS(0,j+1)) : - r[j]*WS(0,j+1);
        }
        dr[NR-1] += (x[NA-1] != 0.0) ? -R*(log(x[NA-1])-log(1.0-x[NA-1])) : 0.0;
    }

    if(mask & THIRD) {
        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dr2[i][j] = -R/x[0] - WS(0,i+1) - WS(0,j+1);
                dr2[i][j] += (i != j) ? WS(i+1,j+1) : 0.0;
            }
            dr2[i][i] += (r[i] != 0.0) ? - R/r[i] : 0.0;
        }
        dr2[NR-1][NR-1] +=
            (x[NA-1] != 0.0) ? -R*(1.0/x[NA-1] + 1.0/(1.0-x[NA-1])) : 0.0;
    }
}

void
cpmixLiq_v34(int mask, double t, double p, double *r,
    double *cpmix, /* Heat capacity of mixing BINARY MASK: 001 */
    double *dt,    /* d(cp)/d(t)              BINARY MASK: 010 */
    double *dr     /* d(cp)/d(x[])            BINARY MASK: 100 */
    )
{
    double x[NA];
    int i;

    /* x[0] --> x[NA-1] is an array of mole fractions of liquid components      */
    for (x[0] = 1.0, i=0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }
    for (i=0; i<NA; i++) if (x[i] < 0.0)
        printf("CAUTION: Liquid component %s has negative mole fraction.\n",
            liquid[i].label);

    if (mask & FIRST) {
        *cpmix = 0.0;
    }

    if(mask & SECOND) {
        *dt = 0.0;
    }

    if(mask & THIRD) {
        for(i=0; i<NR; i++) dr[i] = 0.0;
    }
}

void
vmixLiq_v34(int mask, double t, double p, double *r,
    double *vmix, /* Volume of mixing                BINARY MASK: 0000000001 */
    double *dr,   /* (pointer to dx[]) d(v)/d(x[])   BINARY MASK: 0000000010 */
    double **dr2, /* (point to dx2[][]) d(v)/d(x[])2 BINARY MASK: 0000000100 */
    double *dt,   /* d(v)/d(t)                       BINARY MASK: 0000001000 */
    double *dp,   /* d(v)/d(p)                       BINARY MASK: 0000010000 */
    double *dt2,  /* d2(v)/d(t)2                     BINARY MASK: 0000100000 */
    double *dtdp, /* d2(v)/d(t)d(p)                  BINARY MASK: 0001000000 */
    double *dp2,  /* d2(v)/d(p)2                     BINARY MASK: 0010000000 */
    double *drdt, /* d2(v)/d(x[])d(t)                BINARY MASK: 0100000000 */
    double *drdp  /* d2(v)/d(x[])d(p)                BINARY MASK: 1000000000 */
    )
{
    double x[NA];
    int i, j;

    /* x[0] --> x[NA-1] is an array of mole fractions of liquid components      */
    for (x[0] = 1.0, i=0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }
    for (i=0; i<NA; i++) if (x[i] < 0.0)
        printf("CAUTION: Liquid component %s has negative mole fraction.\n",
            liquid[i].label);

    if (mask & FIRST) {
        for (*vmix = 0.0, i=0; i<NA; i++)
            for (j=i+1; j<NA; j++) *vmix += x[i]*x[j]*WV(i,j);
    }

    if(mask & SECOND) {
        for (i=0; i<NR; i++) {
            dr[i] = x[0]*WV(0,i+1);
            for (j=0; j<NR; j++)
                dr[i] += (i != j) ? r[j]*(WV(i+1,j+1) - WV(0,j+1)) : - r[j]*WV(0,j+1);
        }
    }

    if(mask & THIRD) {
        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dr2[i][j]  = - WV(0,i+1) - WV(0,j+1);
                dr2[i][j] += (i != j) ? WV(i+1,j+1) : 0.0;
            }
        }
    }

    if(mask & FOURTH) {
        *dt = 0.0;
    }

    if(mask & FIFTH) {
        *dp = 0.0;
    }

    if(mask & SIXTH) {
        *dt2 = 0.0;
    }

    if(mask & SEVENTH) {
        *dtdp = 0.0;
    }

    if(mask & EIGHTH) {
        *dp2 = 0.0;
    }

    if(mask & NINTH) {
        for (i=0; i<NR; i++) drdt[i] = 0.0;
    }

    if(mask & TENTH) {
        for (i=0; i<NR; i++) drdp[i] = 0.0;
    }
}

/* ============================================================================
   In the following public routine:
   m = m[i] = moles of the ith component in the liquid, and
   mu*O2    = mu O2 - mu0 O2, defined from the vector o[]
   ==========================================================================*/

void
muO2Liq_v34(int mask, double t, double p, double *m,
    double *muO2, /* muO2      = mu*O2                 BINARY MASK: 0000000001 */
    double *dm,   /* dm[i]     = d mu*O2/dm[i]         BINARY MASK: 0000000010 */
    double *dt,   /* dt        = d mu*O2/d T           BINARY MASK: 0000000100 */
    double *dp,   /* dp        = d mu*O2/d P           BINARY MASK: 0000001000 */
    double **d2m, /* d2m[i][j] = d mu*O2/dm[i][j]      BINARY MASK: 0000010000 */
    double *d2mt, /* d2mt[i]   = d mu*O2/dm[i]dt       BINARY MASK: 0000100000 */
    double *d2mp, /* d2mp[i]   = d mu*O2/dm[i]dp       BINARY MASK: 0001000000 */
    double *d2t2, /* d2t2      = d mu*O2/dt2           BINARY MASK: 0010000000 */
    double *d2tp, /* d2tp      = d mu*O2/dtdp          BINARY MASK: 0100000000 */
    double *d2p2) /* d2p2      = d mu*O2/dp2           BINARY MASK: 1000000000 */
{
    /*-------------------------------------------------------------------------
        The algorithm used is that given by:
        Kress, VC and Carmichael, ISE (1991) The compressibility of silicate
            liquids containing Fe2O3 and the effect of composition, temperature,
            oxygen fugacity and pressure on their redox states.
            Contributions to Mineralogy and Petrology (in press)
        Coefficients for the oxides are initialized in LIQ_STRUCT_DATA.H
    --------------------------------------------------------------------------*/
    static const double t0 = 1673.15,                      /* K       */
                       a =    0.196,
                       b =    1.1492e4,                  /* K       */
                       c =   -6.675,
                       e =   -3.364,
                       f =   -7.01e-7  * 1.0e5,          /* K/bar   */
                       g =   -1.54e-10 * 1.0e5,          /* 1/bar   */
                       h =    3.85e-17 * 1.0e5 * 1.0e5;  /* K/bar^2 */
    double mOx[NA], vTemp[NA], mTemp[NA][NA], sum;
    int indexFeO = -1, indexFe2O3 = -1, i, j;

    for (i=0; i<NA; i++) {
        if (bulkSystem[i].type == FEO)   indexFeO	= i;
        if (bulkSystem[i].type == FE2O3) indexFe2O3 = i;
    }
    if (indexFeO == -1 || indexFe2O3 == -1) {
     printf("Fatal error in muO2Liq_v34 (LIQUID.C)\n");
     printf("The oxides FeO and Fe2O3 cannot be identified.\n");
     return;
    }

    for (i=0; i<NA; i++)
        for (j=0, mOx[i]=0.0; j<nlc; j++) mOx[i] += (liquid[j].liqToOx)[i]*m[j];

    for (i=0, sum=0.0; i<NA; i++) { sum += mOx[i]; } sum += mOx[indexFe2O3];
    if (sum == 0.0 || mOx[indexFeO] == 0.0 || mOx[indexFe2O3] == 0.0) {
        if (mask & FIRST)   *muO2 = 0.0;
        if (mask & SECOND)  for (i=0; i<nlc; i++) dm[i] = 0.0;
        if (mask & THIRD)   *dt = 0.0;
        if (mask & FOURTH)  *dp = 0.0;
        if (mask & FIFTH)   {
            for (i=0; i<nlc; i++) for (j=0; j<nlc; j++) d2m[i][j] = 0.0;
        }
        if (mask & SIXTH)   for (i=0; i<nlc; i++) d2mt[i] = 0.0;
        if (mask & SEVENTH) for (i=0; i<nlc; i++) d2mp[i] = 0.0;
        if (mask & EIGHTH)  *d2t2 = 0.0;
        if (mask & NINTH)   *d2tp = 0.0;
        if (mask & TENTH)   *d2p2 = 0.0;
        return;
    }

    /*-------------------------------------------------------------------------*/

    if (mask & FIRST) {
        double temp;
        temp = b/t + c + e*(1.0-t0/t - log(t/t0)) + f*p/t + g*(t-t0)*p/t
         + h*SQUARE(p)/t;
        for (i=0; i<NA; i++) temp += bulkSystem[i].coeff*mOx[i]/sum;
        temp += 2.0*bulkSystem[indexFeO].coeff*mOx[indexFe2O3]/sum
                    - bulkSystem[indexFe2O3].coeff*mOx[indexFe2O3]/sum;
        *muO2 = R*t*(log(mOx[indexFe2O3]/mOx[indexFeO]) - temp)/a;
    }

    if (mask & SECOND) {
        for (j=0; j<NA; j++) {
            double factor = (j == indexFe2O3) ? 2.0 : 1.0;
            for (i=0, vTemp[j]=0.0; i<NA; i++)
                vTemp[j] -= (i == j) ?
                    bulkSystem[i].coeff*(1.0-factor*mOx[i]/sum)/sum :
                    - bulkSystem[i].coeff*factor*mOx[i]/SQUARE(sum);
            vTemp[j] += - (factor*mOx[indexFe2O3]/SQUARE(sum))
                *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
            if      (j == indexFeO) vTemp[j] += -1.0/mOx[indexFeO];
            else if (j == indexFe2O3) {
                vTemp[j] += 1.0/mOx[indexFe2O3];
                vTemp[j] += (1.0/sum)
                    *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
            }
            vTemp[j] *= R*t/a;
        }
        for (i=0; i<nlc; i++)
            for (j=0, dm[i]=0.0; j<NA; j++) dm[i] += vTemp[j]*(liquid[i].liqToOx)[j];
    }

    if (mask & THIRD) {
        double temp;
        temp = b/t + c + e*(1.0-t0/t - log(t/t0)) + f*p/t + g*(t-t0)*p/t
         + h*SQUARE(p)/t;
        for (i=0; i<NA; i++) temp += bulkSystem[i].coeff*mOx[i]/sum;
        temp += 2.0*bulkSystem[indexFeO].coeff*mOx[indexFe2O3]/sum
                    - bulkSystem[indexFe2O3].coeff*mOx[indexFe2O3]/sum;
        *dt = R*(log(mOx[indexFe2O3]/mOx[indexFeO]) - temp)/a
                + R*t*(b/SQUARE(t) - e*(t0/t-1.0)*(1.0/t) + f*p/SQUARE(t)
                - g*(t0/t)*(p/t) + h*SQUARE(p/t))/a;
    }

    if (mask & FOURTH) {
        *dp = R*t*(-f/t - g*(t-t0)/t - 2.0*h*p/t)/a;
    }

    if (mask & FIFTH) {
        int k, l;
        for (k=0; k<NA; k++) {
            double factorK = (k == indexFe2O3) ? 2.0 : 1.0;
            for (j=0; j<NA; j++) {
                double factorJ = (j == indexFe2O3) ? 2.0 : 1.0;
                for (i=0, mTemp[k][j]=0.0; i<NA; i++) {
                    mTemp[k][j] -=
                        2.0*factorJ*factorK*bulkSystem[i].coeff*mOx[i]/CUBE(sum);
                    if (i == j) mTemp[k][j] -= - factorK*bulkSystem[i].coeff/SQUARE(sum);
                    if (i == k) mTemp[k][j] -= - factorJ*bulkSystem[i].coeff/SQUARE(sum);
                }
                mTemp[k][j] += 2.0*(factorJ*factorK*mOx[indexFe2O3]/CUBE(sum))
                    *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
                if      (j == indexFeO && k == indexFeO)
                    mTemp[k][j] += 1.0/SQUARE(mOx[indexFeO]);
                else if (j == indexFe2O3 && k == indexFe2O3) {
                    mTemp[k][j] += -1.0/SQUARE(mOx[indexFe2O3]);
                    mTemp[k][j] += -((factorJ+factorK)/SQUARE(sum))
                        *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
                }
                else if (j == indexFe2O3) mTemp[k][j] += -(factorK/SQUARE(sum))
                        *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
                else if (k == indexFe2O3) mTemp[k][j] += -(factorJ/SQUARE(sum))
                        *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
                mTemp[k][j] *= R*t/a;
            }
        }
        for (i=0; i<nlc; i++) for (j=0; j<nlc; j++)
            for (k=0, d2m[i][j]=0.0; k<NA; k++) for (l=0; l<NA; l++)
                d2m[i][j] += mTemp[k][l]*(liquid[i].liqToOx)[k]*(liquid[j].liqToOx)[l];
    }

    if (mask & SIXTH) {
        for (j=0; j<NA; j++) {
            double factor = (j == indexFe2O3) ? 2.0 : 1.0;
            for (i=0, vTemp[j]=0.0; i<NA; i++)
                vTemp[j] -= (i == j) ?
                    bulkSystem[i].coeff*(1.0-factor*mOx[i]/sum)/sum :
                    - bulkSystem[i].coeff*factor*mOx[i]/SQUARE(sum);
            vTemp[j] += - (factor*mOx[indexFe2O3]/SQUARE(sum))
                *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
            if      (j == indexFeO) vTemp[j] += -1.0/mOx[indexFeO];
            else if (j == indexFe2O3) {
                vTemp[j] += 1.0/mOx[indexFe2O3];
                vTemp[j] += (1.0/sum)
                    *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
            }
            vTemp[j] *= R/a;
        }
        for (i=0; i<nlc; i++) for (j=0, d2mt[i]=0.0; j<NA; j++)
            d2mt[i] += vTemp[j]*(liquid[i].liqToOx)[j];
    }

    if (mask & SEVENTH) {
        for (i=0; i<nlc; i++) d2mp[i] = 0.0;
    }

    if (mask & EIGHTH) {

        *d2t2 = 2.0*R*(b/SQUARE(t) - e*(t0/t-1.0)*(1.0/t) + f*p/SQUARE(t)
               - g*(t0/t)*(p/t) + h*SQUARE(p/t))/a
                    + R*t*(-2.0*b/CUBE(t) - e*(1.0 - 2.0*t0/t)/SQUARE(t)
                        - 2.0*f*p/CUBE(t) + 2.0*g*t0*p/CUBE(t) - 2.0*h*SQUARE(p/t)/t)/a;
    }

    if (mask & NINTH) {
        *d2tp = -R*(f/t + g*(t-t0)/t + 2.0*h*p/t)/a
                    + R*t*(f/SQUARE(t) - g*(t0/t)/t + 2.0*h*p/SQUARE(t))/a;
    }

    if (mask & TENTH) {
        *d2p2 = - (R*t*2.0*h/t)/a;
    }

}

void
visLiq_v34(int mask, double t, double p, double *r,
    double *viscosity  /* log(10) viscosity            BINARY MASK: 00000001 */
    )
{
    double coeff[NA], factor[NA], m[NA], x[NA], sum;
    int nSiO2 = -1, i, j;

    struct _shawModel {
        char *oxide;
        double     coeff;
        double     factor;
    } shawModel[] = {
        { "TiO2",	4.5, 1.0 }, { "Al2O3",  6.7, 2.0 },
        { "Fe2O3",  3.4, 2.0 }, { "FeO",	3.4, 1.0 },
        { "MgO",	3.4, 1.0 }, { "CaO",	4.5, 1.0 },
        { "Na2O",	2.8, 1.0 }, { "K2O",	2.8, 1.0 },
        { "H2O",	2.0, 1.0 }
    };
    const int nShaw = (sizeof shawModel / sizeof(struct _shawModel));

    for (j=0; j<NA; j++) { coeff[j] = 0.0; factor[j] = 0.0; }
    for (i=0; i<nShaw; i++) {
        for (j=0; j<NA; j++) if (strcmp(shawModel[i].oxide, bulkSystem[j].label) == 0) {
            coeff[j]  = shawModel[i].coeff;
            factor[j] = shawModel[i].factor;
            break;
        }
    }
    for (i=0; i<NA; i++) if (strcmp("SiO2", bulkSystem[i].label) == 0) { nSiO2 = i; break; }

    if (nSiO2 == -1) { *viscosity = 0.0; return; }

    /* m[0] --> m[NA-1] is an array of mole fractions of liquid components      */
    for (m[0]=1.0, i=0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* convert m[] -> x[] : mole fractions of liquid comp -> moles of oxides    */
    for (i=0; i<NA; i++)
        for (j=0, x[i]=0.0; j<NA; j++) x[i] += (liquid[j].liqToOx)[i]*m[j];

    /* Convert to the Shaw mole fractions                                       */
    for (i=0, sum=0.0; i<NA; i++) {
        if (factor[i] > 0.0) x[i] *= factor[i];
        sum += x[i];
    }
    for (i=0; i<NA; i++) x[i] /= (sum != 0.0) ? sum : 1.0;

    if (mask & FIRST) {
        for (i=0, *viscosity=0.0; i<NA; i++) *viscosity += coeff[i]*x[nSiO2]*x[i];
        *viscosity /= (x[nSiO2] < 1.0) ? 1.0 - x[nSiO2] : 1.0;
        *viscosity  = (*viscosity)*(10000.0/t - 1.50)  - 6.40;
        *viscosity /= log(10.0);
    }
}

/* end of file LIQUID_V34.C */
