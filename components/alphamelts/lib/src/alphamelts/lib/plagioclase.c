const char *plagioclase_ver(void) { return "$Id: plagioclase.c,v 1.3 2006/10/20 00:59:22 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: plagioclase.c,v $
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
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:58:12  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
 * Revision 3.7  1997/06/21  22:49:57  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.6  1997/05/03  20:23:36  ghiorso
 * *** empty log message ***
 *
 * Revision 3.5  1997/03/27  17:03:39  ghiorso
 * *** empty log message ***
 *
 * Revision 3.4  1996/09/24  20:33:42  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 3.3  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 3.2  1995/09/09  23:37:35  ghiorso
 * Modifications by Asimow to include new Derivatives of gmix and convert
 * options for liquid-absent fO2 buffering. Also new olivine components
 * implemented.
 *
 * Revision 3.2  1995/09/09  23:37:35  ghiorso
 * Modifications by Asimow to include new Derivatives of gmix and convert
 * options for liquid-absent fO2 buffering. Also new olivine components
 * implemented.
 *
 * Revision 3.1  1995/08/18  17:31:17  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Routines to compute plagioclase solution properties
**      (file: PLAGIOCLASE.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  June 3, 1991   Original Version - started
**      V1.0-2  Mark S. Ghiorso  June 19, 1991
**              Added testPlg and revised conPlg
**      V1.0-3  Mark S. Ghiorso  July 5, 1991
**              Added exclusion criteria option to actPlg
**      V1.0-4  Mark S. Ghiorso  August 26, 1991
**              Extended conPlg to convert indep comp var -> mole frac of comp
**      V1.0-5  Mark S. Ghiorso  September 16, 1991
**              Altered parameter list and structure of conPlg in order
**              to provide derivatives of the independent compositional
**              variables with respect to moles of endmember components.
**              Added several options and tests to conPlg.
**      V1.0-6  Mark S. Ghiorso  September 24, 1991
**              (1) Altered parameter list and structure of conPlg in order
**                  to provide derivatives of the endmember mole fractions with
**                  respect to independent compositional variables.
**              (2) Altered parameter list and structure of testPlg to
**                  provide for bounds testing on both the independent
**                  compositional variables and on endmember moles
**      V1.0-7  Mark S. Ghiorso  October 30, 1991
**              Corrected definition of enthalpy of mixing
**      V1.0-8  Mark S. Ghiorso  January 28, 1992
**              Altered initialization of automatic arrays to conform to
**              ANSI standard (DEC C RISC compiler)
**      V1.0-9  Mark S. Ghiorso  March 3, 1992
**              Corrected local initialization of xab, xan, xor to account for
**              zero concentrations of endmembers
**      V1.0-10 Mark S. Ghiorso  March 23, 1992
**              Removed references to the pow() function and added macros
**              for SQUARE and CUBE
**      V1.1-1  Mark S. Ghiorso  July 13, 1992
**              Added (*display) function
**      V2.0-1  Mark S. Ghiorso  May 10, 1994
**              (1) Began modifications for new isenthalpic, isentropic,
**                  isochoric derivatives
**      V2.0-2  Mark S. Ghiorso  May 17, 1994
**              (1) Changed definition of the enthalpy of mixing
**      V2.0-3  Mark S. Ghiorso  March 2, 1995
**              Implemented suggestions made by Paul D. Asimow (CalTech)
**              (1) changed != to > in xab, xan, xor calculations
**	V3.0-1  Paul D. Asimow  July 31, 1995
**		(1) added d3rdm3 to (*display) and d3gdr3 to (*gmix)
**
*/

#include "silmin.h"  /* Structure definitions foor SILMIN package */

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

/*
 *=============================================================================
 * Plagioclase solution parameters:
 * Elkins, Linda T., Grove, Timothy L.
 * Ternary plagioclase experiments and thermodynamic models
 * American Mineralogist 75, 544-559
 */
#ifndef TESTDYNAMICLIB
// Compare 33818.5 on sanidine
#define PLAG_ADJUSTMENT   -66200.0
#endif

#ifndef TESTDYNAMICLIB
static double whabor   = -6498.5;  /* joules     */
#else
static double whabor   = 18810.0;  /* joules     */
#endif
static double wsabor   = 10.3;     /* joules/K   */
static double wvabor   = 0.4602;   /* joules/bar */
#ifndef TESTDYNAMICLIB
static double whorab   = -6498.5;  /* joules     */
#else
static double whorab   = 27320.0;  /* joules     */
#endif
static double wsorab   = 10.3;     /* joules/K   */
static double wvorab   = 0.3264;   /* joules/bar */
static double whaban   = 7924.0;   /* joules     */
static double whanab   = 0.0;      /* joules     */
#ifndef TESTDYNAMICLIB
static double whoran   =  6498.5;  /* joules     */
static double whanor   =  6498.5;  /* joules     */
#else
static double whoran   = 40317.0;  /* joules     */
static double whanor   = 38974.0;  /* joules     */
#endif
static double wvanor   = -0.1037;  /* joules/bar */
#ifndef TESTDYNAMICLIB
static double whabanor = 12545.0+(23065.0+39645.5+PLAG_ADJUSTMENT);  /* joules     */
#else
static double whabanor = 12545.0;  /* joules     */
#endif
static double wvabanor = -1.095;   /* joules/bar */

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conPlg defines the conversion from m[i], to r[j]
 */
#define NR         2
#define NS         0
#define NA         3
#define FR0(i)     (i == 0) ? 1.0 - xab : - xab
#define FR1(i)     (i == 1) ? 1.0 - xan : - xan
#define DFR0DR0(i) - 1.0
#define DFR1DR1(i) - 1.0

#define LOG(x) ((x > 0.0) ? log(x) : log(DBL_EPSILON))

/*
 * Global (to this file): derivative definitions
 */
#define R       8.3143
#define S       - R*(xab*LOG(xab) + xan*LOG(xan) + xor*LOG(xor)) + \
                                wsabor*xab*xor*(xor+xan/2.0) + wsorab*xab*xor*(xab+xan/2.0)
#define H       whaban*xab*xan*(xan+xor/2.0) + whanab*xab*xan*(xab+xor/2.0) + \
                                whabor*xab*xor*(xor+xan/2.0) + whorab*xab*xor*(xab+xan/2.0) + \
                                whanor*xan*xor*(xor+xab/2.0) + whoran*xan*xor*(xan+xab/2.0) + \
                                whabanor*xab*xan*xor
#define V       wvabor*xab*xor*(xor+xan/2.0) + wvorab*xab*xor*(xab+xan/2.0) + \
                                wvanor*xan*xor*(xor+xab/2.0) + wvabanor*xab*xan*xor
#define G       H - t*(S) + (p-1.0)*(V)

#define DGDR0   R*t*(LOG(xab) - LOG(xor)) + \
                                whaban*(xan*(xan+xor/2.0) - 0.5*xab*xan) + \
                                whanab*(xan*(xab+xor/2.0) + 0.5*xab*xan) + \
                                (whabor-t*wsabor+(p-1.0)*wvabor)* \
                   ((xor-xab)*(xor+xan/2.0) - xab*xor) + \
                                (whorab-t*wsorab+(p-1.0)*wvorab)* \
                   ((xor-xab)*(xab+xan/2.0) + xab*xor) + \
                                (whanor+(p-1.0)*wvanor)*(- xan*(xor+xab/2.0) - 0.5*xan*xor) + \
                                whoran*(- xan*(xan+xab/2.0) + 0.5*xan*xor) + \
                                (whabanor+(p-1.0)*wvabanor)*xan*(xor - xab)
#define DGDR1   R*t*(LOG(xan) - LOG(xor)) + \
                                whaban*(xab*(xan+xor/2.0) + 0.5*xab*xan) + \
                                whanab*(xab*(xab+xor/2.0) - 0.5*xab*xan) + \
                                (whabor-t*wsabor+(p-1.0)*wvabor)* \
                   (- xab*(xor+xan/2.0) - 0.5*xab*xor) + \
                                (whorab-t*wsorab+(p-1.0)*wvorab)* \
                   (- xab*(xab+xan/2.0) + 0.5*xab*xor) + \
                                (whanor+(p-1.0)*wvanor)*((xor-xan)*(xor+xab/2.0) - xan*xor) + \
                                whoran*((xor-xan)*(xan+xab/2.0) + xan*xor) + \
                                (whabanor+(p-1.0)*wvabanor)*xab*(xor - xan)
#define DGDP    (V)

#define D2GDR0R0  R*t*(1.0/xab + 1.0/xor) + (whanab - whaban)*xan \
                                    - 2.0*(whabor-t*wsabor+(p-1.0)*wvabor)* \
                     ((xor+xan/2.0) + (xor-xab)) \
                                    - 2.0*(whorab-t*wsorab+(p-1.0)*wvorab)* \
                     ((xab+xan/2.0) - (xor-xab)) + \
                                    (whanor+(p-1.0)*wvanor - whoran)*xan + \
                                    -2.0*(whabanor+(p-1.0)*wvabanor)*xan
#define D2GDR0R1  R*t/xor + \
                                    whaban*((xan+xor/2.0) + 0.5*xan - 0.5*xab) + \
                                    whanab*((xab+xor/2.0) - 0.5*xan + 0.5*xab) + \
                                    0.5*(whabor-t*wsabor+(p-1.0)*wvabor)*(3.0*xab-xan-3.0*xor) + \
                                    0.5*(whorab-t*wsorab+(p-1.0)*wvorab)*(xor - 5.0*xab - xan) + \
                                    0.5*(whanor+(p-1.0)*wvanor)*(3.0*xan - 3.0*xor - xab) + \
                                    0.5*whoran*(xor - 5.0*xan - xab) + \
                                    (whabanor+(p-1.0)*wvabanor)*(xor - xab - xan)
#define D2GDR0DT  R*(LOG(xab) - LOG(xor)) - \
                                    wsabor*((xor-xab)*(xor+xan/2.0) - xab*xor) - \
                                    wsorab*((xor-xab)*(xab+xan/2.0) + xab*xor)
#define D2GDR0DP  wvabor*((xor-xab)*(xor+xan/2.0) - xab*xor) + \
                                    wvorab*((xor-xab)*(xab+xan/2.0) + xab*xor) + \
                                    wvanor*(- xan*(xor+xab/2.0) - 0.5*xan*xor) + \
                                    wvabanor*xan*(xor - xab)

#define D2GDR1R1  R*t*(1.0/xan + 1.0/xor) + xab*(whaban - whanab) + \
                                    (whabor-t*wsabor+(p-1.0)*wvabor)*xab - \
                                    (whorab-t*wsorab+(p-1.0)*wvorab)*xab + \
                                    (whanor+(p-1.0)*wvanor)*(2.0*xan - 4.0*xor - xab) + \
                                    whoran*(2.0*xor - 4.0*xan - xab) - \
                                    2.0*(whabanor+(p-1.0)*wvabanor)*xab
#define D2GDR1DT  R*(LOG(xan) - LOG(xor)) - \
                                    wsabor*(- xab*(xor+xan/2.0) - 0.5*xab*xor) - \
                                    wsorab*(- xab*(xab+xan/2.0) + 0.5*xab*xor)
#define D2GDR1DP  wvabor*(- xab*(xor+xan/2.0) - 0.5*xab*xor) + \
                                    wvorab*(- xab*(xab+xan/2.0) + 0.5*xab*xor) + \
                                    wvanor*((xor-xan)*(xor+xab/2.0) - xan*xor) + \
                                    wvabanor*xab*(xor - xan)
#define D2GDT2    0.0
#define D2GDTDP   0.0
#define D2GDP2    0.0

#define D3GDR0R0R0 R*t*(1.0/SQUARE(xor)-1.0/SQUARE(xab)) + 6.0* \
		   (whabor-t*wsabor+(p-1)*wvabor) - 6.0*(whorab - \
		   t*wsorab+(p-1)*wvorab)
#define D3GDR0R0R1 R*t*(1.0/SQUARE(xor)) + (whanab-whaban) + 3.0* \
		   (whabor-t*wsabor+(p-1)*wvabor) - 3.0*(whorab - \
		   t*wsorab+(p-1)*wvorab) + (whanor+(p-1)*wvanor - \
		   whoran) - 2*(whabanor + (p-1)*wvabanor)
#define D3GDR0R0DT R*(1.0/xab + 1.0/xor) + \
                   2.0*wsabor*((xor+xan/2.0) + (xor-xab)) + \
                   2.0*wsorab*((xab+xan/2.0) - (xor-xab))
#define D3GDR0R0DP - 2.0*wvabor*((xor+xan/2.0) + (xor-xab)) - \
                   2.0*wvorab*((xab+xan/2.0) - (xor-xab)) + \
                   wvanor*xan - 2.0*wvabanor*xan
#define D3GDR0R1DT R/xor - 0.5*wsabor*(3.0*xab - xan - 3.0*xor) - \
                   0.5*wsorab*(xor - 5.0*xab - xan)
#define D3GDR0R1DP 0.5*wvabor*(3.0*xab - xan - 3.0*xor) + \
                   0.5*wvorab*(xor - 5.0*xab - xan) + \
                   0.5*wvanor*(3.0*xan - 3.0*xor - xab) + \
                   wvabanor*(xor - xab - xan)
#define D3GDR1R1R1 R*t*(1.0/SQUARE(xor)-1.0/SQUARE(xan)) + 6.0* \
		   (whanor+(p-1.0)*wvanor) - 6.0*whoran
#define D3GDR1R1R0 R*t*(1.0/SQUARE(xor)) + (whaban-whanab) + \
		   (whabor-t*wsabor+(p-1)*wvabor) - (whorab - \
		   t*wsorab+(p-1)*wvorab) + 3.0*(whanor+(p-1.0)*wvanor) \
		   -3.0*whoran - 2*(whabanor + (p-1)*wvabanor)
#define D3GDR1R1DT R*(1.0/xan + 1.0/xor) + (wsorab-wsabor)*xab
#define D3GDR1R1DP (wvabor - wvorab)*xab + \
                   wvanor*(2.0*xan - 4.0*xor - xab) - 2.0*wvabanor*xab
#define D3GDT3     0.0
#define D3GDT2DP   0.0
#define D3GDTDP2   0.0
#define D3GDP3     0.0

#define D3GDR0DT2  0.0
#define D3GDR0DTDP 0.0
#define D3GDR0DP2  0.0
#define D3GDR1DT2  0.0
#define D3GDR1DTDP 0.0
#define D3GDR1DP2  0.0

/*
 *=============================================================================
 * Public functions:
 *    mask  -  bitwise mask for selecting output
 *    t     -  Temperature (K)
 *    p     -  Pressure (bars)
 *    *x    -  (pointer to x[]) Array of independent compositional variables
 */
int
testPlg(int mask, double t, double p,
    int na,          /* Expected number of endmember components                 */
    int nr,          /* Expected number of independent compositional variables  */
    char **names,    /* array of strings of names of endmember components       */
    char **formulas, /* array of strings of formulas of endmember components    */
    double *r,       /* array of indepependent compos variables, check bounds   */
    double *m)       /* array of moles of endmember components, check bounds    */
{
    const char *phase = "plagioclase.c";
    const char *NAMES[NA]    = { "albite", "anorthite", "sanidine" };
    const char *FORMULAS[NA] = { "NaAlSi3O8", "CaAl2Si2O8", "KAlSi3O8" };
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
                printf("<<%s>> Component[%d] should be %s not %s.\n",
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
conPlg(int inpMask, int outMask, double t, double p,
    double *e,      /* comp of plagioclase in moles of elements                    */
    double *m,      /* comp of plagioclase in moles of endmember components        */
    double *r,      /* comp of plagioclase in terms of the independent comp var    */
    double *x,      /* comp of plagioclase in mole fractions of endmember comp     */
    double **dm,    /* Jacobian matrix: dm[i][j] = dr[i]/dm[j]                  */
    double ***d2m,  /* vector of matrices: d2m[i][j][k] = d2r[i]/dm[j]dm[k]     */
    double **dr,    /* Jacobian matrix: dr[i][j] = dx[i]/dr[j]                  */
    double ****d3m) /* 3rd deriv matrix: d3m[i][j][k][l]=d3r[i]/dm[j]dm[k]dm[l] */
{
    /*---------------------------------------------------------------------------
    Not all combinations of inpMask and outMask are feasible. Valid
        combinations are:

       inpMask          outMask
    (1)  FIRST            SECOND
    (2)  SECOND           THIRD  | FOURTH  | FIFTH | SIXTH | EIGHTH
    (3)  THIRD            FOURTH | SEVENTH

    (1) converts a vector of moles of elements into a vector of moles of
            endmember plagioclase components.
    (2) calculates from a vector of moles of endmember components, one or
            all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
    (3) calculates from a vector of independent compositional variables
            mole fractions of endmember components and/or the Jacobian matrix
            dx[]/dr[]

    In this routine it is assumed that the elements are in the order of atomic
    numbers and that the order of plagioclase components has been verified as:
                m[0] = albite    (NaAlSi3O8) ,
                m[1] = anorthite (CaAl2Si2O8) and
                m[2] = sanidine  (KAlSi3O8)

    ----------------------------------------------------------------------------*/

    int i, j, k;

    if (inpMask == FIRST && outMask == SECOND) {
        static const int Na = 11;
        static const int K  = 19;
        static const int Ca = 20;

                                    /* Projection into the Na, Ca, K  triangle */
        m[0] = e[Na]; /* moles of NaAlSi3O8                      */
        m[1] = e[Ca]; /* Moles of CaAl2Si2O8                     */
        m[2] = e[K];  /* Moles of KAlSi3O8                       */

    } else if (inpMask == SECOND) {
        double sum;

        if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
            printf("Illegal call to conPlg with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        for (i=0, sum=0.0; i<NA; i++) sum += m[i];

        if (outMask & THIRD) {
            for (i=0; i<NR; i++) r[i] = (sum != 0.0) ? m[i]/sum : 0.0;
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
                        dm[i][j] = (i == j) ? (1.0-m[i]/sum)/sum : - m[i]/SQUARE(sum);
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
                            d2m[i][j][k]  = 2.0*m[i]/CUBE(sum);
                            d2m[i][j][k] -= (i == j) ? 1.0/SQUARE(sum) : 0.0;
                            d2m[i][j][k] -= (i == k) ? 1.0/SQUARE(sum) : 0.0;
                        }
                    }
                }
            }

        }

        if (outMask & EIGHTH) {
            /* calculates the matrix d3r[i]/dm[j]dm[k]dm[l] using m[] as input        */
            int l;

            if (sum == 0.0) {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++)  {
                        for (k=0; k<NA; k++)  {
                            for (l=0; l<NA; l++) d3m[i][j][k][l] = 0.0;
                        }
                    }
                }
            } else {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++)  {
                        for (k=0; k<NA; k++)  {
                            for (l=0; l<NA; l++)  {
                                d3m[i][j][k][l]  = -6.0*m[i]/QUARTIC(sum);
                                d3m[i][j][k][l] += (i == j) ? 2.0/CUBE(sum) : 0.0;
                                d3m[i][j][k][l] += (i == k) ? 2.0/CUBE(sum) : 0.0;
                                d3m[i][j][k][l] += (i == l) ? 2.0/CUBE(sum) : 0.0;
                            }
                        }
                    }
                }
            }

        }

    } else if (inpMask == THIRD) {

        if (outMask & ~(FOURTH | SEVENTH))
            printf("Illegal call to conPlg with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        if (outMask & FOURTH) {
            /* Converts a vector of independent compositional variables (r)
         into a vector of mole fractions of endmember components (x).         */

            for (i=0, x[2]=1.0; i<NR; i++) { x[i] = r[i]; x[2] -= r[i]; }
            if (fabs(x[2]) < sqrt(DBL_EPSILON)) x[2] = 0.0;
        }

        if (outMask & SEVENTH) {
            /* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
            for (i=0; i<NR; i++) for (j=0; j<NR; j++) dr[i][j] = (i == j) ? 1.0 : 0.0;
            for (j=0; j<NR; j++) dr[2][j] = -1.0;
        }

    } else  {
        printf("Illegal call to conPlg with inpMask = %o and outMask = %o\n",
            inpMask, outMask);
    }

}

void
dispPlg(int mask, double t, double p, double *x,
    char **formula            /* Mineral formula for interface display MASK: 1 */
    )
{
    double *r = x;
    static char masterString[] = {
/*             1111111111222222222233333333334444444444555555555566666666667
     01234567890123456789012345678901234567890123456789012345678901234567890 */
        "K_.__Na_.__Ca_.__Al_.__Si_.__O8" };

    if (mask & FIRST) {
        char *string = strcpy((char *) malloc((size_t) (strlen(masterString)+1)*sizeof(char)), masterString);
        double totAl, totCa, totNa, totK, totSi;
        char n[5];
        int i;

        totK   = 1.0 - r[0] - r[1];
        totNa  = r[0];
        totCa  = r[1];
        totAl  = 1.0 + r[1];
        totSi  = 3.0 - r[1];

        (void) snprintf(n, 5, "%4.2f", totK);   for (i=0; i<4; i++) string[ 1+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totNa);  for (i=0; i<4; i++) string[ 7+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totCa);  for (i=0; i<4; i++) string[13+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totAl);  for (i=0; i<4; i++) string[19+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totSi);  for (i=0; i<4; i++) string[25+i] = n[i];

        *formula = string;
    }
}

void
actPlg(int mask, double t, double p, double *x,
    double *a,  /* (pointer to a[]) activities              BINARY MASK: 0001 */
    double *mu, /* (pointer to mu[]) chemical potentials    BINARY MASK: 0010 */
    double **dx /* (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100 */
    )           /* exclusion criteria applied to results if BINARY MASK: 1000 */
{
    double xab = (x[0]          >  DBL_EPSILON) ? x[0]          : DBL_EPSILON;
    double xan = (x[1]          >  DBL_EPSILON) ? x[1]          : DBL_EPSILON;
    double xor = (1.0-x[0]-x[1] >  DBL_EPSILON) ? 1.0-x[0]-x[1] : DBL_EPSILON;

    double g, dgdr[NR], fr[NA][NR];
    int i, j;

    for(i=0; i<NA; i++) {
     fr[i][0] = FR0(i);
     fr[i][1] = FR1(i);
    }

    g       = G;
    dgdr[0] = DGDR0;
    dgdr[1] = DGDR1;

    /* activities for library */
    if (!mask && a != NULL) {
        for(i=0; i<NA; i++) {
       for (a[i]=g, j=0; j<NR; j++) a[i] += fr[i][j]*dgdr[j];
       a[i] = exp(a[i]/(R*t));
        }
    }

    if (mask & FIRST) {
        for(i=0; i<NA; i++) {
       for (a[i]=g, j=0; j<NR; j++) a[i] += fr[i][j]*dgdr[j];
       a[i] = exp(a[i]/(R*t));
        }
    }

    if (mask & SECOND) {
        for(i=0; i<NA; i++) {
       for (mu[i]=g, j=0; j<NR; j++) mu[i] += fr[i][j]*dgdr[j];
        }
    }

    if (mask & THIRD) {
        double d2gdr2[NR][NR], dfrdr[NA][NR], sum;
        int k;

        d2gdr2[0][0] = D2GDR0R0;
        d2gdr2[0][1] = D2GDR0R1;
        d2gdr2[1][0] = d2gdr2[0][1];
        d2gdr2[1][1] = D2GDR1R1;

        for(i=0; i<NA; i++) {
       dfrdr[i][0] = DFR0DR0(i);
       dfrdr[i][1] = DFR1DR1(i);
        }

        for (i=0; i<NA; i++) {
       for (k=0; k<NR; k++) {
                    for (dx[i][k]=g, j=0; j<NR; j++) dx[i][k] += fr[i][j]*dgdr[j];
                    dx[i][k] = exp(dx[i][k]/(R*t));
                    sum = (1.0+dfrdr[i][k])*dgdr[k];
                    for (j=0; j<NR; j++) sum += fr[i][j]*d2gdr2[j][k];
                    dx[i][k] *= sum/(R*t);
       }
        }
    }

    if (mask & FOURTH) {
        /* implement exclusion criteria on quantities for preclb routines  */
        static const double exclusion[NA] = {
       0.05,  /* exclusion criteria on the mole fraction of albite     */
       0.05,  /* exclusion criteria on the mole fraction of anorthite  */
       0.05   /* exclusion criteria on the mole fraction of orthoclase */
        };
        double x[NA];

        x[0] = xab; /* mole fraction of albite                             */
        x[1] = xan; /* mole fraction of anorthite                          */
        x[2] = xor; /* mole fraction of orthoclase                         */

        for (i=0; i<NA; i++) {
            if (x[i] < exclusion[i]) {
                if (mask & FIRST)  a[i]  = 0.0;
                if (mask & SECOND) mu[i] = 0.0;
                if (mask & THIRD)  for (j=0; j<NR; j++) dx[i][j] = 0.0;
            }
        }
    }

}

void
gmixPlg(int mask, double t, double p, double *x,
    double *gmix, /* Gibbs energy of mixing             BINARY MASK: 0001 */
    double *dx,   /* (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010 */
    double **dx2, /* (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100 */
    double ***dx3 /* (pointer to dx3[][][]) d3(g)/d(x[])3 NARY MASK: 1000 */
    )
{
    double xab = (x[0]          >  DBL_EPSILON) ? x[0]          : DBL_EPSILON;
    double xan = (x[1]          >  DBL_EPSILON) ? x[1]          : DBL_EPSILON;
    double xor = (1.0-x[0]-x[1] >  DBL_EPSILON) ? 1.0-x[0]-x[1] : DBL_EPSILON;

    if (mask & FIRST) {
        *gmix = G;
    }

    if(mask & SECOND) {
        dx[0] = DGDR0;
        dx[1] = DGDR1;
    }

    if(mask & THIRD) {
        double d2gdr2[NR][NR];
        int i, j;

        d2gdr2[0][0] = D2GDR0R0;
        d2gdr2[0][1] = D2GDR0R1;
        d2gdr2[1][0] = d2gdr2[0][1];
        d2gdr2[1][1] = D2GDR1R1;

        for (i=0; i<NR; i++) {
       for (j=0; j<NR; j++) dx2[i][j] = d2gdr2[i][j];
        }
    }

    if (mask & FOURTH) {
   double d3gdr3[NR][NR][NR];
   int i, j, k;

        d3gdr3[0][0][0] = D3GDR0R0R0;
        d3gdr3[0][0][1] = D3GDR0R0R1;
        d3gdr3[0][1][0] = d3gdr3[0][0][1];
        d3gdr3[0][1][1] = D3GDR1R1R0;
        d3gdr3[1][0][0] = d3gdr3[0][0][1];
        d3gdr3[1][0][1] = d3gdr3[0][1][1];
        d3gdr3[1][1][0] = d3gdr3[0][1][1];
        d3gdr3[1][1][1] = D3GDR1R1R1;

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                for (k=0; k<NR; k++) dx3[i][j][k] = d3gdr3[i][j][k];
            }
        }
    }

}

void
hmixPlg(int mask, double t, double p, double *x,
    double *hmix /* Enthalpy of mixing BINARY MASK: 1 */
    )
{
    double xab = (x[0]          >  DBL_EPSILON) ? x[0]          : DBL_EPSILON;
    double xan = (x[1]          >  DBL_EPSILON) ? x[1]          : DBL_EPSILON;
    double xor = (1.0-x[0]-x[1] >  DBL_EPSILON) ? 1.0-x[0]-x[1] : DBL_EPSILON;

    *hmix = (G) + t*(S);
}

void
smixPlg(int mask, double t, double p, double *x,
    double *smix, /* Entropy of mixing                  BINARY MASK: 001 */
    double *dx,   /* (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010 */
    double **dx2  /* (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100 */
    )
{
    double xab = (x[0]          >  DBL_EPSILON) ? x[0]          : DBL_EPSILON;
    double xan = (x[1]          >  DBL_EPSILON) ? x[1]          : DBL_EPSILON;
    double xor = (1.0-x[0]-x[1] >  DBL_EPSILON) ? 1.0-x[0]-x[1] : DBL_EPSILON;

    if (mask & FIRST) {
        *smix = S;
    }

    if(mask & SECOND) {
        double d2gdrdt[NR];
        int i;

        d2gdrdt[0] = D2GDR0DT;
        d2gdrdt[1] = D2GDR1DT;

        for (i=0; i<NR; i++) dx[i] = - d2gdrdt[i];
    }

    if(mask & THIRD) {
        double d3gdr2dt[NR][NR];
        int i, j;

        d3gdr2dt[0][0] = D3GDR0R0DT;
        d3gdr2dt[0][1] = D3GDR0R1DT;
        d3gdr2dt[1][0] = d3gdr2dt[0][1];
        d3gdr2dt[1][1] = D3GDR1R1DT;

        for (i=0; i<NR; i++) {
       for (j=0; j<NR; j++) dx2[i][j] = - d3gdr2dt[i][j];
        }
    }

}

void
cpmixPlg(int mask, double t, double p, double *x,
    double *cpmix, /* Heat capacity of mixing               BINARY MASK: 001 */
    double *dt,    /* d(cp)/d(t)                            BINARY MASK: 010 */
    double *dx     /* d(cp)/d(x[])d(t)                      BINARY MASK: 100 */
    )
{
    double d2gdt2 = D2GDT2;

    if (mask & FIRST) {
        *cpmix = - t*d2gdt2;
    }

    if(mask & SECOND) {
        double d3gdt3   = D3GDT3;

        *dt = -t*d3gdt3 - d2gdt2;
    }

    if(mask & THIRD) {
        double d3gdrdt2[NR];
        int i;

        d3gdrdt2[0] = D3GDR0DT2;
        d3gdrdt2[1] = D3GDR1DT2;

        for (i=0; i<NR; i++) dx[i] = -t*d3gdrdt2[i];
    }

}

void
vmixPlg(int mask, double t, double p, double *x,
    double *vmix, /* Volume of mixing                BINARY MASK: 0000000001 */
    double *dx,   /* (pointer to dx[]) d(v)/d(x[])   BINARY MASK: 0000000010 */
    double **dx2, /* (point to dx2[][]) d(v)/d(x[])2 BINARY MASK: 0000000100 */
    double *dt,   /* d(v)/d(t)                       BINARY MASK: 0000001000 */
    double *dp,   /* d(v)/d(p)                       BINARY MASK: 0000010000 */
    double *dt2,  /* d2(v)/d(t)2                     BINARY MASK: 0000100000 */
    double *dtdp, /* d2(v)/d(t)d(p)                  BINARY MASK: 0001000000 */
    double *dp2,  /* d2(v)/d(p)2                     BINARY MASK: 0010000000 */
    double *dxdt, /* d2(v)/d(x[])d(t)                BINARY MASK: 0100000000 */
    double *dxdp  /* d2(v)/d(x[])d(p)                BINARY MASK: 1000000000 */
    )
{
    double xab = (x[0]          >  DBL_EPSILON) ? x[0]          : DBL_EPSILON;
    double xan = (x[1]          >  DBL_EPSILON) ? x[1]          : DBL_EPSILON;
    double xor = (1.0-x[0]-x[1] >  DBL_EPSILON) ? 1.0-x[0]-x[1] : DBL_EPSILON;

    if (mask & FIRST) {
        *vmix = DGDP;
    }

    if(mask & SECOND) {
        double d2gdrdp[NR];
        int i;

        d2gdrdp[0] = D2GDR0DP;
        d2gdrdp[1] = D2GDR1DP;

        for (i=0; i<NR; i++) dx[i] = d2gdrdp[i];
    }

    if(mask & THIRD) {
        double d3gdr2dp[NR][NR];
        int i, j;

        d3gdr2dp[0][0] = D3GDR0R0DP;
        d3gdr2dp[0][1] = D3GDR0R1DP;
        d3gdr2dp[1][0] = d3gdr2dp[0][1];
        d3gdr2dp[1][1] = D3GDR1R1DP;

        for (i=0; i<NR; i++) {
       for (j=0; j<NR; j++) dx2[i][j] = d3gdr2dp[i][j];
        }
    }

    if(mask & FOURTH) {
        *dt = D2GDTDP;
    }

    if(mask & FIFTH) {
        *dp = D2GDP2;
    }

    if(mask & SIXTH) {
        *dt2 = D3GDT2DP;
    }

    if(mask & SEVENTH) {
        *dtdp = D3GDTDP2;
    }

    if(mask & EIGHTH) {
        *dp2 = D3GDP3;
    }

    if(mask & NINTH) {
        double d3gdrdtdp[NR];
        int i;

        d3gdrdtdp[0] = D3GDR0DTDP;
        d3gdrdtdp[1] = D3GDR1DTDP;

        for (i=0; i<NR; i++) dxdt[i] = d3gdrdtdp[i];
    }

    if(mask & TENTH) {
        double d3gdrdp2[NR];
        int i;

        d3gdrdp2[0] = D3GDR0DP2;
        d3gdrdp2[1] = D3GDR1DP2;

        for (i=0; i<NR; i++) dxdp[i] = d3gdrdp2[i];
    }

/* original end of file plagioclase.C */


#if NEVER_DEFINED
static double DNa = 0.39;
/*
    {"K", "plagioclase", 0.18},
    {"Sc", "plagioclase", 0.02},
    {"Ti", "plagioclase", 0.04},
    {"Cr", "plagioclase", 0.05},
    {"Mn", "plagioclase", 0.05},
    {"Co", "plagioclase", 0.05},
    {"Rb", "plagioclase", 0.03},
    {"Sr", "plagioclase", 2.},
    {"Y", "plagioclase", 0.03},
    {"Zr", "plagioclase", 0.01},
    {"Nb", "plagioclase", 0.01},
    {"Cs", "plagioclase", 0.025},
    {"Ba", "plagioclase", 0.33},
    {"La", "plagioclase", 0.27},
    {"Ce", "plagioclase", 0.2},
    {"Pr", "plagioclase", 0.17},
    {"Nd", "plagioclase", 0.14},
    {"Sm", "plagioclase", 0.11},
    {"Eu", "plagioclase", 0.73},
    {"Gd", "plagioclase", 0.066},
    {"Tb", "plagioclase", 0.06},
    {"Dy", "plagioclase", 0.055},
    {"Ho", "plagioclase", 0.048},
    {"Er", "plagioclase", 0.041},
    {"Tm", "plagioclase", 0.036},
    {"Yb", "plagioclase", 0.031},
    {"Lu", "plagioclase", 0.025},
    {"Hf", "plagioclase", 0.01},
    {"Pb", "plagioclase", 0.36},
    {"Ra", "plagioclase", 0.01},
    {"Th", "plagioclase", 0.05},
    {"U", "plagioclase", 0.11},
    {"H2O", "plagioclase", 0.0005},
    {"plagioclase",
*/

void
conPlg(int inpMask, int outMask, double t, double p,
    double *e,      /* comp of plagioclase in moles of elements                    */
    double *m,      /* comp of plagioclase in moles of endmember components        */
    double *r,      /* comp of plagioclase in terms of the independent comp var    */
    double *x,      /* comp of plagioclase in mole fractions of endmember comp     */
    double **dm,    /* Jacobian matrix: dm[i][j] = dr[i]/dm[j]                  */
    double ***d2m,  /* vector of matrices: d2m[i][j][k] = d2r[i]/dm[j]dm[k]     */
    double **dr,    /* Jacobian matrix: dr[i][j] = dx[i]/dr[j]                  */



/* mask is an input mask */
void
getdPlg(int mask, double t, double p,
    double *xSol, /* solid composition  */
    double *xLiq  /* liquid composition */
    double *d0pt, /* partition coefficients */
    )
{
    int i;

    double xab = (x[0]          >  DBL_EPSILON) ? x[0]          : DBL_EPSILON;
    double xan = (x[1]          >  DBL_EPSILON) ? x[1]          : DBL_EPSILON;
    double xor = (1.0-x[0]-x[1] >  DBL_EPSILON) ? 1.0-x[0]-x[1] : DBL_EPSILON;

    //for (i=1; i<elsN - 1; i++) d0pt = 0.0;

    if (mask & FIRST) { /* default or user defined constant */


        /* if incoming is non-zero then use that - otherwise take defaut */
        fillD0PT;
    }
    else if (mask & SECOND) {

            for (i=0, oxSum=0.0; i<nc; i++) {
    	for (j=0, oxVal[i]=0.0; j<nlc; j++) oxVal[i] += (liquid[j].liqToOx)[i]*xliq[j];
    	oxVal[i] *= bulkSystem[i].mw;
    	oxSum	 += oxVal[i];
            }
            if (oxSum != 0.0) for (i=0; i<nc; i++) oxVal[i] /= oxSum;

            /* find Na2O value  of oxVal*/


        totNa  = r[0];
        na20 = totaNa/2.0* mw



        liqNa2O =

        D1_0 =



    }
    else if (mask & THIRD) {





    }


}

#endif

}

/* end of file plagioclase.C */
