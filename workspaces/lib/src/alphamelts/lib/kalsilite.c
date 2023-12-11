const char *kalsilite_ver(void) { return "$Id: kalsilite.c,v 1.3 2006/10/20 00:59:22 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: kalsilite.c,v $
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
 * Revision 1.3  1997/06/21  22:49:50  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 1.2  1997/05/03  20:23:29  ghiorso
 * *** empty log message ***
 *
 * Revision 1.1  1997/03/27  17:05:27  ghiorso
 * MELTS 3.0 (March 1997)
 *
 * Revision 1.3  1996/09/24  20:33:28  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 1.2  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 1.1  1995/11/01  22:40:27  ghiorso
 * Initial revision
 *
 * Revision 1.1  1995/11/01  22:40:27  ghiorso
 * Initial revision
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Routines to compute nepheline solution properties 
**      (file: KALSILITE.C) 
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso October 14, 1995 Original Version - started
**              Modified from SPINEL.C
**--
*/

#include "silmin.h"  /* Structure definitions foor SILMIN package */

#ifdef DEBUG
#undef DEBUG
#endif

#define SQUARE(x)  ((x)*(x))
#define CUBE(x)    ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

/*
 *=============================================================================
 * Nepheline solution parameters:
 * Sack, R.O., Ghiorso, M.S. (1995)
 * Thermodynamics of Nepheline Solutions
 */

                    /* kals - neph */
#define DH1                  32145.67 /* Na4Al4Si4O16 joules */
#define DS1                     10.46 /* joules/T */
#define DV1                       0.0 /* joules/bar */
#define DH2                   -5648.4 /* K4Al4Si4O16 joules */
#define DS2                       0.0 /* joules/T */
#define DV2                     -0.04 /* joules/bar */
#define DH3                   14644.0 /* []Na3Al3Si5O16 joules */
#define DS3                     -18.7 /* joules/T */
#define DV3                       0.0 /* joules/bar */
#define DH4                   35184.0 /* []CaNa2Al2Si2O16 joules */
#define DS4                     -7.17 /* joules/T */
#define DV4                       0.0 /* joules/bar */

                    /* kalsilite structure */
#define WKALS                 29288.0 /* WNa4-K4 joules */
#define WKVKALS               30334.0 /* WK4-[]Na3 joules */
#define WNVKALS               14646.0 /* WNa4-[]Na3 joules */
#define WNCKALS                   0.0 /* WNa4-[]CaNa2 joules */
#define WKCKALS               14644.0 /* WK4-[]CaNa2 joules */
#define WVVKALS                   0.0 /* W[]Na3-[]CaNa2 joules */

/*
 * Global (to this file): variables 
 */

#define R  8.3143
#define NR         3    /* Three independent composition variables */
#define NA         4    /* Four endmember compositions             */
                   /* Penalty constant for approaching the zero
                      concentration of r[1] - see #define G        */
#define PENALTY    sqrt(DBL_EPSILON)

                /* correction to zero point entropy of []CaNa2Al4Si4O16 */
#define S4      -15.8765
#define G4      -t*(S4)

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conSpn defines the conversion from m[i], to r[j]
 */
                   /* Order: X2, X3, X4 */
#define FR2(i)     (i == 1) ? 1.0 - r[0] : - r[0]
#define FR3(i)     (i == 2) ? 1.0 - r[1] : - r[1]
#define FR4(i)     (i == 3) ? 1.0 - r[2] : - r[2]

#define DFR2DR2(i) - 1.0                                 
#define DFR3DR3(i) - 1.0
#define DFR4DR4(i) - 1.0

/*
 * Global (to this file): derivative definitions
 */

#define S (DS1)*(1.0-r[0]-r[1]-r[2]) + (DS2)*r[0] + (DS3)*r[1] + ((DS4)+(S4))*r[2] \
                   - 4.0*R*(xk*log(xk) + (xvc-xca)*log(xvc-xca) + xna*log(xna) + xca*log(xca))
#define H (DH1)*(1.0-r[0]-r[1]-r[2]) + (DH2)*r[0] + (DH3)*r[1] + (DH4)*r[2] \
                   + (WKALS)*r[0]*(1.0-r[0]-r[1]-r[2]) + (WNVKALS)*r[1]*(1.0-r[0]-r[1]-r[2]) \
                   + (WNCKALS)*r[2]*(1.0-r[0]-r[1]-r[2]) + (WKVKALS)*r[0]*r[1] \
                   + (WKCKALS)*r[0]*r[2] + (WVVKALS)*r[1]*r[2]
#define V (DV1)*(1.0-r[0]-r[1]-r[2]) + (DV2)*r[0] + (DV3)*r[1] + (DV4)*r[2]
#define G (H) - t*(S) + (p-1.0)*(V) + ((r[1] > 0.0) ? (PENALTY)/r[1] : DBL_MAX)

/*----------------------------------------------------------------------------*/

#define DGDR0 - ((DH1)-t*(DS1)+(p-1.0)*(DV1)) \
              + ((DH2)-t*(DS2)+(p-1.0)*(DV2)) + R*t*4.0*(log(xk) - log(xna)) \
              + (WKALS)*(1.0-2.0*r[0]-r[1]-r[2]) - (WNVKALS)*r[1] \
              - (WNCKALS)*r[2] + (WKVKALS)*r[1] + (WKCKALS)*r[2]
#define DGDR1 - ((DH1)-t*(DS1)+(p-1.0)*(DV1)) \
              + ((DH3)-t*(DS3)+(p-1.0)*(DV3)) + R*t*(log(xvc-xca) - log(xna)) \
              - (WKALS)*r[0] + (WNVKALS)*(1.0-r[0]-2.0*r[1]-r[2]) \
              - (WNCKALS)*r[2] + (WKVKALS)*r[0] + (WVVKALS)*r[2] \
              - ((r[1] > 0.0) ? (PENALTY)/(r[1]*r[1]) : DBL_MAX)
#define DGDR2 - ((DH1)-t*(DS1)+(p-1.0)*(DV1)) \
              + ((DH4)-t*(DS4)+(p-1.0)*(DV4)) + (G4) \
              + R*t*(log(xca) - 2.0*log(xna) - 1.0) \
              - (WKALS)*r[0] - (WNVKALS)*r[1] \
              + (WNCKALS)*(1.0-r[0]-r[1]-2.0*r[2]) \
              + (WKCKALS)*r[0] + (WVVKALS)*r[1]
#define DGDT  - (S)
#define DGDP  (V)

/*----------------------------------------------------------------------------*/

#define D2GDR0R0 R*t*4.0*(1.0/xk + 1.0/xna) - 2.0*(WKALS)
#define D2GDR0R1 R*t/xna - (WKALS) + (WKVKALS) - (WNVKALS)
#define D2GDR0R2 R*t*2.0/xna - (WKALS) - (WNCKALS) + (WKCKALS)
#define D2GDR0DT (DS1) - (DS2) + R*4.0*(log(xk) - log(xna))
#define D2GDR0DP - (DV1) + (DV2)

#define D2GDR1R1 R*t*0.25*(1.0/(xvc-xca) + 1.0/xna) - 2.0*(WNVKALS) \
                 + ((r[1] > 0.0) ? 2.0*(PENALTY)/(r[1]*r[1]*r[1]) : DBL_MAX)
#define D2GDR1R2 R*t*(0.5/xna) - (WNVKALS) - (WNCKALS) + (WVVKALS)
#define D2GDR1DT (DS1) - (DS3) + R*(log(xvc-xca) - log(xna))
#define D2GDR1DP - (DV1) + (DV3)

#define D2GDR2R2 R*t*(0.25/xca + 1.0/xna) - 2.0*(WNCKALS)
#define D2GDR2DT (DS1) - (DS4) - (S4) + R*(log(xca) - 2.0*log(xna) - 1.0)
#define D2GDR2DP - (DV1) + (DV4)

#define D2GDT2   0.0
#define D2GDTDP  0.0
#define D2GDP2   0.0

/*----------------------------------------------------------------------------*/

#define D3GDR0R0R0 - R*t*4.0*(1.0/(xk*xk) - 1.0/(xna*xna))
#define D3GDR0R0R1 R*t/(xna*xna)

#define D3GDR0R0R2 R*t*2.0/(xna*xna)
#define D3GDR0R0DT R*4.0*(1.0/xk + 1.0/xna)
#define D3GDR0R0DP 0.0

#define D3GDR0R1R1 R*t*0.25/(xna*xna)
#define D3GDR0R1R2 R*t*0.5/(xna*xna)
#define D3GDR0R1DT R/xna
#define D3GDR0R1DP 0.0

#define D3GDR0R2R2 R*t/(xna*xna)
#define D3GDR0R2DT R*2.0/xna
#define D3GDR0R2DP 0.0
  
#define D3GDR1R1R1 - R*t*0.25*0.25*(1.0/((xvc-xca)*(xvc-xca)) - 1.0/(xna*xna)) \
                   - ((r[1] > 0.0) ? 6.0*(PENALTY)/(r[1]*r[1]*r[1]*r[1]) : DBL_MAX)
#define D3GDR1R1R2 R*t*0.25*(0.5/(xna*xna))
#define D3GDR1R1DT R*0.25*(1.0/(xvc-xca) + 1.0/xna)
#define D3GDR1R1DP 0.0

#define D3GDR1R2R2 R*t*(0.5*0.5/(xna*xna))
#define D3GDR1R2DT R*(0.5/xna)
#define D3GDR1R2DP 0.0

#define D3GDR2R2R2 R*t*(- 0.25*0.25/(xca*xca) + 0.5/(xna*xna))
#define D3GDR2R2DT R*(0.25/xca + 1.0/xna)
#define D3GDR2R2DP 0.0

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
#define D3GDR2DT2  0.0
#define D3GDR2DTDP 0.0
#define D3GDR2DP2  0.0

/*
 *=============================================================================
 * Macros for automatic array initialization
 */
#define fillD2GDR2 \
 d2gdr2[0][0] = D2GDR0R0;     d2gdr2[0][1] = D2GDR0R1;     d2gdr2[0][2] = D2GDR0R2; \
 d2gdr2[1][0] = d2gdr2[0][1]; d2gdr2[1][1] = D2GDR1R1;     d2gdr2[1][2] = D2GDR1R2; \
 d2gdr2[2][0] = d2gdr2[0][2]; d2gdr2[2][1] = d2gdr2[1][2]; d2gdr2[2][2] = D2GDR2R2;

#define fillD2GDRDT d2gdrdt[0] = D2GDR0DT;    d2gdrdt[1] = D2GDR1DT;    d2gdrdt[2] = D2GDR2DT;
#define fillD2GDRDP d2gdrdp[0] = D2GDR0DP;    d2gdrdp[1] = D2GDR1DP;    d2gdrdp[2] = D2GDR2DP;

#define fillD3GDR3 \
 d3gdr3[0][0][0] = D3GDR0R0R0;		d3gdr3[0][0][1] = D3GDR0R0R1;      d3gdr3[0][0][2] = D3GDR0R0R2; \
 d3gdr3[0][1][0] = d3gdr3[0][0][1];	d3gdr3[0][1][1] = D3GDR0R1R1;      d3gdr3[0][1][2] = D3GDR0R1R2; \
 d3gdr3[0][2][0] = d3gdr3[0][0][2];     d3gdr3[0][2][1] = d3gdr3[0][1][2]; d3gdr3[0][2][2] = D3GDR0R2R2; \
 d3gdr3[1][0][0] = d3gdr3[0][0][1];	d3gdr3[1][0][1] = d3gdr3[0][1][1]; d3gdr3[1][0][2] = d3gdr3[0][1][2]; \
 d3gdr3[1][1][0] = d3gdr3[0][1][1];	d3gdr3[1][1][1] = D3GDR1R1R1;      d3gdr3[1][1][2] = D3GDR1R1R2; \
 d3gdr3[1][2][0] = d3gdr3[0][1][2];	d3gdr3[1][2][1] = d3gdr3[1][1][2]; d3gdr3[1][2][2] = D3GDR1R2R2; \
 d3gdr3[2][0][0] = d3gdr3[0][0][2];	d3gdr3[2][0][1] = d3gdr3[0][1][2]; d3gdr3[2][0][2] = d3gdr3[0][2][2]; \
 d3gdr3[2][1][0] = d3gdr3[0][1][2];	d3gdr3[2][1][1] = d3gdr3[1][1][2]; d3gdr3[2][1][2] = d3gdr3[1][2][2]; \
 d3gdr3[2][2][0] = d3gdr3[0][2][2];	d3gdr3[2][2][1] = d3gdr3[1][2][2]; d3gdr3[2][2][2] = D3GDR2R2R2; \

#define fillD3GDR2DT \
 d3gdr2dt[0][0] = D3GDR0R0DT;     d3gdr2dt[0][1] = D3GDR0R1DT;     d3gdr2dt[0][2] = D3GDR0R2DT; \
 d3gdr2dt[1][0] = d3gdr2dt[0][1]; d3gdr2dt[1][1] = D3GDR1R1DT;     d3gdr2dt[1][2] = D3GDR1R2DT; \
 d3gdr2dt[2][0] = d3gdr2dt[0][2]; d3gdr2dt[2][1] = d3gdr2dt[1][2]; d3gdr2dt[2][2] = D3GDR2R2DT;

#define fillD3GDR2DP \
 d3gdr2dp[0][0] = D3GDR0R0DP;     d3gdr2dp[0][1] = D3GDR0R1DP;     d3gdr2dp[0][2] = D3GDR0R2DP; \
 d3gdr2dp[1][0] = d3gdr2dp[0][1]; d3gdr2dp[1][1] = D3GDR1R1DP;     d3gdr2dp[1][2] = D3GDR1R2DP; \
 d3gdr2dp[2][0] = d3gdr2dp[0][2]; d3gdr2dp[2][1] = d3gdr2dp[1][2]; d3gdr2dp[2][2] = D3GDR2R2DP;

#define fillD3GDRDT2  d3gdrdt2[0] = D3GDR0DT2;   d3gdrdt2[1] = D3GDR1DT2;   d3gdrdt2[2] = D3GDR2DT2;
#define fillD3GDRDTDP d3gdrdtdp[0] = D3GDR0DTDP; d3gdrdtdp[1] = D3GDR1DTDP; d3gdrdtdp[2] = D3GDR2DTDP;
#define fillD3GDRDP2  d3gdrdp2[0] = D3GDR0DP2;   d3gdrdp2[1] = D3GDR1DP2;   d3gdrdp2[2] = D3GDR2DP2;

/*
 *=============================================================================
 * Public functions:
 *    mask  -  bitwise mask for selecting output
 *    t     -  Temperature (K)
 *    p     -  Pressure (bars)
 *    *x    -  (pointer to x[]) Array of independent compositional variables
 */
int
testKal(int mask, double t, double p,
  int na,          /* Expected number of endmember components                 */
  int nr,          /* Expected number of independent compositional variables  */
  char **names,    /* array of strings of names of endmember components       */
  char **formulas, /* array of strings of formulas of endmember components    */
  double *r,       /* array of indepependent compos variables, check bounds   */
  double *m)       /* array of moles of endmember components, check bounds    */
{
  const char *phase = "kalsilite.c";
  const char *NAMES[NA]    = { "na-nepheline", "k-nepheline", "vc-nepheline", "ca-nepheline" }; 
  const char *FORMULAS[NA] = { "Na4Al4Si4O16", "K4Al4Si4O16", "Na3Al3Si5O16", "CaNa2Al4Si4O16" }; 
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

    result = result && (r[0] >= 0.0) && (r[0] <= 1.0-r[1]/4.0-r[2]/2.0);
    result = result && (r[1] >= 0.0) && (r[1] <= 1.0);
    result = result && (r[2] >= 0.0) && (r[2] <= 1.0);

    result = result && (4.0*r[0] >= 0.0) && (4.0*r[0] <= 4.0+sqrt(DBL_EPSILON));               /* tot K  */
    result = result && (4.0-4.0*r[0]-r[1]-2.0*r[2] >= 0.0)                                     /* tot Na */
                    && (4.0-4.0*r[0]-r[1]-2.0*r[2] <= 4.0+sqrt(DBL_EPSILON));
    result = result && (r[2] >= 0.0) && (r[2] <= 1.0+sqrt(DBL_EPSILON));                       /* tot Ca */
    result = result && (r[1]+r[2] >= 0.0) && (r[1]+r[2] <= 1.0+sqrt(DBL_EPSILON));             /* tot Vc */
    result = result && (4.0-r[1] >= 3.0-sqrt(DBL_EPSILON)) && (4.0-r[1] <= 5.0+sqrt(DBL_EPSILON));   /* tot Al */
    result = result && (4.0+r[1] >= 3.0-sqrt(DBL_EPSILON)) && (4.0+r[1] <= 5.0+sqrt(DBL_EPSILON));   /* tot Si */ 
  }
  /* Check bounds on moles of endmember components */
  if (mask & SIXTH) {
    int pResult;
    for (i=0, sum=0.0; i<NA; i++) sum += m[i];
    result = result && (sum >= 0.0);
    if (sum > 0.0) {
      pResult = (m[0] >= -3.0*m[2]/4.0-m[3]/2.0-DBL_EPSILON) && (m[0] <= sum+DBL_EPSILON);
#ifdef DEBUG
      if (!pResult) printf("Bound check m[0] : %g <= %g <= %g\n", -3.0*m[2]/4.0-m[3]/2.0, m[0], sum);
#endif
      result = result && pResult;

      pResult = (m[1] >= 0.0) && (m[1] <= sum-m[2]/4.0-m[3]/2.0+DBL_EPSILON);
#ifdef DEBUG
      if (!pResult) printf("Bound check m[1] : %g <= %g <= %g\n", 0.0, m[1], sum-m[2]/4.0-m[3]/2.0);
#endif
      result = result && pResult;

      pResult = (m[2] >= 0.0) && (m[2] <= sum+DBL_EPSILON);
#ifdef DEBUG
      if (!pResult) printf("Bound check m[2] : %g <= %g <= %g\n", 0.0, m[2], sum);
#endif
      result = result && pResult;

      pResult = (m[3] >= 0.0) && (m[3] <= sum+DBL_EPSILON);
#ifdef DEBUG
      if (!pResult) printf("Bound check m[3] : %g <= %g <= %g\n", 0.0, m[3], sum);
#endif
      result = result && pResult;

      pResult = (4.0*m[1] >= 0.0) && (4.0*m[1] <= 4.0*sum+sqrt(DBL_EPSILON));       /* tot K  */
#ifdef DEBUG
      if (!pResult) printf("Bound check tot K : %g <= %g <= %g\n", 0.0, 4.0*m[1], 4.0*sum);
#endif
      result = result && pResult;

      pResult = (4.0*m[0]+3.0*m[2]+2.0*m[3] >= 0.0)                                 /* tot Na */
             && (4.0*m[0]+3.0*m[2]+2.0*m[3] <= 4.0*sum+sqrt(DBL_EPSILON));
#ifdef DEBUG
      if (!pResult) printf("Bound check tot Na: %g <= %g <= %g\n", 0.0, 4.0*m[0]+3.0*m[2]+2.0*m[3], 4.0*sum);
#endif
      result = result && pResult;

      pResult = (m[3] >= 0.0) && (m[3] <= sum+sqrt(DBL_EPSILON));                   /* tot Ca */
#ifdef DEBUG
      if (!pResult) printf("Bound check tot Ca: %g <= %g <= %g\n", 0.0, m[3], sum);
#endif
      result = result && pResult;

      pResult = (m[2]+m[3] >= 0.0) && (m[2]+m[3] <= sum+sqrt(DBL_EPSILON));         /* tot Vc */
#ifdef DEBUG
      if (!pResult) printf("Bound check tot Vc: %g <= %g <= %g\n", 0.0, m[2]+m[3], sum);
#endif
      result = result && pResult;

      pResult = (4.0*m[0]+4.0*m[1]+3.0*m[2]+4.0*m[3] >= 3.0*sum-sqrt(DBL_EPSILON))  /* tot Al */
             && (4.0*m[0]+4.0*m[1]+3.0*m[2]+4.0*m[3] <= 5.0*sum+sqrt(DBL_EPSILON));
#ifdef DEBUG
      if (!pResult) printf("Bound check tot Al: %g <= %g <= %g\n", 3.0*sum, 
        4.0*m[0]+4.0*m[1]+3.0*m[2]+4.0*m[3], 5.0*sum);
#endif
      result = result && pResult;

      pResult = (4.0*m[0]+4.0*m[1]+5.0*m[2]+4.0*m[3] >= 3.0*sum-sqrt(DBL_EPSILON))  /* tot Si */
             && (4.0*m[0]+4.0*m[1]+5.0*m[2]+4.0*m[3] <= 5.0*sum+sqrt(DBL_EPSILON));
#ifdef DEBUG
      if (!pResult) printf("Bound check tot Si: %g <= %g <= %g\n", 3.0*sum, 
        4.0*m[0]+4.0*m[1]+5.0*m[2]+4.0*m[3], 5.0*sum);
#endif
      result = result && pResult;
    }
  }

  return result;
}

void
conKal(int inpMask, int outMask, double t, double p,
  double *e,      /* comp of spinel in moles of elements                      */
  double *m,      /* comp of spinel in moles of endmember components          */
  double *r,      /* comp of spinel in terms of the independent comp var      */
  double *x,      /* comp of spinel in mole fractions of endmember comp       */
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
      endmember spinel components.
  (2) calculates from a vector of moles of endmember components, one or
      all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
  (3) calculates from a vector of independent compositional variables
      mole fractions of endmember components and/or the Jacobian matrix
      dx[]/dr[]

  In this routine it is assumed that the elements are in the order of atomic
  numbers and that the order of kalsilite components has been verified as:
      m[0] = na-nepheline   (Na4Al4Si4O16) ,
      m[1] = k-nepheline    (K4Al4Si4O16), 
      m[2] = vc-nepheline   (Na3Al3Si5O16),
      m[3] = ca-nepheline   (CaNa2Al4Si4O16),

  ----------------------------------------------------------------------------*/

  int i, j, k;

  if (inpMask == FIRST && outMask == SECOND) {
    /* Converts a vector of moles of elements into a vector of moles of
       end-member components.                                                 */
    static const int Na = 11;
    static const int Al = 13;
    static const int Si = 14;
    static const int K  = 19;
    static const int Ca = 20;
    double divalent = e[Ca];
    double vacancy  = (e[Al]+e[Si])/2.0 - e[Na] - e[K] - divalent;

    m[0] = (e[Na] - 3.0*vacancy + divalent)/4.0;  /* Moles of Na4Al4Si4O16 */
    m[1] = e[K]/4.0;                              /* Moles of K4Al4Si4O16  */
    m[2] = vacancy-divalent;                      /* Moles of Na3Al3Si5O16 */
    m[3] = divalent;                              /* Moles of CaNa2Al4Si4O16 */

  } else if (inpMask == SECOND) {
    double sum;

    if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
      printf("Illegal call to conNph with inpMask = %o and outMask = %o\n",
        inpMask, outMask);

    for (i=0, sum=0.0; i<NA; i++) sum += m[i];

    if (outMask & THIRD) {
      /* Converts a vector of moles of end-member components (m) into a vector 
         of independent compositional variables (r) required as input for the 
         remaining public functions.                                          */
      r[0] = (sum != 0.0) ? m[1]/sum : 0.0;  /* X2 = X K4Al4Si4O16    */ 
      r[1] = (sum != 0.0) ? m[2]/sum : 0.0;  /* X3 = X Na3Al3Si5O16   */
      r[2] = (sum != 0.0) ? m[3]/sum : 0.0;  /* X4 = X CaNa2Al4Si4O16 */
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
        for (j=0; j<NA; j++) {
          dm[0][j] = (j == 1) ? (1.0 - m[1]/sum)/sum : -m[1]/SQUARE(sum);
          dm[1][j] = (j == 2) ? (1.0 - m[2]/sum)/sum : -m[2]/SQUARE(sum);
          dm[2][j] = (j == 3) ? (1.0 - m[3]/sum)/sum : -m[3]/SQUARE(sum);
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
        for (j=0; j<NA; j++) {
          for (k=0; k<NA; k++) {
            d2m[0][j][k]  = 2.0*m[1]/CUBE(sum);
            d2m[0][j][k] -= (j == 1) ? 1.0/SQUARE(sum) : 0.0;
            d2m[0][j][k] -= (k == 1) ? 1.0/SQUARE(sum) : 0.0;
            d2m[1][j][k]  = 2.0*m[2]/CUBE(sum);
            d2m[1][j][k] -= (j == 2) ? 1.0/SQUARE(sum) : 0.0;
            d2m[1][j][k] -= (k == 2) ? 1.0/SQUARE(sum) : 0.0;
            d2m[2][j][k]  = 2.0*m[3]/CUBE(sum);
            d2m[2][j][k] -= (j == 3) ? 1.0/SQUARE(sum) : 0.0;
            d2m[2][j][k] -= (k == 3) ? 1.0/SQUARE(sum) : 0.0;
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
        for (j=0; j<NA; j++)  {
          for (k=0; k<NA; k++)  {
            for (l=0; l<NA; l++)  {
              d3m[0][j][k][l]  = -6.0*m[1]/QUARTIC(sum);
              d3m[0][j][k][l] += (j == 1) ? 2.0/CUBE(sum) : 0.0;
              d3m[0][j][k][l] += (k == 1) ? 2.0/CUBE(sum) : 0.0;
              d3m[0][j][k][l] += (l == 1) ? 2.0/CUBE(sum) : 0.0;
              d3m[1][j][k][l]  = -6.0*m[2]/QUARTIC(sum);
              d3m[1][j][k][l] += (j == 2) ? 2.0/CUBE(sum) : 0.0;
              d3m[1][j][k][l] += (k == 2) ? 2.0/CUBE(sum) : 0.0;
              d3m[1][j][k][l] += (l == 2) ? 2.0/CUBE(sum) : 0.0;
              d3m[2][j][k][l]  = -6.0*m[3]/QUARTIC(sum);
              d3m[2][j][k][l] += (j == 3) ? 2.0/CUBE(sum) : 0.0;
              d3m[2][j][k][l] += (k == 3) ? 2.0/CUBE(sum) : 0.0;
              d3m[2][j][k][l] += (l == 3) ? 2.0/CUBE(sum) : 0.0;
            }
          }
        }
      }
    }

  } else if (inpMask == THIRD) {

    if (outMask & ~(FOURTH | SEVENTH))
      printf("Illegal call to conNph with inpMask = %o and outMask = %o\n",
        inpMask, outMask);

    if (outMask & FOURTH) {
      /* Converts a vector of independent compositional variables (r) into a
         vector of mole fractions of endmember components (x).                */
      x[0] = 1.0 - r[0] - r[1] - r[2];
      x[1] = r[0];
      x[2] = r[1];
      x[3] = r[2];
    }

    if (outMask & SEVENTH) {
      /* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
      for (i=0; i<NA; i++) for (j=0; j<NR; j++) dr[i][j] = 0.0;
      dr[0][0] = -1.0; dr[0][1] = -1.0; dr[0][2] = -1.0;
      dr[1][0] =  1.0;
      dr[2][1] =  1.0;
      dr[3][2] =  1.0;
    }

  } else  {
    printf("Illegal call to conNph with inpMask = %o and outMask = %o\n",
      inpMask, outMask);
  }

}

void
dispKal(int mask, double t, double p, double *x,
  char **formula            /* Mineral formula for interface display MASK: 1 */
  )
{
  double *r = x;
  static char masterString[] = {
/*             1111111111222222222233333333334444444444555555555566666666667
     01234567890123456789012345678901234567890123456789012345678901234567890 */
    "kals Na_.__K_.__Ca_.__[]_.__Al_.__Si_.__O16" };

  if (mask & FIRST) {
    char *string = strcpy((char *) malloc((size_t) (strlen(masterString)+1)*sizeof(char)), masterString);
    double totNa, totK, totCa, totVc, totAl, totSi;
    char n[5];
    int i;

    totNa = 4.0*(1.0-r[0]-r[1]-r[2]) + 3.0*r[1] + 2.0*r[2];
    totK  = 4.0*r[0];
    totCa = r[2];
    totVc = r[1] + r[2];
    totAl = 4.0*(1.0-r[0]-r[1]-r[2]) + 4.0*r[0] + 3.0*r[1] + 4.0*r[2];
    totSi = 4.0*(1.0-r[0]-r[1]-r[2]) + 4.0*r[0] + 5.0*r[1] + 4.0*r[2];

    (void) snprintf(n, 5, "%4.2f", totNa); for (i=0; i<4; i++) string[ 7+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totK);  for (i=0; i<4; i++) string[12+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totCa); for (i=0; i<4; i++) string[18+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totVc); for (i=0; i<4; i++) string[24+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totAl); for (i=0; i<4; i++) string[30+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totSi); for (i=0; i<4; i++) string[36+i] = n[i];

    *formula = string;
  }
}

void 
actKal(int mask, double t, double p, double *x, 
  double *a,  /* (pointer to a[]) activities              BINARY MASK: 0001 */
  double *mu, /* (pointer to mu[]) chemical potentials    BINARY MASK: 0010 */
  double **dx /* (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100 */
  )           /* exclusion criteria applied to results if BINARY MASK: 1000 */
{
  double *r = x;
  double g, dgdr[NR];
  double fr[NA][NR];
  int i, j;
  double xk  = r[0];
  double xvc = (r[1]+r[2])/4.0;
  double xna = 1.0 - r[0] - r[1]/4.0 - r[2]/2.0;
  double xca = r[2]/4.0;

  if (xk  <= DBL_EPSILON) xk  = DBL_EPSILON;
  if (xvc <= DBL_EPSILON) xvc = DBL_EPSILON;
  if (xna <= DBL_EPSILON) xna = DBL_EPSILON;
  if (xca <= DBL_EPSILON) xca = DBL_EPSILON;
  if (xk  >= 1.0 - DBL_EPSILON) xk  = 1.0 - DBL_EPSILON;
  if (xvc >= 1.0 - DBL_EPSILON) xvc = 1.0 - DBL_EPSILON;
  if (xna >= 1.0 - DBL_EPSILON) xna = 1.0 - DBL_EPSILON;
  if (xca >= 1.0 - DBL_EPSILON) xca = 1.0 - DBL_EPSILON;

  for(i=0; i<NA; i++) {
     fr[i][0] = FR2(i); /* X2 */
     fr[i][1] = FR3(i); /* X3 */
     fr[i][2] = FR4(i); /* X4 */
  }

  g       = G;
  dgdr[0] = DGDR0;          
  dgdr[1] = DGDR1;
  dgdr[2] = DGDR2;

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
    for(i=0; i<NA; i++)
       for (mu[i]=g, j=0; j<NR; j++) mu[i] += fr[i][j]*dgdr[j];
  }

  if (mask & THIRD) {
    double d2gdr2[NR][NR], dfrdr[NA][NR], sum;
    int k;

    fillD2GDR2

    for(i=0; i<NA; i++) {
       dfrdr[i][0] = DFR2DR2(i); /* X2 */
       dfrdr[i][1] = DFR3DR3(i); /* X3 */
       dfrdr[i][2] = DFR4DR4(i); /* X4 */
    }

    for (i=0; i<NA; i++) {
      for (k=0; k<NR; k++) {
        /* compute activity of the i-th component */
        for (dx[i][k]=g, j=0; j<NR; j++) dx[i][k] += fr[i][j]*dgdr[j];
        dx[i][k] = exp(dx[i][k]/(R*t));

        /* compute derivative of i-th activity with respect to r(k) */
        sum = (1.0+dfrdr[i][k])*dgdr[k];
        for (j=0; j<NR; j++) sum += fr[i][j]*d2gdr2[j][k];
        dx[i][k] *= sum/(R*t);
      }
    }
  }

  if (mask & FOURTH) {
    /* implement exclusion criteria on quantities for preclb routines         */
    static const double exclusion[NA] = {
       0.05,  /* exclusion criteria on the mole fraction of Na+  */
       0.05,  /* exclusion criteria on the mole fraction of K+   */
       0.05,  /* exclusion criteria on the mole fraction of Ca++ */
       0.05,  /* exclusion criteria on the mole fraction of Vc   */ 
    };
    double x[NA], totNa, totK, totCa, totVc;

    totNa = 4.0*(1.0-r[0]-r[1]-r[2]) + 3.0*r[1] + 2.0*r[2];
    totK  = 4.0*r[0];
    totCa = r[2];
    totVc = r[1] + r[2];

    x[0] = totNa/4.0;   /* Na+ averaged on both sites */
    x[1] = totK/4.0;    /* K+ averaged on both sites  */
    x[2] = totVc-totCa; /* Vc on large site           */
    x[3] = totCa;       /* Ca on small site           */
    
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
gmixKal(int mask, double t, double p, double *x, 
  double *gmix, /* Gibbs energy of mixing             BINARY MASK: 0001 */
  double *dx,   /* (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010 */
  double **dx2, /* (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100 */
  double ***dx3 /* (pointer to dx3[][][]) d3(g)/d(x[])3 NARY MASK: 1000 */
  )
{
  double *r = x;
  double xk  = r[0];
  double xvc = (r[1]+r[2])/4.0;
  double xna = 1.0 - r[0] - r[1]/4.0 - r[2]/2.0;
  double xca = r[2]/4.0;

  if (xk  <= DBL_EPSILON) xk  = DBL_EPSILON;
  if (xvc <= DBL_EPSILON) xvc = DBL_EPSILON;
  if (xna <= DBL_EPSILON) xna = DBL_EPSILON;
  if (xca <= DBL_EPSILON) xca = DBL_EPSILON;
  if (xk  >= 1.0 - DBL_EPSILON) xk  = 1.0 - DBL_EPSILON;
  if (xvc >= 1.0 - DBL_EPSILON) xvc = 1.0 - DBL_EPSILON;
  if (xna >= 1.0 - DBL_EPSILON) xna = 1.0 - DBL_EPSILON;
  if (xca >= 1.0 - DBL_EPSILON) xca = 1.0 - DBL_EPSILON;

  if (mask & FIRST) {
    *gmix = G;
  }
  
  if(mask & SECOND) {
    dx[0] = DGDR0;
    dx[1] = DGDR1;
    dx[2] = DGDR2;
  }

  if(mask & THIRD) {
    double d2gdr2[NR][NR];
    int i, j;

    fillD2GDR2
    for (i=0; i<NR; i++) for (j=0; j<NR; j++) dx2[i][j] = d2gdr2[i][j];
  }

  if(mask & FOURTH) {
    double d3gdr3[NR][NR][NR];
    int i, j, k;

    fillD3GDR3
    for (i=0; i<NR; i++) for (j=0; j<NR; j++) for (k=0; k<NR; k++) dx3[i][j][k] = d3gdr3[i][j][k];
  }

}

void 
hmixKal(int mask, double t, double p, double *x, 
  double *hmix /* Enthalpy of mixing BINARY MASK: 1 */
  )
{
  double *r = x;
  double xk  = r[0];
  double xvc = (r[1]+r[2])/4.0;
  double xna = 1.0 - r[0] - r[1]/4.0 - r[2]/2.0;
  double xca = r[2]/4.0;

  if (xk  <= DBL_EPSILON) xk  = DBL_EPSILON;
  if (xvc <= DBL_EPSILON) xvc = DBL_EPSILON;
  if (xna <= DBL_EPSILON) xna = DBL_EPSILON;
  if (xca <= DBL_EPSILON) xca = DBL_EPSILON;
  if (xk  >= 1.0 - DBL_EPSILON) xk  = 1.0 - DBL_EPSILON;
  if (xvc >= 1.0 - DBL_EPSILON) xvc = 1.0 - DBL_EPSILON;
  if (xna >= 1.0 - DBL_EPSILON) xna = 1.0 - DBL_EPSILON;
  if (xca >= 1.0 - DBL_EPSILON) xca = 1.0 - DBL_EPSILON;

  *hmix = (G) + t*(S);
}

void 
smixKal(int mask, double t, double p, double *x, 
  double *smix, /* Entropy of mixing                  BINARY MASK: 001 */
  double *dx,   /* (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010 */
  double **dx2  /* (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100 */
  )
{
  double *r = x;
  double xk  = r[0];
  double xvc = (r[1]+r[2])/4.0;
  double xna = 1.0 - r[0] - r[1]/4.0 - r[2]/2.0;
  double xca = r[2]/4.0;

  if (xk  <= DBL_EPSILON) xk  = DBL_EPSILON;
  if (xvc <= DBL_EPSILON) xvc = DBL_EPSILON;
  if (xna <= DBL_EPSILON) xna = DBL_EPSILON;
  if (xca <= DBL_EPSILON) xca = DBL_EPSILON;
  if (xk  >= 1.0 - DBL_EPSILON) xk  = 1.0 - DBL_EPSILON;
  if (xvc >= 1.0 - DBL_EPSILON) xvc = 1.0 - DBL_EPSILON;
  if (xna >= 1.0 - DBL_EPSILON) xna = 1.0 - DBL_EPSILON;
  if (xca >= 1.0 - DBL_EPSILON) xca = 1.0 - DBL_EPSILON;

  if (mask & FIRST) {
    *smix = S; 
  }
  
  if(mask & SECOND) {
    double d2gdrdt[NR];
    int i;

    fillD2GDRDT    
    for (i=0; i<NR; i++) dx[i] = -d2gdrdt[i];
  }

  if(mask & THIRD) {
    double d3gdr2dt[NR][NR];
    int i, j;

    fillD3GDR2DT

    for (i=0; i<NR; i++) for (j=0; j<NR; j++) dx2[i][j] = -d3gdr2dt[i][j];
  }

}

void 
cpmixKal(int mask, double t, double p, double *x, 
  double *cpmix, /* Heat capacity of mixing               BINARY MASK: 001 */
  double *dt,    /* d(cp)/d(t)                            BINARY MASK: 010 */
  double *dx     /* d(cp)/d(x[])d(t)                      BINARY MASK: 100 */
  )
{
  double *r = x;
  double d2gdt2;
  int i;
  double xk  = r[0];
  double xvc = (r[1]+r[2])/4.0;
  double xna = 1.0 - r[0] - r[1]/4.0 - r[2]/2.0;
  double xca = r[2]/4.0;

  if (xk  <= DBL_EPSILON) xk  = DBL_EPSILON;
  if (xvc <= DBL_EPSILON) xvc = DBL_EPSILON;
  if (xna <= DBL_EPSILON) xna = DBL_EPSILON;
  if (xca <= DBL_EPSILON) xca = DBL_EPSILON;
  if (xk  >= 1.0 - DBL_EPSILON) xk  = 1.0 - DBL_EPSILON;
  if (xvc >= 1.0 - DBL_EPSILON) xvc = 1.0 - DBL_EPSILON;
  if (xna >= 1.0 - DBL_EPSILON) xna = 1.0 - DBL_EPSILON;
  if (xca >= 1.0 - DBL_EPSILON) xca = 1.0 - DBL_EPSILON;

  d2gdt2  = D2GDT2;

  if (mask & FIRST) {
    *cpmix = -t*d2gdt2;
  }

  if(mask & SECOND) {
    double d3gdt3 = D3GDT3;

    *dt = -t*d3gdt3 - d2gdt2;
  }

  if(mask & THIRD) {
    double d3gdrdt2[NR];

    fillD3GDRDT2
    for (i=0; i<NR; i++) dx[i] = -t*d3gdrdt2[i];
  }

}

void 
vmixKal(int mask, double t, double p, double *x, 
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
  double *r = x;
  double xk  = r[0];
  double xvc = (r[1]+r[2])/4.0;
  double xna = 1.0 - r[0] - r[1]/4.0 - r[2]/2.0;
  double xca = r[2]/4.0;

  if (xk  <= DBL_EPSILON) xk  = DBL_EPSILON;
  if (xvc <= DBL_EPSILON) xvc = DBL_EPSILON;
  if (xna <= DBL_EPSILON) xna = DBL_EPSILON;
  if (xca <= DBL_EPSILON) xca = DBL_EPSILON;
  if (xk  >= 1.0 - DBL_EPSILON) xk  = 1.0 - DBL_EPSILON;
  if (xvc >= 1.0 - DBL_EPSILON) xvc = 1.0 - DBL_EPSILON;
  if (xna >= 1.0 - DBL_EPSILON) xna = 1.0 - DBL_EPSILON;
  if (xca >= 1.0 - DBL_EPSILON) xca = 1.0 - DBL_EPSILON;
  
  if (mask & FIRST) {
    *vmix = DGDP;
  }

  if(mask & SECOND) {
    double d2gdrdp[NR];
    int i;

    fillD2GDRDP
    for (i=0; i<NR; i++) dx[i] = d2gdrdp[i];
  }

  if(mask & THIRD) {
    double d3gdr2dp[NR][NR];
    int i, j;

    fillD3GDR2DP
    for (i=0; i<NR; i++) for (j=0; j<NR; j++) dx2[i][j] = d3gdr2dp[i][j]; 
  }

  if(mask & FOURTH) {
    double d2gdtdp = D2GDTDP;

    *dt = d2gdtdp;
  }

  if(mask & FIFTH) {
    double d2gdp2 = D2GDP2;

    *dp = d2gdp2;
  }

  if(mask & SIXTH) {
    double d3gdt2dp = D3GDT2DP;

    *dt2 = d3gdt2dp;
  }

  if(mask & SEVENTH) {
    double d3gdtdp2 = D3GDTDP2;

    *dtdp = d3gdtdp2;
  }

  if(mask & EIGHTH) {
    double d3gdp3 = D3GDP3;

    *dp2 = d3gdp3;
  }

  if(mask & NINTH) {
    double d3gdrdtdp[NR];
    int i;

    fillD3GDRDTDP
    for (i=0; i<NR; i++) dxdt[i]=d3gdrdtdp[i];
  }

  if(mask & TENTH) {
    double d3gdrdp2[NR];
    int i;

    fillD3GDRDP2
    for (i=0; i<NR; i++) dxdp[i]=d3gdrdp2[i];
  }

}

/* end of file KALSILITE.C */
