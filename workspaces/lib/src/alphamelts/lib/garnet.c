const char *garnet_ver(void) { return "$Id: garnet.c,v 1.3 2006/10/20 00:59:22 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: garnet.c,v $
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
MELTS Source Code: RCS Revision 1.2  2005/02/13 22:40:46  cvsaccount
MELTS Source Code: RCS *** empty log message ***
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
 * Revision 3.7  1997/06/21  22:49:56  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.6  1997/05/03  20:23:35  ghiorso
 * *** empty log message ***
 *
 * Revision 3.5  1997/03/27  17:03:38  ghiorso
 * *** empty log message ***
 *
 * Revision 3.4  1996/09/24  20:33:41  ghiorso
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
 * Revision 3.1  1995/08/18  17:31:53  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Routines to compute garnet solution properties 
**      (file: GARNET.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  April 7, 1992 Original Version
**      V1.1-1  Mark S. Ghiorso  July 13, 1992 
**              Added (*display) functions
**      V2.0-1  Mark S. Ghiorso  May 10, 1994
**              (1) Began modifications for new isenthalpic, isentropic,
**                  isochoric derivatives
**      V2.0-2  Mark S. Ghiorso  May 17, 1994
**              (1) Changed definition of the enthalpy of mixing
**	V3.0-1  Paul D. Asimow  July 28, 1995
**		(1) Added d3rdm3 to (*convert) and d3gdx3 to (gmix)
**
*/

#include "silmin.h"  /* Structure definitions foor SILMIN package */

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

/*
 *=============================================================================
 * Garnet solution parameters:
 *
 * Berman, R.G.
 * Mixing properties of Ca-Mg-Fe-Mn garnets
 * American Mineralogist 75, 328-344
 *
 * Berman, R.G., Koziol, A.M.
 * Ternary excess properties of grossular-pyrope-almandine garnet and
 * their influence in geothermobarometry
 * American Mineralogist 76, 1223-1231
 *
 * 1 == Grossular, 2 == Pyrope, 3 == Almandine
 */
#define WH112 21560.0   /* joules     */
#define WS112    18.79  /* joules/K   */
#define WV112     0.10  /* joules/bar */
#define WH122 69200.0   /* joules     */
#define WS122    18.79  /* joules/K   */
#define WV122     0.10  /* joules/bar */
#define WH113 20320.0   /* joules     */
#define WS113     5.08  /* joules/K   */
#define WV113     0.17  /* joules/bar */
#define WH133  2620.0   /* joules     */
#define WS133     5.08  /* joules/K   */
#define WV133     0.09  /* joules/bar */
#define WH223   230.0   /* joules     */
#define WS223     0.0   /* joules/K   */
#define WV223     0.01  /* joules/bar */
#define WH233  3720.0   /* joules     */
#define WS233     0.0   /* joules/K   */
#define WV233     0.06  /* joules/bar */
#define WH123     0.0   /* joules     */
#define WS123     0.0   /* joules/K   */
#define WV123     0.00  /* joules/bar */

#define WG112  (WH112-t*WS112+p*WV112)  
#define WG122  (WH122-t*WS122+p*WV122)  
#define WG113  (WH113-t*WS113+p*WV113)  
#define WG133  (WH133-t*WS133+p*WV133)  
#define WG223  (WH223-t*WS223+p*WV223)  
#define WG233  (WH233-t*WS233+p*WV233)  
#define WG123  (WH123-t*WS123+p*WV123)  

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conFld defines the conversion from m[i], to r[j]
 */
#define NR         2
#define NS         0
#define NA         3
#define FR0(i)     (i == 0) ? 1.0 - xal : - xal
#define FR1(i)     (i == 1) ? 1.0 - xgr : - xgr
#define DFR0DR0(i) - 1.0
#define DFR1DR1(i) - 1.0

/*
 * Global (to this file): derivative definitions
 */
#define R         8.3143
#define S         - 3.0*R*(xgr*log(xgr) + xpy*log(xpy) + xal*log(xal)) + \
                  WS112*xgr*xpy*(xgr+xal/2.0) + WS122*xgr*xpy*(xpy+xal/2.0) + \
                  WS113*xgr*xal*(xgr+xpy/2.0) + WS133*xgr*xal*(xal+xpy/2.0) + \
                  WS223*xpy*xal*(xpy+xgr/2.0) + WS233*xpy*xal*(xal+xgr/2.0) + \
                  WS123*xgr*xpy*xal
#define H         WH112*xgr*xpy*(xgr+xal/2.0) + WH122*xgr*xpy*(xpy+xal/2.0) + \
                  WH113*xgr*xal*(xgr+xpy/2.0) + WH133*xgr*xal*(xal+xpy/2.0) + \
                  WH223*xpy*xal*(xpy+xgr/2.0) + WH233*xpy*xal*(xal+xgr/2.0) + \
                  WH123*xgr*xpy*xal
#define V         WV112*xgr*xpy*(xgr+xal/2.0) + WV122*xgr*xpy*(xpy+xal/2.0) + \
                  WV113*xgr*xal*(xgr+xpy/2.0) + WV133*xgr*xal*(xal+xpy/2.0) + \
                  WV223*xpy*xal*(xpy+xgr/2.0) + WV233*xpy*xal*(xal+xgr/2.0) + \
                  WV123*xgr*xpy*xal
#define G         H - t*(S) + p*(V)

#define DGDR0     3.0*R*t*(log(xal) - log(xpy)) + \
                  WG112*( xgr*xpy/2.0 - xgr*(xgr+xal/2.0)) + \
                  WG122*(-xgr*xpy/2.0 - xgr*(xpy+xal/2.0)) + \
                  WG113*(-xgr*xal/2.0 + xgr*(xgr+xpy/2.0)) + \
                  WG133*( xgr*xal/2.0 + xgr*(xal+xpy/2.0)) + \
                  WG223*(-xpy*xal     + (xpy-xal)*(xpy+xgr/2.0)) + \
                  WG233*( xpy*xal     + (xpy-xal)*(xal+xgr/2.0)) + \
                  WG123*xgr*(xpy-xal)
#define DGDR1     3.0*R*t*(log(xgr) - log(xpy)) + \
                  WG112*( xgr*xpy     + (xpy-xgr)*(xgr+xal/2.0)) + \
                  WG122*(-xgr*xpy     + (xpy-xgr)*(xpy+xal/2.0)) + \
                  WG113*( xgr*xal/2.0 + xal*(xgr+xpy/2.0)) + \
                  WG133*(-xgr*xal/2.0 + xal*(xal+xpy/2.0)) + \
                  WG223*(-xpy*xal/2.0 - xal*(xpy+xgr/2.0)) + \
                  WG233*( xpy*xal/2.0 - xal*(xal+xgr/2.0)) + \
                  WG123*(xpy-xgr)*xal
#define DGDP      (V)

#define D2GDR0R0  3.0*R*t*(1.0/xal + 1.0/xpy) + \
                  xgr*(WG122-WG112+WG133-WG113) + 2.0*(xpy-xal)*(WG233-WG223) \
                  - 2.0*WG223*(xpy+xgr/2.0) - 2.0*WG233*(xal+xgr/2.0) \
                  - 2.0*WG123*xgr
#define D2GDR0R1  3.0*R*t*(1.0/xpy) + \
                  WG112*( (xpy-xgr)/2.0 - (xgr+xal/2.0) - xgr) + \
                  WG122*(-(xpy-xgr)/2.0 - (xpy+xal/2.0) + xgr) + \
                  WG113*(-xal/2.0 + (xgr+xpy/2.0) + xgr/2.0) + \
                  WG133*( xal/2.0 + (xal+xpy/2.0) - xgr/2.0) + \
                  WG223*( xal     - (xpy+xgr/2.0) - (xpy-xal)/2.0) + \
                  WG233*(-xal     - (xal+xgr/2.0) + (xpy-xal)/2.0) + \
                  WG123*(xpy-xal-xgr)
#define D2GDR0DT  3.0*R*(log(xal) - log(xpy)) - (\
                  WS112*( xgr*xpy/2.0 - xgr*(xgr+xal/2.0)) + \
                  WS122*(-xgr*xpy/2.0 - xgr*(xpy+xal/2.0)) + \
                  WS113*(-xgr*xal/2.0 + xgr*(xgr+xpy/2.0)) + \
                  WS133*( xgr*xal/2.0 + xgr*(xal+xpy/2.0)) + \
                  WS223*(-xpy*xal     + (xpy-xal)*(xpy+xgr/2.0)) + \
                  WS233*( xpy*xal     + (xpy-xal)*(xal+xgr/2.0)) + \
                  WS123*xgr*(xpy-xal) )
#define D2GDR0DP  WV112*( xgr*xpy/2.0 - xgr*(xgr+xal/2.0)) + \
                  WV122*(-xgr*xpy/2.0 - xgr*(xpy+xal/2.0)) + \
                  WV113*(-xgr*xal/2.0 + xgr*(xgr+xpy/2.0)) + \
                  WV133*( xgr*xal/2.0 + xgr*(xal+xpy/2.0)) + \
                  WV223*(-xpy*xal     + (xpy-xal)*(xpy+xgr/2.0)) + \
                  WV233*( xpy*xal     + (xpy-xal)*(xal+xgr/2.0)) + \
                  WV123*xgr*(xpy-xal)
#define D2GDR1R1  3.0*R*t*(1.0/xgr + 1.0/xpy) + \
                  2.0*(xpy-xgr)*(WG112-WG122) - 2.0*WG112*(xgr+xal/2.0) \
                  - 2.0*WG122*(xpy+xal/2.0) + xal*(WG113-WG133+WG223-WG233) \
                  - 2.0*WG123*xal
#define D2GDR1DT  3.0*R*(log(xgr) - log(xpy)) - ( \
                  WS112*( xgr*xpy     + (xpy-xgr)*(xgr+xal/2.0)) + \
                  WS122*(-xgr*xpy     + (xpy-xgr)*(xpy+xal/2.0)) + \
                  WS113*( xgr*xal/2.0 + xal*(xgr+xpy/2.0)) + \
                  WS133*(-xgr*xal/2.0 + xal*(xal+xpy/2.0)) + \
                  WS223*(-xpy*xal/2.0 - xal*(xpy+xgr/2.0)) + \
                  WS233*( xpy*xal/2.0 - xal*(xal+xgr/2.0)) + \
                  WS123*(xpy-xgr)*xal )
#define D2GDR1DP  WV112*( xgr*xpy     + (xpy-xgr)*(xgr+xal/2.0)) + \
                  WV122*(-xgr*xpy     + (xpy-xgr)*(xpy+xal/2.0)) + \
                  WV113*( xgr*xal/2.0 + xal*(xgr+xpy/2.0)) + \
                  WV133*(-xgr*xal/2.0 + xal*(xal+xpy/2.0)) + \
                  WV223*(-xpy*xal/2.0 - xal*(xpy+xgr/2.0)) + \
                  WV233*( xpy*xal/2.0 - xal*(xal+xgr/2.0)) + \
                  WV123*(xpy-xgr)*xal
#define D2GDT2    0.0
#define D2GDTDP   0.0
#define D2GDP2    0.0

#define D3GDR0R0R0 3.0*R*t*(1.0/SQUARE(xpy)-1.0/SQUARE(xal)) + \
		   6.0*WG223 - 6.0*WG233
#define D3GDR0R0R1 3.0*R*t*(1.0/SQUARE(xpy)) + WG122 - WG112 \
		   + WG133 - WG113 - 3.0*WG233 + 3.0*WG223 - 2.0*WG123
#define D3GDR0R0DT 3.0*R*(1.0/xal + 1.0/xpy) - ( \
                   xgr*(WS122-WS112+WS133-WS113) \
                   + 2.0*(xpy-xal)*(WS233-WS223) - 2.0*WS223*(xpy+xgr/2.0) \
                   - 2.0*WS233*(xal+xgr/2.0) - 2.0*WS123*xgr )
#define D3GDR0R0DP xgr*(WV122-WV112+WV133-WV113) \
                   + 2.0*(xpy-xal)*(WV233-WV223) - 2.0*WV223*(xpy+xgr/2.0) \
                   - 2.0*WV233*(xal+xgr/2.0) - 2.0*WV123*xgr
#define D3GDR0R1R1 3.0*R*t*(1.0/SQUARE(xpy)) + WG113 - WG133 \
		   + WG223 - WG233 - 3.0*WG112 + 3.0*WG122 - 2.0*WG123
#define D3GDR0R1DT 3.0*R*(1.0/xpy) - ( \
                   WS112*( (xpy-xgr)/2.0 - (xgr+xal/2.0) - xgr) + \
                   WS122*(-(xpy-xgr)/2.0 - (xpy+xal/2.0) + xgr) + \
                   WS113*(-xal/2.0 + (xgr+xpy/2.0) + xgr/2.0) + \
                   WS133*( xal/2.0 + (xal+xpy/2.0) - xgr/2.0) + \
                   WS223*( xal     - (xpy+xgr/2.0) - (xpy-xal)/2.0) + \
                   WS233*(-xal     - (xal+xgr/2.0) + (xpy-xal)/2.0) + \
                   WS123*(xpy-xal-xgr) )
#define D3GDR0R1DP WV112*( (xpy-xgr)/2.0 - (xgr+xal/2.0) - xgr) + \
                   WV122*(-(xpy-xgr)/2.0 - (xpy+xal/2.0) + xgr) + \
                   WV113*(-xal/2.0 + (xgr+xpy/2.0) + xgr/2.0) + \
                   WV133*( xal/2.0 + (xal+xpy/2.0) - xgr/2.0) + \
                   WV223*( xal     - (xpy+xgr/2.0) - (xpy-xal)/2.0) + \
                   WV233*(-xal     - (xal+xgr/2.0) + (xpy-xal)/2.0) + \
                   WV123*(xpy-xal-xgr)
#define D3GDR1R1R1 3.0*R*t*(1.0/SQUARE(xpy)-1.0/SQUARE(xgr)) + \
		   6.0*WG122 - 6.0*WG112
#define D3GDR1R1DT 3.0*R*(1.0/xgr + 1.0/xpy) - ( \
                   2.0*(xpy-xgr)*(WS112-WS122) - 2.0*WS112*(xgr+xal/2.0) \
                   - 2.0*WS122*(xpy+xal/2.0) + xal*(WS113-WS133+WS223-WS233) \
                   - 2.0*WS123*xal )
#define D3GDR1R1DP 2.0*(xpy-xgr)*(WV112-WV122) - 2.0*WV112*(xgr+xal/2.0) \
                   - 2.0*WV122*(xpy+xal/2.0) + xal*(WV113-WV133+WV223-WV233) \
                   - 2.0*WV123*xal
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
testGrn(int mask, double t, double p,
  int na,          /* Expected number of endmember components                 */
  int nr,          /* Expected number of independent compositional variables  */
  char **names,    /* array of strings of names of endmember components       */
  char **formulas, /* array of strings of formulas of endmember components    */
  double *r,       /* array of indepependent compos variables, check bounds   */
  double *m)       /* array of moles of endmember components, check bounds    */
{
  const char *phase = "garnet.c";
  const char *NAMES[NA]    = { "almandine", "grossular", "pyrope" };
  const char *FORMULAS[NA] = { "Fe3Al2Si3O12", "Ca3Al2Si3O12", "Mg3Al2Si3O12" };
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
conGrn(int inpMask, int outMask, double t, double p,
  double *e,      /* comp of garnet in moles of elements                      */
  double *m,      /* comp of garnet in moles of endmember components          */
  double *r,      /* comp of garnet in terms of the independent comp var      */
  double *x,      /* comp of garnet in mole fractions of endmember comp       */
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
      endmember garnet components.
  (2) calculates from a vector of moles of endmember components, one or
      all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
  (3) calculates from a vector of independent compositional variables
      mole fractions of endmember components and/or the Jacobian matrix
      dx[]/dr[]

  In this routine it is assumed that the elements are in the order of atomic 
  numbers and that the order of garnet components has been verified as:
        m[0] = almandine (Fe3Al2Si3O12) ,
        m[1] = grossular (Ca3Al2Si3O12) and
        m[2] = pyrope    (Mg3Al2Si3O12) 

  ----------------------------------------------------------------------------*/

  int i, j, k;

  if (inpMask == FIRST && outMask == SECOND) {
    static const int Mg = 12;
    static const int Ca = 20;
    static const int Fe = 26;

                      /* Projection into the Fe, Ca, Mg triangle */
    m[0] = e[Fe]/3.0; /* moles of Fe3Al2Si3O12                   */
    m[1] = e[Ca]/3.0; /* Moles of Ca3Al2Si3O12                   */
    m[2] = e[Mg]/3.0; /* Moles of Mg3Al2Si3O12                   */

  } else if (inpMask == SECOND) {
    double sum;

    if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
      printf("Illegal call to conGrn with inpMask = %o and outMask = %o\n",
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
      /* Calculates the matrix d3r[i]/dm[j]dm[k]dm[l] using m[] as input      */
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
          for (j=0; j<NA; j++) {
            for (k=0; k<NA; k++) {
	      for (l=0; l<NA; l++) {
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
      printf("Illegal call to conGrn with inpMask = %o and outMask = %o\n",
        inpMask, outMask);

    if (outMask & FOURTH) {
      /* Converts a vector of independent compositional variables (r) 
         into a vector of mole fractions of endmember components (x).         */

      for (i=0, x[2]=1.0; i<NR; i++) { x[i] = r[i]; x[2] -= r[i]; }
    }

    if (outMask & SEVENTH) {
      /* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
      for (i=0; i<NR; i++) for (j=0; j<NR; j++) dr[i][j] = (i == j) ? 1.0 : 0.0;
      for (j=0; j<NR; j++) dr[2][j] = -1.0;
    }

  } else  {
    printf("Illegal call to conGrn with inpMask = %o and outMask = %o\n",
      inpMask, outMask);
  }

}

void
dispGrn(int mask, double t, double p, double *x,
  char **formula            /* Mineral formula for interface display MASK: 1 */
  )
{
  double *r = x;
  static char masterString[] = {
/*             1111111111222222222233333333334444444444555555555566666666667
     01234567890123456789012345678901234567890123456789012345678901234567890 */
    "(Ca_.__Fe''_.__Mg_.__)3Al2Si3O12" };

  if (mask & FIRST) {
    char *string = strcpy((char *) malloc((size_t) (strlen(masterString)+1)*sizeof(char)), masterString);
    double totCa, totFe2, totMg;
    char n[5];
    int i;

    totCa  = r[1];
    totFe2 = r[0];
    totMg  = 1.0 - r[0] - r[1];

    (void) snprintf(n, 5, "%4.2f", totCa);  for (i=0; i<4; i++) string[ 3+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totFe2); for (i=0; i<4; i++) string[11+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totMg);  for (i=0; i<4; i++) string[17+i] = n[i];
 
    *formula = string;
  }
}

void 
actGrn(int mask, double t, double p, double *x, 
  double *a,  /* (pointer to a[]) activities              BINARY MASK: 0001 */
  double *mu, /* (pointer to mu[]) chemical potentials    BINARY MASK: 0010 */
  double **dx /* (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100 */
  )           /* exclusion criteria applied to results if BINARY MASK: 1000 */
{
  double xal = (x[0]          > DBL_EPSILON) ? x[0]          : DBL_EPSILON;
  double xgr = (x[1]          > DBL_EPSILON) ? x[1]          : DBL_EPSILON;
  double xpy = (1.0-x[0]-x[1] > DBL_EPSILON) ? 1.0-x[0]-x[1] : DBL_EPSILON;

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
       0.05,  /* exclusion criteria on the mole fraction of almandine  */
       0.05,  /* exclusion criteria on the mole fraction of grossular  */
       0.05   /* exclusion criteria on the mole fraction of pyrope     */
    };
    double x[NA];

    x[0] = xal; /* mole fraction of almandine                          */ 
    x[1] = xgr; /* mole fraction of grossular                          */
    x[2] = xpy; /* mole fraction of pyrope                             */

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
gmixGrn(int mask, double t, double p, double *x, 
  double *gmix, /* Gibbs energy of mixing             BINARY MASK: 0001 */
  double *dx,   /* (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010 */
  double **dx2, /* (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100 */
  double ***dx3 /* (pointer to dx3[][]) d3(g)/d(x[])3 BINARY MASK: 1000 */
  )
{
  double xal = (x[0]          > DBL_EPSILON) ? x[0]          : DBL_EPSILON;
  double xgr = (x[1]          > DBL_EPSILON) ? x[1]          : DBL_EPSILON;
  double xpy = (1.0-x[0]-x[1] > DBL_EPSILON) ? 1.0-x[0]-x[1] : DBL_EPSILON;
  
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

  if(mask & FOURTH) {
    double d3gdr3[NR][NR][NR];
    int i, j, k;

    d3gdr3[0][0][0] = D3GDR0R0R0;
    d3gdr3[0][0][1] = D3GDR0R0R1;
    d3gdr3[0][1][0] = d3gdr3[0][0][1];
    d3gdr3[0][1][1] = D3GDR0R1R1;
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
hmixGrn(int mask, double t, double p, double *x, 
  double *hmix /* Enthalpy of mixing BINARY MASK: 1 */
  )
{
  double xal = (x[0]          > DBL_EPSILON) ? x[0]          : DBL_EPSILON;
  double xgr = (x[1]          > DBL_EPSILON) ? x[1]          : DBL_EPSILON;
  double xpy = (1.0-x[0]-x[1] > DBL_EPSILON) ? 1.0-x[0]-x[1] : DBL_EPSILON;
  
  *hmix = (G) + t*(S);
}

void 
smixGrn(int mask, double t, double p, double *x, 
  double *smix, /* Entropy of mixing                  BINARY MASK: 001 */
  double *dx,   /* (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010 */
  double **dx2  /* (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100 */
  )
{
  double xal = (x[0]          > DBL_EPSILON) ? x[0]          : DBL_EPSILON;
  double xgr = (x[1]          > DBL_EPSILON) ? x[1]          : DBL_EPSILON;
  double xpy = (1.0-x[0]-x[1] > DBL_EPSILON) ? 1.0-x[0]-x[1] : DBL_EPSILON;

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
cpmixGrn(int mask, double t, double p, double *x, 
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
vmixGrn(int mask, double t, double p, double *x, 
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
  double xal = (x[0]          > DBL_EPSILON) ? x[0]          : DBL_EPSILON;
  double xgr = (x[1]          > DBL_EPSILON) ? x[1]          : DBL_EPSILON;
  double xpy = (1.0-x[0]-x[1] > DBL_EPSILON) ? 1.0-x[0]-x[1] : DBL_EPSILON;
  
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

}

/* end of file GARNET.C */
