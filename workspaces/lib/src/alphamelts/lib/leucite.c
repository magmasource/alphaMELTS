const char *leucite_ver(void) { return "$Id: leucite.c,v 1.4 2007/03/12 20:06:36 ghiorso Exp $"; }

/*
MELTS Source Code: RCS $Log: leucite.c,v $
MELTS Source Code: RCS Revision 1.4  2007/03/12 20:06:36  ghiorso
MELTS Source Code: RCS Changed fO2 inclusion criteria and adjusted some mineral endmember inclusion
MELTS Source Code: RCS tolerances.
MELTS Source Code: RCS
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
 * Revision 1.4  1997/06/21  22:49:49  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 1.3  1997/05/03  20:23:28  ghiorso
 * *** empty log message ***
 *
 * Revision 1.2  1997/03/27  17:03:32  ghiorso
 * *** empty log message ***
 *
 * Revision 1.1  1996/09/24  20:38:54  ghiorso
 * New modules created for MELTS 3.0.x
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Routines to compute leucite solution properties 
**      (file: LEUCITE.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  June 19, 1993 Original Version
**      V2.0-1  Mark S. Ghiorso  May 10, 1994
**              (1) Began modifications for new isenthalpic, isentropic,
**                  isochoric derivatives
**      V2.0-2  Mark S. Ghiorso  May 17, 1994
**              (1) Changed definition of the enthalpy of mixing
**
*/

#include "silmin.h"  /* Structure definitions foor SILMIN package */

#define SQUARE(x)  ((x)*(x))
#define CUBE(x)    ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

/*
 *=============================================================================
 * Leucite regular solution for K,Na mixing
 */

#define WNAK   7000.0  /* joules */
#define WNAH   7000.0  /* joules */
#define DGR   53000.0  /* joules */

/* Change nothing below this line */

#define G0   0.0
#define GX   (WNAK)
#define GY   (WNAH)
#define GXX -(WNAK)
#define GYY -(WNAH)
#define GXY  (DGR)

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conLeu defines the conversion from m[i], to r[j]
 */
#define NR         2
#define NS         0
#define NA         3
#define FR0(i)     (i == 0) ? 1.0 - x[0] : - x[0]
#define FR1(i)     (i == 1) ? 1.0 - x[1] : - x[1]
#define DFR0DR0(i) - 1.0
#define DFR1DR1(i) - 1.0

/*
 * Global (to this file): derivative definitions
 */
#define R          8.3143
#define S          -R*(xKL*log(xKL)+xNaL*log(xNaL)+xH2OL*log(xH2OL) \
                       +3.0*xKS*log(xKS)/2.0+3.0*xNaS*log(xNaS)/2.0 \
                       +3.0*xVcS*log(xVcS)/2.0 \
                       -r[1]*((3.0/2.0)*log(2.0/3.0)+(1.0/2.0)*log(1.0/2.0)) )
#define H          (G0)+(GX)*r[0]+(GY)*r[1]+(GXX)*r[0]*r[0] \
                       +(GYY)*r[1]*r[1]+(GXY)*r[0]*r[1]
#define V          0.0
#define G          (H) - t*(S) + (p-1.0)*(V)

#define DGDR0      R*t*((1.0-r[1])*log(xKL) - (1.0-r[1])*log(xNaL) \
                        + r[1]*log(xKS) - r[1]*log(xNaS) ) \
                     + (GX) + 2.0*(GXX)*r[0] + (GXY)*r[1]
#define DGDR1      R*t*(-r[0]*log(xKL) - (1.0-r[0])*log(xNaL) + log(xH2OL) \
                        +r[0]*log(xKS) + (1.0-r[0])*log(xNaS) -log(xVcS)  \
                        - ((3.0/2.0)*log(2.0/3.0)+(1.0/2.0)*log(1.0/2.0)) ) \
                     + (GY) + 2.0*(GYY)*r[1] + (GXY)*r[0]
#define DGDT       -(S)
#define DGDP       (V)

#define D2GDR0R0   R*t*((1.0-r[1])*(1.0-r[1])/xKL + (1.0-r[1])*(1.0-r[1])/xNaL \
                      + 2.0*r[1]*r[1]/(3.0*xKS) + 2.0*r[1]*r[1]/(3.0*xNaS)) \
                   + 2.0*(GXX)
#define D2GDR0R1   (GXY)
#define D2GDR1R1   R*t*(r[0]*r[0]/xKL + (1.0-r[0])*(1.0-r[0])/xNaL + 1.0/xH2OL \
                      + 2.0*r[0]*r[0]/(3.0*xKS) + 2.0*(1.0-r[0])*(1.0-r[0])/(3.0*xNaS) \
                      + 2.0/(3.0*xVcS) ) + 2.0*(GYY)
#define D2GDR0DT   R*((1.0-r[1])*log(xKL) - (1.0-r[1])*log(xNaL) \
                      + r[1]*log(xKS) - r[1]*log(xNaS) )
#define D2GDR1DT   R*(-r[0]*log(xKL) - (1.0-r[0])*log(xNaL) + log(xH2OL) \
                      + r[0]*log(xKS) + (1.0-r[0])*log(xNaS) - log(xVcS) \
                      - ((3.0/2.0)*log(2.0/3.0)+(1.0/2.0)*log(1.0/2.0)) )
#define D2GDR0DP   0.0
#define D2GDR1DP   0.0
#define D2GDT2     0.0
#define D2GDTDP    0.0
#define D2GDP2     0.0

#define D3GDR0R0R0 R*t*(-(1.0-r[1])*(1.0-r[1])*(1.0-r[1])/(xKL*xKL) \
                      + (1.0-r[1])*(1.0-r[1])*(1.0-r[1])/(xNaL*xNaL) \
                      - 4.0*r[1]*r[1]*r[1]/(9.0*xKS*xKS) \
                      + 4.0*r[1]*r[1]*r[1]/(9.0*xNaS*xNaS))
#define D3GDR0R0R1 0.0
#define D3GDR0R0DT R*((1.0-r[1])*(1.0-r[1])/xKL + (1.0-r[1])*(1.0-r[1])/xNaL \
                    + 2.0*r[1]*r[1]/(3.0*xKS) + 2.0*r[1]*r[1]/(3.0*xNaS))
#define D3GDR0R0DP 0.0
#define D3GDR0R1R1 0.0
#define D3GDR0R1DT 0.0
#define D3GDR0R1DP 0.0
#define D3GDR1R1R1 R*t*(r[0]*r[0]*r[0]/(xKL*xKL) \
                      + (1.0-r[0])*(1.0-r[0])*(1.0-r[0])/(xNaL*xNaL) \
                      - 1.0/(xH2OL*xH2OL) - 4.0*r[0]*r[0]*r[0]/(9.0*xKS*xKS) \
                      - 4.0*(1.0-r[0])*(1.0-r[0])*(1.0-r[0])/(9.0*xNaS*xNaS) \
                      + 4.0/(9.0*xVcS*xVcS) )
#define D3GDR1R1DT R*(r[0]*r[0]/xKL + (1.0-r[0])*(1.0-r[0])/xNaL + 1.0/xH2OL \
                    + 2.0*r[0]*r[0]/(3.0*xKS) + 2.0*(1.0-r[0])*(1.0-r[0])/(3.0*xNaS) \
                    + 2.0/(3.0*xVcS) )
#define D3GDR1R1DP 0.0
#define D3GDT3     0.0
#define D3GDT2DP   0.0
#define D3GDTDP2   0.0
#define D3GDP3     0.0

#define D3GDR0DT2  0.0
#define D3GDR1DT2  0.0
#define D3GDR0DTDP 0.0
#define D3GDR1DTDP 0.0
#define D3GDR0DP2  0.0
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
testLeu(int mask, double t, double p,
  int na,          /* Expected number of endmember components                 */
  int nr,          /* Expected number of independent compositional variables  */
  char **names,    /* array of strings of names of endmember components       */
  char **formulas, /* array of strings of formulas of endmember components    */
  double *r,       /* array of indepependent compos variables, check bounds   */
  double *m)       /* array of moles of endmember components, check bounds    */
{
  const char *phase = "leucite.c";
  const char *NAMES[NA]    = { "leucite",  "analcime",       "na-leucite" };
  const char *FORMULAS[NA] = { "KAlSi2O6", "NaAlSi2O5(OH)2", "NaAlSi2O6"  };
  int result = TRUE, i;

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
    for (i=0; i<NR; i++) {
      result = result && (r[i] >= 0.0) && (r[i] <= 1.0);
    }
  }
  /* Check bounds on moles of endmember components */
  if (mask & SIXTH) {
    for (i=0; i<NR; i++) result = result && (m[i] >= 0.0);
    result = result && (m[2] >= -m[1]) && (m[2] >= -m[0]);
  }

  return result;
}

void
conLeu(int inpMask, int outMask, double t, double p,
  double *e,      /* comp of leucite in moles of elements                     */
  double *m,      /* comp of leucite in moles of endmember components         */
  double *r,      /* comp of leucite in terms of the independent comp var     */
  double *x,      /* comp of leucite in mole fractions of endmember comp      */
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
      endmember leucite components.
  (2) calculates from a vector of moles of endmember components, one or
      all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
  (3) calculates from a vector of independent compositional variables
      mole fractions of endmember components and/or the Jacobian matrix
      dx[]/dr[]

  In this routine it is assumed that the elements are in the order of atomic 
  numbers and that the order of leucite components has been verified as:
        m[0] = leucite    (KAlSi2O6),
        m[1] = analcime   (NaAlSi2O5(OH)2 .
        m[2] = na-leucite (NaAlSi2O6) . 

  ----------------------------------------------------------------------------*/

  int i, j, k;

  if (inpMask == FIRST && outMask == SECOND) {
    static const int Hy =  1;
    static const int Na = 11;
    static const int K  = 19;

                         /* Projection into the Na-K-H ternary */
    m[0] = e[K ];        /* moles of KAlSi2O6                  */
    m[1] = e[Hy]/2.0;    /* moles of NaAlSi2O5(OH)2            */
    m[2] = e[Na] - m[1]; /* moles of NaAlSi2O6                 */

  } else if (inpMask == SECOND) {
    double sum;

    if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
      printf("Illegal call to conLeu with inpMask = %o and outMask = %o\n",
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
      /* calculates the matrix d3r[i]/dm[j]dm[k]dm[l] using m[] as input */
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
      printf("Illegal call to conLeu with inpMask = %o and outMask = %o\n",
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
    printf("Illegal call to conLeu with inpMask = %o and outMask = %o\n",
      inpMask, outMask);
  }

}

void
dispLeu(int mask, double t, double p, double *x,
  char **formula            /* Mineral formula for interface display MASK: 1 */
  )
{
  double *r = x;
  static char masterString[] = {
/*             1111111111222222222233333333334444444444555555555566666666667
     01234567890123456789012345678901234567890123456789012345678901234567890 */
    "K_.__Na_.__AlSi2O_.__(OH)_.__" };

  if (mask & FIRST) {
    char *string = strcpy((char *) malloc((size_t) (strlen(masterString)+1)*sizeof(char)), masterString);
    double totK, totNa, totO, totOH;
    char n[5];
    int i;

    totK  = r[0];
    totNa = 1.0-r[0];
    totO  = 6.0-r[1];
    totOH = 2.0*r[1];

    (void) snprintf(n, 5, "%4.2f", totK);  for (i=0; i<4; i++) string[ 1+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totNa); for (i=0; i<4; i++) string[ 7+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totO);  for (i=0; i<4; i++) string[17+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totOH); for (i=0; i<4; i++) string[25+i] = n[i];

    *formula = string;
  }
}

void 
actLeu(int mask, double t, double p, double *x, 
  double *a,  /* (pointer to a[]) activities              BINARY MASK: 0001 */
  double *mu, /* (pointer to mu[]) chemical potentials    BINARY MASK: 0010 */
  double **dx /* (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100 */
  )           /* exclusion criteria applied to results if BINARY MASK: 1000 */
{
  double *r    = x;
  double xKL   = (x[0]*(1.0-x[1])         > DBL_EPSILON) ? x[0]*(1.0-x[1])         : DBL_EPSILON;
  double xNaL  = ((1.0-x[0])*(1.0-x[1])   > DBL_EPSILON) ? (1.0-x[0])*(1.0-x[1])   : DBL_EPSILON;
  double xH2OL = (x[1]                    > DBL_EPSILON) ? x[1]                    : DBL_EPSILON;
  double xKS   = (2.0*x[0]*x[1]/3.0       > DBL_EPSILON) ? 2.0*x[0]*x[1]/3.0       : DBL_EPSILON;
  double xNaS  = (2.0*(1.0-x[0])*x[1]/3.0 > DBL_EPSILON) ? 2.0*(1.0-x[0])*x[1]/3.0 : DBL_EPSILON;
  double xVcS  = ((1.0-2.0*x[1]/3.0)      > DBL_EPSILON) ? (1.0-2.0*x[1]/3.0)      : DBL_EPSILON;

  double g, dgdr[NR], fr[NA][NR];
  int i, j;
  
  for(i=0; i<NA; i++) fr[i][0] = FR0(i);
  for(i=0; i<NA; i++) fr[i][1] = FR1(i);

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
    d2gdr2[0][1] = D2GDR0R1; d2gdr2[1][0] = d2gdr2[0][1];
    d2gdr2[1][1] = D2GDR1R1;

    for(i=0; i<NA; i++) dfrdr[i][0] = DFR0DR0(i);
    for(i=0; i<NA; i++) dfrdr[i][1] = DFR1DR1(i);

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
       0.05,  /* exclusion criteria on the mole fraction of leucite    */
       0.05,  /* exclusion criteria on the mole fraction of analcime   */
       0.05,  /* exclusion criteria on the mole fraction of na-leucite */
    };
    double x[NA];

    x[0] = xKL;  /* mole fraction of leucite    */ 
    x[1] = xNaS; /* mole fraction of analcime   */
    x[2] = xNaL; /* mole fraction of na-leucite */

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
gmixLeu(int mask, double t, double p, double *x, 
  double *gmix, /* Gibbs energy of mixing             BINARY MASK: 0001 */
  double *dx,   /* (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010 */
  double **dx2, /* (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100 */
  double ***dx3 /* (pointer to dx3[][][]) d3(g)/d(x[])3 NARY MASK: 1000 */
  )
{
  double *r    = x;
  double xKL   = (x[0]*(1.0-x[1])         > DBL_EPSILON) ? x[0]*(1.0-x[1])         : DBL_EPSILON;
  double xNaL  = ((1.0-x[0])*(1.0-x[1])   > DBL_EPSILON) ? (1.0-x[0])*(1.0-x[1])   : DBL_EPSILON;
  double xH2OL = (x[1]                    > DBL_EPSILON) ? x[1]                    : DBL_EPSILON;
  double xKS   = (2.0*x[0]*x[1]/3.0       > DBL_EPSILON) ? 2.0*x[0]*x[1]/3.0       : DBL_EPSILON;
  double xNaS  = (2.0*(1.0-x[0])*x[1]/3.0 > DBL_EPSILON) ? 2.0*(1.0-x[0])*x[1]/3.0 : DBL_EPSILON;
  double xVcS  = ((1.0-2.0*x[1]/3.0)      > DBL_EPSILON) ? (1.0-2.0*x[1]/3.0)      : DBL_EPSILON;
  
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
    d2gdr2[0][1] = D2GDR0R1; d2gdr2[1][0] = d2gdr2[0][1];
    d2gdr2[1][1] = D2GDR1R1;

    for (i=0; i<NR; i++) {
       for (j=0; j<NR; j++) dx2[i][j] = d2gdr2[i][j];
    }
  }
  if (mask & FOURTH) {
   double d3gdr3[NR][NR][NR];
   int i, j, k;

    d3gdr3[0][0][0] = D3GDR0R0R0;
    d3gdr3[0][0][1] = D3GDR0R0R1; d3gdr3[0][1][0] = d3gdr3[0][0][1]; d3gdr3[1][0][0] = d3gdr3[0][0][1];
    d3gdr3[0][1][1] = D3GDR0R1R1; d3gdr3[1][0][1] = d3gdr3[0][1][1]; d3gdr3[1][1][0] = d3gdr3[0][1][1];
    d3gdr3[1][1][1] = D3GDR1R1R1;

    for (i=0; i<NR; i++) {
      for (j=0; j<NR; j++) {
        for (k=0; k<NR; k++) dx3[i][j][k] = d3gdr3[i][j][k];
      }
    }
  }

}

void 
hmixLeu(int mask, double t, double p, double *x, 
  double *hmix /* Enthalpy of mixing BINARY MASK: 1 */
  )
{
  double *r    = x;
  double xKL   = (x[0]*(1.0-x[1])         > DBL_EPSILON) ? x[0]*(1.0-x[1])         : DBL_EPSILON;
  double xNaL  = ((1.0-x[0])*(1.0-x[1])   > DBL_EPSILON) ? (1.0-x[0])*(1.0-x[1])   : DBL_EPSILON;
  double xH2OL = (x[1]                    > DBL_EPSILON) ? x[1]                    : DBL_EPSILON;
  double xKS   = (2.0*x[0]*x[1]/3.0       > DBL_EPSILON) ? 2.0*x[0]*x[1]/3.0       : DBL_EPSILON;
  double xNaS  = (2.0*(1.0-x[0])*x[1]/3.0 > DBL_EPSILON) ? 2.0*(1.0-x[0])*x[1]/3.0 : DBL_EPSILON;
  double xVcS  = ((1.0-2.0*x[1]/3.0)      > DBL_EPSILON) ? (1.0-2.0*x[1]/3.0)      : DBL_EPSILON;
  
  *hmix = (G) + t*(S);
}

void 
smixLeu(int mask, double t, double p, double *x, 
  double *smix, /* Entropy of mixing                  BINARY MASK: 001 */
  double *dx,   /* (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010 */
  double **dx2  /* (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100 */
  )
{
  double *r    = x;
  double xKL   = (x[0]*(1.0-x[1])         > DBL_EPSILON) ? x[0]*(1.0-x[1])         : DBL_EPSILON;
  double xNaL  = ((1.0-x[0])*(1.0-x[1])   > DBL_EPSILON) ? (1.0-x[0])*(1.0-x[1])   : DBL_EPSILON;
  double xH2OL = (x[1]                    > DBL_EPSILON) ? x[1]                    : DBL_EPSILON;
  double xKS   = (2.0*x[0]*x[1]/3.0       > DBL_EPSILON) ? 2.0*x[0]*x[1]/3.0       : DBL_EPSILON;
  double xNaS  = (2.0*(1.0-x[0])*x[1]/3.0 > DBL_EPSILON) ? 2.0*(1.0-x[0])*x[1]/3.0 : DBL_EPSILON;
  double xVcS  = ((1.0-2.0*x[1]/3.0)      > DBL_EPSILON) ? (1.0-2.0*x[1]/3.0)      : DBL_EPSILON;

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
    d3gdr2dt[0][1] = D3GDR0R1DT; d3gdr2dt[1][0] = d3gdr2dt[0][1]; 
    d3gdr2dt[1][1] = D3GDR1R1DT; 

    for (i=0; i<NR; i++) {
       for (j=0; j<NR; j++) dx2[i][j] = - d3gdr2dt[i][j];
    }
  }

}

void 
cpmixLeu(int mask, double t, double p, double *x, 
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
vmixLeu(int mask, double t, double p, double *x, 
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
    d3gdr2dp[0][1] = D3GDR0R1DP; d3gdr2dp[1][0] = d3gdr2dp[0][1];
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

/* end of file LEUCITE.C */
