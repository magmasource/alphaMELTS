const char *hornblende_ver(void) { return "$Id: hornblende.c,v 1.3 2006/10/20 00:59:22 ghiorso Exp $"; }

/*
MELTS Source Code: RCS $Log: hornblende.c,v $
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
 * Revision 1.3  1997/06/21  22:49:52  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 1.2  1997/05/03  20:23:31  ghiorso
 * *** empty log message ***
 *
 * Revision 1.1  1997/03/27  17:05:26  ghiorso
 * MELTS 3.0 (March 1997)
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Routines to compute hornblende solution properties 
**      (file: HORNBLENDE.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  March 26, 1997 Original Version
**
*/

#include "silmin.h"  /* Structure definitions foor SILMIN package */

#define SQUARE(x)  ((x)*(x))
#define CUBE(x)    ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

/*
 *=============================================================================
 * hornblende (pargasite-hastingite) reciprocal solution
 */

#define WFEMG   4.0*1.68*4.184*1000.0  /* 4*W12 pyroxenes joules */
#define WFEAL  16.78*1000.0            /* W34   pyroxenes joules */
#define DGR     0.0                    /* joules */

/* Change nothing below this line */

#define G0   0.0
#define GX   (WFEMG)
#define GY   (WFEAL)
#define GXX -(WFEMG)
#define GYY -(WFEAL)
#define GXY  (DGR)

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conHrn defines the conversion from m[i], to r[j]
 */
#define NR         2
#define NS         0
#define NA         3
#define FR0(i)     (i == 1) ? 1.0 - x[0] : - x[0]
#define FR1(i)     (i == 2) ? 1.0 - x[1] : - x[1]
#define DFR0DR0(i) - 1.0
#define DFR1DR1(i) - 1.0

/*
 * Global (to this file): derivative definitions
 */
#define R          8.3143
#define S          -R*(4.0*xMgM12*log(xMgM12) + 4.0*xFe2M12*log(xFe2M12) \
                     + xFe3M3*log(xFe3M3) + xAlM3*log(xAlM3))
#define H          (G0)+(GX)*r[0]+(GY)*r[1]+(GXX)*r[0]*r[0] \
                       +(GYY)*r[1]*r[1]+(GXY)*r[0]*r[1]
#define V          0.0
#define G          (H) - t*(S) + (p-1.0)*(V)

#define DGDR0      R*t*(4.0*log(xFe2M12) - 4.0*log(xMgM12)) \
                     + (GX) + 2.0*(GXX)*r[0] + (GXY)*r[1]
#define DGDR1      R*t*(log(xFe3M3) - log(xAlM3)) \
                     + (GY) + 2.0*(GYY)*r[1] + (GXY)*r[0]
#define DGDT       -(S)
#define DGDP       (V)

#define D2GDR0R0   R*t*(4.0/xFe2M12 + 4.0/xMgM12) + 2.0*(GXX)
#define D2GDR0R1   (GXY)
#define D2GDR1R1   R*t*(1.0/xFe3M3 + 1.0/xAlM3) + 2.0*(GYY)
#define D2GDR0DT   R*(4.0*log(xFe2M12) - 4.0*log(xMgM12))
#define D2GDR1DT   R*(log(xFe3M3) - log(xAlM3))
#define D2GDR0DP   0.0
#define D2GDR1DP   0.0
#define D2GDT2     0.0
#define D2GDTDP    0.0
#define D2GDP2     0.0

#define D3GDR0R0R0 R*t*(-4.0/(xFe2M12*xFe2M12) + 4.0/(xMgM12*xMgM12))
#define D3GDR0R0R1 0.0
#define D3GDR0R0DT R*(4.0/xFe2M12 + 4.0/xMgM12)
#define D3GDR0R0DP 0.0
#define D3GDR0R1R1 0.0
#define D3GDR0R1DT 0.0
#define D3GDR0R1DP 0.0
#define D3GDR1R1R1 R*t*(-1.0/(xFe3M3*xFe3M3) + 1.0/(xAlM3*xAlM3))
#define D3GDR1R1DT R*(1.0/xFe3M3 + 1.0/xAlM3)
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
testHrn(int mask, double t, double p,
  int na,          /* Expected number of endmember components                 */
  int nr,          /* Expected number of independent compositional variables  */
  char **names,    /* array of strings of names of endmember components       */
  char **formulas, /* array of strings of formulas of endmember components    */
  double *r,       /* array of indepependent compos variables, check bounds   */
  double *m)       /* array of moles of endmember components, check bounds    */
{
  const char *phase = "leucite.c";
  const char *NAMES[NA]    = { "pargasite",  "ferropargasite",       "magnesiohastingsite" };
  const char *FORMULAS[NA] = { "NaCa2Mg4AlAl2Si6O22(OH)2", "NaCa2Fe4AlAl2Si6O22(OH)2", 
                         "NaCa2Mg4FeAl2Si6O22(OH)2"  };
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
    for (i=1; i<NA; i++) result = result && (m[i] >= 0.0);
    result = result && (m[0] >= -m[1]) && (m[0] >= -m[2]);
  }

  return result;
}

void
conHrn(int inpMask, int outMask, double t, double p,
  double *e,      /* comp of hornblende in moles of elements                  */
  double *m,      /* comp of hornblende in moles of endmember components      */
  double *r,      /* comp of hornblende in terms of the independent comp var  */
  double *x,      /* comp of hornblende in mole fractions of endmember comp   */
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
        m[0] = pargasite           (NaCa2Mg4AlAl2Si6O22(OH)2),
        m[1] = ferropargasite      (NaCa2Fe4AlAl2Si6O22(OH)2),
        m[2] = magnesiohastingsite (NaCa2Mg4FeAl2Si6O22(OH)2). 

  ----------------------------------------------------------------------------*/

  int i, j, k;

  if (inpMask == FIRST && outMask == SECOND) {
    double sumcat, sumchg, fe2, fe3;
    static const int Na = 11;
    static const int Mg = 12;
    static const int Al = 13;
    static const int Si = 14;
    static const int K  = 19;
    static const int Ca = 20;
    static const int Ti = 22;
    static const int Cr = 24;
    static const int Mn = 25;
    static const int Fe = 26;

    /* Sum the cations and correct the analysis for silica deficiency */
    sumcat  = e[Na] + e[Mg] + e[Al] + e[Si] + e[K] + e[Ca] + e[Ti] + e[Cr] 
            + e[Mn] + e[Fe];
    sumchg  = e[Na] + 2.0*e[Mg] + 3.0*e[Al] + 4.0*e[Si] + e[K] + 2.0*e[Ca] 
            + 4.0*e[Ti] + 3.0*e[Cr] + 2.0*e[Mn];
  
    /* Compute the ferric/ferrous ratio */
    fe3 = 23.0*sumcat/8.0 - sumchg - 2.0*e[Fe];
    fe2 = e[Fe] - fe3;
       
    if (fe3 < 0.01*e[Fe]) { fe3 = 0.01*e[Fe]; fe2 = 0.99*e[Fe]; }
    if (fe2 < 0.01*e[Fe]) { fe2 = 0.01*e[Fe]; fe3 = 0.99*e[Fe]; }

                    /* Projection into the  ternary */
    m[0] = e[Mg];   /* moles of NaCa2Mg4AlAl2Si6O22(OH)2  */
    m[1] = fe2/4.0; /* moles of NaCa2Fe4AlAl2Si6O22(OH)2  */
    m[2] = fe3;     /* moles of NaCa2Mg4FeAl2Si6O22(OH)2  */

  } else if (inpMask == SECOND) {
    double sum;

    if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
      printf("Illegal call to conHrn with inpMask = %o and outMask = %o\n",
        inpMask, outMask);

    for (i=0, sum=0.0; i<NA; i++) sum += m[i];

    if (outMask & THIRD) {
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
                d3m[i][j][k][l]  = -6.0*m[i+1]/QUARTIC(sum);
                d3m[i][j][k][l] += (i+1 == j) ? 2.0/CUBE(sum) : 0.0;
                d3m[i][j][k][l] += (i+1 == k) ? 2.0/CUBE(sum) : 0.0;
                d3m[i][j][k][l] += (i+1 == l) ? 2.0/CUBE(sum) : 0.0;
              }
            }
          }
        }
      }

    }

  } else if (inpMask == THIRD) {

    if (outMask & ~(FOURTH | SEVENTH))
      printf("Illegal call to conHrn with inpMask = %o and outMask = %o\n",
        inpMask, outMask);

    if (outMask & FOURTH) {
      /* Converts a vector of independent compositional variables (r) 
         into a vector of mole fractions of endmember components (x).         */

      for (i=0, x[0]=1.0; i<NR; i++) { x[i+1] = r[i]; x[0] -= r[i]; }
    }

    if (outMask & SEVENTH) {
      /* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
      for (i=0; i<NR; i++) for (j=0; j<NR; j++) dr[i+1][j] = (i == j) ? 1.0 : 0.0;
      for (j=0; j<NR; j++) dr[0][j] = -1.0;
    }

  } else  {
    printf("Illegal call to conHrn with inpMask = %o and outMask = %o\n",
      inpMask, outMask);
  }

}

void
dispHrn(int mask, double t, double p, double *x,
  char **formula            /* Mineral formula for interface display MASK: 1 */
  )
{
  double *r = x;
  static char masterString[] = {
/*             1111111111222222222233333333334444444444555555555566666666667
     01234567890123456789012345678901234567890123456789012345678901234567890 */
    "NaCa2Mg_.__Fe2+_.__Al_.__Fe3+_.__Al2Si6O22(OH)2" };

  if (mask & FIRST) {
    char *string = strcpy((char *) malloc((size_t) (strlen(masterString)+1)*sizeof(char)), masterString);
    double totMg, totFe2, totAl, totFe3;
    char n[5];
    int i;

    totMg  = (1.0-r[0])*4.0;
    totFe2 = r[0]*4.0;
    totAl  = 1.0-r[1];
    totFe3 = r[1];

    (void) snprintf(n, 5, "%4.2f", totMg);  for (i=0; i<4; i++) string[ 7+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totFe2); for (i=0; i<4; i++) string[15+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totAl);  for (i=0; i<4; i++) string[21+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totFe3); for (i=0; i<4; i++) string[29+i] = n[i];

    *formula = string;
  }
}

void 
actHrn(int mask, double t, double p, double *x, 
  double *a,  /* (pointer to a[]) activities              BINARY MASK: 0001 */
  double *mu, /* (pointer to mu[]) chemical potentials    BINARY MASK: 0010 */
  double **dx /* (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100 */
  )           /* exclusion criteria applied to results if BINARY MASK: 1000 */
{
  double *r    = x;
  double xMgM12  = 1.0-r[0];
  double xFe2M12 = r[0];
  double xAlM3   = 1.0-r[1];
  double xFe3M3  = r[1];

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
       0.05,  /* exclusion criteria on the mole fraction of pargasite           */
       0.05,  /* exclusion criteria on the mole fraction of ferropargasite      */
       0.05,  /* exclusion criteria on the mole fraction of magnesiohastingsite */
    };
    double x[NA];

    x[0] = xMgM12;  /* mole fraction of pargasite           */ 
    x[1] = xFe2M12; /* mole fraction of ferropargasite      */
    x[2] = xAlM3;   /* mole fraction of magnesiohastingsite */

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
gmixHrn(int mask, double t, double p, double *x, 
  double *gmix, /* Gibbs energy of mixing             BINARY MASK: 0001 */
  double *dx,   /* (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010 */
  double **dx2, /* (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100 */
  double ***dx3 /* (pointer to dx3[][][]) d3(g)/d(x[])3 NARY MASK: 1000 */
  )
{
  double *r    = x;
  double xMgM12  = 1.0-r[0];
  double xFe2M12 = r[0];
  double xAlM3   = 1.0-r[1];
  double xFe3M3  = r[1];
  
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
hmixHrn(int mask, double t, double p, double *x, 
  double *hmix /* Enthalpy of mixing BINARY MASK: 1 */
  )
{
  double *r    = x;
  double xMgM12  = 1.0-r[0];
  double xFe2M12 = r[0];
  double xAlM3   = 1.0-r[1];
  double xFe3M3  = r[1];
  
  *hmix = (G) + t*(S);
}

void 
smixHrn(int mask, double t, double p, double *x, 
  double *smix, /* Entropy of mixing                  BINARY MASK: 001 */
  double *dx,   /* (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010 */
  double **dx2  /* (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100 */
  )
{
  double *r    = x;
  double xMgM12  = 1.0-r[0];
  double xFe2M12 = r[0];
  double xAlM3   = 1.0-r[1];
  double xFe3M3  = r[1];

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
cpmixHrn(int mask, double t, double p, double *x, 
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
vmixHrn(int mask, double t, double p, double *x, 
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

/* end of file HORNBLENDE.C */
