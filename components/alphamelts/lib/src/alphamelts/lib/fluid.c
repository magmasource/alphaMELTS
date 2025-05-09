const char *fluid_ver(void) { return "$Id: fluid.c,v 1.4 2009/04/16 16:35:23 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: fluid.c,v $
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
 * Revision 1.4  1997/06/21  22:49:57  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 1.3  1997/05/03  20:23:35  ghiorso
 * *** empty log message ***
 *
 * Revision 1.2  1997/03/27  17:03:38  ghiorso
 * *** empty log message ***
 *
 * Revision 1.1  1996/09/24  20:38:56  ghiorso
 * New modules created for MELTS 3.0.x
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Collection of functions to implement the Pitzer equations for H2O and CO2
**      (file: FLUID.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  February 10, 1995 Original Version
**              Canabalized from WATER.C
**--
*/


#ifdef DEBUG
#undef DEBUG
#endif

#ifdef DEBUG

#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define ABS(x)   ((x) < 0 ? -(x) : (x))

void whaar(double p, double t, double *gH2O, double *hH2O, double *sH2O,
  double *cpH2O, double *dcpdtH2O, double *vH2O, double *dvdtH2O,
  double *dvdpH2O, double *d2vdt2H2O, double *d2vdtdpH2O, double *d2vdp2H2O);

void fluidPhase (
  double    t,       /* Input: temperature in kelvins               */
  double    p,       /* Input: pressure in bars                     */
  double   *x,       /* Input: [NA] Composition in mole fractions   */ 
  double   *g,       /* Mask: 000000000000000000000001 [scaler]     */
  double   *dgdx,    /* Mask: 000000000000000000000010 [NA]         */
  double  **d2gdx2,  /* Mask: 000000000000000000000100 [NA][NA]     */
  double ***d3gdx3,  /* Mask: 000000000000000000001000 [NA][NA][NA] */
  double   *h,       /* Mask: 000000000000000000010000 [scaler]     */
  double   *s,       /* Mask: 000000000000000000100000 [scaler]     */
  double   *dsdx,    /* Mask: 000000000000000001000000 [NA]         */
  double  **d2sdx2,  /* Mask: 000000000000000010000000 [NA][NA]     */
  double   *cp,      /* Mask: 000000000000000100000000 [scaler]     */
  double   *dcpdt,   /* Mask: 000000000000001000000000 [scaler]     */
  double   *dcpdx,   /* Mask: 000000000000010000000000 [NA]         */
  double   *v,       /* Mask: 000000000000100000000000 [scaler]     */
  double   *dvdx,    /* Mask: 000000000001000000000000 [NA]         */
  double   **d2vdx2, /* Mask: 000000000010000000000000 [NA]         */
  double   *dvdt,    /* Mask: 000000000100000000000000 [scaler]     */
  double   *dvdp,    /* Mask: 000000001000000000000000 [scaler]     */
  double   *d2vdt2,  /* Mask: 000000010000000000000000 [scaler]     */
  double   *d2vdtdp, /* Mask: 000000100000000000000000 [scaler]     */
  double   *d2vdp2,  /* Mask: 000001000000000000000000 [scaler]     */
  double   *d2vdxdt, /* Mask: 000010000000000000000000 [NA]         */
  double   *d2vdxdp, /* Mask: 000100000000000000000000 [NA]         */
  double   *a,       /* Mask: 001000000000000000000000 [NA]         */
  double   *mu,      /* Mask: 010000000000000000000000 [NA]         */
  double  **dadx);   /* Mask: 100000000000000000000000 [NA][NA]     */

#define TEST(x) \
   flag = ( ABS(temp - (x)) <= sqrt(sqrt(DBL_EPSILON))*ABS(x) ) \
   ? "OK " : "BAD"

int calculationMode = 1;

int main() {
  float ftemp;
  double t, p, x[2], deltaT, deltaP, temp; 
  double g, *dgdx, **d2gdx2, ***d3gdx3, h, s, *dsdx, **d2sdx2, cp, dcpdt, 
    *dcpdx, v, *dvdx, **d2vdx2, dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2, 
    *d2vdxdt, *d2vdxdp, *a, *mu, **dadx;
  double gRef, *dgdxRef, **d2gdx2Ref, ***d3gdx3Ref, hRef, sRef, *dsdxRef, 
    **d2sdx2Ref, cpRef, dcpdtRef, *dcpdxRef, vRef, *dvdxRef, **d2vdx2Ref, 
    dvdtRef, dvdpRef, d2vdt2Ref, d2vdtdpRef, d2vdp2Ref, *d2vdxdtRef, 
    *d2vdxdpRef, *aRef, *muRef, **dadxRef;
  double gHaar, hHaar, sHaar, cpHaar, dcpdtHaar, vHaar, dvdtHaar, dvdpHaar, 
    d2vdt2Haar, d2vdtdpHaar, d2vdp2Haar;
  char *flag;

  printf("Input t in K   : "); scanf("%f", &ftemp); t    = (double) ftemp;
  printf("Input p in bars: "); scanf("%f", &ftemp); p    = (double) ftemp;
  printf("Input x H2O    : "); scanf("%f", &ftemp); x[0] = (double) ftemp;

  x[1] = 1.0 - x[0];

  printf("%s: %f K, %f bars:\n", (x[0] == 1.0) ? "H2O" : "CO2", t, p);
  fluidPhase(t, p, x, &gRef, dgdxRef, d2gdx2Ref, d3gdx3Ref, &hRef, &sRef, 
    dsdxRef, d2sdx2Ref, &cpRef, &dcpdtRef, dcpdxRef, &vRef, dvdxRef, 
    d2vdx2Ref, &dvdtRef, &dvdpRef, &d2vdt2Ref, &d2vdtdpRef, &d2vdp2Ref, 
    d2vdxdtRef, d2vdxdpRef, aRef, muRef, dadxRef);

  if (x[0] == 1.0) whaar(p, t, &gHaar, &hHaar, &sHaar, &cpHaar, &dcpdtHaar, 
    &vHaar, &dvdtHaar, &dvdpHaar, &d2vdt2Haar, &d2vdtdpHaar, &d2vdp2Haar);

  printf("%s %s %4.2f g  = %g", "   ", "x H2O : ", x[0], gRef);
  if (x[0] == 1.0) printf(" [H: %g]\n", gHaar); else printf("\n");
  printf("%s %s %4.2f h  = %g", "   ", "x H2O : ", x[0], hRef);
  if (x[0] == 1.0) printf(" [H: %g]\n", hHaar); else printf("\n");

  deltaT = sqrt(DBL_EPSILON)*(1.0+ABS(t));
  fluidPhase(t + deltaT, p, x, &g, dgdx, d2gdx2, d3gdx3, &h, &s, dsdx, d2sdx2,  
    &cp, &dcpdt, dcpdx, &v, dvdx, d2vdx2, &dvdt, &dvdp, &d2vdt2, 
    &d2vdtdp, &d2vdp2, d2vdxdt, d2vdxdp, a, mu, dadx);

  temp = - (g-gRef)/deltaT;
  TEST(sRef);
  printf("%s %s %4.2f s  = %g (%g, %s)", flag, "x H2O : ", x[0], sRef,
    temp, "from g");
  if (x[0] == 1.0) printf(" [H: %g]\n", sHaar); else printf("\n");

  temp = t*(s-sRef)/deltaT;
  TEST(cpRef);
  printf("%s %s %4.2f cp  = %g (%g, %s)", flag, "x H2O : ", x[0], cpRef, 
    temp, "from s");
  if (x[0] == 1.0) printf(" [H: %g]\n", cpHaar); else printf("\n");

  temp = (h-hRef)/deltaT;
  TEST(cpRef);
  printf("%s %s %4.2f cp  = %g (%g, %s)", flag, "x H2O : ", x[0], cpRef, 
    temp, "from h");
  if (x[0] == 1.0) printf(" [H: %g]\n", cpHaar); else printf("\n");

  temp = (cp-cpRef)/deltaT;
  TEST(dcpdtRef);
  printf("%s %s %4.2f dcpdt  = %g (%g, %s)", flag, "x H2O : ", x[0], dcpdtRef,
    temp, "from cp");
  if (x[0] == 1.0) printf(" [H: %g]\n", dcpdtHaar); else printf("\n");

  temp = (v-vRef)/deltaT;
  TEST(dvdtRef);
  printf("%s %s %4.2f dvdt  = %g (%g, %s)", flag, "x H2O : ", x[0], dvdtRef,
    temp, "from v");
  if (x[0] == 1.0) printf(" [H: %g]\n", dvdtHaar/10.0); else printf("\n");

  temp = (dvdt-dvdtRef)/deltaT;
  TEST(d2vdt2Ref);
  printf("%s %s %4.2f d2vdt2  = %g (%g, %s)", flag, "x H2O : ", x[0], d2vdt2Ref, 
    temp, "from dvdt");
  if (x[0] == 1.0) printf(" [H: %g]\n", d2vdt2Haar/10.0); else printf("\n");

  temp = (dvdp-dvdpRef)/deltaT;
  TEST(d2vdtdpRef);
  printf("%s %s %4.2f d2vdtdp  = %g (%g, %s)", flag, "x H2O : ", x[0], d2vdtdpRef,
     temp, "from dvdp");
  if (x[0] == 1.0) printf(" [H: %g]\n", d2vdtdpHaar/10.0); else printf("\n");

  deltaP = sqrt(DBL_EPSILON)*(1.0+ABS(p));
  fluidPhase(t, p + deltaP, x, &g, dgdx, d2gdx2, d3gdx3, &h, &s, dsdx, d2sdx2,  
    &cp, &dcpdt, dcpdx, &v, dvdx, d2vdx2, &dvdt, &dvdp, &d2vdt2, 
    &d2vdtdp, &d2vdp2, d2vdxdt, d2vdxdp, a, mu, dadx);

  temp = (g-gRef)/deltaP;
  TEST(vRef);
  printf("%s %s %4.2f v  = %g (%g, %s)", flag, "x H2O : ", x[0], vRef, 
    temp, "from g");
  if (x[0] == 1.0) printf(" [H: %g]\n", vHaar/10.0); else printf("\n");

  temp = (v-vRef)/deltaP;
  TEST(dvdpRef);
  printf("%s %s %4.2f dvdp  = %g (%g, %s)", flag, "x H2O : ", x[0], dvdpRef, 
    temp, "from v");
  if (x[0] == 1.0) printf(" [H: %g]\n", dvdpHaar/10.0); else printf("\n");

  temp = (dvdp-dvdpRef)/deltaP;
  TEST(d2vdp2Ref);
  printf("%s %s %4.2f d2vdp2  = %g (%g, %s)", flag, "x H2O : ", x[0], d2vdp2Ref, 
    temp, "from dvdp");
  if (x[0] == 1.0) printf(" [H: %g]\n", d2vdp2Haar/10.0); else printf("\n");

  temp = (dvdt-dvdtRef)/deltaP;
  TEST(d2vdtdpRef);
  printf("%s %s %4.2f d2vdtdp  = %g (%g, %s)", flag, "x H2O : ", x[0], d2vdtdpRef, 
    temp, "from dvdt");
  if (x[0] == 1.0) printf(" [H: %g]\n", d2vdtdpHaar/10.0); else printf("\n");

  return 0;
}
#else

#include "silmin.h"               /*SILMIN structures include file         */
#undef R

#endif

#define SQUARE(x)  ((x)*(x))
#define CUBIC(x)   ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))
#define QUINTIC(x) ((x)*(x)*(x)*(x)*(x))

#define NA 2
#define NR 1

/*****************************************************************************
 H2O : Pitzer KS and Sterner SM (1994) Equations of state valid continuously
       from zero to extreme pressures for H2O and CO2. J Chem Phys 101: 3111-6

 CO2 : Sterner SM and Pitzer KS (1994) An equation of state for carbon dioxide
       valid from zero to extreme pressures. Contr Mineral Petrol 117: 362-74
******************************************************************************/ 

static void redlichKwong(double t, double p, double b, double a2b, double *Z);

static void idealGas(double t, double rho, double *x, double *Ai, double *dAidt, 
  double *dAidrh, double *d2Aidt2, double *d2Aidtdrh, double *d2Aidrh2, 
  double *d3Aidt3, double *d3Aidt2drh, double *d3Aidtdrh2, double *d3Aidrh3);

#define nTERMS 7  /* 1 + number of temperature terms in PS EOS equation */ 
#define nCOEFF 11 /* 1 + number of c[] coefficients in PS EOS equation  */

typedef struct _coeffEOS {
  const double c[nTERMS];
} CoeffEOS;

typedef struct _refProp {
  const double g;
  const double h;
  const double s;
} RefProp;

void fluidPhase(
  double    t,       /* Input: temperature in kelvins               */
  double    p,       /* Input: pressure in bars                     */
  double   *x,       /* Input: [NA] Composition in mole fractions   */ 
  double   *g,       /* Mask: 000000000000000000000001 [scaler]     */
  double   *dgdx,    /* Mask: 000000000000000000000010 [NA]         */
  double  **d2gdx2,  /* Mask: 000000000000000000000100 [NA][NA]     */
  double ***d3gdx3,  /* Mask: 000000000000000000001000 [NA][NA][NA] */
  double   *h,       /* Mask: 000000000000000000010000 [scaler]     */
  double   *s,       /* Mask: 000000000000000000100000 [scaler]     */
  double   *dsdx,    /* Mask: 000000000000000001000000 [NA]         */
  double  **d2sdx2,  /* Mask: 000000000000000010000000 [NA][NA]     */
  double   *cp,      /* Mask: 000000000000000100000000 [scaler]     */
  double   *dcpdt,   /* Mask: 000000000000001000000000 [scaler]     */
  double   *dcpdx,   /* Mask: 000000000000010000000000 [NA]         */
  double   *v,       /* Mask: 000000000000100000000000 [scaler]     */
  double   *dvdx,    /* Mask: 000000000001000000000000 [NA]         */
  double   **d2vdx2, /* Mask: 000000000010000000000000 [NA]         */
  double   *dvdt,    /* Mask: 000000000100000000000000 [scaler]     */
  double   *dvdp,    /* Mask: 000000001000000000000000 [scaler]     */
  double   *d2vdt2,  /* Mask: 000000010000000000000000 [scaler]     */
  double   *d2vdtdp, /* Mask: 000000100000000000000000 [scaler]     */
  double   *d2vdp2,  /* Mask: 000001000000000000000000 [scaler]     */
  double   *d2vdxdt, /* Mask: 000010000000000000000000 [NA]         */
  double   *d2vdxdp, /* Mask: 000100000000000000000000 [NA]         */
  double   *a,       /* Mask: 001000000000000000000000 [NA]         */
  double   *mu,      /* Mask: 010000000000000000000000 [NA]         */
  double  **dadx)    /* Mask: 100000000000000000000000 [NA][NA]     */
{
  CoeffEOS h2o[nCOEFF] = {
    /*               T^-4            T^-2           T^-1                          T              T^2    */
     { { 0.0e0,        0.0e0,          0.0e0,         0.0e0,         0.0e0,         0.0e0,         0.0e0  } },
     { { 0.0e0,        0.0e0,          0.0e0,  0.24657688e+6, 0.51359951e+2,        0.0e0,         0.0e0  } },
     { { 0.0e0,        0.0e0,          0.0e0,  0.58638965e+0, -.28646939e-2, 0.31375577e-4,        0.0e0  } },
     { { 0.0e0,        0.0e0,          0.0e0,  -.62783840e+1, 0.14791599e-1, 0.35779579e-3, 0.15432925e-7 } },
     { { 0.0e0,        0.0e0,          0.0e0,         0.0e0,  -.42719875e+0, -.16325155e-4,        0.0e0  } },
     { { 0.0e0,        0.0e0,          0.0e0,  0.56654978e+4, -.16580167e+2, 0.76560762e-1,        0.0e0  } },
     { { 0.0e0,        0.0e0,          0.0e0,         0.0e0,  0.10917883e+0,        0.0e0,         0.0e0  } },
     { { 0.0e0, 0.38878656e+13, -.13494878e+9, 0.30916564e+6, 0.75591105e+1,        0.0e0,         0.0e0  } },
     { { 0.0e0,        0.0e0,          0.0e0,  -.65537898e+5, 0.18810675e+3,        0.0e0,         0.0e0  } }, 
     { { 0.0e0, -.14182435e+14, 0.18165390e+9, -.19769068e+6, -.23530318e+2,        0.0e0,         0.0e0  } },
     { { 0.0e0,        0.0e0,          0.0e0,  0.92093375e+5, 0.12246777e+3,        0.0e0,         0.0e0  } }
  };

  CoeffEOS co2[nCOEFF] = {
    /*               T^-4            T^-2           T^-1                          T              T^2    */
     { { 0.0e0,        0.0e0,	       0.0e0,	      0.0e0,	     0.0e0,	    0.0e0,	   0.0e0  } },
     { { 0.0e0,        0.0e0,	       0.0e0,  0.18261340e+7, 0.79224365e+2,	    0.0e0,	   0.0e0  } },
     { { 0.0e0,        0.0e0,	       0.0e0,	      0.0e0,  0.66560660e-4, 0.57152798e-5, 0.30222363e-9 } },
     { { 0.0e0,        0.0e0,	       0.0e0,	      0.0e0,  0.59957845e-2, 0.71669631e-4, 0.62416103e-8 } },
     { { 0.0e0,        0.0e0,	       0.0e0,  -.13270279e+1, -.15210731e+0, 0.53654244e-3, -.71115142e-7 } },
     { { 0.0e0,        0.0e0,	       0.0e0,  0.12456776e+0, 0.49045367e+1, 0.98220560e-2, 0.55962121e-5 } },
     { { 0.0e0,        0.0e0,	       0.0e0,	      0.0e0,  0.75522299e+0,	     0.0e0,	   0.0e0  } },
     { { 0.0e0, -.39344644e+12, 0.90918237e+8, 0.42776716e+6, -.22347856e+2,	     0.0e0,	   0.0e0  } },
     { { 0.0e0,        0.0e0,	       0.0e0,  0.40282608e+3, 0.11971627e+3,	     0.0e0,	   0.0e0  } },
     { { 0.0e0,        0.0e0,	0.22995650e+8, -.78971817e+5, -.63376456e+2,	     0.0e0,	   0.0e0  } },
     { { 0.0e0,        0.0e0,	       0.0e0,  0.95029765e+5, 0.18038071e+2,	     0.0e0,	   0.0e0  } }
  };

  RefProp refH2O = { /* Berman (1988)   - H2O gas */
       -228538.00,   /* gRef 25 c and 1 bar            */
       -241816.00,   /* href 25 c and 1 bar            */
           188.72    /* sref 25 c and 1 bar            */
  };
  RefProp refCO2 = { /* Berman (1988)   - CO2 gas */
       -394341.00,   /* gRef 25 c and 1 bar            */
       -393510.00,   /* href 25 c and 1 bar            */
           213.677   /* sref 25 c and 1 bar            */
  };

   double R       =  83.14241; /* cm^3-bar/mol-K */
   
   double dp, rh, rhn, temp, dpdrh, dpdt, drhdt, d2pdt2, d2pdtdrh, d2pdrh2;
   double term1, dterm1dt, dterm1drh, d2term1dt2, d2term1dtdrh, d2term1drh2,
          d3term1dt3, d3term1dt2drh, d3term1dtdrh2, d3term1drh3;
   double term2, dterm2dt, dterm2drh, d2term2dt2, d2term2dtdrh, d2term2drh2,
          d3term2dt3, d3term2dt2drh, d3term2dtdrh2, d3term2drh3;
   double term2a, dterm2adt, d2term2adt2, d3term2adt3;
   double term2b, dterm2bdt, dterm2bdrh, d2term2bdt2, d2term2bdtdrh, d2term2bdrh2,
          d3term2bdt3, d3term2bdt2drh, d3term2bdtdrh2, d3term2bdrh3;
   double term3, dterm3dt, dterm3drh, d2term3dt2, d2term3dtdrh, d2term3drh2,
          d3term3dt3, d3term3dt2drh, d3term3dtdrh2, d3term3drh3;
   double term3a, dterm3adt, d2term3adt2, d3term3adt3;
   double term3b, dterm3bdt, dterm3bdrh, d2term3bdt2, d2term3bdtdrh, d2term3bdrh2,
          d3term3bdt3, d3term3bdt2drh, d3term3bdtdrh2, d3term3bdrh3;
   double c[nCOEFF], dcdt[nCOEFF], d2cdt2[nCOEFF], d3cdt3[nCOEFF];
   double Ar, dArdt, dArdrh, d2Ardt2, d2Ardtdrh, d2Ardrh2, d3Ardt3, d3Ardt2drh, 
          d3Ardtdrh2, d3Ardrh3;
   double Ai, dAidt, dAidrh, d2Aidt2, d2Aidtdrh, d2Aidrh2, d3Aidt3, d3Aidt2drh, 
          d3Aidtdrh2, d3Aidrh3;
   double A, dAdt, dAdrh, d2Adt2, d2Adtdrh, d2Adrh2, d3Adt3, d3Adt2drh, 
          d3Adtdrh2, d3Adrh3;
   int i, count;

   /**************************************************************************
    set initial guess for rho using redlich-kwong equation
    **************************************************************************/ 

   double a0CO2 = 46.0e6;
   double aCO2  = (73.03 - 0.0714*(t-273.15) + 2.157e-5*pow(t-273.15, (double) 2.0))*10.0e5;
   double bCO2  = 29.7;
   double a0H2O = exp(4.881243 + 0.1823047e-2*(t-273.15) - 0.1712269e-4*pow(t-273.15, (double) 2.0)
                      + 6.479419e-8*pow(t-273.15, (double) 3.0))*10.0e5;
   double aH2O  = (111.3057 + 50.70033*exp(-0.982646e-2*(t-273.15)))*10.0e5;
   double bH2O = 14.6;
   double k    = exp(-11.071 + 5953.0/t - 2.746e6/(t*t) + 4.646e8/(t*t*t));
   double aMix = sqrt(a0H2O*a0CO2) + k*pow(82.05, (double) 2.0)*pow(t, (double) 2.5);

   double aRK  = x[0]*x[0]*aH2O + x[1]*x[1]*aCO2 + 2.0*x[0]*x[1]*aMix;
   double bRK  = x[0]*bH2O + x[1]*bCO2;
   double zRK[2] = {0.0, 0.0};

   redlichKwong(t, p, bRK/(82.05*t), aRK/(bRK*82.05*pow(t, (double) 1.5)), zRK);
   rhn = p/(zRK[0]*R*t);
#ifdef DEBUG
   printf("Redlich-Kwong rho: %14.6g", p/(zRK[0]*R*t));
   if (zRK[0] != zRK[1]) printf(" rho: %14.6g", p/(zRK[1]*R*t));
   printf("\n");
#endif

   /**************************************************************************/

   for (i=1; i<nCOEFF; i++) {
     double cH2O[nCOEFF], dcH2Odt[nCOEFF], d2cH2Odt2[nCOEFF], d3cH2Odt3[nCOEFF]; 
     double cCO2[nCOEFF], dcCO2dt[nCOEFF], d2cCO2dt2[nCOEFF], d3cCO2dt3[nCOEFF];

     cH2O[i] = (h2o[i].c)[1]/(t*t*t*t) +  (h2o[i].c)[2]/(t*t) + (h2o[i].c)[3]/t
       + (h2o[i].c)[4] + (h2o[i].c)[5]*t + (h2o[i].c)[6]*t*t;
     dcH2Odt[i]   = -4.0*(h2o[i].c)[1]/(t*t*t*t*t) - 2.0*(h2o[i].c)[2]/(t*t*t) 
       - (h2o[i].c)[3]/(t*t) + (h2o[i].c)[5] + 2.0*(h2o[i].c)[6]*t;
     d2cH2Odt2[i] = 20.0*(h2o[i].c)[1]/(t*t*t*t*t*t) + 6.0*(h2o[i].c)[2]/(t*t*t*t)
       + 2.0*(h2o[i].c)[3]/(t*t*t) + 2.0*(h2o[i].c)[6];
     d3cH2Odt3[i] = -120.0*(h2o[i].c)[1]/(t*t*t*t*t*t*t) 
       - 24.0*(h2o[i].c)[2]/(t*t*t*t*t) - 6.0*(h2o[i].c)[3]/(t*t*t*t); 

     cCO2[i] = (co2[i].c)[1]/(t*t*t*t) +  (co2[i].c)[2]/(t*t) + (co2[i].c)[3]/t
             + (co2[i].c)[4] + (co2[i].c)[5]*t + (co2[i].c)[6]*t*t;
     dcCO2dt[i]   = -4.0*(co2[i].c)[1]/(t*t*t*t*t) - 2.0*(co2[i].c)[2]/(t*t*t) 
       - (co2[i].c)[3]/(t*t) + (co2[i].c)[5] + 2.0*(co2[i].c)[6]*t;
     d2cCO2dt2[i] = 20.0*(co2[i].c)[1]/(t*t*t*t*t*t) + 6.0*(co2[i].c)[2]/(t*t*t*t)
       + 2.0*(co2[i].c)[3]/(t*t*t) + 2.0*(co2[i].c)[6];
     d3cCO2dt3[i] = -120.0*(co2[i].c)[1]/(t*t*t*t*t*t*t) 
       - 24.0*(co2[i].c)[2]/(t*t*t*t*t) - 6.0*(co2[i].c)[3]/(t*t*t*t); 

     /***** Below is an initial guess at the functional form *****/
     c[i]      = x[0]*cH2O[i]      + x[1]*cCO2[i]; 
     dcdt[i]   = x[0]*dcH2Odt[i]   + x[1]*dcCO2dt[i];
     d2cdt2[i] = x[0]*d2cH2Odt2[i] + x[1]*d2cCO2dt2[i];
     d3cdt3[i] = x[0]*d3cH2Odt3[i] + x[1]*d3cCO2dt3[i];
     /***** Above is an initial guess at the functional form *****/
   }

   dp = DBL_MAX;
   for (count=1, rh = rhn; count<=100 && (fabs(dp) > 10.0*DBL_EPSILON); count++) {
      double pr, dpr, temp1, temp2, dtemp1, dtemp2;

      temp1  = c[3] + 2.0*c[4]*rh + 3.0*c[5]*rh*rh + 4.0*c[6]*rh*rh*rh;
      dtemp1 = 2.0*c[4] + 6.0*c[5]*rh + 12.0*c[6]*rh*rh;

      temp2  = c[2] + c[3]*rh + c[4]*rh*rh + c[5]*rh*rh*rh + c[6]*rh*rh*rh*rh;
      dtemp2 = c[3] + 2.0*c[4]*rh + 3.0*c[5]*rh*rh + 4.0*c[6]*rh*rh*rh;

      pr = rh + c[1]*rh*rh - rh*rh*temp1/(temp2*temp2) + c[7]*rh*rh*exp(-c[8]*rh)
         + c[9]*rh*rh*exp(-c[10]*rh);
      pr *= R*t;
      pr -= p;

      dpr = 1.0 + 2.0*c[1]*rh - 2.0*rh*temp1/(temp2*temp2) 
          - rh*rh*(temp2*temp2*dtemp1 - temp1*2.0*temp2*dtemp2)/(temp2*temp2*temp2*temp2)
          + 2.0*c[7]*rh*exp(-c[8]*rh) - c[7]*c[8]*rh*rh*exp(-c[8]*rh)
          + 2.0*c[9]*rh*exp(-c[10]*rh) - c[9]*c[10]*rh*rh*exp(-c[10]*rh);
      dpr *= R*t;
      
      dp  = - pr/dpr;
      rhn = rh;
      rh += dp;
      if (rh < 0.0) rh = 10.0*DBL_EPSILON;
   }
#ifdef DEBUG
   printf("count, rh, dp, Z: %d %14.6g %14.6g %14.6g\n", count, rh, dp, p/(rh*R*t));
#endif

   /* Calculate the Residual Function */

   term1         = c[2] + c[3]*rh + c[4]*rh*rh + c[5]*rh*rh*rh + c[6]*rh*rh*rh*rh;
   dterm1dt      = dcdt[2] + dcdt[3]*rh + dcdt[4]*rh*rh + dcdt[5]*rh*rh*rh 
                 + dcdt[6]*rh*rh*rh*rh;
   dterm1drh     = c[3] + 2.0*c[4]*rh + 3.0*c[5]*rh*rh + 4.0*c[6]*rh*rh*rh;
   d2term1dt2    = d2cdt2[2] + d2cdt2[3]*rh + d2cdt2[4]*rh*rh + d2cdt2[5]*rh*rh*rh
                 + d2cdt2[6]*rh*rh*rh*rh;
   d2term1dtdrh  = dcdt[3] + 2.0*dcdt[4]*rh + 3.0*dcdt[5]*rh*rh + 4.0*dcdt[6]*rh*rh*rh;
   d2term1drh2   = 2.0*c[4] + 6.0*c[5]*rh + 12.0*c[6]*rh*rh;
   d3term1dt3    = d3cdt3[2] + d3cdt3[3]*rh + d3cdt3[4]*rh*rh + d3cdt3[5]*rh*rh*rh
                 + d3cdt3[6]*rh*rh*rh*rh;
   d3term1dt2drh = d2cdt2[3] + 2.0*d2cdt2[4]*rh + 3.0*d2cdt2[5]*rh*rh 
                 + 4.0*d2cdt2[6]*rh*rh*rh;
   d3term1dtdrh2 = 2.0*dcdt[4] + 6.0*dcdt[5]*rh + 12.0*dcdt[6]*rh*rh;
   d3term1drh3   = 6.0*c[5] + 24.0*c[6]*rh;

   term2a         = c[7]/c[8];
   dterm2adt      = dcdt[7]/c[8] - c[7]*dcdt[8]/(c[8]*c[8]);
   d2term2adt2    = d2cdt2[7]/c[8] - 2.0*dcdt[7]*dcdt[8]/(c[8]*c[8])
                  - c[7]*d2cdt2[8]/(c[8]*c[8]) 
                  + 2.0*c[7]*dcdt[8]*dcdt[8]/(c[8]*c[8]*c[8]);
   d3term2adt3    = d3cdt3[7]/c[8] - d2cdt2[7]*dcdt[8]/(c[8]*c[8])
                  - 2.0*(d2cdt2[7]*dcdt[8] + dcdt[7]*d2cdt2[8])/(c[8]*c[8])
                  + 4.0*dcdt[7]*dcdt[8]*dcdt[8]/(c[8]*c[8]*c[8])
                  - (dcdt[7]*d2cdt2[8] + c[7]*d3cdt3[8])/(c[8]*c[8])
                  + 2.0*c[7]*d2cdt2[8]*dcdt[8]/(c[8]*c[8]*c[8])
                  + 2.0*(dcdt[7]*dcdt[8]*dcdt[8] + c[7]*dcdt[8]*d2cdt2[8]
                      + c[7]*d2cdt2[8]*dcdt[8])/(c[8]*c[8]*c[8])
                  - 6.0*c[7]*dcdt[8]*dcdt[8]*dcdt[8]/(c[8]*c[8]*c[8]*c[8]);

   term2b         = exp(-c[8]*rh) - 1.0;
   dterm2bdt      = - dcdt[8]*rh*exp(-c[8]*rh);
   dterm2bdrh     = - c[8]*exp(-c[8]*rh);
   d2term2bdt2    = (dcdt[8]*dcdt[8]*rh*rh - d2cdt2[8]*rh)*exp(-c[8]*rh);
   d2term2bdtdrh  = (c[8]*dcdt[8]*rh - dcdt[8])*exp(-c[8]*rh);
   d2term2bdrh2   = c[8]*c[8]*exp(-c[8]*rh);
   d3term2bdt3    = (2.0*dcdt[8]*d2cdt2[8]*rh*rh - d3cdt3[8]*rh
                      - dcdt[8]*dcdt[8]*dcdt[8]*rh*rh*rh
                      + d2cdt2[8]*dcdt[8]*rh*rh)*exp(-c[8]*rh);
   d3term2bdt2drh = (2.0*dcdt[8]*dcdt[8]*rh - d2cdt2[8])*exp(-c[8]*rh)
                  - c[8]*(dcdt[8]*dcdt[8]*rh*rh - d2cdt2[8]*rh)*exp(-c[8]*rh);
   d3term2bdtdrh2 = (2.0*c[8]*dcdt[8] - c[8]*c[8]*dcdt[8]*rh)*exp(-c[8]*rh);
   d3term2bdrh3   = - c[8]*c[8]*c[8]*exp(-c[8]*rh);
 
   term2          = term2a*term2b;
   dterm2dt       = dterm2adt*term2b + term2a*dterm2bdt;
   dterm2drh      = term2a*dterm2bdrh;
   d2term2dt2     = d2term2adt2*term2b + 2.0*dterm2adt*dterm2bdt
                  + term2a*d2term2bdt2;
   d2term2dtdrh   = dterm2adt*dterm2bdrh + term2a*d2term2bdtdrh;
   d2term2drh2    = term2a*d2term2bdrh2;
   d3term2dt3     = d3term2adt3*term2b + 3.0*d2term2adt2*dterm2bdt
                  + 3.0*dterm2adt*d2term2bdt2 + term2a*d3term2bdt3;
   d3term2dt2drh  = d2term2adt2*dterm2bdrh + 2.0*dterm2adt*d2term2bdtdrh
                  + term2a*d3term2bdt2drh;
   d3term2dtdrh2  = dterm2adt*d2term2bdrh2 + term2a*d3term2bdtdrh2;
   d3term2drh3    = term2a*d3term2bdrh3;

   term3a         = c[9]/c[10];
   dterm3adt      = dcdt[9]/c[10] - c[9]*dcdt[10]/(c[10]*c[10]);
   d2term3adt2    = d2cdt2[9]/c[10] - 2.0*dcdt[9]*dcdt[10]/(c[10]*c[10])
                  - c[9]*d2cdt2[10]/(c[10]*c[10]) 
                  + 2.0*c[9]*dcdt[10]*dcdt[10]/(c[10]*c[10]*c[10]);
   d3term3adt3    = d3cdt3[9]/c[10] - d2cdt2[9]*dcdt[10]/(c[10]*c[10])
                  - 2.0*(d2cdt2[9]*dcdt[10] + dcdt[9]*d2cdt2[10])/(c[10]*c[10])
                  + 4.0*dcdt[9]*dcdt[10]*dcdt[10]/(c[10]*c[10]*c[10])
                  - (dcdt[9]*d2cdt2[10] + c[9]*d3cdt3[10])/(c[10]*c[10])
                  + 2.0*c[9]*d2cdt2[10]*dcdt[10]/(c[10]*c[10]*c[10])
                  + 2.0*(dcdt[9]*dcdt[10]*dcdt[10] + c[9]*dcdt[10]*d2cdt2[10]
                      + c[9]*d2cdt2[10]*dcdt[10])/(c[10]*c[10]*c[10])
                  - 6.0*c[9]*dcdt[10]*dcdt[10]*dcdt[10]/(c[10]*c[10]*c[10]*c[10]);

   term3b         = exp(-c[10]*rh) - 1.0;
   dterm3bdt      = - dcdt[10]*rh*exp(-c[10]*rh);
   dterm3bdrh     = - c[10]*exp(-c[10]*rh);
   d2term3bdt2    = (dcdt[10]*dcdt[10]*rh*rh - d2cdt2[10]*rh)*exp(-c[10]*rh);
   d2term3bdtdrh  = (c[10]*dcdt[10]*rh - dcdt[10])*exp(-c[10]*rh);
   d2term3bdrh2   = c[10]*c[10]*exp(-c[10]*rh);
   d3term3bdt3    = (2.0*dcdt[10]*d2cdt2[10]*rh*rh - d3cdt3[10]*rh
                      - dcdt[10]*dcdt[10]*dcdt[10]*rh*rh*rh
                      + d2cdt2[10]*dcdt[10]*rh*rh)*exp(-c[10]*rh);
   d3term3bdt2drh = (2.0*dcdt[10]*dcdt[10]*rh - d2cdt2[10])*exp(-c[10]*rh)
                  - c[10]*(dcdt[10]*dcdt[10]*rh*rh - d2cdt2[10]*rh)*exp(-c[10]*rh);
   d3term3bdtdrh2 = (2.0*c[10]*dcdt[10] - c[10]*c[10]*dcdt[10]*rh)*exp(-c[10]*rh);
   d3term3bdrh3   = - c[10]*c[10]*c[10]*exp(-c[10]*rh);

   term3          = term3a*term3b;
   dterm3dt       = dterm3adt*term3b + term3a*dterm3bdt;
   dterm3drh      = term3a*dterm3bdrh;
   d2term3dt2     = d2term3adt2*term3b + 2.0*dterm3adt*dterm3bdt
                  + term3a*d2term3bdt2;
   d2term3dtdrh   = dterm3adt*dterm3bdrh + term3a*d2term3bdtdrh;
   d2term3drh2    = term3a*d2term3bdrh2;
   d3term3dt3     = d3term3adt3*term3b + 3.0*d2term3adt2*dterm3bdt
                  + 3.0*dterm3adt*d2term3bdt2 + term3a*d3term3bdt3;
   d3term3dt2drh  = d2term3adt2*dterm3bdrh + 2.0*dterm3adt*d2term3bdtdrh
                  + term3a*d3term3bdt2drh;
   d3term3dtdrh2  = dterm3adt*d2term3bdrh2 + term3a*d3term3bdtdrh2;
   d3term3drh3    = term3a*d3term3bdrh3;

   Ar         = c[1]*rh + 1.0/term1 - 1.0/c[2] - term2 - term3;
#ifdef DEBUG
   printf("Ares/RT = %f\n", Ar);
#endif
   dArdt      = dcdt[1]*rh - dterm1dt/(term1*term1) + dcdt[2]/(c[2]*c[2]) 
              - dterm2dt - dterm3dt;
   dArdrh     = c[1] - dterm1drh/(term1*term1) - dterm2drh - dterm3drh;
   d2Ardt2    = d2cdt2[1]*rh - d2term1dt2/(term1*term1) 
              + 2.0*dterm1dt*dterm1dt/(term1*term1*term1)
              + d2cdt2[2]/(c[2]*c[2]) - 2.0*dcdt[2]*dcdt[2]/(c[2]*c[2]*c[2])
              - d2term2dt2 - d2term3dt2;
   d2Ardtdrh  = dcdt[1] - d2term1dtdrh/(term1*term1) 
              + 2.0*dterm1drh*dterm1dt/(term1*term1*term1) - d2term2dtdrh 
              - d2term3dtdrh;
   d2Ardrh2   = - d2term1drh2/(term1*term1) 
              + 2.0*dterm1drh*dterm1drh/(term1*term1*term1) - d2term2drh2 
              - d2term3drh2;
   d3Ardt3    = d3cdt3[1]*rh - d3term1dt3/(term1*term1) 
              + 2.0*d2term1dt2*dterm1dt/(term1*term1*term1)
              + 4.0*dterm1dt*d2term1dt2/(term1*term1*term1) 
              - 6.0*dterm1dt*dterm1dt*dterm1dt/(term1*term1*term1*term1)
              + d3cdt3[2]/(c[2]*c[2]) - 2.0*d2cdt2[2]*dcdt[2]/(c[2]*c[2]*c[2])
              - 4.0*dcdt[2]*d2cdt2[2]/(c[2]*c[2]*c[2]) 
              + 6.0*dcdt[2]*dcdt[2]*dcdt[2]/(c[2]*c[2]*c[2]*c[2])
              - d3term2dt3 - d3term3dt3;
   d3Ardt2drh = d2cdt2[1] - d3term1dt2drh/(term1*term1) 
              + 2.0*d2term1dt2*dterm1drh/(term1*term1*term1)
              + 4.0*dterm1dt*d2term1dtdrh/(term1*term1*term1) 
              - 6.0*dterm1dt*dterm1dt*dterm1drh/(term1*term1*term1*term1)
              - d3term2dt2drh - d3term3dt2drh;
   d3Ardtdrh2 = - d3term1dtdrh2/(term1*term1) 
              + 2.0*d2term1drh2*dterm1dt/(term1*term1*term1)
              + 4.0*dterm1drh*d2term1dtdrh/(term1*term1*term1) 
              - 6.0*dterm1drh*dterm1drh*dterm1dt/(term1*term1*term1*term1)
              - d3term2dtdrh2 - d3term3dtdrh2;
   d3Ardrh3   = - d3term1drh3/(term1*term1) 
              + 2.0*d2term1drh2*dterm1drh/(term1*term1*term1)
              + 4.0*dterm1drh*d2term1drh2/(term1*term1*term1) 
              - 6.0*dterm1drh*dterm1drh*dterm1drh/(term1*term1*term1*term1) 
              - d3term2drh3 - d3term3drh3;

   /* calculate ideal gas contributions */ 

   idealGas(t, rh, x, &Ai, &dAidt, &dAidrh, &d2Aidt2, &d2Aidtdrh, &d2Aidrh2, 
     &d3Aidt3, &d3Aidt2drh, &d3Aidtdrh2, &d3Aidrh3);
 
   /* Calculate the sum                                                      */

   A         =                   R*t*Ar         + Ai;
   dAdt      =     R*Ar        + R*t*dArdt      + dAidt;
   dAdrh     =                   R*t*dArdrh     + dAidrh;
   d2Adt2    = 2.0*R*dArdt     + R*t*d2Ardt2    + d2Aidt2;
   d2Adtdrh  =     R*dArdrh    + R*t*d2Ardtdrh  + d2Aidtdrh;
   d2Adrh2   =                   R*t*d2Ardrh2   + d2Aidrh2;
   d3Adt3    = 3.0*R*d2Ardt2   + R*t*d3Ardt3    + d3Aidt3;
   d3Adt2drh = 2.0*R*d2Ardtdrh + R*t*d3Ardt2drh + d3Aidt2drh;
   d3Adtdrh2 =     R*d2Ardrh2  + R*t*d3Ardtdrh2 + d3Aidtdrh2;
   d3Adrh3   =                   R*t*d3Ardrh3   + d3Aidrh3;

   /* calculate g = A + p/rh  and  v = 1/rh */

   p           = rh*rh*dAdrh;
#ifdef DEBUG
   printf("Calculated p: %f\n", p);
#endif
   dpdrh       = 2.0*rh*dAdrh + rh*rh*d2Adrh2;
   dpdt        = rh*rh*d2Adtdrh;
   drhdt       = -dpdt/dpdrh;
   d2pdt2      = rh*rh*d3Adt2drh;
   d2pdtdrh    = 2.0*rh*d2Adtdrh + rh*rh*d3Adtdrh2;
   d2pdrh2     = 2.0*dAdrh + 4.0*rh*d2Adrh2 + rh*rh*d3Adrh3;
 
   *g       = A + p/rh;
   *h       = A + p/rh - t*dAdt; /* g + ts */
   *s       = - dAdt;
   *cp      = -t*d2Adt2 + (t/(rh*rh))*SQUARE(dpdt)/(dpdrh);
   *v       = 1.0/rh;
   *dvdt    = -(1.0/SQUARE(rh))*drhdt;
   *dvdp    = -(1.0/SQUARE(rh))/dpdrh;
       temp = (-d2pdt2 - 2.0*d2pdtdrh*drhdt - d2pdrh2*SQUARE(drhdt))/dpdrh;
   *d2vdt2  = 2.0*SQUARE(drhdt)/CUBIC(rh) - temp/SQUARE(rh);
       temp = (-d2pdtdrh/dpdrh - d2pdrh2*drhdt/dpdrh)/dpdrh;
   *d2vdtdp = -2.0*dpdt/(CUBIC(rh)*SQUARE(dpdrh)) - temp/SQUARE(rh);
   *d2vdp2  = (1.0/SQUARE(rh))*d2pdrh2/CUBIC(dpdrh) 
            + (2.0/CUBIC(rh))/SQUARE(dpdrh);
       temp = - d2Adt2 - t*d3Adt3 + SQUARE(dpdt/rh)/dpdrh
            + t*(2.0*dpdrh*dpdt*d2pdt2 - SQUARE(dpdt)*d2pdtdrh)
              /SQUARE(rh*dpdrh);
   *dcpdt   = temp + t*(*d2vdt2)*dpdt;

               /* -> joules/mol */
   *g       /= 10.0; 
   *h       /= 10.0;
   *s       /= 10.0;
   *cp      /= 10.0;
   *v       /= 10.0;
   *dvdt    /= 10.0;
   *dvdp    /= 10.0;
   *d2vdt2  /= 10.0;
   *d2vdtdp /= 10.0;
   *d2vdp2  /= 10.0; 
   *dcpdt   /= 10.0;
   
         /* correct 298.15 K, 1 bar value to Berman */
   *g += x[0]*(refH2O.g - (-46493.8016496949  )) 
         /* - t*x[0]*(refH2O.s - (   187.572582951064)) + 298.15*x[0]*(refH2O.s - (   187.572582951064)) */;
   *h += x[0]*(refH2O.h - (  9430.96281231262 ));
   *s += x[0]*(refH2O.s - (   187.572582951064));
         /* correct 298.15 K, 1 bar value to Berman */
   *g += x[1]*(refCO2.g - (-394450.0)) 
         /* - t*x[1]*(refCO2.s - (213.698)) + 298.15*x[1]*(refCO2.s - (213.698)) */;
   *h += x[1]*(refCO2.h - (-330735.0  ));
   *s += x[1]*(refCO2.s - (    213.698));
}

static void idealGas(double t, double rho, double *x, double *Ai, double *dAidt, 
  double *dAidrh, double *d2Aidt2, double *d2Aidtdrh, double *d2Aidrh2, 
  double *d3Aidt3, double *d3Aidt2drh, double *d3Aidtdrh2, double *d3Aidrh3)
{
  const double P0 =  1.01325; /* reference pressure in bars */
  const double R  = 83.14241; /* cm^3-bar/mol-K */

  if (x[0] == 1.0) { /* water */
    static struct {
      const double b;
      const double beta;
    } H2O[] = { /* Cooper (1982) */
      {  0.0,         0.0 },
      { 0.134865,     0.0 },
      {-5.005143,     0.0 },
      { 4.006320,     0.0 },
      { 0.012436,   833.0 },
      { 0.973150,  2289.0 },
      { 1.279500,  5009.0 },
      { 0.969560,  5982.0 },
      { 0.248730, 17800.0 }
    };
    double f, DfDrho, D2fDrho2, D3fDrho3, DfDt, D2fDt2, D3fDt3;
    int i;
    
    f  = -1.0;
    f += log(rho*R*t/P0);
    f +=  H2O[1].b;
    f += (H2O[2].b)/t;
    f += (H2O[3].b)*(1.0-log(t));

    DfDrho   = 1.0/rho;
    D2fDrho2 = -1.0/(rho*rho);
    D3fDrho3 = 2.0/(rho*rho*rho);
    DfDt   = 1.0/t - (H2O[2].b)/(t*t) - (H2O[3].b)/t;
    D2fDt2 = -1.0/(t*t) + 2.0*(H2O[2].b)/(t*t*t) + (H2O[3].b)/(t*t);
    D3fDt3 = 2.0/(t*t*t) - 6.0*(H2O[2].b)/(t*t*t*t) - 2.0*(H2O[3].b)/(t*t*t);

    for (i=4; i<=8; i++) {
      double b = H2O[i].beta;
      double term = 1.0 - exp(-b/t);
      double DtermDt = - b*exp(-b/t)/(t*t);
      double D2termDt2 = - b*b*exp(-b/t)/(t*t*t*t) 
                        + 2.0*b*exp(-b/t)/(t*t*t);
      double D3termDt3 = - b*b*b*exp(-b/t)/(t*t*t*t*t*t)
                        + 4.0*b*b*exp(-b/t)/(t*t*t*t*t)
                        + 2.0*b*b*exp(-b/t)/(t*t*t*t*t)
                        - 6.0*b*exp(-b/t)/(t*t*t*t);

      f      += (H2O[i].b)*log(term);
      DfDt   += ((H2O[i].b)/term)*DtermDt;
      D2fDt2 += - ((H2O[i].b)/(term*term))*DtermDt*DtermDt
              + ((H2O[i].b)/term)*D2termDt2;
      D3fDt3 += 2.0*((H2O[i].b)/(term*term*term))*DtermDt*DtermDt*DtermDt
              - ((H2O[i].b)/(term*term))*2.0*DtermDt*D2termDt2
              - ((H2O[i].b)/(term*term))*DtermDt*D2termDt2
              + ((H2O[i].b)/term)*D3termDt3;
    }

    *Ai         = R*t*f;
    *dAidt      = R*f + R*t*DfDt;
    *dAidrh     = R*t*DfDrho;
    *d2Aidt2    = 2.0*R*DfDt + R*t*D2fDt2;
    *d2Aidtdrh  = R*DfDrho;
    *d2Aidrh2   = R*t*D2fDrho2;
    *d3Aidt3    = 3.0*R*D2fDt2 + R*t*D3fDt3;
    *d3Aidt2drh = 0.0;
    *d3Aidtdrh2 = R*D2fDrho2;
    *d3Aidrh3   = R*t*D3fDrho3;

  } else if (x[0] == 0.0) { /* CO2 */
    const double k[] = { /* Berman (1988) units changed to cm^3-bar/mol-K */
    	 93.0*10.0, -13.409e2*10.0, 1.238e5*10.0, 0.0, 
    	 -0.002876*10.0, 6336.2*10.0 };
    	                   /* Stull and Prophet */
    const double Gref = -94265.0*4.184*10.0; /* cm^3-bar/mol    */
    const double Sref =   51.072*4.184*10.0; /* cm^3-bar/mol-K */
    const double Tref = 298.15;

    *Ai = -k[2]/(2.0*t) - 2.0*k[1]*t/sqrt(Tref) - k[2]*t/(2.0*Tref*Tref)
        - k[5]*t/Tref - k[0]*t*log(t) + k[5]*log(t) + R*t*log(R*t*rho)
        + 4.0*k[1]*sqrt(t) - k[4]*t*t/2.0 + k[0]*t*log(Tref) - R*t*log(P0)
        - R*t - Sref*t + k[0]*t + k[4]*Tref*t + k[2]/Tref - k[5]*log(Tref)
        + Tref*Sref + Gref - Tref*k[0] - 2.0*k[1]*sqrt(Tref) 
        - k[4]*Tref*Tref/2.0 + k[5];
    *dAidt = 2.0*k[1]/sqrt(t) + k[2]/(2.0*t*t) + k[5]/t - k[0]*log(t)
        + R*log(R*t*rho) - k[4]*t - 2.0*k[1]/sqrt(Tref) 
        - k[2]/(2.0*Tref*Tref) - k[5]/Tref + k[0]*log(Tref)
        - R*log(P0) - Sref + k[4]*Tref;
    *dAidrh = R*t/rho;
    *d2Aidt2 = R/t - k[0]/t - k[1]/pow(t, (double) 3.0/2.0) 
        - k[2]/(t*t*t) - k[5]/(t*t) - k[4];
    *d2Aidtdrh = R/rho;
    *d2Aidrh2 = - R*t/(rho*rho);
    *d3Aidt3 = -R/(t*t) + k[0]/(t*t) + 3.0*k[1]/(2.0*pow(t, (double) 5.0/2.0))
        + 3.0*k[2]/(t*t*t*t) + 2.0*k[5]/(t*t*t);
    *d3Aidt2drh = 0.0;
    *d3Aidtdrh2 = - R/(rho*rho);
    *d3Aidrh3 = 2.0*R*t/(rho*rho*rho);
     
  } else {
  	double xTemp[2], AiTemp, dAidtTemp, dAidrhTemp, d2Aidt2Temp, 
  	  d2AidtdrhTemp, d2Aidrh2Temp, d3Aidt3Temp, d3Aidt2drhTemp, 
  	  d3Aidtdrh2Temp, d3Aidrh3Temp;
  	
  	xTemp[0] = 1.0; xTemp[1] = 0.0; 
  	idealGas(t, rho, xTemp, &AiTemp, &dAidtTemp, &dAidrhTemp, &d2Aidt2Temp, 
  	  &d2AidtdrhTemp, &d2Aidrh2Temp, &d3Aidt3Temp, &d3Aidt2drhTemp, 
  	  &d3Aidtdrh2Temp, &d3Aidrh3Temp);
     
  	*Ai         = x[0]*AiTemp;
    *dAidt      = x[0]*dAidtTemp;
    *dAidrh     = x[0]*dAidrhTemp;
    *d2Aidt2    = x[0]*d2Aidt2Temp;
    *d2Aidtdrh  = x[0]*d2AidtdrhTemp;
    *d2Aidrh2   = x[0]*d2Aidrh2Temp;
    *d3Aidt3    = x[0]*d3Aidt3Temp;
    *d3Aidt2drh = x[0]*d3Aidt2drhTemp;
    *d3Aidtdrh2 = x[0]*d3Aidtdrh2Temp;
    *d3Aidrh3   = x[0]*d3Aidrh3Temp;
  
    xTemp[0] = 0.0; xTemp[1] = 1.0; 
  	idealGas(t, rho, xTemp, &AiTemp, &dAidtTemp, &dAidrhTemp, &d2Aidt2Temp, 
  	  &d2AidtdrhTemp, &d2Aidrh2Temp, &d3Aidt3Temp, &d3Aidt2drhTemp, 
  	  &d3Aidtdrh2Temp, &d3Aidrh3Temp);
     
  	*Ai         += x[1]*AiTemp;
    *dAidt      += x[1]*dAidtTemp;
    *dAidrh     += x[1]*dAidrhTemp;
    *d2Aidt2    += x[1]*d2Aidt2Temp;
    *d2Aidtdrh  += x[1]*d2AidtdrhTemp;
    *d2Aidrh2   += x[1]*d2Aidrh2Temp;
    *d3Aidt3    += x[1]*d3Aidt3Temp;
    *d3Aidt2drh += x[1]*d3Aidt2drhTemp;
    *d3Aidtdrh2 += x[1]*d3Aidtdrh2Temp;
    *d3Aidrh3   += x[1]*d3Aidrh3Temp;
	
  }
}

static void redlichKwong(double t, double p, double b, double a2b, double *Z)
{
  double n, m, arg;

  if (a2b <= 0.0) a2b = 0.001;

  n = b*p*(a2b-b*p-1.0)/3.0 - a2b*b*p*b*p - 2.0/27.0;
  m = b*p*(a2b-b*p-1.0) - 1.0/3.0;

  arg = n*n/4.0 + m*m*m/27.0;
  if        (arg  > 0.0) {     /* Edminster (1968) - CASE I  */
    double term1 = - n/2.0 + sqrt(arg);
    double term2 = - n/2.0 - sqrt(arg);
    arg  = (term1 >= 0.0) ?   pow( term1, (double) (1.0/3.0))
                          : - pow(-term1, (double) (1.0/3.0));
    arg += (term2 >= 0.0) ?   pow( term2, (double) (1.0/3.0))
                          : - pow(-term2, (double) (1.0/3.0));
    arg += 1.0/3.0;
    Z[0] = arg;
    Z[1] = arg;
    return;

  } else if (arg == 0.0) {     /* Edminster (1968) - CASE II  */
#ifdef DEBUG
    printf("Error in redlichKwong(). Case II encountered!\n");
#endif
    Z[0] = (double) 1.0;
    Z[1] = (double) 1.0;
    return;

  } else if (arg  < 0.0) {     /* Edminster (1968) - CASE III */
    double PI = acos((double) -1.0);
    double cosPhi = (n > 0) ? - sqrt(- (n*n/4.0) / (m*m*m/27.0) )
                             :   sqrt(- (n*n/4.0) / (m*m*m/27.0) );
    double phi = acos(cosPhi);
    double r1  = 2.0*sqrt(-m/3.0)*cos(phi/3.0) + 1.0/3.0;
    double r2  = 2.0*sqrt(-m/3.0)*cos(phi/3.0 + 2.0*PI/3.0) + 1.0/3.0;
    double r3  = 2.0*sqrt(-m/3.0)*cos(phi/3.0 + 4.0*PI/3.0) + 1.0/3.0;
    if (r1 > r2 && r1 > r3) Z[0] = r1; /* gas */
    if (r2 > r1 && r2 > r3) Z[0] = r2;
    if (r3 > r1 && r3 > r2) Z[0] = r3;
	if (r1 < r2 && r1 < r3) Z[1] = r1; /* liquid */
    if (r2 < r1 && r2 < r3) Z[1] = r2;
    if (r3 < r1 && r3 < r2) Z[1] = r3;
    return;
  }
  
#ifdef DEBUG
  printf("Error in redlichKwong(). Bad exit value!\n");
#endif

  Z[0] = (double) 1.0;
  Z[1] = (double) 1.0;
  return;;
}

/* END-OF-FILE: FLUID.C */
