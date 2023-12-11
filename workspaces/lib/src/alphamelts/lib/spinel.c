const char *spinel_ver(void) { return "$Id: spinel.c,v 1.4 2007/03/12 20:06:36 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: spinel.c,v $
MELTS Source Code: RCS Revision 1.4  2007/03/12 20:06:36  ghiorso
MELTS Source Code: RCS Changed fO2 inclusion criteria and adjusted some mineral endmember inclusion
MELTS Source Code: RCS tolerances.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.3  2006/10/20 00:59:22  ghiorso
MELTS Source Code: RCS (1) Made initial modifications for thread safe code.
MELTS Source Code: RCS (2) Added support for XML I/O in batch mode
MELTS Source Code: RCS (3) Added support for Melts-batch listener for eventual integration into VIGMCS
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2006/08/17 16:47:19  ghiorso
MELTS Source Code: RCS Made modifications to protect strings.  These modifications allow removal
MELTS Source Code: RCS of the flag -fwritable-strings during gcc compilation.  This brings the
MELTS Source Code: RCS code up to gcc 4.x standards.
MELTS Source Code: RCS
MELTS Source Code: RCS Other minor rearrangements and cleanup.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2006/08/15 16:57:36  ghiorso
MELTS Source Code: RCS xMELTS gcc 3.x sources
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.3  2005/06/10 19:00:16  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2004/08/22 02:13:19  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2004/01/02 19:21:49  cvsaccount
MELTS Source Code: RCS CTserver University of Chicago
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.3  2003/09/27 15:35:22  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2003/06/24 16:42:42  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:58:12  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
 * Revision 3.7  1997/06/21  22:49:25  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.6  1997/05/03  20:23:04  ghiorso
 * *** empty log message ***
 *
 * Revision 3.5  1997/03/27  17:03:09  ghiorso
 * *** empty log message ***
 *
 * Revision 3.4  1996/09/24  20:33:19  ghiorso
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
 * Revision 3.1  1995/08/18  19:15:18  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Routines to compute spinel solution properties
**      (file: SPINEL.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  June 4, 1991   Original Version - started
**                               June 6, 1991                    - compiled
**      V1.0-2  Mark S. Ghiorso  June 7, 1991
**              Removed calls to matrix and free_matrix.
**              Reorganized logic in calls to order to avoid passing
**                double and triple pointers to double.
**              Added call to convert_matrix function to pass a two-
**                dimensional array as a pointer to pointer to double
**                to the numerical recipes routine gaussj.
**      V1.0-3  Mark S. Ghiorso  June 9, 1991
**              Corrected the bounds portion of the ordering algorithm.
**      V1.0-4  Mark S. Ghiorso  June 11, 1991
**              Added endmember correction terms
**      V1.0-5  Mark S. Ghiorso  June 13, 1991
**              Modified endmember algorithm to be more generic
**      V1.0-6  Mark S. Ghiorso  June 19, 1991
**              Added testSpn and revised conSpn
**      V1.0-7  Mark S. Ghiorso  July 5, 1991
**              Added exclusion criteria option to actSpn
**      V1.0-8  Mark S. Ghiorso  July 12, 1991
**              Corrected bound constraints in ORDER and altered projection
**              algorithm to account for negative calculated Fe3+ in some cases
**      V1.0-9  Mark S. Ghiorso  July 13, 1991
**              Finished initial guess algorithm for order routine
**      V1.0-10 Mark S. Ghiorso  August 26, 1991
**              Extended conSpn to convert indep comp var -> mole frac comp
**      V1.0-11 Mark S. Ghiorso  September 16, 1991
**              Altered parameter list and structure of conSpn in order
**              to provide derivatives of the independent compositional
**              variables with respect to moles of endmember components.
**              Added several options and tests to conSpn.
**      V1.0-12 Mark S. Ghiorso  September 24, 1991
**              Added a minimal bound constraint on site mole fractions
**              to prevent division by zero in second and third derivative
**              routines
**      V1.0-13 Mark S. Ghiorso  September 24, 1991
**              (1) Altered parameter list and structure of conSpn in order
**                  to provide derivatives of the endmember mole fractions with
**                  respect to independent compositional variables.
**              (2) Altered parameter list and structure of testSpn to
**                  provide for bounds testing on both the independent
**                  compositional variables and on endmember moles
**      V1.0-14 Mark S. Ghiorso  September 24, 1991
**              Revised convergence criteria in ordering routines to be
**              10.0*DBL_EPSILON
**      V1.0-15 Mark S. Ghiorso  January 31, 1992
**              Altered initialization of automatic arrays to conform to the
**              ANSI standard (DEC C RISC)
**      V1.0-16 Mark S. Ghiorso  February 19, 1992
**              Corrected minor casting violations for ANSI C compliance
**      V1.0-17 Mark S. Ghiorso  March 4, 1992
**              Added a maximum iteration limit in the order() routine
**      V1.0-18 Mark S. Ghiorso  March 23, 1992
**              (1) Corrected indexing error in d2x[][] for vmixSpn
**              (2) Replaced illegal calls to pow() function with references
**                  to macros SQUARE and CUBE
**      V1.0-19 Mark S. Ghiorso  April 1, 1992
**              Added a test for round-off created negative mole fractions
**              in order()
**      V1.0-20 Mark S. Ghiorso  April 17, 1992
**              Modified bound constraint algorith in function order()
**      V1.0-21 Mark S. Ghiorso  July 1, 1992
**              Modified exclusion criteria in actPyx
**      V1.1-1  Mark S. Ghiorso  July 13, 1992
**              Added (*display) function
**      V2.0-1  Mark S. Ghiorso  May 10, 1994
**              (1) Began modifications for new isenthalpic, isentropic,
**                  isochoric derivatives
**      V2.0-2  Mark S. Ghiorso  May 17, 1994
**              (1) Changed definition of the enthalpy of mixing
**	V3.0-1  Paul D. Asimow  July 31, 1995
**		(1) added d3rdm3 to (*convert) and d3gdr3 to (*gmix)
**--
*/

#define CALIBRATE_SPINEL

#ifdef DEBUG
#undef DEBUG
#endif

#include "melts_gsl.h"
#include "silmin.h"  /* Structure definitions foor SILMIN package */

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

/*
 *=============================================================================
 * Spinel solution parameters:
 * Sack, R.O., Ghiorso, M.S. (1991)
 *   An internally consistent model for the thermodynamic properties of
 *   Fe-Mg-titanomagnetite-aluminate spinels
 *   Contributions to Mineralogy and Petrology 106: 474-505
 * Sack, R.O., Ghiorso, M.S. (1991)
 *   Chromian spinels and petrogenetic indicators: Thermodynamics and
 *   petrological applications
 *   American Mineralogist (in press)
 */

#define H11     -8.7 * 1000.0 * 4.184 /* joules */
#define W11      4.5 * 1000.0 * 4.184 /* joules */
#define W14     20.8 * 1000.0 * 4.184 /* joules */
#define W1P4    12.4 * 1000.0 * 4.184 /* joules */
#define W15     10.0 * 1000.0 * 4.184 /* joules */
#define W1P5    14.4 * 1000.0 * 4.184 /* joules */
#define W15P    11.7 * 1000.0 * 4.184 /* joules */
#define W1P5P    7.0 * 1000.0 * 4.184 /* joules */
#define W22      3.6 * 1000.0 * 4.184 /* joules */
#define H24      6.55* 1000.0 * 4.184 /* joules */
#define W24U    12.6 * 1000.0 * 4.184 /* joules */
#define W2P4U   10.9 * 1000.0 * 4.184 /* joules */
#define H25      8.05* 1000.0 * 4.184 /* joules */
#define W25PU   15.3 * 1000.0 * 4.184 /* joules */
#define W2P5U   14.4 * 1000.0 * 4.184 /* joules */
#define W45      6.0 * 1000.0 * 4.184 /* joules */
#define W45P     2.0 * 1000.0 * 4.184 /* joules */
#define W4U5PU   1.8 * 1000.0 * 4.184 /* joules */
#define H55      6.25* 1000.0 * 4.184 /* joules */
#define S55      0.0                  /* joules */
#define W55      0.0 * 1000.0 * 4.184 /* joules */
#define HEX     -3.6 * 1000.0 * 4.184 /* joules */
#define HX       2.4 * 1000.0 * 4.184 /* joules */
#define WOCT     2.0 * 1000.0 * 4.184 /* joules */
#define WTET     2.0 * 1000.0 * 4.184 /* joules */
#define H33    -20.0 * 1000.0 * 4.184 /* joules */
#define H23      0.0 * 1000.0 * 4.184 /* joules */
#define W33      7.0 * 1000.0 * 4.184 /* joules */
#define W13     10.0 * 1000.0 * 4.184 /* joules */
#define W1P3    15.2 * 1000.0 * 4.184 /* joules */
#define W13P    11.3 * 1000.0 * 4.184 /* joules */
#define W1P3P    5.9 * 1000.0 * 4.184 /* joules */
#define W2P3U   10.0 * 1000.0 * 4.184 /* joules */
#define W23PU    9.7 * 1000.0 * 4.184 /* joules */
#define W34     12.5 * 1000.0 * 4.184 /* joules */
#define W3P4    10.0 * 1000.0 * 4.184 /* joules */
#define W3PU4U  10.4 * 1000.0 * 4.184 /* joules */
#define W35      0.0 * 1000.0 * 4.184 /* joules */
#define W3P5    10.0 * 1000.0 * 4.184 /* joules */
#define W35P     8.0 * 1000.0 * 4.184 /* joules */
#define W3P5P    0.0 * 1000.0 * 4.184 /* joules */

#define WV1     -0.1250               /* joules/bar  xmg^2*xuv term */
#define WV2      0.1018               /* joules/bar  xmg*xuv^2 term */

/*
 * Definitions of Taylor expansion coefficients in terms of solution
 * parameters. Independent variables are x2, x3, x4, x5, s1, s2, s4
 * (s3 = x3).
 */

#ifndef CALIBRATE_SPINEL

static const double g0    = 0.0;
static const double gx2   = 0.0;
static const double gx3   = 0.0;
static const double gx4   = 0.0;
static const double gx5   = 0.0;
static const double gs1   = 0.5*(-(WOCT) + (W24U) - (W14) + 0.5*(HX) + 0.5*(HEX)
                                                    + (H24));
static const double gs2   = (W11) + (H11);
static const double gs3   = (W13P) - (W13) + (H33);
static const double hs4   = (W15P) - (W15) + (H55);
static const double ss4   = (S55);
static const double gx2x2 = -0.25*((WTET) + (WOCT) + (HX));
static const double gx2x3 = 0.5*((W2P3U) - (W22) -(W1P3) + (W11) + 2.0*(H23));
static const double gx2x4 = 0.5*((WTET) - (W24U) + (W14) + 0.5*(HX) - 0.5*(HEX)
                                                    + (H24));
static const double gx2x5 = 0.5*(-(W22) + (W11) + (W2P5U) - (W1P5) + 2.0*(H25));
static const double gx2s1 = 0.5*((WOCT) - (WTET));
static const double gx2s2 = 0.5*((WTET) - (W22) + (W11) - (WOCT) + 2.0*(W2P4U)
             		  - 2.0*(W1P4) - (HEX) + 2.0*(H24));
static const double gx2s3 = 0.5*((WTET) - (WOCT) + 2.0*(W3PU4U) - 2.0*(W3P4)
             		  - (W2P3U) - (W23PU) + (W22) + (W1P3) + (W13P) - (W11)
             		  - (HEX) + 2.0*(H24) - 2.0*(H23));
static const double gx2s4 = 0.5*((WTET) - (WOCT) - (W11) + (W22) + 2.0*(W4U5PU)
             		  - 2.0*(W45P) - (W2P5U) - (W25PU) + (W1P5) + (W15P) - (HEX)
             		  - 2.0*(H25) + 2.0*(H24));
static const double gx3x3 = - (W13);
static const double gx3x4 = (W34) - (W14) - (W13);
static const double gx3x5 = (W35) - (W15) -(W13);
static const double gx3s1 = 0.5*((W2P3U) - (W22) - (W1P3) + (W11));
static const double gx3s2 = (W1P3) - (W13) - (W11);
static const double gx3s3 = (W33) - (W13P) + (W13);
static const double gx3s4 = (W35P) - (W35) - (W15P) + (W15);
static const double gx4x4 = -(W14);
static const double gx4x5 = (W45) - (W15) - (W14);
static const double gx4s1 = 0.5*((WTET) - (W24U) + (W14) - 0.5*(HX) + 0.5*(HEX)
             		  - (H24));
static const double gx4s2 = -(W11) + (W1P4) - (W14);
static const double gx4s3 = (W3P4) - (W34) - (W13P) + (W13);
static const double gx4s4 = (W45P) - (W45) - (W15P) + (W15);
static const double gx5x5 = -(W15);
static const double gx5s1 = 0.5*(-(W22) + (W11) + (W2P5U) - (W1P5));
static const double gx5s2 = -(W11) + (W1P5) - (W15);
static const double gx5s3 = (W3P5) - (W35) - (W13P) + (W13);
static const double gx5s4 = (W55) - (W15P) + (W15);
static const double gs1s1 = 0.25*(-(WTET) - (WOCT) + (HX));
static const double gs1s2 = 0.5*((WTET) - (W22) + (W11) + (WOCT) -(HX));
static const double gs1s3 = 0.5*((WTET) + (WOCT) - (W2P3U) - (W23PU) + (W22)
             		  + (W1P3) + (W13P) - (W11) - (HX));
static const double gs1s4 = 0.5*((WTET) + (WOCT) + (W22) - (W11) - (W2P5U)
             		  - (W25PU) + (W1P5) + (W15P) - (HX));
static const double gs2s2 = -(W11);
static const double gs2s3 = (W1P3P) - (W1P3) - (W13P) + (W13);
static const double gs2s4 = (W1P5P) - (W1P5) - (W15P) + (W15);
static const double gs3s3 = -(W33);
static const double gs3s4 = (W3P5P) - (W3P5) - (W35P) + (W35);
static const double gs4s4 = -(W55);

#else

/* The code from here ... */

static double g0    = 0.0;
static double gx2   = 0.0;
static double gx3   = 0.0;
static double gx4   = 0.0;
static double gx5   = 0.0;
static double gs1   = 0.5*(-(WOCT) + (W24U) - (W14) + 0.5*(HX) + 0.5*(HEX)
        	    + (H24));
static double gs2   = (W11) + (H11);
static double gs3   = (W13P) - (W13) + (H33);
static double hs4   = (W15P) - (W15) + (H55);
static double ss4   = (S55);
static double gx2x2 = -0.25*((WTET) + (WOCT) + (HX));
static double gx2x3 = 0.5*((W2P3U) - (W22) -(W1P3) + (W11) + 2.0*(H23));
static double gx2x4 = 0.5*((WTET) - (W24U) + (W14) + 0.5*(HX) - 0.5*(HEX)
        	    + (H24));
static double gx2x5 = 0.5*(-(W22) + (W11) + (W2P5U) - (W1P5) + 2.0*(H25));
static double gx2s1 = 0.5*((WOCT) - (WTET));
static double gx2s2 = 0.5*((WTET) - (W22) + (W11) - (WOCT) + 2.0*(W2P4U)
        	    - 2.0*(W1P4) - (HEX) + 2.0*(H24));
static double gx2s3 = 0.5*((WTET) - (WOCT) + 2.0*(W3PU4U) - 2.0*(W3P4)
        	    - (W2P3U) - (W23PU) + (W22) + (W1P3) + (W13P) - (W11)
        	    - (HEX) + 2.0*(H24) - 2.0*(H23));
static double gx2s4 = 0.5*((WTET) - (WOCT) - (W11) + (W22) + 2.0*(W4U5PU)
        	    - 2.0*(W45P) - (W2P5U) - (W25PU) + (W1P5) + (W15P) - (HEX)
        	    - 2.0*(H25) + 2.0*(H24));
static double gx3x3 = - (W13);
static double gx3x4 = (W34) - (W14) - (W13);
static double gx3x5 = (W35) - (W15) -(W13);
static double gx3s1 = 0.5*((W2P3U) - (W22) - (W1P3) + (W11));
static double gx3s2 = (W1P3) - (W13) - (W11);
static double gx3s3 = (W33) - (W13P) + (W13);
static double gx3s4 = (W35P) - (W35) - (W15P) + (W15);
static double gx4x4 = -(W14);
static double gx4x5 = (W45) - (W15) - (W14);
static double gx4s1 = 0.5*((WTET) - (W24U) + (W14) - 0.5*(HX) + 0.5*(HEX)
        	    - (H24));
static double gx4s2 = -(W11) + (W1P4) - (W14);
static double gx4s3 = (W3P4) - (W34) - (W13P) + (W13);
static double gx4s4 = (W45P) - (W45) - (W15P) + (W15);
static double gx5x5 = -(W15);
static double gx5s1 = 0.5*(-(W22) + (W11) + (W2P5U) - (W1P5));
static double gx5s2 = -(W11) + (W1P5) - (W15);
static double gx5s3 = (W3P5) - (W35) - (W13P) + (W13);
static double gx5s4 = (W55) - (W15P) + (W15);
static double gs1s1 = 0.25*(-(WTET) - (WOCT) + (HX));
static double gs1s2 = 0.5*((WTET) - (W22) + (W11) + (WOCT) -(HX));
static double gs1s3 = 0.5*((WTET) + (WOCT) - (W2P3U) - (W23PU) + (W22)
        	    + (W1P3) + (W13P) - (W11) - (HX));
static double gs1s4 = 0.5*((WTET) + (WOCT) + (W22) - (W11) - (W2P5U)
        	    - (W25PU) + (W1P5) + (W15P) - (HX));
static double gs2s2 = -(W11);
static double gs2s3 = (W1P3P) - (W1P3) - (W13P) + (W13);
static double gs2s4 = (W1P5P) - (W1P5) - (W15P) + (W15);
static double gs3s3 = -(W33);
static double gs3s4 = (W3P5P) - (W3P5) - (W35P) + (W35);
static double gs4s4 = -(W55);


static int noSave = FALSE;

void resetValueOfW45andW45p (double newW45, double newW45p) {
    gx2s4 = 0.5*((WTET) - (WOCT) - (W11) + (W22) + 2.0*(W4U5PU)
                - 2.0*newW45p - (W2P5U) - (W25PU) + (W1P5) + (W15P) - (HEX)
                - 2.0*(H25) + 2.0*(H24));
    gx4x5 = newW45 - (W15) - (W14);
    gx4s4 = newW45p - newW45 - (W15P) + (W15);
    noSave = TRUE;
}

void resetValueOfH55 (double newH55) {
    hs4    = (W15P) - (W15) + newH55;
    noSave = TRUE;
}

void resetValueOfS55 (double newS55) {
    ss4    = newS55;
    noSave = TRUE;
}

void resetValueOfW55 (double newW55) {
    gx5s4 =  newW55 - (W15P) + (W15);
    gs4s4 = -newW55;
    noSave = TRUE;
}

#endif

/* ... to here is NOT thread safe.  It is nerver called by the
       except during calibration of the spinel model.           */

/*
 * Global (to this file): variables
 */

#define R  8.3143
#define NR         4    /* Four independent composition variables */
#define NS         3    /* Three ordering parameters s3 = X3      */
#define NA         5    /* Five endmember compositions            */

#define MAX_ITER 200    /* Maximum number of iterations allowed in order */

/*************************************/
/* Statics for Ordering Calculations */
/*************************************/

static MTHREAD_ONCE_T initThreadOBlock = MTHREAD_ONCE_INIT;

static MTHREAD_KEY_T tOldKey;
static MTHREAD_KEY_T pOldKey;
static MTHREAD_KEY_T rOldKey;
static MTHREAD_KEY_T sOldKey;
static MTHREAD_KEY_T ptToD2gds2Key;
static MTHREAD_KEY_T d2gds2Key;
static MTHREAD_KEY_T indexD2gds2Key;
static MTHREAD_KEY_T tOldPureKey;
static MTHREAD_KEY_T pOldPureKey;
static MTHREAD_KEY_T sOldPureKey;
static MTHREAD_KEY_T d2gds2PureKey;

static void freeNSarray(void *NSarray) {
    gsl_vector_free((gsl_vector *) NSarray);
}

static void freePtToD2gds2(void *ptToD2gds2) {
    gsl_matrix_free((gsl_matrix *) ptToD2gds2);
}

static void freeIndexD2gds2(void *indexD2gds2) {
    gsl_permutation_free((gsl_permutation *) indexD2gds2);
}

static void threadOInit(void) {
    MTHREAD_KEY_CREATE(&tOldKey,       free);
    MTHREAD_KEY_CREATE(&pOldKey,       free);
    MTHREAD_KEY_CREATE(&rOldKey,       freeNSarray);
    MTHREAD_KEY_CREATE(&sOldKey,       freeNSarray);
    MTHREAD_KEY_CREATE(&ptToD2gds2Key, freePtToD2gds2);
    MTHREAD_KEY_CREATE(&d2gds2Key,     free);
    MTHREAD_KEY_CREATE(&indexD2gds2Key, freeIndexD2gds2);
    MTHREAD_KEY_CREATE(&tOldPureKey,   free);
    MTHREAD_KEY_CREATE(&pOldPureKey,   free);
    MTHREAD_KEY_CREATE(&sOldPureKey,   freeNSarray);
    MTHREAD_KEY_CREATE(&d2gds2PureKey, freeNSarray);
}

static double getTOld() {
    double *tOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    tOldPt = (double *) MTHREAD_GETSPECIFIC(tOldKey);
    if (tOldPt == NULL) {
        tOldPt  = (double *) malloc(sizeof(double));
        *tOldPt = -9999.0;
        MTHREAD_SETSPECIFIC(tOldKey, (void *) tOldPt);
    }
    return *tOldPt;
}

static void setTOld(double tOld) {
    double *tOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    tOldPt = (double *) MTHREAD_GETSPECIFIC(tOldKey);
    if (tOldPt == NULL) {
        tOldPt  = (double *) malloc(sizeof(double));
        *tOldPt = -9999.0;
        MTHREAD_SETSPECIFIC(tOldKey, (void *) tOldPt);
    }
    *tOldPt = tOld;
}

static double getPOld() {
    double *pOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    pOldPt = (double *) MTHREAD_GETSPECIFIC(pOldKey);
    if (pOldPt == NULL) {
        pOldPt  = (double *) malloc(sizeof(double));
        *pOldPt = -9999.0;
        MTHREAD_SETSPECIFIC(pOldKey, (void *) pOldPt);
    }
    return *pOldPt;
}

static void setPOld(double pOld) {
    double *pOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    pOldPt = (double *) MTHREAD_GETSPECIFIC(pOldKey);
    if (pOldPt == NULL) {
        pOldPt  = (double *) malloc(sizeof(double));
        *pOldPt = -9999.0;
        MTHREAD_SETSPECIFIC(pOldKey, (void *) pOldPt);
    }
    *pOldPt = pOld;
}

static double *getROld() {
    gsl_vector *rOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    rOldPt = (gsl_vector *) MTHREAD_GETSPECIFIC(rOldKey);
    if (rOldPt == NULL) {
        rOldPt = gsl_vector_alloc((size_t) NR);
        gsl_vector_set_all(rOldPt, -9999.0);
        MTHREAD_SETSPECIFIC(rOldKey, (void *) rOldPt);
    }
    return rOldPt->data;
}

static double *getSOld() {
    gsl_vector *sOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    sOldPt = (gsl_vector *) MTHREAD_GETSPECIFIC(sOldKey);
    if (sOldPt == NULL) {
        sOldPt = gsl_vector_alloc((size_t) NS);
        gsl_vector_set_all(sOldPt, 2.0);
        MTHREAD_SETSPECIFIC(sOldKey, (void *) sOldPt);
    }
    return sOldPt->data;
}

static gsl_matrix *getPtToD2gds2() {
    gsl_matrix *ptToD2gds2Pt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    ptToD2gds2Pt = (gsl_matrix *) MTHREAD_GETSPECIFIC(ptToD2gds2Key);
    if (ptToD2gds2Pt == NULL) {
        ptToD2gds2Pt  = gsl_matrix_alloc((size_t) NS, (size_t) NS);
        gsl_matrix_set_zero(ptToD2gds2Pt);
        MTHREAD_SETSPECIFIC(ptToD2gds2Key, (void *) ptToD2gds2Pt);
    }
    return ptToD2gds2Pt;
}

static double **getD2gds2() {
    double **d2gds2Pt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    d2gds2Pt = (double **) MTHREAD_GETSPECIFIC(d2gds2Key);
    if (d2gds2Pt == NULL) {
        int i;
        gsl_matrix *ptToD2gds2Pt = getPtToD2gds2();
        d2gds2Pt  = (double **) malloc((size_t) NS*sizeof(double *));
        for (i=0; i<NS; i++) d2gds2Pt[i] = &(ptToD2gds2Pt->data[i*NS]);
        MTHREAD_SETSPECIFIC(d2gds2Key, (void *) d2gds2Pt);
    }
    return d2gds2Pt;
}

static gsl_permutation *getIndexD2gds2() {
    gsl_permutation *indexD2gds2Pt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    indexD2gds2Pt = (gsl_permutation *) MTHREAD_GETSPECIFIC(indexD2gds2Key);
    if (indexD2gds2Pt == NULL) {
        indexD2gds2Pt = gsl_permutation_alloc((size_t) NS);
        gsl_permutation_init(indexD2gds2Pt);
        MTHREAD_SETSPECIFIC(indexD2gds2Key, (void *) indexD2gds2Pt);
    }
    return indexD2gds2Pt;
}

static double getTOldPure() {
    double *tOldPurePt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    tOldPurePt = (double *) MTHREAD_GETSPECIFIC(tOldPureKey);
    if (tOldPurePt == NULL) {
        tOldPurePt  = (double *) malloc(sizeof(double));
        *tOldPurePt = -9999.0;
        MTHREAD_SETSPECIFIC(tOldPureKey, (void *) tOldPurePt);
    }
    return *tOldPurePt;
}

static void setTOldPure(double tOldPure) {
    double *tOldPurePt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    tOldPurePt = (double *) MTHREAD_GETSPECIFIC(tOldPureKey);
    if (tOldPurePt == NULL) {
        tOldPurePt  = (double *) malloc(sizeof(double));
        *tOldPurePt = -9999.0;
        MTHREAD_SETSPECIFIC(tOldPureKey, (void *) tOldPurePt);
    }
    *tOldPurePt = tOldPure;
}

static double getPOldPure() {
    double *pOldPurePt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    pOldPurePt = (double *) MTHREAD_GETSPECIFIC(pOldPureKey);
    if (pOldPurePt == NULL) {
        pOldPurePt  = (double *) malloc(sizeof(double));
        *pOldPurePt = -9999.0;
        MTHREAD_SETSPECIFIC(pOldPureKey, (void *) pOldPurePt);
    }
    return *pOldPurePt;
}

static void setPOldPure(double pOldPure) {
    double *pOldPurePt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    pOldPurePt = (double *) MTHREAD_GETSPECIFIC(pOldPureKey);
    if (pOldPurePt == NULL) {
        pOldPurePt  = (double *) malloc(sizeof(double));
        *pOldPurePt = -9999.0;
        MTHREAD_SETSPECIFIC(pOldPureKey, (void *) pOldPurePt);
    }
    *pOldPurePt = pOldPure;
}

static double *getSOldPure() {
    gsl_vector *sOldPurePt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    sOldPurePt = (gsl_vector *) MTHREAD_GETSPECIFIC(sOldPureKey);
    if (sOldPurePt == NULL) {
        sOldPurePt = gsl_vector_alloc((size_t) NS);
        gsl_vector_set_all(sOldPurePt, 2.0);
        MTHREAD_SETSPECIFIC(sOldPureKey, (void *) sOldPurePt);
    }
    return sOldPurePt->data;
}

static double *getD2gds2Pure() {
    gsl_vector *d2gds2PurePt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    d2gds2PurePt = (gsl_vector *) MTHREAD_GETSPECIFIC(d2gds2PureKey);
    if (d2gds2PurePt == NULL) {
        d2gds2PurePt  = gsl_vector_alloc((size_t) NS);
        gsl_vector_set_zero(d2gds2PurePt);
        MTHREAD_SETSPECIFIC(d2gds2PureKey, (void *) d2gds2PurePt);
    }
    return d2gds2PurePt->data;
}

/***********************************/
/* Statics for Site Mole Fractions */
/***********************************/

#define DECLARE_SITE_FRACTIONS \
    double xmg2tet, xfe2tet, xal3tet, xfe3tet, xmg2oct, xfe2oct, xal3oct, xfe3oct, xcr3oct, xti4oct;

#define XMG2TET 0
#define XFE2TET 1
#define XAL3TET 2
#define XFE3TET 3
#define XMG2OCT 4
#define XFE2OCT 5
#define XAL3OCT 6
#define XFE3OCT 7
#define XCR3OCT 8
#define XTI4OCT 9

#define NX     10

static MTHREAD_ONCE_T initThreadXBlock = MTHREAD_ONCE_INIT;

static MTHREAD_KEY_T xKey[NX];

static void threadXInit(void) {
    int i;
    for (i=0; i<NX; i++) MTHREAD_KEY_CREATE(&xKey[i], free);
}

static double getX(int n) {
    double *xPt;
    MTHREAD_ONCE(&initThreadXBlock, threadXInit);

    xPt = (double *) MTHREAD_GETSPECIFIC(xKey[n]);
    if (xPt == NULL) {
        xPt = (double *) malloc(sizeof(double)); *xPt = 0.0;
        MTHREAD_SETSPECIFIC(xKey[n], (void *) xPt);
    }
    return *xPt;
}

static void setX(int n, double x) {
    double *xPt;
    MTHREAD_ONCE(&initThreadXBlock, threadXInit);

    xPt = (double *) MTHREAD_GETSPECIFIC(xKey[n]);
    if (xPt == NULL) {
        xPt = (double *) malloc(sizeof(double));
        MTHREAD_SETSPECIFIC(xKey[n], (void *) xPt);
    }
    *xPt = x;
}

#define GET_SITE_FRACTIONS \
    xmg2tet = getX(XMG2TET); \
    xfe2tet = getX(XFE2TET); \
    xal3tet = getX(XAL3TET); \
    xfe3tet = getX(XFE3TET); \
    xmg2oct = getX(XMG2OCT); \
    xfe2oct = getX(XFE2OCT); \
    xal3oct = getX(XAL3OCT); \
    xfe3oct = getX(XFE3OCT); \
    xcr3oct = getX(XCR3OCT); \
    xti4oct = getX(XTI4OCT);


#define SET_SITE_FRACTIONS \
    setX(XMG2TET, xmg2tet); \
    setX(XFE2TET, xfe2tet); \
    setX(XAL3TET, xal3tet); \
    setX(XFE3TET, xfe3tet); \
    setX(XMG2OCT, xmg2oct); \
    setX(XFE2OCT, xfe2oct); \
    setX(XAL3OCT, xal3oct); \
    setX(XFE3OCT, xfe3oct); \
    setX(XCR3OCT, xcr3oct); \
    setX(XTI4OCT, xti4oct);

/*
 *=============================================================================
 * Local function to compute ordering state and associated derivatives of
 * pure component endmembers
 */

#define HC_S  -R*( s[1]*log(s[1]) + 2.0*(1.0-s[1])*log(1.0-s[1]) \
                                    + (1.0+s[1])*log(1.0+s[1]) - 2.0*log(2.0) )
#define HC_H  g0 + gs2*s[1] + gs2s2*s[1]*s[1]
#define HC_G  HC_H - t*(HC_S)

#define DHC_GDS1 gs2 + 2.0*gs2s2*s[1] \
                 + R*t*(log(s[1]) + log(1.0+s[1]) - 2.0*log(1.0-s[1]))
#define DHC_GDT  -(HC_S)
#define DHC_GDP  0.0

#define D2HC_GDS1S1 2.0*gs2s2 + R*t*(1.0/s[1] + 1.0/(1.0+s[1]) + 2.0/(1.0-s[1]))
#define D2HC_GDS1DT R*(log(s[1]) + log(1.0+s[1]) - 2.0*log(1.0-s[1]))
#define D2HC_GDS1DP 0.0
#define D2HC_GDT2   0.0
#define D2HC_GDTDP  0.0
#define D2HC_GDP2   0.0

#define D3HC_GDS1S1S1 - R*t*(1.0/SQUARE(s[1]) + 1.0/SQUARE(1.0+s[1]) \
                             - 2.0/SQUARE(1.0-s[1]))
#define D3HC_GDS1S1DT R*(1.0/s[1] + 1.0/(1.0+s[1]) + 2.0/(1.0-s[1]))
#define D3HC_GDS1S1DP 0.0
#define D3HC_GDS1DT2  0.0
#define D3HC_GDS1DTDP 0.0
#define D3HC_GDS1DP2  0.0
#define D3HC_GDT3     0.0
#define D3HC_GDT2DP   0.0
#define D3HC_GDTDP2   0.0
#define D3HC_GDP3     0.0

#define SP_S  -R*(  0.5*(1.0+s[0])*log(1.0+s[0]) + (1.0-s[0])*log(1.0-s[0]) \
                                    + 0.5*(3.0+s[0])*log(3.0+s[0]) - 5.0*log(2.0) )
#define SP_H  g0 + gx2 + gs1*s[0] + gs2*(1.0+s[0])/2.0 + gx2x2 \
                 + gx2s1*s[0] + gx2s2*(1.0+s[0])/2.0 + gs1s1*s[0]*s[0] \
                 + gs1s2*s[0]*(1.0+s[0])/2.0 + gs2s2*SQUARE(1.0+s[0])/4.0
#define SP_G  SP_H - t*(SP_S)

#define DSP_GDS0 gs1 + gs2/2.0 + gx2s1 + gx2s2/2.0 + 2.0*gs1s1*s[0] \
                 + gs1s2*(0.5+s[0]) + gs2s2*(1.0+s[0])/2.0 \
                 + R*t*(0.5*log(1.0+s[0]) + 0.5*log(3.0+s[0]) - log(1.0-s[0]))
#define DSP_GDT  -(SP_S)
#define DSP_GDP  0.0

#define D2SP_GDS0S0 2.0*gs1s1 + gs1s2 + gs2s2/2.0 + R*t*(0.5/(1.0 + s[0]) \
                                        + 0.5/(3.0 + s[0]) + 1.0/(1.0 - s[0]))
#define D2SP_GDS0DT R*(0.5*log(1.0+s[0]) + 0.5*log(3.0+s[0]) - log(1.0-s[0]))
#define D2SP_GDS0DP 0.0
#define D2SP_GDT2   0.0
#define D2SP_GDTDP  0.0
#define D2SP_GDP2   0.0

#define D3SP_GDS0S0S0 - R*t*(0.5/SQUARE(1.0 + s[0]) + 0.5/SQUARE(3.0 + s[0]) \
                             - 1.0/SQUARE(1.0 - s[0]))
#define D3SP_GDS0S0DT R*(0.5/(1.0+s[0]) + 0.5/(3.0+s[0]) + 1.0/(1.0-s[0]))
#define D3SP_GDS0S0DP 0.0
#define D3SP_GDS0DT2  0.0
#define D3SP_GDS0DTDP 0.0
#define D3SP_GDS0DP2  0.0
#define D3SP_GDT3     0.0
#define D3SP_GDT2DP   0.0
#define D3SP_GDTDP2   0.0
#define D3SP_GDP3     0.0

#define CR_S        0.0
#define CR_H        g0 + gx3 + gs3 + gx3x3 + gx3s3 + gs3s3
#define CR_G        CR_H - t*(CR_S)

#define DCR_GDT     -(CR_S)
#define DCR_GDP     0.0

#define D2CR_GDT2   0.0
#define D2CR_GDTDP  0.0
#define D2CR_GDP2   0.0

#define D3CR_GDT3   0.0
#define D3CR_GDT2DP 0.0
#define D3CR_GDTDP2 0.0
#define D3CR_GDP3   0.0

#define UV_S        R*2.0*log(2.0)
#define UV_H        g0 + gx4 + gx4x4
#define UV_G        UV_H - t*(UV_S)

#define DUV_GDT     -(UV_S)
#define DUV_GDP     0.0

#define D2UV_GDT2   0.0
#define D2UV_GDTDP  0.0
#define D2UV_GDP2   0.0

#define D3UV_GDT3   0.0
#define D3UV_GDT2DP 0.0
#define D3UV_GDTDP2 0.0
#define D3UV_GDP3   0.0

#define MT_S  -R*( s[2]*log(s[2]) + 2.0*(1.0-s[2])*log(1.0-s[2]) \
                   + (1.0+s[2])*log(1.0+s[2]) - 2.0*log(2.0) ) + ss4*s[2]
#define MT_H  g0 + gx5 + hs4*s[2] + gx5x5 + gx5s4*s[2] + gs4s4*s[2]*s[2]
#define MT_G  MT_H - t*(MT_S)

#define DMT_GDS2 hs4 - t*ss4 + gx5s4 + 2.0*gs4s4*s[2] \
                 + R*t*(log(s[2]) - 2.0*log(1.0-s[2]) + log(1.0+s[2]))
#define DMT_GDT  -(MT_S)
#define DMT_GDP  0.0

#define D2MT_GDS2S2 2.0*gs4s4 + R*t*(1.0/s[2] + 2.0/(1.0-s[2]) + 1.0/(1.0+s[2]))
#define D2MT_GDS2DT R*(log(s[2]) - 2.0*log(1.0-s[2]) + log(1.0+s[2])) - ss4
#define D2MT_GDS2DP 0.0
#define D2MT_GDT2   0.0
#define D2MT_GDTDP  0.0
#define D2MT_GDP2   0.0

#define D3MT_GDS2S2S2 - R*t*(1.0/SQUARE(s[2]) - 2.0/SQUARE(1.0-s[2]) \
                             + 1.0/SQUARE(1.0+s[2]))
#define D3MT_GDS2S2DT R*(1.0/s[2] + 2.0/(1.0-s[2]) + 1.0/(1.0+s[2]))
#define D3MT_GDS2S2DP 0.0
#define D3MT_GDS2DT2  0.0
#define D3MT_GDS2DTDP 0.0
#define D3MT_GDS2DP2  0.0
#define D3MT_GDT3     0.0
#define D3MT_GDT2DP   0.0
#define D3MT_GDTDP2   0.0
#define D3MT_GDP3     0.0

#define fillD2GDSDT \
 d2gdsdt[0] = D2SP_GDS0DT; d2gdsdt[1] = D2HC_GDS1DT; d2gdsdt[2] = D2MT_GDS2DT;

#define fillD2GDSDP \
 d2gdsdp[0] = D2SP_GDS0DP; d2gdsdp[1] = D2HC_GDS1DP; d2gdsdp[2] = D2MT_GDS2DP;

#define fillD3GDS3 \
 d3gds3[0] = D3SP_GDS0S0S0; d3gds3[1] = D3HC_GDS1S1S1; \
 d3gds3[2] = D3MT_GDS2S2S2;

#define fillD3GDS2DT \
 d3gds2dt[0] = D3SP_GDS0S0DT; d3gds2dt[1] = D3HC_GDS1S1DT; \
 d3gds2dt[2] = D3MT_GDS2S2DT;

#define fillD3GDS2DP \
 d3gds2dp[0] = D3SP_GDS0S0DP; d3gds2dp[1] = D3HC_GDS1S1DP; \
 d3gds2dp[2] = D3MT_GDS2S2DP;

#define fillD3GDSDT2 \
 d3gdsdt2[0] = D3SP_GDS0DT2; d3gdsdt2[1] = D3HC_GDS1DT2; \
 d3gdsdt2[2] = D3MT_GDS2DT2;

#define fillD3GDSDTDP \
 d3gdsdtdp[0] = D3SP_GDS0DTDP; d3gdsdtdp[1] = D3HC_GDS1DTDP; \
 d3gdsdtdp[2] = D3MT_GDS2DTDP;

#define fillD3GDSDP2 \
 d3gdsdp2[0] = D3SP_GDS0DP2; d3gdsdp2[1] = D3HC_GDS1DP2; \
 d3gdsdp2[2] = D3MT_GDS2DP2;

static void
pureOrder(int mask, double t, double p,
            double s[NS],   /* s[NS]       BINARY MASK: 000001 */
            double dt[NS],  /* ds[NS]/dt   BINARY MASK: 000010 */
            double dp[NS],  /* ds[NS]/dp   BINARY MASK: 000100 */
            double dt2[NS], /* d2s[NS]/dt2 BINARY MASK: 001000 */
            double dtp[NS], /* d2s[NS]/dtp BINARY MASK: 010000 */
            double dp2[NS]  /* d2s[NS]/dp2 BINARY MASK: 100000 */
            )
{
    double tOld    = getTOldPure();
    double pOld    = getPOldPure();
    double *sOld   = getSOldPure();
    double *d2gds2 = getD2gds2Pure();
    int i;

#ifndef CALIBRATE_SPINEL
    if ( (t != tOld) || (p != pOld) ) {
#else
    if ( (t != tOld) || (p != pOld) || noSave) {
#endif
        double dgds[NS], sNew[NS];
        for (i=0; i<NS; i++) sOld[i] = 2.0;
        sNew[0] = 0.5;
        sNew[1] = 0.9;
        sNew[2] = 0.1;
        while ( (ABS(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON) ||
                        (ABS(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON) ||
                        (ABS(sNew[2]-sOld[2]) > 10.0*DBL_EPSILON) ) {
            double s[NS];

            for (i=0; i<NS; i++) s[i] = sNew[i];

            dgds[0] = DSP_GDS0;
            dgds[1] = DHC_GDS1;
            dgds[2] = DMT_GDS2;

            d2gds2[0] = D2SP_GDS0S0;
            d2gds2[1] = D2HC_GDS1S1;
            d2gds2[2] = D2MT_GDS2S2;

            for (i=0; i<NS; i++) sOld[i] = s[i];

            for (i=0; i<NS; i++) {
         s[i] += - dgds[i]/d2gds2[i];
         s[i] = MIN(s[i], 1.0 - DBL_EPSILON);
            }
            s[0] = MAX(s[0], -1.0 + DBL_EPSILON);
            s[1] = MAX(s[1], DBL_EPSILON);
            s[2] = MAX(s[2], DBL_EPSILON);

            for (i=0; i<NS; i++) sNew[i] = s[i];
        }
        tOld = t;
        pOld = p;

        setTOldPure(tOld);
        setPOldPure(pOld);
        /* sOld and d2gds2 are automatically saved */
    }

    if (mask & FIRST  ) {   /* return s        */
        for (i=0; i<NS; i++) s[i] = sOld[i];
    }

    if (mask & SECOND ) {   /* compute ds/dt:  */
        double *s = sOld;
        double d2gdsdt[NS];

        fillD2GDSDT

        for (i=0; i<NS; i++) dt[i] = - d2gdsdt[i]/d2gds2[i];
    }

    if (mask & THIRD  ) {   /* compute ds/dp:  */
        double d2gdsdp[NS];

        fillD2GDSDP

        for (i=0; i<NS; i++) dp[i] = - d2gdsdp[i]/d2gds2[i];
    }

    if (mask & FOURTH ) {   /* compute d2s/dt2 */
        double *s = sOld;
        double d2gdsdt[NS], d3gds3[NS], d3gds2dt[NS], d3gdsdt2[NS], dsdt[NS];

        fillD2GDSDT
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDSDT2

        for (i=0; i<NS; i++) dsdt[i] = - d2gdsdt[i]/d2gds2[i];
        for (i=0; i<NS; i++) dt2[i] = - (d3gdsdt2[i] + 2.0*d3gds2dt[i]*dsdt[i]
            + d3gds3[i]*dsdt[i]*dsdt[i])/d2gds2[i];
    }

    if (mask & FIFTH  ) {   /* compute d2s/dtp */
        double *s = sOld;
        double d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS], d3gds2dt[NS], d3gds2dp[NS],
            d3gdsdtdp[NS], dsdt[NS], dsdp[NS];

        fillD2GDSDT
        fillD2GDSDP
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDS2DP
        fillD3GDSDTDP

        for (i=0; i<NS; i++) dsdt[i] = - d2gdsdt[i]/d2gds2[i];
        for (i=0; i<NS; i++) dsdp[i] = - d2gdsdp[i]/d2gds2[i];

        for (i=0; i<NS; i++)
            dtp[i] = - (d3gdsdtdp[i] + d3gds2dt[i]*dsdp[i] + d3gds2dp[i]*dsdt[i]
             + d3gds3[i]*dsdt[i]*dsdp[i])/d2gds2[i];

    }

    if (mask & SIXTH  ) {   /* compute d2s/dp2 */
        double *s = sOld;
        double d2gdsdp[NS], d3gds3[NS], d3gds2dp[NS], d3gdsdp2[NS], dsdp[NS];

        fillD2GDSDP
        fillD3GDS3
        fillD3GDS2DP
        fillD3GDSDP2

        for (i=0; i<NS; i++) dsdp[i] = - d2gdsdp[i]/d2gds2[i];

        for (i=0; i<NS; i++) dp2[i] = - (d3gdsdp2[i] + 2.0*d3gds2dp[i]*dsdp[i]
            + d3gds3[i]*dsdp[i]*dsdp[i])/d2gds2[i];
    }
}

#undef fillD2GDSDT
#undef fillD2GDSDP
#undef fillD3GDS3
#undef fillD3GDS2DT
#undef fillD3GDS2DP
#undef fillD3GDSDT2
#undef fillD3GDSDTDP
#undef fillD3GDSDP2

#define fillD2GDS2 \
 d2gds2[0][0] = 0.0;         d2gds2[0][1] = 0.0; \
 d2gds2[0][2] = 0.0; \
 d2gds2[1][0] = 0.0;         d2gds2[1][1] = D2HC_GDS1S1; \
 d2gds2[1][2] = 0.0; \
 d2gds2[2][0] = 0.0;         d2gds2[2][1] = 0.0; \
 d2gds2[2][2] = D2MT_GDS2S2; \
 d2gds2[3][0] = D2SP_GDS0S0; d2gds2[3][1] = 0.0; \
 d2gds2[3][2] = 0.0; \
 d2gds2[4][0] = 0.0;         d2gds2[4][1] = 0.0; \
 d2gds2[4][2] = 0.0;

#define fillD2GDSDT \
 d2gdsdt[0][0] = 0.0;         d2gdsdt[0][1] = 0.0; \
 d2gdsdt[0][2] = 0.0; \
 d2gdsdt[1][0] = 0.0;         d2gdsdt[1][1] = D2HC_GDS1DT; \
 d2gdsdt[1][2] = 0.0; \
 d2gdsdt[2][0] = 0.0;         d2gdsdt[2][1] = 0.0; \
 d2gdsdt[2][2] = D2MT_GDS2DT; \
 d2gdsdt[3][0] = D2SP_GDS0DT; d2gdsdt[3][1] = 0.0; \
 d2gdsdt[3][2] = 0.0; \
 d2gdsdt[4][0] = 0.0;         d2gdsdt[4][1] = 0.0; \
 d2gdsdt[4][2] = 0.0;

#define fillD2GDSDP \
 d2gdsdp[0][0] = 0.0;         d2gdsdp[0][1] = 0.0; \
 d2gdsdp[0][2] = 0.0; \
 d2gdsdp[1][0] = 0.0;         d2gdsdp[1][1] = D2HC_GDS1DP; \
 d2gdsdp[1][2] = 0.0; \
 d2gdsdp[2][0] = 0.0;         d2gdsdp[2][1] = 0.0; \
 d2gdsdp[2][2] = D2MT_GDS2DP; \
 d2gdsdp[3][0] = D2SP_GDS0DP; d2gdsdp[3][1] = 0.0; \
 d2gdsdp[3][2] = 0.0; \
 d2gdsdp[4][0] = 0.0;         d2gdsdp[4][1] = 0.0; \
 d2gdsdp[4][2] = 0.0;

#define fillD2GDT2 \
 d2gdt2[0] = D2CR_GDT2; d2gdt2[1] = D2HC_GDT2; d2gdt2[2] = D2MT_GDT2; \
 d2gdt2[3] = D2SP_GDT2; d2gdt2[4] = D2UV_GDT2;

#define fillD2GDTDP \
 d2gdtdp[0] = D2CR_GDTDP; d2gdtdp[1] = D2HC_GDTDP; d2gdtdp[2] = D2MT_GDTDP; \
 d2gdtdp[3] = D2SP_GDTDP; d2gdtdp[4] = D2UV_GDTDP;

#define fillD2GDP2 \
 d2gdp2[0] = D2CR_GDP2; d2gdp2[1] = D2HC_GDP2; d2gdp2[2] = D2MT_GDP2; \
 d2gdp2[3] = D2SP_GDP2; d2gdp2[4] = D2UV_GDP2;

#define fillD3GDS3 \
 d3gds3[0][0] = 0.0;           d3gds3[0][1] = 0.0; \
 d3gds3[0][2] = 0.0; \
 d3gds3[1][0] = 0.0;           d3gds3[1][1] = D3HC_GDS1S1S1; \
 d3gds3[1][2] = 0.0; \
 d3gds3[2][0] = 0.0;           d3gds3[2][1] = 0.0; \
 d3gds3[2][2] = D3MT_GDS2S2S2; \
 d3gds3[3][0] = D3SP_GDS0S0S0; d3gds3[3][1] = 0.0; \
 d3gds3[3][2] = 0.0; \
 d3gds3[4][0] = 0.0;           d3gds3[4][1] = 0.0; \
 d3gds3[4][2] = 0.0;

#define fillD3GDS2DT \
 d3gds2dt[0][0] = 0.0;           d3gds2dt[0][1] = 0.0; \
 d3gds2dt[0][2] = 0.0; \
 d3gds2dt[1][0] = 0.0;           d3gds2dt[1][1] = D3HC_GDS1S1DT; \
 d3gds2dt[1][2] = 0.0; \
 d3gds2dt[2][0] = 0.0;           d3gds2dt[2][1] = 0.0; \
 d3gds2dt[2][2] = D3MT_GDS2S2DT; \
 d3gds2dt[3][0] = D3SP_GDS0S0DT; d3gds2dt[3][1] = 0.0; \
 d3gds2dt[3][2] = 0.0; \
 d3gds2dt[4][0] = 0.0;           d3gds2dt[4][1] = 0.0; \
 d3gds2dt[4][2] = 0.0;

#define fillD3GDS2DP \
 d3gds2dp[0][0] = 0.0;           d3gds2dp[0][1] = 0.0; \
 d3gds2dp[0][2] = 0.0; \
 d3gds2dp[1][0] = 0.0;           d3gds2dp[1][1] = D3HC_GDS1S1DP; \
 d3gds2dp[1][2] = 0.0; \
 d3gds2dp[2][0] = 0.0;           d3gds2dp[2][1] = 0.0; \
 d3gds2dp[2][2] = D3MT_GDS2S2DP; \
 d3gds2dp[3][0] = D3SP_GDS0S0DP; d3gds2dp[3][1] = 0.0; \
 d3gds2dp[3][2] = 0.0; \
 d3gds2dp[4][0] = 0.0;           d3gds2dp[4][1] = 0.0; \
 d3gds2dp[4][2] = 0.0;

#define fillD3GDSDT2 \
 d3gdsdt2[0][0] = 0.0;          d3gdsdt2[0][1] = 0.0; \
 d3gdsdt2[0][2] = 0.0; \
 d3gdsdt2[1][0] = 0.0;          d3gdsdt2[1][1] = D3HC_GDS1DT2; \
 d3gdsdt2[1][2] = 0.0; \
 d3gdsdt2[2][0] = 0.0;          d3gdsdt2[2][1] = 0.0; \
 d3gdsdt2[2][2] = D3MT_GDS2DT2; \
 d3gdsdt2[3][0] = D3SP_GDS0DT2; d3gdsdt2[3][1] = 0.0; \
 d3gdsdt2[3][2] = 0.0; \
 d3gdsdt2[4][0] = 0.0;          d3gdsdt2[4][1] = 0.0; \
 d3gdsdt2[4][2] = 0.0;

#define fillD3GDSDTDP \
 d3gdsdtdp[0][0] = 0.0;           d3gdsdtdp[0][1] = 0.0; \
 d3gdsdtdp[0][2] = 0.0; \
 d3gdsdtdp[1][0] = 0.0;           d3gdsdtdp[1][1] = D3HC_GDS1DTDP; \
 d3gdsdtdp[1][2] = 0.0; \
 d3gdsdtdp[2][0] = 0.0;           d3gdsdtdp[2][1] = 0.0; \
 d3gdsdtdp[2][2] = D3MT_GDS2DTDP; \
 d3gdsdtdp[3][0] = D3SP_GDS0DTDP; d3gdsdtdp[3][1] = 0.0; \
 d3gdsdtdp[3][2] = 0.0; \
 d3gdsdtdp[4][0] = 0.0;           d3gdsdtdp[4][1] = 0.0; \
 d3gdsdtdp[4][2] = 0.0;

#define fillD3GDSDP2 \
 d3gdsdp2[0][0] = 0.0;          d3gdsdp2[0][1] = 0.0; \
 d3gdsdp2[0][2] = 0.0; \
 d3gdsdp2[1][0] = 0.0;          d3gdsdp2[1][1] = D3HC_GDS1DP2; \
 d3gdsdp2[1][2] = 0.0; \
 d3gdsdp2[2][0] = 0.0;          d3gdsdp2[2][1] = 0.0; \
 d3gdsdp2[2][2] = D3MT_GDS2DP2; \
 d3gdsdp2[3][0] = D3SP_GDS0DP2; d3gdsdp2[3][1] = 0.0; \
 d3gdsdp2[3][2] = 0.0; \
 d3gdsdp2[4][0] = 0.0;          d3gdsdp2[4][1] = 0.0; \
 d3gdsdp2[4][2] = 0.0;

#define fillD3GDT3 \
 d3gdt3[0] = D3CR_GDT3; d3gdt3[1] = D3HC_GDT3; d3gdt3[2] = D3MT_GDT3; \
 d3gdt3[3] = D3SP_GDT3; d3gdt3[4] = D3UV_GDT3;

#define fillD3GDT2DP \
 d3gdt2dp[0] = D3CR_GDT2DP; d3gdt2dp[1] = D3HC_GDT2DP; \
 d3gdt2dp[2] = D3MT_GDT2DP; d3gdt2dp[3] = D3SP_GDT2DP; \
 d3gdt2dp[4] = D3UV_GDT2DP;

#define fillD3GDTDP2 \
 d3gdtdp2[0] = D3CR_GDTDP2; d3gdtdp2[1] = D3HC_GDTDP2; \
 d3gdtdp2[2] = D3MT_GDTDP2; d3gdtdp2[3] = D3SP_GDTDP2; \
 d3gdtdp2[4] = D3UV_GDTDP2;

#define fillD3GDP3 \
 d3gdp3[0] = D3CR_GDP3; d3gdp3[1] = D3HC_GDP3; d3gdp3[2] = D3MT_GDP3; \
 d3gdp3[3] = D3SP_GDP3; d3gdp3[4] = D3UV_GDP3;

static void
pureSpn(int mask, double t, double p,
    double a[NA],        /* activities              BINARY MASK: 0000000000001 */
    double mu[NA],       /* chemical potentials     BINARY MASK: 0000000000010 */
    double gmix[NA],     /* Gibbs energy            BINARY MASK: 0000000000100 */
    double hmix[NA],     /* Enthalpy of mixing      BINARY MASK: 0000000001000 */
    double smix[NA],     /* Entropy of mixing       BINARY MASK: 0000000010000 */
    double cpmix[NA],    /* Heat capacity of mixing BINARY MASK: 0000000100000 */
    double cpmixdt[NA],  /* d(cp)/d(t)              BINARY MASK: 0000001000000 */
    double vmix[NA],     /* Volume of mixing        BINARY MASK: 0000010000000 */
    double vmixdt[NA],   /* d(v)/d(t)               BINARY MASK: 0000100000000 */
    double vmixdp[NA],   /* d(v)/d(p)               BINARY MASK: 0001000000000 */
    double vmixdt2[NA],  /* d2(v)/d(t)2             BINARY MASK: 0010000000000 */
    double vmixdtdp[NA], /* d2(v)/d(t)d(p)          BINARY MASK: 0100000000000 */
    double vmixdp2[NA]   /* d2(v)/d(p)2             BINARY MASK: 1000000000000 */
    )
{
    double s[NS];
    int i, j;

    pureOrder(FIRST, t, p,
                        s,               (double *) NULL, (double *) NULL, (double *) NULL,
                        (double *) NULL, (double *) NULL);

    if (mask & FIRST) {
        a[0] = CR_G;
        a[0] = exp(a[0]/(R*t));
        a[1] = HC_G;
        a[1] = exp(a[1]/(R*t));
        a[2] = MT_G;
        a[2] = exp(a[2]/(R*t));
        a[3] = SP_G;
        a[3] = exp(a[3]/(R*t));
        a[4] = UV_G;
        a[4] = exp(a[4]/(R*t));
    }

    if (mask & SECOND) {
        mu[0] = CR_G;
        mu[1] = HC_G;
        mu[2] = MT_G;
        mu[3] = SP_G;
        mu[4] = UV_G;
    }

    if (mask & THIRD) {
        gmix[0] = CR_G;
        gmix[1] = HC_G;
        gmix[2] = MT_G;
        gmix[3] = SP_G;
        gmix[4] = UV_G;
    }

    if (mask & FOURTH) {
        hmix[0] = (CR_G) + t*(CR_S);
        hmix[1] = (HC_G) + t*(HC_S);
        hmix[2] = (MT_G) + t*(MT_S);
        hmix[3] = (SP_G) + t*(SP_S);
        hmix[4] = (UV_G) + t*(UV_S);
    }

    if (mask & FIFTH) {
        smix[0] = CR_S;
        smix[1] = HC_S;
        smix[2] = MT_S;
        smix[3] = SP_S;
        smix[4] = UV_S;
    }

    if (mask & SIXTH) {
        double d2gdsdt[NA][NS], d2gds2[NA][NS], dsdt[NS];

        fillD2GDS2
        fillD2GDSDT

        pureOrder(SECOND, t, p,
                        (double *) NULL, dsdt,            (double *) NULL, (double *) NULL,
                        (double *) NULL, (double *) NULL);

        cpmix[0] = D2CR_GDT2;
        cpmix[1] = D2HC_GDT2;
        cpmix[2] = D2MT_GDT2;
        cpmix[3] = D2SP_GDT2;
        cpmix[4] = D2UV_GDT2;

        for (i=0; i<NA; i++) {
            for (j=0; j<NS; j++)
                cpmix[i] += 2.0*d2gdsdt[i][j]*dsdt[j] + d2gds2[i][j]*SQUARE(dsdt[j]);
            cpmix[i] *= -t;
        }
    }

    if(mask & SEVENTH) {
        double d2gdsdt[NA][NS], d2gds2[NA][NS], d2gdt2[NA], d3gds3[NA][NS],
            d3gds2dt[NA][NS], d3gdsdt2[NA][NS], d3gdt3[NA], dsdt[NS], d2sdt2[NS],
            temp;

        fillD2GDT2
        fillD2GDSDT
        fillD2GDS2
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDSDT2
        fillD3GDT3

        pureOrder(SECOND | FOURTH, t, p,
                        (double *) NULL, dsdt,            (double *) NULL, d2sdt2,
                        (double *) NULL, (double *) NULL);

        for (i=0; i<NA; i++) {
            temp = d2gdt2[i];
            for (j=0; j<NS; j++)
                temp += 2.0*d2gdsdt[i][j]*dsdt[j] + d2gds2[i][j]*SQUARE(dsdt[j]);

            cpmixdt[i] = d3gdt3[i];
            for (j=0; j<NS; j++)
                cpmixdt[i] += 3.0*d3gdsdt2[i][j]*dsdt[j]
                                        + 3.0*d2gdsdt[i][j]*d2sdt2[j]
                                        + 3.0*d2gds2[i][j]*dsdt[j]*d2sdt2[j]
                                        + 3.0*d3gds2dt[i][j]*dsdt[j]*dsdt[j]
                                        + d3gds3[i][j]*dsdt[j]*dsdt[j]*dsdt[j];

            cpmixdt[i] = -t*cpmixdt[i] - temp;
        }
    }

    if (mask & EIGHTH) {
        vmix[0] = DCR_GDP;
        vmix[1] = DHC_GDP;
        vmix[2] = DMT_GDP;
        vmix[3] = DSP_GDP;
        vmix[4] = DUV_GDP;
    }

    if(mask & NINTH) {
        double d2gdsdt[NA][NS], d2gdsdp[NA][NS], d2gds2[NA][NS], d2gdtdp[NA],
            dsdt[NS], dsdp[NS];

        fillD2GDSDT
        fillD2GDS2
        fillD2GDSDP
        fillD2GDTDP

        pureOrder(SECOND | THIRD, t, p,
                        (double *) NULL, dsdt,            dsdp,            (double *) NULL,
                        (double *) NULL, (double *) NULL);

        for (i=0; i<NA; i++) {
            vmixdt[i] = d2gdtdp[i];
            for (j=0; j<NS; j++)
                vmixdt[i] += d2gdsdt[i][j]*dsdp[j] + d2gdsdp[i][j]*dsdt[j]
                                + d2gds2[i][j]*dsdt[j]*dsdp[j];
        }
    }

    if(mask & TENTH) {
        double d2gdsdp[NA][NS], d2gds2[NA][NS], d2gdp2[NA], dsdp[NS];

        fillD2GDS2
        fillD2GDSDP
        fillD2GDP2

        pureOrder(THIRD, t, p,
                        (double *) NULL, (double *) NULL, dsdp,            (double *) NULL,
                        (double *) NULL, (double *) NULL);

        for (i=0; i<NA; i++) {
            vmixdp[i] = d2gdp2[i];
            for (j=0; j<NS; j++)
                vmixdp[i] += 2.0*d2gdsdp[i][j]*dsdp[j] + d2gds2[i][j]*dsdp[j]*dsdp[j];
        }
    }

    if(mask & ELEVENTH) {
        double d2gdsdt[NA][NS], d2gdsdp[NA][NS], d2gds2[NA][NS], d3gds3[NA][NS],
            d3gds2dt[NA][NS], d3gdsdt2[NA][NS], d3gds2dp[NA][NS], d3gdsdtdp[NA][NS],
            d3gdt2dp[NA], dsdt[NS], dsdp[NS], d2sdt2[NS], d2sdtdp[NS];

        fillD2GDSDT
        fillD2GDS2
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDSDT2
        fillD2GDSDP
        fillD3GDS2DP
        fillD3GDSDTDP
        fillD3GDT2DP

        pureOrder(SECOND | THIRD | FOURTH | FIFTH, t, p,
                        (double *) NULL, dsdt,            dsdp,            d2sdt2,
                        d2sdtdp,         (double *) NULL);

        for (i=0; i<NA; i++) {
            vmixdt2[i] = d3gdt2dp[i];
            for (j=0; j<NS; j++)
                vmixdt2[i] += d3gdsdt2[i][j]*dsdp[j]
                                        + 2.0*d2gdsdt[i][j]*d2sdtdp[j]
                                        + d2gdsdp[i][j]*d2sdt2[j] + 2.0*d3gdsdtdp[i][j]*dsdt[j]
                                        + 2.0*d3gds2dt[i][j]*dsdt[j]*dsdp[j]
                                        + d2gds2[i][j]*d2sdt2[j]*dsdp[j]
                                        + 2.0*d2gds2[i][j]*dsdt[j]*d2sdtdp[j]
                                        + d3gds2dp[i][j]*dsdt[j]*dsdt[j]
                                        + d3gds3[i][j]*dsdt[j]*dsdt[j]*dsdp[j];
        }
    }

    if(mask & TWELFTH) {
        double d2gdsdt[NA][NS], d2gdsdp[NA][NS], d2gds2[NA][NS], d3gds3[NA][NS],
            d3gds2dt[NA][NS], d3gds2dp[NA][NS], d3gdsdtdp[NA][NS], d3gdsdp2[NA][NS],
            d3gdtdp2[NA], dsdt[NS], dsdp[NS], d2sdtdp[NS], d2sdp2[NS];

        fillD2GDSDT
        fillD2GDS2
        fillD3GDS3
        fillD3GDS2DT
        fillD2GDSDP
        fillD3GDS2DP
        fillD3GDSDTDP
        fillD3GDSDP2
        fillD3GDTDP2

        pureOrder(SECOND | THIRD | FIFTH | SIXTH, t, p,
                        (double *) NULL, dsdt,            dsdp,            (double *) NULL,
                        d2sdtdp,         d2sdp2);

        for (i=0; i<NA; i++) {
            vmixdtdp[i] = d3gdtdp2[i];
            for (j=0; j<NS; j++)
                vmixdtdp[i] += 2.0*d3gdsdtdp[i][j]*dsdp[j] + d2gdsdt[i][j]*d2sdp2[j]
                     + 2.0*d2gdsdp[i][j]*d2sdtdp[j] + d3gdsdp2[i][j]*dsdt[j]
                     + 2.0*d3gds2dp[i][j]*dsdt[j]*dsdp[j]
                     + d2gds2[i][j]*dsdt[j]*d2sdp2[j]
                     + 2.0*d2gds2[i][j]*d2sdtdp[j]*dsdp[j]
                     + d3gds2dt[i][j]*dsdp[j]*dsdp[j]
                     + d3gds3[i][j]*dsdt[j]*dsdp[j]*dsdp[j];
        }
    }

    if(mask & THIRTEENTH) {
        double d2gdsdp[NA][NS], d2gds2[NA][NS], d3gds3[NA][NS], d3gds2dp[NA][NS],
            d3gdsdp2[NA][NS], d3gdp3[NA], dsdp[NS], d2sdp2[NS];

        fillD2GDS2
        fillD3GDS3
        fillD2GDSDP
        fillD3GDS2DP
        fillD3GDSDP2
        fillD3GDP3

        pureOrder(THIRD | SIXTH, t, p,
                        (double *) NULL, (double *) NULL, dsdp,            (double *) NULL,
                        (double *) NULL, d2sdp2);

        for (i=0; i<NA; i++) {
            vmixdp2[i] = d3gdp3[i];
            for (j=0; j<NS; j++)
                vmixdp2[i] += 3.0*d3gdsdp2[i][j]*dsdp[j] + 3.0*d2gdsdp[i][j]*d2sdp2[j]
                                        + 3.0*d2gds2[i][j]*dsdp[j]*d2sdp2[j]
                                        + 3.0*d3gds2dp[i][j]*dsdp[j]*dsdp[j]
                                        + d3gds3[i][j]*dsdp[j]*dsdp[j]*dsdp[j];
        }
    }

}

#undef fillD2GDS2
#undef fillD2GDSDT
#undef fillD2GDSDP
#undef fillD2GDT2
#undef fillD2GDTDP
#undef fillD2GDP2
#undef fillD3GDS3
#undef fillD3GDS2DT
#undef fillD3GDS2DP
#undef fillD3GDSDT2
#undef fillD3GDSDTDP
#undef fillD3GDSDP2
#undef fillD3GDT3
#undef fillD3GDT2DP
#undef fillD3GDTDP2
#undef fillD3GDP3

#undef HC_S
#undef HC_H
#undef HC_G
#undef DHC_GDS1
#undef DHC_GDT
#undef DHC_GDP
#undef D2HC_GDS1S1
#undef D2HC_GDS1DT
#undef D2HC_GDS1DP
#undef D2HC_GDT2
#undef D2HC_GDTDP
#undef D2HC_GDP2
#undef D3HC_GDS1S1S1
#undef D3HC_GDS1S1DT
#undef D3HC_GDS1S1DP
#undef D3HC_GDS1DT2
#undef D3HC_GDS1DTDP
#undef D3HC_GDS1DP2
#undef D3HC_GDT3
#undef D3HC_GDT2DP
#undef D3HC_GDTDP2
#undef D3HC_GDP3

#undef SP_S
#undef SP_H
#undef SP_G
#undef DSP_GDS0
#undef DSP_GDT
#undef DSP_GDP
#undef D2SP_GDS0S0
#undef D2SP_GDS0DT
#undef D2SP_GDS0DP
#undef D2SP_GDT2
#undef D2SP_GDTDP
#undef D2SP_GDP2
#undef D3SP_GDS0S0S0
#undef D3SP_GDS0S0DT
#undef D3SP_GDS0S0DP
#undef D3SP_GDS0DT2
#undef D3SP_GDS0DTDP
#undef D3SP_GDS0DP2
#undef D3SP_GDT3
#undef D3SP_GDT2DP
#undef D3SP_GDTDP2
#undef D3SP_GDP3

#undef CR_S
#undef CR_H
#undef CR_G
#undef DCR_GDT
#undef DCR_GDP
#undef D2CR_GDT2
#undef D2CR_GDTDP
#undef D2CR_GDP2
#undef D3CR_GDT3
#undef D3CR_GDT2DP
#undef D3CR_GDTDP2
#undef D3CR_GDP3

#undef UV_S
#undef UV_H
#undef UV_G
#undef DUV_GDT
#undef DUV_GDP
#undef D2UV_GDT2
#undef D2UV_GDTDP
#undef D2UV_GDP2
#undef D3UV_GDT3
#undef D3UV_GDT2DP
#undef D3UV_GDTDP2
#undef D3UV_GDP3

#undef MT_S
#undef MT_H
#undef MT_G
#undef DMT_GDS2
#undef DMT_GDT
#undef DMT_GDP
#undef D2MT_GDS2S2
#undef D2MT_GDS2DT
#undef D2MT_GDS2DP
#undef D2MT_GDT2
#undef D2MT_GDTDP
#undef D2MT_GDP2
#undef D3MT_GDS2S2S2
#undef D3MT_GDS2S2DT
#undef D3MT_GDS2S2DP
#undef D3MT_GDS2DT2
#undef D3MT_GDS2DTDP
#undef D3MT_GDS2DP2
#undef D3MT_GDT3
#undef D3MT_GDT2DP
#undef D3MT_GDTDP2
#undef D3MT_GDP3

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conSpn defines the conversion from m[i], to r[j]
 */
                   /* Order: X2, X3. X4, X5 */
#define FR2(i)     (i == 3) ? 1.0 - r[0] : - r[0]
#define FR3(i)     (i == 0) ? 1.0 - r[1] : - r[1]
#define FR4(i)     (i == 4) ? 1.0 - r[2] : - r[2]
#define FR5(i)     (i == 2) ? 1.0 - r[3] : - r[3]

                   /* Order: S1, S2, S3 = X3, S4 */
#define GS1(i)     (i == 3)           ? 1.0 - s[0] : - s[0]
#define GS2(i)     (i == 1 || i == 3) ? 1.0 - s[1] : - s[1]
#define GS3(i)     (i == 0)           ? 1.0 - r[1] : - r[1]
#define GS4(i)     - s[2]

#define DFR2DR2(i) - 1.0
#define DFR3DR3(i) - 1.0
#define DFR4DR4(i) - 1.0
#define DFR5DR5(i) - 1.0

#define DGS1DS1(i) - 1.0
#define DGS2DS2(i) - 1.0
#define DGS3DS3(i) - 1.0
#define DGS4DS4(i) - 1.0

#define ENDMEMBERS  (  r[1]*ends[0] + (1.0-r[0]-r[1]-r[2]-r[3])*ends[1] \
                     + r[3]*ends[2] + r[0]*ends[3] + r[2]*ends[4] )

#define DENDDR0  (ends[3] - ends[1])
#define DENDDR1  (ends[0] - ends[1])
#define DENDDR2  (ends[4] - ends[1])
#define DENDDR3  (ends[2] - ends[1])

/*
 * Global (to this file): derivative definitions
 */

#define S -R*(    xmg2tet*log(xmg2tet) +     xfe2tet*log(xfe2tet) + \
                                    xal3tet*log(xal3tet) +     xfe3tet*log(xfe3tet) + \
                            2.0*xmg2oct*log(xmg2oct) + 2.0*xfe2oct*log(xfe2oct) + \
                            2.0*xal3oct*log(xal3oct) + 2.0*xfe3oct*log(xfe3oct) + \
                            2.0*xcr3oct*log(xcr3oct) + 2.0*xti4oct*log(xti4oct) ) + ss4*s[2]
#define H     g0 + \
                            gx2*r[0] + gx3*r[1] + gx4*r[2] + gx5*r[3] + \
                            gs1*s[0] + gs2*s[1] + gs3*r[1] + hs4*s[2] + \
                            gx2x2*r[0]*r[0] + gx2x3*r[0]*r[1] + gx2x4*r[0]*r[2] + \
                            gx2x5*r[0]*r[3] + gx2s1*r[0]*s[0] + gx2s2*r[0]*s[1] + \
                            gx2s3*r[0]*r[1] + gx2s4*r[0]*s[2] + gx3x3*r[1]*r[1] + \
                            gx3x4*r[1]*r[2] + gx3x5*r[1]*r[3] + gx3s1*r[1]*s[0] + \
                            gx3s2*r[1]*s[1] + gx3s3*r[1]*r[1] + gx3s4*r[1]*s[2] + \
                            gx4x4*r[2]*r[2] + gx4x5*r[2]*r[3] + gx4s1*r[2]*s[0] + \
                            gx4s2*r[2]*s[1] + gx4s3*r[2]*r[1] + gx4s4*r[2]*s[2] + \
                            gx5x5*r[3]*r[3] + gx5s1*r[3]*s[0] + gx5s2*r[3]*s[1] + \
                            gx5s3*r[3]*r[1] + gx5s4*r[3]*s[2] + gs1s1*s[0]*s[0] + \
                            gs1s2*s[0]*s[1] + gs1s3*s[0]*r[1] + gs1s4*s[0]*s[2] + \
                            gs2s2*s[1]*s[1] + gs2s3*s[1]*r[1] + gs2s4*s[1]*s[2] + \
                            gs3s3*r[1]*r[1] + gs3s4*r[1]*s[2] + gs4s4*s[2]*s[2]
#define G     H - t*(S) + r[2]*r[3]*((WV1)*r[3]+(WV2)*r[2])*(p-1.0)

/*----------------------------------------------------------------------------*/

#define DGDR0 gx2 + 2.0*gx2x2*r[0] + gx2x3*r[1] + gx2x4*r[2] + gx2x5*r[3] + \
                            gx2s1*s[0] + gx2s2*s[1] + gx2s3*r[1] + gx2s4*s[2] + \
                            0.5*R*t*(log(xmg2tet/xfe2tet) + log(xmg2oct/xfe2oct) )
#define DGDR1 gx3 + gx2x3*r[0] + 2.0*gx3x3*r[1] + gx3x4*r[2] + gx3x5*r[3] + \
                            gx3s1*s[0] + gx3s2*s[1] + gx3s3*r[1] + gx3s4*s[2] + \
                            gs3 + gx2s3*r[0] + gx3s3*r[1] + gx4s3*r[2] + gx5s3*r[3] + \
                            gs1s3*s[0] + gs2s3*s[1] + 2.0*gs3s3*r[1] + gs3s4*s[2] + \
                            R*t*(log(xfe2tet/xal3tet) + \
                                    2.0*log(xcr3oct) - log(xfe2oct) - log(xal3oct) )
#define DGDR2 gx4 + gx2x4*r[0] + gx3x4*r[1] + 2.0*gx4x4*r[2] + gx4x5*r[3] + \
                            gx4s1*s[0] + gx4s2*s[1] + gx4s3*r[1] + gx4s4*s[2] + \
                            R*t*(log(xfe2tet/xal3tet) + log(xti4oct/xal3oct) ) + \
	      r[3]*((WV1)*r[3]+(WV2)*r[2])*(p-1.0) + r[2]*r[3]*(WV2)*(p-1.0)
#define DGDR3 gx5 + gx2x5*r[0] + gx3x5*r[1] + gx4x5*r[2] + 2.0*gx5x5*r[3] + \
                            gx5s1*s[0] + gx5s2*s[1] + gx5s3*r[1] + gx5s4*s[2] + \
                            R*t*(log(xfe3tet/xal3tet) + log(xfe3oct/xal3oct) ) + \
	      r[2]*((WV1)*r[3]+(WV2)*r[2])*(p-1.0) + r[2]*r[3]*(WV1)*(p-1.0)
#define DGDS0 gs1 + gx2s1*r[0] + gx3s1*r[1] + gx4s1*r[2] + gx5s1*r[3] + \
                            2.0*gs1s1*s[0] + gs1s2*s[1] + gs1s3*r[1] + gs1s4*s[2] + \
                            0.5*R*t*(log(xmg2tet/xfe2tet) + log(xfe2oct/xmg2oct) )
#define DGDS1 gs2 + gx2s2*r[0] + gx3s2*r[1] + gx4s2*r[2] + gx5s2*r[3] + \
                            gs1s2*s[0] + 2.0*gs2s2*s[1] + gs2s3*r[1] + gs2s4*s[2] + \
                            R*t*(log(xfe2tet/xal3tet) + log(xal3oct/xfe2oct) )
#define DGDS2 hs4 - t*ss4 + gx2s4*r[0] + gx3s4*r[1] + gx4s4*r[2] + gx5s4*r[3] + \
                            gs1s4*s[0] + gs2s4*s[1] + gs3s4*r[1] + 2.0*gs4s4*s[2] + \
                            R*t*(log(xfe2tet/xfe3tet) + log(xfe3oct/xfe2oct) )
#define DGDT  - (S)
#define DGDP  r[2]*r[3]*((WV1)*r[3]+(WV2)*r[2])

/*----------------------------------------------------------------------------*/

#define D2GDR0R0 2.0*gx2x2 + 0.25*R*t*( \
                 1.0/xmg2tet + 1.0/xfe2tet + 0.5/xmg2oct + 0.5/xfe2oct )
#define D2GDR0R1 gx2x3 + gx2s3 + 0.5*R*t*(- 1.0/xfe2tet + 0.5/xfe2oct)
#define D2GDR0R2 gx2x4 + 0.5*R*t*(- 1.0/xfe2tet)
#define D2GDR0R3 gx2x5
#define D2GDR0S0 gx2s1+ 0.25*R*t*( \
                 1.0/xmg2tet + 1.0/xfe2tet - 0.5/xmg2oct - 0.5/xfe2oct )
#define D2GDR0S1 gx2s2 + 0.5*R*t*(- 1.0/xfe2tet + 0.5/xfe2oct)
#define D2GDR0S2 gx2s4 + 0.5*R*t*(- 1.0/xfe2tet + 0.5/xfe2oct)
#define D2GDR0DT 0.5*R*( \
                 log(xmg2tet) - log(xfe2tet) + log(xmg2oct) - log(xfe2oct) )
#define D2GDR0DP 0.0

#define D2GDR1R1 2.0*(gx3x3 + gx3s3 + gs3s3) + R*t*( 1.0/xfe2tet + \
                 1.0/xal3tet + 2.0/xcr3oct + 0.5/xfe2oct + 0.5/xal3oct )
#define D2GDR1R2 gx3x4 + gx4s3 + R*t*(1.0/xfe2tet + 1.0/xal3tet + 0.5/xal3oct)
#define D2GDR1R3 gx3x5 + gx5s3 + R*t*(1.0/xal3tet + 0.5/xal3oct)
#define D2GDR1S0 gx3s1 + gs1s3 + R*t*(-0.5/xfe2tet - 0.25/xfe2oct)
#define D2GDR1S1 gx3s2 + gs2s3 + R*t*( \
                 1.0/xfe2tet + 1.0/xal3tet  + 0.5/xfe2oct - 0.5/xal3oct )
#define D2GDR1S2 gx3s4 + gs3s4 + R*t*(1.0/xfe2tet + 0.5/xfe2oct)
#define D2GDR1DT R*(log(xfe2tet/xal3tet) + \
                 2.0*log(xcr3oct) - log(xfe2oct) + - log(xal3oct) )
#define D2GDR1DP 0.0

#define D2GDR2R2 2.0*gx4x4 + R*t*(1.0/xfe2tet + 1.0/xal3tet + 0.5/xti4oct + 0.5/xal3oct ) + \
                 2.0*r[3]*(WV2)*(p-1.0)
#define D2GDR2R3 gx4x5 + R*t*(1.0/xal3tet + 0.5/xal3oct) + \
                 ((WV1)*r[3]+(WV2)*r[2])*(p-1.0) + r[3]*(WV1)*(p-1.0) + r[2]*(WV2)*(p-1.0)
#define D2GDR2S0 gx4s1 + R*t*(-0.5/xfe2tet)
#define D2GDR2S1 gx4s2 + R*t*(1.0/xfe2tet + 1.0/xal3tet - 0.5/xal3oct)
#define D2GDR2S2 gx4s4 + R*t*(1.0/xfe2tet)
#define D2GDR2DT R*(log(xfe2tet/xal3tet) + log(xti4oct/xal3oct) )
#define D2GDR2DP r[3]*((WV1)*r[3]+(WV2)*r[2]) + r[2]*r[3]*(WV2)

#define D2GDR3R3 2.0*gx5x5 + R*t*(1.0/xfe3tet + 1.0/xal3tet + 0.5/xfe3oct + 0.5/xal3oct ) + \
                 2.0*r[2]*(WV1)*(p-1.0)
#define D2GDR3S0 gx5s1
#define D2GDR3S1 gx5s2 + R*t*(1.0/xal3tet - 0.5/xal3oct)
#define D2GDR3S2 gx5s4 + R*t*(-1.0/xfe3tet + 0.5/xfe3oct)
#define D2GDR3DT R*(log(xfe3tet/xal3tet) + log(xfe3oct/xal3oct) )
#define D2GDR3DP r[2]*((WV1)*r[3]+(WV2)*r[2]) + r[2]*r[3]*(WV1)

#define D2GDS0S0 2.0*gs1s1 + 0.25*R*t*( \
                 1.0/xmg2tet + 1.0/xfe2tet + 0.5/xfe2oct + 0.5/xmg2oct )
#define D2GDS0S1 gs1s2 + 0.5*R*t*(- 1.0/xfe2tet - 0.5/xfe2oct)
#define D2GDS0S2 gs1s4 + 0.5*R*t*(- 1.0/xfe2tet - 0.5/xfe2oct)
#define D2GDS0DT 0.5*R*(log(xmg2tet/xfe2tet) + log(xfe2oct/xmg2oct) )
#define D2GDS0DP 0.0

#define D2GDS1S1 2.0*gs2s2 + R*t*( \
                 1.0/xfe2tet + 1.0/xal3tet + 0.5/xal3oct + 0.5/xfe2oct)
#define D2GDS1S2 gs2s4 + R*t*(1.0/xfe2tet + 0.5/xfe2oct)
#define D2GDS1DT R*(log(xfe2tet/xal3tet) + log(xal3oct/xfe2oct) )
#define D2GDS1DP 0.0

#define D2GDS2S2 2.0*gs4s4 + R*t*( \
                 1.0/xfe2tet + 1.0/xfe3tet + 0.5/xfe3oct + 0.5/xfe2oct )
#define D2GDS2DT R*(log(xfe2tet/xfe3tet) + log(xfe3oct/xfe2oct) ) - ss4
#define D2GDS2DP 0.0

#define D2GDT2   0.0
#define D2GDTDP  0.0
#define D2GDP2   0.0

/*----------------------------------------------------------------------------*/

#define D3GDR0R0R0 - 0.125*R*t*(1.0/SQUARE(xmg2tet) - 1.0/SQUARE(xfe2tet) \
                   + 0.25/SQUARE(xmg2oct) - 0.25/SQUARE(xfe2oct) )
#define D3GDR0R0R1 - 0.25*R*t*(1.0/SQUARE(xfe2tet)-0.25/SQUARE(xfe2oct) )
#define D3GDR0R0R2 - 0.25*R*t*(1.0/SQUARE(xfe2tet) )
#define D3GDR0R0R3 0.0
#define D3GDR0R1R1 + 0.5*R*t*(1.0/SQUARE(xfe2tet)+0.25/SQUARE(xfe2oct) )
#define D3GDR0R1R2 + 0.5*R*t*(1.0/SQUARE(xfe2tet) )
#define D3GDR0R1R3 0.0
#define D3GDR0R2R2 + 0.5*R*t*(1.0/SQUARE(xfe2tet) )
#define D3GDR0R2R3 0.0
#define D3GDR0R3R3 0.0

#define D3GDR0R0S0 - 0.125*R*t*(1.0/SQUARE(xmg2tet) - 1.0/SQUARE(xfe2tet) \
                   - 0.25/SQUARE(xmg2oct) + 0.25/SQUARE(xfe2oct) )
#define D3GDR0R0S1 - 0.25*R*t*(1.0/SQUARE(xfe2tet) - 0.25/SQUARE(xfe2oct))
#define D3GDR0R0S2 - 0.25*R*t*(1.0/SQUARE(xfe2tet) - 0.25/SQUARE(xfe2oct))
#define D3GDR0R0DT 0.25*R*( \
                   1.0/xmg2tet + 1.0/xfe2tet + 0.5/xmg2oct + 0.5/xfe2oct )
#define D3GDR0R0DP 0.0

#define D3GDR0R1S0 - 0.25*R*t*(1.0/SQUARE(xfe2tet) + 0.25/SQUARE(xfe2oct))
#define D3GDR0R1S1 - 0.5*R*t*(- 1.0/SQUARE(xfe2tet) - 0.25/SQUARE(xfe2oct))
#define D3GDR0R1S2 - 0.5*R*t*(- 1.0/SQUARE(xfe2tet) - 0.25/SQUARE(xfe2oct))
#define D3GDR0R1DT 0.5*R*(- 1.0/xfe2tet + 0.5/xfe2oct)
#define D3GDR0R1DP 0.0

#define D3GDR0R2S0 - 0.25*R*t*(1.0/SQUARE(xfe2tet))
#define D3GDR0R2S1 - 0.5*R*t*(- 1.0/SQUARE(xfe2tet))
#define D3GDR0R2S2 - 0.5*R*t*(- 1.0/SQUARE(xfe2tet))
#define D3GDR0R2DT 0.5*R*(- 1.0/xfe2tet)
#define D3GDR0R2DP 0.0

#define D3GDR0R3S0 0.0
#define D3GDR0R3S1 0.0
#define D3GDR0R3S2 0.0
#define D3GDR0R3DT 0.0
#define D3GDR0R3DP 0.0

#define D3GDR0S0S0 - 0.125*R*t*( 1.0/SQUARE(xmg2tet) - 1.0/SQUARE(xfe2tet) \
                   + 0.25/SQUARE(xmg2oct) - 0.25/SQUARE(xfe2oct) )
#define D3GDR0S0S1 - 0.25*R*t*(1.0/SQUARE(xfe2tet) + 0.25/SQUARE(xfe2oct))
#define D3GDR0S0S2 - 0.25*R*t*(1.0/SQUARE(xfe2tet) + 0.25/SQUARE(xfe2oct))
#define D3GDR0S0DT 0.25*R*( \
                   1.0/xmg2tet + 1.0/xfe2tet - 0.5/xmg2oct - 0.5/xfe2oct )
#define D3GDR0S0DP 0.0

#define D3GDR0S1S1 - 0.5*R*t*(- 1.0/SQUARE(xfe2tet) - 0.25/SQUARE(xfe2oct))
#define D3GDR0S1S2 - 0.5*R*t*(- 1.0/SQUARE(xfe2tet) - 0.25/SQUARE(xfe2oct))
#define D3GDR0S1DT 0.5*R*(- 1.0/xfe2tet + 0.5/xfe2oct)
#define D3GDR0S1DP 0.0

#define D3GDR0S2S2 - 0.5*R*t*(- 1.0/SQUARE(xfe2tet) - 0.25/SQUARE(xfe2oct))
#define D3GDR0S2DT 0.5*R*(- 1.0/xfe2tet + 0.5/xfe2oct)
#define D3GDR0S2DP 0.0

#define D3GDR1R1R1 R*t*(-1.0/SQUARE(xfe2tet)+1.0/SQUARE(xal3tet) \
		   -2.0/SQUARE(xcr3oct)+0.25/SQUARE(xfe2oct) \
		   +0.25/SQUARE(xal3oct) )
#define D3GDR1R1R2 R*t*(-1.0/SQUARE(xfe2tet)+1.0/SQUARE(xal3tet) \
		   +0.25/SQUARE(xal3oct) )
#define D3GDR1R1R3 R*t*(1.0/SQUARE(xal3tet) + 0.25/SQUARE(xal3oct) )
#define D3GDR1R2R2 R*t*(-1.0/SQUARE(xfe2tet)+1.0/SQUARE(xal3tet) \
		   +0.25/SQUARE(xal3oct) )
#define D3GDR1R2R3 R*t*(1.0/SQUARE(xal3tet) + 0.25/SQUARE(xal3oct) )
#define D3GDR1R3R3 R*t*(1.0/SQUARE(xal3tet) + 0.25/SQUARE(xal3oct) )
#define D3GDR1R1S0 - 0.5*R*t*(- 1.0/SQUARE(xfe2tet) + 0.25/SQUARE(xfe2oct))
#define D3GDR1R1S1 - R*t*(1.0/SQUARE(xfe2tet) - 1.0/SQUARE(xal3tet) \
                   - 0.25/SQUARE(xfe2oct) + 0.25/SQUARE(xal3oct) )
#define D3GDR1R1S2 - R*t*(1.0/SQUARE(xfe2tet) - 0.25/SQUARE(xfe2oct))
#define D3GDR1R1DT R*(1.0/xfe2tet + \
                   1.0/xal3tet + 2.0/xcr3oct + 0.5/xfe2oct + 0.5/xal3oct )
#define D3GDR1R1DP 0.0

#define D3GDR1R2S0 - R*t*(- 0.5/SQUARE(xfe2tet))
#define D3GDR1R2S1 - R*t*(1.0/SQUARE(xfe2tet) - 1.0/SQUARE(xal3tet) \
                   + 0.25/SQUARE(xal3oct))
#define D3GDR1R2S2 - R*t*(1.0/SQUARE(xfe2tet))
#define D3GDR1R2DT R*(1.0/xfe2tet + 1.0/xal3tet + 0.5/xal3oct)
#define D3GDR1R2DP 0.0

#define D3GDR1R3S0 0.0
#define D3GDR1R3S1 - R*t*(- 1.0/SQUARE(xal3tet) + 0.25/SQUARE(xal3oct))
#define D3GDR1R3S2 0.0
#define D3GDR1R3DT R*(1.0/xal3tet + 0.5/xal3oct)
#define D3GDR1R3DP 0.0

#define D3GDR1S0S0 - 0.25*R*t*(1.0/SQUARE(xfe2tet) - 0.25/SQUARE(xfe2oct))
#define D3GDR1S0S1 - 0.5*R*t*(-1.0/SQUARE(xfe2tet) + 0.25/SQUARE(xfe2oct))
#define D3GDR1S0S2 - 0.5*R*t*(-1.0/SQUARE(xfe2tet) + 0.25/SQUARE(xfe2oct))
#define D3GDR1S0DT R*(-0.5/xfe2tet - 0.25/xfe2oct)
#define D3GDR1S0DP 0.0

#define D3GDR1S1S1 - R*t*(1.0/SQUARE(xfe2tet) - 1.0/SQUARE(xal3tet)  \
                   - 0.25/SQUARE(xfe2oct) - 0.25/SQUARE(xal3oct))
#define D3GDR1S1S2 - R*t*(1.0/SQUARE(xfe2tet) - 0.25/SQUARE(xfe2oct))
#define D3GDR1S1DT R*(1.0/xfe2tet + 1.0/xal3tet  + 0.5/xfe2oct - 0.5/xal3oct)
#define D3GDR1S1DP 0.0

#define D3GDR1S2S2 - R*t*(1.0/SQUARE(xfe2tet) - 0.25/SQUARE(xfe2oct))
#define D3GDR1S2DT R*(1.0/xfe2tet + 0.5/xfe2oct)
#define D3GDR1S2DP 0.0

#define D3GDR2R2R2 R*t*(-1.0/SQUARE(xfe2tet)+1.0/SQUARE(xal3tet) \
		   +0.25/SQUARE(xal3oct) -0.25/SQUARE(xti4oct) )
#define D3GDR2R2R3 R*t*(1.0/SQUARE(xal3tet) + 0.25/SQUARE(xal3oct) ) + 2.0*(WV2)*(p-1.0)
#define D3GDR2R3R3 R*t*(1.0/SQUARE(xal3tet) + 0.25/SQUARE(xal3oct) ) + 2.0*(WV1)*(p-1.0)
#define D3GDR2R2S0 - R*t*(- 0.5/SQUARE(xfe2tet))
#define D3GDR2R2S1 - R*t*(1.0/SQUARE(xfe2tet) - 1.0/SQUARE(xal3tet) \
                   + 0.25/SQUARE(xal3oct))
#define D3GDR2R2S2 - R*t*(1.0/SQUARE(xfe2tet))
#define D3GDR2R2DT R*(1.0/xfe2tet + 1.0/xal3tet + 0.5/xti4oct + 0.5/xal3oct)
#define D3GDR2R2DP 2.0*r[3]*(WV2)

#define D3GDR2R3S0 0.0
#define D3GDR2R3S1 - R*t*(- 1.0/SQUARE(xal3tet) + 0.25/SQUARE(xal3oct))
#define D3GDR2R3S2 0.0
#define D3GDR2R3DT R*(1.0/xal3tet + 0.5/xal3oct)
#define D3GDR2R3DP 2.0*r[3]*(WV1) + 2.0*r[2]*(WV2)

#define D3GDR2S0S0 - R*t*(0.25/SQUARE(xfe2tet))
#define D3GDR2S0S1 - R*t*(-0.5/SQUARE(xfe2tet))
#define D3GDR2S0S2 - R*t*(-0.5/SQUARE(xfe2tet))
#define D3GDR2S0DT R*(-0.5/xfe2tet)
#define D3GDR2S0DP 0.0

#define D3GDR2S1S1 - R*t*(1.0/SQUARE(xfe2tet) - 1.0/SQUARE(xal3tet) \
                   - 0.25/SQUARE(xal3oct))
#define D3GDR2S1S2 - R*t*(1.0/SQUARE(xfe2tet))
#define D3GDR2S1DT R*(1.0/xfe2tet + 1.0/xal3tet - 0.5/xal3oct)
#define D3GDR2S1DP 0.0

#define D3GDR2S2S2 - R*t*(1.0/SQUARE(xfe2tet))
#define D3GDR2S2DT R*(1.0/xfe2tet)
#define D3GDR2S2DP 0.0

#define D3GDR3R3R3 R*t*(-1.0/SQUARE(xfe3tet)+1.0/SQUARE(xal3tet) \
		   -0.25/SQUARE(xfe3oct) + 0.25/SQUARE(xal3oct) )
#define D3GDR3R3S0 0.0
#define D3GDR3R3S1 - R*t*(- 1.0/SQUARE(xal3tet) + 0.25/SQUARE(xal3oct))
#define D3GDR3R3S2 - R*t*(- 1.0/SQUARE(xfe3tet) + 0.25/SQUARE(xfe3oct))
#define D3GDR3R3DT R*(1.0/xfe3tet + 1.0/xal3tet + 0.5/xfe3oct + 0.5/xal3oct)
#define D3GDR3R3DP 2.0*r[2]*(WV1)

#define D3GDR3S0S0 0.0
#define D3GDR3S0S1 0.0
#define D3GDR3S0S2 0.0
#define D3GDR3S0DT 0.0
#define D3GDR3S0DP 0.0

#define D3GDR3S1S1 - R*t*(- 1.0/SQUARE(xal3tet) - 0.25/SQUARE(xal3oct))
#define D3GDR3S1S2 0.0
#define D3GDR3S1DT R*(1.0/xal3tet - 0.5/xal3oct)
#define D3GDR3S1DP 0.0

#define D3GDR3S2S2 - R*t*(1.0/SQUARE(xfe3tet) + 0.25/SQUARE(xfe3oct))
#define D3GDR3S2DT R*(-1.0/xfe3tet + 0.5/xfe3oct)
#define D3GDR3S2DP 0.0

#define D3GDS0S0S0 - 0.125*R*t*(1.0/SQUARE(xmg2tet) - 1.0/SQUARE(xfe2tet) \
                   + 0.25/SQUARE(xfe2oct) - 0.25/SQUARE(xmg2oct))
#define D3GDS0S0S1 - 0.25*R*t*(1.0/SQUARE(xfe2tet) - 0.25/SQUARE(xfe2oct))
#define D3GDS0S0S2 - 0.25*R*t*(1.0/SQUARE(xfe2tet) - 0.25/SQUARE(xfe2oct))
#define D3GDS0S0DT 0.25*R*(1.0/xmg2tet + 1.0/xfe2tet + 0.5/xfe2oct + \
                   0.5/xmg2oct )
#define D3GDS0S0DP 0.0
#define D3GDS0S1S1 - 0.5*R*t*(- 1.0/SQUARE(xfe2tet) + 0.25/SQUARE(xfe2oct))
#define D3GDS0S1S2 - 0.5*R*t*(- 1.0/SQUARE(xfe2tet) + 0.25/SQUARE(xfe2oct))
#define D3GDS0S1DT 0.5*R*(- 1.0/xfe2tet - 0.5/xfe2oct)
#define D3GDS0S1DP 0.0
#define D3GDS0S2S2 - 0.5*R*t*(- 1.0/SQUARE(xfe2tet) + 0.25/SQUARE(xfe2oct))
#define D3GDS0S2DT 0.5*R*(- 1.0/xfe2tet - 0.5/xfe2oct)
#define D3GDS0S2DP 0.0
#define D3GDS0DT2  0.0
#define D3GDS0DTDP 0.0
#define D3GDS0DP2  0.0

#define D3GDS1S1S1 - R*t*(1.0/SQUARE(xfe2tet) - 1.0/SQUARE(xal3tet) \
                   + 0.25/SQUARE(xal3oct) - 0.25/SQUARE(xfe2oct))
#define D3GDS1S1S2 - R*t*(1.0/SQUARE(xfe2tet) - 0.25/SQUARE(xfe2oct))
#define D3GDS1S1DT R*(1.0/xfe2tet + 1.0/xal3tet + 0.5/xal3oct + 0.5/xfe2oct)
#define D3GDS1S1DP 0.0
#define D3GDS1S2S2 - R*t*(1.0/SQUARE(xfe2tet) - 0.25/SQUARE(xfe2oct))
#define D3GDS1S2DT R*(1.0/xfe2tet + 0.5/xfe2oct)
#define D3GDS1S2DP 0.0
#define D3GDS1DT2  0.0
#define D3GDS1DTDP 0.0
#define D3GDS1DP2  0.0

#define D3GDS2S2S2 - R*t*(1.0/SQUARE(xfe2tet) - 1.0/SQUARE(xfe3tet) \
                   + 0.25/SQUARE(xfe3oct) - 0.25/SQUARE(xfe2oct))
#define D3GDS2S2DT R*(1.0/xfe2tet + 1.0/xfe3tet + 0.5/xfe3oct + 0.5/xfe2oct)
#define D3GDS2S2DP 0.0
#define D3GDS2DT2  0.0
#define D3GDS2DTDP 0.0
#define D3GDS2DP2  0.0

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
#define D3GDR3DT2  0.0
#define D3GDR3DTDP 0.0
#define D3GDR3DP2  0.0


/*
 *=============================================================================
 * Macros for automatic array initialization
 */
#define fillD2GDR2 \
 d2gdr2[0][0] = D2GDR0R0;     d2gdr2[0][1] = D2GDR0R1;     \
 d2gdr2[0][2] = D2GDR0R2;     d2gdr2[0][3] = D2GDR0R3;     \
 d2gdr2[1][0] = d2gdr2[0][1]; d2gdr2[1][1] = D2GDR1R1;     \
 d2gdr2[1][2] = D2GDR1R2;     d2gdr2[1][3] = D2GDR1R3;     \
 d2gdr2[2][0] = d2gdr2[0][2]; d2gdr2[2][1] = d2gdr2[1][2]; \
 d2gdr2[2][2] = D2GDR2R2;     d2gdr2[2][3] = D2GDR2R3;     \
 d2gdr2[3][0] = d2gdr2[0][3]; d2gdr2[3][1] = d2gdr2[1][3]; \
 d2gdr2[3][2] = d2gdr2[2][3]; d2gdr2[3][3] = D2GDR3R3;

#define fillD2GDRDS \
 d2gdrds[0][0] = D2GDR0S0; d2gdrds[0][1] = D2GDR0S1; d2gdrds[0][2] = D2GDR0S2; \
 d2gdrds[1][0] = D2GDR1S0; d2gdrds[1][1] = D2GDR1S1; d2gdrds[1][2] = D2GDR1S2; \
 d2gdrds[2][0] = D2GDR2S0; d2gdrds[2][1] = D2GDR2S1; d2gdrds[2][2] = D2GDR2S2; \
 d2gdrds[3][0] = D2GDR3S0; d2gdrds[3][1] = D2GDR3S1; d2gdrds[3][2] = D2GDR3S2;

#define fillD2GDRDT \
 d2gdrdt[0] = D2GDR0DT; d2gdrdt[1] = D2GDR1DT; d2gdrdt[2] = D2GDR2DT; \
 d2gdrdt[3] = D2GDR3DT;

#define fillD2GDRDP \
 d2gdrdp[0] = D2GDR0DP; d2gdrdp[1] = D2GDR1DP; d2gdrdp[2] = D2GDR2DP; \
 d2gdrdp[3] = D2GDR3DP;

#define fillD2GDS2 \
 d2gds2[0][0] = D2GDS0S0;     d2gds2[0][1] = D2GDS0S1; \
 d2gds2[0][2] = D2GDS0S2; \
 d2gds2[1][0] = d2gds2[0][1]; d2gds2[1][1] = D2GDS1S1; \
 d2gds2[1][2] = D2GDS1S2; \
 d2gds2[2][0] = d2gds2[0][2]; d2gds2[2][1] = d2gds2[1][2]; \
 d2gds2[2][2] = D2GDS2S2;

#define fillD2GDSDT \
 d2gdsdt[0] = D2GDS0DT;  d2gdsdt[1] = D2GDS1DT; d2gdsdt[2] = D2GDS2DT;

#define fillD2GDSDP \
 d2gdsdp[0] = D2GDS0DP;  d2gdsdp[1] = D2GDS1DP; d2gdsdp[2] = D2GDS2DP;

#define fillD3GDR3 \
 d3gdr3[0][0][0] = D3GDR0R0R0;		d3gdr3[0][0][1] = D3GDR0R0R1; \
 d3gdr3[0][0][2] = D3GDR0R0R2;		d3gdr3[0][0][3] = D3GDR0R0R3; \
 d3gdr3[0][1][0] = d3gdr3[0][0][1];	d3gdr3[0][1][1] = D3GDR0R1R1; \
 d3gdr3[0][1][2] = D3GDR0R1R2;		d3gdr3[0][1][3] = D3GDR0R1R3; \
 d3gdr3[0][2][0] = d3gdr3[0][0][2];	d3gdr3[0][2][1] = d3gdr3[0][1][2]; \
 d3gdr3[0][2][2] = D3GDR0R2R2;		d3gdr3[0][2][3] = D3GDR0R2R3; \
 d3gdr3[0][3][0] = d3gdr3[0][0][3];	d3gdr3[0][3][1] = d3gdr3[0][1][3]; \
 d3gdr3[0][3][2] = d3gdr3[0][2][3];	d3gdr3[0][3][3] = D3GDR0R3R3; \
 d3gdr3[1][0][0] = d3gdr3[0][0][1];	d3gdr3[1][0][1] = d3gdr3[0][1][1]; \
 d3gdr3[1][0][2] = d3gdr3[0][1][2];	d3gdr3[1][0][3] = d3gdr3[0][1][3]; \
 d3gdr3[1][1][0] = d3gdr3[0][1][1];	d3gdr3[1][1][1] = D3GDR1R1R1; \
 d3gdr3[1][1][2] = D3GDR1R1R2;		d3gdr3[1][1][3] = D3GDR1R1R3; \
 d3gdr3[1][2][0] = d3gdr3[0][1][2];	d3gdr3[1][2][1] = d3gdr3[1][1][2]; \
 d3gdr3[1][2][2] = D3GDR1R2R2;		d3gdr3[1][2][3] = D3GDR1R2R3; \
 d3gdr3[1][3][0] = d3gdr3[0][1][3];	d3gdr3[1][3][1] = d3gdr3[1][1][3]; \
 d3gdr3[1][3][2] = d3gdr3[1][2][3];	d3gdr3[1][3][3] = D3GDR1R3R3; \
 d3gdr3[2][0][0] = d3gdr3[0][0][2];	d3gdr3[2][0][1] = d3gdr3[0][1][2]; \
 d3gdr3[2][0][2] = d3gdr3[0][2][2];	d3gdr3[2][0][3] = d3gdr3[0][2][3]; \
 d3gdr3[2][1][0] = d3gdr3[0][1][2];	d3gdr3[2][1][1] = d3gdr3[1][1][2]; \
 d3gdr3[2][1][2] = d3gdr3[1][2][2];	d3gdr3[2][1][3] = d3gdr3[1][2][3]; \
 d3gdr3[2][2][0] = d3gdr3[0][2][2];	d3gdr3[2][2][1] = d3gdr3[1][2][2]; \
 d3gdr3[2][2][2] = D3GDR2R2R2;		d3gdr3[2][2][3] = D3GDR2R2R3; \
 d3gdr3[2][3][0] = d3gdr3[0][2][3];	d3gdr3[2][3][1] = d3gdr3[1][2][3]; \
 d3gdr3[2][3][2] = d3gdr3[2][2][3];	d3gdr3[2][3][3] = D3GDR2R3R3; \
 d3gdr3[3][0][0] = d3gdr3[0][0][3];	d3gdr3[3][0][1] = d3gdr3[0][1][3]; \
 d3gdr3[3][0][2] = d3gdr3[0][2][3];	d3gdr3[3][0][3] = d3gdr3[0][3][3]; \
 d3gdr3[3][1][0] = d3gdr3[0][1][3];	d3gdr3[3][1][1] = d3gdr3[1][1][3]; \
 d3gdr3[3][1][2] = d3gdr3[1][2][3];	d3gdr3[3][1][3] = d3gdr3[1][3][3]; \
 d3gdr3[3][2][0] = d3gdr3[0][2][3];	d3gdr3[3][2][1] = d3gdr3[1][2][3]; \
 d3gdr3[3][2][2] = d3gdr3[2][2][3];	d3gdr3[3][2][3] = d3gdr3[2][3][3]; \
 d3gdr3[3][3][0] = d3gdr3[0][3][3];	d3gdr3[3][3][1] = d3gdr3[1][3][3]; \
 d3gdr3[3][3][2] = d3gdr3[2][3][3];	d3gdr3[3][3][3] = D3GDR3R3R3;

#define fillD3GDR2DS \
 d3gdr2ds[0][0][0] = D3GDR0R0S0;        d3gdr2ds[0][0][1] = D3GDR0R0S1; \
 d3gdr2ds[0][0][2] = D3GDR0R0S2; \
 d3gdr2ds[0][1][0] = D3GDR0R1S0;        d3gdr2ds[0][1][1] = D3GDR0R1S1; \
 d3gdr2ds[0][1][2] = D3GDR0R1S2; \
 d3gdr2ds[0][2][0] = D3GDR0R2S0;        d3gdr2ds[0][2][1] = D3GDR0R2S1; \
 d3gdr2ds[0][2][2] = D3GDR0R2S2; \
 d3gdr2ds[0][3][0] = D3GDR0R3S0;        d3gdr2ds[0][3][1] = D3GDR0R3S1; \
 d3gdr2ds[0][3][2] = D3GDR0R3S2; \
 d3gdr2ds[1][0][0] = d3gdr2ds[0][1][0]; d3gdr2ds[1][0][1] = d3gdr2ds[0][1][1]; \
 d3gdr2ds[1][0][2] = d3gdr2ds[0][1][2]; \
 d3gdr2ds[1][1][0] = D3GDR1R1S0;        d3gdr2ds[1][1][1] = D3GDR1R1S1; \
 d3gdr2ds[1][1][2] = D3GDR1R1S2; \
 d3gdr2ds[1][2][0] = D3GDR1R2S0;        d3gdr2ds[1][2][1] = D3GDR1R2S1; \
 d3gdr2ds[1][2][2] = D3GDR1R2S2; \
 d3gdr2ds[1][3][0] = D3GDR1R3S0;        d3gdr2ds[1][3][1] = D3GDR1R3S1; \
 d3gdr2ds[1][3][2] = D3GDR1R3S2; \
 d3gdr2ds[2][0][0] = d3gdr2ds[0][2][0]; d3gdr2ds[2][0][1] = d3gdr2ds[0][2][1]; \
 d3gdr2ds[2][0][2] = d3gdr2ds[0][2][2]; \
 d3gdr2ds[2][1][0] = d3gdr2ds[1][2][0]; d3gdr2ds[2][1][1] = d3gdr2ds[1][2][1]; \
 d3gdr2ds[2][1][2] = d3gdr2ds[1][2][2]; \
 d3gdr2ds[2][2][0] = D3GDR2R2S0;        d3gdr2ds[2][2][1] = D3GDR2R2S1; \
 d3gdr2ds[2][2][2] = D3GDR2R2S2; \
 d3gdr2ds[2][3][0] = D3GDR2R3S0;        d3gdr2ds[2][3][1] = D3GDR2R3S1; \
 d3gdr2ds[2][3][2] = D3GDR2R3S2; \
 d3gdr2ds[3][0][0] = d3gdr2ds[0][3][0]; d3gdr2ds[3][0][1] = d3gdr2ds[0][3][1]; \
 d3gdr2ds[3][0][2] = d3gdr2ds[0][3][2]; \
 d3gdr2ds[3][1][0] = d3gdr2ds[1][3][0]; d3gdr2ds[3][1][1] = d3gdr2ds[1][3][1]; \
 d3gdr2ds[3][1][2] = d3gdr2ds[1][3][2]; \
 d3gdr2ds[3][2][0] = d3gdr2ds[2][3][0]; d3gdr2ds[3][2][1] = d3gdr2ds[2][3][1]; \
 d3gdr2ds[3][2][2] = d3gdr2ds[2][3][2]; \
 d3gdr2ds[3][3][0] = D3GDR3R3S0;        d3gdr2ds[3][3][1] = D3GDR3R3S1; \
 d3gdr2ds[3][3][2] = D3GDR3R3S2;

#define fillD3GDR2DT \
 d3gdr2dt[0][0] = D3GDR0R0DT;     d3gdr2dt[0][1] = D3GDR0R1DT; \
 d3gdr2dt[0][2] = D3GDR0R2DT;     d3gdr2dt[0][3] = D3GDR0R3DT; \
 d3gdr2dt[1][0] = d3gdr2dt[0][1]; d3gdr2dt[1][1] = D3GDR1R1DT; \
 d3gdr2dt[1][2] = D3GDR1R2DT;     d3gdr2dt[1][3] = D3GDR1R3DT; \
 d3gdr2dt[2][0] = d3gdr2dt[0][2]; d3gdr2dt[2][1] = d3gdr2dt[1][2]; \
 d3gdr2dt[2][2] = D3GDR2R2DT;     d3gdr2dt[2][3] = D3GDR2R3DT; \
 d3gdr2dt[3][0] = d3gdr2dt[0][3]; d3gdr2dt[3][1] = d3gdr2dt[1][3]; \
 d3gdr2dt[3][2] = d3gdr2dt[2][3]; d3gdr2dt[3][3] = D3GDR3R3DT;

#define fillD3GDR2DP \
 d3gdr2dp[0][0] = D3GDR0R0DP;     d3gdr2dp[0][1] = D3GDR0R1DP; \
 d3gdr2dp[0][2] = D3GDR0R2DP;     d3gdr2dp[0][3] = D3GDR0R3DP; \
 d3gdr2dp[1][0] = d3gdr2dp[0][1]; d3gdr2dp[1][1] = D3GDR1R1DP; \
 d3gdr2dp[1][2] = D3GDR1R2DP;     d3gdr2dp[1][3] = D3GDR1R3DP; \
 d3gdr2dp[2][0] = d3gdr2dp[0][2]; d3gdr2dp[2][1] = d3gdr2dp[1][2]; \
 d3gdr2dp[2][2] = D3GDR2R2DP;     d3gdr2dp[2][3] = D3GDR2R3DP; \
 d3gdr2dp[3][0] = d3gdr2dp[0][3]; d3gdr2dp[3][1] = d3gdr2dp[1][3]; \
 d3gdr2dp[3][2] = d3gdr2dp[2][3]; d3gdr2dp[3][3] = D3GDR3R3DP;

#define fillD3GDRDS2 \
 d3gdrds2[0][0][0] = D3GDR0S0S0;        d3gdrds2[0][0][1] = D3GDR0S0S1; \
 d3gdrds2[0][0][2] = D3GDR0S0S2; \
 d3gdrds2[0][1][0] = d3gdrds2[0][0][1]; d3gdrds2[0][1][1] = D3GDR0S1S1; \
 d3gdrds2[0][1][2] = D3GDR0S1S2; \
 d3gdrds2[0][2][0] = d3gdrds2[0][0][2]; d3gdrds2[0][2][1] = d3gdrds2[0][1][2]; \
 d3gdrds2[0][2][2] = D3GDR0S2S2; \
 d3gdrds2[1][0][0] = D3GDR1S0S0;        d3gdrds2[1][0][1] = D3GDR1S0S1; \
 d3gdrds2[1][0][2] = D3GDR1S0S2; \
 d3gdrds2[1][1][0] = d3gdrds2[1][0][1]; d3gdrds2[1][1][1] = D3GDR1S1S1; \
 d3gdrds2[1][1][2] = D3GDR1S1S2; \
 d3gdrds2[1][2][0] = d3gdrds2[1][0][2]; d3gdrds2[1][2][1] = d3gdrds2[1][1][2]; \
 d3gdrds2[1][2][2] = D3GDR1S2S2; \
 d3gdrds2[2][0][0] = D3GDR2S0S0;        d3gdrds2[2][0][1] = D3GDR2S0S1; \
 d3gdrds2[2][0][2] = D3GDR2S0S2; \
 d3gdrds2[2][1][0] = d3gdrds2[2][0][1]; d3gdrds2[2][1][1] = D3GDR2S1S1; \
 d3gdrds2[2][1][2] = D3GDR2S1S2; \
 d3gdrds2[2][2][0] = d3gdrds2[2][0][2]; d3gdrds2[2][2][1] = d3gdrds2[2][1][2]; \
 d3gdrds2[2][2][2] = D3GDR2S2S2; \
 d3gdrds2[3][0][0] = D3GDR3S0S0;        d3gdrds2[3][0][1] = D3GDR3S0S1; \
 d3gdrds2[3][0][2] = D3GDR3S0S2; \
 d3gdrds2[3][1][0] = d3gdrds2[3][0][1]; d3gdrds2[3][1][1] = D3GDR3S1S1; \
 d3gdrds2[3][1][2] = D3GDR3S1S2; \
 d3gdrds2[3][2][0] = d3gdrds2[3][0][2]; d3gdrds2[3][2][1] = d3gdrds2[3][1][2]; \
 d3gdrds2[3][2][2] = D3GDR3S2S2;

#define fillD3GDRDSDT \
 d3gdrdsdt[0][0] = D3GDR0S0DT;  d3gdrdsdt[0][1] = D3GDR0S1DT; \
 d3gdrdsdt[0][2] = D3GDR0S2DT; \
 d3gdrdsdt[1][0] = D3GDR1S0DT;  d3gdrdsdt[1][1] = D3GDR1S1DT; \
 d3gdrdsdt[1][2] = D3GDR1S2DT; \
 d3gdrdsdt[2][0] = D3GDR2S0DT;  d3gdrdsdt[2][1] = D3GDR2S1DT; \
 d3gdrdsdt[2][2] = D3GDR2S2DT; \
 d3gdrdsdt[3][0] = D3GDR3S0DT;  d3gdrdsdt[3][1] = D3GDR3S1DT; \
 d3gdrdsdt[3][2] = D3GDR3S2DT;

#define fillD3GDRDSDP \
 d3gdrdsdp[0][0] = D3GDR0S0DP; d3gdrdsdp[0][1] = D3GDR0S1DP; \
 d3gdrdsdp[0][2] = D3GDR0S2DP; \
 d3gdrdsdp[1][0] = D3GDR1S0DP; d3gdrdsdp[1][1] = D3GDR1S1DP; \
 d3gdrdsdp[1][2] = D3GDR1S2DP; \
 d3gdrdsdp[2][0] = D3GDR2S0DP; d3gdrdsdp[2][1] = D3GDR2S1DP; \
 d3gdrdsdp[2][2] = D3GDR2S2DP; \
 d3gdrdsdp[3][0] = D3GDR3S0DP; d3gdrdsdp[3][1] = D3GDR3S1DP; \
 d3gdrdsdp[3][2] = D3GDR3S2DP;

#define fillD3GDS3 \
 d3gds3[0][0][0] = D3GDS0S0S0;      d3gds3[0][0][1] = D3GDS0S0S1; \
 d3gds3[0][0][2] = D3GDS0S0S2; \
 d3gds3[0][1][0] = d3gds3[0][0][1]; d3gds3[0][1][1] = D3GDS0S1S1; \
 d3gds3[0][1][2] = D3GDS0S1S2; \
 d3gds3[0][2][0] = d3gds3[0][0][2]; d3gds3[0][2][1] = d3gds3[0][1][2]; \
 d3gds3[0][2][2] = D3GDS0S2S2; \
 d3gds3[1][0][0] = d3gds3[0][0][1]; d3gds3[1][0][1] = d3gds3[0][1][1]; \
 d3gds3[1][0][2] = d3gds3[0][1][2]; \
 d3gds3[1][1][0] = d3gds3[0][1][1]; d3gds3[1][1][1] = D3GDS1S1S1; \
 d3gds3[1][1][2] = D3GDS1S1S2; \
 d3gds3[1][2][0] = d3gds3[0][1][2]; d3gds3[1][2][1] = d3gds3[1][1][2]; \
 d3gds3[1][2][2] = D3GDS1S2S2; \
 d3gds3[2][0][0] = d3gds3[0][0][2]; d3gds3[2][0][1] = d3gds3[0][1][2]; \
 d3gds3[2][0][2] = d3gds3[0][2][2]; \
 d3gds3[2][1][0] = d3gds3[0][1][2]; d3gds3[2][1][1] = d3gds3[1][1][2]; \
 d3gds3[2][1][2] = d3gds3[1][2][2]; \
 d3gds3[2][2][0] = d3gds3[0][2][2]; d3gds3[2][2][1] = d3gds3[1][2][2]; \
 d3gds3[2][2][2] = D3GDS2S2S2;

#define fillD3GDS2DT \
 d3gds2dt[0][0] = D3GDS0S0DT;     d3gds2dt[0][1] = D3GDS0S1DT; \
 d3gds2dt[0][2] = D3GDS0S2DT; \
 d3gds2dt[1][0] = d3gds2dt[0][1]; d3gds2dt[1][1] = D3GDS1S1DT; \
 d3gds2dt[1][2] = D3GDS1S2DT; \
 d3gds2dt[2][0] = d3gds2dt[0][2]; d3gds2dt[2][1] = d3gds2dt[1][2]; \
 d3gds2dt[2][2] = D3GDS2S2DT;

#define fillD3GDS2DP \
 d3gds2dp[0][0] = D3GDS0S0DP;     d3gds2dp[0][1] = D3GDS0S1DP; \
 d3gds2dp[0][2] = D3GDS0S2DP; \
 d3gds2dp[1][0] = d3gds2dp[0][1]; d3gds2dp[1][1] = D3GDS1S1DP; \
 d3gds2dp[1][2] = D3GDS1S2DP; \
 d3gds2dp[2][0] = d3gds2dp[0][2]; d3gds2dp[2][1] = d3gds2dp[1][2]; \
 d3gds2dp[2][2] = D3GDS2S2DP;

#define fillD3GDSDT2 \
 d3gdsdt2[0] = D3GDS0DT2; d3gdsdt2[1] = D3GDS1DT2; d3gdsdt2[2] = D3GDS2DT2;

#define fillD3GDSDTDP \
 d3gdsdtdp[0] = D3GDS0DTDP; d3gdsdtdp[1] = D3GDS1DTDP; \
 d3gdsdtdp[2] = D3GDS2DTDP;

#define fillD3GDSDP2 \
 d3gdsdp2[0] = D3GDS0DP2; d3gdsdp2[1] = D3GDS1DP2; d3gdsdp2[2] = D3GDS2DP2;

#define fillD3GDRDT2 \
 d3gdrdt2[0] = D3GDR0DT2; d3gdrdt2[1] = D3GDR1DT2; \
 d3gdrdt2[2] = D3GDR2DT2; d3gdrdt2[3] = D3GDR3DT2;

#define fillD3GDRDTDP \
 d3gdrdtdp[0] = D3GDR0DTDP; d3gdrdtdp[1] = D3GDR1DTDP; \
 d3gdrdtdp[2] = D3GDR2DTDP; d3gdrdtdp[3] = D3GDR3DTDP;

#define fillD3GDRDP2 \
 d3gdrdp2[0] = D3GDR0DP2; d3gdrdp2[1] = D3GDR1DP2; \
 d3gdrdp2[2] = D3GDR2DP2; d3gdrdp2[3] = D3GDR3DP2;

/*
 *=============================================================================
 * Local function to compute ordering state and associated derivatives
 */

static void
order(int mask, double t, double p, double r[NR],
            double s[NS],           /* s[NS]                BINARY MASK: 0000000001 */
            double dr[NS][NR] ,     /* ds[NS]/dr[NR]        BINARY MASK: 0000000010 */
            double dt[NS],          /* ds[NS]/dt            BINARY MASK: 0000000100 */
            double dp[NS],          /* ds[NS]/dp            BINARY MASK: 0000001000 */
            double dr2[NS][NR][NR], /* d2s[NS]/dr[NR]dr[NR] BINARY MASK: 0000010000 */
            double drt[NS][NR],     /* d2s[NS]/dr[NR]dt     BINARY MASK: 0000100000 */
            double drp[NS][NR],     /* d2s[NS]/dr[NR]dp     BINARY MASK: 0001000000 */
            double dt2[NS],         /* d2s[NS]/dt2          BINARY MASK: 0010000000 */
            double dtp[NS],         /* d2s[NS]/dtp          BINARY MASK: 0100000000 */
            double dp2[NS]          /* d2s[NS]/dp2          BINARY MASK: 1000000000 */
            )
{
    DECLARE_SITE_FRACTIONS
    double tOld         = getTOld();
    double pOld         = getPOld();
    double *rOld        = getROld();
    double *sOld        = getSOld();
    double **d2gds2     = getD2gds2();
    gsl_matrix      *ptToD2gds2  = getPtToD2gds2();
    gsl_permutation *indexD2gds2 = getIndexD2gds2();
    int i, j, iter = 0, signum;

    GET_SITE_FRACTIONS

    /* look-up or compute the current ordering state */
#ifndef CALIBRATE_SPINEL
    if ( (t != tOld)       || (p != pOld) ||
#else
    if ( (t != tOld)       || (p != pOld) || noSave ||
#endif
       (r[0] != rOld[0]) || (r[1] != rOld[1]) || (r[2] != rOld[2]) ||
       (r[3] != rOld[3]) ) {
        double dgds[NS], sNew[NS], ratio;
        double totAl, totCr, totFe2, totFe3, totMg, totTi;
        gsl_vector_view vvToDgds = gsl_vector_view_array(dgds, (size_t) NS);

        for (i=0; i<NS; i++) sOld[i] = 2.0;

        totAl  = 2.0*(1.0 - r[1] - r[2] - r[3]);
        totCr  = 2.0*r[1];
        totFe2 = 1.0 - r[0] + r[2];
        totFe3 = 2.0*r[3];
        totMg  = r[0];
        totTi  = r[2];
        ratio = 2.0 - totCr - totTi;      /* available oct / available tet sites */

        xmg2oct = totMg  * ratio/(1.0+ratio);
        xfe2oct = totFe2 * ratio/(1.0+ratio);
        xal3oct = totAl  * ratio/(1.0+ratio);
        xfe3oct = totFe3 * ratio/(1.0+ratio);

        xmg2tet = totMg  - xmg2oct;
        xfe2tet = totFe2 - xfe2oct;
        xal3tet = totAl  - xal3oct;
        xfe3tet = totFe3 - xfe3oct;

        xmg2oct /= 2.0;
        xfe2oct /= 2.0;
        xal3oct /= 2.0;
        xfe3oct /= 2.0;

        sNew[0] = xmg2tet - 2.0*xmg2oct;
        sNew[1] = xal3oct - xal3tet/2.0;
        sNew[2] = xfe3oct - xfe3tet/2.0;

        while ( ((ABS(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON) ||
             (ABS(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON) ||
             (ABS(sNew[2]-sOld[2]) > 10.0*DBL_EPSILON) ) && (iter < MAX_ITER)) {
            double s[NS], deltaS[NS], lambda;
            gsl_vector_view	vvToDeltaS = gsl_vector_view_array(deltaS, (size_t) NS);

            for (i=0; i<NS; i++) s[i] = sNew[i];

            xmg2tet = (r[0] + s[0])/2.0;
            xfe2tet = r[2] - 0.5*r[0] - 0.5*s[0] + s[1] + r[1] + s[2];
            xal3tet = 1.0 - r[1] - r[2] - r[3] - s[1];
            xfe3tet = r[3] - s[2];
            xmg2oct = (r[0] - s[0])/4.0;
            xfe2oct = (2.0 - r[0] + s[0] - 2.0*s[1] - 2.0*r[1] - 2.0*s[2])/4.0;
            xal3oct = (1.0 - r[1] - r[2] - r[3] + s[1])/2.0;
            xfe3oct = (r[3] + s[2])/2.0;
            xcr3oct = r[1];
            xti4oct = r[2]/2.0;

            if (xmg2tet <= 0.0) xmg2tet = DBL_EPSILON;
            if (xfe2tet <= 0.0) xfe2tet = DBL_EPSILON;
            if (xal3tet <= 0.0) xal3tet = DBL_EPSILON;
            if (xfe3tet <= 0.0) xfe3tet = DBL_EPSILON;
            if (xmg2oct <= 0.0) xmg2oct = DBL_EPSILON;
            if (xfe2oct <= 0.0) xfe2oct = DBL_EPSILON;
            if (xal3oct <= 0.0) xal3oct = DBL_EPSILON;
            if (xfe3oct <= 0.0) xfe3oct = DBL_EPSILON;
            if (xcr3oct <= 0.0) xcr3oct = DBL_EPSILON;
            if (xti4oct <= 0.0) xti4oct = DBL_EPSILON;

            if (xmg2tet >= 1.0) xmg2tet = 1.0 - DBL_EPSILON;
            if (xfe2tet >= 1.0) xfe2tet = 1.0 - DBL_EPSILON;
            if (xal3tet >= 1.0) xal3tet = 1.0 - DBL_EPSILON;
            if (xfe3tet >= 1.0) xfe3tet = 1.0 - DBL_EPSILON;
            if (xmg2oct >= 1.0) xmg2oct = 1.0 - DBL_EPSILON;
            if (xfe2oct >= 1.0) xfe2oct = 1.0 - DBL_EPSILON;
            if (xal3oct >= 1.0) xal3oct = 1.0 - DBL_EPSILON;
            if (xfe3oct >= 1.0) xfe3oct = 1.0 - DBL_EPSILON;
            if (xcr3oct >= 1.0) xcr3oct = 1.0 - DBL_EPSILON;
            if (xti4oct >= 1.0) xti4oct = 1.0 - DBL_EPSILON;

            dgds[0] = DGDS0;
            dgds[1] = DGDS1;
            dgds[2] = DGDS2;

            d2gds2[0][0] = D2GDS0S0;
            d2gds2[0][1] = D2GDS0S1;
            d2gds2[0][2] = D2GDS0S2;
            d2gds2[1][0] = d2gds2[0][1];
            d2gds2[1][1] = D2GDS1S1;
            d2gds2[1][2] = D2GDS1S2;
            d2gds2[2][0] = d2gds2[0][2];
            d2gds2[2][1] = d2gds2[1][2];
            d2gds2[2][2] = D2GDS2S2;

            for (i=0; i<NS; i++) sOld[i] = s[i];

            /* original: gaussj(ptToD2gds2, NS, (double **) NULL, 0);
            for (i=0; i<NS; i++) {
            for(j=0; j<NS; j++) s[i] += - d2gds2[i][j]*dgds[j];
            deltaS[i] = s[i] - sOld[i];
            } */

            gsl_matrix_scale(ptToD2gds2, -1.0);
            melts_LU_decomp(ptToD2gds2, indexD2gds2, &signum);

            melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToDgds.vector, &vvToDeltaS.vector);
            for (i=0; i<NS; i++) s[i] += deltaS[i];

            lambda = 1.0;

            /* xmg2tet = (r[0] + s[0])/2.0; */
            if      ((deltaS[0]/2.0) != 0.0 && xmg2tet+lambda*(deltaS[0]/2.0) < 0.0)
                lambda = -xmg2tet/(deltaS[0]/2.0);
            else if ((deltaS[0]/2.0) != 0.0 && xmg2tet+lambda*(deltaS[0]/2.0) > 1.0)
                lambda = (1.0-xmg2tet)/(deltaS[0]/2.0);

            /* xfe2tet = r[2] - 0.5*r[0] - 0.5*s[0] + s[1] + r[1] + s[2]; */
            if      ((-deltaS[0]/2.0+deltaS[1]+deltaS[2]) != 0.0 &&
                xfe2tet+lambda*(-deltaS[0]/2.0+deltaS[1]+deltaS[2]) < 0.0)
                lambda = -xfe2tet/(-deltaS[0]/2.0+deltaS[1]+deltaS[2]);
            else if ((-deltaS[0]/2.0+deltaS[1]+deltaS[2]) != 0.0 &&
                xfe2tet+lambda*(-deltaS[0]/2.0+deltaS[1]+deltaS[2]) > 1.0)
                lambda = (1.0-xfe2tet)/(-deltaS[0]/2.0+deltaS[1]+deltaS[2]);

            /* xal3tet = 1.0 - r[1] - r[2] - r[3] - s[1]; */
            if      ((-deltaS[1]) != 0.0 && xal3tet+lambda*(-deltaS[1]) < 0.0)
                lambda = -xal3tet/(-deltaS[1]);
            else if ((-deltaS[1]) != 0.0 && xal3tet+lambda*(-deltaS[1]) > 1.0)
                lambda = (1.0-xal3tet)/(-deltaS[1]);

            /* xfe3tet = r[3] - s[2]; */
            if      ((-deltaS[2]) != 0.0 && xfe3tet+lambda*(-deltaS[2]) < 0.0)
                lambda = -xfe3tet/(-deltaS[2]);
            else if ((-deltaS[2]) != 0.0 && xfe3tet+lambda*(-deltaS[2]) > 1.0)
                lambda = (1.0-xfe3tet)/(-deltaS[2]);

            /* xmg2oct = (r[0] - s[0])/4.0; */
            if      ((-deltaS[0]/4.0) != 0.0 && xmg2oct+lambda*(-deltaS[0]/4.0) < 0.0)
                lambda = -xmg2oct/(-deltaS[0]/4.0);
            else if ((-deltaS[0]/4.0) != 0.0 && xmg2oct+lambda*(-deltaS[0]/4.0) > 1.0)
                lambda = (1.0-xmg2oct)/(-deltaS[0]/4.0);

            /* xfe2oct = (2.0 - r[0] + s[0] - 2.0*s[1] - 2.0*r[1] - 2.0*s[2])/4.0; */
            if      ((deltaS[0]/4.0-deltaS[1]/2.0-deltaS[2]/2.0) != 0.0 &&
                xfe2oct+lambda*(deltaS[0]/4.0-deltaS[1]/2.0-deltaS[2]/2.0) < 0.0)
                lambda = -xfe2oct/(deltaS[0]/4.0-deltaS[1]/2.0-deltaS[2]/2.0);
            else if ((deltaS[0]/4.0-deltaS[1]/2.0-deltaS[2]/2.0) != 0.0 &&
                xfe2oct+lambda*(deltaS[0]/4.0-deltaS[1]/2.0-deltaS[2]/2.0) > 1.0)
                lambda = (1.0-xfe2oct)/(deltaS[0]/4.0-deltaS[1]/2.0-deltaS[2]/2.0);

            /* xal3oct = (1.0 - r[1] - r[2] - r[3] + s[1])/2.0; */
            if      ((deltaS[1]/2.0) != 0.0 && xal3oct+lambda*(deltaS[1]/2.0) < 0.0)
                lambda = -xal3oct/(deltaS[1]/2.0);
            else if ((deltaS[1]/2.0) != 0.0 && xal3oct+lambda*(deltaS[1]/2.0) > 1.0)
                lambda = (1.0-xal3oct)/(deltaS[1]/2.0);

            /* xfe3oct = (r[3] + s[2])/2.0; */
            if      ((deltaS[2]/2.0) != 0.0 && xfe3oct+lambda*(deltaS[2]/2.0) < 0.0)
                lambda = -xfe3oct/(deltaS[2]/2.0);
            else if ((deltaS[2]/2.0) != 0.0 && xfe3oct+lambda*(deltaS[2]/2.0) > 1.0)
                lambda = (1.0-xfe3oct)/(deltaS[2]/2.0);

            /* Modify steplength if required to maintain feasibility */
            if (lambda < 1.0) for (i=0; i<NS; i++) s[i] = sOld[i] + lambda*deltaS[i];

            for (i=0; i<NS; i++) sNew[i] = s[i];
            iter++;
        }
        tOld = t;
        pOld = p;
        for (i=0; i<NR; i++) rOld[i] = r[i];

#ifdef DEBUG
        for (i=0; i<NS; i++) {
            if (dgds[i] > sqrt(DBL_EPSILON) && ABS(sOld[i]) > DBL_EPSILON) {
                printf("ERROR in SPINEL.C (function ORDER). Failed to converge!\n");
                if (iter >= MAX_ITER)
                    printf("  Iteration limit (%4d) exceeded.\n", iter);
                printf("  X2    = %13.6g, X3    = %13.6g, X4    = %13.6g\n",
                    r[0], r[1], r[2]);
                printf("  X5    = %13.6g\n", r[3]);
                printf("  s1    = %13.6g, s2    = %13.6g, s4    = %13.6g\n",
                    sOld[0], sOld[1], sOld[2]);
                printf("  dgds1 = %13.6g, dgds2 = %13.6g, dgds4 = %13.6g\n",
                    dgds[0], dgds[1], dgds[2]);
                printf("  X Al   oct: %13.6g  X Al   tet: %13.6g\n", xal3oct, xal3tet);
                printf("  X Cr   oct: %13.6g\n",                     xcr3oct         );
                printf("  X Mg   oct: %13.6g  X Mg   tet: %13.6g\n", xmg2oct, xmg2tet);
                printf("  X Fe2+ oct: %13.6g  X Fe2+ tet: %13.6g\n", xfe2oct, xfe2tet);
                printf("  X Fe3+ oct: %13.6g  X Fe3+ tet: %13.6g\n", xfe3oct, xfe3tet);
                printf("  X Ti   oct: %13.6g\n",                     xti4oct         );
                break;
            }
        }
#endif

        setTOld(tOld);
        setPOld(pOld);
        /* arrays (rOld, sOld, d2gds2, indexD2gds2) should be preserved automatically */

        SET_SITE_FRACTIONS
    }

    if (mask & FIRST  ) {   /* return s        */
        for (i=0; i<NS; i++) s[i] = sOld[i];
    }
    if (mask & SECOND ) {   /* compute ds/dr:  */
        double d2gdrds[NR][NS];
        gsl_matrix_view mvToDr = gsl_matrix_view_array((double *) dr, (size_t) NS, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NS);

        fillD2GDRDS

        /* original: dr[i][j] += - d2gds2[i][k]*d2gdrds[j][k]; */
        for (j=0; j<NR; j++) {
            gsl_vector_view vvToDr = gsl_matrix_column(&mvToDr.matrix, j),
            	vvToD2gdrds = gsl_matrix_row(&mvToD2gdrds.matrix, j);
            melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdrds.vector, &vvToDr.vector);
        }
    }
    if (mask & THIRD  ) {   /* compute ds/dt:  */
        double d2gdsdt[NS];
        gsl_vector_view vvToDt = gsl_vector_view_array(dt, (size_t) NS),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NS);

        fillD2GDSDT

        /* original: dt[i] += - d2gds2[i][j]*d2gdsdt[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdt.vector, &vvToDt.vector);

    }
    if (mask & FOURTH ) {   /* compute ds/dp:  */
        double d2gdsdp[NS];
        gsl_vector_view vvToDp = gsl_vector_view_array(dp, (size_t) NS),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NS);

        fillD2GDSDP

        /* original: dp[i] += - d2gds2[i][j]*d2gdsdp[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdp.vector, &vvToDp.vector);

    }
    if (mask & FIFTH  ) {   /* compute d2s/dr2 */
        double d2gdrds[NR][NS], d3gdr2ds[NR][NR][NS], d3gdrds2[NR][NS][NS],
            d3gds3[NS][NS][NS], dsdr[NS][NR], temp[NS];
        gsl_matrix_view mvToDsdr = gsl_matrix_view_array((double *) dsdr, (size_t) NS, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NS);
        gsl_vector_view vvToTemp = gsl_vector_view_array(temp, (size_t) NS);
        int k, l, m, n;

        fillD2GDRDS
        fillD3GDR2DS
        fillD3GDRDS2
        fillD3GDS3

        /* compute dsdr matrix */
        /* original: dsdr[i][j] += - d2gds2[i][k]*d2gdrds[j][k]; */
        for (j=0; j<NR; j++) {
            gsl_vector_view vvToDsdr = gsl_matrix_column(&mvToDsdr.matrix, j),
            	vvToD2gdrds = gsl_matrix_row(&mvToD2gdrds.matrix, j);
            melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdrds.vector, &vvToDsdr.vector);
        }

        /* compute dsdr2 cube */
        for (j=0; j<NR; j++) {
            for (k=0; k<NR; k++) {
                for (l=0; l<NS; l++) {
                    temp[l] = d3gdr2ds[j][k][l];
                    for (m=0; m<NS; m++) {
                        temp[l] += d3gdrds2[j][l][m]*dsdr[m][k]
                            + d3gdrds2[k][l][m]*dsdr[m][j];
                        for (n=0; n<NS; n++)
                            temp[l] += d3gds3[l][m][n]*dsdr[m][j]*dsdr[n][k];
                    }
                }
                /* original: dr2[i][j][k] += - d2gds2[i][l]*temp[l]; */
                melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
                for (l=0; l<NS; l++) dr2[l][j][k] = temp[l];
            }
        }

    }
    if (mask & SIXTH  ) {   /* compute d2s/drt */
        double d2gdrds[NR][NS], d2gdsdt[NS], d3gdrds2[NR][NS][NS],
            d3gdrdsdt[NR][NS], d3gds3[NS][NS][NS], d3gds2dt[NS][NS], dsdr[NS][NR],
            dsdt[NS], temp[NS];
        gsl_matrix_view mvToDsdr = gsl_matrix_view_array((double *) dsdr, (size_t) NS, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NS);
        gsl_vector_view vvToDsdt = gsl_vector_view_array(dsdt, (size_t) NS),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NS),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NS);
        int k, l, m;

        fillD2GDRDS
        fillD2GDSDT
        fillD3GDRDS2
        fillD3GDRDSDT
        fillD3GDS3
        fillD3GDS2DT

        /* compute dsdr matrix */
        /* original: dsdr[i][j] += - d2gds2[i][k]*d2gdrds[j][k]; */
        for (j=0; j<NR; j++) {
            gsl_vector_view vvToDsdr = gsl_matrix_column(&mvToDsdr.matrix, j),
            	vvToD2gdrds = gsl_matrix_row(&mvToD2gdrds.matrix, j);
            melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdrds.vector, &vvToDsdr.vector);
        }

        /* compute dsdt vector */
        /* original: dsdt[i] += - d2gds2[i][j]*d2gdsdt[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdt.vector, &vvToDsdt.vector);

        /* compute dsdrdt matrix */
        for (j=0; j<NR; j++) {
            for (k=0; k<NS; k++) {
                temp[k] = d3gdrdsdt[j][k];
                for (l=0; l<NS; l++) {
                    temp[k] += d3gdrds2[j][k][l]*dsdt[l] + d3gds2dt[k][l]*dsdr[l][j];
                    for (m=0; m<NS; m++) temp[k] += d3gds3[k][l][m]*dsdr[l][j]*dsdt[m];
                }
            }
            /* original: drt[i][j] += - d2gds2[i][k]*temp[k]; */
            melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
            for (k=0; k<NS; k++) drt[k][j] = temp[k];
        }

    }
    if (mask & SEVENTH) {   /* compute d2s/drp */
        double d2gdrds[NR][NS], d2gdsdp[NS], d3gdrds2[NR][NS][NS],
            d3gdrdsdp[NR][NS], d3gds3[NS][NS][NS], d3gds2dp[NS][NS], dsdr[NS][NR],
            dsdp[NS], temp[NS];
        gsl_matrix_view mvToDsdr = gsl_matrix_view_array((double *) dsdr, (size_t) NS, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NS);
        gsl_vector_view vvToDsdp = gsl_vector_view_array(dsdp, (size_t) NS),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NS),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NS);
        int k, l, m;

        fillD2GDRDS
        fillD2GDSDP
        fillD3GDRDS2
        fillD3GDRDSDP
        fillD3GDS3
        fillD3GDS2DP

        /* compute dsdr matrix */
        /* original: dsdr[i][j] += - d2gds2[i][k]*d2gdrds[j][k]; */
        for (j=0; j<NR; j++) {
            gsl_vector_view vvToDsdr = gsl_matrix_column(&mvToDsdr.matrix, j),
            	vvToD2gdrds = gsl_matrix_row(&mvToD2gdrds.matrix, j);
            melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdrds.vector, &vvToDsdr.vector);
        }

        /* compute dsdp vector */
        /* original: dsdp[i] += - d2gds2[i][j]*d2gdsdp[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdp.vector, &vvToDsdp.vector);

        /* compute dsdrdp matrix */
        for (j=0; j<NR; j++) {
            for (k=0; k<NS; k++) {
                temp[k] = d3gdrdsdp[j][k];
                for (l=0; l<NS; l++) {
                    temp[k] += d3gdrds2[j][k][l]*dsdp[l] + d3gds2dp[k][l]*dsdr[l][j];
                    for (m=0; m<NS; m++) temp[k] += d3gds3[k][l][m]*dsdr[l][j]*dsdp[m];
                }
            }
            /* original: drp[i][j] += - d2gds2[i][k]*temp[k]; */
            melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
            for (k=0; k<NS; k++) drp[k][j] = temp[k];
        }

    }
    if (mask & EIGHTH ) {   /* compute d2s/dt2 */
        double d2gdsdt[NS], d3gds3[NS][NS][NS], d3gds2dt[NS][NS], d3gdsdt2[NS],
            dsdt[NS], temp[NS];
        gsl_vector_view vvToDsdt = gsl_vector_view_array(dsdt, (size_t) NS),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NS),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NS);
        int k, l;

        fillD2GDSDT
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDSDT2

        /* compute dsdt vector */
        /* original: dsdt[i] += - d2gds2[i][j]*d2gdsdt[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdt.vector, &vvToDsdt.vector);

        /* compute dsdt2 vector */
        for (j=0; j<NS; j++) {
            temp[j] = d3gdsdt2[j];
            for (k=0; k<NS; k++) {
                temp[j] +=  2.0*d3gds2dt[j][k]*dsdt[k];
                for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdt[l];
            }
        }
        /* original: dt2[i] += - d2gds2[i][j]*temp[j]; */
        melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
        for (j=0; j<NS; j++) dt2[j] = temp[j];

    }
    if (mask & NINTH  ) {   /* compute d2s/dtp */
        double d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS][NS][NS], d3gds2dt[NS][NS],
            d3gds2dp[NS][NS], d3gdsdtdp[NS], dsdt[NS], dsdp[NS], temp[NS];
        gsl_vector_view vvToDsdt = gsl_vector_view_array(dsdt, (size_t) NS),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NS),
            vvToDsdp = gsl_vector_view_array(dsdp, (size_t) NS),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NS),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NS);
        int k, l;

        fillD2GDSDT
        fillD2GDSDP
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDS2DP
        fillD3GDSDTDP

        /* compute dsdt vector */
        /* original: dsdt[i] += - d2gds2[i][j]*d2gdsdt[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdt.vector, &vvToDsdt.vector);

        /* compute dsdp vector */
        /* original: dsdp[i] += - d2gds2[i][j]*d2gdsdp[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdp.vector, &vvToDsdp.vector);

        /* compute dsdtp vector */
        for (j=0; j<NS; j++) {
            temp[j] = d3gdsdtdp[j];
            for (k=0; k<NS; k++) {
                temp[j] += d3gds2dt[j][k]*dsdp[k] + d3gds2dp[j][k]*dsdt[k];
                for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdp[l];
            }
        }
        /* original: dtp[i] += - d2gds2[i][j]*temp[j]; */
        melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
        for (j=0; j<NS; j++) dtp[j] = temp[j];

    }
    if (mask & TENTH  ) {   /* compute d2s/dp2 */
        double d2gdsdp[NS], d3gds3[NS][NS][NS], d3gds2dp[NS][NS], d3gdsdp2[NS],
            dsdp[NS], temp[NS];
        gsl_vector_view vvToDsdp = gsl_vector_view_array(dsdp, (size_t) NS),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NS),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NS);
        int k, l;

        fillD2GDSDP
        fillD3GDS3
        fillD3GDS2DP
        fillD3GDSDP2

        /* compute dsdp vector */
        /* original: dsdp[i] += - d2gds2[i][j]*d2gdsdp[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdp.vector, &vvToDsdp.vector);

        /* compute dsdp2 vector */
        for (j=0; j<NS; j++) {
            temp[j] = d3gdsdp2[j];
            for (k=0; k<NS; k++) {
                temp[j] +=  2.0*d3gds2dp[j][k]*dsdp[k];
                for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdp[k]*dsdp[l];
            }
        }
        /* original: dp2[i] += - d2gds2[i][j]*temp[j]; */
        melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
        for (j=0; j<NS; j++) dp2[j] = temp[j];

    }

}

/*
 *=============================================================================
 * Public functions:
 *    mask  -  bitwise mask for selecting output
 *    t     -  Temperature (K)
 *    p     -  Pressure (bars)
 *    *x    -  (pointer to x[]) Array of independent compositional variables
 */
int
testSpn(int mask, double t, double p,
    int na,          /* Expected number of endmember components                 */
    int nr,          /* Expected number of independent compositional variables  */
    char **names,    /* array of strings of names of endmember components       */
    char **formulas, /* array of strings of formulas of endmember components    */
    double *r,       /* array of indepependent compos variables, check bounds   */
    double *m)       /* array of moles of endmember components, check bounds    */
{
    const char *phase = "spinel.c";
    const char *NAMES[NA]    = { "chromite", "hercynite", "magnetite", "spinel", "ulvospinel" };
    const char *FORMULAS[NA] = { "FeCr2O4",  "FeAl2O4",   "Fe3O4",     "MgAl2O4","Fe2TiO4"    };
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
        result = result && (r[0] >= 0.0) && (r[0] <= 1.0+r[2]);
        result = result && (r[1] >= 0.0) && (r[1] <= 1.0);
        result = result && (r[2] >= 0.0) && (r[2] <= 1.0);
        result = result && (r[3] >= 0.0) && (r[3] <= 1.0);
        result = result && (1.0-r[0]-r[1]-r[2]-r[3] >= -r[0])
                                        && (1.0-r[0]-r[1]-r[2]-r[3] <= 1.0);
    }
    /* Check bounds on moles of endmember components */
    if (mask & SIXTH) {
        for (i=0, sum=0.0; i<NA; i++) sum += m[i];
        result = result && (sum >= 0.0);
        if (sum > 0.0) {
            result = result && (m[3]/sum >= 0.0) && (m[3]/sum <= 1.0+m[4]/sum);
            result = result && (m[0]/sum >= 0.0) && (m[0]/sum <= 1.0);
            result = result && (m[4]/sum >= 0.0) && (m[4]/sum <= 1.0);
            result = result && (m[2]/sum >= 0.0) && (m[2]/sum <= 1.0);
            result = result && (m[1]/sum >= -m[3]/sum) && (m[1]/sum <= 1.0);
        }
    }

    return result;
}

void
conSpn(int inpMask, int outMask, double t, double p,
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
    (3)  THIRD            FOURTH | SEVENTH | EIGHTH

    (1) converts a vector of moles of elements into a vector of moles of
            endmember spinel components.
    (2) calculates from a vector of moles of endmember components, one or
            all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
    (3) calculates from a vector of independent compositional variables
            mole fractions of endmember components and/or the Jacobian matrix
            dx[]/dr[] or a vector of ordering parameters (returned in *x)

    In this routine it is assumed that the elements are in the order of atomic
    numbers and that the order of spinel components has been verified as:
            m[0] = chromite   (FeCr2O4) ,
            m[1] = hercynite  (FeAl2O4),
            m[2] = magnetite  (Fe3O4),
            m[3] = spinel     (MgAl2O4),
            m[4] = ulvospinel (Fe2TiO4),

    ----------------------------------------------------------------------------*/

    int i, j, k;

    if (inpMask == FIRST && outMask == SECOND) {
        /* Converts a vector of moles of elements into a vector of moles of
       end-member components.                                                 */
        double sumcat, sumchg, fe2, fe3, proj, fTet, fOct;
        static const int Mg = 12;
        static const int Al = 13;
        static const int Ti = 22;
        static const int Cr = 24;
        static const int Mn = 25;
        static const int Fe = 26;
        static const int Co = 27;
        static const int Ni = 28;

        sumcat = e[Mg] + e[Al] + e[Ti] + e[Cr] + e[Mn] + e[Fe] + e[Co] + e[Ni];

        sumchg = 2.0*e[Mg] + 3.0*e[Al] + 4.0*e[Ti] + 3.0*e[Cr] + 2.0*e[Mn]
           + 2.0*e[Co] + 2.0*e[Ni];
        fe3 = 8.0*sumcat/3.0 - sumchg - 2.0*e[Fe];
        fe2 = e[Fe] - fe3;

        if (fe3 < 0.0) {
            fe3 = DBL_EPSILON;
            fe2 = e[Fe];
            fTet = ((fe2+e[Mg]) > 0.0) ?
             ((1.0/3.0)*sumcat + e[Ti])/(fe2 + e[Mg]) : 1.0;
            fOct = ((2.0*e[Ti]+e[Al]+e[Cr]) > 0.0) ?
                            (2.0/3.0)*sumcat/(2.0*e[Ti]+e[Al]+e[Cr]) : 1.0;
        } else {
            fTet = 1.0;
            fOct = 1.0;
        }

        proj = ((fe2 + e[Mg] + e[Mn]) == 0.0) ? 1.0 :
           (fe2 + e[Mn] + e[Mg] + e[Co] + e[Ni])/(fe2 + e[Mg] + e[Mn]);

        m[0] = e[Cr]*fOct/2.0;                               /* Moles of FeCr2O4             */
#ifdef TRUE_xMELTS
        m[1] = (fe2*fTet+e[Mn]*fTet)*proj-e[Cr]*fOct/2.0-fe3/2.0-2.0*e[Ti]*fOct;
#else
        m[1] = (fe2*fTet+e[Mn]*fTet-e[Cr]*fOct/2.0-fe3/2.0-2.0*e[Ti]*fOct)*proj;
#endif
                                                         /* Moles of FeAl2O4 and MnAl2O4 */
        m[2] = fe3/2.0;                                      /* Moles of Fe3O4               */
        m[3] = e[Mg]*proj*fTet;                              /* Moles of MgAl2O4             */
        m[4] = e[Ti]*fOct;                                   /* Moles of Fe2TiO4             */

        if (m[1] < 0.0 && fabs(m[1]) < sqrt(DBL_EPSILON)) m[1] = 0.0;

    } else if (inpMask == SECOND) {
        double sum;

        if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
            printf("Illegal call to conSpn with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        for (i=0, sum=0.0; i<NA; i++) sum += m[i];

        if (outMask & THIRD) {
            /* Converts a vector of moles of end-member components (m) into a vector
         of independent compositional variables (r) required as input for the
         remaining public functions.                                          */
            r[0] = (sum != 0.0) ? m[3]/sum : 0.0;  /* X2 = X MgAl2O4 */
            r[1] = (sum != 0.0) ? m[0]/sum : 0.0;  /* X3 = X FeCr2O4 */
            r[2] = (sum != 0.0) ? m[4]/sum : 0.0;  /* X4 = X Fe2TiO4 */
            r[3] = (sum != 0.0) ? m[2]/sum : 0.0;  /* X5 = X Fe3O4   */
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
                    dm[0][j] = (j == 3) ? (1.0 - m[3]/sum)/sum : -m[3]/SQUARE(sum);
                    dm[1][j] = (j == 0) ? (1.0 - m[0]/sum)/sum : -m[0]/SQUARE(sum);
                    dm[2][j] = (j == 4) ? (1.0 - m[4]/sum)/sum : -m[4]/SQUARE(sum);
                    dm[3][j] = (j == 2) ? (1.0 - m[2]/sum)/sum : -m[2]/SQUARE(sum);
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
                        d2m[0][j][k]  = 2.0*m[3]/CUBE(sum);
                        d2m[0][j][k] -= (j == 3) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[0][j][k] -= (k == 3) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[1][j][k]  = 2.0*m[0]/CUBE(sum);
                        d2m[1][j][k] -= (j == 0) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[1][j][k] -= (k == 0) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[2][j][k]  = 2.0*m[4]/CUBE(sum);
                        d2m[2][j][k] -= (j == 4) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[2][j][k] -= (k == 4) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[3][j][k]  = 2.0*m[2]/CUBE(sum);
                        d2m[3][j][k] -= (j == 2) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[3][j][k] -= (k == 2) ? 1.0/SQUARE(sum) : 0.0;
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
                            d3m[0][j][k][l]  = -6.0*m[3]/QUARTIC(sum);
                            d3m[0][j][k][l] += (j == 3) ? 2.0/CUBE(sum) : 0.0;
                            d3m[0][j][k][l] += (k == 3) ? 2.0/CUBE(sum) : 0.0;
                            d3m[0][j][k][l] += (l == 3) ? 2.0/CUBE(sum) : 0.0;
                            d3m[1][j][k][l]  = -6.0*m[0]/QUARTIC(sum);
                            d3m[1][j][k][l] += (j == 0) ? 2.0/CUBE(sum) : 0.0;
                            d3m[1][j][k][l] += (k == 0) ? 2.0/CUBE(sum) : 0.0;
                            d3m[1][j][k][l] += (l == 0) ? 2.0/CUBE(sum) : 0.0;
                            d3m[2][j][k][l]  = -6.0*m[4]/QUARTIC(sum);
                            d3m[2][j][k][l] += (j == 4) ? 2.0/CUBE(sum) : 0.0;
                            d3m[2][j][k][l] += (k == 4) ? 2.0/CUBE(sum) : 0.0;
                            d3m[2][j][k][l] += (l == 4) ? 2.0/CUBE(sum) : 0.0;
                            d3m[3][j][k][l]  = -6.0*m[2]/QUARTIC(sum);
                            d3m[3][j][k][l] += (j == 2) ? 2.0/CUBE(sum) : 0.0;
                            d3m[3][j][k][l] += (k == 2) ? 2.0/CUBE(sum) : 0.0;
                            d3m[3][j][k][l] += (l == 2) ? 2.0/CUBE(sum) : 0.0;
                        }
                    }
                }
            }
        }

    } else if (inpMask == THIRD) {

        if (outMask & ~(FOURTH | SEVENTH | EIGHTH))
            printf("Illegal call to conFld with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        if (outMask & FOURTH) {
            /* Converts a vector of independent compositional variables (r) into a
         vector of mole fractions of endmember components (x).                */
            x[0] = r[1];
            x[1] = 1.0 - r[0] - r[1] - r[2] - r[3];
            x[2] = r[3];
            x[3] = r[0];
            x[4] = r[2];
        }

        if (outMask & SEVENTH) {
            /* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
            for (i=0; i<NA; i++) for (j=0; j<NR; j++) dr[i][j] = 0.0;
            dr[0][1] =  1.0;
            dr[1][0] = -1.0; dr[1][1] = -1.0; dr[1][2] = -1.0; dr[1][3] = -1.0;
            dr[2][3] =  1.0;
            dr[3][0] =  1.0;
            dr[4][2] =  1.0;
        }

        if (outMask & EIGHTH) {
            /* Computes a vector of ordering parameters and returns it in x[] */
            order(FIRST, t, p, r, x, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
        }

    } else  {
        printf("Illegal call to conSpn with inpMask = %o and outMask = %o\n",
            inpMask, outMask);
    }

}

void
dispSpn(int mask, double t, double p, double *x,
    char **formula            /* Mineral formula for interface display MASK: 1 */
    )
{
    double *r = x;
    static char masterString[] = {
/*             1111111111222222222233333333334444444444555555555566666666667
     01234567890123456789012345678901234567890123456789012345678901234567890 */
        "Fe''_.__Mg_.__Fe'''_.__Al_.__Cr_.__Ti_.__O4" };

    if (mask & FIRST) {
        char *string = strcpy((char *) malloc((size_t) (strlen(masterString)+1)*sizeof(char)), masterString);
        double totAl, totFe2, totFe3, totMg, totTi, totCr;
        char n[5];
        int i;

        totFe2 = 1.0-r[0]+r[2];
        totMg  = r[0];
        totFe3 = 2.0*r[3];
        totAl  = 2.0*(1.0-r[1]-r[2]-r[3]);
        totCr  = 2.0*r[1];
        totTi  = r[2];

        (void) snprintf(n, 5, "%4.2f", totFe2); for (i=0; i<4; i++) string[ 4+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totMg);  for (i=0; i<4; i++) string[10+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totFe3); for (i=0; i<4; i++) string[19+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totAl);  for (i=0; i<4; i++) string[25+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totCr);  for (i=0; i<4; i++) string[31+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totTi);  for (i=0; i<4; i++) string[37+i] = n[i];

        *formula = string;
    }
}

void
actSpn(int mask, double t, double p, double *x,
    double *a,  /* (pointer to a[]) activities              BINARY MASK: 0001 */
    double *mu, /* (pointer to mu[]) chemical potentials    BINARY MASK: 0010 */
    double **dx /* (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100 */
    )           /* exclusion criteria applied to results if BINARY MASK: 1000 */
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS], g, dgdr[NR];
    double fr[NA][NR];
    int i, j;

    for(i=0; i<NA; i++) {
     fr[i][0] = FR2(i); /* X2 */
     fr[i][1] = FR3(i); /* X3 */
     fr[i][2] = FR4(i); /* X4 */
     fr[i][3] = FR5(i); /* X5 */
    }

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    g       = G;
    dgdr[0] = DGDR0;
    dgdr[1] = DGDR1;
    dgdr[2] = DGDR2;
    dgdr[3] = DGDR3;

    /* activities for library */
    if (!mask && a != NULL) {
        for(i=0; i<NA; i++) {
       for (a[i]=g, j=0; j<NR; j++) a[i] += fr[i][j]*dgdr[j];
       a[i] = exp(a[i]/(R*t));
        }
    }

    if (mask & FIRST) {
        double a0[NA];

        pureSpn(FIRST, t, p,
       a0,              (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        for(i=0; i<NA; i++) {
       for (a[i]=g, j=0; j<NR; j++) a[i] += fr[i][j]*dgdr[j];
       a[i] = exp(a[i]/(R*t));
       if (a0[i] != 0.0) a[i] = a[i]/a0[i];
        }

    }

    if (mask & SECOND) {
        double mu0[NA];

        pureSpn(SECOND, t, p,
       (double *) NULL, mu0,             (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        for(i=0; i<NA; i++) {
       for (mu[i]=g-mu0[i], j=0; j<NR; j++) mu[i] += fr[i][j]*dgdr[j];
        }

    }

    if (mask & THIRD) {
        double d2gdr2[NR][NR], d2gdrds[NR][NS], d2gds2[NS][NS], dsdr[NS][NR],
            dfrdr[NA][NR], gs[NA][NS], dgsds[NA][NS], sum, a0[NA];
        int k, l;

        fillD2GDR2
        fillD2GDRDS
        fillD2GDS2

        pureSpn(FIRST, t, p,
       a0,              (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        for(i=0; i<NA; i++) {
       gs[i][0] = GS1(i); /* s1 */
       gs[i][1] = GS2(i); /* s2 */
       gs[i][2] = GS4(i); /* s4 */
       dfrdr[i][0] = DFR2DR2(i); /* X2 */
       dfrdr[i][1] = DFR3DR3(i); /* X3 */
       dfrdr[i][2] = DFR4DR4(i); /* X4 */
       dfrdr[i][3] = DFR5DR5(i); /* X5 */
       dgsds[i][0] = DGS1DS1(i); /* s1 */
       dgsds[i][1] = DGS2DS2(i); /* s2 */
       dgsds[i][2] = DGS4DS4(i); /* s4 */
        }

        order(SECOND, t, p, r,
                    NULL, dsdr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NA; i++) {
            for (k=0; k<NR; k++) {
                /* compute activity of the i-th component */
                for (dx[i][k]=g, j=0; j<NR; j++) dx[i][k] += fr[i][j]*dgdr[j];
                dx[i][k] = exp(dx[i][k]/(R*t));
                if (a0[i] != 0.0) dx[i][k]  = dx[i][k]/a0[i];

                /* compute derivative of i-th activity with respect to r(k) */
                sum = (1.0+dfrdr[i][k])*dgdr[k];
                for (j=0; j<NR; j++) {
                    sum += fr[i][j]*d2gdr2[j][k];
                    for (l=0; l<NS; l++) sum += fr[i][j]*d2gdrds[j][l]*dsdr[l][k];
                }
                for (j=0; j<NS; j++) {
                    sum += gs[i][j]*d2gdrds[k][j];
                    for (l=0; l<NS; l++) sum += gs[i][j]*d2gds2[j][l]*dsdr[l][k];
                }
                dx[i][k] *= sum/(R*t);
            }
        }
    }

    if (mask & FOURTH) {
        /* implement exclusion criteria on quantities for preclb routines         */
        static const double exclusion[NA] = {
       0.05,  /* exclusion criteria on the mole fraction of Cr3+ (chromite)   */
       0.05,  /* exclusion criteria on the mole fraction of Al3+ (hercynite)  */
       0.05,  /* exclusion criteria on the mole fraction of Fe3+ (magnetite)  */
       0.05,  /* exclusion criteria on the mole fraction of Al3+ (spinel)     */
       0.05   /* exclusion criteria on the mole fraction of Ti4+ (ulvospinel) */
        };
        double x[NA], totFe2, totMg, sum;

        totFe2 = 1.0 - r[0] + r[2];
        totMg  = r[0];
        sum    = 1.0 + r[2];

        x[0] = (sum != 0.0) ? (totFe2/sum)*r[1] : 0.0;                 /* Cr3+/2 */
        x[1] = (sum != 0.0) ? (totFe2/sum)*(1.0-r[1]-r[2]-r[3]) : 0.0; /* Al3+/2 */
        x[2] = (sum != 0.0) ? (totFe2/sum)*r[3] : 0.0;                 /* Fe3+/2 */
        x[3] = (sum != 0.0) ? (totMg /sum)*(1.0-r[1]-r[2]-r[3]) : 0.0; /* Al3+/2 */
        x[4] = (sum != 0.0) ? (totFe2/sum)*r[2] : 0.0;                 /* Ti4+   */

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
gmixSpn(int mask, double t, double p, double *x,
    double *gmix, /* Gibbs energy of mixing             BINARY MASK: 0001 */
    double *dx,   /* (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010 */
    double **dx2, /* (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100 */
    double ***dx3 /* (pointer to dx3[][][]) d3(g)/d(x[])3 NARY MASK: 1000 */
    )
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS];

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    if (mask & FIRST) {
        double ends[NA];

        *gmix = G;

        pureSpn(THIRD, t, p,
       (double *) NULL, (double *) NULL, ends,            (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        *gmix -= ENDMEMBERS;
    }

    if(mask & SECOND) {
        double ends[NA];

        dx[0] = DGDR0;
        dx[1] = DGDR1;
        dx[2] = DGDR2;
        dx[3] = DGDR3;

        pureSpn(THIRD, t, p,
       (double *) NULL, (double *) NULL, ends,            (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        dx[0] -= DENDDR0;
        dx[1] -= DENDDR1;
        dx[2] -= DENDDR2;
        dx[3] -= DENDDR3;
    }

    if(mask & THIRD) {
        double d2gdr2[NR][NR], d2gdrds[NR][NS], d2gds2[NS][NS], dsdr[NS][NR];
        int i, j, k, l;

        fillD2GDR2
        fillD2GDRDS
        fillD2GDS2

        order(SECOND, t, p, r,
                    NULL, dsdr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dx2[i][j] = d2gdr2[i][j];
                for (k=0; k<NS; k++) {
                    dx2[i][j] += d2gdrds[i][k]*dsdr[k][j] + d2gdrds[j][k]*dsdr[k][i];
                    for (l=0; l<NS; l++) dx2[i][j] += d2gds2[k][l]*dsdr[k][i]*dsdr[l][j];
                }
            }
        }

    }

    if(mask & FOURTH) {
        double d3gdr3[NR][NR][NR], d3gdr2ds[NR][NR][NS], d3gdrds2[NR][NS][NS];
        double d3gds3[NS][NS][NS], dsdr[NS][NR];
        int i, j, k, l, m, n;

        fillD3GDR3
        fillD3GDR2DS
        fillD3GDRDS2
        fillD3GDS3

        order(SECOND, t, p, r,
                    NULL, dsdr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                for (k=0; k<NR; k++) {
                    dx3[i][j][k] = d3gdr3[i][j][k];
                    for (l=0; l<NS; l++) {
                        dx3[i][j][k] += d3gdr2ds[i][j][l]*dsdr[l][k] +
                            d3gdr2ds[j][k][l]*dsdr[l][i] + d3gdr2ds[k][i][l]*dsdr[l][j];
                        for (m=0; m<NS; m++) {
                            dx3[i][j][k] +=
                                d3gdrds2[i][l][m]*dsdr[l][j]*dsdr[m][k] +
                                d3gdrds2[j][l][m]*dsdr[l][k]*dsdr[m][i] +
                                d3gdrds2[k][l][m]*dsdr[l][i]*dsdr[m][j];
                            for (n=0; n<NS; n++)
                                dx3[i][j][k] +=
                                    d3gds3[l][m][n]*dsdr[l][i]*dsdr[m][j]*dsdr[n][k];
                        }
                    }
                }
            }
        }

    }

}

void
hmixSpn(int mask, double t, double p, double *x,
    double *hmix /* Enthalpy of mixing BINARY MASK: 1 */
    )
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS], ends[NA];

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    *hmix = (G) + t*(S);

    pureSpn(FOURTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, ends,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

    *hmix -= ENDMEMBERS;
}

void
smixSpn(int mask, double t, double p, double *x,
    double *smix, /* Entropy of mixing                  BINARY MASK: 001 */
    double *dx,   /* (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010 */
    double **dx2  /* (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100 */
    )
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS];

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    if (mask & FIRST) {
        double ends[NA];

        *smix = S;

        pureSpn(FIFTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       ends,            (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        *smix -= ENDMEMBERS;
    }

    if(mask & SECOND) {
        double d2gdrds[NR][NS], d2gdrdt[NR], d2gds2[NS][NS], d2gdsdt[NS],
            dsdr[NS][NR], dsdt[NS], ends[NA];
        int i, k, l;

        fillD2GDRDS
        fillD2GDRDT
        fillD2GDS2
        fillD2GDSDT

        order(SECOND | THIRD, t, p, r,
                    NULL, dsdr, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            dx[i] = d2gdrdt[i];
            for (k=0; k<NS; k++) {
                dx[i] += d2gdrds[i][k]*dsdt[k] + d2gdsdt[k]*dsdr[k][i];
                for (l=0; l<NS; l++) dx[i] += d2gds2[k][l]*dsdt[k]*dsdr[l][i] ;
            }
            dx[i] *= -1.0;
        }

        pureSpn(FIFTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       ends,            (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        dx[0] -= DENDDR0;
        dx[1] -= DENDDR1;
        dx[2] -= DENDDR2;
        dx[3] -= DENDDR3;
    }

    if(mask & THIRD) {
        double d2gdrds[NR][NS], d2gds2[NS][NS], d2gdsdt[NS], d3gdr2ds[NR][NR][NS],
            d3gdr2dt[NR][NR], d3gdrds2[NR][NS][NS], d3gdrdsdt[NR][NS],
            d3gds3[NS][NS][NS], d3gds2dt[NS][NS], dsdr[NS][NR], dsdt[NS],
            d2sdr2[NS][NR][NR], d2sdrdt[NS][NR];
        int i, j, k, l, m;

        fillD2GDRDS
        fillD2GDS2
        fillD2GDSDT
        fillD3GDR2DS
        fillD3GDR2DT
        fillD3GDRDS2
        fillD3GDRDSDT
        fillD3GDS3
        fillD3GDS2DT

        order(SECOND | THIRD | FIFTH | SIXTH, t, p, r,
                    NULL, dsdr, dsdt, NULL, d2sdr2, d2sdrdt, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dx2[i][j] = d3gdr2dt[i][j];
                for (k=0; k<NS; k++) {
                    dx2[i][j] += d3gdr2ds[i][j][k]*dsdt[k]
                     + d3gdrdsdt[i][k]*dsdr[k][j]
                     + d3gdrdsdt[j][k]*dsdr[k][i]
                     + d2gdsdt[k]*d2sdr2[k][i][j]
                     + d2gdrds[i][k]*d2sdrdt[k][j]
                     + d2gdrds[j][k]*d2sdrdt[k][i];
                    for (l=0; l<NS; l++) {
                        dx2[i][j] += d3gdrds2[i][k][l]*dsdr[k][j]*dsdt[l]
                       + d3gdrds2[j][k][l]*dsdr[k][i]*dsdt[l]
                       + d2gds2[k][l]*d2sdr2[k][i][j]*dsdt[l]
                       + d3gds2dt[k][l]*dsdr[k][i]*dsdr[l][j]
                       + d2gds2[k][l]*dsdr[k][i]*d2sdrdt[l][j]
                       + d2gds2[k][l]*dsdr[k][j]*d2sdrdt[l][i];
                        for (m=0; m<NS; m++)
                            dx2[i][j] += d3gds3[k][l][m]*dsdr[k][i]*dsdr[l][j]*dsdt[m];
                    }
                }
                dx2[i][j] *= -1.0;
            }
        }

    }

}

void
cpmixSpn(int mask, double t, double p, double *x,
    double *cpmix, /* Heat capacity of mixing               BINARY MASK: 001 */
    double *dt,    /* d(cp)/d(t)                            BINARY MASK: 010 */
    double *dx     /* d(cp)/d(x[])d(t)                      BINARY MASK: 100 */
    )
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS], dsdt[NS], d2gdsdt[NS], d2gds2[NS][NS], d2gdt2;
    int i, j;

    order(FIRST | THIRD, t, p, r,
                s, NULL, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    fillD2GDS2
    fillD2GDSDT
    d2gdt2  = D2GDT2;

    if (mask & FIRST) {
        double ends[NA];

        *cpmix = d2gdt2;
        for (i=0; i<NS; i++) {
            *cpmix += 2.0*d2gdsdt[i]*dsdt[i];
            for (j=0; j<NS; j++) *cpmix += d2gds2[i][j]*dsdt[i]*dsdt[j];
        }
        *cpmix *= -t;

        pureSpn(SIXTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, ends,            (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        *cpmix -= ENDMEMBERS;
    }

    if(mask & SECOND) {
        double d3gds3[NS][NS][NS], d3gds2dt[NS][NS], d3gdsdt2[NS], d2sdt2[NS],
            temp, ends[NA];
        double d3gdt3 = D3GDT3;
        int k;

        fillD3GDS3
        fillD3GDS2DT
        fillD3GDSDT2

        order(EIGHTH, t, p, r,
                    NULL, NULL, NULL, NULL, NULL, NULL, NULL, d2sdt2, NULL, NULL);

        /* compute d2gdt2 */
        temp = d2gdt2;
        for (i=0; i<NS; i++) {
            temp += 2.0*d2gdsdt[i]*dsdt[i];
            for (j=0; j<NS; j++) temp += d2gds2[i][j]*dsdt[i]*dsdt[j];
        }

        *dt = d3gdt3;
        for (i=0; i<NS; i++) {
            *dt += 3.0*d3gdsdt2[i]*dsdt[i] + 3.0*d2gdsdt[i]*d2sdt2[i];
            for (j=0; j<NS; j++) {
                *dt += 3.0*d2gds2[i][j]*dsdt[i]*d2sdt2[j]
             + 3.0*d3gds2dt[i][j]*dsdt[i]*dsdt[j];
                for (k=0; k<NS; k++) *dt += d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdt[k];
            }
        }
        *dt = -t*(*dt) - temp;

        pureSpn(SEVENTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, ends,            (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        *dt -= ENDMEMBERS;
    }

    if(mask & THIRD) {
        double d3gds3[NS][NS][NS], d3gdrds2[NR][NS][NS], d3gdrdsdt[NR][NS],
            d3gds2dt[NS][NS], d2gdrds[NR][NS], d3gdrdt2[NR], d3gdsdt2[NS],
            dsdr[NS][NR], d2sdrdt[NS][NR], d2sdt2[NS], ends[NA];
        int k, l;

        fillD2GDRDS
        fillD3GDRDS2
        fillD3GDRDSDT
        fillD3GDRDT2
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDSDT2

        order(SECOND | SIXTH | EIGHTH, t, p, r,
                    NULL, dsdr, NULL, NULL, NULL, d2sdrdt, NULL, d2sdt2, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0,dx[i]=d3gdrdt2[i]; j<NS; j++) {
                dx[i] += d3gdsdt2[j]*dsdr[j][i] + 2.0*d2gdsdt[j]*d2sdrdt[j][i] +
                 2.0*d3gdrdsdt[i][j]*dsdt[j] + d2gdrds[i][j]*d2sdt2[j];
                for (k=0; k<NS; k++) {
                    dx[i] += d3gdrds2[i][j][k]*dsdt[j]*dsdt[k] +
                   2.0*d2gds2[j][k]*dsdt[j]*d2sdrdt[k][i] +
                   2.0*d3gds2dt[j][k]*dsdr[j][i]*dsdt[k] +
                   d2gds2[j][k]*dsdr[j][i]*d2sdt2[k];
                    for (l=0; l<NS; l++)
                        dx[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdt[k]*dsdt[l];
                }
            }
            dx[i] *= -t;
        }

        pureSpn(SIXTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, ends,            (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        dx[0] -= DENDDR0;
        dx[1] -= DENDDR1;
        dx[2] -= DENDDR2;
        dx[3] -= DENDDR3;
    }

}

void
vmixSpn(int mask, double t, double p, double *x,
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
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS];

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    if (mask & FIRST) {
        double ends[NA];

        *vmix = DGDP;

        pureSpn(EIGHTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, ends,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        *vmix -= ENDMEMBERS;
    }

    if(mask & SECOND) {
        double d2gdrds[NR][NS], d2gdrdp[NR], d2gds2[NS][NS], d2gdsdp[NS],
            dsdr[NS][NR], dsdp[NS], ends[NA];
        int i, j, k;

        fillD2GDRDS
        fillD2GDRDP
        fillD2GDS2
        fillD2GDSDP

        order(SECOND | FOURTH, t, p, r,
                    NULL, dsdr, NULL, dsdp, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            dx[i] = d2gdrdp[i];
            for (j=0; j<NS; j++) {
                dx[i] += d2gdrds[i][j]*dsdp[j] + d2gdsdp[j]*dsdr[j][i];
                for (k=0; k<NS; k++) dx[i] += d2gds2[j][k]*dsdp[j]*dsdr[k][i];
            }
        }

        pureSpn(EIGHTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, ends,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        dx[0] -= DENDDR0;
        dx[1] -= DENDDR1;
        dx[2] -= DENDDR2;
        dx[3] -= DENDDR3;
    }

    if(mask & THIRD) {
        double d2gdrds[NR][NS], d2gds2[NS][NS], d2gdsdp[NS], d3gdr2ds[NR][NR][NS],
            d3gdr2dp[NR][NR], d3gdrds2[NR][NS][NS], d3gdrdsdp[NR][NS],
            d3gds3[NS][NS][NS], d3gds2dp[NS][NS], dsdr[NS][NR], dsdp[NS],
            d2sdr2[NS][NR][NR], d2sdrdp[NS][NR];
        int i, j, k, l, m;

        fillD2GDRDS
        fillD2GDS2
        fillD2GDSDP
        fillD3GDR2DS
        fillD3GDR2DP
        fillD3GDRDS2
        fillD3GDRDSDP
        fillD3GDS3
        fillD3GDS2DP

        order(SECOND | FOURTH | FIFTH | SEVENTH, t, p, r,
                    NULL, dsdr, NULL, dsdp, d2sdr2, NULL, d2sdrdp,  NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dx2[i][j] = d3gdr2dp[i][j];
                for (k=0; k<NS; k++) {
                    dx2[i][j] += d3gdr2ds[i][j][k]*dsdp[k]
                     + d3gdrdsdp[i][k]*dsdr[k][j]
                     + d3gdrdsdp[j][k]*dsdr[k][i]
                     + d2gdsdp[k]*d2sdr2[k][i][j]
                     + d2gdrds[i][k]*d2sdrdp[k][j]
                     + d2gdrds[j][k]*d2sdrdp[k][i];
                    for (l=0; l<NS; l++) {
                        dx2[i][j] += d3gdrds2[i][k][l]*dsdr[k][j]*dsdp[l]
                       + d3gdrds2[j][k][l]*dsdr[k][i]*dsdp[l]
                       + d2gds2[k][l]*d2sdr2[k][i][j]*dsdp[l]
                       + d3gds2dp[k][l]*dsdr[k][i]*dsdr[l][j]
                       + d2gds2[k][l]*dsdr[k][i]*d2sdrdp[l][j]
                       + d2gds2[k][l]*dsdr[k][j]*d2sdrdp[l][i];
                        for (m=0; m<NS; m++)
                            dx2[i][j] += d3gds3[k][l][m]*dsdr[k][i]*dsdr[l][j]*dsdp[m];
                    }
                }
            }
        }

    }

    if(mask & FOURTH) {
        double d2gdtdp = D2GDTDP;
        double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], dsdt[NS], dsdp[NS],
            ends[NA];
        int i, j;

        fillD2GDS2
        fillD2GDSDT
        fillD2GDSDP

        order(THIRD | FOURTH, t, p, r,
                    NULL, NULL, dsdt, dsdp, NULL, NULL, NULL, NULL, NULL, NULL);

        *dt = d2gdtdp;
        for (i=0; i<NS; i++) {
            *dt += d2gdsdt[i]*dsdp[i] + d2gdsdp[i]*dsdt[i];
            for (j=0; j<NS; j++) *dt += d2gds2[i][j]*dsdt[i]*dsdp[j];
        }

        pureSpn(NINTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       ends,            (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        *dt -= ENDMEMBERS;
    }

    if(mask & FIFTH) {
        double d2gdp2 = D2GDP2;
        double d2gds2[NS][NS], d2gdsdp[NS], dsdp[NS], ends[NA];
        int i,j;

        fillD2GDS2
        fillD2GDSDP

        order(FOURTH, t, p, r,
                    NULL, NULL, NULL, dsdp, NULL, NULL, NULL, NULL, NULL, NULL);

        *dp = d2gdp2;
        for (i=0; i<NS; i++) {
            *dp += 2.0*d2gdsdp[i]*dsdp[i];
            for (j=0; j<NS; j++) *dp += d2gds2[i][j]*dsdp[i]*dsdp[j];
        }

        pureSpn(TENTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, ends,            (double *) NULL, (double *) NULL,
       (double *) NULL);

        *dp -= ENDMEMBERS;
    }

    if(mask & SIXTH) {
        double d3gdt2dp = D3GDT2DP;
        double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS][NS][NS],
            d3gds2dp[NS][NS], d3gds2dt[NS][NS], d3gdsdtdp[NS], d3gdsdt2[NS],
            dsdt[NS], dsdp[NS], d2sdt2[NS], d2sdtdp[NS], ends[NA];
        int i, j, k;

        fillD2GDS2
        fillD2GDSDT
        fillD2GDSDP
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDS2DP
        fillD3GDSDT2
        fillD3GDSDTDP

        order(THIRD | FOURTH | EIGHTH | NINTH, t, p, r,
                    NULL, NULL, dsdt, dsdp, NULL, NULL, NULL, d2sdt2, d2sdtdp, NULL);

        *dt2 = d3gdt2dp;
        for (i=0; i<NS; i++) {
            *dt2 += d3gdsdt2[i]*dsdp[i] + 2.0*d2gdsdt[i]*d2sdtdp[i]
                        + d2gdsdp[i]*d2sdt2[i] + 2.0*d3gdsdtdp[i]*dsdt[i];
            for (j=0; j<NS; j++) {
                *dt2 += 2.0*d3gds2dt[i][j]*dsdt[i]*dsdp[j]
                            + d2gds2[i][j]*d2sdt2[i]*dsdp[j]
                            + 2.0*d2gds2[i][j]*dsdt[i]*d2sdtdp[j]
                            + d3gds2dp[i][j]*dsdt[i]*dsdt[j];
                for (k=0; k<NS; k++) *dt2 += d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdp[k];
            }
        }

        pureSpn(ELEVENTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, ends,            (double *) NULL,
       (double *) NULL);

        *dt2 -= ENDMEMBERS;
    }

    if(mask & SEVENTH) {
        double d3gdtdp2 = D3GDTDP2;
        double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS][NS][NS],
            d3gds2dt[NS][NS], d3gds2dp[NS][NS], d3gdsdtdp[NS], d3gdsdp2[NS],
            dsdt[NS], dsdp[NS], d2sdtdp[NS], d2sdp2[NS], ends[NA];
        int i, j, k;

        fillD2GDS2
        fillD2GDSDT
        fillD2GDSDP
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDS2DP
        fillD3GDSDTDP
        fillD3GDSDP2

        order(THIRD | FOURTH | NINTH | TENTH, t, p, r,
                    NULL, NULL, dsdt, dsdp, NULL, NULL, NULL, NULL, d2sdtdp, d2sdp2);

        *dtdp = d3gdtdp2;
        for (i=0; i<NS; i++) {
            *dtdp += 2.0*d3gdsdtdp[i]*dsdp[i] + d2gdsdt[i]*d2sdp2[i]
             + 2.0*d2gdsdp[i]*d2sdtdp[i] + d3gdsdp2[i]*dsdt[i];
            for (j=0; j<NS; j++) {
                *dtdp += 2.0*d3gds2dp[i][j]*dsdt[i]*dsdp[j]
               + d2gds2[i][j]*dsdt[i]*d2sdp2[j]
               + 2.0*d2gds2[i][j]*d2sdtdp[i]*dsdp[j]
               + d3gds2dt[i][j]*dsdp[i]*dsdp[j];
                for (k=0; k<NS; k++) *dtdp += d3gds3[i][j][k]*dsdt[i]*dsdp[j]*dsdp[k];
            }
        }

        pureSpn(TWELFTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, ends,
       (double *) NULL);

        *dtdp -= ENDMEMBERS;
    }

    if(mask & EIGHTH) {
        double d3gdp3 = D3GDP3;
        double d2gds2[NS][NS], d2gdsdp[NS], d3gds3[NS][NS][NS], d3gds2dp[NS][NS],
            d3gdsdp2[NS], dsdp[NS], d2sdp2[NS], ends[NA];
        int i, j, k;

        fillD2GDS2
        fillD2GDSDP
        fillD3GDS3
        fillD3GDS2DP
        fillD3GDSDP2

        order(FOURTH | TENTH, t, p, r,
                    NULL, NULL, NULL, dsdp, NULL, NULL, NULL, NULL, NULL, d2sdp2);

        *dp2 = d3gdp3;
        for (i=0; i<NS; i++) {
            *dp2 += 3.0*d3gdsdp2[i]*dsdp[i] + 3.0*d2gdsdp[i]*d2sdp2[i];
            for (j=0; j<NS; j++) {
                *dp2 += 3.0*d2gds2[i][j]*dsdp[i]*d2sdp2[j]
                            + 3.0*d3gds2dp[i][j]*dsdp[i]*dsdp[j];
                for (k=0; k<NS; k++) *dp2 += d3gds3[i][j][k]*dsdp[i]*dsdp[j]*dsdp[k];
            }
        }

        pureSpn(THIRTEENTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       ends);

        *dp2 -= ENDMEMBERS;
    }

    if(mask & NINTH) {
        double d3gds3[NS][NS][NS], d3gdrds2[NR][NS][NS], d3gdrdsdt[NR][NS],
            d3gds2dp[NS][NS], d2gdrds[NR][NS], d3gdrdtdp[NR], d3gdsdtdp[NS],
            dsdt[NS], dsdp[NS], dsdr[NS][NR], d2sdrdt[NS][NR], d2sdrdp[NS][NR],
            d2gds2[NS][NS], d2gdsdt[NS], d3gdrdsdp[NR][NS], d2gdsdp[NS],
            d2sdtdp[NS], d3gds2dt[NS][NS], ends[NA];
        int i, j, k, l;

        fillD2GDRDS
        fillD2GDS2
        fillD2GDSDT
        fillD2GDSDP
        fillD3GDRDS2
        fillD3GDRDSDT
        fillD3GDRDSDP
        fillD3GDRDTDP
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDSDTDP
        fillD3GDS2DP

        order(SECOND | THIRD | FOURTH | SIXTH | SEVENTH | NINTH, t, p, r,
            NULL, dsdr, dsdt, dsdp, NULL, d2sdrdt, d2sdrdp, NULL, d2sdtdp, NULL);

        for (i=0; i<NR; i++) {
            for (j=0,dxdt[i]=d3gdrdtdp[i]; j<NS; j++) {
                dxdt[i] += d3gdsdtdp[j]*dsdr[j][i] + d2gdsdt[j]*d2sdrdp[j][i] +
                   d3gdrdsdt[i][j]*dsdp[j] + d2gdrds[i][j]*d2sdtdp[j] +
                   d3gdrdsdp[i][j]*dsdt[j] + d2gdsdp[j]*d2sdrdt[j][i];
                for (k=0; k<NS; k++) {
                    dxdt[i] += d3gdrds2[i][j][k]*dsdt[j]*dsdp[k] +
                     d2gds2[j][k]*dsdt[j]*d2sdrdp[k][i] +
                     d2gds2[j][k]*dsdp[j]*d2sdrdt[k][i] +
                     d3gds2dt[j][k]*dsdr[j][i]*dsdp[k] +
                     d3gds2dp[j][k]*dsdr[j][i]*dsdt[k] +
                     d2gds2[j][k]*dsdr[j][i]*d2sdtdp[k];
                    for (l=0; l<NS; l++)
                        dxdt[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdt[k]*dsdp[l];
                }
            }
        }

        pureSpn(NINTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       ends,            (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        dxdt[0] -= DENDDR0;
        dxdt[1] -= DENDDR1;
        dxdt[2] -= DENDDR2;
        dxdt[3] -= DENDDR3;
    }

    if(mask & TENTH) {
        double d3gds3[NS][NS][NS], d3gdrds2[NR][NS][NS], d3gds2dp[NS][NS],
            d2gdrds[NR][NS], dsdp[NS], dsdr[NS][NR], d2sdrdp[NS][NR], d2gds2[NS][NS],
            d3gdrdsdp[NR][NS], d3gdrdp2[NR], d3gdsdp2[NS], d2gdsdp[NS], d2sdp2[NS],
            ends[NA];
        int i, j, k, l;

        fillD2GDRDS
        fillD2GDS2
        fillD2GDSDP
        fillD3GDRDS2
        fillD3GDRDSDP
        fillD3GDRDP2
        fillD3GDS3
        fillD3GDS2DP
        fillD3GDSDP2

        order(SECOND | FOURTH | SEVENTH | TENTH, t, p, r,
            NULL, dsdr, NULL, dsdp, NULL, NULL, d2sdrdp, NULL, NULL, d2sdp2);

        for (i=0; i<NR; i++) {
            for (j=0,dxdp[i]=d3gdrdp2[i]; j<NS; j++) {
                dxdp[i] += d3gdsdp2[j]*dsdr[j][i] + d2gdsdp[j]*d2sdrdp[j][i] +
                   2.0*d3gdrdsdp[i][j]*dsdp[j] + d2gdrds[i][j]*d2sdp2[j] +
                   d2gdsdp[j]*d2sdrdp[j][i];
                for (k=0; k<NS; k++) {
                    dxdp[i] += d3gdrds2[i][j][k]*dsdp[j]*dsdp[k] +
                     2.0*d2gds2[j][k]*dsdp[j]*d2sdrdp[k][i] +
                     2.0*d3gds2dp[j][k]*dsdr[j][i]*dsdp[k] +
                     d2gds2[j][k]*dsdr[j][i]*d2sdp2[k];
                    for (l=0; l<NS; l++)
                        dxdp[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdp[k]*dsdp[l];
                }
            }
        }

        pureSpn(TENTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, ends,            (double *) NULL, (double *) NULL,
       (double *) NULL);

        dxdp[0] -= DENDDR0;
        dxdp[1] -= DENDDR1;
        dxdp[2] -= DENDDR2;
        dxdp[3] -= DENDDR3;
    }

}

/* end of file SPINEL.C */
