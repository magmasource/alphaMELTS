#ifndef _Nash_h
#define _Nash_h

/*
MELTS Source Code: RCS $Log: nash.h,v $
MELTS Source Code: RCS Revision 1.1.1.1  2006/08/15 16:57:36  ghiorso
MELTS Source Code: RCS xMELTS gcc 3.x sources
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2004/01/02 19:21:49  cvsaccount
MELTS Source Code: RCS CTserver University of Chicago
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:48:56  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
 * Revision 3.6  1997/06/21  22:49:37  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.5  1997/05/03  20:23:14  ghiorso
 * *** empty log message ***
 *
 * Revision 3.4  1997/03/27  17:03:21  ghiorso
 * *** empty log message ***
 *
 * Revision 3.3  1996/09/24  20:33:28  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 3.2  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 3.1  1995/08/18  17:51:16  ghiorso
 * MELTS Version 3 - Initial Entry
 *
 * Revision 3.1  1995/08/18  17:51:16  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Include file for macros and definitions of functions translated
**      from PASCAL algorithms published in:
**
**      Nash, J. C. (1990)
**      Compact Numerical Methods for Computers. Linear algebra and
**      function minimisation (second edition)
**      Adam Hilger, New York, 278 pages
**
**      File: NASH.H
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  August 14, 1991  Original Version
**--
*/

#ifdef MINGW
#ifndef DBL_EPSILON
#define DBL_EPSILON __DBL_EPSILON__
#endif
#endif

#define MIN1D_SUCCESS          0
#define MIN1D_BAD_INITIAL      1
#define MIN1D_ITERS_EXCEEDED   2

#define MODMRT_SUCCESS         0
#define MODMRT_BAD_INITIAL     1
#define MODMRT_ITERS_EXCEEDED  2

#define RQMCG_SUCCESS          0
#define RQMCG_ITERS_EXCEEDED   1
#define RQMCG_SINGULAR_B       2
#define RQMCG_ZERO_DETERMINANT 3

#define VMMIN_SUCCESS          0
#define VMMIN_BAD_INITIAL      1

int 
min1d(double *bb, double *st, double reltest, int *ifn, double *fnminval, 
      double (*fn1d)(double bb, int *notcomp));

int 
modmrt(int n, int m, double *Bvec, double *Fmin, double reltest, int maxIter, 
       double (*nlres)(int i, int n, double *Bvec, int *notcomp), 
       void   (*nljac)(int i, int n, double *Bvec, double *X));

int 
rqmcg(int n, double **A, double **B, double *X, int *ipr, double *rq);

int vmmin(int n, double *Bvec, double *Fmin, double intol, 
          double (*fminfn)(int n, double *Bvec, int *notcomp),
          void   (*fmingr)(int n, double *Bvec, double *g));

#endif /* _Nash_h */
