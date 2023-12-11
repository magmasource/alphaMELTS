const char *nash_ver(void) { return "$Id: nash.c,v 1.2 2006/08/17 16:47:19 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: nash.c,v $
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
MELTS Source Code: RCS Revision 1.1.1.1  2004/01/02 19:21:49  cvsaccount
MELTS Source Code: RCS CTserver University of Chicago
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2003/05/03 18:43:56  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:58:12  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
 * Revision 3.7  1997/06/21  22:49:38  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.6  1997/05/03  20:23:14  ghiorso
 * *** empty log message ***
 *
 * Revision 3.5  1997/03/27  17:03:22  ghiorso
 * *** empty log message ***
 *
 * Revision 3.4  1996/09/24  20:33:28  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 3.3  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 3.2  1995/11/01  22:40:27  ghiorso
 * Implementation of subsolidus options after Asimow.
 * Additional implementation of nepheline solid solutions.
 *
 * Revision 3.2  1995/11/01  22:40:27  ghiorso
 * Implementation of subsolidus options after Asimow.
 * Additional implementation of nepheline solid solutions.
 *
 * Revision 3.1  1995/08/18  17:50:56  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Source file for functions translated from PASCAL algorithms
**      published in:
**
**      Nash, J. C. (1990)
**      Compact Numerical Methods for Computers. Linear algebra and
**      function minimisation (second edition)
**      Adam Hilger, New York, 278 pages
**
**      File: NASH.C
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  August 14, 1991  Original Version
**              function modmrt (C version of Marquardt routine for
**                               nonlinear least squares)
**      V1.1-1  Mark S. Ghiorso  August 26, 1991
**              function min1d  (C version of minimisation routine
**                               for a function of one variable)
**      V1.2-1  Mark S. Ghiorso  October 1, 1991
**              function rqmcg  (C version of Rayleigh quotient minimisation
**                               by conjugate gradients)
**      V1.3-1  Mark S. Ghiorso (C version of Variable metric minimiser)
**      V1.3-2  Mark S. Ghiorso February 12, 1992
**              (1) Added definition of BIG to min1d in order to prevent
**                  floating point overflow during parabolic searches
**              (2) Altered division algorithm in parabalic searches in order
**                  to prevent floating point overflow
**      V1.3-3  Mark S. Ghiorso March 23, 1992
**              Replaced illegal calls to pow() function with references to
**              SQUARE macro
**      V1.4-1  Mark S. Ghiorso March 23, 1995
**              Added trivial case return to modmrt
**--
*/

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "nash.h"   /* file of external macro definitions and function decl */

#define SQUARE(x) ((x)*(x))

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE  1
#endif

/******************************************************************************
 * Algorithm 17.
 * One-dimensional minimisation of a function using success-failure search and
 * parabolic inverse interpolation
 ******************************************************************************/

#define A1  1.5
#define A2 -0.25

int min1d(          /* returned value, MODE flag as defined in NASH.H         */
    double *bb,       /* initial guess to minimum, resulting minimum position   */
    double *st,       /* initial and final step-length                          */
    double reltest,   /* maximal |error| allowed in computing value of bb       */
    int    *ifn,      /* input: maximum No. of func eval; output: actual No.    */
    double *fnminval, /* minimum function value on return                       */
    double (*fn1d)(double bb, int *notcomp))   /* computes function value at bb */
{
    double fii, s0, s1, s2, tt0, tt1, tt2, x0, x1, x2, xii;
    int notcomp, tripleok;
    int ifnMax = *ifn;
    double BIG = sqrt(DBL_MAX);

    *ifn = 0;
    x1   = *bb;
    s0 = (*fn1d)(x1, &notcomp); if (notcomp) return MIN1D_BAD_INITIAL; (*ifn)++;

    do {
        x0 = x1;
        *bb = x0;
        x1 = x0 + (*st);
        s1 = (*fn1d)(x1, &notcomp); if (notcomp) s1 = BIG; (*ifn)++;
            if (*ifn > ifnMax) return MIN1D_ITERS_EXCEEDED;
        tripleok = FALSE;
        if (s1 < s0) {
            do {
                *st *= A1;
                x2 = x1 + (*st);
                s2 = (*fn1d)(x2, &notcomp); if (notcomp) s2 = BIG; (*ifn)++;
                    if (*ifn > ifnMax) return MIN1D_ITERS_EXCEEDED;
                if (s2 < s1) {
                    s0 = s1; s1 = s2;
                    x0 = x1; x1 = x2;
                } else {
                    tripleok = TRUE;
                }
            } while (!tripleok);
        } else {
            *st *= A2;
            tt2 = s0; s0 = s1; s1 = tt2;
            tt2 = x0; x0 = x1; x1 = tt2;
            do {
                x2 = x1 + (*st);
                s2 = (*fn1d)(x2, &notcomp); if (notcomp) s2 = BIG; (*ifn)++;
                    if (*ifn > ifnMax) return MIN1D_ITERS_EXCEEDED;
                if (s2 < s1) {
                    s0 = s1; s1 = s2; x0 = x1; x1 = x2;
                    *st *= A1;
                } else {
                    tripleok = TRUE;
                }
            } while (!tripleok);
        }
        tt0 = x0 - x1;
        tt1 = (s0 - s1)*(*st); tt2 = (s2 - s1)*tt0;
        if (tt1 != tt2) {
            /* was: *st = 0.5*(tt2*tt0 - tt1*(*st))/(tt2 - tt1); */
            *st = 0.5*( (tt2/(tt2 - tt1))*tt0 - (tt1/(tt2 - tt1))*(*st) );
            xii = x1 + (*st);
            if (fabs(xii-x1) > reltest) {
                fii = (*fn1d)(xii, &notcomp); if (notcomp) fii = BIG; (*ifn)++;
                    if (*ifn > ifnMax) return MIN1D_ITERS_EXCEEDED;
                if (fii < s1) {
                    s1 = fii; x1 = xii;
                }
            }
        }
        s0 = s1;
    } while (fabs(*bb - x1) > reltest);

    *fnminval = s1;
    (void) (*fn1d)(*bb, &notcomp); /* resets constraints->liquidDelta to correct value */
    return MIN1D_SUCCESS;
}

#undef A1
#undef A2

/******************************************************************************
 * Algorithm 21.
 * Variable metric minimiser
 ******************************************************************************/

#define STEPREDN  0.2
#define ACCTOL    0.0001

int vmmin(         /* returned value, MODE flag as set in NASH.H              */
    int n,           /* number of parameters in function to be minimized        */
    double *Bvec,    /* input: initial guess; output: best estimate of solution */
    double *Fmin,    /* function value at the minimum                           */
    double reltest,  /* user-initialized convergence tolerance for linear search*/
    double (*fminfn)(int n, double *Bvec, int *notcomp), /* function value      */
    void   (*fmingr)(int n, double *Bvec, double *g))    /* gradients           */
{
    double **B, *c, D1, D2, f, *g, gradproj, s, steplength, *t, *X;
    int accpoint, count, funcount, gradcount, i, ilast, j, notcomp, iterations;

    f = fminfn(n, Bvec, &notcomp); if (notcomp) return VMMIN_BAD_INITIAL;

    B = (double **) malloc((unsigned) n*sizeof(double *));
    for (i=0; i<n; i++) B[i] = (double *) malloc((unsigned) n*sizeof(double));
    c = (double *) malloc((unsigned) n*sizeof(double));
    g = (double *) malloc((unsigned) n*sizeof(double));
    t = (double *) malloc((unsigned) n*sizeof(double));
    X = (double *) malloc((unsigned) n*sizeof(double));

    *Fmin = f; funcount = 1; gradcount = 1;
    fmingr(n, Bvec, g); ilast = gradcount;
    iterations = 0;

    do {
            iterations++; if (iterations > 10000.0) {
                    for (i=0; i<n; i++) free(B[i]);
                    free(B); free(c); free(g); free(t); free(X);
                    return VMMIN_BAD_INITIAL;
            }
        if (ilast == gradcount)
            for (i=0; i<n; i++) { for (j=0; j<n; j++) B[i][j] = 0.0; B[i][i] = 1.0; }
        for (i=0; i<n; i++) { X[i] = Bvec[i]; c[i] = g[i]; }
        for (i=0, gradproj=0.0; i<n; i++) {
            for (j=0, s=0.0; j<n; j++) s -= B[i][j]*g[j];
            t[i] = s; gradproj += s*g[i];
        }
        if (gradproj < 0.0) {
            steplength = 1.0; accpoint = FALSE;
            do {
                for (i=0, count=0; i<n; i++) {
                    Bvec[i] = X[i] + steplength*t[i];
                    if (fabs(X[i]-Bvec[i]) < reltest) count++;
                }
                if (count < n) {
                    f = fminfn(n, Bvec, &notcomp); funcount++;
                    accpoint = (!notcomp) && (f <= *Fmin+gradproj*steplength*ACCTOL);
                    if (!accpoint) steplength *= STEPREDN;
                }
            } while (count != n && !accpoint);
            if (count < n) {
                *Fmin = f;
                fmingr(n, Bvec, g); gradcount++;
                for (i=0, D1=0.0; i<n; i++)
                    { t[i] *= steplength; c[i] = g[i] - c[i]; D1 += t[i]*c[i]; }
                if (D1 > 0.0) {
                    for (i=0, D2=0.0; i<n; i++) {
                        for (j=0, s=0.0; j<n; j++) s += B[i][j]*c[j];
                        X[i] = s; D2 += s*c[i];
                    }
                    D2 = 1.0 + D2/D1;
                    for (i=0; i<n; i++) for (j=0; j<n; j++)
                        B[i][j] -= (t[i]*X[j] + X[i]*t[j] - D2*t[i]*t[j])/D1;
                } else ilast = gradcount;
            } else {
                if (ilast < gradcount) { count = 0; ilast = gradcount; }
            }
        } else { count = 0; ilast = gradcount; }
    } while (count != n || ilast != gradcount);

    for (i=0; i<n; i++) free(B[i]);
    free(B); free(c); free(g); free(t); free(X);
    return VMMIN_SUCCESS;
}

#undef STEPREDN
#undef ACCTOL

/******************************************************************************
 * Algorithm 23.
 * Modified Marquardt method for minimising a nonlinear sum-of-squares function
 ******************************************************************************/

#define DEC  0.4
#define INC 10.0
#define PHI  1.0

static void Choldcmp(int n, double *a, int *singmat);
static void Cholback(int n, double *a, double *delta);

int modmrt(        /* returned value, MODE flag as defined in NASH.H          */
    int n,           /* number of parameters in nonlinear least-squares problem */
    int m,           /* number of residuals in nonlinear least-squares problem  */
    double *Bvec,    /* vector of parameters in nonlinear least-squares problem */
    double *Fmin,    /* returned value, minimum value of sum-of-squares function*/
    double reltest,  /* maximal |error| allowed in computing elements of Bvec   */
    int    maxIter,  /* maximum number of times the func (*nlres) may be called */
    double (*nlres)(int i, int n, double *Bvec, int *notcomp), /* residual func */
    void   (*nljac)(int i, int n, double *Bvec, double *X))    /* gradient func */
{
    /* The residual function returns the residual of the ith equation and a
       flag (notcomp) which is set to true if the residual cannot be computed.
     The gradient function returns a vector corresponding to the ith row
       of the Jacobian of the system of nonlinear equations                   */

    double *a, *c, *delta, *res, *v, *X;
    double lambda, p;
    int count, i, ifn, igrad, j, k, nn2, q;
    int notcomp, singmat, calcmat;

    *Fmin   = DBL_MAX;                                               /* STEP  0 */
    lambda  = 0.0001;
    ifn     = 0;
    igrad   = 0;
    calcmat = TRUE;
    nn2     = n*(n+1)/2;

    if (n == 0) return MODMRT_SUCCESS;

    a     = (double *) malloc((unsigned) nn2*sizeof(double));
    c     = (double *) malloc((unsigned) nn2*sizeof(double));
    delta = (double *) malloc((unsigned) n*sizeof(double));
    res   = (double *) malloc((unsigned) m*sizeof(double));
    v     = (double *) malloc((unsigned) n*sizeof(double));
    X     = (double *) malloc((unsigned) n*sizeof(double));

    for (i=0, p=0.0; i<m; i++) {                                     /* STEP  1 */
        res[i] = (*nlres)(i, n, Bvec, &notcomp);
        if (notcomp) {
            free(a); free(c); free(delta); free(res); free(v); free(X);
            return MODMRT_BAD_INITIAL;
        }
        p += SQUARE(res[i]);
    }
    ifn++; *Fmin = p; count = 0;

    while (count < n) {                                              /* STEP  2 */
        if (calcmat) {
            igrad++;                                                     /* STEP  3 */
            for (j=0; j<nn2; j++) a[j] = 0.0;
            for (j=0; j<n;   j++) v[j] = 0.0;
            for (i=0; i<m; i++) {                                        /* STEP  4 */
                (*nljac)(i, n, Bvec, X);
                for (j=1; j<=n; j++) {
                    v[j-1] = v[j-1] + X[j-1]*res[i];
                    q = j*(j-1)/2;
                    for (k=1; k<=j; k++) a[q+k-1] += X[j-1]*X[k-1];
                }
            }
            for (j=0; j<nn2; j++) c[j] = a[j];                           /* STEP  5 */
            for (j=0; j<n;   j++) X[j] = Bvec[j];
        }
        for (j=1; j<=n; j++) {                                          /* STEP  6 */
            q = j*(j+1)/2;
            a[q-1] = c[q-1]*(1.0+lambda) + PHI*lambda;
            delta[j-1] = -v[j-1];
            if (j > 1) for (i=1; i<=(j-1); i++) a[q-i-1] = c[q-i-1];
        }
        notcomp = FALSE;
        Choldcmp(n, a, &singmat);                                      /* STEP  7 */
        if (!singmat) {
            Cholback(n, a, delta);                                       /* STEP  8 */
            for (i=0, count = 0; i<n; i++) {                             /* STEP  9 */
                Bvec[i] = X[i] + delta[i];
                if (fabs(Bvec[i]-X[i]) <= reltest) count++;
            }
            if (count < n) {                                             /* STEP 10 */
                p = 0.0; i = 0;
                do {
                    res[i] = (*nlres)(i, n, Bvec, &notcomp);
                    if (!notcomp) p += SQUARE(res[i]);
                    i++;
                } while (!notcomp && i < m); ifn++;
                if (ifn > maxIter) {
                    free(a); free(c); free(delta); free(res); free(v); free(X);
                    return MODMRT_ITERS_EXCEEDED;
                }
            }
        }
        if (count < n) {
            if (!singmat && !notcomp && p < *Fmin) {
                lambda *= DEC;
                *Fmin = p;
                calcmat = TRUE;
            } else {
                lambda *= INC;
                if (lambda < SQUARE(DBL_EPSILON)) lambda = DBL_EPSILON;
                calcmat = FALSE;
            }
        }
    }

    free(a); free(c); free(delta); free(res); free(v); free(X);
    return MODMRT_SUCCESS;
}

#undef DEC
#undef INC
#undef PHI

/******************************************************************************
 * Algorithm 7.
 * Choleski decomposition in compact storage
 ******************************************************************************/

static void Choldcmp(int n, double *a, int *singmat)
{
    int i, j, k, m, q;
    double s;

    *singmat = FALSE;
    for (j=1; j<=n; j++) {                                           /* STEP  1 */
        q = j*(j+1)/2;                                                 /* STEP  2 */
        if (j > 1) {                                                   /* STEP  3 */
            for (i=j; i<=n; i++) {                                       /* STEP  4 */
                m = i*(i-1)/2 + j; s = a[m-1];
                for (k=1; k<=(j-1); k++) s -= a[m-k-1]*a[q-k-1];
                a[m-1] = s;
            }
        }
        if (a[q-1] <= 0.0) {                                           /* STEP  5 */
            *singmat = TRUE;
            a[q-1] = 0.0;
        }
        s = sqrt(a[q-1]);                                               /* STEP 7 */
        for (i=j; i<=n; i++) {                                          /* STEP 8 */
            m = i*(i-1)/2 + j;
            a[m-1] = (s == 0.0) ? 0.0 : a[m-1]/s;
        }
    }
}

/******************************************************************************
 * Algorithm 8.
 * Choleski back-substitution
 ******************************************************************************/

static void Cholback(int n, double *a, double *x)
{
    int i, j, q;

    x[0] = (a[0] == 0.0) ? 0.0 : x[0]/a[0];                          /* STEP  1 */
    if (n > 1) {                                                     /* STEP  2 */
        q = 1;                                                         /* STEP  3 */
        for (i=2; i<=n; i++) {                                         /* STEP  4 */
            for (j=1; j<=(i-1); j++) {q++; x[i-1] -= a[q-1]*x[j-1];}     /* STEP  5 */
            q++;                                                         /* STEP  6 */
            x[i-1] = (a[q-1] == 0.0) ? 0.0 : x[i-1]/a[q-1];              /* STEP  7 */
        }                                                              /* STEP  8 */
    }
                                                                   /* STEP  9 */
    x[n-1] = (a[n*(n+1)/2-1] == 0.0) ? 0.0 : x[n-1]/a[n*(n+1)/2-1];
    if (n > 1) {                                                     /* STEP 10 */
        for (i=n; i>=2; i--) {                                     /* STEPS 11/12 */
            q = i*(i-1)/2;
            for (j=1; j<=(i-1); j++) x[j-1] -= x[i-1]*a[q+j-1];          /* STEP 13 */
            x[i-2] = (a[q-1] == 0.0) ? 0.0 : x[i-2]/a[q-1];              /* STEP 14 */
        }                                                              /* STEP 15 */
    }                                                                /* STEP 16 */
}

/******************************************************************************
 * Algorithm 25.
 * Rayleigh quotient minimisation by conjugate gradients
 *   The quotient is defined by x(T)Ax / x(T)Bx, where A and B are square
 *   matrices and B is positive definite.
 ******************************************************************************/

int rqmcg(          /* returned value, MODE flag as defined in NASH.H         */
    int    n,         /* input:  order of the square matrices A and B           */
    double **A,       /* input:  square matrix in numerator of the quotient     */
    double **B,       /* input:  square matrix in denominator of the quotient   */
    double *X,        /* i/o:    approximation to minimal eigenvector           */
    int    *ipr,      /* i/o:    limit on number of iterations/ actual number   */
    double *rq)       /* output: Rayleigh quotient at minimum (eigenvalue)      */
{
    int conv, count, fail, i, itn, itlimit, j;
    double *avec, *bvec, *yvec, *zvec, *g, *t;
    double beta, d, gg, pa, pn, step, ta, tabt, tat,
         tbt, tol, u, v, w, xat, xax, xbt, xbx;

    avec = (double *) malloc((unsigned) n*sizeof(double));
    bvec = (double *) malloc((unsigned) n*sizeof(double));
    yvec = (double *) malloc((unsigned) n*sizeof(double));
    zvec = (double *) malloc((unsigned) n*sizeof(double));
    g    = (double *) malloc((unsigned) n*sizeof(double));
    t    = (double *) malloc((unsigned) n*sizeof(double));

    itlimit = *ipr; fail = FALSE; conv = FALSE;
    *ipr = 0; tol = SQUARE(n)*DBL_EPSILON; pa = DBL_MAX;

    while (*ipr <= itlimit && !conv) {                               /* STEP  1 */
        for (i=0; i<n; i++) for (j=0, avec[i]=0.0, bvec[i]=0.0; j<n; j++)
            { avec[i] += A[i][j]*X[j]; bvec[i] += B[i][j]*X[j]; }
        (*ipr)++;

        for (i=0, xax=0.0, xbx=0.0; i<n; i++)                          /* STEP  2 */
            { xax += X[i]*avec[i]; xbx += X[i]*bvec[i]; }
        if (xbx <= tol) {                                              /* STEP  3 */
            free(avec); free(bvec); free(yvec); free(zvec); free(g); free(t);
            return RQMCG_SINGULAR_B;
        }
        *rq = xax/xbx;                                                 /* STEP  4 */

        if (*rq < pa) {                                                /* STEP  5 */
            pa = *rq;
            for (i=0, gg=0.0; i<n; i++)
                { g[i] = 2.0*(avec[i] - (*rq)*bvec[i])/xbx; gg += g[i]*g[i]; }
            if (gg > tol) {                                              /* STEP  7 */
                for (i=0; i<n; i++) t[i] = -g[i];                          /* STEP  8 */
                itn = 0;                                                   /* STEP  9 */
                do {
                    itn++;                                                   /* STEP 10 */
                    for (i=0; i<n; i++) for (j=0, yvec[i]=0.0, zvec[i]=0.0; j<n; j++)
                        { yvec[i] += A[i][j]*t[j]; zvec[i] += B[i][j]*t[j]; }
                    (*ipr)++;
                    for (i=0, tat = 0.0, tbt = 0.0, xat = 0.0, xbt = 0.0; i<n; i++) {
                        xat += X[i]*yvec[i]; tat += t[i]*yvec[i];
                        xbt += X[i]*zvec[i]; tbt += t[i]*zvec[i];
                    }
                    u = tat*xbt - xat*tbt; v = tat*xbx - xax*tbt;            /* STEP 12 */
                    w = xat*xbx - xax*xbt; d = v*v - 4.0*u*w;
                    if (d < 0.0) {                                           /* STEP 13 */
                        free(avec); free(bvec); free(yvec); free(zvec); free(g); free(t);
                        return RQMCG_ZERO_DETERMINANT;
                    }
                    d = sqrt(d);                                             /* STEP 14 */
                    step = (v > 0.0) ? -2.0*w/(v+d) : 0.5*(d-v)/u;
                    for (i=0, count=0, xax=0.0, xbx=0.0; i<n; i++) {         /* STEP 15 */
                        avec[i] += step*yvec[i]; bvec[i] += step*zvec[i];
                        w = X[i]; X[i] = w + step*t[i];
                        if (fabs(X[i]-w) < 10.0*DBL_EPSILON) count++;
                        xax += X[i]*avec[i];  xbx += X[i]*bvec[i];
                    }
                    if (xbx <= tol) {                                        /* STEP 16 */
                        free(avec); free(bvec); free(yvec); free(zvec); free(g); free(t);
                        return RQMCG_SINGULAR_B;
                    } else pn = xax/xbx;
                    if (count < n && pn < *rq) {                     /* STEPS 17 and 18 */
                        *rq = pn;                                              /* STEP 19 */
                        for (i=0, gg=0.0; i<n; i++)
                            { g[i] = 2.0*(avec[i]-pn*bvec[i])/xbx; gg += g[i]*g[i]; }
                        if (gg > tol) {                                        /* STEP 20 */
                            for (i=0, xbt=0.0; i<n; i++) xbt += X[i]*zvec[i];    /* STEP 21 */
                            for (i=0, tabt=0.0, beta=0.0; i<n; i++) {            /* STEP 22 */
                                w = yvec[i] - pn*zvec[i]; tabt += t[i]*w;
                                beta += g[i]*(w - g[i]*xbt);
                            }
                            beta /= tabt;                                        /* STEP 23 */
                            for (i=0; i<n; i++) t[i] = beta*t[i] - g[i];
                        }
                    } else {
                        if (itn == 1) conv = TRUE;
                        itn = n + 1;
                    }
                } while (itn < n && count != n && gg > tol && !conv);
            } else conv = TRUE;
        } else conv = TRUE;
        for (i=0, ta=0.0; i<n; i++) ta += SQUARE(X[i]);
        ta = 1.0/sqrt(ta);
        for (i=0; i<n; i++) X[i] *= ta;
    }

    free(avec); free(bvec); free(yvec); free(zvec); free(g); free(t);
    if (*ipr > itlimit) return RQMCG_ITERS_EXCEEDED;
    else                return RQMCG_SUCCESS;

}

#undef FALSE
#undef TRUE

/* end of file NASH.C */
