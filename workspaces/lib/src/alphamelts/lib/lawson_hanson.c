const char *lawson_hanson_ver(void) { return "$Id: lawson_hanson.c,v 1.2 2006/08/17 16:47:18 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: lawson_hanson.c,v $
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
MELTS Source Code: RCS Revision 1.2  2005/01/08 03:14:02  cvsaccount
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
 * Revision 3.6  1997/06/21  22:49:50  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.5  1997/05/03  20:23:29  ghiorso
 * *** empty log message ***
 *
 * Revision 3.4  1997/03/27  17:03:33  ghiorso
 * *** empty log message ***
 *
 * Revision 3.3  1996/09/24  20:33:37  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 3.2  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 3.1  1995/08/18  17:44:59  ghiorso
 * MELTS Version 3 - Initial Entry
 *
 * Revision 3.1  1995/08/18  17:44:59  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Source file for functions translated from algorithms published in:
**
**      Lawson, Charles L, and Hanson, Richard J. (1974)
**      Solving Least Squares Problems
**      Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 340 pp
**
**      File: LAWSON_HANSON.C
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  August 29, 1991
**              Implementation of H1, H2 and HFTI algorithms
**      V1.0-2  Mark S. Ghiorso  March 23, 1992
**              Removed incorrect calls to pow() function and replaced them
**              with a macro SQUARE
**--
*/

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "lawson_hanson.h"   /* file of external macro defn and function decl */

#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#define SQUARE(x) ((x)*(x))

/* ==========================================================================
   Algorithm HFTI (page 81-82)
   Solution of a least-squares problem with possible rank deficiency.

   The LS problem is specified by designating the matrices a and b
   which are stored as pointers to pointers to double. Note that the
   LS solution can be provided for multiple right-hand-side vectors.

        -----------n-----------    -----nb------      -----nb------
        |                     |    |           |      |           |
        |                     |    |           |      |           |
        |                     |    |           |      |           |
        |                     |    |           |      |           |
        |                     |    |           |      |           |
        |                     |    |           |      |           |
        m         a           |    n     x     |  =   m     b     |
        |                     |    |           |      |           |
        |                     |    |           |      |           |
        |                     |    |           |      |           |
        |                     |    |           |      |           |
        |                     |    |           |      |           |
        |                     |    |           |      |           |
        -----------------------    -------------      -------------

   The Method involves Householder decomposition of the columns of a, 
   accompanied by suitable permutations (P) of the elements to render the 
   resulting matrix (R) upper triangular with strictly decreasing diagonal 
   elements. The resulting transformations are in turn applied to the columns 
   of b, i.e.

                           Q a P = R   Q b  = c.

   The diagonal elements of R are then evaluated against an absolute
   tolerance parameter (tau) to determine the pseudorank (k) of the 
   resulting solution. Thus:

                            ------------             ------
                            | R11  R12 |             | c1 |
                Q a P = R = |          |   Q b = c = |    |
                            | 0    R22 |             | c2 |
                            ------------             ------

   where R11 is k by k; R12 k by n-k, R22 m-k by n-k, c1 is k by nb
   and c2 is m-k by nb. If the problem is rank deficient, Householder 
   tranformations (K) are applied to the rows of [ R11 : R12 ] to yield
   a matrix W:

                       [ R11 : R12 ] K = [ W : 0 ] .

   If the solution is full rank (k = n), W is simply R and K is the
   identity matrix. A solution (y1) is then computed (W y1 = c1) to the stable 
   subproblem, where if the system is full rank, c1 is taken to be c. The 
   final (possibly rank deficient) solution is computed by setting y2 = 0 
   (if k < n, otherwise y2 has zero length), and solving the system:

                                  ------
                                  | y1 |
                          x = P K |    | = P K y .
                                  | y2 |
                                  ------

   Note that the orthogonal transformation (K) is applied to the columns of 
   the solution vector even though it was derived from the rows of R.

   In the following algorithm, we let the solution vector(s) and the 
   right-hand-side vector(s) occupy the same storage, denoted b.

   rnrom is the residual norm (= || c2 ||), p is a vector representation
   of the permutation matrix P, and h and g are vectors storing quantities
   related to the Q and K orthogonal decompositions, respectively. The 
   remainder of the information related to the orthogonal transformations
   Q and K, is returned in the storage for a.

   ========================================================================== */

void hfti(       /* all arrays have an index base of zero                     */
  double **a,    /* i/o: matrix of coefficients of the least squares problem  */
  int m, int n,  /* inp: dimension of a, i.e. a[m][n]                         */
  double **b,    /* i/o: [m] right-hand-side vectors of the least-squares prob*/
  int    nb,     /* inp: number of vectors (columns) in the matrix b          */
  double tau,    /* inp: absolute tolerance param for determ of pseudorank    */
  int    *k,     /* out: pseudorank                                           */
  double *rnorm, /* out: [nb] residual norms for the nb solution vectors      */
  double *h,     /* out: [n] pivot scalars for Householder transformations Qi */
  double *g,     /* out: [n] pivot scalars for Householder transformations Ki */
  int    *p)     /* out: [n] interchange record                               */
{
  double hBar = 0.0, tmp;
  int i, j, l, lambda = 0;

  int mu = MIN(m, n);

  for (j=0; j<mu; j++) {
    if (j > 0) {
      for (l=j; l<n; l++) h[l] -= SQUARE(a[j-1][l]);
      for (l=j, lambda=j; l<n; l++) if (h[l] > h[lambda]) lambda = l;
    }
    if (j == 0 || (hBar + 0.001*h[lambda]) <= hBar) {
      for (l=j; l<n; l++) {
        for (i=j, h[l]=0.0; i<m; i++) h[l] += SQUARE(a[i][l]);
      }
      for (l=j, lambda=j; l<n; l++) if (h[l] > h[lambda]) lambda = l;
      hBar = h[lambda];
    } 
    p[j] = lambda;
    if (p[j] != j) {
      for (i=0; i<m; i++) {
        tmp = a[i][j]; a[i][j] = a[i][lambda]; a[i][lambda] = tmp;
      }
      h[lambda] = h[j];
    }
    householderColCol(
      HOUSEHOLDER_CALC_MODE_H1, j, j+1, m-1, a, j, &h[j], a, j+1, n-1);
    householderColCol(
      HOUSEHOLDER_CALC_MODE_H2, j, j+1, m-1, a, j, &h[j], b, 0, nb-1);
  }

  if (tau < 0.0) (*k) = (int) fabs(tau); /* special case -tau = known pseudorank */
  else {
    for (i=0, *k=0; i<n; i++) if (fabs(a[i][i]) > tau) (*k)++;
  }

  for (l=0; l<nb; l++) {
    for (i=(*k), rnorm[l]=0.0; i<m; i++) rnorm[l] += SQUARE(b[i][l]);
    rnorm[l] = sqrt(rnorm[l]);
  }

  if (*k == 0) { 
    for (i=0; i<n; i++) { for (j=0; j<nb; j++) b[i][j] = 0.0; } return; }

  if (*k < n) for (i=(*k)-1; i>=0; i--) householderRowRow(
    HOUSEHOLDER_CALC_MODE_H1, i, *k, n-1, a, i, &g[i], a, 0, i-1);

  for (l=0; l<nb; l++) {
    b[*k-1][l] /= a[*k-1][*k-1];
    for (i=(*k)-2; i>=0; i--) { 
      for (j=i+1; j<(*k); j++) b[i][l] -= a[i][j]*b[j][l];
      b[i][l] /= a[i][i];
    }  
    if (*k < n) for (i=(*k); i<n; i++) b[i][l] = 0.0;
  }
  if (*k < n) for (i=0; i<(*k); i++) householderRowCol(
    HOUSEHOLDER_CALC_MODE_H2, i, *k, n-1, a, i, &g[i], b, 0, nb-1);

  for (j=mu-1; j>=0; j--) if (p[j] != j) {
    for (l=0; l<nb; l++) {
      tmp = b[j][l]; b[j][l] = b[p[j]][l]; b[p[j]][l] = tmp;
    }
  }

}

/* ==========================================================================
   Algorithm H12 (page 57)
   Householder decomposition of a vector and application of the resulting 
     transformation to column or row vectors of a matrix.

   If mode = HOUSEHOLDER_CALC_MODE_H1
 
   Householder decomposition of a vector v into:

                              -        -
                              |  v[0]  |
                              |    .   |
                              | v[p-1] |
                              |  y[p]  |
                      Q v  =  | v[p+1] |  =  y (stored in v)
                              |    .   |
                              | v[l-1] |
                              |    0   |
                              |    .   |
                              |    0   |
                              -        -

   with y[p] = -s*(v[p]^2 + sum(l<=i<=m) v[i]^2)^1/2. On return, the pth 
   element of v is used to store s, and the storage location h is used to 
   contain the quantity v[p] - s. The vector v may be stored in a column 
   of the input matrix (in which case one uses the householderCol* routines 
   and indicates the designated column, vCol), or as a row of the input matrix 
   (in which case one uses the householderRow* routines and indicates the row, 
   vRow). The elements of v are in either case stored in a pointer to pointer 
   to double matrix. 

   If mode = HOUSEHOLDER_CALC_MODE_H1 or HOUSEHOLDER_CALC_MODE_H2

   The transformation is applied to the set of vectors stored in the matrix c. 
   These may be column vectors beginning with cColStart and ending with 
   cColEnd, for which the householder*Col routines are suitable. Alternatively, 
   the transformation may be applied to row vectors beginning with cRowStart
   and ending with cRowEnd, in which case the householder*Row routines are
   appropriate. Note that if end < start, c is ignored (in effect the
   identity transformation is performed). The information used to tranform
   c is provided in the elements of v and the scaler h.

   For both modes the user must input the pivot element p, and the indices
   of the elements of the vector v that are to be zeroed (H1 mode) or
   transformed in c (H1 mode, if end >= start, or H2 mode).

   Note that all vectors and matrices have zero-based indices and that
   on input 0 <= p < l and that l <= m. The elemnets of v indexed l through
   m are zero in H1 or transformed in H1/2.

   ========================================================================== */

void householderColCol(              /* all arrays have an index base of zero */
  int    mode, /* Calculation mode, i.e. algorithm H1 or H2                   */
  int    p,    /* index of the pivot element of the vector stored in v        */
  int    l,    /* range of indices to zero elements of vector stored in v     */
  int    m, 
  double **v,  /* matrix containing the vector v stored as v[][vCol]          */
  int    vCol, /* index number of column of v that contains vector            */
  double *h,   /* extra storage for the transformed pth element of v (u[p])   */
  double **c,  /* matrix containing the set of column vectors to apply the 
                  transformation, i.e. c[][cColStart] -> c[][cColEnd]         */
  int cColStart, int cColEnd)
{
  double b, s;
  int i, j;
  
  if (0 > p || p >= l || l > m) return;
  switch (mode) {

  case HOUSEHOLDER_CALC_MODE_H1:
    for (i=l, s=SQUARE(v[p][vCol]); i<=m; i++) s += SQUARE(v[i][vCol]);
    s = sqrt(s);
    if (v[p][vCol] > 0.0) s *= -1.0;
    *h = v[p][vCol] - s; v[p][vCol] = s;

  case HOUSEHOLDER_CALC_MODE_H2:
    b = v[p][vCol]*(*h);
    if (b == 0.0) return;

    for (j=cColStart; j<=cColEnd; j++) {
      for (i=l, s=c[p][j]*(*h); i<=m; i++) s += c[i][j]*v[i][vCol];
      s /= b;
      c[p][j] += s*(*h);
      for (i=l; i<=m; i++) c[i][j] += s*v[i][vCol];
    }

  } /* end switch */
}

void householderColRow(              /* all arrays have an index base of zero */
  int    mode, /* Calculation mode, i.e. algorithm H1 or H2                   */
  int    p,    /* index of the pivot element of the vector stored in v        */
  int    l,    /* range of indices to zero elements of vector stored in v     */
  int    m, 
  double **v,  /* matrix containing the vector v stored as v[][vCol]          */
  int    vCol, /* index number of column of v that contains vector            */
  double *h,   /* extra storage for the transformed pth element of v (u[p])   */
  double **c,  /* matrix containing the set of row vectors to apply the 
                  transformation, i.e. c[cRowStart][] -> c[cRowEnd][]         */
  int cRowStart, int cRowEnd)
{
  double b, s;
  int i, j;
  
  if (0 > p || p >= l || l > m) return;
  switch (mode) {

  case HOUSEHOLDER_CALC_MODE_H1:
    for (i=l, s=SQUARE(v[p][vCol]); i<=m; i++) s += SQUARE(v[i][vCol]);
    s = sqrt(s);
    if (v[p][vCol] > 0.0) s *= -1.0;
    *h = v[p][vCol] - s; v[p][vCol] = s;

  case HOUSEHOLDER_CALC_MODE_H2:
    b = v[p][vCol]*(*h);
    if (b == 0.0) return;

    for (j=cRowStart; j<=cRowEnd; j++) {
      for (i=l, s=c[j][p]*(*h); i<=m; i++) s += c[j][i]*v[i][vCol];
      s /= b;
      c[j][p] += s*(*h);
      for (i=l; i<=m; i++) c[j][i] += s*v[i][vCol];
    }

  } /* end switch */
}

void householderRowCol(              /* all arrays have an index base of zero */
  int    mode, /* Calculation mode, i.e. algorithm H1 or H2                   */
  int    p,    /* index of the pivot element of the vector stored in v        */
  int    l,    /* range of indices to zero elements of vector stored in v     */
  int    m, 
  double **v,  /* matrix containing the vector v stored as v[vRow][]          */
  int    vRow, /* index number of row of v that contains vector               */
  double *h,   /* extra storage for the transformed pth element of v (u[p])   */
  double **c,  /* matrix containing the set of column vectors to apply the 
                  transformation, i.e. c[][cColStart] -> c[][cColEnd]         */
  int cColStart, int cColEnd)
{
  double b, s;
  int i, j;
  
  if (0 > p || p >= l || l > m) return;
  switch (mode) {

  case HOUSEHOLDER_CALC_MODE_H1:
    for (i=l, s=SQUARE(v[vRow][p]); i<=m; i++) s += SQUARE(v[vRow][i]);
    s = sqrt(s);
    if (v[vRow][p] > 0.0) s *= -1.0;
    *h = v[vRow][p] - s; v[vRow][p] = s;

  case HOUSEHOLDER_CALC_MODE_H2:
    b = v[vRow][p]*(*h);
    if (b == 0.0) return;

    for (j=cColStart; j<=cColEnd; j++) {
      for (i=l, s=c[p][j]*(*h); i<=m; i++) s += c[i][j]*v[vRow][i];
      s /= b;
      c[p][j] += s*(*h);
      for (i=l; i<=m; i++) c[i][j] += s*v[vRow][i];
    }

  } /* end switch */
}

void householderRowRow(              /* all arrays have an index base of zero */
  int    mode, /* Calculation mode, i.e. algorithm H1 or H2                   */
  int    p,    /* index of the pivot element of the vector stored in v        */
  int    l,    /* range of indices to zero elements of vector stored in v     */
  int    m, 
  double **v,  /* matrix containing the vector v stored as v[vRow][]          */
  int    vRow, /* index number of row of v that contains vector               */
  double *h,   /* extra storage for the transformed pth element of v (u[p])   */
  double **c,  /* matrix containing the set of row vectors to apply the 
                  transformation, i.e. c[cRowStart][] -> c[cRowEnd][]         */
  int cRowStart, int cRowEnd)
{
  double b, s;
  int i, j;
  
  if (0 > p || p >= l || l > m) return;
  switch (mode) {

  case HOUSEHOLDER_CALC_MODE_H1:
    for (i=l, s=SQUARE(v[vRow][p]); i<=m; i++) s += SQUARE(v[vRow][i]);
    s = sqrt(s);
    if (v[vRow][p] > 0.0) s *= -1.0;
    *h = v[vRow][p] - s; v[vRow][p] = s;

  case HOUSEHOLDER_CALC_MODE_H2:
    b = v[vRow][p]*(*h);
    if (b == 0.0) return;

    for (j=cRowStart; j<=cRowEnd; j++) {
      for (i=l, s=c[j][p]*(*h); i<=m; i++) s += c[j][i]*v[vRow][i];
      s /= b;
      c[j][p] += s*(*h);
      for (i=l; i<=m; i++) c[j][i] += s*v[vRow][i];
    }

  } /* end switch */
}

/* end of file LAWSON_HANSON.C */
