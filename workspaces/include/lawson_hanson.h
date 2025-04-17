#ifndef _Lawson_Hanson_h
#define _Lawson_Hanson_h

/*
MELTS Source Code: RCS $Log: lawson_hanson.h,v $
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
 * Revision 3.6  1997/06/21  22:49:50  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.5  1997/05/03  20:23:28  ghiorso
 * *** empty log message ***
 *
 * Revision 3.4  1997/03/27  17:03:32  ghiorso
 * *** empty log message ***
 *
 * Revision 3.3  1996/09/24  20:33:37  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 3.2  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 3.1  1995/08/18  17:45:34  ghiorso
 * MELTS Version 3 - Initial Entry
 *
 * Revision 3.1  1995/08/18  17:45:34  ghiorso
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
**      from FORTRAN algorithms published in:
**
**      Lawson, Charles L, and Hanson, Richard J. (1974)
**      Solving Least Squares Problems
**      Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 340 pp
**
**      File: LAWSON_HANSON.H
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  August 26, 1991  Original Version
**--
*/

#define HOUSEHOLDER_CALC_MODE_H1     1
#define HOUSEHOLDER_CALC_MODE_H2     2

void hfti(double **a, int m, int n, double **b, int nb, double tau, int *k,
  double *rnorm, double *h, double *g, int *p);

void householderColCol(int mode, int p, int l, int m, double **v, int vCol, 
  double *h, double **c, int cColStart, int cColEnd);

void householderColRow(int mode, int p, int l, int m, double **v, int vCol, 
  double *h, double **c, int cRowStart, int cRowEnd);

void householderRowCol(int mode, int p, int l, int m, double **v, int vRow, 
  double *h, double **c, int cColStart, int cColEnd);

void householderRowRow(int mode, int p, int l, int m, double **v, int vRow, 
  double *h, double **c, int cRowStart, int cRowEnd);

#endif /* _Lawson_Hanson_h */
