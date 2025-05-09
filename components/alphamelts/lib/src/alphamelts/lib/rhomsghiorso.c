const char *rhomsghiorso_ver(void) { return "$Id: rhomsghiorso.c,v 1.4 2008/03/06 17:51:23 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: rhomsghiorso.c,v $
MELTS Source Code: RCS Revision 1.4  2008/03/06 17:51:23  ghiorso
MELTS Source Code: RCS New fluid fractionation mode and other enhancements.
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
MELTS Source Code: RCS Revision 1.4  2003/09/27 15:35:22  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.3  2003/06/24 16:42:42  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2003/05/08 18:07:30  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1  2003/05/07 04:39:09  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:58:12  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Routines to compute rhombohedral oxide solution properties
**      (file: RHOMBOHEDRAL.C)
**
**  MODIFICATION HISTORY:
**
**  v 1.0 New model May 2003
**
**--
*/

#ifdef DEBUG
#undef DEBUG
#endif

#include "melts_gsl.h"
#include "silmin.h"  /* Structure definitions foor SILMIN package */

#define SQUARE(x)  ((x)*(x))
#define CUBE(x)    ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

/*
 *=============================================================================
 * Rhombohedral solution parameters:
 * Ghiorso, M.S., Evans, B.W. (2008)
 */

static double dvilm     =  0.010758; /* joules/bar */
static double dvgei     =  0.010758; /* joules/bar */
static double dvpyr     =  0.010758; /* joules/bar */

static double wvilm     =  0.035089; /* joules/bar */
static double wvgei     =  0.035089; /* joules/bar */
static double wvpyr     =  0.035089; /* joules/bar */

static double dwvhmilm  =  0.013701; /* joules/bar */
static double dwvhmgei  =  0.013701; /* joules/bar */
static double dwvhmpyr  =  0.013701; /* joules/bar */

static double wvhmilm2  =  0.000000; /* joules/bar */
static double wvhmgei2  =  0.000000; /* joules/bar */
static double wvhmpyr2  =  0.000000; /* joules/bar */

static double dwvcrnilm =  0.013701; /* joules/bar - guess */
static double dwvcrngei =  0.013701; /* joules/bar - guess */
static double dwvcrnpyr =  0.013701; /* joules/bar - guess */

static double wvhmilm   = -0.11764;  /* joules/bar */
static double wvhmgei   = -0.11764;  /* joules/bar - guess */
static double wvhmpyr   = -0.11764;  /* joules/bar - guess */
static double wvilmgei  =  0.000000; /* joules/bar - guess */
static double wvilmpyr  =  0.000000; /* joules/bar - guess */
static double wvgeipyr  =  0.000000; /* joules/bar - guess */

static double dhilm     =  17477.0; /* joules ordering analysis     	- fixed    */
static double dhgei     =  17477.0; /* joules ordering analysis     	- fixed    */
static double dhpyr     =  17477.0; /* joules ordering analysis     	- fixed    */

static double whilm     =   3189.0; /* joules ordering analysis     	- fixed    */
static double whgei     =   3189.0; /* joules ordering analysis     	- fixed    */
static double whpyr     =   3189.0; /* joules ordering analysis     	- fixed    */

static double dwhhmilm  =  -5626.63; /* -1510.8 joules ordering analysis - dwhhmilm/2+whhmilm2/4 is fixed at -3021.6 */
static double dwhhmgei  =  -5626.63; /* -1510.8 joules ordering analysis - dwhhmgei/2+whhmgei2/4 is fixed at -3021.6 */
static double dwhhmpyr  =  -5626.63; /* -1510.8 joules ordering analysis - dwhhmpyr/2+whhmpyr2/4 is fixed at -3021.6 */

static double whhmilm2  =   -833.14; /* -9064.8 joules ordering analysis - dwhhmilm/2+whhmilm2/4 is fixed at -3021.6 */
static double whhmgei2  =   -833.14; /* -9064.8 joules ordering analysis - dwhhmgei/2+whhmgei2/4 is fixed at -3021.6 */
static double whhmpyr2  =   -833.14; /* -9064.8 joules ordering analysis - dwhhmpyr/2+whhmpyr2/4 is fixed at -3021.6 */

static double dwhcrnilm =      0.0; /* joules - guess                   - fixed    */
static double dwhcrngei =      0.0; /* joules - guess                   - fixed    */
static double dwhcrnpyr =      0.0; /* joules - guess                   - fixed    */

static double whmcrn    =  69000.0; /* 69000.0 joules Majzlan et al., 2002      - fixed   */
static double whhmilm   =  22535.6; /* 21200.0 joules solvus analysis	       - variable */
static double whhmgei   =  22535.6; /* 21200.0 joules guess		       - variable */
static double whhmpyr   =  22535.6; /* 21200.0 joules guess		       - variable */
static double wcrnilm   =  22535.6; /* 21200.0 joules guess		       - variable */
static double wcrngei   =  22535.6; /* 21200.0 joules guess		       - variable */
static double wcrnpyr   =  22535.6; /* 21200.0 joules guess		       - variable */
static double whilmgei  =   2600.0; /*  2600.0 joules Pownceby and Fisher-White - fixed   */
static double whilmpyr  =   2200.0; /*  2200.0 joules O'Neill et al.	       - fixed    */
static double whgeipyr  =   2600.0; /*  2600.0 joules analogy with ilm-gei      - fixed   */

static double whilmgeiT =  88099.9; /* 63694.0 guess, equal to whilmgei zeroes (G*)st	  */
static double whilmpyrT =  30244.0; /* 45214.0 guess, equal to whilmpyr zeroes (G*)su	  */
static double whgeipyrT =   2600.0; /* 54454.0 guess, equal to whilmpyr zeroes (G*)tu	  */

static double whilmilmgei =    0.0; /* joules - guess */
static double whilmgeigei =    0.0; /* joules - guess */
static double whilmilmpyr =    0.0; /* joules - guess */
static double whilmpyrpyr =    0.0; /* joules - guess */
static double whgeigeipyr =    0.0; /* joules - guess */
static double whgeipyrpyr =    0.0; /* joules - guess */

static double SROconst = 0.0730205; /* Scaling factor for SRO */
static double SRO600   = 0.0730205; /* Scaling factor for SRO  600 C */
static double SRO700   = 0.0730205; /* Scaling factor for SRO  700 C */
static double SRO800   = 0.0730205; /* Scaling factor for SRO  800 C */
static double SRO900   = 0.0730205; /* Scaling factor for SRO  900 C */
static double SRO1000  = 0.0730205; /* Scaling factor for SRO 1000 C */
static double SRO1100  = 0.0730205; /* Scaling factor for SRO 1100 C */
static double SRO1200  = 0.0730205; /* Scaling factor for SRO 1200 C */
static double SRO1300  = 0.0730205; /* Scaling factor for SRO 1300 C */

static int makeSpline = TRUE;

/* Short-range Order correction for entropy */

void printParameters() {
    printf("rhomsghiorso: dhilm	   : %g + %g P\n", dhilm      , dvilm	 );
    printf("rhomsghiorso: dhgei	   : %g + %g P\n", dhgei      , dvgei	 );
    printf("rhomsghiorso: dhpyr	   : %g + %g P\n", dhpyr      , dvpyr	 );

    printf("rhomsghiorso: whilm	   : %g + %g P\n", whilm      , wvilm	 );
    printf("rhomsghiorso: whgei	   : %g + %g P\n", whgei      , wvgei	 );
    printf("rhomsghiorso: whpyr	   : %g + %g P\n", whpyr      , wvpyr	 );

    printf("rhomsghiorso: dwhhmilm   : %g + %g P\n", dwhhmilm   , dwvhmilm );
    printf("rhomsghiorso: dwhhmgei   : %g + %g P\n", dwhhmgei   , dwvhmgei );
    printf("rhomsghiorso: dwhhmpyr   : %g + %g P\n", dwhhmpyr   , dwvhmpyr );

    printf("rhomsghiorso: whhmilm2   : %g + %g P\n", whhmilm2   , wvhmilm2 );
    printf("rhomsghiorso: whhmgei2   : %g + %g P\n", whhmgei2   , wvhmgei2 );
    printf("rhomsghiorso: whhmpyr2   : %g + %g P\n", whhmpyr2   , wvhmpyr2 );

    printf("rhomsghiorso: dwhcrnilm  : %g + %g P\n", dwhcrnilm  , dwvcrnilm);
    printf("rhomsghiorso: dwhcrngei  : %g + %g P\n", dwhcrngei  , dwvcrngei);
    printf("rhomsghiorso: dwhcrnpyr  : %g + %g P\n", dwhcrnpyr  , dwvcrnpyr);

    printf("rhomsghiorso: whmcrn     : %g\n",        whmcrn                );
    printf("rhomsghiorso: whhmilm    : %g + %g P\n", whhmilm    , wvhmilm  );
    printf("rhomsghiorso: whhmgei    : %g + %g P\n", whhmgei    , wvhmgei  );
    printf("rhomsghiorso: whhmpyr    : %g + %g P\n", whhmpyr    , wvhmpyr  );
    printf("rhomsghiorso: wcrnilm    : %g\n",        wcrnilm      	 );
    printf("rhomsghiorso: wcrngei    : %g\n",        wcrngei      	 );
    printf("rhomsghiorso: wcrnpyr    : %g\n",        wcrnpyr      	 );
    printf("rhomsghiorso: whilmgei   : %g + %g P\n", whilmgei   , wvilmgei );
    printf("rhomsghiorso: whilmpyr   : %g + %g P\n", whilmpyr   , wvilmpyr );
    printf("rhomsghiorso: whgeipyr   : %g + %g P\n", whgeipyr   , wvgeipyr );

    printf("rhomsghiorso: whilmgeiT  : %g\n",  	   whilmgeiT		 );
    printf("rhomsghiorso: whilmpyrT  : %g\n",  	   whilmpyrT		 );
    printf("rhomsghiorso: whgeipyrT  : %g\n",  	   whgeipyrT		 );

    printf("rhomsghiorso: whilmilmgei: %g\n",  	   whilmilmgei  	 );
    printf("rhomsghiorso: whilmgeigei: %g\n",  	   whilmgeigei  	 );
    printf("rhomsghiorso: whilmilmpyr: %g\n",  	   whilmilmpyr  	 );
    printf("rhomsghiorso: whilmpyrpyr: %g\n",  	   whilmpyrpyr  	 );
    printf("rhomsghiorso: whgeigeipyr: %g\n",  	   whgeigeipyr  	 );
    printf("rhomsghiorso: whgeipyrpyr: %g\n",  	   whgeipyrpyr  	 );

    printf("rhomsghiorso: SROconst   : %g\n",  	   SROconst		 );
    printf("rhomsghiorso: SRO600     : %g\n",  	   SRO600		 );
    printf("rhomsghiorso: SRO700     : %g\n",  	   SRO700		 );
    printf("rhomsghiorso: SRO800     : %g\n",  	   SRO800		 );
    printf("rhomsghiorso: SRO900     : %g\n",  	   SRO900		 );
    printf("rhomsghiorso: SRO1000    : %g\n",  	   SRO1000		 );
    printf("rhomsghiorso: SRO1100    : %g\n",  	   SRO1100		 );
    printf("rhomsghiorso: SRO1200    : %g\n",  	   SRO1200		 );
    printf("rhomsghiorso: SRO1300    : %g\n",  	   SRO1300		 );
}

static double fSRO(double tk, int order) {
    static double x[8], y[8], y2[8];
    static double tkOld = -9999.0;
    static double value = 0.0, dvalue = 0.0, d2value = 0.0, d3value = 0.0;
    int i, n=8;

    if (makeSpline) {
        double u[7];
        x[0] =  873.15;  y[0] = SRO600;
        x[1] =  973.15;  y[1] = SRO700;
        x[2] = 1073.15;  y[2] = SRO800;
        x[3] = 1173.15;  y[3] = SRO900;
        x[4] = 1273.15;  y[4] = SRO1000;
        x[5] = 1373.15;  y[5] = SRO1100;
        x[6] = 1473.15;  y[6] = SRO1200;
        x[7] = 1573.15;  y[7] = SRO1300;

        y2[0] = 0.0;
        u[0]  = 0.0;
        for (i=1; i<(n-1); i++) {
            double sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
            double p = sig*y2[i-1]+2.0;
            y2[i] = (sig-1.0)/p;
            u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
            u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
        }
        y2[n-1] = 0.0;
        for (i=(n-2); i>=0; i--) y2[i] = y2[i]*y2[i+1] + u[i];

        tkOld = -9999.0;
        makeSpline = FALSE;
    }

    if (tk != tkOld) {
        int k, klo = 0, khi = n-1;
        double h, b, a;
        while ((khi-klo) > 1) {
            k = (khi+klo) >> 1;
            if (x[k] > tk) khi = k;
            else klo = k;
        }
        h = x[khi] - x[klo];
        a = (x[khi] - tk)/h;
        b = (tk - x[klo])/h;
        value   = a*y[klo] + b*y[khi] + ((a*a*a - a)*y2[klo] + (b*b*b - b)*y2[khi])*(h*h)/6.0;
        dvalue  = (y[khi] - y[klo])/h - ((3.0*a*a - 1.0)*y2[klo] - (3.0*b*b - 1.0)*y2[khi])*h/6.0;
        d2value = a*y2[klo] + b*y2[khi];
        d3value = (y2[khi]-y2[klo])/h;
        tkOld = tk;
    }

    if      (order == 0) return value;
    else if (order == 1) return dvalue;
    else if (order == 2) return d2value;
    else if (order == 3) return d3value;
    else                 return 0.0;
}

/*
 * Global (to this file): variables
 */

#define R  8.3143
#define NR 4      /* Four independent composition variables */
#define NS 3      /* Three ordering parameters              */
#define NA 5      /* Five endmember compositions            */

#define MAX_ITER 100  /* Max number of iterations allowed in ordering alg */

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

/*static double getTOldPure() {
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
    }*/

/***********************************/
/* Statics for Site Mole Fractions */
/***********************************/

static double xmg2a; /* site mole fractions */
static double xfe2a;
static double xmn2a;
static double xti4a;
static double xfe3a;
static double xal3a;

static double xmg2b;
static double xfe2b;
static double xmn2b;
static double xti4b;
static double xfe3b;
static double xal3b;

static double xmg2ID; /* ideal mixing site mole fractions */
static double xfe2ID;
static double xmn2ID;
static double xti4ID;
static double xfe3ID;
static double xal3ID;

static int noSave = FALSE;

void resetValueOfDhrnd   (double newDhrnd )   { dhilm     = newDhrnd;	 dhgei    = newDhrnd;	 dhpyr    = newDhrnd;	 noSave = TRUE; }
void resetValueOfDhilm   (double newDhilm )   { dhilm     = newDhilm;							 noSave = TRUE; }
void resetValueOfDhgei   (double newDhgei )   { dhgei     = newDhgei;							 noSave = TRUE; }
void resetValueOfDhpyr   (double newDhpyr )   { dhpyr     = newDhpyr;							 noSave = TRUE; }

void resetValueOfWhs     (double newWhs   )   { whilm     = newWhs;	 whgei    = newWhs;	 whpyr    = newWhs;	 noSave = TRUE; }
void resetValueOfWhilm   (double newWhilm )   { whilm     = newWhilm;							 noSave = TRUE; }
void resetValueOfWhgei   (double newWhgei )   { whgei	  = newWhgei;							 noSave = TRUE; }
void resetValueOfWhpyr   (double newWhpyr )   { whpyr	  = newWhpyr;							 noSave = TRUE; }

void resetValueOfWhmcrn  (double newWhmcrn)   { whmcrn    = newWhmcrn;  						 noSave = TRUE; }
void resetValueOfWhmilm  (double newWhmilm)   { whhmilm   = newWhmilm;  						 noSave = TRUE; }
void resetValueOfWhmgei  (double newWhmgei)   { whhmgei   = newWhmgei;  						 noSave = TRUE; }
void resetValueOfWhmpyr  (double newWhmpyr)   { whhmpyr   = newWhmpyr;  						 noSave = TRUE; }
void resetValueOfWcrnilm (double newWcrnilm)  { wcrnilm   = newWcrnilm; 						 noSave = TRUE; }
void resetValueOfWcrngei (double newWcrngei)  { wcrngei   = newWcrngei; 						 noSave = TRUE; }
void resetValueOfWcrnpyr (double newWcrnpyr)  { wcrnpyr   = newWcrnpyr; 						 noSave = TRUE; }
void resetValueOfWilmgei (double newWilmgei)  { whilmgei  = newWilmgei; 						 noSave = TRUE; }
void resetValueOfWilmpyr (double newWilmpyr)  { whilmpyr  = newWilmpyr; 						 noSave = TRUE; }
void resetValueOfWgeipyr (double newWgeipyr)  { whgeipyr  = newWgeipyr; 						 noSave = TRUE; }
void resetValueOfWilmgeiT(double newWilmgeiT) { whilmgeiT = newWilmgeiT;   						 noSave = TRUE; }
void resetValueOfWilmpyrT(double newWilmpyrT) { whilmpyrT = newWilmpyrT;   						 noSave = TRUE; }
void resetValueOfWgeipyrT(double newWgeipyrT) { whgeipyrT = newWgeipyrT;   						 noSave = TRUE; }

void resetValueOfWilmilmgei(double newWilmilmgei) { whilmilmgei = newWilmilmgei;   				         noSave = TRUE; }
void resetValueOfWilmgeigei(double newWilmgeigei) { whilmgeigei = newWilmgeigei;   				         noSave = TRUE; }
void resetValueOfWilmilmpyr(double newWilmilmpyr) { whilmilmpyr = newWilmilmpyr;   				         noSave = TRUE; }
void resetValueOfWilmpyrpyr(double newWilmpyrpyr) { whilmpyrpyr = newWilmpyrpyr;   				         noSave = TRUE; }
void resetValueOfWgeigeipyr(double newWgeigeipyr) { whgeigeipyr = newWgeigeipyr;   				         noSave = TRUE; }
void resetValueOfWgeipyrpyr(double newWgeipyrpyr) { whgeipyrpyr = newWgeipyrpyr;   				         noSave = TRUE; }

void resetValueOfDwhhmrnd(double newDwhhmrnd) { dwhhmilm  = newDwhhmrnd; dwhhmgei = newDwhhmrnd; dwhhmpyr = newDwhhmrnd; noSave = TRUE; }
void resetValueOfDwhhmilm(double newDwhhmilm) { dwhhmilm  = newDwhhmilm;						 noSave = TRUE; }
void resetValueOfDwhhmgei(double newDwhhmgei) { dwhhmgei  = newDwhhmgei;						 noSave = TRUE; }
void resetValueOfDwhhmpyr(double newDwhhmpyr) { dwhhmpyr  = newDwhhmpyr;						 noSave = TRUE; }

void resetValueOfWhhmrnd2(double newWhhmrnd2) { whhmilm2  = newWhhmrnd2; whhmgei2 = newWhhmrnd2; whhmpyr2 = newWhhmrnd2; noSave = TRUE; }
void resetValueOfWhhmilm2(double newWhhmilm2) { whhmilm2  = newWhhmilm2;						 noSave = TRUE; }
void resetValueOfWhhmgei2(double newWhhmgei2) { whhmgei2  = newWhhmgei2;						 noSave = TRUE; }
void resetValueOfWhhmpyr2(double newWhhmpyr2) { whhmpyr2  = newWhhmpyr2;						 noSave = TRUE; }

void resetValueOfSROconst (double newSROconst) { SROconst = newSROconst;
                                                 SRO600   = newSROconst;
  						 SRO700   = newSROconst;
  						 SRO800   = newSROconst;
  						 SRO900   = newSROconst;
  						 SRO1000  = newSROconst;
  						 SRO1100  = newSROconst;
  						 SRO1200  = newSROconst;
  						 SRO1300  = newSROconst; makeSpline = TRUE;  				 noSave = TRUE; }
void resetValueOfSRO600   (double newSRO600  ) { SRO600   = newSRO600;   makeSpline = TRUE;  				 noSave = TRUE; }
void resetValueOfSRO700   (double newSRO700  ) { SRO700   = newSRO700;   makeSpline = TRUE;  				 noSave = TRUE; }
void resetValueOfSRO800   (double newSRO800  ) { SRO800   = newSRO800;   makeSpline = TRUE;  				 noSave = TRUE; }
void resetValueOfSRO900   (double newSRO900  ) { SRO900   = newSRO900;   makeSpline = TRUE;  				 noSave = TRUE; }
void resetValueOfSRO1000  (double newSRO1000 ) { SRO1000  = newSRO1000;  makeSpline = TRUE;  				 noSave = TRUE; }
void resetValueOfSRO1100  (double newSRO1100 ) { SRO1100  = newSRO1100;  makeSpline = TRUE;  				 noSave = TRUE; }
void resetValueOfSRO1200  (double newSRO1200 ) { SRO1200  = newSRO1200;  makeSpline = TRUE;  				 noSave = TRUE; }
void resetValueOfSRO1300  (double newSRO1300 ) { SRO1300  = newSRO1300;  makeSpline = TRUE;  				 noSave = TRUE; }

/*
 *=============================================================================
 * Local function to compute ordering state and associated derivatives of
 * pure component endmembers
 */

#define IL_S          -R*((1.0+s[0])*log(1.0+s[0]) + (1.0-s[0])*log(1.0-s[0]) - 2.0*log(2.0))
#define IL_H          (dhilm+(p-1.0)*dvilm)*(1.0 - s[0]*s[0]) + (whilm+(p-1.0)*wvilm)*s[0]*s[0]*(1.0 - s[0]*s[0])
#define IL_G          IL_H - t*(IL_S)

#define DIL_GDS0      R*t*(log(1.0+s[0]) - log(1.0-s[0])) \
                       - 2.0*(dhilm+(p-1.0)*dvilm)*s[0] + (whilm+(p-1.0)*wvilm)*(2.0*s[0] - 4.0*s[0]*s[0]*s[0])
#define DIL_GDT       -(IL_S)
#define DIL_GDP       dvilm*(1.0 - s[0]*s[0]) + wvilm*s[0]*s[0]*(1.0 - s[0]*s[0])

#define D2IL_GDS0S0   R*t*(1.0/(1.0+s[0]) + 1.0/(1.0-s[0])) - 2.0*(dhilm+(p-1.0)*dvilm) + (whilm+(p-1.0)*wvilm)*(2.0 - 12.0*s[0]*s[0])
#define D2IL_GDS0DT   R*(log(1.0+s[0]) - log(1.0-s[0]))
#define D2IL_GDS0DP   - 2.0*dvilm*s[0] + wvilm*(2.0*s[0] - 4.0*s[0]*s[0]*s[0])
#define D2IL_GDT2     0.0
#define D2IL_GDTP     0.0
#define D2IL_GDP2     0.0

#define D3IL_GDS0S0S0 R*t*(-1.0/SQUARE(1.0+s[0]) + 1.0/SQUARE(1.0-s[0])) - 24.0*(whilm+(p-1.0)*wvilm)*s[0]
#define D3IL_GDS0S0DT R*(1.0/(1.0+s[0]) + 1.0/(1.0-s[0]))
#define D3IL_GDS0S0DP - 2.0*dvilm + wvilm*(2.0 - 12.0*s[0]*s[0])
#define D3IL_GDS0DT2  0.0
#define D3IL_GDS0DTDP 0.0
#define D3IL_GDS0DP2  0.0
#define D3IL_GDT3     0.0
#define D3IL_GDT2DP   0.0
#define D3IL_GDTDP2   0.0
#define D3IL_GDP3     0.0

#define GK_S          -R*((1.0+s[1])*log(1.0+s[1]) + (1.0-s[1])*log(1.0-s[1]) - 2.0*log(2.0))
#define GK_H          (dhgei+(p-1.0)*dvgei)*(1.0 - s[1]*s[1]) + (whgei+(p-1.0)*wvgei)*s[1]*s[1]*(1.0 - s[1]*s[1])
#define GK_G          GK_H - t*(GK_S)

#define DGK_GDS1      R*t*(log(1.0+s[1]) - log(1.0-s[1])) \
                       - 2.0*(dhgei+(p-1.0)*dvgei)*s[1] + (whgei+(p-1.0)*wvgei)*(2.0*s[1] - 4.0*s[1]*s[1]*s[1])
#define DGK_GDT       -(GK_S)
#define DGK_GDP       dvgei*(1.0 - s[1]*s[1]) + wvgei*s[1]*s[1]*(1.0 - s[1]*s[1])

#define D2GK_GDS1S1   R*t*(1.0/(1.0+s[1]) + 1.0/(1.0-s[1])) - 2.0*(dhgei+(p-1.0)*dvgei) + (whgei+(p-1.0)*wvgei)*(2.0 - 12.0*s[1]*s[1])
#define D2GK_GDS1DT   R*(log(1.0+s[1]) - log(1.0-s[1]))
#define D2GK_GDS1DP   - 2.0*dvgei*s[1] + wvgei*(2.0*s[1] - 4.0*s[1]*s[1]*s[1])
#define D2GK_GDT2     0.0
#define D2GK_GDTP     0.0
#define D2GK_GDP2     0.0

#define D3GK_GDS1S1S1 R*t*(-1.0/SQUARE(1.0+s[1]) + 1.0/SQUARE(1.0-s[1])) - 24.0*(whgei+(p-1.0)*wvgei)*s[1]
#define D3GK_GDS1S1DT R*(1.0/(1.0+s[1]) + 1.0/(1.0-s[1]))
#define D3GK_GDS1S1DP - 2.0*dvgei + wvgei*(2.0 - 12.0*s[1]*s[1])
#define D3GK_GDS1DT2  0.0
#define D3GK_GDS1DTDP 0.0
#define D3GK_GDS1DP2  0.0
#define D3GK_GDT3     0.0
#define D3GK_GDT2DP   0.0
#define D3GK_GDTDP2   0.0
#define D3GK_GDP3     0.0

#define PY_S          -R*((1.0+s[2])*log(1.0+s[2]) + (1.0-s[2])*log(1.0-s[2]) - 2.0*log(2.0))
#define PY_H          (dhpyr+(p-1.0)*dvpyr)*(1.0 - s[2]*s[2]) + (whpyr+(p-1.0)*wvpyr)*s[2]*s[2]*(1.0 - s[2]*s[2])
#define PY_G          PY_H - t*(PY_S)

#define DPY_GDS2      R*t*(log(1.0+s[2]) - log(1.0-s[2])) \
                       - 2.0*(dhpyr+(p-1.0)*dvpyr)*s[2] + (whpyr+(p-1.0)*wvpyr)*(2.0*s[2] - 4.0*s[2]*s[2]*s[2])
#define DPY_GDT       -(PY_S)
#define DPY_GDP       dvpyr*(1.0 - s[2]*s[2]) + wvpyr*s[2]*s[2]*(1.0 - s[2]*s[2])

#define D2PY_GDS2S2   R*t*(1.0/(1.0+s[2]) + 1.0/(1.0-s[2])) - 2.0*(dhpyr+(p-1.0)*dvpyr) + (whpyr+(p-1.0)*wvpyr)*(2.0 - 12.0*s[2]*s[2])
#define D2PY_GDS2DT   R*(log(1.0+s[2]) - log(1.0-s[2]))
#define D2PY_GDS2DP   - 2.0*dvpyr*s[2] + wvpyr*(2.0*s[2] - 4.0*s[2]*s[2]*s[2])
#define D2PY_GDT2     0.0
#define D2PY_GDTP     0.0
#define D2PY_GDP2     0.0

#define D3PY_GDS2S2S2 R*t*(-1.0/SQUARE(1.0+s[2]) + 1.0/SQUARE(1.0-s[2])) - 24.0*(whpyr+(p-1.0)*wvpyr)*s[2]
#define D3PY_GDS2S2DT R*(1.0/(1.0+s[2]) + 1.0/(1.0-s[2]))
#define D3PY_GDS2S2DP - 2.0*dvpyr + wvpyr*(2.0 - 12.0*s[2]*s[2])
#define D3PY_GDS2DT2  0.0
#define D3PY_GDS2DTDP 0.0
#define D3PY_GDS2DP2  0.0
#define D3PY_GDT3     0.0
#define D3PY_GDT2DP   0.0
#define D3PY_GDTDP2   0.0
#define D3PY_GDP3     0.0

#define HM_S          fSRO(t,0)*2.0*R*log(2.0)
#define HM_H          0.0
#define HM_G          HM_H - t*(HM_S)
#define DHM_GDT       - fSRO(t,0)*2.0*R*log(2.0) - fSRO(t,1)*2.0*R*t*log(2.0)
#define DHM_GDP       0.0
#define D2HM_GDT2     - fSRO(t,1)*4.0*R*log(2.0) - fSRO(t,2)*2.0*R*t*log(2.0)
#define D2HM_GDTP     0.0
#define D2HM_GDP2     0.0
#define D3HM_GDT3     - fSRO(t,2)*6.0*R*log(2.0) - fSRO(t,3)*2.0*R*t*log(2.0)
#define D3HM_GDT2DP   0.0
#define D3HM_GDTDP2   0.0
#define D3HM_GDP3     0.0

#define CR_S          fSRO(t,0)*2.0*R*log(2.0)
#define CR_H          0.0
#define CR_G          CR_H - t*(CR_S)
#define DCR_GDT       - fSRO(t,0)*2.0*R*log(2.0) - fSRO(t,1)*2.0*R*t*log(2.0)
#define DCR_GDP       0.0
#define D2CR_GDT2     - fSRO(t,1)*4.0*R*log(2.0) - fSRO(t,2)*2.0*R*t*log(2.0)
#define D2CR_GDTP     0.0
#define D2CR_GDP2     0.0
#define D3CR_GDT3     - fSRO(t,2)*6.0*R*log(2.0) - fSRO(t,3)*2.0*R*t*log(2.0)
#define D3CR_GDT2DP   0.0
#define D3CR_GDTDP2   0.0
#define D3CR_GDP3     0.0

#define fillD2GDSDT   d2gdsdt[0]   = D2IL_GDS0DT;   d2gdsdt[1]   = D2GK_GDS1DT;   d2gdsdt[2]   = D2PY_GDS2DT;
#define fillD2GDSDP   d2gdsdp[0]   = D2IL_GDS0DP;   d2gdsdp[1]   = D2GK_GDS1DP;   d2gdsdp[2]   = D2PY_GDS2DP;
#define fillD3GDS3    d3gds3[0]    = D3IL_GDS0S0S0; d3gds3[1]    = D3GK_GDS1S1S1; d3gds3[2]    = D3PY_GDS2S2S2;
#define fillD3GDS2DT  d3gds2dt[0]  = D3IL_GDS0S0DT; d3gds2dt[1]  = D3GK_GDS1S1DT; d3gds2dt[2]  = D3PY_GDS2S2DT;
#define fillD3GDS2DP  d3gds2dp[0]  = D3IL_GDS0S0DP; d3gds2dp[1]  = D3GK_GDS1S1DP; d3gds2dp[2]  = D3PY_GDS2S2DP;
#define fillD3GDSDT2  d3gdsdt2[0]  = D3IL_GDS0DT2;  d3gdsdt2[1]  = D3GK_GDS1DT2;  d3gdsdt2[2]  = D3PY_GDS2DT2;
#define fillD3GDSDTDP d3gdsdtdp[0] = D3IL_GDS0DTDP; d3gdsdtdp[1] = D3GK_GDS1DTDP; d3gdsdtdp[2] = D3PY_GDS2DTDP;
#define fillD3GDSDP2  d3gdsdp2[0]  = D3IL_GDS0DP2;  d3gdsdp2[1]  = D3GK_GDS1DP2;  d3gdsdp2[2]  = D3PY_GDS2DP2;

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
    static double tOld = -9999.0;
    static double pOld = -9999.0;
    static double sOld[NS];
    static double d2gds2[NS];
    static const int iterMax = 1000;
    int i;

    if ((t != tOld) || (p != pOld) || noSave) {
        double dgds[NS], sNew[NS];
        int iter = 0;
        for (i=0; i<NS; i++) { sOld[i] = 2.0; sNew[i] = 0.98; }
        while ((iter < iterMax) && (
           (ABS(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON) ||
           (ABS(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON) ||
           (ABS(sNew[2]-sOld[2]) > 10.0*DBL_EPSILON) )) {
            double s[NS];

            for (i=0; i<NS; i++) s[i] = sNew[i];

            dgds[0] = DIL_GDS0;
            dgds[1] = DGK_GDS1;
            dgds[2] = DPY_GDS2;

            d2gds2[0] = D2IL_GDS0S0;
            d2gds2[1] = D2GK_GDS1S1;
            d2gds2[2] = D2PY_GDS2S2;

            for (i=0; i<NS; i++) sOld[i] = s[i];

            for (i=0; i<NS; i++) {
         s[i] += - dgds[i]/d2gds2[i];
         s[i] = MIN(s[i],  1.0 - 10.0*DBL_EPSILON);
         s[i] = MAX(s[i],  0.0);
            }

            for (i=0; i<NS; i++) sNew[i] = s[i];
            iter++;
        }
        tOld = t;
        pOld = p;
        if (iter == iterMax) {
#ifdef DEBUG
            printf("Iteration limit exceeded in pureOrder (rhomsghiorso.c)\n");
            printf("  t,p    = %13.6g, %13.6g\n",         t,         p);
            printf("  s      = %13.6g, %13.6g, %13.6g\n", sOld[0],   sOld[1],   sOld[2]);
            printf("  dgds   = %13.6g, %13.6g, %13.6g\n", dgds[0],   dgds[1],   dgds[2]);
            printf("  d2gds2 = %13.6g, %13.6g, %13.6g\n", d2gds2[0], d2gds2[1], d2gds2[2]);
#endif
            /* Condition arises at a T above the 1st order transition when Newton's method has
         arrived at a local minimum.  In this case, the global minimum is at zero.       */
            for (i=0; i<NS; i++) sOld[i] = 0.0;
        }
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
        double *s = sOld;
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
 d2gds2[0][0] = 0.0;         d2gds2[0][1] = D2GK_GDS1S1;  d2gds2[0][2] = 0.0;         \
 d2gds2[1][0] = 0.0;         d2gds2[1][1] = 0.0;          d2gds2[1][2] = 0.0;         \
 d2gds2[2][0] = D2IL_GDS0S0; d2gds2[2][1] = 0.0;          d2gds2[2][2] = 0.0;         \
 d2gds2[3][0] = 0.0;         d2gds2[3][1] = 0.0;          d2gds2[3][2] = D2PY_GDS2S2; \
 d2gds2[4][0] = 0.0;         d2gds2[4][1] = 0.0;          d2gds2[4][2] = 0.0;

#define fillD2GDSDT \
 d2gdsdt[0][0] = 0.0;         d2gdsdt[0][1] = D2GK_GDS1DT;  d2gdsdt[0][2] =  0.0;        \
 d2gdsdt[1][0] = 0.0;         d2gdsdt[1][1] = 0.0;          d2gdsdt[1][2] =  0.0;        \
 d2gdsdt[2][0] = D2IL_GDS0DT; d2gdsdt[2][1] = 0.0;          d2gdsdt[2][2] =  0.0;        \
 d2gdsdt[3][0] = 0.0;         d2gdsdt[3][1] = 0.0;          d2gdsdt[3][2] = D2PY_GDS2DT; \
 d2gdsdt[4][0] = 0.0;         d2gdsdt[4][1] = 0.0;          d2gdsdt[4][2] =  0.0;

#define fillD2GDSDP \
 d2gdsdp[0][0] = 0.0;         d2gdsdp[0][1] = D2GK_GDS1DP;  d2gdsdp[0][2] = 0.0;         \
 d2gdsdp[1][0] = 0.0;         d2gdsdp[1][1] = 0.0;          d2gdsdp[1][2] = 0.0;         \
 d2gdsdp[2][0] = D2IL_GDS0DP; d2gdsdp[2][1] = 0.0;          d2gdsdp[2][2] = 0.0;         \
 d2gdsdp[3][0] = 0.0;         d2gdsdp[3][1] = 0.0;          d2gdsdp[3][2] = D2PY_GDS2DP; \
 d2gdsdp[4][0] = 0.0;         d2gdsdp[4][1] = 0.0;          d2gdsdp[4][2] = 0.0;

#define fillD2GDT2 \
 d2gdt2[0]  = D2GK_GDT2; d2gdt2[1]  = D2HM_GDT2; d2gdt2[2]  = D2IL_GDT2; d2gdt2[3]  = D2PY_GDT2; d2gdt2[4]  = D2CR_GDT2;

#define fillD2GDTDP \
 d2gdtdp[0] = D2GK_GDTP; d2gdtdp[1] = D2HM_GDTP; d2gdtdp[2] = D2IL_GDTP; d2gdtdp[3] = D2PY_GDTP; d2gdtdp[4] = D2CR_GDTP;

#define fillD2GDP2 \
 d2gdp2[0]  = D2GK_GDP2; d2gdp2[1]  = D2HM_GDP2; d2gdp2[2]  = D2IL_GDP2; d2gdp2[3]  = D2PY_GDP2; d2gdp2[4]  = D2CR_GDP2;

#define fillD3GDS3 \
 d3gds3[0][0] = 0.0;           d3gds3[0][1] = D3GK_GDS1S1S1; d3gds3[0][2] = 0.0;           \
 d3gds3[1][0] = 0.0;           d3gds3[1][1] = 0.0;           d3gds3[1][2] = 0.0;           \
 d3gds3[2][0] = D3IL_GDS0S0S0; d3gds3[2][1] = 0.0;           d3gds3[2][2] = 0.0;           \
 d3gds3[3][0] = 0.0;           d3gds3[3][1] = 0.0;           d3gds3[3][2] = D3PY_GDS2S2S2; \
 d3gds3[4][0] = 0.0;           d3gds3[4][1] = 0.0;           d3gds3[4][2] = 0.0;

#define fillD3GDS2DT \
 d3gds2dt[0][0] = 0.0;           d3gds2dt[0][1] = D3GK_GDS1S1DT; d3gds2dt[0][2] = 0.0;           \
 d3gds2dt[1][0] = 0.0;           d3gds2dt[1][1] = 0.0;           d3gds2dt[1][2] = 0.0;           \
 d3gds2dt[2][0] = D3IL_GDS0S0DT; d3gds2dt[2][1] = 0.0;           d3gds2dt[2][2] = 0.0;           \
 d3gds2dt[3][0] = 0.0;           d3gds2dt[3][1] = 0.0;           d3gds2dt[3][2] = D3PY_GDS2S2DT; \
 d3gds2dt[4][0] = 0.0;           d3gds2dt[4][1] = 0.0;           d3gds2dt[4][2] = 0.0;

#define fillD3GDS2DP \
 d3gds2dp[0][0] = 0.0;           d3gds2dp[0][1] = D3GK_GDS1S1DP; d3gds2dp[0][2] = 0.0;           \
 d3gds2dp[1][0] = 0.0;           d3gds2dp[1][1] = 0.0;           d3gds2dp[1][2] = 0.0;           \
 d3gds2dp[2][0] = D3IL_GDS0S0DP; d3gds2dp[2][1] = 0.0;           d3gds2dp[2][2] = 0.0;           \
 d3gds2dp[3][0] = 0.0;           d3gds2dp[3][1] = 0.0;           d3gds2dp[3][2] = D3PY_GDS2S2DP; \
 d3gds2dp[4][0] = 0.0;           d3gds2dp[4][1] = 0.0;           d3gds2dp[4][2] = 0.0;

#define fillD3GDSDT2 \
 d3gdsdt2[0][0] = 0.0;          d3gdsdt2[0][1] = D3GK_GDS1DT2; d3gdsdt2[0][2] = 0.0;           \
 d3gdsdt2[1][0] = 0.0;          d3gdsdt2[1][1] = 0.0;          d3gdsdt2[1][2] = 0.0;           \
 d3gdsdt2[2][0] = D3IL_GDS0DT2; d3gdsdt2[2][1] = 0.0;          d3gdsdt2[2][2] = 0.0;           \
 d3gdsdt2[3][0] = 0.0;          d3gdsdt2[3][1] = 0.0;          d3gdsdt2[3][2] =  D3PY_GDS2DT2; \
 d3gdsdt2[4][0] = 0.0;          d3gdsdt2[4][1] = 0.0;          d3gdsdt2[4][2] = 0.0;

#define fillD3GDSDTDP \
 d3gdsdtdp[0][0] = 0.0;           d3gdsdtdp[0][1] = D3GK_GDS1DTDP; d3gdsdtdp[0][2] = 0.0;           \
 d3gdsdtdp[1][0] = 0.0;           d3gdsdtdp[1][1] = 0.0;           d3gdsdtdp[1][2] = 0.0;           \
 d3gdsdtdp[2][0] = D3IL_GDS0DTDP; d3gdsdtdp[2][1] = 0.0;           d3gdsdtdp[2][2] = 0.0;           \
 d3gdsdtdp[3][0] = 0.0;           d3gdsdtdp[3][1] = 0.0;           d3gdsdtdp[3][2] = D3PY_GDS2DTDP; \
 d3gdsdtdp[4][0] = 0.0;           d3gdsdtdp[4][1] = 0.0;           d3gdsdtdp[4][2] = 0.0;

#define fillD3GDSDP2 \
 d3gdsdp2[0][0] = 0.0;          d3gdsdp2[0][1] = D3GK_GDS1DP2; d3gdsdp2[0][2] = 0.0;          \
 d3gdsdp2[1][0] = 0.0;          d3gdsdp2[1][1] = 0.0;          d3gdsdp2[1][2] = 0.0;          \
 d3gdsdp2[2][0] = D3IL_GDS0DP2; d3gdsdp2[2][1] = 0.0;          d3gdsdp2[2][2] = 0.0;          \
 d3gdsdp2[3][0] = 0.0;          d3gdsdp2[3][1] = 0.0;          d3gdsdp2[3][2] = D3PY_GDS2DP2; \
 d3gdsdp2[4][0] = 0.0;          d3gdsdp2[4][1] = 0.0;          d3gdsdp2[4][2] = 0.0;

#define fillD3GDT3 \
 d3gdt3[0]   = D3GK_GDT3;   d3gdt3[1]   = D3HM_GDT3;   d3gdt3[2]   = D3IL_GDT3;   d3gdt3[3]   = D3PY_GDT3;   d3gdt3[4]   = D3CR_GDT3;

#define fillD3GDT2DP \
 d3gdt2dp[0] = D3GK_GDT2DP; d3gdt2dp[1] = D3HM_GDT2DP; d3gdt2dp[2] = D3IL_GDT2DP; d3gdt2dp[3] = D3PY_GDT2DP; d3gdt2dp[4] = D3CR_GDT2DP;

#define fillD3GDTDP2 \
 d3gdtdp2[0] = D3GK_GDTDP2; d3gdtdp2[1] = D3HM_GDTDP2; d3gdtdp2[2] = D3IL_GDTDP2; d3gdtdp2[3] = D3PY_GDTDP2; d3gdtdp2[4] = D3CR_GDTDP2;

#define fillD3GDP3 \
 d3gdp3[0]   = D3GK_GDP3;   d3gdp3[1]   = D3HM_GDP3;   d3gdp3[2]   = D3IL_GDP3;   d3gdp3[3]   = D3PY_GDP3;   d3gdp3[4]   = D3CR_GDP3;

static void
pureRhm(int mask, double t, double p,
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

    pureOrder(FIRST, t, p, s, NULL, NULL, NULL, NULL, NULL);

    if (mask & FIRST) {
        a[0] = GK_G;
        a[0] = exp(a[0]/(R*t));
        a[1] = HM_G;
        a[1] = exp(a[1]/(R*t));
        a[2] = IL_G;
        a[2] = exp(a[2]/(R*t));
        a[3] = PY_G;
        a[3] = exp(a[3]/(R*t));
        a[4] = CR_G;
        a[4] = exp(a[4]/(R*t));
    }

    if (mask & SECOND) {
        mu[0] = GK_G;
        mu[1] = HM_G;
        mu[2] = IL_G;
        mu[3] = PY_G;
        mu[4] = CR_G;
    }

    if (mask & THIRD) {
        gmix[0] = GK_G;
        gmix[1] = HM_G;
        gmix[2] = IL_G;
        gmix[3] = PY_G;
        gmix[4] = CR_G;
    }

    if (mask & FOURTH) {
        hmix[0] = (GK_G) + t*(GK_S);
        hmix[1] = (HM_G) - t*(DHM_GDT);
        hmix[2] = (IL_G) + t*(IL_S);
        hmix[3] = (PY_G) + t*(PY_S);
        hmix[4] = (CR_G) - t*(DCR_GDT);
    }

    if (mask & FIFTH) {
        smix[0] = GK_S;
        smix[1] = -(DHM_GDT);
        smix[2] = IL_S;
        smix[3] = PY_S;
        smix[4] = -(DCR_GDT);
    }

    if (mask & SIXTH) {
        double d2gdsdt[NA][NS], d2gds2[NA][NS], dsdt[NS];

        fillD2GDS2
        fillD2GDSDT

        pureOrder(SECOND, t, p, NULL, dsdt, NULL, NULL, NULL, NULL);

        cpmix[0] = D2GK_GDT2;
        cpmix[1] = D2HM_GDT2;
        cpmix[2] = D2IL_GDT2;
        cpmix[3] = D2PY_GDT2;
        cpmix[4] = D2CR_GDT2;

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

        pureOrder(SECOND | FOURTH, t, p, NULL, dsdt, NULL, d2sdt2, NULL, NULL);

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
        vmix[0] = DGK_GDP;
        vmix[1] = DHM_GDP;
        vmix[2] = DIL_GDP;
        vmix[3] = DPY_GDP;
        vmix[4] = DCR_GDP;
    }

    if(mask & NINTH) {
        double d2gdsdt[NA][NS], d2gdsdp[NA][NS], d2gds2[NA][NS], d2gdtdp[NA],
            dsdt[NS], dsdp[NS];

        fillD2GDSDT
        fillD2GDS2
        fillD2GDSDP
        fillD2GDTDP

        pureOrder(SECOND | THIRD, t, p, NULL, dsdt, dsdp, NULL, NULL, NULL);

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

        pureOrder(THIRD, t, p, NULL, NULL, dsdp, NULL, NULL, NULL);

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

        pureOrder(SECOND | THIRD | FOURTH | FIFTH, t, p, NULL, dsdt, dsdp, d2sdt2, d2sdtdp, NULL);

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

        pureOrder(SECOND | THIRD | FIFTH | SIXTH, t, p, NULL, dsdt, dsdp, NULL, d2sdtdp, d2sdp2);

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

        pureOrder(THIRD | SIXTH, t, p, NULL, NULL, dsdp, NULL, NULL, d2sdp2);

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

#undef IL_S
#undef IL_H
#undef IL_G
#undef DIL_GDS0
#undef DIL_GDT
#undef DIL_GDP
#undef D2IL_GDS0S0
#undef D2IL_GDS0DT
#undef D2IL_GDS0DP
#undef D2IL_GDT2
#undef D2IL_GDTP
#undef D2IL_GDP2
#undef D3IL_GDS0S0S0
#undef D3IL_GDS0S0DT
#undef D3IL_GDS0S0DP
#undef D3IL_GDS0DT2
#undef D3IL_GDS0DTDP
#undef D3IL_GDS0DP2
#undef D3IL_GDT3
#undef D3IL_GDT2DP
#undef D3IL_GDTDP2
#undef D3IL_GDP3

#undef GK_S
#undef GK_H
#undef GK_G
#undef DGK_GDS1
#undef DGK_GDT
#undef DGK_GDP
#undef D2GK_GDS1S1
#undef D2GK_GDS1DT
#undef D2GK_GDS1DP
#undef D2GK_GDT2
#undef D2GK_GDTP
#undef D2GK_GDP2
#undef D3GK_GDS1S1S1
#undef D3GK_GDS1S1DT
#undef D3GK_GDS1S1DP
#undef D3GK_GDS1DT2
#undef D3GK_GDS1DTDP
#undef D3GK_GDS1DP2
#undef D3GK_GDT3
#undef D3GK_GDT2DP
#undef D3GK_GDTDP2
#undef D3GK_GDP3

#undef PY_S
#undef PY_H
#undef PY_G
#undef DPY_GDS2
#undef DPY_GDT
#undef DPY_GDP
#undef D2PY_GDS2S2
#undef D2PY_GDS2DT
#undef D2PY_GDS2DP
#undef D2PY_GDT2
#undef D2PY_GDTP
#undef D2PY_GDP2
#undef D3PY_GDS2S2S2
#undef D3PY_GDS2S2DT
#undef D3PY_GDS2S2DP
#undef D3PY_GDS2DT2
#undef D3PY_GDS2DTDP
#undef D3PY_GDS2DP2
#undef D3PY_GDT3
#undef D3PY_GDT2DP
#undef D3PY_GDTDP2
#undef D3PY_GDP3

#undef HM_S
#undef HM_H
#undef HM_G
#undef DHM_GDT
#undef DHM_GDP
#undef D2HM_GDT2
#undef D2HM_GDTP
#undef D2HM_GDP2
#undef D3HM_GDT3
#undef D3HM_GDT2DP
#undef D3HM_GDTDP2
#undef D3HM_GDP3

#undef CR_S
#undef CR_H
#undef CR_G
#undef DCR_GDT
#undef DCR_GDP
#undef D2CR_GDT2
#undef D2CR_GDTP
#undef D2CR_GDP2
#undef D3CR_GDT3
#undef D3CR_GDT2DP
#undef D3CR_GDTDP2
#undef D3CR_GDP3

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conRhm defines the conversion from m[i], to r[j]
 */
                   /* Order: Xil, Xgk, Xpy, Xcn */
#define FR0(i)     (i == 2) ? 1.0 - r[0] : - r[0]
#define FR1(i)     (i == 0) ? 1.0 - r[1] : - r[1]
#define FR2(i)     (i == 3) ? 1.0 - r[2] : - r[2]
#define FR3(i)     (i == 4) ? 1.0 - r[3] : - r[3]

                   /* Order: s, t, u */
#define GS0(i)     (i == 2) ? 1.0 - s[0] : - s[0]
#define GS1(i)     (i == 0) ? 1.0 - s[1] : - s[1]
#define GS2(i)     (i == 3) ? 1.0 - s[2] : - s[2]

#define DFR0DR0(i) - 1.0
#define DFR1DR1(i) - 1.0
#define DFR2DR2(i) - 1.0
#define DFR3DR3(i) - 1.0

#define DGS0DS0(i) - 1.0
#define DGS1DS1(i) - 1.0
#define DGS2DS2(i) - 1.0

#define ENDMEMBERS (r[1]*ends[0] + (1.0-r[0]-r[1]-r[2]-r[3])*ends[1] + r[0]*ends[2] + r[2]*ends[3] + r[3]*ends[4])

#define DENDDR0    (ends[2] - ends[1])
#define DENDDR1    (ends[0] - ends[1])
#define DENDDR2    (ends[3] - ends[1])
#define DENDDR3    (ends[4] - ends[1])

/*
 * Global (to this file): derivative definitions
 */

#define S  -R*(  xfe2a*log(xfe2a) + xmg2a*log(xmg2a) + xmn2a*log(xmn2a) + xti4a*log(xti4a) \
               + xfe2b*log(xfe2b) + xmg2b*log(xmg2b) + xmn2b*log(xmn2b) + xti4b*log(xti4b) \
                - 2.0*xfe2ID*log(xfe2ID) - 2.0*xmg2ID*log(xmg2ID) - 2.0*xmn2ID*log(xmn2ID) - 2.0*xti4ID*log(xti4ID) \
	      ) \
           -(1.0-fSRO(t,0))*2.0*R*(  xfe3ID*log(xfe3ID) + xal3ID*log(xal3ID) + xfe2ID*log(xfe2ID) \
	                      + xmg2ID*log(xmg2ID) + xmn2ID*log(xmn2ID) + xti4ID*log(xti4ID) ) + fSRO(t,0)*2.0*R*log(2.0)

#define H    ((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
        + ((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
           + ((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
        +  whmcrn*r[3]*(1.0-r[0]-r[1]-r[2]-r[3]) \
        + (whilmgei+(p-1.0)*wvilmgei+whilmgeiT)*r[0]*r[1]/2.0 + (whilmgei-whilmgeiT)*s[0]*s[1]/2.0 \
        + (whilmpyr+(p-1.0)*wvilmpyr+whilmpyrT)*r[0]*r[2]/2.0 + (whilmpyr-whilmpyrT)*s[0]*s[2]/2.0 \
        + (whgeipyr+(p-1.0)*wvgeipyr+whgeipyrT)*r[1]*r[2]/2.0 + (whgeipyr-whgeipyrT)*s[1]*s[2]/2.0 \
        + (wcrnilm+(dwhcrnilm+(p-1.0)*dwvcrnilm)*(r[0]-r[3]))*r[0]*r[3] \
        + (wcrngei+(dwhcrngei+(p-1.0)*dwvcrngei)*(r[1]-r[3]))*r[1]*r[3] \
        + (wcrnpyr+(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*(r[2]-r[3]))*r[2]*r[3] \
        + (dhilm+(p-1.0)*dvilm+((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
                                                    -(dwhcrnilm+(p-1.0)*dwvcrnilm)*r[3]/2.0)*(r[0]*r[0]-s[0]*s[0]) \
        + (dhgei+(p-1.0)*dvgei+((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
                                                    -(dwhcrngei+(p-1.0)*dwvcrngei)*r[3]/2.0)*(r[1]*r[1]-s[1]*s[1]) \
        + (dhpyr+(p-1.0)*dvpyr+((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
                                                    -(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[3]/2.0)*(r[2]*r[2]-s[2]*s[2]) \
        + (whilm+(p-1.0)*wvilm)*s[0]*s[0]*(r[0]*r[0]-s[0]*s[0]) \
        + (whgei+(p-1.0)*wvgei)*s[1]*s[1]*(r[1]*r[1]-s[1]*s[1]) \
        + (whpyr+(p-1.0)*wvpyr)*s[2]*s[2]*(r[2]*r[2]-s[2]*s[2]) \
        + (whilmilmgei)*(r[0]*r[0]-s[0]*s[0])*r[1]/4.0 + (whilmgeigei)*(r[1]*r[1]-s[1]*s[1])*r[0]/4.0 \
        + (whilmilmpyr)*(r[0]*r[0]-s[0]*s[0])*r[2]/4.0 + (whilmpyrpyr)*(r[2]*r[2]-s[2]*s[2])*r[0]/4.0 \
        + (whgeigeipyr)*(r[1]*r[1]-s[1]*s[1])*r[2]/4.0 + (whgeipyrpyr)*(r[2]*r[2]-s[2]*s[2])*r[1]/4.0

#define G  H - t*(S)

#define DGDR0  R*t*(  0.5*log(xfe2a) + 0.5*log(xti4a) + 0.5*log(xfe2b) + 0.5*log(xti4b) - log(xfe2ID) - log(xti4ID) ) \
               + (1.0-fSRO(t,0))*2.0*R*t*(- log(xfe3ID) + 0.5*log(xfe2ID) + 0.5*log(xti4ID) ) \
               + ((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
                - ((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*r[1] \
               - ((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*r[2] \
                - 2.0*(dwhhmilm+(p-1.0)*dwvhmilm)*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                -     (dwhhmgei+(p-1.0)*dwvhmgei)*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
               -     (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                -  whmcrn*r[3] \
                + (whilmgei+(p-1.0)*wvilmgei+whilmgeiT)*r[1]/2.0 \
                + (whilmpyr+(p-1.0)*wvilmpyr+whilmpyrT)*r[2]/2.0 \
                + (wcrnilm+(dwhcrnilm+(p-1.0)*dwvcrnilm)*(r[0]-r[3]))*r[3] + (dwhcrnilm+(p-1.0)*dwvcrnilm)*r[0]*r[3] \
                + 2.0*(dhilm+(p-1.0)*dvilm+((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
	    			         -(dwhcrnilm+(p-1.0)*dwvcrnilm)*r[3]/2.0)*r[0] \
                - ((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*(r[0]*r[0]-s[0]*s[0]) \
                - ((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*(r[1]*r[1]-s[1]*s[1]) \
                - ((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*(r[2]*r[2]-s[2]*s[2]) \
                + 2.0*(whilm+(p-1.0)*wvilm)*s[0]*s[0]*r[0] \
                + (whilmilmgei)*r[0]*r[1]/2.0 + (whilmgeigei)*(r[1]*r[1]-s[1]*s[1])/4.0 \
                + (whilmilmpyr)*r[0]*r[2]/2.0 + (whilmpyrpyr)*(r[2]*r[2]-s[2]*s[2])/4.0

#define DGDR1  R*t*(  0.5*log(xmg2a) + 0.5*log(xti4a) + 0.5*log(xmg2b) + 0.5*log(xti4b) - log(xmg2ID) - log(xti4ID) ) \
               + (1.0-fSRO(t,0))*2.0*R*t*(- log(xfe3ID) + 0.5*log(xmg2ID) + 0.5*log(xti4ID) ) \
               - ((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*r[0] \
                + ((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
               - ((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*r[2] \
                -     (dwhhmilm+(p-1.0)*dwvhmilm)*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                - 2.0*(dwhhmgei+(p-1.0)*dwvhmgei)*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
               -     (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                -  whmcrn*r[3] \
                + (whilmgei+(p-1.0)*wvilmgei+whilmgeiT)*r[0]/2.0 \
                + (whgeipyr+(p-1.0)*wvgeipyr+whgeipyrT)*r[2]/2.0 \
                + (wcrngei+(dwhcrngei+(p-1.0)*dwvcrngei)*(r[1]-r[3]))*r[3] +(dwhcrngei+(p-1.0)*dwvcrngei)*r[1]*r[3] \
                - ((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*(r[0]*r[0]-s[0]*s[0]) \
                + 2.0*(dhgei+(p-1.0)*dvgei+((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
	    			         -(dwhcrngei+(p-1.0)*dwvcrngei)*r[3]/2.0)*r[1] \
                - ((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*(r[1]*r[1]-s[1]*s[1]) \
                - ((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*(r[2]*r[2]-s[2]*s[2]) \
                + 2.0*(whgei+(p-1.0)*wvgei)*s[1]*s[1]*r[1] \
                + (whilmilmgei)*(r[0]*r[0]-s[0]*s[0])/4.0 + (whilmgeigei)*r[1]*r[0]/2.0 \
                + (whgeigeipyr)*r[1]*r[2]/2.0 + (whgeipyrpyr)*(r[2]*r[2]-s[2]*s[2])/4.0

#define DGDR2  R*t*(  0.5*log(xmn2a) + 0.5*log(xti4a) + 0.5*log(xmn2b) + 0.5*log(xti4b) - log(xmn2ID) - log(xti4ID) ) \
               + (1.0-fSRO(t,0))*2.0*R*t*(- log(xfe3ID) + 0.5*log(xmn2ID) + 0.5*log(xti4ID) ) \
               - ((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*r[0] \
                - ((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*r[1] \
               + ((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
               -     (dwhhmilm+(p-1.0)*dwvhmilm)*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                -     (dwhhmgei+(p-1.0)*dwvhmgei)*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
               - 2.0*(dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                -  whmcrn*r[3] \
                + (whilmpyr+(p-1.0)*wvilmpyr+whilmpyrT)*r[0]/2.0 \
                + (whgeipyr+(p-1.0)*wvgeipyr+whgeipyrT)*r[1]/2.0 \
                + (wcrnpyr+(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*(r[2]-r[3]))*r[3] + (dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[2]*r[3] \
                - ((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*(r[0]*r[0]-s[0]*s[0]) \
                - ((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*(r[1]*r[1]-s[1]*s[1]) \
                + 2.0*(dhpyr+(p-1.0)*dvpyr+((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
                                        -(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[3]/2.0)*r[2] \
                - ((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*(r[2]*r[2]-s[2]*s[2]) \
                + 2.0*(whpyr+(p-1.0)*wvpyr)*s[2]*s[2]*r[2] \
                + (whilmilmpyr)*(r[0]*r[0]-s[0]*s[0])/4.0 + (whilmpyrpyr)*r[2]*r[0]/2.0 \
                + (whgeigeipyr)*(r[1]*r[1]-s[1]*s[1])/4.0 + (whgeipyrpyr)*r[2]*r[1]/2.0

#define DGDR3  (1.0-fSRO(t,0))*2.0*R*t*(- log(xfe3ID) + log(xal3ID) ) \
               - ((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*r[0] \
                - ((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*r[1] \
               - ((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*r[2] \
               - (dwhhmilm+(p-1.0)*dwvhmilm)*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                - (dwhhmgei+(p-1.0)*dwvhmgei)*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
               - (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                +  whmcrn*(1.0-r[0]-r[1]-r[2]-2.0*r[3]) \
                + (wcrnilm+(dwhcrnilm+(p-1.0)*dwvcrnilm)*(r[0]-r[3]))*r[0] - (dwhcrnilm+(p-1.0)*dwvcrnilm)*r[0]*r[3] \
                + (wcrngei+(dwhcrngei+(p-1.0)*dwvcrngei)*(r[1]-r[3]))*r[1] - (dwhcrngei+(p-1.0)*dwvcrngei)*r[1]*r[3] \
                + (wcrnpyr+(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*(r[2]-r[3]))*r[2] - (dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[2]*r[3] \
                - ((dwhhmilm+(p-1.0)*dwvhmilm)/2.0 + (whhmilm2+(p-1.0)*wvhmilm2)/4.0 + (dwhcrnilm+(p-1.0)*dwvcrnilm)/2.0)*(r[0]*r[0]-s[0]*s[0]) \
                - ((dwhhmgei+(p-1.0)*dwvhmgei)/2.0 + (whhmgei2+(p-1.0)*wvhmgei2)/4.0 + (dwhcrngei+(p-1.0)*dwvcrngei)/2.0)*(r[1]*r[1]-s[1]*s[1]) \
                - ((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0 + (whhmpyr2+(p-1.0)*wvhmpyr2)/4.0 + (dwhcrnpyr+(p-1.0)*dwvcrnpyr)/2.0)*(r[2]*r[2]-s[2]*s[2])

#define DGDS0  0.5*R*t*( log(xfe2a) - log(xti4a) - log(xfe2b) + log(xti4b) ) \
                + (whilmgei-whilmgeiT)*s[1]/2.0 + (whilmpyr-whilmpyrT)*s[2]/2.0 \
                - 2.0*(dhilm+(p-1.0)*dvilm+((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
	    			         -(dwhcrnilm+(p-1.0)*dwvcrnilm)*r[3]/2.0)*s[0] \
                + 2.0*(whilm+(p-1.0)*wvilm)*s[0]*(r[0]*r[0]-2.0*s[0]*s[0]) \
                - (whilmilmgei)*s[0]*r[1]/2.0 - (whilmilmpyr)*s[0]*r[2]/2.0

#define DGDS1  0.5*R*t*( log(xmg2a) - log(xti4a) - log(xmg2b) + log(xti4b) ) \
                + (whilmgei-whilmgeiT)*s[0]/2.0 + (whgeipyr-whgeipyrT)*s[2]/2.0 \
                - 2.0*(dhgei+(p-1.0)*dvgei+((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
	    			         -(dwhcrngei+(p-1.0)*dwvcrngei)*r[3]/2.0)*s[1] \
                + 2.0*(whgei+(p-1.0)*wvgei)*s[1]*(r[1]*r[1]-2.0*s[1]*s[1]) \
                - (whilmgeigei)*s[1]*r[0]/2.0 - (whgeigeipyr)*s[1]*r[2]/2.0

#define DGDS2  0.5*R*t*( log(xmn2a) - log(xti4a) - log(xmn2b) + log(xti4b) ) \
                + (whilmpyr-whilmpyrT)*s[0]/2.0 + (whgeipyr-whgeipyrT)*s[1]/2.0 \
                - 2.0*(dhpyr+(p-1.0)*dvpyr+((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
	    			         -(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[3]/2.0)*s[2] \
                + 2.0*(whpyr+(p-1.0)*wvpyr)*s[2]*(r[2]*r[2]-2.0*s[2]*s[2]) \
                - (whilmpyrpyr)*s[2]*r[0]/2.0 - (whgeipyrpyr)*s[2]*r[1]/2.0

#define DGDT   - (S) - fSRO(t,1)*2.0*R*t*(  xfe3ID*log(xfe3ID) + xal3ID*log(xal3ID) + xfe2ID*log(xfe2ID) \
                + xmg2ID*log(xmg2ID) + xmn2ID*log(xmn2ID) + xti4ID*log(xti4ID) ) - fSRO(t,1)*2.0*R*t*log(2.0)

#define DGDP     (wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                + (wvhmgei+dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
               + (wvhmpyr+dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                +  wvilmgei*r[0]*r[1] \
                +  wvilmpyr*r[0]*r[2] \
                +  wvgeipyr*r[1]*r[2] \
                +  dwvcrnilm*(r[0]-r[3])*r[0]*r[3] \
                +  dwvcrngei*(r[1]-r[3])*r[1]*r[3] \
                +  dwvcrnpyr*(r[2]-r[3])*r[2]*r[3] \
                + (dvilm + (dwvhmilm/2.0+wvhmilm2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) - dwvcrnilm*r[3]/2.0)*(r[0]*r[0]-s[0]*s[0]) \
                + (dvgei + (dwvhmgei/2.0+wvhmgei2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) - dwvcrngei*r[3]/2.0)*(r[1]*r[1]-s[1]*s[1]) \
                + (dvpyr + (dwvhmpyr/2.0+wvhmpyr2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) - dwvcrnpyr*r[3]/2.0)*(r[2]*r[2]-s[2]*s[2]) \
                +  wvilm*s[0]*s[0]*(r[0]*r[0]-s[0]*s[0]) \
                +  wvgei*s[1]*s[1]*(r[1]*r[1]-s[1]*s[1]) \
                +  wvpyr*s[2]*s[2]*(r[2]*r[2]-s[2]*s[2])

#define D2GDR0R0  R*t*(  0.25/xfe2a + 0.25/xti4a + 0.25/xfe2b + 0.25/xti4b - 0.5/xfe2ID - 0.5/xti4ID ) \
                                    + (1.0-fSRO(t,0))*2.0*R*t*( 1.0/xfe3ID + 0.25/xfe2ID + 0.25/xti4ID ) \
                                    - 2.0*((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3])) \
                                    - 2.0*(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
	          + (dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
                                    + (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
	          - 2.0*(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
	          +	(dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
                                    +	(dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
	          + 2.0*(dwhcrnilm+(p-1.0)*dwvcrnilm)*r[3] \
	          + 2.0*(dhilm+(p-1.0)*dvilm+((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
                                -(dwhcrnilm+(p-1.0)*dwvcrnilm)*r[3]/2.0) \
	          - 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*r[0] \
	          + 2.0*(whilm+(p-1.0)*wvilm)*s[0]*s[0] \
	          + (whilmilmgei)*r[1]/2.0 + (whilmilmpyr)*r[2]/2.0

#define D2GDR0R1  R*t*( 0.25/xti4a + 0.25/xti4b - 0.5/xti4ID ) \
                                    + (1.0-fSRO(t,0))*2.0*R*t*( 1.0/xfe3ID + 0.25/xti4ID ) \
                                    - ((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3])) \
                                    - (dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
	          - ((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-4.0*r[1]-r[2]-r[3])) \
                                    + (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
	          + 2.0*(dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
	          -	(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
                                    +	(dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
	          + (whilmgei+(p-1.0)*wvilmgei+whilmgeiT)/2.0 \
	          - 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*r[0] \
	          - 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*r[1] \
	          + (whilmilmgei)*r[0]/2.0 + (whilmgeigei)*r[1]/2.0

#define D2GDR0R2  R*t*( 0.25/xti4a + 0.25/xti4b - 0.5/xti4ID ) \
                                    + (1.0-fSRO(t,0))*2.0*R*t*( 1.0/xfe3ID + 0.25/xti4ID ) \
                                    - ((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3])) \
                                    - (dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
	          + (dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
                                    - ((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-4.0*r[2]-r[3])) \
	          + 2.0*(dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
	          +	(dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
                                    -	(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
	          + (whilmpyr+(p-1.0)*wvilmpyr+whilmpyrT)/2.0 \
	          - 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*r[0] \
	          - 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*r[2] \
	          + (whilmilmpyr)*r[0]/2.0 + (whilmpyrpyr)*r[2]/2.0

#define D2GDR0R3  (1.0-fSRO(t,0))*2.0*R*t*( 1.0/xfe3ID ) \
                                    - ((whhmilm+(p-1.0)*wvhmilm)+(dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3])) \
                                    - (dwhhmilm+(p-1.0)*dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
	          + (dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
                                    + (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
	          + 2.0*(dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
	          +	(dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
                                    +	(dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
	          -  whmcrn \
	          + (wcrnilm+(dwhcrnilm+(p-1.0)*dwvcrnilm)*(r[0]-2.0*r[3])) + (dwhcrnilm+(p-1.0)*dwvcrnilm)*r[0] \
	          - 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0+(dwhcrnilm+(p-1.0)*dwvcrnilm)/2.0)*r[0]

#define D2GDR0S0  0.25*R*t*( 1.0/xfe2a - 1.0/xti4a - 1.0/xfe2b + 1.0/xti4b ) \
	          + 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*s[0] \
	          + 4.0*(whilm+(p-1.0)*wvilm)*s[0]*r[0]

#define D2GDR0S1  0.25*R*t*( - 1.0/xti4a + 1.0/xti4b ) \
	          + 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*s[1] \
	          - (whilmgeigei)*s[1]/2.0

#define D2GDR0S2  0.25*R*t*( - 1.0/xti4a + 1.0/xti4b ) \
	          + 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*s[2] \
	          - (whilmpyrpyr)*s[2]/2.0

#define D2GDR0DT  R*(  0.5*log(xfe2a) + 0.5*log(xti4a) + 0.5*log(xfe2b) + 0.5*log(xti4b) - log(xfe2ID) - log(xti4ID) ) \
                                    + (1.0-fSRO(t,0))*2.0*R*(- log(xfe3ID) + 0.5*log(xfe2ID) + 0.5*log(xti4ID) ) \
        - fSRO(t,1)*2.0*R*t*(- log(xfe3ID) + 0.5*log(xfe2ID) + 0.5*log(xti4ID) )

#define D2GDR0DP (wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
                - (wvhmgei+dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*r[1] \
               - (wvhmpyr+dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*r[2] \
                - 2.0*dwvhmilm*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                -     dwvhmgei*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
               -     dwvhmpyr*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                + wvilmgei*r[1] \
                + wvilmpyr*r[2] \
                + dwvcrnilm*(r[0]-r[3])*r[3] + dwvcrnilm*r[0]*r[3] \
                + 2.0*(dvilm+(dwvhmilm/2.0+wvhmilm2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrnilm*r[3]/2.0)*r[0] \
                - (dwvhmilm/2.0+wvhmilm2/4.0)*(r[0]*r[0]-s[0]*s[0]) \
                - (dwvhmgei/2.0+wvhmgei2/4.0)*(r[1]*r[1]-s[1]*s[1]) \
                - (dwvhmpyr/2.0+wvhmpyr2/4.0)*(r[2]*r[2]-s[2]*s[2]) \
                + 2.0*wvilm*s[0]*s[0]*r[0]

#define D2GDR1R1  R*t*(  0.25/xmg2a + 0.25/xti4a + 0.25/xmg2b + 0.25/xti4b - 0.5/xmg2ID - 0.5/xti4ID ) \
               + (1.0-fSRO(t,0))*2.0*R*t*( 1.0/xfe3ID + 0.25/xmg2ID + 0.25/xti4ID ) \
               + (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
                - 2.0*((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3])) \
                - 2.0*(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
               + (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
                +     (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
                - 2.0*(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
               +     (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
                + 2.0*(dwhcrngei+(p-1.0)*dwvcrngei)*r[3] \
                + 2.0*(dhgei+(p-1.0)*dvgei+((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
	    			         -(dwhcrngei+(p-1.0)*dwvcrngei)*r[3]/2.0) \
                - 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*r[1] \
                + 2.0*(whgei+(p-1.0)*wvgei)*s[1]*s[1] \
                + (whilmgeigei)*r[0]/2.0 + (whgeigeipyr)*r[2]/2.0

#define D2GDR1R2  R*t*(  0.25/xti4a + 0.25/xti4b - 0.5/xti4ID ) \
               + (1.0-fSRO(t,0))*2.0*R*t*( 1.0/xfe3ID + 0.25/xti4ID ) \
               + (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
                - ((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3])) \
                - (dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
               - ((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-4.0*r[2]-r[3])) \
                +     (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
                + 2.0*(dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
               -     (dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
                + (whgeipyr+(p-1.0)*wvgeipyr+whgeipyrT)/2.0 \
                - 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*r[1] \
                - 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*r[2] \
                + (whgeigeipyr)*r[1]/2.0 + (whgeipyrpyr)*r[2]/2.0

#define D2GDR1R3  (1.0-fSRO(t,0))*2.0*R*t*( 1.0/xfe3ID ) \
               + (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
                - ((whhmgei+(p-1.0)*wvhmgei)+(dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3])) \
                - (dwhhmgei+(p-1.0)*dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
               + (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
                +     (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
                + 2.0*(dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
               +     (dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
                -  whmcrn \
                + (wcrngei+(dwhcrngei+(p-1.0)*dwvcrngei)*(r[1]-2.0*r[3])) +(dwhcrngei+(p-1.0)*dwvcrngei)*r[1] \
                - 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0+(dwhcrngei+(p-1.0)*dwvcrngei)/2.0)*r[1]

#define D2GDR1S0  0.25*R*t*( - 1.0/xti4a + 1.0/xti4b ) \
                + 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*s[0] \
                - (whilmilmgei)*s[0]/2.0

#define D2GDR1S1  0.25*R*t*( 1.0/xmg2a - 1.0/xti4a - 1.0/xmg2b + 1.0/xti4b ) \
                + 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*s[1] \
                + 4.0*(whgei+(p-1.0)*wvgei)*s[1]*r[1]

#define D2GDR1S2  0.25*R*t*( - 1.0/xti4a + 1.0/xti4b ) \
                + 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*s[2] \
                - (whgeipyrpyr)*s[2]/2.0

#define D2GDR1DT  R*(  0.5*log(xmg2a) + 0.5*log(xti4a) + 0.5*log(xmg2b) + 0.5*log(xti4b) - log(xmg2ID) - log(xti4ID) ) \
               + (1.0-fSRO(t,0))*2.0*R*(- log(xfe3ID) + 0.5*log(xmg2ID) + 0.5*log(xti4ID) ) \
                - fSRO(t,1)*2.0*R*t*(- log(xfe3ID) + 0.5*log(xmg2ID) + 0.5*log(xti4ID) )

#define D2GDR1DP - (wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*r[0] \
                + (wvhmgei+dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
               - (wvhmpyr+dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*r[2] \
                -     dwvhmilm*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                - 2.0*dwvhmgei*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
               -     dwvhmpyr*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                + wvilmgei*r[0] \
                + wvgeipyr*r[2] \
                + dwvcrngei*(r[1]-r[3])*r[3] +dwvcrngei*r[1]*r[3] \
                - (dwvhmilm/2.0+wvhmilm2/4.0)*(r[0]*r[0]-s[0]*s[0]) \
                + 2.0*(dvgei+(dwvhmgei/2.0+wvhmgei2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrngei*r[3]/2.0)*r[1] \
                - (dwvhmgei/2.0+wvhmgei2/4.0)*(r[1]*r[1]-s[1]*s[1]) \
                - (dwvhmpyr/2.0+wvhmpyr2/4.0)*(r[2]*r[2]-s[2]*s[2]) \
                + 2.0*wvgei*s[1]*s[1]*r[1]

#define D2GDR2R2  R*t*(  0.25/xmn2a + 0.25/xti4a + 0.25/xmn2b + 0.25/xti4b - 0.5/xmn2ID - 0.5/xti4ID ) \
               + (1.0-fSRO(t,0))*2.0*R*t*( 1.0/xfe3ID + 0.25/xmn2ID + 0.25/xti4ID ) \
               + (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
                + (dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
               - 2.0*((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3])) \
               - 2.0*(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
               +     (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
                +     (dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
               - 2.0*(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
                + 2.0*(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[3] \
                + 2.0*(dhpyr+(p-1.0)*dvpyr+((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
                                        -(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[3]/2.0) \
                - 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*r[2] \
                + 2.0*(whpyr+(p-1.0)*wvpyr)*s[2]*s[2] \
                + (whilmpyrpyr)*r[0]/2.0 + (whgeipyrpyr)*r[1]/2.0

#define D2GDR2R3  (1.0-fSRO(t,0))*2.0*R*t*( 1.0/xfe3ID ) \
               + (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
                + (dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
               - ((whhmpyr+(p-1.0)*wvhmpyr)+(dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3])) \
               - (dwhhmpyr+(p-1.0)*dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
               +     (dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
                +     (dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
               + 2.0*(dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
                -  whmcrn \
                + (wcrnpyr+(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*(r[2]-2.0*r[3])) + (dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[2] \
                - 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0+(dwhcrnpyr+(p-1.0)*dwvcrnpyr)/2.0)*r[2]

#define D2GDR2S0  0.25*R*t*( - 1.0/xti4a + 1.0/xti4b ) \
                + 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*s[0] \
                - (whilmilmpyr)*s[0]/2.0

#define D2GDR2S1  0.25*R*t*( - 1.0/xti4a + 1.0/xti4b ) \
                + 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*s[1] \
                - (whgeigeipyr)*s[1]/2.0

#define D2GDR2S2  0.25*R*t*( 1.0/xmn2a - 1.0/xti4a - 1.0/xmn2b + 1.0/xti4b ) \
                + 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*s[2] \
                + 4.0*(whpyr+(p-1.0)*wvpyr)*s[2]*r[2]

#define D2GDR2DT  R*(  0.5*log(xmn2a) + 0.5*log(xti4a) + 0.5*log(xmn2b) + 0.5*log(xti4b) - log(xmn2ID) - log(xti4ID) ) \
               + (1.0-fSRO(t,0))*2.0*R*(- log(xfe3ID) + 0.5*log(xmn2ID) + 0.5*log(xti4ID) ) \
                - fSRO(t,1)*2.0*R*t*(- log(xfe3ID) + 0.5*log(xmn2ID) + 0.5*log(xti4ID) )

#define D2GDR2DP  - (wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*r[0] \
                - (wvhmgei+dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*r[1] \
               + (wvhmpyr+dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
               -     dwvhmilm*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                -     dwvhmgei*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
               - 2.0*dwvhmpyr*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                + wvilmpyr*r[0] \
                + wvgeipyr*r[1] \
                + dwvcrnpyr*(r[2]-r[3])*r[3] + dwvcrnpyr*r[2]*r[3] \
                - (dwvhmilm/2.0+wvhmilm2/4.0)*(r[0]*r[0]-s[0]*s[0]) \
                - (dwvhmgei/2.0+wvhmgei2/4.0)*(r[1]*r[1]-s[1]*s[1]) \
                + 2.0*(dvpyr+(dwvhmpyr/2.0+wvhmgei2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrnpyr*r[3]/2.0)*r[2] \
                - (dwvhmpyr/2.0+wvhmpyr2/4.0)*(r[2]*r[2]-s[2]*s[2]) \
                + 2.0*wvpyr*s[2]*s[2]*r[2]

#define D2GDR3R3  (1.0-fSRO(t,0))*2.0*R*t*( 1.0/xfe3ID + 1.0/xal3ID ) \
               + 2.0*(dwhhmilm+(p-1.0)*dwvhmilm)*r[0] \
                + 2.0*(dwhhmgei+(p-1.0)*dwvhmgei)*r[1] \
               + 2.0*(dwhhmpyr+(p-1.0)*dwvhmpyr)*r[2] \
                - 2.0*whmcrn \
                - 2.0*(dwhcrnilm+(p-1.0)*dwvcrnilm)*r[0] \
                - 2.0*(dwhcrngei+(p-1.0)*dwvcrngei)*r[1] \
                - 2.0*(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[2]

#define D2GDR3S0 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0 + (whhmilm2+(p-1.0)*wvhmilm2)/4.0 + (dwhcrnilm+(p-1.0)*dwvcrnilm)/2.0)*s[0]

#define D2GDR3S1 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0 + (whhmgei2+(p-1.0)*wvhmgei2)/4.0 + (dwhcrngei+(p-1.0)*dwvcrngei)/2.0)*s[1]

#define D2GDR3S2 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0 + (whhmpyr2+(p-1.0)*wvhmpyr2)/4.0 + (dwhcrnpyr+(p-1.0)*dwvcrnpyr)/2.0)*s[2]

#define D2GDR3DT  (1.0-fSRO(t,0))*2.0*R*(- log(xfe3ID) + log(xal3ID) ) - fSRO(t,1)*2.0*R*t*(- log(xfe3ID) + log(xal3ID) )

#define D2GDR3DP - (wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]))*r[0] \
                - (wvhmgei+dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]))*r[1] \
               - (wvhmpyr+dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]))*r[2] \
               - dwvhmilm*r[0]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                - dwvhmgei*r[1]*(1.0-r[0]-r[1]-r[2]-r[3]) \
               - dwvhmpyr*r[2]*(1.0-r[0]-r[1]-r[2]-r[3]) \
                + (dwvcrnilm*(r[0]-r[3]))*r[0] - dwvcrnilm*r[0]*r[3] \
                + (dwvcrngei*(r[1]-r[3]))*r[1] - dwvcrngei*r[1]*r[3] \
                + (dwvcrnpyr*(r[2]-r[3]))*r[2] - dwvcrnpyr*r[2]*r[3] \
                - (dwvhmilm/2.0+wvhmilm2/4.0+dwvcrnilm/2.0)*(r[0]*r[0]-s[0]*s[0]) \
                - (dwvhmgei/2.0+wvhmgei2/4.0+dwvcrngei/2.0)*(r[1]*r[1]-s[1]*s[1]) \
                - (dwvhmpyr/2.0+wvhmpyr2/4.0+dwvcrnpyr/2.0)*(r[2]*r[2]-s[2]*s[2])

#define D2GDS0S0  0.25*R*t*( 1.0/xfe2a + 1.0/xti4a + 1.0/xfe2b + 1.0/xti4b ) \
	          - 2.0*(dhilm+(p-1.0)*dvilm+((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
	    		   	            -(dwhcrnilm+(p-1.0)*dwvcrnilm)*r[3]/2.0) \
	          + 2.0*(whilm+(p-1.0)*wvilm)*(r[0]*r[0]-6.0*s[0]*s[0]) \
	          - (whilmilmgei)*r[1]/2.0 - (whilmilmpyr)*r[2]/2.0

#define D2GDS0S1  0.25*R*t*( 1.0/xti4a + 1.0/xti4b ) + (whilmgei-whilmgeiT)/2.0

#define D2GDS0S2  0.25*R*t*( 1.0/xti4a + 1.0/xti4b ) + (whilmpyr-whilmpyrT)/2.0

#define D2GDS0DT  0.5*R*( log(xfe2a) - log(xti4a) - log(xfe2b) + log(xti4b) )

#define D2GDS0DP - 2.0*(dvilm+(dwvhmilm/2.0+wvhmilm2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrnilm*r[3]/2.0)*s[0] \
                    + 2.0*wvilm*s[0]*(r[0]*r[0]-2.0*s[0]*s[0])

#define D2GDS1S1  0.25*R*t*( 1.0/xmg2a + 1.0/xti4a + 1.0/xmg2b + 1.0/xti4b ) \
	          - 2.0*(dhgei+(p-1.0)*dvgei+((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
                                        -(dwhcrngei+(p-1.0)*dwvcrngei)*r[3]/2.0) \
	          + 2.0*(whgei+(p-1.0)*wvgei)*(r[1]*r[1]-6.0*s[1]*s[1]) \
	          - (whilmgeigei)*r[0]/2.0 - (whgeigeipyr)*r[2]/2.0

#define D2GDS1S2  0.25*R*t*( 1.0/xti4a + 1.0/xti4b ) + (whgeipyr-whgeipyrT)/2.0

#define D2GDS1DT  0.5*R*( log(xmg2a) - log(xti4a) - log(xmg2b) + log(xti4b) )

#define D2GDS1DP - 2.0*(dvgei+(dwvhmgei/2.0+wvhmgei2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrngei*r[3]/2.0)*s[1] \
                    + 2.0*wvgei*s[1]*(r[1]*r[1]-2.0*s[1]*s[1])

#define D2GDS2S2  0.25*R*t*( 1.0/xmn2a + 1.0/xti4a + 1.0/xmn2b + 1.0/xti4b ) \
	          - 2.0*(dhpyr+(p-1.0)*dvpyr+((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0)*(1.0-r[0]-r[1]-r[2]-r[3]) \
                                        -(dwhcrnpyr+(p-1.0)*dwvcrnpyr)*r[3]/2.0) \
	          + 2.0*(whpyr+(p-1.0)*wvpyr)*(r[2]*r[2]-6.0*s[2]*s[2]) \
	          - (whilmpyrpyr)*r[0]/2.0 - (whgeipyrpyr)*r[1]/2.0

#define D2GDS2DT  0.5*R*( log(xmn2a) - log(xti4a) - log(xmn2b) + log(xti4b) )

#define D2GDS2DP - 2.0*(dvpyr+(dwvhmpyr/2.0+wvhmpyr2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrnpyr*r[3]/2.0)*s[2] \
                    + 2.0*wvpyr*s[2]*(r[2]*r[2]-2.0*s[2]*s[2])

#define D2GDT2     - fSRO(t,1)*4.0*R*(  xfe3ID*log(xfe3ID) + xal3ID*log(xal3ID) + xfe2ID*log(xfe2ID) \
	                              + xmg2ID*log(xmg2ID) + xmn2ID*log(xmn2ID) + xti4ID*log(xti4ID) ) \
		   - fSRO(t,1)*4.0*R*log(2.0) \
                        - fSRO(t,2)*2.0*R*t*(  xfe3ID*log(xfe3ID) + xal3ID*log(xal3ID) + xfe2ID*log(xfe2ID) \
                        + xmg2ID*log(xmg2ID) + xmn2ID*log(xmn2ID) + xti4ID*log(xti4ID) ) - fSRO(t,2)*2.0*R*t*log(2.0)
#define D2GDTDP   0.0
#define D2GDP2    0.0

/*----------------------------------------------------------------------------*/
/* Excess free energy modifications incomplete below this point               */
/*----------------------------------------------------------------------------*/


#define D3GDR0R0R0 R*t*(- 0.125/SQUARE(xfe2a) - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xfe2b) - 0.125/SQUARE(xti4b) \
                                                + 0.25/SQUARE(xfe2ID) + 0.25/SQUARE(xti4ID) ) \
                   + (1.0-fSRO(t,0))*2.0*R*t*(1.0/SQUARE(xfe3ID) - 0.125/SQUARE(xfe2ID) - 0.125/SQUARE(xti4ID) ) \
                   + 12.0*(dwhhmilm+(p-1.0)*dwvhmilm) \
                        - 6.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0)

#define D3GDR0R0R1 R*t*(  0.125/SQUARE(xti4a) + 0.125/SQUARE(xti4b) - 0.25/SQUARE(xti4ID) ) \
                   + (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 0.125/SQUARE(xti4ID) )
#define D3GDR0R0R2 R*t*(  0.125/SQUARE(xti4a) + 0.125/SQUARE(xti4b) - 0.25/SQUARE(xti4ID) ) \
                   + (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 0.125/SQUARE(xti4ID) )
#define D3GDR0R0R3 (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR0R0S0 R*t*( - 0.125/SQUARE(xfe2a) + 0.125/SQUARE(xti4a) + 0.125/SQUARE(xfe2b) - 0.125/SQUARE(xti4b) )  \
                        + 4.0*(whilm+(p-1.0)*wvilm)*s[0]
#define D3GDR0R0S1 R*t*( 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR0R0S2 R*t*( 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR0R0DT R*(  0.25/xfe2a + 0.25/xti4a + 0.25/xfe2b + 0.25/xti4b - 0.5/xfe2ID - 0.5/xti4ID ) \
                   + (1.0-fSRO(t,0))*2.0*R*( 1.0/xfe3ID + 0.25/xfe2ID + 0.25/xti4ID ) \
		   - fSRO(t,1)*2.0*R*t*( 1.0/xfe3ID + 0.25/xfe2ID + 0.25/xti4ID )
#define D3GDR0R0DP - 2.0*(wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3])) \
                   - 2.0*dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
                        + dwvhmgei*r[1] \
                   + dwvhmpyr*r[2] \
                        - 2.0*dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
                        +	dwvhmgei*r[1] \
                   +	dwvhmpyr*r[2] \
                        + 2.0*dwvcrnilm*r[3] \
                        + 2.0*(dvilm+(dwvhmilm/2.0+wvhmilm2/4.0)*(1.0-2.0*r[0]-r[1]-r[2]-r[3])-dwvcrnilm*r[3]/2.0) \
                        - 2.0*(dwvhmilm/2.0+wvhmilm2/4.0)*r[0] \
                        + 2.0*wvilm*s[0]*s[0]
/*end-OK*/
#define D3GDR0R1R1 R*t*( 0.125/SQUARE(xti4a) + 0.125/SQUARE(xti4b) - 0.25/SQUARE(xti4ID) ) \
                   + (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 0.125/SQUARE(xti4ID) )
#define D3GDR0R1R2 R*t*( 0.125/SQUARE(xti4a) + 0.125/SQUARE(xti4b) - 0.25/SQUARE(xti4ID) ) \
                   + (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 0.125/SQUARE(xti4ID) )
#define D3GDR0R1R3 (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR0R1S0 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR0R1S1 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR0R1S2 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR0R1DT R*( 0.25/xti4a + 0.25/xti4b - 0.5/xti4ID ) \
                   + (1.0-fSRO(t,0))*2.0*R*( 1.0/xfe3ID + 0.25/xti4ID ) \
		   - fSRO(t,1)*2.0*R*t*( 1.0/xfe3ID + 0.25/xti4ID )
#define D3GDR0R1DP - (wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3])) \
                   - (dwvhmilm)*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
                        - (wvhmgei+dwvhmgei*(1.0-r[0]-4.0*r[1]-r[2]-r[3])) \
                   + dwvhmpyr*r[2] \
                        + 2.0*dwvhmilm*r[0] \
                        -	dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
                   +	dwvhmpyr*r[2] \
                        + wvilmgei \
                        - 2.0*(dwvhmilm/2.0+wvhmilm2/4.0)*r[0] \
                        - 2.0*(dwvhmgei/2.0+wvhmgei2/4.0)*r[1]
/*end-OK*/
#define D3GDR0R2R2 R*t*( 0.125/SQUARE(xti4a) + 0.125/SQUARE(xti4b) - 0.25/SQUARE(xti4ID) ) \
                   + (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 0.125/SQUARE(xti4ID) )
#define D3GDR0R2R3 (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR0R2S0 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR0R2S1 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR0R2S2 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR0R2DT R*( 0.25/xti4a + 0.25/xti4b - 0.5/xti4ID ) \
                   + (1.0-fSRO(t,0))*2.0*R*( 1.0/xfe3ID + 0.25/xti4ID ) \
		   - fSRO(t,1)*2.0*R*t*( 1.0/xfe3ID + 0.25/xti4ID )
#define D3GDR0R2DP - (wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3])) \
                   - dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
                        + dwvhmgei*r[1] \
                   - (wvhmpyr+dwvhmpyr*(1.0-r[0]-r[1]-4.0*r[2]-r[3])) \
                        + 2.0*dwvhmilm*r[0] \
                        +	dwvhmgei*r[1] \
                   -	dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
                        + wvilmpyr \
                        - 2.0*(dwvhmilm/2.0+wvhmilm2/4.0)*r[0] \
                        - 2.0*(dwvhmpyr/2.0+wvhmpyr2/4.0)*r[2]
/*end-OK*/
#define D3GDR0R3R3 (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR0R3S0 0.0
#define D3GDR0R3S1 0.0
#define D3GDR0R3S2 0.0
#define D3GDR0R3DT (1.0-fSRO(t,0))*2.0*R*( 1.0/xfe3ID ) -fSRO(t,1)*2.0*R*t*( 1.0/xfe3ID )
#define D3GDR0R3DP - (wvhmilm+dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3])) \
                   - dwvhmilm*(1.0-2.0*r[0]-r[1]-r[2]-r[3]) \
                        + dwvhmgei*r[1] \
                   + dwvhmpyr*r[2] \
                        + 2.0*dwvhmilm*r[0] \
	          +	dwvhmgei*r[1] \
                                    +	dwvhmpyr*r[2] \
	          + (dwvcrnilm*(r[0]-2.0*r[3])) + dwvcrnilm*r[0] \
	          - 2.0*(dwvhmilm/2.0+wvhmilm2/4.0+dwvcrnilm/2.0)*r[0]

#define D3GDR0S0S0 R*t*( - 0.125/SQUARE(xfe2a) - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xfe2b) - 0.125/SQUARE(xti4b) ) \
                   + 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0) + 4.0*(whilm+(p-1.0)*wvilm)*r[0]
#define D3GDR0S0S1 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR0S0S2 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR0S0DT R*( 0.25/xfe2a - 0.25/xti4a - 0.25/xfe2b + 0.25/xti4b )
#define D3GDR0S0DP 2.0*(dwvhmilm/2.0+wvhmilm2/4.0)*s[0] + 4.0*wvilm*s[0]*r[0]

#define D3GDR0S1S1 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) ) \
                   + 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0) - (whilmgeigei)/2.0
#define D3GDR0S1S2 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR0S1DT R*( - 0.25/xti4a + 0.25/xti4b )
#define D3GDR0S1DP 2.0*(dwvhmgei/2.0+wvhmgei2/4.0)*s[1]

#define D3GDR0S2S2 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) ) \
                   + 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0) - (whilmpyrpyr)/2.0
#define D3GDR0S2DT R*( - 0.25/xti4a + 0.25/xti4b )
#define D3GDR0S2DP 2.0*(dwvhmpyr/2.0+wvhmpyr2/4.0)*s[2]

#define D3GDR0DT2  -fSRO(t,1)*4.0*R*(- log(xfe3ID) + 0.5*log(xfe2ID) + 0.5*log(xti4ID) ) \
		   -fSRO(t,2)*2.0*R*t*(- log(xfe3ID) + 0.5*log(xfe2ID) + 0.5*log(xti4ID) )
#define D3GDR0DTDP 0.0
#define D3GDR0DP2  0.0
/*end-OK*/
#define D3GDR1R1R1 R*t*(  0.125/SQUARE(xmg2a) + 0.125/SQUARE(xti4a) + 0.125/SQUARE(xmg2b) + 0.125/SQUARE(xti4b) \
                                                - 0.25/SQUARE(xmg2ID) - 0.25/SQUARE(xti4ID) ) \
                   + (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 0.125/SQUARE(xmg2ID) + 0.125/SQUARE(xti4ID) )
#define D3GDR1R1R2 R*t*(  0.125/SQUARE(xti4a) + 0.125/SQUARE(xti4b) - 0.25/SQUARE(xti4ID) ) \
                   + (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID)  + 0.125/SQUARE(xti4ID) )
#define D3GDR1R1R3 (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR1R1S0 R*t*(   0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR1R1S1 R*t*( - 0.125/SQUARE(xmg2a) + 0.125/SQUARE(xti4a) + 0.125/SQUARE(xmg2b) - 0.125/SQUARE(xti4b) ) \
                   + 4.0*(whgei+(p-1.0)*wvgei)*s[1]
#define D3GDR1R1S2 R*t*(   0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR1R1DT R*(  0.25/xmg2a + 0.25/xti4a + 0.25/xmg2b + 0.25/xti4b - 0.5/xmg2ID - 0.5/xti4ID ) \
                   + (1.0-fSRO(t,0))*2.0*R*( 1.0/xfe3ID + 0.25/xmg2ID + 0.25/xti4ID ) \
		   - fSRO(t,1)*2.0*R*t*( 1.0/xfe3ID + 0.25/xmg2ID + 0.25/xti4ID )
#define D3GDR1R1DP dwvhmilm*r[0] \
                        - 2.0*(wvhmgei+dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3])) \
                        - 2.0*dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
                   + dwvhmpyr*r[2] \
                        + dwvhmilm*r[0] \
                        - 2.0*dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
                   + dwvhmpyr*r[2] \
                        + 2.0*dwvcrngei*r[3] \
                        + 2.0*(dvgei+(dwvhmgei/2.0+wvhmgei2/4.0)*(1.0-r[0]-2.0*r[1]-r[2]-r[3])-dwvcrngei*r[3]/2.0) \
                        - 2.0*(dwvhmgei/2.0+wvhmgei2/4.0)*r[1] \
                        + 2.0*wvgei*s[1]*s[1]
/*end-OK*/
#define D3GDR1R2R2 R*t*(  0.125/SQUARE(xti4a) + 0.125/SQUARE(xti4b) - 0.25/SQUARE(xti4ID) ) \
                   + (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 0.125/SQUARE(xti4ID) )
#define D3GDR1R2R3 (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR1R2S0 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR1R2S1 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR1R2S2 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR1R2DT R*(  0.25/xti4a + 0.25/xti4b - 0.5/xti4ID ) \
                   + (1.0-fSRO(t,0))*2.0*R*( 1.0/xfe3ID + 0.25/xti4ID ) \
		   - fSRO(t,1)*2.0*R*t*( 1.0/xfe3ID + 0.25/xti4ID )
#define D3GDR1R2DP (dwvhmilm)*r[0] \
                - ((wvhmgei)+(dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3])) \
                - (dwvhmgei)*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
               - ((wvhmpyr)+(dwvhmpyr)*(1.0-r[0]-r[1]-4.0*r[2]-r[3])) \
                +     (dwvhmilm)*r[0] \
                + 2.0*(dwvhmgei)*r[1] \
               -     (dwvhmpyr)*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
                + (wvgeipyr) \
                - 2.0*((dwvhmgei)/2.0+(wvhmgei2)/4.0)*r[1] \
                - 2.0*((dwvhmpyr)/2.0+(wvhmpyr2)/4.0)*r[2]
/*end-OK*/
#define D3GDR1R3R3 (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR1R3S0 0.0
#define D3GDR1R3S1 0.0
#define D3GDR1R3S2 0.0
#define D3GDR1R3DT (1.0-fSRO(t,0))*2.0*R*( 1.0/xfe3ID ) - fSRO(t,1)*2.0*R*t*( 1.0/xfe3ID )
#define D3GDR1R3DP dwvhmilm*r[0] \
                        - (wvhmgei+dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3])) \
                        - dwvhmgei*(1.0-r[0]-2.0*r[1]-r[2]-r[3]) \
                   + dwvhmpyr*r[2] \
                        +     dwvhmilm*r[0] \
                        + 2.0*dwvhmgei*r[1] \
                   +     dwvhmpyr*r[2] \
                        + dwvcrngei*(r[1]-2.0*r[3])+dwvcrngei*r[1] \
                        - 2.0*(dwvhmgei/2.0 + wvhmgei2/4.0 + dwvcrngei/2.0)*r[1]

#define D3GDR1S0S0 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) ) \
                   + 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0) - (whilmilmgei)/2.0
#define D3GDR1S0S1 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR1S0S2 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR1S0DT 0.25*R*( - 1.0/xti4a + 1.0/xti4b )
#define D3GDR1S0DP 2.0*(dwvhmilm/2.0+wvhmilm2/4.0)*s[0]

#define D3GDR1S1S1 R*t*( - 0.125/SQUARE(xmg2a) - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xmg2b) - 0.125/SQUARE(xti4b) ) \
                   + 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0) + 4.0*(whgei+(p-1.0)*wvgei)*r[1]
#define D3GDR1S1S2 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) )
#define D3GDR1S1DT R*( 0.25/xmg2a - 0.25/xti4a - 0.25/xmg2b + 0.25/xti4b )
#define D3GDR1S1DP 2.0*(dwvhmgei/2.0+wvhmgei2/4.0)*s[1] + 4.0*wvgei*s[1]*r[1]

#define D3GDR1S2S2 R*t*( - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xti4b) ) \
                   + 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0) - (whgeipyrpyr)/2.0
#define D3GDR1S2DT 0.25*R*( - 1.0/xti4a + 1.0/xti4b )
#define D3GDR1S2DP 2.0*(dwvhmpyr/2.0+wvhmpyr2/4.0)*s[2]

#define D3GDR1DT2  -fSRO(t,1)*4.0*R*(- log(xfe3ID) + 0.5*log(xmg2ID) + 0.5*log(xti4ID) ) \
                        -fSRO(t,2)*2.0*R*t*(- log(xfe3ID) + 0.5*log(xmg2ID) + 0.5*log(xti4ID) )
#define D3GDR1DTDP 0.0
#define D3GDR1DP2  0.0
/*end-OK*/
#define D3GDR2R2R2 R*t*(  0.125/SQUARE(xmn2a) + 0.125/SQUARE(xti4a) + 0.125/SQUARE(xmn2b) + 0.125/SQUARE(xti4b) \
                                                - 0.25/SQUARE(xmn2ID) - 0.25/SQUARE(xti4ID) ) \
                   + (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 0.125/SQUARE(xmn2ID) + 0.125/SQUARE(xti4ID) )
#define D3GDR2R2R3 (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR2R2S0 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR2R2S1 - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR2R2S2 - 0.125*R*t*( 1.0/SQUARE(xmn2a) - 1.0/SQUARE(xti4a) - 1.0/SQUARE(xmn2b) + 1.0/SQUARE(xti4b) ) \
                   + 4.0*(whpyr+(p-1.0)*wvpyr)*s[2]
#define D3GDR2R2DT R*(  0.25/xmn2a + 0.25/xti4a + 0.25/xmn2b + 0.25/xti4b - 0.5/xmn2ID - 0.5/xti4ID ) \
                   + (1.0-fSRO(t,0))*2.0*R*( 1.0/xfe3ID + 0.25/xmn2ID + 0.25/xti4ID ) \
		   - fSRO(t,1)*2.0*R*t*( 1.0/xfe3ID + 0.25/xmn2ID + 0.25/xti4ID )
#define D3GDR2R2DP dwvhmilm*r[0] \
                        + dwvhmgei*r[1] \
                   - 2.0*(wvhmpyr+dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3])) \
                   - 2.0*dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
                   +     dwvhmilm*r[0] \
                        +     dwvhmgei*r[1] \
                   - 2.0*dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
                        + 2.0*dwvcrnpyr*r[3] \
                        + 2.0*(dvpyr+(dwvhmpyr/2.0+wvhmpyr2/4.0)*(1.0-r[0]-r[1]-2.0*r[2]-r[3])-dwvcrnpyr*r[3]/2.0) \
                        - 2.0*(dwvhmpyr/2.0+wvhmpyr2/4.0)*r[2] \
                        + 2.0*wvpyr*s[2]*s[2]
/*end-OK*/
#define D3GDR2R3R3 (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) )
/*OK*/
#define D3GDR2R3S0 0.0
#define D3GDR2R3S1 0.0
#define D3GDR2R3S2 0.0
#define D3GDR2R3DT (1.0-fSRO(t,0))*2.0*R*( 1.0/xfe3ID ) - fSRO(t,1)*2.0*R*t*( 1.0/xfe3ID )
#define D3GDR2R3DP dwvhmilm*r[0] \
                        + dwvhmgei*r[1] \
                   - (wvhmpyr+dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3])) \
                   - dwvhmpyr*(1.0-r[0]-r[1]-2.0*r[2]-r[3]) \
                   +     dwvhmilm*r[0] \
                        +     dwvhmgei*r[1] \
                   + 2.0*dwvhmpyr*r[2] \
                        + dwvcrnpyr*(r[2]-2.0*r[3]) + dwvcrnpyr*r[2] \
                        - 2.0*(dwvhmpyr/2.0+wvhmpyr2/4.0+dwvcrnpyr/2.0)*r[2]

#define D3GDR2S0S0 - 0.125*R*t*( 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) ) \
                   + 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0+(whhmilm2+(p-1.0)*wvhmilm2)/4.0) - (whilmilmpyr)/2.0
#define D3GDR2S0S1 - 0.125*R*t*( 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR2S0S2 - 0.125*R*t*( 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR2S0DT R*( - 0.25/xti4a + 0.25/xti4b )
#define D3GDR2S0DP 2.0*(dwvhmilm/2.0+wvhmilm2/4.0)*s[0]

#define D3GDR2S1S1 - 0.125*R*t*( 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) ) \
                   + 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0+(whhmgei2+(p-1.0)*wvhmgei2)/4.0) - (whgeigeipyr)/2.0
#define D3GDR2S1S2 - 0.125*R*t*( 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDR2S1DT 0.25*R*( - 1.0/xti4a + 1.0/xti4b )
#define D3GDR2S1DP 2.0*(dwvhmgei/2.0+wvhmgei2/4.0)*s[1]

#define D3GDR2S2S2 R*t*( - 0.125/SQUARE(xmn2a) - 0.125/SQUARE(xti4a) - 0.125/SQUARE(xmn2b) - 0.125/SQUARE(xti4b) ) \
                   + 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0+(whhmpyr2+(p-1.0)*wvhmpyr2)/4.0) + 4.0*(whpyr+(p-1.0)*wvpyr)*r[2]
#define D3GDR2S2DT 0.25*R*( 1.0/xmn2a - 1.0/xti4a - 1.0/xmn2b + 1.0/xti4b )
#define D3GDR2S2DP 2.0*(dwvhmpyr/2.0+wvhmpyr2/4.0)*s[2] + 4.0*wvpyr*s[2]*r[2]

#define D3GDR2DT2  -fSRO(t,1)*4.0*R*(- log(xfe3ID) + 0.5*log(xmn2ID) + 0.5*log(xti4ID) ) \
                        -fSRO(t,2)*2.0*R*t*(- log(xfe3ID) + 0.5*log(xmn2ID) + 0.5*log(xti4ID) )
#define D3GDR2DTDP 0.0
#define D3GDR2DP2  0.0

#define D3GDR3R3R3 (1.0-fSRO(t,0))*2.0*R*t*( -1.0/SQUARE(xfe3ID) + 1.0/SQUARE(xal3ID) )
#define D3GDR3R3S0 0.0
#define D3GDR3R3S1 0.0
#define D3GDR3R3S2 0.0
#define D3GDR3R3DT (1.0-fSRO(t,0))*2.0*R*( 1.0/xfe3ID + 1.0/xal3ID ) - fSRO(t,1)*2.0*R*t*( 1.0/xfe3ID + 1.0/xal3ID )
#define D3GDR3R3DP   2.0*dwvhmilm*r[0] \
                        + 2.0*dwvhmgei*r[1] \
                   + 2.0*dwvhmpyr*r[2] \
                        - 2.0*dwvcrnilm*r[0] \
                        - 2.0*dwvcrngei*r[1] \
                        - 2.0*dwvcrnpyr*r[2]

#define D3GDR3S0S0 2.0*((dwhhmilm+(p-1.0)*dwvhmilm)/2.0 + (whhmilm2+(p-1.0)*wvhmilm2)/4.0 + (dwhcrnilm+(p-1.0)*dwvcrnilm)/2.0)
#define D3GDR3S0S1 0.0
#define D3GDR3S0S2 0.0
#define D3GDR3S0DT 0.0
#define D3GDR3S0DP 2.0*(dwvhmilm/2.0 + wvhmilm2/4.0 + dwvcrnilm/2.0)*s[0]

#define D3GDR3S1S1 2.0*((dwhhmgei+(p-1.0)*dwvhmgei)/2.0 + (whhmgei2+(p-1.0)*wvhmgei2)/4.0 + (dwhcrngei+(p-1.0)*dwvcrngei)/2.0)
#define D3GDR3S1S2 0.0
#define D3GDR3S1DT 0.0
#define D3GDR3S1DP 2.0*(dwvhmgei/2.0 + wvhmgei2/4.0 + dwvcrngei/2.0)*s[1]

#define D3GDR3S2S2 2.0*((dwhhmpyr+(p-1.0)*dwvhmpyr)/2.0 + (whhmpyr2+(p-1.0)*wvhmpyr2)/4.0 + (dwhcrnpyr+(p-1.0)*dwvcrnpyr)/2.0)
#define D3GDR3S2DT 0.0
#define D3GDR3S2DP 2.0*(dwvhmpyr/2.0 + wvhmpyr2/4.0 + dwvcrnpyr/2.0)*s[2]

#define D3GDR3DT2  -fSRO(t,1)*4.0*R*(- log(xfe3ID) + log(xal3ID) ) - fSRO(t,2)*2.0*R*t*(- log(xfe3ID) + log(xal3ID) )
#define D3GDR3DTDP 0.0
#define D3GDR3DP2  0.0

/*----------------------------------------------------------------------------*/
/* Excess free energy modifications are complete below this point             */
/*----------------------------------------------------------------------------*/

#define D3GDS0S0S0  - 0.125*R*t*(   1.0/SQUARE(xfe2a) - 1.0/SQUARE(xti4a) - 1.0/SQUARE(xfe2b) + 1.0/SQUARE(xti4b) ) \
	            - 24.0*(whilm+(p-1.0)*wvilm)*s[0]
#define D3GDS0S0S1  - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDS0S0S2  - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDS0S0DT  0.25*R*( 1.0/xfe2a + 1.0/xti4a + 1.0/xfe2b + 1.0/xti4b )
#define D3GDS0S0DP  - 2.0*(dvilm+(dwvhmilm/2.0+wvhmilm2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrnilm*r[3]/2.0) + 2.0*wvilm*(r[0]*r[0]-6.0*s[0]*s[0])
#define D3GDS0S1S1  - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDS0S1S2  - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDS0S1DT  0.25*R*( 1.0/xti4a + 1.0/xti4b )
#define D3GDS0S1DP  0.0
#define D3GDS0S2S2  - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDS0S2DT  0.25*R*( 1.0/xti4a + 1.0/xti4b )
#define D3GDS0S2DP  0.0
#define D3GDS0DT2   0.0
#define D3GDS0DTDP  0.0
#define D3GDS0DP2   0.0

#define D3GDS1S1S1  - 0.125*R*t*( 1.0/SQUARE(xmg2a) - 1.0/SQUARE(xti4a) - 1.0/SQUARE(xmg2b) + 1.0/SQUARE(xti4b) ) \
	            - 24.0*(whgei+(p-1.0)*wvgei)*s[1]
#define D3GDS1S1S2  - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDS1S1DT  0.25*R*( 1.0/xmg2a + 1.0/xti4a + 1.0/xmg2b + 1.0/xti4b )
#define D3GDS1S1DP  - 2.0*(dvgei+(dwvhmgei/2.0+wvhmgei2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrngei*r[3]/2.0) + 2.0*wvgei*(r[1]*r[1]-6.0*s[1]*s[1])
#define D3GDS1S2S2  - 0.125*R*t*( - 1.0/SQUARE(xti4a) + 1.0/SQUARE(xti4b) )
#define D3GDS1S2DT  0.25*R*( 1.0/xti4a + 1.0/xti4b )
#define D3GDS1S2DP  0.0
#define D3GDS1DT2   0.0
#define D3GDS1DTDP  0.0
#define D3GDS1DP2   0.0

#define D3GDS2S2S2  - 0.125*R*t*( 1.0/SQUARE(xmn2a) - 1.0/SQUARE(xti4a) - 1.0/SQUARE(xmn2b) + 1.0/SQUARE(xti4b) ) \
	            - 24.0*(whpyr+(p-1.0)*wvpyr)*s[2]
#define D3GDS2S2DT  0.25*R*( 1.0/xmn2a + 1.0/xti4a + 1.0/xmn2b + 1.0/xti4b )
#define D3GDS2S2DP  - 2.0*(dvpyr+(dwvhmpyr/2.0+wvhmpyr2/4.0)*(1.0-r[0]-r[1]-r[2]-r[3])-dwvcrnpyr*r[3]/2.0) + 2.0*wvpyr*(r[2]*r[2]-6.0*s[2]*s[2])
#define D3GDS2DT2   0.0
#define D3GDS2DTDP  0.0
#define D3GDS2DP2   0.0

#define D3GDT3      - fSRO(t,2)*6.0*R*(  xfe3ID*log(xfe3ID) + xal3ID*log(xal3ID) + xfe2ID*log(xfe2ID) \
                                                                + xmg2ID*log(xmg2ID) + xmn2ID*log(xmn2ID) + xti4ID*log(xti4ID) ) \
            - fSRO(t,2)*6.0*R*log(2.0) \
	            - fSRO(t,3)*2.0*R*t*(  xfe3ID*log(xfe3ID) + xal3ID*log(xal3ID) + xfe2ID*log(xfe2ID) \
                                                                    + xmg2ID*log(xmg2ID) + xmn2ID*log(xmn2ID) + xti4ID*log(xti4ID) ) \
            - fSRO(t,3)*2.0*R*t*log(2.0)
#define D3GDT2DP    0.0
#define D3GDTDP2    0.0
#define D3GDP3      0.0

/*
 *=============================================================================
 * Macros for automatic array initialization
 */
#define fillD2GDR2 \
 d2gdr2[0][0] = D2GDR0R0;     d2gdr2[0][1] = D2GDR0R1;      d2gdr2[0][2] = D2GDR0R2;     d2gdr2[0][3] = D2GDR0R3; \
 d2gdr2[1][0] = d2gdr2[0][1]; d2gdr2[1][1] = D2GDR1R1;      d2gdr2[1][2] = D2GDR1R2;     d2gdr2[1][3] = D2GDR1R3; \
 d2gdr2[2][0] = d2gdr2[0][2]; d2gdr2[2][1] = d2gdr2[1][2];  d2gdr2[2][2] = D2GDR2R2;     d2gdr2[2][3] = D2GDR2R3; \
 d2gdr2[3][0] = d2gdr2[0][3]; d2gdr2[3][1] = d2gdr2[1][3];  d2gdr2[3][2] = d2gdr2[2][3]; d2gdr2[3][3] = D2GDR3R3;

#define fillD2GDRDS \
 d2gdrds[0][0] = D2GDR0S0; d2gdrds[0][1] = D2GDR0S1; d2gdrds[0][2] = D2GDR0S2; \
 d2gdrds[1][0] = D2GDR1S0; d2gdrds[1][1] = D2GDR1S1; d2gdrds[1][2] = D2GDR1S2; \
 d2gdrds[2][0] = D2GDR2S0; d2gdrds[2][1] = D2GDR2S1; d2gdrds[2][2] = D2GDR2S2; \
 d2gdrds[3][0] = D2GDR3S0; d2gdrds[3][1] = D2GDR3S1; d2gdrds[3][2] = D2GDR3S2;

#define fillD2GDRDT \
 d2gdrdt[0] = D2GDR0DT; d2gdrdt[1] = D2GDR1DT; d2gdrdt[2] = D2GDR2DT; d2gdrdt[3] = D2GDR3DT;

#define fillD2GDRDP \
 d2gdrdp[0] = D2GDR0DP; d2gdrdp[1] = D2GDR1DP; d2gdrdp[2] = D2GDR2DP; d2gdrdp[3] = D2GDR3DP;

#define fillD2GDS2 \
 d2gds2[0][0] = D2GDS0S0;     d2gds2[0][1] = D2GDS0S1;     d2gds2[0][2] = D2GDS0S2; \
 d2gds2[1][0] = d2gds2[0][1]; d2gds2[1][1] = D2GDS1S1;     d2gds2[1][2] = D2GDS1S2; \
 d2gds2[2][0] = d2gds2[0][2]; d2gds2[2][1] = d2gds2[1][2]; d2gds2[2][2] = D2GDS2S2;

#define fillD2GDSDT \
 d2gdsdt[0] = D2GDS0DT;  d2gdsdt[1] = D2GDS1DT; d2gdsdt[2] = D2GDS2DT;

#define fillD2GDSDP \
 d2gdsdp[0] = D2GDS0DP;  d2gdsdp[1] = D2GDS1DP; d2gdsdp[2] = D2GDS2DP;

#define fillD3GDR3 \
 d3gdr3[0][0][0] = D3GDR0R0R0;	     d3gdr3[0][0][1] = D3GDR0R0R1;       d3gdr3[0][0][2] = D3GDR0R0R2;      d3gdr3[0][0][3] = D3GDR0R0R3;      \
 d3gdr3[0][1][0] = d3gdr3[0][0][1];  d3gdr3[0][1][1] = D3GDR0R1R1;       d3gdr3[0][1][2] = D3GDR0R1R2;      d3gdr3[0][1][3] = D3GDR0R1R3;      \
 d3gdr3[0][2][0] = d3gdr3[0][0][2];  d3gdr3[0][2][1] = d3gdr3[0][1][2];  d3gdr3[0][2][2] = D3GDR0R2R2;      d3gdr3[0][2][3] = D3GDR0R2R3;      \
 d3gdr3[0][3][0] = d3gdr3[0][0][3];  d3gdr3[0][3][1] = d3gdr3[0][1][3];  d3gdr3[0][3][2] = d3gdr3[0][2][3]; d3gdr3[0][3][3] = D3GDR0R3R3;      \
 d3gdr3[1][0][0] = d3gdr3[0][0][1];  d3gdr3[1][0][1] = d3gdr3[0][1][1];	 d3gdr3[1][0][2] = d3gdr3[0][1][2]; d3gdr3[1][0][3] = d3gdr3[0][1][3]; \
 d3gdr3[1][1][0] = d3gdr3[0][1][1];  d3gdr3[1][1][1] = D3GDR1R1R1;       d3gdr3[1][1][2] = D3GDR1R1R2;      d3gdr3[1][1][3] = D3GDR1R1R3;      \
 d3gdr3[1][2][0] = d3gdr3[0][1][2];  d3gdr3[1][2][1] = d3gdr3[1][1][2];  d3gdr3[1][2][2] = D3GDR1R2R2;      d3gdr3[1][2][3] = D3GDR1R2R3;      \
 d3gdr3[1][3][0] = d3gdr3[0][1][3];  d3gdr3[1][3][1] = d3gdr3[1][1][3];  d3gdr3[1][3][2] = d3gdr3[1][2][3]; d3gdr3[1][3][3] = D3GDR1R3R3;      \
 d3gdr3[2][0][0] = d3gdr3[0][0][2];  d3gdr3[2][0][1] = d3gdr3[0][1][2];  d3gdr3[2][0][2] = d3gdr3[0][2][2]; d3gdr3[2][0][3] = d3gdr3[0][2][3]; \
 d3gdr3[2][1][0] = d3gdr3[0][1][2];  d3gdr3[2][1][1] = d3gdr3[1][1][2];	 d3gdr3[2][1][2] = d3gdr3[1][2][2]; d3gdr3[2][1][3] = d3gdr3[1][2][3]; \
 d3gdr3[2][2][0] = d3gdr3[0][2][2];  d3gdr3[2][2][1] = d3gdr3[1][2][2];  d3gdr3[2][2][2] = D3GDR2R2R2;      d3gdr3[2][2][3] = D3GDR2R2R3;      \
 d3gdr3[2][3][0] = d3gdr3[0][2][3];  d3gdr3[2][3][1] = d3gdr3[1][2][3];  d3gdr3[2][3][2] = d3gdr3[2][2][3]; d3gdr3[2][3][3] = D3GDR2R3R3;      \
 d3gdr3[3][0][0] = d3gdr3[0][0][3];  d3gdr3[3][0][1] = d3gdr3[0][1][3];  d3gdr3[3][0][2] = d3gdr3[0][2][3]; d3gdr3[3][0][3] = d3gdr3[0][3][3]; \
 d3gdr3[3][1][0] = d3gdr3[0][1][3];  d3gdr3[3][1][1] = d3gdr3[1][1][3];  d3gdr3[3][1][2] = d3gdr3[1][2][3]; d3gdr3[3][1][3] = d3gdr3[1][3][3]; \
 d3gdr3[3][2][0] = d3gdr3[0][2][3];  d3gdr3[3][2][1] = d3gdr3[1][2][3];  d3gdr3[3][2][2] = d3gdr3[2][2][3]; d3gdr3[3][2][3] = d3gdr3[2][3][3]; \
 d3gdr3[3][3][0] = d3gdr3[0][3][3];  d3gdr3[3][3][1] = d3gdr3[1][3][3];  d3gdr3[3][3][2] = d3gdr3[2][3][3]; d3gdr3[3][3][3] = D3GDR3R3R3;

#define fillD3GDR2DS \
 d3gdr2ds[0][0][0] = D3GDR0R0S0;        d3gdr2ds[0][0][1] = D3GDR0R0S1;        d3gdr2ds[0][0][2] = D3GDR0R0S2;	      \
 d3gdr2ds[0][1][0] = D3GDR0R1S0;        d3gdr2ds[0][1][1] = D3GDR0R1S1;        d3gdr2ds[0][1][2] = D3GDR0R1S2;	      \
 d3gdr2ds[0][2][0] = D3GDR0R2S0;        d3gdr2ds[0][2][1] = D3GDR0R2S1;        d3gdr2ds[0][2][2] = D3GDR0R2S2;	      \
 d3gdr2ds[0][3][0] = D3GDR0R3S0;        d3gdr2ds[0][3][1] = D3GDR0R3S1;        d3gdr2ds[0][3][2] = D3GDR0R3S2;	      \
 d3gdr2ds[1][0][0] = d3gdr2ds[0][1][0]; d3gdr2ds[1][0][1] = d3gdr2ds[0][1][1]; d3gdr2ds[1][0][2] = d3gdr2ds[0][1][2]; \
 d3gdr2ds[1][1][0] = D3GDR1R1S0;        d3gdr2ds[1][1][1] = D3GDR1R1S1;        d3gdr2ds[1][1][2] = D3GDR1R1S2;	      \
 d3gdr2ds[1][2][0] = D3GDR1R2S0;        d3gdr2ds[1][2][1] = D3GDR1R2S1;        d3gdr2ds[1][2][2] = D3GDR1R2S2;	      \
 d3gdr2ds[1][3][0] = D3GDR1R3S0;        d3gdr2ds[1][3][1] = D3GDR1R3S1;        d3gdr2ds[1][3][2] = D3GDR1R3S2;	      \
 d3gdr2ds[2][0][0] = d3gdr2ds[0][2][0]; d3gdr2ds[2][0][1] = d3gdr2ds[0][2][1]; d3gdr2ds[2][0][2] = d3gdr2ds[0][2][2]; \
 d3gdr2ds[2][1][0] = d3gdr2ds[1][2][0]; d3gdr2ds[2][1][1] = d3gdr2ds[1][2][1]; d3gdr2ds[2][1][2] = d3gdr2ds[1][2][2]; \
 d3gdr2ds[2][2][0] = D3GDR2R2S0;        d3gdr2ds[2][2][1] = D3GDR2R2S1;        d3gdr2ds[2][2][2] = D3GDR2R2S2;	      \
 d3gdr2ds[2][3][0] = D3GDR2R3S0;        d3gdr2ds[2][3][1] = D3GDR2R3S1;        d3gdr2ds[2][3][2] = D3GDR2R3S2;	      \
 d3gdr2ds[3][0][0] = d3gdr2ds[0][3][0]; d3gdr2ds[3][0][1] = d3gdr2ds[0][3][1]; d3gdr2ds[3][0][2] = d3gdr2ds[0][3][2]; \
 d3gdr2ds[3][1][0] = d3gdr2ds[1][3][0]; d3gdr2ds[3][1][1] = d3gdr2ds[1][3][1]; d3gdr2ds[3][1][2] = d3gdr2ds[1][3][2]; \
 d3gdr2ds[3][2][0] = d3gdr2ds[2][3][0]; d3gdr2ds[3][2][1] = d3gdr2ds[2][3][1]; d3gdr2ds[3][2][2] = d3gdr2ds[2][3][2]; \
 d3gdr2ds[3][3][0] = D3GDR3R3S0;	d3gdr2ds[3][3][1] = D3GDR3R3S1;        d3gdr2ds[3][3][2] = D3GDR3R3S2;

#define fillD3GDR2DT \
 d3gdr2dt[0][0] = D3GDR0R0DT;     d3gdr2dt[0][1] = D3GDR0R1DT;      d3gdr2dt[0][2] = D3GDR0R2DT;     d3gdr2dt[0][3] = D3GDR0R3DT; \
 d3gdr2dt[1][0] = d3gdr2dt[0][1]; d3gdr2dt[1][1] = D3GDR1R1DT;      d3gdr2dt[1][2] = D3GDR1R2DT;     d3gdr2dt[1][3] = D3GDR1R3DT; \
 d3gdr2dt[2][0] = d3gdr2dt[0][2]; d3gdr2dt[2][1] = d3gdr2dt[1][2];  d3gdr2dt[2][2] = D3GDR2R2DT;     d3gdr2dt[2][3] = D3GDR2R3DT; \
 d3gdr2dt[3][0] = d3gdr2dt[0][3]; d3gdr2dt[3][1] = d3gdr2dt[1][3];  d3gdr2dt[3][2] = d3gdr2dt[2][3]; d3gdr2dt[3][3] = D3GDR3R3DT;

#define fillD3GDR2DP \
 d3gdr2dp[0][0] = D3GDR0R0DP;     d3gdr2dp[0][1] = D3GDR0R1DP;      d3gdr2dp[0][2] = D3GDR0R2DP;     d3gdr2dp[0][3] = D3GDR0R3DP; \
 d3gdr2dp[1][0] = d3gdr2dp[0][1]; d3gdr2dp[1][1] = D3GDR1R1DP;      d3gdr2dp[1][2] = D3GDR1R2DP;     d3gdr2dp[1][3] = D3GDR1R3DP; \
 d3gdr2dp[2][0] = d3gdr2dp[0][2]; d3gdr2dp[2][1] = d3gdr2dp[1][2];  d3gdr2dp[2][2] = D3GDR2R2DP;     d3gdr2dp[2][3] = D3GDR2R3DP; \
 d3gdr2dp[3][0] = d3gdr2dp[0][3]; d3gdr2dp[3][1] = d3gdr2dp[1][3];  d3gdr2dp[3][2] = d3gdr2dp[2][3]; d3gdr2dp[3][3] = D3GDR3R3DP;

#define fillD3GDRDS2 \
 d3gdrds2[0][0][0] = D3GDR0S0S0;        d3gdrds2[0][0][1] = D3GDR0S0S1;        d3gdrds2[0][0][2] = D3GDR0S0S2; \
 d3gdrds2[0][1][0] = d3gdrds2[0][0][1]; d3gdrds2[0][1][1] = D3GDR0S1S1;        d3gdrds2[0][1][2] = D3GDR0S1S2; \
 d3gdrds2[0][2][0] = d3gdrds2[0][0][2]; d3gdrds2[0][2][1] = d3gdrds2[0][1][2]; d3gdrds2[0][2][2] = D3GDR0S2S2; \
 d3gdrds2[1][0][0] = D3GDR1S0S0;        d3gdrds2[1][0][1] = D3GDR1S0S1;        d3gdrds2[1][0][2] = D3GDR1S0S2; \
 d3gdrds2[1][1][0] = d3gdrds2[1][0][1]; d3gdrds2[1][1][1] = D3GDR1S1S1;        d3gdrds2[1][1][2] = D3GDR1S1S2; \
 d3gdrds2[1][2][0] = d3gdrds2[1][0][2]; d3gdrds2[1][2][1] = d3gdrds2[1][1][2]; d3gdrds2[1][2][2] = D3GDR1S2S2; \
 d3gdrds2[2][0][0] = D3GDR2S0S0;        d3gdrds2[2][0][1] = D3GDR2S0S1;        d3gdrds2[2][0][2] = D3GDR2S0S2; \
 d3gdrds2[2][1][0] = d3gdrds2[2][0][1]; d3gdrds2[2][1][1] = D3GDR2S1S1;        d3gdrds2[2][1][2] = D3GDR2S1S2; \
 d3gdrds2[2][2][0] = d3gdrds2[2][0][2]; d3gdrds2[2][2][1] = d3gdrds2[2][1][2]; d3gdrds2[2][2][2] = D3GDR2S2S2; \
 d3gdrds2[3][0][0] = D3GDR3S0S0;        d3gdrds2[3][0][1] = D3GDR3S0S1;        d3gdrds2[3][0][2] = D3GDR3S0S2; \
 d3gdrds2[3][1][0] = d3gdrds2[3][0][1]; d3gdrds2[3][1][1] = D3GDR3S1S1;        d3gdrds2[3][1][2] = D3GDR3S1S2; \
 d3gdrds2[3][2][0] = d3gdrds2[3][0][2]; d3gdrds2[3][2][1] = d3gdrds2[3][1][2]; d3gdrds2[3][2][2] = D3GDR3S2S2;

#define fillD3GDRDSDT \
 d3gdrdsdt[0][0] = D3GDR0S0DT; d3gdrdsdt[0][1] = D3GDR0S1DT; d3gdrdsdt[0][2] = D3GDR0S2DT; \
 d3gdrdsdt[1][0] = D3GDR1S0DT; d3gdrdsdt[1][1] = D3GDR1S1DT; d3gdrdsdt[1][2] = D3GDR1S2DT; \
 d3gdrdsdt[2][0] = D3GDR2S0DT; d3gdrdsdt[2][1] = D3GDR2S1DT; d3gdrdsdt[2][2] = D3GDR2S2DT; \
 d3gdrdsdt[3][0] = D3GDR3S0DT; d3gdrdsdt[3][1] = D3GDR3S1DT; d3gdrdsdt[3][2] = D3GDR3S2DT;

#define fillD3GDRDSDP \
 d3gdrdsdp[0][0] = D3GDR0S0DP; d3gdrdsdp[0][1] = D3GDR0S1DP; d3gdrdsdp[0][2] = D3GDR0S2DP; \
 d3gdrdsdp[1][0] = D3GDR1S0DP; d3gdrdsdp[1][1] = D3GDR1S1DP; d3gdrdsdp[1][2] = D3GDR1S2DP; \
 d3gdrdsdp[2][0] = D3GDR2S0DP; d3gdrdsdp[2][1] = D3GDR2S1DP; d3gdrdsdp[2][2] = D3GDR2S2DP; \
 d3gdrdsdp[3][0] = D3GDR3S0DP; d3gdrdsdp[3][1] = D3GDR3S1DP; d3gdrdsdp[3][2] = D3GDR3S2DP;

#define fillD3GDS3 \
 d3gds3[0][0][0] = D3GDS0S0S0;      d3gds3[0][0][1] = D3GDS0S0S1;      d3gds3[0][0][2] = D3GDS0S0S2; \
 d3gds3[0][1][0] = d3gds3[0][0][1]; d3gds3[0][1][1] = D3GDS0S1S1;      d3gds3[0][1][2] = D3GDS0S1S2; \
 d3gds3[0][2][0] = d3gds3[0][0][2]; d3gds3[0][2][1] = d3gds3[0][1][2]; d3gds3[0][2][2] = D3GDS0S2S2; \
 d3gds3[1][0][0] = d3gds3[0][0][1]; d3gds3[1][0][1] = d3gds3[0][1][1]; d3gds3[1][0][2] = d3gds3[0][1][2]; \
 d3gds3[1][1][0] = d3gds3[0][1][1]; d3gds3[1][1][1] = D3GDS1S1S1;      d3gds3[1][1][2] = D3GDS1S1S2; \
 d3gds3[1][2][0] = d3gds3[0][1][2]; d3gds3[1][2][1] = d3gds3[1][1][2]; d3gds3[1][2][2] = D3GDS1S2S2; \
 d3gds3[2][0][0] = d3gds3[0][0][2]; d3gds3[2][0][1] = d3gds3[0][1][2]; d3gds3[2][0][2] = d3gds3[0][2][2]; \
 d3gds3[2][1][0] = d3gds3[0][1][2]; d3gds3[2][1][1] = d3gds3[1][1][2]; d3gds3[2][1][2] = d3gds3[1][2][2]; \
 d3gds3[2][2][0] = d3gds3[0][2][2]; d3gds3[2][2][1] = d3gds3[1][2][2]; d3gds3[2][2][2] = D3GDS2S2S2;

#define fillD3GDS2DT \
 d3gds2dt[0][0] = D3GDS0S0DT;     d3gds2dt[0][1] = D3GDS0S1DT;     d3gds2dt[0][2] = D3GDS0S2DT; \
 d3gds2dt[1][0] = d3gds2dt[0][1]; d3gds2dt[1][1] = D3GDS1S1DT;     d3gds2dt[1][2] = D3GDS1S2DT; \
 d3gds2dt[2][0] = d3gds2dt[0][2]; d3gds2dt[2][1] = d3gds2dt[1][2]; d3gds2dt[2][2] = D3GDS2S2DT;

#define fillD3GDS2DP \
 d3gds2dp[0][0] = D3GDS0S0DP;     d3gds2dp[0][1] = D3GDS0S1DP;     d3gds2dp[0][2] = D3GDS0S2DP; \
 d3gds2dp[1][0] = d3gds2dp[0][1]; d3gds2dp[1][1] = D3GDS1S1DP;     d3gds2dp[1][2] = D3GDS1S2DP; \
 d3gds2dp[2][0] = d3gds2dp[0][2]; d3gds2dp[2][1] = d3gds2dp[1][2]; d3gds2dp[2][2] = D3GDS2S2DP;

#define fillD3GDSDT2 \
 d3gdsdt2[0] = D3GDS0DT2; d3gdsdt2[1] = D3GDS1DT2; d3gdsdt2[2] = D3GDS2DT2;

#define fillD3GDSDTDP \
 d3gdsdtdp[0] = D3GDS0DTDP; d3gdsdtdp[1] = D3GDS1DTDP; d3gdsdtdp[2] = D3GDS2DTDP;

#define fillD3GDSDP2 \
 d3gdsdp2[0] = D3GDS0DP2; d3gdsdp2[1] = D3GDS1DP2; d3gdsdp2[2] = D3GDS2DP2;

#define fillD3GDRDT2 \
 d3gdrdt2[0] = D3GDR0DT2; d3gdrdt2[1] = D3GDR1DT2; d3gdrdt2[2] = D3GDR2DT2; d3gdrdt2[3] = D3GDR3DT2;

#define fillD3GDRDTDP \
 d3gdrdtdp[0] = D3GDR0DTDP; d3gdrdtdp[1] = D3GDR1DTDP; d3gdrdtdp[2] = D3GDR2DTDP; d3gdrdtdp[3] = D3GDR3DTDP;

#define fillD3GDRDP2 \
 d3gdrdp2[0] = D3GDR0DP2; d3gdrdp2[1] = D3GDR1DP2; d3gdrdp2[2] = D3GDR2DP2; d3gdrdp2[3] = D3GDR3DP2;

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
    //DECLARE_SITE_FRACTIONS
    double tOld         = getTOld();
    double pOld         = getPOld();
    double *rOld        = getROld();
    double *sOld        = getSOld();
    double **d2gds2     = getD2gds2();
    gsl_matrix      *ptToD2gds2  = getPtToD2gds2();
    gsl_permutation *indexD2gds2 = getIndexD2gds2();
    int i, j, iter = 0, signum;

    //GET_SITE_FRACTIONS

    xfe2ID = (r[0]/2.0 > 0.0)		   ? r[0]/2.0		     : DBL_EPSILON;
    xmg2ID = (r[1]/2.0 > 0.0)		   ? r[1]/2.0		     : DBL_EPSILON;
    xmn2ID = (r[2]/2.0 > 0.0)		   ? r[2]/2.0		     : DBL_EPSILON;
    xti4ID = ((r[0]+r[1]+r[2])/2.0 > 0.0)    ? (r[0]+r[1]+r[2])/2.0    : DBL_EPSILON;
    xal3ID = (r[3] > 0.0) 		   ? r[3]		     : DBL_EPSILON;
    xfe3ID = (1.0-r[0]-r[1]-r[2]-r[3] > 0.0) ? 1.0-r[0]-r[1]-r[2]-r[3] : DBL_EPSILON;

    /* look-up or compute the current ordering state */
    if ( (t != tOld)       || (p != pOld) || noSave ||
       (r[0] != rOld[0]) || (r[1] != rOld[1]) || (r[2] != rOld[2]) || (r[3] != rOld[3]) ) {
        double dgds[NS], sNew[NS];
        gsl_vector_view vvToDgds = gsl_vector_view_array(dgds, (size_t) NS);

        for (i=0; i<NS; i++) { sOld[i] = 2.0; sNew[i] = 0.9*r[i]; }

        while ( ((ABS(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON) ||
             (ABS(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON) ||
             (ABS(sNew[2]-sOld[2]) > 10.0*DBL_EPSILON) ) && (iter < MAX_ITER)) {
            double s[NS], deltaS[NS];
            gsl_vector_view vvToDeltaS = gsl_vector_view_array(deltaS, (size_t) NS);

            for (i=0; i<NS; i++) s[i] = sNew[i];

            xfe2a = (r[0] + s[0])/2.0;
            xmg2a = (r[1] + s[1])/2.0;
            xmn2a = (r[2] + s[2])/2.0;
            xti4a = (r[0] - s[0] + r[1] - s[1] + r[2] - s[2])/2.0;
            xal3a = r[3];
            xfe3a = 1.0 - r[0] - r[1] - r[2] - r[3];

            xfe2b = (r[0] - s[0])/2.0;
            xmg2b = (r[1] - s[1])/2.0;
            xmn2b = (r[2] - s[2])/2.0;
            xti4b = (r[0] + s[0] + r[1] + s[1] + r[2] + s[2])/2.0;
            xal3b = r[3];
            xfe3b = 1.0 - r[0] - r[1] - r[2] - r[3];

            if (xfe2a <= 0.0) xfe2a = DBL_EPSILON; /* added in V1.0-7 */
            if (xmg2a <= 0.0) xmg2a = DBL_EPSILON; /* added in V1.0-7 */
            if (xmn2a <= 0.0) xmn2a = DBL_EPSILON; /* added in V1.0-7 */
            if (xti4a <= 0.0) xti4a = DBL_EPSILON; /* added in V1.0-7 */
            if (xal3a <= 0.0) xal3a = DBL_EPSILON; /* added in V1.0-7 */
            if (xfe3a <= 0.0) xfe3a = DBL_EPSILON; /* added in V1.0-7 */

            if (xfe2b <= 0.0) xfe2b = DBL_EPSILON; /* added in V1.0-7 */
            if (xmg2b <= 0.0) xmg2b = DBL_EPSILON; /* added in V1.0-7 */
            if (xmn2b <= 0.0) xmn2b = DBL_EPSILON; /* added in V1.0-7 */
            if (xti4b <= 0.0) xti4b = DBL_EPSILON; /* added in V1.0-7 */
            if (xal3b <= 0.0) xal3b = DBL_EPSILON; /* added in V1.0-7 */
            if (xfe3b <= 0.0) xfe3b = DBL_EPSILON; /* added in V1.0-7 */

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
                for(j=0; j<NS; j++) s[i] += - d2gds2[i][j]*dgds[j];
	} */

            gsl_matrix_scale(ptToD2gds2, -1.0);
            melts_LU_decomp(ptToD2gds2, indexD2gds2, &signum);

            melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToDgds.vector, &vvToDeltaS.vector);
            for (i=0; i<NS; i++) {
                s[i] += deltaS[i];
                s[i] = MIN(s[i],  r[i] - DBL_EPSILON);
                s[i] = MAX(s[i], -r[i] + DBL_EPSILON);
            }

            for (i=0; i<NS; i++) sNew[i] = s[i];
            iter++;
        }
        tOld = t;
        pOld = p;
        for (i=0; i<NR; i++) rOld[i] = r[i];

#ifdef DEBUG
        for (i=0; i<NS; i++) {
            if (dgds[i]*r[i] > sqrt(DBL_EPSILON) && ABS(sOld[i]) > 10.0*DBL_EPSILON) {
                printf("ERROR in RHOMBOHEDRAL.C (function ORDER). Failed to converge!\n");
                printf("  t,p  = %13.6g, %13.6g\n", t, p);
                printf("  X Il = %13.6g, X Gk = %13.6g, X Py = %13.6g, X Crn = %13.6g\n", r[0], r[1], r[2], r[3]);
                printf("  s    = %13.6g, %13.6g, %13.6g\n", sOld[0], sOld[1], sOld[2]);
                printf("  dgds = %13.6g, %13.6g, %13.6g\n", dgds[0], dgds[1], dgds[2]);
                printf("  X Fe2+ A: %13.6g  X Fe2+ B: %13.6g\n", xfe2a, xfe2b);
                printf("  X Mg   A: %13.6g  X MG   B: %13.6g\n", xmg2a, xmg2b);
                printf("  X Mn   A: %13.6g  X Mn   B: %13.6g\n", xmn2a, xmn2b);
                printf("  X Ti   A: %13.6g  X Ti   B: %13.6g\n", xti4a, xti4b);
                printf("  X Al3+ A: %13.6g  X Al3+ B: %13.6g\n", xal3a, xal3b);
                printf("  X Fe3+ A: %13.6g  X Fe3+ B: %13.6g\n", xfe3a, xfe3b);
                break;
            }
        }
#endif

        setTOld(tOld);
        setPOld(pOld);
        /* arrays (rOld, sOld, d2gds2, indexD2gds2) should be preserved automatically */

        //SET_SITE_FRACTIONS
    }

    if (mask & FIRST  ) {   /* return s        */
        for (i=0; i<NS; i++) s[i] = sOld[i];
    }
    if (mask & SECOND ) {   /* compute ds/dr:  */
        double *s = sOld;
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
        double *s = sOld;
        double d2gdsdp[NS];
        gsl_vector_view vvToDp = gsl_vector_view_array(dp, (size_t) NS),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NS);

        fillD2GDSDP

        /* original: dp[i] += - d2gds2[i][j]*d2gdsdp[j]; */
        melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdp.vector, &vvToDp.vector);

    }
    if (mask & FIFTH  ) {   /* compute d2s/dr2 */
        double *s = sOld;
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
        double *s = sOld;
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
        double *s = sOld;
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
        double *s = sOld;
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
        double *s = sOld;
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
        double *s = sOld;
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
testMsg(int mask, double t, double p,
    int na,          /* Expected number of endmember components                 */
    int nr,          /* Expected number of independent compositional variables  */
    char **names,    /* array of strings of names of endmember components       */
    char **formulas, /* array of strings of formulas of endmember components    */
    double *r,       /* array of indepependent compos variables, check bounds   */
    double *m)       /* array of moles of endmember components, check bounds    */
{
    const char *phase = "rhombohedral.c";
    const char *NAMES[NA]    = { "geikielite", "hematite", "ilmenite", "pyrophanite", "corundum" };
    const char *FORMULAS[NA] = { "MgTiO3", "Fe2O3", "FeTiO3", "MnTiO3", "Al2O3" };
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
conMsg(int inpMask, int outMask, double t, double p,
    double *e,      /* comp of rhm oxides in moles of elements                  */
    double *m,      /* comp of rhm oxides in moles of endmember components      */
    double *r,      /* comp of rhm oxides in terms of the independent comp var  */
    double *x,      /* comp of rhm oxides in mole fractions of endmember comp   */
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
    (2)  SECOND           THIRD | FOURTH | FIFTH | SIXTH | EIGHTH
    (3)  THIRD            FOURTH | SEVENTH | EIGHTH

    (1) converts a vector of moles of elements into a vector of moles of
            endmember rhm oxides components.
    (2) calculates from a vector of moles of endmember components, one or
            all of: r[], x[], dr[]/dm[] d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
    (3) calculates from a vector of independent compositional variables
            mole fractions of endmember components and/or the Jacobian matrix
            dx[]/dr[] or a vector of ordering parameters (returned in *x)

    In this routine it is assumed that the elements are in the order of atomic
    numbers and that the order of rhm oxides components has been verified as:
            m[0] = geikielite  (MgTiO3) ,
            m[1] = hematite    (Fe2O3),
            m[2] = ilmenite    (FeTiO3),
            m[3] = pyrophanite (MnTiO3),
            m[4] = corundum    (Al2O3)

    ----------------------------------------------------------------------------*/

    int i, j, k;

    if (inpMask == FIRST && outMask == SECOND) {
        /* Converts a vector of moles of elements into a vector of moles of
       end-member components.                                                 */
        double sumcat, sumchg, fe2, fe3;
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
        fe3 = 6.0*sumcat/2.0 - sumchg - 2.0*e[Fe];
        fe2 = e[Fe] - fe3;

        m[0] = e[Mg];                     /* Moles of MgTiO3 */
        m[1] = fe3/2.0;                   /* Moles of Fe2O3  */
        m[2] = fe2;                       /* Moles of FeTiO3 */
        m[3] = e[Mn] + e[Co] + e[Ni];     /* Moles of MnTiO3 */
        m[4] = (e[Al] + e[Cr])/2.0;       /* Moles of Al2O3 */

        if (m[1] < 0.0) m[1] = 0.0;
        if (m[2] < 0.0) m[2] = 0.0;

    } else if (inpMask == SECOND) {
        double sum;

        if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
            printf("Illegal call to conRhm with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        for (i=0, sum=0.0; i<NA; i++) sum += m[i];

        if (outMask & THIRD) {
            /* Converts a vector of moles of end-member components (m) into a vector
         of independent compositional variables (r) required as input for the
         remaining public functions.                                          */
            r[0] = (sum != 0.0) ? m[2]/sum : 0.0;  /* Xil = X FeTiO3 */
            r[1] = (sum != 0.0) ? m[0]/sum : 0.0;  /* Xgk = X MgTiO3 */
            r[2] = (sum != 0.0) ? m[3]/sum : 0.0;  /* Xpy = X MnTiO3 */
            r[3] = (sum != 0.0) ? m[4]/sum : 0.0;  /* XCn = X Al2O3  */
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
                    dm[0][j] = (j == 2) ? (1.0-m[2]/sum)/sum : -m[2]/SQUARE(sum);
                    dm[1][j] = (j == 0) ? (1.0-m[0]/sum)/sum : -m[0]/SQUARE(sum);
                    dm[2][j] = (j == 3) ? (1.0-m[3]/sum)/sum : -m[3]/SQUARE(sum);
                    dm[3][j] = (j == 4) ? (1.0-m[4]/sum)/sum : -m[4]/SQUARE(sum);
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
                        d2m[0][j][k]  = 2.0*m[2]/CUBE(sum);
                        d2m[0][j][k] -= (j == 2) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[0][j][k] -= (k == 2) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[1][j][k]  = 2.0*m[0]/CUBE(sum);
                        d2m[1][j][k] -= (j == 0) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[1][j][k] -= (k == 0) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[2][j][k]  = 2.0*m[3]/CUBE(sum);
                        d2m[2][j][k] -= (j == 3) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[2][j][k] -= (k == 3) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[3][j][k]  = 2.0*m[4]/CUBE(sum);
                        d2m[3][j][k] -= (j == 4) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[3][j][k] -= (k == 4) ? 1.0/SQUARE(sum) : 0.0;
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
                            d3m[0][j][k][l]  = -6.0*m[2]/QUARTIC(sum);
                            d3m[0][j][k][l] += (j == 2) ? 2.0/CUBE(sum) : 0.0;
                            d3m[0][j][k][l] += (k == 2) ? 2.0/CUBE(sum) : 0.0;
                            d3m[0][j][k][l] += (l == 2) ? 2.0/CUBE(sum) : 0.0;
                            d3m[1][j][k][l]  = -6.0*m[0]/QUARTIC(sum);
                            d3m[1][j][k][l] += (j == 0) ? 2.0/CUBE(sum) : 0.0;
                            d3m[1][j][k][l] += (k == 0) ? 2.0/CUBE(sum) : 0.0;
                            d3m[1][j][k][l] += (l == 0) ? 2.0/CUBE(sum) : 0.0;
                            d3m[2][j][k][l]  = -6.0*m[3]/QUARTIC(sum);
                            d3m[2][j][k][l] += (j == 3) ? 2.0/CUBE(sum) : 0.0;
                            d3m[2][j][k][l] += (k == 3) ? 2.0/CUBE(sum) : 0.0;
                            d3m[2][j][k][l] += (l == 3) ? 2.0/CUBE(sum) : 0.0;
                            d3m[3][j][k][l]  = -6.0*m[4]/QUARTIC(sum);
                            d3m[3][j][k][l] += (j == 4) ? 2.0/CUBE(sum) : 0.0;
                            d3m[3][j][k][l] += (k == 4) ? 2.0/CUBE(sum) : 0.0;
                            d3m[3][j][k][l] += (l == 4) ? 2.0/CUBE(sum) : 0.0;
	    }
                    }
                }
            }
        }

    } else if (inpMask == THIRD) {

        if (outMask & ~(FOURTH | SEVENTH | EIGHTH))
            printf("Illegal call to conRhm with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        if (outMask & FOURTH) {
            /* Converts a vector of independent compositional variables (r) into a
         vector of mole fractions of endmember components (x).                */
            x[0] = r[1];
            x[1] = 1.0 - r[0] - r[1] - r[2] - r[3];
            x[2] = r[0];
            x[3] = r[2];
            x[4] = r[3];
        }

        if (outMask & SEVENTH) {
            /* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
            for (i=0; i<NA; i++) for (j=0; j<NR; j++) dr[i][j] = 0.0;
            dr[0][1] =  1.0;
            dr[1][0] = -1.0; dr[1][1] = -1.0; dr[1][2] = -1.0; dr[1][3] = -1.0;
            dr[2][0] = 1.0;
            dr[3][2] = 1.0;
            dr[4][3] = 1.0;
        }

        if (outMask & EIGHTH) {
            /* Computes a vector of ordering parameters and returns it in x[] */
            order(FIRST, t, p, r, x, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
        }

    } else  {
        printf("Illegal call to conRhm with inpMask = %o and outMask = %o\n",
            inpMask, outMask);
    }

}

void
dispMsg(int mask, double t, double p, double *x,
    char **formula            /* Mineral formula for interface display MASK: 1 */
    )
{
    double *r = x;
    static char masterString[] = {
/*             1111111111222222222233333333334444444444555555555566666666667
     01234567890123456789012345678901234567890123456789012345678901234567890 */
        "Mn_.__Fe''_.__Mg_.__Fe'''_.__Al_.__Ti_.__O3" };

    if (mask & FIRST) {
        char *string = strcpy((char *) malloc((size_t) (strlen(masterString)+1)*sizeof(char)), masterString);
        double totMn, totFe2, totMg, totFe3, totAl, totTi;
        char n[5];
        int i;

        totMn  = r[2];
        totFe2 = r[0];
        totMg  = r[1];
        totFe3 = 2.0*(1.0-r[0]-r[1]-r[2]-r[3]);
        totAl  = 2.0*r[3];
        totTi  = r[0]+r[1]+r[2];

        (void) snprintf(n, 5, "%4.2f", totMn);  for (i=0; i<4; i++) string[ 2+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totFe2); for (i=0; i<4; i++) string[10+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totMg);  for (i=0; i<4; i++) string[16+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totFe3); for (i=0; i<4; i++) string[25+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totAl);  for (i=0; i<4; i++) string[31+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totTi);  for (i=0; i<4; i++) string[37+i] = n[i];

        *formula = string;
    }
}

void
actMsg(int mask, double t, double p, double *x,
    double *a,  /* (pointer to a[]) activities              BINARY MASK: 0001 */
    double *mu, /* (pointer to mu[]) chemical potentials    BINARY MASK: 0010 */
    double **dx /* (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100 */
    )           /* exclusion criteria applied to results if BINARY MASK: 1000 */
{
    double *r = x;
    double s[NS], g, dgdr[NR];
    double fr[NA][NR];
    int i, j;

    for(i=0; i<NA; i++) {
     fr[i][0] = FR0(i); /* Xil */
     fr[i][1] = FR1(i); /* Xgk */
     fr[i][2] = FR2(i); /* Xpy */
     fr[i][3] = FR3(i); /* Xcn */
    }

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

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

        pureRhm(FIRST, t, p,
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

        pureRhm(SECOND, t, p,
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

        pureRhm(FIRST, t, p,
       a0,              (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        for(i=0; i<NA; i++) {
       gs[i][0] = GS0(i);        /* s   */
       gs[i][1] = GS1(i);        /* t   */
       gs[i][2] = GS2(i);        /* u   */
       dfrdr[i][0] = DFR0DR0(i); /* Xil */
       dfrdr[i][1] = DFR1DR1(i); /* Xgk */
       dfrdr[i][2] = DFR2DR2(i); /* Xpy */
       dfrdr[i][3] = DFR3DR3(i); /* Xcn */
       dgsds[i][0] = DGS0DS0(i); /* s   */
       dgsds[i][1] = DGS1DS1(i); /* t   */
       dgsds[i][2] = DGS2DS2(i); /* u   */
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
        /* implement exclusion criteria on quantities for preclb routines   */
        static const double exclusion[NA] = {
       0.05,  /* exclusion criteria on the mole fraction of geikielite  */
       0.05,  /* exclusion criteria on the mole fraction of hematite    */
       0.05,  /* exclusion criteria on the mole fraction of ilmenite    */
       0.05,  /* exclusion criteria on the mole fraction of pyrophanite */
       0.05   /* exclusion criteria on the mole fraction of corundum    */
        };
        double x[NA];

        x[0] = r[1];                            /* mole fraction of geikielite  */
        x[1] = 1.0 - r[0] - r[1] - r[2] - r[3]; /* mole fraction of hematite    */
        x[2] = r[0];                            /* mole fraction of ilmenite    */
        x[3] = r[2];                            /* mole fraction of pyrophanite */
        x[4] = r[3];                            /* mole fraction of pyrophanite */

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
gmixMsg(int mask, double t, double p, double *x,
    double *gmix, /* Gibbs energy of mixing             BINARY MASK: 0001 */
    double *dx,   /* (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010 */
    double **dx2, /* (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100 */
    double ***dx3 /* (pointer to dx3[][][]) d3(g)/d(x[])3 NARY MASK: 1000 */
    )
{
    double *r = x;
    double s[NS];

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    if (mask & FIRST) {
        double ends[NA];

        *gmix = G;

        pureRhm(THIRD, t, p,
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

        pureRhm(THIRD, t, p,
       (double *) NULL, (double *) NULL, ends,           (double *) NULL,
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
hmixMsg(int mask, double t, double p, double *x,
    double *hmix /* Enthalpy of mixing BINARY MASK: 1 */
    )
{
    double *r = x;
    double s[NS], ends[NA];

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    *hmix = (G) - t*(DGDT);

    pureRhm(FOURTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, ends,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

    *hmix -= ENDMEMBERS;
}

void
smixMsg(int mask, double t, double p, double *x,
    double *smix, /* Entropy of mixing                  BINARY MASK: 001 */
    double *dx,   /* (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010 */
    double **dx2  /* (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100 */
    )
{
    double *r = x;
    double s[NS];

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    if (mask & FIRST) {
        double ends[NA];

        *smix = -(DGDT);

        pureRhm(FIFTH, t, p,
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

        pureRhm(FIFTH, t, p,
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
cpmixMsg(int mask, double t, double p, double *x,
    double *cpmix, /* Heat capacity of mixing               BINARY MASK: 001 */
    double *dt,    /* d(cp)/d(t)                            BINARY MASK: 010 */
    double *dx     /* d(cp)/d(x[])d(t)                      BINARY MASK: 100 */
    )
{
    double *r = x;
    double s[NS], dsdt[NS], d2gdsdt[NS], d2gds2[NS][NS], d2gdt2;
    int i, j;

    order(FIRST | THIRD, t, p, r,
                s, NULL, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

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

        pureRhm(SIXTH, t, p,
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

        pureRhm(SEVENTH, t, p,
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

        pureRhm(SIXTH, t, p,
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
vmixMsg(int mask, double t, double p, double *x,
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
    double s[NS];

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    if (mask & FIRST) {
        double ends[NA];

        *vmix = DGDP;

        pureRhm(EIGHTH, t, p,
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

        pureRhm(EIGHTH, t, p,
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

        pureRhm(NINTH, t, p,
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

        pureRhm(TENTH, t, p,
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

        pureRhm(ELEVENTH, t, p,
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

        pureRhm(TWELFTH, t, p,
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

        pureRhm(THIRTEENTH, t, p,
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

        pureRhm(NINTH, t, p,
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


        pureRhm(TENTH, t, p,
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

/* end of file RHOMBOHEDRAL.C */
