const char *gibbs_ver(void) { return "$Id: gibbs.c,v 1.11 2008/10/15 00:40:06 ghiorso Exp $"; }
/*
 MELTS Source Code: RCS $Log: gibbs.c,v $
 MELTS Source Code: RCS Revision 1.10  2008/05/23 17:07:31  ghiorso
 MELTS Source Code: RCS Added more hydroxyl species.
 MELTS Source Code: RCS Improved preclb interface for larger number of species.
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.9  2008/05/02 19:03:52  ghiorso
 MELTS Source Code: RCS Revised liquid speciation model.
 MELTS Source Code: RCS Created new test routine for homogeneous equilibrium fO2 at P.
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.8  2007/12/22 22:43:30  ghiorso
 MELTS Source Code: RCS Fixed error in BM integration in Gibbs.c
 MELTS Source Code: RCS Updated param_struct_data.h file for AGU 2007 xMELTS parameters
 MELTS Source Code: RCS Added support for status file production in MELTS-batch (XML)
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.7  2007/11/29 05:32:12  ghiorso
 MELTS Source Code: RCS Majorite testMaj corrections.  Fixed "O2" in sol_struct_data.h to "o2".
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.6  2007/10/03 21:33:47  ghiorso
 MELTS Source Code: RCS Updated liquid eos thermodynamics.
 MELTS Source Code: RCS Added regression of ferric/ferrous parameters from data file.
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.5  2007/09/13 16:12:02  ghiorso
 MELTS Source Code: RCS (1) Revised standard state liquid properties.
 MELTS Source Code: RCS (2) Revised standard state solid properties (removed non-Berman) Cp, and
 MELTS Source Code: RCS     removed Saxena EOS treatment.  All EOS parameterizations are Vinet.
 MELTS Source Code: RCS     Updated K, K', alpha to conform to Knittle (1995) and Fei (1995)
 MELTS Source Code: RCS     except where refitted Berman (1988) makes more sense.
 MELTS Source Code: RCS (3) Updated code to allow for fusion entropies of liquid components to
 MELTS Source Code: RCS     be adjusted (fusion enthalpies are dependent).
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.4  2006/10/20 00:59:22  ghiorso
 MELTS Source Code: RCS (1) Made initial modifications for thread safe code.
 MELTS Source Code: RCS (2) Added support for XML I/O in batch mode
 MELTS Source Code: RCS (3) Added support for Melts-batch listener for eventual integration into VIGMCS
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.3  2006/08/17 20:47:54  ghiorso
 MELTS Source Code: RCS Clarified variable initialization issues in routines.  Problems discovered
 MELTS Source Code: RCS when compiler optimization is turned on.
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
 MELTS Source Code: RCS Revision 1.4  2005/06/10 19:00:16  cvsaccount
 MELTS Source Code: RCS *** empty log message ***
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.3  2004/12/04 19:10:36  cvsaccount
 MELTS Source Code: RCS *** empty log message ***
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.2  2004/10/02 19:41:47  cvsaccount
 MELTS Source Code: RCS *** empty log message ***
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.1.1.1  2004/01/02 19:21:49  cvsaccount
 MELTS Source Code: RCS CTserver University of Chicago
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.7  2003/09/27 15:35:22  ghiorso
 MELTS Source Code: RCS *** empty log message ***
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.6  2002/07/22 03:02:40  ghiorso
 MELTS Source Code: RCS *** empty log message ***
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.5  2002/06/28 03:33:12  ghiorso
 MELTS Source Code: RCS *** empty log message ***
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.4  2002/01/21 19:48:03  ghiorso
 MELTS Source Code: RCS *** empty log message ***
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.3  2002/01/13 05:48:56  ghiorso
 MELTS Source Code: RCS *** empty log message ***
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.2  2002/01/10 02:28:04  ghiorso
 MELTS Source Code: RCS *** empty log message ***
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
 MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
 MELTS Source Code: RCS
 MELTS Source Code: RCS Revision 5.1  2000/02/15 17:58:12  ghiorso
 MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
 MELTS Source Code: RCS
 * Revision 3.7  1997/06/21  22:49:55  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.6  1997/05/03  20:23:34  ghiorso
 * *** empty log message ***
 *
 * Revision 3.5  1997/03/27  17:03:37  ghiorso
 * *** empty log message ***
 *
 * Revision 3.4  1996/09/24  20:33:40  ghiorso
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
 * Revision 3.1  1995/08/18  17:32:47  ghiorso
 * MELTS Version 3 - Initial Entry
 *
 */

/*
 **++
 **  FACILITY:  Silicate Melts Regression/Crystallization Package
 **
 **  MODULE DESCRIPTION:
 **
 **      Collection of functions to implement Berman (1988) thermodynamic
 **         database (file: GIBBS.C)
 **
 **  MODIFICATION HISTORY:
 **
 **      V1.0-1  Mark S. Ghiorso  April 28, 1991 Original Version
 **      V1.1-1  Mark S. Ghiorso  September 7, 1991
 **              Added global function getlog10fo2
 **      V1.1-2  Mark S. Ghiorso  October 10, 1991
 **              Corrected error in convergence criteria for albite
 **      V1.1-3  Mark S. Ghiorso  January 22, 1992
 **              (1) Altered liquid water properties to conform to Nicholls
 **              (2) Corrected error in liquid properties for extrapolation
 **                  below the glass transition
 **      V1.1-4  Mark S. Ghiorso  March 23, 1992
 **              Corrected calls to pow() functions for (double) arguments and
 **              defined integer power macros
 **      V1.1-5  Mark S. Ghiorso  January 12, 1992
 **              (1) Removed Richet algorithm for liquid SiO2
 **      V1.2-1  Mark S. Ghiorso  May 31, 1993
 **              Corrected error in array initialization in wdh78 (caused
 **              bad G for water above 10 kb)
 **      V1.2-2  Mark S. Ghiorso  June 14, 1993
 **              Added new fo2 path options to getlog10fo2().
 **      V2.0-1  Mark S. Ghiorso  May 12, 1994
 **              (1) Updated water calculations to obtain needed endmember
 **                  properties
 **              (2) Updated liquid H2O calculations to obtain needed
 **                  endmember properties
 **              (3) Updated albite calculations to obtain endmember
 **                  properties
 **              (4) Changed analytical calculation of H, S and Cp to output
 **                  properties at t and p
 **      V2.0-2  Mark S. Ghiorso  June 17, 1994
 **              (1) Added new fo2 path functions:
 **                  getdlog10fo2dt(), getdlog10fo2dp(), getd2log10fo2dt2() and
 **                  getd2log10fo2dp2().
 **      V2.0-3  Mark S. Ghiorso  July 11, 1994
 **              (1) Added macros flag SPECIAL_O2 to change standard state of O2
 **                  to that of the melt (new function getO2properties)
 **      V2.1-1  Mark S. Ghiorso  February 4, 1995
 **              (1) Corrected error in cps, dcpsdt and d2vdt2 equation for
 **                  gehlenite and sanidine (coeff of 0.5 in d1 term)
 **              (2) Corrected error in computation of s for albite (added
 **                  pressure correction from albite() call
 **              (3) Created new external albite() routine
 **              (4) Corrected error in generic solids dcpsdt computation
 **                  (neglected pressure correction term)
 **                               February 15, 1995
 **              (5) Replaced water routines
 **      V2.1-2  Mark S. Ghiorso  February 25, 1995
 **              (1) Altered return parameters for wdh78
 **              (2) Corrected errors in volume and derivatives in quartz and
 **                  cristobalite (below the lambda transition)
 **      V2.1-3  Mark S. Ghiorso  February 28, 1995
 **              (1) Transferred Berman's reference state correction for
 **                  water and H2O to the water.c file
 **      V2.2-1  Mark S. Ghiorso  March 24, 1995
 **              (1) Experiments with O2 properties
 **--
 */

#include "silmin.h"  /* SILMIN structures include file */

/*
 #define SPECIAL_O2
 #define ZERO_O2
 */
#define CONSTANT_P_O2

#if defined(RHYOLITE_ADJUSTMENTS) && ! defined(TRUE_xMELTS)
#define QUARTZ_ADJUSTMENT ((calculationMode != MODE_pMELTS) ? -1291.0 : 0.0)
#define TRIDYMITE_ADJUSTMENT ((calculationMode != MODE_pMELTS) ? -2625.0 : 0.0)
#define CRISTOBALITE_ADJUSTMENT ((calculationMode != MODE_pMELTS) ? -450.0 : 0.0)
#else
#define QUARTZ_ADJUSTMENT 0.0
#define TRIDYMITE_ADJUSTMENT 0.0
#define CRISTOBALITE_ADJUSTMENT 0.0
#endif

#define NO_HIGH_STRUCTURAL_STATE_FELDSPAR

#define WETTING_ANGLE_CORR 0.0

#define SQUARE(x)  ((x)*(x))
#define CUBE(x)    ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))
#define QUINTIC(x) ((x)*(x)*(x)*(x)*(x))

extern void propertiesOfPureH2O(double t, double p,
                                double *g, double *h, double *s, double *cp, double *dcpdt, double *v, double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2);
extern void propertiesOfPureCO2(double t, double p,
                                double *g, double *h, double *s, double *cp, double *dcpdt, double *v, double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2);

static double getlog10fo2COH(int mask, double t, double p) {
    double pLim = (p < 67201.0) ? p : 67201.0;
    double result = 0.0;

    if     (mask == FIRST) { /* log fO2 */
        result = (-22324.0 + 189.0*pLim/1000.0 - 1.41*pLim*pLim/1000000.0)/t +  4.62;
        if (p > 67201.0) result += 608.966*(p/10000.0)/t - 608.966*(67201.0/10000.0)/t;
    } else if (mask == SECOND) { /* /dT */
        result = -(-22324.0 + 189.0*pLim/1000.0 - 1.41*pLim*pLim/1000000.0)/SQUARE(t);
        if (p > 67201.0) result += -608.966*(p/10000.0)/SQUARE(t) + 608.966*(67201.0/10000.0)/SQUARE(t);
    } else if (mask == THIRD) { /* /dP */
        result = (189.0/1000.0 - 2.0*1.41*pLim/1000000.0)/t;
        if (p > 67201.0) result += 608.966*(1.0/10000.0)/t;
    } else if (mask == FOURTH) { /* /dT2 */
        result = 2.0*(-22324.0 + 189.0*pLim/1000.0 - 1.41*pLim*pLim/1000000.0)/CUBE(t);
        if (p > 67201.0) result += 2.0*608.966*(p/10000.0)/CUBE(t) - 2.0*608.966*(67201.0/10000.0)/CUBE(t);
    } else if (mask == FIFTH) { /* /dP2 */
        result = (-2.0*1.41/1000000.0)/t;
    }

    return result;
}

/*
 *=============================================================================
 * Public function (These calculations should be reviewed for internal
 *                  consistency with Berman 1988 and for the pressure
 *                  effect)
 * The silminState->fo2Delta has been added to implement non-integer offsets
 * to implement non-integer offsets to fO2 buffer in Java MELTS
 */
double getlog10fo2(double t, double p, int buffer)
{
    double fo2Delta = (silminState != NULL) ? silminState->fo2Delta : 0.0;
    if      (buffer == FO2_HM )    return -23847.6/t + 13.480 + fo2Delta;
    else if (buffer == FO2_NNO)    return -24930.0/t +  9.360 + fo2Delta;
    else if (buffer == FO2_QFM)    return -24441.9/t + 0.110*(p-1.0)/t +  8.290 + fo2Delta;
    else if (buffer == FO2_COH)    return getlog10fo2COH(FIRST, t, p)  + fo2Delta;
    else if (buffer == FO2_IW )    return -26834.7/t +  6.471 + fo2Delta;
    else if (buffer == FO2_QFM_P3) return -24441.9/t + 0.110*(p-1.0)/t +  8.290 + 3.0;
    else if (buffer == FO2_QFM_P2) return -24441.9/t + 0.110*(p-1.0)/t +  8.290 + 2.0;
    else if (buffer == FO2_QFM_P1) return -24441.9/t + 0.110*(p-1.0)/t +  8.290 + 1.0;
    else if (buffer == FO2_QFM_M1) return -24441.9/t + 0.110*(p-1.0)/t +  8.290 - 1.0;
    else if (buffer == FO2_QFM_M2) return -24441.9/t + 0.110*(p-1.0)/t +  8.290 - 2.0;
    else if (buffer == FO2_QFM_M3) return -24441.9/t + 0.110*(p-1.0)/t +  8.290 - 3.0;
    else if (buffer == FO2_QFM_M4) return -24441.9/t + 0.110*(p-1.0)/t +  8.290 - 4.0;
    else if (buffer == FO2_QFM_M5) return -24441.9/t + 0.110*(p-1.0)/t +  8.290 - 5.0;
    else if (buffer == FO2_QFM_M6) return -24441.9/t + 0.110*(p-1.0)/t +  8.290 - 6.0;
    else if (buffer == FO2_QFM_M7) return -24441.9/t + 0.110*(p-1.0)/t +  8.290 - 7.0;
    else if (buffer == FO2_QFM_M8) return -24441.9/t + 0.110*(p-1.0)/t +  8.290 - 8.0;
    else if (buffer == FO2_QFM_M9) return -24441.9/t + 0.110*(p-1.0)/t +  8.290 - 9.0;
    else if (buffer == FO2_QFM_P0_5) return -24441.9/t + 0.110*(p-1.0)/t +  8.290 + 0.5;
    else if (buffer == FO2_QFM_P1_5) return -24441.9/t + 0.110*(p-1.0)/t +  8.290 + 1.5;
#ifdef PHMELTS_ADJUSTMENTS
    else if (buffer == FO2_ABS)    return fo2Delta;
#endif
    else                           return 0.0;
}

double getdlog10fo2dt(double t, double p, int buffer)
{
    if      (buffer == FO2_HM )    return 23847.6/SQUARE(t);
    else if (buffer == FO2_NNO)    return 24930.0/SQUARE(t);
    else if (buffer == FO2_QFM)    return 24441.9/SQUARE(t) - 0.110*(p-1.0)/SQUARE(t);
    else if (buffer == FO2_COH)    return getlog10fo2COH(SECOND, t, p);
    else if (buffer == FO2_IW )    return 26834.7/SQUARE(t);
    else if (buffer == FO2_QFM_P3) return 24441.9/SQUARE(t) - 0.110*(p-1.0)/SQUARE(t);
    else if (buffer == FO2_QFM_P2) return 24441.9/SQUARE(t) - 0.110*(p-1.0)/SQUARE(t);
    else if (buffer == FO2_QFM_P1) return 24441.9/SQUARE(t) - 0.110*(p-1.0)/SQUARE(t);
    else if (buffer == FO2_QFM_M1) return 24441.9/SQUARE(t) - 0.110*(p-1.0)/SQUARE(t);
    else if (buffer == FO2_QFM_M2) return 24441.9/SQUARE(t) - 0.110*(p-1.0)/SQUARE(t);
    else if (buffer == FO2_QFM_M3) return 24441.9/SQUARE(t) - 0.110*(p-1.0)/SQUARE(t);
    else if (buffer == FO2_QFM_M4) return 24441.9/SQUARE(t) - 0.110*(p-1.0)/SQUARE(t);
    else if (buffer == FO2_QFM_M5) return 24441.9/SQUARE(t) - 0.110*(p-1.0)/SQUARE(t);
    else if (buffer == FO2_QFM_M6) return 24441.9/SQUARE(t) - 0.110*(p-1.0)/SQUARE(t);
    else if (buffer == FO2_QFM_M7) return 24441.9/SQUARE(t) - 0.110*(p-1.0)/SQUARE(t);
    else if (buffer == FO2_QFM_M8) return 24441.9/SQUARE(t) - 0.110*(p-1.0)/SQUARE(t);
    else if (buffer == FO2_QFM_M9) return 24441.9/SQUARE(t) - 0.110*(p-1.0)/SQUARE(t);
    else if (buffer == FO2_QFM_P0_5) return 24441.9/SQUARE(t) - 0.110*(p-1.0)/SQUARE(t);
    else if (buffer == FO2_QFM_P1_5) return 24441.9/SQUARE(t) - 0.110*(p-1.0)/SQUARE(t);
    else                           return 0.0;
}

double getdlog10fo2dp(double t, double p, int buffer)
{
    if      (buffer == FO2_HM )    return 0.0;
    else if (buffer == FO2_NNO)    return 0.0;
    else if (buffer == FO2_QFM)    return 0.110/t;
    else if (buffer == FO2_COH)    return getlog10fo2COH(THIRD, t, p);
    else if (buffer == FO2_IW )    return 0.0;
    else if (buffer == FO2_QFM_P3) return 0.110/t;
    else if (buffer == FO2_QFM_P2) return 0.110/t;
    else if (buffer == FO2_QFM_P1) return 0.110/t;
    else if (buffer == FO2_QFM_M1) return 0.110/t;
    else if (buffer == FO2_QFM_M2) return 0.110/t;
    else if (buffer == FO2_QFM_M3) return 0.110/t;
    else if (buffer == FO2_QFM_M4) return 0.110/t;
    else if (buffer == FO2_QFM_M5) return 0.110/t;
    else if (buffer == FO2_QFM_M6) return 0.110/t;
    else if (buffer == FO2_QFM_M7) return 0.110/t;
    else if (buffer == FO2_QFM_M8) return 0.110/t;
    else if (buffer == FO2_QFM_M9) return 0.110/t;
    else if (buffer == FO2_QFM_P0_5) return 0.110/t;
    else if (buffer == FO2_QFM_P1_5) return 0.110/t;
    else                           return 0.110/t;
}

double getd2log10fo2dt2(double t, double p, int buffer)
{
    if      (buffer == FO2_HM )    return -2.0*23847.6/CUBE(t);
    else if (buffer == FO2_NNO)    return -2.0*24930.0/CUBE(t);
    else if (buffer == FO2_QFM)    return -2.0*24441.9/CUBE(t) + 2.0*0.110*(p-1.0)/CUBE(t);
    else if (buffer == FO2_COH)    return getlog10fo2COH(FOURTH, t, p);
    else if (buffer == FO2_IW )    return -2.0*26834.7/CUBE(t);
    else if (buffer == FO2_QFM_P3) return -2.0*24441.9/CUBE(t) + 2.0*0.110*(p-1.0)/CUBE(t);
    else if (buffer == FO2_QFM_P2) return -2.0*24441.9/CUBE(t) + 2.0*0.110*(p-1.0)/CUBE(t);
    else if (buffer == FO2_QFM_P1) return -2.0*24441.9/CUBE(t) + 2.0*0.110*(p-1.0)/CUBE(t);
    else if (buffer == FO2_QFM_M1) return -2.0*24441.9/CUBE(t) + 2.0*0.110*(p-1.0)/CUBE(t);
    else if (buffer == FO2_QFM_M2) return -2.0*24441.9/CUBE(t) + 2.0*0.110*(p-1.0)/CUBE(t);
    else if (buffer == FO2_QFM_M3) return -2.0*24441.9/CUBE(t) + 2.0*0.110*(p-1.0)/CUBE(t);
    else if (buffer == FO2_QFM_M4) return -2.0*24441.9/CUBE(t) + 2.0*0.110*(p-1.0)/CUBE(t);
    else if (buffer == FO2_QFM_M5) return -2.0*24441.9/CUBE(t) + 2.0*0.110*(p-1.0)/CUBE(t);
    else if (buffer == FO2_QFM_M6) return -2.0*24441.9/CUBE(t) + 2.0*0.110*(p-1.0)/CUBE(t);
    else if (buffer == FO2_QFM_M7) return -2.0*24441.9/CUBE(t) + 2.0*0.110*(p-1.0)/CUBE(t);
    else if (buffer == FO2_QFM_M8) return -2.0*24441.9/CUBE(t) + 2.0*0.110*(p-1.0)/CUBE(t);
    else if (buffer == FO2_QFM_M9) return -2.0*24441.9/CUBE(t) + 2.0*0.110*(p-1.0)/CUBE(t);
    else if (buffer == FO2_QFM_P0_5) return -2.0*24441.9/CUBE(t) + 2.0*0.110*(p-1.0)/CUBE(t);
    else if (buffer == FO2_QFM_P1_5) return -2.0*24441.9/CUBE(t) + 2.0*0.110*(p-1.0)/CUBE(t);
    else                           return 0.0;
}

double getd2log10fo2dp2(double t, double p, int buffer)
{
    if      (buffer == FO2_HM )    return 0.0;
    else if (buffer == FO2_NNO)    return 0.0;
    else if (buffer == FO2_QFM)    return 0.0;
    else if (buffer == FO2_COH)    return getlog10fo2COH(FIFTH, t, p);
    else if (buffer == FO2_IW )    return 0.0;
    else if (buffer == FO2_QFM_P3) return 0.0;
    else if (buffer == FO2_QFM_P2) return 0.0;
    else if (buffer == FO2_QFM_P1) return 0.0;
    else if (buffer == FO2_QFM_M1) return 0.0;
    else if (buffer == FO2_QFM_M2) return 0.0;
    else if (buffer == FO2_QFM_M3) return 0.0;
    else if (buffer == FO2_QFM_M4) return 0.0;
    else if (buffer == FO2_QFM_M5) return 0.0;
    else if (buffer == FO2_QFM_M6) return 0.0;
    else if (buffer == FO2_QFM_M7) return 0.0;
    else if (buffer == FO2_QFM_M8) return 0.0;
    else if (buffer == FO2_QFM_M9) return 0.0;
    else if (buffer == FO2_QFM_P0_5) return 0.0;
    else if (buffer == FO2_QFM_P1_5) return 0.0;
    else                           return 0.0;
}

/*
 *=============================================================================
 * Private functions:
 */

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
void whaar(double p, double t, double *gH2O, double *hH2O, double *sH2O,
           double *cpH2O, double *dcpdtH2O, double *vH2O, double *dvdtH2O,
           double *dvdpH2O, double *d2vdt2H2O, double *d2vdtdpH2O, double *d2vdp2H2O);
void wdh78(double p, double tk, double *gDELTA, double *hDELTA, double *sDELTA,
           double *cpDELTA, double *dcpdtDELTA, double *vTOTAL, double *dvdtTOTAL,
           double *dvdpTOTAL, double *d2vdt2TOTAL, double *d2vdtdpTOTAL,
           double *d2vdp2TOTAL);
void albite(double p, double t, double *gDis, double *hDis, double *sDis,
            double *cpDis, double *dcpdtDis, double *vDis, double *dvdtDis,
            double *dvdpDis, double *d2vdt2Dis, double *d2vdtdpDis, double *d2vdp2Dis);
#ifdef SPECIAL_O2
static void getO2properties(double t, double p, double *g, double *h,
                            double *s, double *cp, double *dcpdt, double *v, double *dvdt, double *dvdp,
                            double *d2vdt2, double *d2vdp2, double *d2vdtdp);
#endif
static void fe_metal(double t, double p, double *gs, double *hs,
                     double *ss, double *cps, double *dcpdts, double *vs, double *dvsdt,
                     double *dvsdp, double *d2vsdt2, double *d2vsdp2, double *d2vsdtdp);
static void fe_liquid(double t, double p, double *gl, double *hl,
                      double *sl, double *cpl, double *dcpldt, double *vl, double *dvldt,
                      double *dvldp, double *d2vldt2, double *d2vldp2, double *d2vldtdp);
static void ni_metal(double t, double p, double *gs, double *hs,
                     double *ss, double *cps, double *dcpdts, double *vs, double *dvsdt,
                     double *dvsdp, double *d2vsdt2, double *d2vsdp2, double *d2vsdtdp);
static void ni_liquid(double t, double p, double *gl, double *hl,
                      double *sl, double *cpl, double *dcpldt, double *vl, double *dvldt,
                      double *dvldp, double *d2vldt2, double *d2vldp2, double *d2vldtdp);
static void getSiOHproperties(double t, double *g, double *h, double *s,
                              double *cp, double *dcpdt);
static void getAlOHproperties(double t, double *g, double *h, double *s,
                              double *cp, double *dcpdt);
static void getFe3OHproperties(double t, double *g, double *h, double *s,
                               double *cp, double *dcpdt);
static void getFe2OHproperties(double t, double *g, double *h, double *s,
                               double *cp, double *dcpdt);
static void getMgOHproperties(double t, double *g, double *h, double *s,
                              double *cp, double *dcpdt);
static void getCaOHproperties(double t, double *g, double *h, double *s,
                              double *cp, double *dcpdt);
static void getNaOHproperties(double t, double *g, double *h, double *s,
                              double *cp, double *dcpdt);
static void getKOHproperties(double t, double *g, double *h, double *s,
                             double *cp, double *dcpdt);
static void getFe2AlO4_5properties(double t, double *g, double *h, double *s,
                                   double *cp, double *dcpdt);
static void getFe2AlO4_1properties(double t, double *g, double *h, double *s,
                                   double *cp, double *dcpdt);
static void getCaCO3properties(double t, double p, double *g, double *h, double *s,
                               double *cp, double *dcpdt, double *v, double *dvdt, double *dvdp,
                               double *d2vdt2, double *d2vdtdp, double *d2vdp2);

static void intEOSsolid(ThermoRef *phase, double t, double p, double *g,
                        double *h, double *s, double *cp, double *dcpdt, double *v, double *dvdt,
                        double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2)
{
    const double tr = 298.15; /* K    */
    const double pr = 1.0;    /* bars */

    if (phase->eos_type == EOS_BERMAN) {
        double v0 = phase->v;
        double v1 = phase->eos.Berman.v1;
        double v2 = phase->eos.Berman.v2;
        double v3 = phase->eos.Berman.v3;
        double v4 = phase->eos.Berman.v4;

        *g       += v0*((v1/2.0-v2)*(p*p-pr*pr)+v2*(p*p*p-pr*pr*pr)/3.0 + (1.0-v1+v2+v3*(t-tr)+v4*SQUARE(t-tr))*(p-pr));
        *h       += v0*((v1/2.0-v2)*(p*p-pr*pr)+v2*(p*p*p-pr*pr*pr)/3.0 + (1.0-v1+v2+v3*(t-tr)+v4*SQUARE(t-tr))*(p-pr)) - t*v0*(v3 + 2.0*(t-tr)*v4)*(p-pr);
        *s       += -v0*(v3 + 2.0*(t-tr)*v4)*(p-pr);
        *cp      += -t*v0*2.0*v4*(p-pr);
        *dcpdt   += -v0*2.0*v4*(p-pr);
        *v       += v0*(1.0 + (p-pr)*v1 + (p-pr)*(p-pr)*v2 + (t-tr)*v3 + (t-tr)*(t-tr)*v4);
        *dvdt    += v0*(v3 + 2.0*(t-tr)*v4);
        *dvdp    += v0*(v1 + 2.0*(p-pr)*v2);
        *d2vdt2  += v0*2.0*v4;
        *d2vdtdp += 0.0;
        *d2vdp2  += v0*2.0*v2;
    } else if(phase->eos_type == EOS_VINET) {
        double v0    = phase->v;
        double alpha = phase->eos.Vinet.alpha;
        double K     = phase->eos.Vinet.K;
        double Kp    = phase->eos.Vinet.Kp;
        double eta   = 3.0*(Kp-1.0)/2.0;
        double x     = 1.0;
        double x0    = 1.0;
        double fn, dfn, a, dxdt, dxdp, d2xdt2, d2xdtdp, d2xdp2, dx0dt, d2x0dt2;
        int iter;

        iter = 0;
        do {
            fn = x*x*(p/10000.0) - 3.0*K*(1.0-x)*exp(eta*(1.0-x)) - x*x*alpha*K*(t-tr);
            dfn = 2.0*x*(p/10000.0) + 3.0*K*(1.0+eta*(1.0-x))*exp(eta*(1.0-x)) - 2.0*x*alpha*K*(t-tr);
            x = x - fn/dfn;
            iter++;
        } while ((iter < 500) && (fn*fn > DBL_EPSILON));
        dxdt    = -(1.0/3.0)*x*x*x*alpha*K/(K*exp(eta*(1.0-x))*(-2.0+x-eta*x+eta*x*x));
        d2xdt2  = -dxdt*(
                         2.0*(p/10000.0)*dxdt-6.0*K*eta*dxdt*exp(eta*(1.0-x))-3.0*K*eta*eta*dxdt*exp(eta*(1.0-x))+3.0*K*eta*eta*x*dxdt*exp(eta*(1.0-x))
                         -4.0*alpha*K*x-2.0*alpha*K*dxdt*(t-tr))/(2.0*(p/10000.0)*x+3.0*K*exp(eta*(1.0-x))+3.0*K*eta*exp(eta*(1.0-x))
                                                                  -3.0*K*eta*x*exp(eta*(1.0-x))-2.0*alpha*K*x*(t-tr));
        dxdp    = -(1.0/3.0)*x*x*x/(K*exp(eta*(1.0-x))*(2.0-x+eta*x-eta*x*x));
        d2xdtdp = -(2.0*x*dxdt+2.0*(p/10000.0)*dxdp*dxdt-2.0*alpha*K*dxdp*dxdt*(t-tr)-3.0*K*eta*eta*dxdt*dxdp*exp(eta*(1.0-x))
                    -6.0*K*dxdt*eta*dxdp*exp(eta*(1.0-x))-2.0*alpha*K*x*dxdp+3.0*K*eta*eta*dxdt*dxdp*exp(eta*(1.0-x))*x)/
        (2.0*(p/10000.0)*x+3.0*K*exp(eta*(1.0-x))+3.0*K*eta*exp(eta*(1.0-x))-3.0*K*eta*exp(eta*(1.0-x))*x-2.0*alpha*K*x*(t-tr));
        d2xdp2  = - dxdp*dxdp*((-6.0+2.0*x-4.0*eta*x+2.0*eta*x*x-eta*eta*x*x+eta*eta*x*x*x)*exp(eta*(1.0-x)))/(x*exp(eta*(1.0-x))*(2.0-x+eta*x-eta*x*x));

        iter = 0;
        do {
            fn = x0*x0*(pr/10000.0) - 3.0*K*(1.0-x0)*exp(eta*(1.0-x0)) - x0*x0*alpha*K*(t-tr);
            dfn = 2.0*x0*(pr/10000.0) + 3.0*K*(1.0+eta*(1.0-x0))*exp(eta*(1.0-x0)) - 2.0*x0*alpha*K*(t-tr);
            x0 = x0 - fn/dfn;
            iter++;
        } while ((iter < 500) && (fn*fn > DBL_EPSILON));
        dx0dt    = -(1.0/3.0)*x0*x0*x0*alpha*K/(K*exp(eta*(1.0-x0))*(-2.0+x0-eta*x0+eta*x0*x0));
        d2x0dt2  = -dx0dt*(
                           2.0*(pr/10000.0)*dx0dt-6.0*K*eta*dx0dt*exp(eta*(1.0-x0))-3.0*K*eta*eta*dx0dt*exp(eta*(1.0-x0))+3.0*K*eta*eta*x0*dx0dt*exp(eta*(1.0-x0))
                           -4.0*alpha*K*x0-2.0*alpha*K*dx0dt*(t-tr))/(2.0*(pr/10000.0)*x0+3.0*K*exp(eta*(1.0-x0))+3.0*K*eta*exp(eta*(1.0-x0))
                                                                      -3.0*K*eta*x0*exp(eta*(1.0-x0))-2.0*alpha*K*x0*(t-tr));

        a  = (9.0*v0*K/(eta*eta))*(1.0 - eta*(1.0-x))*exp(eta*(1.0-x));
        a += v0*(t-tr)*K*alpha*(x*x*x - 1.0) - 9.0*v0*K/(eta*eta);

        a -= (9.0*v0*K/(eta*eta))*(1.0 - eta*(1.0-x0))*exp(eta*(1.0-x0));
        a -= v0*(t-tr)*K*alpha*(x0*x0*x0 - 1.0) - 9.0*v0*K/(eta*eta);

        *g       += -a*10000.0 + p*v0*x*x*x - pr*v0*x0*x0*x0;
        *h       += -a*10000.0 + p*v0*x*x*x - pr*v0*x0*x0*x0 + 10000.0*t*alpha*K*v0*(x*x*x-x0*x0*x0);
        *s       += 10000.0*alpha*K*v0*(x*x*x-x0*x0*x0);
        *cp      += 10000.0*t*alpha*K*v0*3.0*(x*x*dxdt - x0*x0*dx0dt);
        *dcpdt   += 10000.0*alpha*K*v0*3.0*(x*x*dxdt - x0*x0*dx0dt)
        + 10000.0*t*alpha*K*v0*3.0*(2.0*x*dxdt*dxdt + x*x*d2xdt2 - 2.0*x0*dx0dt*dx0dt - x0*x0*d2x0dt2);
        *v       += v0*x*x*x;
        *dvdt    += 3.0*v0*x*x*dxdt;
        *dvdp    += 3.0*v0*x*x*dxdp/10000.0;
        *d2vdt2  += 3.0*v0*(2.0*x*dxdt*dxdt + x*x*d2xdt2);
        *d2vdtdp += 3.0*v0*(2.0*x*dxdt*dxdp + x*x*d2xdtdp)/10000.0;
        *d2vdp2  += 3.0*v0*(2.0*x*dxdp*dxdp + x*x*d2xdp2)/(10000.0*10000.0);
    } else if(phase->eos_type == EOS_SAXENA) {
        double a0 = phase->eos.Saxena.a0;
        double a1 = phase->eos.Saxena.a1;
        double a2 = phase->eos.Saxena.a2;
        double a3 = phase->eos.Saxena.a3;
        double b0 = phase->eos.Saxena.b0;
        double b1 = phase->eos.Saxena.b1;
        double b2 = phase->eos.Saxena.b2;
        double b3 = phase->eos.Saxena.b3;

        double alpha      = a0 + a1*t + a2/t + a3/(t*t);
        double dalphadt   = a1 - a2/(t*t) - 2.0*a3/(t*t*t);
        double d2alphadt2 = 2.0*a2/(t*t*t) + 6.0*a3/(t*t*t*t);
        double beta       = b0 + b1*t + b2*t*t + b3*t*t*t;
        double dbetaDt    = b1 + 2.0*b2*t + 3.0*b3*t*t;
        double d2betaDt2  = 2.0*b2 + 6.0*b3*t;
        double d3betaDt3  = 6.0*b3;
        double K          = 1.0/beta;
        double dKdt       = -dbetaDt/(beta*beta);
        double d2Kdt2     = 2.0*dbetaDt*dbetaDt/(beta*beta*beta) - d2betaDt2/(beta*beta);
        double d3Kdt3     = -6.0*dbetaDt*dbetaDt*dbetaDt/(beta*beta*beta*beta) + 6.0*dbetaDt*d2betaDt2/(beta*beta*beta) - d3betaDt3/(beta*beta);
        double Kp         = phase->eos.Saxena.dKdP + (phase->eos.Saxena.d2KdTdP)*(t-300.0)*log(t/300.0);
        double dKpdt      = (phase->eos.Saxena.d2KdTdP)*(log(t/300.0) + 1.0 - 300.0/t);
        double d2Kpdt2    = (phase->eos.Saxena.d2KdTdP)*(1.0/t + 300.0/(t*t));
        double d3Kpdt3    = (phase->eos.Saxena.d2KdTdP)*(-1.0/(t*t) - 600.0/(t*t*t));
        double vTPr       = (phase->v)*exp(a0*(t-tr)+(a1/2.0)*(t*t-tr*tr)+a2*log(t/tr)-a3*(1.0/t-1.0/tr));
        double dvTPrdt    =  vTPr*alpha;
        double d2vTPrdt2  =  vTPr*(dalphadt + alpha*alpha);
        double d3vTPrdt3  =  vTPr*alpha*(dalphadt + alpha*alpha) + vTPr*(d2alphadt2 + 2.0*alpha*dalphadt);

        int iter          = 0;
        const int iterMax = 1000;
        double vTP, vLast, fn , dfn;
        double x, dxdt, d2xdt2, d3xdt3, y, dydp, d2ydp2, dydt, d2ydt2, d3ydt3, d2ydtdp;
        double IntVdP, dIntVdPdt, d2IntVdPdt2, d3IntVdPdt3;
        double s1, s2, s3, s4, s5, s6, s7;

        /* Newton's method for volume from BM */
        vLast  = vTPr;
        vTP    = vTPr*0.99;
        while ((fabs(vTP- vLast) > 100.0*DBL_EPSILON) && (iter < iterMax)) {
            fn      = (3.0/2.0)*K*(pow(vTPr/vTP, (double) 7.0/3.0) - pow(vTPr/vTP, (double) 5.0/3.0))
            *(1.0-(3.0/4.0)*(4.0-Kp)*(pow(vTPr/vTP, (double) 2.0/3.0) - 1.0))
            - p;
            dfn     = (3.0/2.0)*K*((7.0/3.0)*pow(vTPr/vTP, (double) 4.0/3.0) - (5.0/3.0)*pow(vTPr/vTP, (double) 2.0/3.0))*(-vTPr/(vTP*vTP))
            *(1.0-(3.0/4.0)*(4.0-Kp)*(pow(vTPr/vTP, (double) 2.0/3.0) - 1.0))
            + (3.0/2.0)*K*(pow(vTPr/vTP, (double) 7.0/3.0) - pow(vTPr/vTP, (double) 5.0/3.0))
            *(-(3.0/4.0)*(4.0-Kp)*(2.0/3.0)*pow(vTPr/vTP, (double) -1.0/3.0))*(-vTPr/(vTP*vTP));
            vLast   = vTP;
            vTP    += -fn/dfn;
            if (vTP > 2.00*vTPr) vTP = 2.00*vTPr;
            if (vTP < 0.01*vTPr) vTP = 0.01*vTPr;
            iter++;
        }
        if (iter >= iterMax) {
            printf("Convergence error in Saxena Volume (BM) routine in function gibbs().\n");
            printf("  v = %g, dV = %g, f = %g at T = %g and P = %g in %d iterations.\n", vTP, fabs(vTP - vLast), fn, t, p, iter);
        }

        x      =  (3.0/4.0)*(4.0-Kp);
        dxdt   = -(3.0/4.0)*dKpdt;
        d2xdt2 = -(3.0/4.0)*d2Kpdt2;
        d3xdt3 = -(3.0/4.0)*d3Kpdt3;

        y = vTPr/vTP;
        /* This code is taken directly from Maple output */
        dydp    = -2.0/pow(y,(double) (2.0/3.0))/K/(-7.0*pow(y,(double) (2.0/3.0))+9.0*pow(y,(double) (4.0/3.0))*x-14.0*x*pow(y,(double) (2.0/3.0))+5.0+5.0*x);
        d2ydp2  = -2.0/3.0*dydp*dydp*(27.0*pow(y,(double) (4.0/3.0))*x-28.0*x*pow(y,(double) (2.0/3.0))-14.0*pow(y,(double) (2.0/3.0))+5.0+5.0*x)/y/(-7.0*pow(y,
                                                                                                                                                              (double) (2.0/3.0))+9.0*pow(y,(double) (4.0/3.0))*x-14.0*x*pow(y, (double) (2.0/3.0))+5.0+5.0*x);
        dydt    = -3.0*y*(-dKdt*pow(y,(double) (2.0/3.0))+dKdt*pow(y,(double) (4.0/3.0))*x-2.0*dKdt*pow(y,(double) (2.0/3.0))*x+dKdt+dKdt*x+pow(y,
                                                                                                                                                (double) (4.0/3.0))*dxdt*K-2.0*pow(y,(double) (2.0/3.0))*dxdt*K+dxdt*K)/K/(-7.0*pow(y,(double) (2.0/3.0))+9.0*pow(y,(double) (4.0/3.0))*x-14.0*x*pow(y,
                                                                                                                                                                                                                                                                                                     (double) (2.0/3.0))+5.0+5.0*x);
        d2ydt2  = -(9.0*d2Kdt2*y*y*x-18.0*d2Kdt2*pow(y,(double) (8.0/3.0))*x-9.0*pow(y,(double) (8.0/3.0))*d2Kdt2+10.0*K*dydt*dydt*x+18.0*dKdt*y*y*dxdt-36.0*dKdt
                    *pow(y,(double) (8.0/3.0))*dxdt+30.0*dKdt*y*dydt+54.0*pow(y, (double) (4.0/3.0))*dydt*dydt*x*K+18.0*pow(y,(double) (10.0/3.0))*dKdt*dxdt
                    -42.0*pow(y,(double) (5.0/3.0))*dKdt*dydt+9.0*pow(y,(double) (10.0/3.0))*d2Kdt2*x-18.0*pow(y,(double) (8.0/3.0))*d2xdt2*K+9.0*pow(y,
                                                                                                                                                      (double) (10.0/3.0))*d2xdt2*K+10.0*dydt*dydt*K-28.0*pow(y,(double) (2.0/3.0))*dydt*dydt*K+30.0*dKdt*y*dydt*x-84.0*dKdt*pow(y,(double) (5.0/3.0))*dydt*x
                    -56.0*K*dydt*dydt*pow(y,(double) (2.0/3.0))*x+30.0*K*dydt*dxdt*y-84.0*K*dydt*pow(y,(double) (5.0/3.0))*dxdt+54.0*pow(y,(double) (7.0/3.0))*dKdt*dydt*x+
                    9.0*d2xdt2*y*y*K+54.0*pow(y,(double) (7.0/3.0))*dxdt*dydt*K+9.0*d2Kdt2*y*y)/K/y/(-7.0*pow(y,(double) (2.0/3.0))+9.0*pow(y,(double) (4.0/3.0))*x-14.0*x*pow
                                                                                                     (y,(double) (2.0/3.0))+5.0+5.0*x)/3.0;
        d2ydtdp = -dydp*(-21.0*pow(y,(double) (5.0/3.0))*dKdt+27.0*pow(y,(double) (7.0/3.0))*dKdt*x-42.0*dKdt*pow(y,(double) (5.0/3.0))*x+15.0*dKdt*
                         y+15.0*dKdt*y*x+54.0*pow(y,(double) (4.0/3.0))*dydt*x*K-56.0*K*dydt*x*pow(y, (double) (2.0/3.0))+27.0*pow(y,(double) (7.0/3.0))*dxdt*K-42.0*K*pow(y,
                                                                                                                                                                           (double) (5.0/3.0))*dxdt+15.0*K*dxdt*y-28.0*pow(y,(double) (2.0/3.0))*dydt*K+10.0*dydt*K+10.0*K*dydt*x)/K/y/(-7.0*pow(y,(double) (2.0/3.0))+9.0*pow(y,
                                                                                                                                                                                                                                                                                                                               (double) (4.0/3.0))*x-14.0*x*pow(y,(double) (2.0/3.0))+5.0+5.0*x)/3.0;

        s2 = -10.0*dKdt*y*dydt*dydt+18.0*d2Kdt2*pow(y,(double) (11.0/3.0))*dxdt+28.0*dKdt*pow(y,(double) (5.0/3.0))*dydt*dydt
        +28.0/9.0*K*dydt*dydt*dydt*pow(y,(double) (2.0/3.0))+6.0*d3Kdt3*pow(y,(double) (11.0/3.0))*x-15.0*d2Kdt2*y*y*dydt
        +10.0/9.0*K*x*dydt*dydt*dydt+21.0*d2Kdt2*pow(y,(double) (8.0/3.0))*dydt
        +21.0*dKdt*pow(y,(double) (8.0/3.0))*d2ydt2-3.0*d3xdt3*y*y*y*K-3.0*pow(y,(double) (14.0/3.0))*d3Kdt3*x-9.0*d2Kdt2*y*y*y*dxdt
        -9.0*dKdt*y*y*y*d2xdt2+56.0*K*pow(y,(double) (5.0/3.0))*dxdt*dydt*dydt+18.0*dKdt*pow(y,(double) (11.0/3.0))*d2xdt2
        +56.0/9.0*K*dydt*dydt*dydt*pow(y,(double) (2.0/3.0))*x+6.0*pow(y,(double) (11.0/3.0))*d3xdt3*K-10.0*K*y*dydt*dydt*dxdt-15.0*K*y*y*d2ydt2*dxdt
        -3.0*d3Kdt3*y*y*y*x-10.0*dKdt*y*dydt*dydt*x+42.0*K*pow(y,(double) (8.0/3.0))*dxdt*d2ydt2-9.0*pow(y,(double) (14.0/3.0))*d2Kdt2*dxdt
        +56.0*dKdt*pow(y,(double) (5.0/3.0))*dydt*dydt*x-15.0*d2Kdt2*y*y*dydt*x;
        s1 = 42.0*d2Kdt2*pow(y,(double) (8.0/3.0))*dydt*x-15.0*K*dydt*d2xdt2*y*y-9.0*pow(y,(double) (14.0/3.0))*dKdt*d2xdt2
        -3.0*pow(y,(double) (14.0/3.0))*d3xdt3*K-30.0*dKdt*y*y*dydt*dxdt-10.0*dydt*d2ydt2*y*K+28.0*K*pow(y,(double) (5.0/3.0))*dydt*d2ydt2
        -27.0*pow(y,(double) (10.0/3.0))*d2xdt2*dydt*K-54.0*pow(y,(double) (7.0/3.0))*dxdt*dydt*dydt*K-27.0*pow(y,(double) (10.0/3.0))*dxdt*d2ydt2*K
        -18.0*x*dydt*dydt*dydt*pow(y,(double) (4.0/3.0))*K-15.0*dKdt*y*y*d2ydt2*x+42.0*dKdt*pow(y,(double) (8.0/3.0))*d2ydt2*x
        -27.0*pow(y,(double) (10.0/3.0))*dKdt*d2ydt2*x-54.0*pow(y,(double) (7.0/3.0))*dKdt*dydt*dydt*x-27.0*pow(y,(double) (10.0/3.0))*d2Kdt2*dydt*x
        +84.0*dKdt*pow(y,(double) (8.0/3.0))*dydt*dxdt-54.0*pow(y,(double) (10.0/3.0))*dKdt*dydt*dxdt-10.0*x*dydt*d2ydt2*y*K
        +56.0*K*pow(y,(double) (5.0/3.0))*dydt*d2ydt2*x-54.0*pow(y,(double) (7.0/3.0))*dydt*d2ydt2*x*K+42.0*K*dydt*pow(y,(double) (8.0/3.0))*d2xdt2
        -3.0*d3Kdt3*y*y*y+10.0/9.0*K*dydt*dydt*dydt+3.0*d3Kdt3*pow(y,(double) (11.0/3.0))-15.0*dKdt*y*y*d2ydt2+s2;
        s2 = 1/K/(y*y)/(-7.0*pow(y,(double) (2.0/3.0))+9.0*pow(y,(double) (4.0/3.0))*x-14.0*x*pow(y,(double) (2.0/3.0))+5.0+5.0*x);
        d3ydt3 = s1*s2;

        /* p = (3.0/2.0)*K*(pow(y, (double) 7.0/3.0) - pow(y, (double) 5.0/3.0))*(1.0-x*(pow(y, (double) 2.0/3.0) - 1.0)) */

        IntVdP = vTP*(p-pr) + (3.0/2.0)*K*vTPr*(  (3.0/4.0)*(1.0+2.0*x)*(pow(y, (double) 4.0/3.0) - 1.0)
                                                - (3.0/2.0)*(1.0+x)    *(pow(y, (double) 2.0/3.0) - 1.0)
                                                - (1.0/2.0)*x*(y*y - 1.0) );

        /* This code is taken directly from Maple output */
        dIntVdPdt = dvTPrdt*(p-1.0)/y-vTPr*(p-1.0)/(y*y)*dydt+3.0/2.0*dKdt*vTPr*(3.0/4.0*(1.0+2.0*x)*(pow(y,(double) (4.0/3.0))-1.0)
                                                                                 -3.0/2.0*(1.0+x)*(pow(y,(double) (2.0/3.0))-1.0)-x*(y*y-1.0)/2.0)+3.0/2.0*K*dvTPrdt*(3.0/4.0*(1.0+2.0*x)*(pow(y,(double) (4.0/3.0))-1.0)
                                                                                                                                                                      -3.0/2.0*(1.0+x)*(pow(y,(double) (2.0/3.0))-1.0)-x*(y*y-1.0)/2.0)+3.0/2.0*K*vTPr*(3.0/2.0*dxdt*(pow(y,(double) (4.0/3.0))-1.0)
                                                                                                                                                                                                                                                        +(1.0+2.0*x)*pow(y,(double) (1.0/3.0))*dydt-3.0/2.0*dxdt*(pow(y,(double) (2.0/3.0))-1.0)-(1.0+x)/pow(y,(double) (1.0/3.0))*dydt
                                                                                                                                                                                                                                                        -dxdt*(y*y-1.0)/2.0-x*y*dydt);

        s1 = d2vTPrdt2*(p-1.0)/y-2.0*dvTPrdt*(p-1.0)/(y*y)*dydt+2.0*vTPr*(p-1.0)/(y*y*y)*dydt*dydt-vTPr*(p-1.0)/(y*y)*d2ydt2
        +3.0/2.0*d2Kdt2*vTPr*(3.0/4.0*(1.0+2.0*x)*(pow(y,(double) (4.0/3.0))-1.0)-3.0/2.0*(1.0+x)*(pow(y,(double) (2.0/3.0))-1.0)-x*(y*y-1.0)/2.0);
        s2 = s1+3.0*dKdt*dvTPrdt*(3.0/4.0*(1.0+2.0*x)*(pow(y,(double) (4.0/3.0))-1.0)-3.0/2.0*(1.0+x)*(pow(y,(double) (2.0/3.0))-1.0)-x*(y*y-1.0)/2.0)
        +3.0*dKdt*vTPr*(3.0/2.0*dxdt*(pow(y,(double) (4.0/3.0))-1.0)+(1.0+2.0*x)*pow(y,(double) (1.0/3.0))*dydt-3.0/2.0*dxdt*(pow(y,(double) (2.0/3.0))-1.0)
                        -(1.0+x)/pow(y,(double) (1.0/3.0))*dydt-dxdt*(y*y-1.0)/2.0-x*y*dydt);
        d2IntVdPdt2 = s2+3.0/2.0*K*d2vTPrdt2*(3.0/4.0*(1.0+2.0*x)*(pow(y,(double) (4.0/3.0))-1.0)-3.0/2.0*(1.0+x)*(pow(y,(double) (2.0/3.0))-1.0)
                                              -x*(y*y-1.0)/2.0)+3.0*K*dvTPrdt*(3.0/2.0*dxdt*(pow(y,(double) (4.0/3.0))-1.0)+(1.0+2.0*x)*pow(y,(double) (1.0/3.0))*dydt
                                                                               -3.0/2.0*dxdt*(pow(y,(double) (2.0/3.0))-1.0)-(1.0+x)/pow(y,(double) (1.0/3.0))*dydt-dxdt*(y*y-1.0)/2.0-x*y*dydt)
        +3.0/2.0*K*vTPr*(3.0/2.0*d2xdt2*(pow(y,(double) (4.0/3.0))-1.0)+4.0*dxdt*pow(y,(double) (1.0/3.0))*dydt
                         +(1.0+2.0*x)/pow(y,(double) (2.0/3.0))*dydt*dydt/3.0+(1.0+2.0*x)*pow(y,(double) (1.0/3.0))*d2ydt2-3.0/2.0*d2xdt2*(pow(y,(double) (2.0/3.0))-1.0)
                         -2.0*dxdt/pow(y,(double) (1.0/3.0))*dydt+(1.0+x)/pow(y,(double) (4.0/3.0))*dydt*dydt/3.0-(1.0+x)/pow(y,(double) (1.0/3.0))*d2ydt2
                         -d2xdt2*(y*y-1.0)/2.0-2.0*dxdt*y*dydt-x*dydt*dydt-x*y*d2ydt2);

        s1 = -3.0*d2vTPrdt2*(p-1.0)/(y*y)*dydt+d3vTPrdt3*(p-1.0)/y-3.0*dvTPrdt*(p-1.0)/(y*y)*d2ydt2-vTPr*(p-1.0)/(y*y)*d3ydt3
        +9.0/2.0*d2Kdt2*vTPr*(3.0/2.0*dxdt*(pow(y,(double) (4.0/3.0))-1.0)+(1.0+2.0*x)*pow(y,(double) (1.0/3.0))*dydt
                              -3.0/2.0*dxdt*(pow(y,(double) (2.0/3.0))-1.0)-(1.0+x)/pow(y,(double) (1.0/3.0))*dydt-dxdt*(y*y-1.0)/2.0-x*y*dydt)
        +9.0/2.0*d2Kdt2*dvTPrdt*(3.0/4.0*(1.0+2.0*x)*(pow(y,(double) (4.0/3.0))-1.0)-3.0/2.0*(1.0+x)*(pow(y,(double) (2.0/3.0))-1.0)
                                 -x*(y*y-1.0)/2.0)+3.0/2.0*d3Kdt3*vTPr*(3.0/4.0*(1.0+2.0*x)*(pow(y,(double) (4.0/3.0))-1.0)-3.0/2.0*(1.0+x)*(pow(y,(double) (2.0/3.0))-1.0)
                                                                        -x*(y*y-1.0)/2.0)+9.0/2.0*K*d2vTPrdt2*(3.0/2.0*dxdt*(pow(y,(double) (4.0/3.0))-1.0)+(1.0+2.0*x)*pow(y,(double) (1.0/3.0))*dydt
                                                                                                               -3.0/2.0*dxdt*(pow(y,(double) (2.0/3.0))-1.0)-(1.0+x)/pow(y,(double) (1.0/3.0))*dydt-dxdt*(y*y-1.0)/2.0-x*y*dydt);
        s2 = s1+3.0/2.0*K*d3vTPrdt3*(3.0/4.0*(1.0+2.0*x)*(pow(y,(double) (4.0/3.0))-1.0)-3.0/2.0*(1.0+x)*(pow(y,(double) (2.0/3.0))-1.0)
                                     -x*(y*y-1.0)/2.0)+9.0/2.0*dKdt*vTPr*(3.0/2.0*d2xdt2*(pow(y,(double) (4.0/3.0))-1.0)+4.0*dxdt*pow(y,(double) (1.0/3.0))*dydt
                                                                          +(1.0+2.0*x)/pow(y,(double) (2.0/3.0))*dydt*dydt/3.0+(1.0+2.0*x)*pow(y,(double) (1.0/3.0))*d2ydt2-3.0/2.0*d2xdt2*(pow(y,(double) (2.0/3.0))-1.0)
                                                                          -2.0*dxdt/pow(y,(double) (1.0/3.0))*dydt+(1.0+x)/pow(y,(double) (4.0/3.0))*dydt*dydt/3.0-(1.0+x)/pow(y,(double) (1.0/3.0))*d2ydt2
                                                                          -d2xdt2*(y*y-1.0)/2.0-2.0*dxdt*y*dydt-x*dydt*dydt-x*y*d2ydt2)+9.0*dKdt*dvTPrdt*(3.0/2.0*dxdt*(pow(y,(double) (4.0/3.0))-1.0)
                                                                                                                                                          +(1.0+2.0*x)*pow(y,(double) (1.0/3.0))*dydt-3.0/2.0*dxdt*(pow(y,(double) (2.0/3.0))-1.0)-(1.0+x)/pow(y,(double) (1.0/3.0))*dydt
                                                                                                                                                          -dxdt*(y*y-1.0)/2.0-x*y*dydt)+9.0/2.0*dKdt*d2vTPrdt2*(3.0/4.0*(1.0+2.0*x)*(pow(y,(double) (4.0/3.0))-1.0)
                                                                                                                                                                                                                -3.0/2.0*(1.0+x)*(pow(y,(double) (2.0/3.0))-1.0)-x*(y*y-1.0)/2.0);
        s4 = s2;
        s6 = 9.0/2.0*K*dvTPrdt*(3.0/2.0*d2xdt2*(pow(y,(double) (4.0/3.0))-1.0)+4.0*dxdt*pow(y,(double) (1.0/3.0))*dydt
                                +(1.0+2.0*x)/pow(y,(double) (2.0/3.0))*dydt*dydt/3.0+(1.0+2.0*x)*pow(y,(double) (1.0/3.0))*d2ydt2-3.0/2.0*d2xdt2*(pow(y,(double) (2.0/3.0))-1.0)
                                -2.0*dxdt/pow(y,(double) (1.0/3.0))*dydt+(1.0+x)/pow(y,(double) (4.0/3.0))*dydt*dydt/3.0-(1.0+x)/pow(y,(double) (1.0/3.0))*d2ydt2
                                -d2xdt2*(y*y-1.0)/2.0-2.0*dxdt*y*dydt-x*dydt*dydt-x*y*d2ydt2);
        s7 = 3.0/2.0*K*vTPr*(-2.0/9.0*(1.0+2.0*x)/pow(y,(double) (5.0/3.0))*dydt*dydt*dydt+2.0*dxdt/pow(y,(double) (2.0/3.0))*dydt*dydt
                             -4.0/9.0*(1.0+x)/pow(y,(double) (7.0/3.0))*dydt*dydt*dydt+3.0/2.0*d3xdt3*(pow(y,(double) (4.0/3.0))-1.0)
                             +dxdt/pow(y,(double) (4.0/3.0))*dydt*dydt+(1.0+2.0*x)/pow(y,(double) (2.0/3.0))*dydt*d2ydt2+(1.0+x)/pow(y,(double) (4.0/3.0))*dydt*d2ydt2
                             -d3xdt3*(y*y-1.0)/2.0-3.0/2.0*d3xdt3*(pow(y,(double) (2.0/3.0))-1.0)+6.0*dxdt*pow(y,(double) (1.0/3.0))*d2ydt2
                             +6.0*d2xdt2*pow(y,(double) (1.0/3.0))*dydt-3.0*dxdt/pow(y,(double) (1.0/3.0))*d2ydt2-3.0*d2xdt2/pow(y,(double) (1.0/3.0))*dydt
                             +(1.0+2.0*x)*pow(y,(double) (1.0/3.0))*d3ydt3-3.0*d2xdt2*y*dydt-(1.0+x)/pow(y,(double) (1.0/3.0))*d3ydt3-x*y*d3ydt3-3.0*x*dydt*d2ydt2
                             -3.0*dxdt*y*d2ydt2-3.0*dxdt*dydt*dydt);
        s5 = s6+s7;
        s3 = s4+s5;
        d3IntVdPdt3 = s3+6.0*dvTPrdt*(p-1.0)/(y*y*y)*dydt*dydt+6.0*vTPr*(p-1.0)/(y*y*y)*dydt*d2ydt2-6.0*vTPr*(p-1.0)/(y*y*y*y)*dydt*dydt*dydt;

        *g       += IntVdP;
        *h       += IntVdP -t*dIntVdPdt;
        *s       += -dIntVdPdt;
        *cp      += -t*d2IntVdPdt2;
        *dcpdt   += -d2IntVdPdt2 - t*d3IntVdPdt3;
        *v       += vTP;
        *dvdt    += (dvTPrdt - vTP*dydt)/y;
        *dvdp    += -vTP*dydp/y;
        *d2vdt2  += (y*d2vTPrdt2 - vTPr*d2ydt2)/(y*y) - 2.0*dydt*(y*dvTPrdt - vTPr*dydt)/(y*y*y);
        *d2vdtdp += -dvTPrdt*dydp/(y*y) + 2.0*vTPr*dydt*dydp/(y*y*y) - vTPr*d2ydtdp/(y*y);
        *d2vdp2  += vTPr*(2.0*dydp*dydp/(y*y*y) - d2ydp2/(y*y));
    }
}

static void HPproperties(double t, double p, double h298, double s298, double v298, double a, double b, double c, double d,
                         double a0, double K298, double Tc0, double Smax, double Vmax, double *gs, double *ss, double *hs, double *cps,
                         double *dcpsdt, double *vs, double *dvsdt, double *dvsdp, double *d2vsdt2, double *d2vsdtdp, double *d2vsdp2) {
    const double tr = 298.15;
    double v1T      = v298*(1.0+a0*(t-tr) - 20.0*a0*(sqrt(t)-sqrt(tr)));
    double dv1Tdt   = v298*(a0 - 10.0*a0/sqrt(t));
    double d2v1Tdt2 = v298*5.0*a0/pow(t, 3.0/2.0);
    double d3v1Tdt3 = -3.0*v298*5.0*a0/pow(t, 5.0/2.0)/2.0;
    double KT    = K298*(1.0 - 1.5e-4*(t-tr));
    double dKTdt = -1.5e-4*K298;

    double vInt      = pow(1.0 + 4.0*p/KT, 3.0/4.0)/3.0 - 1.0/3.0;
    double dvIntdt   = -p*dKTdt/KT/KT/pow(1.0 + 4.0*p/KT, 1.0/4.0);
    double d2vIntdt2 =  2.0*p*dKTdt*dKTdt/KT/KT/KT/pow(1.0 + 4.0*p/KT, 1.0/4.0) - p*p*dKTdt*dKTdt/KT/KT/KT/KT/pow(1.0 + 4.0*p/KT, 5.0/4.0);
    double d3vIntdt3 =  - 6.0*p*dKTdt*dKTdt*dKTdt/KT/KT/KT/KT/pow(1.0 + 4.0*p/KT, 1.0/4.0)
    + 2.0*p*p*dKTdt*dKTdt*dKTdt/KT/KT/KT/KT/KT/pow(1.0 + 4.0*p/KT, 5.0/4.0)
    + 4.0*p*p*dKTdt*dKTdt*dKTdt/KT/KT/KT/KT/KT/pow(1.0 + 4.0*p/KT, 5.0/4.0)
    - 5.0*p*p*p*dKTdt*dKTdt*dKTdt/KT/KT/KT/KT/KT/KT/pow(1.0 + 4.0*p/KT, 9.0/4.0);

    *vs = v1T*pow(1.0 - 4.0*p/(KT + 4.0*p), 1.0/4.0);
    *dvsdt = dv1Tdt*pow(1.0 - 4.0*p/(KT + 4.0*p), 1.0/4.0)
    + v1T*(4.0*p/(KT+4.0*p)/(KT+4.0*p))*dKTdt/4.0/pow(1.0 - 4.0*p/(KT + 4.0*p), 3.0/4.0);
    *dvsdp = v1T*(-4.0*KT/(KT+4.0*p)/(KT+4.0*p))/4.0/pow(1.0 - 4.0*p/(KT + 4.0*p), 3.0/4.0);
    *d2vsdt2 = d2v1Tdt2*pow(1.0 - 4.0*p/(KT + 4.0*p), 1.0/4.0)
    + 2.0*dv1Tdt*(4.0*p/(KT+4.0*p)/(KT+4.0*p))*dKTdt/4.0/pow(1.0 - 4.0*p/(KT + 4.0*p), 3.0/4.0)
    - 3.0*v1T*pow(4.0*p*dKTdt/(KT+4.0*p)/(KT+4.0*p), 2.0)/4.0/4.0/pow(1.0 - 4.0*p/(KT + 4.0*p), 7.0/4.0)
    + v1T*(-8.0*p/(KT+4.0*p)/(KT+4.0*p)/(KT+4.0*p))*dKTdt*dKTdt/4.0/pow(1.0 - 4.0*p/(KT + 4.0*p), 3.0/4.0);
    *d2vsdtdp = dv1Tdt*(-4.0*KT/(KT+4.0*p)/(KT+4.0*p))/4.0/pow(1.0 - 4.0*p/(KT + 4.0*p), 3.0/4.0)
    - 3.0*v1T*(-4.0*KT/(KT+4.0*p)/(KT+4.0*p))*(4.0*p/(KT+4.0*p)/(KT+4.0*p))*dKTdt/4.0/4.0
    /pow(1.0 - 4.0*p/(KT + 4.0*p), 7.0/4.0)
    + v1T*(4.0/(KT+4.0*p)/(KT+4.0*p) - 32.0*p/(KT+4.0*p)/(KT+4.0*p)/(KT+4.0*p))*dKTdt/4.0
    /pow(1.0 - 4.0*p/(KT + 4.0*p), 3.0/4.0);
    *d2vsdp2 = - 3.0*v1T*pow(-4.0*KT/(KT+4.0*p)/(KT+4.0*p), 2.0)/4.0/4.0/pow(1.0 - 4.0*p/(KT + 4.0*p), 7.0/4.0)
    + v1T*(32.0*KT/(KT+4.0*p)/(KT+4.0*p)/(KT+4.0*p))/4.0/pow(1.0 - 4.0*p/(KT + 4.0*p), 3.0/4.0);

    *gs = v1T*KT*vInt;
    *ss = dv1Tdt*KT*vInt + v1T*dKTdt*vInt + v1T*KT*dvIntdt;
    *ss *= -1.0;
    *hs = *gs + t*(*ss);
    *cps = d2v1Tdt2*KT*vInt + 2.0*dv1Tdt*dKTdt*vInt + 2.0*dv1Tdt*KT*dvIntdt + 2.0*v1T*dKTdt*dvIntdt + v1T*KT*d2vIntdt2;
    *dcpsdt = d3v1Tdt3*KT*vInt + 3.0*d2v1Tdt2*dKTdt*vInt + 6.0*dv1Tdt*dKTdt*dvIntdt + 3.0*d2v1Tdt2*KT*dvIntdt  + 3.0*dv1Tdt*KT*d2vIntdt2
    + 3.0*v1T*dKTdt*d2vIntdt2 + v1T*KT*d3vIntdt3;
    *dcpsdt *= -t;
    *dcpsdt -= *cps;
    *cps *= -t;

    *cps    += a + b*t + c/t/t + d/sqrt(t);
    *dcpsdt += b - 2.0*c/t/t/t - d/2.0/pow(t, 3.0/2.0);
    *hs     += h298 + a*(t-tr) + b*(t*t-tr*tr)/2.0 - c*(1.0/t-1.0/tr) + 2.0*d*(sqrt(t)-sqrt(tr));
    *ss     += s298 + a*log(t/tr) + b*(t-tr) - c*(1.0/t/t-1.0/tr/tr)/2.0 - 2.0*d*(1.0/sqrt(t)-1.0/sqrt(tr));
    *gs = *hs - t*(*ss);

    if ((Tc0 > 0.0) && (Smax > 0.0)) {
        double Tc = Tc0 + Vmax*p/Smax;
        if (t < Tc) {
            double dTcdp     = Vmax/Smax;
            double Q2at298   = sqrt(1.0 - tr/Tc0);
            double Q2        = sqrt(1.0 - t/Tc);
            double dQ2dt     = -1.0/sqrt(1.0 - t/Tc)/2.0/Tc;
            double d2Q2dt2   = -1.0/pow(1.0 - t/Tc, 3.0/2.0)/2.0/2.0/Tc/Tc;
            double d3Q2dt3   = -3.0/pow(1.0 - t/Tc, 5.0/2.0)/2.0/2.0/2.0/Tc/Tc/Tc;
            double dQ2dp     = t*dTcdp/sqrt(1.0 - t/Tc)/2.0/Tc/Tc;
            double d2Q2dtdp  = -t*dTcdp/pow(1.0 - t/Tc, 3.0/2.0)/2.0/2.0/Tc/Tc + dTcdp/pow(1.0 - t/Tc, 1.0/2.0)/2.0/Tc/Tc;
            double d2Q2dp2   = -t*t*dTcdp*dTcdp/pow(1.0 - t/Tc, 3.0/2.0)/2.0/2.0/Tc/Tc/Tc/Tc - 2.0*t*dTcdp*dTcdp/sqrt(1.0 - t/Tc)/2.0/Tc/Tc/Tc;
            double d3Q2dp3   = 3.0*t*t*t*dTcdp*dTcdp*dTcdp/pow(1.0 - t/Tc, 5.0/2.0)/2.0/2.0/2.0/Tc/Tc/Tc/Tc/Tc/Tc
            + 4.0*t*t*dTcdp*dTcdp*dTcdp/pow(1.0 - t/Tc, 3.0/2.0)/2.0/2.0/Tc/Tc/Tc/Tc/Tc
            + t*t*dTcdp*dTcdp*dTcdp/pow(1.0 - t/Tc, 3.0/2.0)/2.0/Tc/Tc/Tc/Tc/Tc
            + 6.0*t*dTcdp*dTcdp*dTcdp/sqrt(1.0 - t/Tc)/2.0/Tc/Tc/Tc/Tc;
            double d3Q2dt2dp = -dTcdp/pow(1.0 - t/Tc, 3.0/2.0)/2.0/2.0/Tc/Tc - 3.0*t*dTcdp/pow(1.0 - t/Tc, 5.0/2.0)/2.0/2.0/2.0/Tc/Tc/Tc
            + dTcdp/pow(1.0 - t/Tc, 3.0/2.0)/2.0/2.0/Tc/Tc/Tc;
            double d3Q2dtdp2 = 3.0*t*t*dTcdp*dTcdp/pow(1.0 - t/Tc, 5.0/2.0)/2.0/2.0/2.0/Tc/Tc/Tc/Tc + 2.0*t*dTcdp*dTcdp/pow(1.0 - t/Tc, 3.0/2.0)/2.0/2.0/Tc/Tc/Tc
            - t*dTcdp*dTcdp/pow(1.0 - t/Tc, 3.0/2.0)/2.0/2.0/Tc/Tc/Tc/Tc - 2.0*dTcdp*dTcdp/pow(1.0 - t/Tc, 1.0/2.0)/2.0/Tc/Tc/Tc;

            double hP298     = Smax*Tc0*(Q2at298 - Q2at298*Q2at298*Q2at298/3.0);
            double sP298     = Smax*Q2at298;

            double vPT       = Vmax*Q2at298*(1.0+a0*(t-tr) - 20.0*a0*(sqrt(t)-sqrt(tr)));
            double dvPTdt    = Vmax*Q2at298*(a0 - 10.0*a0/sqrt(t));
            double d2vPTdt2  = Vmax*Q2at298*5.0*a0/pow(t, 3.0/2.0);
            double d3vPTdt3  = -3.0*Vmax*Q2at298*5.0*a0/pow(t, 5.0/2.0)/2.0;

            double vIntegral    = vPT*KT*vInt;
            double dvIntegraldt = dvPTdt*KT*vInt + vPT*dKTdt*vInt + vPT*KT*dvIntdt;
            double d2vIntegraldt2 = d2vPTdt2*KT*vInt + 2.0*dvPTdt*dKTdt*vInt + 2.0*dvPTdt*KT*dvIntdt + 2.0*vPT*dKTdt*dvIntdt + vPT*KT*d2vIntdt2;
            double d3vIntegraldt3 = d3vPTdt3*KT*vInt + 3.0*d2vPTdt2*dKTdt*vInt + 6.0*dvPTdt*dKTdt*dvIntdt + 3.0*d2vPTdt2*KT*dvIntdt  + 3.0*dvPTdt*KT*d2vIntdt2
            + 3.0*vPT*dKTdt*d2vIntdt2 + vPT*KT*d3vIntdt3;

            double Glandau      = Smax*((t-Tc)*Q2 + Tc*Q2*Q2*Q2/3.0);
            double dGlandaudt   = Smax*(Q2 + (t-Tc)*dQ2dt + Tc*Q2*Q2*dQ2dt);
            double d2Glandaudt2 = Smax*(2.0*dQ2dt + (t-Tc)*d2Q2dt2 + 2.0*Tc*Q2*dQ2dt*dQ2dt + Tc*Q2*Q2*d2Q2dt2);
            double d3Glandaudt3 = Smax*(3.0*d2Q2dt2 + (t-Tc)*d3Q2dt3 + 2.0*Tc*dQ2dt*dQ2dt*dQ2dt + 6.0*Tc*Q2*dQ2dt*d2Q2dt2 + Tc*Q2*Q2*d3Q2dt3);

            double dGlandaudp     = Smax*(-dTcdp*Q2 + (t-Tc)*dQ2dp + dTcdp*Q2*Q2*Q2/3.0 + Tc*Q2*Q2*dQ2dp);
            double d2Glandaudtdp  = Smax*(-dTcdp*dQ2dt + dQ2dp + (t-Tc)*d2Q2dtdp + dTcdp*Q2*Q2*dQ2dt + 2.0*Tc*Q2*dQ2dt*dQ2dp + Tc*Q2*Q2*d2Q2dtdp);
            double d2Glandaudp2   = Smax*(-dTcdp*dQ2dp -dTcdp*dQ2dp + (t-Tc)*d2Q2dtdp + dTcdp*Q2*Q2*dQ2dp/3.0 + dTcdp*Q2*Q2*dQ2dp + 2.0*Tc*Q2*dQ2dp*dQ2dp + Tc*Q2*Q2*d2Q2dp2);
            double d3Glandaudt2dp = Smax*(-dTcdp*d2Q2dt2 + d2Q2dtdp + d2Q2dtdp + (t-Tc)*d3Q2dt2dp + 2.0*dTcdp*Q2*dQ2dt*dQ2dt + dTcdp*Q2*Q2*d2Q2dt2
                                          + 2.0*Tc*dQ2dt*dQ2dt*dQ2dp + 2.0*Tc*Q2*d2Q2dt2*dQ2dp + 2.0*Tc*Q2*dQ2dt*d2Q2dtdp + 2.0*Tc*Q2*dQ2dt*d2Q2dtdp + + Tc*Q2*Q2*d3Q2dt2dp);
            double d3Glandaudtdp2 = Smax*(-dTcdp*d2Q2dtdp + d2Q2dp2 + (t-Tc)*d3Q2dtdp2 + 2.0*dTcdp*Q2*dQ2dp*dQ2dt + dTcdp*Q2*Q2*d2Q2dtdp
                                          + 2.0*dTcdp*Q2*dQ2dt*dQ2dp + 2.0*Tc*dQ2dp*dQ2dt*dQ2dp + 2.0*Tc*Q2*d2Q2dtdp*dQ2dp + 2.0*Tc*Q2*dQ2dt*d2Q2dp2
                                          + dTcdp*Q2*Q2*d2Q2dtdp + 2.0*Tc*Q2*dQ2dp*d2Q2dtdp + Tc*Q2*Q2*d3Q2dtdp2);
            double d3Glandaudp3   = Smax*(-dTcdp*d2Q2dp2 -dTcdp*d2Q2dp2 + (t-Tc)*d3Q2dtdp2 + 2.0*dTcdp*Q2*dQ2dp*dQ2dp/3.0 + dTcdp*Q2*Q2*d2Q2dp2/3.0
                                          + 2.0*dTcdp*Q2*dQ2dp*dQ2dp + dTcdp*Q2*Q2*d2Q2dp2 + 2.0*dTcdp*Q2*dQ2dp*dQ2dp + 2.0*Tc*dQ2dp*dQ2dp*dQ2dp + 4.0*Tc*Q2*dQ2dp*d2Q2dp2
                                          + dTcdp*Q2*Q2*d2Q2dp2 + 2.0*Tc*Q2*dQ2dp*d2Q2dp2 + Tc*Q2*Q2*d3Q2dp3);

            *gs     += hP298 - t*sP298 + vIntegral + Glandau;
            *ss     += sP298 -(dvIntegraldt + dGlandaudt);
            *hs     += hP298 + vIntegral + Glandau -t*(dvIntegraldt + dGlandaudt);
            *cps    += -t*(d2vIntegraldt2 + d2Glandaudt2);
            *dcpsdt += -(d2vIntegraldt2 + d2Glandaudt2) - t*(d3vIntegraldt3 + d3Glandaudt3);

            *vs       += vPT*pow(1.0 - 4.0*p/(KT + 4.0*p), 1.0/4.0)
            + dGlandaudp;
            *dvsdt    += dvPTdt*pow(1.0 - 4.0*p/(KT + 4.0*p), 1.0/4.0)
            + vPT*(4.0*p/(KT+4.0*p)/(KT+4.0*p))*dKTdt/4.0/pow(1.0 - 4.0*p/(KT + 4.0*p), 3.0/4.0)
            + d2Glandaudtdp;
            *dvsdp    += vPT*(-4.0*KT/(KT+4.0*p)/(KT+4.0*p))/4.0/pow(1.0 - 4.0*p/(KT + 4.0*p), 3.0/4.0)
            + d2Glandaudp2;
            *d2vsdt2  += d2vPTdt2*pow(1.0 - 4.0*p/(KT + 4.0*p), 1.0/4.0)
            + 2.0*dvPTdt*(4.0*p/(KT+4.0*p)/(KT+4.0*p))*dKTdt/4.0/pow(1.0 - 4.0*p/(KT + 4.0*p), 3.0/4.0)
            - 3.0*vPT*pow(4.0*p*dKTdt/(KT+4.0*p)/(KT+4.0*p), 2.0)/4.0/4.0/pow(1.0 - 4.0*p/(KT + 4.0*p), 7.0/4.0)
            + vPT*(-8.0*p/(KT+4.0*p)/(KT+4.0*p)/(KT+4.0*p))*dKTdt*dKTdt/4.0/pow(1.0 - 4.0*p/(KT + 4.0*p), 3.0/4.0)
            + d3Glandaudt2dp;
            *d2vsdtdp += dvPTdt*(-4.0*KT/(KT+4.0*p)/(KT+4.0*p))/4.0/pow(1.0 - 4.0*p/(KT + 4.0*p), 3.0/4.0)
            - 3.0*vPT*(-4.0*KT/(KT+4.0*p)/(KT+4.0*p))*(4.0*p/(KT+4.0*p)/(KT+4.0*p))*dKTdt/4.0/4.0
            /pow(1.0 - 4.0*p/(KT + 4.0*p), 7.0/4.0)
            + vPT*(4.0/(KT+4.0*p)/(KT+4.0*p) - 32.0*p/(KT+4.0*p)/(KT+4.0*p)/(KT+4.0*p))*dKTdt/4.0
            /pow(1.0 - 4.0*p/(KT + 4.0*p), 3.0/4.0)
            + d3Glandaudtdp2;
            *d2vsdp2  += - 3.0*vPT*pow(-4.0*KT/(KT+4.0*p)/(KT+4.0*p), 2.0)/4.0/4.0/pow(1.0 - 4.0*p/(KT + 4.0*p), 7.0/4.0)
            + vPT*(32.0*KT/(KT+4.0*p)/(KT+4.0*p)/(KT+4.0*p))/4.0/pow(1.0 - 4.0*p/(KT + 4.0*p), 3.0/4.0)
            + d3Glandaudp3;

        }
    }
}

/*
 *=============================================================================
 * Public function
 */

void gibbs(double t, double p, char *name, ThermoRef *phase,
           ThermoLiq *liquid, ThermoData *fusion, ThermoData *result)
{
    static const double r = 8.3143;
    static const double tr = 298.15;
    static const double pr = 1.0;
    static const double trl = 1673.0;

    double k0, k1, k2, k3, v1 = 0.0, v2 = 0.0, v3 = 0.0, v4 = 0.0, vR = 0.0, cp_t, cp_h, l1, l2;
    double gs       = 0.0,
    hs       = 0.0,
    ss       = 0.0,
    vs       = 0.0,
    cps      = 0.0,
    dcpsdt   = 0.0,
    dvsdt    = 0.0,
    dvsdp    = 0.0,
    d2vsdt2  = 0.0,
    d2vsdp2  = 0.0,
    d2vsdtdp = 0.0;

    /* Liquid : Special cases handled first!                                  */
    if (liquid != NULL) {
        double hl, sl, cpl, gl = 0.0, dcpldt, vl = 0.0, dvldt = 0.0, dvldp = 0.0, d2vldt2 = 0.0,
        d2vldp2 = 0.0, d2vldtdp = 0.0;

        /* special case - no EOS option */
        if ((strcmp(name, "SiO2") == 0) && (
#ifndef TRUE_xMELTS
                                            (calculationMode == MODE_xMELTS) ||
#endif
                                            (calculationMode == MODE__MELTS) ||
                                            (calculationMode == MODE__MELTSandCO2) ||
                                            (calculationMode == MODE__MELTSandCO2_H2O)) ) {
            double h0_sio2 = -901554.0;
            double s0_sio2 = 48.475;
            double al_sio2 = 127.200;
            double bl_sio2 = -10.777e-3;
            double cl_sio2 = 4.3127e5;
            double dl_sio2 = -1463.8;
            double tg_sio2 = 1480.0;
            double cp_sio2 = 81.373;

            if (t >= tg_sio2) {
                hl = h0_sio2
                + al_sio2*(tg_sio2 - tr)
                + bl_sio2*(SQUARE(tg_sio2)-tr*tr)/2.0
                - cl_sio2*(1.0/tg_sio2-1.0/tr)
                + 2.0 *dl_sio2*(sqrt(tg_sio2)-sqrt(tr));
                sl = s0_sio2
                + al_sio2*log(tg_sio2/tr)
                + bl_sio2*(tg_sio2 - tr)
                - (cl_sio2/2.0)*(1.0/SQUARE(tg_sio2) - 1.0/SQUARE(tr))
                - 2.0*dl_sio2*(1.0/sqrt(tg_sio2)-1.0/sqrt(tr));
                hl = hl + (t - tg_sio2)*cp_sio2;
                sl = sl + cp_sio2*log(t/tg_sio2);
                cpl = cp_sio2;
                dcpldt = 0.0;
            } else {
                hl = h0_sio2 + al_sio2*(t-tr) + bl_sio2*(t*t-tr*tr)/2.0
                - cl_sio2*(1.0/t - 1.0/tr) + 2.0*dl_sio2*(sqrt(t)-sqrt(tr));
                sl = s0_sio2 + al_sio2*log(t/tr) + bl_sio2*(t - tr)
                - (cl_sio2/2.0)*(1.0/SQUARE(t) - 1.0/SQUARE(tr))
                - 2.0*dl_sio2*(1.0/sqrt(t)-1.0/sqrt(tr));
                cpl = al_sio2 + bl_sio2*t + cl_sio2/(t*t) + dl_sio2/sqrt(t);
                dcpldt = bl_sio2 - 2.0*cl_sio2/CUBE(t)
                - 0.5*dl_sio2/pow(t, (double) 1.5);
            }

            gl = hl - t*sl
            + (liquid->v + liquid->eos.Kress.dvdt*(t-trl))*(p-pr)
            + 0.5*(liquid->eos.Kress.dvdp + (t-trl)*liquid->eos.Kress.d2vdtp)*(p*p-pr*pr)
            - (liquid->eos.Kress.dvdp + (t-trl)*liquid->eos.Kress.d2vdtp)*pr*(p-pr)
            + liquid->eos.Kress.d2vdp2*((CUBE(p)-CUBE(pr))/6.0 - pr*(p*p-pr*pr)/2.0 + pr*pr*(p-pr)/2.0 );
            hl += (liquid->v + liquid->eos.Kress.dvdt*(t-trl))*(p-pr)
            + 0.5*(liquid->eos.Kress.dvdp + (t-trl)*liquid->eos.Kress.d2vdtp)*(p*p-pr*pr)
            - (liquid->eos.Kress.dvdp + (t-trl)*liquid->eos.Kress.d2vdtp)*pr*(p-pr)
            + liquid->eos.Kress.d2vdp2*( (CUBE(p)-CUBE(pr))/6.0 - pr*(p*p-pr*pr)/2.0 + pr*pr*(p-pr)/2.0 )
            - t*(liquid->eos.Kress.dvdt*(p-pr) + 0.5*liquid->eos.Kress.d2vdtp*(p-pr)*(p-pr));
            sl += -(liquid->eos.Kress.dvdt*(p-pr) + 0.5*liquid->eos.Kress.d2vdtp*(p-pr)*(p-pr));

            vl	  = liquid->v + liquid->eos.Kress.dvdt*(t-trl)
            + (liquid->eos.Kress.dvdp + liquid->eos.Kress.d2vdtp*(t-trl))*(p-pr)
            + liquid->eos.Kress.d2vdp2*(0.5*p*p - pr*(p-pr));
            dvldt    = liquid->eos.Kress.dvdt + liquid->eos.Kress.d2vdtp*(p-pr);
            dvldp    = liquid->eos.Kress.dvdp + liquid->eos.Kress.d2vdtp*(t-trl)
            + liquid->eos.Kress.d2vdp2*(p-pr);
            d2vldt2  = 0.0;
            d2vldp2  = liquid->eos.Kress.d2vdp2;
            d2vldtdp = liquid->eos.Kress.d2vdtp;

            /* special case - no EOS option */
        } else if( (strcmp(name, "H2O") == 0) && (
#ifdef TRUE_xMELTS
						  (calculationMode == MODE_xMELTS) ||
#endif
						  (calculationMode == MODE_pMELTS)) ) {
            double a = phase->h/r;
            double b = -phase->s/r;
            double gH2O, hH2O, sH2O, cpH2O, dcpdtH2O, vH2O, dvdtH2O, dvdpH2O,
            d2vdt2H2O, d2vdtdpH2O, d2vdp2H2O;
            double x[2] = { 1.0, 0.0};

            fluidPhase(t, 9550.0, x, &gH2O, NULL, NULL, NULL, &hH2O, &sH2O,
                       NULL, NULL, &cpH2O, &dcpdtH2O, NULL, &vH2O, NULL, NULL, &dvdtH2O,
                       &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O, NULL, NULL, NULL, NULL,
                       NULL);

            gl     = gH2O + r*t*(a/t + b);
            hl     = hH2O + r*a;
            sl     = sH2O - r*b;
            cpl    = cpH2O;
            dcpldt = dcpdtH2O;

            /* special case - no EOS option */
        } else if( (strcmp(name, "H2O") == 0) && ((calculationMode == MODE__MELTS) ||
                                                  (calculationMode == MODE__MELTSandCO2) )) {
            double a = -33676.0 + phase->h/r, b = 18.3527 - phase->s/r;
            double phiP = (0.110/t + 4.432e-5 + 1.405e-7*t - 2.394e-11*t*t)*p
            + (7.337e-8/t - 1.170e-8 - 9.502e-13*t)*p*p
            + (1.876e-10/t + 4.586e-13)*CUBE(p) - 1.191e-14*QUARTIC(p)/t;
            double dphiPdt = (-0.110/SQUARE(t) + 1.405e-7 - 2.0*2.394e-11*t)*p
            + (-7.337e-8/SQUARE(t) - 9.502e-13)*p*p
            - 1.876e-10*CUBE(p)/SQUARE(t) + 1.191e-14*QUARTIC(p)/SQUARE(t);
            double d2phiPdt2 = (2.0*0.110/CUBE(t) - 2.0*2.394e-11)*p
            + 2.0*7.337e-8*p*p/CUBE(t) + 2.0*1.876e-10*CUBE(p)/CUBE(t)
            - 2.0*1.191e-14*QUARTIC(p)/CUBE(t);
            double d3phiPdt3 = -6.0*0.110*p/QUARTIC(t)
            - 6.0*7.337e-8*p*p/QUARTIC(t) - 6.0*1.876e-10*CUBE(p)/QUARTIC(t)
            + 6.0*1.191e-14*QUARTIC(p)/QUARTIC(t);
            double gRobie = r*t*(2.9147*log(t) - 9.6863e-4*t + 6.8593e-8*t*t
                                 + 77.8899/sqrt(t) - 28954.8/t - 2263.27/SQUARE(t) - 15.8997);
            double dgRobiedt = r*(2.9147*log(t) - 9.6863e-4*t + 6.8593e-8*t*t
                                  + 77.8899/sqrt(t) - 28954.8/t - 2263.27/SQUARE(t) - 15.8997) +
            r*t*(2.9147/t - 9.6863e-4 + 2.0*6.8593e-8*t - 0.5*77.8899/pow(t,
                                                                          (double) 1.5) + 28954.8/SQUARE(t) + 2.0*2263.27/CUBE(t));
            double d2gRobiedt2 = 2.0*r*(2.9147/t - 9.6863e-4 + 2.0*6.8593e-8*t
                                        - 0.5*77.8899/pow(t, (double) 1.5) + 28954.8/SQUARE(t)
                                        + 2.0*2263.27/CUBE(t)) + r*t*(-2.9147/SQUARE(t) + 2.0*6.8593e-8
                                                                      + 1.5*0.5*77.8899/pow(t, (double) 2.5) - 2.0*28954.8/CUBE(t)
                                                                      - 6.0*2263.27/QUARTIC(t));
            double d3gRobiedt3 = 3.0*r*(-2.9147/SQUARE(t) + 2.0*6.8593e-8
                                        + 1.5*0.5*77.8899/pow(t, (double) 2.5) - 2.0*28954.8/CUBE(t)
                                        - 6.0*2263.27/QUARTIC(t)) + r*t*(2.0*2.9147/CUBE(t)
                                                                         - 2.5*1.5*0.5*77.8899/pow(t, (double) 3.5) + 6.0*28954.8/QUARTIC(t)
                                                                         + 24.0*2263.27/QUINTIC(t));
            double gH2O, hH2O, sH2O, cpH2O, dcpdtH2O, vH2O, dvdtH2O, dvdpH2O,
            d2vdt2H2O, d2vdtdpH2O, d2vdp2H2O, dgdt, d2gdt2, d3gdt3;

            whaar(1.0, t, &gH2O, &hH2O, &sH2O, &cpH2O, &dcpdtH2O, &vH2O, &dvdtH2O,
                  &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O);

            gl     = r*t*(a/t + b + phiP) + gH2O - gRobie;
            dgdt   = r*(a/t + b + phiP) + r*t*(-a/SQUARE(t) + dphiPdt) - dgRobiedt;
            sl     = sH2O - dgdt;
            hl     = gl + t*sl;
            d2gdt2 = 2.0*r*(-a/SQUARE(t) + dphiPdt) + r*t*(2.0*a/CUBE(t) + d2phiPdt2) - d2gRobiedt2;
            d3gdt3 = 3.0*r*(2.0*a/CUBE(t) + d2phiPdt2) + r*t*(-6.0*a/QUARTIC(t) + d3phiPdt3) - d3gRobiedt3;

            cpl    = cpH2O - t*d2gdt2;
            dcpldt = dcpdtH2O - d2gdt2 - t*d3gdt3;

            vl       = r*(0.110 + 4.432e-5*t + 1.405e-7*t*t - 2.394e-11*CUBE(t))
            + 2.0*r*(7.337e-8 - 1.170e-8*t - 9.502e-13*t*t)*p
            + 3.0*r*(1.876e-10 + 4.586e-13*t)*p*p
            - 4.0*r*1.191e-14*CUBE(p);
            dvldt    = r*(4.432e-5 + 2.0*1.405e-7*t - 3.0*2.394e-11*t*t)
            - 2.0*r*(1.170e-8 + 2.0*9.502e-13*t)*p + 3.0*r*4.586e-13*p*p;
            dvldp    = 2.0*r*(7.337e-8 - 1.170e-8*t - 9.502e-13*t*t)
            + 6.0*r*(1.876e-10 + 4.586e-13*t)*p - 12.0*r*1.191e-14*p*p;
            d2vldt2  = r*(2.0*1.405e-7 - 6.0*2.394e-11*t) - 4.0*r*9.502e-13*p;
            d2vldp2  = 6.0*r*(1.876e-10 + 4.586e-13*t) - 24.0*r*1.191e-14*p;
            d2vldtdp = - 2.0*r*(1.170e-8 + 2.0*9.502e-13*t) + 6.0*r*4.586e-13*p;

            /* special case - no EOS option */
        } else if( (strcmp(name, "H2O") == 0) && (
#ifndef TRUE_xMELTS
						  (calculationMode == MODE_xMELTS) ||
#endif
						  (calculationMode == MODE__MELTSandCO2_H2O) )) {
            double a = -33676.0 + phase->h/r, b = 18.3527 - phase->s/r;
            double gRobie  = r*t*(2.9147*log(t) - 9.6863e-4*t + 6.8593e-8*t*t + 77.8899/sqrt(t) - 28954.8/t - 2263.27/SQUARE(t) - 15.8997);
            double dgRobiedt = r*(2.9147*log(t) - 9.6863e-4*t + 6.8593e-8*t*t + 77.8899/sqrt(t) - 28954.8/t - 2263.27/SQUARE(t) - 15.8997)
                             + r*t*(2.9147/t - 9.6863e-4 + 2.0*6.8593e-8*t - 0.5*77.8899/pow(t, (double) 1.5) + 28954.8/SQUARE(t) + 2.0*2263.27/CUBE(t));
            double d2gRobiedt2 = 2.0*r*(2.9147/t - 9.6863e-4 + 2.0*6.8593e-8*t
                                        - 0.5*77.8899/pow(t, (double) 1.5) + 28954.8/SQUARE(t)
                                        + 2.0*2263.27/CUBE(t)) + r*t*(-2.9147/SQUARE(t) + 2.0*6.8593e-8
                                                                      + 1.5*0.5*77.8899/pow(t, (double) 2.5) - 2.0*28954.8/CUBE(t)
                                                                      - 6.0*2263.27/QUARTIC(t));
            double d3gRobiedt3 = 3.0*r*(-2.9147/SQUARE(t) + 2.0*6.8593e-8
                                        + 1.5*0.5*77.8899/pow(t, (double) 2.5) - 2.0*28954.8/CUBE(t)
                                        - 6.0*2263.27/QUARTIC(t)) + r*t*(2.0*2.9147/CUBE(t)
                                                                         - 2.5*1.5*0.5*77.8899/pow(t, (double) 3.5) + 6.0*28954.8/QUARTIC(t)
                                                                         + 24.0*2263.27/QUINTIC(t));
            double gH2O, hH2O, sH2O, cpH2O, dcpdtH2O, vH2O, dvdtH2O, dvdpH2O, d2vdt2H2O, d2vdtdpH2O, d2vdp2H2O, dgdt, d2gdt2, d3gdt3;

            double vOaksLange    = 2.775;
            double dvdtOaksLange = 1.086e-3;
            double dvdpOaksLange = -0.382e-4;

            a +=  2783.6851512128/r;
            b -=     2.3838467967178/r;

            whaar(1.0, t, &gH2O, &hH2O, &sH2O, &cpH2O, &dcpdtH2O, &vH2O, &dvdtH2O, &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O);

            gl     = r*t*(a/t + b) + gH2O - gRobie + vOaksLange*(p-1.0) + dvdtOaksLange*(t-1673.15)*(p-1.0) + 0.5*dvdpOaksLange*(p-1.0)*(p-1.0);
            dgdt   = r*(a/t + b) + r*t*(-a/SQUARE(t)) - dgRobiedt + dvdtOaksLange*(p-1.0);
            sl     = sH2O - dgdt;
            hl     = r*t*(a/t + b) + gH2O - gRobie  + t*(sH2O - dgdt) + vOaksLange*(p-1.0) + dvdtOaksLange*(t-1673.15)*(p-1.0) + 0.5*dvdpOaksLange*(p-1.0)*(p-1.0);
            d2gdt2 = 2.0*r*(-a/SQUARE(t)) + r*t*(2.0*a/CUBE(t)) - d2gRobiedt2;
            d3gdt3 = 3.0*r*(2.0*a/CUBE(t)) + r*t*(-6.0*a/QUARTIC(t)) - d3gRobiedt3;

            cpl    = cpH2O - t*d2gdt2;
            dcpldt = dcpdtH2O - d2gdt2 - t*d3gdt3;

            vl       = vOaksLange + dvdtOaksLange*(t-1673.15) + dvdpOaksLange*(p-1.0);
            dvldt    = dvdtOaksLange;
            dvldp    = dvdpOaksLange;
            d2vldt2  = 0.0;
            d2vldp2  = 0.0;
            d2vldtdp = 0.0;

            /* special case - no EOS option */
        } else if( (strcmp(name, "CO2") == 0) && (
#ifndef TRUE_xMELTS
                                                  (calculationMode == MODE_xMELTS) ||
#endif
                                                  (calculationMode == MODE__MELTSandCO2) ||
                                                  (calculationMode == MODE__MELTSandCO2_H2O) )) {
            double hCO2       =  - 630.93193811701;
            double sCO2       =  - 109.39331414050;
            double vCO2       =  4.0157994267547;
            double dvCO2dt    =  1.213189e-3;
            double dvCO2dp    = -0.4267387e-4;

            double gDuan1bar, hDuan1bar, sDuan1bar, cpDuan1bar, dcpdtDuan1bar, vDuan1bar, dvdtDuan1bar, dvdpDuan1bar, d2vdt2Duan1bar, d2vdtdpDuan1bar, d2vdp2Duan1bar;
            propertiesOfPureCO2(t, 1.0, &gDuan1bar, &hDuan1bar, &sDuan1bar, &cpDuan1bar, &dcpdtDuan1bar,
                                &vDuan1bar, &dvdtDuan1bar, &dvdpDuan1bar, &d2vdt2Duan1bar, &d2vdtdpDuan1bar, &d2vdp2Duan1bar);

            vl       = vCO2 + dvCO2dt*(t-trl) + dvCO2dp*(p-pr);
            dvldt    = dvCO2dt;
            dvldp    = dvCO2dp;
            d2vldt2  = 0.0;
            d2vldp2  = 0.0;
            d2vldtdp = 0.0;

            gl     = gDuan1bar + hCO2 - t*sCO2 + vCO2*(p-pr) + dvCO2dt*(t-trl)*(p-pr) + dvCO2dp*(p*p/2.0 - pr*pr/2.0 -pr*(p-pr));
            sl     = sDuan1bar + sCO2 - dvCO2dt*(p-pr);
            hl     = hDuan1bar + hCO2 + vCO2*(p-pr) + dvCO2dt*(t-trl)*(p-pr) + dvCO2dp*(p*p/2.0 - pr*pr/2.0 -pr*(p-pr)) - t*dvCO2dt*(p-pr);;
            cpl    = cpDuan1bar;
            dcpldt = dcpdtDuan1bar;

            /* special case - no EOS option */
        } else if( (strcmp(name, "CaCO3") == 0) && ((calculationMode == MODE__MELTSandCO2) ||
                                                    (calculationMode == MODE__MELTSandCO2_H2O) ||
                                                    (calculationMode == MODE_xMELTS) ) ) {
            double hCorr = - 17574.497522747;
            double vCorr = - 1.9034060173857;

            getCaCO3properties(t, p, &gl, &hl, &sl, &cpl, &dcpldt, &vl, &dvldt, &dvldp, &d2vldt2, &d2vldtdp, &d2vldp2);

            hl += hCorr;
            vl += vCorr;
            gl += hCorr + vCorr*(p-pr);

            /* special case - no EOS option */
        } else if ( (strcmp(name, "Si0.25OH") == 0) && (calculationMode == MODE_xMELTS) ) {
            hl     = phase->h;
            sl     = phase->s;
            cpl    = 0.0;
            dcpldt = 0.0;
            gl     = hl - t*sl;
            getSiOHproperties(t, &gl, &hl, &sl, &cpl, &dcpldt);

            /* special case - no EOS option */
        } else if ( (strcmp(name, "Al0.33OH") == 0) && (calculationMode == MODE_xMELTS) ) {
            hl     = phase->h;
            sl     = phase->s;
            cpl    = 0.0;
            dcpldt = 0.0;
            gl     = hl - t*sl;
            getAlOHproperties(t, &gl, &hl, &sl, &cpl, &dcpldt);

            /* special case - no EOS option */
        } else if ( (strcmp(name, "Fe0.33OH") == 0) && (calculationMode == MODE_xMELTS) ) {
            hl     = phase->h;
            sl     = phase->s;
            cpl    = 0.0;
            dcpldt = 0.0;
            gl     = hl - t*sl;
            getFe3OHproperties(t, &gl, &hl, &sl, &cpl, &dcpldt);

            /* special case - no EOS option */
        } else if ( (strcmp(name, "Fe0.5OH") == 0) && (calculationMode == MODE_xMELTS) ) {
            hl     = phase->h;
            sl     = phase->s;
            cpl    = 0.0;
            dcpldt = 0.0;
            gl     = hl - t*sl;
            getFe2OHproperties(t, &gl, &hl, &sl, &cpl, &dcpldt);

            /* special case - no EOS option */
        } else if ( (strcmp(name, "Mg0.5OH") == 0) && (calculationMode == MODE_xMELTS) ) {
            hl     = phase->h;
            sl     = phase->s;
            cpl    = 0.0;
            dcpldt = 0.0;
            gl     = hl - t*sl;
            getMgOHproperties(t, &gl, &hl, &sl, &cpl, &dcpldt);

            /* special case - no EOS option */
        } else if ( (strcmp(name, "Ca0.5OH") == 0) && (calculationMode == MODE_xMELTS) ) {
            hl     = phase->h;
            sl     = phase->s;
            cpl    = 0.0;
            dcpldt = 0.0;
            gl     = hl - t*sl;
            getCaOHproperties(t, &gl, &hl, &sl, &cpl, &dcpldt);

            /* special case - no EOS option */
        } else if ( (strcmp(name, "NaOH") == 0) && (calculationMode == MODE_xMELTS) ) {
            hl     = phase->h;
            sl     = phase->s;
            cpl    = 0.0;
            dcpldt = 0.0;
            gl     = hl - t*sl;
            getNaOHproperties(t, &gl, &hl, &sl, &cpl, &dcpldt);

            /* special case - no EOS option */
        } else if ( (strcmp(name, "KOH") == 0) && (calculationMode == MODE_xMELTS) ) {
            hl     = phase->h;
            sl     = phase->s;
            cpl    = 0.0;
            dcpldt = 0.0;
            gl     = hl - t*sl;
            getKOHproperties(t, &gl, &hl, &sl, &cpl, &dcpldt);

            /* special case - no EOS option */
        } else if ( (strcmp(name, "Fe2SiO5") == 0) && (calculationMode == MODE_xMELTS) ) {
            double hO2, sO2, cpO2, dcpdtO2;
            /* tabulated properties of fayalite solid */
            ThermoRef  tempRef = {
                -1479360.0, 150.930, 0.0,
                CP_BERMAN,  {{248.93, -19.239E2, 0.0, -13.910E7, 0.0, 0.0, 0.0, 0.0}},
                EOS_BERMAN, {{0.0, 0.0, 0.0, 0.0}}
            };
            /* Fusion properties of fayalite */
            static const double tFus = 1490.0; /* K */
            static const double sFus = 59.9;   /* J/K-mol */
            /* Heat capacity of fayalite liquid (after LN) */
            static const double cpLiq = 82.6 + 2.0*78.8;
            /* Speciation properties from Kress and Carmichael (1991) */
            static const double deltaH  = -106000.0; /* J/mol   */
            static const double deltaS  = -55.1;     /* J/K-mol */
            static const double deltaCp = 31.86;     /* J/K-mol */
            static const double T0      = 1673.15;   /* K       */
            /* Execute first time through - thereafter fusion->g is set */
            if (fusion->g == 0.0) {
                gibbs(tFus, pr, "dummy", &tempRef, NULL, NULL, fusion);
                fusion->h += tFus*sFus;
                fusion->s += sFus;
            }
            /* These are now the properties of fayalite liquid at t, pr */
            hl  = fusion->h + cpLiq*(t - tFus);
            sl  = fusion->s + cpLiq*log(t/tFus);
            cpl = cpLiq;

            /* These are the properties of O2 gas at pr */
            hO2 = 23.10248*(t-tr) + 2.0*804.8876*(sqrt(t)-sqrt(tr)) - 1762835.0*(1.0/t-1.0/tr)
            - 18172.91960*log(t/tr) + 0.5*0.002676*(t*t-tr*tr);
            sO2 = 205.15 + 23.10248*log(t/tr) - 2.0*804.8876*(1.0/sqrt(t)-1.0/sqrt(tr))
            - 0.5*1762835.0*(1.0/(t*t)-1.0/(tr*tr)) + 18172.91960*(1.0/t-1.0/tr) + 0.002676*(t-tr);
            cpO2 = 23.10248 + 804.8876/sqrt(t) + 1762835.0/SQUARE(t) - 18172.91960/t + 0.002676*t;
            dcpdtO2 = - 0.5*804.8876/pow(t, (double) 1.5) - 2.0*1762835.0/CUBE(t) + 18172.91960/SQUARE(t) + 0.002676;

            hl     += 2.0*(deltaH + deltaCp*(t-T0))    + hO2/2.0;
            sl     += 2.0*(deltaS + deltaCp*log(t/T0)) + sO2/2.0;
            cpl    += 2.0*deltaCp                      + cpO2/2.0;
            dcpldt  = dcpdtO2/2.0;

            /* special case - no EOS option */
        } else if ( (strcmp(name, "Fe2SiO4.6") == 0) && (calculationMode == MODE_xMELTS) ) {
            double hO2, sO2, cpO2, dcpdtO2;
            /* tabulated properties of fayalite solid */
            ThermoRef  tempRef = {
                -1479360.0, 150.930, 0.0,
                CP_BERMAN,  {{248.93, -19.239E2, 0.0, -13.910E7, 0.0, 0.0, 0.0, 0.0}},
                EOS_BERMAN, {{0.0, 0.0, 0.0, 0.0}}
            };
            /* Fusion properties of fayalite */
            static const double tFus = 1490.0; /* K */
            static const double sFus = 59.9;   /* J/K-mol */
            /* Heat capacity of fayalite liquid (after LN) */
            static const double cpLiq = 82.6 + 2.0*78.8;
            /* Speciation properties from Kress and Carmichael (1991) */
            static const double deltaH  = -106000.0; /* J/mol   */
            static const double deltaS  = -55.1;     /* J/K-mol */
            static const double deltaCp = 31.86;     /* J/K-mol */
            static const double K2      = 0.4;
            static const double T0      = 1673.15;   /* K       */
            static const double y       = 0.3;
            /* Execute first time through - thereafter fusion->g is set */
            if (fusion->g == 0.0) {
                gibbs(tFus, pr, "dummy", &tempRef, NULL, NULL, fusion);
                fusion->h += tFus*sFus;
                fusion->s += sFus;
            }
            /* These are now the properties of fayalite liquid at t, pr */
            hl  = fusion->h + cpLiq*(t - tFus);
            sl  = fusion->s + cpLiq*log(t/tFus);
            cpl = cpLiq;

            /* These are the properties of O2 gas at pr */
            hO2 = 23.10248*(t-tr) + 2.0*804.8876*(sqrt(t)-sqrt(tr)) - 1762835.0*(1.0/t-1.0/tr)
            - 18172.91960*log(t/tr) + 0.5*0.002676*(t*t-tr*tr);
            sO2 = 205.15 + 23.10248*log(t/tr) - 2.0*804.8876*(1.0/sqrt(t)-1.0/sqrt(tr))
            - 0.5*1762835.0*(1.0/(t*t)-1.0/(tr*tr)) + 18172.91960*(1.0/t-1.0/tr) + 0.002676*(t-tr);
            cpO2 = 23.10248 + 804.8876/sqrt(t) + 1762835.0/SQUARE(t) - 18172.91960/t + 0.002676*t;
            dcpdtO2 = - 0.5*804.8876/pow(t, (double) 1.5) - 2.0*1762835.0/CUBE(t) + 18172.91960/SQUARE(t) + 0.002676;

            hl     += 4.0*y*(deltaH + deltaCp*(t-T0))    + y*hO2;
            sl     += 4.0*y*(deltaS + deltaCp*log(t/T0)) + y*sO2 - 2.0*r*log(K2);
            cpl    += 4.0*y*deltaCp                      + y*cpO2;
            dcpldt  = y*dcpdtO2;


            /* special case - no EOS option */
        } else if ( (strcmp(name, "Fe2AlO4.5") == 0) && (calculationMode == MODE_xMELTS) ) {
            hl     = phase->h;
            sl     = phase->s;
            cpl    = 0.0;
            dcpldt = 0.0;
            gl     = hl - t*sl;
            getFe2AlO4_5properties(t, &gl, &hl, &sl, &cpl, &dcpldt);

            /* special case - no EOS option */
        } else if ( (strcmp(name, "Fe2AlO4.1") == 0) && (calculationMode == MODE_xMELTS) ) {
            hl     = phase->h;
            sl     = phase->s;
            cpl    = 0.0;
            dcpldt = 0.0;
            gl     = hl - t*sl;
            getFe2AlO4_1properties(t, &gl, &hl, &sl, &cpl, &dcpldt);

            /* all other components, MELTS, pMELTS, XMELTS */
        } else {
            if (fusion->g == 0.0 && phase->h != 0.0) {
                gibbs(liquid->tfus, pr, "dummy", phase, NULL, NULL, fusion);
                fusion->h += liquid->tfus * liquid->sfus;
                fusion->s += liquid->sfus;
            }
            if (t > liquid->tglass) {
                hl = fusion->h + liquid->cp*(t - liquid->tfus);
                sl = fusion->s + liquid->cp*log(t/liquid->tfus);
                cpl = liquid->cp;
                dcpldt = 0.0;
            } else {
                hl = fusion->h + liquid->cp*(liquid->tglass - liquid->tfus);
                sl = fusion->s + liquid->cp*log(liquid->tglass/liquid->tfus);

                k0 =   phase->cp.Berman.k0;
                k1 =   phase->cp.Berman.k1;
                k2 =   phase->cp.Berman.k2;
                k3 =   phase->cp.Berman.k3;
                cp_t = phase->cp.Berman.Tt;
                cp_h = phase->cp.Berman.deltah;
                l1 =   phase->cp.Berman.l1;
                l2 =   phase->cp.Berman.l2;
                k0 += liquid->cp - (k0 + k1/sqrt(liquid->tglass)
                                    + k2/SQUARE(liquid->tglass) + k3/CUBE(liquid->tglass));
                hl = hl + k0*(t-liquid->tglass)
                + 2.0*k1*(sqrt(t)-sqrt(liquid->tglass))
                - k2*(1.0/t-1.0/liquid->tglass)
                - 0.5*k3*(1.0/SQUARE(t) - 1.0/SQUARE(liquid->tglass));
                sl = sl + k0*log(t/liquid->tglass)
                - 2.0*k1*(1.0/sqrt(t) - 1.0/sqrt(liquid->tglass))
                - 0.5*k2*(1.0/SQUARE(t) - 1.0/SQUARE(liquid->tglass))
                - (1.0/3.0)*k3*(1.0/CUBE(t) - 1.0/CUBE(liquid->tglass));
                cpl = k0 + k1/sqrt(t) + k2/(t*t) + k3/CUBE(t);
                dcpldt = - 0.5*k1/pow(t, (double) 1.5) - 2.0*k2/CUBE(t)
                - 3.0*k3/QUARTIC(t);
                if(t < cp_t ) {
                    hl = hl + 0.5*l1*l1*(t*t-cp_t*cp_t)
                    + (2.0/3.0)*l1*l2*(CUBE(t) - CUBE(cp_t))
                    + 0.25*l2*l2*(QUARTIC(t) - QUARTIC(cp_t));
                    sl = sl + l1*l1*(t - cp_t) + l1*l2*(t*t - cp_t*cp_t)
                    + (1.0/3.0)*l2*l2*(CUBE(t) - CUBE(cp_t));
                    cpl += t*SQUARE(l1+l2*t);
                    dcpldt += SQUARE(l1+l2*t) + t*2.0*(l1+l2*t)*l2;
                }
            }

            /* MELTS case - always assume Kress (Lange/polynomial) EOS                  */
            /* Placed here because SiO2, H2O and CO2 special cases are dealt with above */
            if ((calculationMode == MODE__MELTS) ||
                (calculationMode == MODE__MELTSandCO2) ||
                (calculationMode == MODE__MELTSandCO2_H2O) ||
                ((calculationMode == MODE_xMELTS) && (liquid->eos_type != EOS_GHIORSO) )) {
                gl = hl - t*sl
                + (liquid->v + liquid->eos.Kress.dvdt*(t-trl))*(p-pr)
                + 0.5*(liquid->eos.Kress.dvdp + (t-trl)*liquid->eos.Kress.d2vdtp)*(p*p-pr*pr)
                - (liquid->eos.Kress.dvdp + (t-trl)*liquid->eos.Kress.d2vdtp)*pr*(p-pr)
                + liquid->eos.Kress.d2vdp2*( (CUBE(p)-CUBE(pr))/6.0 - pr*(p*p-pr*pr)/2.0 + pr*pr*(p-pr)/2.0 );
                hl += (liquid->v + liquid->eos.Kress.dvdt*(t-trl))*(p-pr)
                + 0.5*(liquid->eos.Kress.dvdp + (t-trl)*liquid->eos.Kress.d2vdtp)*(p*p-pr*pr)
                - (liquid->eos.Kress.dvdp + (t-trl)*liquid->eos.Kress.d2vdtp)*pr*(p-pr)
                + liquid->eos.Kress.d2vdp2*( (CUBE(p)-CUBE(pr))/6.0 - pr*(p*p-pr*pr)/2.0 + pr*pr*(p-pr)/2.0 )
                - t*(liquid->eos.Kress.dvdt*(p-pr) + 0.5*liquid->eos.Kress.d2vdtp*(p-pr)*(p-pr));
                sl += -(liquid->eos.Kress.dvdt*(p-pr) + 0.5*liquid->eos.Kress.d2vdtp*(p-pr)*(p-pr));

                vl       = liquid->v + liquid->eos.Kress.dvdt*(t-trl)
                + (liquid->eos.Kress.dvdp + liquid->eos.Kress.d2vdtp*(t-trl))*(p-pr)
                + liquid->eos.Kress.d2vdp2*(0.5*p*p - pr*(p-pr));
                dvldt    = liquid->eos.Kress.dvdt + liquid->eos.Kress.d2vdtp*(p-pr);
                dvldp    = liquid->eos.Kress.dvdp + liquid->eos.Kress.d2vdtp*(t-trl) + liquid->eos.Kress.d2vdp2*(p-pr);
                d2vldt2  = 0.0;
                d2vldp2  = liquid->eos.Kress.d2vdp2;
                d2vldtdp = liquid->eos.Kress.d2vdtp;
            }
        }

        /* xMELTS case - option for Ghiorso EOS */
        if (calculationMode == MODE_xMELTS && liquid->eos_type == EOS_GHIORSO) {
            double alpha     = liquid->eos.Ghiorso.alpha;
            double v0        = liquid->v*exp(alpha*(t-trl));
            double cRef      = liquid->eos.Ghiorso.c;
            double dcdt      = liquid->eos.Ghiorso.dcdt;
            double c         = cRef + (t-trl)*dcdt;
            double v2        = liquid->eos.Ghiorso.d2vdp2;
            double v3        = liquid->eos.Ghiorso.d3vdp3;
            double v4        = liquid->eos.Ghiorso.d4vdp4;
            double mw        = liquid->eos.Ghiorso.mw;
            double cp        = liquid->eos.Ghiorso.cp;
            double gInt      = 0.0;
            double dgIntdt   = 0.0;
            double d2gIntdt2 = 0.0;
            double d3gIntdt3 = 0.0;
            double v1Ref, v1, dv1dt, d2v1dt2, d3v1dt3, a, b, sum;

            v1Ref = -(liquid->v)*(liquid->v)*(1000.0/(mw*cRef*cRef) + trl*alpha*alpha/(10.0*cp));
            v1    = -v0*v0*(1000.0/(mw*c*c)       +   t*alpha*alpha/(10.0*cp));
            a     = ((2.0*v1Ref*v3 - 3.0*v2*v2) != 0.0) ? (v2*v3 - v1Ref*v4/2.0)/(2.0*v1Ref*v3 - 3.0*v2*v2)  : 0.0;
            b     = ((2.0*v1Ref*v3 - 3.0*v2*v2) != 0.0) ? (v2*v4/4.0 - v3*v3/3.0)/(2.0*v1Ref*v3 - 3.0*v2*v2) : 0.0;
            sum   = a*a - 4.0*b;

            dv1dt   = -2.0*v0*v0*alpha*(1000.0/(mw*c*c) + t*alpha*alpha/(10.0*cp))
            - v0*v0*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(10.0*cp));
            d2v1dt2 = -4.0*v0*v0*alpha*alpha*(1000.0/(mw*c*c) + t*alpha*alpha/(10.0*cp))
            - 4.0*v0*v0*alpha*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(10.0*cp))
            - v0*v0*(6000.0*dcdt*dcdt/(mw*c*c*c*c));
            d3v1dt3 = -8.0*v0*v0*alpha*alpha*alpha*(1000.0/(mw*c*c) + t*alpha*alpha/(10.0*cp))
            - 12.0*v0*v0*alpha*alpha*(-2000.0*dcdt/(mw*c*c*c) + alpha*alpha/(10.0*cp))
            - 6.0*v0*v0*alpha*(6000.0*dcdt*dcdt/(mw*c*c*c*c)) - v0*v0*(-24000.0*dcdt*dcdt*dcdt/(mw*c*c*c*c*c));

            vl	  = (v0 + (v1+v0*a)*(p-pr) + (v2/2.0+v1*a+v0*b)*(p-pr)*(p-pr))/(1.0 + a*(p-pr) + b*(p-pr)*(p-pr));
            dvldt    = (alpha*v0 + (dv1dt+alpha*v0*a)*(p-pr) + (dv1dt*a+alpha*v0*b)*(p-pr)*(p-pr))/(1.0 + a*(p-pr) + b*(p-pr)*(p-pr));
            dvldp    = (v1 + v0*a + 2.0*(v2/2.0+v1*a+v0*b)*(p-pr) - vl*(a + 2.0*b*(p-pr)))/(1.0 + a*(p-pr) + b*(p-pr)*(p-pr));
            d2vldt2  = (alpha*alpha*v0 + (d2v1dt2+alpha*alpha*v0*a)*(p-pr) + (d2v1dt2*a+alpha*alpha*v0*b)*(p-pr)*(p-pr))
            /(1.0 + a*(p-pr) + b*(p-pr)*(p-pr));
            d2vldp2  = (v2 + 2.0*v1*a + 2.0*v0*b - 2.0*dvldp*(a + 2.0*b*(p-pr)) - 2.0*vl*b)/(1.0 + a*(p-pr) + b*(p-pr)*(p-pr));
            d2vldtdp = (dv1dt + alpha*v0*a + 2.0*(dv1dt*a+alpha*v0*b)*(p-pr) - dvldt*(a + 2.0*b*(p-pr)))/(1.0 + a*(p-pr) + b*(p-pr)*(p-pr));

            if ((a == 0.0) && (b == 0.0)) {
                gInt      += v0*(p-pr) + v1*(p-pr)*(p-pr)/2.0 + v2*(p-pr)*(p-pr)*(p-pr)/6.0;
                dgIntdt   += alpha*v0*(p-pr) + dv1dt*(p-pr)*(p-pr)/2.0;
                d2gIntdt2 += alpha*alpha*v0*(p-pr) + d2v1dt2*(p-pr)*(p-pr)/2.0;
                d3gIntdt3 += alpha*alpha*alpha*v0*(p-pr) + d3v1dt3*(p-pr)*(p-pr)/2.0;

            } else if ((a != 0.0) && (b == 0.0)) {
                gInt      += (v0 - v2/(2.0*a*a))*(p-pr) + (v1 + v2/(2.0*a))*(p-pr)*(p-pr)/2.0 + v2*log(1.0+a*(p-pr))/(2.0*a*a*a);
                dgIntdt   += alpha*v0*(p-pr) + dv1dt*(p-pr)*(p-pr)/2.0;
                d2gIntdt2 += alpha*alpha*v0*(p-pr) + d2v1dt2*(p-pr)*(p-pr)/2.0;
                d3gIntdt3 += alpha*alpha*alpha*v0*(p-pr) + d3v1dt3*(p-pr)*(p-pr)/2.0;

            } else if ((a == 0.0) && (b != 0.0)) {
                gInt      += (v0 + v2/(2.0*b))*(p-pr) + v1*log(1.0 + b*(p-pr)*(p-pr))/(2.0*b);
                gInt      += (b > 0.0) ? -v2*atan(sqrt(b)*(p-pr))/(2.0*b*sqrt(b)) : -v2*log((1.0+sqrt(-b)*(p-pr))/(1.0-sqrt(-b)*(p-pr)))/(4.0*b*sqrt(-b));
                dgIntdt   += alpha*v0*(p-pr) + dv1dt*log(1.0 + b*(p-pr)*(p-pr))/(2.0*b);
                d2gIntdt2 += alpha*alpha*v0*(p-pr) + d2v1dt2*log(1.0 + b*(p-pr)*(p-pr))/(2.0*b);
                d3gIntdt3 += alpha*alpha*alpha*v0*(p-pr) + d3v1dt3*log(1.0 + b*(p-pr)*(p-pr))/(2.0*b);

            } else if (sum > 0.0) {
                double x = sqrt(sum);
                double y = (a + 2.0*b*(p-pr))/x;
                double z = a/x;
                double PcA = (2.0*b*pr - a + x)/(2.0*b);
                double PcB = (2.0*b*pr - a - x)/(2.0*b);
                double arg = (1.0 + y)*(1.0 - z)/((1.0 - y)*(1.0 + z));

                gInt      += (v0 + a*v1/b + v2/(2.0*b))*(p-pr);
                gInt      += (((pr < PcA) && (PcA < p)) || ((pr < PcB) && (PcB < p))) ? (v1*(1.0-a*a/b) - v2*a/(2.0*b))*(-sqrt(DBL_MAX))/(2.0*b)
                : (v1*(1.0-a*a/b) - v2*a/(2.0*b))*log(1.0+a*(p-pr)+b*(p-pr)*(p-pr))/(2.0*b);
                gInt      += (arg > 0) ? (v1*a*(3.0-a*a/b) + v2*(1.0 - a*a/(2.0*b)))*log(arg)/(2.0*b*x) : (v1*a*(3.0-a*a/b) + v2*(1.0 - a*a/(2.0*b)))*(-sqrt(DBL_MAX))/(2.0*b*x);

                dgIntdt   += (alpha*v0 + a*dv1dt/b)*(p-pr);
                dgIntdt   += (((pr < PcA) && (PcA < p)) || ((pr < PcB) && (PcB < p))) ? dv1dt*(1.0-a*a/b)*(-sqrt(DBL_MAX))/(2.0*b)
                : dv1dt*(1.0-a*a/b)*log(1.0+a*(p-pr)+b*(p-pr)*(p-pr))/(2.0*b);
                dgIntdt   += (arg > 0) ? dv1dt*a*(3.0-a*a/b)*log(arg)/(2.0*b*x) : dv1dt*a*(3.0-a*a/b)*(-sqrt(DBL_MAX))/(2.0*b*x);

                d2gIntdt2 += (alpha*alpha*v0 + a*d2v1dt2/b)*(p-pr);
                d2gIntdt2 += (((pr < PcA) && (PcA < p)) || ((pr < PcB) && (PcB < p))) ? d2v1dt2*(1.0-a*a/b)*(-sqrt(DBL_MAX))/(2.0*b)
                : d2v1dt2*(1.0-a*a/b)*log(1.0+a*(p-pr)+b*(p-pr)*(p-pr))/(2.0*b);
                d2gIntdt2 += (arg > 0) ? d2v1dt2*a*(3.0-a*a/b)*log(arg)/(2.0*b*x) : d2v1dt2*a*(3.0-a*a/b)*(-sqrt(DBL_MAX))/(2.0*b*x);

                d3gIntdt3 += (alpha*alpha*alpha*v0 + a*d3v1dt3/b)*(p-pr);
                d3gIntdt3 += (((pr < PcA) && (PcA < p)) || ((pr < PcB) && (PcB < p))) ? d3v1dt3*(1.0-a*a/b)*(-sqrt(DBL_MAX))/(2.0*b)
                : d3v1dt3*(1.0-a*a/b)*log(1.0+a*(p-pr)+b*(p-pr)*(p-pr))/(2.0*b);
                d3gIntdt3 += (arg > 0) ? d3v1dt3*a*(3.0-a*a/b)*log(arg)/(2.0*b*x) : d3v1dt3*a*(3.0-a*a/b)*(-sqrt(DBL_MAX))/(2.0*b*x);

            } else if (sum == 0.0) {
                gInt      += (v0 + 4.0*v1/a + 2.0*v2/(a*a))*(p-pr);
                gInt      += -8.0*(v2/a + v1)/(a*a*(2.0+a*(p-pr))) + 4.0*(v2/a + v1)/(a*a);
                gInt      += -4.0*(3.0*v1 + 2.0*v2/a)*log(1.0 + a*(p-pr)/2.0)/(a*a);

                dgIntdt   += (alpha*v0 + 4.0*dv1dt/a)*(p-pr);
                dgIntdt   += -8.0*dv1dt/(a*a*(2.0+a*(p-pr))) + 4.0*dv1dt/(a*a);
                dgIntdt   += -4.0*3.0*dv1dt*log(1.0 + a*(p-pr)/2.0)/(a*a);

                d2gIntdt2 += (alpha*alpha*v0 + 4.0*d2v1dt2/a)*(p-pr);
                d2gIntdt2 += -8.0*d2v1dt2/(a*a*(2.0+a*(p-pr))) + 4.0*d2v1dt2/(a*a);
                d2gIntdt2 += -4.0*3.0*d2v1dt2*log(1.0 + a*(p-pr)/2.0)/(a*a);

                d3gIntdt3 += (alpha*alpha*alpha*v0 + 4.0*d3v1dt3/a)*(p-pr);
                d3gIntdt3 += -8.0*d3v1dt3/(a*a*(2.0+a*(p-pr))) + 4.0*d3v1dt3/(a*a);
                d3gIntdt3 += -4.0*3.0*d3v1dt3*log(1.0 + a*(p-pr)/2.0)/(a*a);

            } else if(sum < 0.0) {
                double x = sqrt(-sum);
                double y = (a + 2.0*b*(p-pr))/x;
                double z = a/x;

                gInt      += (v0 + a*v1/b + v2/(2.0*b))*(p-pr);
                gInt      += (v1*(1.0-a*a/b) - v2*a/(2.0*b))*log(1.0+a*(p-pr)+b*(p-pr)*(p-pr))/(2.0*b);
                gInt      += (v1*a*(3.0-a*a/b) + v2*(1.0 - a*a/(2.0*b)))*atan((z-y)/(1.0+z*y))/(b*x);

                dgIntdt   += (alpha*v0 + a*dv1dt/b)*(p-pr);
                dgIntdt   += dv1dt*(1.0-a*a/b)*log(1.0+a*(p-pr)+b*(p-pr)*(p-pr))/(2.0*b);
                dgIntdt   += dv1dt*a*(3.0-a*a/b)*atan((z-y)/(1.0+z*y))/(b*x);

                d2gIntdt2 += (alpha*alpha*v0 + a*d2v1dt2/b)*(p-pr);
                d2gIntdt2 += d2v1dt2*(1.0-a*a/b)*log(1.0+a*(p-pr)+b*(p-pr)*(p-pr))/(2.0*b);
                d2gIntdt2 += d2v1dt2*a*(3.0-a*a/b)*atan((z-y)/(1.0+z*y))/(b*x);

                d3gIntdt3 += (alpha*alpha*alpha*v0 + a*d3v1dt3/b)*(p-pr);
                d3gIntdt3 += d3v1dt3*(1.0-a*a/b)*log(1.0+a*(p-pr)+b*(p-pr)*(p-pr))/(2.0*b);
                d3gIntdt3 += d3v1dt3*a*(3.0-a*a/b)*atan((z-y)/(1.0+z*y))/(b*x);

            }

            gl      = hl - t*sl + gInt;
            sl     += -dgIntdt;
            hl     += gInt - t*dgIntdt;
            cpl    += -t*d2gIntdt2;
            dcpldt += -d2gIntdt2 -t*d3gIntdt3;
        }

        /* pMELTS case; xMELTS case - option for Birch-Murnaghan EOS */

#ifdef TRUE_xMELTS
        if ((calculationMode == MODE_pMELTS) || (calculationMode == MODE_xMELTS && liquid->eos_type == EOS_KRESS)) {
#else
        if (calculationMode == MODE_pMELTS) {
#endif
            if (p > pr && (liquid->eos.Kress.dvdp + liquid->eos.Kress.d2vdtp*(t-trl)) != 0.0) {
                double d2v0dtdp   = liquid->eos.Kress.d2vdtp;
                double v0	      = liquid->v + liquid->eos.Kress.dvdt*(t-trl);
                double dv0dt      = liquid->eos.Kress.dvdt;
                double dv0dp      = liquid->eos.Kress.dvdp + d2v0dtdp*(t-trl);
                double K	      = -v0/dv0dp;
                double dKdt       =      d2v0dtdp*v0/(dv0dp*dv0dp)                               -     dv0dt/dv0dp;
                double d2Kdt2     = -2.0*d2v0dtdp*d2v0dtdp*v0/(dv0dp*dv0dp*dv0dp)                + 2.0*d2v0dtdp*dv0dt/(dv0dp*dv0dp);
                double d3Kdt3     =  6.0*d2v0dtdp*d2v0dtdp*d2v0dtdp*v0/(dv0dp*dv0dp*dv0dp*dv0dp) - 6.0*d2v0dtdp*d2v0dtdp*dv0dt/(dv0dp*dv0dp*dv0dp);
                double Kp	      = (calculationMode == MODE_pMELTS) ? 5.0 : liquid->eos.Kress.d2vdp2;
                int iter          = 0;
                const int iterMax = 1000;
                double v, vLast, fn , dfn;
                double f, dfdt, d2fdt2, d3fdt3, A, B, C, D, E, F, dDdt, dEdt, dFdt, d3vldt3, d2gdt2, d3gdt3;
                double minusIntPdV, d_minusIntPdV_df, d2_minusIntPdV_df2, d3_minusIntPdV_df3;

                /* Newton's method for volume from BM */
                vLast  = v0;
                v      = v0*0.99;
                while ((fabs(v - vLast) > 10.0*DBL_EPSILON) && (iter < iterMax)) {
                    fn      = (3.0/2.0)*K*(pow(v0/v, (double) 7.0/3.0) - pow(v0/v, (double) 5.0/3.0))
                    *(1.0-(3.0/4.0)*(4.0-Kp)*(pow(v0/v, (double) 2.0/3.0) - 1.0))
                    - p;
                    dfn     = (3.0/2.0)*K*((7.0/3.0)*pow(v0/v, (double) 4.0/3.0) - (5.0/3.0)*pow(v0/v, (double) 2.0/3.0))*(-v0/(v*v))
                    *(1.0-(3.0/4.0)*(4.0-Kp)*(pow(v0/v, (double) 2.0/3.0) - 1.0))
                    + (3.0/2.0)*K*(pow(v0/v, (double) 7.0/3.0) - pow(v0/v, (double) 5.0/3.0))
                    *(-(3.0/4.0)*(4.0-Kp)*(2.0/3.0)*pow(v0/v, (double) -1.0/3.0))*(-v0/(v*v));
                    vLast   = v;
                    v      += -fn/dfn;
                    if (v > v0     ) v = v0;
                    if (v < 0.01*v0) v = 0.01*v0;
                    iter++;
                }
                if (iter >= iterMax) {
#ifndef ALPHAMELTS_UPDATE_SYSTEM
                    printf("Convergence error in Birch-Murnaghan Volume routine in function gibbs().\n");
                    printf("  For %s,  v = %g, dV = %g, f = %g at T = %g and P = %g in %d iterations.\n", name, v, fabs(v - vLast), fn, t, p, iter);
#else
                    fprintf(stderr, "Convergence error in Birch-Murnaghan Volume routine in function gibbs().\n");
                    fprintf(stderr, "  For %s,  v = %g, dV = %g, f = %g at T = %g and P = %g in %d iterations.\n", name, v, fabs(v - vLast), fn, t, p, iter);
#endif
                }

                f = (pow(v0/v, (double) 2.0/3.0) - 1.0)/2.0;
                A = 1.0/f + 5.0/(1.0+2.0*f) + (3.0/2.0)*(Kp-4.0)/(1.0+(3.0/2.0)*f*(Kp-4.0));
                B = 1.0/(f*f) + 10.0/((1.0+2.0*f)*(1.0+2.0*f)) + pow((3.0/2.0)*(Kp-4.0)/(1.0+(3.0/2.0)*f*(Kp-4.0)), (double) 2.0);
                C = -2.0/(f*f*f) - 40.0/((1.0+2.0*f)*(1.0+2.0*f)*(1.0+2.0*f)) - 2.0*pow((3.0/2.0)*(Kp-4.0)/(1.0+(3.0/2.0)*f*(Kp-4.0)), (double) 3.0);

                D    = dKdt*dv0dt/K;
                dDdt = dv0dt*(d2Kdt2/K - dKdt*dKdt/(K*K));
                E    = (dv0dt*dKdt + v0*d2Kdt2)/K;
                dEdt = 2.0*dv0dt*d2Kdt2/K + v0*d3Kdt3/K - dv0dt*dKdt*dKdt/(K*K) - v0*dKdt*d2Kdt2/(K*K);
                F    = - v0*dKdt*dKdt*(A - 5.0/(1.0+2.0*f) + B/A)/(A*K*K);
                dFdt = (-dv0dt*dKdt*dKdt - 2.0*v0*dKdt*d2Kdt2 + 2.0*v0*dKdt*dKdt*dKdt/K + v0*dKdt*dKdt*dKdt*B/(A*A*K))
                *(A - 5.0/(1.0+2.0*f) + B/A)/(A*K*K)
                - v0*dKdt*dKdt*dKdt*(B - 10.0/((1.0+2.0*f)*(1.0+2.0*f)) + C/A - B*B/(A*A))/(A*A*K*K*K);

                minusIntPdV        = (9.0/2.0)*K*(Kp-4.0)*v0*f*f*f + (9.0/2.0)*K*v0*f*f;
                d_minusIntPdV_df   = (27.0/2.0)*K*(Kp-4.0)*v0*f*f + 9.0*K*v0*f;
                d2_minusIntPdV_df2 = 27.0*K*(Kp-4.0)*v0*f + 9.0*K*v0;
                d3_minusIntPdV_df3 = 27.0*K*(Kp-4.0)*v0;

                dfdt = -dKdt/(A*K);
                d2fdt2 = -d2Kdt2/(A*K) + dKdt*dKdt/(A*K*K) + dKdt*dKdt*B/(A*A*A*K*K);
                d3fdt3 = -d3Kdt3/(A*K) + 3.0*dKdt*d2Kdt2*(B/(A*A)+1.0)/(A*K*K)
                + dKdt*dKdt*dKdt*(C/(A*A*A)-3.0*B/(A*A)-2.0-3.0*B*B/(A*A*A*A))/(A*K*K*K);

                vl       = v;
                dvldt    = 3.0*v0*dKdt/(A*pow(1.0+2.0*f, (double) 5.0/2.0)*K) + dv0dt/pow(1.0+2.0*f, (double) 3.0/2.0);
                dvldp    = -3.0*v0/(p*pow(1.0+2.0*f, (double) 5.0/2.0)*A);
                d2vldt2  = (3.0/(A*pow(1.0+2.0*f, (double) 5.0/2.0)))
                *(dv0dt*dKdt/K + (d2Kdt2*v0+dKdt*dv0dt)/K - v0*dKdt*dKdt*(A-5.0/(1.0+2.0*f)+B/A)/(A*K*K));
                d2vldp2  = (3.0*v0/(p*p*pow(1.0+2.0*f, (double) 5.0/2.0)*A*A))*(A-B/A+5.0/(1.0+2.0*f));
                d2vldtdp = (3.0/(p*pow(1.0+2.0*f, (double) 5.0/2.0)*A))*(-dv0dt + v0*dKdt*(B/A-5.0/(1.0+2.0*f))/(K*A));

                gl      = hl - t*sl + p*v - pr*v0 + minusIntPdV;
                sl     += -(p*dvldt - pr*dv0dt + (minusIntPdV/v0)*dv0dt + (minusIntPdV/K)*dKdt + d_minusIntPdV_df*dfdt);
                hl      = gl + t*sl;

                d2gdt2  = p*d2vldt2 + 2.0*(minusIntPdV/(v0*K))*dKdt*dv0dt + 2.0*(d_minusIntPdV_df/v0)*dfdt*dv0dt
                + 2.0*(d_minusIntPdV_df/K)*dfdt*dKdt + d2_minusIntPdV_df2*dfdt*dfdt
                + (minusIntPdV/K)*d2Kdt2 + d_minusIntPdV_df*d2fdt2;

                cpl    += - t*d2gdt2;

                d3vldt3 = 3.0*(dDdt + dEdt + dFdt)/(A*pow(1.0+2.0*f, (double) 5.0/2.0))
                - 3.0*(D + E + F)*dKdt*(B/A - 5.0/(1.0+2.0*f))/(K*A*A*pow(1.0+2.0*f, (double) 5.0/2.0));
                d3gdt3  = p*d3vldt3 + 3.0*d2_minusIntPdV_df2*dfdt*dfdt*dKdt/K + 6.0*d_minusIntPdV_df*dv0dt*dfdt*dKdt/(K*v0)
                + 3.0*d_minusIntPdV_df*d2fdt2*dKdt/K + d3_minusIntPdV_df3*dfdt*dfdt*dfdt + 3.0*d2_minusIntPdV_df2*dfdt*dfdt*dv0dt/v0
                + 3.0*d_minusIntPdV_df*dfdt*d2Kdt2/K + 3.0*d2_minusIntPdV_df2*dfdt*d2fdt2 + 3.0*minusIntPdV*dv0dt*d2Kdt2/(v0*K)
                + 3.0*d_minusIntPdV_df*dv0dt*d2fdt2/v0 + minusIntPdV*d3Kdt3/K + d_minusIntPdV_df*d3fdt3;


                dcpldt += -t*d3gdt3 - d2gdt2;
            } else {
                vl       = liquid->v;
                dvldt    = 0.0;
                dvldp    = 0.0;
                d2vldt2  = 0.0;
                d2vldp2  = 0.0;
                d2vldtdp = 0.0;

                gl      = hl - t*sl;
                hl      = gl + t*sl;
            }
        }

        result->g       = gl;
        result->h       = hl;
        result->s       = sl;
        result->v       = vl;
        result->cp      = cpl;
        result->dcpdt   = dcpldt;
        result->dvdt    = dvldt;
        result->dvdp    = dvldp;
        result->d2vdt2  = d2vldt2;
        result->d2vdp2  = d2vldp2;
        result->d2vdtdp = d2vldtdp;

        return;
    }

    /* Solids : Special cases handled first! (Oxygen gas and Water are included in this section) */

    if(strcmp(name, "o2") == 0) {
#ifdef SPECIAL_O2
        getO2properties(t, p, &gs, &hs, &ss, &cps, &dcpsdt, &vs, &dvsdt, &dvsdp, &d2vsdt2, &d2vsdp2, &d2vsdtdp);
#elif defined(CONSTANT_P_O2)
        hs = 23.10248*(t-tr) + 2.0*804.8876*(sqrt(t)-sqrt(tr))
        - 1762835.0*(1.0/t-1.0/tr)
        - 18172.91960*log(t/tr) + 0.5*0.002676*(t*t-tr*tr);
        ss = 205.15 + 23.10248*log(t/tr)
        - 2.0*804.8876*(1.0/sqrt(t)-1.0/sqrt(tr))
        - 0.5*1762835.0*(1.0/(t*t)-1.0/(tr*tr))
        + 18172.91960*(1.0/t-1.0/tr) + 0.002676*(t-tr);
        vs = r*t/pr;
        cps = 23.10248 + 804.8876/sqrt(t) + 1762835.0/SQUARE(t)
        - 18172.91960/t + 0.002676*t;
        dcpsdt = - 0.5*804.8876/pow(t, (double) 1.5) - 2.0*1762835.0/CUBE(t)
        + 18172.91960/SQUARE(t) + 0.002676;
        dvsdt = r/pr;
        dvsdp = - r*t/SQUARE(pr);
        d2vsdt2 = 0.0;
        d2vsdp2 = 2.0*r*t/CUBE(pr);
        d2vsdtdp = - r/SQUARE(pr);
        gs = hs - t*ss;         /* standard state is any t and 1 bar for gases */
#elif defined(ZERO_O2)
        hs = 23.10248*(t-tr) + 2.0*804.8876*(sqrt(t)-sqrt(tr))
        - 1762835.0*(1.0/t-1.0/tr)
        - 18172.91960*log(t/tr) + 0.5*0.002676*(t*t-tr*tr);
        ss = 205.15 + 23.10248*log(t/tr)
        - 2.0*804.8876*(1.0/sqrt(t)-1.0/sqrt(tr))
        - 0.5*1762835.0*(1.0/(t*t)-1.0/(tr*tr))
        + 18172.91960*(1.0/t-1.0/tr) + 0.002676*(t-tr);
        vs = r*t/pr;
        cps = 23.10248 + 804.8876/sqrt(t) + 1762835.0/SQUARE(t)
        - 18172.91960/t + 0.002676*t;
        dcpsdt = - 0.5*804.8876/pow(t, (double) 1.5) - 2.0*1762835.0/CUBE(t)
        + 18172.91960/SQUARE(t) + 0.002676;
        dvsdt = r/pr;
        dvsdp = - r*t/SQUARE(pr);
        d2vsdt2 = 0.0;
        d2vsdp2 = 2.0*r*t/CUBE(pr);
        d2vsdtdp = - r/SQUARE(pr);
        gs = hs - t*ss;         /* standard state is any t and 1 bar for gases */
        hs       = 0.0;
        ss       = 0.0;
        vs       = 0.0;
        cps      = 0.0;
        dcpsdt   = 0.0;
        dvsdt    = 0.0;
        dvsdp    = 0.0;
        d2vsdt2  = 0.0;
        d2vsdp2  = 0.0;
        d2vsdtdp = 0.0;
#else
        /* This is an ideal gas treatment of the ss */
        hs = 23.10248*(t-tr) + 2.0*804.8876*(sqrt(t)-sqrt(tr))
        - 1762835.0*(1.0/t-1.0/tr)
        - 18172.91960*log(t/tr) + 0.5*0.002676*(t*t-tr*tr);
        ss = 205.15 + 23.10248*log(t/tr)
        - 2.0*804.8876*(1.0/sqrt(t)-1.0/sqrt(tr))
        - 0.5*1762835.0*(1.0/(t*t)-1.0/(tr*tr))
        + 18172.91960*(1.0/t-1.0/tr) + 0.002676*(t-tr) - r*log(p/pr);
        vs = r*t/p;
        cps = 23.10248 + 804.8876/sqrt(t) + 1762835.0/SQUARE(t)
        - 18172.91960/t + 0.002676*t;
        dcpsdt = - 0.5*804.8876/pow(t, (double) 1.5) - 2.0*1762835.0/CUBE(t)
        + 18172.91960/SQUARE(t) + 0.002676;
        dvsdt = r/p;
        dvsdp = - r*t/SQUARE(p);
        d2vsdt2 = 0.0;
        d2vsdp2 = 2.0*r*t/CUBE(p);
        d2vsdtdp = - r/SQUARE(p);
        gs = hs - t*ss;         /* standard state is any t and 1 bar for gases */
#endif

    } else if(strcmp(name, "s2") == 0) {
        /* This is an ideal gas treatment of the ss Data are from the JANAF tables */
        hs = 128600.0 + 32.2768105*(t-tr) + 4.35324992e-3*(t*t-tr*tr)/2.0 + 3.40018077e5*(1.0/t-1.0/tr)
        - 4.64141632e-7*(t*t*t-tr*tr*tr)/3.0 + 2.0*47.1832248*(sqrt(t)-sqrt(tr));
        ss = 228.165 + 32.2768105*log(t/tr) + 4.35324992e-3*(t-tr) + (3.40018077e5/2.0)*(1.0/(t*t) - 1.0/(tr*tr))
        - 4.64141632e-7*(t*t-tr*tr)/2.0 - (47.1832248/2.0)*(1.0/sqrt(t) - 1.0/sqrt(tr)) - r*log(p/pr);
        vs = r*t/p;
        cps = 32.2768105 + 4.35324992e-3*t - 3.40018077e5/(t*t) - 4.64141632e-7*t*t + 47.1832248/sqrt(t);
        dcpsdt = 4.35324992e-3 + 2.0*3.40018077e5/(t*t*t) - 2.0*4.64141632e-7*t - 0.5*47.1832248/(t*sqrt(t));;
        dvsdt = r/p;
        dvsdp = - r*t/SQUARE(p);
        d2vsdt2 = 0.0;
        d2vsdp2 = 2.0*r*t/CUBE(p);
        d2vsdtdp = - r/SQUARE(p);
        gs = hs - t*ss;         /* standard state is any t and 1 bar for gases */

    } else if(strcmp(name, "quartz") == 0) {
        double lambdaV = 0.0, lambdadVdt = 0.0, lambdadVdp = 0.0,
        lambdad2Vdt2 = 0.0, lambdad2Vdtdp = 0.0, lambdad2Vdp2 = 0.0;

        hs   = phase->h;
        ss   = phase->s;
        k0   = phase->cp.Berman.k0;
        k1   = phase->cp.Berman.k1;
        k2   = phase->cp.Berman.k2;
        k3   = phase->cp.Berman.k3;
        cp_t = phase->cp.Berman.Tt + 0.0237*(p-1.0); /*  Berman (1988) dt/dp */
        cp_h = phase->cp.Berman.deltah;
        l1   = phase->cp.Berman.l1;
        l2   = phase->cp.Berman.l2;
        if (phase->eos_type == EOS_BERMAN) {
            vR = phase->v;
            v1 = phase->eos.Berman.v1;
            v2 = phase->eos.Berman.v2;
            v3 = phase->eos.Berman.v3;
            v4 = phase->eos.Berman.v4;
        }

        cps = k0 + k1/sqrt(t) + k2/SQUARE(t) + k3/CUBE(t);
        dcpsdt = - 0.5*k1/pow(t, (double) 1.5) - 2.0*k2/CUBE(t) - 3.0*k3/QUARTIC(t);

        if(t > cp_t) {
            hs = -908627.0;         /* Properties of beta-quartz (Berman, 1988) */
            ss = 44.207;
            if (phase->eos_type == EOS_BERMAN) {
                phase->v             = 2.370;
                phase->eos.Berman.v1 = -1.238e-6;
                phase->eos.Berman.v2 = 7.087e-13;
                phase->eos.Berman.v3 = 0.0;
                phase->eos.Berman.v4 = 0.0;
            }
            hs = hs + k0*(t-tr) + 2.0*k1*(sqrt(t)-sqrt(tr))
            - k2*(1.0/t-1.0/tr) - 0.5*k3*(1.0/(t*t)-1.0/(tr*tr));
            ss = ss + k0*log(t/tr)
            - 2.0*k1*(1.0/sqrt(t)-1.0/sqrt(tr))
            - 0.5*k2*(1.0/(t*t)-1.0/(tr*tr))
            - (1.0/3.0)*k3*(1.0/(t*t*t)-1.0/(tr*tr*tr));
        } else {
            double delt, x1, x2, x3, x4;
            double DdeltDp = - 0.0237;
            double Dx1Dp, D2x1Dp2, D3x1Dp3, Dx2Dp, D2x2Dp2, D3x2Dp3,
            Dx3Dp, D2x3Dp2, D3x3Dp3, Dx4Dp, D2x4Dp2, D3x4Dp3;

            hs = hs + k0*(t-tr) + 2.0*k1*(sqrt(t)-sqrt(tr))
            - k2*(1.0/t-1.0/tr) - 0.5*k3*(1.0/(t*t)-1.0/(tr*tr));
            ss = ss + k0*log(t/tr)
            - 2.0*k1*(1.0/sqrt(t)-1.0/sqrt(tr))
            - 0.5*k2*(1.0/(t*t)-1.0/(tr*tr))
            - (1.0/3.0)*k3*(1.0/(t*t*t)-1.0/(tr*tr*tr));
            delt = phase->cp.Berman.Tt - cp_t;

            x1 = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
            Dx1Dp   = l1*l1*DdeltDp + 4.0*l1*l2*delt*DdeltDp
            + 3.0*l2*l2*delt*delt*DdeltDp;
            D2x1Dp2 = 4.0*l1*l2*DdeltDp*DdeltDp + 6.0*l2*l2*delt*DdeltDp*DdeltDp;
            D3x1Dp3 = 6.0*l2*l2*DdeltDp*DdeltDp*DdeltDp;

            x2 = l1*l1 + 4.0*l1*l2*delt + 3.0*l2*l2*delt*delt;
            Dx2Dp   = 4.0*l1*l2*DdeltDp + 6.0*l2*l2*delt*DdeltDp;
            D2x2Dp2 = 6.0*l2*l2*DdeltDp*DdeltDp;
            D3x2Dp3 = 0.0;

            x3 = 2.0*l1*l2 + 3.0*l2*l2*delt;
            Dx3Dp   = 3.0*l2*l2*DdeltDp;
            D2x3Dp2 = 0.0;
            D3x3Dp3 = 0.0;

            x4 = l2*l2;
            Dx4Dp   = 0.0;
            D2x4Dp2 = 0.0;
            D3x4Dp3 = 0.0;

            hs = hs + x1*(t-(373.0-delt))
            + x2*(t*t-SQUARE(373.0-delt))/2.0
            + x3*(t*t*t-CUBE(373.0-delt))/3.0
            + x4*(t*t*t*t-QUARTIC(373.0-delt))/4.0;
            ss = ss + x1*(log(t)-log(373.0-delt))
            + x2*(t-(373.0-delt))
            + x3*(t*t-SQUARE(373.0-delt))/2.0
            + x4*(t*t*t-CUBE(373.0-delt))/3.0;
            lambdaV = Dx1Dp*(t-(373.0-delt)) + Dx2Dp*(t*t-SQUARE(373.0-delt))/2.0
            + Dx3Dp*(t*t*t-CUBE(373.0-delt))/3.0
            + Dx4Dp*(t*t*t*t-QUARTIC(373.0-delt))/4.0
            - t*( Dx1Dp*(log(t)-log(373.0-delt)) + Dx2Dp*(t-(373.0-delt))
                 + Dx3Dp*(t*t-SQUARE(373.0-delt))/2.0
                 + Dx4Dp*(t*t*t-CUBE(373.0-delt))/3.0
                 )
            + x1*DdeltDp + x2*(373.0-delt)*DdeltDp
            + x3*SQUARE(373.0-delt)*DdeltDp
            + x4*CUBE(373.0-delt)*DdeltDp
            - t*(x1*DdeltDp/(373.0-delt) + x2*DdeltDp + x3*(373.0-delt)*DdeltDp
                 + x4*SQUARE(373.0-delt)*DdeltDp
                 );
            lambdadVdt = Dx1Dp + Dx2Dp*t + Dx3Dp*t*t + Dx4Dp*t*t*t
            - ( Dx1Dp*(log(t)-log(373.0-delt)) + Dx2Dp*(t-(373.0-delt))
               + Dx3Dp*(t*t-SQUARE(373.0-delt))/2.0
               + Dx4Dp*(t*t*t-CUBE(373.0-delt))/3.0
               ) - t*(Dx1Dp/t + Dx2Dp + Dx3Dp*t + Dx4Dp*t*t)
            - (x1*DdeltDp/(373.0-delt) + x2*DdeltDp + x3*(373.0-delt)*DdeltDp
               + x4*SQUARE(373.0-delt)*DdeltDp
               );
            lambdadVdp = D2x1Dp2*(t-(373.0-delt)) + Dx1Dp*DdeltDp
            + D2x2Dp2*(t*t-SQUARE(373.0-delt))/2.0 + Dx2Dp*(373.0-delt)*DdeltDp
            + D2x3Dp2*(t*t*t-CUBE(373.0-delt))/3.0
            + Dx3Dp*SQUARE(373.0-delt)*DdeltDp
            + D2x4Dp2*(t*t*t*t-QUARTIC(373.0-delt))/4.0
            + Dx4Dp*CUBE(373.0-delt)*DdeltDp
            - t*(  D2x1Dp2*(log(t)-log(373.0-delt)) + Dx1Dp*DdeltDp/(373.0-delt)
                 + D2x2Dp2*(t-(373.0-delt)) + Dx2Dp*DdeltDp
                 + D2x3Dp2*(t*t-SQUARE(373.0-delt))/2.0 + Dx3Dp*(373.0-delt)*DdeltDp
                 + D2x4Dp2*(t*t*t-CUBE(373.0-delt))/3.0
                 + Dx4Dp*SQUARE(373.0-delt)*DdeltDp
                 )
            + Dx1Dp*DdeltDp
            + Dx2Dp*(373.0-delt)*DdeltDp - x2*DdeltDp*DdeltDp
            + Dx3Dp*SQUARE(373.0-delt)*DdeltDp
            - x3*2.0*(373.0-delt)*SQUARE(DdeltDp)
            + Dx4Dp*CUBE(373.0-delt)*DdeltDp
            - x4*3.0*SQUARE(373.0-delt)*SQUARE(DdeltDp)
            - t*(  Dx1Dp*DdeltDp/(373.0-delt)
                 + x1*SQUARE(DdeltDp)/SQUARE(373.0-delt)
                 + Dx2Dp*DdeltDp
                 + Dx3Dp*(373.0-delt)*DdeltDp - x3*SQUARE(DdeltDp)
                 + Dx4Dp*SQUARE(373.0-delt)*DdeltDp
                 - x4*2.0*(373.0-delt)*SQUARE(DdeltDp)
                 );
            lambdad2Vdt2 = Dx2Dp + 2.0*Dx3Dp*t + 3.0*Dx4Dp*t*t
            - 2.0*(Dx1Dp/t + Dx2Dp + Dx3Dp*t + Dx4Dp*t*t)
            - t*(-Dx1Dp/(t*t) + Dx3Dp + 2.0*Dx4Dp*t);
            lambdad2Vdtdp = D2x1Dp2 + D2x2Dp2*t + D2x3Dp2*t*t + D2x4Dp2*t*t*t
            - (  D2x1Dp2*(log(t)-log(373.0-delt)) + Dx1Dp*DdeltDp/(373.0-delt)
               + D2x2Dp2*(t-(373.0-delt)) + Dx2Dp*DdeltDp
               + D2x3Dp2*(t*t-SQUARE(373.0-delt))/2.0 + Dx3Dp*(373.0-delt)*DdeltDp
               + D2x4Dp2*(t*t*t-CUBE(373.0-delt))/3.0
               + Dx4Dp*SQUARE(373.0-delt)*DdeltDp
               )
            - t*(D2x1Dp2/t + D2x2Dp2 + D2x3Dp2*t + D2x4Dp2*t*t)
            - (  Dx1Dp*DdeltDp/(373.0-delt)
               + x1*SQUARE(DdeltDp)/SQUARE(373.0-delt)
               + Dx2Dp*DdeltDp
               + Dx3Dp*(373.0-delt)*DdeltDp - x3*SQUARE(DdeltDp)
               + Dx4Dp*SQUARE(373.0-delt)*DdeltDp
               - x4*2.0*(373.0-delt)*SQUARE(DdeltDp)
               );
            lambdad2Vdp2 = D3x1Dp3*(t-(373.0-delt)) + D2x1Dp2*DdeltDp
            + D2x1Dp2*DdeltDp
            + D3x2Dp3*(t*t-SQUARE(373.0-delt))/2.0 + D2x2Dp2*(373.0-delt)*DdeltDp
            + D2x2Dp2*(373.0-delt)*DdeltDp - Dx2Dp*SQUARE(DdeltDp)
            + D3x3Dp3*(t*t*t-CUBE(373.0-delt))/3.0
            + D2x3Dp2*SQUARE(373.0-delt)*DdeltDp
            + D2x3Dp2*SQUARE(373.0-delt)*DdeltDp
            - Dx3Dp*2.0*(373.0-delt)*SQUARE(DdeltDp)
            + D3x4Dp3*(t*t*t*t-QUARTIC(373.0-delt))/4.0
            + D2x4Dp2*CUBE(373.0-delt)*DdeltDp
            + D2x4Dp2*CUBE(373.0-delt)*DdeltDp
            - Dx4Dp*3.0*SQUARE(373.0-delt)*SQUARE(DdeltDp)
            - t*(  D3x1Dp3*(log(t)-log(373.0-delt)) + D2x1Dp2*DdeltDp/(373.0-delt)
                 + D2x1Dp2*DdeltDp/(373.0-delt)
                 + Dx1Dp*SQUARE(DdeltDp/(373.0-delt))
                 + D3x2Dp3*(t-(373.0-delt)) + D2x2Dp2*DdeltDp
                 + D2x2Dp2*DdeltDp
                 + D3x3Dp3*(t*t-SQUARE(373.0-delt))/2.0
                 + D2x3Dp2*(373.0-delt)*DdeltDp
                 + D2x3Dp2*(373.0-delt)*DdeltDp - Dx3Dp*SQUARE(DdeltDp)
                 + D3x4Dp3*(t*t*t-CUBE(373.0-delt))/3.0
                 + D2x4Dp2*SQUARE(373.0-delt)*DdeltDp
                 + D2x4Dp2*SQUARE(373.0-delt)*DdeltDp
                 - Dx4Dp*2.0*(373.0-delt)*SQUARE(DdeltDp)
                 )
            + D2x1Dp2*DdeltDp
            + D2x2Dp2*(373.0-delt)*DdeltDp - Dx2Dp*SQUARE(DdeltDp)
            - Dx2Dp*DdeltDp*DdeltDp
            + D2x3Dp2*SQUARE(373.0-delt)*DdeltDp
            - Dx3Dp*2.0*(373.0-delt)*SQUARE(DdeltDp)
            - Dx3Dp*2.0*(373.0-delt)*SQUARE(DdeltDp) + x3*2.0*CUBE(DdeltDp)
            + D2x4Dp2*CUBE(373.0-delt)*DdeltDp
            - Dx4Dp*3.0*SQUARE(373.0-delt)*SQUARE(DdeltDp)
            - Dx4Dp*3.0*SQUARE(373.0-delt)*SQUARE(DdeltDp)
            + x4*6.0*(373.0-delt)*CUBE(DdeltDp)
            - t*(  D2x1Dp2*DdeltDp/(373.0-delt)
                 + Dx1Dp*SQUARE(DdeltDp/(373.0-delt))
                 + Dx1Dp*SQUARE(DdeltDp/(373.0-delt))
                 + x1*2.0*CUBE(DdeltDp/(373.0-delt))
                 + D2x2Dp2*DdeltDp
                 + D2x3Dp2*(373.0-delt)*DdeltDp - Dx3Dp*SQUARE(DdeltDp)
                 - Dx3Dp*SQUARE(DdeltDp)
                 + D2x4Dp2*SQUARE(373.0-delt)*DdeltDp
                 - Dx4Dp*2.0*(373.0-delt)*SQUARE(DdeltDp)
                 - Dx4Dp*2.0*(373.0-delt)*SQUARE(DdeltDp)
                 + x4*2.0*CUBE(DdeltDp)
                 );
            cps += (t+delt)*SQUARE(l1+l2*(t+delt));
            dcpsdt += SQUARE(l1+l2*(t+delt)) + (t+delt)*2.0*(l1+l2*(t+delt))*l2;
        }

        /* Do equation of state integral */
        gs = hs -t*ss + (QUARTZ_ADJUSTMENT)+(WETTING_ANGLE_CORR);
        hs += (QUARTZ_ADJUSTMENT)+(WETTING_ANGLE_CORR);
        intEOSsolid(phase, t, p, &gs, &hs, &ss, &cps, &dcpsdt, &vs, &dvsdt, &dvsdp, &d2vsdt2, &d2vsdtdp, &d2vsdp2);
        if (phase->eos_type == EOS_BERMAN) {
            phase->v             = vR;
            phase->eos.Berman.v1 = v1;
            phase->eos.Berman.v2 = v2;
            phase->eos.Berman.v3 = v3;
            phase->eos.Berman.v4 = v4;
        }

        vs       += lambdaV;
        dvsdt    += lambdadVdt;
        dvsdp    += lambdadVdp;
        d2vsdt2  += lambdad2Vdt2;
        d2vsdtdp += lambdad2Vdtdp;
        d2vsdp2  += lambdad2Vdp2;

    } else if(strcmp(name, "tridymite") == 0) {
        hs   = phase->h;
        ss   = phase->s;
        k0   = phase->cp.Berman.k0;
        k1   = phase->cp.Berman.k1;
        k2   = phase->cp.Berman.k2;
        k3   = phase->cp.Berman.k3;
        cp_t = phase->cp.Berman.Tt;
        cp_h = phase->cp.Berman.deltah;
        l1   = phase->cp.Berman.l1;
        l2   = phase->cp.Berman.l2;
        if (phase->eos_type == EOS_BERMAN) {
            vR = phase->v;
            v1 = phase->eos.Berman.v1;
            v2 = phase->eos.Berman.v2;
            v3 = phase->eos.Berman.v3;
            v4 = phase->eos.Berman.v4;
        }

        cps = k0 + k1/sqrt(t) + k2/SQUARE(t) + k3/CUBE(t);
        dcpsdt = - 0.5*k1/pow(t, (double) 1.5) - 2.0*k2/CUBE(t) - 3.0*k3/QUARTIC(t);

        if(t > cp_t) {
            hs = -907045.0;      /* Properties of beta-tridymite (Berman, 1988) */
            ss = 45.524;
            if (phase->eos_type == EOS_BERMAN) {
                phase->v             = 2.737;
                phase->eos.Berman.v1 = -0.740e-6;
                phase->eos.Berman.v2 = 3.735e-12;
                phase->eos.Berman.v3 = 4.829e-6;
                phase->eos.Berman.v4 = 0.0;
            }
            hs = hs + k0*(t-tr) + 2.0*k1*(sqrt(t)-sqrt(tr))
            - k2*(1.0/t-1.0/tr) - 0.5*k3*(1.0/(t*t)-1.0/(tr*tr));
            ss = ss + k0*log(t/tr)
            - 2.0*k1*(1.0/sqrt(t)-1.0/sqrt(tr))
            - 0.5*k2*(1.0/(t*t)-1.0/(tr*tr))
            - (1.0/3.0)*k3*(1.0/(t*t*t)-1.0/(tr*tr*tr));
        } else {
            double delt, x1, x2, x3, x4;
            hs = hs + k0*(t-tr) + 2.0*k1*(sqrt(t)-sqrt(tr))
            - k2*(1.0/t-1.0/tr) - 0.5*k3*(1.0/(t*t)-1.0/(tr*tr));
            ss = ss + k0*log(t/tr)
            - 2.0*k1*(1.0/sqrt(t)-1.0/sqrt(tr))
            - 0.5*k2*(1.0/(t*t)-1.0/(tr*tr))
            - (1.0/3.0)*k3*(1.0/(t*t*t)-1.0/(tr*tr*tr));
            delt = phase->cp.Berman.Tt - cp_t;
            x1 = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
            x2 = l1*l1 + 4.0*l1*l2*delt + 3.0*l2*l2*delt*delt;
            x3 = 2.0*l1*l2 + 3.0*l2*l2*delt;
            x4 = l2*l2;
            hs = hs + x1*(t-(tr-delt))
            + x2*(t*t-SQUARE(tr-delt))/2.0
            + x3*(t*t*t-CUBE(tr-delt))/3.0
            + x4*(t*t*t*t-QUARTIC(tr-delt))/4.0;
            ss = ss + x1*(log(t)-log(tr-delt))
            + x2*(t-(tr-delt))
            + x3*(t*t-SQUARE(tr-delt))/2.0
            + x4*(t*t*t-CUBE(tr-delt))/3.0;
            cps += (t+delt)*SQUARE(l1+l2*(t+delt));
            dcpsdt += SQUARE(l1+l2*(t+delt)) + (t+delt)*2.0*(l1+l2*(t+delt))*l2;
        }

        /* Do equation of state integral */
        gs = hs - t*ss + (TRIDYMITE_ADJUSTMENT);
        hs += (TRIDYMITE_ADJUSTMENT);
        intEOSsolid(phase, t, p, &gs, &hs, &ss, &cps, &dcpsdt, &vs, &dvsdt, &dvsdp, &d2vsdt2, &d2vsdtdp, &d2vsdp2);
        if (phase->eos_type == EOS_BERMAN) {
            phase->v             = vR;
            phase->eos.Berman.v1 = v1;
            phase->eos.Berman.v2 = v2;
            phase->eos.Berman.v3 = v3;
            phase->eos.Berman.v4 = v4;
        }

    } else if(strcmp(name, "cristobalite") == 0) {
        double lambdaV = 0.0, lambdadVdt = 0.0, lambdadVdp = 0.0,
        lambdad2Vdt2 = 0.0, lambdad2Vdtdp = 0.0, lambdad2Vdp2 = 0.0;

        hs   = phase->h;
        ss   = phase->s;
        k0   = phase->cp.Berman.k0;
        k1   = phase->cp.Berman.k1;
        k2   = phase->cp.Berman.k2;
        k3   = phase->cp.Berman.k3;
        cp_t = phase->cp.Berman.Tt + 0.0480*(p-1.0); /*  Berman (1988) dt/dp */
        cp_h = phase->cp.Berman.deltah;
        l1   = phase->cp.Berman.l1;
        l2   = phase->cp.Berman.l2;
        if (phase->eos_type == EOS_BERMAN) {
            vR = phase->v;
            v1 = phase->eos.Berman.v1;
            v2 = phase->eos.Berman.v2;
            v3 = phase->eos.Berman.v3;
            v4 = phase->eos.Berman.v4;
        }

        cps = k0 + k1/sqrt(t) + k2/SQUARE(t) + k3/CUBE(t);
        dcpsdt = - 0.5*k1/pow(t, (double) 1.5) - 2.0*k2/CUBE(t) - 3.0*k3/QUARTIC(t);

        if(t > cp_t) {
            hs = -906377.0;   /* Properties of beta-cristobalite (Berman, 1988) */
            ss = 46.029;
            if (phase->eos_type == EOS_BERMAN) {
                phase->v             = 2.730;
                phase->eos.Berman.v1 = -1.100e-6;
                phase->eos.Berman.v2 = 5.535e-12;
                phase->eos.Berman.v3 = 3.189e-6;
                phase->eos.Berman.v4 = 0.0;
            }
            hs = hs + k0*(t-tr) + 2.0*k1*(sqrt(t)-sqrt(tr))
            - k2*(1.0/t-1.0/tr) - 0.5*k3*(1.0/(t*t)-1.0/(tr*tr));
            ss = ss + k0*log(t/tr)
            - 2.0*k1*(1.0/sqrt(t)-1.0/sqrt(tr))
            - 0.5*k2*(1.0/(t*t)-1.0/(tr*tr))
            - (1.0/3.0)*k3*(1.0/(t*t*t)-1.0/(tr*tr*tr));
        } else {
            double delt, x1, x2, x3, x4;
            double DdeltDp = - 0.0480;
            double Dx1Dp, D2x1Dp2, D3x1Dp3, Dx2Dp, D2x2Dp2, D3x2Dp3,
            Dx3Dp, D2x3Dp2, D3x3Dp3, Dx4Dp, D2x4Dp2, D3x4Dp3;
            hs = hs + k0*(t-tr) + 2.0*k1*(sqrt(t)-sqrt(tr))
            - k2*(1.0/t-1.0/tr) - 0.5*k3*(1.0/(t*t)-1.0/(tr*tr));
            ss = ss + k0*log(t/tr)
            - 2.0*k1*(1.0/sqrt(t)-1.0/sqrt(tr))
            - 0.5*k2*(1.0/(t*t)-1.0/(tr*tr))
            - (1.0/3.0)*k3*(1.0/(t*t*t)-1.0/(tr*tr*tr));

            delt = phase->cp.Berman.Tt - cp_t;

            x1      = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
            Dx1Dp   = l1*l1*DdeltDp + 4.0*l1*l2*delt*DdeltDp
            + 3.0*l2*l2*delt*delt*DdeltDp;
            D2x1Dp2 = 4.0*l1*l2*DdeltDp*DdeltDp + 6.0*l2*l2*delt*DdeltDp*DdeltDp;
            D3x1Dp3 = 6.0*l2*l2*DdeltDp*DdeltDp*DdeltDp;

            x2      = l1*l1 + 4.0*l1*l2*delt + 3.0*l2*l2*delt*delt;
            Dx2Dp   = 4.0*l1*l2*DdeltDp + 6.0*l2*l2*delt*DdeltDp;
            D2x2Dp2 = 6.0*l2*l2*DdeltDp*DdeltDp;
            D3x2Dp3 = 0.0;

            x3      = 2.0*l1*l2 + 3.0*l2*l2*delt;
            Dx3Dp   = 3.0*l2*l2*DdeltDp;
            D2x3Dp2 = 0.0;
            D3x3Dp3 = 0.0;

            x4      = l2*l2;
            Dx4Dp   = 0.0;
            D2x4Dp2 = 0.0;
            D3x4Dp3 = 0.0;

            hs = hs + x1*(t-(tr-delt))
            + x2*(t*t-SQUARE(tr-delt))/2.0
            + x3*(t*t*t-CUBE(tr-delt))/3.0
            + x4*(t*t*t*t-QUARTIC(tr-delt))/4.0;
            ss = ss + x1*(log(t)-log(tr-delt))
            + x2*(t-(tr-delt))
            + x3*(t*t-SQUARE(tr-delt))/2.0
            + x4*(t*t*t-CUBE(tr-delt))/3.0;
            lambdaV = Dx1Dp*(t-(tr-delt)) + Dx2Dp*(t*t-SQUARE(tr-delt))/2.0
            + Dx3Dp*(t*t*t-CUBE(tr-delt))/3.0
            + Dx4Dp*(t*t*t*t-QUARTIC(tr-delt))/4.0
            - t*( Dx1Dp*(log(t)-log(tr-delt)) + Dx2Dp*(t-(tr-delt))
                 + Dx3Dp*(t*t-SQUARE(tr-delt))/2.0
                 + Dx4Dp*(t*t*t-CUBE(tr-delt))/3.0
                 )
            + x1*DdeltDp + x2*(tr-delt)*DdeltDp + x3*SQUARE(tr-delt)*DdeltDp
            + x4*CUBE(tr-delt)*DdeltDp
            - t*(x1*DdeltDp/(tr-delt) + x2*DdeltDp + x3*(tr-delt)*DdeltDp
                 + x4*SQUARE(tr-delt)*DdeltDp
                 );
            lambdadVdt = Dx1Dp + Dx2Dp*t + Dx3Dp*t*t + Dx4Dp*t*t*t
            - ( Dx1Dp*(log(t)-log(tr-delt)) + Dx2Dp*(t-(tr-delt))
               + Dx3Dp*(t*t-SQUARE(tr-delt))/2.0
               + Dx4Dp*(t*t*t-CUBE(tr-delt))/3.0
               ) - t*(Dx1Dp/t + Dx2Dp + Dx3Dp*t + Dx4Dp*t*t)
            - (x1*DdeltDp/(tr-delt) + x2*DdeltDp + x3*(tr-delt)*DdeltDp
               + x4*SQUARE(tr-delt)*DdeltDp
               );
            lambdadVdp = D2x1Dp2*(t-(tr-delt)) + Dx1Dp*DdeltDp
            + D2x2Dp2*(t*t-SQUARE(tr-delt))/2.0 + Dx2Dp*(tr-delt)*DdeltDp
            + D2x3Dp2*(t*t*t-CUBE(tr-delt))/3.0 + Dx3Dp*SQUARE(tr-delt)*DdeltDp
            + D2x4Dp2*(t*t*t*t-QUARTIC(tr-delt))/4.0 + Dx4Dp*CUBE(tr-delt)*DdeltDp
            - t*(  D2x1Dp2*(log(t)-log(tr-delt)) + Dx1Dp*DdeltDp/(tr-delt)
                 + D2x2Dp2*(t-(tr-delt)) + Dx2Dp*DdeltDp
                 + D2x3Dp2*(t*t-SQUARE(tr-delt))/2.0 + Dx3Dp*(tr-delt)*DdeltDp
                 + D2x4Dp2*(t*t*t-CUBE(tr-delt))/3.0 + Dx4Dp*SQUARE(tr-delt)*DdeltDp
                 )
            + Dx1Dp*DdeltDp
            + Dx2Dp*(tr-delt)*DdeltDp - x2*DdeltDp*DdeltDp
            + Dx3Dp*SQUARE(tr-delt)*DdeltDp - x3*2.0*(tr-delt)*SQUARE(DdeltDp)
            + Dx4Dp*CUBE(tr-delt)*DdeltDp - x4*3.0*SQUARE(tr-delt)*SQUARE(DdeltDp)
            - t*(  Dx1Dp*DdeltDp/(tr-delt) + x1*SQUARE(DdeltDp)/SQUARE(tr-delt)
                 + Dx2Dp*DdeltDp
                 + Dx3Dp*(tr-delt)*DdeltDp - x3*SQUARE(DdeltDp)
                 + Dx4Dp*SQUARE(tr-delt)*DdeltDp - x4*2.0*(tr-delt)*SQUARE(DdeltDp)
                 );
            lambdad2Vdt2 = Dx2Dp + 2.0*Dx3Dp*t + 3.0*Dx4Dp*t*t
            - 2.0*(Dx1Dp/t + Dx2Dp + Dx3Dp*t + Dx4Dp*t*t)
            - t*(-Dx1Dp/(t*t) + Dx3Dp + 2.0*Dx4Dp*t);
            lambdad2Vdtdp = D2x1Dp2 + D2x2Dp2*t + D2x3Dp2*t*t + D2x4Dp2*t*t*t
            - (  D2x1Dp2*(log(t)-log(tr-delt)) + Dx1Dp*DdeltDp/(tr-delt)
               + D2x2Dp2*(t-(tr-delt)) + Dx2Dp*DdeltDp
               + D2x3Dp2*(t*t-SQUARE(tr-delt))/2.0 + Dx3Dp*(tr-delt)*DdeltDp
               + D2x4Dp2*(t*t*t-CUBE(tr-delt))/3.0 + Dx4Dp*SQUARE(tr-delt)*DdeltDp
               )
            - t*(D2x1Dp2/t + D2x2Dp2 + D2x3Dp2*t + D2x4Dp2*t*t)
            - (  Dx1Dp*DdeltDp/(tr-delt) + x1*SQUARE(DdeltDp)/SQUARE(tr-delt)
               + Dx2Dp*DdeltDp
               + Dx3Dp*(tr-delt)*DdeltDp - x3*SQUARE(DdeltDp)
               + Dx4Dp*SQUARE(tr-delt)*DdeltDp - x4*2.0*(tr-delt)*SQUARE(DdeltDp)
               );
            lambdad2Vdp2 = D3x1Dp3*(t-(tr-delt)) + D2x1Dp2*DdeltDp
            + D2x1Dp2*DdeltDp
            + D3x2Dp3*(t*t-SQUARE(tr-delt))/2.0 + D2x2Dp2*(tr-delt)*DdeltDp
            + D2x2Dp2*(tr-delt)*DdeltDp - Dx2Dp*SQUARE(DdeltDp)
            + D3x3Dp3*(t*t*t-CUBE(tr-delt))/3.0 + D2x3Dp2*SQUARE(tr-delt)*DdeltDp
            + D2x3Dp2*SQUARE(tr-delt)*DdeltDp
            - Dx3Dp*2.0*(tr-delt)*SQUARE(DdeltDp)
            + D3x4Dp3*(t*t*t*t-QUARTIC(tr-delt))/4.0
            + D2x4Dp2*CUBE(tr-delt)*DdeltDp
            + D2x4Dp2*CUBE(tr-delt)*DdeltDp
            - Dx4Dp*3.0*SQUARE(tr-delt)*SQUARE(DdeltDp)
            - t*(  D3x1Dp3*(log(t)-log(tr-delt)) + D2x1Dp2*DdeltDp/(tr-delt)
                 + D2x1Dp2*DdeltDp/(tr-delt) + Dx1Dp*SQUARE(DdeltDp/(tr-delt))
                 + D3x2Dp3*(t-(tr-delt)) + D2x2Dp2*DdeltDp
                 + D2x2Dp2*DdeltDp
                 + D3x3Dp3*(t*t-SQUARE(tr-delt))/2.0 + D2x3Dp2*(tr-delt)*DdeltDp
                 + D2x3Dp2*(tr-delt)*DdeltDp - Dx3Dp*SQUARE(DdeltDp)
                 + D3x4Dp3*(t*t*t-CUBE(tr-delt))/3.0
                 + D2x4Dp2*SQUARE(tr-delt)*DdeltDp
                 + D2x4Dp2*SQUARE(tr-delt)*DdeltDp
                 - Dx4Dp*2.0*(tr-delt)*SQUARE(DdeltDp)
                 )
            + D2x1Dp2*DdeltDp
            + D2x2Dp2*(tr-delt)*DdeltDp - Dx2Dp*SQUARE(DdeltDp)
            - Dx2Dp*DdeltDp*DdeltDp
            + D2x3Dp2*SQUARE(tr-delt)*DdeltDp
            - Dx3Dp*2.0*(tr-delt)*SQUARE(DdeltDp)
            - Dx3Dp*2.0*(tr-delt)*SQUARE(DdeltDp) + x3*2.0*CUBE(DdeltDp)
            + D2x4Dp2*CUBE(tr-delt)*DdeltDp
            - Dx4Dp*3.0*SQUARE(tr-delt)*SQUARE(DdeltDp)
            - Dx4Dp*3.0*SQUARE(tr-delt)*SQUARE(DdeltDp)
            + x4*6.0*(tr-delt)*CUBE(DdeltDp)
            - t*(  D2x1Dp2*DdeltDp/(tr-delt) + Dx1Dp*SQUARE(DdeltDp/(tr-delt))
                 + Dx1Dp*SQUARE(DdeltDp/(tr-delt))
                 + x1*2.0*CUBE(DdeltDp/(tr-delt))
                 + D2x2Dp2*DdeltDp
                 + D2x3Dp2*(tr-delt)*DdeltDp - Dx3Dp*SQUARE(DdeltDp)
                 - Dx3Dp*SQUARE(DdeltDp)
                 + D2x4Dp2*SQUARE(tr-delt)*DdeltDp
                 - Dx4Dp*2.0*(tr-delt)*SQUARE(DdeltDp)
                 - Dx4Dp*2.0*(tr-delt)*SQUARE(DdeltDp)
                 + x4*2.0*CUBE(DdeltDp)
                 );
            cps += (t+delt)*SQUARE(l1+l2*(t+delt));
            dcpsdt += SQUARE(l1+l2*(t+delt)) + (t+delt)*2.0*(l1+l2*(t+delt))*l2;
        }

        /* Do equation of state integral */
        gs = hs - t*ss + (CRISTOBALITE_ADJUSTMENT);
        hs += (CRISTOBALITE_ADJUSTMENT);
        intEOSsolid(phase, t, p, &gs, &hs, &ss, &cps, &dcpsdt, &vs, &dvsdt, &dvsdp, &d2vsdt2, &d2vsdtdp, &d2vsdp2);
        if (phase->eos_type == EOS_BERMAN) {
            phase->v             = vR;
            phase->eos.Berman.v1 = v1;
            phase->eos.Berman.v2 = v2;
            phase->eos.Berman.v3 = v3;
            phase->eos.Berman.v4 = v4;
        }

        vs       += lambdaV;
        dvsdt    += lambdadVdt;
        dvsdp    += lambdadVdp;
        d2vsdt2  += lambdad2Vdt2;
        d2vsdtdp += lambdad2Vdtdp;
        d2vsdp2  += lambdad2Vdp2;

    } else if(strcmp(name, "enstatite") == 0) {
        double h_cen, s_cen, v_cen, g_cen, cp_cen, dcpdt_cen, dvdt_cen,
        dvdp_cen, d2vdt2_cen, d2vdp2_cen, d2vdtdp_cen;
        double h_oen, s_oen, v_oen, g_oen, cp_oen, dcpdt_oen, dvdt_oen,
        dvdp_oen, d2vdt2_oen, d2vdp2_oen, d2vdtdp_oen;
        double h_pen, s_pen, v_pen, g_pen, cp_pen, dcpdt_pen, dvdt_pen,
        dvdp_pen, d2vdt2_pen, d2vdp2_pen, d2vdtdp_pen;

        hs = -1545926.0*2.0;                 /* clino-enstatite (Berman, 1988) */
        ss = 66.325*2.0;
        k0 = 139.96*2.0;
        k1 = -4.970e2*2.0;
        k2 = -44.002e5*2.0;
        k3 = 53.571e7*2.0;
        if (phase->eos_type == EOS_BERMAN) {
            vR = phase->v;
            v1 = phase->eos.Berman.v1;
            v2 = phase->eos.Berman.v2;
            v3 = phase->eos.Berman.v3;
            v4 = phase->eos.Berman.v4;
            phase->v	     = 3.131*2.0;
            phase->eos.Berman.v1 = -0.750e-6;
            phase->eos.Berman.v2 = 0.448e-12;
            phase->eos.Berman.v3 = 21.915e-6;
            phase->eos.Berman.v4 = 74.920e-10;
        }

        h_cen = hs + k0*(t-tr) + 2.0*k1*(sqrt(t)-sqrt(tr)) - k2*(1.0/t-1.0/tr) - 0.5*k3*(1.0/(t*t)-1.0/(tr*tr));
        s_cen = ss + k0*log(t/tr) - 2.0*k1*(1.0/sqrt(t)-1.0/sqrt(tr)) - 0.5*k2*(1.0/(t*t)-1.0/(tr*tr))
        - (1.0/3.0)*k3*(1.0/(t*t*t)-1.0/(tr*tr*tr));
        cp_cen      = k0 + k1/sqrt(t) + k2/SQUARE(t) + k3/CUBE(t);
        dcpdt_cen   = - 0.5*k1/pow(t, (double) 1.5) - 2.0*k2/CUBE(t) - 3.0*k3/QUARTIC(t);

        g_cen = h_cen - t*s_cen;
        intEOSsolid(phase, t, p, &g_cen, &h_cen, &s_cen, &cp_cen, &dcpdt_cen, &v_cen, &dvdt_cen, &dvdp_cen, &d2vdt2_cen, &d2vdtdp_cen, &d2vdp2_cen);
        if (phase->eos_type == EOS_BERMAN) {
            phase->v             = vR;
            phase->eos.Berman.v1 = v1;
            phase->eos.Berman.v2 = v2;
            phase->eos.Berman.v3 = v3;
            phase->eos.Berman.v4 = v4;
        }

        hs = -1545552.0*2.0;                 /* ortho-enstatite (Berman, 1988) */
        ss = 66.170*2.0;
        k0 = 166.58*2.0;
        k1 = -12.006e2*2.0;
        k2 = -22.706e5*2.0;
        k3 = 27.915e7*2.0;
        if (phase->eos_type == EOS_BERMAN) {
            vR = phase->v;
            v1 = phase->eos.Berman.v1;
            v2 = phase->eos.Berman.v2;
            v3 = phase->eos.Berman.v3;
            v4 = phase->eos.Berman.v4;
            phase->v	      = 3.133*2.0;
            phase->eos.Berman.v1 = -0.749e-6;
            phase->eos.Berman.v2 = 0.447e-12;
            phase->eos.Berman.v3 = 24.656e-6;
            phase->eos.Berman.v4 = 74.670e-10;
        }

        h_oen = hs + k0*(t-tr) + 2.0*k1*(sqrt(t)-sqrt(tr)) - k2*(1.0/t-1.0/tr) - 0.5*k3*(1.0/(t*t)-1.0/(tr*tr));
        s_oen = ss + k0*log(t/tr) - 2.0*k1*(1.0/sqrt(t)-1.0/sqrt(tr)) - 0.5*k2*(1.0/(t*t)-1.0/(tr*tr))
        - (1.0/3.0)*k3*(1.0/(t*t*t)-1.0/(tr*tr*tr));
        cp_oen      = k0 + k1/sqrt(t) + k2/SQUARE(t) + k3/CUBE(t);
        dcpdt_oen   = - 0.5*k1/pow(t, (double) 1.5) - 2.0*k2/CUBE(t) - 3.0*k3/QUARTIC(t);

        g_oen = h_oen - t*s_oen;
        intEOSsolid(phase, t, p, &g_oen, &h_oen, &s_oen, &cp_oen, &dcpdt_oen, &v_oen, &dvdt_oen, &dvdp_oen, &d2vdt2_oen, &d2vdtdp_oen, &d2vdp2_oen);
        if (phase->eos_type == EOS_BERMAN) {
            phase->v             = vR;
            phase->eos.Berman.v1 = v1;
            phase->eos.Berman.v2 = v2;
            phase->eos.Berman.v3 = v3;
            phase->eos.Berman.v4 = v4;
        }

        hs = -1543959.0*2.0;                 /* proto-enstatite (Berman, 1988) */
        ss = 67.438*2.0;
        k0 = 166.58*2.0;
        k1 = -12.006e2*2.0;
        k2 = -22.706e5*2.0;
        k3 = 27.915e7*2.0;
        if (phase->eos_type == EOS_BERMAN) {
            vR = phase->v;
            v1 = phase->eos.Berman.v1;
            v2 = phase->eos.Berman.v2;
            v3 = phase->eos.Berman.v3;
            v4 = phase->eos.Berman.v4;
            phase->v	     = 3.242*2.0;
            phase->eos.Berman.v1 = -0.750e-6;
            phase->eos.Berman.v2 = 0.448e-12;
            phase->eos.Berman.v3 = 16.832e-6;
            phase->eos.Berman.v4 = 116.650e-10;
        }

        h_pen = hs + k0*(t-tr) + 2.0*k1*(sqrt(t)-sqrt(tr)) - k2*(1.0/t-1.0/tr) - 0.5*k3*(1.0/(t*t)-1.0/(tr*tr));
        s_pen = ss + k0*log(t/tr) - 2.0*k1*(1.0/sqrt(t)-1.0/sqrt(tr)) - 0.5*k2*(1.0/(t*t)-1.0/(tr*tr))
        - (1.0/3.0)*k3*(1.0/(t*t*t)-1.0/(tr*tr*tr));
        cp_pen      = k0 + k1/sqrt(t) + k2/SQUARE(t) + k3/CUBE(t);
        dcpdt_pen   = - 0.5*k1/pow(t, (double) 1.5) - 2.0*k2/CUBE(t) - 3.0*k3/QUARTIC(t);

        g_pen = h_pen - t*s_pen;
        intEOSsolid(phase, t, p, &g_pen, &h_pen, &s_pen, &cp_pen, &dcpdt_pen, &v_pen, &dvdt_pen, &dvdp_pen, &d2vdt2_pen, &d2vdtdp_pen, &d2vdp2_pen);
        if (phase->eos_type == EOS_BERMAN) {
            phase->v             = vR;
            phase->eos.Berman.v1 = v1;
            phase->eos.Berman.v2 = v2;
            phase->eos.Berman.v3 = v3;
            phase->eos.Berman.v4 = v4;
        }

        if((g_cen <= g_oen) && (g_cen <= g_pen)) {
            hs = h_cen; ss = s_cen; gs = g_cen; vs = v_cen; dvsdt = dvdt_cen;
            dvsdp = dvdp_cen; d2vsdt2 = d2vdt2_cen; d2vsdp2 = d2vdp2_cen;
            d2vsdtdp = d2vdtdp_cen; cps = cp_cen; dcpsdt = dcpdt_cen;
        } else if((g_oen <= g_cen) && (g_oen <= g_pen)) {
            hs = h_oen; ss = s_oen; gs = g_oen; vs = v_oen; dvsdt = dvdt_oen;
            dvsdp = dvdp_oen; d2vsdt2 = d2vdt2_oen; d2vsdp2 = d2vdp2_oen;
            d2vsdtdp = d2vdtdp_oen; cps = cp_oen; dcpsdt = dcpdt_oen;
        } else if((g_pen <= g_oen) && (g_pen <= g_cen)) {
            hs = h_pen; ss = s_pen; gs = g_pen; vs = v_pen; dvsdt = dvdt_pen;
            dvsdp = dvdp_pen; d2vsdt2 = d2vdt2_pen; d2vsdp2 = d2vdp2_pen;
            d2vsdtdp = d2vdtdp_pen; cps = cp_pen; dcpsdt = dcpdt_pen;
        } else {
            printf("Absurd logic problem in GIBBS\n");
            exit(1);
        }
    } else if(strcmp(name, "fluid") == 0) { // i.e. "water"
        double gH2O=0.0, hH2O=0.0, sH2O=0.0, cpH2O=0.0, dcpdtH2O=0.0, vH2O=0.0, dvdtH2O=0.0, dvdpH2O=0.0,
        d2vdt2H2O=0.0, d2vdtdpH2O=0.0, d2vdp2H2O=0.0;

#ifdef TRUE_xMELTS
	if ((calculationMode == MODE_xMELTS) || (calculationMode == MODE_pMELTS)) {
#else
        if (calculationMode == MODE_pMELTS) {
#endif

            double x[2] = { 1.0, 0.0};
            fluidPhase(t, p, x, &gH2O, NULL, NULL, NULL, &hH2O, &sH2O,
                       NULL, NULL, &cpH2O, &dcpdtH2O, NULL, &vH2O, NULL, NULL, &dvdtH2O,
                       &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O, NULL, NULL, NULL, NULL,
                       NULL);
            gH2O        = hH2O - t*sH2O; /* used gH2O for calibration */
            vH2O       *=10.0;
            dvdtH2O    *=10.0;
            dvdpH2O    *=10.0;
            d2vdt2H2O  *=10.0;
            d2vdtdpH2O *=10.0;
            d2vdp2H2O  *=10.0;
        } else if ((calculationMode == MODE__MELTS) ||
#ifndef TRUE_xMELTS
                   (calculationMode == MODE_xMELTS) ||
#endif
                   (calculationMode == MODE__MELTSandCO2) ||
                   (calculationMode == MODE__MELTSandCO2_H2O)) {
            double pHaar = (p <= 10000.0) ? p : 10000.0;
            whaar(pHaar, t, &gH2O, &hH2O, &sH2O, &cpH2O, &dcpdtH2O, &vH2O, &dvdtH2O,
                  &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O);

            if(p > 10000.0) {
                double gDELTA, hDELTA, sDELTA, cpDELTA, dcpdtDELTA;

                wdh78(p, t, &gDELTA, &hDELTA, &sDELTA, &cpDELTA, &dcpdtDELTA, &vH2O,
                      &dvdtH2O, &dvdpH2O, &d2vdt2H2O, &d2vdtdpH2O, &d2vdp2H2O);

                gH2O       += gDELTA;
                hH2O       += hDELTA;
                sH2O       += sDELTA;
                cpH2O      += cpDELTA;
                dcpdtH2O   += dcpdtDELTA;
            }
        }

        gs       = gH2O;
        ss       = sH2O;
        hs       = hH2O; /* gH2O + t*sH2O */;
        vs       = vH2O/10.0;
        cps      = cpH2O;
        dcpsdt   = dcpdtH2O;
        vs       = vH2O/10.0;
        dvsdt    = dvdtH2O/10.0;
        dvsdp    = dvdpH2O/10.0;
        d2vsdt2  = d2vdt2H2O/10.0;
        d2vsdtdp = d2vdtdpH2O/10.0;
        d2vsdp2  = d2vdp2H2O/10.0;

    } else if(strcmp(name, "h2oduan") == 0) {
        propertiesOfPureH2O(t, p, &gs, &hs, &ss, &cps, &dcpsdt, &vs, &dvsdt, &dvsdp, &d2vsdt2, &d2vsdtdp, &d2vsdp2);

    } else if(strcmp(name, "co2duan") == 0) {
        propertiesOfPureCO2(t, p, &gs, &hs, &ss, &cps, &dcpsdt, &vs, &dvsdt, &dvsdp, &d2vsdt2, &d2vsdtdp, &d2vsdp2);

    } else if(strcmp(name, "Fe-metal") == 0) {
        fe_metal(t, p, &gs, &hs, &ss, &cps, &dcpsdt, &vs, &dvsdt,
                 &dvsdp, &d2vsdt2, &d2vsdp2, &d2vsdtdp);

    } else if(strcmp(name, "Fe-liquid") == 0) {
        fe_liquid(t, p, &gs, &hs, &ss, &cps, &dcpsdt, &vs, &dvsdt,
                  &dvsdp, &d2vsdt2, &d2vsdp2, &d2vsdtdp);

    } else if(strcmp(name, "Ni-metal") == 0) {
        ni_metal(t, p, &gs, &hs, &ss, &cps, &dcpsdt, &vs, &dvsdt,
                 &dvsdp, &d2vsdt2, &d2vsdp2, &d2vsdtdp);

    } else if(strcmp(name, "Ni-liquid") == 0) {
        ni_liquid(t, p, &gs, &hs, &ss, &cps, &dcpsdt, &vs, &dvsdt,
                  &dvsdp, &d2vsdt2, &d2vsdp2, &d2vsdtdp);

        /* Holland and Powell style properties for use by Tajcmanova biotite model */
    } else if (strcmp(name, "fbiTaj") == 0) {
        double gsEAST, ssEAST, hsEAST, cpsEAST, dcpsdtEAST, vsEAST, dvsdtEAST, dvsdpEAST, d2vsdt2EAST, d2vsdtdpEAST, d2vsdp2EAST;
        double gsCOR, ssCOR, hsCOR, cpsCOR, dcpsdtCOR, vsCOR, dvsdtCOR, dvsdpCOR, d2vsdt2COR, d2vsdtdpCOR, d2vsdp2COR;
        double gsHEM, ssHEM, hsHEM, cpsHEM, dcpsdtHEM, vsHEM, dvsdtHEM, dvsdpHEM, d2vsdt2HEM, d2vsdtdpHEM, d2vsdp2HEM;

        HPproperties(t, p, -6338170.0, 318.0, 14.738, 785.5, -0.038031, -2130300.0, -6893.7, 5.79e-5, 513000.0, 0.0, 0.0, 0.0,
                     &gsEAST, &ssEAST, &hsEAST, &cpsEAST, &dcpsdtEAST, &vsEAST, &dvsdtEAST, &dvsdpEAST, &d2vsdt2EAST, &d2vsdtdpEAST,
                     &d2vsdp2EAST);
        HPproperties(t, p, -1675140.0, 50.90, 2.558,  139.5,  0.005890, -2460600.0,  -589.2, 4.19e-5, 2520000.0, 0.0, 0.0, 0.0,
                     &gsCOR, &ssCOR, &hsCOR, &cpsCOR, &dcpsdtCOR, &vsCOR, &dvsdtCOR, &dvsdpCOR, &d2vsdt2COR, &d2vsdtdpCOR,
                     &d2vsdp2COR);
        HPproperties(t, p,  -825690.0, 87.4,  3.0274, 163.9, 0.0,       -2257200.0,  -657.6, 5.99e-5, 1996000.0, 955.0, 15.6, 0.0,
                     &gsHEM, &ssHEM, &hsHEM, &cpsHEM, &dcpsdtHEM, &vsHEM, &dvsdtHEM, &dvsdpHEM, &d2vsdt2HEM, &d2vsdtdpHEM,
                     &d2vsdp2HEM);
        gs       =  gsEAST	      - 0.5*gsCOR	+ 0.5*gsHEM;
        ss       =  ssEAST	      - 0.5*ssCOR	+ 0.5*ssHEM;
        hs       =  hsEAST	      - 0.5*hsCOR	+ 0.5*hsHEM;
        cps      =  cpsEAST      - 0.5*cpsCOR	+ 0.5*cpsHEM;
        dcpsdt   =  dcpsdtEAST   - 0.5*dcpsdtCOR	+ 0.5*dcpsdtHEM;
        vs       =  vsEAST	      - 0.5*vsCOR	+ 0.5*vsHEM;
        dvsdt    =  dvsdtEAST    - 0.5*dvsdtCOR	+ 0.5*dvsdtHEM;
        dvsdp    =  dvsdpEAST    - 0.5*dvsdpCOR	+ 0.5*dvsdpHEM;
        d2vsdt2  =  d2vsdt2EAST  - 0.5*d2vsdt2COR  + 0.5*d2vsdt2HEM;
        d2vsdtdp =  d2vsdtdpEAST - 0.5*d2vsdtdpCOR + 0.5*d2vsdtdpHEM;
        d2vsdp2  =  d2vsdp2EAST  - 0.5*d2vsdp2COR  + 0.5*d2vsdp2HEM;

        gs += 6000.0; /* Tajcmanova correction */
        hs += 6000.0;

    } else if (strcmp(name, "tbiTaj") == 0) {
        double gsPHL, ssPHL, hsPHL, cpsPHL, dcpsdtPHL, vsPHL, dvsdtPHL, dvsdpPHL, d2vsdt2PHL, d2vsdtdpPHL, d2vsdp2PHL;
        double gsBRU, ssBRU, hsBRU, cpsBRU, dcpsdtBRU, vsBRU, dvsdtBRU, dvsdpBRU, d2vsdt2BRU, d2vsdtdpBRU, d2vsdp2BRU;
        double gsRUT, ssRUT, hsRUT, cpsRUT, dcpsdtRUT, vsRUT, dvsdtRUT, dvsdpRUT, d2vsdt2RUT, d2vsdtdpRUT, d2vsdp2RUT;

        HPproperties(t, p, -6219160.0, 328.0, 14.964, 770.3, -0.036939, -2328900.0, -6531.6, 5.79e-5, 513000.0, 0.0, 0.0, 0.0,
                     &gsPHL, &ssPHL, &hsPHL, &cpsPHL, &dcpsdtPHL, &vsPHL, &dvsdtPHL, &dvsdpPHL, &d2vsdt2PHL, &d2vsdtdpPHL,
                     &d2vsdp2PHL);
        HPproperties(t, p,  -924940.0, 64.50,  2.463, 158.4, -0.004076, -1052300.0, -1171.3, 13.0e-5, 485000.0, 0.0, 0.0, 0.0,
                     &gsBRU, &ssBRU, &hsBRU, &cpsBRU, &dcpsdtBRU, &vsBRU, &dvsdtBRU, &dvsdpBRU, &d2vsdt2BRU, &d2vsdtdpBRU,
                     &d2vsdp2BRU);
        HPproperties(t, p,  -944180.0, 50.60,  1.882,  90.4,  0.002900,        0.0,  -623.8, 4.43e-5, 2225000.0, 0.0, 0.0, 0.0,
                     &gsRUT, &ssRUT, &hsRUT, &cpsRUT, &dcpsdtRUT, &vsRUT, &dvsdtRUT, &dvsdpRUT, &d2vsdt2RUT, &d2vsdtdpRUT,
                     &d2vsdp2RUT);
        gs       =  gsPHL	     - gsBRU	   + gsRUT;
        ss       =  ssPHL	     - ssBRU	   + ssRUT;
        hs       =  hsPHL	     - hsBRU	   + hsRUT;
        cps      =  cpsPHL      - cpsBRU	   + cpsRUT;
        dcpsdt   =  dcpsdtPHL   - dcpsdtBRU   + dcpsdtRUT;
        vs       =  vsPHL	     - vsBRU	   + vsRUT;
        dvsdt    =  dvsdtPHL    - dvsdtBRU    + dvsdtRUT;
        dvsdp    =  dvsdpPHL    - dvsdpBRU    + dvsdpRUT;
        d2vsdt2  =  d2vsdt2PHL  - d2vsdt2BRU  + d2vsdt2RUT;
        d2vsdtdp =  d2vsdtdpPHL - d2vsdtdpBRU + d2vsdtdpRUT;
        d2vsdp2  =  d2vsdp2PHL  - d2vsdp2BRU  + d2vsdp2RUT;

        gs += 84000.0 - t*11.5; /* Tajcmanova correction */
        hs += 84000.0;
        ss += 11.5;

    } else if (strcmp(name, "eastTaj") == 0) {
        /* H 298       S 298  V 298   a      b          c           d        a0       K 298     Tc   Smax Vmax */
        HPproperties(t, p, -6338170.0, 318.0, 14.738, 785.5, -0.038031, -2130300.0, -6893.7, 5.79e-5, 513000.0, 0.0, 0.0, 0.0,
                     &gs, &ss, &hs, &cps, &dcpsdt, &vs, &dvsdt, &dvsdp, &d2vsdt2, &d2vsdtdp, &d2vsdp2);
    } else if (strcmp(name, "annTaj") == 0) {
        /* H 298       S 298  V 298   a      b          c           d        a0       K 298     Tc   Smax Vmax */
        HPproperties(t, p, -5151670.0, 418.0, 15.432, 815.7, -0.034861,    19800.0, -7466.7, 5.79e-5, 513000.0, 0.0, 0.0, 0.0,
                     &gs, &ss, &hs, &cps, &dcpsdt, &vs, &dvsdt, &dvsdp, &d2vsdt2, &d2vsdtdp, &d2vsdp2);
    } else if (strcmp(name, "phlTaj") == 0) {
        /* H 298       S 298  V 298   a      b          c           d        a0       K 298     Tc   Smax Vmax */
        HPproperties(t, p, -6219160.0, 328.0, 14.964, 770.3, -0.036939, -2328900.0, -6531.6, 5.79e-5, 513000.0, 0.0, 0.0, 0.0,
                     &gs, &ss, &hs, &cps, &dcpsdt, &vs, &dvsdt, &dvsdp, &d2vsdt2, &d2vsdtdp, &d2vsdp2);

        /* Solids: General equations from Berman (1988) or Saxena (1993) */

    } else {
        hs   = phase->h;
        ss   = phase->s;

        if (phase->cp_type == CP_BERMAN) {
            k0   = phase->cp.Berman.k0;
            k1   = phase->cp.Berman.k1;
            k2   = phase->cp.Berman.k2;
            k3   = phase->cp.Berman.k3;
            cp_t = phase->cp.Berman.Tt;
            cp_h = phase->cp.Berman.deltah;
            l1   = phase->cp.Berman.l1;
            l2   = phase->cp.Berman.l2;

            cps = k0 + k1/sqrt(t) + k2/SQUARE(t) + k3/CUBE(t);
            dcpsdt = - 0.5*k1/pow(t, (double) 1.5) - 2.0*k2/CUBE(t)
            - 3.0*k3/QUARTIC(t);

            hs = hs + k0*(t-tr) + 2.0*k1*(sqrt(t)-sqrt(tr))
            - k2*(1.0/t-1.0/tr) - 0.5*k3*(1.0/(t*t)-1.0/(tr*tr));
            ss = ss + k0*log(t/tr)
            - 2.0*k1*(1.0/sqrt(t)-1.0/sqrt(tr))
            - 0.5*k2*(1.0/(t*t)-1.0/(tr*tr))
            - (1.0/3.0)*k3*(1.0/(t*t*t)-1.0/(tr*tr*tr));
            if(cp_t != 0.0) {
                if(t > cp_t) {
                    hs = hs + cp_h
                    + 0.5*l1*l1*(cp_t*cp_t-tr*tr)
                    + (2.0/3.0)*l1*l2*(cp_t*cp_t*cp_t-tr*tr*tr)
                    + 0.25*l2*l2*(cp_t*cp_t*cp_t*cp_t-tr*tr*tr*tr);
                    ss = ss + cp_h/cp_t
                    + l1*l1*(cp_t-tr)
                    + l1*l2*(cp_t*cp_t-tr*tr)
                    + (1.0/3.0)*l2*l2*(cp_t*cp_t*cp_t-tr*tr*tr);
                } else {
                    hs = hs + 0.5*l1*l1*(t*t-tr*tr)
                    + (2.0/3.0)*l1*l2*(t*t*t-tr*tr*tr)
                    + 0.25*l2*l2*(t*t*t*t-tr*tr*tr*tr);
                    ss = ss + l1*l1*(t-tr)
                    + l1*l2*(t*t-tr*tr)
                    + (1.0/3.0)*l2*l2*(t*t*t-tr*tr*tr);
                    cps += t*SQUARE(l1+l2*t);
                    dcpsdt += SQUARE(l1+l2*t) + t*2.0*(l1+l2*t)*l2;
                }
            }

        } else if (phase->cp_type == CP_SAXENA) {
            double aCp = phase->cp.Saxena.a;
            double bCp = phase->cp.Saxena.b;
            double cCp = phase->cp.Saxena.c;
            double dCp = phase->cp.Saxena.d;
            double eCp = phase->cp.Saxena.e;
            double gCp = phase->cp.Saxena.g;
            double hCp = phase->cp.Saxena.h;

            cps    = aCp + bCp*t + cCp/SQUARE(t) + dCp*SQUARE(t) + eCp/CUBE(t) + gCp/sqrt(t) + hCp/t;
            hs     = hs + aCp*(t-tr) + (bCp/2.0)*(SQUARE(t)-SQUARE(tr)) - cCp*(1.0/t-1.0/tr) + (dCp/3.0)*(CUBE(t)-CUBE(tr))
            - (eCp/2.0)*(1.0/SQUARE(t)-1.0/SQUARE(tr)) + 2.0*gCp*(sqrt(t)-sqrt(tr)) + hCp*log(t/tr);
            ss     = ss + aCp*log(t/tr) + bCp*(t-tr) - (cCp/2.0)*(1.0/SQUARE(t)-1.0/SQUARE(tr)) + (dCp/2.0)*(SQUARE(t)-SQUARE(tr))
            - (eCp/3.0)*(1.0/CUBE(t)-1.0/CUBE(tr)) - 2.0*gCp*(1.0/sqrt(t)-1.0/sqrt(tr)) - hCp*(1.0/t-1.0/tr);
            dcpsdt = bCp - 2.0*cCp/CUBE(t) + 2.0*dCp*t - 3.0*eCp/QUARTIC(t) - (gCp/2.0)/pow(t, (double) 3.0/2.0) - hCp/SQUARE(t);
        }

        gs       = hs - t*ss;
        intEOSsolid(phase, t, p, &gs, &hs, &ss, &cps, &dcpsdt, &vs, &dvsdt, &dvsdp, &d2vsdt2, &d2vsdtdp, &d2vsdp2);

        /* Special terms for minerals that involve order-disorder functions
         or additional first-order phase transitions                         */

        if(strcmp(name, "albite") == 0) {
#ifdef HIGH_STRUCTURAL_STATE_FELDSPAR
#else
            double gDis, hDis, sDis, cpDis, dcpdtDis, vDis, dvdtDis, dvdpDis,
            d2vdt2Dis, d2vdtdpDis, d2vdp2Dis;

            albite(p, t, &gDis, &hDis, &sDis, &cpDis, &dcpdtDis, &vDis,
                   &dvdtDis, &dvdpDis, &d2vdt2Dis, &d2vdtdpDis, &d2vdp2Dis);
            gs       += gDis;
            hs       += hDis;
            ss       += sDis;
            cps      += cpDis;
            dcpsdt   += dcpdtDis;
            vs       += vDis;
            dvsdt    += dvdtDis;
            dvsdp    += dvdpDis;
            d2vsdt2  += d2vdt2Dis;
            d2vsdtdp += d2vdtdpDis;
            d2vsdp2  += d2vdp2Dis;
#endif
        } else if(strcmp(name, "gehlenite") == 0) {

            double td, d0, d1, d2, d3, d4, d5, dhdis, dsdis, dvdis;

            td = MIN(1600.0,t);
            d0 = -221.74;
            d1 = 0.0;
            d2 = 172.91e5;
            d3 = 36.950e-2;
            d4 = -146.900e-6;
            d5 = 0.0;

            dhdis = d0*(td-698.0) + 2.0*d1*(sqrt(td)-sqrt(698.0))
            - d2*(1.0/td-1.0/698.0) + d3*(td*td-SQUARE(698.0))/2.0
            + d4*(CUBE(td)-CUBE(698.0))/3.0;
            dsdis = d0*(log(td)-log(698.0))
            - 2.0*d1*(1.0/sqrt(td)-1.0/sqrt(698.0))
            - d2*(1.0/SQUARE(td)-1.0/SQUARE(698.0))/2.0
            + d3*(td-698.0) + d4*(SQUARE(td)-SQUARE(698.0))/2.0;
            dvdis = (d5 != 0.0) ? dhdis/d5 : 0.0;

            gs       += (t > 698.0) ? dhdis - t*dsdis + dvdis*(p-pr) : 0.0;
            hs       += (t > 698.0) ? dhdis + dvdis*(p-pr) : 0.0;
            ss       += (t > 698.0) ? dsdis : 0.0;
            vs       += (t > 698.0) ? dvdis : 0.0;

            cps      += (t > 698.0 && t < 1600.0) ?
            d0+d1/sqrt(t)+d2/SQUARE(t)+d3*t+d4*SQUARE(t) : 0.0;
            dcpsdt   += (t > 698.0 && t < 1600.0) ?
            -d1*0.5/pow(t, (double) 1.5)-2.0*d2/CUBE(t)+d3+2.0*d4*t :
            0.0;
            dvsdt    += (t > 698.0 && t < 1600.0 && d5 != 0.0) ?
            (d0+d1/sqrt(t)+d2/SQUARE(t)+d3*t+d4*SQUARE(t))/d5 : 0.0;
            d2vsdt2  += (t > 698.0 && t < 1600.0 && d5 != 0.0) ?
            (-d1*0.5/pow(t, (double) 1.5)-2.0*d2/CUBE(t)+d3+2.0*d4*t)/
            d5 : 0.0;
            hs       += (t > tr && t < 1436.0 && d5 != 0.0) ? -(p-pr)*t*
            (d0+d1/sqrt(t)+d2/SQUARE(t)+d3*t+d4*SQUARE(t))/d5 : 0.0;
            ss       += (t > tr && t < 1436.0 && d5 != 0.0) ? -(p-pr)*
            (d0+d1/sqrt(t)+d2/SQUARE(t)+d3*t+d4*SQUARE(t))/d5 : 0.0;
            cps      += (t > tr && t < 1436.0 && d5 != 0.0) ? -(p-pr)*t*
            (-d1*0.5/pow(t, (double) 1.5)-2.0*d2/CUBE(t)+d3+2.0*d4*t)/
            d5 : 0.0;
            dcpsdt   += (t > tr && t < 1436.0 && d5 != 0.0) ? -(p-pr)*
            (-d1*0.5/pow(t, (double) 1.5)-2.0*d2/CUBE(t)+d3+2.0*d4*t)/
            d5
            -(p-pr)*t*(1.5*0.5*d1/pow(t,(double) 2.5)+6.0*d2/QUARTIC(t)
                       +2.0*d4)/d5 : 0.0;

        } else if(strcmp(name, "vc-nepheline") == 0) {
            ThermoRef  tempRef = {    /* prop are na-nepheline, unless noted */
                -2093004.0*4.0,                                      /* H ref (J)     */
                124.641*4.0,                                         /* S ref (J/K)   */
                5.434181*8.0,          /* V ref (J/bar) This is the vacancy endmember */
                /* 8 times because of factor of two below      */
                CP_BERMAN,  {{205.24*4.0, -7.599E2*4.0, -108.383E5*4.0, 208.182E7*4.0,
                    467.0, 241.0*4.0, -50.249E-2*2.0, 165.95E-5*2.0}},
                EOS_VINET,  {{31.802e-6, 48.7805, 1.4747}}
            };
            ThermoData tempRes;

            /* Get properties of Na4Al4Si4O16 at p and t */
            gibbs(t, p, "na-nepheline", &tempRef, NULL, NULL, &tempRes);

            /* [2*NaAlSi3O8 (high albite) + Na4Al4Si4O16 (beta-nepheline)]/2
             gs, hs, etc already contain properties of high-albite (Berman 1988)
             Note that the volumetric properties of high-albite have been zeroed */
            gs       = (2.0*gs       + tempRes.g)/2.0;
            hs       = (2.0*hs       + tempRes.h)/2.0;
            ss       = (2.0*ss       + tempRes.s)/2.0;
            vs       = (2.0*vs       + tempRes.v)/2.0;
            cps      = (2.0*cps      + tempRes.cp)/2.0;
            dcpsdt   = (2.0*dcpsdt   + tempRes.dcpdt)/2.0;
            dvsdt    = (2.0*dvsdt    + tempRes.dvdt)/2.0;
            dvsdp    = (2.0*dvsdp    + tempRes.dvdp)/2.0;
            d2vsdt2  = (2.0*d2vsdt2  + tempRes.d2vdt2)/2.0;
            d2vsdp2  = (2.0*d2vsdp2  + tempRes.d2vdp2)/2.0;
            d2vsdtdp = (2.0*d2vsdtdp + tempRes.d2vdtdp)/2.0;

        } else if(strcmp(name, "ca-nepheline") == 0) {
            ThermoRef  tempRef = {    /* prop are na-nepheline, unless noted */
                -2093004.0*4.0,                                      /* H ref (J)     */
                124.641*4.0,                                         /* S ref (J/K)   */
                5.433181*8.0,          /* V ref (J/bar) This is the vacancy endmember */
                /* 8 times because of factor of two below      */
                CP_BERMAN,  {{205.24*4.0, -7.599E2*4.0, -108.383E5*4.0, 208.182E7*4.0,
                    467.0, 241.0*4.0, -50.249E-2*2.0, 165.95E-5*2.0}},
                EOS_VINET,  {{31.802e-6, 48.7805, 1.4747}}
            };
            ThermoData tempRes;

            /* Get properties of Na4Al4Si4O16 at p and t */
            gibbs(t, p, "na-nepheline", &tempRef, NULL, NULL, &tempRes);

            /* Enthalpy correction */
            hs += 23096.0;
            gs += 23096.0;

            /* zero point correction for []CaNa2Al4Si4O16 */
            ss +=    15.8765;
            gs += -t*15.8765;

            /* [2*CaAl2Si2O8 (high anorthite) + Na4Al4Si4O16 (beta-nepheline)]/2
             gs, hs, etc already contain properties of high-anorthite (Berman 1988)
             Note that the volumetric properties of anorthite have been zeroed.   */
            gs       = (2.0*gs       + tempRes.g)/2.0;
            hs       = (2.0*hs       + tempRes.h)/2.0;
            ss       = (2.0*ss       + tempRes.s)/2.0;
            vs       = (2.0*vs       + tempRes.v)/2.0;
            cps      = (2.0*cps      + tempRes.cp)/2.0;
            dcpsdt   = (2.0*dcpsdt   + tempRes.dcpdt)/2.0;
            dvsdt    = (2.0*dvsdt    + tempRes.dvdt)/2.0;
            dvsdp    = (2.0*dvsdp    + tempRes.dvdp)/2.0;
            d2vsdt2  = (2.0*d2vsdt2  + tempRes.d2vdt2)/2.0;
            d2vsdp2  = (2.0*d2vsdp2  + tempRes.d2vdp2)/2.0;
            d2vsdtdp = (2.0*d2vsdtdp + tempRes.d2vdtdp)/2.0;

        } else if(strcmp(name, "iron-akermanite") == 0) {
            ThermoRef faRef = {
                -1479360.0,       /* H ref (J)                        Berman (1988) */
                150.930,          /* S ref (J/K)                      Berman (1988) */
                4.630,            /* V ref (J/bar)                    Berman (1988) */
                CP_BERMAN,  {{248.93, -19.239E2, 0.0, -13.910E7, 0.0, 0.0, 0.0, 0.0}},
                EOS_BERMAN, {{-0.730E-6, 0.0, 26.546E-6, 79.482E-10}}
            };
            ThermoRef foRef = {
                -2174420.0,       /* H ref (J)                        Berman (1988) */
                94.010,           /* S ref (J/K)                      Berman (1988) */
                4.366,            /* V ref (J/bar)                    Berman (1988) */
                CP_BERMAN,  {{238.64, -20.013E2, 0.0, -11.624E7, 0.0, 0.0, 0.0, 0.0}},
                EOS_BERMAN, {{-0.791E-6, 1.351E-12, 29.464E-6, 88.633E-10}}
            };

            ThermoData tempRes;

            gibbs(t, p, "fayalite", &faRef, NULL, NULL, &tempRes);

            /* Ca2MgSi2O7 + 1/2 Fe2SiO4 - 1/2 Mg2SiO4 in Gibbs.c */
            gs       = gs       + tempRes.g/2.0;
            hs       = hs       + tempRes.h/2.0;
            ss       = ss       + tempRes.s/2.0;
            vs       = vs       + tempRes.v/2.0;
            cps      = cps      + tempRes.cp/2.0;
            dcpsdt   = dcpsdt   + tempRes.dcpdt/2.0;
            dvsdt    = dvsdt    + tempRes.dvdt/2.0;
            dvsdp    = dvsdp    + tempRes.dvdp/2.0;
            d2vsdt2  = d2vsdt2  + tempRes.d2vdt2/2.0;
            d2vsdp2  = d2vsdp2  + tempRes.d2vdp2/2.0;
            d2vsdtdp = d2vsdtdp + tempRes.d2vdtdp/2.0;

            gibbs(t, p, "forsterite", &foRef, NULL, NULL, &tempRes);

            /* Ca2MgSi2O7 + 1/2 Fe2SiO4 - 1/2 Mg2SiO4 in Gibbs.c */
            gs       = gs       - tempRes.g/2.0;
            hs       = hs       - tempRes.h/2.0;
            ss       = ss       - tempRes.s/2.0;
            vs       = vs       - tempRes.v/2.0;
            cps      = cps      - tempRes.cp/2.0;
            dcpsdt   = dcpsdt   - tempRes.dcpdt/2.0;
            dvsdt    = dvsdt    - tempRes.dvdt/2.0;
            dvsdp    = dvsdp    - tempRes.dvdp/2.0;
            d2vsdt2  = d2vsdt2  - tempRes.d2vdt2/2.0;
            d2vsdp2  = d2vsdp2  - tempRes.d2vdp2/2.0;
            d2vsdtdp = d2vsdtdp - tempRes.d2vdtdp/2.0;

        } else if(strcmp(name, "soda-melilite") == 0) {
            ThermoRef abRef = {
                -3921618.0,      /* H ref (J)             monalbite - Berman (1988) */
                224.412,         /* S ref (J/K)                       Berman (1988) */
                10.083,          /* V ref (J/bar)                     Berman (1988) */
                CP_BERMAN,  {{393.64, -24.155E2, -78.928E5, 107.064E7, 0.0, 0.0, 0.0, 0.0}},
                EOS_BERMAN, {{-1.945E-6, 4.861E-12, 26.307E-6, 32.407E-10}}
            };
            ThermoRef anRef = {
                -4228730.0+3.7*4184.0,     /* H ref (J)   Berman (1988)  Carpenter I1->C1 */
                200.186+3.7*4184.0/2200.0, /* S ref (J/K) Berman (1988)  Carpenter I1->C1 */
                10.075,           /* V ref (J/bar)                          Berman (1988) */
                CP_BERMAN,  {{439.37, -37.341E2, 0.0, -31.702E7, 0.0, 0.0, 0.0, 0.0}},
                EOS_BERMAN, {{-1.272E-6, 3.176E-12, 10.918E-6, 41.985E-10}}
            };

            ThermoData tempRes;

            gibbs(t, p, "albite", &abRef, NULL, NULL, &tempRes);

            /* Ca2Al2SiO7 + 2NaAlSi3O8 - 2CaAl2Si2O8 in Gibbs.c */
            gs       = gs       + 2.0*tempRes.g;
            hs       = hs       + 2.0*tempRes.h;
            ss       = ss       + 2.0*tempRes.s;
            vs       = vs       + 2.0*tempRes.v;
            cps      = cps      + 2.0*tempRes.cp;
            dcpsdt   = dcpsdt   + 2.0*tempRes.dcpdt;
            dvsdt    = dvsdt    + 2.0*tempRes.dvdt;
            dvsdp    = dvsdp    + 2.0*tempRes.dvdp;
            d2vsdt2  = d2vsdt2  + 2.0*tempRes.d2vdt2;
            d2vsdp2  = d2vsdp2  + 2.0*tempRes.d2vdp2;
            d2vsdtdp = d2vsdtdp + 2.0*tempRes.d2vdtdp;

            gibbs(t, p, "anorthite", &anRef, NULL, NULL, &tempRes);

            /* Ca2Al2SiO7 + 2NaAlSi3O8 - 2CaAl2Si2O8 in Gibbs.c */
            gs       = gs       - 2.0*tempRes.g;
            hs       = hs       - 2.0*tempRes.h;
            ss       = ss       - 2.0*tempRes.s;
            vs       = vs       - 2.0*tempRes.v;
            cps      = cps      - 2.0*tempRes.cp;
            dcpsdt   = dcpsdt   - 2.0*tempRes.dcpdt;
            dvsdt    = dvsdt    - 2.0*tempRes.dvdt;
            dvsdp    = dvsdp    - 2.0*tempRes.dvdp;
            d2vsdt2  = d2vsdt2  - 2.0*tempRes.d2vdt2;
            d2vsdp2  = d2vsdp2  - 2.0*tempRes.d2vdp2;
            d2vsdtdp = d2vsdtdp - 2.0*tempRes.d2vdtdp;

        } else if((strcmp(name, "nepheline") == 0) && (t >= 1180.15) ) {
            hs += 2393.0;
            ss += 2393.0/1180.15;
            gs += 2393.0 - t*2393.0/1180.15;
        } else if(strcmp(name, "sanidine") == 0) {
            double td, d0, d1, d2, d3, d4, d5, dhdis, dsdis, dvdis;

#ifdef HIGH_STRUCTURAL_STATE_FELDSPAR
            td = 1436.0;
#else
            td = MIN(1436.0,t);
#endif
            d0 = 282.98;
            d1 = -4.83e3;
            d2 = 36.21e5;
            d3 = -15.733e-2;
            d4 = 34.770e-6;
            d5 = 41.063e4;

            dhdis = d0*(td-tr) + 2.0*d1*(sqrt(td)-sqrt(tr)) - d2*(1.0/td-1.0/tr)
            + d3*(td*td-SQUARE(tr))/2.0 + d4*(CUBE(td)-CUBE(tr))/3.0;
            dsdis = d0*(log(td)-log(tr)) - 2.0*d1*(1.0/sqrt(td)-1.0/sqrt(tr))
            - d2*(1.0/SQUARE(td)-1.0/SQUARE(tr))/2.0 + d3*(td-tr)
            + d4*(SQUARE(td)-SQUARE(tr))/2.0;
            dvdis = (d5 != 0.0) ? dhdis/d5 : 0.0;

            gs       += (t > tr) ? dhdis - t*dsdis + dvdis*(p-pr) : 0.0;
            hs       += (t > tr) ? dhdis + dvdis*(p-pr) : 0.0;
            ss       += (t > tr) ? dsdis : 0.0;
            vs       += (t > tr) ? dvdis : 0.0;

            cps      += (t > tr && t < 1436.0) ?
            d0+d1/sqrt(t)+d2/SQUARE(t)+d3*t+d4*SQUARE(t) : 0.0;
            dcpsdt   += (t > tr && t < 1436.0) ?
            -d1*0.5/pow(t, (double) 1.5)-2.0*d2/CUBE(t)+d3+2.0*d4*t :
            0.0;
            dvsdt    += (t > tr && t < 1436.0 && d5 != 0.0) ?
            (d0+d1/sqrt(t)+d2/SQUARE(t)+d3*t+d4*SQUARE(t))/d5 : 0.0;
            d2vsdt2  += (t > tr && t < 1436.0 && d5 != 0.0) ?
            (-d1*0.5/pow(t, (double) 1.5)-2.0*d2/CUBE(t)+d3+2.0*d4*t)/
            d5 : 0.0;
            hs       += (t > tr && t < 1436.0 && d5 != 0.0) ? -(p-pr)*t*
            (d0+d1/sqrt(t)+d2/SQUARE(t)+d3*t+d4*SQUARE(t))/d5 : 0.0;
            ss       += (t > tr && t < 1436.0 && d5 != 0.0) ? -(p-pr)*
            (d0+d1/sqrt(t)+d2/SQUARE(t)+d3*t+d4*SQUARE(t))/d5 : 0.0;
            cps      += (t > tr && t < 1436.0 && d5 != 0.0) ? -(p-pr)*t*
            (-d1*0.5/pow(t, (double) 1.5)-2.0*d2/CUBE(t)+d3+2.0*d4*t)/
            d5 : 0.0;
            dcpsdt   += (t > tr && t < 1436.0 && d5 != 0.0) ? -(p-pr)*
            (-d1*0.5/pow(t, (double) 1.5)-2.0*d2/CUBE(t)+d3+2.0*d4*t)/
            d5
            -(p-pr)*t*(1.5*0.5*d1/pow(t,(double) 2.5)+6.0*d2/QUARTIC(t)
                       +2.0*d4)/d5
            : 0.0;
        }
    }

    result->g       = gs;
    result->h       = hs;
    result->s       = ss;
    result->v       = vs;
    result->cp      = cps;
    result->dcpdt   = dcpsdt;
    result->dvdt    = dvsdt;
    result->dvdp    = dvsdp;
    result->d2vdt2  = d2vsdt2;
    result->d2vdp2  = d2vsdp2;
    result->d2vdtdp = d2vsdtdp;
}

#ifdef SPECIAL_O2
static void getO2properties(double t, double p, double *g, double *h,
                            double *s, double *cp, double *dcpdt, double *v, double *dvdt, double *dvdp,
                            double *d2vdt2, double *d2vdp2, double *d2vdtdp)
{
    int indFeO = -1, indFe2O3 = -1;
    int i;

    if (indFeO == -1 || indFe2O3 == -1) {
        for (i=0; i<nc; i++) {
            if (bulkSystem[i].type == FEO)   indFeO   = i;
            if (bulkSystem[i].type == FE2O3) indFe2O3 = i;
        }
        if (indFeO == -1 || indFe2O3 == -1) {
            printf("Internal Error in function gibbs. Cannot find FeO or Fe2O3\n");
            exit(0);
        }
    }

    for (i=0, *g=0.0, *h=0.0, *s=0.0, *cp=0.0, *dcpdt=0.0, *v=0.0, *dvdt=0.0,
         *dvdp=0.0, *d2vdt2=0.0, *d2vdp2=0.0, *d2vdtdp=0.0; i<nlc; i++) {
        double coeff = 2.0*(bulkSystem[indFe2O3].oxToLiq)[i]
        - 4.0*(bulkSystem[indFeO].oxToLiq)[i];
        if( coeff != 0.0) {
            gibbs(t, p, (char *) liquid[i].label, &(liquid[i].ref), &(liquid[i].liq),
                  &(liquid[i].fus), &(liquid[i].cur));
            *g       += coeff*(liquid[i].cur).g;
            *h       += coeff*(liquid[i].cur).h;
            *s       += coeff*(liquid[i].cur).s;
            *cp      += coeff*(liquid[i].cur).cp;
            *dcpdt   += coeff*(liquid[i].cur).dcpdt;
            *v       += coeff*(liquid[i].cur).v;
            *dvdt    += coeff*(liquid[i].cur).dvdt;
            *dvdp    += coeff*(liquid[i].cur).dvdp;
            *d2vdt2  += coeff*(liquid[i].cur).d2vdt2;
            *d2vdp2  += coeff*(liquid[i].cur).d2vdp2;
            *d2vdtdp += coeff*(liquid[i].cur).d2vdtdp;
        }
    }
}
#endif

static void getCaCO3properties(double t, double p, double *g, double *h, double *s, double *cp, double *dcpdt,
                               double *v, double *dvdt, double *dvdp, double *d2vdt2, double *d2vdtdp, double *d2vdp2) {
    int nSiO2 = -1, nCaSiO3 = -1, nCO2 = -1;
    ThermoData result;

    if ( (nSiO2 == -1) || (nCaSiO3 == -1) || (nCO2 == -1) ) {
        int i;
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "SiO2")   == 0) { nSiO2   = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "CaSiO3") == 0) { nCaSiO3 = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "CO2")    == 0) { nCO2    = i; break; }
        if ( (nSiO2 == -1) || (nCaSiO3 == -1) || (nCO2 == -1)) {
            printf("FATAL ERROR in Gibbs.c. Request for CaCO3 properties when nSiO2 = %d, nCaSiO3 = %d and nCO2 = %d.\n", nSiO2, nCaSiO3, nCO2);
            exit(0);
        }
    }

    /* reaction: CaSiO3 + CO2 - SiO2 = CaCO3 */

    gibbs(t, p, (char *) liquid[nCaSiO3].label, &(liquid[nCaSiO3].ref), &(liquid[nCaSiO3].liq), &(liquid[nCaSiO3].fus), &result);
    *g       = result.g;
    *h       = result.h;
    *s       = result.s;
    *cp      = result.cp;
    *dcpdt   = result.dcpdt;
    *v       = result.v;
    *dvdt    = result.dvdt;
    *dvdp    = result.dvdp;
    *d2vdp2  = result.d2vdt2;
    *d2vdtdp = result.d2vdtdp;
    *d2vdp2  = result.d2vdp2;

    gibbs(t, p, (char *) liquid[nCO2].label, &(liquid[nCO2].ref), &(liquid[nCO2].liq), &(liquid[nCO2].fus), &result);
    *g       += result.g;
    *h       += result.h;
    *s       += result.s;
    *cp      += result.cp;
    *dcpdt   += result.dcpdt;
    *v       += result.v;
    *dvdt    += result.dvdt;
    *dvdp    += result.dvdp;
    *d2vdp2  += result.d2vdt2;
    *d2vdtdp += result.d2vdtdp;
    *d2vdp2  += result.d2vdp2;

    gibbs(t, p, (char *) liquid[nSiO2].label, &(liquid[nSiO2].ref), &(liquid[nSiO2].liq), &(liquid[nSiO2].fus), &result);
    *g       -= result.g;
    *h       -= result.h;
    *s       -= result.s;
    *cp      -= result.cp;
    *dcpdt   -= result.dcpdt;
    *v       -= result.v;
    *dvdt    -= result.dvdt;
    *dvdp    -= result.dvdp;
    *d2vdp2  -= result.d2vdt2;
    *d2vdtdp -= result.d2vdtdp;
    *d2vdp2  -= result.d2vdp2;
}

static void getSiOHproperties(double t, double *g, double *h, double *s, double *cp, double *dcpdt) {
    int nSiO2 = -1, nH2O = -1;
    ThermoData result;

    if ( (nSiO2 == -1) || (nH2O == -1) ) {
        int i;
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "SiO2") == 0) { nSiO2 = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "H2O")  == 0) { nH2O  = i; break; }
        if ( (nSiO2 == -1) || (nH2O == -1) ) {
            printf("FATAL ERROR in Gibbs.c. Request for Si0.25OH properties when nSiO2 = %d and nH2O = %d.\n", nSiO2, nH2O);
            exit(0);
        }
    }

    /* reaction: (1/4) SiO2 + 1/2 H2O = Si(1/4)OH */

    gibbs(t, 1.0, (char *) liquid[nSiO2].label, &(liquid[nSiO2].ref), &(liquid[nSiO2].liq), &(liquid[nSiO2].fus), &result);
    *g     += 0.25*result.g;
    *h     += 0.25*result.h;
    *s     += 0.25*result.s;
    *cp    += 0.25*result.cp;
    *dcpdt += 0.25*result.dcpdt;

    gibbs(t, 1.0, (char *) liquid[nH2O].label, &(liquid[nH2O].ref), &(liquid[nH2O].liq), &(liquid[nH2O].fus), &result);
    *g     += 0.5*result.g;
    *h     += 0.5*result.h;
    *s     += 0.5*result.s;
    *cp    += 0.5*result.cp;
    *dcpdt += 0.5*result.dcpdt;
}

static void getAlOHproperties(double t, double *g, double *h, double *s, double *cp, double *dcpdt) {
    int nAl2O3 = -1, nH2O = -1;
    ThermoData result;
    const double FRAC = 1.0/6.0;

    if ( (nAl2O3 == -1) || (nH2O == -1) ) {
        int i;
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "Al2O3") == 0) { nAl2O3 = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "H2O")   == 0) { nH2O   = i; break; }
        if ( (nAl2O3 == -1) || (nH2O == -1) ) {
            printf("FATAL ERROR in Gibbs.c. Request for Al0.33OH properties when nAl2O3 = %d and nH2O = %d.\n", nAl2O3, nH2O);
            exit(0);
        }
    }

    /* reaction: (1/6) Al2O3 + 1/2 H2O = Al(1/3)OH */

    gibbs(t, 1.0, (char *) liquid[nAl2O3].label, &(liquid[nAl2O3].ref), &(liquid[nAl2O3].liq), &(liquid[nAl2O3].fus), &result);
    *g     += FRAC*result.g;
    *h     += FRAC*result.h;
    *s     += FRAC*result.s;
    *cp    += FRAC*result.cp;
    *dcpdt += FRAC*result.dcpdt;

    gibbs(t, 1.0, (char *) liquid[nH2O].label, &(liquid[nH2O].ref), &(liquid[nH2O].liq), &(liquid[nH2O].fus), &result);
    *g     += 0.5*result.g;
    *h     += 0.5*result.h;
    *s     += 0.5*result.s;
    *cp    += 0.5*result.cp;
    *dcpdt += 0.5*result.dcpdt;
}

static void getFe3OHproperties(double t, double *g, double *h, double *s, double *cp, double *dcpdt) {
    int nSiO2 = -1, nFe2SiO5 = -1, nH2O = -1;
    ThermoData result;
    const double FRAC = 1.0/6.0;

    if ( (nSiO2 == -1) || (nFe2SiO5 == -1) || (nH2O == -1) ) {
        int i;
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "SiO2")    == 0) { nSiO2    = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "Fe2SiO5") == 0) { nFe2SiO5 = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "H2O")     == 0) { nH2O     = i; break; }
        if ( (nSiO2 == -1) || (nFe2SiO5 == -1) || (nH2O == -1) ) {
            printf("FATAL ERROR in Gibbs.c. Request for Fe0.33OH properties when nSiO2 = %d, nFe2SiO5 = %d and nH2O = %d.\n", nSiO2, nFe2SiO5, nH2O);
            exit(0);
        }
    }

    /* reaction: (1/6) Fe2SiO5 + 1/2 H2O - (1/6) SiO2 = Fe(1/3)OH */

    gibbs(t, 1.0, (char *) liquid[nFe2SiO5].label, &(liquid[nFe2SiO5].ref), &(liquid[nFe2SiO5].liq), &(liquid[nFe2SiO5].fus), &result);
    *g     += FRAC*result.g;
    *h     += FRAC*result.h;
    *s     += FRAC*result.s;
    *cp    += FRAC*result.cp;
    *dcpdt += FRAC*result.dcpdt;

    gibbs(t, 1.0, (char *) liquid[nH2O].label, &(liquid[nH2O].ref), &(liquid[nH2O].liq), &(liquid[nH2O].fus), &result);
    *g     += 0.5*result.g;
    *h     += 0.5*result.h;
    *s     += 0.5*result.s;
    *cp    += 0.5*result.cp;
    *dcpdt += 0.5*result.dcpdt;

    gibbs(t, 1.0, (char *) liquid[nSiO2].label, &(liquid[nSiO2].ref), &(liquid[nSiO2].liq), &(liquid[nSiO2].fus), &result);
    *g     -= FRAC*result.g;
    *h     -= FRAC*result.h;
    *s     -= FRAC*result.s;
    *cp    -= FRAC*result.cp;
    *dcpdt -= FRAC*result.dcpdt;
}

static void getFe2OHproperties(double t, double *g, double *h, double *s, double *cp, double *dcpdt) {
    int nSiO2 = -1, nFe2SiO4 = -1, nH2O = -1;
    ThermoData result;
    const double FRAC = 1.0/4.0;

    if ( (nSiO2 == -1) || (nFe2SiO4 == -1) || (nH2O == -1) ) {
        int i;
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "SiO2")    == 0) { nSiO2    = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "Fe2SiO4") == 0) { nFe2SiO4 = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "H2O")     == 0) { nH2O     = i; break; }
        if ( (nSiO2 == -1) || (nFe2SiO4 == -1) || (nH2O == -1) ) {
            printf("FATAL ERROR in Gibbs.c. Request for Fe0.5OH properties when nSiO2 = %d, nFe2SiO4 = %d and nH2O = %d.\n", nSiO2, nFe2SiO4, nH2O);
            exit(0);
        }
    }

    /* reaction: (1/4) Fe2SiO4 + 1/2 H2O - (1/4) SiO2 = Fe(1/2)OH */

    gibbs(t, 1.0, (char *) liquid[nFe2SiO4].label, &(liquid[nFe2SiO4].ref), &(liquid[nFe2SiO4].liq), &(liquid[nFe2SiO4].fus), &result);
    *g     += FRAC*result.g;
    *h     += FRAC*result.h;
    *s     += FRAC*result.s;
    *cp    += FRAC*result.cp;
    *dcpdt += FRAC*result.dcpdt;

    gibbs(t, 1.0, (char *) liquid[nH2O].label, &(liquid[nH2O].ref), &(liquid[nH2O].liq), &(liquid[nH2O].fus), &result);
    *g     += 0.5*result.g;
    *h     += 0.5*result.h;
    *s     += 0.5*result.s;
    *cp    += 0.5*result.cp;
    *dcpdt += 0.5*result.dcpdt;

    gibbs(t, 1.0, (char *) liquid[nSiO2].label, &(liquid[nSiO2].ref), &(liquid[nSiO2].liq), &(liquid[nSiO2].fus), &result);
    *g     -= FRAC*result.g;
    *h     -= FRAC*result.h;
    *s     -= FRAC*result.s;
    *cp    -= FRAC*result.cp;
    *dcpdt -= FRAC*result.dcpdt;
}

static void getMgOHproperties(double t, double *g, double *h, double *s, double *cp, double *dcpdt) {
    int nSiO2 = -1, nMg2SiO4 = -1, nH2O = -1;
    ThermoData result;
    const double FRAC = 1.0/4.0;

    if ( (nSiO2 == -1) || (nMg2SiO4 == -1) || (nH2O == -1) ) {
        int i;
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "SiO2")    == 0) { nSiO2    = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "Mg2SiO4") == 0) { nMg2SiO4 = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "H2O")     == 0) { nH2O     = i; break; }
        if ( (nSiO2 == -1) || (nMg2SiO4 == -1) || (nH2O == -1) ) {
            printf("FATAL ERROR in Gibbs.c. Request for Mg0.5OH properties when nSiO2 = %d, nMg2SiO4 = %d and nH2O = %d.\n", nSiO2, nMg2SiO4, nH2O);
            exit(0);
        }
    }

    /* reaction: (1/4) Mg2SiO4 + 1/2 H2O - (1/4) SiO2 = Mg(1/2)OH */

    gibbs(t, 1.0, (char *) liquid[nMg2SiO4].label, &(liquid[nMg2SiO4].ref), &(liquid[nMg2SiO4].liq), &(liquid[nMg2SiO4].fus), &result);
    *g     += FRAC*result.g;
    *h     += FRAC*result.h;
    *s     += FRAC*result.s;
    *cp    += FRAC*result.cp;
    *dcpdt += FRAC*result.dcpdt;

    gibbs(t, 1.0, (char *) liquid[nH2O].label, &(liquid[nH2O].ref), &(liquid[nH2O].liq), &(liquid[nH2O].fus), &result);
    *g     += 0.5*result.g;
    *h     += 0.5*result.h;
    *s     += 0.5*result.s;
    *cp    += 0.5*result.cp;
    *dcpdt += 0.5*result.dcpdt;

    gibbs(t, 1.0, (char *) liquid[nSiO2].label, &(liquid[nSiO2].ref), &(liquid[nSiO2].liq), &(liquid[nSiO2].fus), &result);
    *g     -= FRAC*result.g;
    *h     -= FRAC*result.h;
    *s     -= FRAC*result.s;
    *cp    -= FRAC*result.cp;
    *dcpdt -= FRAC*result.dcpdt;
}

static void getCaOHproperties(double t, double *g, double *h, double *s, double *cp, double *dcpdt) {
    int nSiO2 = -1, nCa2SiO4 = -1, nH2O = -1;
    ThermoData result;
    const double FRAC = 1.0/4.0;

    if ( (nSiO2 == -1) || (nCa2SiO4 == -1) || (nH2O == -1) ) {
        int i;
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "SiO2")    == 0) { nSiO2    = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "Ca2SiO4") == 0) { nCa2SiO4 = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "H2O")     == 0) { nH2O     = i; break; }
        if ( (nSiO2 == -1) || (nCa2SiO4 == -1) || (nH2O == -1) ) {
            printf("FATAL ERROR in Gibbs.c. Request for Ca0.5OH properties when nSiO2 = %d, nCa2SiO4 = %d and nH2O = %d.\n", nSiO2, nCa2SiO4, nH2O);
            exit(0);
        }
    }

    /* reaction: (1/4) Ca2SiO4 + 1/2 H2O - (1/4) SiO2 = Ca(1/2)OH */

    gibbs(t, 1.0, (char *) liquid[nCa2SiO4].label, &(liquid[nCa2SiO4].ref), &(liquid[nCa2SiO4].liq), &(liquid[nCa2SiO4].fus), &result);
    *g     += FRAC*result.g;
    *h     += FRAC*result.h;
    *s     += FRAC*result.s;
    *cp    += FRAC*result.cp;
    *dcpdt += FRAC*result.dcpdt;

    gibbs(t, 1.0, (char *) liquid[nH2O].label, &(liquid[nH2O].ref), &(liquid[nH2O].liq), &(liquid[nH2O].fus), &result);
    *g     += 0.5*result.g;
    *h     += 0.5*result.h;
    *s     += 0.5*result.s;
    *cp    += 0.5*result.cp;
    *dcpdt += 0.5*result.dcpdt;

    gibbs(t, 1.0, (char *) liquid[nSiO2].label, &(liquid[nSiO2].ref), &(liquid[nSiO2].liq), &(liquid[nSiO2].fus), &result);
    *g     -= FRAC*result.g;
    *h     -= FRAC*result.h;
    *s     -= FRAC*result.s;
    *cp    -= FRAC*result.cp;
    *dcpdt -= FRAC*result.dcpdt;
}

static void getNaOHproperties(double t, double *g, double *h, double *s, double *cp, double *dcpdt) {
    int nSiO2 = -1, nNa2SiO3 = -1, nNa4SiO4 = -1, nH2O = -1;
    ThermoData result;
    const double FRAC2 = 1.0/2.0;
    const double FRAC4 = 1.0/4.0;

    if ( (nSiO2 == -1) || ((nNa2SiO3 == -1) && (nNa4SiO4 == -1)) || (nH2O == -1) ) {
        int i;
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "SiO2")    == 0) { nSiO2    = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "Na2SiO3") == 0) { nNa2SiO3 = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "Na4SiO4") == 0) { nNa4SiO4 = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "H2O")     == 0) { nH2O     = i; break; }
        if ( (nSiO2 == -1) || ((nNa2SiO3 == -1) && (nNa4SiO4 == -1)) || (nH2O == -1) ) {
            printf("FATAL ERROR in Gibbs.c. Request for NaOH properties when nSiO2 = %d, nNa2SiO3 = %d, nNa4SiO4 = %d and nH2O = %d.\n",
                   nSiO2, nNa2SiO3, nNa4SiO4, nH2O);
            exit(0);
        }
    }

    if (nNa2SiO3 != -1) {
        /* reaction: (1/2) Na2SiO3 + 1/2 H2O - (1/2) SiO2 = NaOH */

        gibbs(t, 1.0, (char *) liquid[nNa2SiO3].label, &(liquid[nNa2SiO3].ref), &(liquid[nNa2SiO3].liq), &(liquid[nNa2SiO3].fus), &result);
        *g     += FRAC2*result.g;
        *h     += FRAC2*result.h;
        *s     += FRAC2*result.s;
        *cp    += FRAC2*result.cp;
        *dcpdt += FRAC2*result.dcpdt;

        gibbs(t, 1.0, (char *) liquid[nH2O].label, &(liquid[nH2O].ref), &(liquid[nH2O].liq), &(liquid[nH2O].fus), &result);
        *g     += FRAC2*result.g;
        *h     += FRAC2*result.h;
        *s     += FRAC2*result.s;
        *cp    += FRAC2*result.cp;
        *dcpdt += FRAC2*result.dcpdt;

        gibbs(t, 1.0, (char *) liquid[nSiO2].label, &(liquid[nSiO2].ref), &(liquid[nSiO2].liq), &(liquid[nSiO2].fus), &result);
        *g     -= FRAC2*result.g;
        *h     -= FRAC2*result.h;
        *s     -= FRAC2*result.s;
        *cp    -= FRAC2*result.cp;
        *dcpdt -= FRAC2*result.dcpdt;
    } else {
        /* reaction: (1/4) Na4SiO4 + 1/2 H2O - (1/4) SiO2 = NaOH */

        gibbs(t, 1.0, (char *) liquid[nNa4SiO4].label, &(liquid[nNa4SiO4].ref), &(liquid[nNa4SiO4].liq), &(liquid[nNa4SiO4].fus), &result);
        *g     += FRAC4*result.g;
        *h     += FRAC4*result.h;
        *s     += FRAC4*result.s;
        *cp    += FRAC4*result.cp;
        *dcpdt += FRAC4*result.dcpdt;

        gibbs(t, 1.0, (char *) liquid[nH2O].label, &(liquid[nH2O].ref), &(liquid[nH2O].liq), &(liquid[nH2O].fus), &result);
        *g     += FRAC2*result.g;
        *h     += FRAC2*result.h;
        *s     += FRAC2*result.s;
        *cp    += FRAC2*result.cp;
        *dcpdt += FRAC2*result.dcpdt;

        gibbs(t, 1.0, (char *) liquid[nSiO2].label, &(liquid[nSiO2].ref), &(liquid[nSiO2].liq), &(liquid[nSiO2].fus), &result);
        *g     -= FRAC4*result.g;
        *h     -= FRAC4*result.h;
        *s     -= FRAC4*result.s;
        *cp    -= FRAC4*result.cp;
        *dcpdt -= FRAC4*result.dcpdt;
    }
}

static void getKOHproperties(double t, double *g, double *h, double *s, double *cp, double *dcpdt) {
    int nSiO2 = -1, nK2SiO3 = -1, nH2O = -1;
    ThermoData result;
    const double FRAC = 1.0/2.0;

    if ( (nSiO2 == -1) || (nK2SiO3 == -1) || (nH2O == -1) ) {
        int i;
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "SiO2")   == 0) { nSiO2    = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "K2SiO3") == 0) { nK2SiO3 = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "H2O")    == 0) { nH2O     = i; break; }
        if ( (nSiO2 == -1) || (nK2SiO3 == -1) || (nH2O == -1) ) {
            printf("FATAL ERROR in Gibbs.c. Request for KOH properties when nSiO2 = %d, nK2SiO3 = %d and nH2O = %d.\n", nSiO2, nK2SiO3, nH2O);
            exit(0);
        }
    }

    /* reaction: (1/2) K2SiO3 + 1/2 H2O - (1/2) SiO2 = KOH */

    gibbs(t, 1.0, (char *) liquid[nK2SiO3].label, &(liquid[nK2SiO3].ref), &(liquid[nK2SiO3].liq), &(liquid[nK2SiO3].fus), &result);
    *g     += FRAC*result.g;
    *h     += FRAC*result.h;
    *s     += FRAC*result.s;
    *cp    += FRAC*result.cp;
    *dcpdt += FRAC*result.dcpdt;

    gibbs(t, 1.0, (char *) liquid[nH2O].label, &(liquid[nH2O].ref), &(liquid[nH2O].liq), &(liquid[nH2O].fus), &result);
    *g     += 0.5*result.g;
    *h     += 0.5*result.h;
    *s     += 0.5*result.s;
    *cp    += 0.5*result.cp;
    *dcpdt += 0.5*result.dcpdt;

    gibbs(t, 1.0, (char *) liquid[nSiO2].label, &(liquid[nSiO2].ref), &(liquid[nSiO2].liq), &(liquid[nSiO2].fus), &result);
    *g     -= FRAC*result.g;
    *h     -= FRAC*result.h;
    *s     -= FRAC*result.s;
    *cp    -= FRAC*result.cp;
    *dcpdt -= FRAC*result.dcpdt;
}

static void getFe2AlO4_5properties(double t, double *g, double *h, double *s, double *cp, double *dcpdt) {
    int nSiO2 = -1, nAl2O3 = -1, nFe2SiO5 = -1;
    ThermoData result;

    if ( (nSiO2 == -1) || (nAl2O3 == -1) || (nFe2SiO5 == -1) ) {
        int i;
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "SiO2")	  == 0) { nSiO2    = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "Al2O3")    == 0) { nAl2O3   = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "Fe2SiO5")  == 0) { nFe2SiO5 = i; break; }
        if ( (nSiO2 == -1) || (nAl2O3 == -1) || (nFe2SiO5 == -1) ) {
            printf("FATAL ERROR in Gibbs.c. Request for Fe2AlO4.5 properties when nSiO2 = %d, nAl2O3 = %d, and nFe2SiO5 = %d.\n",
                   nSiO2, nAl2O3, nFe2SiO5);
            exit(0);
        }
    }

    /* reaction: Fe2SiO5 - SiO2 + 1/2 Al2O3 = Fe2AlSiO4.5 */

    gibbs(t, 1.0, (char *) liquid[nFe2SiO5].label, &(liquid[nFe2SiO5].ref), &(liquid[nFe2SiO5].liq),
          &(liquid[nFe2SiO5].fus), &result);
    *g	 += result.g;
    *h	 += result.h;
    *s	 += result.s;
    *cp	 += result.cp;
    *dcpdt += result.dcpdt;

    gibbs(t, 1.0, (char *) liquid[nSiO2].label, &(liquid[nSiO2].ref), &(liquid[nSiO2].liq),
          &(liquid[nSiO2].fus), &result);
    *g	 += -result.g;
    *h	 += -result.h;
    *s	 += -result.s;
    *cp	 += -result.cp;
    *dcpdt += -result.dcpdt;

    gibbs(t, 1.0, (char *) liquid[nAl2O3].label, &(liquid[nAl2O3].ref), &(liquid[nAl2O3].liq),
          &(liquid[nAl2O3].fus), &result);
    *g	 += 0.5*result.g;
    *h	 += 0.5*result.h;
    *s	 += 0.5*result.s;
    *cp	 += 0.5*result.cp;
    *dcpdt += 0.5*result.dcpdt;
}

static void getFe2AlO4_1properties(double t, double *g, double *h, double *s, double *cp, double *dcpdt) {
    int nSiO2 = -1, nAl2O3 = -1, nFe2SiO4_6 = -1;
    ThermoData result;

    if ( (nSiO2 == -1) || (nAl2O3 == -1) || (nFe2SiO4_6 == -1) ) {
        int i;
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "SiO2")	   == 0) { nSiO2      = i; break; }
        for (i=0; i<nlc; i++) if (strcmp(liquid[i].label, "Al2O3")     == 0) { nAl2O3     = i; break; }
        for (i=0; i<nls; i++) if (strcmp(liquid[i].label, "Fe2SiO4.6") == 0) { nFe2SiO4_6 = i; break; }
        if ( (nSiO2 == -1) || (nAl2O3 == -1) || (nFe2SiO4_6 == -1) ) {
            printf("FATAL ERROR in Gibbs.c. Request for Fe2AlO4.1 properties when nSiO2 = %d, nAl2O3 = %d, and nFe2SiO4_6 = %d.\n",
                   nSiO2, nAl2O3, nFe2SiO4_6);
            exit(0);
        }
    }

    /* reaction: Fe2SiO4.6 - SiO2 + 1/2 Al2O3 = Fe2AlSiO4.1 */

    gibbs(t, 1.0, (char *) liquid[nFe2SiO4_6].label, &(liquid[nFe2SiO4_6].ref), &(liquid[nFe2SiO4_6].liq),
          &(liquid[nFe2SiO4_6].fus), &result);
    *g	 += result.g;
    *h	 += result.h;
    *s	 += result.s;
    *cp	 += result.cp;
    *dcpdt += result.dcpdt;

    gibbs(t, 1.0, (char *) liquid[nSiO2].label, &(liquid[nSiO2].ref), &(liquid[nSiO2].liq),
          &(liquid[nSiO2].fus), &result);
    *g	 += -result.g;
    *h	 += -result.h;
    *s	 += -result.s;
    *cp	 += -result.cp;
    *dcpdt += -result.dcpdt;

    gibbs(t, 1.0, (char *) liquid[nAl2O3].label, &(liquid[nAl2O3].ref), &(liquid[nAl2O3].liq),
          &(liquid[nAl2O3].fus), &result);
    *g	 += 0.5*result.g;
    *h	 += 0.5*result.h;
    *s	 += 0.5*result.s;
    *cp	 += 0.5*result.cp;
    *dcpdt += 0.5*result.dcpdt;
}

static void fe_metal(double t, double p, double *gs, double *hs,
                     double *ss, double *cps, double *dcpsdt, double *vs, double *dvsdt,
                     double *dvsdp, double *d2vsdt2, double *d2vsdp2, double *d2vsdtdp)
{
    double tr =     298.15;
    double pr =       1.0;
    double k0 =      23.9788;
    double k1 =      -0.12875;
    double k2 =     662.0648;
    double k3 = -111298.07;
    double ct =       0.00836712;

    double v1 = -0.594e-6;
    double v2 =  0.0;
    double v3 = 75.8e-06;
    double v4 =  0.0;

    double h0 = 7788.0;
    double s0 =   35.545;
    double v0 =    0.709258*0.964596;

    *cps = k0 + k1/sqrt(t) + ct*t + k2/SQUARE(t) + k3/CUBE(t);
    *dcpsdt = - 0.5*k1/pow(t, (double) 1.5) + ct - 2.0*k2/CUBE(t)
    - 3.0*k3/QUARTIC(t);

    *hs = h0 + k0*(t-tr) + 0.5*ct*(t*t-tr*tr) + 2.0*k1*(sqrt(t)-sqrt(tr))
    - k2*(1.0/t-1.0/tr) - 0.5*k3*(1.0/(t*t)-1.0/(tr*tr));
    *ss = s0 + k0*log(t/tr)
    - 2.0*k1*(1.0/sqrt(t)-1.0/sqrt(tr)) + ct*(t-tr)
    - 0.5*k2*(1.0/(t*t)-1.0/(tr*tr))
    - (1.0/3.0)*k3*(1.0/(t*t*t)-1.0/(tr*tr*tr));

    *gs       = (*hs) - t*(*ss)
    + v0*((v1/2.0-v2)*(p*p-pr*pr)+v2*(p*p*p-pr*pr*pr)/3.0
          + (1.0-v1+v2+v3*(t-tr)+v4*SQUARE(t-tr))*(p-pr));
    *hs      += v0*((v1/2.0-v2)*(p*p-pr*pr)+v2*(p*p*p-pr*pr*pr)/3.0
                    + (1.0-v1+v2+v3*(t-tr)+v4*SQUARE(t-tr))*(p-pr))
    - t*v0*(v3 + 2.0*(t-tr)*v4)*(p-pr);
    *ss      += -v0*(v3 + 2.0*(t-tr)*v4)*(p-pr);
    *cps     += -t*v0*2.0*v4*(p-pr);
    *dcpsdt  += -v0*2.0*v4*(p-pr);
    *dvsdt    = v0*(v3 + 2.0*(t-tr)*v4);
    *dvsdp    = v0*(v1 + 2.0*(p-pr)*v2);
    *d2vsdt2  = v0*2.0*v4;
    *d2vsdp2  = v0*2.0*v2;
    *d2vsdtdp = 0.0;
    *vs       = v0*(1.0 + (p-pr)*v1 + (p-pr)*(p-pr)*v2 + (t-tr)*v3
                    + (t-tr)*(t-tr)*v4);
}


static void fe_liquid(double t, double p, double *gl, double *hl,
                      double *sl, double *cpl, double *dcpldt, double *vl, double *dvldt,
                      double *dvldp, double *d2vldt2, double *d2vldp2, double *d2vldtdp)
{
    double tr  =  298.15;
    double tr1 = 1200.0;
    double pr  =    1.0;

    double h0 = 12395.0+32028.0;
    double s0 = 80.932;

    double v0 =  0.652697;
    double v1 = -8.87646e-7;
    double v2 =  1.306567e-12;
    double v3 =  1.31136e-4;
    double v4 =  7.2557e-9;

    double cp = 46.024;

    *cpl    = cp;
    *dcpldt = 0.0;

    *hl = h0 + cp*(t-tr1);
    *sl = s0 + cp*log(t/tr1);

    *gl       = (*hl) - t*(*sl)
    + v0*((v1/2.0-v2)*(p*p-pr*pr)+v2*(p*p*p-pr*pr*pr)/3.0
          + (1.0-v1+v2+v3*(t-tr)+v4*SQUARE(t-tr))*(p-pr));
    *hl      += v0*((v1/2.0-v2)*(p*p-pr*pr)+v2*(p*p*p-pr*pr*pr)/3.0
                    + (1.0-v1+v2+v3*(t-tr)+v4*SQUARE(t-tr))*(p-pr))
    - t*v0*(v3 + 2.0*(t-tr)*v4)*(p-pr);
    *sl      += -v0*(v3 + 2.0*(t-tr)*v4)*(p-pr);
    *cpl     += -t*v0*2.0*v4*(p-pr);
    *dcpldt  += -v0*2.0*v4*(p-pr);
    *dvldt    = v0*(v3 + 2.0*(t-tr)*v4);
    *dvldp    = v0*(v1 + 2.0*(p-pr)*v2);
    *d2vldt2  = v0*2.0*v4;
    *d2vldp2  = v0*2.0*v2;
    *d2vldtdp = 0.0;
    *vl       = v0*(1.0 + (p-pr)*v1 + (p-pr)*(p-pr)*v2 + (t-tr)*v3
                    + (t-tr)*(t-tr)*v4);
}

static void ni_metal(double t, double p, double *gs, double *hs,
                     double *ss, double *cps, double *dcpsdt, double *vs, double *dvsdt,
                     double *dvsdp, double *d2vsdt2, double *d2vsdp2, double *d2vsdtdp)
{
    double tr1   = 298.15;
    double tr    = 750.0;
    double pr    =   1.0;

    double h0    = 3.32124*4184.0;
    double s0    = 0.013795*4184.0;

    double v0    =  0.659;
    double vni1  = -0.536e-6;
    double vni3a = 47.900e-6;
    double vni3b = 53.900e-6;

    double cpa   =     5.071e-3*4184.0;
    double cpb   = 2.0*1.157e-6*4184.0;
    double cpc   =    -3.136e-2*4184.0;

    *cps      = cpa + cpb*t - cpc/(t*t);
    *dcpsdt   = cpb + 2.0*cpc/(t*t*t);

    *hs       = h0 + cpa*(t-tr) + 0.5*cpb*(t*t-tr*tr) + cpc*(1.0/t-1.0/tr);
    *ss       = s0 + cpa*(log(t)-log(tr)) + cpb*(t-tr) + 0.5*cpc*(1/t/t-1/tr/tr);

    *gs       = (*hs) - t*(*ss)
    + v0*(1.0+vni3a*(631.0-tr1))*((vni1/2.0)*(p*p-pr*pr) + (1.0-vni1+vni3b*(t-631.0))*(p-pr));

    *hs      += v0*(1.0+vni3a*(631.0-tr1))*(0.5*vni1*(p*p-pr*pr) + (1.0-vni1+vni3b*(t-631.0))*(p-pr))
    - t*v0*(1.0+vni3a*(631.0-tr1))*vni3b*(p-pr);
    *ss      += -v0*(1.0+vni3a*(631.0-tr1))*vni3b*(p-pr);
    *cps     += 0.0;
    *dcpsdt  += 0.0;
    *dvsdt    = v0*(1.0+vni3a*(631.0-tr1))*vni3b;
    *dvsdp    = v0*(1.0+vni3a*(631.0-tr1))*vni1;
    *d2vsdt2  = 0.0;
    *d2vsdp2  = 0.0;
    *d2vsdtdp = 0.0;
    *vs       = v0*(1.0+vni3a*(631.0-tr1))*(vni1*p + (1.0-vni1+vni3b*(t-631.0)));
}


static void ni_liquid(double t, double p, double *gl, double *hl,
                      double *sl, double *cpl, double *dcpldt, double *vl, double *dvldt,
                      double *dvldp, double *d2vldt2, double *d2vldp2, double *d2vldtdp)
{
    double tr  =  298.15;
    double tr1 = 1300.0;
    double pr  =    1.0;

    double h0 = 17479.0+30382.0;
    double s0 = 84.680;

    double v0 =  0.6286951;
    double v1 = -8.87646e-7;
    double v2 =  1.306567e-12;
    double v3 =  1.2486e-4;
    double v4 =  3.2325e-9;

    double cp = 38.911;

    *cpl    = cp;
    *dcpldt = 0.0;

    *hl = h0 + cp*(t-tr1);
    *sl = s0 + cp*log(t/tr1);

    *gl       = (*hl) - t*(*sl)
    + v0*((v1/2.0-v2)*(p*p-pr*pr)+v2*(p*p*p-pr*pr*pr)/3.0
          + (1.0-v1+v2+v3*(t-tr)+v4*SQUARE(t-tr))*(p-pr));
    *hl      += v0*((v1/2.0-v2)*(p*p-pr*pr)+v2*(p*p*p-pr*pr*pr)/3.0
                    + (1.0-v1+v2+v3*(t-tr)+v4*SQUARE(t-tr))*(p-pr))
    - t*v0*(v3 + 2.0*(t-tr)*v4)*(p-pr);
    *sl      += -v0*(v3 + 2.0*(t-tr)*v4)*(p-pr);
    *cpl     += -t*v0*2.0*v4*(p-pr);
    *dcpldt  += -v0*2.0*v4*(p-pr);
    *dvldt    = v0*(v3 + 2.0*(t-tr)*v4);
    *dvldp    = v0*(v1 + 2.0*(p-pr)*v2);
    *d2vldt2  = v0*2.0*v4;
    *d2vldp2  = v0*2.0*v2;
    *d2vldtdp = 0.0;
    *vl       = v0*(1.0 + (p-pr)*v1 + (p-pr)*(p-pr)*v2 + (t-tr)*v3
                    + (t-tr)*(t-tr)*v4);
}

/* END-OF-FILE: Gibbs.c */
