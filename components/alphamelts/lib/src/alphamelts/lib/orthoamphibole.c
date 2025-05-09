const char *orthoamphibole_ver(void) { return "$Id: orthoamphibole.c,v 1.3 2006/10/20 00:59:22 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: amphibole.c,v $
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
MELTS Source Code: RCS Revision 1.1.1.1  2006/08/15 16:57:35  ghiorso
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
 * Revision 1.4  1997/06/21  22:50:10  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 1.3  1997/05/03  20:23:45  ghiorso
 * *** empty log message ***
 *
 * Revision 1.2  1997/03/27  17:03:53  ghiorso
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
**      Routines to compute amphibole solution properties
**      (file: ORTHOAMPHIBOLE.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  June 1, 1996 Original Version
**              Canabalized  PYROXENE.C file
**--
*/

#ifdef ISCLINO
#undef ISCLINO
#endif

#ifdef DEBUG
#undef DEBUG
#endif

#define MAX_ITER 200  /* Maximum number of iterations allowed in order */

#include "melts_gsl.h"
#include "silmin.h"   /* Structure definitions for SILMIN package      */

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

/*
 *=============================================================================
 * Amphibole solution parameters:
 * Ghiorso, M.S. and B.W. Evans (1996)
 *   Work in progress.
 */

#define  R       8.3143
#define  NR      2       /* Two independent composition variables */
#define  NS      2       /* Two ordering parameters               */
#define  NA      3       /* Seven endmember compositions          */

#define  DH1      (-12073132.08) - (-12067920.38) /* Ortho-Clino: Mg7Si8O22(OH)2     (joules) */
#define  DS1      (535.2587) - (540.2587)         /*                               (joules/K) */
#define  DV1      (26.310) - (26.340)             /*                             (joules/bar) */
#define  DH2      (-9627014.85) - (-9623550.0)    /* Ortho-Clino: Fe7Si8O22(OH)2     (joules) */
#define  DS2      (720.0) - (725.0)               /*                               (joules/K) */
#define  DV2      (27.810) - (27.840)             /*                             (joules/bar) */
#define  DH3      (-12307863.0) - (-12307863.0)   /* Ortho-Clino: Ca2Mg5Si8O22(OH)2  (joules) */
#define  DS3      (544.500) - (549.500)           /*                               (joules/K) */
#define  DV3      (27.312) - (27.312)             /*                             (joules/bar) */

#define  GEX        -25.8556    * 1000.0 /* joules */
#define  GDEX         5.90595   * 1000.0 /* joules */

#define  cGX123_4    12.1027    * 1000.0 /* joules */
#define  cGX13_24     3.29327   * 1000.0 /* joules */
#define  cGX134_2    15.1135    * 1000.0 /* joules */
#define  cW123        5.11564   * 1000.0 /* joules */
#define  cW4          3.55954   * 1000.0 /* joules */
#define  cDZ         -0.781732  * 1000.0 /* joules */
#define  cGXCA_13    17.801     * 1000.0 /* joules */
#define  cGXCA_2    -11.395     * 1000.0 /* joules */
#define  cWCAMG      28.100     * 1000.0 /* joules */
#define  cWCAFE      19.500     * 1000.0 /* joules */
#define  cDWCAMG      0.0000    * 1000.0 /* joules */
#define  cDWCAFE      0.0000    * 1000.0 /* joules */

#define  oGX123_4    13.7354    * 1000.0 /* joules */
#define  oGX13_24     5.25246   * 1000.0 /* joules */
#define  oGX134_2    16.7461    * 1000.0 /* joules */
#define  oW123        5.52380   * 1000.0 /* joules */
#define  oW4          3.72280   * 1000.0 /* joules */
#define  oDZ         -0.210304  * 1000.0 /* joules */
#define  oGXCA_13    17.801     * 1000.0 /* joules */
#define  oGXCA_2    -11.395     * 1000.0 /* joules */
#define  oWCAMG      30.100     * 1000.0 /* joules */
#define  oWCAFE      19.500     * 1000.0 /* joules */
#define  oDWCAMG      0.0000    * 1000.0 /* joules */
#define  oDWCAFE      0.0000    * 1000.0 /* joules */

/*
 *=============================================================================
 * Dependent parameters (DO NOT change definitions below this line:
 */

#define  DG1     (DH1)-t*(DS1)+(p-1.0)*DV1
#define  DG2     (DH2)-t*(DS2)+(p-1.0)*DV2
#define  DG3     (DH3)-t*(DS3)+(p-1.0)*DV3

#define  GX123_4  (clino) ? (cGX123_4) : (oGX123_4)
#define  GX13_24  (clino) ? (cGX13_24) : (oGX13_24)
#define  GX134_2  (clino) ? (cGX134_2) : (oGX134_2)
#define  W123     (clino) ? (cW123)    : (oW123)
#define  W4       (clino) ? (cW4)      : (oW4)
#define  DZ       (clino) ? (cDZ)      : (oDZ)
#define  GXCA_13  (clino) ? (cGXCA_13) : (oGXCA_13)
#define  GXCA_2   (clino) ? (cGXCA_2)  : (oGXCA_2)
#define  WCAMG    (clino) ? (cWCAMG)   : (oWCAMG)
#define  WCAFE    (clino) ? (cWCAFE)   : (oWCAFE)
#define  DWCAMG   (clino) ? (cDWCAMG)  : (oDWCAMG)
#define  DWCAFE   (clino) ? (cDWCAFE)  : (oDWCAFE)

/*
 * Vertices of composition space
 */

         /* Mg7Si8O22(OH)2 */
#define  H1      (clino) ? 0.0 : (DH1)
#define  S1      (clino) ? 0.0 : (DS1)
#define  V1      (clino) ? 0.0 : (DV1)
#define  G1      (clino) ? 0.0 : (DH1)-t*(DS1)+(p-1.0)*(DV1)
         /* Fe7Si8O22(OH)2 */
#define  H2      (clino) ? 0.0 : (DH2)
#define  S2      (clino) ? 0.0 : (DS2)
#define  V2      (clino) ? 0.0 : (DV2)
#define  G2      (clino) ? 0.0 : (DH2)-t*(DS2)+(p-1.0)*(DV2)
         /* Ca2Mg5Si8O22(OH)2 */
#define  H3      (clino) ? 0.0 : (DH3)
#define  S3      (clino) ? 0.0 : (DS3)
#define  V3      (clino) ? 0.0 : (DV3)
#define  G3      (clino) ? 0.0 : (DH3)-t*(DS3)+(p-1.0)*(DV3)

/*
 * Definitions of Taylor expansion coefficients in terms of solution
 * parameters. Independent variables are r, p, s13, s2
 */

#define  H0         (H1)/2.0 + (H2)/2.0 + 5.0*(W123)/4.0 + (W4)/2.0 + (GX123_4)/4.0
#define  HR        -(H1)/2.0 + (H2)/2.0
#define  HS13       3.0*(GEX)/5.0 + (GDEX)
#define  HS2        2.0*(GEX)/5.0 - (GDEX)
#define  HRR       -5.0*(W123)/4.0 - (W4)/2.0 - (GX123_4)/4.0
#define  HRS13      5.0*(W123)/14.0 - 6.0*(W4)/7.0 - 5.0*(GX123_4)/28.0 + (DZ)/2.0
#define  HRS2       15.0*(W123)/14.0 - 4.0*(W4)/7.0 - (GX123_4)/28.0 - (DZ)/2.0
#define  HS13S13   -125.0*(W123)/98.0 - 18.0*(W4)/49.0 - 43.0*(GX123_4)/196.0 \
                                        + (GX13_24)/2.0 - (DZ)/14.0
#define  HS13S2     115.0*(W123)/49.0 - 24.0*(W4)/49.0 + 30.0*(GX123_4)/49.0 \
                                        - (GX13_24)/2.0 - (GX134_2)/2.0 - (DZ)/7.0
#define  HS2S2     -145.0*(W123)/98.0 - 8.0*(W4)/49.0 - 37.0*(GX123_4)/196.0 \
                                        + (GX134_2)/2.0 + 3.0*(DZ)/14.0

#define  HP        -(H1) + (H3) - (W4) + (WCAFE) + (WCAMG) + (GXCA_13)/2.0 \
                                        + (GXCA_2) - (DWCAFE)/2.0 - (DWCAMG)/2.0
#define  HRP       -(W4) + (WCAFE) - (WCAMG) + (GXCA_13)/2.0 + (GXCA_2) \
                                        - (DWCAFE) + (DWCAMG)
#define  HPS13     -6.0*(W4)/7.0 + 6.0*(WCAFE)/7.0 - 6.0*(WCAMG)/7.0 \
                                        - 4.0*(GXCA_13)/7.0 - (GXCA_2)/7.0 - 6.0*(DWCAFE)/7.0 + 6.0*(DWCAMG)/7.0
#define  HPS2      -4.0*(W4)/7.0 + 4.0*(WCAFE)/7.0 - 4.0*(WCAMG)/7.0 \
                                        + 2.0*(GXCA_13)/7.0 - 3.0*(GXCA_2)/7.0 - 4.0*(DWCAFE)/7.0 + 4.0*(DWCAMG)/7.0
#define  HPP       -2.0*(WCAMG) + (DWCAFE) + 3.0*(DWCAMG)

#define  HPS13S13  -18.0*(DWCAFE)/49.0 - 18.0*(DWCAMG)/49.0
#define  HPPS13     6.0*(DWCAFE)/7.0 - 18.0*(DWCAMG)/7.0
#define  HPPS2      4.0*(DWCAFE)/7.0 - 12.0*(DWCAMG)/7.0
#define  HRPS13    -6.0*(DWCAFE)/7.0 - 6.0*(DWCAMG)/7.0
#define  HRPS2     -4.0*(DWCAFE)/7.0 - 4.0*(DWCAMG)/7.0
#define  HPS2S2    -8.0*(DWCAFE)/49.0 - 8.0*(DWCAMG)/49.0
#define  HPS13S2   -24.0*(DWCAFE)/49.0 - 24.0*(DWCAMG)/49.0
#define  HRRP      -(DWCAFE)/2.0 - (DWCAMG)/2.0
#define  HRPP       (DWCAFE) - 3.0*(DWCAMG)
#define  HPPP      -4.0*(DWCAMG)

#define  S0         (S1)/2.0 + (S2)/2.0
#define  SR        -(S1)/2.0 + (S2)/2.0
#define  SS13       0.0
#define  SS2        0.0
#define  SRR        0.0
#define  SRS13      0.0
#define  SRS2       0.0
#define  SS13S13    0.0
#define  SS13S2     0.0
#define  SS2S2      0.0

#define  SP        -(S1) + (S3)
#define  SRP        0.0
#define  SPS13      0.0
#define  SPS2       0.0
#define  SPP        0.0

#define  SPS13S13   0.0
#define  SPPS13     0.0
#define  SPPS2      0.0
#define  SRPS13     0.0
#define  SRPS2      0.0
#define  SPS2S2     0.0
#define  SPS13S2    0.0
#define  SRRP       0.0
#define  SRPP       0.0
#define  SPPP       0.0

#define  V0         (V1)/2.0 + (V2)/2.0
#define  VR        -(V1)/2.0 + (V2)/2.0
#define  VS13       0.0
#define  VS2        0.0
#define  VRR        0.0
#define  VRS13      0.0
#define  VRS2       0.0
#define  VS13S13    0.0
#define  VS13S2     0.0
#define  VS2S2      0.0

#define  VP        -(V1) + (V3)
#define  VRP        0.0
#define  VPS13      0.0
#define  VPS2       0.0
#define  VPP        0.0

#define  VPS13S13   0.0
#define  VPPS13     0.0
#define  VPPS2      0.0
#define  VRPS13     0.0
#define  VRPS2      0.0
#define  VPS2S2     0.0
#define  VPS13S2    0.0
#define  VRRP       0.0
#define  VRPP       0.0
#define  VPPP       0.0

#define  G0         (G1)/2.0 + (G2)/2.0 + 5.0*(W123)/4.0 + (W4)/2.0 + (GX123_4)/4.0
#define  GR        -(G1)/2.0 + (G2)/2.0
#define  GS13       3.0*(GEX)/5.0 + (GDEX)
#define  GS2        2.0*(GEX)/5.0 - (GDEX)
#define  GRR       -5.0*(W123)/4.0 - (W4)/2.0 - (GX123_4)/4.0
#define  GRS13      5.0*(W123)/14.0 - 6.0*(W4)/7.0 - 5.0*(GX123_4)/28.0 + (DZ)/2.0
#define  GRS2       15.0*(W123)/14.0 - 4.0*(W4)/7.0 - (GX123_4)/28.0 - (DZ)/2.0
#define  GS13S13   -125.0*(W123)/98.0 - 18.0*(W4)/49.0 - 43.0*(GX123_4)/196.0 \
                                        + (GX13_24)/2.0 - (DZ)/14.0
#define  GS13S2     115.0*(W123)/49.0 - 24.0*(W4)/49.0 + 30.0*(GX123_4)/49.0 \
                                        - (GX13_24)/2.0 - (GX134_2)/2.0 - (DZ)/7.0
#define  GS2S2     -145.0*(W123)/98.0 - 8.0*(W4)/49.0 - 37.0*(GX123_4)/196.0 \
                                        + (GX134_2)/2.0 + 3.0*(DZ)/14.0

#define  GP        -(G1) + (G3) - (W4) + (WCAFE) + (WCAMG) + (GXCA_13)/2.0 \
                                        + (GXCA_2) - (DWCAFE)/2.0 - (DWCAMG)/2.0
#define  GRP       -(W4) + (WCAFE) - (WCAMG) + (GXCA_13)/2.0 + (GXCA_2) \
                                        - (DWCAFE) + (DWCAMG)
#define  GPS13     -6.0*(W4)/7.0 + 6.0*(WCAFE)/7.0 - 6.0*(WCAMG)/7.0 \
                                        - 4.0*(GXCA_13)/7.0 - (GXCA_2)/7.0 - 6.0*(DWCAFE)/7.0 + 6.0*(DWCAMG)/7.0
#define  GPS2      -4.0*(W4)/7.0 + 4.0*(WCAFE)/7.0 - 4.0*(WCAMG)/7.0 \
                                        + 2.0*(GXCA_13)/7.0 - 3.0*(GXCA_2)/7.0 - 4.0*(DWCAFE)/7.0 + 4.0*(DWCAMG)/7.0
#define  GPP       -2.0*(WCAMG) + (DWCAFE) + 3.0*(DWCAMG)

#define  GPS13S13  -18.0*(DWCAFE)/49.0 - 18.0*(DWCAMG)/49.0
#define  GPPS13     6.0*(DWCAFE)/7.0 - 18.0*(DWCAMG)/7.0
#define  GPPS2      4.0*(DWCAFE)/7.0 - 12.0*(DWCAMG)/7.0
#define  GRPS13    -6.0*(DWCAFE)/7.0 - 6.0*(DWCAMG)/7.0
#define  GRPS2     -4.0*(DWCAFE)/7.0 - 4.0*(DWCAMG)/7.0
#define  GPS2S2    -8.0*(DWCAFE)/49.0 - 8.0*(DWCAMG)/49.0
#define  GPS13S2   -24.0*(DWCAFE)/49.0 - 24.0*(DWCAMG)/49.0
#define  GRRP      -(DWCAFE)/2.0 - (DWCAMG)/2.0
#define  GRPP       (DWCAFE) - 3.0*(DWCAMG)
#define  GPPP      -4.0*(DWCAMG)

/*
 * Global (to this file): variables
 */

/*******************************/
/* Static for Structural State */
/*******************************/

static MTHREAD_ONCE_T initThreadCBlock = MTHREAD_ONCE_INIT;

static MTHREAD_KEY_T clinoKey;

static void threadCInit(void) {
    MTHREAD_KEY_CREATE(&clinoKey, free);
}

static int getClino() {
    int *clinoPt;
    MTHREAD_ONCE(&initThreadCBlock, threadCInit);

    clinoPt = (int *) MTHREAD_GETSPECIFIC(clinoKey);
    if (clinoPt == NULL) {
        clinoPt  = (int *) malloc(sizeof(int));
        *clinoPt = FALSE;
        MTHREAD_SETSPECIFIC(clinoKey, (void *) clinoPt);
    }
    return *clinoPt;
}

static void setClino(int clino) {
    int *clinoPt;
    MTHREAD_ONCE(&initThreadCBlock, threadCInit);

    clinoPt = (int *) MTHREAD_GETSPECIFIC(clinoKey);
    if (clinoPt == NULL) {
        clinoPt  = (int *) malloc(sizeof(int));
        MTHREAD_SETSPECIFIC(clinoKey, (void *) clinoPt);
    }
    *clinoPt = clino;
}

/*************************************/
/* Statics for Ordering Calculations */
/*************************************/

static MTHREAD_ONCE_T initThreadOBlock = MTHREAD_ONCE_INIT;

static MTHREAD_KEY_T tOldKey;
static MTHREAD_KEY_T pOldKey;
static MTHREAD_KEY_T rOldKey;
static MTHREAD_KEY_T structOldKey;
static MTHREAD_KEY_T sOldKey;
static MTHREAD_KEY_T ptToD2gds2Key;
static MTHREAD_KEY_T d2gds2Key;
static MTHREAD_KEY_T indexD2gds2Key;

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
    MTHREAD_KEY_CREATE(&structOldKey,  free);
    MTHREAD_KEY_CREATE(&sOldKey,       freeNSarray);
    MTHREAD_KEY_CREATE(&ptToD2gds2Key, freePtToD2gds2);
    MTHREAD_KEY_CREATE(&d2gds2Key,     free);
    MTHREAD_KEY_CREATE(&indexD2gds2Key, freeIndexD2gds2);
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

static int getStructOld() {
    int *structOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    structOldPt = (int *) MTHREAD_GETSPECIFIC(structOldKey);
    if (structOldPt == NULL) {
        structOldPt  = (int *) malloc(sizeof(int));
        *structOldPt = FALSE;
        MTHREAD_SETSPECIFIC(structOldKey, (void *) structOldPt);
    }
    return *structOldPt;
}

static void setStructOld(int structOld) {
    int *structOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    structOldPt = (int *) MTHREAD_GETSPECIFIC(structOldKey);
    if (structOldPt == NULL) {
        structOldPt  = (int *) malloc(sizeof(int));
        *structOldPt = FALSE;
        MTHREAD_SETSPECIFIC(structOldKey, (void *) structOldPt);
    }
    *structOldPt = structOld;
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

/***********************************/
/* Statics for Site Mole Fractions */
/***********************************/

#define DECLARE_SITE_FRACTIONS \
    double xfe2m13, xmg2m13, xfe2m2, xmg2m2, xfe2m4, xmg2m4, xca2m4;

#define XFE2M13 0
#define XMG2M13 1
#define XFE2M2  2
#define XMG2M2  3
#define XFE2M4  4
#define XMG2M4  5
#define XCA2M4  6

#define NX      7

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
    xfe2m13 = getX(XFE2M13); \
    xmg2m13 = getX(XMG2M13); \
    xfe2m2  = getX(XFE2M2); \
    xmg2m2  = getX(XMG2M2); \
    xfe2m4  = getX(XFE2M4); \
    xmg2m4  = getX(XMG2M4); \
    xca2m4  = getX(XCA2M4);

#define SET_SITE_FRACTIONS \
    setX(XFE2M13, xfe2m13); \
    setX(XMG2M13, xmg2m13); \
    setX(XFE2M2,  xfe2m2); \
    setX(XMG2M2,  xmg2m2); \
    setX(XFE2M4,  xfe2m4); \
    setX(XMG2M4,  xmg2m4); \
    setX(XCA2M4,  xca2m4);

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conOph defines the conversion from m[i], to r[j]
 */
                   /* Order: r, p */
#define FR(i)     (i == 0 || i == 2) ? -(1.0+r[0]) : (1.0-r[0])
#define FP(i)     (i == 2)           ?   1.0-r[1]  :     -r[1]

                   /* Order: S13, S2 */
#define FS13(i)    - s[0]
#define FS2(i)     - s[1]

#define DFRDR(i)     - 1.0
#define DFPDP(i)     - 1.0

#define DFS13DS13(i) - 1.0
#define DFS2DS2(i)   - 1.0

/*
 * Global (to this file): derivative definitions
 */

#define SIC -R*(2.0*xfe2m4*log(xfe2m4) + 2.0*xmg2m4*log(xmg2m4) \
                                + 2.0*xca2m4*log(xca2m4) + 3.0*xfe2m13*log(xfe2m13) \
                                + 3.0*xmg2m13*log(xmg2m13) + 2.0*xfe2m2*log(xfe2m2) \
                                + 2.0*xmg2m2*log(xmg2m2) )
#define S   (SIC) + \
                        (S0) + (SR)*r[0] + (SS13)*s[0] + (SS2)*s[1] + (SRR)*r[0]*r[0] \
                    + (SRS13)*r[0]*s[0] + (SRS2)*r[0]*s[1] + (SS13S13)*s[0]*s[0] \
                    + (SS13S2)*s[0]*s[1] + (SS2S2)*s[1]*s[1] + (SP)*r[1] + (SRP)*r[0]*r[1] \
                    + (SPS13)*r[1]*s[0] + (SPS2)*r[1]*s[1] + (SPP)*r[1]*r[1] \
                    + (SPS13S13)*r[1]*s[0]*s[0] + (SPPS13)*r[1]*r[1]*s[0] \
                    + (SPPS2)*r[1]*r[1]*s[1] + (SRPS13)*r[0]*r[1]*s[0] + (SRPS2)*r[0]*r[1]*s[1] \
                    + (SPS2S2)*r[1]*s[1]*s[1] + (SPS13S2)*r[1]*s[0]*s[1] + (SRRP)*r[0]*r[0]*r[1] \
                    + (SRPP)*r[0]*r[1]*r[1] + (SPPP)*r[1]*r[1]*r[1]
#define H   (H0) + (HR)*r[0] + (HS13)*s[0] + (HS2)*s[1] + (HRR)*r[0]*r[0] \
                    + (HRS13)*r[0]*s[0] + (HRS2)*r[0]*s[1] + (HS13S13)*s[0]*s[0] \
                    + (HS13S2)*s[0]*s[1] + (HS2S2)*s[1]*s[1] + (HP)*r[1] + (HRP)*r[0]*r[1] \
                    + (HPS13)*r[1]*s[0] + (HPS2)*r[1]*s[1] + (HPP)*r[1]*r[1] \
                    + (HPS13S13)*r[1]*s[0]*s[0] + (HPPS13)*r[1]*r[1]*s[0] \
                    + (HPPS2)*r[1]*r[1]*s[1] + (HRPS13)*r[0]*r[1]*s[0] + (HRPS2)*r[0]*r[1]*s[1] \
                    + (HPS2S2)*r[1]*s[1]*s[1] + (HPS13S2)*r[1]*s[0]*s[1] + (HRRP)*r[0]*r[0]*r[1] \
                    + (HRPP)*r[0]*r[1]*r[1] + (HPPP)*r[1]*r[1]*r[1]
#define V   (V0) + (VR)*r[0] + (VS13)*s[0] + (VS2)*s[1] + (VRR)*r[0]*r[0] \
                    + (VRS13)*r[0]*s[0] + (VRS2)*r[0]*s[1] + (VS13S13)*s[0]*s[0] \
                    + (VS13S2)*s[0]*s[1] + (VS2S2)*s[1]*s[1] + (VP)*r[1] + (VRP)*r[0]*r[1] \
                    + (VPS13)*r[1]*s[0] + (VPS2)*r[1]*s[1] + (VPP)*r[1]*r[1] \
                    + (VPS13S13)*r[1]*s[0]*s[0] + (VPPS13)*r[1]*r[1]*s[0] \
                    + (VPPS2)*r[1]*r[1]*s[1] + (VRPS13)*r[0]*r[1]*s[0] + (VRPS2)*r[0]*r[1]*s[1] \
                    + (VPS2S2)*r[1]*s[1]*s[1] + (VPS13S2)*r[1]*s[0]*s[1] + (VRRP)*r[0]*r[0]*r[1] \
                    + (VRPP)*r[0]*r[1]*r[1] + (VPPP)*r[1]*r[1]*r[1]
#define G   -t*(SIC) + \
                        (G0) + (GR)*r[0] + (GS13)*s[0] + (GS2)*s[1] + (GRR)*r[0]*r[0] \
                    + (GRS13)*r[0]*s[0] + (GRS2)*r[0]*s[1] + (GS13S13)*s[0]*s[0] \
                    + (GS13S2)*s[0]*s[1] + (GS2S2)*s[1]*s[1] + (GP)*r[1] + (GRP)*r[0]*r[1] \
                    + (GPS13)*r[1]*s[0] + (GPS2)*r[1]*s[1] + (GPP)*r[1]*r[1] \
                    + (GPS13S13)*r[1]*s[0]*s[0] + (GPPS13)*r[1]*r[1]*s[0] \
                    + (GPPS2)*r[1]*r[1]*s[1] + (GRPS13)*r[0]*r[1]*s[0] + (GRPS2)*r[0]*r[1]*s[1] \
                    + (GPS2S2)*r[1]*s[1]*s[1] + (GPS13S2)*r[1]*s[0]*s[1] + (GRRP)*r[0]*r[0]*r[1] \
                    + (GRPP)*r[0]*r[1]*r[1] + (GPPP)*r[1]*r[1]*r[1]

/*----------------------------------------------------------------------------*/

#define DGDR0 R*t*(log(xfe2m4) - log(xmg2m4) + 3.0*log(xfe2m13)/2.0 \
                 - 3.0*log(xmg2m13)/2.0 + log(xfe2m2) - log(xmg2m2) ) \
                        + (GR) + 2.0*(GRR)*r[0] + (GRS13)*s[0] + (GRS2)*s[1] \
                        + (GRP)*r[1] + (GRPS13)*r[1]*s[0] + (GRPS2)*r[1]*s[1] \
                        + 2.0*(GRRP)*r[0]*r[1] + (GRPP)*r[1]*r[1]
#define DGDR1 R*t*(-2.0*log(xmg2m4) + 2.0*log(xca2m4) ) \
                        + (GP) + (GRP)*r[0] + (GPS13)*s[0] + (GPS2)*s[1] + 2.0*(GPP)*r[1] \
                        + (GPS13S13)*s[0]*s[0] + 2.0*(GPPS13)*r[1]*s[0] \
                        + 2.0*(GPPS2)*r[1]*s[1] + (GRPS13)*r[0]*s[0] + (GRPS2)*r[0]*s[1] \
                        + (GPS2S2)*s[1]*s[1] + (GPS13S2)*s[0]*s[1] + (GRRP)*r[0]*r[0] \
                        + 2.0*(GRPP)*r[0]*r[1] + 3.0*(GPPP)*r[1]*r[1]
#define DGDS0 R*t*(6.0*log(xfe2m4)/7.0 - 6.0*log(xmg2m4)/7.0 \
                 - 12.0*log(xfe2m13)/7.0 + 12.0*log(xmg2m13)/7.0 \
                 + 6.0*log(xfe2m2)/7.0 - 6.0*log(xmg2m2)/7.0 ) \
                        + (GS13) + (GRS13)*r[0] + 2.0*(GS13S13)*s[0] + (GS13S2)*s[1] \
                        + (GPS13)*r[1] + 2.0*(GPS13S13)*r[1]*s[0] + (GPPS13)*r[1]*r[1] \
                        + (GRPS13)*r[0]*r[1] + (GPS13S2)*r[1]*s[1]
#define DGDS1 R*t*(4.0*log(xfe2m4)/7.0 - 4.0*log(xmg2m4)/7.0 \
                 + 6.0*log(xfe2m13)/7.0 - 6.0*log(xmg2m13)/7.0 \
                 - 10.0*log(xfe2m2)/7.0 + 10.0*log(xmg2m2)/7.0 ) \
                        + (GS2) + (GRS2)*r[0] + (GS13S2)*s[0] + 2.0*(GS2S2)*s[1] \
                        + (GPS2)*r[1] + (GPPS2)*r[1]*r[1] + (GRPS2)*r[0]*r[1] \
                        + 2.0*(GPS2S2)*r[1]*s[1] + (GPS13S2)*r[1]*s[0]
#define DGDT  -(S)
#define DGDP   (V)

/*----------------------------------------------------------------------------*/

#define D2GDR0R0 R*t*(1.0/(2.0*xfe2m4) + 1.0/(2.0*xmg2m4) + 3.0/(4.0*xfe2m13) \
                                        + 3.0/(4.0*xmg2m13) + 1.0/(2.0*xfe2m2) + 1.0/(2.0*xmg2m2) ) \
               + 2.0*(GRR) + 2.0*(GRRP)*r[1]
#define D2GDR0R1 R*t*(1.0/xmg2m4) + (GRP) + (GRPS13)*s[0] + (GRPS2)*s[1] \
               + 2.0*(GRRP)*r[0] + 2.0*(GRPP)*r[1]
#define D2GDR0S0 R*t*(3.0/(7.0*xfe2m4) + 3.0/(7.0*xmg2m4) - 6.0/(7.0*xfe2m13) \
                                        - 6.0/(7.0*xmg2m13) + 3.0/(7.0*xfe2m2) + 3.0/(7.0*xmg2m2) ) \
               + (GRS13) + (GRPS13)*r[1]
#define D2GDR0S1 R*t*(2.0/(7.0*xfe2m4) + 2.0/(7.0*xmg2m4) + 3.0/(7.0*xfe2m13) \
                                        + 3.0/(7.0*xmg2m13) - 5.0/(7.0*xfe2m2) - 5.0/(7.0*xmg2m2) ) \
               + (GRS2) + (GRPS2)*r[1]
#define D2GDR0DT R*(log(xfe2m4) - log(xmg2m4) + 3.0*log(xfe2m13)/2.0 \
                                    - 3.0*log(xmg2m13)/2.0 + log(xfe2m2) - log(xmg2m2) ) \
               - (SR) - 2.0*(SRR)*r[0] - (SRS13)*s[0] - (SRS2)*s[1] \
               - (SRP)*r[1] - (SRPS13)*r[1]*s[0] - (SRPS2)*r[1]*s[1] \
               - 2.0*(SRRP)*r[0]*r[1] - (SRPP)*r[1]*r[1]
#define D2GDR0DP (VR) + 2.0*(VRR)*r[0] + (VRS13)*s[0] + (VRS2)*s[1] \
               + (VRP)*r[1] + (VRPS13)*r[1]*s[0] + (VRPS2)*r[1]*s[1] \
               + 2.0*(VRRP)*r[0]*r[1] + (VRPP)*r[1]*r[1]

#define D2GDR1R1 R*t*(2.0/(xmg2m4) + 2.0/(xca2m4) ) \
               + 2.0*(GPP) + 2.0*(GPPS13)*s[0] + 2.0*(GPPS2)*s[1] \
               + 2.0*(GRPP)*r[0] + 6.0*(GPPP)*r[1]
#define D2GDR1S0 R*t*(6.0/(7.0*xmg2m4) ) \
               + (GPS13) + 2.0*(GPS13S13)*s[0] + 2.0*(GPPS13)*r[1] \
               + (GRPS13)*r[0] + (GPS13S2)*s[1]
#define D2GDR1S1 R*t*(4.0/(7.0*xmg2m4) ) \
               + (GPS2) + 2.0*(GPPS2)*r[1] + (GRPS2)*r[0] \
               + 2.0*(GPS2S2)*s[1] + (GPS13S2)*s[0]
#define D2GDR1DT R*(-2.0*log(xmg2m4) + 2.0*log(xca2m4) ) \
               - (SP) - (SRP)*r[0] - (SPS13)*s[0] - (SPS2)*s[1] - 2.0*(SPP)*r[1] \
               - (SPS13S13)*s[0]*s[0] - 2.0*(SPPS13)*r[1]*s[0] \
               - 2.0*(SPPS2)*r[1]*s[1] - (SRPS13)*r[0]*s[0] - (SRPS2)*r[0]*s[1] \
               - (SPS2S2)*s[1]*s[1] - (SPS13S2)*s[0]*s[1] - (SRRP)*r[0]*r[0] \
               - 2.0*(SRPP)*r[0]*r[1] - 3.0*(SPPP)*r[1]*r[1]
#define D2GDR1DP (VP) + (VRP)*r[0] + (VPS13)*s[0] + (VPS2)*s[1] + 2.0*(VPP)*r[1] \
               + (VPS13S13)*s[0]*s[0] + 2.0*(VPPS13)*r[1]*s[0] \
               + 2.0*(VPPS2)*r[1]*s[1] + (VRPS13)*r[0]*s[0] + (VRPS2)*r[0]*s[1] \
               + (VPS2S2)*s[1]*s[1] + (VPS13S2)*s[0]*s[1] + (VRRP)*r[0]*r[0] \
               + 2.0*(VRPP)*r[0]*r[1] + 3.0*(VPPP)*r[1]*r[1]

#define D2GDS0S0 R*t*(6.0*3.0/(7.0*7.0*xfe2m4) + 6.0*3.0/(7.0*7.0*xmg2m4) \
                                        + 12.0*4.0/(7.0*7.0*xfe2m13) + 12.0*4.0/(7.0*7.0*xmg2m13) \
                                        + 6.0*3.0/(7.0*7.0*xfe2m2) + 6.0*3.0/(7.0*7.0*xmg2m2) ) \
               + 2.0*(GS13S13) + 2.0*(GPS13S13)*r[1]
#define D2GDS0S1 R*t*(6.0*2.0/(7.0*7.0*xfe2m4) + 6.0*2.0/(7.0*7.0*xmg2m4) \
                                        - 12.0*2.0/(7.0*7.0*xfe2m13) - 12.0*2.0/(7.0*7.0*xmg2m13) \
                                        - 6.0*5.0/(7.0*7.0*xfe2m2) - 6.0*5.0/(7.0*7.0*xmg2m2) ) \
               + (GS13S2) + (GPS13S2)*r[1]
#define D2GDS0DT R*(6.0*log(xfe2m4)/7.0 - 6.0*log(xmg2m4)/7.0 \
                                    - 12.0*log(xfe2m13)/7.0 + 12.0*log(xmg2m13)/7.0 \
                                    + 6.0*log(xfe2m2)/7.0 - 6.0*log(xmg2m2)/7.0 ) \
               - (SS13) - (SRS13)*r[0] - 2.0*(SS13S13)*s[0] - (SS13S2)*s[1] \
               - (SPS13)*r[1] - 2.0*(SPS13S13)*r[1]*s[0] - (SPPS13)*r[1]*r[1] \
               - (SRPS13)*r[0]*r[1] - (SPS13S2)*r[1]*s[1]
#define D2GDS0DP (VS13) + (VRS13)*r[0] + 2.0*(VS13S13)*s[0] + (VS13S2)*s[1] \
                                + (VPS13)*r[1] + 2.0*(VPS13S13)*r[1]*s[0] + (VPPS13)*r[1]*r[1] \
                                + (VRPS13)*r[0]*r[1] + (VPS13S2)*r[1]*s[1]

#define D2GDS1S1 R*t*(4.0*2.0/(7.0*7.0*xfe2m4) + 4.0*2.0/(7.0*7.0*xmg2m4) \
                                        + 6.0*2.0/(7.0*7.0*xfe2m13) + 6.0*2.0/(7.0*7.0*xmg2m13) \
                                        + 10.0*5.0/(7.0*7.0*xfe2m2) + 10.0*5.0/(7.0*7.0*xmg2m2) ) \
               + 2.0*(GS2S2) + 2.0*(GPS2S2)*r[1]
#define D2GDS1DT R*(4.0*log(xfe2m4)/7.0 - 4.0*log(xmg2m4)/7.0 \
                                    + 6.0*log(xfe2m13)/7.0 - 6.0*log(xmg2m13)/7.0 \
                                    - 10.0*log(xfe2m2)/7.0 + 10.0*log(xmg2m2)/7.0 ) \
               - (SS2) - (SRS2)*r[0] - (SS13S2)*s[0] - 2.0*(SS2S2)*s[1] \
               - (SPS2)*r[1] - (SPPS2)*r[1]*r[1] - (SRPS2)*r[0]*r[1] \
               - 2.0*(SPS2S2)*r[1]*s[1] - (SPS13S2)*r[1]*s[0]
#define D2GDS1DP (VS2) + (VRS2)*r[0] + (VS13S2)*s[0] + 2.0*(VS2S2)*s[1] \
               + (VPS2)*r[1] + (VPPS2)*r[1]*r[1] + (VRPS2)*r[0]*r[1] \
               + 2.0*(VPS2S2)*r[1]*s[1] + (VPS13S2)*r[1]*s[0]

#define D2GDT2   0.0
#define D2GDTDP  0.0
#define D2GDP2   0.0

/*----------------------------------------------------------------------------*/

#define D3GDR0R0R0 R*t*(-1.0/(4.0*xfe2m4*xfe2m4) + 1.0/(4.0*xmg2m4*xmg2m4) \
                                            - 3.0/(8.0*xfe2m13*xfe2m13) + 3.0/(8.0*xmg2m13*xmg2m13) \
                                            - 1.0/(4.0*xfe2m2*xfe2m2) + 1.0/(4.0*xmg2m2*xmg2m2) )
#define D3GDR0R0R1 R*t*(1.0/(2.0*xmg2m4*xmg2m4) ) + 2.0*(GRRP)
#define D3GDR0R0S0 R*t*(-3.0/(2.0*7.0*xfe2m4*xfe2m4) + 3.0/(2.0*7.0*xmg2m4*xmg2m4) \
                                            + 3.0*4.0/(4.0*7.0*xfe2m13*xfe2m13) - 3.0*4.0/(4.0*7.0*xmg2m13*xmg2m13) \
                                            - 3.0/(2.0*7.0*xfe2m2*xfe2m2) + 3.0/(2.0*7.0*xmg2m2*xmg2m2) )
#define D3GDR0R0S1 R*t*(-2.0/(2.0*7.0*xfe2m4*xfe2m4) + 2.0/(2.0*7.0*xmg2m4*xmg2m4) \
                                            - 3.0*2.0/(4.0*7.0*xfe2m13*xfe2m13) + 3.0*2.0/(4.0*7.0*xmg2m13*xmg2m13) \
                                            + 5.0/(2.0*7.0*xfe2m2*xfe2m2) - 5.0/(2.0*7.0*xmg2m2*xmg2m2) )
#define D3GDR0R0DT R*(1.0/(2.0*xfe2m4) + 1.0/(2.0*xmg2m4) + 3.0/(4.0*xfe2m13) \
                                            + 3.0/(4.0*xmg2m13) + 1.0/(2.0*xfe2m2) + 1.0/(2.0*xmg2m2) ) \
                                            - 2.0*(SRR) - 2.0*(SRRP)*r[1]
#define D3GDR0R0DP 2.0*(VRR) + 2.0*(VRRP)*r[1]

#define D3GDR0R1R1 R*t*(1.0/(xmg2m4*xmg2m4)) + 2.0*(GRPP)
#define D3GDR0R1S0 R*t*(3.0/(7.0*xmg2m4*xmg2m4)) + (GRPS13)
#define D3GDR0R1S1 R*t*(2.0/(7.0*xmg2m4*xmg2m4)) + (GRPS2)
#define D3GDR0R1DT R*(1.0/xmg2m4) - (SRP) - (SRPS13)*s[0] - (SRPS2)*s[1] \
                   - 2.0*(SRRP)*r[0] - 2.0*(SRPP)*r[1]
#define D3GDR0R1DP (VRP) + (VRPS13)*s[0] + (VRPS2)*s[1] + 2.0*(VRRP)*r[0] \
                   + 2.0*(VRPP)*r[1]

#define D3GDR0S0S0 R*t*(-3.0*3.0/(7.0*7.0*xfe2m4*xfe2m4) + 3.0*3.0/(7.0*7.0*xmg2m4*xmg2m4) \
                                            - 6.0*4.0/(7.0*7.0*xfe2m13*xfe2m13) + 6.0*4.0/(7.0*7.0*xmg2m13*xmg2m13) \
                                            - 3.0*3.0/(7.0*7.0*xfe2m2*xfe2m2) + 3.0*3.0/(7.0*7.0*xmg2m2*xmg2m2) )
#define D3GDR0S0S1 R*t*(-3.0*2.0/(7.0*7.0*xfe2m4*xfe2m4) + 3.0*2.0/(7.0*7.0*xmg2m4*xmg2m4) \
                                            + 6.0*2.0/(7.0*7.0*xfe2m13*xfe2m13) - 6.0*2.0/(7.0*7.0*xmg2m13*xmg2m13) \
                                            + 3.0*5.0/(7.0*7.0*xfe2m2*xfe2m2) - 3.0*5.0/(7.0*7.0*xmg2m2*xmg2m2) )
#define D3GDR0S0DT R*(3.0/(7.0*xfe2m4) + 3.0/(7.0*xmg2m4) - 6.0/(7.0*xfe2m13) \
                                        - 6.0/(7.0*xmg2m13) + 3.0/(7.0*xfe2m2) + 3.0/(7.0*xmg2m2) ) \
                   - (SRS13) - (SRPS13)*r[1]
#define D3GDR0S0DP (VRS13) + (VRPS13)*r[1]

#define D3GDR0S1S1 R*t*(-2.0*2.0/(7.0*7.0*xfe2m4*xfe2m4) + 2.0*2.0/(7.0*7.0*xmg2m4*xmg2m4) \
                                            - 3.0*2.0/(7.0*7.0*xfe2m13*xfe2m13) + 3.0*2.0/(7.0*7.0*xmg2m13*xmg2m13) \
                                            - 5.0*5.0/(7.0*7.0*xfe2m2*xfe2m2) + 5.0*5.0/(7.0*7.0*xmg2m2*xmg2m2) )
#define D3GDR0S1DT R*(2.0/(7.0*xfe2m4) + 2.0/(7.0*xmg2m4) + 3.0/(7.0*xfe2m13) \
                                        + 3.0/(7.0*xmg2m13) - 5.0/(7.0*xfe2m2) - 5.0/(7.0*xmg2m2) ) \
                   - (SRS2) - (SRPS2)*r[1]
#define D3GDR0S1DP (VRS2) + (VRPS2)*r[1]

#define D3GDR1R1R1 R*t*(2.0/(xmg2m4*xmg2m4) - 2.0/(xca2m4*xca2m4) ) + 6.0*(GPPP)
#define D3GDR1R1S0 R*t*(2.0*3.0/(7.0*xmg2m4*xmg2m4) ) + 2.0*(GPPS13)
#define D3GDR1R1S1 R*t*(2.0*2.0/(7.0*xmg2m4*xmg2m4) ) + 2.0*(GPPS2)
#define D3GDR1R1DT R*(2.0/(xmg2m4) + 2.0/(xca2m4) ) - 2.0*(SPP) \
                   - 2.0*(SPPS13)*s[0] - 2.0*(SPPS2)*s[1] - 2.0*(SRPP)*r[0] - 6.0*(SPPP)*r[1]
#define D3GDR1R1DP 2.0*(VPP) + 2.0*(VPPS13)*s[0] + 2.0*(VPPS2)*s[1] \
                   + 2.0*(VRPP)*r[0] + 6.0*(VPPP)*r[1]

#define D3GDR1S0S0 R*t*(6.0*3.0/(7.0*7.0*xmg2m4*xmg2m4) ) + 2.0*(GPS13S13)
#define D3GDR1S0S1 R*t*(6.0*2.0/(7.0*7.0*xmg2m4*xmg2m4) ) + (GPS13S2)
#define D3GDR1S0DT R*(6.0/(7.0*xmg2m4) ) \
                   - (SPS13) - 2.0*(SPS13S13)*s[0] - 2.0*(SPPS13)*r[1] \
                   - (SRPS13)*r[0] - (SPS13S2)*s[1]
#define D3GDR1S0DP (VPS13) + 2.0*(VPS13S13)*s[0] + 2.0*(VPPS13)*r[1] \
                   + (VRPS13)*r[0] + (VPS13S2)*s[1]

#define D3GDR1S1S1 R*t*(4.0*2.0/(7.0*7.0*xmg2m4*xmg2m4) ) + 2.0*(GPS2S2)
#define D3GDR1S1DT R*(4.0/(7.0*xmg2m4) ) \
                   - (SPS2) - 2.0*(SPPS2)*r[1] - (SRPS2)*r[0] \
                   - 2.0*(SPS2S2)*s[1] - (SPS13S2)*s[0]
#define D3GDR1S1DP (VPS2) + 2.0*(VPPS2)*r[1] + (VRPS2)*r[0] \
                   + 2.0*(VPS2S2)*s[1] + (VPS13S2)*s[0]

#define D3GDS0S0S0 R*t*(-6.0*3.0*3.0/(7.0*7.0*7.0*xfe2m4*xfe2m4) + 6.0*3.0*3.0/(7.0*7.0*7.0*xmg2m4*xmg2m4) \
                                            + 12.0*4.0*4.0/(7.0*7.0*7.0*xfe2m13*xfe2m13) - 12.0*4.0*4.0/(7.0*7.0*7.0*xmg2m13*xmg2m13) \
                                            - 6.0*3.0*3.0/(7.0*7.0*7.0*xfe2m2*xfe2m2) + 6.0*3.0*3.0/(7.0*7.0*7.0*xmg2m2*xmg2m2) )
#define D3GDS0S0S1 R*t*(-6.0*3.0*2.0/(7.0*7.0*7.0*xfe2m4*xfe2m4) + 6.0*3.0*2.0/(7.0*7.0*7.0*xmg2m4*xmg2m4) \
                                            - 12.0*4.0*2.0/(7.0*7.0*7.0*xfe2m13*xfe2m13) + 12.0*4.0*2.0/(7.0*7.0*7.0*xmg2m13*xmg2m13) \
                                            + 6.0*3.0*5.0/(7.0*7.0*7.0*xfe2m2*xfe2m2) - 6.0*3.0*5.0/(7.0*7.0*7.0*xmg2m2*xmg2m2) )
#define D3GDS0S0DT R*(6.0*3.0/(7.0*7.0*xfe2m4) + 6.0*3.0/(7.0*7.0*xmg2m4) \
                                        + 12.0*4.0/(7.0*7.0*xfe2m13) + 12.0*4.0/(7.0*7.0*xmg2m13) \
                                        + 6.0*3.0/(7.0*7.0*xfe2m2) + 6.0*3.0/(7.0*7.0*xmg2m2) ) \
                                        - 2.0*(SS13S13) - 2.0*(SPS13S13)*r[1]
#define D3GDS0S0DP 2.0*(VS13S13) + 2.0*(VPS13S13)*r[1]

#define D3GDS0S1S1 R*t*(-6.0*2.0*2.0/(7.0*7.0*7.0*xfe2m4*xfe2m4) + 6.0*2.0*2.0/(7.0*7.0*7.0*xmg2m4*xmg2m4) \
                                            + 12.0*2.0*2.0/(7.0*7.0*7.0*xfe2m13*xfe2m13) - 12.0*2.0*2.0/(7.0*7.0*7.0*xmg2m13*xmg2m13) \
                                            - 6.0*5.0*5.0/(7.0*7.0*7.0*xfe2m2*xfe2m2) + 6.0*5.0*5.0/(7.0*7.0*7.0*xmg2m2*xmg2m2) )
#define D3GDS0S1DT R*(6.0*2.0/(7.0*7.0*xfe2m4) + 6.0*2.0/(7.0*7.0*xmg2m4) \
                                        - 12.0*2.0/(7.0*7.0*xfe2m13) - 12.0*2.0/(7.0*7.0*xmg2m13) \
                                        - 6.0*5.0/(7.0*7.0*xfe2m2) - 6.0*5.0/(7.0*7.0*xmg2m2) ) \
                                        - (SS13S2) - (SPS13S2)*r[1]
#define D3GDS0S1DP (VS13S2) + (VPS13S2)*r[1]

#define D3GDS0DT2  0.0
#define D3GDS0DTDP 0.0
#define D3GDS0DP2  0.0

#define D3GDS1S1S1 R*t*(-4.0*2.0*2.0/(7.0*7.0*7.0*xfe2m4*xfe2m4) + 4.0*2.0*2.0/(7.0*7.0*7.0*xmg2m4*xmg2m4) \
                                            - 6.0*2.0*2.0/(7.0*7.0*7.0*xfe2m13*xfe2m13) + 6.0*2.0*2.0/(7.0*7.0*7.0*xmg2m13*xmg2m13) \
                                            + 10.0*5.0*5.0/(7.0*7.0*7.0*xfe2m2*xfe2m2) - 10.0*5.0*5.0/(7.0*7.0*7.0*xmg2m2*xmg2m2) )
#define D3GDS1S1DT R*(4.0*2.0/(7.0*7.0*xfe2m4) + 4.0*2.0/(7.0*7.0*xmg2m4) \
                                        + 6.0*2.0/(7.0*7.0*xfe2m13) + 6.0*2.0/(7.0*7.0*xmg2m13) \
                                        + 10.0*5.0/(7.0*7.0*xfe2m2) + 10.0*5.0/(7.0*7.0*xmg2m2) ) \
                                        - 2.0*(SS2S2) - 2.0*(SPS2S2)*r[1]
#define D3GDS1S1DP 2.0*(VS2S2) + 2.0*(VPS2S2)*r[1]

#define D3GDS1DT2  0.0
#define D3GDS1DTDP 0.0
#define D3GDS1DP2  0.0

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

/*
 *=============================================================================
 * Macros for automatic array initialization
 */
#define fillD2GDR2 \
 d2gdr2[0][0] = D2GDR0R0;     d2gdr2[0][1] = D2GDR0R1; \
 d2gdr2[1][0] = d2gdr2[0][1]; d2gdr2[1][1] = D2GDR1R1;

#define fillD2GDRDS \
 d2gdrds[0][0] = D2GDR0S0; d2gdrds[0][1] = D2GDR0S1; \
 d2gdrds[1][0] = D2GDR1S0; d2gdrds[1][1] = D2GDR1S1;

#define fillD2GDRDT \
 d2gdrdt[0] = D2GDR0DT; d2gdrdt[1] = D2GDR1DT;

#define fillD2GDRDP \
 d2gdrdp[0] = D2GDR0DP; d2gdrdp[1] = D2GDR1DP;

#define fillD2GDS2 \
 d2gds2[0][0] = D2GDS0S0;     d2gds2[0][1] = D2GDS0S1; \
 d2gds2[1][0] = d2gds2[0][1]; d2gds2[1][1] = D2GDS1S1;

#define fillD2GDSDT \
 d2gdsdt[0] = D2GDS0DT;  d2gdsdt[1] = D2GDS1DT;

#define fillD2GDSDP \
 d2gdsdp[0] = D2GDS0DP;  d2gdsdp[1] = D2GDS1DP;

#define fillD3GDR3 \
 d3gdr3[0][0][0] = D3GDR0R0R0;		d3gdr3[0][0][1] = D3GDR0R0R1; \
 d3gdr3[0][1][0] = d3gdr3[0][0][1];	d3gdr3[0][1][1] = D3GDR0R1R1; \
 d3gdr3[1][0][0] = d3gdr3[0][0][1];	d3gdr3[1][0][1] = d3gdr3[0][1][1]; \
 d3gdr3[1][1][0] = d3gdr3[0][1][1];	d3gdr3[1][1][1] = D3GDR1R1R1;

#define fillD3GDR2DS \
 d3gdr2ds[0][0][0] = D3GDR0R0S0;        d3gdr2ds[0][0][1] = D3GDR0R0S1; \
 d3gdr2ds[0][1][0] = D3GDR0R1S0;        d3gdr2ds[0][1][1] = D3GDR0R1S1; \
 d3gdr2ds[1][0][0] = d3gdr2ds[0][1][0]; d3gdr2ds[1][0][1] = d3gdr2ds[0][1][1]; \
 d3gdr2ds[1][1][0] = D3GDR1R1S0;        d3gdr2ds[1][1][1] = D3GDR1R1S1;

#define fillD3GDR2DT \
 d3gdr2dt[0][0] = D3GDR0R0DT;     d3gdr2dt[0][1] = D3GDR0R1DT; \
 d3gdr2dt[1][0] = d3gdr2dt[0][1]; d3gdr2dt[1][1] = D3GDR1R1DT;

#define fillD3GDR2DP \
 d3gdr2dp[0][0] = D3GDR0R0DP;     d3gdr2dp[0][1] = D3GDR0R1DP; \
 d3gdr2dp[1][0] = d3gdr2dp[0][1]; d3gdr2dp[1][1] = D3GDR1R1DP;

#define fillD3GDRDS2 \
 d3gdrds2[0][0][0] = D3GDR0S0S0;        d3gdrds2[0][0][1] = D3GDR0S0S1; \
 d3gdrds2[0][1][0] = d3gdrds2[0][0][1]; d3gdrds2[0][1][1] = D3GDR0S1S1; \
 d3gdrds2[1][0][0] = D3GDR1S0S0;        d3gdrds2[1][0][1] = D3GDR1S0S1; \
 d3gdrds2[1][1][0] = d3gdrds2[1][0][1]; d3gdrds2[1][1][1] = D3GDR1S1S1;

#define fillD3GDRDSDT \
 d3gdrdsdt[0][0] = D3GDR0S0DT;  d3gdrdsdt[0][1] = D3GDR0S1DT; \
 d3gdrdsdt[1][0] = D3GDR1S0DT;  d3gdrdsdt[1][1] = D3GDR1S1DT;

#define fillD3GDRDSDP \
 d3gdrdsdp[0][0] = D3GDR0S0DP; d3gdrdsdp[0][1] = D3GDR0S1DP; \
 d3gdrdsdp[1][0] = D3GDR1S0DP; d3gdrdsdp[1][1] = D3GDR1S1DP;

#define fillD3GDS3 \
 d3gds3[0][0][0] = D3GDS0S0S0;      d3gds3[0][0][1] = D3GDS0S0S1; \
 d3gds3[0][1][0] = d3gds3[0][0][1]; d3gds3[0][1][1] = D3GDS0S1S1; \
 d3gds3[1][0][0] = d3gds3[0][0][1]; d3gds3[1][0][1] = d3gds3[0][1][1]; \
 d3gds3[1][1][0] = d3gds3[0][1][1]; d3gds3[1][1][1] = D3GDS1S1S1;

#define fillD3GDS2DT \
 d3gds2dt[0][0] = D3GDS0S0DT;     d3gds2dt[0][1] = D3GDS0S1DT; \
 d3gds2dt[1][0] = d3gds2dt[0][1]; d3gds2dt[1][1] = D3GDS1S1DT;

#define fillD3GDS2DP \
 d3gds2dp[0][0] = D3GDS0S0DP;     d3gds2dp[0][1] = D3GDS0S1DP; \
 d3gds2dp[1][0] = d3gds2dp[0][1]; d3gds2dp[1][1] = D3GDS1S1DP;

#define fillD3GDSDT2 \
 d3gdsdt2[0] = D3GDS0DT2; d3gdsdt2[1] = D3GDS1DT2;

#define fillD3GDSDTDP \
 d3gdsdtdp[0] = D3GDS0DTDP; d3gdsdtdp[1] = D3GDS1DTDP;

#define fillD3GDSDP2 \
 d3gdsdp2[0] = D3GDS0DP2; d3gdsdp2[1] = D3GDS1DP2;

#define fillD3GDRDT2 \
 d3gdrdt2[0] = D3GDR0DT2; d3gdrdt2[1] = D3GDR1DT2;

#define fillD3GDRDTDP \
 d3gdrdtdp[0] = D3GDR0DTDP; d3gdrdtdp[1] = D3GDR1DTDP;

#define fillD3GDRDP2 \
 d3gdrdp2[0] = D3GDR0DP2; d3gdrdp2[1] = D3GDR1DP2;

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
    int clino           = getClino();
    double tOld         = getTOld();
    double pOld         = getPOld();
    double *rOld        = getROld();
    double *sOld        = getSOld();
    double **d2gds2     = getD2gds2();
    int    structOld    = getStructOld();
    gsl_matrix      *ptToD2gds2  = getPtToD2gds2();
    gsl_permutation *indexD2gds2 = getIndexD2gds2();

    int i, j, iter=0, signum;

    GET_SITE_FRACTIONS

    /* look-up or compute the current ordering state */
    if ( (t != tOld)       || (p != pOld)    || (structOld != clino) ||
       (r[0] != rOld[0]) || (r[1] != rOld[1]) ) {
        double dgds[NS], sNew[NS], xfe, xmg;
        gsl_vector_view vvToDgds = gsl_vector_view_array(dgds, (size_t) NS);
        structOld = clino;
        for (i=0; i<NS; i++) sOld[i] = 2.0;

        xca2m4  = r[1];
        xfe     = (1.0+r[0])/2.0;
        xmg     = (7.0-2.0*xca2m4-7.0*xfe)/7.0;
        sNew[0] = ((1.0-xca2m4)*xfe - xfe)/(xfe+xmg);
        sNew[1] = ((1.0-xca2m4)*xfe - xfe)/(xfe+xmg);

        if (xfe < 10.0*DBL_EPSILON) {
            xfe2m4  = DBL_EPSILON;
            xfe2m13 = DBL_EPSILON;
            xfe2m2  = DBL_EPSILON;
            xmg2m4  = 1.0 - xca2m4;
            xmg2m13 = 1.0;
            xmg2m2  = 1.0;
            sOld[0] = sNew[0];
            sOld[1] = sNew[1];
        }

        while ( ((ABS(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON) ||
             (ABS(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON)    ) &&
                        (iter < MAX_ITER)) {
            double s[NS], sCorr[NS], lambda = 1.0;
            gsl_vector_view vvToSCorr = gsl_vector_view_array(sCorr, (size_t) NS);

            for (i=0; i<NS; i++) s[i] = sNew[i];

            xfe2m4  = (1.0 + r[0] + 6.0*s[0]/7.0 +  4.0*s[1]/7.0)/2.0;
            xfe2m13 = (1.0 + r[0] - 8.0*s[0]/7.0 +  4.0*s[1]/7.0)/2.0;
            xfe2m2  = (1.0 + r[0] + 6.0*s[0]/7.0 - 10.0*s[1]/7.0)/2.0;
            xmg2m4  = 1.0 - xfe2m4  - xca2m4;
            xmg2m13 = 1.0 - xfe2m13;
            xmg2m2  = 1.0 - xfe2m2;

            if (xfe2m4  <= DBL_EPSILON) xfe2m4  = DBL_EPSILON;
            if (xfe2m13 <= DBL_EPSILON) xfe2m13 = DBL_EPSILON;
            if (xfe2m2  <= DBL_EPSILON) xfe2m2  = DBL_EPSILON;
            if (xmg2m4  <= DBL_EPSILON) xmg2m4  = DBL_EPSILON;
            if (xmg2m13 <= DBL_EPSILON) xmg2m13 = DBL_EPSILON;
            if (xmg2m2  <= DBL_EPSILON) xmg2m2  = DBL_EPSILON;

            dgds[0] = DGDS0;
            dgds[1] = DGDS1;

            fillD2GDS2

            for (i=0; i<NS; i++) sOld[i] = s[i];

            /* original: gaussj(ptToD2gds2, NS, (double **) NULL, 0);
            	for (i=0; i<NS; i++) {
                    for(j=0, sCorr[i]=0.0; j<NS; j++) sCorr[i] += - d2gds2[i][j]*dgds[j];
                    s[i] = sOld[i] + lambda*sCorr[i];
        	  }*/

            gsl_matrix_scale(ptToD2gds2, -1.0);
            melts_LU_decomp(ptToD2gds2, indexD2gds2, &signum);

            melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToDgds.vector, &vvToSCorr.vector);
            for (i=0; i<NS; i++) s[i] += lambda*sCorr[i];

            xfe2m4  = (1.0 + r[0] + 6.0*s[0]/7.0 +  4.0*s[1]/7.0)/2.0;
            xfe2m13 = (1.0 + r[0] - 8.0*s[0]/7.0 +  4.0*s[1]/7.0)/2.0;
            xfe2m2  = (1.0 + r[0] + 6.0*s[0]/7.0 - 10.0*s[1]/7.0)/2.0;
            xmg2m4  = 1.0 - xfe2m4  - xca2m4;
            xmg2m13 = 1.0 - xfe2m13;
            xmg2m2  = 1.0 - xfe2m2;

            while ((   xfe2m4  < 0.0 || xfe2m4  > 1.0 || xmg2m4  < 0.0 || xmg2m4  > 1.0
                            || xfe2m13 < 0.0 || xfe2m13 > 1.0 || xmg2m13 < 0.0 || xmg2m13 > 1.0
                            || xfe2m2  < 0.0 || xfe2m2  > 1.0 || xmg2m2  < 0.0 || xmg2m2  > 1.0)
             && lambda > DBL_EPSILON ) {
                lambda /= 2.0;
                for (i=0; i<NS; i++) s[i] = sOld[i] + lambda*sCorr[i];
                xfe2m4  = (1.0 + r[0] + 6.0*s[0]/7.0 +  4.0*s[1]/7.0)/2.0;
                xfe2m13 = (1.0 + r[0] - 8.0*s[0]/7.0 +  4.0*s[1]/7.0)/2.0;
                xfe2m2  = (1.0 + r[0] + 6.0*s[0]/7.0 - 10.0*s[1]/7.0)/2.0;
                xmg2m4  = 1.0 - xfe2m4  - xca2m4;
                xmg2m13 = 1.0 - xfe2m13;
                xmg2m2  = 1.0 - xfe2m2;
            }

            for (i=0; i<NS; i++) sNew[i] = s[i];
            iter++;
        }
        tOld = t;
        pOld = p;
        for (i=0; i<NR; i++) rOld[i] = r[i];

#ifdef DEBUG
        for (i=0; i<NS; i++) {
            if (dgds[i] > sqrt(DBL_EPSILON) && ABS(sOld[i]) > DBL_EPSILON) {
                printf("ERROR in ORTHOAMPHIBOLE.C (function ORDER). Failed to converge!\n");
                printf("  X2    = %13.6g, X3    = %13.6g\n", r[0], r[1]);
                printf("  s1    = %13.6g, s2    = %13.6g\n", sOld[0], sOld[1]);
                printf("  dgds1 = %13.6g, dgds2 = %13.6g\n", dgds[0], dgds[1]);
                printf("  X Ca2+ m4:  %13.6g\n", xca2m4);
                printf("  X Mg2+ m4:  %13.6g  X Mg2+ m13: %13.6g  X Mg2+ m2: %13.6g\n", xmg2m4, xmg2m13, xmg2m2);
                printf("  X Fe2+ m4:  %13.6g  X Fe2+ m13: %13.6g  X Fe2+ m2: %13.6g\n", xfe2m4, xfe2m13, xfe2m2);
                break;
            }
        }
#endif

        setTOld(tOld);
        setPOld(pOld);
        setStructOld(structOld);
        /* arrays (rOld, sOld, d2gds2, indexD2gds2) should be preserved automatically */

        SET_SITE_FRACTIONS
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
        double *s = sOld;
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
 * Private function to determine structural state of amphibole
 */

static int isClino(double t, double p, double r[NR])
{
#ifdef ISCLINO
    DECLARE_SITE_FRACTIONS
    double gOrtho, gClino, s[NS];
    int clino;

    clino = TRUE;
    setClino(clino);

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    GET_SITE_FRACTIONS
    gClino = G;

    clino = FALSE;
    setClino(clino);

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    GET_SITE_FRACTIONS
    gOrtho = G;

    if (gOrtho < gClino) return FALSE; else return TRUE;
#else
    return FALSE;
#endif
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
testOph(int mask, double t, double p,
    int na,          /* Expected number of endmember components                 */
    int nr,          /* Expected number of independent compositional variables  */
    char **names,    /* array of strings of names of endmember components       */
    char **formulas, /* array of strings of formulas of endmember components    */
    double *r,       /* array of indepependent compos variables, check bounds   */
    double *m)       /* array of moles of endmember components, check bounds    */
{
    const char *phase = "orthoamphibole.c";
    const char *NAMES[NA]    = { "cummingtonite", "grunerite", "tremolite" };
    const char *FORMULAS[NA] = { "Mg7Si8O22(OH)2", "Fe7Si8O22(OH)2", "Ca2Mg5Si8O22(OH)2" };
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
        result = result && ((r[0]+1.0)/2.0 >= 0.0)
                                        && ((r[0]+1.0)/2.0 <= 1.0-2.0*r[1]/7.0 );         /* Fe2+ */
        result = result && (r[1] >= 0.0) && (r[1] <= 1.0);                /* Ca2+ */
        result = result && (1.0-(r[0]+1.0)/2.0-r[1]+5.0*r[1]/7.0 >= 0.0)  /* Mg2+ */
                                        && (1.0-(r[0]+1.0)/2.0-r[1]+5.0*r[1]/7.0 <= 1.0-2.0*r[1]/7.0);
    }
    /* Check bounds on moles of endmember components */
    if (mask & SIXTH) {
        for (i=0, sum=0.0; i<NA; i++) sum += m[i];
        result = result && (sum >= 0.0);
        if (sum > 0.0) {
            result = result && (m[1]/sum >= 0.0)
                                            && (m[1]/sum <= 1.0-2.0*m[2]/(7.0*sum));        /* Fe2+ */
            result = result && (m[2]/sum >= 0.0) && (m[2]/sum <= 1.0);      /* Ca2+ */
            result = result && (m[0]/sum+5.0*m[2]/(7.0*sum) >= 0.0)         /* Mg2+ */
                                            && (m[0]/sum+5.0*m[2]/(7.0*sum) <= 1.0-2.0*m[2]/(7.0*sum));
        }
    }

    return result;
}

void
conOph(int inpMask, int outMask, double t, double p,
    double *e,      /* comp of amphibole in moles of elements                   */
    double *m,      /* comp of amphibole in moles of endmember components       */
    double *r,      /* comp of amphibole in terms of the independent comp var   */
    double *x,      /* comp of amphibole in mole fractions of endmember comp    */
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
            endmember amphibole components.
    (2) calculates from a vector of moles of endmember components, one or
            all of: r[], x[], dr[]/dm[] d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
    (3) calculates from a vector of independent compositional variables
            mole fractions of endmember components and/or the Jacobian matrix
            dx[]/dr[]

    In this routine it is assumed that the elements are in the order of atomic
    numbers and that the order of amphibole components has been verified as:
            m[0] = cummingtonite  (Mg7Si8O22(OH)2),
            m[1] = grunerite      (Fe7Si8O22(OH)2),
            m[2] = tremolite      (Ca2Mg5Si8O22(OH)2),

    ----------------------------------------------------------------------------*/

    int i, j, k;

    if (inpMask == FIRST && outMask == SECOND) {
        /* Converts a vector of moles of elements into a vector of moles of
       end-member components.                                                 */
        double sumcat, sumchg, fe2, fe3;
        static const int Hy =  1;
        static const int O  =  8;
        static const int F  =  9;
        static const int Na = 11;
        static const int Mg = 12;
        static const int Al = 13;
        static const int Si = 14;
        static const int Cl = 17;
        static const int Ca = 20;
        static const int Ti = 22;
        static const int Mn = 25;
        static const int Fe = 26;

        /* Sum the cations and correct the analysis for silica deficiency */
        sumcat  = e[Na] + e[Mg] + e[Al] + e[Si] + e[Ca] + e[Ti] + e[Mn] + e[Fe];
        sumchg  = e[Na] + 2.0*e[Mg] + 3.0*e[Al] + 4.0*e[Si] + 2.0*e[Ca] + 4.0*e[Ti]
                        + 2.0*e[Mn];

        /* Compute the ferric/ferrous ratio */
        fe3 = 2.0*(e[O]-e[Hy]) + e[Hy] + e[F] + e[Cl] - sumchg - 2.0*e[Fe];
        fe2 = e[Fe] - fe3;

/*
        if (fe3 < 0.01*e[Fe]) { fe3 = 0.01*e[Fe]; fe2 = 0.99*e[Fe]; }
        if (fe2 < 0.01*e[Fe]) { fe2 = 0.01*e[Fe]; fe3 = 0.99*e[Fe]; }
*/
        if (fe3 < 0.0) { fe3 = 0.01*e[Fe]; fe2 = 0.99*e[Fe]; }
        if (fe2 < 0.0) { fe2 = 0.01*e[Fe]; fe3 = 0.99*e[Fe]; }

        /* Assign moles of endmembers */
        m[0] = (e[Mg]  - 5.0*e[Ca]/2.0)/7.0;
        m[1] =  fe2/7.0;
        m[2] =  e[Ca]/2.0;

    } else if (inpMask == SECOND) {
        double sum;

        if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
            printf("Illegal call to conOph with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        for (i=0, sum=0.0; i<NA; i++) sum += m[i];

        if (outMask & THIRD) {
            /* Converts a vector of moles of end-member components (m) into a vector
         of independent compositional variables (r) required as input for the
         remaining public functions.                                          */
            r[0] = (sum != 0.0) ? 2.0*m[1]/sum - 1.0 : 0.0;
            r[1] = (sum != 0.0) ?     m[2]/sum : 0.0;
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
                    dm[0][j]  = -2.0*m[1]/SQUARE(sum);
                    dm[0][j] += (j == 1) ? 2.0/sum : 0.0;
                    dm[1][j]  = -(m[2])/SQUARE(sum);
                    dm[1][j] += (j == 2) ? 1.0/sum : 0.0;
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
                        d2m[0][j][k]  = 4.0*m[1]/CUBE(sum);
                        d2m[0][j][k] -= (j == 1) ? 2.0/SQUARE(sum) : 0.0;
                        d2m[0][j][k] -= (k == 1) ? 2.0/SQUARE(sum) : 0.0;

                        d2m[1][j][k]  = 2.0*m[2]/CUBE(sum);
                        d2m[1][j][k] -= (j == 2) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[1][j][k] -= (k == 2) ? 1.0/SQUARE(sum) : 0.0;
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
                            d3m[0][j][k][l]  = -12.0*m[1]/QUARTIC(sum);
                            d3m[0][j][k][l] += (j == 1) ? 4.0/CUBE(sum) : 0.0;
                            d3m[0][j][k][l] += (k == 1) ? 4.0/CUBE(sum) : 0.0;
                            d3m[0][j][k][l] += (l == 1) ? 4.0/CUBE(sum) : 0.0;

                            d3m[1][j][k][l]  = -6.0*m[2]/QUARTIC(sum);
                            d3m[1][j][k][l] += (j == 2) ? 2.0/CUBE(sum) : 0.0;
                            d3m[1][j][k][l] += (k == 2) ? 2.0/CUBE(sum) : 0.0;
                            d3m[1][j][k][l] += (l == 2) ? 2.0/CUBE(sum) : 0.0;
                        }
                    }
                }
            }
        }

    } else if (inpMask == THIRD) {

        if (outMask & ~(FOURTH | SEVENTH))
            printf("Illegal call to conOph with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        if (outMask & FOURTH) {
            /* Converts a vector of independent compositional variables (r) into a
         vector of mole fractions of endmember components (x).                */
            x[0] = 1.0 - (1.0+r[0])/2.0 - r[1];
            x[1] = (1.0+r[0])/2.0;
            x[2] = r[1];
        }

        if (outMask & SEVENTH) {
            /* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
            for (i=0; i<NA; i++) for (j=0; j<NR; j++) dr[i][j] = 0.0;
            dr[0][0] = -1.0/2.0; dr[0][1] = -1.0;
            dr[1][0] =  1.0/2.0;
            dr[2][1] =  1.0;
        }

    } else  {
        printf("Illegal call to conOph with inpMask = %o and outMask = %o\n",
            inpMask, outMask);
    }

}

/**************************************************************************/
/* This routine is not thread-safe and should be synchronized when called */
/**************************************************************************/

void
dispOph(int mask, double t, double p, double *x,
    char **formula            /* Mineral formula for interface display MASK: 1 */
    )
{
    double *r = x;
    static char masterString[] = {
/*             1111111111222222222233333333334444444444555555555566666666667
     01234567890123456789012345678901234567890123456789012345678901234567890 */
#ifdef PHMELTS_ADJUSTMENTS
        "     Ca_.__Fe_.__Mg_.__Si8O22(OH)2" };
#else
        "_Aph Ca_.__Fe_.__Mg_.__Si8O22(OH)2" };
#endif

    if (mask & FIRST) {
        char *string = strcpy((char *) malloc((size_t) (strlen(masterString)+1)*sizeof(char)), masterString);
        double totCa, totFe2, totMg;
        char n[5];
        int i;
#ifndef PHMELTS_ADJUSTMENTS
        if (isClino(t, p, r) == TRUE) string[0] = 'c'; else string[0] = 'o';
#endif

        totCa  = 2.0*r[1];
        totFe2 = 7.0*(r[0]+1.0)/2.0;
        totMg  = (1.0 - (r[0]+1.0)/2.0 - r[1])*7.0 +  5.0*r[1];

        (void) snprintf(n, 5, "%4.2f", totCa);  for (i=0; i<4; i++) string[ 7+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totFe2); for (i=0; i<4; i++) string[13+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totMg);  for (i=0; i<4; i++) string[19+i] = n[i];

        *formula = string;
    }
}

void
actOph(int mask, double t, double p, double *x,
    double *a,  /* (pointer to a[]) activities              BINARY MASK: 0001 */
    double *mu, /* (pointer to mu[]) chemical potentials    BINARY MASK: 0010 */
    double **dx /* (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100 */
    )           /* exclusion criteria applied to results if BINARY MASK: 1000 */
{
    double *r = x;
    double s[NS], g, dgdr[NR], fr[NA][NR];
    DECLARE_SITE_FRACTIONS
    int i, j, clino;

    for(i=0; i<NA; i++) {
     fr[i][0] = FR(i); /* R */
     fr[i][1] = FP(i); /* P */
    }

    clino = isClino(t, p, r);
    setClino(clino);

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

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
        double d2gdr2[NR][NR], d2gdrds[NR][NS], d2gds2[NS][NS], dsdr[NS][NR],
            dfrdr[NA][NR], gs[NA][NS], dgsds[NA][NS], sum;
        int k, l;

        fillD2GDR2
        fillD2GDRDS
        fillD2GDS2

        for(i=0; i<NA; i++) {
       gs[i][0] = FS13(i);         /* s13 */
       gs[i][1] = FS2(i);          /* s2  */
       dfrdr[i][0] = DFRDR(i);     /* R   */
       dfrdr[i][1] = DFPDP(i);     /* P   */
       dgsds[i][0] = DFS13DS13(i); /* s13 */
       dgsds[i][1] = DFS2DS2(i);   /* s2  */
        }

        order(SECOND, t, p, r,
                NULL, dsdr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NA; i++) {
            for (k=0; k<NR; k++) {
                /* compute activity of the i-th component */
                for (dx[i][k]=g, j=0; j<NR; j++) dx[i][k] += fr[i][j]*dgdr[j];
                dx[i][k] = exp(dx[i][k]/(R*t));

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
        int flags[NA];

        flags[0] = (xmg2m4 < 0.05) | (xmg2m13 < 0.05) | (xmg2m2 < 0.05);
        flags[1] = (xfe2m4 < 0.05) | (xfe2m13 < 0.05) | (xfe2m2 < 0.05);
        flags[2] = (xca2m4 < 0.05) | (xmg2m13 < 0.05) | (xmg2m2 < 0.05);

        for (i=0; i<NA; i++) {
            if (flags[i]) {
                if (mask & FIRST)  a[i]  = 0.0;
                if (mask & SECOND) mu[i] = 0.0;
                if (mask & THIRD)  for (j=0; j<NR; j++) dx[i][j] = 0.0;
            }
        }

    }

}

void
gmixOph(int mask, double t, double p, double *x,
    double *gmix, /* Gibbs energy of mixing             BINARY MASK: 0001 */
    double *dx,   /* (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010 */
    double **dx2, /* (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100 */
    double ***dx3 /* (pointer to dx3[][][]) d3(g)/d(x[])3 NARY MASK: 1000 */
    )
{
    double *r = x;
    double s[NS];
    DECLARE_SITE_FRACTIONS
    int clino;

    clino = isClino(t, p, r);
    setClino(clino);

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    if (mask & FIRST) {
        *gmix = G;
    }

    if(mask & SECOND) {
        dx[0] = DGDR0;
        dx[1] = DGDR1;
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
hmixOph(int mask, double t, double p, double *x,
    double *hmix /* Enthalpy of mixing BINARY MASK: 1 */
    )
{
    double *r = x;
    double s[NS];
    DECLARE_SITE_FRACTIONS
    int clino;

    clino = isClino(t, p, r);
    setClino(clino);

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    *hmix = (G) + t*(S);
}

void
smixOph(int mask, double t, double p, double *x,
    double *smix, /* Entropy of mixing                  BINARY MASK: 001 */
    double *dx,   /* (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010 */
    double **dx2  /* (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100 */
    )
{
    double *r = x;
    double s[NS];
    DECLARE_SITE_FRACTIONS
    int clino;

    clino = isClino(t, p, r);
    setClino(clino);

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    if (mask & FIRST) {
        *smix = S;
    }

    if(mask & SECOND) {
        double d2gdrds[NR][NS], d2gdrdt[NR], d2gds2[NS][NS], d2gdsdt[NS],
            dsdr[NS][NR], dsdt[NS];
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
                       + d2gds2[k][l]*dsdr[l][j]*d2sdrdt[k][i];
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
cpmixOph(int mask, double t, double p, double *x,
    double *cpmix, /* Heat capacity of mixing               BINARY MASK: 001 */
    double *dt,    /* d(cp)/d(t)                            BINARY MASK: 010 */
    double *dx     /* d(cp)/d(x[])d(t)                      BINARY MASK: 100 */
    )
{
    double *r = x;
    double s[NS], dsdt[NS], d2gdsdt[NS], d2gds2[NS][NS], d2gdt2;
    DECLARE_SITE_FRACTIONS
    int clino, i, j;

    clino = isClino(t, p, r);
    setClino(clino);

    order(FIRST | THIRD, t, p, r,
                s, NULL, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    fillD2GDS2
    fillD2GDSDT
    d2gdt2 = D2GDT2;

    if (mask & FIRST) {
        *cpmix = d2gdt2;
        for (i=0; i<NS; i++) {
            *cpmix += 2.0*d2gdsdt[i]*dsdt[i];
            for (j=0; j<NS; j++) *cpmix += d2gds2[i][j]*dsdt[i]*dsdt[j];
        }
        *cpmix *= -t;
    }

    if(mask & SECOND) {
        double d3gds3[NS][NS][NS], d3gds2dt[NS][NS], d3gdsdt2[NS], d2sdt2[NS],
            temp;
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
    }

    if(mask & THIRD) {
        double d3gds3[NS][NS][NS], d3gdrds2[NR][NS][NS], d3gdrdsdt[NR][NS],
            d3gds2dt[NS][NS], d2gdrds[NR][NS], d3gdrdt2[NR], d3gdsdt2[NS],
            dsdr[NS][NR], d2sdrdt[NS][NR], d2sdt2[NS];
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
    }

}

void
vmixOph(int mask, double t, double p, double *x,
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
    DECLARE_SITE_FRACTIONS
    int clino;

    clino = isClino(t, p, r);
    setClino(clino);

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    if (mask & FIRST) {
        *vmix = DGDP;
    }

    if(mask & SECOND) {
        double d2gdrds[NR][NS], d2gdrdp[NR], d2gds2[NS][NS], d2gdsdp[NS],
            dsdr[NS][NR], dsdp[NS];
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
                    NULL, dsdr, NULL, dsdp, d2sdr2, NULL, d2sdrdp, NULL, NULL, NULL);

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
        double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], dsdt[NS], dsdp[NS];
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
    }

    if(mask & FIFTH) {
        double d2gdp2 = D2GDP2;
        double d2gds2[NS][NS], d2gdsdp[NS], dsdp[NS];
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
    }

    if(mask & SIXTH) {
        double d3gdt2dp = D3GDT2DP;
        double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS][NS][NS],
            d3gds2dp[NS][NS], d3gds2dt[NS][NS], d3gdsdtdp[NS], d3gdsdt2[NS],
            dsdt[NS], dsdp[NS], d2sdt2[NS], d2sdtdp[NS];
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
    }

    if(mask & SEVENTH) {
        double d3gdtdp2 = D3GDTDP2;
        double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS][NS][NS],
            d3gds2dt[NS][NS], d3gds2dp[NS][NS], d3gdsdtdp[NS], d3gdsdp2[NS],
            dsdt[NS], dsdp[NS], d2sdtdp[NS], d2sdp2[NS];
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
    }

    if(mask & EIGHTH) {
        double d3gdp3 = D3GDP3;
        double d2gds2[NS][NS], d2gdsdp[NS], d3gds3[NS][NS][NS], d3gds2dp[NS][NS],
            d3gdsdp2[NS], dsdp[NS], d2sdp2[NS];
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
    }

    if(mask & NINTH) {
        double d3gds3[NS][NS][NS], d3gdrds2[NR][NS][NS], d3gdrdsdt[NR][NS],
            d3gds2dp[NS][NS], d2gdrds[NR][NS], d3gdrdtdp[NR], d3gdsdtdp[NS],
            dsdt[NS], dsdp[NS], dsdr[NS][NR], d2sdrdt[NS][NR], d2sdrdp[NS][NR],
            d2gds2[NS][NS], d2gdsdt[NS], d3gdrdsdp[NR][NS], d2gdsdp[NS],
            d2sdtdp[NS], d3gds2dt[NS][NS];
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
    }

    if(mask & TENTH) {
        double d3gds3[NS][NS][NS], d3gdrds2[NR][NS][NS], d3gds2dp[NS][NS],
            d2gdrds[NR][NS], dsdp[NS], dsdr[NS][NR], d2sdrdp[NS][NR], d2gds2[NS][NS],
            d3gdrdsdp[NR][NS], d3gdrdp2[NR], d3gdsdp2[NS], d2gdsdp[NS], d2sdp2[NS];
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
    }

}

/* end of file ORTHOAMPHIBOLE.C */
