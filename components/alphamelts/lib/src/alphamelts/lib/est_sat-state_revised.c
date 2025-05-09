const char *est_satState_revised_ver(void) { return "January 25, 2011"; }

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Routines to estimate the composition and chemical affinity for
**      saturation with a solid solution at a particular T and P.
**      (file: EST_SATSTATE_REVISED.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  January 25, 2011  Original Version
**--
*/

#include "silmin.h"

#ifdef DEBUG
#undef DEBUG
#endif

#define SCALE   1000.0  /* Scaling factor for chemical affinities           */

#define SUCCESS TRUE
#define FAILURE FALSE

/***************************************************************************
 * Consider an ideal solution of n+1 endmembers:
 * - A + RT ln x[1] + RT ln g[1]                    = - mu[1],
 * - A + RT ln x[2] + RT ln g[2]                    = - mu[2],
 *   ...
 * - A + RT ln x[n] + RT ln g[n]                    = - mu[n],
 * - A + RT ln(1-x[1]-x[2]-...-x[n]) + RT ln g[n+1] = - mu[n+1].
 *
 * Setting: f[1] = exp(-(mu[1]-RTlng[1]-mu[n]+RTlng[n])/(R*t)),
 *          f[2] = exp(-(mu[2]-RTlng[2]-mu[n]+RTlng[n])/(R*t)),
 *          ...
 *          f[n] = exp(-(mu[n]-RTlng[n]-mu[n+1]+RTlng[n+1])/(R*t))
 *
 * A solution is given by:
 * x[n]   = f[n]/(1+f[n]*(1+f[1]+f[2]+...+f[n-1])
 * x[n-1] = f[n-1]*x[n]
 * x[n-2] = f[n-2]*x[n]
 * ....
 * x[2]   = f[2]*x[n]
 * x[1]   = f[1]*x[n]
 *
 * Any of the original equations may be used to compute A
 ***************************************************************************/

int getAffinityAndCompositionGeneric( /* Returns a MODE flag for success or failure */
    double t,              /* temperature (K)					*/
    double p,              /* pressure (bars)					*/
    int index,             /* index of solid phase in the solids[] structure	*/
                         /* -1 indicates liquid 				*/
    int    *zeroX,         /* TRUE if endmember component has zero mole fraction  */
    double *muMinusMu0,    /* vector of end-member mu - mu0 i.e. A + RTln(a)	*/
    double *affinity,      /* returned value, chemical affinity (J)		*/
    double *indepVar)      /* returned vector, composition of phase (length nr)	*/
{
    int i, j, iter = 0, foundSolution = FALSE;
    static double *moleFrac = NULL, *bVec = NULL, *gVec = NULL, *activity = NULL, *mu = NULL;
    static int    *nullComp = NULL, *nullList = NULL;
    int    hasNull, na, nr, nz, solidID = -1, liquidMode;

    if (activity == NULL) {
        activity  = (double *)  malloc((unsigned) nlc*sizeof (double));
        bVec      = (double *)  malloc((unsigned) nlc*sizeof (double));
        gVec      = (double *)  malloc((unsigned) nlc*sizeof (double));
        moleFrac  = (double *)  malloc((unsigned) nlc*sizeof (double));
        mu        = (double *)  malloc((unsigned) nlc*sizeof (double));
        nullComp  = (int *)     malloc((unsigned) nlc*sizeof (int));
        nullList  = (int *)     malloc((unsigned) nlc*sizeof (int));
    }

    /* Test input parameters */
    if (t <= 0.0) 		  return FAILURE;
    if (p <  0.0) 		  return FAILURE;
    if (index < -1 || index >= npc) return FAILURE;

    /* Check parameters for algorithmic assumptions */
    liquidMode = (index == -1);
    if (!liquidMode) {
        na = solids[index].na;
        nr = solids[index].nr;
        if (na != (nr+1)) return FAILURE;
        solidID = index;
    } else {
        na = nlc;
        nr = nlc-1;
    }

    /* remove excluded components for computation of mole fraction estimates.
     the static nz is initialized here.                                    */
    for (i=0, hasNull=FALSE, nz=0; i<na; i++) {
        nullComp[i] = zeroX[i]; hasNull |= zeroX[i];
        if (!nullComp[i]) { nullList[nz] = i; mu[nz++] = muMinusMu0[i]; }
    }
    if (nz == 0) return FAILURE;

    for (i=0; i<na; i++) gVec[i] = 1.0;
    *affinity = 100000.0;

    while (!foundSolution && (iter < 100)) {
        double sum;

        /* take appropriate action for the number of non-zero components */
         if (nz == 0)                      bVec[na-1] = 0.0;
        else if (nz == 1) { moleFrac[0] = 1.0; bVec[na-1] = (mu[0]-gVec[0])/SCALE; }
        else {
            /* condense the activity cofficient vector */
            for (i=0, j=0; i<na; i++) if (!nullComp[i]) gVec[j++] = R*t*log(gVec[i]);

            /* Compute the f[n] terms and store them temporarily in moleFrac[] */
            sum = 1.0;
            if (nz > 2) for (i=0; i<(nz-2); i++) {
                moleFrac[i] = exp(-(mu[i]+gVec[i]-mu[nz-2]-gVec[nz-2])/(R*t)); sum += moleFrac[i];
            }
            moleFrac[nz-2] = exp(-(mu[nz-2]+gVec[nz-2]-mu[nz-1]-gVec[nz-1])/(R*t));

            /* Solve for the composition variables (mole fractions) */
            moleFrac[nz-2] /= 1.0 + moleFrac[nz-2]*sum;
            moleFrac[nz-1]  = 1.0 - moleFrac[nz-2];
            if (nz > 2) for (i=0; i<(nz-2); i++) {
                moleFrac[i] *= moleFrac[nz-2]; moleFrac[nz-1] -= moleFrac[i];
            }
            /* This fix is taken from Thermoengine c246b6c56ebddfd170fc2b36262fa4abb9327bad */
            for (i=0; i<nz; i++) if (moleFrac[i] < DBL_EPSILON) moleFrac[i] = DBL_EPSILON;

            /* compute the chemical affinity (choice of mu[] is arbitrary) */
            bVec[na-1] = (mu[0] + R*t*log(moleFrac[0]) + gVec[0])/SCALE;
        }

        /* Reassemble the mole fraction and chemical potential vectors with zeros for the absent endmembers */
        if (hasNull) for (i=na-1, j=nz; i>=0; i--) {
            if(!nullComp[i]) moleFrac[i] = moleFrac[--j];
            else             moleFrac[i] = 0.0;
        }

        /* convert mole fractions of endmembers into independent compos var */
        if (!liquidMode) (*solids[solidID].convert)(SECOND, THIRD, t, p, NULL, moleFrac, bVec, NULL, NULL, NULL, NULL, NULL);
        else                                 conLiq(SECOND, THIRD, t, p, NULL, moleFrac, bVec, NULL, NULL, NULL, NULL);

        if (fabs(*affinity-bVec[nr]) < 0.1/SCALE) foundSolution = TRUE;
        else {
            if (!liquidMode) (*solids[solidID].activity)(FIRST, t, p, bVec, activity, NULL, NULL);
            else                                  actLiq(FIRST, t, p, bVec, activity, NULL, NULL, NULL);
            for (i=0; i<na; i++) if (moleFrac[i] != 0.0) gVec[i] = activity[i]/moleFrac[i];
        }
        *affinity = bVec[nr];
#ifdef DEBUG
        printf("SS iter = %2.2d X,g:", iter);
        for (i=0; i<na; i++) printf(" %13.6g (%13.6g)", moleFrac[i], gVec[i]);
        printf(" A: %13.6g\n", bVec[nr]);
#endif
        iter++;
    }
    for (i=0; i<nr; i++) indepVar[i] = bVec[i];
    *affinity = bVec[nr]*SCALE;

    return (iter < 100) ? SUCCESS : FAILURE;
}

#define  DH1      0.5573     * 1000.0 * 4.184 /* joules     */
#define  DS1      0.00033    * 1000.0 * 4.184 /* joules/K   */
#define  DV1      0.018               * 4.184 /* joules/bar */
#define  DH2      1.200      * 1000.0 * 4.184 /* joules     */
#define  DS2      0.000675   * 1000.0 * 4.184 /* joules/K   */
#define  DV2      0.0148              * 4.184 /* joules/bar */
#define  DH3      0.6991252  * 1000.0 * 4.184 /* joules     */
#define  DS3      0.0005926  * 1000.0 * 4.184 /* joules/K   */
#define  DV3      0.00865             * 4.184 /* joules/bar */
#define  DH4     -0.49       * 1000.0 * 4.184 /* joules     */
#define  DS4     -0.000132   * 1000.0 * 4.184 /* joules/K   */

#define cV0DI     6.620                       /* joules/bar */  /* reference */
#define cV0HD     6.7894                      /* joules/bar */  /* reference */
#define oV0HD     1.64620             * 4.184 /* joules/bar */
#define oV0EN     3.133 * 2.0                 /* joules/bar */  /* reference */
#define oV0FS     3.296 * 2.0                 /* joules/bar */  /* reference */

#define oH027    -3.30       * 1000.0 * 4.184 /* joules     */
#define oS027    -0.00055428 * 1000.0 * 4.184 /* joules/K   */
#define pH027    -1.97547480 * 1000.0 * 4.184 /* joules     */
#define pS027     0.00071134 * 1000.0 * 4.184 /* joules/K   */
#define pV027    -0.00920929          * 4.184 /* joules/bar */

#define oHEX     -1.87       * 1000.0 * 4.184 /* joules     */
#define oVEX     -0.029               * 4.184 /* joules/bar */
#define cHEX     -2.2        * 1000.0 * 4.184 /* joules     */
#define cVEX      0.01                * 4.184 /* joules/bar */
#define pHEX     -0.65       * 1000.0 * 4.184 /* joules     */
#define pVEX      0.0                 * 4.184 /* joules/bar */

#define oHX      -0.45       * 1000.0 * 4.184 /* joules     */
#define oVX       0.00675             * 4.184 /* joules/bar */
#define cHX      -0.10       * 1000.0 * 4.184 /* joules     */
#define cVX       0.0                 * 4.184 /* joules/bar */
#define pHX      -0.10       * 1000.0 * 4.184 /* joules     */
#define pVX       0.0                 * 4.184 /* joules/bar */

#define oWHFEMG   2.0        * 1000.0 * 4.184 /* joules     */  /* W Fe-Mg M2 */
#define oWVFEMG   0.003375            * 4.184 /* joules/bar */
#define cWHFEMG   1.15       * 1000.0 * 4.184 /* joules     */
#define cWVFEMG  -0.0035              * 4.184 /* joules/bar */

#define oWH12     2.0        * 1000.0 * 4.184 /* joules     */  /* W12 */
#define oWV12     0.003375            * 4.184 /* joules/bar */
#define cWH12     1.68       * 1000.0 * 4.184 /* joules     */
#define cWV12     0.0                 * 4.184 /* joules/bar */

#define oWHCAMG   7.56       * 1000.0 * 4.184 /* joules     */  /* W17 */
#define oWVCAMG   0.008               * 4.184 /* joules/bar */
#define cWHCAMG   6.72       * 1000.0 * 4.184 /* joules     */
#define cWVCAMG  -0.009               * 4.184 /* joules/bar */

#define oWHCAFE   4.12       * 1000.0 * 4.184 /* joules     */  /* W17U */
#define oWVCAFE   0.011               * 4.184 /* joules/bar */
#define cWHCAFE   4.515      * 1000.0 * 4.184 /* joules     */
#define cWVCAFE   0.003               * 4.184 /* joules/bar */

#define oDWHCAMG  -1.3       * 1000.0 * 4.184 /* joules     */  /* delta W17 */
#define oDWVCAMG   0.012              * 4.184 /* joules/bar */
#define cDWHCAMG  -0.7       * 1000.0 * 4.184 /* joules     */
#define cDWVCAMG   0.008              * 4.184 /* joules/bar */

#define oDWHCAFE  -1.1       * 1000.0 * 4.184 /* joules     */  /* delta W17U */
#define oDWVCAFE   0.005              * 4.184 /* joules/bar */
#define cDWHCAFE  -0.5375    * 1000.0 * 4.184 /* joules     */
#define cDWVCAFE   0.0025             * 4.184 /* joules/bar */

#define  W13     16.318     * 1000.0 /* joules     */
#define  W14     18.82800   * 1000.0 /* joules     */
#define  W15     20.920     * 1000.0 /* joules     */
#define  W15P    16.78900   * 1000.0 /* joules     */
#define  W16      0.00000   * 1000.0 /* joules     */
#define  W25      3.98176   * 1000.0 /* joules     */
#define  W25P     7.49005   * 1000.0 /* joules     */
#define  W26      0.00000   * 1000.0 /* joules     */
#define  W34     16.78900   * 1000.0 /* joules     */
#define  W35     12.02812   * 1000.0 /* joules     */
#define  W35P    47.38600   * 1000.0 /* joules     */
#define  W36     34.35062   * 1000.0 /* joules     */
#define  cW37    28.65091   * 1000.0 /* joules     */
#define  cW3U7U  34.90070   * 1000.0 /* joules     */
#define  W45     27.14312   * 1000.0 /* joules     */
#define  W45P    16.318     * 1000.0 /* joules     */
#define  W46     40.40662   * 1000.0 /* joules     */
#define  cW47    35.38104   * 1000.0 /* joules     */
#define  cW4U7U  28.54527   * 1000.0 /* joules     */
#define  cW55    35.670     * 1000.0 /* joules     */
#define  W56      0.00000   * 1000.0 /* joules     */
#define  cW57    33.25291   * 1000.0 /* joules     */
#define  cW57U   15.29737   * 1000.0 /* joules     */
#define  W5P6     0.00000   * 1000.0 /* joules     */
#define  cW5P7   36.24989   * 1000.0 /* joules     */
#define  cW5P7U   2.35321   * 1000.0 /* joules     */
#define  cW67     0.00000   * 1000.0 /* joules     */
#define  cW67U    0.00000   * 1000.0 /* joules     */
#define  oW37    28.65091   * 1000.0 /* joules     */
#define  oW3U7U  34.90070   * 1000.0 /* joules     */
#define  oW47    35.38104   * 1000.0 /* joules     */
#define  oW4U7U  28.54527   * 1000.0 /* joules     */
#define  oW55    35.670     * 1000.0 /* joules     */
#define  oW57    33.25291   * 1000.0 /* joules     */
#define  oW57U   15.29737   * 1000.0 /* joules     */
#define  oW5P7   36.24989   * 1000.0 /* joules     */
#define  oW5P7U   2.35321   * 1000.0 /* joules     */
#define  oW67     0.00000   * 1000.0 /* joules     */
#define  oW67U    0.00000   * 1000.0 /* joules     */

#define  H23     -2.71700   * 1000.0 /* joules     */
#define  S23      0.0                /* joules/K   */
#define  H24     -7.36648   * 1000.0 /* joules     */
#define  S24      0.0                /* joules/K   */
#define  cH55     9.48524   * 1000.0 /* joules     */
#define  S55      0.0                /* joules/K   */
#define  oH55     9.48524   * 1000.0 /* joules     */

#define  DcTOoH3    25.11924    * 1000.0 /* joules     */
#define  DcTOoS3    -2.00000             /* joules/K   */
#define  DcTOoV3    -0.05129             /* joules/bar */
#define  DcTOoH4    25.11924    * 1000.0 /* joules     */
#define  DcTOoS4    -2.00000             /* joules/K   */
#define  DcTOoV4    -0.05129             /* joules/bar */
#define  DcTOoH5    25.11924    * 1000.0 /* joules     */
#define  DcTOoS5    -2.00000             /* joules/K   */
#define  DcTOoV5    -0.05129             /* joules/bar */
#define  DcTOoH6    25.11924    * 1000.0 /* joules     */
#define  DcTOoS6    -2.00000             /* joules/K   */
#define  DcTOoV6    -0.05129             /* joules/bar */

#define  DcTOpH1    -0.9        * 1000.0 * 4.184 /* joules     */
#define  DcTOpS1    -0.0006920  * 1000.0 * 4.184 /* joules/K   */
#define  DcTOpV1     0.0                 * 4.184 /* joules/bar */
#define  DcTOpH2    -1.88       * 1000.0 * 4.184 /* joules     */
#define  DcTOpS2    -0.00156401 * 1000.0 * 4.184 /* joules/K   */
#define  DcTOpV2     0.0                 * 4.184 /* joules/bar */
#define  DcTOpH3    21.11924    * 1000.0 /* joules     */
#define  DcTOpS3    -2.00000             /* joules/K   */
#define  DcTOpV3    -0.05129             /* joules/bar */
#define  DcTOpH4    21.11924    * 1000.0 /* joules     */
#define  DcTOpS4    -2.00000             /* joules/K   */
#define  DcTOpV4    -0.05129             /* joules/bar */
#define  DcTOpH5    21.11924    * 1000.0 /* joules     */
#define  DcTOpS5    -2.00000             /* joules/K   */
#define  DcTOpV5    -0.05129             /* joules/bar */
#define  DcTOpHj    21.11924    * 1000.0 /* joules     */ /* Jadeite */
#define  DcTOpSj    -2.00000             /* joules/K   */ /* Jadeite */
#define  DcTOpVj    -0.05129             /* joules/bar */ /* Jadeite */
#define  DcTOpH7     0.0        * 1000.0 * 4.184 /* joules     */
#define  DcTOpS7     0.0        * 1000.0 * 4.184 /* joules/K   */
#define  DcTOpV7     0.0                 * 4.184 /* joules/bar */

/*
 *=============================================================================
 * Dependent parameters (DO NOT change definitions below this line:
 */

#define oV0DI    (cV0DI)+(DV1)
#define cV0EN    (oV0EN)+(DV2)
#define cV0FS    (oV0FS)+(DV3)

#define DV4      (cV0HD)-(oV0HD)

#define oV027    (oV0FS)-(oV0EN)+2.0*(oV0DI)-2.0*(oV0HD)
#define cV027    (cV0FS)-(cV0EN)+2.0*(cV0DI)-2.0*(cV0HD)
#define cH027    (oH027)+(DH3)-(DH2)-2.0*(DH1)-2.0*(DH4)
#define cS027    (oS027)+(DS3)-(DS2)-2.0*(DS1)-2.0*(DS4)

#define  DG1     (DH1)-t*(DS1)+(p-1.0)*DV1
#define  DG2     (DH2)-t*(DS2)+(p-1.0)*DV2
#define  DG3     (DH3)-t*(DS3)+(p-1.0)*DV3
#define  DG4     (DH4)-t*(DS4)+(p-1.0)*DV4

#define cG027    (cH027)-t*(cS027)+(p-1.0)*(cV027)
#define pG027    (pH027)-t*(pS027)+(p-1.0)*(pV027)
#define oG027    (oH027)-t*(oS027)+(p-1.0)*(oV027)
#define  H027    (clino) ? (cH027) : (oH027)
#define  S027    (clino) ? (cS027) : (oS027)
#define  V027    (clino) ? (cV027) : (oV027)
#define  G027    (clino) ? (cG027) : (oG027)
#define dH027    (clino) ? (pH027)-(cH027) : 0.0
#define dS027    (clino) ? (pS027)-(cS027) : 0.0
#define dV027    (clino) ? (pV027)-(cV027) : 0.0
#define dG027    (clino) ? (pG027)-(cG027) : 0.0

#define cGEX     (cHEX)+(p-1.0)*(cVEX)
#define pGEX     (pHEX)+(p-1.0)*(pVEX)
#define oGEX     (oHEX)+(p-1.0)*(oVEX)
#define  HEX     (clino) ? (cHEX) : (oHEX)
#define  VEX     (clino) ? (cVEX) : (oVEX)
#define  GEX     (clino) ? (cGEX) : (oGEX)
#define dHEX     (clino) ? (pHEX)-(cHEX) : 0.0
#define dVEX     (clino) ? (pVEX)-(cVEX) : 0.0
#define dGEX     (clino) ? (pGEX)-(cGEX) : 0.0

#define cGX      (cHX)+(p-1.0)*(cVX)
#define pGX      (pHX)+(p-1.0)*(pVX)
#define oGX      (oHX)+(p-1.0)*(oVX)
#define  HX      (clino) ? (cHX) : (oHX)
#define  VX      (clino) ? (cVX) : (oVX)
#define  GX      (clino) ? (cGX) : (oGX)
#define dHX      (clino) ? (pHX)-(cHX) : 0.0
#define dVX      (clino) ? (pVX)-(cVX) : 0.0
#define dGX      (clino) ? (pGX)-(cGX) : 0.0

#define oWFEMG   (oWHFEMG)+(p-1.0)*(oWVFEMG)
#define cWFEMG   (cWHFEMG)+(p-1.0)*(cWVFEMG)
#define  WHFEMG  (clino) ? (cWHFEMG) : (oWHFEMG)
#define  WVFEMG  (clino) ? (cWVFEMG) : (oWVFEMG)
#define  WFEMG   (clino) ? (cWFEMG) : (oWFEMG)

#define oW12     (oWH12)+(p-1.0)*(oWV12)
#define cW12     (cWH12)+(p-1.0)*(cWV12)
#define  WH12    (clino) ? (cWH12) : (oWH12)
#define  WV12    (clino) ? (cWV12) : (oWV12)
#define  W12     (clino) ? (cW12) : (oW12)

#define oWCAMG   (oWHCAMG)+(p-1.0)*(oWVCAMG)
#define cWCAMG   (cWHCAMG)+(p-1.0)*(cWVCAMG)
#define  WHCAMG  (clino) ? (cWHCAMG) : (oWHCAMG)
#define  WVCAMG  (clino) ? (cWVCAMG) : (oWVCAMG)
#define  WCAMG   (clino) ? (cWCAMG) : (oWCAMG)

#define oDWCAMG   (oDWHCAMG)+(p-1.0)*(oDWVCAMG)
#define cDWCAMG   (cDWHCAMG)+(p-1.0)*(cDWVCAMG)
#define  DWHCAMG  (clino) ? (cDWHCAMG) : (oDWHCAMG)
#define  DWVCAMG  (clino) ? (cDWVCAMG) : (oDWVCAMG)
#define  DWCAMG   (clino) ? (cDWCAMG) : (oDWCAMG)

#define oWCAFE   (oWHCAFE)+(p-1.0)*(oWVCAFE)
#define cWCAFE   (cWHCAFE)+(p-1.0)*(cWVCAFE)
#define  WHCAFE  (clino) ? (cWHCAFE) : (oWHCAFE)
#define  WVCAFE  (clino) ? (cWVCAFE) : (oWVCAFE)
#define  WCAFE   (clino) ? (cWCAFE) : (oWCAFE)

#define oDWCAFE   (oDWHCAFE)+(p-1.0)*(oDWVCAFE)
#define cDWCAFE   (cDWHCAFE)+(p-1.0)*(cDWVCAFE)
#define  DWHCAFE  (clino) ? (cDWHCAFE) : (oDWHCAFE)
#define  DWVCAFE  (clino) ? (cDWVCAFE) : (oDWVCAFE)
#define  DWCAFE   (clino) ? (cDWCAFE) : (oDWCAFE)

#define  W37      (clino) ? (cW37)   : (oW37)
#define  W3U7U    (clino) ? (cW3U7U) : (oW3U7U)
#define  W47      (clino) ? (cW47)   : (oW47)
#define  W4U7U    (clino) ? (cW4U7U) : (oW4U7U)
#define  W55      (clino) ? (cW55)   : (oW55)
#define  W57      (clino) ? (cW57)   : (oW57)
#define  W57U     (clino) ? (cW57U)  : (oW57U)
#define  W5P7     (clino) ? (cW5P7)  : (oW5P7)
#define  W5P7U    (clino) ? (cW5P7U) : (oW5P7U)
#define  W67      (clino) ? (cW67)   : (oW67)
#define  W67U     (clino) ? (cW67U)  : (oW67U)
#define  H55      (clino) ? (cH55)   : (oH55)

#define  G23     (H23)-t*(S23)
#define  G24     (H24)-t*(S24)
#define  G55     (H55)-t*(S55)

#define  DcTOoG3 (DcTOoH3)-t*(DcTOoS3)+(p-1.0)*(DcTOoV3)
#define  DcTOoG4 (DcTOoH4)-t*(DcTOoS4)+(p-1.0)*(DcTOoV4)
#define  DcTOoG5 (DcTOoH5)-t*(DcTOoS5)+(p-1.0)*(DcTOoV5)
#define  DcTOoG6 (DcTOoH6)-t*(DcTOoS6)+(p-1.0)*(DcTOoV6)

#define  DcTOpG1 (DcTOpH1)-t*(DcTOpS1)+(p-1.0)*(DcTOpV1)
#define  DcTOpG2 (DcTOpH2)-t*(DcTOpS2)+(p-1.0)*(DcTOpV2)
#define  DcTOpG3 (DcTOpH3)-t*(DcTOpS3)+(p-1.0)*(DcTOpV3)
#define  DcTOpG4 (DcTOpH4)-t*(DcTOpS4)+(p-1.0)*(DcTOpV4)
#define  DcTOpG5 (DcTOpH5)-t*(DcTOpS5)+(p-1.0)*(DcTOpV5)
#define  DcTOpGj (DcTOpHj)-t*(DcTOpSj)+(p-1.0)*(DcTOpVj)
#define  DcTOpG7 (DcTOpH7)-t*(DcTOpS7)+(p-1.0)*(DcTOpV7)

/*
 * Vertices of composition space
 */

/* (Ca)(Mg)(Si)2 O6 */
#define  H1      (clino) ? 0.0 :  (DH1)
#define  S1      (clino) ? 0.0 :  (DS1)
#define  V1      (clino) ? 0.0 :  (DV1)
#define  G1      (clino) ? 0.0 :  (DH1)-t*(DS1)+(p-1.0)*(DV1)
/* (Ca) (Fe2+) (Si)2 O6 */
#define  H2      (clino) ? 0.0 : -(DH4)
#define  S2      (clino) ? 0.0 : -(DS4)
#define  V2      (clino) ? 0.0 : -(DV4)
#define  G2      (clino) ? 0.0 : -(DH4)+t*(DS4)-(p-1.0)*(DV4)
/* (Ca) (Ti,Mg) (Al,Si)2 O6 */
#define  H3      (clino) ? 0.0 : (DcTOoH3)
#define  S3      (clino) ? 0.0 : (DcTOoS3)
#define  V3      (clino) ? 0.0 : (DcTOoV3)
#define  G3      (clino) ? 0.0 : (DcTOoG3)
/* (Ca) (Ti,Mg) (Fe3+,Si)2 O6 */
#define  H4      (clino) ? 0.0 : (DcTOoH4)
#define  S4      (clino) ? 0.0 : (DcTOoS4)
#define  V4      (clino) ? 0.0 : (DcTOoV4)
#define  G4      (clino) ? 0.0 : (DcTOoG4)
/* (Ca) (Fe3+,Al) [(Fe3+,Al),Si]2 O6 */
#define  H5      (clino) ? 0.0 : (DcTOoH5)
#define  S5      (clino) ? 0.0 : (DcTOoS5)
#define  V5      (clino) ? 0.0 : (DcTOoV5)
#define  G5      (clino) ? 0.0 : (DcTOoG5)
/* (Na) (Al) (Si)2 O6 */
#define  H6      (clino) ? 0.0 : (DcTOoH6)
#define  S6      (clino) ? 0.0 : (DcTOoS6)
#define  V6      (clino) ? 0.0 : (DcTOoV6)
#define  G6      (clino) ? 0.0 : (DcTOoG6)
/* (Mg) (Mg) (Si)2 O6 */
#define  H7      (clino) ? 0.0 : -(DH2)
#define  S7      (clino) ? 0.0 : -(DS2)
#define  V7      (clino) ? 0.0 : -(DV2)
#define  G7      (clino) ? 0.0 : -(DH2)+t*(DS2)-(p-1.0)*(DV2)
/* Note that ferrosilite (and DG3 terms) are dependent */

/* (Ca)(Mg)(Si)2 O6 */
#define  pH1     (clino) ? (DcTOpH1) : 0.0
#define  pS1     (clino) ? (DcTOpS1) : 0.0
#define  pV1     (clino) ? (DcTOpV1) : 0.0
#define  pG1     (clino) ? (DcTOpG1) : 0.0
/* (Ca) (Fe2+) (Si)2 O6 */
#define  pH2     (clino) ? (DcTOpH2) : 0.0
#define  pS2     (clino) ? (DcTOpS2) : 0.0
#define  pV2     (clino) ? (DcTOpV2) : 0.0
#define  pG2     (clino) ? (DcTOpG2) : 0.0
/* (Ca) (Ti,Mg) (Al,Si)2 O6 */
#define  pH3     (clino) ? (DcTOpH3) : 0.0
#define  pS3     (clino) ? (DcTOpS3) : 0.0
#define  pV3     (clino) ? (DcTOpV3) : 0.0
#define  pG3     (clino) ? (DcTOpG3) : 0.0
/* (Ca) (Ti,Mg) (Fe3+,Si)2 O6 */
#define  pH4     (clino) ? (DcTOpH4) : 0.0
#define  pS4     (clino) ? (DcTOpS4) : 0.0
#define  pV4     (clino) ? (DcTOpV4) : 0.0
#define  pG4     (clino) ? (DcTOpG4) : 0.0
/* (Ca) (Fe3+,Al) [(Fe3+,Al),Si]2 O6 */
#define  pH5     (clino) ? (DcTOpH5) : 0.0
#define  pS5     (clino) ? (DcTOpS5) : 0.0
#define  pV5     (clino) ? (DcTOpV5) : 0.0
#define  pG5     (clino) ? (DcTOpG5) : 0.0
/* (Na) (Al) (Si)2 O6 */
#define  pHj     (clino) ? (DcTOpHj) : 0.0
#define  pSj     (clino) ? (DcTOpSj) : 0.0
#define  pVj     (clino) ? (DcTOpVj) : 0.0
#define  pGj     (clino) ? (DcTOpGj) : 0.0
/* (Mg) (Mg) (Si)2 O6 */
#define  pH7     (clino) ? (DcTOpH7) : 0.0
#define  pS7     (clino) ? (DcTOpS7) : 0.0
#define  pV7     (clino) ? (DcTOpV7) : 0.0
#define  pG7     (clino) ? (DcTOpG7) : 0.0


int getAffinityAndCompositionPyroxene( /* Returns a MODE flag for success or failure */
    double t,              /* temperature (K)					*/
    double p,              /* pressure (bars)					*/
    int    index,          /* index of solid phase in the solids[] structure	*/
    int    *zeroX,         /* TRUE if endmember component has zero mole fraction  */
    double *muMinusMu0,    /* vector of end-member mu - mu0 i.e. A + RTln(a)	*/
    double *affinity,      /* returned value, chemical affinity (J)		*/
    double *indepVar)      /* returned vector, composition of phase (length nr)	*/
{
    int i, j, iter = 0, foundSolution = FALSE;
    static double *mF = NULL, *mFR = NULL, *bVec = NULL, *gVec = NULL, *activity = NULL, *mu = NULL, *mu0 = NULL, *muR = NULL;
    static int    *nullComp = NULL, *nullList = NULL;
    int    hasNull, na, nr, nz, solidID = -1, ns;

    /* Test input parameters */
    if (t <= 0.0) 		return FAILURE;
    if (p <  0.0) 		return FAILURE;
    if (index < 0 || index > npc) return FAILURE;

    /* Check parameters for algorithmic assumptions */
    na = solids[index].na;
    nr = solids[index].nr;
    ns = na + 7;
    solidID = index;

    if (activity == NULL) {
        activity  = (double *)  malloc((unsigned) ns*sizeof (double));
        bVec      = (double *)  malloc((unsigned) ns*sizeof (double));
        gVec      = (double *)  malloc((unsigned) ns*sizeof (double));
        mF        = (double *)  malloc((unsigned) ns*sizeof (double));
        mFR       = (double *)  malloc((unsigned) na*sizeof (double));
        mu        = (double *)  malloc((unsigned) ns*sizeof (double));
        mu0       = (double *)  malloc((unsigned) ns*sizeof (double));
        muR       = (double *)  malloc((unsigned) na*sizeof (double));
        nullComp  = (int *)     malloc((unsigned) ns*sizeof (int));
        nullList  = (int *)     malloc((unsigned) ns*sizeof (int));
    }

    for (i=0, hasNull=FALSE, nz=0; i<na; i++) {
        nullComp[i] = zeroX[i]; hasNull |= zeroX[i];
        if (!nullComp[i]) { nullList[nz] = i; mu[nz++] = muMinusMu0[i]; }
        mu0[i] = (solids[solidID+1+i].cur).g;
    }
    if (nz == 0) return FAILURE;

    /* compute dependent species mu0 for Ca(Fe,Ti)AlSiO6 */
    if (!nullComp[0] && !nullComp[2] && !nullComp[3]) {
        double correction = (G23);
        mu0[7] = mu0[3] + mu0[2]/2.0 - mu0[0]/2.0 + correction;
        nullComp[7] = FALSE;
        nullList[nz] = 7;
        mu[nz++] = muMinusMu0[3] + muMinusMu0[2]/2.0 - muMinusMu0[0]/2.0 + correction;
    } else { mu0[7] = 0.0; nullComp[7] = TRUE; }

    /* compute dependent species mu0 for Ca(Fe,Ti)FeSiO6 */
    if (!nullComp[0] && !nullComp[2] && !nullComp[4]) {
        double correction = (G24);
        mu0[8] = mu0[4] + mu0[2]/2.0 - mu0[0]/2.0 + correction;
        nullComp[8] = FALSE;
        nullList[nz] = 8;
        mu[nz++] = muMinusMu0[4] + muMinusMu0[2]/2.0 - muMinusMu0[0]/2.0 + correction;
    } else { mu0[8] = 0.0; nullComp[8] = TRUE; }

    /* compute dependent species mu0 for CaAl2SiO6 */
    if (!nullComp[3] && !nullComp[4] && !nullComp[5]) {
        static const int clino = TRUE;
        double correction = (G55) - (W34) - (W45P) + (W35P);
        mu0[9] = mu0[5] - mu0[4] + mu0[3] + correction;
        nullComp[9] = FALSE;
        nullList[nz] = 9;
        mu[nz++] = muMinusMu0[5] - muMinusMu0[4] + muMinusMu0[3] + correction;
    } else { mu0[9] = 0.0; nullComp[9] = TRUE; }

    /* compute dependent species mu0 for MgAl2SiO6 */
    if (!nullComp[0] && !nullComp[1] && !nullComp[3] && !nullComp[4] && !nullComp[5]) {
        static const int clino = TRUE;
        double correction = (pG5) - (pG4) + (pG3) - (pG1) + (pG7)
            + (G55) - (W47) + (W37) + (W5P7) - (W34) + (W14) - (W13) - (W45P) + (W35P) - (W15P) - (WCAMG) + 2.0*(DWCAMG) - t*0.25*(S027);
        mu0[10] = mu0[5] - mu0[4] + mu0[3] + mu0[1] - mu0[0] + correction;
        nullComp[10] = FALSE;
        nullList[nz] = 10;
        mu[nz++] = muMinusMu0[5] - muMinusMu0[4] + muMinusMu0[3] + muMinusMu0[1] - muMinusMu0[0] + correction;
    } else { mu0[10] = 0.0; nullComp[10] = TRUE; }

    /* compute dependent species mu0 for FeAl2SiO6 */
    if (!nullComp[0] && !nullComp[1] && !nullComp[2] && !nullComp[3] && !nullComp[4] && !nullComp[5]) {
        static const int clino = TRUE;
        double correction = (pG5) - (pG4) + (pG3) - 2.0*(pG1) + (pG7) + (pG2)
            + (G55) + (G24) - (G23) + (G027) + (dG027) - (W4U7U) + (W3U7U) + (W5P7U) - (W34) + (W14) - (W13) - (W45P)
            + (W35P) - (W25P) - (WCAFE) + 3.0*(DWCAFE)/2.0  + t*0.25*(S027);
        mu0[11] = mu0[5] - mu0[4] + mu0[3] + mu0[2] + mu0[1] - 2.0*mu0[0] + correction;
        nullComp[11] = FALSE;
        nullList[nz] = 11;
        mu[nz++] = muMinusMu0[5] - muMinusMu0[4] + muMinusMu0[3] + muMinusMu0[2] + muMinusMu0[1] - 2.0*muMinusMu0[0] + correction;
    } else { mu0[11] = 0.0; nullComp[11] = TRUE; }

    /* compute dependent species mu0 for NaFeSi2O6 */
    if (!nullComp[3] && !nullComp[4] && !nullComp[6]) {
        static const int clino = TRUE;
        double correction = - (G55) + (W56) - (W5P6) + (W46) - (W36) + (W45) - (W45P) + (W35P) - (W35) - (W34) - (W55);
        mu0[12] = mu0[6] + mu0[4] - mu0[3] + correction;
        nullComp[12] = FALSE;
        nullList[nz] = 12;
        mu[nz++] = muMinusMu0[6] + muMinusMu0[4] - muMinusMu0[3] + correction;
    } else { mu0[12] = 0.0; nullComp[12] = TRUE; }

    /* compute dependent species mu0 for Fe2Si2O6 */
    if (!nullComp[0] && !nullComp[1] && !nullComp[2]) {
        static const int clino = TRUE;
        double correction = (pG7) - 2.0*(pG1) + 2.0*(pG2) + (G027) + (dG027)  + t*0.25*(S027);
        mu0[13] = mu0[1] + 2.0*mu0[2] - 2.0*mu0[0] + correction;
        nullComp[13] = FALSE;
        nullList[nz] = 13;
        mu[nz++] = muMinusMu0[1] + 2.0*muMinusMu0[2] - 2.0*muMinusMu0[0] + correction;
    } else { mu0[13] = 0.0; nullComp[13] = TRUE; }

#ifdef DEBUG_EXTRA
    for (i=0, j=0; i<ns; i++) if (!nullComp[i]) printf("...mu0[%2.2d] = %13.6g mu = %13.6g\n", i, mu0[i], mu[j++]);
#endif

    for (i=0; i<ns; i++) gVec[i] = 1.0;
    *affinity = 100000.0;

    while (!foundSolution && (iter < 100)) {
        double sum;

        /* take appropriate action for the number of non-zero components */
         if (nz == 0)                      bVec[na-1] = 0.0;
        else if (nz == 1) { mF[0] = 1.0; bVec[na-1] = (mu[0]-gVec[0])/SCALE; }
        else {
            /* condense the activity cofficient vector */
            for (i=0, j=0; i<ns; i++) if (!nullComp[i]) gVec[j++] = R*t*log(gVec[i]);

            /* Compute the f[n] terms and store them temporarily in mF[] */
            sum = 1.0;
            if (nz > 2) for (i=0; i<(nz-2); i++) {
                mF[i] = exp(-(mu[i]+gVec[i]-mu[nz-2]-gVec[nz-2])/(R*t)); sum += mF[i];
            }
            mF[nz-2] = exp(-(mu[nz-2]+gVec[nz-2]-mu[nz-1]-gVec[nz-1])/(R*t));

            /* Solve for the composition variables (mole fractions) */
            mF[nz-2] /= 1.0 + mF[nz-2]*sum;
            mF[nz-1]  = 1.0 - mF[nz-2];
            if (nz > 2) for (i=0; i<(nz-2); i++) {
                mF[i] *= mF[nz-2]; mF[nz-1] -= mF[i];
            }
            /* This fix is taken from Thermoengine c246b6c56ebddfd170fc2b36262fa4abb9327bad */
            for (i=0; i<nz; i++) if (mF[i] < DBL_EPSILON) mF[i] = DBL_EPSILON;

            /* compute the chemical affinity (choice of mu[] is arbitrary) */
            bVec[na-1] = (mu[0] + R*t*log(mF[0]) + gVec[0])/SCALE;
        }

        /* Reassemble the mole fraction and chemical potential vectors with zeros for the absent endmembers */
        if (hasNull) for (i=ns-1, j=nz; i>=0; i--) {
            if(!nullComp[i]) mF[i] = mF[--j];
            else             mF[i] = 0.0;
        }

        mFR[0] = mF[0] - mF[7]/2.0 - mF[8]/2.0 - mF[10] - 2.0*mF[11] - 2.0*mF[13]; /* CaMgSi2O6       */
        mFR[1] = mF[1] + mF[10] + mF[11] + mF[13];  			       /* Mg2Si2O6        */
        mFR[2] = mF[2] + mF[7]/2.0 + mF[8]/2.0 + mF[11] + 2.0*mF[13];	       /* CaFeSi2O6       */
        mFR[3] = mF[3] + mF[7] + mF[9] + mF[10] + mF[11] - mF[12];  	       /* Ca(Mg,Ti)AlSiO6 */
        mFR[4] = mF[4] + mF[8] - mF[9] - mF[10] - mF[11] + mF[12];  	       /* Ca(Mg,Ti)FeSiO6 */
        mFR[5] = mF[5] + mF[9] + mF[10] + mF[11];				       /* CaFeAlSiO6      */
        mFR[6] = mF[6] + mF[12];						       /* NaAlSi2O6       */

        /* convert mole fractions of endmembers into independent compos var */
        (*solids[solidID].convert)(SECOND, THIRD, t, p, NULL, mFR, bVec, NULL, NULL, NULL, NULL, NULL);

        if (fabs(*affinity-bVec[nr]) < 0.1/SCALE) foundSolution = TRUE;
        else {
            (*solids[solidID].activity)(FIRST | SECOND, t, p, bVec, activity, muR, NULL);
            for (i=0; i<na; i++) if (mFR[i] != 0.0) gVec[i] = activity[i]/mF[i];

            if (mF[7] != 0.0) {
                activity[7] = exp((muR[3] - muR[0]/2.0 + muR[2]/2.0 + mu0[3] - mu0[0]/2.0 + mu0[2]/2.0 - mu0[7])/(R*t));
                gVec[7] = activity[7]/mF[7];
            } else gVec[7] = 1.0;
            if (mF[8] != 0.0) {
                activity[8] = exp((muR[4] - muR[0]/2.0 + muR[2]/2.0 + mu0[4] - mu0[0]/2.0 + mu0[2]/2.0 - mu0[8])/(R*t));
                gVec[8] = activity[8]/mF[8];
            } else gVec[8] = 1.0;
            if (mF[9] != 0.0) {
                activity[9] = exp((muR[5] - muR[4] + muR[3] + mu0[5] - mu0[4] + mu0[3] - mu0[9])/(R*t));
                gVec[9] = activity[9]/mF[9];
            } else gVec[9] = 1.0;
            if (mF[10] != 0.0) {
                activity[10] = exp((muR[5] - muR[4] + muR[3] - muR[0] + muR[1] + mu0[5] - mu0[4] + mu0[3] - mu0[0] + mu0[1] - mu0[10])/(R*t));
                gVec[10] = activity[10]/mF[10];
            } else gVec[10] = 1.0;
            if (mF[11] != 0.0) {
                activity[11] = exp((muR[5] - muR[4] + muR[3] - 2.0*muR[0] + muR[1] + muR[2] + mu0[5] - mu0[4] + mu0[3] - 2.0*mu0[0] + mu0[1] + mu0[2] - mu0[11])/(R*t));
                gVec[11] = activity[11]/mF[11];
            } else gVec[11] = 1.0;
            if (mF[12] != 0.0) {
                activity[12] = exp((muR[6] - muR[3] + muR[4] + mu0[6] - mu0[3] + mu0[4] - mu0[12])/(R*t));
                gVec[12] = activity[12]/mF[12];
            } else gVec[12] = 1.0;
            if (mF[13] != 0.0) {
                activity[13] = exp((muR[1] - 2.0*muR[0] + 2.0*muR[2] + mu0[1] - 2.0*mu0[0] + 2.0*mu0[2] - mu0[13])/(R*t));
                gVec[13] = activity[13]/mF[13];
            } else gVec[13] = 1.0;
        }
        *affinity = bVec[nr];
#ifdef DEBUG
        {
            char *string;
            (*solids[solidID].display)(FIRST, t, p, bVec, &string);
            printf("SS iter = %2.2d X,g:", iter);
            for (i=0; i<ns; i++) {
                printf(" %13.6g (%13.6g)", mF[i], gVec[i]);
                if ((i == 3) || (i == 7) || (i == 11)) printf("\n                 ");
            }
            printf("\n A: %13.6g Formula: %s\n", bVec[nr], string);
        }
#endif
        iter++;
    }
    for (i=0; i<nr; i++) indepVar[i] = bVec[i];
    *affinity = bVec[nr]*SCALE;

    return (iter < 100) ? SUCCESS : FAILURE;
}

#undef W14
#undef H24
#undef H23

#define W14     20.8 * 1000.0 * 4.184 /* joules */
#define H24      6.55* 1000.0 * 4.184 /* joules */
#define H23      0.0 * 1000.0 * 4.184 /* joules */
#define W13P    11.3 * 1000.0 * 4.184 /* joules */
#define W23PU    9.7 * 1000.0 * 4.184 /* joules */
#define W3P4    10.0 * 1000.0 * 4.184 /* joules */
#define W3PU4U  10.4 * 1000.0 * 4.184 /* joules */
#define W24U    12.6 * 1000.0 * 4.184 /* joules */
#define H25      8.05* 1000.0 * 4.184 /* joules */

int getAffinityAndCompositionSpinel( /* Returns a MODE flag for success or failure */
    double t,              /* temperature (K)					*/
    double p,              /* pressure (bars)					*/
    int    index,          /* index of solid phase in the solids[] structure	*/
    int    *zeroX,         /* TRUE if endmember component has zero mole fraction  */
    double *muMinusMu0,    /* vector of end-member mu - mu0 i.e. A + RTln(a)	*/
    double *affinity,      /* returned value, chemical affinity (J)		*/
    double *indepVar)      /* returned vector, composition of phase (length nr)	*/
{
    int i, j, iter = 0, foundSolution = FALSE;
    static double *mF = NULL, *mFR = NULL, *bVec = NULL, *gVec = NULL, *activity = NULL, *mu = NULL, *mu0 = NULL, *muR = NULL;
    static int    *nullComp = NULL, *nullList = NULL;
    int    hasNull, na, nr, nz, solidID = -1, ns;

    /* Test input parameters */
    if (t <= 0.0) 		return FAILURE;
    if (p <  0.0) 		return FAILURE;
    if (index < 0 || index > npc) return FAILURE;

    /* Check parameters for algorithmic assumptions */
    na = solids[index].na;
    nr = solids[index].nr;
    ns = na + 3;
    solidID = index;

    if (activity == NULL) {
        activity  = (double *)  malloc((unsigned) ns*sizeof (double));
        bVec      = (double *)  malloc((unsigned) ns*sizeof (double));
        gVec      = (double *)  malloc((unsigned) ns*sizeof (double));
        mF        = (double *)  malloc((unsigned) ns*sizeof (double));
        mFR       = (double *)  malloc((unsigned) na*sizeof (double));
        mu        = (double *)  malloc((unsigned) ns*sizeof (double));
        mu0       = (double *)  malloc((unsigned) ns*sizeof (double));
        muR       = (double *)  malloc((unsigned) na*sizeof (double));
        nullComp  = (int *)     malloc((unsigned) ns*sizeof (int));
        nullList  = (int *)     malloc((unsigned) ns*sizeof (int));
    }

    for (i=0, hasNull=FALSE, nz=0; i<na; i++) {
        nullComp[i] = zeroX[i]; hasNull |= zeroX[i];
        if (!nullComp[i]) { nullList[nz] = i; mu[nz++] = muMinusMu0[i]; }
        mu0[i] = (solids[solidID+1+i].cur).g;
    }
    if (nz == 0) return FAILURE;

    /* compute dependent species mu0 for MgCr2O4 */
    if (!nullComp[0] && !nullComp[1] && !nullComp[3]) {
        double correction = 2.0*(H24) - ((W14)+(W23PU)+(W3P4)) + ((W13P)+(W24U)+(W3PU4U));
        mu0[5] = mu0[0] + mu0[3] - mu0[1] + correction;
        nullComp[5] = FALSE;
        nullList[nz] = 5;
        mu[nz++] = muMinusMu0[0] + muMinusMu0[3] - muMinusMu0[1] + correction;
    } else { mu0[5] = 0.0; nullComp[5] = TRUE; }

    /* compute dependent species mu0 for MgFe2O3 */
    if (!nullComp[1] && !nullComp[2] && !nullComp[3]) {
        double correction = (H25);
        mu0[6] = mu0[2] + mu0[3] - mu0[1] + correction;
        nullComp[6] = FALSE;
        nullList[nz] = 6;
        mu[nz++] = muMinusMu0[2] + muMinusMu0[3] - muMinusMu0[1] + correction;
    } else { mu0[6] = 0.0; nullComp[6] = TRUE; }

    /* compute dependent species mu0 for Mg2TiO4 */
    if (!nullComp[1] && !nullComp[3] && !nullComp[4]) {
        double correction = 2.0*(H24);
        mu0[7] = mu0[4] + 2.0*(mu0[3] - mu0[1]) + correction;
        nullComp[7] = FALSE;
        nullList[nz] = 7;
        mu[nz++] = muMinusMu0[4] + 2.0*(muMinusMu0[3] - muMinusMu0[1]) + correction;
    } else { mu0[7] = 0.0; nullComp[7] = TRUE; }

#ifdef DEBUG_EXTRA
    for (i=0, j=0; i<ns; i++) if (!nullComp[i]) printf("...mu0[%2.2d] = %13.6g mu = %13.6g\n", i, mu0[i], mu[j++]);
#endif

    for (i=0; i<ns; i++) gVec[i] = 1.0;
    *affinity = 100000.0;

    while (!foundSolution && (iter < 100)) {
        double sum;

        /* take appropriate action for the number of non-zero components */
         if (nz == 0)                      bVec[na-1] = 0.0;
        else if (nz == 1) { mF[0] = 1.0; bVec[na-1] = (mu[0]-gVec[0])/SCALE; }
        else {
            /* condense the activity cofficient vector */
            for (i=0, j=0; i<ns; i++) if (!nullComp[i]) gVec[j++] = R*t*log(gVec[i]);

            /* Compute the f[n] terms and store them temporarily in mF[] */
            sum = 1.0;
            if (nz > 2) for (i=0; i<(nz-2); i++) {
                mF[i] = exp(-(mu[i]+gVec[i]-mu[nz-2]-gVec[nz-2])/(R*t)); sum += mF[i];
            }
            mF[nz-2] = exp(-(mu[nz-2]+gVec[nz-2]-mu[nz-1]-gVec[nz-1])/(R*t));

            /* Solve for the composition variables (mole fractions) */
            mF[nz-2] /= 1.0 + mF[nz-2]*sum;
            mF[nz-1]  = 1.0 - mF[nz-2];
            if (nz > 2) for (i=0; i<(nz-2); i++) {
                mF[i] *= mF[nz-2]; mF[nz-1] -= mF[i];
            }
            /* This fix is taken from Thermoengine c246b6c56ebddfd170fc2b36262fa4abb9327bad */
            for (i=0; i<nz; i++) if (mF[i] < DBL_EPSILON) mF[i] = DBL_EPSILON;

            /* compute the chemical affinity (choice of mu[] is arbitrary) */
            bVec[na-1] = (mu[0] + R*t*log(mF[0]) + gVec[0])/SCALE;
        }

        /* Reassemble the mole fraction and chemical potential vectors with zeros for the absent endmembers */
        if (hasNull) for (i=ns-1, j=nz; i>=0; i--) {
            if(!nullComp[i]) mF[i] = mF[--j];
            else             mF[i] = 0.0;
        }

        mFR[0] = mF[0] + mF[5];                                       /* FeCr2O4 */
        mFR[1] = mF[1] + mF[3] - (mF[3] + mF[5] + mF[6] + 2.0*mF[7]); /* FeAl2O4 */
        mFR[2] = mF[2] + mF[6];	                                      /* Fe3O4   */
        mFR[3] = mF[3] + mF[5] + mF[6] + 2.0*mF[7];    	              /* MgAl2O4 */
        mFR[4] = mF[4] + mF[7];  	                                  /* Fe2TiO4 */

        /* convert mole fractions of endmembers into independent compos var */
        (*solids[solidID].convert)(SECOND, THIRD, t, p, NULL, mFR, bVec, NULL, NULL, NULL, NULL, NULL);

        if (fabs(*affinity-bVec[nr]) < 0.1/SCALE) foundSolution = TRUE;
        else {
            (*solids[solidID].activity)(FIRST | SECOND, t, p, bVec, activity, muR, NULL);
            for (i=0; i<na; i++) if (mFR[i] != 0.0) gVec[i] = activity[i]/mF[i];

            if (mF[5] != 0.0) {
                activity[5] = exp((muR[0] + muR[3] - muR[1] + mu0[0] + mu0[3] - mu0[1] - mu0[5])/(R*t));
                gVec[5] = activity[5]/mF[5];
            } else gVec[5] = 1.0;
            if (mF[6] != 0.0) {
                activity[6] = exp((muR[2] + muR[3] - muR[1] + mu0[2] + mu0[3] - mu0[1] - mu0[6])/(R*t));
                gVec[6] = activity[6]/mF[6];
            } else gVec[6] = 1.0;
            if (mF[7] != 0.0) {
                activity[7] = exp((muR[4] + 2.0*(muR[3] - muR[1]) + mu0[4] + 2.0*(mu0[3] - mu0[1]) - mu0[7])/(R*t));
                gVec[7] = activity[7]/mF[7];
            } else gVec[7] = 1.0;
        }
        *affinity = bVec[nr];
#ifdef DEBUG
        {
            char *string;
            (*solids[solidID].display)(FIRST, t, p, bVec, &string);
            printf("SS iter = %2.2d X,g:", iter);
            for (i=0; i<ns; i++) {
                printf(" %13.6g (%13.6g)", mF[i], gVec[i]);
                if ((i == 3) || (i == 7) || (i == 11)) printf("\n                 ");
            }
            printf("\n A: %13.6g Formula: %s\n", bVec[nr], string);
        }
#endif
        iter++;
    }
    for (i=0; i<nr; i++) indepVar[i] = bVec[i];
    *affinity = bVec[nr]*SCALE;

    return (iter < 100) ? SUCCESS : FAILURE;
}


/* end of file EST_SATSTATE_REVISED.C */
