const char *orthopyroxene_ver(void) { return "$Id: orthopyroxene.c,v 1.4 2007/05/07 18:23:21 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: orthopyroxene.c,v $
MELTS Source Code: RCS Revision 1.4  2007/05/07 18:23:21  ghiorso
MELTS Source Code: RCS Modifications to LEPR and calibration algorithms following visit by
MELTS Source Code: RCS Marc and Tim.  Mostly work on cpx and opx inclusion and reclassification.
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
MELTS Source Code: RCS Revision 1.1.1.1  2004/01/02 19:21:49  cvsaccount
MELTS Source Code: RCS CTserver University of Chicago
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2001/12/21 03:04:26  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:58:12  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
 * Revision 3.7  1997/06/21  22:49:32  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.6  1997/05/03  20:23:11  ghiorso
 * *** empty log message ***
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Routines to compute orthopyroxene solution properties
**      (file: ORTHOPYROXENE.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0  Mark S. Ghiorso  May 3, 1997 Original Version
**            Canalbalized PYROXENE.C (now CLINOPYROXENE.C)
**--
*/

#ifdef ISCLINO /* If defined,  structural state is calculated          */
#undef ISCLINO /* If undefine, structural state is always orthorhombic */
#endif

#ifdef REGRESS_DI_EN_ONLY
#undef REGRESS_DI_EN_ONLY
#endif

#ifdef REGRESS_LOW_P_ONLY
#undef REGRESS_LOW_P_ONLY
#endif

#ifdef DEBUG
#undef DEBUG
#endif

#include "melts_gsl.h"
#include "silmin.h"  /* Structure definitions foor SILMIN package */

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

#define MAX_ITER 100  /* Max number of iterations allowed in ordering alg */

/*
 *=============================================================================
 * Pyroxene solution parameters:
 * Sack, R.O., Ghiorso, M.S. (1993)
 *   Work in progress.
 */

#define  R       8.3143
#define  NR      6       /* Six independent composition variables */
#define  NS      2       /* Two ordering parameters               */
#define  NA      7       /* Seven endmember compositions          */

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

/*
 * Definitions of Taylor expansion coefficients in terms of solution
 * parameters. Independent variables are x2, x3, x4, x5, x6, x7, s1, s2
 */

#define H0    (H1)
#define HX2   (H2)-(H1) + (WH12)
#define HX3   (H3)-(H1) + (W13)
#define HX4   (H4)-(H1) + (W14)
#define HX5   (H5)-(H1) + 0.5*((W15)+(W15P)+(H55))
#define HX6   (H6) + 0.5*((H4)-(H1)-(H3)-(H5)) - 0.5*(W13) + 0.5*(W14) \
                            + 0.25*(W15) - 0.75*(W15P) + (W16) - 0.75*(H55)
#define HX7   (H7)-(H1) + (pH1) + 0.5*((WHCAFE)+(WHCAMG)-(WH12)) \
                            + 0.25*((HEX)+(HX)+(H027)) + 0.5*(DWHCAMG) - 1.5*(DWHCAFE)

#define HS1   0.5*((W15)-(W15P)-(H55))
#define HS2   0.5*((WHCAFE)-(WHCAMG)-(WH12)) + 0.25*((HEX)+(HX)+(H027)) \
                            - 0.5*(DWHCAMG) - 1.5*(DWHCAFE)

#define HX2X2 -(WH12)
#define HX2X3 2.0*(H23) - 0.5*(WH12)
#define HX2X4 2.0*(H24) - 0.5*(WH12)
#define HX2X5 0.5*((W25)+(W25P)-(W15)-(W15P)) - (WH12)
#define HX2X6 (W26)-(W16) + 0.25*((W25)-(W15)) - 0.75*((W25P)-(W15P)) \
                            + (H24) - (H23) - 0.5*(WH12)
#define HX2X7 (pH2)-(pH1) + (WH12) + 0.5*((H027)-(HEX)) - 2.0*(DWHCAMG) \
                            + 2.0*(DWHCAFE)
#define HX2S1 0.5*((W25)-(W25P)-(W15)+(W15P))
#define HX2S2 (WH12) - 0.5*(HX) + 2.0*(DWHCAMG) + 2.0*(DWHCAFE)

#define HX3X3 -(W13)
#define HX3X4 (W34) - (W13) - (W14)
#define HX3X5 0.5*((W35)+(W35P)-(W15)-(W15P)) - (W13)
#define HX3X6 (W36)-(W16) + 0.25*((W35)-(W15)) + 0.5*((W34)-(W14)) \
                            - 0.75*((W35P)-(W15P))
#define HX3X7 0.125*((H027)-(HEX)-(HX))+0.5*((W37)+(W3U7U)-(WHCAFE)-(WHCAMG)) \
                            - 1.5*(H23) - (W13) + 0.25*(WH12) - (DWHCAMG) + (DWHCAFE) \
                            + (pH3)-(pH1)
#define HX3S1 0.5*((W35)-(W35P)-(W15)+(W15P))
#define HX3S2 0.125*((H027)-(HEX)-(HX))+0.5*((W3U7U)-(W37)+(WHCAMG)-(WHCAFE)) \
                            - 1.5*(H23) + 0.25*(WH12) + (DWHCAMG) + (DWHCAFE)

#define HX4X4 -(W14)
#define HX4X5 0.5*((W45)+(W45P)-(W15)-(W15P)) - (W14)
#define HX4X6 (W46)-(W16)-(W14) + 0.25*((W45)-(W15)) - 0.5*((W34)-(W13)) \
                            - 0.75*((W45P)-(W15P))
#define HX4X7 0.125*((H027)-(HEX)-(HX))+0.5*((W47)+(W4U7U)-(WHCAFE)-(WHCAMG)) \
                            - 1.5*(H24) - (W14) + 0.25*(WH12) - (DWHCAMG) + (DWHCAFE) \
                            + (pH4)-(pH1)
#define HX4S1 0.5*((W45)-(W45P)-(W15)+(W15P))
#define HX4S2 0.125*((H027)-(HEX)-(HX))+0.5*((W4U7U)-(W47)+(WHCAMG)-(WHCAFE)) \
                            - 1.5*(H24) + 0.25*(WH12) + (DWHCAMG) + (DWHCAFE)

#define HX5X5 0.25*(W55) - 0.5*((W15)+(W15P))
#define HX5X6 0.5*((W56)+(W5P6)-(W15)+(W15P)-2.0*(W16)) - 0.25*(W55) \
                            + 0.25*((W45)+(W45P)-2.0*(W14)) - 0.25*((W35)+(W35P)-2.0*(W13))
#define HX5X7 0.25*((H027)-(HEX)-(HX)) + 0.25*((W57U)+(W5P7U)+(W57)+(W5P7)) \
                            - 0.5*((WHCAFE)+(WHCAMG)+(W25)+(W25P)-(WH12)) - (DWHCAMG) \
                            + 2.0*(DWHCAFE) + (pH5)-(pH1)
#define HX5S1 - 0.5*((W15)-(W15P))
#define HX5S2 0.25*((H027)-(HEX)-(HX)) + 0.25*((W57U)+(W5P7U)-(W57)-(W5P7)) \
                            + 0.5*((W15)+(W15P)-(W25)-(W25P)+(WHCAMG)-(WHCAFE)+(WH12)) \
                            + (DWHCAMG) + 2.0*(DWHCAFE)

#define HX6X6 0.25*(W56) - 0.75*(W5P6) + 0.5*((W46)-(W36)-(W16)) \
                            - 0.1875*(W55) + 0.125*((W45)-(W35)-(W15)) \
                            - 0.375*((W45P)-(W35P)-(W15P)) - 0.25*((W34)+(W14)-(W13))
#define HX6X7 0.5*((W67)+(W67U)) - 0.375*((W5P7)+(W5P7U)) \
                            + 0.125*((W57)+(W57U)) + 0.25*((W47)+(W4U7U)) \
                            - 0.25*((W37)+(W3U7U)) - (W26) - 0.25*(W25) + 0.75*(W25P) \
                            - 0.5*((W14)-(W13)) + 0.25*(WH12) - 0.25*((WHCAFE)+(WHCAMG)) \
                            - 0.125*((HEX)+(HX)-(H027)) - 0.75*((H24)-(H23)) \
                            - 0.50*(DWHCAMG) + (DWHCAFE) + (pHj)-0.5*((pH1)+(pH3)-(pH4)+(pH5))
#define HX6S1 0.5*((W56)-(W5P6)) + 0.25*((W45)-(W45P)) - 0.25*((W35)-(W35P)) \
                            - 0.25*((W15)-(W15P)) - 0.5*(W55)
#define HX6S2 0.5*((W67U)-(W67)) + 0.125*((W57U)-(W57)) \
                            + 0.375*((W5P7)-(W5P7U)) + 0.25*((W4U7U)-(W47)) \
                            + 0.25*((W37)-(W3U7U)) + (W16)-(W26) + 0.75*((W25P)-(W15P)) \
                            + 0.25*((W15)-(W25)) + 0.25*(WH12) - 0.25*((WHCAFE)-(WHCAMG)) \
                            - 0.125*((HEX)+(HX)-(H027)) - 0.75*((H24)-(H23)) \
                            + 0.50*(DWHCAMG) + (DWHCAFE)

#define HX7X7 0.25*((HEX)-(H027)) - 0.5*((WHCAFE)+(WHCAMG)) \
                            + 0.25*((WHFEMG)-(WH12)) + 0.75*(DWHCAMG) + 0.25*(DWHCAFE) \
                            + (pH7)-(pH1)+0.25*((dH027)+(dHEX)+(dHX))
#define HX7S1 0.25*((W57U)+(W57)-(W5P7U)-(W5P7)) - 0.5*((W25)-(W25P))
#define HX7S2 0.5*((WHCAMG)-(WHCAFE)-(WH12)) + 0.25*((HEX)+(HX)-(H027)) \
                            + 1.5*(DWHCAMG) - 1.5*(DWHCAFE) + 0.25*((dH027)+(dHEX)+(dHX))

#define HS1S1 -0.25*(W55)
#define HS1S2 0.25*((W57U)+(W5P7)-(W5P7U)-(W57)) \
                            + 0.5*((W25P)-(W25)+(W15)-(W15P))

#define HS2S2 0.25*((HX)-(WHFEMG)-(WH12)) - 2.25*(DWHCAMG) - 1.75*(DWHCAFE)

#define HX2X2X7 - 0.5*(DWHCAMG) + 0.5*(DWHCAFE)
#define HX2X2S2 0.5*(DWHCAMG) - 0.5*(DWHCAFE)
#define HX2X3X7 (DWHCAMG) - (DWHCAFE)
#define HX2X3S2 - (DWHCAMG) - (DWHCAFE)
#define HX2X4X7 (DWHCAMG) - (DWHCAFE)
#define HX2X4S2 - (DWHCAMG) - (DWHCAFE)
#define HX2X5X7 (DWHCAMG) - (DWHCAFE)
#define HX2X5S2 - (DWHCAMG) - (DWHCAFE)
#define HX2X6X7 0.5*(DWHCAMG) - 0.5*(DWHCAFE)
#define HX2X6S2 - 0.5*(DWHCAMG) - 0.5*(DWHCAFE)
#define HX2X7X7 2.75*(DWHCAMG) - 2.75*(DWHCAFE) + 0.5*((dH027)-(dHEX))
#define HX2X7S2 - 2.5*(DWHCAMG) - 1.5*(DWHCAFE) - 0.5*(dHX)
#define HX2S2S2 - 0.25*(DWHCAMG) + 0.25*(DWHCAFE)
#define HX3X3X7 0.5*(DWHCAMG)
#define HX3X3S2 - 0.5*(DWHCAMG)
#define HX3X4X7 (DWHCAMG) + 0.25*(DWHCAFE)
#define HX3X4S2 - (DWHCAMG) + 0.25*(DWHCAFE)
#define HX3X5X7 (DWHCAMG) - 0.25*(DWHCAFE)
#define HX3X5S2 - (DWHCAMG) - 0.25*(DWHCAFE)
#define HX3X6X7 0.5*(DWHCAMG)
#define HX3X6S2 - 0.5*(DWHCAMG)
#define HX3X7X7 0.25*(DWHCAMG) + 0.625*(DWHCAFE) + 0.125*((dH027)-(dHEX)-(dHX))
#define HX3X7S2 - 1.5*(DWHCAMG) + 1.5*(DWHCAFE) + 0.125*((dH027)-(dHEX)-(dHX))
#define HX3S2S2 1.25*(DWHCAMG) + 0.875*(DWHCAFE)
#define HX4X4X7 0.5*(DWHCAMG)
#define HX4X4S2 - 0.5*(DWHCAMG)
#define HX4X5X7 (DWHCAMG) - 0.25*(DWHCAFE)
#define HX4X5S2 - (DWHCAMG) - 0.25*(DWHCAFE)
#define HX4X6X7 0.5*(DWHCAMG) - 0.25*(DWHCAFE)
#define HX4X6S2 - 0.5*(DWHCAMG) - 0.25*(DWHCAFE)
#define HX4X7X7 0.25*(DWHCAMG) + 0.625*(DWHCAFE) + 0.125*((dH027)-(dHEX)-(dHX))
#define HX4X7S2 - 1.5*(DWHCAMG) + 1.5*(DWHCAFE) + 0.125*((dH027)-(dHEX)-(dHX))
#define HX4S2S2 1.25*(DWHCAMG) + 0.875*(DWHCAFE)
#define HX5X5X7 0.5*(DWHCAMG) - 0.5*(DWHCAFE)
#define HX5X5S2 - 0.5*(DWHCAMG) - 0.5*(DWHCAFE)
#define HX5X6X7 0.5*(DWHCAMG) - 0.5*(DWHCAFE)
#define HX5X6S2 - 0.5*(DWHCAMG) - 0.5*(DWHCAFE)
#define HX5X7X7 0.25*(DWHCAMG) - 0.25*(DWHCAFE) + 0.25*((dH027)-(dHEX)-(dHX))
#define HX5X7S2 - 1.5*(DWHCAMG) + 0.5*(DWHCAFE) + 0.25*((dH027)-(dHEX)-(dHX))
#define HX5S2S2 1.25*(DWHCAMG) + 0.75*(DWHCAFE)
#define HX6X6X7 0.125*(DWHCAMG) - 0.1875*(DWHCAFE)
#define HX6X6S2 - 0.125*(DWHCAMG) - 0.1875*(DWHCAFE)
#define HX6X7X7 0.125*(DWHCAMG) - 0.125*(DWHCAFE) + 0.125*((dH027)-(dHEX)-(dHX))
#define HX6X7S2 - 0.75*(DWHCAMG) + 0.25*(DWHCAFE) + 0.125*((dH027)-(dHEX)-(dHX))
#define HX6S2S2 0.625*(DWHCAMG) + 0.375*(DWHCAFE)
#define HX7X7X7 - 1.5*(DWHCAMG) + 1.5*(DWHCAFE) + 0.25*((dHEX)-(dH027))
#define HX7X7S2 - (DWHCAMG) + 3.0*(DWHCAFE) + 0.25*((dHEX)+(dHX)-(dH027))
#define HX7S2S2 2.5*(DWHCAMG) + 1.5*(DWHCAFE) + 0.25*(dHX)

#define S0    (S1)
#define SX2   (S2) - (S1)
#define SX3   (S3) - (S1)
#define SX4   (S4) - (S1)
#define SX5   (S5) - (S1) + 0.5*(S55)
#define SX6   (S6) + 0.5*((S4)-(S1)-(S3)-(S5)) - 0.75*(S55)
#define SX7   (S7) - (S1) + (pS1) + 0.25*(S027)

#define SS1   -0.5*(S55)
#ifdef TRUE_xMELTS
#define SS2   ((calculationMode == MODE_xMELTS) ? 0.25*(S027) : 0.0)
#else
#define SS2   0.0
#endif

#define SX2X2 0.0
#define SX2X3 2.0*(S23)
#define SX2X4 2.0*(S24)
#define SX2X5 0.0
#define SX2X6 (S24) - (S23)
#define SX2X7 (pS2)-(pS1) + 0.5*(S027)
#define SX2S1 0.0
#define SX2S2 0.0

#define SX3X3 0.0
#define SX3X4 0.0
#define SX3X5 0.0
#define SX3X6 0.0
#define SX3X7 0.125*(S027) - 1.5*(S23) + (pS3)-(pS1)
#define SX3S1 0.0
#define SX3S2 0.125*(S027) - 1.5*(S23)

#define SX4X4 0.0
#define SX4X5 0.0
#define SX4X6 0.0
#define SX4X7 0.125*(S027) - 1.5*(S24) + (pS4)-(pS1)
#define SX4S1 0.0
#define SX4S2 0.125*(S027) - 1.5*(S24)

#define SX5X5 0.0
#define SX5X6 0.0
#define SX5X7 0.25*(S027) + (pS5)-(pS1)
#define SX5S1 0.0
#define SX5S2 0.25*(S027)

#define SX6X6 0.0
#define SX6X7 0.125*(S027) - 0.75*((S24)-(S23)) \
                            + (pSj)-0.5*((pS1)+(pS3)-(pS4)+(pS5))
#define SX6S1 0.0
#define SX6S2 0.125*(S027) - 0.75*((S24)-(S23))

#define SX7X7 -0.25*(S027) + (pS7)-(pS1)+0.25*(dS027)
#define SX7S1 0.0
#define SX7S2 -0.25*(S027) + 0.25*(dS027)

#define SS1S1 0.0
#define SS1S2 0.0

#define SS2S2 0.0

#define SX2X7X7 0.5*(dS027)
#define SX3X7X7 0.125*(dS027)
#define SX3X7S2 0.125*(dS027)
#define SX4X7X7 0.125*(dS027)
#define SX4X7S2 0.125*(dS027)
#define SX5X7X7 0.25*(dS027)
#define SX5X7S2 0.25*(dS027)
#define SX6X7X7 0.125*(dS027)
#define SX6X7S2 0.125*(dS027)
#define SX7X7X7 -0.25*(dS027)
#define SX7X7S2 -0.25*(dS027)

#define V0    (V1)
#define VX2   (V2)-(V1) + (WV12)
#define VX3   (V3)-(V1)
#define VX4   (V4)-(V1)
#define VX5   (V5)-(V1)
#define VX6   (V6) + 0.5*((V4)-(V1)-(V3)-(V5))
#define VX7   (V7)-(V1) + (pV1) + 0.5*((WVCAFE)+(WVCAMG)-(WV12)) \
                            + 0.25*((VEX)+(VX)+(V027)) + 0.5*(DWVCAMG) - 1.5*(DWVCAFE)

#define VS1   0.0
#define VS2   0.5*((WVCAFE)-(WVCAMG)-(WV12)) + 0.25*((VEX)+(VX)+(V027)) \
                            - 0.5*(DWVCAMG) - 1.5*(DWVCAFE)

#define VX2X2 -(WV12)
#define VX2X3 -0.5*(WV12)
#define VX2X4 -0.5*(WV12)
#define VX2X5 -(WV12)
#define VX2X6 -0.5*(WV12)
#define VX2X7 (pV2)-(pV1) + (WV12) + 0.5*((V027)-(VEX)) - 2.0*(DWVCAMG) \
                            + 2.0*(DWVCAFE)
#define VX2S1 0.0
#define VX2S2 (WV12) - 0.5*(VX) + 2.0*(DWVCAMG) + 2.0*(DWVCAFE)

#define VX3X3 0.0
#define VX3X4 0.0
#define VX3X5 0.0
#define VX3X6 0.0
#define VX3X7 0.125*((V027)-(VEX)-(VX)) - 0.5*((WVCAFE)+(WVCAMG)) \
                            + 0.25*(WV12) - (DWVCAMG) + (DWVCAFE) + (pV3)-(pV1)
#define VX3S1 0.0
#define VX3S2 0.125*((V027)-(VEX)-(VX)) + 0.5*((WVCAMG)-(WVCAFE)) \
                            + 0.25*(WV12) + (DWVCAMG) + (DWVCAFE)

#define VX4X4 0.0
#define VX4X5 0.0
#define VX4X6 0.0
#define VX4X7 0.125*((V027)-(VEX)-(VX)) - 0.5*((WVCAFE)+(WVCAMG)) \
                            + 0.25*(WV12) - (DWVCAMG) + (DWVCAFE) + (pV4)-(pV1)
#define VX4S1 0.0
#define VX4S2 0.125*((V027)-(VEX)-(VX)) + 0.5*((WVCAMG)-(WVCAFE)) \
                            + 0.25*(WV12) + (DWVCAMG) + (DWVCAFE)

#define VX5X5 0.0
#define VX5X6 0.0
#define VX5X7 0.25*((V027)-(VEX)-(VX)) - 0.5*((WVCAFE)+(WVCAMG)-(WV12)) \
                            - (DWVCAMG) + 2.0*(DWVCAFE) + (pV5)-(pV1)
#define VX5S1 0.0
#define VX5S2 0.25*((V027)-(VEX)-(VX)) + 0.5*((WVCAMG)-(WVCAFE)+(WV12)) \
                            + (DWVCAMG) + 2.0*(DWVCAFE)

#define VX6X6 0.0
#define VX6X7 0.25*(WV12) - 0.25*((WVCAFE)+(WVCAMG)) \
                            - 0.125*((VEX)+(VX)-(V027)) - 0.5*(DWVCAMG) + (DWVCAFE) \
                            + (pVj)-0.5*((pV1)+(pV3)-(pV4)+(pV5))
#define VX6S1 0.0
#define VX6S2 0.25*(WV12) - 0.25*((WVCAFE)-(WVCAMG)) \
                            - 0.125*((VEX)+(VX)-(V027)) + 0.5*(DWVCAMG) + (DWVCAFE)

#define VX7X7 0.25*((VEX)-(V027)) - 0.5*((WVCAFE)+(WVCAMG)) \
                            + 0.25*((WVFEMG)-(WV12)) + 0.75*(DWVCAMG) + 0.25*(DWVCAFE) \
                            + (pV7)-(pV1)+0.25*((dV027)+(dVEX)+(dVX))
#define VX7S1 0.0
#define VX7S2 0.5*((WVCAMG)-(WVCAFE)-(WV12)) + 0.25*((VEX)+(VX)-(V027)) \
                            + 1.5*(DWVCAMG) - 1.5*(DWVCAFE) + 0.25*((dV027)+(dVEX)+(dVX))

#define VS1S1 0.0
#define VS1S2 0.0

#define VS2S2 0.25*((VX)-(WVFEMG)-(WV12)) - 2.25*(DWVCAMG) - 1.75*(DWVCAFE)

#define VX2X2X7 - 0.5*(DWVCAMG) + 0.5*(DWVCAFE)
#define VX2X2S2 0.5*(DWVCAMG) - 0.5*(DWVCAFE)
#define VX2X3X7 (DWVCAMG) - (DWVCAFE)
#define VX2X3S2 - (DWVCAMG) - (DWVCAFE)
#define VX2X4X7 (DWVCAMG) - (DWVCAFE)
#define VX2X4S2 - (DWVCAMG) - (DWVCAFE)
#define VX2X5X7 (DWVCAMG) - (DWVCAFE)
#define VX2X5S2 - (DWVCAMG) - (DWVCAFE)
#define VX2X6X7 0.5*(DWVCAMG) - 0.5*(DWVCAFE)
#define VX2X6S2 - 0.5*(DWVCAMG) - 0.5*(DWVCAFE)
#define VX2X7X7 2.75*(DWVCAMG) - 2.75*(DWVCAFE) + 0.5*((dV027)-(dVEX))
#define VX2X7S2 - 2.5*(DWVCAMG) - 1.5*(DWVCAFE) - 0.5*(dVX)
#define VX2S2S2 - 0.25*(DWVCAMG) + 0.25*(DWVCAFE)
#define VX3X3X7 0.5*(DWVCAMG)
#define VX3X3S2 - 0.5*(DWVCAMG)
#define VX3X4X7 (DWVCAMG) + 0.25*(DWVCAFE)
#define VX3X4S2 - (DWVCAMG) + 0.25*(DWVCAFE)
#define VX3X5X7 (DWVCAMG) - 0.25*(DWVCAFE)
#define VX3X5S2 - (DWVCAMG) - 0.25*(DWVCAFE)
#define VX3X6X7 0.5*(DWVCAMG)
#define VX3X6S2 - 0.5*(DWVCAMG)
#define VX3X7X7 0.25*(DWVCAMG) + 0.625*(DWVCAFE) + 0.125*((dV027)-(dVEX)-(dVX))
#define VX3X7S2 - 1.5*(DWVCAMG) + 1.5*(DWVCAFE) + 0.125*((dV027)-(dVEX)-(dVX))
#define VX3S2S2 1.25*(DWVCAMG) + 0.875*(DWVCAFE)
#define VX4X4X7 0.5*(DWVCAMG)
#define VX4X4S2 - 0.5*(DWVCAMG)
#define VX4X5X7 (DWVCAMG) - 0.25*(DWVCAFE)
#define VX4X5S2 - (DWVCAMG) - 0.25*(DWVCAFE)
#define VX4X6X7 0.5*(DWVCAMG) - 0.25*(DWVCAFE)
#define VX4X6S2 - 0.5*(DWVCAMG) - 0.25*(DWVCAFE)
#define VX4X7X7 0.25*(DWVCAMG) + 0.625*(DWVCAFE) + 0.125*((dV027)-(dVEX)-(dVX))
#define VX4X7S2 - 1.5*(DWVCAMG) + 1.5*(DWVCAFE) + 0.125*((dV027)-(dVEX)-(dVX))
#define VX4S2S2 1.25*(DWVCAMG) + 0.875*(DWVCAFE)
#define VX5X5X7 0.5*(DWVCAMG) - 0.5*(DWVCAFE)
#define VX5X5S2 - 0.5*(DWVCAMG) - 0.5*(DWVCAFE)
#define VX5X6X7 0.5*(DWVCAMG) - 0.5*(DWVCAFE)
#define VX5X6S2 - 0.5*(DWVCAMG) - 0.5*(DWVCAFE)
#define VX5X7X7 0.25*(DWVCAMG) - 0.25*(DWVCAFE) + 0.25*((dV027)-(dVEX)-(dVX))
#define VX5X7S2 - 1.5*(DWVCAMG) + 0.5*(DWVCAFE) + 0.25*((dV027)-(dVEX)-(dVX))
#define VX5S2S2 1.25*(DWVCAMG) + 0.75*(DWVCAFE)
#define VX6X6X7 0.125*(DWVCAMG) - 0.1875*(DWVCAFE)
#define VX6X6S2 - 0.125*(DWVCAMG) - 0.1875*(DWVCAFE)
#define VX6X7X7 0.125*(DWVCAMG) - 0.125*(DWVCAFE) + 0.125*((dV027)-(dVEX)-(dVX))
#define VX6X7S2 - 0.75*(DWVCAMG) + 0.25*(DWVCAFE) + 0.125*((dV027)-(dVEX)-(dVX))
#define VX6S2S2 0.625*(DWVCAMG) + 0.375*(DWVCAFE)
#define VX7X7X7 - 1.5*(DWVCAMG) + 1.5*(DWVCAFE) + 0.25*((dVEX)-(dV027))
#define VX7X7S2 - (DWVCAMG) + 3.0*(DWVCAFE) + 0.25*((dVEX)+(dVX)-(dV027))
#define VX7S2S2 2.5*(DWVCAMG) + 1.5*(DWVCAFE) + 0.25*(dVX)

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
    MTHREAD_KEY_CREATE(&sOldPureKey,   free);
    MTHREAD_KEY_CREATE(&d2gds2PureKey, free);
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

static double getSOldPure() {
    double *sOldPurePt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    sOldPurePt = (double *) MTHREAD_GETSPECIFIC(sOldPureKey);
    if (sOldPurePt == NULL) {
        sOldPurePt  = (double *) malloc(sizeof(double));
        *sOldPurePt = -9999.0;
        MTHREAD_SETSPECIFIC(sOldPureKey, (void *) sOldPurePt);
    }
    return *sOldPurePt;
}

static void setSOldPure(double sOldPure) {
    double *sOldPurePt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    sOldPurePt = (double *) MTHREAD_GETSPECIFIC(sOldPureKey);
    if (sOldPurePt == NULL) {
        sOldPurePt  = (double *) malloc(sizeof(double));
        *sOldPurePt = -9999.0;
        MTHREAD_SETSPECIFIC(sOldPureKey, (void *) sOldPurePt);
    }
    *sOldPurePt = sOldPure;
}

static double getD2gds2Pure() {
    double *d2gds2PurePt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    d2gds2PurePt = (double *) MTHREAD_GETSPECIFIC(d2gds2PureKey);
    if (d2gds2PurePt == NULL) {
        d2gds2PurePt  = (double *) malloc(sizeof(double));
        *d2gds2PurePt = -9999.0;
        MTHREAD_SETSPECIFIC(d2gds2PureKey, (void *) d2gds2PurePt);
    }
    return *d2gds2PurePt;
}

static void setD2gds2Pure(double d2gds2Pure) {
    double *d2gds2PurePt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    d2gds2PurePt = (double *) MTHREAD_GETSPECIFIC(d2gds2PureKey);
    if (d2gds2PurePt == NULL) {
        d2gds2PurePt  = (double *) malloc(sizeof(double));
        *d2gds2PurePt = -9999.0;
        MTHREAD_SETSPECIFIC(d2gds2PureKey, (void *) d2gds2PurePt);
    }
    *d2gds2PurePt = d2gds2Pure;
}

/***********************************/
/* Statics for Site Mole Fractions */
/***********************************/

#define DECLARE_SITE_FRACTIONS \
    double xal3m1, xfe2m1, xfe3m1, xmg2m1, xti4m1, xca2m2, xfe2m2, xmg2m2, xna1m2, xal3tet, xfe3tet, xsi4tet;

#define XAL3M1   0
#define XFE2M1   1
#define XFE3M1   2
#define XMG2M1   3
#define XTI4M1   4
#define XCA2M2   5
#define XFE2M2   6
#define XMG2M2   7
#define XNA1M2   8
#define XAL3TET  9
#define XFE3TET 10
#define XSI4TET 11

#define NX      12

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
    xal3m1  = getX(XAL3M1); \
    xfe2m1  = getX(XFE2M1); \
    xfe3m1  = getX(XFE3M1); \
    xmg2m1  = getX(XMG2M1); \
    xti4m1  = getX(XTI4M1); \
    xca2m2  = getX(XCA2M2); \
    xfe2m2  = getX(XFE2M2); \
    xmg2m2  = getX(XMG2M2); \
    xna1m2  = getX(XNA1M2); \
    xal3tet = getX(XAL3TET); \
    xfe3tet = getX(XFE3TET); \
    xsi4tet = getX(XSI4TET);

#define SET_SITE_FRACTIONS \
    setX(XAL3M1,  xal3m1); \
    setX(XFE2M1,  xfe2m1); \
    setX(XFE3M1,  xfe3m1); \
    setX(XMG2M1,  xmg2m1); \
    setX(XTI4M1,  xti4m1); \
    setX(XCA2M2,  xca2m2); \
    setX(XFE2M2,  xfe2m2); \
    setX(XMG2M2,  xmg2m2); \
    setX(XNA1M2,  xna1m2); \
    setX(XAL3TET, xal3tet); \
    setX(XFE3TET, xfe3tet); \
    setX(XSI4TET, xsi4tet);

/*
 *=============================================================================
 * Local function to compute ordering state and associated derivatives of
 * pure component endmembers
 */

#define DI_S          (S0)
#define DI_H          (H0)
#define DI_G          (H0)-t*(S0)+(p-1.0)*(V0)

#define DDI_GDT       -(S0)
#define DDI_GDP       (V0)

#define D2DI_GDT2     0.0
#define D2DI_GDTP     0.0
#define D2DI_GDP2     0.0

#define D3DI_GDT3     0.0
#define D3DI_GDT2DP   0.0
#define D3DI_GDTDP2   0.0
#define D3DI_GDP3     0.0

#define EN_S          (S0)+(SX7)-(SS2)+(SX7X7)-(SX7S2)+(SS2S2)+(SX7X7X7) \
                                            -(SX7X7S2)
#define EN_H          (H0)+(HX7)-(HS2)+(HX7X7)-(HX7S2)+(HS2S2)+(HX7X7X7) \
                                            -(HX7X7S2)+(HX7S2S2)
#define EN_G          (H0)+(HX7)-(HS2)+(HX7X7)-(HX7S2)+(HS2S2)+(HX7X7X7) \
                                            -(HX7X7S2)+(HX7S2S2) -t*((S0)+(SX7)-(SS2)+(SX7X7)- \
                                            (SX7S2)+(SS2S2)+(SX7X7X7)-(SX7X7S2)) +(p-1.0)*((V0)+ \
                                            (VX7)-(VS2)+(VX7X7)-(VX7S2)+(VS2S2)+(VX7X7X7)-(VX7X7S2)+ \
                                            (VX7S2S2))

#define DEN_GDT       -(EN_S)
#define DEN_GDP       (V0)+(VX7)-(VS2)+(VX7X7)-(VX7S2)+(VS2S2)+(VX7X7X7) \
                                            -(VX7X7S2)+(VX7S2S2)

#define D2EN_GDT2     0.0
#define D2EN_GDTP     0.0
#define D2EN_GDP2     0.0

#define D3EN_GDT3     0.0
#define D3EN_GDT2DP   0.0
#define D3EN_GDTDP2   0.0
#define D3EN_GDP3     0.0

#define HD_S          (S0)+(SX2)+(SX2X2)
#define HD_H          (H0)+(HX2)+(HX2X2)
#define HD_G          (H0)+(HX2)+(HX2X2) -t*((S0)+(SX2)+(SX2X2)) \
                                            + (p-1.0)*((V0)+(VX2)+(VX2X2))

#define DHD_GDT       -(HD_S)
#define DHD_GDP       (V0)+(VX2)+(VX2X2)

#define D2HD_GDT2     0.0
#define D2HD_GDTP     0.0
#define D2HD_GDP2     0.0

#define D3HD_GDT3     0.0
#define D3HD_GDT2DP   0.0
#define D3HD_GDTDP2   0.0
#define D3HD_GDP3     0.0

#define CA_S          (S0)+(SX3)+(SX3X3)
#define CA_H          (H0)+(HX3)+(HX3X3)
#define CA_G          (H0)+(HX3)+(HX3X3) -t*((S0)+(SX3)+(SX3X3)) \
                                            + (p-1.0)*((V0)+(VX3)+(VX3X3))

#define DCA_GDT       -(CA_S)
#define DCA_GDP       (V0)+(VX3)+(VX3X3)

#define D2CA_GDT2     0.0
#define D2CA_GDTP     0.0
#define D2CA_GDP2     0.0

#define D3CA_GDT3     0.0
#define D3CA_GDT2DP   0.0
#define D3CA_GDTDP2   0.0
#define D3CA_GDP3     0.0

#define CF_S          (S0)+(SX4)+(SX4X4)
#define CF_H          (H0)+(HX4)+(HX4X4)
#define CF_G          (H0)+(HX4)+(HX4X4) - t*((S0)+(SX4)+(SX4X4)) \
                                            + (p-1.0)*((V0)+(VX4)+(VX4X4))

#define DCF_GDT       -(CF_S)
#define DCF_GDP       (V0)+(VX4)+(VX4X4)

#define D2CF_GDT2     0.0
#define D2CF_GDTP     0.0
#define D2CF_GDP2     0.0

#define D3CF_GDT3     0.0
#define D3CF_GDT2DP   0.0
#define D3CF_GDTDP2   0.0
#define D3CF_GDP3     0.0

#define ES_S          -R*((1.0-s)*log(1.0-s) + (1.0+s)*log(1.0+s) \
                                                - 2.0*log(2.0)) + \
                                            (S0)+(SX5)+(SS1)*s+(SX5X5)+(SX5S1)*s+(SS1S1)*s*s
#define ES_H          (H0)+(HX5)+(HS1)*s+(HX5X5)+(HX5S1)*s+(HS1S1)*s*s
#define ES_G          R*t*((1.0-s)*log(1.0-s) + (1.0+s)*log(1.0+s) \
                                                - 2.0*log(2.0)) + \
                                            (H0)+(HX5)+(HS1)*s+(HX5X5)+(HX5S1)*s+(HS1S1)*s*s -t*( \
                                            (S0)+(SX5)+(SS1)*s+(SX5X5)+(SX5S1)*s+(SS1S1)*s*s) + \
                                            (p-1.0)*((V0)+(VX5)+(VS1)*s+(VX5X5)+(VX5S1)*s+(VS1S1)*s*s)

#define DES_GDS1      R*t*(log(1.0+s)-log(1.0-s)) + (HS1) + (HX5S1) \
                                            + (HS1S1)*s*2.0 - t*((SS1)+(SX5S1)+(SS1S1)*s*2.0) \
                                            + (p-1.0)*((VS1)+(VX5S1)+(VS1S1)*s*2.0)
#define DES_GDT       -(ES_S)
#define DES_GDP       (V0)+(VX5)+(VS1)*s+(VX5X5)+(VX5S1)*s+(VS1S1)*s*s

#define D2ES_GDS1S1   R*t*(1.0/(1.0+s) + 1.0/(1.0-s)) + (HS1S1)*2.0 \
                                            - t*(SS1S1)*2.0 + (p-1.0)*(VS1S1)*2.0
#define D2ES_GDS1DT   R*(log(1.0+s)-log(1.0-s)) - ((SS1)+(SX5S1)+(SS1S1)*s*2.0)
#define D2ES_GDS1DP   (VS1)+(VX5S1)+(VS1S1)*s*2.0
#define D2ES_GDT2     0.0
#define D2ES_GDTP     0.0
#define D2ES_GDP2     0.0

#define D3ES_GDS1S1S1 R*t*(1.0/SQUARE(1.0-s) - 1.0/SQUARE(1.0+s))
#define D3ES_GDS1S1DT R*(1.0/(1.0+s) + 1.0/(1.0-s)) - (SS1S1)*2.0
#define D3ES_GDS1S1DP (VS1S1)*2.0
#define D3ES_GDS1DT2  0.0
#define D3ES_GDS1DTDP 0.0
#define D3ES_GDS1DP2  0.0
#define D3ES_GDT3     0.0
#define D3ES_GDT2DP   0.0
#define D3ES_GDTDP2   0.0
#define D3ES_GDP3     0.0

#define JD_S          (S0) + 0.5*(SX3) - 0.5*(SX4) + 0.5*(SX5) + (SX6) \
                                            - (SS1) + 0.25*(SX3X3) - 0.25*(SX3X4) + 0.25*(SX3X5) \
                                            + 0.5*(SX3X6) - 0.5*(SX3S1) + 0.25*(SX4X4) \
                                            - 0.25*(SX4X5) - 0.5*(SX4X6) + 0.5*(SX4S1) \
                                            + 0.25*(SX5X5) + 0.5*(SX5X6) - 0.5*(SX5S1) \
                                            + (SX6X6) - (SX6S1) + (SS1S1)
#define JD_H          (H0) + 0.5*(HX3) - 0.5*(HX4) + 0.5*(HX5) + (HX6) \
                                            - (HS1) + 0.25*(HX3X3) - 0.25*(HX3X4) + 0.25*(HX3X5) \
                                            + 0.5*(HX3X6) - 0.5*(HX3S1) + 0.25*(HX4X4) \
                                            - 0.25*(HX4X5) - 0.5*(HX4X6) + 0.5*(HX4S1) \
                                            + 0.25*(HX5X5) + 0.5*(HX5X6) - 0.5*(HX5S1) \
                                            + (HX6X6) - (HX6S1) + (HS1S1)
#define JD_G          (H0)-t*(S0)+(p-1.0)*(V0) \
                                + 0.5*((HX3)-t*(SX3)+(p-1.0)*(VX3)) \
                                - 0.5*((HX4)-t*(SX4)+(p-1.0)*(VX4)) \
                                + 0.5*((HX5)-t*(SX5)+(p-1.0)*(VX5)) \
                                        + ((HX6)-t*(SX6)+(p-1.0)*(VX6)) \
                                        - ((HS1)-t*(SS1)+(p-1.0)*(VS1)) \
               + 0.25*((HX3X3)-t*(SX3X3)+(p-1.0)*(VX3X3)) \
               - 0.25*((HX3X4)-t*(SX3X4)+(p-1.0)*(VX3X4)) \
               + 0.25*((HX3X5)-t*(SX3X5)+(p-1.0)*(VX3X5)) \
                                + 0.5*((HX3X6)-t*(SX3X6)+(p-1.0)*(VX3X6)) \
                                - 0.5*((HX3S1)-t*(SX3S1)+(p-1.0)*(VX3S1)) \
               + 0.25*((HX4X4)-t*(SX4X4)+(p-1.0)*(VX4X4)) \
               - 0.25*((HX4X5)-t*(SX4X5)+(p-1.0)*(VX4X5)) \
                                - 0.5*((HX4X6)-t*(SX4X6)+(p-1.0)*(VX4X6)) \
                                + 0.5*((HX4S1)-t*(SX4S1)+(p-1.0)*(VX4S1)) \
               + 0.25*((HX5X5)-t*(SX5X5)+(p-1.0)*(VX5X5)) \
                                + 0.5*((HX5X6)-t*(SX5X6)+(p-1.0)*(VX5X6)) \
                                - 0.5*((HX5S1)-t*(SX5S1)+(p-1.0)*(VX5S1)) \
                                        + ((HX6X6)-t*(SX6X6)+(p-1.0)*(VX6X6)) \
                                        - ((HX6S1)-t*(SX6S1)+(p-1.0)*(VX6S1)) \
                                        + ((HS1S1)-t*(SS1S1)+(p-1.0)*(VS1S1))

#define DJD_GDT       -(JD_S)
#define DJD_GDP       (V0) + 0.5*(VX3) - 0.5*(VX4) + 0.5*(VX5) + (VX6) \
                                            - (VS1) + 0.25*(VX3X3) - 0.25*(VX3X4) + 0.25*(VX3X5) \
                                            + 0.5*(VX3X6) - 0.5*(VX3S1) + 0.25*(VX4X4) \
                                            - 0.25*(VX4X5) - 0.5*(VX4X6) + 0.5*(VX4S1) \
                                            + 0.25*(VX5X5) + 0.5*(VX5X6) - 0.5*(VX5S1) \
                                            + (VX6X6) - (VX6S1) + (VS1S1)

#define D2JD_GDT2     0.0
#define D2JD_GDTP     0.0
#define D2JD_GDP2     0.0

#define D3JD_GDT3     0.0
#define D3JD_GDT2DP   0.0
#define D3JD_GDTDP2   0.0
#define D3JD_GDP3     0.0

/*
 In the endmember properties routines, a static variable clino overrides the
 current global variable clino. This insures that the solution properties
 returned by the global functions in this file are always reference to a
 monoclinic state (i.e. the thermodynamic constants stored
 SOLID_STRUCT_DATA.H and used in GIBBS.C refer to the monoclinic structure).
 If the C 2/c -> Pbam transition occurs, the solution properties will
 contain the standard state contribution.
*/

static void
pureOrder(int mask, double t, double p,
            double *s,   /* s1       BINARY MASK: 000001 */
            double *dt,  /* ds1/dt   BINARY MASK: 000010 */
            double *dp,  /* ds1/dp   BINARY MASK: 000100 */
            double *dt2, /* d2s1/dt2 BINARY MASK: 001000 */
            double *dtp, /* d2s1/dtp BINARY MASK: 010000 */
            double *dp2  /* d2s1/dp2 BINARY MASK: 100000 */
            )
{
    double tOld   = getTOldPure();
    double pOld   = getPOldPure();
    double sOld   = getSOldPure();
    double d2gds2 = getD2gds2Pure();

    static const int clino = TRUE;

    if ( (t != tOld) || (p != pOld) ) {
        double dgds, sNew;
        int iter = 0;
        sOld = 2.0;
        sNew = 0.5;
        while ((ABS(sNew-sOld) > 10.0*DBL_EPSILON) && (iter < MAX_ITER)) {
            double s;
            s      = sNew;
            dgds   = DES_GDS1;
            d2gds2 = D2ES_GDS1S1;
            sOld   = s;
            s     += - dgds/d2gds2;
            s      = MIN(s,  1.0 - DBL_EPSILON);
            s      = MAX(s, -1.0 + DBL_EPSILON);
            sNew   = s;
                iter++;
        }
        tOld = t;
        pOld = p;

        setTOldPure(tOld);
        setPOldPure(pOld);
        setSOldPure(sOld);
        setD2gds2Pure(d2gds2);
    }

    if (mask & FIRST  ) {   /* return s        */
        *s = sOld;
    }

    if (mask & SECOND ) {   /* compute ds/dt:  */
        double s = sOld;
        *dt = - (D2ES_GDS1DT)/d2gds2;
    }

    if (mask & THIRD  ) {   /* compute ds/dp:  */
        double s = sOld;
        *dp = - (D2ES_GDS1DP)/d2gds2;
    }

    if (mask & FOURTH ) {   /* compute d2s/dt2 */
        double s = sOld;
        double dsdt = - (D2ES_GDS1DT)/d2gds2;
        *dt2  = - ((D3ES_GDS1DT2) + 2.0*(D3ES_GDS1S1DT)*dsdt
                    + (D3ES_GDS1S1S1)*dsdt*dsdt)/d2gds2;
    }

    if (mask & FIFTH  ) {   /* compute d2s/dtp */
        double s = sOld;
        double dsdt = - (D2ES_GDS1DT)/d2gds2;
        double dsdp = - (D2ES_GDS1DP)/d2gds2;
        *dtp = - ((D3ES_GDS1DTDP) + (D3ES_GDS1S1DT)*dsdp + (D3ES_GDS1S1DP)*dsdt
         + (D3ES_GDS1S1S1)*dsdt*dsdp)/d2gds2;

    }

    if (mask & SIXTH  ) {   /* compute d2s/dp2 */
        double s = sOld;
        double dsdp = - (D2ES_GDS1DP)/d2gds2;
        *dp2 = - ((D3ES_GDS1DP2) + 2.0*(D3ES_GDS1S1DP)*dsdp
         + (D3ES_GDS1S1S1)*dsdp*dsdp)/d2gds2;
    }
}

static void
purePyx(int mask, double t, double p,
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
    double s;
    int i;
    static const int clino = TRUE;

    pureOrder(FIRST, t, p,
                        &s,              (double *) NULL, (double *) NULL, (double *) NULL,
                        (double *) NULL, (double *) NULL);

    if (mask & FIRST) {
        a[0] = DI_G;
        a[0] = exp(a[0]/(R*t));
        a[1] = EN_G;
        a[1] = exp(a[1]/(R*t));
        a[2] = HD_G;
        a[2] = exp(a[2]/(R*t));
        a[3] = CA_G;
        a[3] = exp(a[3]/(R*t));
        a[4] = CF_G;
        a[4] = exp(a[4]/(R*t));
        a[5] = ES_G;
        a[5] = exp(a[5]/(R*t));
        a[6] = JD_G;
        a[6] = exp(a[6]/(R*t));
    }

    if (mask & SECOND) {
        mu[0] = DI_G;
        mu[1] = EN_G;
        mu[2] = HD_G;
        mu[3] = CA_G;
        mu[4] = CF_G;
        mu[5] = ES_G;
        mu[6] = JD_G;
    }

    if (mask & THIRD) {
        gmix[0] = DI_G;
        gmix[1] = EN_G;
        gmix[2] = HD_G;
        gmix[3] = CA_G;
        gmix[4] = CF_G;
        gmix[5] = ES_G;
        gmix[6] = JD_G;
    }

    if (mask & FOURTH) {
        hmix[0] = (DI_G) + t*(DI_S);
        hmix[1] = (EN_G) + t*(EN_S);
        hmix[2] = (HD_G) + t*(HD_S);
        hmix[3] = (CA_G) + t*(CA_S);
        hmix[4] = (CF_G) + t*(CF_S);
        hmix[5] = (ES_G) + t*(ES_S);
        hmix[6] = (JD_G) + t*(JD_S);
    }

    if (mask & FIFTH) {
        smix[0] = DI_S;
        smix[1] = EN_S;
        smix[2] = HD_S;
        smix[3] = CA_S;
        smix[4] = CF_S;
        smix[5] = ES_S;
        smix[6] = JD_S;
    }

    if (mask & SIXTH) {
        double dsdt;

        pureOrder(SECOND, t, p,
                        (double *) NULL, &dsdt,           (double *) NULL, (double *) NULL,
                        (double *) NULL, (double *) NULL);

        cpmix[0]  = D2DI_GDT2;
        cpmix[1]  = D2EN_GDT2;
        cpmix[2]  = D2HD_GDT2;
        cpmix[3]  = D2CA_GDT2;
        cpmix[4]  = D2CF_GDT2;
        cpmix[5]  = D2ES_GDT2;
        cpmix[5] += 2.0*(D2ES_GDS1DT)*dsdt;
        cpmix[5] += (D2ES_GDS1S1)*SQUARE(dsdt);
        cpmix[6]  = D2JD_GDT2;

        for (i=0; i<NA; i++) cpmix[i] *= -t;
    }

    if(mask & SEVENTH) {
        double d2gdt2[NA], d3gdt3[NA], dsdt, d2sdt2;

        d2gdt2[0] = D2DI_GDT2; d2gdt2[1] = D2EN_GDT2; d2gdt2[2] = D2HD_GDT2;
        d2gdt2[3] = D2CA_GDT2; d2gdt2[4] = D2CF_GDT2; d2gdt2[5] = D2ES_GDT2;
        d2gdt2[6] = D2JD_GDT2;

        d3gdt3[0] = D3DI_GDT3; d3gdt3[1] = D3EN_GDT3; d3gdt3[2] = D3HD_GDT3;
        d3gdt3[3] = D3CA_GDT3; d3gdt3[4] = D3CF_GDT3; d3gdt3[5] = D3ES_GDT3;
        d3gdt3[6] = D3JD_GDT3;

        pureOrder(SECOND | FOURTH, t, p,
                        (double *) NULL, &dsdt,           (double *) NULL, &d2sdt2,
                        (double *) NULL, (double *) NULL);

        d2gdt2[5] += 2.0*(D2ES_GDS1DT)*dsdt;
        d2gdt2[5] += (D2ES_GDS1S1)*SQUARE(dsdt);

        d3gdt3[5] += 3.0*(D3ES_GDS1DT2)*dsdt;
        d3gdt3[5] += 3.0*(D2ES_GDS1DT)*d2sdt2;
        d3gdt3[5] += 3.0*(D2ES_GDS1S1)*dsdt*d2sdt2;
        d3gdt3[5] += 3.0*(D3ES_GDS1S1DT)*dsdt*dsdt;
        d3gdt3[5] += (D3ES_GDS1S1S1)*CUBE(dsdt);

        for (i=0; i<NA; i++) cpmixdt[i] = -t*d3gdt3[i] - d2gdt2[i];
    }

    if (mask & EIGHTH) {
        vmix[0] = DDI_GDP;
        vmix[1] = DEN_GDP;
        vmix[2] = DHD_GDP;
        vmix[3] = DCA_GDP;
        vmix[4] = DCF_GDP;
        vmix[5] = DES_GDP;
        vmix[6] = DJD_GDP;
    }

    if(mask & NINTH) {
        double dsdt, dsdp;

        pureOrder(SECOND | THIRD, t, p,
                        (double *) NULL, &dsdt,           &dsdp,           (double *) NULL,
                        (double *) NULL, (double *) NULL);

        vmixdt[0]  = D2DI_GDTP;
        vmixdt[1]  = D2EN_GDTP;
        vmixdt[2]  = D2HD_GDTP;
        vmixdt[3]  = D2CA_GDTP;
        vmixdt[4]  = D2CF_GDTP;
        vmixdt[5]  = D2ES_GDTP;
        vmixdt[5] += (D2ES_GDS1DT)*dsdp;
        vmixdt[5] += (D2ES_GDS1DP)*dsdt;
        vmixdt[5] += (D2ES_GDS1S1)*dsdt*dsdp;
        vmixdt[6]  = D2JD_GDTP;
    }

    if(mask & TENTH) {
        double dsdp;

        pureOrder(THIRD, t, p,
                        (double *) NULL, (double *) NULL, &dsdp,           (double *) NULL,
                        (double *) NULL, (double *) NULL);

        vmixdp[0]  = D2DI_GDP2;
        vmixdp[1]  = D2EN_GDP2;
        vmixdp[2]  = D2HD_GDP2;
        vmixdp[3]  = D2CA_GDP2;
        vmixdp[4]  = D2CF_GDP2;
        vmixdp[5]  = D2ES_GDP2;
        vmixdp[5] += 2.0*(D2ES_GDS1DP)*dsdp;
        vmixdp[5] += (D2ES_GDS1S1)*dsdp*dsdp;
        vmixdp[6]  = D2JD_GDP2;
    }

    if(mask & ELEVENTH) {
        double dsdt, dsdp, d2sdt2, d2sdtdp;

        pureOrder(SECOND | THIRD | FOURTH | FIFTH, t, p,
                        (double *) NULL, &dsdt,           &dsdp,           &d2sdt2,
                        &d2sdtdp,        (double *) NULL);

        vmixdt2[0]  = D3DI_GDT2DP;
        vmixdt2[1]  = D3EN_GDT2DP;
        vmixdt2[2]  = D3HD_GDT2DP;
        vmixdt2[3]  = D3CA_GDT2DP;
        vmixdt2[4]  = D3CF_GDT2DP;
        vmixdt2[5]  = D3ES_GDT2DP;
        vmixdt2[5] += (D3ES_GDS1DT2)*dsdp;
        vmixdt2[5] += 2.0*(D2ES_GDS1DT)*d2sdtdp;
        vmixdt2[5] += (D2ES_GDS1DP)*d2sdt2;
        vmixdt2[5] += 2.0*(D3ES_GDS1DTDP)*dsdt;
        vmixdt2[5] += 2.0*(D3ES_GDS1S1DT)*dsdt*dsdp;
        vmixdt2[5] += (D2ES_GDS1S1)*(d2sdt2*dsdp + 2.0*dsdt*d2sdtdp);
        vmixdt2[5] += (D3ES_GDS1S1DP)*dsdt*dsdt;
        vmixdt2[5] += (D3ES_GDS1S1S1)*dsdt*dsdt*dsdp;
        vmixdt2[6]  = D3JD_GDT2DP;
    }

    if(mask & TWELFTH) {
        double dsdt, dsdp, d2sdtdp, d2sdp2;

        pureOrder(SECOND | THIRD | FIFTH | SIXTH, t, p,
                        (double *) NULL, &dsdt,           &dsdp,           (double *) NULL,
                        &d2sdtdp,        &d2sdp2);

        vmixdtdp[0]  = D3DI_GDTDP2;
        vmixdtdp[1]  = D3EN_GDTDP2;
        vmixdtdp[2]  = D3HD_GDTDP2;
        vmixdtdp[3]  = D3CA_GDTDP2;
        vmixdtdp[4]  = D3CF_GDTDP2;
        vmixdtdp[5]  = D3ES_GDTDP2;
        vmixdtdp[5] += 2.0*(D3ES_GDS1DTDP)*dsdp;
        vmixdtdp[5] += (D2ES_GDS1DT)*d2sdp2;
        vmixdtdp[5] += 2.0*(D2ES_GDS1DP)*d2sdtdp;
        vmixdtdp[5] += (D3ES_GDS1DP2)*dsdt;
        vmixdtdp[5] += 2.0*(D3ES_GDS1S1DP)*dsdt*dsdp;
        vmixdtdp[5] += (D2ES_GDS1S1)*(dsdt*d2sdp2 + 2.0*d2sdtdp*dsdp);
        vmixdtdp[5] += (D3ES_GDS1S1DT)*dsdp*dsdp;
        vmixdtdp[5] += (D3ES_GDS1S1S1)*dsdt*dsdp*dsdp;
        vmixdtdp[6]  = D3JD_GDTDP2;
    }

    if(mask & THIRTEENTH) {
        double dsdp, d2sdp2;

        pureOrder(THIRD | SIXTH, t, p,
                        (double *) NULL, (double *) NULL, &dsdp,           (double *) NULL,
                        (double *) NULL, &d2sdp2);

        vmixdp2[0]  = D3DI_GDP3;
        vmixdp2[1]  = D3EN_GDP3;
        vmixdp2[2]  = D3HD_GDP3;
        vmixdp2[3]  = D3CA_GDP3;
        vmixdp2[4]  = D3CF_GDP3;
        vmixdp2[5]  = D3ES_GDP3;
        vmixdp2[5] += 3.0*(D3ES_GDS1DP2)*dsdp;
        vmixdp2[5] += 3.0*(D2ES_GDS1DP)*d2sdp2;
        vmixdp2[5] += 3.0*(D2ES_GDS1S1)*dsdp*d2sdp2;
        vmixdp2[5] += 3.0*(D3ES_GDS1S1DP)*dsdp*dsdp;
        vmixdp2[5] += (D3ES_GDS1S1S1)*dsdp*dsdp*dsdp;
        vmixdp2[6]  = D3JD_GDP3;
    }

}

#undef DI_S
#undef DI_H
#undef DI_G
#undef DDI_GDT
#undef DDI_GDP
#undef D2DI_GDT2
#undef D2DI_GDTP
#undef D2DI_GDP2
#undef D3DI_GDT3
#undef D3DI_GDT2DP
#undef D3DI_GDTDP2
#undef D3DI_GDP3

#undef EN_S
#undef EN_H
#undef EN_G
#undef DEN_GDT
#undef DEN_GDP
#undef D2EN_GDT2
#undef D2EN_GDTP
#undef D2EN_GDP2
#undef D3EN_GDT3
#undef D3EN_GDT2DP
#undef D3EN_GDTDP2
#undef D3EN_GDP3

#undef HD_S
#undef HD_H
#undef HD_G
#undef DHD_GDT
#undef DHD_GDP
#undef D2HD_GDT2
#undef D2HD_GDTP
#undef D2HD_GDP2
#undef D3HD_GDT3
#undef D3HD_GDT2DP
#undef D3HD_GDTDP2
#undef D3HD_GDP3

#undef CA_S
#undef CA_H
#undef CA_G
#undef DCA_GDT
#undef DCA_GDP
#undef D2CA_GDT2
#undef D2CA_GDTP
#undef D2CA_GDP2
#undef D3CA_GDT3
#undef D3CA_GDT2DP
#undef D3CA_GDTDP2
#undef D3CA_GDP3

#undef CF_S
#undef CF_H
#undef CF_G
#undef DCF_GDT
#undef DCF_GDP
#undef D2CF_GDT2
#undef D2CF_GDTP
#undef D2CF_GDP2
#undef D3CF_GDT3
#undef D3CF_GDT2DP
#undef D3CF_GDTDP2
#undef D3CF_GDP3

#undef ES_S
#undef ES_H
#undef ES_G
#undef DES_GDS1
#undef DES_GDT
#undef DES_GDP
#undef D2ES_GDS1S1
#undef D2ES_GDS1DT
#undef D2ES_GDS1DP
#undef D2ES_GDT2
#undef D2ES_GDTP
#undef D2ES_GDP2
#undef D3ES_GDS1S1S1
#undef D3ES_GDS1S1DT
#undef D3ES_GDS1S1DP
#undef D3ES_GDS1DT2
#undef D3ES_GDS1DTDP
#undef D3ES_GDS1DP2
#undef D3ES_GDT3
#undef D3ES_GDT2DP
#undef D3ES_GDTDP2
#undef D3ES_GDP3

#undef JD_S
#undef JD_H
#undef JD_G
#undef DJD_GDT
#undef DJD_GDP
#undef D2JD_GDT2
#undef D2JD_GDTP
#undef D2JD_GDP2
#undef D3JD_GDT3
#undef D3JD_GDTDP
#undef D3JD_GDTDP2
#undef D3JD_GDP3

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conOpx defines the conversion from m[i], to r[j]
 */
                   /* Order: X2, X3. X4, X5, X6, X7 */
#define FR2(i)     (i == 2) ? 1.0 - r[0] : - r[0]
#define FR3(i)     (i == 3) ? 1.0 - r[1] : ((i == 6) ?   0.5 - r[1] : - r[1])
#define FR4(i)     (i == 4) ? 1.0 - r[2] : ((i == 6) ? - 0.5 - r[2] : - r[2])
#define FR5(i)     (i == 5) ? 1.0 - r[3] : ((i == 6) ?   0.5 - r[3] : - r[3])
#define FR6(i)     (i == 6) ? 1.0 - r[4] : - r[4]
#define FR7(i)     (i == 1) ? 1.0 - r[5] : - r[5]

                   /* Order: S1, S2 */
#define FS1(i)     (i == 5) ?   1.0 - s[0] : ((i == 6) ? - 1.0 - s[0] : - s[0])
#define FS2(i)     (i == 1) ? - 1.0 - s[1] : - s[1]

#define DFR2DR2(i) - 1.0
#define DFR3DR3(i) - 1.0
#define DFR4DR4(i) - 1.0
#define DFR5DR5(i) - 1.0
#define DFR6DR6(i) - 1.0
#define DFR7DR7(i) - 1.0

#define DFS1DS1(i) - 1.0
#define DFS2DS2(i) - 1.0

#define ENDMEMBERS  (  (1.0-r[0]-r[1]-r[2]-r[3]-r[4]/2.0-r[5])*ends[0] \
                     + r[5]*ends[1] + r[0]*ends[2] + (r[1]-r[4]/2.0)*ends[3] \
                     + (r[2]+r[4]/2.0)*ends[4] + (r[3]-r[4]/2.0)*ends[5] \
                     + r[4]*ends[6] )

#define DENDDR0  (ends[2] - ends[0])
#define DENDDR1  (ends[3] - ends[0])
#define DENDDR2  (ends[4] - ends[0])
#define DENDDR3  (ends[5] - ends[0])
#define DENDDR4  (ends[6] + 0.5*(ends[4] - ends[0] - ends[3] - ends[5]))
#define DENDDR5  (ends[1] - ends[0])

/*
 * Global (to this file): derivative definitions
 */

#define SIC -R*(xmg2m1*log(xmg2m1) \
                            + xfe2m1*log(xfe2m1) \
                            + xal3m1*log(xal3m1) \
                            + xfe3m1*log(xfe3m1) \
                            + xti4m1*log(xti4m1) \
                            + (xmg2m1+xfe2m1-xti4m1)*log(xmg2m1+xfe2m1-xti4m1) \
                            - (1.0-xti4m1)*log(1.0-xti4m1) \
                            - (xmg2m1+xfe2m1)*log(xmg2m1+xfe2m1) \
                            - 2.0*(1.0-xsi4tet)*log(1.0-xsi4tet) \
                            + 2.0*xal3tet*log(xal3tet) \
                            + 2.0*xfe3tet*log(xfe3tet) \
                            + xca2m2*log(xca2m2) \
                            + xna1m2*log(xna1m2) \
                            + xmg2m2*log(xmg2m2) \
                            + xfe2m2*log(xfe2m2) \
                            + (1.0-xmg2m1-xfe2m1-xna1m2)*log(1.0-xmg2m1-xfe2m1-xna1m2) \
                            - (1.0-xmg2m1-xfe2m1)*log(1.0-xmg2m1-xfe2m1) \
                            - (1.0-xna1m2)*log(1.0-xna1m2) )
#define S   (SIC) + (S0) + \
                        (SX2)*r[0] + (SX3)*r[1] + (SX4)*r[2] + (SX5)*r[3] + \
                        (SX6)*r[4] + (SX7)*r[5] + (SS1)*s[0] + (SS2)*s[1] + \
                        (SX2X2)*r[0]*r[0] + (SX2X3)*r[0]*r[1] + (SX2X4)*r[0]*r[2] + \
                        (SX2X5)*r[0]*r[3] + (SX2X6)*r[0]*r[4] + (SX2X7)*r[0]*r[5] + \
                        (SX2S1)*r[0]*s[0] + (SX2S2)*r[0]*s[1] + (SX3X3)*r[1]*r[1] + \
                        (SX3X4)*r[1]*r[2] + (SX3X5)*r[1]*r[3] + (SX3X6)*r[1]*r[4] + \
                        (SX3X7)*r[1]*r[5] + (SX3S1)*r[1]*s[0] + (SX3S2)*r[1]*s[1] + \
                        (SX4X4)*r[2]*r[2] + (SX4X5)*r[2]*r[3] + (SX4X6)*r[2]*r[4] + \
                        (SX4X7)*r[2]*r[5] + (SX4S1)*r[2]*s[0] + (SX4S2)*r[2]*s[1] + \
                        (SX5X5)*r[3]*r[3] + (SX5X6)*r[3]*r[4] + (SX5X7)*r[3]*r[5] + \
                        (SX5S1)*r[3]*s[0] + (SX5S2)*r[3]*s[1] + (SX6X6)*r[4]*r[4] + \
                        (SX6X7)*r[4]*r[5] + (SX6S1)*r[4]*s[0] + (SX6S2)*r[4]*s[1] + \
                        (SX7X7)*r[5]*r[5] + (SX7S1)*r[5]*s[0] + (SX7S2)*r[5]*s[1] + \
                        (SS1S1)*s[0]*s[0] + (SS1S2)*s[0]*s[1] + (SS2S2)*s[1]*s[1] + \
                        (SX2X7X7)*r[0]*r[5]*r[5] + (SX3X7X7)*r[1]*r[5]*r[5] + \
                        (SX3X7S2)*r[1]*r[5]*s[1] + (SX4X7X7)*r[2]*r[5]*r[5] + \
                        (SX4X7S2)*r[2]*r[5]*s[1] + (SX5X7X7)*r[3]*r[5]*r[5] + \
                        (SX5X7S2)*r[3]*r[5]*s[1] + (SX6X7X7)*r[4]*r[5]*r[5] + \
                        (SX6X7S2)*r[4]*r[5]*s[1] + (SX7X7X7)*r[5]*r[5]*r[5] + \
                        (SX7X7S2)*r[5]*r[5]*s[1]
#define H   (H0) + \
                        (HX2)*r[0] + (HX3)*r[1] + (HX4)*r[2] + (HX5)*r[3] + \
                        (HX6)*r[4] + (HX7)*r[5] + (HS1)*s[0] + (HS2)*s[1] + \
                        (HX2X2)*r[0]*r[0] + (HX2X3)*r[0]*r[1] + (HX2X4)*r[0]*r[2] + \
                        (HX2X5)*r[0]*r[3] + (HX2X6)*r[0]*r[4] + (HX2X7)*r[0]*r[5] + \
                        (HX2S1)*r[0]*s[0] + (HX2S2)*r[0]*s[1] + (HX3X3)*r[1]*r[1] + \
                        (HX3X4)*r[1]*r[2] + (HX3X5)*r[1]*r[3] + (HX3X6)*r[1]*r[4] + \
                        (HX3X7)*r[1]*r[5] + (HX3S1)*r[1]*s[0] + (HX3S2)*r[1]*s[1] + \
                        (HX4X4)*r[2]*r[2] + (HX4X5)*r[2]*r[3] + (HX4X6)*r[2]*r[4] + \
                        (HX4X7)*r[2]*r[5] + (HX4S1)*r[2]*s[0] + (HX4S2)*r[2]*s[1] + \
                        (HX5X5)*r[3]*r[3] + (HX5X6)*r[3]*r[4] + (HX5X7)*r[3]*r[5] + \
                        (HX5S1)*r[3]*s[0] + (HX5S2)*r[3]*s[1] + (HX6X6)*r[4]*r[4] + \
                        (HX6X7)*r[4]*r[5] + (HX6S1)*r[4]*s[0] + (HX6S2)*r[4]*s[1] + \
                        (HX7X7)*r[5]*r[5] + (HX7S1)*r[5]*s[0] + (HX7S2)*r[5]*s[1] + \
                        (HS1S1)*s[0]*s[0] + (HS1S2)*s[0]*s[1] + (HS2S2)*s[1]*s[1] + \
                        (HX2X2X7)*r[0]*r[0]*r[5] + (HX2X2S2)*r[0]*r[0]*s[1] + \
                        (HX2X3X7)*r[0]*r[1]*r[5] + (HX2X3S2)*r[0]*r[1]*s[1] + \
                        (HX2X4X7)*r[0]*r[2]*r[5] + (HX2X4S2)*r[0]*r[2]*s[1] + \
                        (HX2X5X7)*r[0]*r[3]*r[5] + (HX2X5S2)*r[0]*r[3]*s[1] + \
                        (HX2X6X7)*r[0]*r[4]*r[5] + (HX2X6S2)*r[0]*r[4]*s[1] + \
                        (HX2X7X7)*r[0]*r[5]*r[5] + (HX2X7S2)*r[0]*r[5]*s[1] + \
                        (HX2S2S2)*r[0]*s[1]*s[1] + (HX3X3X7)*r[1]*r[1]*r[5] + \
                        (HX3X3S2)*r[1]*r[1]*s[1] + (HX3X4X7)*r[1]*r[2]*r[5] + \
                        (HX3X4S2)*r[1]*r[2]*s[1] + (HX3X5X7)*r[1]*r[3]*r[5] + \
                        (HX3X5S2)*r[1]*r[3]*s[1] + (HX3X6X7)*r[1]*r[4]*r[5] + \
                        (HX3X6S2)*r[1]*r[4]*s[1] + (HX3X7X7)*r[1]*r[5]*r[5] + \
                        (HX3X7S2)*r[1]*r[5]*s[1] + (HX3S2S2)*r[1]*s[1]*s[1] + \
                        (HX4X4X7)*r[2]*r[2]*r[5] + (HX4X4S2)*r[2]*r[2]*s[1] + \
                        (HX4X5X7)*r[2]*r[3]*r[5] + (HX4X5S2)*r[2]*r[3]*s[1] + \
                        (HX4X6X7)*r[2]*r[4]*r[5] + (HX4X6S2)*r[2]*r[4]*s[1] + \
                        (HX4X7X7)*r[2]*r[5]*r[5] + (HX4X7S2)*r[2]*r[5]*s[1] + \
                        (HX4S2S2)*r[2]*s[1]*s[1] + (HX5X5X7)*r[3]*r[3]*r[5] + \
                        (HX5X5S2)*r[3]*r[3]*s[1] + (HX5X6X7)*r[3]*r[4]*r[5] + \
                        (HX5X6S2)*r[3]*r[4]*s[1] + (HX5X7X7)*r[3]*r[5]*r[5] + \
                        (HX5X7S2)*r[3]*r[5]*s[1] + (HX5S2S2)*r[3]*s[1]*s[1] + \
                        (HX6X6X7)*r[4]*r[4]*r[5] + (HX6X6S2)*r[4]*r[4]*s[1] + \
                        (HX6X7X7)*r[4]*r[5]*r[5] + (HX6X7S2)*r[4]*r[5]*s[1] + \
                        (HX6S2S2)*r[4]*s[1]*s[1] + (HX7X7X7)*r[5]*r[5]*r[5] + \
                        (HX7X7S2)*r[5]*r[5]*s[1] + (HX7S2S2)*r[5]*s[1]*s[1]
#define V   (V0) + \
                        (VX2)*r[0] + (VX3)*r[1] + (VX4)*r[2] + (VX5)*r[3] + \
                        (VX6)*r[4] + (VX7)*r[5] + (VS1)*s[0] + (VS2)*s[1] + \
                        (VX2X2)*r[0]*r[0] + (VX2X3)*r[0]*r[1] + (VX2X4)*r[0]*r[2] + \
                        (VX2X5)*r[0]*r[3] + (VX2X6)*r[0]*r[4] + (VX2X7)*r[0]*r[5] + \
                        (VX2S1)*r[0]*s[0] + (VX2S2)*r[0]*s[1] + (VX3X3)*r[1]*r[1] + \
                        (VX3X4)*r[1]*r[2] + (VX3X5)*r[1]*r[3] + (VX3X6)*r[1]*r[4] + \
                        (VX3X7)*r[1]*r[5] + (VX3S1)*r[1]*s[0] + (VX3S2)*r[1]*s[1] + \
                        (VX4X4)*r[2]*r[2] + (VX4X5)*r[2]*r[3] + (VX4X6)*r[2]*r[4] + \
                        (VX4X7)*r[2]*r[5] + (VX4S1)*r[2]*s[0] + (VX4S2)*r[2]*s[1] + \
                        (VX5X5)*r[3]*r[3] + (VX5X6)*r[3]*r[4] + (VX5X7)*r[3]*r[5] + \
                        (VX5S1)*r[3]*s[0] + (VX5S2)*r[3]*s[1] + (VX6X6)*r[4]*r[4] + \
                        (VX6X7)*r[4]*r[5] + (VX6S1)*r[4]*s[0] + (VX6S2)*r[4]*s[1] + \
                        (VX7X7)*r[5]*r[5] + (VX7S1)*r[5]*s[0] + (VX7S2)*r[5]*s[1] + \
                        (VS1S1)*s[0]*s[0] + (VS1S2)*s[0]*s[1] + (VS2S2)*s[1]*s[1] + \
                        (VX2X2X7)*r[0]*r[0]*r[5] + (VX2X2S2)*r[0]*r[0]*s[1] + \
                        (VX2X3X7)*r[0]*r[1]*r[5] + (VX2X3S2)*r[0]*r[1]*s[1] + \
                        (VX2X4X7)*r[0]*r[2]*r[5] + (VX2X4S2)*r[0]*r[2]*s[1] + \
                        (VX2X5X7)*r[0]*r[3]*r[5] + (VX2X5S2)*r[0]*r[3]*s[1] + \
                        (VX2X6X7)*r[0]*r[4]*r[5] + (VX2X6S2)*r[0]*r[4]*s[1] + \
                        (VX2X7X7)*r[0]*r[5]*r[5] + (VX2X7S2)*r[0]*r[5]*s[1] + \
                        (VX2S2S2)*r[0]*s[1]*s[1] + (VX3X3X7)*r[1]*r[1]*r[5] + \
                        (VX3X3S2)*r[1]*r[1]*s[1] + (VX3X4X7)*r[1]*r[2]*r[5] + \
                        (VX3X4S2)*r[1]*r[2]*s[1] + (VX3X5X7)*r[1]*r[3]*r[5] + \
                        (VX3X5S2)*r[1]*r[3]*s[1] + (VX3X6X7)*r[1]*r[4]*r[5] + \
                        (VX3X6S2)*r[1]*r[4]*s[1] + (VX3X7X7)*r[1]*r[5]*r[5] + \
                        (VX3X7S2)*r[1]*r[5]*s[1] + (VX3S2S2)*r[1]*s[1]*s[1] + \
                        (VX4X4X7)*r[2]*r[2]*r[5] + (VX4X4S2)*r[2]*r[2]*s[1] + \
                        (VX4X5X7)*r[2]*r[3]*r[5] + (VX4X5S2)*r[2]*r[3]*s[1] + \
                        (VX4X6X7)*r[2]*r[4]*r[5] + (VX4X6S2)*r[2]*r[4]*s[1] + \
                        (VX4X7X7)*r[2]*r[5]*r[5] + (VX4X7S2)*r[2]*r[5]*s[1] + \
                        (VX4S2S2)*r[2]*s[1]*s[1] + (VX5X5X7)*r[3]*r[3]*r[5] + \
                        (VX5X5S2)*r[3]*r[3]*s[1] + (VX5X6X7)*r[3]*r[4]*r[5] + \
                        (VX5X6S2)*r[3]*r[4]*s[1] + (VX5X7X7)*r[3]*r[5]*r[5] + \
                        (VX5X7S2)*r[3]*r[5]*s[1] + (VX5S2S2)*r[3]*s[1]*s[1] + \
                        (VX6X6X7)*r[4]*r[4]*r[5] + (VX6X6S2)*r[4]*r[4]*s[1] + \
                        (VX6X7X7)*r[4]*r[5]*r[5] + (VX6X7S2)*r[4]*r[5]*s[1] + \
                        (VX6S2S2)*r[4]*s[1]*s[1] + (VX7X7X7)*r[5]*r[5]*r[5] + \
                        (VX7X7S2)*r[5]*r[5]*s[1] + (VX7S2S2)*r[5]*s[1]*s[1]
#define G   -t*(SIC) + (H0)-t*(S0)+(p-1.0)*(V0) + \
                        ((HX2)-t*(SX2)+(p-1.0)*(VX2))*r[0] + \
                        ((HX3)-t*(SX3)+(p-1.0)*(VX3))*r[1] + \
                        ((HX4)-t*(SX4)+(p-1.0)*(VX4))*r[2] + \
                        ((HX5)-t*(SX5)+(p-1.0)*(VX5))*r[3] + \
                        ((HX6)-t*(SX6)+(p-1.0)*(VX6))*r[4] + \
                        ((HX7)-t*(SX7)+(p-1.0)*(VX7))*r[5] + \
                        ((HS1)-t*(SS1)+(p-1.0)*(VS1))*s[0] + \
                        ((HS2)-t*(SS2)+(p-1.0)*(VS2))*s[1] + \
                        ((HX2X2)-t*(SX2X2)+(p-1.0)*(VX2X2))*r[0]*r[0] + \
                        ((HX2X3)-t*(SX2X3)+(p-1.0)*(VX2X3))*r[0]*r[1] + \
                        ((HX2X4)-t*(SX2X4)+(p-1.0)*(VX2X4))*r[0]*r[2] + \
                        ((HX2X5)-t*(SX2X5)+(p-1.0)*(VX2X5))*r[0]*r[3] + \
                        ((HX2X6)-t*(SX2X6)+(p-1.0)*(VX2X6))*r[0]*r[4] + \
                        ((HX2X7)-t*(SX2X7)+(p-1.0)*(VX2X7))*r[0]*r[5] + \
                        ((HX2S1)-t*(SX2S1)+(p-1.0)*(VX2S1))*r[0]*s[0] + \
                        ((HX2S2)-t*(SX2S2)+(p-1.0)*(VX2S2))*r[0]*s[1] + \
                        ((HX3X3)-t*(SX3X3)+(p-1.0)*(VX3X3))*r[1]*r[1] + \
                        ((HX3X4)-t*(SX3X4)+(p-1.0)*(VX3X4))*r[1]*r[2] + \
                        ((HX3X5)-t*(SX3X5)+(p-1.0)*(VX3X5))*r[1]*r[3] + \
                        ((HX3X6)-t*(SX3X6)+(p-1.0)*(VX3X6))*r[1]*r[4] + \
                        ((HX3X7)-t*(SX3X7)+(p-1.0)*(VX3X7))*r[1]*r[5] + \
                        ((HX3S1)-t*(SX3S1)+(p-1.0)*(VX3S1))*r[1]*s[0] + \
                        ((HX3S2)-t*(SX3S2)+(p-1.0)*(VX3S2))*r[1]*s[1] + \
                        ((HX4X4)-t*(SX4X4)+(p-1.0)*(VX4X4))*r[2]*r[2] + \
                        ((HX4X5)-t*(SX4X5)+(p-1.0)*(VX4X5))*r[2]*r[3] + \
                        ((HX4X6)-t*(SX4X6)+(p-1.0)*(VX4X6))*r[2]*r[4] + \
                        ((HX4X7)-t*(SX4X7)+(p-1.0)*(VX4X7))*r[2]*r[5] + \
                        ((HX4S1)-t*(SX4S1)+(p-1.0)*(VX4S1))*r[2]*s[0] + \
                        ((HX4S2)-t*(SX4S2)+(p-1.0)*(VX4S2))*r[2]*s[1] + \
                        ((HX5X5)-t*(SX5X5)+(p-1.0)*(VX5X5))*r[3]*r[3] + \
                        ((HX5X6)-t*(SX5X6)+(p-1.0)*(VX5X6))*r[3]*r[4] + \
                        ((HX5X7)-t*(SX5X7)+(p-1.0)*(VX5X7))*r[3]*r[5] + \
                        ((HX5S1)-t*(SX5S1)+(p-1.0)*(VX5S1))*r[3]*s[0] + \
                        ((HX5S2)-t*(SX5S2)+(p-1.0)*(VX5S2))*r[3]*s[1] + \
                        ((HX6X6)-t*(SX6X6)+(p-1.0)*(VX6X6))*r[4]*r[4] + \
                        ((HX6X7)-t*(SX6X7)+(p-1.0)*(VX6X7))*r[4]*r[5] + \
                        ((HX6S1)-t*(SX6S1)+(p-1.0)*(VX6S1))*r[4]*s[0] + \
                        ((HX6S2)-t*(SX6S2)+(p-1.0)*(VX6S2))*r[4]*s[1] + \
                        ((HX7X7)-t*(SX7X7)+(p-1.0)*(VX7X7))*r[5]*r[5] + \
                        ((HX7S1)-t*(SX7S1)+(p-1.0)*(VX7S1))*r[5]*s[0] + \
                        ((HX7S2)-t*(SX7S2)+(p-1.0)*(VX7S2))*r[5]*s[1] + \
                        ((HS1S1)-t*(SS1S1)+(p-1.0)*(VS1S1))*s[0]*s[0] + \
                        ((HS1S2)-t*(SS1S2)+(p-1.0)*(VS1S2))*s[0]*s[1] + \
                        ((HS2S2)-t*(SS2S2)+(p-1.0)*(VS2S2))*s[1]*s[1] + \
                        ((HX2X2X7)+(p-1.0)*(VX2X2X7))*r[0]*r[0]*r[5] + \
                        ((HX2X2S2)+(p-1.0)*(VX2X2S2))*r[0]*r[0]*s[1] + \
                        ((HX2X3X7)+(p-1.0)*(VX2X3X7))*r[0]*r[1]*r[5] + \
                        ((HX2X3S2)+(p-1.0)*(VX2X3S2))*r[0]*r[1]*s[1] + \
                        ((HX2X4X7)+(p-1.0)*(VX2X4X7))*r[0]*r[2]*r[5] + \
                        ((HX2X4S2)+(p-1.0)*(VX2X4S2))*r[0]*r[2]*s[1] + \
                        ((HX2X5X7)+(p-1.0)*(VX2X5X7))*r[0]*r[3]*r[5] + \
                        ((HX2X5S2)+(p-1.0)*(VX2X5S2))*r[0]*r[3]*s[1] + \
                        ((HX2X6X7)+(p-1.0)*(VX2X6X7))*r[0]*r[4]*r[5] + \
                        ((HX2X6S2)+(p-1.0)*(VX2X6S2))*r[0]*r[4]*s[1] + \
                        ((HX2X7X7)-t*(SX2X7X7)+(p-1.0)*(VX2X7X7))*r[0]*r[5]*r[5] + \
                        ((HX2X7S2)+(p-1.0)*(VX2X7S2))*r[0]*r[5]*s[1] + \
                        ((HX2S2S2)+(p-1.0)*(VX2S2S2))*r[0]*s[1]*s[1] + \
                        ((HX3X3X7)+(p-1.0)*(VX3X3X7))*r[1]*r[1]*r[5] + \
                        ((HX3X3S2)+(p-1.0)*(VX3X3S2))*r[1]*r[1]*s[1] + \
                        ((HX3X4X7)+(p-1.0)*(VX3X4X7))*r[1]*r[2]*r[5] + \
                        ((HX3X4S2)+(p-1.0)*(VX3X4S2))*r[1]*r[2]*s[1] + \
                        ((HX3X5X7)+(p-1.0)*(VX3X5X7))*r[1]*r[3]*r[5] + \
                        ((HX3X5S2)+(p-1.0)*(VX3X5S2))*r[1]*r[3]*s[1] + \
                        ((HX3X6X7)+(p-1.0)*(VX3X6X7))*r[1]*r[4]*r[5] + \
                        ((HX3X6S2)+(p-1.0)*(VX3X6S2))*r[1]*r[4]*s[1] + \
                        ((HX3X7X7)-t*(SX3X7X7)+(p-1.0)*(VX3X7X7))*r[1]*r[5]*r[5] + \
                        ((HX3X7S2)-t*(SX3X7S2)+(p-1.0)*(VX3X7S2))*r[1]*r[5]*s[1] + \
                        ((HX3S2S2)+(p-1.0)*(VX3S2S2))*r[1]*s[1]*s[1] + \
                        ((HX4X4X7)+(p-1.0)*(VX4X4X7))*r[2]*r[2]*r[5] + \
                        ((HX4X4S2)+(p-1.0)*(VX4X4S2))*r[2]*r[2]*s[1] + \
                        ((HX4X5X7)+(p-1.0)*(VX4X5X7))*r[2]*r[3]*r[5] + \
                        ((HX4X5S2)+(p-1.0)*(VX4X5S2))*r[2]*r[3]*s[1] + \
                        ((HX4X6X7)+(p-1.0)*(VX4X6X7))*r[2]*r[4]*r[5] + \
                        ((HX4X6S2)+(p-1.0)*(VX4X6S2))*r[2]*r[4]*s[1] + \
                        ((HX4X7X7)-t*(SX4X7X7)+(p-1.0)*(VX4X7X7))*r[2]*r[5]*r[5] + \
                        ((HX4X7S2)-t*(SX4X7S2)+(p-1.0)*(VX4X7S2))*r[2]*r[5]*s[1] + \
                        ((HX4S2S2)+(p-1.0)*(VX4S2S2))*r[2]*s[1]*s[1] + \
                        ((HX5X5X7)+(p-1.0)*(VX5X5X7))*r[3]*r[3]*r[5] + \
                        ((HX5X5S2)+(p-1.0)*(VX5X5S2))*r[3]*r[3]*s[1] + \
                        ((HX5X6X7)+(p-1.0)*(VX5X6X7))*r[3]*r[4]*r[5] + \
                        ((HX5X6S2)+(p-1.0)*(VX5X6S2))*r[3]*r[4]*s[1] + \
                        ((HX5X7X7)-t*(SX5X7X7)+(p-1.0)*(VX5X7X7))*r[3]*r[5]*r[5] + \
                        ((HX5X7S2)-t*(SX5X7S2)+(p-1.0)*(VX5X7S2))*r[3]*r[5]*s[1] + \
                        ((HX5S2S2)+(p-1.0)*(VX5S2S2))*r[3]*s[1]*s[1] + \
                        ((HX6X6X7)+(p-1.0)*(VX6X6X7))*r[4]*r[4]*r[5] + \
                        ((HX6X6S2)+(p-1.0)*(VX6X6S2))*r[4]*r[4]*s[1] + \
                        ((HX6X7X7)-t*(SX6X7X7)+(p-1.0)*(VX6X7X7))*r[4]*r[5]*r[5] + \
                        ((HX6X7S2)-t*(SX6X7S2)+(p-1.0)*(VX6X7S2))*r[4]*r[5]*s[1] + \
                        ((HX6S2S2)+(p-1.0)*(VX6S2S2))*r[4]*s[1]*s[1] + \
                        ((HX7X7X7)-t*(SX7X7X7)+(p-1.0)*(VX7X7X7))*r[5]*r[5]*r[5] + \
                        ((HX7X7S2)-t*(SX7X7S2)+(p-1.0)*(VX7X7S2))*r[5]*r[5]*s[1] + \
                        ((HX7S2S2)+(p-1.0)*(VX7S2S2))*r[5]*s[1]*s[1]

/*----------------------------------------------------------------------------*/

#define DGDR0 R*t*(log(xfe2m1) - log(xmg2m1)) + \
                            ((HX2)-t*(SX2)+(p-1.0)*(VX2)) + \
                            ((HX2X2)-t*(SX2X2)+(p-1.0)*(VX2X2))*r[0]*2.0 + \
                            ((HX2X3)-t*(SX2X3)+(p-1.0)*(VX2X3))*r[1] + \
                            ((HX2X4)-t*(SX2X4)+(p-1.0)*(VX2X4))*r[2] + \
                            ((HX2X5)-t*(SX2X5)+(p-1.0)*(VX2X5))*r[3] + \
                            ((HX2X6)-t*(SX2X6)+(p-1.0)*(VX2X6))*r[4] + \
                            ((HX2X7)-t*(SX2X7)+(p-1.0)*(VX2X7))*r[5] + \
                            ((HX2S1)-t*(SX2S1)+(p-1.0)*(VX2S1))*s[0] + \
                            ((HX2S2)-t*(SX2S2)+(p-1.0)*(VX2S2))*s[1] + \
                            ((HX2X2X7)+(p-1.0)*(VX2X2X7))*r[0]*r[5]*2.0 + \
                            ((HX2X2S2)+(p-1.0)*(VX2X2S2))*r[0]*s[1]*2.0 + \
                            ((HX2X3X7)+(p-1.0)*(VX2X3X7))*r[1]*r[5] + \
                            ((HX2X3S2)+(p-1.0)*(VX2X3S2))*r[1]*s[1] + \
                            ((HX2X4X7)+(p-1.0)*(VX2X4X7))*r[2]*r[5] + \
                            ((HX2X4S2)+(p-1.0)*(VX2X4S2))*r[2]*s[1] + \
                            ((HX2X5X7)+(p-1.0)*(VX2X5X7))*r[3]*r[5] + \
                            ((HX2X5S2)+(p-1.0)*(VX2X5S2))*r[3]*s[1] + \
                            ((HX2X6X7)+(p-1.0)*(VX2X6X7))*r[4]*r[5] + \
                            ((HX2X6S2)+(p-1.0)*(VX2X6S2))*r[4]*s[1] + \
                            ((HX2X7X7)-t*(SX2X7X7)+(p-1.0)*(VX2X7X7))*r[5]*r[5] + \
                            ((HX2X7S2)+(p-1.0)*(VX2X7S2))*r[5]*s[1] + \
                            ((HX2S2S2)+(p-1.0)*(VX2S2S2))*s[1]*s[1]
#define DGDR1 R*t*(0.5*log(xti4m1) - 0.5*log(xmg2m1) + 0.5*log(1.0-xti4m1) \
                                - log(xmg2m1+xfe2m1-xti4m1) + 0.5*log(xmg2m1+xfe2m1) \
                                + 0.5*log(1.0-xmg2m1-xfe2m1-xna1m2) - log(1.0-xsi4tet) \
                                - 0.5*log(1.0-xmg2m1-xfe2m1) + log(xal3tet)) + \
                            ((HX3)-t*(SX3)+(p-1.0)*(VX3)) + \
                            ((HX2X3)-t*(SX2X3)+(p-1.0)*(VX2X3))*r[0] + \
                            ((HX3X3)-t*(SX3X3)+(p-1.0)*(VX3X3))*r[1]*2.0 + \
                            ((HX3X4)-t*(SX3X4)+(p-1.0)*(VX3X4))*r[2] + \
                            ((HX3X5)-t*(SX3X5)+(p-1.0)*(VX3X5))*r[3] + \
                            ((HX3X6)-t*(SX3X6)+(p-1.0)*(VX3X6))*r[4] + \
                            ((HX3X7)-t*(SX3X7)+(p-1.0)*(VX3X7))*r[5] + \
                            ((HX3S1)-t*(SX3S1)+(p-1.0)*(VX3S1))*s[0] + \
                            ((HX3S2)-t*(SX3S2)+(p-1.0)*(VX3S2))*s[1] + \
                            ((HX2X3X7)+(p-1.0)*(VX2X3X7))*r[0]*r[5] + \
                            ((HX2X3S2)+(p-1.0)*(VX2X3S2))*r[0]*s[1] + \
                            ((HX3X3X7)+(p-1.0)*(VX3X3X7))*r[1]*r[5]*2.0 + \
                            ((HX3X3S2)+(p-1.0)*(VX3X3S2))*r[1]*s[1]*2.0 + \
                            ((HX3X4X7)+(p-1.0)*(VX3X4X7))*r[2]*r[5] + \
                            ((HX3X4S2)+(p-1.0)*(VX3X4S2))*r[2]*s[1] + \
                            ((HX3X5X7)+(p-1.0)*(VX3X5X7))*r[3]*r[5] + \
                            ((HX3X5S2)+(p-1.0)*(VX3X5S2))*r[3]*s[1] + \
                            ((HX3X6X7)+(p-1.0)*(VX3X6X7))*r[4]*r[5] + \
                            ((HX3X6S2)+(p-1.0)*(VX3X6S2))*r[4]*s[1] + \
                            ((HX3X7X7)-t*(SX3X7X7)+(p-1.0)*(VX3X7X7))*r[5]*r[5] + \
                            ((HX3X7S2)-t*(SX3X7S2)+(p-1.0)*(VX3X7S2))*r[5]*s[1] + \
                            ((HX3S2S2)+(p-1.0)*(VX3S2S2))*s[1]*s[1]
#define DGDR2 R*t*(0.5*log(xti4m1) - 0.5*log(xmg2m1) + 0.5*log(1.0-xti4m1) \
                                - log(xmg2m1+xfe2m1-xti4m1) + 0.5*log(xmg2m1+xfe2m1) \
                                + 0.5*log(1.0-xmg2m1-xfe2m1-xna1m2) - log(1.0-xsi4tet) \
                                - 0.5*log(1.0-xmg2m1-xfe2m1) + log(xfe3tet)) + \
                            ((HX4)-t*(SX4)+(p-1.0)*(VX4)) + \
                            ((HX2X4)-t*(SX2X4)+(p-1.0)*(VX2X4))*r[0] + \
                            ((HX3X4)-t*(SX3X4)+(p-1.0)*(VX3X4))*r[1] + \
                            ((HX4X4)-t*(SX4X4)+(p-1.0)*(VX4X4))*r[2]*2.0 + \
                            ((HX4X5)-t*(SX4X5)+(p-1.0)*(VX4X5))*r[3] + \
                            ((HX4X6)-t*(SX4X6)+(p-1.0)*(VX4X6))*r[4] + \
                            ((HX4X7)-t*(SX4X7)+(p-1.0)*(VX4X7))*r[5] + \
                            ((HX4S1)-t*(SX4S1)+(p-1.0)*(VX4S1))*s[0] + \
                            ((HX4S2)-t*(SX4S2)+(p-1.0)*(VX4S2))*s[1] + \
                            ((HX2X4X7)+(p-1.0)*(VX2X4X7))*r[0]*r[5] + \
                            ((HX2X4S2)+(p-1.0)*(VX2X4S2))*r[0]*s[1] + \
                            ((HX3X4X7)+(p-1.0)*(VX3X4X7))*r[1]*r[5] + \
                            ((HX3X4S2)+(p-1.0)*(VX3X4S2))*r[1]*s[1] + \
                            ((HX4X4X7)+(p-1.0)*(VX4X4X7))*r[2]*r[5]*2.0 + \
                            ((HX4X4S2)+(p-1.0)*(VX4X4S2))*r[2]*s[1]*2.0 + \
                            ((HX4X5X7)+(p-1.0)*(VX4X5X7))*r[3]*r[5] + \
                            ((HX4X5S2)+(p-1.0)*(VX4X5S2))*r[3]*s[1] + \
                            ((HX4X6X7)+(p-1.0)*(VX4X6X7))*r[4]*r[5] + \
                            ((HX4X6S2)+(p-1.0)*(VX4X6S2))*r[4]*s[1] + \
                            ((HX4X7X7)-t*(SX4X7X7)+(p-1.0)*(VX4X7X7))*r[5]*r[5] + \
                            ((HX4X7S2)-t*(SX4X7S2)+(p-1.0)*(VX4X7S2))*r[5]*s[1] + \
                            ((HX4S2S2)+(p-1.0)*(VX4S2S2))*s[1]*s[1]
#define DGDR3 R*t*(0.5*log(xal3m1) + 0.5*log(xfe3m1) - log(xmg2m1) \
                                - log(xmg2m1+xfe2m1-xti4m1) + log(xmg2m1+xfe2m1) \
                                - log(1.0-xsi4tet) + 0.5*log(xal3tet) + 0.5*log(xfe3tet) \
                                + log(1.0-xmg2m1-xfe2m1-xna1m2) - log(1.0-xmg2m1-xfe2m1)) + \
                            ((HX5)-t*(SX5)+(p-1.0)*(VX5)) + \
                            ((HX2X5)-t*(SX2X5)+(p-1.0)*(VX2X5))*r[0] + \
                            ((HX3X5)-t*(SX3X5)+(p-1.0)*(VX3X5))*r[1] + \
                            ((HX4X5)-t*(SX4X5)+(p-1.0)*(VX4X5))*r[2] + \
                            ((HX5X5)-t*(SX5X5)+(p-1.0)*(VX5X5))*r[3]*2.0 + \
                            ((HX5X6)-t*(SX5X6)+(p-1.0)*(VX5X6))*r[4] + \
                            ((HX5X7)-t*(SX5X7)+(p-1.0)*(VX5X7))*r[5] + \
                            ((HX5S1)-t*(SX5S1)+(p-1.0)*(VX5S1))*s[0] + \
                            ((HX5S2)-t*(SX5S2)+(p-1.0)*(VX5S2))*s[1] + \
                            ((HX2X5X7)+(p-1.0)*(VX2X5X7))*r[0]*r[5] + \
                            ((HX2X5S2)+(p-1.0)*(VX2X5S2))*r[0]*s[1] + \
                            ((HX3X5X7)+(p-1.0)*(VX3X5X7))*r[1]*r[5] + \
                            ((HX3X5S2)+(p-1.0)*(VX3X5S2))*r[1]*s[1] + \
                            ((HX4X5X7)+(p-1.0)*(VX4X5X7))*r[2]*r[5] + \
                            ((HX4X5S2)+(p-1.0)*(VX4X5S2))*r[2]*s[1] + \
                            ((HX5X5X7)+(p-1.0)*(VX5X5X7))*r[3]*r[5]*2.0 + \
                            ((HX5X5S2)+(p-1.0)*(VX5X5S2))*r[3]*s[1]*2.0 + \
                            ((HX5X6X7)+(p-1.0)*(VX5X6X7))*r[4]*r[5] + \
                            ((HX5X6S2)+(p-1.0)*(VX5X6S2))*r[4]*s[1] + \
                            ((HX5X7X7)-t*(SX5X7X7)+(p-1.0)*(VX5X7X7))*r[5]*r[5] + \
                            ((HX5X7S2)-t*(SX5X7S2)+(p-1.0)*(VX5X7S2))*r[5]*s[1] + \
                            ((HX5S2S2)+(p-1.0)*(VX5S2S2))*s[1]*s[1]
#define DGDR4 R*t*(0.25*log(xal3m1) + 0.25*log(xfe3m1) - 0.5*log(xmg2m1) \
                                - 0.5*log(xmg2m1+xfe2m1-xti4m1) + 0.5*log(xmg2m1+xfe2m1) \
                                + 0.5*log(1.0-xsi4tet) - 0.25*log(xal3tet) \
                                - 0.25*log(xfe3tet) - log(xca2m2) + log(xna1m2) \
                                - 0.5*log(1.0-xmg2m1-xfe2m1-xna1m2) \
                                - 0.5*log(1.0-xmg2m1-xfe2m1) + log(1.0-xna1m2)) + \
                            ((HX6)-t*(SX6)+(p-1.0)*(VX6)) + \
                            ((HX2X6)-t*(SX2X6)+(p-1.0)*(VX2X6))*r[0] + \
                            ((HX3X6)-t*(SX3X6)+(p-1.0)*(VX3X6))*r[1] + \
                            ((HX4X6)-t*(SX4X6)+(p-1.0)*(VX4X6))*r[2] + \
                            ((HX5X6)-t*(SX5X6)+(p-1.0)*(VX5X6))*r[3] + \
                            ((HX6X6)-t*(SX6X6)+(p-1.0)*(VX6X6))*r[4]*2.0 + \
                            ((HX6X7)-t*(SX6X7)+(p-1.0)*(VX6X7))*r[5] + \
                            ((HX6S1)-t*(SX6S1)+(p-1.0)*(VX6S1))*s[0] + \
                            ((HX6S2)-t*(SX6S2)+(p-1.0)*(VX6S2))*s[1] + \
                            ((HX2X6X7)+(p-1.0)*(VX2X6X7))*r[0]*r[5] + \
                            ((HX2X6S2)+(p-1.0)*(VX2X6S2))*r[0]*s[1] + \
                            ((HX3X6X7)+(p-1.0)*(VX3X6X7))*r[1]*r[5] + \
                            ((HX3X6S2)+(p-1.0)*(VX3X6S2))*r[1]*s[1] + \
                            ((HX4X6X7)+(p-1.0)*(VX4X6X7))*r[2]*r[5] + \
                            ((HX4X6S2)+(p-1.0)*(VX4X6S2))*r[2]*s[1] + \
                            ((HX5X6X7)+(p-1.0)*(VX5X6X7))*r[3]*r[5] + \
                            ((HX5X6S2)+(p-1.0)*(VX5X6S2))*r[3]*s[1] + \
                            ((HX6X6X7)+(p-1.0)*(VX6X6X7))*r[4]*r[5]*2.0 + \
                            ((HX6X6S2)+(p-1.0)*(VX6X6S2))*r[4]*s[1]*2.0 + \
                            ((HX6X7X7)-t*(SX6X7X7)+(p-1.0)*(VX6X7X7))*r[5]*r[5] + \
                            ((HX6X7S2)-t*(SX6X7S2)+(p-1.0)*(VX6X7S2))*r[5]*s[1] + \
                            ((HX6S2S2)+(p-1.0)*(VX6S2S2))*s[1]*s[1]
#define DGDR5 0.5*R*t*(log(xmg2m1) - 2.0*log(xca2m2) \
                                + log(xmg2m2) + log(xfe2m2/xfe2m1)) + \
                            ((HX7)-t*(SX7)+(p-1.0)*(VX7)) + \
                            ((HX2X7)-t*(SX2X7)+(p-1.0)*(VX2X7))*r[0] + \
                            ((HX3X7)-t*(SX3X7)+(p-1.0)*(VX3X7))*r[1] + \
                            ((HX4X7)-t*(SX4X7)+(p-1.0)*(VX4X7))*r[2] + \
                            ((HX5X7)-t*(SX5X7)+(p-1.0)*(VX5X7))*r[3] + \
                            ((HX6X7)-t*(SX6X7)+(p-1.0)*(VX6X7))*r[4] + \
                            ((HX7X7)-t*(SX7X7)+(p-1.0)*(VX7X7))*r[5]*2.0 + \
                            ((HX7S1)-t*(SX7S1)+(p-1.0)*(VX7S1))*s[0] + \
                            ((HX7S2)-t*(SX7S2)+(p-1.0)*(VX7S2))*s[1] + \
                            ((HX2X2X7)+(p-1.0)*(VX2X2X7))*r[0]*r[0] + \
                            ((HX2X3X7)+(p-1.0)*(VX2X3X7))*r[0]*r[1] + \
                            ((HX2X4X7)+(p-1.0)*(VX2X4X7))*r[0]*r[2] + \
                            ((HX2X5X7)+(p-1.0)*(VX2X5X7))*r[0]*r[3] + \
                            ((HX2X6X7)+(p-1.0)*(VX2X6X7))*r[0]*r[4] + \
                            ((HX2X7X7)-t*(SX2X7X7)+(p-1.0)*(VX2X7X7))*r[0]*r[5]*2.0 + \
                            ((HX2X7S2)+(p-1.0)*(VX2X7S2))*r[0]*s[1] + \
                            ((HX3X3X7)+(p-1.0)*(VX3X3X7))*r[1]*r[1] + \
                            ((HX3X4X7)+(p-1.0)*(VX3X4X7))*r[1]*r[2] + \
                            ((HX3X5X7)+(p-1.0)*(VX3X5X7))*r[1]*r[3] + \
                            ((HX3X6X7)+(p-1.0)*(VX3X6X7))*r[1]*r[4] + \
                            ((HX3X7X7)-t*(SX3X7X7)+(p-1.0)*(VX3X7X7))*r[1]*r[5]*2.0 + \
                            ((HX3X7S2)-t*(SX3X7S2)+(p-1.0)*(VX3X7S2))*r[1]*s[1] + \
                            ((HX4X4X7)+(p-1.0)*(VX4X4X7))*r[2]*r[2] + \
                            ((HX4X5X7)+(p-1.0)*(VX4X5X7))*r[2]*r[3] + \
                            ((HX4X6X7)+(p-1.0)*(VX4X6X7))*r[2]*r[4] + \
                            ((HX4X7X7)-t*(SX4X7X7)+(p-1.0)*(VX4X7X7))*r[2]*r[5]*2.0 + \
                            ((HX4X7S2)-t*(SX4X7S2)+(p-1.0)*(VX4X7S2))*r[2]*s[1] + \
                            ((HX5X5X7)+(p-1.0)*(VX5X5X7))*r[3]*r[3] + \
                            ((HX5X6X7)+(p-1.0)*(VX5X6X7))*r[3]*r[4] + \
                            ((HX5X7X7)-t*(SX5X7X7)+(p-1.0)*(VX5X7X7))*r[3]*r[5]*2.0 + \
                            ((HX5X7S2)-t*(SX5X7S2)+(p-1.0)*(VX5X7S2))*r[3]*s[1] + \
                            ((HX6X6X7)+(p-1.0)*(VX6X6X7))*r[4]*r[4] + \
                            ((HX6X7X7)-t*(SX6X7X7)+(p-1.0)*(VX6X7X7))*r[4]*r[5]*2.0 + \
                            ((HX6X7S2)-t*(SX6X7S2)+(p-1.0)*(VX6X7S2))*r[4]*s[1] + \
                            ((HX7X7X7)-t*(SX7X7X7)+(p-1.0)*(VX7X7X7))*r[5]*r[5]*3.0 + \
                            ((HX7X7S2)-t*(SX7X7S2)+(p-1.0)*(VX7X7S2))*r[5]*s[1]*2.0 + \
                            ((HX7S2S2)+(p-1.0)*(VX7S2S2))*s[1]*s[1]
#define DGDS0 0.5*R*t*(log(xfe3m1)-log(xal3m1)+log(xal3tet)-log(xfe3tet)) + \
                            ((HS1)-t*(SS1)+(p-1.0)*(VS1)) + \
                            ((HX2S1)-t*(SX2S1)+(p-1.0)*(VX2S1))*r[0] + \
                            ((HX3S1)-t*(SX3S1)+(p-1.0)*(VX3S1))*r[1] + \
                            ((HX4S1)-t*(SX4S1)+(p-1.0)*(VX4S1))*r[2] + \
                            ((HX5S1)-t*(SX5S1)+(p-1.0)*(VX5S1))*r[3] + \
                            ((HX6S1)-t*(SX6S1)+(p-1.0)*(VX6S1))*r[4] + \
                            ((HX7S1)-t*(SX7S1)+(p-1.0)*(VX7S1))*r[5] + \
                            ((HS1S1)-t*(SS1S1)+(p-1.0)*(VS1S1))*s[0]*2.0 + \
                            ((HS1S2)-t*(SS1S2)+(p-1.0)*(VS1S2))*s[1]
#define DGDS1 0.5*R*t*(log(xmg2m1)-log(xfe2m1)+log(xfe2m2)-log(xmg2m2)) + \
                            ((HS2)-t*(SS2)+(p-1.0)*(VS2)) + \
                            ((HX2S2)-t*(SX2S2)+(p-1.0)*(VX2S2))*r[0] + \
                            ((HX3S2)-t*(SX3S2)+(p-1.0)*(VX3S2))*r[1] + \
                            ((HX4S2)-t*(SX4S2)+(p-1.0)*(VX4S2))*r[2] + \
                            ((HX5S2)-t*(SX5S2)+(p-1.0)*(VX5S2))*r[3] + \
                            ((HX6S2)-t*(SX6S2)+(p-1.0)*(VX6S2))*r[4] + \
                            ((HX7S2)-t*(SX7S2)+(p-1.0)*(VX7S2))*r[5] + \
                            ((HS1S2)-t*(SS1S2)+(p-1.0)*(VS1S2))*s[0] + \
                            ((HS2S2)-t*(SS2S2)+(p-1.0)*(VS2S2))*s[1]*2.0 + \
                            ((HX2X2S2)+(p-1.0)*(VX2X2S2))*r[0]*r[0] + \
                            ((HX2X3S2)+(p-1.0)*(VX2X3S2))*r[0]*r[1] + \
                            ((HX2X4S2)+(p-1.0)*(VX2X4S2))*r[0]*r[2] + \
                            ((HX2X5S2)+(p-1.0)*(VX2X5S2))*r[0]*r[3] + \
                            ((HX2X6S2)+(p-1.0)*(VX2X6S2))*r[0]*r[4] + \
                            ((HX2X7S2)+(p-1.0)*(VX2X7S2))*r[0]*r[5] + \
                            ((HX2S2S2)+(p-1.0)*(VX2S2S2))*r[0]*s[1]*2.0 + \
                            ((HX3X3S2)+(p-1.0)*(VX3X3S2))*r[1]*r[1] + \
                            ((HX3X4S2)+(p-1.0)*(VX3X4S2))*r[1]*r[2] + \
                            ((HX3X5S2)+(p-1.0)*(VX3X5S2))*r[1]*r[3] + \
                            ((HX3X6S2)+(p-1.0)*(VX3X6S2))*r[1]*r[4] + \
                            ((HX3X7S2)-t*(SX3X7S2)+(p-1.0)*(VX3X7S2))*r[1]*r[5] + \
                            ((HX3S2S2)+(p-1.0)*(VX3S2S2))*r[1]*s[1]*2.0 + \
                            ((HX4X4S2)+(p-1.0)*(VX4X4S2))*r[2]*r[2] + \
                            ((HX4X5S2)+(p-1.0)*(VX4X5S2))*r[2]*r[3] + \
                            ((HX4X6S2)+(p-1.0)*(VX4X6S2))*r[2]*r[4] + \
                            ((HX4X7S2)-t*(SX4X7S2)+(p-1.0)*(VX4X7S2))*r[2]*r[5] + \
                            ((HX4S2S2)+(p-1.0)*(VX4S2S2))*r[2]*s[1]*2.0 + \
                            ((HX5X5S2)+(p-1.0)*(VX5X5S2))*r[3]*r[3] + \
                            ((HX5X6S2)+(p-1.0)*(VX5X6S2))*r[3]*r[4] + \
                            ((HX5X7S2)-t*(SX5X7S2)+(p-1.0)*(VX5X7S2))*r[3]*r[5] + \
                            ((HX5S2S2)+(p-1.0)*(VX5S2S2))*r[3]*s[1]*2.0 + \
                            ((HX6X6S2)+(p-1.0)*(VX6X6S2))*r[4]*r[4] + \
                            ((HX6X7S2)-t*(SX6X7S2)+(p-1.0)*(VX6X7S2))*r[4]*r[5] + \
                            ((HX6S2S2)+(p-1.0)*(VX6S2S2))*r[4]*s[1]*2.0 + \
                            ((HX7X7S2)-t*(SX7X7S2)+(p-1.0)*(VX7X7S2))*r[5]*r[5] + \
                            ((HX7S2S2)+(p-1.0)*(VX7S2S2))*r[5]*s[1]*2.0
#define DGDT  -(S)
#define DGDP  (V)

/*----------------------------------------------------------------------------*/

#define D2GDR0R0 R*t*(1.0/xfe2m1 + 1.0/xmg2m1) + \
                 ((HX2X2)-t*(SX2X2)+(p-1.0)*(VX2X2))*2.0 + \
                 ((HX2X2X7)+(p-1.0)*(VX2X2X7))*r[5]*2.0 + \
                 ((HX2X2S2)+(p-1.0)*(VX2X2S2))*s[1]*2.0
#define D2GDR0R1 R*t*0.5/xmg2m1 + ((HX2X3)-t*(SX2X3)+(p-1.0)*(VX2X3)) + \
                 ((HX2X3X7)+(p-1.0)*(VX2X3X7))*r[5] + \
                 ((HX2X3S2)+(p-1.0)*(VX2X3S2))*s[1]
#define D2GDR0R2 R*t*0.5/xmg2m1 + ((HX2X4)-t*(SX2X4)+(p-1.0)*(VX2X4)) + \
                 ((HX2X4X7)+(p-1.0)*(VX2X4X7))*r[5] + \
                 ((HX2X4S2)+(p-1.0)*(VX2X4S2))*s[1]
#define D2GDR0R3 R*t/xmg2m1 + ((HX2X5)-t*(SX2X5)+(p-1.0)*(VX2X5)) + \
                 ((HX2X5X7)+(p-1.0)*(VX2X5X7))*r[5] + \
                 ((HX2X5S2)+(p-1.0)*(VX2X5S2))*s[1]
#define D2GDR0R4 R*t*0.5/xmg2m1 + ((HX2X6)-t*(SX2X6)+(p-1.0)*(VX2X6)) + \
                 ((HX2X6X7)+(p-1.0)*(VX2X6X7))*r[5] + \
                 ((HX2X6S2)+(p-1.0)*(VX2X6S2))*s[1]
#define D2GDR0R5 -R*t*0.5*(1.0/xfe2m1 + 1.0/xmg2m1) + \
                 ((HX2X7)-t*(SX2X7)+(p-1.0)*(VX2X7)) + \
                 ((HX2X2X7)+(p-1.0)*(VX2X2X7))*r[0]*2.0 + \
                 ((HX2X3X7)+(p-1.0)*(VX2X3X7))*r[1] + \
                 ((HX2X4X7)+(p-1.0)*(VX2X4X7))*r[2] + \
                 ((HX2X5X7)+(p-1.0)*(VX2X5X7))*r[3] + \
                 ((HX2X6X7)+(p-1.0)*(VX2X6X7))*r[4] + \
                 ((HX2X7X7)-t*(SX2X7X7)+(p-1.0)*(VX2X7X7))*r[5]*2.0 + \
                 ((HX2X7S2)+(p-1.0)*(VX2X7S2))*s[1]
#define D2GDR0S0 ((HX2S1)-t*(SX2S1)+(p-1.0)*(VX2S1))
#define D2GDR0S1 -R*t*0.5*(1.0/xfe2m1 + 1.0/xmg2m1) + \
                 ((HX2S2)-t*(SX2S2)+(p-1.0)*(VX2S2)) + \
                 ((HX2X2S2)+(p-1.0)*(VX2X2S2))*r[0]*2.0 + \
                 ((HX2X3S2)+(p-1.0)*(VX2X3S2))*r[1] + \
                 ((HX2X4S2)+(p-1.0)*(VX2X4S2))*r[2] + \
                 ((HX2X5S2)+(p-1.0)*(VX2X5S2))*r[3] + \
                 ((HX2X6S2)+(p-1.0)*(VX2X6S2))*r[4] + \
                 ((HX2X7S2)+(p-1.0)*(VX2X7S2))*r[5] + \
                 ((HX2S2S2)+(p-1.0)*(VX2S2S2))*s[1]*2.0
#define D2GDR0DT R*(log(xfe2m1) - log(xmg2m1)) - (SX2) - (SX2X2)*r[0]*2.0 - \
                 (SX2X3)*r[1] - (SX2X4)*r[2] - (SX2X5)*r[3] - (SX2X6)*r[4] - \
                 (SX2X7)*r[5] - (SX2S1)*s[0] - (SX2S2)*s[1] - \
                 (SX2X7X7)*r[5]*r[5]
#define D2GDR0DP (VX2) + (VX2X2)*r[0]*2.0 + (VX2X3)*r[1] + (VX2X4)*r[2] + \
                 (VX2X5)*r[3] + (VX2X6)*r[4] + (VX2X7)*r[5] + (VX2S1)*s[0] + \
                 (VX2S2)*s[1] + (VX2X2X7)*r[0]*r[5]*2.0 + \
                 (VX2X2S2)*r[0]*s[1]*2.0 + (VX2X3X7)*r[1]*r[5] + \
                 (VX2X3S2)*r[1]*s[1] + (VX2X4X7)*r[2]*r[5] + \
                 (VX2X4S2)*r[2]*s[1] + (VX2X5X7)*r[3]*r[5] + \
                 (VX2X5S2)*r[3]*s[1] + (VX2X6X7)*r[4]*r[5] + \
                 (VX2X6S2)*r[4]*s[1] + (VX2X7X7)*r[5]*r[5] + \
                 (VX2X7S2)*r[5]*s[1] + (VX2S2S2)*s[1]*s[1]

#define D2GDR1R1 R*t*(0.25/xti4m1 + 0.25/xmg2m1 - 0.25/(1.0-xti4m1) \
                 + 1.0/(xmg2m1+xfe2m1-xti4m1) - 0.25/(xmg2m1+xfe2m1) \
                 + 0.25/(1.0-xmg2m1-xfe2m1-xna1m2) - 0.5/(1.0-xsi4tet) \
                 - 0.25/(1.0-xmg2m1-xfe2m1) + 0.5/xal3tet) + \
                 ((HX3X3)-t*(SX3X3)+(p-1.0)*(VX3X3))*2.0 + \
                 ((HX3X3X7)+(p-1.0)*(VX3X3X7))*r[5]*2.0 + \
                 ((HX3X3S2)+(p-1.0)*(VX3X3S2))*s[1]*2.0
#define D2GDR1R2 R*t*(0.25/xti4m1 + 0.25/xmg2m1 - 0.25/(1.0-xti4m1) \
                 + 1.0/(xmg2m1+xfe2m1-xti4m1) - 0.25/(xmg2m1+xfe2m1) \
                 + 0.25/(1.0-xmg2m1-xfe2m1-xna1m2) - 0.5/(1.0-xsi4tet) \
                 - 0.25/(1.0-xmg2m1-xfe2m1)) + \
                 ((HX3X4)-t*(SX3X4)+(p-1.0)*(VX3X4)) + \
                 ((HX3X4X7)+(p-1.0)*(VX3X4X7))*r[5] + \
                 ((HX3X4S2)+(p-1.0)*(VX3X4S2))*s[1]
#define D2GDR1R3 R*t*(0.5/xmg2m1 + 1.0/(xmg2m1+xfe2m1-xti4m1) \
                 - 0.5/(xmg2m1+xfe2m1) + 0.5/(1.0-xmg2m1-xfe2m1-xna1m2) \
                 - 0.5/(1.0-xsi4tet) \
                 - 0.5/(1.0-xmg2m1-xfe2m1) + 0.25/xal3tet) + \
                 ((HX3X5)-t*(SX3X5)+(p-1.0)*(VX3X5)) + \
                 ((HX3X5X7)+(p-1.0)*(VX3X5X7))*r[5] + \
                 ((HX3X5S2)+(p-1.0)*(VX3X5S2))*s[1]
#define D2GDR1R4 R*t*(0.25/xmg2m1 + 0.5/(xmg2m1+xfe2m1-xti4m1) \
                 - 0.25/(xmg2m1+xfe2m1) - 0.25/(1.0-xmg2m1-xfe2m1-xna1m2) \
                 + 0.25/(1.0-xsi4tet) - 0.25/(1.0-xmg2m1-xfe2m1) \
                 - 0.125/xal3tet) + ((HX3X6)-t*(SX3X6)+(p-1.0)*(VX3X6)) + \
                 ((HX3X6X7)+(p-1.0)*(VX3X6X7))*r[5] + \
                 ((HX3X6S2)+(p-1.0)*(VX3X6S2))*s[1]
#define D2GDR1R5 -R*t*0.25/xmg2m1 + ((HX3X7)-t*(SX3X7)+(p-1.0)*(VX3X7)) + \
                 ((HX2X3X7)+(p-1.0)*(VX2X3X7))*r[0] + \
                 ((HX3X3X7)+(p-1.0)*(VX3X3X7))*r[1]*2.0 + \
                 ((HX3X4X7)+(p-1.0)*(VX3X4X7))*r[2] + \
                 ((HX3X5X7)+(p-1.0)*(VX3X5X7))*r[3] + \
                 ((HX3X6X7)+(p-1.0)*(VX3X6X7))*r[4] + \
                 ((HX3X7X7)-t*(SX3X7X7)+(p-1.0)*(VX3X7X7))*r[5]*2.0 + \
                 ((HX3X7S2)-t*(SX3X7S2)+(p-1.0)*(VX3X7S2))*s[1]
#define D2GDR1S0 R*t*0.25/xal3tet + ((HX3S1)-t*(SX3S1)+(p-1.0)*(VX3S1))
#define D2GDR1S1 -R*t*0.25/xmg2m1 + ((HX3S2)-t*(SX3S2)+(p-1.0)*(VX3S2)) + \
                 ((HX2X3S2)+(p-1.0)*(VX2X3S2))*r[0] + \
                 ((HX3X3S2)+(p-1.0)*(VX3X3S2))*r[1]*2.0 + \
                 ((HX3X4S2)+(p-1.0)*(VX3X4S2))*r[2] + \
                 ((HX3X5S2)+(p-1.0)*(VX3X5S2))*r[3] + \
                 ((HX3X6S2)+(p-1.0)*(VX3X6S2))*r[4] + \
                 ((HX3X7S2)-t*(SX3X7S2)+(p-1.0)*(VX3X7S2))*r[5] + \
                 ((HX3S2S2)+(p-1.0)*(VX3S2S2))*s[1]*2.0
#define D2GDR1DT R*(0.5*log(xti4m1) - 0.5*log(xmg2m1) + 0.5*log(1.0-xti4m1) \
                 - log(xmg2m1+xfe2m1-xti4m1) + 0.5*log(xmg2m1+xfe2m1) \
                 + 0.5*log(1.0-xmg2m1-xfe2m1-xna1m2) - log(1.0-xsi4tet) \
                 - 0.5*log(1.0-xmg2m1-xfe2m1) + log(xal3tet)) - \
                 (SX3) - (SX2X3)*r[0] - (SX3X3)*r[1]*2.0 - (SX3X4)*r[2] - \
                 (SX3X5)*r[3] - (SX3X6)*r[4] - (SX3X7)*r[5] - (SX3S1)*s[0] - \
                 (SX3S2)*s[1] - (SX3X7X7)*r[5]*r[5] - (SX3X7S2)*r[5]*s[1]
#define D2GDR1DP (VX3) + (VX2X3)*r[0] + (VX3X3)*r[1]*2.0 + (VX3X4)*r[2] + \
                 (VX3X5)*r[3] + (VX3X6)*r[4] + (VX3X7)*r[5] + (VX3S1)*s[0] + \
                 (VX3S2)*s[1] + (VX2X3X7)*r[0]*r[5] + \
                 (VX2X3S2)*r[0]*s[1] + (VX3X3X7)*r[1]*r[5]*2.0 + \
                 (VX3X3S2)*r[1]*s[1]*2.0 + (VX3X4X7)*r[2]*r[5] + \
                 (VX3X4S2)*r[2]*s[1] + (VX3X5X7)*r[3]*r[5] + \
                 (VX3X5S2)*r[3]*s[1] + (VX3X6X7)*r[4]*r[5] + \
                 (VX3X6S2)*r[4]*s[1] + (VX3X7X7)*r[5]*r[5] + \
                 (VX3X7S2)*r[5]*s[1] + (VX3S2S2)*s[1]*s[1]

#define D2GDR2R2 R*t*(0.25/xti4m1 + 0.25/xmg2m1 - 0.25/(1.0-xti4m1) \
                 + 1.0/(xmg2m1+xfe2m1-xti4m1) - 0.25/(xmg2m1+xfe2m1) \
                 + 0.25/(1.0-xmg2m1-xfe2m1-xna1m2) - 0.5/(1.0-xsi4tet) \
                 - 0.25/(1.0-xmg2m1-xfe2m1) + 0.5/xfe3tet) + \
                 ((HX4X4)-t*(SX4X4)+(p-1.0)*(VX4X4))*2.0 + \
                 ((HX4X4X7)+(p-1.0)*(VX4X4X7))*r[5]*2.0 + \
                 ((HX4X4S2)+(p-1.0)*(VX4X4S2))*s[1]*2.0
#define D2GDR2R3 R*t*(0.5/xmg2m1 + 1.0/(xmg2m1+xfe2m1-xti4m1) \
                 - 0.5/(xmg2m1+xfe2m1) + 0.5/(1.0-xmg2m1-xfe2m1-xna1m2) \
                 - 0.5/(1.0-xsi4tet) \
                 - 0.5/(1.0-xmg2m1-xfe2m1) + 0.25/xfe3tet) + \
                 ((HX4X5)-t*(SX4X5)+(p-1.0)*(VX4X5)) + \
                 ((HX4X5X7)+(p-1.0)*(VX4X5X7))*r[5] + \
                 ((HX4X5S2)+(p-1.0)*(VX4X5S2))*s[1]
#define D2GDR2R4 R*t*(0.25/xmg2m1 + 0.5/(xmg2m1+xfe2m1-xti4m1) \
                 - 0.25/(xmg2m1+xfe2m1) - 0.25/(1.0-xmg2m1-xfe2m1-xna1m2) \
                 + 0.25/(1.0-xsi4tet) - 0.25/(1.0-xmg2m1-xfe2m1) \
                 - 0.125/xfe3tet) + ((HX4X6)-t*(SX4X6)+(p-1.0)*(VX4X6)) + \
                 ((HX4X6X7)+(p-1.0)*(VX4X6X7))*r[5] + \
                 ((HX4X6S2)+(p-1.0)*(VX4X6S2))*s[1]
#define D2GDR2R5 -R*t*0.25/xmg2m1 + ((HX4X7)-t*(SX4X7)+(p-1.0)*(VX4X7)) + \
                 ((HX2X4X7)+(p-1.0)*(VX2X4X7))*r[0] + \
                 ((HX3X4X7)+(p-1.0)*(VX3X4X7))*r[1] + \
                 ((HX4X4X7)+(p-1.0)*(VX4X4X7))*r[2]*2.0 + \
                 ((HX4X5X7)+(p-1.0)*(VX4X5X7))*r[3] + \
                 ((HX4X6X7)+(p-1.0)*(VX4X6X7))*r[4] + \
                 ((HX4X7X7)-t*(SX4X7X7)+(p-1.0)*(VX4X7X7))*r[5]*2.0 + \
                 ((HX4X7S2)-t*(SX4X7S2)+(p-1.0)*(VX4X7S2))*s[1]
#define D2GDR2S0 -R*t*0.25/xfe3tet + ((HX4S1)-t*(SX4S1)+(p-1.0)*(VX4S1))
#define D2GDR2S1 -R*t*0.25/xmg2m1 + ((HX4S2)-t*(SX4S2)+(p-1.0)*(VX4S2)) + \
                 ((HX2X4S2)+(p-1.0)*(VX2X4S2))*r[0] + \
                 ((HX3X4S2)+(p-1.0)*(VX3X4S2))*r[1] + \
                 ((HX4X4S2)+(p-1.0)*(VX4X4S2))*r[2]*2.0 + \
                 ((HX4X5S2)+(p-1.0)*(VX4X5S2))*r[3] + \
                 ((HX4X6S2)+(p-1.0)*(VX4X6S2))*r[4] + \
                 ((HX4X7S2)-t*(SX4X7S2)+(p-1.0)*(VX4X7S2))*r[5] + \
                 ((HX4S2S2)+(p-1.0)*(VX4S2S2))*s[1]*2.0
#define D2GDR2DT R*(0.5*log(xti4m1) - 0.5*log(xmg2m1) + 0.5*log(1.0-xti4m1) \
                 - log(xmg2m1+xfe2m1-xti4m1) + 0.5*log(xmg2m1+xfe2m1) \
                 + 0.5*log(1.0-xmg2m1-xfe2m1-xna1m2) - log(1.0-xsi4tet) \
                 - 0.5*log(1.0-xmg2m1-xfe2m1) + log(xfe3tet)) - \
                 (SX4) - (SX2X4)*r[0] - (SX3X4)*r[1] - (SX4X4)*r[2]*2.0 - \
                 (SX4X5)*r[3] - (SX4X6)*r[4] - (SX4X7)*r[5] - (SX4S1)*s[0] - \
                 (SX4S2)*s[1] - (SX4X7X7)*r[5]*r[5] - (SX4X7S2)*r[5]*s[1]
#define D2GDR2DP (VX4) + (VX2X4)*r[0] + (VX3X4)*r[1] + (VX4X4)*r[2]*2.0 + \
                 (VX4X5)*r[3] + (VX4X6)*r[4] + (VX4X7)*r[5] + (VX4S1)*s[0] + \
                 (VX4S2)*s[1] + (VX2X4X7)*r[0]*r[5] + (VX2X4S2)*r[0]*s[1] + \
                 (VX3X4X7)*r[1]*r[5] + (VX3X4S2)*r[1]*s[1] + \
                 (VX4X4X7)*r[2]*r[5]*2.0 + (VX4X4S2)*r[2]*s[1]*2.0 + \
                 (VX4X5X7)*r[3]*r[5] + (VX4X5S2)*r[3]*s[1] + \
                 (VX4X6X7)*r[4]*r[5] + (VX4X6S2)*r[4]*s[1] + \
                 (VX4X7X7)*r[5]*r[5] + (VX4X7S2)*r[5]*s[1] + \
                 (VX4S2S2)*s[1]*s[1]

#define D2GDR3R3 R*t*(0.25/xal3m1 + 0.25/xfe3m1 + 1.0/xmg2m1 \
                 + 1.0/(xmg2m1+xfe2m1-xti4m1) - 1.0/(xmg2m1+xfe2m1) \
                 - 0.5/(1.0-xsi4tet) + 0.125/xal3tet + 0.125/xfe3tet \
                 + 1.0/(1.0-xmg2m1-xfe2m1-xna1m2) - 1.0/(1.0-xmg2m1-xfe2m1)) + \
                 ((HX5X5)-t*(SX5X5)+(p-1.0)*(VX5X5))*2.0 + \
                 ((HX5X5X7)+(p-1.0)*(VX5X5X7))*r[5]*2.0 + \
                 ((HX5X5S2)+(p-1.0)*(VX5X5S2))*s[1]*2.0
#define D2GDR3R4 R*t*(0.125/xal3m1 + 0.125/xfe3m1 + 0.5/xmg2m1 \
                 + 0.5/(xmg2m1+xfe2m1-xti4m1) - 0.5/(xmg2m1+xfe2m1) \
                 + 0.25/(1.0-xsi4tet) - 0.0625/xal3tet - 0.0625/xfe3tet \
                 - 0.5/(1.0-xmg2m1-xfe2m1-xna1m2) - 0.5/(1.0-xmg2m1-xfe2m1)) + \
                 ((HX5X6)-t*(SX5X6)+(p-1.0)*(VX5X6)) + \
                 ((HX5X6X7)+(p-1.0)*(VX5X6X7))*r[5] + \
                 ((HX5X6S2)+(p-1.0)*(VX5X6S2))*s[1]
#define D2GDR3R5 -R*t*0.5/xmg2m1 + ((HX5X7)-t*(SX5X7)+(p-1.0)*(VX5X7)) + \
                 ((HX2X5X7)+(p-1.0)*(VX2X5X7))*r[0] + \
                 ((HX3X5X7)+(p-1.0)*(VX3X5X7))*r[1] + \
                 ((HX4X5X7)+(p-1.0)*(VX4X5X7))*r[2] + \
                 ((HX5X5X7)+(p-1.0)*(VX5X5X7))*r[3]*2.0 + \
                 ((HX5X6X7)+(p-1.0)*(VX5X6X7))*r[4] + \
                 ((HX5X7X7)-t*(SX5X7X7)+(p-1.0)*(VX5X7X7))*r[5]*2.0 + \
                 ((HX5X7S2)-t*(SX5X7S2)+(p-1.0)*(VX5X7S2))*s[1]
#define D2GDR3S0 R*t*0.25*(1.0/xfe3m1-1.0/xal3m1+0.5/xal3tet-0.5/xfe3tet) + \
                 ((HX5S1)-t*(SX5S1)+(p-1.0)*(VX5S1))
#define D2GDR3S1 -R*t*0.5/xmg2m1 + ((HX5S2)-t*(SX5S2)+(p-1.0)*(VX5S2)) + \
                 ((HX2X5S2)+(p-1.0)*(VX2X5S2))*r[0] + \
                 ((HX3X5S2)+(p-1.0)*(VX3X5S2))*r[1] + \
                 ((HX4X5S2)+(p-1.0)*(VX4X5S2))*r[2] + \
                 ((HX5X5S2)+(p-1.0)*(VX5X5S2))*r[3]*2.0 + \
                 ((HX5X6S2)+(p-1.0)*(VX5X6S2))*r[4] + \
                 ((HX5X7S2)-t*(SX5X7S2)+(p-1.0)*(VX5X7S2))*r[5] + \
                 ((HX5S2S2)+(p-1.0)*(VX5S2S2))*s[1]*2.0
#define D2GDR3DT R*(0.5*log(xal3m1) + 0.5*log(xfe3m1) - log(xmg2m1) \
                 - log(xmg2m1+xfe2m1-xti4m1) + log(xmg2m1+xfe2m1) \
                 - log(1.0-xsi4tet) + 0.5*log(xal3tet) + 0.5*log(xfe3tet) \
                 + log(1.0-xmg2m1-xfe2m1-xna1m2) - log(1.0-xmg2m1-xfe2m1)) - \
                 (SX5) - (SX2X5)*r[0] - (SX3X5)*r[1] - (SX4X5)*r[2] - \
                 (SX5X5)*r[3]*2.0 - (SX5X6)*r[4] - (SX5X7)*r[5] - \
                 (SX5S1)*s[0] - (SX5S2)*s[1] - (SX5X7X7)*r[5]*r[5] - \
                 (SX5X7S2)*r[5]*s[1]
#define D2GDR3DP (VX5) + (VX2X5)*r[0] + (VX3X5)*r[1] + (VX4X5)*r[2] + \
                 (VX5X5)*r[3]*2.0 + (VX5X6)*r[4] + (VX5X7)*r[5] + \
                 (VX5S1)*s[0] + (VX5S2)*s[1] + (VX2X5X7)*r[0]*r[5] + \
                 (VX2X5S2)*r[0]*s[1] + (VX3X5X7)*r[1]*r[5] + \
                 (VX3X5S2)*r[1]*s[1] + (VX4X5X7)*r[2]*r[5] + \
                 (VX4X5S2)*r[2]*s[1] + (VX5X5X7)*r[3]*r[5]*2.0 + \
                 (VX5X5S2)*r[3]*s[1]*2.0 + (VX5X6X7)*r[4]*r[5] + \
                 (VX5X6S2)*r[4]*s[1] + (VX5X7X7)*r[5]*r[5] + \
                 (VX5X7S2)*r[5]*s[1] + (VX5S2S2)*s[1]*s[1]

#define D2GDR4R4 R*t*(0.0625/xal3m1 + 0.0625/xfe3m1 + 0.25/xmg2m1 \
                 + 0.25/(xmg2m1+xfe2m1-xti4m1) - 0.25/(xmg2m1+xfe2m1) \
                 - 0.125/(1.0-xsi4tet) + 0.03125/xal3tet + 0.03125/xfe3tet \
                 + 1.0/xca2m2 + 1.0/xna1m2 + 0.25/(1.0-xmg2m1-xfe2m1-xna1m2) \
                 - 0.25/(1.0-xmg2m1-xfe2m1) - 1.0/(1.0-xna1m2)) + \
                 ((HX6X6)-t*(SX6X6)+(p-1.0)*(VX6X6))*2.0 + \
                 ((HX6X6X7)+(p-1.0)*(VX6X6X7))*r[5]*2.0 + \
                 ((HX6X6S2)+(p-1.0)*(VX6X6S2))*s[1]*2.0
#define D2GDR4R5 R*t*(-0.25/xmg2m1 + 1.0/xca2m2) + \
                 ((HX6X7)-t*(SX6X7)+(p-1.0)*(VX6X7)) + \
                 ((HX2X6X7)+(p-1.0)*(VX2X6X7))*r[0] + \
                 ((HX3X6X7)+(p-1.0)*(VX3X6X7))*r[1] + \
                 ((HX4X6X7)+(p-1.0)*(VX4X6X7))*r[2] + \
                 ((HX5X6X7)+(p-1.0)*(VX5X6X7))*r[3] + \
                 ((HX6X6X7)+(p-1.0)*(VX6X6X7))*r[4]*2.0 + \
                 ((HX6X7X7)-t*(SX6X7X7)+(p-1.0)*(VX6X7X7))*r[5]*2.0 + \
                 ((HX6X7S2)-t*(SX6X7S2)+(p-1.0)*(VX6X7S2))*s[1]
#define D2GDR4S0 R*t*(-0.125/xal3m1 + 0.125/xfe3m1 - 0.0625/xal3tet \
                 + 0.0625/xfe3tet) + ((HX6S1)-t*(SX6S1)+(p-1.0)*(VX6S1))
#define D2GDR4S1 -R*t*0.25/xmg2m1 + ((HX6S2)-t*(SX6S2)+(p-1.0)*(VX6S2)) + \
                 ((HX2X6S2)+(p-1.0)*(VX2X6S2))*r[0] + \
                 ((HX3X6S2)+(p-1.0)*(VX3X6S2))*r[1] + \
                 ((HX4X6S2)+(p-1.0)*(VX4X6S2))*r[2] + \
                 ((HX5X6S2)+(p-1.0)*(VX5X6S2))*r[3] + \
                 ((HX6X6S2)+(p-1.0)*(VX6X6S2))*r[4]*2.0 + \
                 ((HX6X7S2)-t*(SX6X7S2)+(p-1.0)*(VX6X7S2))*r[5] + \
                 ((HX6S2S2)+(p-1.0)*(VX6S2S2))*s[1]*2.0
#define D2GDR4DT R*(0.25*log(xal3m1) + 0.25*log(xfe3m1) - 0.5*log(xmg2m1) \
                 - 0.5*log(xmg2m1+xfe2m1-xti4m1) + 0.5*log(xmg2m1+xfe2m1) \
                 + 0.5*log(1.0-xsi4tet) - 0.25*log(xal3tet) \
                 - 0.25*log(xfe3tet) - log(xca2m2) + log(xna1m2) \
                 - 0.5*log(1.0-xmg2m1-xfe2m1-xna1m2) \
                 - 0.5*log(1.0-xmg2m1-xfe2m1) + log(1.0-xna1m2)) - \
                 (SX6) - (SX2X6)*r[0] - (SX3X6)*r[1] - (SX4X6)*r[2] - \
                 (SX5X6)*r[3] - (SX6X6)*r[4]*2.0 - (SX6X7)*r[5] - \
                 (SX6S1)*s[0] - (SX6S2)*s[1] - (SX6X7X7)*r[5]*r[5] - \
                 (SX6X7S2)*r[5]*s[1]
#define D2GDR4DP (VX6) + (VX2X6)*r[0] + (VX3X6)*r[1] + (VX4X6)*r[2] + \
                 (VX5X6)*r[3] + (VX6X6)*r[4]*2.0 + (VX6X7)*r[5] + \
                 (VX6S1)*s[0] + (VX6S2)*s[1] + (VX2X6X7)*r[0]*r[5] + \
                 (VX2X6S2)*r[0]*s[1] + (VX3X6X7)*r[1]*r[5] + \
                 (VX3X6S2)*r[1]*s[1] + (VX4X6X7)*r[2]*r[5] + \
                 (VX4X6S2)*r[2]*s[1] + (VX5X6X7)*r[3]*r[5] + \
                 (VX5X6S2)*r[3]*s[1] + (VX6X6X7)*r[4]*r[5]*2.0 + \
                 (VX6X6S2)*r[4]*s[1]*2.0 + (VX6X7X7)*r[5]*r[5] + \
                 (VX6X7S2)*r[5]*s[1] + (VX6S2S2)*s[1]*s[1]

#define D2GDR5R5 R*t*(0.25/xmg2m1 + 0.25/xfe2m1 + 1.0/xca2m2 + 0.25/xmg2m2 \
                 + 0.25/xfe2m2) + ((HX7X7)-t*(SX7X7)+(p-1.0)*(VX7X7))*2.0 + \
                 ((HX2X7X7)-t*(SX2X7X7)+(p-1.0)*(VX2X7X7))*r[0]*2.0 + \
                 ((HX3X7X7)-t*(SX3X7X7)+(p-1.0)*(VX3X7X7))*r[1]*2.0 + \
                 ((HX4X7X7)-t*(SX4X7X7)+(p-1.0)*(VX4X7X7))*r[2]*2.0 + \
                 ((HX5X7X7)-t*(SX5X7X7)+(p-1.0)*(VX5X7X7))*r[3]*2.0 + \
                 ((HX6X7X7)-t*(SX6X7X7)+(p-1.0)*(VX6X7X7))*r[4]*2.0 + \
                 ((HX7X7X7)-t*(SX7X7X7)+(p-1.0)*(VX7X7X7))*r[5]*6.0 + \
                 ((HX7X7S2)-t*(SX7X7S2)+(p-1.0)*(VX7X7S2))*s[1]*2.0
#define D2GDR5S0 ((HX7S1)-t*(SX7S1)+(p-1.0)*(VX7S1))
#define D2GDR5S1 R*t*(0.25/xmg2m1 + 0.25/xfe2m1 - 0.25/xmg2m2 + 0.25/xfe2m2) + \
                 ((HX7S2)-t*(SX7S2)+(p-1.0)*(VX7S2)) + \
                 ((HX2X7S2)+(p-1.0)*(VX2X7S2))*r[0] + \
                 ((HX3X7S2)-t*(SX3X7S2)+(p-1.0)*(VX3X7S2))*r[1] + \
                 ((HX4X7S2)-t*(SX4X7S2)+(p-1.0)*(VX4X7S2))*r[2] + \
                 ((HX5X7S2)-t*(SX5X7S2)+(p-1.0)*(VX5X7S2))*r[3] + \
                 ((HX6X7S2)-t*(SX6X7S2)+(p-1.0)*(VX6X7S2))*r[4] + \
                 ((HX7X7S2)-t*(SX7X7S2)+(p-1.0)*(VX7X7S2))*r[5]*2.0 + \
                 ((HX7S2S2)+(p-1.0)*(VX7S2S2))*s[1]*2.0
#define D2GDR5DT R*(0.5*log(xmg2m1) - 0.5*log(xfe2m1) - log(xca2m2) \
                 + 0.5*log(xmg2m2) + 0.5*log(xfe2m2)) - \
                 (SX7) - (SX2X7)*r[0] - (SX3X7)*r[1] - (SX4X7)*r[2] - \
                 (SX5X7)*r[3] - (SX6X7)*r[4] - (SX7X7)*r[5]*2.0 - \
                 (SX7S1)*s[0] - (SX7S2)*s[1] - (SX2X7X7)*r[0]*r[5]*2.0 - \
                 (SX3X7X7)*r[1]*r[5]*2.0 - (SX3X7S2)*r[1]*s[1] - \
                 (SX4X7X7)*r[2]*r[5]*2.0 - (SX4X7S2)*r[2]*s[1] - \
                 (SX5X7X7)*r[3]*r[5]*2.0 - (SX5X7S2)*r[3]*s[1] - \
                 (SX6X7X7)*r[4]*r[5]*2.0 - (SX6X7S2)*r[4]*s[1] - \
                 (SX7X7X7)*r[5]*r[5]*3.0 - (SX7X7S2)*r[5]*s[1]*2.0
#define D2GDR5DP (VX7) + (VX2X7)*r[0] + (VX3X7)*r[1] + (VX4X7)*r[2] + \
                 (VX5X7)*r[3] + (VX6X7)*r[4] + (VX7X7)*r[5]*2.0 + \
                 (VX7S1)*s[0] + (VX7S2)*s[1] + (VX2X2X7)*r[0]*r[0] + \
                 (VX2X3X7)*r[0]*r[1] + (VX2X4X7)*r[0]*r[2] + \
                 (VX2X5X7)*r[0]*r[3] + (VX2X6X7)*r[0]*r[4] + \
                 (VX2X7X7)*r[0]*r[5]*2.0 + (VX2X7S2)*r[0]*s[1] + \
                 (VX3X3X7)*r[1]*r[1] + (VX3X4X7)*r[1]*r[2] + \
                 (VX3X5X7)*r[1]*r[3] + (VX3X6X7)*r[1]*r[4] + \
                 (VX3X7X7)*r[1]*r[5]*2.0 + (VX3X7S2)*r[1]*s[1] + \
                 (VX4X4X7)*r[2]*r[2] + (VX4X5X7)*r[2]*r[3] + \
                 (VX4X6X7)*r[2]*r[4] + (VX4X7X7)*r[2]*r[5]*2.0 + \
                 (VX4X7S2)*r[2]*s[1] + (VX5X5X7)*r[3]*r[3] + \
                 (VX5X6X7)*r[3]*r[4] + (VX5X7X7)*r[3]*r[5]*2.0 + \
                 (VX5X7S2)*r[3]*s[1] + (VX6X6X7)*r[4]*r[4] + \
                 (VX6X7X7)*r[4]*r[5]*2.0 + (VX6X7S2)*r[4]*s[1] + \
                 (VX7X7X7)*r[5]*r[5]*3.0 + (VX7X7S2)*r[5]*s[1]*2.0 + \
                 (VX7S2S2)*s[1]*s[1]

#define D2GDS0S0 R*t*(0.25/xfe3m1+0.25/xal3m1+0.125/xal3tet+0.125/xfe3tet) + \
                 ((HS1S1)-t*(SS1S1)+(p-1.0)*(VS1S1))*2.0
#define D2GDS0S1 ((HS1S2)-t*(SS1S2)+(p-1.0)*(VS1S2))
#define D2GDS0DT R*0.5*(log(xfe3m1)-log(xal3m1)+log(xal3tet)-log(xfe3tet)) - \
                 (SS1) - (SX2S1)*r[0] - (SX3S1)*r[1] -(SX4S1)*r[2] - \
                 (SX5S1)*r[3] - (SX6S1)*r[4] - (SX7S1)*r[5] - \
                 (SS1S1)*s[0]*2.0 - (SS1S2)*s[1]
#define D2GDS0DP (VS1) + (VX2S1)*r[0] + (VX3S1)*r[1] + (VX4S1)*r[2] + \
                 (VX5S1)*r[3] + (VX6S1)*r[4] + (VX7S1)*r[5] + \
                 (VS1S1)*s[0]*2.0 + (VS1S2)*s[1]

#define D2GDS1S1 R*t*(0.25/xmg2m1+0.25/xfe2m1+0.25/xfe2m2+0.25/xmg2m2) + \
                 ((HS2S2)-t*(SS2S2)+(p-1.0)*(VS2S2))*2.0 + \
                 ((HX2S2S2)+(p-1.0)*(VX2S2S2))*r[0]*2.0 + \
                 ((HX3S2S2)+(p-1.0)*(VX3S2S2))*r[1]*2.0 + \
                 ((HX4S2S2)+(p-1.0)*(VX4S2S2))*r[2]*2.0 + \
                 ((HX5S2S2)+(p-1.0)*(VX5S2S2))*r[3]*2.0 + \
                 ((HX6S2S2)+(p-1.0)*(VX6S2S2))*r[4]*2.0 + \
                 ((HX7S2S2)+(p-1.0)*(VX7S2S2))*r[5]*2.0
#define D2GDS1DT R*0.5*(log(xmg2m1)-log(xfe2m1)+log(xfe2m2)-log(xmg2m2)) - \
                 (SS2) - (SX2S2)*r[0] - (SX3S2)*r[1] - (SX4S2)*r[2] - \
                 (SX5S2)*r[3] - (SX6S2)*r[4] - (SX7S2)*r[5] - (SS1S2)*s[0] - \
                 (SS2S2)*s[1]*2.0 - (SX3X7S2)*r[1]*r[5] - \
                 (SX4X7S2)*r[2]*r[5] - (SX5X7S2)*r[3]*r[5] - \
                 (SX6X7S2)*r[4]*r[5] - (SX7X7S2)*r[5]*r[5]
#define D2GDS1DP (VS2) + (VX2S2)*r[0] + (VX3S2)*r[1] + (VX4S2)*r[2] + \
                 (VX5S2)*r[3] + (VX6S2)*r[4] + (VX7S2)*r[5] + (VS1S2)*s[0] + \
                 (VS2S2)*s[1]*2.0 + (VX2X2S2)*r[0]*r[0] + \
                 (VX2X3S2)*r[0]*r[1] + (VX2X4S2)*r[0]*r[2] + \
                 (VX2X5S2)*r[0]*r[3] + (VX2X6S2)*r[0]*r[4] + \
                 (VX2X7S2)*r[0]*r[5] + (VX2S2S2)*r[0]*s[1]*2.0 + \
                 (VX3X3S2)*r[1]*r[1] + (VX3X4S2)*r[1]*r[2] + \
                 (VX3X5S2)*r[1]*r[3] + (VX3X6S2)*r[1]*r[4] + \
                 (VX3X7S2)*r[1]*r[5] + (VX3S2S2)*r[1]*s[1]*2.0 + \
                 (VX4X4S2)*r[2]*r[2] + (VX4X5S2)*r[2]*r[3] + \
                 (VX4X6S2)*r[2]*r[4] + (VX4X7S2)*r[2]*r[5] + \
                 (VX4S2S2)*r[2]*s[1]*2.0 + (VX5X5S2)*r[3]*r[3] + \
                 (VX5X6S2)*r[3]*r[4] + (VX5X7S2)*r[3]*r[5] + \
                 (VX5S2S2)*r[3]*s[1]*2.0 + (VX6X6S2)*r[4]*r[4] + \
                 (VX6X7S2)*r[4]*r[5] + (VX6S2S2)*r[4]*s[1]*2.0 + \
                 (VX7X7S2)*r[5]*r[5] + (VX7S2S2)*r[5]*s[1]*2.0

#define D2GDT2   0.0
#define D2GDTDP  0.0
#define D2GDP2   0.0

/*----------------------------------------------------------------------------*/

#define D3GDR0R0R0 R*t*(-1.0/SQUARE(xfe2m1) + 1.0/SQUARE(xmg2m1) )
#define D3GDR0R0R1 R*t*(0.5/SQUARE(xmg2m1) )
#define D3GDR0R0R2 R*t*(0.5/SQUARE(xmg2m1) )
#define D3GDR0R0R3 R*t*(1.0/SQUARE(xmg2m1) )
#define D3GDR0R0R4 R*t*(0.5/SQUARE(xmg2m1) )
#define D3GDR0R0R5 R*t*(0.5/SQUARE(xfe2m1) - 0.5/SQUARE(xmg2m1) ) \
		   + 2.0*((HX2X2X7) + (p-1.0)*(VX2X2X7) )
#define D3GDR0R0S0 0.0
#define D3GDR0R0S1 -R*t*(-0.5/SQUARE(xfe2m1) + 0.5/SQUARE(xmg2m1)) + \
                   ((HX2X2S2)+(p-1.0)*(VX2X2S2))*2.0
#define D3GDR0R0DT R*(1.0/xfe2m1 + 1.0/xmg2m1) - (SX2X2)*2.0
#define D3GDR0R0DP (VX2X2)*2.0 + (VX2X2X7)*r[5]*2.0 + (VX2X2S2)*s[1]*2.0

#define D3GDR0R1R1 R*t*(0.25/SQUARE(xmg2m1) )
#define D3GDR0R1R2 R*t*(0.25/SQUARE(xmg2m1) )
#define D3GDR0R1R3 R*t*(0.5/SQUARE(xmg2m1) )
#define D3GDR0R1R4 R*t*(0.25/SQUARE(xmg2m1) )
#define D3GDR0R1R5 -R*t*0.25/SQUARE(xmg2m1) + ((HX2X3X7)+(p-1.0)*(VX2X3X7))
#define D3GDR0R1S0 0.0
#define D3GDR0R1S1 -R*t*0.25/SQUARE(xmg2m1) + ((HX2X3S2)+(p-1.0)*(VX2X3S2))
#define D3GDR0R1DT R*(0.5/xmg2m1) - (SX2X3)
#define D3GDR0R1DP (VX2X3) + (VX2X3X7)*r[5] + (VX2X3S2)*s[1]

#define D3GDR0R2R2 R*t*(0.25/SQUARE(xmg2m1) )
#define D3GDR0R2R3 R*t*(0.5/SQUARE(xmg2m1) )
#define D3GDR0R2R4 R*t*(0.25/SQUARE(xmg2m1) )
#define D3GDR0R2R5 -R*t*0.25/SQUARE(xmg2m1) + ((HX2X4X7)+(p-1.0)*(VX2X4X7))
#define D3GDR0R2S0 0.0
#define D3GDR0R2S1 -R*t*0.25/SQUARE(xmg2m1) + ((HX2X4S2)+(p-1.0)*(VX2X4S2))
#define D3GDR0R2DT R*(0.5/xmg2m1) - (SX2X4)
#define D3GDR0R2DP (VX2X4) + (VX2X4X7)*r[5] + (VX2X4S2)*s[1]

#define D3GDR0R3R3 R*t*(1.0/SQUARE(xmg2m1) )
#define D3GDR0R3R4 R*t*(0.5/SQUARE(xmg2m1) )
#define D3GDR0R3R5 -R*t*0.5/SQUARE(xmg2m1) + ((HX2X5X7)+(p-1.0)*(VX2X5X7))
#define D3GDR0R3S0 0.0
#define D3GDR0R3S1 -R*t*0.5/SQUARE(xmg2m1) + ((HX2X5S2)+(p-1.0)*(VX2X5S2))
#define D3GDR0R3DT R*(1.0/xmg2m1) - (SX2X5)
#define D3GDR0R3DP (VX2X5) + (VX2X5X7)*r[5] + (VX2X5S2)*s[1]

#define D3GDR0R4R4 R*t*(0.25/SQUARE(xmg2m1) )
#define D3GDR0R4R5 -R*t*0.25/SQUARE(xmg2m1) + ((HX2X6X7)+(p-1.0)*(VX2X6X7))
#define D3GDR0R4S0 0.0
#define D3GDR0R4S1 -R*t*0.25/SQUARE(xmg2m1) + ((HX2X6S2)+(p-1.0)*(VX2X6S2))
#define D3GDR0R4DT R*(0.5/xmg2m1) - (SX2X6)
#define D3GDR0R4DP (VX2X6) + (VX2X6X7)*r[5] + (VX2X6S2)*s[1]

#define D3GDR0R5R5 R*t*(-0.25/SQUARE(xfe2m1) + 0.25/SQUARE(xmg2m1) ) \
		   + 2.0*((HX2X7X7)-t*(SX2X7X7)+(p-1.0)*(VX2X7X7) )
#define D3GDR0R5S0 0.0
#define D3GDR0R5S1 -R*t*(0.25/SQUARE(xfe2m1) - 0.25/SQUARE(xmg2m1)) + \
                   ((HX2X7S2)+(p-1.0)*(VX2X7S2))
#define D3GDR0R5DT -R*(0.5/xfe2m1 + 0.5/xmg2m1) - (SX2X7) - (SX2X7X7)*r[5]*2.0
#define D3GDR0R5DP (VX2X7) + (VX2X2X7)*r[0]*2.0 + (VX2X3X7)*r[1] + \
                   (VX2X4X7)*r[2] + (VX2X5X7)*r[3] + (VX2X6X7)*r[4] + \
                   (VX2X7X7)*r[5]*2.0 + (VX2X7S2)*s[1]

#define D3GDR0S0S0 0.0
#define D3GDR0S0S1 0.0
#define D3GDR0S0DT -(SX2S1)
#define D3GDR0S0DP (VX2S1)

#define D3GDR0S1S1 -R*t*(0.25/SQUARE(xfe2m1) - 0.25/SQUARE(xmg2m1)) + \
                   ((HX2S2S2)+(p-1.0)*(VX2S2S2))*2.0
#define D3GDR0S1DT -R*(0.5/xfe2m1 + 0.5/xmg2m1) - (SX2S2)
#define D3GDR0S1DP (VX2S2) + (VX2X2S2)*r[0]*2.0 + (VX2X3S2)*r[1] + \
                   (VX2X4S2)*r[2] + (VX2X5S2)*r[3] + (VX2X6S2)*r[4] + \
                   (VX2X7S2)*r[5] + (VX2S2S2)*s[1]*2.0

#define D3GDR1R1R1 R*t*(-0.125/SQUARE(xti4m1) + 0.125/SQUARE(xmg2m1) \
		   -0.125/SQUARE(1.0-xti4m1)+1.0/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.125/SQUARE(xmg2m1+xfe2m1) + 0.25/SQUARE(1.0-xsi4tet) \
                   - 0.125/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.125/SQUARE(1.0-xmg2m1-xfe2m1) - 0.25/SQUARE(xal3tet) )
#define D3GDR1R1R2 R*t*(-0.125/SQUARE(xti4m1) + 0.125/SQUARE(xmg2m1) \
		   -0.125/SQUARE(1.0-xti4m1)+1.0/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.125/SQUARE(xmg2m1+xfe2m1) + 0.25/SQUARE(1.0-xsi4tet) \
                   - 0.125/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.125/SQUARE(1.0-xmg2m1-xfe2m1) )
#define D3GDR1R1R3 R*t*(0.25/SQUARE(xmg2m1)+1.0/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.25/SQUARE(xmg2m1+xfe2m1) + 0.25/SQUARE(1.0-xsi4tet) \
                   - 0.25/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.25/SQUARE(1.0-xmg2m1-xfe2m1) - 0.125/SQUARE(xal3tet) )
#define D3GDR1R1R4 R*t*(0.125/SQUARE(xmg2m1)+0.5/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.125/SQUARE(xmg2m1+xfe2m1) - 0.125/SQUARE(1.0-xsi4tet) \
                   + 0.125/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.125/SQUARE(1.0-xmg2m1-xfe2m1) + 0.0625/SQUARE(xal3tet) )
#define D3GDR1R1R5 -R*t*0.125/SQUARE(xmg2m1)+((HX3X3X7)+(p-1.0)*(VX3X3X7))*2.0
#define D3GDR1R1S0 -R*t*0.125/SQUARE(xal3tet)
#define D3GDR1R1S1 -R*t*0.125/SQUARE(xmg2m1) + ((HX3X3S2)+(p-1.0)*(VX3X3S2))*2.0
#define D3GDR1R1DT R*(0.25/xti4m1 + 0.25/xmg2m1 - 0.25/(1.0-xti4m1) \
                   + 1.0/(xmg2m1+xfe2m1-xti4m1) - 0.25/(xmg2m1+xfe2m1) \
                   + 0.25/(1.0-xmg2m1-xfe2m1-xna1m2) - 0.5/(1.0-xsi4tet) \
                   - 0.25/(1.0-xmg2m1-xfe2m1) + 0.5/xal3tet) - (SX3X3)*2.0
#define D3GDR1R1DP (VX3X3)*2.0 + (VX3X3X7)*r[5]*2.0 + (VX3X3S2)*s[1]*2.0

#define D3GDR1R2R2 R*t*(-0.125/SQUARE(xti4m1) + 0.125/SQUARE(xmg2m1) \
		   -0.125/SQUARE(1.0-xti4m1)+1.0/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.125/SQUARE(xmg2m1+xfe2m1) + 0.25/SQUARE(1.0-xsi4tet) \
                   - 0.125/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.125/SQUARE(1.0-xmg2m1-xfe2m1) )
#define D3GDR1R2R3 R*t*(0.25/SQUARE(xmg2m1)+1.0/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.25/SQUARE(xmg2m1+xfe2m1) + 0.25/SQUARE(1.0-xsi4tet) \
                   - 0.25/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.25/SQUARE(1.0-xmg2m1-xfe2m1) )
#define D3GDR1R2R4 R*t*(0.125/SQUARE(xmg2m1)+0.5/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.125/SQUARE(xmg2m1+xfe2m1) - 0.125/SQUARE(1.0-xsi4tet) \
                   + 0.125/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.125/SQUARE(1.0-xmg2m1-xfe2m1) )
#define D3GDR1R2R5 -R*t*0.125/SQUARE(xmg2m1)+((HX3X4X7)+(p-1.0)*(VX3X4X7))
#define D3GDR1R2S0 0.0
#define D3GDR1R2S1 -R*t*0.125/SQUARE(xmg2m1) + ((HX3X4S2)+(p-1.0)*(VX3X4S2))
#define D3GDR1R2DT R*(0.25/xti4m1 + 0.25/xmg2m1 - 0.25/(1.0-xti4m1) \
                   + 1.0/(xmg2m1+xfe2m1-xti4m1) - 0.25/(xmg2m1+xfe2m1) \
                   + 0.25/(1.0-xmg2m1-xfe2m1-xna1m2) - 0.5/(1.0-xsi4tet) \
                   - 0.25/(1.0-xmg2m1-xfe2m1)) - (SX3X4)
#define D3GDR1R2DP (VX3X4) + (VX3X4X7)*r[5] + (VX3X4S2)*s[1]

#define D3GDR1R3R3 R*t*(0.5/SQUARE(xmg2m1)+1.0/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.5/SQUARE(xmg2m1+xfe2m1) + 0.25/SQUARE(1.0-xsi4tet) \
                   - 0.5/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.5/SQUARE(1.0-xmg2m1-xfe2m1) - 0.0625/SQUARE(xal3tet) )
#define D3GDR1R3R4 R*t*(0.25/SQUARE(xmg2m1)+0.5/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.25/SQUARE(xmg2m1+xfe2m1) - 0.125/SQUARE(1.0-xsi4tet) \
                   + 0.25/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.25/SQUARE(1.0-xmg2m1-xfe2m1) + 0.03125/SQUARE(xal3tet) )
#define D3GDR1R3R5 -R*t*0.25/SQUARE(xmg2m1)+((HX3X5X7)+(p-1.0)*(VX3X5X7))
#define D3GDR1R3S0 -R*t*0.0625/SQUARE(xal3tet)
#define D3GDR1R3S1 -R*t*0.25/SQUARE(xmg2m1) + ((HX3X5S2)+(p-1.0)*(VX3X5S2))
#define D3GDR1R3DT R*(0.5/xmg2m1 + 1.0/(xmg2m1+xfe2m1-xti4m1) \
                   - 0.5/(xmg2m1+xfe2m1) + 0.5/(1.0-xmg2m1-xfe2m1-xna1m2) \
                   - 0.5/(1.0-xsi4tet) - 0.5/(1.0-xmg2m1-xfe2m1) \
                   + 0.25/xal3tet) - (SX3X5)
#define D3GDR1R3DP (VX3X5) + (VX3X5X7)*r[5] + (VX3X5S2)*s[1]

#define D3GDR1R4R4 R*t*(0.125/SQUARE(xmg2m1)+0.25/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.125/SQUARE(xmg2m1+xfe2m1) + 0.0625/SQUARE(1.0-xsi4tet) \
                   - 0.125/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.125/SQUARE(1.0-xmg2m1-xfe2m1) - 0.015625/SQUARE(xal3tet) )
#define D3GDR1R4R5 -R*t*0.125/SQUARE(xmg2m1)+((HX3X6X7)+(p-1.0)*(VX3X6X7))
#define D3GDR1R4S0 R*t*0.03125/SQUARE(xal3tet)
#define D3GDR1R4S1 -R*t*0.125/SQUARE(xmg2m1) + ((HX3X6S2)+(p-1.0)*(VX3X6S2))
#define D3GDR1R4DT R*(0.25/xmg2m1 + 0.5/(xmg2m1+xfe2m1-xti4m1) \
                   - 0.25/(xmg2m1+xfe2m1) - 0.25/(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.25/(1.0-xsi4tet) - 0.25/(1.0-xmg2m1-xfe2m1) \
                   - 0.125/xal3tet) - (SX3X6)
#define D3GDR1R4DP (VX3X6) + (VX3X6X7)*r[5] + (VX3X6S2)*s[1]

#define D3GDR1R5R5 R*t*0.125/SQUARE(xmg2m1)+ \
		   ((HX3X7X7)-t*(SX3X7X7)+(p-1.0)*(VX3X7X7))*2.0
#define D3GDR1R5S0 0.0
#define D3GDR1R5S1 R*t*0.125/SQUARE(xmg2m1) + \
                   ((HX3X7S2)-t*(SX3X7S2)+(p-1.0)*(VX3X7S2))
#define D3GDR1R5DT -R*0.25/xmg2m1 - (SX3X7) - (SX3X7X7)*r[5]*2.0 - \
                   (SX3X7S2)*s[1]
#define D3GDR1R5DP (VX3X7) + (VX2X3X7)*r[0] + (VX3X3X7)*r[1]*2.0 + \
                   (VX3X4X7)*r[2] + (VX3X5X7)*r[3] + (VX3X6X7)*r[4] + \
                   (VX3X7X7)*r[5]*2.0 + (VX3X7S2)*s[1]

#define D3GDR1S0S0 -R*t*0.0625/SQUARE(xal3tet)
#define D3GDR1S0S1 0.0
#define D3GDR1S0DT R*0.25/xal3tet - (SX3S1)
#define D3GDR1S0DP (VX3S1)

#define D3GDR1S1S1 R*t*0.125/SQUARE(xmg2m1) + ((HX3S2S2)+(p-1.0)*(VX3S2S2))*2.0
#define D3GDR1S1DT -R*0.25/xmg2m1 - (SX3S2) - (SX3X7S2)*r[5]
#define D3GDR1S1DP (VX3S2) + (VX2X3S2)*r[0] + (VX3X3S2)*r[1]*2.0 + \
                   (VX3X4S2)*r[2] + (VX3X5S2)*r[3] + (VX3X6S2)*r[4] + \
                   (VX3X7S2)*r[5] + (VX3S2S2)*s[1]*2.0

#define D3GDR2R2R2 R*t*(-0.125/SQUARE(xti4m1) + 0.125/SQUARE(xmg2m1) \
		   -0.125/SQUARE(1.0-xti4m1)+1.0/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.125/SQUARE(xmg2m1+xfe2m1) + 0.25/SQUARE(1.0-xsi4tet) \
                   - 0.125/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.125/SQUARE(1.0-xmg2m1-xfe2m1) - 0.25/SQUARE(xfe3tet) )
#define D3GDR2R2R3 R*t*(0.25/SQUARE(xmg2m1)+1.0/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.25/SQUARE(xmg2m1+xfe2m1) + 0.25/SQUARE(1.0-xsi4tet) \
                   - 0.25/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.25/SQUARE(1.0-xmg2m1-xfe2m1) - 0.125/SQUARE(xfe3tet) )
#define D3GDR2R2R4 R*t*(0.125/SQUARE(xmg2m1)+0.5/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.125/SQUARE(xmg2m1+xfe2m1) - 0.125/SQUARE(1.0-xsi4tet) \
                   + 0.125/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.125/SQUARE(1.0-xmg2m1-xfe2m1) + 0.0625/SQUARE(xfe3tet) )
#define D3GDR2R2R5 -R*t*0.125/SQUARE(xmg2m1)+((HX4X4X7)+(p-1.0)*(VX4X4X7))*2.0
#define D3GDR2R2S0 R*t*0.125/SQUARE(xfe3tet)
#define D3GDR2R2S1 -R*t*0.125/SQUARE(xmg2m1) + ((HX4X4S2)+(p-1.0)*(VX4X4S2))*2.0
#define D3GDR2R2DT R*(0.25/xti4m1 + 0.25/xmg2m1 - 0.25/(1.0-xti4m1) \
                   + 1.0/(xmg2m1+xfe2m1-xti4m1) - 0.25/(xmg2m1+xfe2m1) \
                   + 0.25/(1.0-xmg2m1-xfe2m1-xna1m2) - 0.5/(1.0-xsi4tet) \
                   - 0.25/(1.0-xmg2m1-xfe2m1) + 0.5/xfe3tet) - (SX4X4)*2.0
#define D3GDR2R2DP (VX4X4)*2.0 + (VX4X4X7)*r[5]*2.0 + (VX4X4S2)*s[1]*2.0

#define D3GDR2R3R3 R*t*(0.5/SQUARE(xmg2m1)+1.0/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.5/SQUARE(xmg2m1+xfe2m1) + 0.25/SQUARE(1.0-xsi4tet) \
                   - 0.5/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.5/SQUARE(1.0-xmg2m1-xfe2m1) - 0.0625/SQUARE(xfe3tet) )
#define D3GDR2R3R4 R*t*(0.25/SQUARE(xmg2m1)+0.5/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.25/SQUARE(xmg2m1+xfe2m1) - 0.125/SQUARE(1.0-xsi4tet) \
                   + 0.25/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.25/SQUARE(1.0-xmg2m1-xfe2m1) + 0.03125/SQUARE(xfe3tet) )
#define D3GDR2R3R5 -R*t*0.25/SQUARE(xmg2m1)+((HX4X5X7)+(p-1.0)*(VX4X5X7))
#define D3GDR2R3S0 R*t*0.0625/SQUARE(xfe3tet)
#define D3GDR2R3S1 -R*t*0.25/SQUARE(xmg2m1) + ((HX4X5S2)+(p-1.0)*(VX4X5S2))
#define D3GDR2R3DT R*(0.5/xmg2m1 + 1.0/(xmg2m1+xfe2m1-xti4m1) \
                   - 0.5/(xmg2m1+xfe2m1) + 0.5/(1.0-xmg2m1-xfe2m1-xna1m2) \
                   - 0.5/(1.0-xsi4tet) - 0.5/(1.0-xmg2m1-xfe2m1) \
                   + 0.25/xfe3tet) - (SX4X5)
#define D3GDR2R3DP (VX4X5) + (VX4X5X7)*r[5] + (VX4X5S2)*s[1]

#define D3GDR2R4R4 R*t*(0.125/SQUARE(xmg2m1)+0.25/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.125/SQUARE(xmg2m1+xfe2m1) + 0.0625/SQUARE(1.0-xsi4tet) \
                   - 0.125/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.125/SQUARE(1.0-xmg2m1-xfe2m1) - 0.015625/SQUARE(xfe3tet) )
#define D3GDR2R4R5 -R*t*0.125/SQUARE(xmg2m1)+((HX4X6X7)+(p-1.0)*(VX4X6X7))
#define D3GDR2R4S0 -R*t*0.03125/SQUARE(xfe3tet)
#define D3GDR2R4S1 -R*t*0.125/SQUARE(xmg2m1) + ((HX4X6S2)+(p-1.0)*(VX4X6S2))
#define D3GDR2R4DT R*(0.25/xmg2m1 + 0.5/(xmg2m1+xfe2m1-xti4m1) \
                   - 0.25/(xmg2m1+xfe2m1) - 0.25/(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.25/(1.0-xsi4tet) - 0.25/(1.0-xmg2m1-xfe2m1) \
                   - 0.125/xfe3tet) - (SX4X6)
#define D3GDR2R4DP (VX4X6) + (VX4X6X7)*r[5] + (VX4X6S2)*s[1]

#define D3GDR2R5R5 R*t*0.125/SQUARE(xmg2m1)+ \
		   ((HX4X7X7)-t*(SX4X7X7)+(p-1.0)*(VX4X7X7))*2.0
#define D3GDR2R5S0 0.0
#define D3GDR2R5S1 R*t*0.125/SQUARE(xmg2m1) + \
                   ((HX4X7S2)-t*(SX4X7S2)+(p-1.0)*(VX4X7S2))
#define D3GDR2R5DT -R*0.25/xmg2m1 - (SX4X7) - (SX4X7X7)*r[5]*2.0 - \
                   (SX4X7S2)*s[1]
#define D3GDR2R5DP (VX4X7) + (VX2X4X7)*r[0] + (VX3X4X7)*r[1] + \
                   (VX4X4X7)*r[2]*2.0 + (VX4X5X7)*r[3] + (VX4X6X7)*r[4] + \
                   (VX4X7X7)*r[5]*2.0 + (VX4X7S2)*s[1]

#define D3GDR2S0S0 -R*t*0.0625/SQUARE(xfe3tet)
#define D3GDR2S0S1 0.0
#define D3GDR2S0DT -R*0.25/xfe3tet - (SX4S1)
#define D3GDR2S0DP (VX4S1)

#define D3GDR2S1S1 R*t*0.125/SQUARE(xmg2m1) + ((HX4S2S2)+(p-1.0)*(VX4S2S2))*2.0
#define D3GDR2S1DT -R*0.25/xmg2m1 - (SX4S2) - (SX4X7S2)*r[5]
#define D3GDR2S1DP (VX4S2) + (VX2X4S2)*r[0] + (VX3X4S2)*r[1] + \
                   (VX4X4S2)*r[2]*2.0 + (VX4X5S2)*r[3] + (VX4X6S2)*r[4] + \
                   (VX4X7S2)*r[5] + (VX4S2S2)*s[1]*2.0

#define D3GDR3R3R3 R*t*(-.125/SQUARE(xal3m1)-.125/SQUARE(xfe3m1)+ \
		   1.0/SQUARE(xmg2m1)+1.0/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -1.0/SQUARE(xmg2m1+xfe2m1) + 0.25/SQUARE(1.0-xsi4tet) \
		   -0.03125/SQUARE(xal3tet)-0.03125/SQUARE(xfe3tet) \
                   - 1.0/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 1.0/SQUARE(1.0-xmg2m1-xfe2m1) )
#define D3GDR3R3R4 R*t*(-.0625/SQUARE(xal3m1)-.0625/SQUARE(xfe3m1)+ \
		   0.5/SQUARE(xmg2m1)+0.5/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.5/SQUARE(xmg2m1+xfe2m1) - 0.125/SQUARE(1.0-xsi4tet) \
		   +0.015625/SQUARE(xal3tet)+0.015625/SQUARE(xfe3tet) \
                   + 0.5/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.5/SQUARE(1.0-xmg2m1-xfe2m1) )
#define D3GDR3R3R5 -R*t*0.5/SQUARE(xmg2m1)+((HX5X5X7)+(p-1.0)*(VX5X5X7))*2.0
#define D3GDR3R3S0 -R*t*(-0.125/SQUARE(xal3m1) + 0.125/SQUARE(xfe3m1) \
                   + 0.03125/SQUARE(xal3tet) - 0.03125/SQUARE(xfe3tet))
#define D3GDR3R3S1 -R*t*0.5/SQUARE(xmg2m1) + ((HX5X5S2)+(p-1.0)*(VX5X5S2))*2.0
#define D3GDR3R3DT R*(0.25/xal3m1 + 0.25/xfe3m1 + 1.0/xmg2m1 \
                   + 1.0/(xmg2m1+xfe2m1-xti4m1) - 1.0/(xmg2m1+xfe2m1) \
                   - 0.5/(1.0-xsi4tet) + 0.125/xal3tet + 0.125/xfe3tet \
                   + 1.0/(1.0-xmg2m1-xfe2m1-xna1m2) \
                   - 1.0/(1.0-xmg2m1-xfe2m1)) - (SX5X5)*2.0
#define D3GDR3R3DP (VX5X5)*2.0 + (VX5X5X7)*r[5]*2.0 + (VX5X5S2)*s[1]*2.0

#define D3GDR3R4R4 R*t*(-.03125/SQUARE(xal3m1)-.03125/SQUARE(xfe3m1)+ \
		   0.25/SQUARE(xmg2m1)+0.25/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.25/SQUARE(xmg2m1+xfe2m1) + 0.0625/SQUARE(1.0-xsi4tet) \
		   -0.0078125/SQUARE(xal3tet)-0.0078125/SQUARE(xfe3tet) \
                   - 0.25/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.25/SQUARE(1.0-xmg2m1-xfe2m1) )
#define D3GDR3R4R5 -R*t*0.25/SQUARE(xmg2m1)+((HX5X6X7)+(p-1.0)*(VX5X6X7))
#define D3GDR3R4S0 -R*t*(-0.0625/SQUARE(xal3m1) + 0.0625/SQUARE(xfe3m1) \
                   - 0.015625/SQUARE(xal3tet) + 0.015625/SQUARE(xfe3tet))
#define D3GDR3R4S1 -R*t*0.25/SQUARE(xmg2m1) + ((HX5X6S2)+(p-1.0)*(VX5X6S2))
#define D3GDR3R4DT R*(0.125/xal3m1 + 0.125/xfe3m1 + 0.5/xmg2m1 \
                   + 0.5/(xmg2m1+xfe2m1-xti4m1) - 0.5/(xmg2m1+xfe2m1) \
                   + 0.25/(1.0-xsi4tet) - 0.0625/xal3tet - 0.0625/xfe3tet \
                   - 0.5/(1.0-xmg2m1-xfe2m1-xna1m2) \
                   - 0.5/(1.0-xmg2m1-xfe2m1)) - (SX5X6)
#define D3GDR3R4DP (VX5X6) + (VX5X6X7)*r[5] + (VX5X6S2)*s[1]

#define D3GDR3R5R5 R*t*0.25/SQUARE(xmg2m1)+ \
		   ((HX5X7X7)-t*(SX5X7X7)+(p-1.0)*(VX5X7X7))*2.0
#define D3GDR3R5S0 0.0
#define D3GDR3R5S1 R*t*0.25/SQUARE(xmg2m1) + \
                   ((HX5X7S2)-t*(SX5X7S2)+(p-1.0)*(VX5X7S2))
#define D3GDR3R5DT -R*0.5/xmg2m1 - (SX5X7) - (SX5X7X7)*r[5]*2.0 - \
                   (SX5X7S2)*s[1]
#define D3GDR3R5DP (VX5X7) + (VX2X5X7)*r[0] + (VX3X5X7)*r[1] + \
                   (VX4X5X7)*r[2] + (VX5X5X7)*r[3]*2.0 + (VX5X6X7)*r[4] + \
                   (VX5X7X7)*r[5]*2.0 + (VX5X7S2)*s[1]

#define D3GDR3S0S0 -R*t*(0.125/SQUARE(xal3m1) + 0.125/SQUARE(xfe3m1) \
                   + 0.03125/SQUARE(xal3tet) + 0.03125/SQUARE(xfe3tet))
#define D3GDR3S0S1 0.0
#define D3GDR3S0DT R*(-0.25/xal3m1 + 0.25/xfe3m1 + 0.125/xal3tet \
                   - 0.125/xfe3tet) - (SX5S1)
#define D3GDR3S0DP (VX5S1)

#define D3GDR3S1S1 R*t*0.25/SQUARE(xmg2m1) + ((HX5S2S2)+(p-1.0)*(VX5S2S2))*2.0
#define D3GDR3S1DT -R*0.5/xmg2m1 - (SX5S2) - (SX5X7S2)*r[5]
#define D3GDR3S1DP (VX5S2) + (VX2X5S2)*r[0] + (VX3X5S2)*r[1] + \
                   (VX4X5S2)*r[2] + (VX5X5S2)*r[3]*2.0 + (VX5X6S2)*r[4] + \
                   (VX5X7S2)*r[5] + (VX5S2S2)*s[1]*2.0

#define D3GDR4R4R4 R*t*(-.015625/SQUARE(xal3m1)-.015625/SQUARE(xfe3m1)+ \
		   0.125/SQUARE(xmg2m1)+0.125/SQUARE(xmg2m1+xfe2m1-xti4m1) \
            -0.125/SQUARE(xmg2m1+xfe2m1) - 0.03125/SQUARE(1.0-xsi4tet) \
		   +0.00390625/SQUARE(xal3tet)+0.00390625/SQUARE(xfe3tet) \
                   +1.0/SQUARE(xca2m2)-1.0/SQUARE(xna1m2) \
		   + 0.125/SQUARE(1.0-xmg2m1-xfe2m1-xna1m2) \
                   + 0.125/SQUARE(1.0-xmg2m1-xfe2m1) - 1.0/SQUARE(1.0-xna1m2) )
#define D3GDR4R4R5 R*t*(-0.125/SQUARE(xmg2m1)+1.0/SQUARE(xca2m2) ) + \
		   2.0*((HX6X6X7)+(p-1.0)*(VX6X6X7))
#define D3GDR4R4S0 -R*t*(-0.03125/SQUARE(xal3m1) + 0.03125/SQUARE(xfe3m1) \
                   + 0.0078125/SQUARE(xal3tet) - 0.0078125/SQUARE(xfe3tet))
#define D3GDR4R4S1 -R*t*0.125/SQUARE(xmg2m1) + ((HX6X6S2)+(p-1.0)*(VX6X6S2))*2.0
#define D3GDR4R4DT R*(0.0625/xal3m1 + 0.0625/xfe3m1 + 0.25/xmg2m1 \
                   + 0.25/(xmg2m1+xfe2m1-xti4m1) - 0.25/(xmg2m1+xfe2m1) \
                   - 0.125/(1.0-xsi4tet) + 0.03125/xal3tet + 0.03125/xfe3tet \
                   + 1.0/xca2m2 + 1.0/xna1m2 + 0.25/(1.0-xmg2m1-xfe2m1-xna1m2) \
                   - 0.25/(1.0-xmg2m1-xfe2m1) - 1.0/(1.0-xna1m2)) - (SX6X6)*2.0
#define D3GDR4R4DP (VX6X6)*2.0 + (VX6X6X7)*r[5]*2.0 + (VX6X6S2)*s[1]*2.0

#define D3GDR4R5R5 R*t*(0.125/SQUARE(xmg2m1)+1.0/SQUARE(xca2m2) ) + \
		   2.0*((HX6X7X7)-t*(SX6X7X7)+(p-1.0)*(VX6X7X7))
#define D3GDR4R5S0 0.0
#define D3GDR4R5S1 R*t*0.125/SQUARE(xmg2m1) + \
                   ((HX6X7S2)-t*(SX6X7S2)+(p-1.0)*(VX6X7S2))
#define D3GDR4R5DT R*(-0.25/xmg2m1 + 1.0/xca2m2) - (SX6X7) - \
                   (SX6X7X7)*r[5]*2.0 - (SX6X7S2)*s[1]
#define D3GDR4R5DP (VX6X7) + (VX2X6X7)*r[0] + (VX3X6X7)*r[1] + \
                   (VX4X6X7)*r[2] + (VX5X6X7)*r[3] + (VX6X6X7)*r[4]*2.0 + \
                   (VX6X7X7)*r[5]*2.0 + (VX6X7S2)*s[1]

#define D3GDR4S0S0 -R*t*(0.0625/SQUARE(xal3m1) + 0.0625/SQUARE(xfe3m1) \
                   - 0.015625/SQUARE(xal3tet) - 0.015625/SQUARE(xfe3tet))
#define D3GDR4S0S1 0.0
#define D3GDR4S0DT R*(-0.125/xal3m1 + 0.125/xfe3m1 - 0.0625/xal3tet \
                   + 0.0625/xfe3tet) - (SX6S1)
#define D3GDR4S0DP (VX6S1)

#define D3GDR4S1S1 R*t*0.125/SQUARE(xmg2m1) + ((HX6S2S2)+(p-1.0)*(VX6S2S2))*2.0
#define D3GDR4S1DT -R*0.25/xmg2m1 - (SX6S2) - (SX6X7S2)*r[5]
#define D3GDR4S1DP (VX6S2) + (VX2X6S2)*r[0] + (VX3X6S2)*r[1] + \
                   (VX4X6S2)*r[2] + (VX5X6S2)*r[3] + (VX6X6S2)*r[4]*2.0 + \
                   (VX6X7S2)*r[5] + (VX6S2S2)*s[1]*2.0

#define D3GDR5R5R5 R*t*(-0.125/SQUARE(xmg2m1)+0.125/SQUARE(xfe2m1) \
		   +1.0/SQUARE(xca2m2)-0.125/SQUARE(xmg2m2) \
		   -0.125/SQUARE(xfe2m2) ) + 6.0*((HX7X7X7)- \
		   t*(SX7X7X7)+(p-1.0)*(VX7X7X7) )
#define D3GDR5R5S0 0.0
#define D3GDR5R5S1 -R*t*(0.125/SQUARE(xmg2m1) - 0.125/SQUARE(xfe2m1) \
                   - 0.125/SQUARE(xmg2m2) + 0.125/SQUARE(xfe2m2)) + \
                   ((HX7X7S2)-t*(SX7X7S2)+(p-1.0)*(VX7X7S2))*2.0
#define D3GDR5R5DT R*(0.25/xmg2m1 + 0.25/xfe2m1 + 1.0/xca2m2 + 0.25/xmg2m2 \
                   + 0.25/xfe2m2) - (SX7X7)*2.0 - (SX2X7X7)*r[0]*2.0 \
                   - (SX3X7X7)*r[1]*2.0 - (SX4X7X7)*r[2]*2.0 \
                   - (SX5X7X7)*r[3]*2.0 - (SX6X7X7)*r[4]*2.0 \
                   - (SX7X7X7)*r[5]*6.0 - (SX7X7S2)*s[1]*2.0
#define D3GDR5R5DP (VX7X7)*2.0 + (VX2X7X7)*r[0]*2.0 + (VX3X7X7)*r[1]*2.0 + \
                   (VX4X7X7)*r[2]*2.0 + (VX5X7X7)*r[3]*2.0 + \
                   (VX6X7X7)*r[4]*2.0 + (VX7X7X7)*r[5]*6.0 + \
                   (VX7X7S2)*s[1]*2.0

#define D3GDR5S0S0 0.0
#define D3GDR5S0S1 0.0
#define D3GDR5S0DT -(SX7S1)
#define D3GDR5S0DP (VX7S1)

#define D3GDR5S1S1 -R*t*(0.125/SQUARE(xmg2m1) - 0.125/SQUARE(xfe2m1) \
                   + 0.125/SQUARE(xmg2m2) + 0.125/SQUARE(xfe2m2)) + \
                   ((HX7S2S2)+(p-1.0)*(VX7S2S2))*2.0
#define D3GDR5S1DT R*(0.25/xmg2m1+0.25/xfe2m1-0.25/xmg2m2+0.25/xfe2m2) \
                   - (SX7S2) - (SX3X7S2)*r[1] - (SX4X7S2)*r[2] - \
                   (SX5X7S2)*r[3] - (SX6X7S2)*r[4] - (SX7X7S2)*r[5]*2.0
#define D3GDR5S1DP (VX7S2) + (VX2X7S2)*r[0] + (VX3X7S2)*r[1] + \
                   (VX4X7S2)*r[2] + (VX5X7S2)*r[3] + (VX6X7S2)*r[4] + \
                   (VX7X7S2)*r[5]*2.0 + (VX7S2S2)*s[1]*2.0

#define D3GDS0S0S0 -R*t*(0.125/SQUARE(xfe3m1) - 0.125/SQUARE(xal3m1) \
                   + 0.03125/SQUARE(xal3tet) - 0.03125/SQUARE(xfe3tet))
#define D3GDS0S0S1 0.0
#define D3GDS0S0DT R*(0.25/xfe3m1+0.25/xal3m1+0.125/xal3tet+0.125/xfe3tet) \
                   - (SS1S1)*2.0
#define D3GDS0S0DP (VS1S1)*2.0
#define D3GDS0S1S1 0.0
#define D3GDS0S1DT -(SS1S2)
#define D3GDS0S1DP (VS1S2)
#define D3GDS0DT2  0.0
#define D3GDS0DTDP 0.0
#define D3GDS0DP2  0.0

#define D3GDS1S1S1 -R*t*(0.125/SQUARE(xmg2m1) - 0.125/SQUARE(xfe2m1) \
                   + 0.125/SQUARE(xfe2m2) - 0.125/SQUARE(xmg2m2))
#define D3GDS1S1DT R*(0.25/xmg2m1+0.25/xfe2m1+0.25/xfe2m2+0.25/xmg2m2) \
                   - (SS2S2)*2.0
#define D3GDS1S1DP (VS2S2)*2.0 + (VX2S2S2)*r[0]*2.0 + (VX3S2S2)*r[1]*2.0 + \
                   (VX4S2S2)*r[2]*2.0 + (VX5S2S2)*r[3]*2.0 + \
                   (VX6S2S2)*r[4]*2.0 + (VX7S2S2)*r[5]*2.0
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
#define D3GDR2DT2  0.0
#define D3GDR2DTDP 0.0
#define D3GDR2DP2  0.0
#define D3GDR3DT2  0.0
#define D3GDR3DTDP 0.0
#define D3GDR3DP2  0.0
#define D3GDR4DT2  0.0
#define D3GDR4DTDP 0.0
#define D3GDR4DP2  0.0
#define D3GDR5DT2  0.0
#define D3GDR5DTDP 0.0
#define D3GDR5DP2  0.0

/*
 *=============================================================================
 * Macros for automatic array initialization
 */
#define fillD2GDR2 \
 d2gdr2[0][0] = D2GDR0R0;     d2gdr2[0][1] = D2GDR0R1; \
 d2gdr2[0][2] = D2GDR0R2;     d2gdr2[0][3] = D2GDR0R3; \
 d2gdr2[0][4] = D2GDR0R4;     d2gdr2[0][5] = D2GDR0R5; \
 d2gdr2[1][0] = d2gdr2[0][1]; d2gdr2[1][1] = D2GDR1R1; \
 d2gdr2[1][2] = D2GDR1R2;     d2gdr2[1][3] = D2GDR1R3; \
 d2gdr2[1][4] = D2GDR1R4;     d2gdr2[1][5] = D2GDR1R5; \
 d2gdr2[2][0] = d2gdr2[0][2]; d2gdr2[2][1] = d2gdr2[1][2]; \
 d2gdr2[2][2] = D2GDR2R2;     d2gdr2[2][3] = D2GDR2R3; \
 d2gdr2[2][4] = D2GDR2R4;     d2gdr2[2][5] = D2GDR2R5; \
 d2gdr2[3][0] = d2gdr2[0][3]; d2gdr2[3][1] = d2gdr2[1][3]; \
 d2gdr2[3][2] = d2gdr2[2][3]; d2gdr2[3][3] = D2GDR3R3; \
 d2gdr2[3][4] = D2GDR3R4;     d2gdr2[3][5] = D2GDR3R5; \
 d2gdr2[4][0] = d2gdr2[0][4]; d2gdr2[4][1] = d2gdr2[1][4]; \
 d2gdr2[4][2] = d2gdr2[2][4]; d2gdr2[4][3] = d2gdr2[3][4]; \
 d2gdr2[4][4] = D2GDR4R4;     d2gdr2[4][5] = D2GDR4R5; \
 d2gdr2[5][0] = d2gdr2[0][5]; d2gdr2[5][1] = d2gdr2[1][5]; \
 d2gdr2[5][2] = d2gdr2[2][5]; d2gdr2[5][3] = d2gdr2[3][5]; \
 d2gdr2[5][4] = d2gdr2[4][5]; d2gdr2[5][5] = D2GDR5R5;

#define fillD2GDRDS \
 d2gdrds[0][0] = D2GDR0S0; d2gdrds[0][1] = D2GDR0S1; \
 d2gdrds[1][0] = D2GDR1S0; d2gdrds[1][1] = D2GDR1S1; \
 d2gdrds[2][0] = D2GDR2S0; d2gdrds[2][1] = D2GDR2S1; \
 d2gdrds[3][0] = D2GDR3S0; d2gdrds[3][1] = D2GDR3S1; \
 d2gdrds[4][0] = D2GDR4S0; d2gdrds[4][1] = D2GDR4S1; \
 d2gdrds[5][0] = D2GDR5S0; d2gdrds[5][1] = D2GDR5S1;

#define fillD2GDRDT \
 d2gdrdt[0] = D2GDR0DT; d2gdrdt[1] = D2GDR1DT; d2gdrdt[2] = D2GDR2DT; \
 d2gdrdt[3] = D2GDR3DT; d2gdrdt[4] = D2GDR4DT; d2gdrdt[5] = D2GDR5DT;

#define fillD2GDRDP \
 d2gdrdp[0] = D2GDR0DP; d2gdrdp[1] = D2GDR1DP; d2gdrdp[2] = D2GDR2DP; \
 d2gdrdp[3] = D2GDR3DP; d2gdrdp[4] = D2GDR4DP; d2gdrdp[5] = D2GDR5DP;

#define fillD2GDS2 \
 d2gds2[0][0] = D2GDS0S0;     d2gds2[0][1] = D2GDS0S1; \
 d2gds2[1][0] = d2gds2[0][1]; d2gds2[1][1] = D2GDS1S1;

#define fillD2GDSDT \
 d2gdsdt[0] = D2GDS0DT;  d2gdsdt[1] = D2GDS1DT;

#define fillD2GDSDP \
 d2gdsdp[0] = D2GDS0DP;  d2gdsdp[1] = D2GDS1DP;

#define fillD3GDR3 \
 d3gdr3[0][0][0] = D3GDR0R0R0;		d3gdr3[0][0][1] = D3GDR0R0R1; \
 d3gdr3[0][0][2] = D3GDR0R0R2;		d3gdr3[0][0][3] = D3GDR0R0R3; \
 d3gdr3[0][0][4] = D3GDR0R0R4;		d3gdr3[0][0][5] = D3GDR0R0R5; \
 d3gdr3[0][1][0] = d3gdr3[0][0][1];	d3gdr3[0][1][1] = D3GDR0R1R1; \
 d3gdr3[0][1][2] = D3GDR0R1R2;		d3gdr3[0][1][3] = D3GDR0R1R3; \
 d3gdr3[0][1][4] = D3GDR0R1R4;		d3gdr3[0][1][5] = D3GDR0R1R5; \
 d3gdr3[0][2][0] = d3gdr3[0][0][2];	d3gdr3[0][2][1] = d3gdr3[0][1][2]; \
 d3gdr3[0][2][2] = D3GDR0R2R2;		d3gdr3[0][2][3] = D3GDR0R2R3; \
 d3gdr3[0][2][4] = D3GDR0R2R4;		d3gdr3[0][2][5] = D3GDR0R2R5; \
 d3gdr3[0][3][0] = d3gdr3[0][0][3];	d3gdr3[0][3][1] = d3gdr3[0][1][3]; \
 d3gdr3[0][3][2] = d3gdr3[0][2][3];	d3gdr3[0][3][3] = D3GDR0R3R3; \
 d3gdr3[0][3][4] = D3GDR0R3R4;		d3gdr3[0][3][5] = D3GDR0R3R5; \
 d3gdr3[0][4][0] = d3gdr3[0][0][4];	d3gdr3[0][4][1] = d3gdr3[0][1][4]; \
 d3gdr3[0][4][2] = d3gdr3[0][2][4];	d3gdr3[0][4][3] = d3gdr3[0][3][4]; \
 d3gdr3[0][4][4] = D3GDR0R4R4;		d3gdr3[0][4][5] = D3GDR0R4R5; \
 d3gdr3[0][5][0] = d3gdr3[0][0][5];	d3gdr3[0][5][1] = d3gdr3[0][1][5]; \
 d3gdr3[0][5][2] = d3gdr3[0][2][5];	d3gdr3[0][5][3] = d3gdr3[0][3][5]; \
 d3gdr3[0][5][4] = d3gdr3[0][4][5];	d3gdr3[0][5][5] = D3GDR0R5R5; \
 d3gdr3[1][0][0] = d3gdr3[0][0][1];	d3gdr3[1][0][1] = d3gdr3[0][1][1]; \
 d3gdr3[1][0][2] = d3gdr3[0][1][2];	d3gdr3[1][0][3] = d3gdr3[0][1][3]; \
 d3gdr3[1][0][4] = d3gdr3[0][1][4];	d3gdr3[1][0][5] = d3gdr3[0][1][5]; \
 d3gdr3[1][1][0] = d3gdr3[0][1][1];	d3gdr3[1][1][1] = D3GDR1R1R1; \
 d3gdr3[1][1][2] = D3GDR1R1R2;		d3gdr3[1][1][3] = D3GDR1R1R3; \
 d3gdr3[1][1][4] = D3GDR1R1R4;		d3gdr3[1][1][5] = D3GDR1R1R5; \
 d3gdr3[1][2][0] = d3gdr3[0][1][2];	d3gdr3[1][2][1] = d3gdr3[1][1][2]; \
 d3gdr3[1][2][2] = D3GDR1R2R2;		d3gdr3[1][2][3] = D3GDR1R2R3; \
 d3gdr3[1][2][4] = D3GDR1R2R4;		d3gdr3[1][2][5] = D3GDR1R2R5; \
 d3gdr3[1][3][0] = d3gdr3[0][1][3];	d3gdr3[1][3][1] = d3gdr3[1][1][3]; \
 d3gdr3[1][3][2] = d3gdr3[1][2][3];	d3gdr3[1][3][3] = D3GDR1R3R3; \
 d3gdr3[1][3][4] = D3GDR1R3R4;		d3gdr3[1][3][5] = D3GDR1R3R5; \
 d3gdr3[1][4][0] = d3gdr3[0][1][4];	d3gdr3[1][4][1] = d3gdr3[1][1][4]; \
 d3gdr3[1][4][2] = d3gdr3[1][2][4];	d3gdr3[1][4][3] = d3gdr3[1][3][4]; \
 d3gdr3[1][4][4] = D3GDR1R4R4;		d3gdr3[1][4][5] = D3GDR1R4R5; \
 d3gdr3[1][5][0] = d3gdr3[0][1][5];	d3gdr3[1][5][1] = d3gdr3[1][1][5]; \
 d3gdr3[1][5][2] = d3gdr3[1][2][5];	d3gdr3[1][5][3] = d3gdr3[1][3][5]; \
 d3gdr3[1][5][4] = d3gdr3[1][4][5];	d3gdr3[1][5][5] = D3GDR1R5R5; \
 d3gdr3[2][0][0] = d3gdr3[0][0][2];	d3gdr3[2][0][1] = d3gdr3[0][1][2]; \
 d3gdr3[2][0][2] = d3gdr3[0][2][2];	d3gdr3[2][0][3] = d3gdr3[0][2][3]; \
 d3gdr3[2][0][4] = d3gdr3[0][2][4];	d3gdr3[2][0][5] = d3gdr3[0][2][5]; \
 d3gdr3[2][1][0] = d3gdr3[0][1][2];	d3gdr3[2][1][1] = d3gdr3[1][1][2]; \
 d3gdr3[2][1][2] = d3gdr3[1][2][2];	d3gdr3[2][1][3] = d3gdr3[1][2][3]; \
 d3gdr3[2][1][4] = d3gdr3[1][2][4];	d3gdr3[2][1][5] = d3gdr3[1][2][5]; \
 d3gdr3[2][2][0] = d3gdr3[0][2][2];	d3gdr3[2][2][1] = d3gdr3[1][2][2]; \
 d3gdr3[2][2][2] = D3GDR2R2R2;		d3gdr3[2][2][3] = D3GDR2R2R3; \
 d3gdr3[2][2][4] = D3GDR2R2R4;		d3gdr3[2][2][5] = D3GDR2R2R5; \
 d3gdr3[2][3][0] = d3gdr3[0][2][3];	d3gdr3[2][3][1] = d3gdr3[1][2][3]; \
 d3gdr3[2][3][2] = d3gdr3[2][2][3];	d3gdr3[2][3][3] = D3GDR2R3R3; \
 d3gdr3[2][3][4] = D3GDR2R3R4;		d3gdr3[2][3][5] = D3GDR2R3R5; \
 d3gdr3[2][4][0] = d3gdr3[0][2][4];	d3gdr3[2][4][1] = d3gdr3[1][2][4]; \
 d3gdr3[2][4][2] = d3gdr3[2][2][4];	d3gdr3[2][4][3] = d3gdr3[2][3][4]; \
 d3gdr3[2][4][4] = D3GDR2R4R4;		d3gdr3[2][4][5] = D3GDR2R4R5; \
 d3gdr3[2][5][0] = d3gdr3[0][2][5];	d3gdr3[2][5][1] = d3gdr3[1][2][5]; \
 d3gdr3[2][5][2] = d3gdr3[2][2][5];	d3gdr3[2][5][3] = d3gdr3[2][3][5]; \
 d3gdr3[2][5][4] = d3gdr3[2][4][5];	d3gdr3[2][5][5] = D3GDR2R5R5; \
 d3gdr3[3][0][0] = d3gdr3[0][0][3];	d3gdr3[3][0][1] = d3gdr3[0][1][3]; \
 d3gdr3[3][0][2] = d3gdr3[0][2][3];	d3gdr3[3][0][3] = d3gdr3[0][3][3]; \
 d3gdr3[3][0][4] = d3gdr3[0][3][4];	d3gdr3[3][0][5] = d3gdr3[0][3][5]; \
 d3gdr3[3][1][0] = d3gdr3[0][1][3];	d3gdr3[3][1][1] = d3gdr3[1][1][3]; \
 d3gdr3[3][1][2] = d3gdr3[1][2][3];	d3gdr3[3][1][3] = d3gdr3[1][3][3]; \
 d3gdr3[3][1][4] = d3gdr3[1][3][4];	d3gdr3[3][1][5] = d3gdr3[1][3][5]; \
 d3gdr3[3][2][0] = d3gdr3[0][2][3];	d3gdr3[3][2][1] = d3gdr3[1][2][3]; \
 d3gdr3[3][2][2] = d3gdr3[2][2][3];	d3gdr3[3][2][3] = d3gdr3[2][3][3]; \
 d3gdr3[3][2][4] = d3gdr3[2][3][4];	d3gdr3[3][2][5] = d3gdr3[2][3][5]; \
 d3gdr3[3][3][0] = d3gdr3[0][3][3];	d3gdr3[3][3][1] = d3gdr3[1][3][3]; \
 d3gdr3[3][3][2] = d3gdr3[2][3][3];	d3gdr3[3][3][3] = D3GDR3R3R3; \
 d3gdr3[3][3][4] = D3GDR3R3R4;		d3gdr3[3][3][5] = D3GDR3R3R5; \
 d3gdr3[3][4][0] = d3gdr3[0][3][4];	d3gdr3[3][4][1] = d3gdr3[1][3][4]; \
 d3gdr3[3][4][2] = d3gdr3[2][3][4];	d3gdr3[3][4][3] = d3gdr3[3][3][4]; \
 d3gdr3[3][4][4] = D3GDR3R4R4;		d3gdr3[3][4][5] = D3GDR3R4R5; \
 d3gdr3[3][5][0] = d3gdr3[0][3][5];	d3gdr3[3][5][1] = d3gdr3[1][3][5]; \
 d3gdr3[3][5][2] = d3gdr3[2][3][5];	d3gdr3[3][5][3] = d3gdr3[3][3][5]; \
 d3gdr3[3][5][4] = d3gdr3[3][4][5];	d3gdr3[3][5][5] = D3GDR3R5R5; \
 d3gdr3[4][0][0] = d3gdr3[0][0][4];	d3gdr3[4][0][1] = d3gdr3[0][1][4]; \
 d3gdr3[4][0][2] = d3gdr3[0][2][4];	d3gdr3[4][0][3] = d3gdr3[0][3][4]; \
 d3gdr3[4][0][4] = d3gdr3[0][4][4];	d3gdr3[4][0][5] = d3gdr3[0][4][5]; \
 d3gdr3[4][1][0] = d3gdr3[0][1][4];	d3gdr3[4][1][1] = d3gdr3[1][1][4]; \
 d3gdr3[4][1][2] = d3gdr3[1][2][4];	d3gdr3[4][1][3] = d3gdr3[1][3][4]; \
 d3gdr3[4][1][4] = d3gdr3[1][4][4];	d3gdr3[4][1][5] = d3gdr3[1][4][5]; \
 d3gdr3[4][2][0] = d3gdr3[0][2][4];	d3gdr3[4][2][1] = d3gdr3[1][2][4]; \
 d3gdr3[4][2][2] = d3gdr3[2][2][4];	d3gdr3[4][2][3] = d3gdr3[2][3][4]; \
 d3gdr3[4][2][4] = d3gdr3[2][4][4];	d3gdr3[4][2][5] = d3gdr3[2][4][5]; \
 d3gdr3[4][3][0] = d3gdr3[0][3][4];	d3gdr3[4][3][1] = d3gdr3[1][3][4]; \
 d3gdr3[4][3][2] = d3gdr3[2][3][4];	d3gdr3[4][3][3] = d3gdr3[3][3][4]; \
 d3gdr3[4][3][4] = d3gdr3[3][4][4];	d3gdr3[4][3][5] = d3gdr3[3][4][5]; \
 d3gdr3[4][4][0] = d3gdr3[0][4][4];	d3gdr3[4][4][1] = d3gdr3[1][4][4]; \
 d3gdr3[4][4][2] = d3gdr3[2][4][4];	d3gdr3[4][4][3] = d3gdr3[3][4][4]; \
 d3gdr3[4][4][4] = D3GDR4R4R4;		d3gdr3[4][4][5] = D3GDR4R4R5; \
 d3gdr3[4][5][0] = d3gdr3[0][4][5];	d3gdr3[4][5][1] = d3gdr3[1][4][5]; \
 d3gdr3[4][5][2] = d3gdr3[2][4][5];	d3gdr3[4][5][3] = d3gdr3[3][4][5]; \
 d3gdr3[4][5][4] = d3gdr3[4][4][5];	d3gdr3[4][5][5] = D3GDR4R5R5; \
 d3gdr3[5][0][0] = d3gdr3[0][0][5];	d3gdr3[5][0][1] = d3gdr3[0][1][5]; \
 d3gdr3[5][0][2] = d3gdr3[0][2][5];	d3gdr3[5][0][3] = d3gdr3[0][3][5]; \
 d3gdr3[5][0][4] = d3gdr3[0][4][5];	d3gdr3[5][0][5] = d3gdr3[0][5][5]; \
 d3gdr3[5][1][0] = d3gdr3[0][1][5];	d3gdr3[5][1][1] = d3gdr3[1][1][5]; \
 d3gdr3[5][1][2] = d3gdr3[1][2][5];	d3gdr3[5][1][3] = d3gdr3[1][3][5]; \
 d3gdr3[5][1][4] = d3gdr3[1][4][5];	d3gdr3[5][1][5] = d3gdr3[1][5][5]; \
 d3gdr3[5][2][0] = d3gdr3[0][2][5];	d3gdr3[5][2][1] = d3gdr3[1][2][5]; \
 d3gdr3[5][2][2] = d3gdr3[2][2][5];	d3gdr3[5][2][3] = d3gdr3[2][3][5]; \
 d3gdr3[5][2][4] = d3gdr3[2][4][5];	d3gdr3[5][2][5] = d3gdr3[2][5][5]; \
 d3gdr3[5][3][0] = d3gdr3[0][3][5];	d3gdr3[5][3][1] = d3gdr3[1][3][5]; \
 d3gdr3[5][3][2] = d3gdr3[2][3][5];	d3gdr3[5][3][3] = d3gdr3[3][3][5]; \
 d3gdr3[5][3][4] = d3gdr3[3][4][5];	d3gdr3[5][3][5] = d3gdr3[3][5][5]; \
 d3gdr3[5][4][0] = d3gdr3[0][4][5];	d3gdr3[5][4][1] = d3gdr3[1][4][5]; \
 d3gdr3[5][4][2] = d3gdr3[2][4][5];	d3gdr3[5][4][3] = d3gdr3[3][4][5]; \
 d3gdr3[5][4][4] = d3gdr3[4][4][5];	d3gdr3[5][4][5] = d3gdr3[4][5][5]; \
 d3gdr3[5][5][0] = d3gdr3[0][5][5];	d3gdr3[5][5][1] = d3gdr3[1][5][5]; \
 d3gdr3[5][5][2] = d3gdr3[2][5][5];	d3gdr3[5][5][3] = d3gdr3[3][5][5]; \
 d3gdr3[5][5][4] = d3gdr3[4][5][5];	d3gdr3[5][5][5] = D3GDR5R5R5; \

#define fillD3GDR2DS \
 d3gdr2ds[0][0][0] = D3GDR0R0S0;        d3gdr2ds[0][0][1] = D3GDR0R0S1; \
 d3gdr2ds[0][1][0] = D3GDR0R1S0;        d3gdr2ds[0][1][1] = D3GDR0R1S1; \
 d3gdr2ds[0][2][0] = D3GDR0R2S0;        d3gdr2ds[0][2][1] = D3GDR0R2S1; \
 d3gdr2ds[0][3][0] = D3GDR0R3S0;        d3gdr2ds[0][3][1] = D3GDR0R3S1; \
 d3gdr2ds[0][4][0] = D3GDR0R4S0;        d3gdr2ds[0][4][1] = D3GDR0R4S1; \
 d3gdr2ds[0][5][0] = D3GDR0R5S0;        d3gdr2ds[0][5][1] = D3GDR0R5S1; \
 d3gdr2ds[1][0][0] = d3gdr2ds[0][1][0]; d3gdr2ds[1][0][1] = d3gdr2ds[0][1][1]; \
 d3gdr2ds[1][1][0] = D3GDR1R1S0;        d3gdr2ds[1][1][1] = D3GDR1R1S1; \
 d3gdr2ds[1][2][0] = D3GDR1R2S0;        d3gdr2ds[1][2][1] = D3GDR1R2S1; \
 d3gdr2ds[1][3][0] = D3GDR1R3S0;        d3gdr2ds[1][3][1] = D3GDR1R3S1; \
 d3gdr2ds[1][4][0] = D3GDR1R4S0;        d3gdr2ds[1][4][1] = D3GDR1R4S1; \
 d3gdr2ds[1][5][0] = D3GDR1R5S0;        d3gdr2ds[1][5][1] = D3GDR1R5S1; \
 d3gdr2ds[2][0][0] = d3gdr2ds[0][2][0]; d3gdr2ds[2][0][1] = d3gdr2ds[0][2][1]; \
 d3gdr2ds[2][1][0] = d3gdr2ds[1][2][0]; d3gdr2ds[2][1][1] = d3gdr2ds[1][2][1]; \
 d3gdr2ds[2][2][0] = D3GDR2R2S0;        d3gdr2ds[2][2][1] = D3GDR2R2S1; \
 d3gdr2ds[2][3][0] = D3GDR2R3S0;        d3gdr2ds[2][3][1] = D3GDR2R3S1; \
 d3gdr2ds[2][4][0] = D3GDR2R4S0;        d3gdr2ds[2][4][1] = D3GDR2R4S1; \
 d3gdr2ds[2][5][0] = D3GDR2R5S0;        d3gdr2ds[2][5][1] = D3GDR2R5S1; \
 d3gdr2ds[3][0][0] = d3gdr2ds[0][3][0]; d3gdr2ds[3][0][1] = d3gdr2ds[0][3][1]; \
 d3gdr2ds[3][1][0] = d3gdr2ds[1][3][0]; d3gdr2ds[3][1][1] = d3gdr2ds[1][3][1]; \
 d3gdr2ds[3][2][0] = d3gdr2ds[2][3][0]; d3gdr2ds[3][2][1] = d3gdr2ds[2][3][1]; \
 d3gdr2ds[3][3][0] = D3GDR3R3S0;        d3gdr2ds[3][3][1] = D3GDR3R3S1; \
 d3gdr2ds[3][4][0] = D3GDR3R4S0;        d3gdr2ds[3][4][1] = D3GDR3R4S1; \
 d3gdr2ds[3][5][0] = D3GDR3R5S0;        d3gdr2ds[3][5][1] = D3GDR3R5S1; \
 d3gdr2ds[4][0][0] = d3gdr2ds[0][4][0]; d3gdr2ds[4][0][1] = d3gdr2ds[0][4][1]; \
 d3gdr2ds[4][1][0] = d3gdr2ds[1][4][0]; d3gdr2ds[4][1][1] = d3gdr2ds[1][4][1]; \
 d3gdr2ds[4][2][0] = d3gdr2ds[2][4][0]; d3gdr2ds[4][2][1] = d3gdr2ds[2][4][1]; \
 d3gdr2ds[4][3][0] = d3gdr2ds[3][4][0]; d3gdr2ds[4][3][1] = d3gdr2ds[3][4][1]; \
 d3gdr2ds[4][4][0] = D3GDR4R4S0;        d3gdr2ds[4][4][1] = D3GDR4R4S1; \
 d3gdr2ds[4][5][0] = D3GDR4R5S0;        d3gdr2ds[4][5][1] = D3GDR4R5S1; \
 d3gdr2ds[5][0][0] = d3gdr2ds[0][5][0]; d3gdr2ds[5][0][1] = d3gdr2ds[0][5][1]; \
 d3gdr2ds[5][1][0] = d3gdr2ds[1][5][0]; d3gdr2ds[5][1][1] = d3gdr2ds[1][5][1]; \
 d3gdr2ds[5][2][0] = d3gdr2ds[2][5][0]; d3gdr2ds[5][2][1] = d3gdr2ds[2][5][1]; \
 d3gdr2ds[5][3][0] = d3gdr2ds[3][5][0]; d3gdr2ds[5][3][1] = d3gdr2ds[3][5][1]; \
 d3gdr2ds[5][4][0] = d3gdr2ds[4][5][0]; d3gdr2ds[5][4][1] = d3gdr2ds[4][5][1]; \
 d3gdr2ds[5][5][0] = D3GDR5R5S0;        d3gdr2ds[5][5][1] = D3GDR5R5S1;

#define fillD3GDR2DT \
 d3gdr2dt[0][0] = D3GDR0R0DT;     d3gdr2dt[0][1] = D3GDR0R1DT; \
 d3gdr2dt[0][2] = D3GDR0R2DT;     d3gdr2dt[0][3] = D3GDR0R3DT; \
 d3gdr2dt[0][4] = D3GDR0R4DT;     d3gdr2dt[0][5] = D3GDR0R5DT; \
 d3gdr2dt[1][0] = d3gdr2dt[0][1]; d3gdr2dt[1][1] = D3GDR1R1DT; \
 d3gdr2dt[1][2] = D3GDR1R2DT;     d3gdr2dt[1][3] = D3GDR1R3DT; \
 d3gdr2dt[1][4] = D3GDR1R4DT;     d3gdr2dt[1][5] = D3GDR1R5DT; \
 d3gdr2dt[2][0] = d3gdr2dt[0][2]; d3gdr2dt[2][1] = d3gdr2dt[1][2]; \
 d3gdr2dt[2][2] = D3GDR2R2DT;     d3gdr2dt[2][3] = D3GDR2R3DT; \
 d3gdr2dt[2][4] = D3GDR2R4DT;     d3gdr2dt[2][5] = D3GDR2R5DT; \
 d3gdr2dt[3][0] = d3gdr2dt[0][3]; d3gdr2dt[3][1] = d3gdr2dt[1][3]; \
 d3gdr2dt[3][2] = d3gdr2dt[2][3]; d3gdr2dt[3][3] = D3GDR3R3DT; \
 d3gdr2dt[3][4] = D3GDR3R4DT;     d3gdr2dt[3][5] = D3GDR3R5DT; \
 d3gdr2dt[4][0] = d3gdr2dt[0][4]; d3gdr2dt[4][1] = d3gdr2dt[1][4]; \
 d3gdr2dt[4][2] = d3gdr2dt[2][4]; d3gdr2dt[4][3] = d3gdr2dt[3][4]; \
 d3gdr2dt[4][4] = D3GDR4R4DT;     d3gdr2dt[4][5] = D3GDR4R5DT; \
 d3gdr2dt[5][0] = d3gdr2dt[0][5]; d3gdr2dt[5][1] = d3gdr2dt[1][5]; \
 d3gdr2dt[5][2] = d3gdr2dt[2][5]; d3gdr2dt[5][3] = d3gdr2dt[3][5]; \
 d3gdr2dt[5][4] = d3gdr2dt[4][5]; d3gdr2dt[5][5] = D3GDR5R5DT;

#define fillD3GDR2DP \
 d3gdr2dp[0][0] = D3GDR0R0DP;     d3gdr2dp[0][1] = D3GDR0R1DP; \
 d3gdr2dp[0][2] = D3GDR0R2DP;     d3gdr2dp[0][3] = D3GDR0R3DP; \
 d3gdr2dp[0][4] = D3GDR0R4DP;     d3gdr2dp[0][5] = D3GDR0R5DP; \
 d3gdr2dp[1][0] = d3gdr2dp[0][1]; d3gdr2dp[1][1] = D3GDR1R1DP; \
 d3gdr2dp[1][2] = D3GDR1R2DP;     d3gdr2dp[1][3] = D3GDR1R3DP; \
 d3gdr2dp[1][4] = D3GDR1R4DP;     d3gdr2dp[1][5] = D3GDR1R5DP; \
 d3gdr2dp[2][0] = d3gdr2dp[0][2]; d3gdr2dp[2][1] = d3gdr2dp[1][2]; \
 d3gdr2dp[2][2] = D3GDR2R2DP;     d3gdr2dp[2][3] = D3GDR2R3DP; \
 d3gdr2dp[2][4] = D3GDR2R4DP;     d3gdr2dp[2][5] = D3GDR2R5DP; \
 d3gdr2dp[3][0] = d3gdr2dp[0][3]; d3gdr2dp[3][1] = d3gdr2dp[1][3]; \
 d3gdr2dp[3][2] = d3gdr2dp[2][3]; d3gdr2dp[3][3] = D3GDR3R3DP; \
 d3gdr2dp[3][4] = D3GDR3R4DP;     d3gdr2dp[3][5] = D3GDR3R5DP; \
 d3gdr2dp[4][0] = d3gdr2dp[0][4]; d3gdr2dp[4][1] = d3gdr2dp[1][4]; \
 d3gdr2dp[4][2] = d3gdr2dp[2][4]; d3gdr2dp[4][3] = d3gdr2dp[3][4]; \
 d3gdr2dp[4][4] = D3GDR4R4DP;     d3gdr2dp[4][5] = D3GDR4R5DP; \
 d3gdr2dp[5][0] = d3gdr2dp[0][5]; d3gdr2dp[5][1] = d3gdr2dp[1][5]; \
 d3gdr2dp[5][2] = d3gdr2dp[2][5]; d3gdr2dp[5][3] = d3gdr2dp[3][5]; \
 d3gdr2dp[5][4] = d3gdr2dp[4][5]; d3gdr2dp[5][5] = D3GDR5R5DP;

#define fillD3GDRDS2 \
 d3gdrds2[0][0][0] = D3GDR0S0S0;        d3gdrds2[0][0][1] = D3GDR0S0S1; \
 d3gdrds2[0][1][0] = d3gdrds2[0][0][1]; d3gdrds2[0][1][1] = D3GDR0S1S1; \
 d3gdrds2[1][0][0] = D3GDR1S0S0;        d3gdrds2[1][0][1] = D3GDR1S0S1; \
 d3gdrds2[1][1][0] = d3gdrds2[1][0][1]; d3gdrds2[1][1][1] = D3GDR1S1S1; \
 d3gdrds2[2][0][0] = D3GDR2S0S0;        d3gdrds2[2][0][1] = D3GDR2S0S1; \
 d3gdrds2[2][1][0] = d3gdrds2[2][0][1]; d3gdrds2[2][1][1] = D3GDR2S1S1; \
 d3gdrds2[3][0][0] = D3GDR3S0S0;        d3gdrds2[3][0][1] = D3GDR3S0S1; \
 d3gdrds2[3][1][0] = d3gdrds2[3][0][1]; d3gdrds2[3][1][1] = D3GDR3S1S1; \
 d3gdrds2[4][0][0] = D3GDR4S0S0;        d3gdrds2[4][0][1] = D3GDR4S0S1; \
 d3gdrds2[4][1][0] = d3gdrds2[4][0][1]; d3gdrds2[4][1][1] = D3GDR4S1S1; \
 d3gdrds2[5][0][0] = D3GDR5S0S0;        d3gdrds2[5][0][1] = D3GDR5S0S1; \
 d3gdrds2[5][1][0] = d3gdrds2[5][0][1]; d3gdrds2[5][1][1] = D3GDR5S1S1;

#define fillD3GDRDSDT \
 d3gdrdsdt[0][0] = D3GDR0S0DT;  d3gdrdsdt[0][1] = D3GDR0S1DT; \
 d3gdrdsdt[1][0] = D3GDR1S0DT;  d3gdrdsdt[1][1] = D3GDR1S1DT; \
 d3gdrdsdt[2][0] = D3GDR2S0DT;  d3gdrdsdt[2][1] = D3GDR2S1DT; \
 d3gdrdsdt[3][0] = D3GDR3S0DT;  d3gdrdsdt[3][1] = D3GDR3S1DT; \
 d3gdrdsdt[4][0] = D3GDR4S0DT;  d3gdrdsdt[4][1] = D3GDR4S1DT; \
 d3gdrdsdt[5][0] = D3GDR5S0DT;  d3gdrdsdt[5][1] = D3GDR5S1DT;

#define fillD3GDRDSDP \
 d3gdrdsdp[0][0] = D3GDR0S0DP; d3gdrdsdp[0][1] = D3GDR0S1DP; \
 d3gdrdsdp[1][0] = D3GDR1S0DP; d3gdrdsdp[1][1] = D3GDR1S1DP; \
 d3gdrdsdp[2][0] = D3GDR2S0DP; d3gdrdsdp[2][1] = D3GDR2S1DP; \
 d3gdrdsdp[3][0] = D3GDR3S0DP; d3gdrdsdp[3][1] = D3GDR3S1DP; \
 d3gdrdsdp[4][0] = D3GDR4S0DP; d3gdrdsdp[4][1] = D3GDR4S1DP; \
 d3gdrdsdp[5][0] = D3GDR5S0DP; d3gdrdsdp[5][1] = D3GDR5S1DP;

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
 d3gdrdt2[0] = D3GDR0DT2; d3gdrdt2[1] = D3GDR1DT2; \
 d3gdrdt2[2] = D3GDR2DT2; d3gdrdt2[3] = D3GDR3DT2; \
 d3gdrdt2[4] = D3GDR4DT2; d3gdrdt2[5] = D3GDR5DT2;

#define fillD3GDRDTDP \
 d3gdrdtdp[0] = D3GDR0DTDP; d3gdrdtdp[1] = D3GDR1DTDP; \
 d3gdrdtdp[2] = D3GDR2DTDP; d3gdrdtdp[3] = D3GDR3DTDP; \
 d3gdrdtdp[4] = D3GDR4DTDP; d3gdrdtdp[5] = D3GDR5DTDP;

#define fillD3GDRDP2 \
 d3gdrdp2[0] = D3GDR0DP2; d3gdrdp2[1] = D3GDR1DP2; \
 d3gdrdp2[2] = D3GDR2DP2; d3gdrdp2[3] = D3GDR3DP2; \
 d3gdrdp2[4] = D3GDR4DP2; d3gdrdp2[5] = D3GDR5DP2;

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
    gsl_matrix      *ptToD2gds2  = getPtToD2gds2();
    gsl_permutation *indexD2gds2 = getIndexD2gds2();
    int    structOld    = getStructOld();

    int i, j, iter=0, signum;

    GET_SITE_FRACTIONS

    /* look-up or compute the current ordering state */
    if ( (t != tOld)       || (p != pOld)    || (structOld != clino) ||
       (r[0] != rOld[0]) || (r[1] != rOld[1]) || (r[2] != rOld[2]) ||
       (r[3] != rOld[3]) || (r[4] != rOld[4]) || (r[5] != rOld[5]) ) {
        double dgds[NS], sNew[NS];
        double totAl, totCa, totFe2, totFe3, totMg, totNa, totTi, totSi;
        double totM1, totM2;
        gsl_vector_view vvToDgds = gsl_vector_view_array(dgds, (size_t) NS);

        totAl  = r[1] + r[3];
        totCa  = 1.0 - r[4] - r[5];
        totFe2 = r[0];
        totFe3 = r[2] + r[3];
        totMg  = 1.0 - r[0] - 0.5*r[1] - 0.5*r[2] - r[3] - 0.5*r[4] + r[5];
        totNa  = r[4];
        totTi  = 0.5*(r[1]+r[2]);
        totSi  = 2.0 - r[1] - r[2] - r[3] + r[4]/2.0;

             /* M1: Ti + extra Fe2+ and Mg not needed to fill M2 site */
        totM2  = totCa + totNa;
        totM1  = totTi + (totMg+totFe2-(1.0-totM2));

        /* Test here for whether to use ab initio initial guess or previous solution */
        if ( (fabs(t-tOld) > 5) || (fabs(p-pOld) > 100) || (structOld != clino) ||
                (fabs((r[0]-rOld[0])/((rOld[0] != 0.0) ? rOld[0] : 1.0)) > .02) ||
                (fabs((r[1]-rOld[1])/((rOld[1] != 0.0) ? rOld[1] : 1.0)) > .02) ||
                (fabs((r[2]-rOld[2])/((rOld[2] != 0.0) ? rOld[2] : 1.0)) > .02) ||
                (fabs((r[3]-rOld[3])/((rOld[3] != 0.0) ? rOld[3] : 1.0)) > .02) ||
                (fabs((r[4]-rOld[4])/((rOld[4] != 0.0) ? rOld[4] : 1.0)) > .02) ||
                (fabs((r[5]-rOld[5])/((rOld[5] != 0.0) ? rOld[5] : 1.0)) > .02) ) {
            sNew[0] = (totFe3+totAl != 0.0) ?
                            (1.0-totM1)*(totFe3-totAl)/(totFe3+totAl) : 0.0;
            sNew[1] = (totFe2+totMg != 0.0) ?
                            (1.0-totM2)*(totFe2-totMg)/(totFe2+totMg) : 0.0;
        } else {
            sNew[0] = sOld[0]; sNew[1] = sOld[1];
        }

        structOld = clino;
        for (i=0; i<NS; i++) sOld[i] = 2.0;

        xti4m1  = (r[1]+r[2])/2.0;
        xca2m2  = 1.0 - r[4] - r[5];
        xna1m2  = r[4];
        xsi4tet = (4.0-2.0*r[1]-2.0*r[2]-2.0*r[3]+r[4])/4.0;

        if (xti4m1  <= DBL_EPSILON) xti4m1  = DBL_EPSILON;
        if (xca2m2  <= DBL_EPSILON) xca2m2  = DBL_EPSILON;
        if (xna1m2  <= DBL_EPSILON) xna1m2  = DBL_EPSILON;
        if (xsi4tet <= DBL_EPSILON) xsi4tet = DBL_EPSILON;

        if (xti4m1  >= 1.0-DBL_EPSILON) xti4m1  = 1.0 - DBL_EPSILON;
        if (xna1m2  >= 1.0-DBL_EPSILON) xna1m2  = 1.0 - DBL_EPSILON;
        if (xsi4tet >= 1.0-DBL_EPSILON) xsi4tet = 1.0 - DBL_EPSILON;

        if ((totAl+totFe3) < 10.0*DBL_EPSILON) {
            xal3m1  = DBL_EPSILON;
            xfe3m1  = DBL_EPSILON;
            xal3tet = DBL_EPSILON;
            xfe3tet = DBL_EPSILON;

            sNew[0] = 0.0;
            sOld[0] = sNew[0];
        }

        if (totFe2 < 10.0*DBL_EPSILON &&
                totAl  < 10.0*DBL_EPSILON &&
                totFe3 < 10.0*DBL_EPSILON &&
                totNa  < 10.0*DBL_EPSILON &&
                totTi  < 10.0*DBL_EPSILON   ) { /* Di-En join */
            xfe2m1 = DBL_EPSILON;
            xmg2m1 = 1.0 - 3.0*DBL_EPSILON;

            xfe2m2 = DBL_EPSILON;
            xmg2m2 = 1.0 - xca2m2;

            sNew[1] = xfe2m2 - xmg2m2;
            sOld[1] = sNew[1];
        }

        if (totMg < 10.0*DBL_EPSILON &&
            totAl  < 10.0*DBL_EPSILON &&
            totFe3 < 10.0*DBL_EPSILON &&
            totNa  < 10.0*DBL_EPSILON &&
            totTi  < 10.0*DBL_EPSILON   ) { /* Hd-Fs join (PMA 05/19/16) */
            xmg2m1 = DBL_EPSILON;
            xfe2m1 = 1.0 - 3.0*DBL_EPSILON;

            xmg2m2 = DBL_EPSILON;
            xfe2m2 = 1.0 - xca2m2;

            sNew[1] = xfe2m2 - xmg2m2;
            sOld[1] = sNew[1];
        }

        while ( ((ABS(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON) ||
             (ABS(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON) ) && (iter < MAX_ITER)) {
            double s[NS], sCorr[NS];
            gsl_vector_view vvToSCorr = gsl_vector_view_array(sCorr, (size_t) NS);

            for (i=0; i<NS; i++) s[i] = sNew[i];

            xal3m1 = (2.0*r[3]+r[4]-2.0*s[0])/4.0;
            xfe2m1 = (2.0*r[0]-r[5]-s[1])/2.0;
            xfe3m1 = (2.0*r[3]+r[4]+2.0*s[0])/4.0;
            xmg2m1 = 1.0 - r[0] - r[3] - 0.5*(r[1]+r[2]+r[4]-r[5]-s[1]);

            xfe2m2 = (r[5]+s[1])/2.0;
            xmg2m2 = (r[5]-s[1])/2.0;

            xal3tet = (4.0*r[1]+2.0*r[3]-r[4]+2.0*s[0])/8.0;
            xfe3tet = (4.0*r[2]+2.0*r[3]-r[4]-2.0*s[0])/8.0;

            if (xal3m1  <= DBL_EPSILON) xal3m1  = DBL_EPSILON;
            if (xfe2m1  <= DBL_EPSILON) xfe2m1  = DBL_EPSILON;
            if (xfe3m1  <= DBL_EPSILON) xfe3m1  = DBL_EPSILON;
            if (xmg2m1  <= DBL_EPSILON) xmg2m1  = DBL_EPSILON;
            if (xfe2m2  <= DBL_EPSILON) xfe2m2  = DBL_EPSILON;
            if (xmg2m2  <= DBL_EPSILON) xmg2m2  = DBL_EPSILON;
            if (xal3tet <= DBL_EPSILON) xal3tet = DBL_EPSILON;
            if (xfe3tet <= DBL_EPSILON) xfe3tet = DBL_EPSILON;

            if (xmg2m1+xfe2m1 >= 1.0 || xmg2m1+xfe2m1+xna1m2 >= 1.0) {
                if (xmg2m1 > 2.0*DBL_EPSILON) xmg2m1 -= 2.0*DBL_EPSILON;
                if (xfe2m1 > 2.0*DBL_EPSILON) xfe2m1 -= 2.0*DBL_EPSILON;
            }

            if (xal3m1+xfe3m1+2.0*xti4m1 >= 1.0) { /* roundoff e.g. if TiO2=0.0 */
                xmg2m1 += 2.0*DBL_EPSILON;
                xfe2m1 += 2.0*DBL_EPSILON;
            }

            dgds[0] = (totFe3 != 0.0 && totAl != 0.0) ? DGDS0 : 0.0;
            dgds[1] = (totFe2 != 0.0 && totMg != 0.0) ? DGDS1 : 0.0;

            d2gds2[0][0] = (totFe3 != 0.0 && totAl != 0.0) ? D2GDS0S0 : 0.0; // was 1.0 with gaussj
            d2gds2[0][1] = (totFe3 != 0.0 && totAl != 0.0 &&
                                            totFe2 != 0.0 && totMg != 0.0) ? D2GDS0S1 : 0.0;
            d2gds2[1][0] = d2gds2[0][1];
            d2gds2[1][1] = (totFe2 != 0.0 && totMg != 0.0) ? D2GDS1S1 : 0.0; // was 1.0 with gaussj

            for (i=0; i<NS; i++) sOld[i] = s[i];

            /* original: gaussj(ptToD2gds2, NS, (double **) NULL, 0);
            if (totFe3 == 0.0 || totAl == 0.0) d2gds2[0][0] = 0.0;
            if (totFe2 == 0.0 || totMg == 0.0) d2gds2[1][1] = 0.0;
            for (i=0; i<NS; i++) {
            for(j=0; j<NS; j++) s[i] += - d2gds2[i][j]*dgds[j];
            } */

            gsl_matrix_scale(ptToD2gds2, -1.0);
            melts_LU_decomp(ptToD2gds2, indexD2gds2, &signum);

            melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToDgds.vector, &vvToSCorr.vector);
            for (i=0; i<NS; i++) s[i] += sCorr[i];

            if (totFe3+totAl != 0.0) {
                s[0] = MAX(s[0],0.5*( r[4]-2.0*r[3]-4.0*r[1]             )+DBL_EPSILON);
                s[0] = MAX(s[0],0.5*(-r[4]+2.0*r[3]         +4.0*r[2]-4.0)+DBL_EPSILON);
                s[0] = MAX(s[0],0.5*( r[4]+2.0*r[3]                  -4.0)+DBL_EPSILON);
                s[0] = MAX(s[0],0.5*(-r[4]-2.0*r[3]                      )+DBL_EPSILON);

                s[0] = MIN(s[0],0.5*( r[4]-2.0*r[3]-4.0*r[1]         +4.0)-DBL_EPSILON);
                s[0] = MIN(s[0],0.5*(-r[4]+2.0*r[3]         +4.0*r[2]    )-DBL_EPSILON);
                s[0] = MIN(s[0],0.5*( r[4]+2.0*r[3]                      )-DBL_EPSILON);
                s[0] = MIN(s[0],0.5*(-r[4]-2.0*r[3]                  +4.0)-DBL_EPSILON);
            }

            if (totFe2+totMg != 0.0) {
                s[1] = MAX(s[1],-r[5]                        +2.0*r[0]-2.0+DBL_EPSILON);
                s[1] = MAX(s[1],-r[5]+r[4]+2.0*r[3]+r[2]+r[1]+2.0*r[0]-2.0+DBL_EPSILON);
                s[1] = MAX(s[1],-r[5]                                     +DBL_EPSILON);
                s[1] = MAX(s[1], r[5]                                 -2.0+DBL_EPSILON);

                s[1] = MIN(s[1],-r[5]                        +2.0*r[0]    -DBL_EPSILON);
                s[1] = MIN(s[1],-r[5]+r[4]+2.0*r[3]+r[2]+r[1]+2.0*r[0]    -DBL_EPSILON);
                s[1] = MIN(s[1],-r[5]                                 +2.0-DBL_EPSILON);
                s[1] = MIN(s[1], r[5]                                     -DBL_EPSILON);
            }

            for (i=0; i<NS; i++) sNew[i] = s[i];
            iter++;
        }
        tOld = t;
        pOld = p;
        for (i=0; i<NR; i++) rOld[i] = r[i];

        if (fabs(totFe3) < sqrt(DBL_EPSILON) || fabs(totAl) < sqrt(DBL_EPSILON))
            dgds[0] = 0.0;
        if (fabs(totFe2) < sqrt(DBL_EPSILON) || fabs(totMg) < sqrt(DBL_EPSILON))
            dgds[1] = 0.0;

#ifdef DEBUG
        for (i=0; i<NS; i++) {
            if (dgds[i] > sqrt(DBL_EPSILON) && ABS(sOld[i]) > DBL_EPSILON) {
                printf("ERROR in ORTHOPYROXENE.C (function ORDER). Failed to converge!\n");
                printf("  X2    = %13.6g, X3    = %13.6g, X4    = %13.6g\n",
                    r[0], r[1], r[2]);
                printf("  X5    = %13.6g, X6    = %13.6g, X7    = %13.6g\n",
                    r[3], r[4], r[5]);
                printf("  s1    = %13.6g, s2    = %13.6g\n", sOld[0], sOld[1]);
                printf("  dgds1 = %13.6g, dgds2 = %13.6g\n", dgds[0], dgds[1]);
                printf("  X Na+  m2:  %13.6g\n",                    xna1m2         );
                printf("  X Ca2+ m2:  %13.6g\n",                    xca2m2         );
                printf("  X Mg2+ m2:  %13.6g  X Mg2+ m1: %13.6g\n", xmg2m2 , xmg2m1);
                printf("  X Fe2+ m2:  %13.6g  X Fe2+ m1: %13.6g\n", xfe2m2 , xfe2m1);
                printf("                             X Ti4+ m1: %13.6g\n",   xti4m1);
                printf("  X Fe3+ tet: %13.6g  X Fe3+ m1: %13.6g\n", xfe3tet, xfe3m1);
                printf("  X Al3+ tet: %13.6g  X Al3+ m1: %13.6g\n", xal3tet, xal3m1);
                printf("  X Si4+ tet: %13.6g\n",                    xsi4tet        );
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
 * Private function to determine structural state of pyroxene
 */

static int isClino(double t, double p, double r[NR])
{
#ifdef ISCLINO
    DECLARE_SITE_FRACTIONS
    double gOrtho, gClino, s[NS];
    int clino;

#ifndef DEBUG
    if ((1.0-r[5]) > 0.20) return TRUE;  /* Ca + Na on M2 > 20 % */
#endif

    clino = TRUE;
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

#ifdef DEBUG
    {
        static char string[] = { "_px Na_.__Ca_.__Fe''_.__Mg_.__Fe'''_.__Ti_.__Al_.__Si_.__O6" };
        double totAl, totCa, totFe2, totFe3, totMg, totNa, totTi, totSi;
        char n[5];
        int i;
        if (gOrtho < gClino) string[0] = 'o'; else string[0] = 'c';

        totAl  = r[1] + r[3];
        totCa  = 1.0 - r[4] - r[5];
        totFe2 = r[0];
        totFe3 = r[2] + r[3];
        totMg  = 1.0 - r[0] - 0.5*r[1] - 0.5*r[2] - r[3] - 0.5*r[4] + r[5];
        totNa  = r[4];
        totTi  = 0.5*(r[1]+r[2]);
        totSi  = 2.0 - r[1] - r[2] - r[3] + r[4]/2.0;

        (void) snprintf(n, 5, "%4.2f", totNa);  for (i=0; i<4; i++) string[ 6+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totCa);  for (i=0; i<4; i++) string[12+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totFe2); for (i=0; i<4; i++) string[20+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totMg);  for (i=0; i<4; i++) string[26+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totFe3); for (i=0; i<4; i++) string[35+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totTi);  for (i=0; i<4; i++) string[41+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totAl);  for (i=0; i<4; i++) string[47+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totSi);  for (i=0; i<4; i++) string[53+i] = n[i];

        printf("ISCLINO: G(c->o) = %f %s\n", gOrtho-gClino, string);
    }
#endif
    if (gOrtho < gClino) return FALSE; else return TRUE;
#else
    return FALSE;
#endif
}

/*
 *=========================================================
 * Unique public function for testing OPXs in LEPR database
 */

int isPigeonite(double t, double p, double r[NR])
{
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
testOpx(int mask, double t, double p,
    int na,          /* Expected number of endmember components                 */
    int nr,          /* Expected number of independent compositional variables  */
    char **names,    /* array of strings of names of endmember components       */
    char **formulas, /* array of strings of formulas of endmember components    */
    double *r,       /* array of indepependent compos variables, check bounds   */
    double *m)       /* array of moles of endmember components, check bounds    */
{
    const char *phase = "orthopyroxene.c";
    const char *NAMES[NA]    = { "diopside", "clinoenstatite", "hedenbergite",
                               "alumino-buffonite", "buffonite",
                               "essenite", "jadeite" };
    const char *FORMULAS[NA] = { "CaMgSi2O6", "Mg2Si2O6", "CaFeSi2O6",
                               "CaTi0.5Mg0.5AlSiO6", "CaTi0.5Mg0.5FeSiO6",
                               "CaFeAlSiO6", "NaAlSi2O6" };
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
        result = result && (r[0] >= 0.0)					&& (r[0] <= 2.0);				     /* Fe2+ */
        result = result && (r[4] >= 0.0)					&& (r[4] <= 1.0);				     /* Na+  */
        result = result && (r[1]+r[3] >= 0.0)				&& (r[1]+r[3] <= 2.0); 				     /* Al3+ */
        result = result && (r[2]+r[3] >= 0.0)				&& (r[2]+r[3] <= 2.0); 				     /* Fe3+ */
        result = result && (r[1]+r[2] >= 0.0)				&& (r[1]+r[2] <= 1.0); 				     /* Ti4+ */
        result = result && (r[4]+r[5] >= 0.0)				&& (r[4]+r[5] <= 1.0); 				     /* Ca2+ */
        result = result && (1.0-r[0]-r[3]+r[5]-0.5*(r[1]+r[2])-r[4] >= 0.0) && (1.0-r[0]-r[3]+r[5]-0.5*(r[1]+r[2])-r[4] <= 2.0); /* Mg2+ */
        /*
        result = result && (2.0-r[1]-r[2]-r[3]+r[4]/2.0 >= 1.0)		&& (2.0-r[1]-r[2]-r[3]+r[4]/2.0 <= 2.0);	     // Si4+
        */
        result = result && (r[5] >= 0.0);				     /* XMgM2 + XFe2+M2 > 0             */
        result = result && (1.0-r[1]-r[2]-r[3]-r[4]/2.0 >= 0.0);	     /* XMgM1 + XFe2+M1 - XTiM1 > 0     */
        result = result && (r[1]/2.0+r[2]/2.0+r[3]-r[4]/2.0 >= 0.0);     /* 1 - XMgM1 - XFe2+M1 - XNaM2 > 0 */
        result = result && (r[1]/2.0+r[2]/2.0+r[3]+r[4]/2.0 >= 0.0);     /* 1 - XMgM1 - XFe2+M1             */
        result = result && (r[1]/2.0+r[2]/2.0+r[3]/2.0+r[4]/2.0 >= 0.0); /* 1 - XSiTet                      */

    }
    /* Check bounds on moles of endmember components */
    if (mask & SIXTH) {
        for (i=0, sum=0.0; i<NA; i++) sum += m[i];
        result = result && (sum > 0.0);
        if (sum > 0.0) {
            /* abundance constraints */
            result = result && (m[2]/sum >= 0.0)			      && (m[2]/sum <= 2.0);			       /* Fe2+ */
            result = result && (m[6]/sum >= 0.0)			      && (m[6]/sum <= 1.0);			       /* Na+  */
            result = result && ((m[3]+m[5]+m[6])/sum >= 0.0)  	      && ((m[3]+m[5]+m[6])/sum <= 2.0); 	       /* Al3+ */
            result = result && ((m[4]+m[5])/sum >= 0.0)		      && ((m[4]+m[5])/sum <= 2.0);		       /* Fe3+ */
            result = result && ((m[3]+m[4])/sum >= 0.0)		      && ((m[3]+m[4])/sum <= 1.0);		       /* Ti4+ */
            result = result && ((m[0]+m[2]+m[3]+m[4]+m[5])/sum >= 0.0)      && ((m[0]+m[2]+m[3]+m[4]+m[5])/sum <= 1.0);      /* Ca2+ */
            result = result && ((m[0]+2.0*m[1]+(m[3]+m[4])/2.0)/sum >= 0.0) && ((m[0]+2.0*m[1]+(m[3]+m[4])/2.0)/sum <= 2.0); /* Mg2+ */

            /* special entropic constraints */
            result = result && (m[1]/sum >= 0.0);					     /* XMgM2 + XFe2+M2 > 0		*/
            result = result && ((1.0 - m[3]/sum - m[4]/sum - m[5]/sum - m[6]/sum) >= 0.0); /* XMgM1 + XFe2+M1 - XTiM1 > 0	*/
            result = result && ((m[5]/sum + (m[3]/sum + m[4]/sum)/2.0) >= 0.0);	     /* 1 - XMgM1 - XFe2+M1 - XNaM2 > 0 */
            result = result && ((m[5]/sum + m[6]/sum + (m[3]/sum + m[4]/sum)/2.0) >= 0.0); /* 1 - XMgM1 - XFe2+M1		*/
            result = result && ((1.0-(sum+m[0]+m[1]+m[2]+m[6])/(2.0*sum)) >= 0.0);	     /* 1 - XSiTet			*/

        }
    }

    return result;
}

void
conOpx(int inpMask, int outMask, double t, double p,
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
    (0)  FIRST            FIRST
    (1)  FIRST            SECOND
    (2)  SECOND           THIRD  | FOURTH  | FIFTH | SIXTH | EIGHTH
    (3)  THIRD            FOURTH | SEVENTH

    (0) use input Fe2O3 for Opx/Cpx as entered (see marc.c)
    (1) converts a vector of moles of elements into a vector of moles of
            endmember pyroxene components.
    (2) calculates from a vector of moles of endmember components, one or
            all of: r[], x[], dr[]/dm[] d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
    (3) calculates from a vector of independent compositional variables
            mole fractions of endmember components and/or the Jacobian matrix
            dx[]/dr[]

    In this routine it is assumed that the elements are in the order of atomic
    numbers and that the order of pyroxene components has been verified as:
            m[0] = diopside              (CaMgSi2O6) ,
            m[1] = enstatite             (Mg2Si2O6),
            m[2] = hedenbergite          (CaFeSi2O6),
            m[3] = Ca(Ti,Mg)(Al,Si)2O6   (CaTi0.5Mg0.5AlSiO6),
            m[4] = Ca(Ti,Mg)(Fe3+,Si)2O6 (CaTi0.5Mg0.5FeSiO6),
            m[5] = essenite              (CaFeAlSiO6)
            m[6] = jadeite               (NaAlSi2O6)

    ----------------------------------------------------------------------------*/

    int i, j, k;

    if (inpMask == FIRST && outMask == SECOND) {
        /* Converts a vector of moles of elements into a vector of moles of
       end-member components.                                                 */
        double sumcat, sumchg, fe2, fe3, corrSi;
        static const int Na = 11;
        static const int Mg = 12;
        static const int Al = 13;
        static const int Si = 14;
        static const int Ca = 20;
        static const int Ti = 22;
        static const int Cr = 24;
        static const int Mn = 25;
        static const int Fe = 26;

        /* Sum the cations and correct the analysis for silica deficiency */
        sumcat  = e[Na] +     e[Mg] +     e[Al] +     e[Si] +     e[Ca] +     e[Ti] +     e[Cr] +     e[Mn] + e[Fe];
        sumchg  = e[Na] + 2.0*e[Mg] + 3.0*e[Al] + 4.0*e[Si] + 2.0*e[Ca] + 4.0*e[Ti] + 3.0*e[Cr] + 2.0*e[Mn];
        corrSi  = (e[Na]+e[Ca] > 0.25*sumcat) ? 4.0*(e[Na]+e[Ca]) - sumcat : 0.0;
        sumcat += corrSi;

#ifdef NEVER_DEFINED
        if (p < 1000.0) {
            /* Compute the ferric/ferrous ratio */
            fe3 = 3.0*sumcat - sumchg - 2.0*e[Fe];
            fe2 = e[Fe] - fe3;
            if (fe3 < 0.0) { fe3 = 0.01*e[Fe]; fe2 = 0.99*e[Fe]; }
            if (fe2 < 0.0) { fe2 = 0.01*e[Fe]; fe3 = 0.99*e[Fe]; }
        } else {
            fe2 = e[Fe];
            fe3 = 0.0;
        }
#else
        fe2 = e[Fe];
        fe3 = 0.0;
#endif

        /* Assign moles of endmembers */
        m[0] = -fe3/2.0 - fe2 - e[Mn] - e[Al]/2.0 - e[Cr]/2.0 + e[Ca] + e[Na]/2.0 - e[Ti];
        m[1] =  fe3/4.0 + fe2/2.0 + e[Mn]/2.0 + e[Al]/4.0 + e[Cr]/4.0 - e[Ca]/2.0 + e[Mg]/2.0 - e[Na]/4.0;
        m[2] =  fe2 + e[Mn];
        m[3] = -fe3/2.0 + e[Al]/2.0 + e[Cr]/2.0 - e[Na]/2.0 + e[Ti];
        m[4] =  fe3/2.0 - e[Al]/2.0 - e[Cr]/2.0 + e[Na]/2.0 + e[Ti];
        m[5] =  fe3/2.0 + e[Al]/2.0 + e[Cr]/2.0 - e[Na]/2.0 - e[Ti];
        m[6] =  e[Na];

#ifdef PHMELTS_ADJUSTMENTS
    }   if (inpMask == FIRST && outMask == FIRST) {
                // use input Fe2O3 for Opx/Cpx as entered (see marc.c)
                double sumcat, sumchg, fe2, fe3;
                static const int Na = 11;
                static const int Mg = 12;
                static const int Al = 13;
                static const int Si = 14;
                static const int Ca = 20;
                static const int Ti = 22;
                static const int Cr = 24;
                static const int Mn = 25;
                static const int Fe = 26;

                /* Sum the cations and correct the analysis for silica deficiency */
                sumcat  = e[Na] +     e[Mg] +     e[Al] +     e[Si] +     e[Ca] +     e[Ti] +     e[Cr] +     e[Mn] + e[Fe];
                sumchg  = e[Na] + 2.0*e[Mg] + 3.0*e[Al] + 4.0*e[Si] + 2.0*e[Ca] + 4.0*e[Ti] + 3.0*e[Cr] + 2.0*e[Mn];

                /* Compute the ferric/ferrous ratio */
                fe3 = 3.0*sumcat - sumchg - 2.0*e[Fe];
                fe2 = e[Fe] - fe3;
                if (fe3 < 0.0) { fe3 = 0.0; fe2 = e[Fe]; }
                if (fe2 < 0.0) { fe2 = 0.0; fe3 = e[Fe]; }

                /* Assign moles of endmembers */
                m[0] = -fe3/2.0 - fe2 - e[Mn] - e[Al]/2.0 - e[Cr]/2.0 + e[Ca] + e[Na]/2.0 - e[Ti];
                m[1] =  fe3/4.0 + fe2/2.0 + e[Mn]/2.0 + e[Al]/4.0 + e[Cr]/4.0 - e[Ca]/2.0 + e[Mg]/2.0 - e[Na]/4.0;
                m[2] =  fe2 + e[Mn];
                m[3] = -fe3/2.0 + e[Al]/2.0 + e[Cr]/2.0 - e[Na]/2.0 + e[Ti];
                m[4] =  fe3/2.0 - e[Al]/2.0 - e[Cr]/2.0 + e[Na]/2.0 + e[Ti];
                m[5] =  fe3/2.0 + e[Al]/2.0 + e[Cr]/2.0 - e[Na]/2.0 - e[Ti];
                m[6] =  e[Na];

#endif
    } else if (inpMask == SECOND) {
        double sum;

        if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
            printf("Illegal call to conOpx with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        for (i=0, sum=0.0; i<NA; i++) sum += m[i];

        if (outMask & THIRD) {
            /* Converts a vector of moles of end-member components (m) into a vector
         of independent compositional variables (r) required as input for the
         remaining public functions.                                          */
            r[0] = (sum != 0.0) ? m[2]/sum            : 0.0;
            r[1] = (sum != 0.0) ? (m[3]+m[6]/2.0)/sum : 0.0;
            r[2] = (sum != 0.0) ? (m[4]-m[6]/2.0)/sum : 0.0;
            r[3] = (sum != 0.0) ? (m[5]+m[6]/2.0)/sum : 0.0;
            r[4] = (sum != 0.0) ? m[6]/sum            : 0.0;
            r[5] = (sum != 0.0) ? m[1]/sum            : 0.0;
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
                    dm[0][j]  = -m[2]/SQUARE(sum);
                    dm[0][j] += (j == 2) ? 1.0/sum : 0.0;
                    dm[1][j]  = -(m[3]+m[6]/2.0)/SQUARE(sum);
                    dm[1][j] += (j == 3) ? 1.0/sum : 0.0;
                    dm[1][j] += (j == 6) ? 0.5/sum : 0.0;
                    dm[2][j]  = -(m[4]-m[6]/2.0)/SQUARE(sum);
                    dm[2][j] += (j == 4) ? 1.0/sum : 0.0;
                    dm[2][j] -= (j == 6) ? 0.5/sum : 0.0;
                    dm[3][j]  = -(m[5]+m[6]/2.0)/SQUARE(sum);
                    dm[3][j] += (j == 5) ? 1.0/sum : 0.0;
                    dm[3][j] += (j == 6) ? 0.5/sum : 0.0;
                    dm[4][j]  = -m[6]/SQUARE(sum);
                    dm[4][j] += (j == 6) ? 1.0/sum : 0.0;
                    dm[5][j]  = -m[1]/SQUARE(sum);
                    dm[5][j] += (j == 1) ? 1.0/sum : 0.0;
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

                        d2m[1][j][k]  = 2.0*(m[3]+m[6]/2.0)/CUBE(sum);
                        d2m[1][j][k] -= (j == 3) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[1][j][k] -= (k == 3) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[1][j][k] -= (j == 6) ? 0.5/SQUARE(sum) : 0.0;
                        d2m[1][j][k] -= (k == 6) ? 0.5/SQUARE(sum) : 0.0;

                        d2m[2][j][k]  = 2.0*(m[4]-m[6]/2.0)/CUBE(sum);
                        d2m[2][j][k] -= (j == 4) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[2][j][k] -= (k == 4) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[2][j][k] += (j == 6) ? 0.5/SQUARE(sum) : 0.0;
                        d2m[2][j][k] += (k == 6) ? 0.5/SQUARE(sum) : 0.0;

                        d2m[3][j][k]  = 2.0*(m[5]+m[6]/2.0)/CUBE(sum);
                        d2m[3][j][k] -= (j == 5) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[3][j][k] -= (k == 5) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[3][j][k] -= (j == 6) ? 0.5/SQUARE(sum) : 0.0;
                        d2m[3][j][k] -= (k == 6) ? 0.5/SQUARE(sum) : 0.0;

                        d2m[4][j][k]  = 2.0*m[6]/CUBE(sum);
                        d2m[4][j][k] -= (j == 6) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[4][j][k] -= (k == 6) ? 1.0/SQUARE(sum) : 0.0;

                        d2m[5][j][k]  = 2.0*m[1]/CUBE(sum);
                        d2m[5][j][k] -= (j == 1) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[5][j][k] -= (k == 1) ? 1.0/SQUARE(sum) : 0.0;
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

                            d3m[1][j][k][l]  = -6.0*(m[3]+m[6]/2.0)/QUARTIC(sum);
                            d3m[1][j][k][l] += (j == 3) ? 2.0/CUBE(sum) : 0.0;
                            d3m[1][j][k][l] += (k == 3) ? 2.0/CUBE(sum) : 0.0;
                            d3m[1][j][k][l] += (l == 3) ? 2.0/CUBE(sum) : 0.0;
                            d3m[1][j][k][l] += (j == 6) ? 1.0/CUBE(sum) : 0.0;
                            d3m[1][j][k][l] += (k == 6) ? 1.0/CUBE(sum) : 0.0;
                            d3m[1][j][k][l] += (l == 6) ? 1.0/CUBE(sum) : 0.0;

                            d3m[2][j][k][l]  = -6.0*(m[4]-m[6]/2.0)/QUARTIC(sum);
                            d3m[2][j][k][l] += (j == 4) ? 2.0/CUBE(sum) : 0.0;
                            d3m[2][j][k][l] += (k == 4) ? 2.0/CUBE(sum) : 0.0;
                            d3m[2][j][k][l] += (l == 4) ? 2.0/CUBE(sum) : 0.0;
                            d3m[2][j][k][l] -= (j == 6) ? 1.0/CUBE(sum) : 0.0;
                            d3m[2][j][k][l] -= (k == 6) ? 1.0/CUBE(sum) : 0.0;
                            d3m[2][j][k][l] -= (l == 6) ? 1.0/CUBE(sum) : 0.0;

                            d3m[3][j][k][l]  = -6.0*(m[5]+m[6]/2.0)/QUARTIC(sum);
                            d3m[3][j][k][l] += (j == 5) ? 2.0/CUBE(sum) : 0.0;
                            d3m[3][j][k][l] += (k == 5) ? 2.0/CUBE(sum) : 0.0;
                            d3m[3][j][k][l] += (l == 5) ? 2.0/CUBE(sum) : 0.0;
                            d3m[3][j][k][l] += (j == 6) ? 1.0/CUBE(sum) : 0.0;
                            d3m[3][j][k][l] += (k == 6) ? 1.0/CUBE(sum) : 0.0;
                            d3m[3][j][k][l] += (l == 6) ? 1.0/CUBE(sum) : 0.0;

                            d3m[4][j][k][l]  = -6.0*m[6]/QUARTIC(sum);
                            d3m[4][j][k][l] += (j == 6) ? 2.0/CUBE(sum) : 0.0;
                            d3m[4][j][k][l] += (k == 6) ? 2.0/CUBE(sum) : 0.0;
                            d3m[4][j][k][l] += (l == 6) ? 2.0/CUBE(sum) : 0.0;

                            d3m[5][j][k][l]  = -6.0*m[1]/QUARTIC(sum);
                            d3m[5][j][k][l] += (j == 1) ? 2.0/CUBE(sum) : 0.0;
                            d3m[5][j][k][l] += (k == 1) ? 2.0/CUBE(sum) : 0.0;
                            d3m[5][j][k][l] += (l == 1) ? 2.0/CUBE(sum) : 0.0;
                        }
                    }
                }
            }
        }

    } else if (inpMask == THIRD) {

        if (outMask & ~(FOURTH | SEVENTH))
            printf("Illegal call to conOpx with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        if (outMask & FOURTH) {
            /* Converts a vector of independent compositional variables (r) into a
         vector of mole fractions of endmember components (x).                */
            x[0] = 1.0 - r[0] - r[1] - r[2] - r[3] - r[4]/2.0 - r[5];
            x[1] = r[5];
            x[2] = r[0];
            x[3] = r[1] - r[4]/2.0;
            x[4] = r[2] + r[4]/2.0;
            x[5] = r[3] - r[4]/2.0;
            x[6] = r[4];
        }

        if (outMask & SEVENTH) {
            /* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
            for (i=0; i<NA; i++) for (j=0; j<NR; j++) dr[i][j] = 0.0;
            dr[0][0] = -1.0; dr[0][1] = -1.0; dr[0][2] = -1.0; dr[0][3] = -1.0;
                       dr[0][4] = -0.5; dr[0][5] = -1.0;
            dr[1][5] =  1.0;
            dr[2][0] =  1.0;
            dr[3][1] =  1.0; dr[3][4] = -0.5;
            dr[4][2] =  1.0; dr[4][4] =  0.5;
            dr[5][3] =  1.0; dr[5][4] = -0.5;
            dr[6][4] =  1.0;
        }

    } else  {
        printf("Illegal call to conOpx with inpMask = %o and outMask = %o\n",
            inpMask, outMask);
    }

}

void
dispOpx(int mask, double t, double p, double *x,
    char **formula            /* Mineral formula for interface display MASK: 1 */
    )
{
    double *r = x;
    static char masterString[] = {
/*             1111111111222222222233333333334444444444555555555566666666667
     01234567890123456789012345678901234567890123456789012345678901234567890 */
#ifdef PHMELTS_ADJUSTMENTS
        "    Na_.__Ca_.__Fe''_.__Mg_.__Fe'''_.__Ti_.__Al_.__Si_.__O6" };
#else
        "_px Na_.__Ca_.__Fe''_.__Mg_.__Fe'''_.__Ti_.__Al_.__Si_.__O6" };
#endif

    if (mask & FIRST) {
        char *string = strcpy((char *) malloc((size_t) (strlen(masterString)+1)*sizeof(char)), masterString);
        double totAl, totCa, totFe2, totFe3, totMg, totNa, totTi, totSi;
        char n[5];
        int i;
#ifndef PHMELTS_ADJUSTMENTS
        if (isClino(t, p, r) == TRUE) string[0] = 'c'; else string[0] = 'o';
#endif

        totAl  = r[1] + r[3];
        totCa  = 1.0 - r[4] - r[5];
        totFe2 = r[0];
        totFe3 = r[2] + r[3];
        totMg  = 1.0 - r[0] - 0.5*r[1] - 0.5*r[2] - r[3] - 0.5*r[4] + r[5];
        totNa  = r[4];
        totTi  = 0.5*(r[1]+r[2]);
        totSi  = 2.0 - r[1] - r[2] - r[3] + r[4]/2.0;

        (void) snprintf(n, 5, "%4.2f", totNa);  for (i=0; i<4; i++) string[ 6+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totCa);  for (i=0; i<4; i++) string[12+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totFe2); for (i=0; i<4; i++) string[20+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totMg);  for (i=0; i<4; i++) string[26+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totFe3); for (i=0; i<4; i++) string[35+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totTi);  for (i=0; i<4; i++) string[41+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totAl);  for (i=0; i<4; i++) string[47+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totSi);  for (i=0; i<4; i++) string[53+i] = n[i];

        *formula = string;
    }
}

void
actOpx(int mask, double t, double p, double *x,
    double *a,  /* (pointer to a[]) activities              BINARY MASK: 0001 */
    double *mu, /* (pointer to mu[]) chemical potentials    BINARY MASK: 0010 */
    double **dx /* (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100 */
    )           /* exclusion criteria applied to results if BINARY MASK: 1000 */
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS], g, dgdr[NR];
    double fr[NA][NR];
    int i, j, clino;

    for(i=0; i<NA; i++) {
     fr[i][0] = FR2(i); /* X2 */
     fr[i][1] = FR3(i); /* X3 */
     fr[i][2] = FR4(i); /* X4 */
     fr[i][3] = FR5(i); /* X5 */
     fr[i][4] = FR6(i); /* X6 */
     fr[i][5] = FR7(i); /* X7 */
    }

    clino = isClino(t, p, r);
    setClino(clino);

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    g       = G;
    dgdr[0] = DGDR0;
    dgdr[1] = DGDR1;
    dgdr[2] = DGDR2;
    dgdr[3] = DGDR3;
    dgdr[4] = DGDR4;
    dgdr[5] = DGDR5;

    /* activities for library */
    if (!mask && a != NULL) {
        for(i=0; i<NA; i++) {
       for (a[i]=g, j=0; j<NR; j++) a[i] += fr[i][j]*dgdr[j];
       a[i] = exp(a[i]/(R*t));
        }
    }

    if (mask & FIRST) {
        double a0[NA];

        purePyx(FIRST, t, p,
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

        purePyx(SECOND, t, p,
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

        purePyx(FIRST, t, p,
       a0,              (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        for(i=0; i<NA; i++) {
       gs[i][0] = FS1(i); /* s1 */
       gs[i][1] = FS2(i); /* s2 */
       dfrdr[i][0] = DFR2DR2(i); /* X2 */
       dfrdr[i][1] = DFR3DR3(i); /* X3 */
       dfrdr[i][2] = DFR4DR4(i); /* X4 */
       dfrdr[i][3] = DFR5DR5(i); /* X5 */
       dfrdr[i][4] = DFR6DR6(i); /* X6 */
       dfrdr[i][5] = DFR7DR7(i); /* X7 */
       dgsds[i][0] = DFS1DS1(i); /* s1 */
       dgsds[i][1] = DFS2DS2(i); /* s2 */
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
        int flags[NA];

        flags[0] = (xca2m2  < 0.05) | (xmg2m1  < 0.05) | (xsi4tet < 0.05);
        flags[1] = (xmg2m2  < 0.05) | (xmg2m1  < 0.05) | (xsi4tet < 0.05);
        flags[2] = (xca2m2  < 0.05) | (xfe2m1  < 0.05) | (xsi4tet < 0.05);
        flags[3] = (xca2m2  < 0.05) | (xmg2m1  < 0.05) | (xti4m1  < 0.05) |
               (xal3tet < 0.05) | (xsi4tet < 0.05);
        flags[4] = (xca2m2  < 0.05) | (xmg2m1  < 0.05) | (xti4m1  < 0.05) |
               (xfe3tet < 0.05) | (xsi4tet < 0.05);
        flags[5] = (xca2m2  < 0.05) | (xfe3m1  < 0.05) | (xal3tet < 0.05) |
               (xsi4tet < 0.05);
        flags[6] = (xna1m2  < 0.05) | (xal3m1 < 0.05)  | (xsi4tet < 0.05);

#ifdef REGRESS_DI_EN_ONLY
        flags[2] = 1;
        flags[3] = 1;
        flags[4] = 1;
        flags[5] = 1;
        flags[6] = 1;
#endif

#ifdef REGRESS_LOW_P_ONLY
        if (p > 2.0) {
            flags[0] = 1;
            flags[1] = 1;
            flags[2] = 1;
            flags[3] = 1;
            flags[4] = 1;
            flags[5] = 1;
            flags[6] = 1;
        }
#endif

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
gmixOpx(int mask, double t, double p, double *x,
    double *gmix, /* Gibbs energy of mixing             BINARY MASK: 0001 */
    double *dx,   /* (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010 */
    double **dx2, /* (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100 */
    double ***dx3 /* (pointer to dx3[][][]) d3(g)/d(x[])3 NARY MASK: 1000 */
    )
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS];
    int clino;

    clino = isClino(t, p, r);
    setClino(clino);

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    if (mask & FIRST) {
        double ends[NA];

        *gmix = G;

        purePyx(THIRD, t, p,
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
        dx[4] = DGDR4;
        dx[5] = DGDR5;

        purePyx(THIRD, t, p,
       (double *) NULL, (double *) NULL, ends,            (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        dx[0] -= DENDDR0;
        dx[1] -= DENDDR1;
        dx[2] -= DENDDR2;
        dx[3] -= DENDDR3;
        dx[4] -= DENDDR4;
        dx[5] -= DENDDR5;
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
hmixOpx(int mask, double t, double p, double *x,
    double *hmix /* Enthalpy of mixing BINARY MASK: 1 */
    )
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS], ends[NA];
    int clino;

    clino = isClino(t, p, r);
    setClino(clino);

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    *hmix = (G) + t*(S);

    purePyx(FOURTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, ends,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

    *hmix -= ENDMEMBERS;
}

void
smixOpx(int mask, double t, double p, double *x,
    double *smix, /* Entropy of mixing                  BINARY MASK: 001 */
    double *dx,   /* (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010 */
    double **dx2  /* (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100 */
    )
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS];
    int clino;

    clino = isClino(t, p, r);
    setClino(clino);

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    if (mask & FIRST) {
        double ends[NA];

        *smix = S;

        purePyx(FIFTH, t, p,
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

        purePyx(FIFTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       ends,            (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        dx[0] -= DENDDR0;
        dx[1] -= DENDDR1;
        dx[2] -= DENDDR2;
        dx[3] -= DENDDR3;
        dx[4] -= DENDDR4;
        dx[5] -= DENDDR5;
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
cpmixOpx(int mask, double t, double p, double *x,
    double *cpmix, /* Heat capacity of mixing               BINARY MASK: 001 */
    double *dt,    /* d(cp)/d(t)                            BINARY MASK: 010 */
    double *dx     /* d(cp)/d(x[])d(t)                      BINARY MASK: 100 */
    )
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS], dsdt[NS], d2gdsdt[NS], d2gds2[NS][NS], d2gdt2;
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
        double ends[NA];

        *cpmix = d2gdt2;
        for (i=0; i<NS; i++) {
            *cpmix += 2.0*d2gdsdt[i]*dsdt[i];
            for (j=0; j<NS; j++) *cpmix += d2gds2[i][j]*dsdt[i]*dsdt[j];
        }
        *cpmix *= -t;

        purePyx(SIXTH, t, p,
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

        purePyx(SEVENTH, t, p,
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

        purePyx(SIXTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, ends,            (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        dx[0] -= DENDDR0;
        dx[1] -= DENDDR1;
        dx[2] -= DENDDR2;
        dx[3] -= DENDDR3;
        dx[4] -= DENDDR4;
        dx[5] -= DENDDR5;
    }

}

void
vmixOpx(int mask, double t, double p, double *x,
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
    int clino;

    clino = isClino(t, p, r);
    setClino(clino);

    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    GET_SITE_FRACTIONS

    if (mask & FIRST) {
        double ends[NA];

        *vmix = DGDP;

        purePyx(EIGHTH, t, p,
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

        purePyx(EIGHTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, ends,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        dx[0] -= DENDDR0;
        dx[1] -= DENDDR1;
        dx[2] -= DENDDR2;
        dx[3] -= DENDDR3;
        dx[4] -= DENDDR4;
        dx[5] -= DENDDR5;
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

        purePyx(NINTH, t, p,
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

        purePyx(TENTH, t, p,
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

        purePyx(ELEVENTH, t, p,
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

        purePyx(TWELFTH, t, p,
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

        purePyx(THIRTEENTH, t, p,
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

        purePyx(NINTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       ends,            (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL);

        dxdt[0] -= DENDDR0;
        dxdt[1] -= DENDDR1;
        dxdt[2] -= DENDDR2;
        dxdt[3] -= DENDDR3;
        dxdt[4] -= DENDDR4;
        dxdt[5] -= DENDDR5;
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

        purePyx(TENTH, t, p,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, (double *) NULL, (double *) NULL, (double *) NULL,
       (double *) NULL, ends,            (double *) NULL, (double *) NULL,
       (double *) NULL);

        dxdp[0] -= DENDDR0;
        dxdp[1] -= DENDDR1;
        dxdp[2] -= DENDDR2;
        dxdp[3] -= DENDDR3;
        dxdp[4] -= DENDDR4;
        dxdp[5] -= DENDDR5;
    }

}

/* end of file ORTHOPYROXENE.C */


