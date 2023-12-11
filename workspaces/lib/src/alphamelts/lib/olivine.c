const char *olivine_ver(void) { return "$Id: olivine.c,v 1.4 2007/03/12 20:06:36 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: olivine.c,v $
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
 * Revision 3.8  1997/06/21  22:49:35  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.7  1997/05/03  20:23:12  ghiorso
 * *** empty log message ***
 *
 * Revision 3.6  1997/03/27  17:03:19  ghiorso
 * *** empty log message ***
 *
 * Revision 3.5  1996/09/24  20:33:26  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 3.4  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 3.3  1995/11/01  22:40:27  ghiorso
 * Implementation of subsolidus options after Asimow.
 * Additional implementation of nepheline solid solutions.
 *
 * Revision 3.3  1995/11/01  22:40:27  ghiorso
 * Implementation of subsolidus options after Asimow.
 * Additional implementation of nepheline solid solutions.
 *
 * Revision 3.2  1995/09/09  23:37:35  ghiorso
 * Modifications by Asimow to include new Derivatives of gmix and convert
 * options for liquid-absent fO2 buffering. Also new olivine components
 * implemented.
 *
 * Revision 3.1  1995/08/18  18:04:10  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Routines to compute (Mg,Fe,Ca,Mn,Co,Ni)2SiO4 olivine activities
**      (file: OLIVINE.C)
**
**  MODIFICATION HISTORY:
**      V1.0-1  Modified by Marc Hirschmann January, 1992
**              from SPINEL.C V1.0-14 by Mark Ghiorso
**      V2.0-1  Paul D. Asimow May 8, 1995
**              Import changes for MELTS2 from spinel.c V2.0-2
**      V3.0-1  Paul D. Asimow  July 27, 1995
**              Add d3rdm3 to (*convert) and d3gdx3 to (*gmix)
**
**--
*/

#ifdef DEBUG
#undef DEBUG
#endif

#include "melts_gsl.h"
#include "silmin.h"  /* Structure definitions for SILMIN package */

#define SQUARE(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

//#define USE_SVD (calculationMode != MODE_pMELTS)
#define USE_SVD 1

/*
 *=============================================================================
 * Olivine solution parameters:
 *  Hirschmann, M. (1991)  Thermodynamics of multicomponent olivines and the
 *    solution properties of (Ni,Mg,Fe)2SiO4 and (Ca,Mg,Fe)SiO4 olivines
 *    American Mineralogist 77:1232-1248.
 *
 * Sack, R.O., Ghiorso, M.S. (1989)
 *   Importance of considerations of mixing properties in establishing
 *   an internally consistent thermodyanmic database:  Thermochemistry of
 *   minerals in the system Mg2SIO4-Fe2SiO4-SiO2 Contributions to Mineralogy
 *   and Petrology 102: 41-68
 *
 *  Properties of Mn- and Co- bearing olivines - Hirschmann and Ghiorso
 *  1994 GCA
 */


#define HEXMGMN  15.80 *1000.0 /* joules */
#define HEXMGFE  00.00 *1000.0 /* joules */
#define HEXMGCO -15.00 *1000.0 /* joules */
#define HEXMGNI -19.75 *1000.0 /* joules */
#define HEXMNFE -11.80 *1000.0 /* joules */
#define HEXMNCO  00.00 *1000.0 /* joules */
#define HEXMNNI  00.00 *1000.0 /* joules */
#define HEXFECO  00.00 *1000.0 /* joules */
#define HEXFENI -20.00 *1000.0 /* joules */
#define HEXCONI   0.00 *1000.0 /* joules */

#define VEXMGMN  00.00         /* joules/bar */
#define VEXMGFE  00.00         /* joules/bar */
#define VEXMGCO  00.00         /* joules/bar */
#define VEXMGNI  00.000        /* joules/bar */
#define VEXMNFE  00.00         /* joules/bar */
#define VEXMNCO  00.00         /* joules/bar */
#define VEXMNNI  00.00         /* joules/bar */
#define VEXFECO  00.00         /* joules/bar */
#define VEXFENI  00.000        /* joules/bar */
#define VEXCONI  00.00         /* joules/bar */

#define HXMGMN    8.75 *1000.0 /* joules */
#define HXMGFE   10.15 *1000.0 /* joules */
#define HXMGCO   03.00 *1000.0 /* joules */
#define HXMGNI    2.20 *1000.0 /* joules */
#define HXMNFE    0.50 *1000.0 /* joules */
#define HXMNCO    0.00 *1000.0 /* joules */
#define HXMNNI   00.00 *1000.0 /* joules */
#define HXFECO   03.00 *1000.0 /* joules */
#define HXFENI   10.00 *1000.0 /* joules */
#define HXCONI    0.00 *1000.0 /* joules */

#define VXMGMN   00.00         /* joules/bar */
#define VXMGFE   00.015        /* joules/bar */
#define VXMGCO   00.00         /* joules/bar */
#define VXMGNI   00.000        /* joules/bar */
#define VXMNFE   00.00         /* joules/bar */
#define VXMNCO   00.00         /* joules/bar */
#define VXMNNI   00.00         /* joules/bar */
#define VXFECO   00.00         /* joules/bar */
#define VXFENI   00.000        /* joules/bar */
#define VXCONI   00.00         /* joules/bar */

#define WH1MGMN    6.625*1000.0 /* joules */
#define WH2MGMN    6.625*1000.0 /* joules */
#define WH1MGFE    5.075*1000.0 /* joules */
#define WH2MGFE    5.075*1000.0 /* joules */
#define WH1MGCO    1.50 *1000.0 /* joules */
#define WH2MGCO    1.50 *1000.0 /* joules */
#define WH1MGNI    -.600*1000.0 /* joules */
#define WH2MGNI    2.800*1000.0 /* joules */
#define WH1MNFE    1.75 *1000.0 /* joules */
#define WH2MNFE    1.75 *1000.0 /* joules */
#define WH1MNCO    0.00 *1000.0 /* joules */
#define WH2MNCO    0.00 *1000.0 /* joules */
#define WH1MNNI    0.00 *1000.0 /* joules */
#define WH2MNNI    0.00 *1000.0 /* joules */
#define WH1FECO    1.50 *1000.0 /* joules */
#define WH2FECO    1.50 *1000.0 /* joules */
#define WH1FENI    5.000*1000.0 /* joules */
#define WH2FENI    5.000*1000.0 /* joules */
#define WH1CONI    0.00 *1000.0 /* joules */
#define WH2CONI    0.00 *1000.0 /* joules */

#define WV1MGMN   00.00         /* joules/bar */
#define WV2MGMN   00.00         /* joules/bar */
#define WV1MGFE   00.0000       /* joules/bar */
#define WV2MGFE   00.0000       /* joules/bar */
#define WV1MGCO   00.00         /* joules/bar */
#define WV2MGCO   00.00         /* joules/bar */
#define WV1MGNI   00.0000       /* joules/bar */
#define WV2MGNI   00.0000       /* joules/bar */
#define WV1MNFE   00.00         /* joules/bar */
#define WV2MNFE   00.00         /* joules/bar */
#define WV1MNCO   00.00         /* joules/bar */
#define WV2MNCO   00.00         /* joules/bar */
#define WV1MNNI   00.00         /* joules/bar */
#define WV2MNNI   00.00         /* joules/bar */
#define WV1FECO   00.00         /* joules/bar */
#define WV2FECO   00.00         /* joules/bar */
#define WV1FENI   00.0000       /* joules/bar */
#define WV2FENI   00.0000       /* joules/bar */
#define WV1CONI   00.00         /* joules/bar */
#define WV2CONI   00.00         /* joules/bar */

#define WH2CAMG   34.50 *1000.0 /* joules */
#define WH2CAMN   16.00 *1000.0 /* joules */
#define WH2CAFE   21.90 *1000.0 /* joules */
#define WH2CACO   30.00 *1000.0 /* joules */
#define WH2CANI   40.00 *1000.0 /* joules */

#define WV2CAMG   00.35         /* joules/bar */
#define WV2CAMN   00.00         /* joules/bar */
#define WV2CAFE   00.00         /* joules/bar */
#define WV2CACO   00.00         /* joules/bar */
#define WV2CANI   00.00         /* joules/bar */

#define F_MN     09.50 *1000.0 /* joules */
#define F_FE     09.50 *1000.0 /* joules */
#define F_CO     00.00 *1000.0 /* joules */
#define F_NI     00.00 *1000.0 /* joules */

#define GEXMGMN  (HEXMGMN) + (p-1.0)*(VEXMGMN)
#define GEXMGFE  (HEXMGFE) + (p-1.0)*(VEXMGFE)
#define GEXMGCO  (HEXMGCO) + (p-1.0)*(VEXMGCO)
#define GEXMGNI  (HEXMGNI) + (p-1.0)*(VEXMGNI)
#define GEXMNFE  (HEXMNFE) + (p-1.0)*(VEXMNFE)
#define GEXMNCO  (HEXMNCO) + (p-1.0)*(VEXMNCO)
#define GEXMNNI  (HEXMNNI) + (p-1.0)*(VEXMNNI)
#define GEXFECO  (HEXFECO) + (p-1.0)*(VEXFECO)
#define GEXFENI  (HEXFENI) + (p-1.0)*(VEXFENI)
#define GEXCONI  (HEXCONI) + (p-1.0)*(VEXCONI)

#define GXMGMN  (HXMGMN) + (p-1.0)*(VXMGMN)
#define GXMGFE  (HXMGFE) + (p-1.0)*(VXMGFE)
#define GXMGCO  (HXMGCO) + (p-1.0)*(VXMGCO)
#define GXMGNI  (HXMGNI) + (p-1.0)*(VXMGNI)
#define GXMNFE  (HXMNFE) + (p-1.0)*(VXMNFE)
#define GXMNCO  (HXMNCO) + (p-1.0)*(VXMNCO)
#define GXMNNI  (HXMNNI) + (p-1.0)*(VXMNNI)
#define GXFECO  (HXFECO) + (p-1.0)*(VXFECO)
#define GXFENI  (HXFENI) + (p-1.0)*(VXFENI)
#define GXCONI  (HXCONI) + (p-1.0)*(VXCONI)

#define W1MGMN  (WH1MGMN) + (p-1.0)*(WV1MGMN)
#define W2MGMN  (WH2MGMN) + (p-1.0)*(WV2MGMN)
#define W1MGFE  (WH1MGFE) + (p-1.0)*(WV1MGFE)
#define W2MGFE  (WH2MGFE) + (p-1.0)*(WV2MGFE)
#define W1MGCO  (WH1MGCO) + (p-1.0)*(WV1MGCO)
#define W2MGCO  (WH2MGCO) + (p-1.0)*(WV2MGCO)
#define W1MGNI  (WH1MGNI) + (p-1.0)*(WV1MGNI)
#define W2MGNI  (WH2MGNI) + (p-1.0)*(WV2MGNI)
#define W1MNFE  (WH1MNFE) + (p-1.0)*(WV1MNFE)
#define W2MNFE  (WH2MNFE) + (p-1.0)*(WV2MNFE)
#define W1MNCO  (WH1MNCO) + (p-1.0)*(WV1MNCO)
#define W2MNCO  (WH2MNCO) + (p-1.0)*(WV2MNCO)
#define W1MNNI  (WH1MNNI) + (p-1.0)*(WV1MNNI)
#define W2MNNI  (WH2MNNI) + (p-1.0)*(WV2MNNI)
#define W1FECO  (WH1FECO) + (p-1.0)*(WV1FECO)
#define W2FECO  (WH2FECO) + (p-1.0)*(WV2FECO)
#define W1FENI  (WH1FENI) + (p-1.0)*(WV1FENI)
#define W2FENI  (WH2FENI) + (p-1.0)*(WV2FENI)
#define W1CONI  (WH1CONI) + (p-1.0)*(WV1CONI)
#define W2CONI  (WH2CONI) + (p-1.0)*(WV2CONI)

#define W2CAMG  (WH2CAMG) + (p-1.0)*(WV2CAMG)
#define W2CAMN  (WH2CAMN) + (p-1.0)*(WV2CAMN)
#define W2CAFE  (WH2CAFE) + (p-1.0)*(WV2CAFE)
#define W2CACO  (WH2CACO) + (p-1.0)*(WV2CACO)
#define W2CANI  (WH2CANI) + (p-1.0)*(WV2CANI)


 /* Definitions of Taylor expansion coefficients in terms of solution
    * parameters. Independent variables are r1,r2,r3,r4,r5,s1,s2,s3,s4
    */

#define G0  0.25*(     ((GXMNFE)+(W1MNFE)+(W2MNFE)) \
                                            +((GXMNCO)+(W1MNCO)+(W2MNCO)) \
                                            +((GXMNNI)+(W1MNNI)+(W2MNNI)) \
                                            +((GXFECO)+(W1FECO)+(W2FECO)) \
                                            +((GXFENI)+(W1FENI)+(W2FENI)) \
                                            +((GXCONI)+(W1CONI)+(W2CONI)) \
                                    -2.0*((GXMGMN)+(W1MGMN)+(W2MGMN)) \
                                    -2.0*((GXMGFE)+(W1MGFE)+(W2MGFE)) \
                                    -2.0*((GXMGCO)+(W1MGCO)+(W2MGCO)) \
                                    -2.0*((GXMGNI)+(W1MGNI)+(W2MGNI)))
#define GR1 0.25*(-3.0*((GXMGMN)+(W1MGMN)+(W2MGMN)) \
                                            -((GXMGFE)+(W1MGFE)+(W2MGFE)) \
                                            -((GXMGCO)+(W1MGCO)+(W2MGCO)) \
                                            -((GXMGNI)+(W1MGNI)+(W2MGNI)) \
                                            +((GXMNFE)+(W1MNFE)+(W2MNFE)) \
                                            +((GXMNCO)+(W1MNCO)+(W2MNCO)) \
                                            +((GXMNNI)+(W1MNNI)+(W2MNNI)))
#define GR2  0.25*(   -((GXMGMN)+(W1MGMN)+(W2MGMN)) \
                                    -3.0*((GXMGFE)+(W1MGFE)+(W2MGFE)) \
                                            -((GXMGCO)+(W1MGCO)+(W2MGCO)) \
                                            -((GXMGNI)+(W1MGNI)+(W2MGNI)) \
                                            +((GXMNFE)+(W1MNFE)+(W2MNFE)) \
                                            +((GXFECO)+(W1FECO)+(W2FECO)) \
                                            +((GXFENI)+(W1FENI)+(W2FENI)))
#define GR3  0.25*(   -((GXMGMN)+(W1MGMN)+(W2MGMN)) \
                                            -((GXMGFE)+(W1MGFE)+(W2MGFE)) \
                                    -3.0*((GXMGCO)+(W1MGCO)+(W2MGCO)) \
                                            -((GXMGNI)+(W1MGNI)+(W2MGNI)) \
                                            +((GXMNCO)+(W1MNCO)+(W2MNCO)) \
                                            +((GXFECO)+(W1FECO)+(W2FECO)) \
                                            +((GXCONI)+(W1CONI)+(W2CONI)))
#define GR4 0.25*(    -((GXMGMN)+(W1MGMN)+(W2MGMN)) \
                                            -((GXMGFE)+(W1MGFE)+(W2MGFE)) \
                                            -((GXMGCO)+(W1MGCO)+(W2MGCO)) \
                                    -3.0*((GXMGNI)+(W1MGNI)+(W2MGNI)) \
                                            +((GXMNNI)+(W1MNNI)+(W2MNNI)) \
                                            +((GXFENI)+(W1FENI)+(W2FENI)) \
                                            +((GXCONI)+(W1CONI)+(W2CONI)))
#define GR5        -(W2CAMG) \
                                    -0.25*(((F_MN)+(GEXMGMN)+(GXMGMN) \
                                                -2.0*(W2CAMN)+2.0*(W2MGMN)) \
                                                +((F_FE)+(GEXMGFE)+(GXMGFE) \
                                                -2.0*(W2CAFE)+2.0*(W2MGFE)) \
                                                +((F_CO)+(GEXMGCO)+(GXMGCO) \
                                                -2.0*(W2CACO)+2.0*(W2MGCO)) \
                                                +((F_NI)+(GEXMGNI)+(GXMGNI) \
                                                -2.0*(W2CANI)+2.0*(W2MGNI)))
#define GR1R5             -0.25*((F_MN)+(GEXMGMN)+(GXMGMN) \
                   +2.0*(W2CAMG)-2.0*(W2CAMN)+2.0*(W2MGMN))
#define GR2R5             -0.25*((F_FE)+(GEXMGFE)+(GXMGFE) \
                   +2.0*(W2CAMG)-2.0*(W2CAFE)+2.0*(W2MGFE))
#define GR3R5             -0.25*((F_CO)+(GEXMGCO)+(GXMGCO) \
                   +2.0*(W2CAMG)-2.0*(W2CACO)+2.0*(W2MGCO))
#define GS1   0.25*(( (GEXMGMN)+3.0*(W1MGMN)-3.0*(W2MGMN)) \
                                                        -((GEXMGFE)-(W1MGFE)+(W2MGFE)) \
                                                        -((GEXMGCO)-(W1MGCO)+(W2MGCO)) \
                                                        -((GEXMGNI)-(W1MGNI)+(W2MGNI)) \
                                                        +((GEXMNFE)-(W1MNFE)+(W2MNFE)) \
                                                        +((GEXMNCO)-(W1MNCO)+(W2MNCO)) \
                                                        +((GEXMNNI)-(W1MNNI)+(W2MNNI)))
#define GS2           0.25*((-(GEXMGMN)+(W1MGMN)-(W2MGMN)) \
                                        +((GEXMGFE)+3.0*(W1MGFE)-3.0*(W2MGFE)) \
                                                        -((GEXMGCO)-(W1MGCO)+(W2MGCO)) \
                                                        -((GEXMGNI)-(W1MGNI)+(W2MGNI)) \
                                                        -((GEXMNFE)+(W1MNFE)-(W2MNFE)) \
                                                        +((GEXFECO)-(W1FECO)+(W2FECO)) \
                                                        +((GEXFENI)-(W1FENI)+(W2FENI)))
#define GS3           0.25*(( (GEXMGMN)-(W1MGMN)+(W2MGMN)) \
                                                        +((GEXMGFE)-(W1MGFE)+(W2MGFE)) \
                                        -((GEXMGCO)+3.0*(W1MGCO)-3.0*(W2MGCO)) \
                                                        +((GEXMGNI)-(W1MGNI)+(W2MGNI)) \
                                                        +((GEXMNCO)+(W1MNCO)-(W2MNCO)) \
                                                        +((GEXFECO)+(W1FECO)-(W2FECO)) \
                                                        -((GEXCONI)-(W1CONI)+(W2CONI)))
#define GS4           0.25*(( (GEXMGMN)-(W1MGMN)+(W2MGMN)) \
                                                        +((GEXMGFE)-(W1MGFE)+(W2MGFE)) \
                                                        +((GEXMGCO)-(W1MGCO)+(W2MGCO)) \
                                        -((GEXMGNI)+3.0*(W1MGNI)-3.0*(W2MGNI)) \
                                                        +((GEXMNNI)+(W1MNNI)-(W2MNNI)) \
                                                        +((GEXFENI)+(W1FENI)-(W2FENI)) \
                                                        +((GEXCONI)+(W1CONI)-(W2CONI)))
#define GR1R1         -0.25*(  (GXMGMN)+(W1MGMN)+(W2MGMN))
#define GR1R2          0.25*(( (GXMNFE)+(W1MNFE)+(W2MNFE)) \
                             -((GXMGMN)+(W1MGMN)+(W2MGMN)) \
                             -((GXMGFE)+(W1MGFE)+(W2MGFE)))
#define GR1R3          0.25*(( (GXMNCO)+(W1MNCO)+(W2MNCO)) \
                             -((GXMGMN)+(W1MGMN)+(W2MGMN)) \
                             -((GXMGCO)+(W1MGCO)+(W2MGCO)))
#define GR1R4          0.25*(( (GXMNNI)+(W1MNNI)+(W2MNNI)) \
                             -((GXMGMN)+(W1MGMN)+(W2MGMN)) \
                             -((GXMGNI)+(W1MGNI)+(W2MGNI)))
#define GR1R5             -0.25*((F_MN)+(GEXMGMN)+(GXMGMN) \
                   +2.0*(W2CAMG)-2.0*(W2CAMN)+2.0*(W2MGMN))
#define GR1S1                       0.5*((W1MGMN)-(W2MGMN))
#define GR1S2         0.25*((-(GEXMGMN)-(W1MGMN)+(W2MGMN)) \
                                                        +((GEXMGFE)+(W1MGFE)-(W2MGFE)) \
                                                        -((GEXMNFE)+(W1MNFE)-(W2MNFE)))
#define GR1S3         0.25*(( (GEXMGMN)-(W1MGMN)+(W2MGMN)) \
                                                        -((GEXMGCO)+(W1MGCO)-(W2MGCO)) \
                                                        +((GEXMNCO)+(W1MNCO)-(W2MNCO)))
#define GR1S4         0.25*(( (GEXMGMN)-(W1MGMN)+(W2MGMN)) \
                                                        -((GEXMGNI)+(W1MGNI)-(W2MGNI)) \
                                                        +((GEXMNNI)+(W1MNNI)-(W2MNNI)))
#define GR2R2         -0.25*(( (GXMGFE)+(W1MGFE)+(W2MGFE)))
#define GR2R3          0.25*(( (GXFECO)+(W1FECO)+(W2FECO)) \
                             -((GXMGFE)+(W1MGFE)+(W2MGFE)) \
                             -((GXMGCO)+(W1MGCO)+(W2MGCO)))
#define GR2R4          0.25*(( (GXFENI)+(W1FENI)+(W2FENI)) \
                             -((GXMGFE)+(W1MGFE)+(W2MGFE)) \
                             -((GXMGNI)+(W1MGNI)+(W2MGNI)))
#define GR2R5             -0.25*((F_FE)+(GEXMGFE)+(GXMGFE) \
                   +2.0*(W2CAMG)-2.0*(W2CAFE)+2.0*(W2MGFE))
#define GR2S1         0.25*(( (GEXMGMN)+(W1MGMN)-(W2MGMN)) \
                                                        -((GEXMGFE)-(W1MGFE)+(W2MGFE)) \
                                                        +((GEXMNFE)-(W1MNFE)+(W2MNFE)))
#define GR2S2                       0.5*((W1MGFE)-(W2MGFE))
#define GR2S3         0.25*(( (GEXMGFE)-(W1MGFE)+(W2MGFE)) \
                                                        -((GEXMGCO)+(W1MGCO)-(W2MGCO)) \
                                                        +((GEXFECO)+(W1FECO)-(W2FECO)))
#define GR2S4         0.25*(( (GEXMGFE)-(W1MGFE)+(W2MGFE)) \
                                                        -((GEXMGNI)+(W1MGNI)-(W2MGNI)) \
                                                        +((GEXFENI)+(W1FENI)-(W2FENI)))
#define GR3R3          -0.25*(  (GXMGCO)+(W1MGCO)+(W2MGCO))
#define GR3R4          0.25*(( (GXCONI)+(W1CONI)+(W2CONI)) \
                             -((GXMGCO)+(W1MGCO)+(W2MGCO)) \
                             -((GXMGNI)+(W1MGNI)+(W2MGNI)))
#define GR3R5             -0.25*((F_CO)+(GEXMGCO)+(GXMGCO) \
                   +2.0*(W2CAMG)-2.0*(W2CACO)+2.0*(W2MGCO))
#define GR3S1         0.25*(( (GEXMGMN)+(W1MGMN)-(W2MGMN)) \
                                                        -((GEXMGCO)-(W1MGCO)+(W2MGCO)) \
                                                        +((GEXMNCO)-(W1MNCO)+(W2MNCO)))
#define GR3S2         0.25*(( (GEXMGFE)+(W1MGFE)-(W2MGFE)) \
                                                        -((GEXMGCO)-(W1MGCO)+(W2MGCO)) \
                                                        +((GEXFECO)-(W1FECO)+(W2FECO)))
#define GR3S3                      0.5*(-(W1MGCO)+(W2MGCO))
#define GR3S4          0.25*(( (GEXMGCO)+(W1MGCO)-(W2MGCO)) \
                             -((GEXMGNI)-(W1MGNI)+(W2MGNI)) \
                             +((GEXCONI)-(W1CONI)+(W2CONI)))
#define GR4R4           -0.25*(  (GXMGNI)+(W1MGNI)+(W2MGNI))
#define GR4R5             -0.25*((F_NI)+(GEXMGNI)+(GXMGNI) \
                   +2.0*(W2CAMG)-2.0*(W2CANI)+2.0*(W2MGNI))
#define GR4S1          0.25*(( (GEXMGMN)+(W1MGMN)-(W2MGMN)) \
                             -((GEXMGNI)-(W1MGNI)+(W2MGNI)) \
                             +((GEXMNNI)-(W1MNNI)+(W2MNNI)))
#define GR4S2          0.25*(( (GEXMGFE)+(W1MGFE)-(W2MGFE)) \
                             -((GEXMGNI)-(W1MGNI)+(W2MGNI)) \
                             +((GEXFENI)-(W1FENI)+(W2FENI)))
#define GR4S3          0.25*((-(GEXMGCO)+(W1MGCO)-(W2MGCO)) \
                             +((GEXMGNI)-(W1MGNI)+(W2MGNI)) \
                             -((GEXCONI)-(W1CONI)+(W2CONI)))
#define GR4S4                       0.5*(-(W1MGNI)+(W2MGNI))
#define GR5R5        -1.0*(W2CAMG)
#define GR5S1                0.25*((F_MN)+(GEXMGMN)+(GXMGMN) \
                     -2.0*(W2CAMG)+2.0*(W2CAMN)-2.0*(W2MGMN))
#define GR5S2                0.25*((F_FE)+(GEXMGFE)+(GXMGFE) \
                     -2.0*(W2CAMG)+2.0*(W2CAFE)-2.0*(W2MGFE))
#define GR5S3               -0.25*((F_CO)+(GEXMGCO)+(GXMGCO) \
                     -2.0*(W2CAMG)+2.0*(W2CACO)-2.0*(W2MGCO))
#define GR5S4               -0.25*((F_NI)+(GEXMGNI)+(GXMGNI) \
                     -2.0*(W2CAMG)+2.0*(W2CANI)-2.0*(W2MGNI))
#define GS1S1               0.25*((GXMGMN)-(W1MGMN)-(W2MGMN))
#define GS1S2            0.25*(( (GXMGMN)-(W1MGMN)-(W2MGMN)) \
                               +((GXMGFE)-(W1MGFE)-(W2MGFE)) \
                               -((GXMNFE)-(W1MNFE)-(W2MNFE)))
#define GS1S3            0.25*(-((GXMGMN)-(W1MGMN)-(W2MGMN)) \
                               -((GXMGCO)-(W1MGCO)-(W2MGCO)) \
                               +((GXMNCO)-(W1MNCO)-(W2MNCO)))
#define GS1S4            0.25*(-((GXMGMN)-(W1MGMN)-(W2MGMN)) \
                               -((GXMGNI)-(W1MGNI)-(W2MGNI)) \
                               +((GXMNNI)-(W1MNNI)-(W2MNNI)))
#define GS2S2             0.25*( (GXMGFE)-(W1MGFE)-(W2MGFE))
#define GS2S3            0.25*(-((GXMGFE)-(W1MGFE)-(W2MGFE)) \
                               -((GXMGCO)-(W1MGCO)-(W2MGCO)) \
                               +((GXFECO)-(W1FECO)-(W2FECO)))
#define GS2S4            0.25*(-((GXMGFE)-(W1MGFE)-(W2MGFE)) \
                               -((GXMGNI)-(W1MGNI)-(W2MGNI)) \
                               +((GXFENI)-(W1FENI)-(W2FENI)))
#define GS3S3             0.25*( (GXMGCO)-(W1MGCO)-(W2MGCO))
#define GS3S4            0.25*( ((GXMGCO)-(W1MGCO)-(W2MGCO)) \
                               +((GXMGNI)-(W1MGNI)-(W2MGNI)) \
                               -((GXCONI)-(W1CONI)-(W2CONI)))
#define GS4S4             0.25*( (GXMGNI)-(W1MGNI)-(W2MGNI))


#define V0  0.25*(     ((VXMNFE)+(WV1MNFE)+(WV2MNFE)) \
                                            +((VXMNCO)+(WV1MNCO)+(WV2MNCO)) \
                                            +((VXMNNI)+(WV1MNNI)+(WV2MNNI)) \
                                            +((VXFECO)+(WV1FECO)+(WV2FECO)) \
                                            +((VXFENI)+(WV1FENI)+(WV2FENI)) \
                                            +((VXCONI)+(WV1CONI)+(WV2CONI)) \
                                    -2.0*((VXMGMN)+(WV1MGMN)+(WV2MGMN)) \
                                    -2.0*((VXMGFE)+(WV1MGFE)+(WV2MGFE)) \
                                    -2.0*((VXMGCO)+(WV1MGCO)+(WV2MGCO)) \
                                    -2.0*((VXMGNI)+(WV1MGNI)+(WV2MGNI)))
#define VR1 0.25*(-3.0*((VXMGMN)+(WV1MGMN)+(WV2MGMN)) \
                                            -((VXMGFE)+(WV1MGFE)+(WV2MGFE)) \
                                            -((VXMGCO)+(WV1MGCO)+(WV2MGCO)) \
                                            -((VXMGNI)+(WV1MGNI)+(WV2MGNI)) \
                                            +((VXMNFE)+(WV1MNFE)+(WV2MNFE)) \
                                            +((VXMNCO)+(WV1MNCO)+(WV2MNCO)) \
                                            +((VXMNNI)+(WV1MNNI)+(WV2MNNI)))
#define VR2  0.25*(   -((VXMGMN)+(WV1MGMN)+(WV2MGMN)) \
                                    -3.0*((VXMGFE)+(WV1MGFE)+(WV2MGFE)) \
                                            -((VXMGCO)+(WV1MGCO)+(WV2MGCO)) \
                                            -((VXMGNI)+(WV1MGNI)+(WV2MGNI)) \
                                            +((VXMNFE)+(WV1MNFE)+(WV2MNFE)) \
                                            +((VXFECO)+(WV1FECO)+(WV2FECO)) \
                                            +((VXFENI)+(WV1FENI)+(WV2FENI)))
#define VR3  0.25*(   -((VXMGMN)+(WV1MGMN)+(WV2MGMN)) \
                                            -((VXMGFE)+(WV1MGFE)+(WV2MGFE)) \
                                    -3.0*((VXMGCO)+(WV1MGCO)+(WV2MGCO)) \
                                            -((VXMGNI)+(WV1MGNI)+(WV2MGNI)) \
                                            +((VXMNCO)+(WV1MNCO)+(WV2MNCO)) \
                                            +((VXFECO)+(WV1FECO)+(WV2FECO)) \
                                            +((VXCONI)+(WV1CONI)+(WV2CONI)))
#define VR4 0.25*(    -((VXMGMN)+(WV1MGMN)+(WV2MGMN)) \
                                            -((VXMGFE)+(WV1MGFE)+(WV2MGFE)) \
                                            -((VXMGCO)+(WV1MGCO)+(WV2MGCO)) \
                                    -3.0*((VXMGNI)+(WV1MGNI)+(WV2MGNI)) \
                                            +((VXMNNI)+(WV1MNNI)+(WV2MNNI)) \
                                            +((VXFENI)+(WV1FENI)+(WV2FENI)) \
                                            +((VXCONI)+(WV1CONI)+(WV2CONI)))
#define VR5        -(WV2CAMG) \
                                    -0.25*(((VEXMGMN)+(VXMGMN) \
                                                -2.0*(WV2CAMN)+2.0*(WV2MGMN)) \
                                                +((VEXMGFE)+(VXMGFE) \
                                                -2.0*(WV2CAFE)+2.0*(WV2MGFE)) \
                                                +((VEXMGCO)+(VXMGCO) \
                                                -2.0*(WV2CACO)+2.0*(WV2MGCO)) \
                                                +((VEXMGNI)+(VXMGNI) \
                                                -2.0*(WV2CANI)+2.0*(WV2MGNI)))
#define VS1   0.25*(( (VEXMGMN)+3.0*(WV1MGMN)-3.0*(WV2MGMN)) \
                                                        -((VEXMGFE)-(WV1MGFE)+(WV2MGFE)) \
                                                        -((VEXMGCO)-(WV1MGCO)+(WV2MGCO)) \
                                                        -((VEXMGNI)-(WV1MGNI)+(WV2MGNI)) \
                                                        +((VEXMNFE)-(WV1MNFE)+(WV2MNFE)) \
                                                        +((VEXMNCO)-(WV1MNCO)+(WV2MNCO)) \
                                                        +((VEXMNNI)-(WV1MNNI)+(WV2MNNI)))
#define VS2           0.25*((-(VEXMGMN)+(WV1MGMN)-(WV2MGMN)) \
                                        +((VEXMGFE)+3.0*(WV1MGFE)-3.0*(WV2MGFE)) \
                                                        -((VEXMGCO)-(WV1MGCO)+(WV2MGCO)) \
                                                        -((VEXMGNI)-(WV1MGNI)+(WV2MGNI)) \
                                                        -((VEXMNFE)+(WV1MNFE)-(WV2MNFE)) \
                                                        +((VEXFECO)-(WV1FECO)+(WV2FECO)) \
                                                        +((VEXFENI)-(WV1FENI)+(WV2FENI)))
#define VS3           0.25*(( (VEXMGMN)-(WV1MGMN)+(WV2MGMN)) \
                                                        +((VEXMGFE)-(WV1MGFE)+(WV2MGFE)) \
                                        -((VEXMGCO)+3.0*(WV1MGCO)-3.0*(WV2MGCO)) \
                                                        +((VEXMGNI)-(WV1MGNI)+(WV2MGNI)) \
                                                        +((VEXMNCO)+(WV1MNCO)-(WV2MNCO)) \
                                                        +((VEXFECO)+(WV1FECO)-(WV2FECO)) \
                                                        -((VEXCONI)-(WV1CONI)+(WV2CONI)))
#define VS4           0.25*(( (VEXMGMN)-(WV1MGMN)+(WV2MGMN)) \
                                                        +((VEXMGFE)-(WV1MGFE)+(WV2MGFE)) \
                                                        +((VEXMGCO)-(WV1MGCO)+(WV2MGCO)) \
                                        -((VEXMGNI)+3.0*(WV1MGNI)-3.0*(WV2MGNI)) \
                                                        +((VEXMNNI)+(WV1MNNI)-(WV2MNNI)) \
                                                        +((VEXFENI)+(WV1FENI)-(WV2FENI)) \
                                                        +((VEXCONI)+(WV1CONI)-(WV2CONI)))
#define VR1R1         -0.25*(  (VXMGMN)+(WV1MGMN)+(WV2MGMN))
#define VR1R2          0.25*(( (VXMNFE)+(WV1MNFE)+(WV2MNFE)) \
                             -((VXMGMN)+(WV1MGMN)+(WV2MGMN)) \
                             -(( VXMGFE)+(WV1MGFE)+(WV2MGFE)))
#define VR1R3          0.25*(( (VXMNCO)+(WV1MNCO)+(WV2MNCO)) \
                             -((VXMGMN)+(WV1MGMN)+(WV2MGMN)) \
                             -((VXMGCO)+(WV1MGCO)+(WV2MGCO)))
#define VR1R4          0.25*(( (VXMNNI)+(WV1MNNI)+(WV2MNNI)) \
                             -((VXMGMN)+(WV1MGMN)+(WV2MGMN)) \
                             -((VXMGNI)+(WV1MGNI)+(WV2MGNI)))
#define VR1R5             -0.25*((VEXMGMN)+(VXMGMN) \
                   +2.0*(WV2CAMG)-2.0*(WV2CAMN)+2.0*(WV2MGMN))
#define VR1S1                       0.5*((WV1MGMN)-(WV2MGMN))
#define VR1S2         0.25*((-(VEXMGMN)-(WV1MGMN)+(WV2MGMN)) \
                                                        +((VEXMGFE)+(WV1MGFE)-(WV2MGFE)) \
                                                        -((VEXMNFE)+(WV1MNFE)-(WV2MNFE)))
#define VR1S3         0.25*(( (VEXMGMN)-(WV1MGMN)+(WV2MGMN)) \
                                                        -((VEXMGCO)+(WV1MGCO)-(WV2MGCO)) \
                                                        +((VEXMNCO)+(WV1MNCO)-(WV2MNCO)))
#define VR1S4         0.25*(( (VEXMGMN)-(WV1MGMN)+(WV2MGMN)) \
                                                        -((VEXMGNI)+(WV1MGNI)-(WV2MGNI)) \
                                                        +((VEXMNNI)+(WV1MNNI)-(WV2MNNI)))
#define VR2R2         -0.25*(( (VXMGFE)+(WV1MGFE)+(WV2MGFE)))
#define VR2R3          0.25*(( (VXFECO)+(WV1FECO)+(WV2FECO)) \
                             -((VXMGFE)+(WV1MGFE)+(WV2MGFE)) \
                             -((VXMGCO)+(WV1MGCO)+(WV2MGCO)))
#define VR2R4          0.25*(( (VXFENI)+(WV1FENI)+(WV2FENI)) \
                             -((VXMGFE)+(WV1MGFE)+(WV2MGFE)) \
                             -((VXMGNI)+(WV1MGNI)+(WV2MGNI)))
#define VR2R5             -0.25*((VEXMGFE)+(VXMGFE) \
                   +2.0*(WV2CAMG)-2.0*(WV2CAFE)+2.0*(WV2MGFE))
#define VR2S1         0.25*(( (VEXMGMN)+(WV1MGMN)-(WV2MGMN)) \
                                                        -((VEXMGFE)-(WV1MGFE)+(WV2MGFE)) \
                                                        +((VEXMNFE)-(WV1MNFE)+(WV2MNFE)))
#define VR2S2                       0.5*((WV1MGFE)-(WV2MGFE))
#define VR2S3         0.25*(( (VEXMGFE)-(WV1MGFE)+(WV2MGFE)) \
                                                        -((VEXMGCO)+(WV1MGCO)-(WV2MGCO)) \
                                                        +((VEXFECO)+(WV1FECO)-(WV2FECO)))
#define VR2S4         0.25*(( (VEXMGFE)-(WV1MGFE)+(WV2MGFE)) \
                                                        -((VEXMGNI)+(WV1MGNI)-(WV2MGNI)) \
                                                        +((VEXFENI)+(WV1FENI)-(WV2FENI)))
#define VR3R3          -0.25*(  (VXMGCO)+(WV1MGCO)+(WV2MGCO))
#define VR3R4          0.25*(( (VXCONI)+(WV1CONI)+(WV2CONI)) \
                             -((VXMGCO)+(WV1MGCO)+(WV2MGCO)) \
                             -((VXMGNI)+(WV1MGNI)+(WV2MGNI)))
#define VR3R5             -0.25*((VEXMGCO)+(VXMGCO) \
                   +2.0*(WV2CAMG)-2.0*(WV2CACO)+2.0*(WV2MGCO))
#define VR3S1         0.25*(( (VEXMGMN)+(WV1MGMN)-(WV2MGMN)) \
                                                        -((VEXMGCO)-(WV1MGCO)+(WV2MGCO)) \
                                                        +((VEXMNCO)-(WV1MNCO)+(WV2MNCO)))
#define VR3S2         0.25*(( (VEXMGFE)+(WV1MGFE)-(WV2MGFE)) \
                                                        -((VEXMGCO)-(WV1MGCO)+(WV2MGCO)) \
                                                        +((VEXFECO)-(WV1FECO)+(WV2FECO)))
#define VR3S3                      0.5*(-(WV1MGCO)+(WV2MGCO))
#define VR3S4          0.25*(( (VEXMGCO)+(WV1MGCO)-(WV2MGCO)) \
                             -((VEXMGNI)-(WV1MGNI)+(WV2MGNI)) \
                             +((VEXCONI)-(WV1CONI)+(WV2CONI)))
#define VR4R4           -0.25*(  (VXMGNI)+(WV1MGNI)+(WV2MGNI))
#define VR4R5              -0.25*((VEXMGNI)+(VXMGNI) \
                                        +2.0*(WV2CAMG)-2.0*(WV2CANI)+2.0*(WV2MGNI))
#define VR4S1          0.25*(( (VEXMGMN)+(WV1MGMN)-(WV2MGMN)) \
                             -((VEXMGNI)-(WV1MGNI)+(WV2MGNI)) \
                             +((VEXMNNI)-(WV1MNNI)+(WV2MNNI)))
#define VR4S2          0.25*(( (VEXMGFE)+(WV1MGFE)-(WV2MGFE)) \
                             -((VEXMGNI)-(WV1MGNI)+(WV2MGNI)) \
                             +((VEXFENI)-(WV1FENI)+(WV2FENI)))
#define VR4S3          0.25*((-(VEXMGCO)+(WV1MGCO)-(WV2MGCO)) \
                             +((VEXMGNI)-(WV1MGNI)+(WV2MGNI)) \
                             -((VEXCONI)-(WV1CONI)+(WV2CONI)))
#define VR4S4                       0.5*(-(WV1MGNI)+(WV2MGNI))
#define VR5R5                                      -(WV2CAMG)
#define VR5S1                0.25*((VEXMGMN)+(VXMGMN) \
                     -2.0*(WV2CAMG)+2.0*(WV2CAMN)-2.0*(WV2MGMN))
#define VR5S2                0.25*((VEXMGFE)+(VXMGFE) \
                     -2.0*(WV2CAMG)+2.0*(WV2CAFE)-2.0*(WV2MGFE))
#define VR5S3               -0.25*((VEXMGCO)+(VXMGCO) \
                     -2.0*(WV2CAMG)+2.0*(WV2CACO)-2.0*(WV2MGCO))
#define VR5S4               -0.25*((VEXMGNI)+(VXMGNI) \
                     -2.0*(WV2CAMG)+2.0*(WV2CANI)-2.0*(WV2MGNI))
#define VS1S1               0.25*((VXMGMN)-(WV1MGMN)-(WV2MGMN))
#define VS1S2            0.25*(( (VXMGMN)-(WV1MGMN)-(WV2MGMN)) \
                               +((VXMGFE)-(WV1MGFE)-(WV2MGFE)) \
                               -((VXMNFE)-(WV1MNFE)-(WV2MNFE)))
#define VS1S3            0.25*(-((VXMGMN)-(WV1MGMN)-(WV2MGMN)) \
                               -((VXMGCO)-(WV1MGCO)-(WV2MGCO)) \
                               +((VXMNCO)-(WV1MNCO)-(WV2MNCO)))
#define VS1S4            0.25*(-((VXMGMN)-(WV1MGMN)-(WV2MGMN)) \
                               -((VXMGNI)-(WV1MGNI)-(WV2MGNI)) \
                               +((VXMNNI)-(WV1MNNI)-(WV2MNNI)))
#define VS2S2             0.25*(  (VXMGFE)-(WV1MGFE)-(WV2MGFE))
#define VS2S3            0.25*(-((VXMGFE)-(WV1MGFE)-(WV2MGFE)) \
                               -((VXMGCO)-(WV1MGCO)-(WV2MGCO)) \
                               +((VXFECO)-(WV1FECO)-(WV2FECO)))
#define VS2S4            0.25*(-((VXMGFE)-(WV1MGFE)-(WV2MGFE)) \
                               -((VXMGNI)-(WV1MGNI)-(WV2MGNI)) \
                               +((VXFENI)-(WV1FENI)-(WV2FENI)))
#define VS3S3             0.25*(  (VXMGCO)-(WV1MGCO)-(WV2MGCO))
#define VS3S4            0.25*( ((VXMGCO)-(WV1MGCO)-(WV2MGCO)) \
                               +((VXMGNI)-(WV1MGNI)-(WV2MGNI)) \
                               -((VXCONI)-(WV1CONI)-(WV2CONI)))
#define VS4S4             0.25*(  (VXMGNI)-(WV1MGNI)-(WV2MGNI))

/*
 * Global (to this file): variables
 */

#define R  8.3143
#define NR         5    /* Five independent composition variables */
#define NS         4    /* Four ordering parameters              */
#define NA         6    /* Six endmember compositions            */
                                                /* site mole fractions */
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
static MTHREAD_KEY_T ptToVKey;
static MTHREAD_KEY_T ptToSKey;
static MTHREAD_KEY_T indexD2gds2Key;

static void freeNSvector(void *NSarray) {
    gsl_vector_free((gsl_vector *) NSarray);
}

static void freeNSmatrix(void *NSarray) {
    gsl_matrix_free((gsl_matrix *) NSarray);
}

static void freeIndexD2gds2(void *indexD2gds2) {
    gsl_permutation_free((gsl_permutation *) indexD2gds2);
}

static void threadOInit(void) {
    MTHREAD_KEY_CREATE(&tOldKey,       free);
    MTHREAD_KEY_CREATE(&pOldKey,       free);
    MTHREAD_KEY_CREATE(&rOldKey,       freeNSvector);
    MTHREAD_KEY_CREATE(&sOldKey,       freeNSvector);
    MTHREAD_KEY_CREATE(&ptToD2gds2Key, freeNSmatrix);
    MTHREAD_KEY_CREATE(&d2gds2Key,     free);
    MTHREAD_KEY_CREATE(&ptToVKey,      freeNSmatrix);
    MTHREAD_KEY_CREATE(&ptToSKey,      freeNSvector);
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

static gsl_matrix *getPtToV() {
    gsl_matrix *ptToVPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    ptToVPt = (gsl_matrix *) MTHREAD_GETSPECIFIC(ptToVKey);
    if (ptToVPt == NULL) {
        ptToVPt  = gsl_matrix_alloc((size_t) NS, (size_t) NS);
        gsl_matrix_set_zero(ptToVPt);
        MTHREAD_SETSPECIFIC(ptToVKey, (void *) ptToVPt);
    }
    return ptToVPt;
}

static gsl_vector *getPtToS() {
    gsl_vector *ptToSPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    ptToSPt = (gsl_vector *) MTHREAD_GETSPECIFIC(ptToSKey);
    if (ptToSPt == NULL) {
        ptToSPt  = gsl_vector_alloc((size_t) NS);
        gsl_vector_set_zero(ptToSPt);
        MTHREAD_SETSPECIFIC(ptToSKey, (void *) ptToSPt);
    }
    return ptToSPt;
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
    double xm1mg, xm1mn, xm1fe, xm1co, xm1ni, xm2mg, xm2mn, xm2fe, xm2co, xm2ni, xm2ca;

#define XM1MG  0
#define XM1MN  1
#define XM1FE  2
#define XM1CO  3
#define XM1NI  4
#define XM2MG  5
#define XM2MN  6
#define XM2FE  7
#define XM2CO  8
#define XM2NI  9
#define XM2CA 10

#define NX    11

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
    xm1mg = getX(XM1MG); \
    xm1mn = getX(XM1MN); \
    xm1fe = getX(XM1FE); \
    xm1co = getX(XM1CO); \
    xm1ni = getX(XM1NI); \
    xm2mg = getX(XM2MG); \
    xm2mn = getX(XM2MN); \
    xm2fe = getX(XM2FE); \
    xm2co = getX(XM2CO); \
    xm2ni = getX(XM2NI); \
    xm2ca = getX(XM2CA);

#define SET_SITE_FRACTIONS \
    setX(XM1MG, xm1mg); \
    setX(XM1MN, xm1mn); \
    setX(XM1FE, xm1fe); \
    setX(XM1CO, xm1co); \
    setX(XM1NI, xm1ni); \
    setX(XM2MG, xm2mg); \
    setX(XM2MN, xm2mn); \
    setX(XM2FE, xm2fe); \
    setX(XM2CO, xm2co); \
    setX(XM2NI, xm2ni); \
    setX(XM2CA, xm2ca);


/*
 * "Darken Equation" coefficients -
 * Global (to this file): activity definitions and component transforms
 *    The function conOlv defines the conversion from m[i], to r[j]
 */
                   /* Order R1,R2,R3,R4,R5 */
#define FR0(i)     (i == 0) ? 1.0 - r[0] : -(1.0 + r[0])
#define FR1(i)     (i == 1) ? 1.0 - r[1] : -(1.0 + r[1])
#define FR2(i)     (i == 2) ? 1.0 - r[2] : -(1.0 + r[2])
#define FR3(i)     (i == 3) ? 1.0 - r[3] : -(1.0 + r[3])
#define FR4(i)     (i == 4) ? 1.0 - r[4] : -r[4]

                                        /* Order: S1, S2, S3, S4 */
#define GGS0(i)     -s[0]
#define GGS1(i)     -s[1]
#define GGS2(i)     -s[2]
#define GGS3(i)     -s[3]

#define DFR0DR0(i) - 1.0
#define DFR1DR1(i) - 1.0
#define DFR2DR2(i) - 1.0
#define DFR3DR3(i) - 1.0
#define DFR4DR4(i) - 1.0

#define DGS0DS0(i) - 1.0
#define DGS1DS1(i) - 1.0
#define DGS2DS2(i) - 1.0
#define DGS3DS3(i) - 1.0

/*
 * Global (to this file): derivative definitions
 */

#define S -R*( xm1mn*log(xm1mn) + xm2mn*log(xm2mn) + xm1fe*log(xm1fe)+ \
               xm2fe*log(xm2fe) + xm1co*log(xm1co) + xm2co*log(xm2co)+ \
               xm1ni*log(xm1ni) + xm2ni*log(xm2ni) + xm2ca*log(xm2ca)+ \
               xm1mg*log(xm1mg) + xm2mg*log(xm2mg))

/*  enthalpy here is enthalpy at P of interest */
#define H     (G0) + \
             (GR1)*r[0] + (GR2)*r[1] + (GR3)*r[2] + (GR4)*r[3] + \
             (GR5)*r[4] + (GS1)*s[0] + (GS2)*s[1] + (GS3)*s[2] + \
             (GS4)*s[3] + \
             (GR1R1)*r[0]*r[0] + (GR1R2)*r[0]*r[1] + (GR1R3)*r[0]*r[2] + \
             (GR1R4)*r[0]*r[3] + (GR1R5)*r[0]*r[4] + (GR1S1)*r[0]*s[0] + \
             (GR1S2)*r[0]*s[1] + (GR1S3)*r[0]*s[2] + (GR1S4)*r[0]*s[3] + \
             (GR2R2)*r[1]*r[1] + (GR2R3)*r[1]*r[2] + (GR2R4)*r[1]*r[3] + \
             (GR2R5)*r[1]*r[4] + (GR2S1)*r[1]*s[0] + (GR2S2)*r[1]*s[1] + \
             (GR2S3)*r[1]*s[2] + (GR2S4)*r[1]*s[3] + (GR3R3)*r[2]*r[2] + \
             (GR3R4)*r[2]*r[3] + (GR3R5)*r[2]*r[4] + (GR3S1)*r[2]*s[0] + \
             (GR3S2)*r[2]*s[1] + (GR3S3)*r[2]*s[2] + (GR3S4)*r[2]*s[3] + \
             (GR4R4)*r[3]*r[3] + (GR4R5)*r[3]*r[4] + (GR4S1)*r[3]*s[0] + \
             (GR4S2)*r[3]*s[1] + (GR4S3)*r[3]*s[2] + (GR4S4)*r[3]*s[3] + \
             (GR5R5)*r[4]*r[4] + (GR5S1)*r[4]*s[0] + (GR5S2)*r[4]*s[1] + \
             (GR5S3)*r[4]*s[2] + (GR5S4)*r[4]*s[3] + (GS1S1)*s[0]*s[0] + \
             (GS1S2)*s[0]*s[1] + (GS1S3)*s[0]*s[2] + (GS1S4)*s[0]*s[3] + \
             (GS2S2)*s[1]*s[1] + (GS2S3)*s[1]*s[2] + (GS2S4)*s[1]*s[3] + \
             (GS3S3)*s[2]*s[2] + (GS3S4)*s[2]*s[3] + (GS4S4)*s[3]*s[3]
#define V    (V0) + \
             (VR1)*r[0] + (VR2)*r[1] + (VR3)*r[2] + (VR4)*r[3] + \
             (VR5)*r[4] + (VS1)*s[0] + (VS2)*s[1] + (VS3)*s[2] + \
             (VS4)*s[3] + \
             (VR1R1)*r[0]*r[0] + (VR1R2)*r[0]*r[1] + (VR1R3)*r[0]*r[2] + \
             (VR1R4)*r[0]*r[3] + (VR1R5)*r[0]*r[4] + (VR1S1)*r[0]*s[0] + \
             (VR1S2)*r[0]*s[1] + (VR1S3)*r[0]*s[2] + (VR1S4)*r[0]*s[3] + \
             (VR2R2)*r[1]*r[1] + (VR2R3)*r[1]*r[2] + (VR2R4)*r[1]*r[3] + \
             (VR2R5)*r[1]*r[4] + (VR2S1)*r[1]*s[0] + (VR2S2)*r[1]*s[1] + \
             (VR2S3)*r[1]*s[2] + (VR2S4)*r[1]*s[3] + (VR3R3)*r[2]*r[2] + \
             (VR3R4)*r[2]*r[3] + (VR3R5)*r[2]*r[4] + (VR3S1)*r[2]*s[0] + \
             (VR3S2)*r[2]*s[1] + (VR3S3)*r[2]*s[2] + (VR3S4)*r[2]*s[3] + \
             (VR4R4)*r[3]*r[3] + (VR4R5)*r[3]*r[4] + (VR4S1)*r[3]*s[0] + \
             (VR4S2)*r[3]*s[1] + (VR4S3)*r[3]*s[2] + (VR4S4)*r[3]*s[3] + \
             (VR5R5)*r[4]*r[4] + (VR5S1)*r[4]*s[0] + (VR5S2)*r[4]*s[1] + \
             (VR5S3)*r[4]*s[2] + (VR5S4)*r[4]*s[3] + (VS1S1)*s[0]*s[0] + \
             (VS1S2)*s[0]*s[1] + (VS1S3)*s[0]*s[2] + (VS1S4)*s[0]*s[3] + \
             (VS2S2)*s[1]*s[1] + (VS2S3)*s[1]*s[2] + (VS2S4)*s[1]*s[3] + \
             (VS3S3)*s[2]*s[2] + (VS3S4)*s[2]*s[3] + (VS4S4)*s[3]*s[3]

#define G    (H) - t*(S)

/*----------------------------------------------------------------------------*/

#define DGDR0  (GR1) + 2.0*(GR1R1)*r[0] + \
                (GR1R2)*r[1] + (GR1R3)*r[2] + (GR1R4)*r[3] + (GR1R5)*r[4] + \
                (GR1S1)*s[0] + (GR1S2)*s[1] + (GR1S3)*s[2] + (GR1S4)*s[3] + \
                         0.5*R*t*(log(xm1mn*xm2mn/xm1mg/xm2mg))
#define DGDR1  (GR2) + 2.0*(GR2R2)*r[1] + \
                (GR1R2)*r[0] + (GR2R3)*r[2] + (GR2R4)*r[3] + (GR2R5)*r[4] + \
                (GR2S1)*s[0] + (GR2S2)*s[1] + (GR2S3)*s[2] + (GR2S4)*s[3] + \
                         0.5*R*t*(log(xm1fe*xm2fe/xm1mg/xm2mg))
#define DGDR2  (GR3) + 2.0*(GR3R3)*r[2] + \
                (GR1R3)*r[0] + (GR2R3)*r[1] + (GR3R4)*r[3] +(GR3R5)*r[4] + \
                (GR3S1)*s[0] + (GR3S2)*s[1] + (GR3S3)*s[2] +(GR3S4)*s[3] + \
                         0.5*R*t*(log(xm1co*xm2co/xm1mg/xm2mg))
#define DGDR3  (GR4) + 2.0*(GR4R4)*r[3] + \
                (GR1R4)*r[0] + (GR2R4)*r[1] + (GR3R4)*r[2] +(GR4R5)*r[4] + \
                (GR4S1)*s[0] + (GR4S2)*s[1] + (GR4S3)*s[2] +(GR4S4)*s[3] + \
                         0.5*R*t*(log(xm1ni*xm2ni/xm1mg/xm2mg))
#define DGDR4  (GR5) + 2.0*(GR5R5)*r[4] + \
                (GR1R5)*r[0] + (GR2R5)*r[1] + (GR3R5)*r[2] + (GR4R5)*r[3] + \
                (GR5S1)*s[0] + (GR5S2)*s[1] + (GR5S3)*s[2] + (GR5S4)*s[3] + \
                                                        R*t*(log(xm2ca/xm2mg))
#define DGDS0  (GS1) + 2.0*(GS1S1)*s[0] + \
                (GR1S1)*r[0] + (GR2S1)*r[1] + (GR3S1)*r[2] + (GR4S1)*r[3] + \
                (GR5S1)*r[4] + (GS1S2)*s[1] + (GS1S3)*s[2] + (GS1S4)*s[3] + \
                         0.5*R*t*(log(xm2mn*xm1mg/xm1mn/xm2mg))
#define DGDS1  (GS2) + 2.0*(GS2S2)*s[1] + \
                (GR1S2)*r[0] + (GR2S2)*r[1] + (GR3S2)*r[2] + (GR4S2)*r[3] + \
                (GR5S2)*r[4] + (GS1S2)*s[0] + (GS2S3)*s[2] + (GS2S4)*s[3] + \
                         0.5*R*t*(log(xm2fe*xm1mg/xm1fe/xm2mg))
#define DGDS2  (GS3) + 2.0*(GS3S3)*s[2] + \
                (GR1S3)*r[0] + (GR2S3)*r[1] + (GR3S3)*r[2] + (GR4S3)*r[3] + \
                (GR5S3)*r[4] + (GS1S3)*s[0] + (GS2S3)*s[1] + (GS3S4)*s[3] + \
                         0.5*R*t*(log(xm1co*xm2mg/xm2co/xm1mg))
#define DGDS3  (GS4) + 2.0*(GS4S4)*s[3] + \
                (GR1S4)*r[0] + (GR2S4)*r[1] + (GR3S4)*r[2] + (GR4S4)*r[3] + \
                (GR5S4)*r[4] + (GS1S4)*s[0] + (GS2S4)*s[1] + (GS3S4)*s[2] + \
                         0.5*R*t*(log(xm1ni*xm2mg/xm2ni/xm1mg))
#define DGDT  (S)
#define DGDP  (V)

/*----------------------------------------------------------------------------*/

#define D2GDR0R0 2.0*(GR1R1) + 0.25*R*t*(\
                 1.0/xm1mn + 1.0/xm2mn + 1.0/xm1mg + 1.0/xm2mg)
#define D2GDR0R1 (GR1R2) + 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDR0R2 (GR1R3) + 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDR0R3 (GR1R4) + 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDR0R4 (GR1R5) + 0.50*R*t*(1.0/xm2mg)
#define D2GDR0S0 (GR1S1) +  0.25*R*t*( \
                 -1.0/xm1mn + 1.0/xm2mn - 1.0/xm1mg + 1.0/xm2mg)
#define D2GDR0S1 (GR1S2) + 0.25*R*t*(- 1.0/xm1mg + 1.0/xm2mg)
#define D2GDR0S2 (GR1S3) + 0.25*R*t*(  1.0/xm1mg - 1.0/xm2mg)
#define D2GDR0S3 (GR1S4) + 0.25*R*t*(  1.0/xm1mg - 1.0/xm2mg)
#define D2GDR0DT 0.5*R*(log(xm1mn*xm2mn/xm1mg/xm2mg))
#define D2GDR0DP (VR1) + 2.0*(VR1R1)*r[0] + \
                (VR1R2)*r[1] + (VR1R3)*r[2] + (VR1R4)*r[3] + (VR1R5)*r[4] + \
                (VR1S1)*s[0] + (VR1S2)*s[1] + (VR1S3)*s[2] + (VR1S4)*s[3]

#define D2GDR1R1 2.0*(GR2R2) + 0.25*R*t* \
                                (1.0/xm1fe + 1.0/xm2fe + 1.0/xm1mg + 1.0/xm2mg)
#define D2GDR1R2 (GR2R3) + 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDR1R3 (GR2R4) + 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDR1R4 (GR2R5) + 0.5*R*t*(1.0/xm2mg)
#define D2GDR1S0 (GR2S1) + 0.25*R*t*(1.0/xm2mg - 1.0/xm1mg)
#define D2GDR1S1 (GR2S2) + 0.25*R*t*( \
                 -1.0/xm1fe + 1.0/xm2fe - 1.0/xm1mg + 1.0/xm2mg)
#define D2GDR1S2 (GR2S3) + 0.25*R*t*(  1.0/xm1mg - 1.0/xm2mg)
#define D2GDR1S3 (GR2S4) + 0.25*R*t*(  1.0/xm1mg - 1.0/xm2mg)
#define D2GDR1DT 0.5*R*(log(xm1fe*xm2fe/xm1mg/xm2mg))
#define D2GDR1DP (VR2) + 2.0*(VR2R2)*r[1] + \
                (VR1R2)*r[0] + (VR2R3)*r[2] + (VR2R4)*r[3] + (VR2R5)*r[4] + \
                (VR2S1)*s[0] + (VR2S2)*s[1] + (VR2S3)*s[2] + (VR2S4)*s[3]

#define D2GDR2R2 2.0*(GR3R3) + 0.25*R*t*( \
                 1.0/xm1co + 1.0/xm2co + 1.0/xm1mg + 1.0/xm2mg)
#define D2GDR2R3 (GR3R4) + 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDR2R4 (GR3R5) + 0.5*R*t*(1.0/xm2mg)
#define D2GDR2S0 (GR3S1) + 0.25*R*t*(1.0/xm2mg-1.0/xm1mg)
#define D2GDR2S1 (GR3S2) + 0.25*R*t*(1.0/xm2mg-1.0/xm1mg)
#define D2GDR2S2 (GR3S3) + 0.25*R*t*( \
                 1.0/xm1co - 1.0/xm2co + 1.0/xm1mg - 1.0/xm2mg)
#define D2GDR2S3 (GR3S4) + 0.25*R*t*(1.0/xm1mg-1.0/xm2mg)
#define D2GDR2DT 0.5*R*(log(xm1co*xm2co/xm1mg/xm2mg))
#define D2GDR2DP (VR3) + 2.0*(VR3R3)*r[2] + \
                (VR1R3)*r[0] + (VR2R3)*r[1] + (VR3R4)*r[3] +(VR3R5)*r[4] + \
                (VR3S1)*s[0] + (VR3S2)*s[1] + (VR3S3)*s[2] +(VR3S4)*s[3]

#define D2GDR3R3 2.0*(GR4R4) + 0.25*R*t*( \
                 1.0/xm1mg + 1.0/xm2mg + 1.0/xm1ni + 1.0/xm2ni)
#define D2GDR3R4 (GR4R5) + 0.5*R*t*(1.0/xm2mg)
#define D2GDR3S0 (GR4S1) + 0.25*R*t*(1.0/xm2mg-1.0/xm1mg)
#define D2GDR3S1 (GR4S2) + 0.25*R*t*(1.0/xm2mg-1.0/xm1mg)
#define D2GDR3S2 (GR4S3) + 0.25*R*t*(1.0/xm1mg-1.0/xm2mg)
#define D2GDR3S3 (GR4S4) + 0.25*R*t*( \
                 1.0/xm1ni - 1.0/xm2ni + 1.0/xm1mg - 1.0/xm2mg)
#define D2GDR3DT 0.5*R*(log(xm1ni*xm2ni/xm1mg/xm2mg))
#define D2GDR3DP (VR4) + 2.0*(VR4R4)*r[3] + \
                (VR1R4)*r[0] + (VR2R4)*r[1] + (VR3R4)*r[2] +(VR4R5)*r[4] + \
                (VR4S1)*s[0] + (VR4S2)*s[1] + (VR4S3)*s[2] +(VR4S4)*s[3]

#define D2GDR4R4 2.0*(GR5R5) + R*t*(1.0/xm2mg + 1.0/xm2ca)
#define D2GDR4S0 (GR5S1) + 0.5*R*t*(1.0/xm2mg)
#define D2GDR4S1 (GR5S2) + 0.5*R*t*(1.0/xm2mg)
#define D2GDR4S2 (GR5S3) - 0.5*R*t*(1.0/xm2mg)
#define D2GDR4S3 (GR5S4) - 0.5*R*t*(1.0/xm2mg)
#define D2GDR4DT R*log(xm2ca/xm2mg)
#define D2GDR4DP (VR5) + 2.0*(VR5R5)*r[4] + \
                (VR1R5)*r[0] + (VR2R5)*r[1] + (VR3R5)*r[2] + (VR4R5)*r[3] + \
                (VR5S1)*s[0] + (VR5S2)*s[1] + (VR5S3)*s[2] + (VR5S4)*s[3]

#define D2GDS0S0 2.0*(GS1S1) + \
                 0.25*R*t*(1.0/xm1mn + 1.0/xm2mn + 1.0/xm1mg + 1.0/xm2mg)
#define D2GDS0S1 (GS1S2) + 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDS0S2 (GS1S3) - 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDS0S3 (GS1S4) - 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDS0DT 0.5*R*(log(xm2mn*xm1mg/xm1mn/xm2mg))
#define D2GDS0DP (VS1) + 2.0*(VS1S1)*s[0] + \
                (VR1S1)*r[0] + (VR2S1)*r[1] + (VR3S1)*r[2] + (VR4S1)*r[3] + \
                (VR5S1)*r[4] + (VS1S2)*s[1] + (VS1S3)*s[2] + (VS1S4)*s[3]

#define D2GDS1S1 2.0*(GS2S2) + \
                 0.25*R*t*(1.0/xm1fe + 1.0/xm2fe + 1.0/xm1mg +1.0/xm2mg)
#define D2GDS1S2 (GS2S3) - 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDS1S3 (GS2S4) - 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDS1DT 0.5*R*(log(xm2fe*xm1mg/xm1fe/xm2mg))
#define D2GDS1DP (VS2) + 2.0*(VS2S2)*s[1] + \
                (VR1S2)*r[0] + (VR2S2)*r[1] + (VR3S2)*r[2] + (VR4S2)*r[3] + \
                (VR5S2)*r[4] + (VS1S2)*s[0] + (VS2S3)*s[2] + (VS2S4)*s[3]

#define D2GDS2S2 2.0*(GS3S3) + \
                0.25*R*t*(1.0/xm1co + 1.0/xm2co + 1.0/xm1mg +1.0/xm2mg)
#define D2GDS2S3 (GS3S4) + 0.25*R*t*(1.0/xm1mg + 1.0/xm2mg)
#define D2GDS2DT 0.5*R*(log(xm1co*xm2mg/xm2co/xm1mg))
#define D2GDS2DP (VS3) + 2.0*(VS3S3)*s[2] + \
                (VR1S3)*r[0] + (VR2S3)*r[1] + (VR3S3)*r[2] + (VR4S3)*r[3] + \
                (VR5S3)*r[4] + (VS1S3)*s[0] + (VS2S3)*s[1] + (VS3S4)*s[3]

#define D2GDS3S3 2.0*(GS4S4) + \
                0.25*R*t*(1.0/xm1ni + 1.0/xm2ni + 1.0/xm1mg +1.0/xm2mg)
#define D2GDS3DT 0.5*R*(log(xm1ni*xm2mg/xm2ni/xm1mg))
#define D2GDS3DP (VS4) + 2.0*(VS4S4)*s[3] + \
                (VR1S4)*r[0] + (VR2S4)*r[1] + (VR3S4)*r[2] + (VR4S4)*r[3] + \
                (VR5S4)*r[4] + (VS1S4)*s[0] + (VS2S4)*s[1] + (VS3S4)*s[2]

#define D2GDT2   0.0
#define D2GDTDP  0.0
#define D2GDP2   0.0

/*----------------------------------------------------------------------------*/
#define D3GDR0R0R0   0.125*R*t*(-1.0/SQUARE(xm1mn)+1.0/SQUARE(xm2mg) \
                       -1.0/SQUARE(xm2mn)+1.0/SQUARE(xm1mg))
#define D3GDR0R0R1   0.125*R*t*(1.0/SQUARE(xm1mg)+1.0/SQUARE(xm2mg))
#define D3GDR0R0S0   0.125*R*t*(1.0/SQUARE(xm1mn)+1.0/SQUARE(xm2mg) \
                       -1.0/SQUARE(xm2mn)-1.0/SQUARE(xm1mg))
#define D3GDR0R0S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR0R0S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR0R0S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR0R0DT   0.25*R*( \
                     1.0/xm1mn + 1.0/xm2mn + 1.0/xm1mg + 1.0/xm2mg)
#define D3GDR0R0DP   2.0*(VR1R1)

#define D3GDR0R1S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR0R1S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR0R1S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR0R1S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR0R1DT   0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDR0R1DP   (VR1R2)

#define D3GDR0R2S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR0R2S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR0R2S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR0R2S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR0R2DT   0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDR0R2DP   (VR1R3)

#define D3GDR0R3S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR0R3S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR0R3S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR0R3S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR0R3DT   0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDR0R3DP   (VR1R4)

#define D3GDR0R4S0   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR0R4S1   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR0R4S2  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR0R4S3  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR0R4DT   0.5 *R*(1.0/xm2mg)
#define D3GDR0R4DP   (VR1R5)

#define D3GDR0S0S0   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg) \
                               -1.0/SQUARE(xm2mn)-1.0/SQUARE(xm1mn))
#define D3GDR0S0S1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S0S2  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S0S3  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S0DT   0.25*R*( \
                     -1.0/xm1mn + 1.0/xm2mn - 1.0/xm1mg + 1.0/xm2mg)
#define D3GDR0S0DP   (VR1S1)

#define D3GDR0S1S1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S1S2  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S1S3  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S1DT   0.25*R*(1.0/xm2mg-1.0/xm1mg)
#define D3GDR0S1DP   (VR1S2)

#define D3GDR0S2S2   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S2S3   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S2DT   0.25*R*(1.0/xm1mg-1.0/xm2mg)
#define D3GDR0S2DP   (VR1S3)

#define D3GDR0S3S3   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR0S3DT   0.25*R*(  1.0/xm1mg - 1.0/xm2mg)

#define D3GDR0S3DP   (VR1S4)

#define D3GDR1R1R1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg) \
                               -1.0/SQUARE(xm1fe)-1.0/SQUARE(xm2fe))
#define D3GDR1R1S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR1R1S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg) \
                               +1.0/SQUARE(xm1fe)-1.0/SQUARE(xm2fe))
#define D3GDR1R1S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR1R1S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR1R1DT   0.25*R*( \
                     1.0/xm1fe + 1.0/xm2fe + 1.0/xm1mg + 1.0/xm2mg)
#define D3GDR1R1DP   2.0*(VR2R2)

#define D3GDR1R2S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR1R2S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR1R2S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR1R2S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR1R2DT   0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDR1R2DP   (VR2R3)

#define D3GDR1R3S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR1R3S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR1R3S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR1R3S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR1R3DT   0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDR1R3DP   (VR2R4)

#define D3GDR1R4S0   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR1R4S1   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR1R4S2  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR1R4S3  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR1R4DT   0.5*R* (1.0/xm2mg)
#define D3GDR1R4DP   (VR2R5)

#define D3GDR1S0S0   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S0S1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S0S2  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S0S3  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S0DT   0.25*R*(1.0/xm2mg - 1.0/xm1mg)
#define D3GDR1S0DP   (VR2S1)

#define D3GDR1S1S1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg) \
                               -1.0/SQUARE(xm2fe)-1.0/SQUARE(xm1fe))
#define D3GDR1S1S2  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S1S3  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S1DT   0.25*R*( \
                     -1.0/xm1fe + 1.0/xm2fe - 1.0/xm1mg + 1.0/xm2mg)
#define D3GDR1S1DP   (VR2S2)

#define D3GDR1S2S2   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S2S3   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S2DT   0.25*R*(  1.0/xm1mg - 1.0/xm2mg)
#define D3GDR1S2DP   (VR2S3)

#define D3GDR1S3S3   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR1S3DT   0.25*R*(  1.0/xm1mg - 1.0/xm2mg)
#define D3GDR1S3DP   (VR2S4)

#define D3GDR2R2R2   0.125*R*t*(1.0/SQUARE(xm1mg)+1.0/SQUARE(xm2mg) \
                               -1.0/SQUARE(xm2co)-1.0/SQUARE(xm1co))
#define D3GDR2R2S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR2R2S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR2R2S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg) \
                               +1.0/SQUARE(xm2co)-1.0/SQUARE(xm1co))
#define D3GDR2R2S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR2R2DT   0.25*R*( \
                     1.0/xm1co + 1.0/xm2co + 1.0/xm1mg + 1.0/xm2mg)
#define D3GDR2R2DP   2.0*(VR3R3)

#define D3GDR2R3S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR2R3S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR2R3S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR2R3S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR2R3DT   0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDR2R3DP   (VR3R4)

#define D3GDR2R4S0   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR2R4S1   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR2R4S2  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR2R4S3  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR2R4DT   0.5*R* (1.0/xm2mg)
#define D3GDR2R4DP   (VR3R5)

#define D3GDR2S0S0   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S0S1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S0S2  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S0S3  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S0DT   0.25*R*(1.0/xm2mg-1.0/xm1mg)
#define D3GDR2S0DP   (VR3S1)

#define D3GDR2S1S1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S1S2  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S1S3  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S1DT   0.25*R*(1.0/xm2mg-1.0/xm1mg)
#define D3GDR2S1DP   (VR3S2)

#define D3GDR2S2S2   0.125*R*t*(1.0/SQUARE(xm1mg)+1.0/SQUARE(xm2mg) \
                               -1.0/SQUARE(xm1co)-1.0/SQUARE(xm2co))
#define D3GDR2S2S3   0.125*R*t*(1.0/SQUARE(xm1mg)+1.0/SQUARE(xm2mg))
#define D3GDR2S2DT   0.25*R*( \
                     1.0/xm1co - 1.0/xm2co + 1.0/xm1mg - 1.0/xm2mg)
#define D3GDR2S2DP   (VR3S3)

#define D3GDR2S3S2   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S3S3   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR2S3DT   0.25*R*(1.0/xm1mg-1.0/xm2mg)
#define D3GDR2S3DP   (VR3S4)

#define D3GDR3R3R3   0.125*R*t*(1.0/SQUARE(xm1mg)+1.0/SQUARE(xm2mg) \
                               -1.0/SQUARE(xm2ni)-1.0/SQUARE(xm1ni))
#define D3GDR3R3S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR3R3S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDR3R3S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDR3R3S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg) \
                               +1.0/SQUARE(xm2ni)-1.0/SQUARE(xm1ni))
#define D3GDR3R3DT   0.25*R*( \
                     1.0/xm1mg + 1.0/xm2mg + 1.0/xm1ni + 1.0/xm2ni)
#define D3GDR3R3DP   2.0*(VR4R4)

#define D3GDR3R4S0   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR3R4S1   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR3R4S2  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR3R4S3  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR3R4DT   0.5*R* (1.0/xm2mg)
#define D3GDR3R4DP   (VR4R5)

#define D3GDR3S0S0   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S0S1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S0S2  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S0S3  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S0DT   0.25*R*(1.0/xm2mg-1.0/xm1mg)
#define D3GDR3S0DP   (VR4S1)

#define D3GDR3S1S1   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S1S2  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S1S3  -0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S1DT   0.25*R*(1.0/xm2mg-1.0/xm1mg)
#define D3GDR3S1DP   (VR4S2)

#define D3GDR3S2S2   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S2S3   0.125*R*t*(1.0/SQUARE(xm2mg)+1.0/SQUARE(xm1mg))
#define D3GDR3S2DT   0.25*R*(1.0/xm1mg-1.0/xm2mg)
#define D3GDR3S2DP   (VR4S3)


#define D3GDR3S3S3   0.125*R*t*(1.0/SQUARE(xm1mg)+1.0/SQUARE(xm2mg) \
                               -1.0/SQUARE(xm1ni)-1.0/SQUARE(xm2ni))
#define D3GDR3S3DT   0.25*R*( \
                     1.0/xm1ni - 1.0/xm2ni + 1.0/xm1mg - 1.0/xm2mg)
#define D3GDR3S3DP   (VR4S4)

#define D3GDR4R4R4   R*t*(1.0/SQUARE(xm2mg) - 1.0/SQUARE(xm2ca))
#define D3GDR4R4R0   0.50*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4R0R0   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4R4S0   0.50*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4R4S1   0.50*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4R4S2  -0.50*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4R4S3  -0.50*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4R4DT   R*(1.0/xm2mg + 1.0/xm2ca)
#define D3GDR4R4DP   2.0*(VR5R5)

#define D3GDR4S0S0   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S0S1   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S0S2  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S0S3  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S0DT   0.5*R*(1.0/xm2mg)
#define D3GDR4S0DP   (VR5S1)

#define D3GDR4S1S1   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S1S2  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S1S3  -0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S1DT   0.5*R*(1.0/xm2mg)
#define D3GDR4S1DP   (VR5S2)

#define D3GDR4S2S2   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S2S3   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S2DT  -0.5*R*(1.0/xm2mg)
#define D3GDR4S2DP   (VR5S3)

#define D3GDR4S3S3   0.25*R*t*(1.0/SQUARE(xm2mg))
#define D3GDR4S3DT  -0.5*R*(1.0/xm2mg)
#define D3GDR4S3DP   (VR5S4)

#define D3GDS0S0S0   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg) \
                               +1.0/SQUARE(xm1mn)-1.0/SQUARE(xm2mn))
#define D3GDS0S0S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDS0S0S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDS0S0S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDS0S0DT   0.25*R*(1.0/xm1mn + 1.0/xm2mn + 1.0/xm1mg +1.0/xm2mg)
#define D3GDS0S0DP   2.0*(VS1S1)

#define D3GDS0S1S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDS0S1S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDS0S1S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDS0S1DT   0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDS0S1DP   (VS1S2)

#define D3GDS0S2S2   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDS0S2S3   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDS0S2DT  -0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDS0S2DP   (VS1S3)

#define D3GDS0S3S3   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDS0S3DT  -0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDS0S3DP   (VS1S4)
#define D3GDS0DT2    0.0
#define D3GDS0DTDP   0.0
#define D3GDS0DP2    0.0

#define D3GDS1S1S1   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg) \
                               +1.0/SQUARE(xm1fe)-1.0/SQUARE(xm2fe))
#define D3GDS1S1S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDS1S1S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDS1S1DT   0.25*R*(1.0/xm1fe + 1.0/xm2fe + 1.0/xm1mg +1.0/xm2mg)
#define D3GDS1S1DP   2.0*(VS2S2)

#define D3GDS1S2S2   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDS1S2S3   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDS1S2DT  -0.25*R*   (1.0/xm1mg + 1.0/xm2mg)
#define D3GDS1S2DP   (VS2S3)

#define D3GDS1S3S3   0.125*R*t*(1.0/SQUARE(xm2mg)-1.0/SQUARE(xm1mg))
#define D3GDS1S3DT  -0.25*R*   (1.0/xm1mg + 1.0/xm2mg)
#define D3GDS1S3DP   (VS2S4)
#define D3GDS1DT2    0.0
#define D3GDS1DTDP   0.0
#define D3GDS1DP2    0.0

#define D3GDS2S2S2   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg) \
                               +1.0/SQUARE(xm2co)-1.0/SQUARE(xm1co))
#define D3GDS2S2S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDS2S2DT   0.25*R*(1.0/xm1co + 1.0/xm2co + 1.0/xm1mg +1.0/xm2mg)
#define D3GDS2S2DP   2.0*(VS3S3)

#define D3GDS2S3S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg))
#define D3GDS2S3DT   0.25*R*(1.0/xm1mg + 1.0/xm2mg)
#define D3GDS2S3DP   (VS3S4)
#define D3GDS2DT2    0.0
#define D3GDS2DTDP   0.0
#define D3GDS2DP2    0.0

#define D3GDS3S3S3   0.125*R*t*(1.0/SQUARE(xm1mg)-1.0/SQUARE(xm2mg) \
                               +1.0/SQUARE(xm2ni)-1.0/SQUARE(xm1ni))
#define D3GDS3S3DT   0.25*R*(1.0/xm1ni + 1.0/xm2ni + 1.0/xm1mg +1.0/xm2mg)
#define D3GDS3S3DP   2.0*(VS4S4)
#define D3GDS3DT2    0.0
#define D3GDS3DTDP   0.0
#define D3GDS3DP2    0.0

#define D3GDT3       0.0
#define D3GDT2DP     0.0
#define D3GDTDP2     0.0
#define D3GDP3       0.0

#define D3GDR0DT2  0.0  /* imported from spinel.c v2.0-2 */
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
/*
 *=============================================================================
 * Macros for automatic array initialization
 */
#define fillD2GDR2 \
 d2gdr2[0][0] = (D2GDR0R0);     d2gdr2[0][1] = (D2GDR0R1); \
 d2gdr2[0][2] = (D2GDR0R2);     d2gdr2[0][3] = (D2GDR0R3); \
 d2gdr2[0][4] = (D2GDR0R4); \
 d2gdr2[1][0] = d2gdr2[0][1]; d2gdr2[1][1] = (D2GDR1R1); \
 d2gdr2[1][2] = (D2GDR1R2);     d2gdr2[1][3] = (D2GDR1R3); \
 d2gdr2[1][4] = (D2GDR1R4); \
 d2gdr2[2][0] = d2gdr2[0][2]; d2gdr2[2][1] = d2gdr2[1][2]; \
 d2gdr2[2][2] = (D2GDR2R2);     d2gdr2[2][3] = (D2GDR2R3); \
 d2gdr2[2][4] = (D2GDR2R4); \
 d2gdr2[3][0] = d2gdr2[0][3]; d2gdr2[3][1] = d2gdr2[1][3]; \
 d2gdr2[3][2] = d2gdr2[2][3]; d2gdr2[3][3] = (D2GDR3R3); \
 d2gdr2[3][4] = (D2GDR3R4); \
 d2gdr2[4][0] = d2gdr2[0][4]; d2gdr2[4][1] = d2gdr2[1][4]; \
 d2gdr2[4][2] = d2gdr2[2][4]; d2gdr2[4][3] = d2gdr2[3][4]; \
 d2gdr2[4][4] = (D2GDR4R4);

#define fillD2GDRDS \
 d2gdrds[0][0] = (D2GDR0S0); d2gdrds[0][1] = (D2GDR0S1); \
 d2gdrds[0][2] = (D2GDR0S2); d2gdrds[0][3] = (D2GDR0S3); \
 d2gdrds[1][0] = (D2GDR1S0); d2gdrds[1][1] = (D2GDR1S1); \
 d2gdrds[1][2] = (D2GDR1S2); d2gdrds[1][3] = (D2GDR1S3); \
 d2gdrds[2][0] = (D2GDR2S0); d2gdrds[2][1] = (D2GDR2S1); \
 d2gdrds[2][2] = (D2GDR2S2); d2gdrds[2][3] = (D2GDR2S3); \
 d2gdrds[3][0] = (D2GDR3S0); d2gdrds[3][1] = (D2GDR3S1); \
 d2gdrds[3][2] = (D2GDR3S2); d2gdrds[3][3] = (D2GDR3S3); \
 d2gdrds[4][0] = (D2GDR4S0); d2gdrds[4][1] = (D2GDR4S1); \
 d2gdrds[4][2] = (D2GDR4S2); d2gdrds[4][3] = (D2GDR4S3);

#define fillD2GDRDT \
 d2gdrdt[0] = D2GDR0DT; d2gdrdt[1] = D2GDR1DT; d2gdrdt[2] = D2GDR2DT; \
 d2gdrdt[3] = D2GDR3DT; d2gdrdt[4] = D2GDR4DT;

#define fillD2GDRDP \
 d2gdrdp[0] = D2GDR0DP; d2gdrdp[1] = D2GDR1DP; d2gdrdp[2] = D2GDR2DP; \
 d2gdrdp[3] = D2GDR3DP; d2gdrdp[4] = D2GDR4DP;

#define fillD2GDS2 \
 d2gds2[0][0] = (D2GDS0S0);     d2gds2[0][1] = (D2GDS0S1); \
 d2gds2[0][2] = (D2GDS0S2);     d2gds2[0][3] = (D2GDS0S3); \
 d2gds2[1][0] = d2gds2[0][1]; d2gds2[1][1] = (D2GDS1S1); \
 d2gds2[1][2] = (D2GDS1S2);     d2gds2[1][3] = (D2GDS1S3); \
 d2gds2[2][0] = d2gds2[0][2]; d2gds2[2][1] = d2gds2[1][2]; \
 d2gds2[2][2] = (D2GDS2S2);     d2gds2[2][3] = (D2GDS2S3); \
 d2gds2[3][0] = d2gds2[0][3]; d2gds2[3][1] = d2gds2[1][3]; \
 d2gds2[3][2] = d2gds2[2][3]; d2gds2[3][3] = (D2GDS3S3);

#define fillD2GDSDT \
 d2gdsdt[0] = D2GDS0DT;  d2gdsdt[1] = D2GDS1DT; \
 d2gdsdt[2] = D2GDS2DT;  d2gdsdt[3] = D2GDS3DT;

#define fillD2GDSDP \
 d2gdsdp[0] = D2GDS0DP;  d2gdsdp[1] = D2GDS1DP; \
 d2gdsdp[2] = D2GDS2DP;  d2gdsdp[3] = D2GDS3DP;

#define fillD3GDR3 \
 d3gdr3[0][0][0] = D3GDR0R0R0;          d3gdr3[0][0][1] = D3GDR0R0R1; \
 d3gdr3[0][0][2] = d3gdr3[0][0][1];     d3gdr3[0][0][3] = d3gdr3[0][0][1]; \
 d3gdr3[0][0][4] = D3GDR4R0R0;          d3gdr3[0][1][0] = d3gdr3[0][0][1]; \
 d3gdr3[0][1][1] = d3gdr3[0][0][1];     d3gdr3[0][1][2] = d3gdr3[0][0][1]; \
 d3gdr3[0][1][3] = d3gdr3[0][0][1];     d3gdr3[0][1][4] = d3gdr3[0][0][4]; \
 d3gdr3[0][2][0] = d3gdr3[0][0][1];     d3gdr3[0][2][1] = d3gdr3[0][0][1]; \
 d3gdr3[0][2][2] = d3gdr3[0][0][1];     d3gdr3[0][2][3] = d3gdr3[0][0][1]; \
 d3gdr3[0][2][4] = d3gdr3[0][0][4];     d3gdr3[0][3][0] = d3gdr3[0][0][1]; \
 d3gdr3[0][3][1] = d3gdr3[0][0][1];     d3gdr3[0][3][2] = d3gdr3[0][0][1]; \
 d3gdr3[0][3][3] = d3gdr3[0][0][1];     d3gdr3[0][3][4] = d3gdr3[0][0][4]; \
 d3gdr3[0][4][0] = d3gdr3[0][0][4];     d3gdr3[0][4][1] = d3gdr3[0][0][4]; \
 d3gdr3[0][4][2] = d3gdr3[0][0][4];     d3gdr3[0][4][3] = d3gdr3[0][0][4]; \
 d3gdr3[0][4][4] = D3GDR4R4R0;          d3gdr3[1][0][0] = d3gdr3[0][0][1]; \
 d3gdr3[1][0][1] = d3gdr3[0][0][1];     d3gdr3[1][0][2] = d3gdr3[0][0][1]; \
 d3gdr3[1][0][3] = d3gdr3[0][0][1];     d3gdr3[1][0][4] = d3gdr3[0][0][4]; \
 d3gdr3[1][1][0] = d3gdr3[0][0][1];     d3gdr3[1][1][1] = D3GDR1R1R1; \
 d3gdr3[1][1][2] = d3gdr3[0][0][1];     d3gdr3[1][1][3] = d3gdr3[0][0][1]; \
 d3gdr3[1][1][4] = d3gdr3[0][0][4];     d3gdr3[1][2][0] = d3gdr3[0][0][1]; \
 d3gdr3[1][2][1] = d3gdr3[0][0][1];     d3gdr3[1][2][2] = d3gdr3[0][0][1]; \
 d3gdr3[1][2][3] = d3gdr3[0][0][1];     d3gdr3[1][2][4] = d3gdr3[0][0][4]; \
 d3gdr3[1][3][0] = d3gdr3[0][0][1];     d3gdr3[1][3][1] = d3gdr3[0][0][1]; \
 d3gdr3[1][3][2] = d3gdr3[0][0][1];     d3gdr3[1][3][3] = d3gdr3[0][0][1]; \
 d3gdr3[1][3][4] = d3gdr3[0][0][4];     d3gdr3[1][4][0] = d3gdr3[0][0][4]; \
 d3gdr3[1][4][1] = d3gdr3[0][0][4];     d3gdr3[1][4][2] = d3gdr3[0][0][4]; \
 d3gdr3[1][4][3] = d3gdr3[0][0][4];     d3gdr3[1][4][4] = d3gdr3[0][4][4]; \
 d3gdr3[2][0][0] = d3gdr3[0][0][1];     d3gdr3[2][0][1] = d3gdr3[0][0][1]; \
 d3gdr3[2][0][2] = d3gdr3[0][0][1];     d3gdr3[2][0][3] = d3gdr3[0][0][1]; \
 d3gdr3[2][0][4] = d3gdr3[0][0][4];     d3gdr3[2][1][0] = d3gdr3[0][0][1]; \
 d3gdr3[2][1][1] = d3gdr3[0][0][1];     d3gdr3[2][1][2] = d3gdr3[0][0][1]; \
 d3gdr3[2][1][3] = d3gdr3[0][0][1];     d3gdr3[2][1][4] = d3gdr3[0][0][4]; \
 d3gdr3[2][2][0] = d3gdr3[0][0][1];     d3gdr3[2][2][1] = d3gdr3[0][0][1]; \
 d3gdr3[2][2][2] = D3GDR2R2R2;          d3gdr3[2][2][3] = d3gdr3[0][0][1]; \
 d3gdr3[2][2][4] = d3gdr3[0][0][4];     d3gdr3[2][3][0] = d3gdr3[0][0][1]; \
 d3gdr3[2][3][1] = d3gdr3[0][0][1];     d3gdr3[2][3][2] = d3gdr3[0][0][1]; \
 d3gdr3[2][3][3] = d3gdr3[0][0][1];     d3gdr3[2][3][4] = d3gdr3[0][0][4]; \
 d3gdr3[2][4][0] = d3gdr3[0][0][4];     d3gdr3[2][4][1] = d3gdr3[0][0][4]; \
 d3gdr3[2][4][2] = d3gdr3[0][0][4];     d3gdr3[2][4][3] = d3gdr3[0][0][4]; \
 d3gdr3[2][4][4] = d3gdr3[0][4][4];     d3gdr3[3][0][0] = d3gdr3[0][0][1]; \
 d3gdr3[3][0][1] = d3gdr3[0][0][1];     d3gdr3[3][0][2] = d3gdr3[0][0][1]; \
 d3gdr3[3][0][3] = d3gdr3[0][0][1];     d3gdr3[3][0][4] = d3gdr3[0][0][4]; \
 d3gdr3[3][1][0] = d3gdr3[0][0][1];     d3gdr3[3][1][1] = d3gdr3[0][0][1]; \
 d3gdr3[3][1][2] = d3gdr3[0][0][1];     d3gdr3[3][1][3] = d3gdr3[0][0][1]; \
 d3gdr3[3][1][4] = d3gdr3[0][0][4];     d3gdr3[3][2][0] = d3gdr3[0][0][1]; \
 d3gdr3[3][2][1] = d3gdr3[0][0][1];     d3gdr3[3][2][2] = d3gdr3[0][0][1]; \
 d3gdr3[3][2][3] = d3gdr3[0][0][1];     d3gdr3[3][2][4] = d3gdr3[0][0][4]; \
 d3gdr3[3][3][0] = d3gdr3[0][0][1];     d3gdr3[3][3][1] = d3gdr3[0][0][1]; \
 d3gdr3[3][3][2] = d3gdr3[0][0][1];     d3gdr3[3][3][3] = D3GDR3R3R3; \
 d3gdr3[3][3][4] = d3gdr3[0][0][4];     d3gdr3[3][4][0] = d3gdr3[0][0][4]; \
 d3gdr3[3][4][1] = d3gdr3[0][0][4];     d3gdr3[3][4][2] = d3gdr3[0][0][4]; \
 d3gdr3[3][4][3] = d3gdr3[0][0][4];     d3gdr3[3][4][4] = d3gdr3[0][4][4]; \
 d3gdr3[4][0][0] = d3gdr3[0][0][4];     d3gdr3[4][0][1] = d3gdr3[0][0][4]; \
 d3gdr3[4][0][2] = d3gdr3[0][0][4];     d3gdr3[4][0][3] = d3gdr3[0][0][4]; \
 d3gdr3[4][0][4] = d3gdr3[0][4][4];     d3gdr3[4][1][0] = d3gdr3[0][0][4]; \
 d3gdr3[4][1][1] = d3gdr3[0][0][4];     d3gdr3[4][1][2] = d3gdr3[0][0][4]; \
 d3gdr3[4][1][3] = d3gdr3[0][0][4];     d3gdr3[4][1][4] = d3gdr3[0][4][4]; \
 d3gdr3[4][2][0] = d3gdr3[0][0][4];     d3gdr3[4][2][1] = d3gdr3[0][0][4]; \
 d3gdr3[4][2][2] = d3gdr3[0][0][4];     d3gdr3[4][2][3] = d3gdr3[0][0][4]; \
 d3gdr3[4][2][4] = d3gdr3[0][4][4];     d3gdr3[4][3][0] = d3gdr3[0][0][4]; \
 d3gdr3[4][3][1] = d3gdr3[0][0][4];     d3gdr3[4][3][2] = d3gdr3[0][0][4]; \
 d3gdr3[4][3][3] = d3gdr3[0][0][4];     d3gdr3[4][3][4] = d3gdr3[0][4][4]; \
 d3gdr3[4][4][0] = d3gdr3[0][4][4];     d3gdr3[4][4][1] = d3gdr3[0][4][4]; \
 d3gdr3[4][4][2] = d3gdr3[0][4][4];     d3gdr3[4][4][3] = d3gdr3[0][4][4]; \
 d3gdr3[4][4][4] = D3GDR4R4R4;

#define fillD3GDR2DS \
 d3gdr2ds[0][0][0] = D3GDR0R0S0;        d3gdr2ds[0][0][1] = D3GDR0R0S1; \
 d3gdr2ds[0][0][2] = D3GDR0R0S2;        d3gdr2ds[0][0][3] = D3GDR0R0S3; \
 d3gdr2ds[0][1][0] = D3GDR0R1S0;        d3gdr2ds[0][1][1] = D3GDR0R1S1; \
 d3gdr2ds[0][1][2] = D3GDR0R1S2;        d3gdr2ds[0][1][3] = D3GDR0R1S3; \
 d3gdr2ds[0][2][0] = D3GDR0R2S0;        d3gdr2ds[0][2][1] = D3GDR0R2S1; \
 d3gdr2ds[0][2][2] = D3GDR0R2S2;        d3gdr2ds[0][2][3] = D3GDR0R2S3; \
 d3gdr2ds[0][3][0] = D3GDR0R3S0;        d3gdr2ds[0][3][1] = D3GDR0R3S1; \
 d3gdr2ds[0][3][2] = D3GDR0R3S2;        d3gdr2ds[0][3][3] = D3GDR0R3S3; \
 d3gdr2ds[0][4][0] = D3GDR0R4S0;        d3gdr2ds[0][4][1] = D3GDR0R4S1; \
 d3gdr2ds[0][4][2] = D3GDR0R4S2;        d3gdr2ds[0][4][3] = D3GDR0R4S3; \
 d3gdr2ds[1][0][0] = d3gdr2ds[0][1][0]; d3gdr2ds[1][0][1] = d3gdr2ds[0][1][1]; \
 d3gdr2ds[1][0][2] = d3gdr2ds[0][1][2]; d3gdr2ds[1][0][3] = d3gdr2ds[0][1][3]; \
 d3gdr2ds[1][1][0] = D3GDR1R1S0;        d3gdr2ds[1][1][1] = D3GDR1R1S1; \
 d3gdr2ds[1][1][2] = D3GDR1R1S2;        d3gdr2ds[1][1][3] = D3GDR1R1S3; \
 d3gdr2ds[1][2][0] = D3GDR1R2S0;        d3gdr2ds[1][2][1] = D3GDR1R2S1; \
 d3gdr2ds[1][2][2] = D3GDR1R2S2;        d3gdr2ds[1][2][3] = D3GDR1R2S3; \
 d3gdr2ds[1][3][0] = D3GDR1R3S0;        d3gdr2ds[1][3][1] = D3GDR1R3S1; \
 d3gdr2ds[1][3][2] = D3GDR1R3S2;        d3gdr2ds[1][3][3] = D3GDR1R3S3; \
 d3gdr2ds[1][4][0] = D3GDR1R4S0;        d3gdr2ds[1][4][1] = D3GDR1R4S1; \
 d3gdr2ds[1][4][2] = D3GDR1R4S2;        d3gdr2ds[1][4][3] = D3GDR1R4S3; \
 d3gdr2ds[2][0][0] = d3gdr2ds[0][2][0]; d3gdr2ds[2][0][1] = d3gdr2ds[0][2][1]; \
 d3gdr2ds[2][0][2] = d3gdr2ds[0][2][2]; d3gdr2ds[2][0][3] = d3gdr2ds[0][2][3]; \
 d3gdr2ds[2][1][0] = d3gdr2ds[1][2][0]; d3gdr2ds[2][1][1] = d3gdr2ds[1][2][1]; \
 d3gdr2ds[2][1][2] = d3gdr2ds[1][2][2]; d3gdr2ds[2][1][3] = d3gdr2ds[1][2][3]; \
 d3gdr2ds[2][2][0] = D3GDR2R2S0;        d3gdr2ds[2][2][1] = D3GDR2R2S1; \
 d3gdr2ds[2][2][2] = D3GDR2R2S2;        d3gdr2ds[2][2][3] = D3GDR2R2S3; \
 d3gdr2ds[2][3][0] = D3GDR2R3S0;        d3gdr2ds[2][3][1] = D3GDR2R3S1; \
 d3gdr2ds[2][3][2] = D3GDR2R3S2;        d3gdr2ds[2][3][3] = D3GDR2R3S3; \
 d3gdr2ds[2][4][0] = D3GDR2R4S0;        d3gdr2ds[2][4][1] = D3GDR2R4S1; \
 d3gdr2ds[2][4][2] = D3GDR2R4S2;        d3gdr2ds[2][4][3] = D3GDR2R4S3; \
 d3gdr2ds[3][0][0] = d3gdr2ds[0][3][0]; d3gdr2ds[3][0][1] = d3gdr2ds[0][3][1]; \
 d3gdr2ds[3][0][2] = d3gdr2ds[0][3][2]; d3gdr2ds[3][0][3] = d3gdr2ds[0][3][3]; \
 d3gdr2ds[3][1][0] = d3gdr2ds[1][3][0]; d3gdr2ds[3][1][1] = d3gdr2ds[1][3][1]; \
 d3gdr2ds[3][1][2] = d3gdr2ds[1][3][2]; d3gdr2ds[3][1][3] = d3gdr2ds[1][3][3]; \
 d3gdr2ds[3][2][0] = d3gdr2ds[2][3][0]; d3gdr2ds[3][2][1] = d3gdr2ds[2][3][1]; \
 d3gdr2ds[3][2][2] = d3gdr2ds[2][3][2]; d3gdr2ds[3][2][3] = d3gdr2ds[2][3][3]; \
 d3gdr2ds[3][3][0] = D3GDR3R3S0;        d3gdr2ds[3][3][1] = D3GDR3R3S1; \
 d3gdr2ds[3][3][2] = D3GDR3R3S2;        d3gdr2ds[3][3][3] = D3GDR3R3S3; \
 d3gdr2ds[3][4][0] = D3GDR3R4S0;        d3gdr2ds[3][4][1] = D3GDR3R4S1; \
 d3gdr2ds[3][4][2] = D3GDR3R4S2;        d3gdr2ds[3][4][3] = D3GDR3R4S3; \
 d3gdr2ds[4][0][0] = d3gdr2ds[0][4][0]; d3gdr2ds[4][0][1] = d3gdr2ds[0][4][1]; \
 d3gdr2ds[4][0][2] = d3gdr2ds[0][4][2]; d3gdr2ds[4][0][3] = d3gdr2ds[0][4][3]; \
 d3gdr2ds[4][1][0] = d3gdr2ds[1][4][0]; d3gdr2ds[4][1][1] = d3gdr2ds[1][4][1]; \
 d3gdr2ds[4][1][2] = d3gdr2ds[1][4][2]; d3gdr2ds[4][1][3] = d3gdr2ds[1][4][3]; \
 d3gdr2ds[4][2][0] = d3gdr2ds[2][4][0]; d3gdr2ds[4][2][1] = d3gdr2ds[2][4][1]; \
 d3gdr2ds[4][2][2] = d3gdr2ds[2][4][2]; d3gdr2ds[4][2][3] = d3gdr2ds[2][4][3]; \
 d3gdr2ds[4][3][0] = d3gdr2ds[3][4][0]; d3gdr2ds[4][3][1] = d3gdr2ds[3][4][1]; \
 d3gdr2ds[4][3][2] = d3gdr2ds[3][4][2]; d3gdr2ds[4][3][3] = d3gdr2ds[3][4][3]; \
 d3gdr2ds[4][4][0] = D3GDR4R4S0;        d3gdr2ds[4][4][1] = D3GDR4R4S1; \
 d3gdr2ds[4][4][2] = D3GDR4R4S2;        d3gdr2ds[4][4][3] = D3GDR4R4S3;

#define fillD3GDR2DT \
 d3gdr2dt[0][0] = D3GDR0R0DT;     d3gdr2dt[0][1] = D3GDR0R1DT; \
 d3gdr2dt[0][2] = D3GDR0R2DT;     d3gdr2dt[0][3] = D3GDR0R3DT; \
 d3gdr2dt[0][4] = D3GDR0R4DT; \
 d3gdr2dt[1][0] = d3gdr2dt[0][1]; d3gdr2dt[1][1] = D3GDR1R1DT; \
 d3gdr2dt[1][2] = D3GDR1R2DT;     d3gdr2dt[1][3] = D3GDR1R3DT; \
 d3gdr2dt[1][4] = D3GDR1R4DT; \
 d3gdr2dt[2][0] = d3gdr2dt[0][2]; d3gdr2dt[2][1] = d3gdr2dt[1][2]; \
 d3gdr2dt[2][2] = D3GDR2R2DT;     d3gdr2dt[2][3] = D3GDR2R3DT; \
 d3gdr2dt[2][4] = D3GDR2R4DT; \
 d3gdr2dt[3][0] = d3gdr2dt[0][3]; d3gdr2dt[3][1] = d3gdr2dt[1][3]; \
 d3gdr2dt[3][2] = d3gdr2dt[2][3]; d3gdr2dt[3][3] = D3GDR3R3DT;     \
 d3gdr2dt[3][4] = D3GDR3R4DT; \
 d3gdr2dt[4][0] = d3gdr2dt[0][4]; d3gdr2dt[4][1] = d3gdr2dt[1][4]; \
 d3gdr2dt[4][2] = d3gdr2dt[2][4]; d3gdr2dt[4][3] = d3gdr2dt[3][4]; \
 d3gdr2dt[4][4] = D3GDR4R4DT;

#define fillD3GDR2DP \
 d3gdr2dp[0][0] = (D3GDR0R0DP);     d3gdr2dp[0][1] = (D3GDR0R1DP); \
 d3gdr2dp[0][2] = (D3GDR0R2DP);     d3gdr2dp[0][3] = (D3GDR0R3DP); \
 d3gdr2dp[0][4] = (D3GDR0R4DP); \
 d3gdr2dp[1][0] = d3gdr2dp[0][1]; d3gdr2dp[1][1] = (D3GDR1R1DP); \
 d3gdr2dp[1][2] = (D3GDR1R2DP);     d3gdr2dp[1][3] = (D3GDR1R3DP); \
 d3gdr2dp[1][4] = (D3GDR1R4DP); \
 d3gdr2dp[2][0] = d3gdr2dp[0][2]; d3gdr2dp[2][1] = d3gdr2dp[1][2]; \
 d3gdr2dp[2][2] = (D3GDR2R2DP);     d3gdr2dp[2][3] = (D3GDR2R3DP); \
 d3gdr2dp[2][4] = (D3GDR2R4DP); \
 d3gdr2dp[3][0] = d3gdr2dp[0][3]; d3gdr2dp[3][1] = d3gdr2dp[1][3]; \
 d3gdr2dp[3][2] = d3gdr2dp[2][3]; d3gdr2dp[3][3] = (D3GDR3R3DP); \
 d3gdr2dp[3][4] = (D3GDR3R4DP); \
 d3gdr2dp[4][0] = d3gdr2dp[0][4]; d3gdr2dp[4][1] = d3gdr2dp[1][4]; \
 d3gdr2dp[4][2] = d3gdr2dp[2][4]; d3gdr2dp[4][3] = d3gdr2dp[3][4]; \
 d3gdr2dp[4][4] = (D3GDR4R4DP);

#define fillD3GDRDS2 \
 d3gdrds2[0][0][0] = D3GDR0S0S0;        d3gdrds2[0][0][1] = D3GDR0S0S1; \
 d3gdrds2[0][0][2] = D3GDR0S0S2;        d3gdrds2[0][0][3] = D3GDR0S0S3; \
 d3gdrds2[0][1][0] = d3gdrds2[0][0][1]; d3gdrds2[0][1][1] = D3GDR0S1S1; \
 d3gdrds2[0][1][2] = D3GDR0S1S2;        d3gdrds2[0][1][3] = D3GDR0S1S3; \
 d3gdrds2[0][2][0] = d3gdrds2[0][0][2]; d3gdrds2[0][2][1] = d3gdrds2[0][1][2]; \
 d3gdrds2[0][2][2] = D3GDR0S2S2;        d3gdrds2[0][2][3] = D3GDR0S2S3; \
 d3gdrds2[0][3][0] = d3gdrds2[0][0][3]; d3gdrds2[0][3][1] = d3gdrds2[0][1][3]; \
 d3gdrds2[0][3][2] = d3gdrds2[0][2][3]; d3gdrds2[0][3][3] = D3GDR0S3S3; \
 d3gdrds2[1][0][0] = D3GDR1S0S0;        d3gdrds2[1][0][1] = D3GDR1S0S1; \
 d3gdrds2[1][0][2] = D3GDR1S0S2;        d3gdrds2[1][0][3] = D3GDR1S0S3; \
 d3gdrds2[1][1][0] = d3gdrds2[1][0][1]; d3gdrds2[1][1][1] = D3GDR1S1S1; \
 d3gdrds2[1][1][2] = D3GDR1S1S2;        d3gdrds2[1][1][3] = D3GDR1S1S3; \
 d3gdrds2[1][2][0] = d3gdrds2[1][0][2]; d3gdrds2[1][2][1] = d3gdrds2[1][1][2]; \
 d3gdrds2[1][2][2] = D3GDR1S2S2;        d3gdrds2[1][2][3] = D3GDR1S2S3; \
 d3gdrds2[1][3][0] = d3gdrds2[1][0][3]; d3gdrds2[1][3][1] = d3gdrds2[1][1][3]; \
 d3gdrds2[1][3][2] = d3gdrds2[1][2][3]; d3gdrds2[1][3][3] = D3GDR1S3S3; \
 d3gdrds2[2][0][0] = D3GDR2S0S0;        d3gdrds2[2][0][1] = D3GDR2S0S1; \
 d3gdrds2[2][0][2] = D3GDR2S0S2;        d3gdrds2[2][0][3] = D3GDR2S0S3; \
 d3gdrds2[2][1][0] = d3gdrds2[2][0][1]; d3gdrds2[2][1][1] = D3GDR2S1S1; \
 d3gdrds2[2][1][2] = D3GDR2S1S2;        d3gdrds2[2][1][3] = D3GDR2S1S3; \
 d3gdrds2[2][2][0] = d3gdrds2[2][0][2]; d3gdrds2[2][2][1] = d3gdrds2[2][1][2]; \
 d3gdrds2[2][2][2] = D3GDR2S2S2;        d3gdrds2[2][2][3] = D3GDR2S2S3; \
 d3gdrds2[2][3][0] = d3gdrds2[2][0][3]; d3gdrds2[2][3][1] = d3gdrds2[2][1][3]; \
 d3gdrds2[2][3][2] = d3gdrds2[2][2][3]; d3gdrds2[2][3][3] = D3GDR2S3S3; \
 d3gdrds2[3][0][0] = D3GDR3S0S0;        d3gdrds2[3][0][1] = D3GDR3S0S1; \
 d3gdrds2[3][0][2] = D3GDR3S0S2;        d3gdrds2[3][0][3] = D3GDR3S0S3; \
 d3gdrds2[3][1][0] = d3gdrds2[3][0][1]; d3gdrds2[3][1][1] = D3GDR3S1S1; \
 d3gdrds2[3][1][2] = D3GDR3S1S2;        d3gdrds2[3][1][3] = D3GDR3S1S3; \
 d3gdrds2[3][2][0] = d3gdrds2[3][0][2]; d3gdrds2[3][2][1] = d3gdrds2[3][1][2]; \
 d3gdrds2[3][2][2] = D3GDR3S2S2;        d3gdrds2[3][2][3] = D3GDR3S2S3; \
 d3gdrds2[3][3][0] = d3gdrds2[3][0][3]; d3gdrds2[3][3][1] = d3gdrds2[3][1][3]; \
 d3gdrds2[3][3][2] = d3gdrds2[3][2][3]; d3gdrds2[3][3][3] = D3GDR3S3S3; \
 d3gdrds2[4][0][0] = D3GDR4S0S0;        d3gdrds2[4][0][1] = D3GDR4S0S1; \
 d3gdrds2[4][0][2] = D3GDR4S0S2;        d3gdrds2[4][0][3] = D3GDR4S0S3; \
 d3gdrds2[4][1][0] = d3gdrds2[4][0][1]; d3gdrds2[4][1][1] = D3GDR4S1S1; \
 d3gdrds2[4][1][2] = D3GDR4S1S2;        d3gdrds2[4][1][3] = D3GDR4S1S3; \
 d3gdrds2[4][2][0] = d3gdrds2[4][0][2]; d3gdrds2[4][2][1] = d3gdrds2[4][1][2]; \
 d3gdrds2[4][2][2] = D3GDR4S2S2;        d3gdrds2[4][2][3] = D3GDR4S2S3; \
 d3gdrds2[4][3][0] = d3gdrds2[4][0][3]; d3gdrds2[4][3][1] = d3gdrds2[4][1][3]; \
 d3gdrds2[4][3][2] = d3gdrds2[4][2][3]; d3gdrds2[4][3][3] = D3GDR4S3S3; \

#define fillD3GDRDSDT \
 d3gdrdsdt[0][0] = D3GDR0S0DT;  d3gdrdsdt[0][1] = D3GDR0S1DT; \
 d3gdrdsdt[0][2] = D3GDR0S2DT;  d3gdrdsdt[0][3] = D3GDR0S3DT; \
 d3gdrdsdt[1][0] = D3GDR1S0DT;  d3gdrdsdt[1][1] = D3GDR1S1DT; \
 d3gdrdsdt[1][2] = D3GDR1S2DT;  d3gdrdsdt[1][3] = D3GDR1S3DT; \
 d3gdrdsdt[2][0] = D3GDR2S0DT;  d3gdrdsdt[2][1] = D3GDR2S1DT; \
 d3gdrdsdt[2][2] = D3GDR2S2DT;  d3gdrdsdt[2][3] = D3GDR2S3DT; \
 d3gdrdsdt[3][0] = D3GDR3S0DT;  d3gdrdsdt[3][1] = D3GDR3S1DT; \
 d3gdrdsdt[3][2] = D3GDR3S2DT;  d3gdrdsdt[3][3] = D3GDR3S3DT; \
 d3gdrdsdt[4][0] = D3GDR4S0DT;  d3gdrdsdt[4][1] = D3GDR4S1DT; \
 d3gdrdsdt[4][2] = D3GDR4S2DT;  d3gdrdsdt[4][3] = D3GDR4S3DT;

#define fillD3GDRDSDP \
 d3gdrdsdp[0][0] = D3GDR0S0DP;  d3gdrdsdp[0][1] = D3GDR0S1DP; \
 d3gdrdsdp[0][2] = D3GDR0S2DP;  d3gdrdsdp[0][3] = D3GDR0S3DP; \
 d3gdrdsdp[1][0] = D3GDR1S0DP;  d3gdrdsdp[1][1] = D3GDR1S1DP; \
 d3gdrdsdp[1][2] = D3GDR1S2DP;  d3gdrdsdp[1][3] = D3GDR1S3DP; \
 d3gdrdsdp[2][0] = D3GDR2S0DP;  d3gdrdsdp[2][1] = D3GDR2S1DP; \
 d3gdrdsdp[2][2] = D3GDR2S2DP;  d3gdrdsdp[2][3] = D3GDR2S3DP; \
 d3gdrdsdp[3][0] = D3GDR3S0DP;  d3gdrdsdp[3][1] = D3GDR3S1DP; \
 d3gdrdsdp[3][2] = D3GDR3S2DP;  d3gdrdsdp[3][3] = D3GDR3S3DP; \
 d3gdrdsdp[4][0] = D3GDR4S0DP;  d3gdrdsdp[4][1] = D3GDR4S1DP; \
 d3gdrdsdp[4][2] = D3GDR4S2DP;  d3gdrdsdp[4][3] = D3GDR4S3DP;

#define fillD3GDS3 \
 d3gds3[0][0][0] = D3GDS0S0S0;      d3gds3[0][0][1] = D3GDS0S0S1; \
 d3gds3[0][0][2] = D3GDS0S0S2;      d3gds3[0][0][3] = D3GDS0S0S3; \
 d3gds3[0][1][0] = d3gds3[0][0][1]; d3gds3[0][1][1] = D3GDS0S1S1; \
 d3gds3[0][1][2] = D3GDS0S1S2;      d3gds3[0][1][3] = D3GDS0S1S3; \
 d3gds3[0][2][0] = d3gds3[0][0][2]; d3gds3[0][2][1] = d3gds3[0][1][2]; \
 d3gds3[0][2][2] = D3GDS0S2S2;      d3gds3[0][2][3] = D3GDS0S2S3; \
 d3gds3[0][3][0] = d3gds3[0][0][3]; d3gds3[0][3][1] = d3gds3[0][1][3]; \
 d3gds3[0][3][2] = d3gds3[0][2][3]; d3gds3[0][3][3] = D3GDS0S3S3; \
 d3gds3[1][0][0] = d3gds3[0][0][1]; d3gds3[1][0][1] = d3gds3[0][1][1]; \
 d3gds3[1][0][2] = d3gds3[0][1][2]; d3gds3[1][0][3] = d3gds3[0][1][3]; \
 d3gds3[1][1][0] = d3gds3[0][1][1]; d3gds3[1][1][1] = D3GDS1S1S1; \
 d3gds3[1][1][2] = D3GDS1S1S2;      d3gds3[1][1][3] = D3GDS1S1S3; \
 d3gds3[1][2][0] = d3gds3[0][1][2]; d3gds3[1][2][1] = d3gds3[1][1][2]; \
 d3gds3[1][2][2] = D3GDS1S2S2;      d3gds3[1][2][3] = D3GDS1S2S3; \
 d3gds3[1][3][0] = d3gds3[0][1][3]; d3gds3[1][3][1] = d3gds3[1][1][3]; \
 d3gds3[1][3][2] = d3gds3[1][2][3]; d3gds3[1][3][3] = D3GDS1S3S3; \
 d3gds3[2][0][0] = d3gds3[0][0][2]; d3gds3[2][0][1] = d3gds3[0][1][2]; \
 d3gds3[2][0][2] = d3gds3[0][2][2]; d3gds3[2][0][3] = d3gds3[0][2][3]; \
 d3gds3[2][1][0] = d3gds3[0][1][2]; d3gds3[2][1][1] = d3gds3[1][1][2]; \
 d3gds3[2][1][2] = d3gds3[1][2][2]; d3gds3[2][1][3] = d3gds3[1][2][3]; \
 d3gds3[2][2][0] = d3gds3[0][2][2]; d3gds3[2][2][1] = d3gds3[1][2][2]; \
 d3gds3[2][2][2] = D3GDS2S2S2;      d3gds3[2][2][3] = D3GDS2S2S3; \
 d3gds3[2][3][0] = d3gds3[0][2][3]; d3gds3[2][3][1] = d3gds3[1][2][3]; \
 d3gds3[2][3][2] = d3gds3[2][2][3]; d3gds3[2][3][3] = D3GDS2S3S3; \
 d3gds3[3][0][0] = d3gds3[0][0][3]; d3gds3[3][0][1] = d3gds3[0][1][3]; \
 d3gds3[3][0][2] = d3gds3[0][2][3]; d3gds3[3][0][3] = d3gds3[0][3][3]; \
 d3gds3[3][1][0] = d3gds3[0][1][3]; d3gds3[3][1][1] = d3gds3[1][1][3]; \
 d3gds3[3][1][2] = d3gds3[1][2][3]; d3gds3[3][1][3] = d3gds3[1][3][3]; \
 d3gds3[3][2][0] = d3gds3[0][2][3]; d3gds3[3][2][1] = d3gds3[1][2][3]; \
 d3gds3[3][2][2] = d3gds3[2][2][3]; d3gds3[3][2][3] = d3gds3[2][3][3]; \
 d3gds3[3][3][0] = d3gds3[0][3][3]; d3gds3[3][3][1] = d3gds3[1][3][3]; \
 d3gds3[3][3][2] = d3gds3[2][3][3]; d3gds3[3][3][3] = D3GDS3S3S3;

#define fillD3GDS2DT \
 d3gds2dt[0][0] = D3GDS0S0DT;     d3gds2dt[0][1] = D3GDS0S1DT; \
 d3gds2dt[0][2] = D3GDS0S2DT;     d3gds2dt[0][3] = D3GDS0S3DT; \
 d3gds2dt[1][0] = d3gds2dt[0][1]; d3gds2dt[1][1] = D3GDS1S1DT; \
 d3gds2dt[1][2] = D3GDS1S2DT;     d3gds2dt[1][3] = D3GDS1S3DT; \
 d3gds2dt[2][0] = d3gds2dt[0][2]; d3gds2dt[2][1] = d3gds2dt[1][2]; \
 d3gds2dt[2][2] = D3GDS2S2DT;     d3gds2dt[2][3] = D3GDS2S3DT; \
 d3gds2dt[3][0] = d3gds2dt[0][3]; d3gds2dt[3][1] = d3gds2dt[1][3]; \
 d3gds2dt[3][2] = d3gds2dt[2][3]; d3gds2dt[3][3] = D3GDS3S3DT;

#define fillD3GDS2DP \
 d3gds2dp[0][0] = D3GDS0S0DP;     d3gds2dp[0][1] = D3GDS0S1DP; \
 d3gds2dp[0][2] = D3GDS0S2DP;     d3gds2dp[0][3] = D3GDS0S3DP; \
 d3gds2dp[1][0] = d3gds2dp[0][1]; d3gds2dp[1][1] = D3GDS1S1DP; \
 d3gds2dp[1][2] = D3GDS1S2DP;     d3gds2dp[1][3] = D3GDS1S3DP; \
 d3gds2dp[2][0] = d3gds2dp[0][2]; d3gds2dp[2][1] = d3gds2dp[1][2]; \
 d3gds2dp[2][2] = D3GDS2S2DP;     d3gds2dp[2][3] = D3GDS2S3DP; \
 d3gds2dp[3][0] = d3gds2dp[0][3]; d3gds2dp[3][1] = d3gds2dp[1][3]; \
 d3gds2dp[3][2] = d3gds2dp[2][3]; d3gds2dp[3][3] = D3GDS3S3DP;

#define fillD3GDSDT2 \
 d3gdsdt2[0] = D3GDS0DT2; d3gdsdt2[1] = D3GDS1DT2; \
 d3gdsdt2[2] = D3GDS2DT2; d3gdsdt2[3] = D3GDS3DT2;

#define fillD3GDSDTDP \
 d3gdsdtdp[0] = D3GDS0DTDP; d3gdsdtdp[1] = D3GDS1DTDP; \
 d3gdsdtdp[2] = D3GDS2DTDP; d3gdsdtdp[3] = D3GDS3DTDP;

#define fillD3GDSDP2 \
 d3gdsdp2[0] = D3GDS0DP2; d3gdsdp2[1] = D3GDS1DP2; \
 d3gdsdp2[2] = D3GDS2DP2; d3gdsdp2[3] = D3GDS3DP2;

#define fillD3GDRDT2 \
 d3gdrdt2[0] = D3GDR0DT2; d3gdrdt2[1] = D3GDR1DT2; \
 d3gdrdt2[2] = D3GDR2DT2; d3gdrdt2[3] = D3GDR3DT2; \
 d3gdrdt2[4] = D3GDR4DT2;

#define fillD3GDRDTDP \
 d3gdrdtdp[0] = D3GDR0DTDP; d3gdrdtdp[1] = D3GDR1DTDP; \
 d3gdrdtdp[2] = D3GDR2DTDP; d3gdrdtdp[3] = D3GDR3DTDP; \
 d3gdrdtdp[4] = D3GDR4DTDP;

#define fillD3GDRDP2 \
 d3gdrdp2[0] = D3GDR0DP2; d3gdrdp2[1] = D3GDR1DP2; \
 d3gdrdp2[2] = D3GDR2DP2; d3gdrdp2[3] = D3GDR3DP2; \
 d3gdrdp2[4] = D3GDR4DP2;
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
    gsl_matrix       *ptToD2gds2 = getPtToD2gds2();
    gsl_matrix       *ptToV = getPtToV();
    gsl_vector       *ptToS = getPtToS();
    gsl_permutation  *indexD2gds2 = getIndexD2gds2();
    int i, j, iter = 0, signum;

    GET_SITE_FRACTIONS

    /* look-up or compute the current ordering state */
    if ( (t != tOld)       || (p != pOld) ||
       (r[0] != rOld[0]) || (r[1] != rOld[1]) || (r[2] != rOld[2]) ||
       (r[3] != rOld[3]) || (r[4] != rOld[4]) ) {
        double dgds[NS], sNew[NS], sMax[NS], sMin[NS];
        gsl_vector_view vvToDgds = gsl_vector_view_array(dgds, (size_t) NS);

        for (i=0; i<NS; i++) sOld[i] = 2.0;

        /* calculate initial guesses for ordering variables  */

        sMax[0]  = MIN(1.0-r[0],1.0+r[0]);
        sMax[0]  = MIN(sMax[0],1.0-r[0]-2.0*r[4]);
        sMin[0]  = MAX(r[0]-1.0,-r[0]-1.0);

        sMax[1]  = MIN(1.0-r[1],1.0+r[1]);
        sMax[1]  = MIN(sMax[1],1.0-r[1]-2.0*r[4]);
        sMin[1]  = MAX(r[1]-1.0,-r[1]-1.0);

        sMax[2]  = MIN(1.0-r[2],1.0+r[2]);
        sMin[2]  = MAX(r[2]-1.0,-r[2]-1.0);
        sMin[2]  = MAX(sMin[2],r[2]-1.0-2.0*r[4]);

        sMax[3]  = MIN(1.0-r[3],1.0+r[3]);
        sMin[3]  = MAX(r[3]-1.0,-r[3]-1.0);
        sMin[3]  = MAX(sMin[3],r[3]-1.0-2.0*r[4]);

        sNew[0] = (r[0]>0.0) ? sMax[0]/3.0  :0.35*sMin[0]+0.65*sMax[0];
        sNew[1] = (r[1]>0.0) ? sMax[1]/15.0 :0.45*sMin[1]+0.55*sMax[1];
        sNew[2] = (r[2]>0.0) ? sMax[2]/3.0  :0.25*sMin[2]+0.75*sMax[2];
        sNew[3] = (r[3]>0.0) ? sMax[3]/2.0  :0.25*sMin[3]+0.75*sMax[3];

        /* END OF INITIAL GUESSES    */

        while ( ((ABS(sNew[0]-sOld[0]) > 100.0*DBL_EPSILON) ||
             (ABS(sNew[1]-sOld[1]) > 100.0*DBL_EPSILON) ||
             (ABS(sNew[2]-sOld[2]) > 100.0*DBL_EPSILON) ||
             (ABS(sNew[3]-sOld[3]) > 100.0*DBL_EPSILON) ) && (iter < MAX_ITER)) {
            double s[NS], sCorr[NS], lambda;
            gsl_vector_view vvToSCorr = gsl_vector_view_array(sCorr, (size_t) NS);

            for (i=0; i<NS; i++) s[i] = sNew[i];

            if (r[0] == -1.0) s[0]=0.0;
            if (r[1] == -1.0) s[1]=0.0;
            if (r[2] == -1.0) s[2]=0.0;
            if (r[3] == -1.0) s[3]=0.0;

            xm1mn= (r[0]-s[0]+1.0)/2.0;
            xm2mn= (r[0]+s[0]+1.0)/2.0;
            xm1fe= (r[1]-s[1]+1.0)/2.0;
            xm2fe= (r[1]+s[1]+1.0)/2.0;
            xm1co= (r[2]+s[2]+1.0)/2.0;
            xm2co= (r[2]-s[2]+1.0)/2.0;
            xm1ni= (r[3]+s[3]+1.0)/2.0;
            xm2ni= (r[3]-s[3]+1.0)/2.0;
            xm2ca=  r[4];
            xm1mg = (1.0 - xm1mn - xm1fe - xm1co - xm1ni);
            xm2mg = (1.0 - xm2mn - xm2fe - xm2co - xm2ni - xm2ca);

            if (xm1mn   <= DBL_EPSILON) xm1mn   = DBL_EPSILON;
            if (xm1fe   <= DBL_EPSILON) xm1fe   = DBL_EPSILON;
            if (xm1co   <= DBL_EPSILON) xm1co   = DBL_EPSILON;
            if (xm1ni   <= DBL_EPSILON) xm1ni   = DBL_EPSILON;
            if (xm1mg   <= DBL_EPSILON) xm1mg   = DBL_EPSILON;
            if (xm2mn   <= DBL_EPSILON) xm2mn   = DBL_EPSILON;
            if (xm2fe   <= DBL_EPSILON) xm2fe   = DBL_EPSILON;
            if (xm2co   <= DBL_EPSILON) xm2co   = DBL_EPSILON;
            if (xm2ni   <= DBL_EPSILON) xm2ni   = DBL_EPSILON;
            if (xm2mg   <= DBL_EPSILON) xm2mg   = DBL_EPSILON;
            if (xm2ca   <= DBL_EPSILON) xm2ca   = DBL_EPSILON;

            dgds[0] = (fabs(r[0]+1.0)>10.0*DBL_EPSILON) ? DGDS0 : 0.0;
            dgds[1] = (fabs(r[1]+1.0)>10.0*DBL_EPSILON) ? DGDS1 : 0.0;
            dgds[2] = (fabs(r[2]+1.0)>10.0*DBL_EPSILON) ? DGDS2 : 0.0;
            dgds[3] = (fabs(r[3]+1.0)>10.0*DBL_EPSILON) ? DGDS3 : 0.0;

            d2gds2[0][0] = (fabs(r[0]+1.0)>10.0*DBL_EPSILON) ? D2GDS0S0 : 0.0;
            d2gds2[0][1] = (fabs(r[0]+1.0)>10.0*DBL_EPSILON && fabs(r[1]+1.0)>10.0*DBL_EPSILON) ? D2GDS0S1 : 0.0;
            d2gds2[1][0] = d2gds2[0][1];
            d2gds2[0][2] = (fabs(r[0]+1.0)>10.0*DBL_EPSILON && fabs(r[2]+1.0)>10.0*DBL_EPSILON) ? D2GDS0S2 : 0.0;
            d2gds2[2][0] = d2gds2[0][2];
            d2gds2[0][3] = (fabs(r[0]+1.0)>10.0*DBL_EPSILON && fabs(r[3]+1.0)>10.0*DBL_EPSILON) ? D2GDS0S3 : 0.0;
            d2gds2[3][0] = d2gds2[0][3];
            d2gds2[1][1] = (fabs(r[1]+1.0)>10.0*DBL_EPSILON) ? D2GDS1S1 : 0.0;
            d2gds2[1][2] = (fabs(r[1]+1.0)>10.0*DBL_EPSILON && fabs(r[2]+1.0)>10.0*DBL_EPSILON) ? D2GDS1S2 : 0.0;
            d2gds2[2][1] = d2gds2[1][2];
            d2gds2[1][3] = (fabs(r[1]+1.0)>10.0*DBL_EPSILON && fabs(r[3]+1.0)>10.0*DBL_EPSILON) ? D2GDS1S3 : 0.0;
            d2gds2[3][1] = d2gds2[1][3];
            d2gds2[2][2] = (fabs(r[2]+1.0)>10.0*DBL_EPSILON) ? D2GDS2S2 : 0.0;
            d2gds2[2][3] = (fabs(r[2]+1.0)>10.0*DBL_EPSILON && fabs(r[3]+1.0)>10.0*DBL_EPSILON) ? D2GDS2S3 : 0.0;
            d2gds2[3][2] = d2gds2[2][3];
            d2gds2[3][3] = (fabs(r[3]+1.0)>10.0*DBL_EPSILON) ? D2GDS3S3 : 0.0;

            for (i=0; i<NS; i++) sOld[i] = s[i];

            /* original: gaussj(ptToD2gds2, NS, (double **) NULL, 0);
            if (fabs(r[0]+1.0)<10.0*DBL_EPSILON)d2gds2[0][0]=0.0;
            if (fabs(r[1]+1.0)<10.0*DBL_EPSILON)d2gds2[1][1]=0.0;
            if (fabs(r[2]+1.0)<10.0*DBL_EPSILON)d2gds2[2][2]=0.0;
            if (fabs(r[3]+1.0)<10.0*DBL_EPSILON)d2gds2[3][3]=0.0;
            for (i=0; i<NS; i++) {
                for(j=0; j<NS; j++) s[i] += - d2gds2[i][j]*dgds[j];
            }
            */

            gsl_matrix_scale(ptToD2gds2, -1.0);

            if (USE_SVD) {
                gsl_linalg_SV_decomp(ptToD2gds2, ptToV, ptToS, &vvToSCorr.vector);
                gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToDgds.vector, &vvToSCorr.vector);
            }
            else {
                melts_LU_decomp(ptToD2gds2, indexD2gds2, &signum);
                melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToDgds.vector, &vvToSCorr.vector);
            }

            if (fabs(r[0]+1.0)<10.0*DBL_EPSILON) sCorr[0]=0.0;
            if (fabs(r[1]+1.0)<10.0*DBL_EPSILON) sCorr[1]=0.0;
            if (fabs(r[2]+1.0)<10.0*DBL_EPSILON) sCorr[2]=0.0;
            if (fabs(r[3]+1.0)<10.0*DBL_EPSILON) sCorr[3]=0.0;

            for (i=0; i<NS; i++) s[i] += sCorr[i];

            /* Canabalised from spinel.c */
            lambda = 1.0;

            /* xm1mn= (r[0]-s[0]+1.0)/2.0; */
            if      ((sCorr[0]/2.0) != 0.0 && xm1mn-lambda*(sCorr[0]/2.0) < 0.0)
                lambda = xm1mn/(sCorr[0]/2.0);
            else if ((sCorr[0]/2.0) != 0.0 && xm1mn-lambda*(sCorr[0]/2.0) > 1.0)
                lambda = (xm1mn-1.0)/(sCorr[0]/2.0);

            /* xm2mn= (r[0]+s[0]+1.0)/2.0; */
            if      ((sCorr[0]/2.0) != 0.0 && xm2mn+lambda*(sCorr[0]/2.0) < 0.0)
                lambda = -xm2mn/(sCorr[0]/2.0);
            else if ((sCorr[0]/2.0) != 0.0 && xm2mn+lambda*(sCorr[0]/2.0) > 1.0)
                lambda = (1.0-xm2mn)/(sCorr[0]/2.0);

            /* xm1fe= (r[1]-s[1]+1.0)/2.0; */
            if      ((sCorr[1]/2.0) != 0.0 && xm1fe-lambda*(sCorr[1]/2.0) < 0.0)
                lambda = xm1fe/(sCorr[1]/2.0);
            else if ((sCorr[1]/2.0) != 0.0 && xm1fe-lambda*(sCorr[1]/2.0) > 1.0)
                lambda = (xm1fe-1.0)/(sCorr[1]/2.0);

            /* xm2fe= (r[1]+s[1]+1.0)/2.0; */
            if      ((sCorr[1]/2.0) != 0.0 && xm2fe+lambda*(sCorr[1]/2.0) < 0.0)
                lambda = -xm2fe/(sCorr[1]/2.0);
            else if ((sCorr[1]/2.0) != 0.0 && xm2fe+lambda*(sCorr[1]/2.0) > 1.0)
                lambda = (1.0-xm2fe)/(sCorr[1]/2.0);

            /* xm1co= (r[2]+s[2]+1.0)/2.0; */
            if      ((sCorr[2]/2.0) != 0.0 && xm1co+lambda*(sCorr[2]/2.0) < 0.0)
                lambda = -xm1co/(sCorr[2]/2.0);
            else if ((sCorr[2]/2.0) != 0.0 && xm1co+lambda*(sCorr[2]/2.0) > 1.0)
                lambda = (1.0-xm1co)/(sCorr[2]/2.0);

            /* xm2co= (r[2]-s[2]+1.0)/2.0; */
            if      ((sCorr[2]/2.0) != 0.0 && xm2co-lambda*(sCorr[2]/2.0) < 0.0)
                lambda = xm2co/(sCorr[2]/2.0);
            else if ((sCorr[2]/2.0) != 0.0 && xm2co-lambda*(sCorr[2]/2.0) > 1.0)
                lambda = (xm2co-1.0)/(sCorr[2]/2.0);

            /* xm1ni= (r[3]+s[3]+1.0)/2.0; */
            if      ((sCorr[3]/2.0) != 0.0 && xm1ni+lambda*(sCorr[3]/2.0) < 0.0)
                lambda = -xm1ni/(sCorr[3]/2.0);
            else if ((sCorr[3]/2.0) != 0.0 && xm1ni+lambda*(sCorr[3]/2.0) > 1.0)
                lambda = (1.0-xm1ni)/(sCorr[3]/2.0);

            /* xm2ni= (r[3]-s[3]+1.0)/2.0; */
            if      ((sCorr[3]/2.0) != 0.0 && xm2ni-lambda*(sCorr[3]/2.0) < 0.0)
                lambda = xm2ni/(sCorr[3]/2.0);
            else if ((sCorr[3]/2.0) != 0.0 && xm2ni-lambda*(sCorr[3]/2.0) > 1.0)
                lambda = (xm2ni-1.0)/(sCorr[3]/2.0);

            /* xm2ca=  r[4]; */

            /* xm1mg = (1.0 - xm1mn - xm1fe - xm1co - xm1ni);
             xm1mg = 1.0 - (r[0]-s[0]+1.0)/2.0 - (r[1]-s[1]+1.0)/2.0 - (r[2]+s[2]+1.0)/2.0 - (r[3]+s[3]+1.0)/2.0; */
            if      (((sCorr[0] + sCorr[1] - sCorr[2] - sCorr[3])/2.0) != 0.0 && xm1mg+lambda*((sCorr[0] + sCorr[1] - sCorr[2] - sCorr[3])/2.0) < 0.0)
                lambda = -xm1mg/((sCorr[0] + sCorr[1] - sCorr[2] - sCorr[3])/2.0);
            else if (((sCorr[0] + sCorr[1] - sCorr[2] - sCorr[3])/2.0) != 0.0 && xm1mg+lambda*((sCorr[0] + sCorr[1] - sCorr[2] - sCorr[3])/2.0) > 1.0)
                lambda = (1.0-xm1mg)/((sCorr[0] + sCorr[1] - sCorr[2] - sCorr[3])/2.0);

            /* xm2mg = (1.0 - xm2mn - xm2fe - xm2co - xm2ni - xm2ca);
             xm2mg = 1.0 - (r[0]+s[0]+1.0)/2.0 - (r[1]+s[1]+1.0)/2.0 - (r[2]-s[2]+1.0)/2.0 - (r[3]-s[3]+1.0)/2.0 - r[4]; */
            if      (((sCorr[2] + sCorr[3] - sCorr[0] - sCorr[1])/2.0) != 0.0 && xm2mg+lambda*((sCorr[2] + sCorr[3] - sCorr[0] - sCorr[1])/2.0) < 0.0)
                lambda = -xm2mg/((sCorr[2] + sCorr[3] - sCorr[0] - sCorr[1])/2.0);
            else if (((sCorr[2] + sCorr[3] - sCorr[0] - sCorr[1])/2.0) != 0.0 && xm2mg+lambda*((sCorr[2] + sCorr[3] - sCorr[0] - sCorr[1])/2.0) > 1.0)
                lambda = (1.0-xm2mg)/((sCorr[2] + sCorr[3] - sCorr[0] - sCorr[1])/2.0);

            /* Modify steplength if required to maintain feasibility */
            if (lambda < 1.0) for (i=0; i<NS; i++) s[i] = sOld[i] + lambda*sCorr[i];

            for (i=0; i<NS; i++) sNew[i] = s[i];
            iter++;
        }
        tOld = t;
        pOld = p;
        for (i=0; i<NR; i++) rOld[i] = r[i];

#ifdef DEBUG
        for (i=0; i<NS; i++) {
            if (dgds[i] > sqrt(DBL_EPSILON) && ABS(sNew[i]) > sqrt(DBL_EPSILON)) {
                printf("ERROR in OLIVINE.C (function ORDER). Failed to converge!\n");
                printf("  r1    = %13.6g, r2    = %13.6g, r3    = %13.6g\n", r[0], r[1], r[2]);
                printf("  r4    = %13.6g, r5    = %13.6g\n",   r[3], r[4]);
                printf("  s1    = %13.6g, s2    = %13.6g\n",   sOld[0], sOld[1]);
                printf("  s3    = %13.6g, s4    = %13.6g\n",   sOld[2], sOld[3]);
                printf("  dgds1 = %13.6g, dgds2 = %13.6g\n",   dgds[0], dgds[1]);
                printf("  dgds3 = %13.6g,   dgds4 = %13.6g\n", dgds[2], dgds[3]);

                printf("xm1 mn = %13.6g, xm2 mn = %13.6g\n", xm1mn, xm2mn);
                printf("xm1 fe = %13.6g, xm2 fe = %13.6g\n", xm1fe, xm2fe);
                printf("xm1 co = %13.6g, xm2 co = %13.6g\n", xm1co, xm2co);
                printf("xm1 ni = %13.6g, xm2 ni = %13.6g\n", xm1ni, xm2ni);
                printf("xm1 mg = %13.6g, xm2 mg = %13.6g\n", xm1mg, xm2mg);
                printf("xm2 ca = %13.6g\n",                  xm2ca);
                break;
            }
        }
#endif

        setTOld(tOld);
        setPOld(pOld);
        /* arrays (rOld, sOld, d2gds2, etc) should be preserved automatically */

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
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdrds.vector, &vvToDr.vector);
            else melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdrds.vector, &vvToDr.vector);
        }
    }
    if (mask & THIRD  ) {   /* compute ds/dt:  */
        double d2gdsdt[NS];
        gsl_vector_view vvToDt = gsl_vector_view_array(dt, (size_t) NS),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NS);

        fillD2GDSDT

        /* original: dt[i] += - d2gds2[i][j]*d2gdsdt[j]; */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdt.vector, &vvToDt.vector);
        else melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdt.vector, &vvToDt.vector);
    }
    if (mask & FOURTH ) {   /* compute ds/dp:  */
        double *s = sOld;
        double d2gdsdp[NS];
        gsl_vector_view vvToDp = gsl_vector_view_array(dp, (size_t) NS),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NS);

        fillD2GDSDP

        /* original: dp[i] += - d2gds2[i][j]*d2gdsdp[j]; */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdp.vector, &vvToDp.vector);
        else melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdp.vector, &vvToDp.vector);
    }
    if (mask & FIFTH  ) {   /* compute d2s/dr2 */
        double d2gdrds[NR][NS];
        double d3gdr2ds[NR][NR][NS];
        double d3gdrds2[NR][NS][NS];
        double d3gds3[NS][NS][NS];
        double dsdr[NS][NR], temp[NS], temp2[NS];
        gsl_matrix_view mvToDsdr = gsl_matrix_view_array((double *) dsdr, (size_t) NS, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NS);
        gsl_vector_view vvToTemp = gsl_vector_view_array(temp, (size_t) NS), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NS);
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
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdrds.vector, &vvToDsdr.vector);
            else melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdrds.vector, &vvToDsdr.vector);
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
                if (USE_SVD) {
                    gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
                    for (l=0; l<NS; l++) dr2[l][j][k] = temp2[l];
                }
                else {
                	melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
                	for (l=0; l<NS; l++) dr2[l][j][k] = temp[l];
                }
            }
        }

    }
    if (mask & SIXTH  ) {   /* compute d2s/drt */
        double d2gdrds[NR][NS];
        double d2gdsdt[NS];
        double d3gdrds2[NR][NS][NS];
        double d3gdrdsdt[NR][NS];
        double d3gds2dt[NS][NS];
        double d3gds3[NS][NS][NS];
        double dsdr[NS][NR], dsdt[NS], temp[NS], temp2[NS];
        gsl_matrix_view mvToDsdr = gsl_matrix_view_array((double *) dsdr, (size_t) NS, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NS);
        gsl_vector_view vvToDsdt = gsl_vector_view_array(dsdt, (size_t) NS),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NS),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NS), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NS);
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
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdrds.vector, &vvToDsdr.vector);
            else melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdrds.vector, &vvToDsdr.vector);
        }

        /* compute dsdt vector */
        /* original: dsdt[i] += - d2gds2[i][j]*d2gdsdt[j]; */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        else melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdt.vector, &vvToDsdt.vector);

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
            if (USE_SVD){
                gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
                for (k=0; k<NS; k++) drt[k][j] = temp2[k];
            }
            else {
                melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
                for (k=0; k<NS; k++) drt[k][j] = temp[k];
            }
        }

    }

    if (mask & SEVENTH) {   /* compute d2s/drp */
        double *s = sOld;
        double d2gdrds[NR][NS], d2gdsdp[NS], d3gdrds2[NR][NS][NS],
            d3gdrdsdp[NR][NS], d3gds3[NS][NS][NS], d3gds2dp[NS][NS], dsdr[NS][NR],
            dsdp[NS], temp[NS], temp2[NS];
        gsl_matrix_view mvToDsdr = gsl_matrix_view_array((double *) dsdr, (size_t) NS, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NS);
        gsl_vector_view vvToDsdp = gsl_vector_view_array(dsdp, (size_t) NS),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NS),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NS), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NS);
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
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdrds.vector, &vvToDsdr.vector);
            else melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdrds.vector, &vvToDsdr.vector);
        }

        /* compute dsdp vector */
        /* original: dsdp[i] += - d2gds2[i][j]*d2gdsdp[j]; */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdp.vector, &vvToDsdp.vector);
        else melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdp.vector, &vvToDsdp.vector);

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
            if (USE_SVD) {
                gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
                for (k=0; k<NS; k++) drp[k][j] = temp2[k];
            }
            else {
                melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
                for (k=0; k<NS; k++) drp[k][j] = temp[k];
            }
        }

    }
    if (mask & EIGHTH ) {   /* compute d2s/dt2 */
        double d2gdsdt[NS];
        double d3gds3[NS][NS][NS];
        double d3gds2dt[NS][NS];
        double d3gdsdt2[NS];
        double dsdt[NS], temp[NS], temp2[NS];
        gsl_vector_view vvToDsdt = gsl_vector_view_array(dsdt, (size_t) NS),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NS),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NS), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NS);
        int k, l;

        fillD2GDSDT
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDSDT2

        /* compute dsdt vector */
        /* original: dsdt[i] += - d2gds2[i][j]*d2gdsdt[j]; */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        else melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdt.vector, &vvToDsdt.vector);

        /* compute dsdt2 vector */
        for (j=0; j<NS; j++) {
            temp[j] = d3gdsdt2[j];
            for (k=0; k<NS; k++) {
                temp[j] +=  2.0*d3gds2dt[j][k]*dsdt[k];
                for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdt[l];
            }
        }
        /* original: dt2[i] += - d2gds2[i][j]*temp[j]; */
        if (USE_SVD) {
            gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
            for (j=0; j<NS; j++) dt2[j] = temp2[j];
        }
        else {
            melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
            for (j=0; j<NS; j++) dt2[j] = temp[j];
        }

    }
    if (mask & NINTH  ) {   /* compute d2s/dtp */
        double *s = sOld;
        double d2gdsdt[NS];
        double d2gdsdp[NS];
        double d3gds3[NS][NS][NS];
        double d3gds2dt[NS][NS];
        double d3gds2dp[NS][NS];
        double d3gdsdtdp[NS];
        double dsdt[NS], dsdp[NS], temp[NS], temp2[NS];
        gsl_vector_view vvToDsdt = gsl_vector_view_array(dsdt, (size_t) NS),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NS),
            vvToDsdp = gsl_vector_view_array(dsdp, (size_t) NS),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NS),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NS), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NS);
        int k, l;

        fillD2GDSDT
        fillD2GDSDP
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDSDTDP
        fillD3GDS2DP

        /* compute dsdt vector */
        /* original: dsdt[i] += - d2gds2[i][j]*d2gdsdt[j]; */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        else melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdt.vector, &vvToDsdt.vector);

        /* compute dsdp vector */
        /* original: dsdp[i] += - d2gds2[i][j]*d2gdsdp[j]; */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdp.vector, &vvToDsdp.vector);
        else melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdp.vector, &vvToDsdp.vector);

        /* compute dsdtp vector */
        for (j=0; j<NS; j++) {
            temp[j] = d3gdsdtdp[j];
            for (k=0; k<NS; k++) {
                temp[j] += d3gds2dt[j][k]*dsdp[k] + d3gds2dp[j][k]*dsdt[k];
                for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdp[l];
            }
        }
        /* original: dtp[i] += - d2gds2[i][j]*temp[j]; */
        if (USE_SVD) {
            gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
            for (j=0; j<NS; j++) dtp[j] = temp2[j];
        }
        else {
            melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
            for (j=0; j<NS; j++) dtp[j] = temp[j];
        }

    }
    if (mask & TENTH  ) {   /* compute d2s/dp2 */
        double *s = sOld;
        double d2gdsdp[NS], d3gds3[NS][NS][NS], d3gds2dp[NS][NS], d3gdsdp2[NS],
            dsdp[NS], temp[NS], temp2[NS];
        gsl_vector_view vvToDsdp = gsl_vector_view_array(dsdp, (size_t) NS),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NS),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NS), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NS);
        int k, l;

        fillD2GDSDP
        fillD3GDS3
        fillD3GDS2DP
        fillD3GDSDP2

        /* compute dsdp vector */
        /* original: dsdp[i] += - d2gds2[i][j]*d2gdsdp[j]; */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdp.vector, &vvToDsdp.vector);
        else melts_LU_solve(ptToD2gds2, indexD2gds2, &vvToD2gdsdp.vector, &vvToDsdp.vector);

        /* compute dsdp2 vector */
        for (j=0; j<NS; j++) {
            temp[j] = d3gdsdp2[j];
            for (k=0; k<NS; k++) {
                temp[j] +=  2.0*d3gds2dp[j][k]*dsdp[k];
                for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdp[k]*dsdp[l];
            }
        }
        /* original: dp2[i] += - d2gds2[i][j]*temp[j]; */
        if (USE_SVD) {
            gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
            for (j=0; j<NS; j++) dp2[j] = temp2[j];
        }
        else {
            melts_LU_svx(ptToD2gds2, indexD2gds2, &vvToTemp.vector);
            for (j=0; j<NS; j++) dp2[j] = temp[j];
        }

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
testOlv(int mask, double t, double p,
    int na,          /* Expected number of endmember components                 */
    int nr,          /* Expected number of independent compositional variables  */
    char **names,    /* array of strings of names of endmember components       */
    char **formulas, /* array of strings of formulas of endmember components    */
    double *r,       /* array of indepependent compos variables, check bounds   */
    double *m)       /* array of moles of endmember components, check bounds    */
{
    const char *phase = "olivine.c";
    const char *NAMES[NA]    = { "tephroite", "fayalite", "co-olivine", "ni-olivine",
                               "monticellite","forsterite" };
    const char *FORMULAS[NA] = { "Mn2SiO4", "Fe2SiO4", "Co2SiO4", "Ni2SiO4",
                               "CaMgSiO4", "Mg2SiO4" };
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
        result = result && (r[0] >= -1.0) && (r[0] <= 1.0-r[4]);
        result = result && (r[1] >= -1.0) && (r[1] <= 1.0-r[4]);
        result = result && (r[2] >= -1.0) && (r[2] <= 1.0-r[4]);
        result = result && (r[3] >= -1.0) && (r[3] <= 1.0-r[4]);
        result = result && (r[4] >= 0.0) && (r[4] <= 1.0);
        result = result && (r[0]+r[1]+r[2]+r[3]+2.0*r[4] <= -2.0-r[4]);
        result = result && (r[0]+r[1]+r[2]+r[3]+2.0*r[4] >= -4.0);
    }
    /* Check bounds on moles of endmember components */
    if (mask & SIXTH) {
        for (i=0, sum=0.0; i<NA; i++) sum += m[i];
        result = result && (sum >= 0.0);

        if (sum > 0.0) {
            result = result && (m[0]/sum >= 0.0) && (m[0]/sum <=  1.0);
            result = result && (m[1]/sum >= 0.0) && (m[1]/sum <=  1.0);
            result = result && (m[2]/sum >= 0.0) && (m[2]/sum <=  1.0);
            result = result && (m[3]/sum >= 0.0) && (m[3]/sum <=  1.0);
            result = result && (m[4]/sum >= 0.0) && (m[4]/sum <=  1.0);
            result = result && (m[5]/sum >= -0.5*m[4]/sum) && (m[5]/sum <= 1.0);
        }
    }

    return result;
}

void
conOlv(int inpMask, int outMask, double t, double p,
    double *e,      /* comp of olivine in moles of elements                     */
    double *m,      /* comp of olivine in moles of endmember components         */
    double *r,      /* comp of olivine in terms of the independent comp var     */
    double *x,      /* comp of olivine in mole fractions of endmember comp      */
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
    (3)  THIRD            FOURTH | SEVENTH

    (1) converts a vector of moles of elements into a vector of moles of
            endmember olivine components.
    (2) calculates from a vector of moles of endmember components, one or
            all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
    (3) calculates from a vector of independent compositional variables
            mole fractions of endmember components and/or the Jacobian matrix
            dx[]/dr[]

    In this routine it is assumed that the elements are in the order of atomic
    numbers and that the order of olivine components has been verified as:
       m[0] = tephroite    (Mn2SiO4) ,
       m[1] = fayalite     (Fe2SiO4) ,
       m[2] = Co-olivine   (Co2SiO4) ,
       m[3] = Ni-olivine   (Ni2SiO4) ,
       m[4] = monticellite (CaMgSiO4) and
       m[5] = forsterite   (Mg2SiO4)

    ----------------------------------------------------------------------------*/

    int i, j, k;

    if (inpMask == FIRST && outMask == SECOND) {
        /* Converts a vector of moles of elements into a vector of moles of
       end-member components.                                                 */
        static const int Mg = 12;
        static const int Ca = 20;
        static const int Mn = 25;
        static const int Fe = 26;
        static const int Co = 27;
        static const int Ni = 28;

        m[0] =  e[Mn]/2.0;        /* moles of Mn2SiO4                        */
        m[1] =  e[Fe]/2.0;        /* moles of Fe2SiO4                        */
        m[2] =  e[Co]/2.0;        /* moles of Co2SiO4                        */
        m[3] =  e[Ni]/2.0;        /* moles of Ni2SiO4                        */
        m[4] =  e[Ca];            /* Moles of CaMgSiO4                       */
        m[5] = (e[Mg]-e[Ca])/2.0; /* Moles of Mg2SiO4                        */

    } else if (inpMask == SECOND) {
        double sum;

        if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
            printf("Illegal call to conOlv with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        for (i=0, sum=0.0; i<NA; i++) sum += m[i];

        if (outMask & THIRD) {
            /* Converts a vector of moles of end-member components (m) into a vector
         of independent compositional variables                               */
            r[0] = (sum != 0.0) ? 2.0*m[0]/sum - 1.0 : 0.0;
            r[1] = (sum != 0.0) ? 2.0*m[1]/sum - 1.0 : 0.0;
            r[2] = (sum != 0.0) ? 2.0*m[2]/sum - 1.0 : 0.0;
            r[3] = (sum != 0.0) ? 2.0*m[3]/sum - 1.0 : 0.0;
            r[4] = (sum != 0.0) ? m[4]/sum           : 0.0;
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
         dm[0][j] = (j == 0) ? 2.0*(1.0-m[0]/sum)/sum : -2.0*m[0]/SQUARE(sum);
         dm[1][j] = (j == 1) ? 2.0*(1.0-m[1]/sum)/sum : -2.0*m[1]/SQUARE(sum);
         dm[2][j] = (j == 2) ? 2.0*(1.0-m[2]/sum)/sum : -2.0*m[2]/SQUARE(sum);
         dm[3][j] = (j == 3) ? 2.0*(1.0-m[3]/sum)/sum : -2.0*m[3]/SQUARE(sum);
         dm[4][j] = (j == 4) ?     (1.0-m[4]/sum)/sum : -1.0*m[4]/SQUARE(sum);
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
                        d2m[0][j][k]  = 4.0*m[0]/CUBE(sum);
                        d2m[0][j][k] -= (j == 0) ? 2.0/SQUARE(sum) : 0.0;
                        d2m[0][j][k] -= (k == 0) ? 2.0/SQUARE(sum) : 0.0;
                        d2m[1][j][k]  = 4.0*m[1]/CUBE(sum);
                        d2m[1][j][k] -= (j == 1) ? 2.0/SQUARE(sum) : 0.0;
                        d2m[1][j][k] -= (k == 1) ? 2.0/SQUARE(sum) : 0.0;
                        d2m[2][j][k]  = 4.0*m[2]/CUBE(sum);
                        d2m[2][j][k] -= (j == 2) ? 2.0/SQUARE(sum) : 0.0;
                        d2m[2][j][k] -= (k == 2) ? 2.0/SQUARE(sum) : 0.0;
                        d2m[3][j][k]  = 4.0*m[3]/CUBE(sum);
                        d2m[3][j][k] -= (j == 3) ? 2.0/SQUARE(sum) : 0.0;
                        d2m[3][j][k] -= (k == 3) ? 2.0/SQUARE(sum) : 0.0;
                        d2m[4][j][k]  = 2.0*m[4]/CUBE(sum);
                        d2m[4][j][k] -= (j == 4) ? 1.0/SQUARE(sum) : 0.0;
                        d2m[4][j][k] -= (k == 4) ? 1.0/SQUARE(sum) : 0.0;
                    }
                }
            }

        }

        if (outMask & EIGHTH) {
            /* Calculates the matrix d3r[i]/dm[j]dm[k]dm[l] using m[] as input      */
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
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++) {
                        for (k=0; k<NA; k++) {
                            for (l=0; l<NA; l++) {
                                d3m[i][j][k][l] = -12.0*m[i]/QUARTIC(sum);
                                d3m[i][j][k][l] += (i == j) ? 4.0/CUBE(sum) : 0.0;
                                d3m[i][j][k][l] += (i == k) ? 4.0/CUBE(sum) : 0.0;
                                d3m[i][j][k][l] += (i == l) ? 4.0/CUBE(sum) : 0.0;
                                if (i == 4) d3m[i][j][k][l] /= 2.0;
                            }
                        }
                    }
                }
            }
        }

    } else if (inpMask == THIRD) {

        if (outMask & ~(FOURTH | SEVENTH))
            printf("Illegal call to conOlv with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        if (outMask & FOURTH) {
            /* Converts a vector of independent compositional variables (r) into
                a vector of mole fractions of endmember components (x).              */

            x[0] = (1.0+r[0])/2.0;
            x[1] = (1.0+r[1])/2.0;
            x[2] = (1.0+r[2])/2.0;
            x[3] = (1.0+r[3])/2.0;
            x[4] = r[4];
            x[5] = -1.0 - (r[0]+r[1]+r[2]+r[3])/2.0 - r[4];

            for (i=0; i<NA; i++) if (fabs(x[i]) < DBL_EPSILON) x[i] = 0.0;
        }

        if (outMask & SEVENTH) {
            /* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j]                 */
            for (i=0; i<NA; i++) for (j=0; j<NR; j++) dr[i][j] = 0.0;
            dr[0][0] =  0.5;
            dr[1][1] =  0.5;
            dr[2][2] =  0.5;
            dr[3][3] =  0.5;
            dr[4][4] =  1.0;

            dr[5][0] = -0.5;
            dr[5][1] = -0.5;
            dr[5][2] = -0.5;
            dr[5][3] = -0.5;
            dr[5][4] = -1.0;
        }

    } else  {
        printf("Illegal call to conOlv with inpMask = %o and outMask = %o\n",
            inpMask, outMask);
    }

}
void
dispOlv(int mask, double t, double p, double *x,
    char **formula            /* Mineral formula for interface display MASK: 1 */
    )
{
    double *r = x;
    static char masterString[] = {
/*             1111111111222222222233333333334444444444555555555566666666667
     01234567890123456789012345678901234567890123456789012345678901234567890 */
        "(Ca_.__Mg_.__Fe''_.__Mn_.__Co_.__Ni_.__)2SiO4" };

    if (mask & FIRST) {
        char *string = strcpy((char *) malloc((size_t) (strlen(masterString)+1)*sizeof(char)), masterString);
        double totCa, totFe2, totMg,totMn,totCo,totNi;
        char n[5];
        int i;

        totCa  = 0.5*r[4];
        totFe2 = (1.0+r[1])/2.0;
        totMn  = (1.0+r[0])/2.0;
        totCo  = (1.0+r[2])/2.0;
        totNi  = (1.0+r[3])/2.0;

        totMg  = (1.0-totFe2-totMn-totCo-totNi);

        (void) snprintf(n, 5, "%4.2f", totCa);  for (i=0; i<4; i++) string[ 3+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totMg);  for (i=0; i<4; i++) string[ 9+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totFe2); for (i=0; i<4; i++) string[17+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totMn);  for (i=0; i<4; i++) string[23+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totCo);  for (i=0; i<4; i++) string[29+i] = n[i];
        (void) snprintf(n, 5, "%4.2f", totNi);  for (i=0; i<4; i++) string[35+i] = n[i];

        *formula = string;
    }
}

void
actOlv(int mask, double t, double p, double *x,
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
     fr[i][0] = FR0(i);
     fr[i][1] = FR1(i);
     fr[i][2] = FR2(i);
     fr[i][3] = FR3(i);
     fr[i][4] = FR4(i);

    }

    order(FIRST, t, p, r, s,NULL,NULL,NULL,
                    NULL,NULL,NULL,NULL,NULL,NULL);
    GET_SITE_FRACTIONS

    g       = G;
    dgdr[0] = DGDR0;
    dgdr[1] = DGDR1;
    dgdr[2] = DGDR2;
    dgdr[3] = DGDR3;
    dgdr[4] = DGDR4;

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
        double d2gdr2[NR][NR];
        double d2gdrds[NR][NS];
        double d2gds2[NS][NS];
        double dsdr[NS][NR], dfrdr[NA][NR], gs[NA][NS], dgsds[NA][NS], sum;
        int k, l;

        fillD2GDR2
        fillD2GDRDS
        fillD2GDS2

        for(i=0; i<NA; i++) {
       gs[i][0] = GGS0(i); /* s1 */
       gs[i][1] = GGS1(i); /* s2 */
       gs[i][2] = GGS2(i); /* s3 */
       gs[i][3] = GGS3(i); /* s4 */
       dfrdr[i][0] = DFR0DR0(i); /* r1 */
       dfrdr[i][1] = DFR1DR1(i); /* r2 */
       dfrdr[i][2] = DFR2DR2(i); /* r3 */
       dfrdr[i][3] = DFR3DR3(i); /* r4 */
       dfrdr[i][4] = DFR4DR4(i); /* r5 */

       dgsds[i][0] = DGS0DS0(i); /* s1 */
       dgsds[i][1] = DGS1DS1(i); /* s2 */
       dgsds[i][2] = DGS2DS2(i); /* s3 */
       dgsds[i][3] = DGS3DS3(i); /* s4 */
        }

        order(SECOND, t, p, r,NULL, dsdr,NULL,NULL,
                        NULL,NULL,NULL,NULL,NULL,NULL);

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
        /* implement exclusion criteria on quantities for preclb routines         */
        static const double exclusion[NA] = {
       0.05,  /* 0.0010 exclusion criteria on the mole fraction of Mn2SiO4  */
       0.05,  /* 0.9900 exclusion criteria on the mole fraction of Fe2SiO4  */
       0.05,  /* 0.0003 exclusion criteria on the mole fraction of Co2SiO4  */
       0.05,  /* 0.0003 exclusion criteria on the mole fraction of Ni2SiO4  */
       0.05,  /* 0.9900 exclusion criteria on the mole fraction of CaMgSiO4 */
       0.05   /* 0.9900 exclusion criteria on the mole fraction of Mg2SiO4  */
        };

        double x[NA];

        x[0] = (r[0]+1.0)/2.0;                       /* total MN2SIO4  */
        x[1] = (r[1]+1.0)/2.0;                       /* total FE2SIO4  */
        x[2] = (r[2]+1.0)/2.0;                       /* total CO2SIO4  */
        x[3] = (r[3]+1.0)/2.0;                       /* total NI2SIO4  */
        x[4] = r[4];                                 /* total CAMGSIO4 */
        x[5] = 1.0-((r[0]+1.0)+(r[1]+1.0)+(r[2]+1.0)+(r[3]+1.0))/2.0-r[4]; /* total MG2SIO4 */

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
gmixOlv(int mask, double t, double p, double *x,
    double *gmix, /* Gibbs energy of mixing             BINARY MASK: 0001 */
    double *dx,   /* (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010 */
    double **dx2, /* (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100 */
    double ***dx3 /* (pointer to dx3[][][]) d3(g)/d(x[])3 NARY MASK: 1000 */
    )
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS];

    order(FIRST, t, p, r,s,NULL, NULL, NULL,
                NULL, NULL, NULL, NULL,NULL, NULL);
    GET_SITE_FRACTIONS

    if (mask & FIRST) {

        *gmix = G;
    }

    if(mask & SECOND) {

        dx[0] = DGDR0;
        dx[1] = DGDR1;
        dx[2] = DGDR2;
        dx[3] = DGDR3;
        dx[4] = DGDR4;
    }

    if(mask & THIRD) {
        double d2gdr2[NR][NR];
        double d2gdrds[NR][NS];
        double d2gds2[NS][NS];
        double dsdr[NS][NR];
        int i, j, k, l;

        fillD2GDR2
        fillD2GDRDS
        fillD2GDS2

        order(SECOND, t, p, r,
                NULL,dsdr,     NULL, NULL,
                NULL, NULL, NULL, NULL,
                NULL, NULL);

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
hmixOlv(int mask, double t, double p, double *x,
    double *hmix /* Enthalpy of mixing BINARY MASK: 1 */
    )
/*  This function calculates enthalpy of mixing corrected
        to 1 bar. */
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS];

    order(FIRST, t, p, r,
                s,               NULL, NULL, NULL,
                NULL, NULL, NULL, NULL,
                NULL, NULL);
    GET_SITE_FRACTIONS

    *hmix = (G) + t*(S);
}

void
smixOlv(int mask, double t, double p, double *x,
    double *smix, /* Entropy of mixing                  BINARY MASK: 001 */
    double *dx,   /* (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010 */
    double **dx2  /* (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100 */
    )
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS];

    order(FIRST, t, p, r,s,NULL, NULL, NULL,
                NULL, NULL, NULL, NULL,NULL, NULL);
    GET_SITE_FRACTIONS

    if (mask & FIRST) {
       *smix = (S);
    }

    if(mask & SECOND) {
        double d2gdrds[NR][NS];
        double d2gdrdt[NR];
        double d2gds2[NS][NS];
        double d2gdsdt[NS];
        double dsdr[NS][NR], dsdt[NS];
        int i, k, l;

        fillD2GDRDS
        fillD2GDRDT
        fillD2GDS2
        fillD2GDSDT

        order(SECOND | THIRD, t, p, r,NULL, dsdr,dsdt,
                    NULL,NULL, NULL, NULL, NULL,NULL, NULL);

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

        double d2gdrds[NR][NS];
        double d2gds2[NS][NS];
        double d2gdsdt[NS];
        double d3gdr2ds[NR][NR][NS];
        double d3gdr2dt[NR][NR];
        double d3gdrds2[NR][NS][NS];
        double d3gdrdsdt[NR][NS];
        double d3gds3[NS][NS][NS];
        double d3gds2dt[NS][NS];

        double dsdr[NS][NR], dsdt[NS], d2sdr2[NS][NR][NR], d2sdrdt[NS][NR];
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
                NULL, dsdr,dsdt,NULL,d2sdr2, d2sdrdt, NULL, NULL,
                NULL, NULL);

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
cpmixOlv(int mask, double t, double p, double *x,
    double *cpmix, /* Heat capacity of mixing         BINARY MASK: 001 */
    double *dt,    /* d(cp)/d(t)                      BINARY MASK: 010 */
    double *dx     /* d(cp)/d(x[])d(t)                BINARY MASK: 100 */
    )
{
    DECLARE_SITE_FRACTIONS
    double *r = x;
    double s[NS], dsdt[NS], d2gdsdt[NS], d2gds2[NS][NS], d2gdt2;
    int i, j;

    order(FIRST | THIRD, t, p, r, s, NULL, dsdt, NULL, NULL, NULL, NULL,
        NULL,NULL, NULL);
    GET_SITE_FRACTIONS

    fillD2GDS2
    fillD2GDSDT
    d2gdt2  = D2GDT2;

    if (mask & FIRST) {

        *cpmix = d2gdt2;
        for (i=0; i<NS; i++) {
            *cpmix += 2.0*d2gdsdt[i]*dsdt[i];
            for (j=0; j<NS; j++) *cpmix += d2gds2[i][j]*dsdt[i]*dsdt[j];
        }
        *cpmix *= -t;

    }

    if(mask & SECOND) {
        double d3gds3[NS][NS][NS];
        double d3gds2dt[NS][NS];
        double d3gdsdt2[NS];
        double d3gdt3 = D3GDT3;
        double d2sdt2[NS], temp;
        int k;

        fillD3GDS3
        fillD3GDS2DT
        fillD3GDSDT2

        order(EIGHTH, t, p, r,NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, d2sdt2,NULL, NULL);

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
vmixOlv(int mask, double t, double p, double *x,
    double *vmix, /* Volume of mixing                BINARY MASK: 0000000001 */
    double *dx,   /* pointer to dx[]) d(v)/d(x[])    BINARY MASK: 0000000010 */
    double **dx2, /* pointer to dx2[][]) d(v)/d(x[])2BINARY MASK: 0000000100 */
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

    order(FIRST, t, p, r,s,NULL, NULL, NULL,
                NULL, NULL, NULL, NULL,NULL, NULL);
    GET_SITE_FRACTIONS

    if (mask & FIRST) {
        *vmix = (DGDP);

    }

    if(mask & SECOND) {
        double d2gdrds[NR][NS];
        double d2gdrdp[NR];
        double d2gds2[NS][NS];
        double d2gdsdp[NS];
        double dsdr[NS][NR], dsdp[NS];
        int i, j, k;

        fillD2GDRDS
        fillD2GDRDP
        fillD2GDS2
        fillD2GDSDP

        order(SECOND | FOURTH, t, p, r,NULL, dsdr,NULL, dsdp,
                NULL, NULL, NULL, NULL,NULL, NULL);

        for (i=0; i<NR; i++) {
            dx[i] = d2gdrdp[i];
            for (j=0; j<NS; j++) {
                dx[i] += d2gdrds[i][j]*dsdp[j] + d2gdsdp[j]*dsdr[j][i];
                for (k=0; k<NS; k++) dx[i] += d2gds2[j][k]*dsdp[j]*dsdr[k][i];
            }
        }

    }

    if(mask & THIRD) {
        double d2gdrds[NR][NS];
        double d2gds2[NS][NS];
        double d2gdsdp[NS];
        double d3gdr2ds[NR][NR][NS];
        double d3gdr2dp[NR][NR];
        double d3gdrds2[NR][NS][NS];
        double d3gdrdsdp[NR][NS];
        double d3gds3[NS][NS][NS];
        double d3gds2dp[NS][NS];
        double dsdr[NS][NR],dsdp[NS],d2sdr2[NS][NR][NR],d2sdrdp[NS][NR];
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
        double d2gds2[NS][NS];
        double d2gdsdt[NS];
        double d2gdsdp[NS];
        double d2gdtdp = D2GDTDP;
        double dsdt[NS], dsdp[NS];
        int i, j;

        fillD2GDS2
        fillD2GDSDT
        fillD2GDSDP

        order(THIRD | FOURTH, t, p, r,NULL, NULL, dsdt,dsdp,
                NULL, NULL, NULL, NULL,NULL, NULL);

        *dt = d2gdtdp;
        for (i=0; i<NS; i++) {
            *dt += d2gdsdt[i]*dsdp[i] + d2gdsdp[i]*dsdt[i];
            for (j=0; j<NS; j++) *dt += d2gds2[i][j]*dsdt[i]*dsdp[j];
        }

    }

    if(mask & FIFTH) {
        double d2gds2[NS][NS];
        double d2gdsdp[NS];
        double d2gdp2 = D2GDP2;
        double dsdp[NS];
        int i,j;

        fillD2GDS2
        fillD2GDSDP

        order(FOURTH, t, p, r,
                NULL, NULL, NULL, dsdp,
                NULL, NULL, NULL, NULL,
                NULL, NULL);

        *dp = d2gdp2;
        for (i=0; i<NS; i++) {
            *dp += 2.0*d2gdsdp[i]*dsdp[i];
            for (j=0; j<NS; j++) *dp += d2gds2[i][j]*dsdp[i]*dsdp[j];
        }

    }

    if(mask & SIXTH) {
        double d2gds2[NS][NS];
        double d2gdsdt[NS];
        double d2gdsdp[NS];
        double d3gds3[NS][NS][NS];
        double d3gds2dt[NS][NS];
        double d3gdsdt2[NS];
        double d3gds2dp[NS][NS];
        double d3gdsdtdp[NS];
        double d3gdt2dp        = D3GDT2DP;
        double dsdt[NS], dsdp[NS], d2sdt2[NS], d2sdtdp[NS];
        int i, j, k;

        fillD2GDS2
        fillD2GDSDT
        fillD2GDSDP
        fillD3GDS3
        fillD3GDS2DT
        fillD3GDSDT2
        fillD3GDS2DP
        fillD3GDSDTDP


        order(THIRD | FOURTH | EIGHTH | NINTH, t, p, r,
                NULL, NULL, dsdt,dsdp,
                NULL, NULL, NULL, d2sdt2,
                d2sdtdp,NULL);

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
        double d2gds2[NS][NS];
        double d2gdsdt[NS];
        double d2gdsdp[NS];
        double d3gds3[NS][NS][NS];
        double d3gds2dt[NS][NS];
        double d3gds2dp[NS][NS];
        double d3gdsdtdp[NS];
        double d3gdsdp2[NS];
        double d3gdtdp2             = D3GDTDP2;
        double dsdt[NS], dsdp[NS], d2sdtdp[NS], d2sdp2[NS];
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
                NULL, NULL, dsdt,dsdp,
                NULL, NULL, NULL, NULL,
                d2sdtdp,d2sdp2);

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
        double d2gds2[NS][NS];
        double d2gdsdp[NS];
        double d3gds3[NS][NS][NS];
        double d3gds2dp[NS][NS];
        double d3gdsdp2[NS];
        double d3gdp3          = D3GDP3;
        double dsdp[NS], d2sdp2[NS];
        int i, j, k;

        fillD2GDS2
        fillD2GDSDP
        fillD3GDS3
        fillD3GDS2DP
        fillD3GDSDP2

        order(FOURTH | TENTH, t, p, r,
                NULL, NULL, NULL, dsdp,
                NULL, NULL, NULL, NULL,
                NULL, d2sdp2);

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

/* end of file OLIVINE.C */
