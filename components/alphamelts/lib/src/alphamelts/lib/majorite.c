const char *majorite_ver(void) { return "$Id: majorite.c,v 1.9 2007/11/29 05:32:13 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: majorite.c,v $
MELTS Source Code: RCS Revision 1.9  2007/11/29 05:32:13  ghiorso
MELTS Source Code: RCS Majorite testMaj corrections.  Fixed "O2" in sol_struct_data.h to "o2".
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.8  2007/09/18 01:42:23  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.7  2007/05/21 21:08:45  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.6  2007/05/14 16:11:37  ghiorso
MELTS Source Code: RCS Reformulation of majorite and change of ss properties of majorite to Berman
MELTS Source Code: RCS and Vinet forms.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.5  2007/05/07 18:23:21  ghiorso
MELTS Source Code: RCS Modifications to LEPR and calibration algorithms following visit by
MELTS Source Code: RCS Marc and Tim.  Mostly work on cpx and opx inclusion and reclassification.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.4  2007/03/07 21:21:58  ghiorso
MELTS Source Code: RCS Revised majorite model and the way garnets are treated during calibration.
MELTS Source Code: RCS Revised calibration XML file to include LEPER + older MELTS/pMELTS data.
MELTS Source Code: RCS
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
MELTS Source Code: RCS Revision 1.1.1.1  2006/08/15 16:57:36  ghiorso
MELTS Source Code: RCS xMELTS gcc 3.x sources
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2004/01/02 19:21:49  cvsaccount
MELTS Source Code: RCS CTserver University of Chicago
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.12  2002/04/17 20:46:06  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.11  2002/04/17 19:43:07  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.10  2002/04/16 04:34:43  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.9  2002/04/16 03:47:50  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.8  2002/04/15 20:47:18  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.7  2002/04/15 04:45:59  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.6  2002/04/13 00:02:45  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.5  2002/04/05 01:40:20  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.3  2002/03/03 22:54:59  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2001/12/31 05:52:19  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Routines to compute garnet solution properties 
**      (file: MAJORITE.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  December 27, 2001 Original Version
**      V2.0-1  Mark S. Ghiorso  April    12, 2002 Revised Model
**      V3.0-1  Mark S. Ghiorso  March     6, 2007 Revised Model
**
*/

#include "silmin.h"  /* Structure definitions foor SILMIN package */

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

/*
 *=============================================================================
 * Majorite solution parameters:
 * Berman and Koziol (see garnet.c) + guesses
 *
 *
 * 1 == Ca3Al2Si3O12, 2 == Mg3Al2Si3O12, 3 == Fe3Al2Si3O12 4 == Mg3MgSiSi3O12 
 *
 */

#define  WH12 ( 45380.0) /* joules W Grossular - Pyrope  	        */
#define dWH12 (-23820.0) /* joules W112 = W12 + dW12, dW12 X1 X2 (X1-X2) */
#define  WH13 ( 11470.0) /* joules W Grossular - Almandine               */
#define dWH13 (  8850.0) /* joules W113 = W13 + dW13, dW13 X1 X3 (X1-X3) */
#define  WH14 (     0.0) /* guess grossular - majorite                   */
#define  WH23 (  1975.0) /* joules W Pyrope    - Almandine               */
#define dWH23 ( -1750.0) /* joules W223 = W23 + dW23, dW23 X2 X3 (X2-X3) */
#define  WH24 (     0.0) /* guess pyrope    - majorite                   */
#define  WH34 (     0.0) /* guess almandine - majorite                   */

#define WS12  ( 18.79 )  /* joules/K */
#define WS13  (  5.08 )  /* joules/K */
#define WS23  (  0.0  )  /* joules/K */

#define  WV12 (  0.10 ) /* (  0.10 ) joules/bar */
#define dWV12 (  0.0  ) /* (  0.0  ) joules/bar */
#define  WV13 (  0.13 ) /* (  0.13 ) joules/bar */
#define dWV13 (  0.04 ) /* (  0.04 ) joules/bar */			     
#define  WV23 (  0.035) /* (  0.035) joules/bar */
#define dWV23 ( -0.025) /* ( -0.025) joules/bar */

#define  WG12  (WH12-t*WS12+p*WV12)
#define  WG13  (WH13-t*WS13+p*WV13)
#define  WG14  (WH14)
#define  WG23  (WH23-t*WS23+p*WV23)
#define  WG24  (WH24) 
#define  WG34  (WH34) 

#define dWG12  (dWH12+p*dWV12)
#define dWG13  (dWH13+p*dWV13)
#define dWG23  (dWH23+p*dWV23)

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conFld defines the conversion from m[i], to r[j]
 */
#define NR         3
#define NS         0
#define NA         4
#define DFR0DR0(i) - 1.0
#define DFR1DR1(i) - 1.0
#define DFR2DR2(i) - 1.0

/*
 * Global (to this file): derivative definitions
 */
#define R         8.3143

#define S         - 3.0*R*(xca*log(xca) + xmg*log(xmg) + xfe*log(xfe)) - 2.0*R*(xsi*log(xsi) + 0.5*xal*log(xal)) \
                  + (1.0-x[0]-x[1]-x[2])*x[0]*WS12 \
		  + (1.0-x[0]-x[1]-x[2])*x[1]*WS13 \
		  + x[0]*x[1]*WS23		  
#define H           (1.0-x[0]-x[1]-x[2])*x[0]*(WH12 + (1.0-2.0*x[0]-x[1]-x[2])*dWH12) \
		  + (1.0-x[0]-x[1]-x[2])*x[1]*(WH13 + (1.0-x[0]-2.0*x[1]-x[2])*dWH13) \
		  + (1.0-x[0]-x[1]-x[2])*x[2]*WH14 \
		  + x[0]*x[1]*(WH23 + (x[0]-x[1])*dWH23) \
		  + x[0]*x[2]*WH24 \
		  + x[1]*x[2]*WH34		  
#define V           (1.0-x[0]-x[1]-x[2])*x[0]*(WV12 + (1.0-2.0*x[0]-x[1]-x[2])*dWV12) \
		  + (1.0-x[0]-x[1]-x[2])*x[1]*(WV13 + (1.0-x[0]-2.0*x[1]-x[2])*dWV13) \
		  + x[0]*x[1]*(WV23 + (x[0]-x[1])*dWV23)		  
#define G         3.0*R*t*(xca*log(xca) + xmg*log(xmg) + xfe*log(xfe)) + 2.0*R*t*(xsi*log(xsi) + 0.5*xal*log(xal)) + \
                  + (1.0-x[0]-x[1]-x[2])*x[0]*(WG12 + (1.0-2.0*x[0]-x[1]-x[2])*dWG12) \
		  + (1.0-x[0]-x[1]-x[2])*x[1]*(WG13 + (1.0-x[0]-2.0*x[1]-x[2])*dWG13) \
		  + (1.0-x[0]-x[1]-x[2])*x[2]*WG14 \
		  + x[0]*x[1]*(WG23 + (x[0]-x[1])*dWG23) \
		  + x[0]*x[2]*WG24 \
		  + x[1]*x[2]*WG34 
		  

#define DGDR0     R*t*(-3.0*log(xca) + 3.0*log(xmg)) \
                  + (1.0-2.0*x[0]-x[1]-x[2])*(WG12 + (1.0-2.0*x[0]-x[1]-x[2])*dWG12) \
		    + (1.0-x[0]-x[1]-x[2])*x[0]*(-2.0*dWG12) \
		  - x[1]*(WG13 + (1.0-x[0]-2.0*x[1]-x[2])*dWG13) + (1.0-x[0]-x[1]-x[2])*x[1]*(-dWG13) \
		  - x[2]*WG14 \
		  + x[1]*(WG23 + (x[0]-x[1])*dWG23) + x[0]*x[1]*dWG23 \
		  + x[2]*WG24
		  
#define DGDR1     R*t*(-3.0*log(xca)+3.0*log(xfe)) \
                  - x[0]*(WG12 + (1.0-2.0*x[0]-x[1]-x[2])*dWG12) + (1.0-x[0]-x[1]-x[2])*x[0]*(-dWG12) \
		  + (1.0-x[0]-2.0*x[1]-x[2])*(WG13 + (1.0-x[0]-2.0*x[1]-x[2])*dWG13) \
		    + (1.0-x[0]-x[1]-x[2])*x[1]*(-2.0*dWG13) \
		  - x[2]*WG14 \
		  + x[0]*(WG23 + (x[0]-x[1])*dWG23) + x[0]*x[1]*(-dWG23) \
		  + x[2]*WG34		    
		  
#define DGDR2     R*t*(-3.0*log(xca) + 3.0*log(xmg) + log(xsi) - log(xal)) \
                  - x[0]*(WG12 + (1.0-2.0*x[0]-x[1]-x[2])*dWG12) + (1.0-x[0]-x[1]-x[2])*x[0]*(-dWG12) \
		  - x[1]*(WG13 + (1.0-x[0]-2.0*x[1]-x[2])*dWG13) + (1.0-x[0]-x[1]-x[2])*x[1]*(-dWG13) \
		  + (1.0-x[0]-x[1]-2.0*x[2])*WG14 \
		  + x[0]*WG24 \
		  + x[1]*WG34		  

#define DGDP      (V)

#define D2GDR0R0  R*t*(3.0/xca + 3.0/xmg) \
                  - 2.0*(WG12 + (1.0-2.0*x[0]-x[1]-x[2])*dWG12) + (1.0-2.0*x[0]-x[1]-x[2])*(-2.0*dWG12) \
		    + (1.0-2.0*x[0]-x[1]-x[2])*(-2.0*dWG12) \
		  - x[1]*(-dWG13) - x[1]*(-dWG13) + x[1]*dWG23 + x[1]*dWG23
#define D2GDR0R1  R*t*(3.0/xca) \
                  - (WG12 + (1.0-2.0*x[0]-x[1]-x[2])*dWG12) + (1.0-2.0*x[0]-x[1]-x[2])*(-dWG12) \
		    - x[0]*(-2.0*dWG12) \
		  - (WG13 + (1.0-x[0]-2.0*x[1]-x[2])*dWG13) - x[1]*(-2.0*dWG13) \
		    - x[1]*(-dWG13) + (1.0-x[0]-x[1]-x[2])*(-dWG13) \
		  + (WG23 + (x[0]-x[1])*dWG23) + x[1]*(-dWG23) + x[0]*dWG23  
#define D2GDR0R2  R*t*(3.0/xca + 3.0/xmg) \
                  - (WG12 + (1.0-2.0*x[0]-x[1]-x[2])*dWG12) + (1.0-2.0*x[0]-x[1]-x[2])*(-dWG12) \
		    - x[0]*(-2.0*dWG12) \
		  - x[1]*(-dWG13) - x[1]*(-dWG13) \
		  - WG14 + WG24		  	  
#define D2GDR0DT  R*(-3.0*log(xca) + 3.0*log(xmg)) \
                  - (1.0-2.0*x[0]-x[1]-x[2])*WS12 + x[1]*WS13 - x[1]*WS23		                  
#define D2GDR0DP  (1.0-2.0*x[0]-x[1]-x[2])*(WV12 + (1.0-2.0*x[0]-x[1]-x[2])*dWV12) \
		    + (1.0-x[0]-x[1]-x[2])*x[0]*(-2.0*dWV12) \
		  - x[1]*(WV13 + (1.0-x[0]-2.0*x[1]-x[2])*dWV13) + (1.0-x[0]-x[1]-x[2])*x[1]*(-dWV13) \
		  + x[1]*(WV23 + (x[0]-x[1])*dWV23) + x[0]*x[1]*dWV23
		  
#define D2GDR1R1  R*t*(3.0/xca + 3.0/xfe) \
                  - x[0]*(-dWG12) - x[0]*(-dWG12) \
		  - 2.0*(WG13 + (1.0-x[0]-2.0*x[1]-x[2])*dWG13) \
		    + (1.0-x[0]-2.0*x[1]-x[2])*(-2.0*dWG13) + (1.0-x[0]-2.0*x[1]-x[2])*(-2.0*dWG13) \
		  + x[0]*(-dWG23) + x[0]*(-dWG23)  
#define D2GDR1R2  R*t*(3.0/xca) \
                  - x[0]*(-dWG12) - x[0]*(-dWG12) \
		  - (WG13 + (1.0-x[0]-2.0*x[1]-x[2])*dWG13) \
		    + (1.0-x[0]-2.0*x[1]-x[2])*(-dWG13) - x[1]*(-2.0*dWG13) \
		  - WG14 + WG34  
#define D2GDR1DT  R*(-3.0*log(xca) + 3.0*log(xfe)) \
                  + x[0]*WS12 \
		  + x[1]*WS13 - (1.0-x[0]-x[1]-x[2])*WS13 \
		  - x[0]*WS23		                   
#define D2GDR1DP  - x[0]*(WV12 + (1.0-2.0*x[0]-x[1]-x[2])*dWV12) + (1.0-x[0]-x[1]-x[2])*x[0]*(-dWV12) \
		  + (1.0-x[0]-2.0*x[1]-x[2])*(WV13 + (1.0-x[0]-2.0*x[1]-x[2])*dWV13) \
		    + (1.0-x[0]-x[1]-x[2])*x[1]*(-2.0*dWV13) \
		  + x[0]*(WV23 + (x[0]-x[1])*dWV23) + x[0]*x[1]*(-dWV23)

#define D2GDR2R2  R*t*(3.0/xca + 3.0/xmg + 0.5/xsi + 1.0/xal) \
                  - x[0]*(-dWG12) - x[0]*(-dWG12) \
		  - x[1]*(-dWG13) - x[1]*(-dWG13) \
		  - 2.0*WG14
#define D2GDR2DT  R*(-3.0*log(xca) + 3.0*log(xmg) + log(xsi) - log(xal)) + x[0]*WS12 + x[1]*WS13     		                               
#define D2GDR2DP  - x[0]*(WV12 + (1.0-2.0*x[0]-x[1]-x[2])*dWV12) + (1.0-x[0]-x[1]-x[2])*x[0]*(-dWV12) \
		  - x[1]*(WV13 + (1.0-x[0]-2.0*x[1]-x[2])*dWV13) + (1.0-x[0]-x[1]-x[2])*x[1]*(-dWV13)

#define D2GDT2    0.0
#define D2GDTDP   0.0
#define D2GDP2    0.0

#define D3GDR0R0R0 R*t*(3.0/SQUARE(xca) - 3.0/SQUARE(xmg)) \
                   - 2.0*(-2.0*dWG12) - 2.0*(-2.0*dWG12) - 2.0*(-2.0*dWG12)
#define D3GDR0R0R1 R*t*(3.0/SQUARE(xca)) \
                   - 2.0*(-dWG12) - (-2.0*dWG12) - (-2.0*dWG12) - (-dWG13) - (-dWG13) + dWG23 + dWG23		  
#define D3GDR0R0R2 R*t*(3.0/SQUARE(xca) - 3.0/SQUARE(xmg)) \
                   - 2.0*(-dWG12) - (-2.0*dWG12) - (-2.0*dWG12)
#define D3GDR0R0DT R*(3.0/xca + 3.0/xmg) + 2.0*WS12
#define D3GDR0R0DP - 2.0*(WV12 + (1.0-2.0*x[0]-x[1]-x[2])*dWV12) + (1.0-2.0*x[0]-x[1]-x[2])*(-2.0*dWV12) \
		    + (1.0-2.0*x[0]-x[1]-x[2])*(-2.0*dWV12) \
		   - x[1]*(-dWV13) - x[1]*(-dWV13) + x[1]*dWV23 + x[1]*dWV23
		   
#define D3GDR0R1R1 R*t*(3.0/SQUARE(xca)) \
                   - (-dWG12) - (-dWG12) - (-2.0*dWG13) - (-2.0*dWG13) - (-dWG13) - (-dWG13) + (-dWG23) + (-dWG23)  	
#define D3GDR0R1R2 R*t*(3.0/SQUARE(xca)) - (-dWG12) - (-dWG12) - (-dWG13) - (-dWG13)
#define D3GDR0R1DT R*(3.0/xca) + WS12 + WS13 - WS23
#define D3GDR0R1DP - (WV12 + (1.0-2.0*x[0]-x[1]-x[2])*dWV12) + (1.0-2.0*x[0]-x[1]-x[2])*(-dWV12) \
		    - x[0]*(-2.0*dWV12) \
		   - (WV13 + (1.0-x[0]-2.0*x[1]-x[2])*dWV13) - x[1]*(-2.0*dWV13) \
		    - x[1]*(-dWV13) + (1.0-x[0]-x[1]-x[2])*(-dWV13) \
		   + (WV23 + (x[0]-x[1])*dWV23) + x[1]*(-dWV23) + x[0]*dWV23	

#define D3GDR0R2R2 R*t*(3.0/SQUARE(xca) - 3.0/SQUARE(xmg)) - (-dWG12) - (-dWG12)
#define D3GDR0R2DT R*(3.0/xca + 3.0/xmg) + WS12
#define D3GDR0R2DP - (WV12 + (1.0-2.0*x[0]-x[1]-x[2])*dWV12) + (1.0-2.0*x[0]-x[1]-x[2])*(-dWV12) \
		    - x[0]*(-2.0*dWV12) - x[1]*(-dWV13) - x[1]*(-dWV13)	

#define D3GDR1R1R1 R*t*(3.0/SQUARE(xca) - 3.0/SQUARE(xfe)) \
		   - 2.0*(-2.0*dWG13) - 2.0*(-2.0*dWG13) - 2.0*(-2.0*dWG13)
#define D3GDR1R1R2 R*t*(3.0/SQUARE(xca)) \
		   - 2.0*(-dWG13) - (-2.0*dWG13) - (-2.0*dWG13)
#define D3GDR1R1DT R*(3.0/xca + 3.0/xfe) + WS13 + WS13
#define D3GDR1R1DP - x[0]*(-dWV12) - x[0]*(-dWV12) \
		   - 2.0*(WV13 + (1.0-x[0]-2.0*x[1]-x[2])*dWV13) \
		    + (1.0-x[0]-2.0*x[1]-x[2])*(-2.0*dWV13) + (1.0-x[0]-2.0*x[1]-x[2])*(-2.0*dWV13) \
		   + x[0]*(-dWV23) + x[0]*(-dWV23) 	

#define D3GDR1R2R2 R*t*(3.0/SQUARE(xca)) + dWG13 - (-dWG13)
#define D3GDR1R2DT R*(3.0/xca) + WS13	
#define D3GDR1R2DP - x[0]*(-dWV12) - x[0]*(-dWV12) \
		   - (WV13 + (1.0-x[0]-2.0*x[1]-x[2])*dWV13) + (1.0-x[0]-2.0*x[1]-x[2])*(-dWV13) - x[1]*(-2.0*dWV13)	

#define D3GDR2R2R2 R*t*(3.0/SQUARE(xca) - 3.0/SQUARE(xmg) - 0.25/SQUARE(xsi) + 1.0/SQUARE(xal))
#define D3GDR2R2DT R*(3.0/xca + 3.0/xmg + 0.5/xsi + 1.0/xal)
#define D3GDR2R2DP - x[0]*(-dWV12) - x[0]*(-dWV12) - x[1]*(-dWV13) - x[1]*(-dWV13)

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

/*
 *=============================================================================
 * Public functions:
 *    mask  -  bitwise mask for selecting output
 *    t     -  Temperature (K)
 *    p     -  Pressure (bars)
 *    *x    -  (pointer to x[]) Array of independent compositional variables
 */
int
testMaj(int mask, double t, double p,
  int na,          /* Expected number of endmember components                 */
  int nr,          /* Expected number of independent compositional variables  */
  char **names,    /* array of strings of names of endmember components       */
  char **formulas, /* array of strings of formulas of endmember components    */
  double *r,       /* array of indepependent compos variables, check bounds   */
  double *m)       /* array of moles of endmember components, check bounds    */
{
  const char *phase = "majorite.c";
  const char *NAMES[NA]    = { "grossular",    "pyrope",       "almandine",    "majorite" };
  const char *FORMULAS[NA] = { "Ca3Al2Si3O12", "Mg3Al2Si3O12", "Fe3Al2Si3O12", "Mg4Si4O12" };
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
    for (i=0, sum=0.0; i<NR; i++) sum += r[i];
    
    result = result && (sum                <= 1.0);
    result = result && (3.0*r[0]+4.0*r[2]  >= 0.0); /* Mg   >= 0.0 */
    result = result && (r[1]               >= 0.0); /* Fe2+ >= 0.0 */
    result = result && (1.0-r[0]-r[1]-r[2] >= 0.0); /* Ca   >= 0.0 */
    result = result && (1.0-r[2]           >= 0.0); /* Al   >= 0.0 */
    result = result && (3.0+r[2]           >= 3.0); /* Si   >= 3.0 */
  }
  /* Check bounds on moles of endmember components */
  if (mask & SIXTH) {
    for (i=0, sum=0.0; i<NA; i++) sum += m[i];

    result = result && (sum               >= 0.0);
    result = result && (3.0*m[1]+4.0*m[3] >= 0.0);     /* Mg   >= 0.0  */
    result = result && (m[2]              >= 0.0);     /* Fe2+ >= 0.0  */
    result = result && (m[0]              >= 0.0);     /* Ca   >= 0.0  */ 
    result = result && (m[0]+m[1]+m[2]    >= 0.0);     /* Al   >= 0.0  */
    result = result && (m[3]              >= 0.0);     /* Si   >= 0.0  */
  }

  return result;
}

void
conMaj(int inpMask, int outMask, double t, double p,
  double *e,      /* comp of garnet in moles of elements                      */
  double *m,      /* comp of garnet in moles of endmember components          */
  double *r,      /* comp of garnet in terms of the independent comp var      */
  double *x,      /* comp of garnet in mole fractions of endmember comp       */
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
      endmember garnet components.
  (2) calculates from a vector of moles of endmember components, one or
      all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
  (3) calculates from a vector of independent compositional variables
      mole fractions of endmember components and/or the Jacobian matrix
      dx[]/dr[]

  In this routine it is assumed that the elements are in the order of atomic 
  numbers and that the order of garnet components has been verified as:
        m[0] = grossular   (Ca3Al2Si3O12),
        m[1] = pyrope      (Mg3Al2Si3O12), 
        m[2] = almandine   (Fe3Al2Si3O12), and
        m[3] = majorite    (Mg3MgSiSi3O12) 

  ----------------------------------------------------------------------------*/

  int i, j, k;

  if (inpMask == FIRST && outMask == SECOND) {
    static const int O  =  8;
    static const int Na = 11;
    static const int Mg = 12;
    static const int Al = 13;
    static const int Si = 14;
    static const int K  = 19;
    static const int Ca = 20;
    static const int Ti = 22;
    static const int Cr = 24;
    static const int Mn = 25;
    static const int Fe = 26;
    static const int Ni = 28;
    
    double sumCat = e[Na] + e[ K] 
                  + e[Mg] + e[Ca] + e[Mn] + e[Fe] + e[Ni] 
                  + e[Al] + e[Cr]
		  + e[Si] + e[Ti];
    double mFe3 = 2.0*(3.0*sumCat/2.0 - (e[O]-e[Fe]) - e[Fe]);
    double mFe2;
    double excessSi = e[Si] - 3.0*sumCat/8.0; /* excess Si = (sumCat/8)(8*e[Si]/sumCat - 3) = e[Si] - 3 sumCat / 8 */

    if (excessSi < 0.0) excessSi = 0.0;
        
    if (mFe3 < 0.0) mFe3 = 0.0;
    mFe2 = e[Fe] - mFe3;
    if (mFe2 < 0.0) { mFe2 = 0.0; mFe3 = e[Fe]; }
                                                	     
    m[0] = (e[Al]/2.0)*e[Ca]/(e[Na]+e[K]+e[Ca]+e[Mg]+e[Mn]+mFe2+e[Ni]); /* moles of Ca3Al2Si3O12 */
    m[1] = (e[Al]/2.0)*e[Mg]/(e[Na]+e[K]+e[Ca]+e[Mg]+e[Mn]+mFe2+e[Ni]); /* Moles of Mg3Al2Si3O12 */
    m[2] = (e[Al]/2.0)* mFe2/(e[Na]+e[K]+e[Ca]+e[Mg]+e[Mn]+mFe2+e[Ni]); /* Moles of Fe3Al2Si3O12 */
    m[3] = excessSi;                                                    /* Moles of *3*SiSi3O12  */

  } else if (inpMask == SECOND) {
    double sum;

    if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
      printf("Illegal call to conMaj with inpMask = %o and outMask = %o\n",
        inpMask, outMask);

    for (i=0, sum=0.0; i<NA; i++) sum += m[i];

    if (outMask & THIRD) {
      if (sum != 0.0) {
        r[0] = m[1]/sum;
	r[1] = m[2]/sum;
	r[2] = m[3]/sum;
      } else for (i=0; i<NR; i++) r[i] = 0.0; 
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
        for (j=0; j<NA; j++) dm[0][j] = (j == 1) ? (1.0-m[1]/sum)/sum : -m[1]/SQUARE(sum);
	for (i=1; i<NR; i++) {
          for (j=0; j<NA; j++) dm[i][j] = (j == (i+1)) ? (1.0-m[i+1]/sum)/sum : -m[i+1]/SQUARE(sum); 
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
	for (j=0; j<NA; j++)  {
          for (k=0; k<NA; k++) {
            d2m[0][j][k]  = 2.0*m[1]/CUBE(sum);
            d2m[0][j][k] -= (1 == j) ? 1.0/SQUARE(sum) : 0.0;
            d2m[0][j][k] -= (1 == k) ? 1.0/SQUARE(sum) : 0.0;
          }
        }
      
        for (i=1; i<NR; i++) {
          for (j=0; j<NA; j++)  {
            for (k=0; k<NA; k++) {
              d2m[i][j][k]  = 2.0*m[i+1]/CUBE(sum);
              d2m[i][j][k] -= ((i+1) == j) ? 1.0/SQUARE(sum) : 0.0;
              d2m[i][j][k] -= ((i+1) == k) ? 1.0/SQUARE(sum) : 0.0;
            }
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
        for (j=0; j<NA; j++) {
          for (k=0; k<NA; k++) {
	    for (l=0; l<NA; l++) {
              d3m[0][j][k][l]  = -6.0*m[1]/QUARTIC(sum);
              d3m[0][j][k][l] += (1 == j) ? 2.0/CUBE(sum) : 0.0;
              d3m[0][j][k][l] += (1 == k) ? 2.0/CUBE(sum) : 0.0;
              d3m[0][j][k][l] += (1 == l) ? 2.0/CUBE(sum) : 0.0;
	    }
          }
        }
      
        for (i=1; i<NR; i++) {
          for (j=0; j<NA; j++) {
            for (k=0; k<NA; k++) {
	      for (l=0; l<NA; l++) {
                d3m[i][j][k][l]  = -6.0*m[i+1]/QUARTIC(sum);
                d3m[i][j][k][l] += ((i+1) == j) ? 2.0/CUBE(sum) : 0.0;
                d3m[i][j][k][l] += ((i+1) == k) ? 2.0/CUBE(sum) : 0.0;
                d3m[i][j][k][l] += ((i+1) == l) ? 2.0/CUBE(sum) : 0.0;
	      }
            }
          }
        }
      }

    }
  } else if (inpMask == THIRD) {

    if (outMask & ~(FOURTH | SEVENTH))
      printf("Illegal call to conMaj with inpMask = %o and outMask = %o\n",
        inpMask, outMask);

    if (outMask & FOURTH) {
      /* Converts a vector of independent compositional variables (r) 
         into a vector of mole fractions of endmember components (x).         */

      x[0] = 1.0 - r[0] - r[1] - r[2];
      x[1] = r[0];
      x[2] = r[1];
      x[3] = r[2];
    }

    if (outMask & SEVENTH) {
      /* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
      dr[0][0] = -1.0; dr[0][1] = -1.0; dr[0][2] = -1.0;
      dr[1][0] =  1.0; dr[1][1] =  0.0; dr[1][2] =  0.0;
      dr[2][0] =  0.0; dr[2][1] =  1.0; dr[2][2] =  0.0;
      dr[3][0] =  0.0; dr[3][1] =  0.0; dr[3][2] =  1.0;
    }

  } else  {
    printf("Illegal call to conMaj with inpMask = %o and outMask = %o\n",
      inpMask, outMask);
  }

}

void
dispMaj(int mask, double t, double p, double *x,
  char **formula            /* Mineral formula for interface display MASK: 1 */
  )
{
  double *r = x;
  static char masterString[] = {
/*             1111111111222222222233333333334444444444555555555566666666667
     01234567890123456789012345678901234567890123456789012345678901234567890 */
    "Ca_.__Fe''_.__Mg_.__Al_.__Si_.__O12" };

  if (mask & FIRST) {
    char *string = strcpy((char *) malloc((size_t) (strlen(masterString)+1)*sizeof(char)), masterString);
    double totalCa, totalFe2, totalMg, totalSi, totalAl;
    char n[5];
    int i;

    totalMg  = 3.0*r[0] + 4.0*r[2];
    totalFe2 = 3.0*r[1];
    totalCa  = 3.0*(1.0 - r[0] - r[1] - r[2]);
    totalAl  = 2.0*(1.0 - r[2]);
    totalSi  = 3.0 + r[2];

    (void) snprintf(n, 5, "%4.2f", totalCa);  for (i=0; i<4; i++) string[ 2+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totalFe2); for (i=0; i<4; i++) string[10+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totalMg);  for (i=0; i<4; i++) string[16+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totalAl);  for (i=0; i<4; i++) string[22+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totalSi);  for (i=0; i<4; i++) string[28+i] = n[i];
 
    *formula = string;
  }
}

void 
actMaj(int mask, double t, double p, double *x, 
  double *a,  /* (pointer to a[]) activities              BINARY MASK: 0001 */
  double *mu, /* (pointer to mu[]) chemical potentials    BINARY MASK: 0010 */
  double **dx /* (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 0100 */
  )           /* exclusion criteria applied to results if BINARY MASK: 1000 */
{
  double xmg = (x[0]+x[2]          > DBL_EPSILON) ? x[0]+x[2]          : DBL_EPSILON;
  double xfe = (x[1]               > DBL_EPSILON) ? x[1]               : DBL_EPSILON;
  double xca = (1.0-x[0]-x[1]-x[2] > DBL_EPSILON) ? 1.0-x[0]-x[1]-x[2] : DBL_EPSILON;
  double xsi = (x[2]/2.0           > DBL_EPSILON) ? x[2]/2.0           : DBL_EPSILON;
  double xal = (1.0-x[2]           > DBL_EPSILON) ? 1.0-x[2]           : DBL_EPSILON;

  double g, dgdr[NR], fr[NA][NR];
  int i, j;

  fr[0][0] =      - x[0];  fr[0][1] =      - x[1];  fr[0][2] = -x[2];
  fr[1][0] = (1.0 - x[0]); fr[1][1] =      - x[1];  fr[1][2] = -x[2];
  fr[2][0] =      - x[0];  fr[2][1] = (1.0 - x[1]); fr[2][2] = -x[2];
  fr[3][0] =      - x[0];  fr[3][1] =      - x[1];  fr[3][2] = (1.0 - x[2]);

  g       = G;
  dgdr[0] = DGDR0;
  dgdr[1] = DGDR1;
  dgdr[2] = DGDR2;

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
    double d2gdr2[NR][NR], dfrdr[NA][NR], sum;
    int k;

    d2gdr2[0][0] = D2GDR0R0;     d2gdr2[0][1] = D2GDR0R1;     d2gdr2[0][2] = D2GDR0R2;
    d2gdr2[1][0] = d2gdr2[0][1]; d2gdr2[1][1] = D2GDR1R1;     d2gdr2[1][2] = D2GDR1R2;
    d2gdr2[2][0] = d2gdr2[0][2]; d2gdr2[2][1] = d2gdr2[1][2]; d2gdr2[2][2] = D2GDR2R2;

    for(i=0; i<NA; i++) {
       dfrdr[i][0] = DFR0DR0(i);
       dfrdr[i][1] = DFR1DR1(i);
       dfrdr[i][2] = DFR2DR2(i);
    }

    for (i=0; i<NA; i++) {
       for (k=0; k<NR; k++) {
          for (dx[i][k]=g, j=0; j<NR; j++) dx[i][k] += fr[i][j]*dgdr[j];
          dx[i][k] = exp(dx[i][k]/(R*t));
          sum = (1.0+dfrdr[i][k])*dgdr[k];
          for (j=0; j<NR; j++) sum += fr[i][j]*d2gdr2[j][k];
          dx[i][k] *= sum/(R*t);
       }
    }  
  }

  if (mask & FOURTH) {
    /* implement exclusion criteria on quantities for preclb routines  */
    double xCa = 3.0*(1.0 - x[0] - x[1] - x[2]);
    double xMg = 3.0*x[0] + 4.0*x[2];
    double xFe = 3.0*x[1];
    double xAl = 2.0*(1.0 - x[0] - x[1] - x[2]) + 2.0*x[0] + 2.0*x[1];
    double xSi = 3.0*(1.0 - x[0] - x[1] - x[2]) + 3.0*x[0] + 3.0*x[1] + 4.0*x[2];
    int temp[NA];
    temp[0] = (xCa > 0.15) && (xAl > 0.10);     /* Ca3Al2Si3O12 > 5% Ca on X, 5% Al on Y */
    temp[1] = (xMg > 0.15) && (xAl > 0.10);     /* Mg3Al2Si3O12 > 5% Mg on X, 5% Al on Y */
    temp[2] = (xFe > 0.15) && (xAl > 0.10);     /* Fe3Al2Si3O12 > 5% Fe on X, 5% Al on Y */
    temp[3] = (xMg > 0.15) && (xSi-3.0 > 0.05); /* Mg4Si4O12    > 5% Mg on X, 5% Si on Y */
    
    temp[3] = temp[3] && (p >= 50000.0); /* pressure must exceed 5 GPa */

    for (i=0; i<NA; i++) {
      if (!temp[i]) {
        if (mask & FIRST)  a[i]  = 0.0;
        if (mask & SECOND) mu[i] = 0.0;
        if (mask & THIRD)  for (j=0; j<NR; j++) dx[i][j] = 0.0;
      }
    }
  }

}

void 
gmixMaj(int mask, double t, double p, double *x, 
  double *gmix, /* Gibbs energy of mixing             BINARY MASK: 0001 */
  double *dx,   /* (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010 */
  double **dx2, /* (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100 */
  double ***dx3 /* (pointer to dx3[][]) d3(g)/d(x[])3 BINARY MASK: 1000 */
  )
{
  double xmg = (x[0]+x[2]          > DBL_EPSILON) ? x[0]+x[2]          : DBL_EPSILON;
  double xfe = (x[1]               > DBL_EPSILON) ? x[1]               : DBL_EPSILON;
  double xca = (1.0-x[0]-x[1]-x[2] > DBL_EPSILON) ? 1.0-x[0]-x[1]-x[2] : DBL_EPSILON;
  double xsi = (x[2]/2.0           > DBL_EPSILON) ? x[2]/2.0           : DBL_EPSILON;
  double xal = (1.0-x[2]           > DBL_EPSILON) ? 1.0-x[2]           : DBL_EPSILON;

  if (mask & FIRST) {
    *gmix = G;
  }
  
  if(mask & SECOND) {
    dx[0] = DGDR0;
    dx[1] = DGDR1;
    dx[2] = DGDR2;
  }

  if(mask & THIRD) {
    double d2gdr2[NR][NR];
    int i, j;

    d2gdr2[0][0] = D2GDR0R0;     d2gdr2[0][1] = D2GDR0R1;     d2gdr2[0][2] = D2GDR0R2;
    d2gdr2[1][0] = d2gdr2[0][1]; d2gdr2[1][1] = D2GDR1R1;     d2gdr2[1][2] = D2GDR1R2;
    d2gdr2[2][0] = d2gdr2[0][2]; d2gdr2[2][1] = d2gdr2[1][2]; d2gdr2[2][2] = D2GDR2R2;

    for (i=0; i<NR; i++) {
       for (j=0; j<NR; j++) dx2[i][j] = d2gdr2[i][j];
    }
  }

  if(mask & FOURTH) {
    double d3gdr3[NR][NR][NR];
    int i, j, k;

    d3gdr3[0][0][0] = D3GDR0R0R0;      d3gdr3[0][0][1] = D3GDR0R0R1;      d3gdr3[0][0][2] = D3GDR0R0R2;
    d3gdr3[0][1][0] = d3gdr3[0][0][1]; d3gdr3[0][1][1] = D3GDR0R1R1;      d3gdr3[0][1][2] = D3GDR0R1R2;
    d3gdr3[0][2][0] = d3gdr3[0][0][2]; d3gdr3[0][2][1] = d3gdr3[0][1][2]; d3gdr3[0][2][2] = D3GDR0R2R2;
    d3gdr3[1][0][0] = d3gdr3[0][0][1]; d3gdr3[1][0][1] = d3gdr3[0][1][1]; d3gdr3[1][0][2] = d3gdr3[0][1][2];
    d3gdr3[1][1][0] = d3gdr3[0][1][1]; d3gdr3[1][1][1] = D3GDR1R1R1;      d3gdr3[1][1][2] = D3GDR1R1R2;
    d3gdr3[1][2][0] = d3gdr3[0][1][2]; d3gdr3[1][2][1] = d3gdr3[1][1][2]; d3gdr3[1][2][2] = D3GDR1R2R2;
    d3gdr3[2][0][0] = d3gdr3[0][0][2]; d3gdr3[2][0][1] = d3gdr3[0][1][2]; d3gdr3[2][0][2] = d3gdr3[0][2][2];
    d3gdr3[2][1][0] = d3gdr3[0][1][2]; d3gdr3[2][1][1] = d3gdr3[1][1][2]; d3gdr3[2][1][2] = d3gdr3[1][2][2];
    d3gdr3[2][2][0] = d3gdr3[0][2][2]; d3gdr3[2][2][1] = d3gdr3[1][2][2]; d3gdr3[2][2][2] = D3GDR2R2R2;

    for (i=0; i<NR; i++) {
      for (j=0; j<NR; j++) {
	for (k=0; k<NR; k++) dx3[i][j][k] = d3gdr3[i][j][k];
      }
    }
  }
}

void 
hmixMaj(int mask, double t, double p, double *x, 
  double *hmix /* Enthalpy of mixing BINARY MASK: 1 */
  )
{
  double xmg = (x[0]+x[2]          > DBL_EPSILON) ? x[0]+x[2]          : DBL_EPSILON;
  double xfe = (x[1]               > DBL_EPSILON) ? x[1]               : DBL_EPSILON;
  double xca = (1.0-x[0]-x[1]-x[2] > DBL_EPSILON) ? 1.0-x[0]-x[1]-x[2] : DBL_EPSILON;
  double xsi = (x[2]/2.0           > DBL_EPSILON) ? x[2]/2.0           : DBL_EPSILON;
  double xal = (1.0-x[2]           > DBL_EPSILON) ? 1.0-x[2]           : DBL_EPSILON;

  *hmix = (G) + t*(S);
}

void 
smixMaj(int mask, double t, double p, double *x, 
  double *smix, /* Entropy of mixing                  BINARY MASK: 001 */
  double *dx,   /* (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 010 */
  double **dx2  /* (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 100 */
  )
{
  double xmg = (x[0]+x[2]          > DBL_EPSILON) ? x[0]+x[2]          : DBL_EPSILON;
  double xfe = (x[1]               > DBL_EPSILON) ? x[1]               : DBL_EPSILON;
  double xca = (1.0-x[0]-x[1]-x[2] > DBL_EPSILON) ? 1.0-x[0]-x[1]-x[2] : DBL_EPSILON;
  double xsi = (x[2]/2.0           > DBL_EPSILON) ? x[2]/2.0           : DBL_EPSILON;
  double xal = (1.0-x[2]           > DBL_EPSILON) ? 1.0-x[2]           : DBL_EPSILON;

  if (mask & FIRST) {
    *smix = S; 
  }
  
  if(mask & SECOND) {
    double d2gdrdt[NR];
    int i;

    d2gdrdt[0] = D2GDR0DT;
    d2gdrdt[1] = D2GDR1DT;
    d2gdrdt[2] = D2GDR2DT;

    for (i=0; i<NR; i++) dx[i] = - d2gdrdt[i];
  }

  if(mask & THIRD) {
    double d3gdr2dt[NR][NR];
    int i, j;

    d3gdr2dt[0][0] = D3GDR0R0DT;     d3gdr2dt[0][1] = D3GDR0R1DT;     d3gdr2dt[0][2] = D3GDR0R2DT;
    d3gdr2dt[1][0] = d3gdr2dt[0][1]; d3gdr2dt[1][1] = D3GDR1R1DT;     d3gdr2dt[1][2] = D3GDR1R2DT;
    d3gdr2dt[2][0] = d3gdr2dt[0][2]; d3gdr2dt[2][1] = d3gdr2dt[1][2]; d3gdr2dt[2][2] = D3GDR2R2DT;

    for (i=0; i<NR; i++) {
       for (j=0; j<NR; j++) dx2[i][j] = - d3gdr2dt[i][j];
    }
  }

}

void 
cpmixMaj(int mask, double t, double p, double *x, 
  double *cpmix, /* Heat capacity of mixing               BINARY MASK: 001 */
  double *dt,    /* d(cp)/d(t)                            BINARY MASK: 010 */
  double *dx     /* d(cp)/d(x[])d(t)                      BINARY MASK: 100 */
  )
{
  double d2gdt2;

  d2gdt2 = D2GDT2;

  if (mask & FIRST) {
    *cpmix = - t*d2gdt2;
  }

  if(mask & SECOND) {
    double d3gdt3   = D3GDT3;

    *dt = -t*d3gdt3 - d2gdt2;
  }

  if(mask & THIRD) {
    double d3gdrdt2[NR];
    int i;

    d3gdrdt2[0] = D3GDR0DT2;
    d3gdrdt2[1] = D3GDR1DT2;
    d3gdrdt2[2] = D3GDR2DT2;

    for (i=0; i<NR; i++) dx[i] = -t*d3gdrdt2[i];
  }

}

void 
vmixMaj(int mask, double t, double p, double *x, 
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

  if (mask & FIRST) {
    *vmix = DGDP;
  }

  if(mask & SECOND) {
    double d2gdrdp[NR];
    int i;

    d2gdrdp[0] = D2GDR0DP; 
    d2gdrdp[1] = D2GDR1DP;
    d2gdrdp[2] = D2GDR2DP;

    for (i=0; i<NR; i++) dx[i] = d2gdrdp[i]; 
  }

  if(mask & THIRD) {
    double d3gdr2dp[NR][NR];
    int i, j;

    d3gdr2dp[0][0] = D3GDR0R0DP;     d3gdr2dp[0][1] = D3GDR0R1DP;     d3gdr2dp[0][2] = D3GDR0R2DP;
    d3gdr2dp[1][0] = d3gdr2dp[0][1]; d3gdr2dp[1][1] = D3GDR1R1DP;     d3gdr2dp[1][2] = D3GDR1R2DP;
    d3gdr2dp[2][0] = d3gdr2dp[0][2]; d3gdr2dp[2][1] = d3gdr2dp[1][2]; d3gdr2dp[2][2] = D3GDR2R2DP;

    for (i=0; i<NR; i++) {
       for (j=0; j<NR; j++) dx2[i][j] = d3gdr2dp[i][j];
    }
  }

  if(mask & FOURTH) {
    *dt = D2GDTDP;
  }

  if(mask & FIFTH) {
    *dp = D2GDP2;
  }

  if(mask & SIXTH) {
    *dt2 = D3GDT2DP;
  }

  if(mask & SEVENTH) {
    *dtdp = D3GDTDP2;
  }

  if(mask & EIGHTH) {
    *dp2 = D3GDP3;
  }

  if(mask & NINTH) {
    double d3gdrdtdp[NR];
    int i;

    d3gdrdtdp[0] = D3GDR0DTDP;
    d3gdrdtdp[1] = D3GDR1DTDP;
    d3gdrdtdp[2] = D3GDR2DTDP;

    for (i=0; i<NR; i++) dxdt[i] = d3gdrdtdp[i];
  }

  if(mask & TENTH) {
    double d3gdrdp2[NR];
    int i;

    d3gdrdp2[0] = D3GDR0DP2;
    d3gdrdp2[1] = D3GDR1DP2;
    d3gdrdp2[2] = D3GDR2DP2;

    for (i=0; i<NR; i++) dxdp[i] = d3gdrdp2[i];
  }

}

/* end of file MAJORITE.C */
