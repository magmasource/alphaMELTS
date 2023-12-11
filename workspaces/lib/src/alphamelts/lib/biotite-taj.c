#include "silmin.h"  /* Structure definitions foor SILMIN package */

#define FILE_TIME_STAMP_BUFFER_SIZE 30

static char fileTimeStamp[FILE_TIME_STAMP_BUFFER_SIZE];

const char *biotiteTaj_ver(void) { 
  strncpy(fileTimeStamp, __BASE_FILE__, FILE_TIME_STAMP_BUFFER_SIZE);
  strncat(fileTimeStamp, __TIMESTAMP__, FILE_TIME_STAMP_BUFFER_SIZE - strlen(fileTimeStamp) - 1);
  return fileTimeStamp; 
}

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))
#define QUARTIC(x) ((x)*(x)*(x)*(x))

/*
 *=============================================================================
 * Biotite solution parameters (Binary regular solution theory):
 * Tajcmanova, Connolly, Cesare
 * J Metamorphic Geology, 27, 153-165
 */
#define WPHLANN  12000.0  /* joules */
#define WPHLEAST 10000.0  /* joules */
#define WPHLOBI   4000.0  /* joules */
#define WANNEAST  3000.0  /* joules */
#define WANNOBI   8000.0  /* joules */
#define WEASTOBI  7000.0  /* joules */

#define MU0FBI  0.0
#define MU0TBI  0.0
#define MU0EAST 0.0
#define MU0ANN  0.0
#define MU0PHL  0.0
#define MU0OBI  -6800.0

/*
 * Global (to this file): activity definitions and component transforms
 *    The function conTaj defines the conversion from m[i], to r[j]
 */
#define NR         4
#define NS         1
#define NA         5

#define FR0(i)     (i == 1) ? 1.0 - r[0] : - r[0]
#define FR1(i)     (i == 2) ? 1.0 - r[1] : - r[1]
#define FR2(i)     (i == 3) ? 1.0 - r[2] : - r[2]
#define FR3(i)     (i == 4) ? 1.0 - r[3] : - r[3]

#define DFR0DR0(i) - 1.0
#define DFR1DR1(i) - 1.0
#define DFR2DR2(i) - 1.0
#define DFR3DR3(i) - 1.0
#define DGS0DS0(i) - 1.0

#define MAX_ITER 200    /* Maximum number of iterations allowed in order */

/*
 * Global (to this file): derivative definitions
 */
#define R          8.3143
#define S          -R*(XM1Fe3*log(XM1Fe3) + XM1Al*log(XM1Al) + XM1Mg*log(XM1Mg) + XM1Fe2*log(XM1Fe2) \
                       + 2.0*XM2Ti*log(XM2Ti) + 2.0*XM2Mg*log(XM2Mg) + 2.0*XM2Fe2*log(XM2Fe2) \
		       + 2.0*XT1Al*log(XT1Al) + 2.0*XT1Si*log(XT1Si) + 2.0*XAOH*log(XAOH) + 2.0*XAO*log(XAO))
#define H          (1.0-r[0]-r[1]-r[2]-r[3])*MU0FBI + r[0]*MU0TBI + r[1]*MU0EAST \
                   + (2.0 - 3.0*r[0] + 4.0*r[2] - 2.0*r[3] - 2.0*s[0])*MU0ANN/6.0 \
		   + (2.0 - 3.0*r[0] - 2.0*r[2] + r[3] - 2.0*s[0])*MU0PHL/3.0 \
		   + (3.0*r[0]/2.0 + r[2] + r[3] - 1.0 + s[0])*MU0OBI \
		   + WPHLANN*(2.0 - 3.0*r[0] - 2.0*r[2] + r[3] - 2.0*s[0])*(2.0 - 3.0*r[0] + 4.0*r[2] - 2.0*r[3] - 2.0*s[0])/3.0/6.0 \
		   + WPHLEAST*(2.0 - 3.0*r[0] - 2.0*r[2] + r[3] - 2.0*s[0])*r[1]/3.0 \
		   + WPHLOBI*(2.0 - 3.0*r[0] - 2.0*r[2] + r[3] - 2.0*s[0])*(3.0*r[0]/2.0 + r[2] + r[3] - 1.0 + s[0])/3.0 \
		   + WANNEAST*(2.0 - 3.0*r[0] + 4.0*r[2] - 2.0*r[3] - 2.0*s[0])*r[1]/6.0 \
		   + WANNOBI*(2.0 - 3.0*r[0] + 4.0*r[2] - 2.0*r[3] - 2.0*s[0])*(3.0*r[0]/2.0 + r[2] + r[3] - 1.0 + s[0])/6.0 \
		   + WEASTOBI*r[1]*(3.0*r[0]/2.0 + r[2] + r[3] - 1.0 + s[0])
#define V          0.0
#define G          H - t*(S) + (p-1.0)*(V)

#define DGDR0      R*t*(-log(XM1Fe3) + log(XM1Fe2) + log(XM2Ti) - log(XM2Fe2) - log(XT1Al) + log(XT1Si) - 2.0*log(XAOH) + 2.0*log(XAO)) \
		   - MU0FBI + MU0TBI - MU0ANN/2.0 - MU0PHL + 3.0*MU0OBI/2.0 \
		   - 3.0*WPHLANN*(2.0 - 3.0*r[0] + 4.0*r[2] - 2.0*r[3] - 2.0*s[0])/3.0/6.0 \
		   - 3.0*WPHLANN*(2.0 - 3.0*r[0] - 2.0*r[2] + r[3] - 2.0*s[0])/3.0/6.0 \
		   - 3.0*WPHLEAST*r[1]/3.0 \
		   - 3.0*WPHLOBI*(3.0*r[0]/2.0 + r[2] + r[3] - 1.0 + s[0])/3.0 \
		   + 3.0*WPHLOBI*(2.0 - 3.0*r[0] - 2.0*r[2] + r[3] - 2.0*s[0])/3.0/2.0 \
		   - 3.0*WANNEAST*r[1]/6.0 \
		   - 3.0*WANNOBI*(3.0*r[0]/2.0 + r[2] + r[3] - 1.0 + s[0])/6.0 \
		   + 3.0*WANNOBI*(2.0 - 3.0*r[0] + 4.0*r[2] - 2.0*r[3] - 2.0*s[0])/6.0/2.0 \
		   + 3.0*WEASTOBI*r[1]/2.0		   
#define DGDR1      R*t*(-log(XM1Fe3) + log(XM1Al)) \
		   - MU0FBI + MU0EAST \
		   + WPHLEAST*(2.0 - 3.0*r[0] - 2.0*r[2] + r[3] - 2.0*s[0])/3.0 \
		   + WANNEAST*(2.0 - 3.0*r[0] + 4.0*r[2] - 2.0*r[3] - 2.0*s[0])/6.0 \
		   + WEASTOBI*(3.0*r[0]/2.0 + r[2] + r[3] - 1.0 + s[0])		   
#define DGDR2      R*t*(-log(XM1Fe3) - 2.0*log(XM1Mg)/3.0 +5.0*log(XM1Fe2)/3.0 - 4.0*log(XM2Mg)/3.0 + 4.0*log(XM2Fe2)/3.0 - log(XT1Al) + log(XT1Si) ) \
		   - MU0FBI + 4.0*MU0ANN/6.0 - 2.0*MU0PHL/3.0 + MU0OBI \
		   - 2.0*WPHLANN*(2.0 - 3.0*r[0] + 4.0*r[2] - 2.0*r[3] - 2.0*s[0])/3.0/6.0 \
		   + 4.0*WPHLANN*(2.0 - 3.0*r[0] - 2.0*r[2] + r[3] - 2.0*s[0])/3.0/6.0 \
		   - 2.0*WPHLEAST*r[1]/3.0 \
		   - 2.0*WPHLOBI*(3.0*r[0]/2.0 + r[2] + r[3] - 1.0 + s[0])/3.0 \
		   + WPHLOBI*(2.0 - 3.0*r[0] - 2.0*r[2] + r[3] - 2.0*s[0])/3.0 \
		   + 4.0*WANNEAST*r[1]/6.0 \
		   + 4.0*WANNOBI*(3.0*r[0]/2.0 + r[2] + r[3] - 1.0 + s[0])/6.0 \
		   + WANNOBI*(2.0 - 3.0*r[0] + 4.0*r[2] - 2.0*r[3] - 2.0*s[0])/6.0 \
		   + WEASTOBI*r[1]		   
#define DGDR3      R*t*(-log(XM1Fe3) + log(XM1Mg)/3.0 + 2.0*log(XM1Fe2)/3.0 + 2.0*log(XM2Mg)/3.0 - 2.0*log(XM2Fe2)/3.0 - log(XT1Al) + log(XT1Si)) \
		   - MU0FBI - 2.0*MU0ANN/6.0 + MU0PHL/3.0 + MU0OBI \
		   + WPHLANN*(2.0 - 3.0*r[0] + 4.0*r[2] - 2.0*r[3] - 2.0*s[0])/3.0/6.0 \
		   - 2.0*WPHLANN*(2.0 - 3.0*r[0] - 2.0*r[2] + r[3] - 2.0*s[0])/3.0/6.0 \
		   + WPHLEAST*r[1]/3.0 \
		   + WPHLOBI*(3.0*r[0]/2.0 + r[2] + r[3] - 1.0 + s[0])/3.0 \
		   + WPHLOBI*(2.0 - 3.0*r[0] - 2.0*r[2] + r[3] - 2.0*s[0])/3.0 \
		   - 2.0*WANNEAST*r[1]/6.0 \
		   - 2.0*WANNOBI*(3.0*r[0]/2.0 + r[2] + r[3] - 1.0 + s[0])/6.0 \
		   + WANNOBI*(2.0 - 3.0*r[0] + 4.0*r[2] - 2.0*r[3] - 2.0*s[0])/6.0 \
		   + WEASTOBI*r[1]   
#define DGDS0      R*t*(- 2.0*log(XM1Mg)/3.0 + 2.0*log(XM1Fe2)/3.0 + 2.0*log(XM2Mg)/3.0 - 2.0*log(XM2Fe2)/3.0) \
		   - 2.0*MU0ANN/6.0 - 2.0*MU0PHL/3.0 + MU0OBI \
		   - 2.0*WPHLANN*(2.0 - 3.0*r[0] + 4.0*r[2] - 2.0*r[3] - 2.0*s[0])/3.0/6.0 \
		   - 2.0*WPHLANN*(2.0 - 3.0*r[0] - 2.0*r[2] + r[3] - 2.0*s[0])/3.0/6.0 \
		   - 2.0*WPHLEAST*r[1]/3.0 \
		   - 2.0*WPHLOBI*(3.0*r[0]/2.0 + r[2] + r[3] - 1.0 + s[0])/3.0 \
		   + WPHLOBI*(2.0 - 3.0*r[0] - 2.0*r[2] + r[3] - 2.0*s[0])/3.0 \
		   - 2.0*WANNEAST*r[1]/6.0 \
		   - 2.0*WANNOBI*(3.0*r[0]/2.0 + r[2] + r[3] - 1.0 + s[0])/6.0 \
		   + WANNOBI*(2.0 - 3.0*r[0] + 4.0*r[2] - 2.0*r[3] - 2.0*s[0])/6.0 \
		   + WEASTOBI*r[1]
#define DGDT       -(S)
#define DGDP       (V)

#define D2GDR0R0   R*t*(1.0/XM1Fe3 + 1.0/XM1Fe2 + 1.0/XM2Ti/2.0 + 1.0/XM2Fe2/2.0 + 1.0/XT1Al/2.0 + 1.0/XT1Si/2.0 + 2.0/XAOH + 2.0/XAO) \
		   + 9.0*WPHLANN/3.0/6.0 + 9.0*WPHLANN/3.0/6.0 - 9.0*WPHLOBI/3.0/2.0 - 9.0*WPHLOBI/3.0/2.0 - 9.0*WANNOBI/6.0/2.0 - 9.0*WANNOBI/6.0/2.0   
#define D2GDR0R1   R*t*(1.0/XM1Fe3) - 3.0*WPHLEAST/3.0 - 3.0*WANNEAST/6.0 + 3.0*WEASTOBI/2.0
		   
#define D2GDR0R2   R*t*(1.0/XM1Fe3 + 5.0/XM1Fe2/3.0 - 2.0/XM2Fe2/3.0 + 1.0/XT1Al/2.0 + 1.0/XT1Si/2.0) \
		   - 12.0*WPHLANN/3.0/6.0 + 6.0*WPHLANN/3.0/6.0 - 3.0*WPHLOBI/3.0 - 6.0*WPHLOBI/3.0/2.0 - 3.0*WANNOBI/6.0 + 12.0*WANNOBI/6.0/2.0 		   		   
#define D2GDR0R3   R*t*(1.0/XM1Fe3 + 2.0/XM1Fe2/3.0 + 1.0/XM2Fe2/3.0 + 1.0/XT1Al/2.0 + 1.0/XT1Si/2.0) \
		   + 6.0*WPHLANN/3.0/6.0 - 3.0*WPHLANN/3.0/6.0 - 3.0*WPHLOBI/3.0 + 3.0*WPHLOBI/3.0/2.0 - 3.0*WANNOBI/6.0 - 6.0*WANNOBI/6.0/2.0
#define D2GDR0S0   R*t*(2.0/XM1Fe2/3.0 + 1.0/XM2Fe2/3.0) \
		   + 6.0*WPHLANN/3.0/6.0 + 6.0*WPHLANN/3.0/6.0 - 3.0*WPHLOBI/3.0 - 6.0*WPHLOBI/3.0/2.0 - 3.0*WANNOBI/6.0 - 6.0*WANNOBI/6.0/2.0
#define D2GDR0DT   R*(-log(XM1Fe3) + log(XM1Fe2) + log(XM2Ti) - log(XM2Fe2) - log(XT1Al) + log(XT1Si) - 2.0*log(XAOH) + 2.0*log(XAO))
#define D2GDR0DP   0.0

#define D2GDR1R1   R*t*(1.0/XM1Fe3 + 1.0/XM1Al)
#define D2GDR1R2   R*t*(1.0/XM1Fe3) - 2.0*WPHLEAST/3.0 + 4.0*WANNEAST/6.0 + WEASTOBI
#define D2GDR1R3   R*t*(1.0/XM1Fe3) + WPHLEAST/3.0 - 2.0*WANNEAST/6.0 + WEASTOBI
#define D2GDR1S0   - 2.0*WPHLEAST/3.0 - 2.0*WANNEAST/6.0 + WEASTOBI
#define D2GDR1DT   R*(-log(XM1Fe3) + log(XM1Al))
#define D2GDR1DP   0.0

#define D2GDR2R2   R*t*(1.0/XM1Fe3 + 4.0/XM1Mg/9.0 + 25.0/XM1Fe2/9.0 + 8.0/XM2Mg/9.0 + 8.0/XM2Fe2/9.0 + 1.0/XT1Al/2.0 + 1.0/XT1Si/2.0 ) \
		   - 8.0*WPHLANN/3.0/6.0 - 8.0*WPHLANN/3.0/6.0 - 2.0*WPHLOBI/3.0 - 2.0*WPHLOBI/3.0 + 4.0*WANNOBI/6.0 + 4.0*WANNOBI/6.0		   
#define D2GDR2R3   R*t*(1.0/XM1Fe3 - 2.0/XM1Mg/9.0 + 10.0/XM1Fe2/9.0 - 4.0/XM2Mg/9.0 - 4.0/XM2Fe2/9.0 + 1.0/XT1Al/2.0 + 1.0/XT1Si/2.0 ) \
		   + 4.0*WPHLANN/3.0/6.0 + 4.0*WPHLANN/3.0/6.0 - 2.0*WPHLOBI/3.0 + WPHLOBI/3.0 + 4.0*WANNOBI/6.0 - 2.0*WANNOBI/6.0		   
#define D2GDR2S0   R*t*(4.0/XM1Mg/9.0 + 10.0/XM1Fe2/9.0 - 4.0/XM2Mg/9.0 - 4.0/XM2Fe2/9.0) \
		   + 4.0*WPHLANN/3.0/6.0 - 8.0*WPHLANN/3.0/6.0 - 2.0*WPHLOBI/3.0 - 2.0*WPHLOBI/3.0 + 4.0*WANNOBI/6.0 - 2.0*WANNOBI/6.0		   
#define D2GDR2DT   R*(-log(XM1Fe3) - 2.0*log(XM1Mg)/3.0 +5.0*log(XM1Fe2)/3.0 - 4.0*log(XM2Mg)/3.0 + 4.0*log(XM2Fe2)/3.0 - log(XT1Al) + log(XT1Si) ) 
#define D2GDR2DP   0.0
#define D2GDR3R3   R*t*(1.0/XM1Fe3 + 1.0/XM1Mg/9.0 + 4.0/XM1Fe2/9.0 + 2.0/XM2Mg/9.0 + 2.0/XM2Fe2/9.0 + 1.0/XT1Al/2.0 + 1.0/XT1Si/2.0) \
		   - 2.0*WPHLANN/3.0/6.0 - 2.0*WPHLANN/3.0/6.0 + WPHLOBI/3.0 + WPHLOBI/3.0 - 2.0*WANNOBI/6.0 - 2.0*WANNOBI/6.0
		   
#define D2GDR3S0   R*t*(-2.0/XM1Mg/9.0 + 4.0/XM1Fe2/9.0 + 2.0/XM2Mg/9.0 + 2.0/XM2Fe2/9.0) \
		   - 2.0*WPHLANN/3.0/6.0 + 4.0*WPHLANN/3.0/6.0 + WPHLOBI/3.0 - 2.0*WPHLOBI/3.0 - 2.0*WANNOBI/6.0 - 2.0*WANNOBI/6.0 		   
#define D2GDR3DT   R*(-log(XM1Fe3) + log(XM1Mg)/3.0 + 2.0*log(XM1Fe2)/3.0 + 2.0*log(XM2Mg)/3.0 - 2.0*log(XM2Fe2)/3.0 - log(XT1Al) + log(XT1Si))
#define D2GDR3DP   0.0

#define D2GDS0S0   R*t*(4.0/XM1Mg/9.0 + 4.0/XM1Fe2/9.0 + 2.0/XM2Mg/9.0 + 2.0/XM2Fe2/9.0) \
		   + 4.0*WPHLANN/3.0/6.0 + 4.0*WPHLANN/3.0/6.0 - 2.0*WPHLOBI/3.0 - 2.0*WPHLOBI/3.0 - 2.0*WANNOBI/6.0 - 2.0*WANNOBI/6.0
#define D2GDS0DT   R*(- 2.0*log(XM1Mg)/3.0 + 2.0*log(XM1Fe2)/3.0 + 2.0*log(XM2Mg)/3.0 - 2.0*log(XM2Fe2)/3.0)
#define D2GDS0DP   0.0

#define D2GDT2     0.0
#define D2GDTDP    0.0
#define D2GDP2     0.0

#define D3GDR0R0R0 -R*t*(-1.0/XM1Fe3/XM1Fe3 + 1.0/XM1Fe2/XM1Fe2 + 1.0/XM2Ti/XM2Ti/4.0 - 1.0/XM2Fe2/XM2Fe2/4.0 - 1.0/XT1Al/XT1Al/4.0 \
                         + 1.0/XT1Si/XT1Si/4.0 - 2.0/XAOH/XAOH + 2.0/XAO/XAO)
#define D3GDR0R0R1 -R*t*(-1.0/XM1Fe3/XM1Fe3)
#define D3GDR0R0R2 -R*t*(-1.0/XM1Fe3/XM1Fe3 + 5.0/XM1Fe2/XM1Fe2/3.0 + 2.0/XM2Fe2/XM2Fe2/6.0 - 1.0/XT1Al/XT1Al/4.0 + 1.0/XT1Si/XT1Si/4.0)
#define D3GDR0R0R3 -R*t*(-1.0/XM1Fe3/XM1Fe3 + 2.0/XM1Fe2/XM1Fe2/3.0 - 1.0/XM2Fe2/XM2Fe2/6.0 - 1.0/XT1Al/XT1Al/4.0 + 1.0/XT1Si/XT1Si/4.0)
#define D3GDR0R0S0 -R*t*(2.0/XM1Fe2/XM1Fe2/3.0 - 1.0/XM2Fe2/XM2Fe2/6.0)
#define D3GDR0R0DT R*(1.0/XM1Fe3 + 1.0/XM1Fe2 + 1.0/XM2Ti/2.0 + 1.0/XM2Fe2/2.0 + 1.0/XT1Al/2.0 + 1.0/XT1Si/2.0 + 2.0/XAOH + 2.0/XAO)
#define D3GDR0R0DP 0.0

#define D3GDR0R1R1 -R*t*(-1.0/XM1Fe3/XM1Fe3)
#define D3GDR0R1R2 -R*t*(-1.0/XM1Fe3/XM1Fe3)
#define D3GDR0R1R3 -R*t*(-1.0/XM1Fe3/XM1Fe3)
#define D3GDR0R1S0 0.0
#define D3GDR0R1DT R*(1.0/XM1Fe3)
#define D3GDR0R1DP 0.0

#define D3GDR0R2R2 -R*t*(-1.0/XM1Fe3/XM1Fe3 + 25.0/XM1Fe2/XM1Fe2/9.0 - 4.0/XM2Fe2/XM2Fe2/9.0 - 1.0/XT1Al/XT1Al/4.0 + 1.0/XT1Si/XT1Si/4.0)
#define D3GDR0R2R3 -R*t*(-1.0/XM1Fe3/XM1Fe3 + 10.0/XM1Fe2/XM1Fe2/9.0 + 2.0/XM2Fe2/XM2Fe2/9.0 - 1.0/XT1Al/XT1Al/4.0 + 1.0/XT1Si/XT1Si/4.0)
#define D3GDR0R2S0 -R*t*(10.0/XM1Fe2/XM1Fe2/9.0 + 2.0/XM2Fe2/XM2Fe2/9.0)
#define D3GDR0R2DT R*(1.0/XM1Fe3 + 5.0/XM1Fe2/3.0 - 2.0/XM2Fe2/3.0 + 1.0/XT1Al/2.0 + 1.0/XT1Si/2.0)
#define D3GDR0R2DP 0.0

#define D3GDR0R3R3 -R*t*(-1.0/XM1Fe3/XM1Fe3 + 4.0/XM1Fe2/XM1Fe2/9.0 - 1.0/XM2Fe2/XM2Fe2/9.0 - 1.0/XT1Al/XT1Al/4.0 + 1.0/XT1Si/XT1Si/4.0)
#define D3GDR0R3S0 -R*t*(4.0/XM1Fe2/XM1Fe2/9.0 - 1.0/XM2Fe2/XM2Fe2/9.0)
#define D3GDR0R3DT R*(1.0/XM1Fe3 + 2.0/XM1Fe2/3.0 + 1.0/XM2Fe2/3.0 + 1.0/XT1Al/2.0 + 1.0/XT1Si/2.0)
#define D3GDR0R3DP 0.0

#define D3GDR0S0S0 -R*t*(4.0/XM1Fe2/XM1Fe2/9.0 - 1.0/XM2Fe2/XM2Fe2/9.0)
#define D3GDR0S0DT R*(2.0/XM1Fe2/3.0 + 1.0/XM2Fe2/3.0)
#define D3GDR0S0DP 0.0

#define D3GDR0DT2  0.0
#define D3GDR0DTDP 0.0
#define D3GDR0DP2  0.0

#define D3GDR1R1R1 -R*t*(-1.0/XM1Fe3/XM1Fe3 + 1.0/XM1Al/XM1Al)
#define D3GDR1R1R2 -R*t*(-1.0/XM1Fe3/XM1Fe3)
#define D3GDR1R1R3 -R*t*(-1.0/XM1Fe3/XM1Fe3)
#define D3GDR1R1S0 0.0
#define D3GDR1R1DT R*(1.0/XM1Fe3 + 1.0/XM1Al)
#define D3GDR1R1DP 0.0

#define D3GDR1R2R2 -R*t*(-1.0/XM1Fe3/XM1Fe3)
#define D3GDR1R2R3 -R*t*(-1.0/XM1Fe3/XM1Fe3)
#define D3GDR1R2S0 0.0
#define D3GDR1R2DT R*(1.0/XM1Fe3)
#define D3GDR1R2DP 0.0

#define D3GDR1R3R3 -R*t*(-1.0/XM1Fe3/XM1Fe3)
#define D3GDR1R3S0 0.0
#define D3GDR1R3DT R*(1.0/XM1Fe3)
#define D3GDR1R3DP 0.0

#define D3GDR1S0S0 0.0
#define D3GDR1S0DT 0.0
#define D3GDR1S0DP 0.0

#define D3GDR1DT2  0.0
#define D3GDR1DTDP 0.0
#define D3GDR1DP2  0.0

#define D3GDR2R2R2 -R*t*(-1.0/XM1Fe3/XM1Fe3 - 8.0/XM1Mg/XM1Mg/27.0 + 125.0/XM1Fe2/XM1Fe2/27.0 - 16.0/XM2Mg/XM2Mg/27.0 + 16.0/XM2Fe2/XM2Fe2/27.0 \
                         - 1.0/XT1Al/XT1Al/4.0 + 1.0/XT1Si/XT1Si/4.0 )
#define D3GDR2R2R3 -R*t*(-1.0/XM1Fe3/XM1Fe3 + 4.0/XM1Mg/XM1Mg/27.0 + 50.0/XM1Fe2/XM1Fe2/27.0 + 8.0/XM2Mg/XM2Mg/27.0 - 8.0/XM2Fe2/XM2Fe2/27.0 \
                         - 1.0/XT1Al/XT1Al/4.0 + 1.0/XT1Si/XT1Si/4.0 )
#define D3GDR2R2S0 -R*t*(-8.0/XM1Mg/XM1Mg/27.0 + 50.0/XM1Fe2/XM1Fe2/27.0 + 8.0/XM2Mg/XM2Mg/27.0 - 8.0/XM2Fe2/XM2Fe2/27.0)
#define D3GDR2R2DT R*(1.0/XM1Fe3 + 4.0/XM1Mg/9.0 + 25.0/XM1Fe2/9.0 + 8.0/XM2Mg/9.0 + 8.0/XM2Fe2/9.0 + 1.0/XT1Al/2.0 + 1.0/XT1Si/2.0 )
#define D3GDR2R2DP 0.0

#define D3GDR2R3R3 -R*t*(-1.0/XM1Fe3/XM1Fe3 - 2.0/XM1Mg/XM1Mg/27.0 + 20.0/XM1Fe2/XM1Fe2/27.0 - 4.0/XM2Mg/XM2Mg/27.0 + 4.0/XM2Fe2/XM2Fe2/27.0 \
                         - 1.0/XT1Al/XT1Al/4.0 + 1.0/XT1Si/XT1Si/4.0 )
#define D3GDR2R3S0 -R*t*(4.0/XM1Mg/XM1Mg/27.0 + 20.0/XM1Fe2/XM1Fe2/27.0 - 4.0/XM2Mg/XM2Mg/27.0 + 4.0/XM2Fe2/XM2Fe2/27.0)
#define D3GDR2R3DT R*(1.0/XM1Fe3 - 2.0/XM1Mg/9.0 + 10.0/XM1Fe2/9.0 - 4.0/XM2Mg/9.0 - 4.0/XM2Fe2/9.0 + 1.0/XT1Al/2.0 + 1.0/XT1Si/2.0 )
#define D3GDR2R3DP 0.0

#define D3GDR2S0S0 -R*t*(-8.0/XM1Mg/XM1Mg/27.0 + 20.0/XM1Fe2/XM1Fe2/27.0 - 4.0/XM2Mg/XM2Mg/27.0 + 4.0/XM2Fe2/XM2Fe2/27.0)
#define D3GDR2S0DT R*(4.0/XM1Mg/9.0 + 10.0/XM1Fe2/9.0 - 4.0/XM2Mg/9.0 - 4.0/XM2Fe2/9.0)
#define D3GDR2S0DP 0.0

#define D3GDR2DT2  0.0
#define D3GDR2DTDP 0.0
#define D3GDR2DP2  0.0

#define D3GDR3R3R3 -R*t*(-1.0/XM1Fe3/XM1Fe3 + 1.0/XM1Mg/XM1Mg/27.0 + 8.0/XM1Fe2/XM1Fe2/27.0 + 2.0/XM2Mg/XM2Mg/27.0 - 2.0/XM2Fe2/XM2Fe2/27.0 \
                         - 1.0/XT1Al/XT1Al/4.0 + 1.0/XT1Si/XT1Si/4.0)
#define D3GDR3R3S0 -R*t*(-2.0/XM1Mg/XM1Mg/27.0 + 8.0/XM1Fe2/XM1Fe2/27.0 + 2.0/XM2Mg/XM2Mg/27.0 - 2.0/XM2Fe2/XM2Fe2/27.0)
#define D3GDR3R3DT R*(1.0/XM1Fe3 + 1.0/XM1Mg/9.0 + 4.0/XM1Fe2/9.0 + 2.0/XM2Mg/9.0 + 2.0/XM2Fe2/9.0 + 1.0/XT1Al/2.0 + 1.0/XT1Si/2.0)
#define D3GDR3R3DP 0.0

#define D3GDR3S0S0 -R*t*(4.0/XM1Mg/XM1Mg/27.0 + 8.0/XM1Fe2/XM1Fe2/27.0 + 2.0/XM2Mg/XM2Mg/27.0 - 2.0/XM2Fe2/XM2Fe2/27.0)
#define D3GDR3S0DT R*(-2.0/XM1Mg/9.0 + 4.0/XM1Fe2/9.0 + 2.0/XM2Mg/9.0 + 2.0/XM2Fe2/9.0)
#define D3GDR3S0DP 0.0

#define D3GDR3DT2  0.0
#define D3GDR3DTDP 0.0
#define D3GDR3DP2  0.0

#define D3GDS0S0S0 -R*t*(-8.0/XM1Mg/XM1Mg/27.0 + 8.0/XM1Fe2/XM1Fe2/27.0 + 2.0/XM2Mg/XM2Mg/27.0 - 2.0/XM2Fe2/XM2Fe2/27.0)
#define D3GDS0S0DT R*(4.0/XM1Mg/9.0 + 4.0/XM1Fe2/9.0 + 2.0/XM2Mg/9.0 + 2.0/XM2Fe2/9.0)
#define D3GDS0S0DP 0.0

#define D3GDS0DT2  0.0
#define D3GDS0DTDP 0.0
#define D3GDS0DP2  0.0

#define D3GDT3     0.0
#define D3GDT2DP   0.0
#define D3GDTDP2   0.0
#define D3GDP3     0.0


/*
 *=============================================================================
 * Macros for automatic array initialization
 */
#define fillD2GDR2 \
 d2gdr2[0][0] = D2GDR0R0;     d2gdr2[0][1] = D2GDR0R1;     \
 d2gdr2[0][2] = D2GDR0R2;     d2gdr2[0][3] = D2GDR0R3;     \
 d2gdr2[1][0] = d2gdr2[0][1]; d2gdr2[1][1] = D2GDR1R1;     \
 d2gdr2[1][2] = D2GDR1R2;     d2gdr2[1][3] = D2GDR1R3;     \
 d2gdr2[2][0] = d2gdr2[0][2]; d2gdr2[2][1] = d2gdr2[1][2]; \
 d2gdr2[2][2] = D2GDR2R2;     d2gdr2[2][3] = D2GDR2R3;     \
 d2gdr2[3][0] = d2gdr2[0][3]; d2gdr2[3][1] = d2gdr2[1][3]; \
 d2gdr2[3][2] = d2gdr2[2][3]; d2gdr2[3][3] = D2GDR3R3;

#define fillD2GDRDS \
 d2gdrds[0][0] = D2GDR0S0; \
 d2gdrds[1][0] = D2GDR1S0; \
 d2gdrds[2][0] = D2GDR2S0; \
 d2gdrds[3][0] = D2GDR3S0;

#define fillD2GDRDT \
 d2gdrdt[0] = D2GDR0DT; d2gdrdt[1] = D2GDR1DT; d2gdrdt[2] = D2GDR2DT; \
 d2gdrdt[3] = D2GDR3DT;

#define fillD2GDRDP \
 d2gdrdp[0] = D2GDR0DP; d2gdrdp[1] = D2GDR1DP; d2gdrdp[2] = D2GDR2DP; \
 d2gdrdp[3] = D2GDR3DP;

#define fillD2GDS2 \
 d2gds2[0][0] = D2GDS0S0;

#define fillD2GDSDT \
 d2gdsdt[0] = D2GDS0DT;

#define fillD2GDSDP \
 d2gdsdp[0] = D2GDS0DP;

#define fillD3GDR3 \
 d3gdr3[0][0][0] = D3GDR0R0R0;		d3gdr3[0][0][1] = D3GDR0R0R1; \
 d3gdr3[0][0][2] = D3GDR0R0R2;		d3gdr3[0][0][3] = D3GDR0R0R3; \
 d3gdr3[0][1][0] = d3gdr3[0][0][1];	d3gdr3[0][1][1] = D3GDR0R1R1; \
 d3gdr3[0][1][2] = D3GDR0R1R2;		d3gdr3[0][1][3] = D3GDR0R1R3; \
 d3gdr3[0][2][0] = d3gdr3[0][0][2];	d3gdr3[0][2][1] = d3gdr3[0][1][2]; \
 d3gdr3[0][2][2] = D3GDR0R2R2;		d3gdr3[0][2][3] = D3GDR0R2R3; \
 d3gdr3[0][3][0] = d3gdr3[0][0][3];	d3gdr3[0][3][1] = d3gdr3[0][1][3]; \
 d3gdr3[0][3][2] = d3gdr3[0][2][3];	d3gdr3[0][3][3] = D3GDR0R3R3; \
 d3gdr3[1][0][0] = d3gdr3[0][0][1];	d3gdr3[1][0][1] = d3gdr3[0][1][1]; \
 d3gdr3[1][0][2] = d3gdr3[0][1][2];	d3gdr3[1][0][3] = d3gdr3[0][1][3]; \
 d3gdr3[1][1][0] = d3gdr3[0][1][1];	d3gdr3[1][1][1] = D3GDR1R1R1; \
 d3gdr3[1][1][2] = D3GDR1R1R2;		d3gdr3[1][1][3] = D3GDR1R1R3; \
 d3gdr3[1][2][0] = d3gdr3[0][1][2];	d3gdr3[1][2][1] = d3gdr3[1][1][2]; \
 d3gdr3[1][2][2] = D3GDR1R2R2;		d3gdr3[1][2][3] = D3GDR1R2R3; \
 d3gdr3[1][3][0] = d3gdr3[0][1][3];	d3gdr3[1][3][1] = d3gdr3[1][1][3]; \
 d3gdr3[1][3][2] = d3gdr3[1][2][3];	d3gdr3[1][3][3] = D3GDR1R3R3; \
 d3gdr3[2][0][0] = d3gdr3[0][0][2];	d3gdr3[2][0][1] = d3gdr3[0][1][2]; \
 d3gdr3[2][0][2] = d3gdr3[0][2][2];	d3gdr3[2][0][3] = d3gdr3[0][2][3]; \
 d3gdr3[2][1][0] = d3gdr3[0][1][2];	d3gdr3[2][1][1] = d3gdr3[1][1][2]; \
 d3gdr3[2][1][2] = d3gdr3[1][2][2];	d3gdr3[2][1][3] = d3gdr3[1][2][3]; \
 d3gdr3[2][2][0] = d3gdr3[0][2][2];	d3gdr3[2][2][1] = d3gdr3[1][2][2]; \
 d3gdr3[2][2][2] = D3GDR2R2R2;		d3gdr3[2][2][3] = D3GDR2R2R3; \
 d3gdr3[2][3][0] = d3gdr3[0][2][3];	d3gdr3[2][3][1] = d3gdr3[1][2][3]; \
 d3gdr3[2][3][2] = d3gdr3[2][2][3];	d3gdr3[2][3][3] = D3GDR2R3R3; \
 d3gdr3[3][0][0] = d3gdr3[0][0][3];	d3gdr3[3][0][1] = d3gdr3[0][1][3]; \
 d3gdr3[3][0][2] = d3gdr3[0][2][3];	d3gdr3[3][0][3] = d3gdr3[0][3][3]; \
 d3gdr3[3][1][0] = d3gdr3[0][1][3];	d3gdr3[3][1][1] = d3gdr3[1][1][3]; \
 d3gdr3[3][1][2] = d3gdr3[1][2][3];	d3gdr3[3][1][3] = d3gdr3[1][3][3]; \
 d3gdr3[3][2][0] = d3gdr3[0][2][3];	d3gdr3[3][2][1] = d3gdr3[1][2][3]; \
 d3gdr3[3][2][2] = d3gdr3[2][2][3];	d3gdr3[3][2][3] = d3gdr3[2][3][3]; \
 d3gdr3[3][3][0] = d3gdr3[0][3][3];	d3gdr3[3][3][1] = d3gdr3[1][3][3]; \
 d3gdr3[3][3][2] = d3gdr3[2][3][3];	d3gdr3[3][3][3] = D3GDR3R3R3;

#define fillD3GDR2DS \
 d3gdr2ds[0][0][0] = D3GDR0R0S0; \
 d3gdr2ds[0][1][0] = D3GDR0R1S0; \
 d3gdr2ds[0][2][0] = D3GDR0R2S0; \
 d3gdr2ds[0][3][0] = D3GDR0R3S0; \
 d3gdr2ds[1][0][0] = d3gdr2ds[0][1][0]; \
 d3gdr2ds[1][1][0] = D3GDR1R1S0; \
 d3gdr2ds[1][2][0] = D3GDR1R2S0; \
 d3gdr2ds[1][3][0] = D3GDR1R3S0; \
 d3gdr2ds[2][0][0] = d3gdr2ds[0][2][0]; \
 d3gdr2ds[2][1][0] = d3gdr2ds[1][2][0]; \
 d3gdr2ds[2][2][0] = D3GDR2R2S0; \
 d3gdr2ds[2][3][0] = D3GDR2R3S0; \
 d3gdr2ds[3][0][0] = d3gdr2ds[0][3][0]; \
 d3gdr2ds[3][1][0] = d3gdr2ds[1][3][0]; \
 d3gdr2ds[3][2][0] = d3gdr2ds[2][3][0]; \
 d3gdr2ds[3][3][0] = D3GDR3R3S0;

#define fillD3GDR2DT \
 d3gdr2dt[0][0] = D3GDR0R0DT;     d3gdr2dt[0][1] = D3GDR0R1DT; \
 d3gdr2dt[0][2] = D3GDR0R2DT;     d3gdr2dt[0][3] = D3GDR0R3DT; \
 d3gdr2dt[1][0] = d3gdr2dt[0][1]; d3gdr2dt[1][1] = D3GDR1R1DT; \
 d3gdr2dt[1][2] = D3GDR1R2DT;     d3gdr2dt[1][3] = D3GDR1R3DT; \
 d3gdr2dt[2][0] = d3gdr2dt[0][2]; d3gdr2dt[2][1] = d3gdr2dt[1][2]; \
 d3gdr2dt[2][2] = D3GDR2R2DT;     d3gdr2dt[2][3] = D3GDR2R3DT; \
 d3gdr2dt[3][0] = d3gdr2dt[0][3]; d3gdr2dt[3][1] = d3gdr2dt[1][3]; \
 d3gdr2dt[3][2] = d3gdr2dt[2][3]; d3gdr2dt[3][3] = D3GDR3R3DT;

#define fillD3GDR2DP \
 d3gdr2dp[0][0] = D3GDR0R0DP;     d3gdr2dp[0][1] = D3GDR0R1DP; \
 d3gdr2dp[0][2] = D3GDR0R2DP;     d3gdr2dp[0][3] = D3GDR0R3DP; \
 d3gdr2dp[1][0] = d3gdr2dp[0][1]; d3gdr2dp[1][1] = D3GDR1R1DP; \
 d3gdr2dp[1][2] = D3GDR1R2DP;     d3gdr2dp[1][3] = D3GDR1R3DP; \
 d3gdr2dp[2][0] = d3gdr2dp[0][2]; d3gdr2dp[2][1] = d3gdr2dp[1][2]; \
 d3gdr2dp[2][2] = D3GDR2R2DP;     d3gdr2dp[2][3] = D3GDR2R3DP; \
 d3gdr2dp[3][0] = d3gdr2dp[0][3]; d3gdr2dp[3][1] = d3gdr2dp[1][3]; \
 d3gdr2dp[3][2] = d3gdr2dp[2][3]; d3gdr2dp[3][3] = D3GDR3R3DP;

#define fillD3GDRDS2 \
 d3gdrds2[0][0][0] = D3GDR0S0S0; \
 d3gdrds2[1][0][0] = D3GDR1S0S0; \
 d3gdrds2[2][0][0] = D3GDR2S0S0; \
 d3gdrds2[3][0][0] = D3GDR3S0S0;

#define fillD3GDRDSDT \
 d3gdrdsdt[0][0] = D3GDR0S0DT; \
 d3gdrdsdt[1][0] = D3GDR1S0DT; \
 d3gdrdsdt[2][0] = D3GDR2S0DT; \
 d3gdrdsdt[3][0] = D3GDR3S0DT;

#define fillD3GDRDSDP \
 d3gdrdsdp[0][0] = D3GDR0S0DP; \
 d3gdrdsdp[1][0] = D3GDR1S0DP; \
 d3gdrdsdp[2][0] = D3GDR2S0DP; \
 d3gdrdsdp[3][0] = D3GDR3S0DP;

#define fillD3GDS3 \
 d3gds3[0][0][0] = D3GDS0S0S0;

#define fillD3GDS2DT \
 d3gds2dt[0][0] = D3GDS0S0DT;

#define fillD3GDS2DP \
 d3gds2dp[0][0] = D3GDS0S0DP;

#define fillD3GDSDT2 \
 d3gdsdt2[0] = D3GDS0DT2;

#define fillD3GDSDTDP \
 d3gdsdtdp[0] = D3GDS0DTDP;

#define fillD3GDSDP2 \
 d3gdsdp2[0] = D3GDS0DP2;

#define fillD3GDRDT2 \
 d3gdrdt2[0] = D3GDR0DT2; d3gdrdt2[1] = D3GDR1DT2; \
 d3gdrdt2[2] = D3GDR2DT2; d3gdrdt2[3] = D3GDR3DT2;

#define fillD3GDRDTDP \
 d3gdrdtdp[0] = D3GDR0DTDP; d3gdrdtdp[1] = D3GDR1DTDP; \
 d3gdrdtdp[2] = D3GDR2DTDP; d3gdrdtdp[3] = D3GDR3DTDP;

#define fillD3GDRDP2 \
 d3gdrdp2[0] = D3GDR0DP2; d3gdrdp2[1] = D3GDR1DP2; \
 d3gdrdp2[2] = D3GDR2DP2; d3gdrdp2[3] = D3GDR3DP2;

static double XM1Fe3, XM1Al, XM1Mg, XM1Fe2, XM2Ti, XM2Mg, XM2Fe2, XT1Al, XT1Si, XAOH, XAO;

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
  static double tOld           = 9999.0;
  static double pOld           = 9999.0;
  static double rOld[NR]       = { -2.0, -2.0, -2.0, -2.0 };
  static double sOld[NS]       = { -2.0 };
  static double d2gds2[NS][NS] = { { 0.0 } };
  int i, j, iter = 0;

  if ( (t != tOld)       || (p != pOld) ||
       (r[0] != rOld[0]) || (r[1] != rOld[1]) || (r[2] != rOld[2]) ||
       (r[3] != rOld[3]) ) {
    double dgds[NS], sNew[NS], Qmin, Qmax, lambda=1.0;
    
                                          	      Qmin = 2.0*r[2]-r[3]-2.0;
    if (Qmin < (-r[2]+r[3]/2.0-1.0/2.0))  	      Qmin = -r[2]+r[3]/2.0-1.0/2.0;
    if (Qmin < (-3.0*r[0]/2.0+2.0*r[2]-r[3]-2.0))     Qmin = -3.0*r[0]/2.0+2.0*r[2]-r[3]-2.0;
    if (Qmin < (-3.0*r[0]/2.0-5.0*r[2]/2.0-r[3]+1.0)) Qmin = -3.0*r[0]/2.0-5.0*r[2]/2.0-r[3]+1.0;
    
                                                  	Qmax = 2.0*r[2] - r[3] + 1.0;
    if (Qmax > (1.0-r[2]+r[3]/2.0))               	Qmax = 1.0-r[2]+r[3]/2.0;
    if (Qmax > (1.0-3.0*r[0]/2.0+2.0*r[2]-r[3]))  	Qmax = 1.0-3.0*r[0]/2.0+2.0*r[2]-r[3];
    if (Qmax > (5.0/2.0-3.0*r[0]/2.-5.0*r[2]/2.0-r[3])) Qmax = 5.0/2.0-3.0*r[0]/2.-5.0*r[2]/2.0-r[3];

    for (i=0; i<NS; i++) sOld[i] = 2.0;
    dgds[0] = 0.0;
    sNew[0] = (Qmin+Qmax)/2.0;

    while ( ((ABS(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON)) && (iter < MAX_ITER)) {
      double s[NS], deltaS[NS];

      for (i=0; i<NS; i++) s[i] = sNew[i];

      XM1Fe3 = 1.0 - r[0] - r[1] -r[2] - r[3]; 
      XM1Al  = r[1]; 
      XM1Mg  = (2.0 - 2.0*r[2] + r[3] - 2.0*s[0])/3.0; 
      XM1Fe2 = (3.0*r[0] + 5.0*r[2] + 2.0*r[3] + 2.0*s[0] - 2.0)/3.0; 
      XM2Ti  = r[0]/2.0; 
      XM2Mg  = (2.0 - 2.0*r[2] + r[3] + s[0])/3.0;
      XM2Fe2 = (2.0 - 3.0*r[0] + 4.0*r[2] - 2.0*r[3] - 2.0*s[0])/6.0;
      XT1Al  = 1.0 - r[0]/2.0 - r[2]/2.0 -r[3]/2.0; 
      XT1Si  = r[0]/2.0 + r[2]/2.0 + r[3]/2.0;
      XAOH   = 1.0 - r[0]; 
      XAO    = r[0];

      if (XM1Fe3 <= 0.0) XM1Fe3 = DBL_EPSILON;
      if (XM1Al  <= 0.0) XM1Al  = DBL_EPSILON;
      if (XM1Mg  <= 0.0) XM1Mg  = DBL_EPSILON;
      if (XM1Fe2 <= 0.0) XM1Fe2 = DBL_EPSILON;
      if (XM2Ti  <= 0.0) XM2Ti  = DBL_EPSILON;
      if (XM2Mg  <= 0.0) XM2Mg  = DBL_EPSILON;
      if (XM2Fe2 <= 0.0) XM2Fe2 = DBL_EPSILON;
      if (XT1Al  <= 0.0) XT1Al  = DBL_EPSILON;
      if (XT1Si  <= 0.0) XT1Si  = DBL_EPSILON;
      if (XAOH   <= 0.0) XAOH   = DBL_EPSILON;
      if (XAO	 <= 0.0) XAO    = DBL_EPSILON;

      if (XM1Fe3 >= 1.0) XM1Fe3 = 1.0 - DBL_EPSILON;
      if (XM1Al  >= 1.0) XM1Al  = 1.0 - DBL_EPSILON;
      if (XM1Mg  >= 1.0) XM1Mg  = 1.0 - DBL_EPSILON;
      if (XM1Fe2 >= 1.0) XM1Fe2 = 1.0 - DBL_EPSILON;
      if (XM2Ti  >= 1.0) XM2Ti  = 1.0 - DBL_EPSILON;
      if (XM2Mg  >= 1.0) XM2Mg  = 1.0 - DBL_EPSILON;
      if (XM2Fe2 >= 1.0) XM2Fe2 = 1.0 - DBL_EPSILON;
      if (XT1Al  >= 1.0) XT1Al  = 1.0 - DBL_EPSILON;
      if (XT1Si  >= 1.0) XT1Si  = 1.0 - DBL_EPSILON;
      if (XAOH   >= 1.0) XAOH   = 1.0 - DBL_EPSILON;
      if (XAO	 >= 1.0) XAO    = 1.0 - DBL_EPSILON;

      dgds[0]      = DGDS0;
      d2gds2[0][0] = D2GDS0S0;

      for (i=0; i<NS; i++) sOld[i] = s[i];

      d2gds2[0][0] = 1.0/d2gds2[0][0];

      for (i=0; i<NS; i++) {
        for(j=0; j<NS; j++) s[i] += - d2gds2[i][j]*dgds[j];
        deltaS[i] = s[i] - sOld[i]; 
      }
      
      lambda = 1.0;
      while (lambda > DBL_EPSILON) {
        int OK = TRUE;
        for (i=0; i<NS; i++) s[i] = sOld[i] + lambda*deltaS[i];
	if ((2.0-2.0*r[2]+r[3]) > 0.0) { /* Mg */
          if(((2.0 - 2.0*r[2] + r[3] - 2.0*s[0])/3.0) < 100.0*DBL_EPSILON) OK = FALSE; 
          if(((2.0 - 2.0*r[2] + r[3] + s[0])/3.0)     < 100.0*DBL_EPSILON) OK = FALSE;
	}
	if (r[2] > 0.0) { /* Fe2+ */
	  if(((3.0*r[0] + 5.0*r[2] + 2.0*r[3] + 2.0*s[0] - 2.0)/3.0) < 100.0*DBL_EPSILON) OK = FALSE; 
          if(((2.0 - 3.0*r[0] + 4.0*r[2] - 2.0*r[3] - 2.0*s[0])/6.0) < 100.0*DBL_EPSILON) OK = FALSE;
	}
	if (OK) break;
	lambda /= 2.0;
      }

      for (i=0; i<NS; i++) sNew[i] = s[i];
      iter++;
    }
    tOld = t;
    pOld = p;
    for (i=0; i<NR; i++) rOld[i] = r[i];

    for (i=0; i<NS; i++) {
      if (dgds[i] > sqrt(DBL_EPSILON) && ABS(sOld[i]) > DBL_EPSILON) {
        printf("ERROR in biotiteTaj.c (function ORDER). Failed to converge!\n");
        if (iter >= MAX_ITER) printf("  Iteration limit (%4d) exceeded.\n", iter);
        printf("  r0    = %13.6g, r1    = %13.6g, r2    = %13.6g,  r3    = %13.6g\n", r[0], r[1], r[2], r[3]);
        printf("  s0    = %13.6g, lam   = %13.6g\n", sOld[0], lambda);
        printf("  dgds0 = %13.6g\n", dgds[0]);
        printf("  X Fe3+ M1: %13.6g  X Al M1: %13.6g  X Mg M1: %13.6g  X Fe2+ M1: %13.6g\n", XM1Fe3, XM1Al, XM1Mg, XM1Fe2); 
        printf("  X Ti   M2: %13.6g                   X Mg M2: %13.6g  X Fe2+ M2: %13.6g\n", XM2Ti, XM2Mg, XM2Fe2); 
        printf("  X Si   T1: %13.6g  X Al T1: %13.6g\n", XT1Si, XT1Al); 
        printf("  X OH   A:  %13.6g  X O  A:  %13.6g\n", XAOH, XAO); 
        break;
      }
    }

  }

  if (mask & FIRST  ) {   /* return s        */
    for (i=0; i<NS; i++) s[i] = sOld[i];
  }   
  if (mask & SECOND ) {   /* compute ds/dr:  */
    double d2gdrds[NR][NS];
    int k;                    

    fillD2GDRDS

    for (i=0; i<NS; i++) {
       for (j=0; j<NR; j++) {
          dr[i][j] = 0.0; 
          for (k=0; k<NS; k++) dr[i][j] += - d2gds2[i][k]*d2gdrds[j][k];
       }
    }
  }
  if (mask & THIRD  ) {   /* compute ds/dt:  */
    double d2gdsdt[NS];

    fillD2GDSDT

    for (i=0; i<NS; i++) {
       dt[i] = 0.0;
       for (j=0; j<NS; j++) dt[i] += - d2gds2[i][j]*d2gdsdt[j];
    }
  }
  if (mask & FOURTH ) {   /* compute ds/dp:  */
    double d2gdsdp[NS];

    fillD2GDSDP

    for (i=0; i<NS; i++) {
      dp[i] = 0.0; 
      for (j=0; j<NS; j++) dp[i] += - d2gds2[i][j]*d2gdsdp[j];
    }
  }
  if (mask & FIFTH  ) {   /* compute d2s/dr2 */
    double d2gdrds[NR][NS], d3gdr2ds[NR][NR][NS], d3gdrds2[NR][NS][NS],
      d3gds3[NS][NS][NS], dsdr[NS][NR], temp[NS];
    int k, l, m, n;                    

    fillD2GDRDS
    fillD3GDR2DS
    fillD3GDRDS2
    fillD3GDS3

    /* compute dsdr matrix */
    for (i=0; i<NS; i++) {
      for (j=0; j<NR; j++) {
        dsdr[i][j] = 0.0; 
        for (k=0; k<NS; k++) dsdr[i][j] += - d2gds2[i][k]*d2gdrds[j][k];
      }
    }

    /* compute dsdr2 cube */
    for (i=0; i<NS; i++) {
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
          dr2[i][j][k] = 0.0;
          for (l=0; l<NS; l++) dr2[i][j][k] += - d2gds2[i][l]*temp[l];
        }
      }
    }

  }
  if (mask & SIXTH  ) {   /* compute d2s/drt */
    double d2gdrds[NR][NS], d2gdsdt[NS], d3gdrds2[NR][NS][NS],
      d3gdrdsdt[NR][NS], d3gds3[NS][NS][NS], d3gds2dt[NS][NS], dsdr[NS][NR],
      dsdt[NS], temp[NS];
    int k, l, m;

    fillD2GDRDS
    fillD2GDSDT
    fillD3GDRDS2
    fillD3GDRDSDT
    fillD3GDS3
    fillD3GDS2DT

    /* compute dsdr matrix */
    for (i=0; i<NS; i++) {
      for (j=0; j<NR; j++) {
        dsdr[i][j] = 0.0; 
        for (k=0; k<NS; k++) dsdr[i][j] += - d2gds2[i][k]*d2gdrds[j][k];
      }
    }

    /* compute dsdt vector */
    for (i=0; i<NS; i++) {
      dsdt[i] = 0.0;
      for (j=0; j<NS; j++) dsdt[i] += - d2gds2[i][j]*d2gdsdt[j];
    }

    /* compute dsdrdt matrix */
    for (i=0; i<NS; i++) {
      for (j=0; j<NR; j++) {
        for (k=0; k<NS; k++) {
          temp[k] = d3gdrdsdt[j][k];
          for (l=0; l<NS; l++) {
             temp[k] += d3gdrds2[j][k][l]*dsdt[l] + d3gds2dt[k][l]*dsdr[l][j];
             for (m=0; m<NS; m++) temp[k] += d3gds3[k][l][m]*dsdr[l][j]*dsdt[m];
          }
        }
        drt[i][j] = 0.0;
        for (k=0; k<NS; k++) drt[i][j] += - d2gds2[i][k]*temp[k];
      }
    }

  }
  if (mask & SEVENTH) {   /* compute d2s/drp */
    double d2gdrds[NR][NS], d2gdsdp[NS], d3gdrds2[NR][NS][NS],
      d3gdrdsdp[NR][NS], d3gds3[NS][NS][NS], d3gds2dp[NS][NS], dsdr[NS][NR],
      dsdp[NS], temp[NS];
    int k, l, m;

    fillD2GDRDS
    fillD2GDSDP
    fillD3GDRDS2
    fillD3GDRDSDP
    fillD3GDS3
    fillD3GDS2DP

    /* compute dsdr matrix */
    for (i=0; i<NS; i++) {
      for (j=0; j<NR; j++) {
        dsdr[i][j] = 0.0; 
        for (k=0; k<NS; k++) dsdr[i][j] += - d2gds2[i][k]*d2gdrds[j][k];
      }
    }

    /* compute dsdp vector */
    for (i=0; i<NS; i++) {
      dsdp[i] = 0.0;
      for (j=0; j<NS; j++) dsdp[i] += - d2gds2[i][j]*d2gdsdp[j];
    }

    /* compute dsdrdp matrix */
    for (i=0; i<NS; i++) {
      for (j=0; j<NR; j++) {
        for (k=0; k<NS; k++) {
          temp[k] = d3gdrdsdp[j][k];
          for (l=0; l<NS; l++) {
             temp[k] += d3gdrds2[j][k][l]*dsdp[l] + d3gds2dp[k][l]*dsdr[l][j];
             for (m=0; m<NS; m++) temp[k] += d3gds3[k][l][m]*dsdr[l][j]*dsdp[m];
          }
        }
        drp[i][j] = 0.0;
        for (k=0; k<NS; k++) drp[i][j] += - d2gds2[i][k]*temp[k];
      }
    }

  }
  if (mask & EIGHTH ) {   /* compute d2s/dt2 */
    double d2gdsdt[NS], d3gds3[NS][NS][NS], d3gds2dt[NS][NS], d3gdsdt2[NS],
      dsdt[NS], temp[NS];
    int k, l;

    fillD2GDSDT
    fillD3GDS3
    fillD3GDS2DT
    fillD3GDSDT2

    /* compute dsdt vector */
    for (i=0; i<NS; i++) {
      dsdt[i] = 0.0;
      for (j=0; j<NS; j++) dsdt[i] += - d2gds2[i][j]*d2gdsdt[j];
    }

    /* compute dsdt2 vector */
    for (i=0; i<NS; i++) {
      for (j=0; j<NS; j++) { 
        temp[j] = d3gdsdt2[j];
        for (k=0; k<NS; k++) {
          temp[j] +=  2.0*d3gds2dt[j][k]*dsdt[k];
          for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdt[l];
        }
      }
      dt2[i] = 0.0;
      for (j=0; j<NS; j++) dt2[i] += - d2gds2[i][j]*temp[j];
    } 

  }
  if (mask & NINTH  ) {   /* compute d2s/dtp */
    double d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS][NS][NS], d3gds2dt[NS][NS],
      d3gds2dp[NS][NS], d3gdsdtdp[NS], dsdt[NS], dsdp[NS], temp[NS];
    int k, l;

    fillD2GDSDT
    fillD2GDSDP
    fillD3GDS3
    fillD3GDS2DT
    fillD3GDS2DP
    fillD3GDSDTDP

    /* compute dsdt vector */
    for (i=0; i<NS; i++) {
      dsdt[i] = 0.0;
      for (j=0; j<NS; j++) dsdt[i] += - d2gds2[i][j]*d2gdsdt[j];
    }

    /* compute dsdp vector */
    for (i=0; i<NS; i++) {
      dsdp[i] = 0.0;
      for (j=0; j<NS; j++) dsdp[i] += - d2gds2[i][j]*d2gdsdp[j];
    }

    /* compute dsdtp vector */
    for (i=0; i<NS; i++) {
      for (j=0; j<NS; j++) {
        temp[j] = d3gdsdtdp[j];
        for (k=0; k<NS; k++) {
          temp[j] += d3gds2dt[j][k]*dsdp[k] + d3gds2dp[j][k]*dsdt[k];
          for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdp[l];
        }
      }
      dtp[i] = 0.0;
      for (j=0; j<NS; j++) dtp[i] += - d2gds2[i][j]*temp[j];
    }

  }
  if (mask & TENTH  ) {   /* compute d2s/dp2 */
    double d2gdsdp[NS], d3gds3[NS][NS][NS], d3gds2dp[NS][NS], d3gdsdp2[NS],
      dsdp[NS], temp[NS];
    int k, l;

    fillD2GDSDP
    fillD3GDS3
    fillD3GDS2DP
    fillD3GDSDP2

    /* compute dsdp vector */
    for (i=0; i<NS; i++) {
      dsdp[i] = 0.0;
      for (j=0; j<NS; j++) dsdp[i] += - d2gds2[i][j]*d2gdsdp[j];
    }

    /* compute dsdp2 vector */
    for (i=0; i<NS; i++) {
      for (j=0; j<NS; j++) { 
        temp[j] = d3gdsdp2[j];
        for (k=0; k<NS; k++) {
          temp[j] +=  2.0*d3gds2dp[j][k]*dsdp[k];
          for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dsdp[k]*dsdp[l];
        }
      }
      dp2[i] = 0.0;
      for (j=0; j<NS; j++) dp2[i] += - d2gds2[i][j]*temp[j];
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
testTaj(int mask, double t, double p,
  int na,          /* Expected number of endmember components                 */
  int nr,          /* Expected number of independent compositional variables  */
  char **names,    /* array of strings of names of endmember components       */
  char **formulas, /* array of strings of formulas of endmember components    */
  double *r,       /* array of indepependent compos variables, check bounds   */
  double *m)       /* array of moles of endmember components, check bounds    */
{
  const char *phase = "biotiteTaj.c";
  const char *NAMES[NA]    = { "fbiTaj", "tbiTaj", "eastTaj", "annTaj", "phlTaj" };
  const char *FORMULAS[NA] = { "KFeMg2Si2Al2O10(OH)2", "KTiMg2Si3AlO10(O)2", "KAlMg2Si2Al2O10(OH)2",
                               "KFe3Si3AlO10(OH)2",    "KMg3Si3AlO10(OH)2" };
  int result = TRUE, i;

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
    double totalFe3  = 1.0 - r[0] - r[1] - r[2] - r[3];
    double totalFe2  = 3.0*r[2];
    double totalAlM1 = r[1]; 
    double totalTi   = r[0];
    double totalMg   = 2.0 - 2.0*r[2] + r[3];
    double totalOH   = 2.0 - 2.0*r[0];
    double totalO    = 2.0*r[0];
    double totalAlT1 = 2.0 - r[0] - r[2] - r[3];
    double totalSi   = r[0] + r[2] + r[3];
    
    result = result && (totalFe3  >= 0.0) && (totalFe3  <= 1.0);
    result = result && (totalFe2  >= 0.0) && (totalFe2  <= 3.0);
    result = result && (totalAlM1 >= 0.0) && (totalAlM1 <= 1.0);
    result = result && (totalTi   >= 0.0) && (totalTi   <= 1.0);
    result = result && (totalMg   >= 0.0) && (totalMg   <= 3.0);
    result = result && (totalOH   >= 0.0) && (totalOH   <= 2.0);
    result = result && (totalO    >= 0.0) && (totalO    <= 2.0);
    result = result && (totalAlT1 >= 1.0) && (totalAlT1 <= 2.0);
    result = result && (totalSi   >= 0.0) && (totalSi   <= 1.0);
  }
  /* Check bounds on moles of endmember components */
  if (mask & SIXTH) {
    double sum, r[NR];
    for (i=0, sum=0.0; i<NA; i++) sum += m[i];
    for (i=0; i<NR; i++) r[i] = (sum != 0.0) ? m[i+1]/sum : 0.0;
    {
      double totalFe3  = 1.0 - r[0] - r[1] - r[2] - r[3];
      double totalFe2  = 3.0*r[2];
      double totalAlM1 = r[1]; 
      double totalTi   = r[0];
      double totalMg   = 2.0 - 2.0*r[2] + r[3];
      double totalOH   = 2.0 - 2.0*r[0];
      double totalO    = 2.0*r[0];
      double totalAlT1 = 2.0 - r[0] - r[2] - r[3];
      double totalSi   = r[0] + r[2] + r[3];
      
      result = result && (totalFe3  >= 0.0) && (totalFe3  <= 1.0);
      result = result && (totalFe2  >= 0.0) && (totalFe2  <= 3.0);
      result = result && (totalAlM1 >= 0.0) && (totalAlM1 <= 1.0);
      result = result && (totalTi   >= 0.0) && (totalTi   <= 1.0);
      result = result && (totalMg   >= 0.0) && (totalMg   <= 3.0);
      result = result && (totalOH   >= 0.0) && (totalOH   <= 2.0);
      result = result && (totalO    >= 0.0) && (totalO    <= 2.0);
      result = result && (totalAlT1 >= 1.0) && (totalAlT1 <= 2.0);
      result = result && (totalSi   >= 0.0) && (totalSi   <= 1.0);
    }
  }

  return result;
}

void
conTaj(int inpMask, int outMask, double t, double p,
  double *e,      /* comp of biotite in moles of elements                     */
  double *m,      /* comp of biotite in moles of endmember components         */
  double *r,      /* comp of biotite in terms of the independent comp var     */
  double *x,      /* comp of biotite in mole fractions of endmember comp      */
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
      endmember biotite components.
  (2) calculates from a vector of moles of endmember components, one or
      all of: r[], x[], dr[]/dm[], d2r[]/dm[]dm[], or d3r[]/dm[]dm[]dm[]
  (3) calculates from a vector of independent compositional variables
      mole fractions of endmember components and/or the Jacobian matrix
      dx[]/dr[]

  ----------------------------------------------------------------------------*/

  int i, j, k;

  if (inpMask == FIRST && outMask == SECOND) {
    static const int Hy =  1;
    static const int Li =  3;
    static const int O  =  8;
    static const int F  =  9;
    static const int Na = 11;
    static const int Mg = 12;
    static const int Al = 13;
    static const int Cl = 17;
    static const int K  = 19;
    static const int Ca = 20;
    static const int Ti = 22;
    static const int Mn = 25;
    static const int Fe = 26;
    static const int Rb = 37;
    static const int Sr = 38;
    static const int Cs = 55;
    static const int Ba = 56;
    
    double hydrogen, oxygen, cationCharge, anionCharge, fe3, fe2, silicon;
    
    if (e[Hy] != 0.0) {
      hydrogen = e[Hy];
      oxygen   = e[O];
    } else {
      hydrogen =  2.0*(e[Li]+e[Na]+e[K]+e[Ca]+e[Rb]+e[Sr]+e[Cs]+e[Ba]) - 2.0*e[Ti] - e[F] - e[Cl];
      oxygen   = 12.0*(e[Li]+e[Na]+e[K]+e[Ca]+e[Rb]+e[Sr]+e[Cs]+e[Ba]);
    }
    silicon = 2.0*oxygen/3.0 
            - (e[Li]+e[Na]+e[K]+e[Ca]+e[Rb]+e[Sr]+e[Cs]+e[Ba] + e[Mg]+e[Al]+e[Ti]+e[Mn]+e[Fe]);
    
    cationCharge = hydrogen + e[Li] + 2.0*e[Mg] + 3.0*e[Al] + 4.0*silicon + e[Na] + e[K] + 2.0*e[Ca] 
                 + 4.0*e[Ti] + 2.0*e[Mn] + e[Rb] + 2.0*e[Sr] + e[Cs] + 2.0*e[Ba];
    anionCharge  = 2.0*oxygen + e[F] + e[Cl];
    fe3 = anionCharge - cationCharge - 2.0*e[Fe];
    fe2 = e[Fe] - fe3;
                        				     			    /* Projection into the Fe-Mg binary */
    m[0] = fe3;         				     			    /* KFeMg2Si2Al2O10(OH)2 */
    m[1] = e[Ti];       				     			    /* KTiMg2Si3AlO10(O)2   */
    m[2] = (3.0*e[Al] - e[Mg] - 4.0*fe3 - e[Ti] - fe2)/7.0;  			    /* KAlMg2Si2Al2O10(OH)2 */
    m[3] = (fe2+e[Mn])/3.0;                                  			    /* KFe3Si3AlO10(OH)2    */        
    m[4] = (3.0*e[Mg] - 2.0*e[Al] - 2.0*fe3 - 4.0*e[Ti] + 2.0*(fe2+e[Mn])/3.0)/7.0; /* KMg3Si3AlO10(OH)2    */
    
  } else if (inpMask == SECOND) {
    double sum;

    if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH | EIGHTH))
      printf("Illegal call to conTaj with inpMask = %o and outMask = %o\n",
        inpMask, outMask);

    for (i=0, sum=0.0; i<NA; i++) sum += m[i];

    if (outMask & THIRD) {
      for (i=0; i<NR; i++) r[i] = (sum != 0.0) ? m[i+1]/sum : 0.0; 
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
        for (i=0; i<NR; i++) {
          for (j=0; j<NA; j++)
            dm[i][j] = (i+1 == j) ? (1.0-m[i+1]/sum)/sum : - m[i+1]/SQUARE(sum);
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
        for (i=0; i<NR; i++) {
          for (j=0; j<NA; j++)  {
            for (k=0; k<NA; k++) {
              d2m[i][j][k]  = 2.0*m[i+1]/CUBE(sum);
              d2m[i][j][k] -= (i+1 == j) ? 1.0/SQUARE(sum) : 0.0;
              d2m[i][j][k] -= (i+1 == k) ? 1.0/SQUARE(sum) : 0.0;
            }
          }
        }
      }
    }

    if (outMask & EIGHTH) {
      /* calculates the matrix d3r[i]/dm[j]dm[k]dm[l] using m[] as input	*/
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
          for (j=0; j<NA; j++)  {
            for (k=0; k<NA; k++)  {
	      for (l=0; l<NA; l++)  {
                d3m[i][j][k][l]  = -6.0*m[i+1]/QUARTIC(sum);
                d3m[i][j][k][l] += (i+1 == j) ? 2.0/CUBE(sum) : 0.0;
                d3m[i][j][k][l] += (i+1 == k) ? 2.0/CUBE(sum) : 0.0;
                d3m[i][j][k][l] += (i+1 == l) ? 2.0/CUBE(sum) : 0.0;
	      }
            }
          }
        }
      }

    }

  } else if (inpMask == THIRD) {

    if (outMask & ~(FOURTH | SEVENTH))
      printf("Illegal call to conTaj with inpMask = %o and outMask = %o\n",
        inpMask, outMask);

    if (outMask & FOURTH) {
      /* Converts a vector of independent compositional variables (r) 
         into a vector of mole fractions of endmember components (x).         */

      for (i=0, x[0]=1.0; i<NR; i++) { x[i+1] = r[i]; x[0] -= r[i]; }
    }

    if (outMask & SEVENTH) {
      /* computes the Jacobian matrix dr[i][j] = dx[i]/dr[j] */
      for (i=0; i<NR; i++) for (j=0; j<NR; j++) dr[i+1][j] = (i == j) ? 1.0 : 0.0;
      for (j=0; j<NR; j++) dr[0][j] = -1.0;
    }

  } else  {
    printf("Illegal call to conTaj with inpMask = %o and outMask = %o\n",
      inpMask, outMask);
  }

}

void
dispTaj(int mask, double t, double p, double *x,
  char **formula            /* Mineral formula for interface display MASK: 1 */
  )
{
  double *r = x;
  static char masterString[] = {
/*             111111111122222222223333333333444444444455555555556
     0123456789012345678901234567890123456789012345678901234567890 */
    "KFe''_.__Mg_.__Fe'''_.__Ti_.__Al_.__Si_.__O10(OH)_.__O_.__" };


  if (mask & FIRST) {
    char *string = strcpy((char *) malloc((size_t) (strlen(masterString)+1)*sizeof(char)), masterString);
    double totFe2, totMg, totFe3, totAl, totTi, totSi, totOH, totO;
    char n[5];
    int i;

    totFe2 = 3.0*r[2];
    totMg  = 2.0*(1.0-r[0]-r[1]-r[2]-r[3]) + 2.0*r[0] + 2.0*r[1] + 3.0*r[3];
    totFe3 = 1.0 - r[0] - r[1] - r[2] - r[3];
    totAl  = 2.0*(1.0-r[0]-r[1]-r[2]-r[3]) + r[0] + 3.0*r[1] + r[2] + r[3];
    totTi  = r[0];
    totSi  = 2.0 + r[0] + r[2] + r[3];
    totOH  = 2.0*(1.0-r[0]-r[1]-r[2]-r[3]) + 2.0*r[1] + 2.0*r[2] + 2.0*r[3];
    totO   = 2.0*r[0]; 

    (void) snprintf(n, 5, "%4.2f", totFe2); for (i=0; i<4; i++) string[ 5+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totMg);  for (i=0; i<4; i++) string[11+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totFe3); for (i=0; i<4; i++) string[20+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totTi);  for (i=0; i<4; i++) string[26+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totAl);  for (i=0; i<4; i++) string[32+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totSi);  for (i=0; i<4; i++) string[38+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totOH);  for (i=0; i<4; i++) string[49+i] = n[i];
    (void) snprintf(n, 5, "%4.2f", totO);   for (i=0; i<4; i++) string[54+i] = n[i];

    *formula = string;
  }
}

void 
actTaj(int mask, double t, double p, double *x, 
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
     fr[i][0] = FR0(i); 
     fr[i][1] = FR1(i); 
     fr[i][2] = FR2(i); 
     fr[i][3] = FR3(i); 
  }

  order(FIRST, t, p, r, 
        s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

  g       = G;
  dgdr[0] = DGDR0;
  dgdr[1] = DGDR1;
  dgdr[2] = DGDR2;
  dgdr[3] = DGDR3;

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
       dfrdr[i][0] = DFR0DR0(i); 
       dfrdr[i][1] = DFR1DR1(i); 
       dfrdr[i][2] = DFR2DR2(i); 
       dfrdr[i][3] = DFR3DR3(i); 
       dgsds[i][0] = DGS0DS0(i); 
    }
    gs[0][0] =    1.0     - s[0];
    gs[1][0] = - (1.0/2.0 + s[0]);
    gs[2][0] =    1.0     - s[0];
    gs[3][0] =            - s[0];
    gs[4][0] =            - s[0];

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
    /* implement exclusion criteria on quantities for preclb routines  */
    static const double exclusion[NA] = {
       0.05, 0.05, 0.05, 0.05, 0.05
    };
    double x[NA];

    x[0] = 1.0 - r[0] - r[1] - r[2] - r[3];  
    x[1] = r[0]; 
    x[2] = r[1]; 
    x[3] = r[2]; 
    x[4] = r[3]; 

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
gmixTaj(int mask, double t, double p, double *x, 
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
    *gmix = G;
  }
  
  if(mask & SECOND) {
    dx[0] = DGDR0;
    dx[1] = DGDR1;
    dx[2] = DGDR2;
    dx[3] = DGDR3;
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

  if (mask & FOURTH) {
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
hmixTaj(int mask, double t, double p, double *x, 
  double *hmix /* Enthalpy of mixing BINARY MASK: 1 */
  )
{
  double *r = x;
  double s[NS];
  
  order(FIRST, t, p, r, 
        s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  
  *hmix = (G) + t*(S);
}

void 
smixTaj(int mask, double t, double p, double *x, 
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
cpmixTaj(int mask, double t, double p, double *x, 
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
vmixTaj(int mask, double t, double p, double *x, 
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

/* end of file BIOTITE.C */
