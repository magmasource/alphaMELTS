/*
**  function prototypes for adiabatic path finding routines
*/

#ifndef _Adiabat_H
#define _Adiabat_H

/*
 *=============================================================================
 *      Global text tools (see subbasic_read_write.c):
 */
extern int *outsw;
extern int guessFlag;

void setoutput();
char *fgetstring(char *s, int size, FILE *stream);
char *scanfilename(char *line);
void putMultipleDataToFile(char *fileName, char *fileName2, SilminState sBlock[], int z, int writeorappend);

double getDspLiquidProperties(SilminState *state, int nl, double *value);
double getDspSolidProperties(SilminState *state, int index, int ns, double *value);

double fractionateSolids();
double extractMelt();

extern int modeFlag, ptpathFlag, guessFlag;

extern double Pmax, Pmin, Tmax, Tmin;
extern SilminState *startState, *oldState, *states;

int adiabat_0ph(int equilibriumGuess);
int adiabat_1ph(int equilibriumGuess, int saveAll, int iterMax, int noprint);

void startingSolution(void);
int adiabatFunc();

void printPhases();

#define ISENTROPIC 0
#define ISOTHERMAL 1
#define ISOBARIC 2
#define GEOTHERMAL 3
#define PTPATH 4
#define ISOCHORIC 5
#define ISENTHALPIC 6
#define PTGRID 7
#define PTFILE 8
#define PSEUDOSECTION 9

#endif
