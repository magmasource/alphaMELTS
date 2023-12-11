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

void setoutput();
char *fgetstring(char *s, int size, FILE *stream);
char *scanfilename(char *line);
void putMultipleDataToFile(char *fileName, char *fileName2, SilminState sBlock[], int z, int writeorappend);

double getDspLiquidProperties(SilminState *state, int nl, double *value);
double getDspSolidProperties(SilminState *state, int index, int ns, double *value);

#endif
