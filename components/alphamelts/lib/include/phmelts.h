/*
**  function prototypes for adiabatic path finding routines
*/

#ifndef _pHMELTS_H
#define _pHMELTS_H

#define THIS_VERSION 2.0301

void initializeLibrary();
int setCalculationMode(int mode);
SilminState *createSilminState(void);

SilminState *reAllocSilminStatePointer(SilminState *q, size_t zOld, size_t z);
void copyStateInfo(SilminState *target, SilminState *source);
void copyThermoData(ThermoData *target, ThermoData *source);
void multiplyThermoData(ThermoData *target, double factor);
void addThermoData(ThermoData *value, ThermoData temp);
void freeSilminStatePointer(SilminState *p);
void printThermoData(SilminState *state);

#endif
