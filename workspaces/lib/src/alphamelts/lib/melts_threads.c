const char *melts_threads_ver(void) { return "$Id: melts_threads.c,v 1.1 2006/10/20 00:59:22 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: melts_threads.c,v $
MELTS Source Code: RCS Revision 1.1  2006/10/20 00:59:22  ghiorso
MELTS Source Code: RCS (1) Made initial modifications for thread safe code.
MELTS Source Code: RCS (2) Added support for XML I/O in batch mode
MELTS Source Code: RCS (3) Added support for Melts-batch listener for eventual integration into VIGMCS
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 4.1  1999/11/13 22:09:00  ghiorso
MELTS Source Code: RCS Master Server Version (MELTS/Calc/pMELTS).
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 4.0  1999/06/18 17:25:43  ghiorso
MELTS Source Code: RCS Java MELTS v 1.1.0 Initial Check in
MELTS Source Code: RCS
*/

/*
**++
**  FACILITY:  Thread Safe Silicate Melts Crystallization Package
**
**  MODULE DESCRIPTION:
**
**  (file: MELTS_THREADS.C)
**--
*/

#include "silmin.h"

static MTHREAD_ONCE_T  meltsThreadControl = MTHREAD_ONCE_INIT;

static MTHREAD_KEY_T liquidFusKey;
static MTHREAD_KEY_T liquidCurKey;
static MTHREAD_KEY_T solidsCurKey;
static MTHREAD_KEY_T oxygenCurKey;
static MTHREAD_KEY_T silminStateKey;
static MTHREAD_KEY_T constraintsKey;

static void meltsThreadInit(void) {
  MTHREAD_KEY_CREATE(&liquidFusKey,   free);
  MTHREAD_KEY_CREATE(&liquidCurKey,   free);
  MTHREAD_KEY_CREATE(&solidsCurKey,   free);
  MTHREAD_KEY_CREATE(&oxygenCurKey,   free);
  MTHREAD_KEY_CREATE(&silminStateKey, destroySilminStateStructure);
  MTHREAD_KEY_CREATE(&constraintsKey, destroyConstraintsStructure);
}

ThermoData *getLiquidFus() {
  ThermoData *liquidFusPt;
  MTHREAD_ONCE(&meltsThreadControl, meltsThreadInit);
  
  liquidFusPt = (ThermoData *) MTHREAD_GETSPECIFIC(liquidFusKey);   
  if (liquidFusPt == NULL) {
    liquidFusPt = (ThermoData *) calloc((size_t) nlc, sizeof(ThermoData));
    MTHREAD_SETSPECIFIC(liquidFusKey, (void *) liquidFusPt);
  }
  return liquidFusPt; 
}

ThermoData *getLiquidCur() {
  ThermoData *liquidCurPt;
  MTHREAD_ONCE(&meltsThreadControl, meltsThreadInit);
  
  liquidCurPt = (ThermoData *) MTHREAD_GETSPECIFIC(liquidCurKey);   
  if (liquidCurPt == NULL) {
    liquidCurPt = (ThermoData *) calloc((size_t) nlc, sizeof(ThermoData));
    MTHREAD_SETSPECIFIC(liquidCurKey, (void *) liquidCurPt);
  }
  return liquidCurPt; 
}

ThermoData *getSolidsCur() {
  ThermoData *solidsCurPt;
  MTHREAD_ONCE(&meltsThreadControl, meltsThreadInit);
  
  solidsCurPt = (ThermoData *) MTHREAD_GETSPECIFIC(solidsCurKey);   
  if (solidsCurPt == NULL) {
    solidsCurPt = (ThermoData *) calloc((size_t) npc, sizeof(ThermoData));
    MTHREAD_SETSPECIFIC(solidsCurKey, (void *) solidsCurPt);
  }
  return solidsCurPt; 
}

ThermoData *getOxygenCur() {
  ThermoData *oxygenCurPt;
  MTHREAD_ONCE(&meltsThreadControl, meltsThreadInit);
  
  oxygenCurPt = (ThermoData *) MTHREAD_GETSPECIFIC(oxygenCurKey);   
  if (oxygenCurPt == NULL) {
    oxygenCurPt = (ThermoData *) calloc((size_t) 1, sizeof(ThermoData));
    MTHREAD_SETSPECIFIC(oxygenCurKey, (void *) oxygenCurPt);
  }
  return oxygenCurPt; 
}

SilminState *getSilminState() {
  SilminState *silminStatePt;
  MTHREAD_ONCE(&meltsThreadControl, meltsThreadInit);
  
  silminStatePt = (SilminState *) MTHREAD_GETSPECIFIC(silminStateKey);   
  if (silminStatePt == NULL) {
    silminStatePt = allocSilminStatePointer();
    MTHREAD_SETSPECIFIC(silminStateKey, (void *) silminStatePt);
  }
  return silminStatePt; 
}

void setSilminState(SilminState *silminState) {
  SilminState *silminStatePt;
  MTHREAD_ONCE(&meltsThreadControl, meltsThreadInit);
  
  silminStatePt = (SilminState *) MTHREAD_GETSPECIFIC(silminStateKey);   
  if (silminStatePt != NULL) destroySilminStateStructure((void *) silminStatePt);
  silminStatePt = silminState;
  MTHREAD_SETSPECIFIC(silminStateKey, (void *) silminStatePt);
}

Constraints *getConstraints() {
  Constraints *constraintsPt;
  MTHREAD_ONCE(&meltsThreadControl, meltsThreadInit);
  
  constraintsPt = (Constraints *) MTHREAD_GETSPECIFIC(constraintsKey);   
  if (constraintsPt == NULL) {
    constraintsPt = allocConstraintsPointer();
    MTHREAD_SETSPECIFIC(constraintsKey, (void *) constraintsPt);
  }
  return constraintsPt; 
}

void setConstraints(Constraints *constraints) {
  Constraints *constraintsPt;
  MTHREAD_ONCE(&meltsThreadControl, meltsThreadInit);
  
  constraintsPt = (Constraints *) MTHREAD_GETSPECIFIC(constraintsKey);   
  if (constraintsPt != NULL) destroyConstraintsStructure((void *) constraintsPt);
  constraintsPt = constraints;
  MTHREAD_SETSPECIFIC(constraintsKey, (void *) constraintsPt);
}
