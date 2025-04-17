#ifndef _Status_h
#define _Status_h

/*
MELTS Source Code: RCS $Log: status.h,v $
MELTS Source Code: RCS Revision 1.1  2007/12/22 22:43:30  ghiorso
MELTS Source Code: RCS Fixed error in BM integration in Gibbs.c
MELTS Source Code: RCS Updated param_struct_data.h file for AGU 2007 xMELTS parameters
MELTS Source Code: RCS Added support for status file production in MELTS-batch (XML)
MELTS Source Code: RCS
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Include file for intercommunication of status variables
**      (file: STATUS.H)
**
**--
*/

#define LIQUIDUS_SUCCESS        0
#define LIQUIDUS_MAX_T          1
#define LIQUIDUS_MIN_T          2
#define LIQUIDUS_TIME           3
#define LIQUIDUS_MULTIPLE       4

#define SILMIN_SUCCESS          5
#define SILMIN_QUAD_MAX         6
#define SILMIN_LIN_ZERO         7
#define SILMIN_LIN_MAX          8
#define SILMIN_ADD_LIQUID_1     9
#define SILMIN_ADD_LIQUID_2    10
#define SILMIN_ADD_LIQUID_3    11
#define SILMIN_RANK            12
#define SILMIN_TIME            13

#define GENERIC_INTERNAL_ERROR 14

#ifdef PHMELTS_ADJUSTMENTS
#define LIQUIDUS_SILMIN_ERROR  15
#endif

typedef struct _meltsStatus {
  int status;
} MeltsStatus;
extern MeltsStatus meltsStatus;

#endif /* _Status_h */
