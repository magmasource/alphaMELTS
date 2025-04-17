#ifndef _alphaMelts_h
#define _alphaMelts_h

#include "adiabat.h"
#include "phmelts.h"

void splashScreen(double newestversion);
int assignGlobalStatics(int mode);
int getEnvironmentSettings(void);
void assignProblemStatics(void);

void menu_option0(void);
int menu_option1(char *filename);
int menu_option2(double dspTstart, double dspTstop, double dspTinc,
  double dspPstart, double dspPstop, double dspPinc);

int menu_option3(int equilibriumGuess);
int menu_option4(int equilibriumGuess, int saveAll, int iterMax);
int menu_option5(int fo2Path, double fo2Delta);

int menu_option7(int constraint, double newRef, double newStop, double newInc);
int menu_option8(int i, int j, int incSolids, int fracSolids, int minType, double newMin);
int menu_option9(int i, int j, int incSolids, int fracSolids, int minType, double newMin);

int menu_option13(char *filename);

#endif
