/* 
To compile the ENKI xMELTS version of the shared library for use with Matlab:

  make realclean Melts-dynamicPrivate "BATCH=-DBATCH_VERSION -DRHYOLITE_ADJUSTMENTS" DYNAMIC=true

To compile the ENKI xMELTS version of Test_dynamicLib but linked to the shared library, as above plus:

  export LD_LIBRARY_PATH=.

*/

void getMeltsVersionString(int *failure, char *errorString, int *nCharInString);
void addConsole(void);
void closeConsole(void);
int getCalculationMode(void);
int setCalculationMode(int mode);
void getMeltsOxideNames(int *failure, char *oxideNames, int *nCharInName, int *numberOxides);
void getMeltsOxideWeights(int *failure, double *oxideWeights, int *numberOxides);
void getMeltsPhaseNames(int *failure, char *phaseNames, int *nCharInName, int *numberPhases, int phaseIndices[]);
void getMeltsWeightsAndFormulas(int *failure, char *phaseName, double *endMemberWeights, char *formulaPtr,
				 int *nCharInName, int *numberEndMembers);
void driveMeltsProcess(int *failure, int *mode, double *pressure, double *bulkComposition,
		       double *entropy, double *temperature,
		       char *phaseNames, int *nCharInName, int *numberPhases, int *output, 
		       char *errorString, int *nCharInString, double *phaseProperties, int phaseIndices[]);
void setMeltsSystemProperties(int *failure, char *strings, int *nCharInString, int *numberStrings);
void getMeltsPhaseProperties(int *failure, char *phaseName, double *temperature, double *pressure,
			      double *bulkComposition, double *phaseProperties);
void getMeltsViscosity(int *failure, char *model, double *temperature, double *bulkComposition, double *viscosity);
void getMeltsMolarProperties(int *failure, char *phaseName, double *temperature, double *pressure,
			      double *bulkComposition, double *phaseProperties);
void getMeltsEndMemberProperties(int *failure, char *phaseName, double *temperature, 
				 double *pressure, double *bulkComposition,
				 char *endMemberFormulas, int *nCharInName, int *numberEndMembers, 
				 double *endMemberProperties);
void getMeltsOxideProperties(int *failure, char *phaseName, double *temperature, 
				 double *pressure, double *bulkComposition,
				 char *oxideNames, int *nCharInName, int *numberOxides, 
				 double *oxideProperties);
void getMeltsSaturationState(int *failure, double *pressure, double *bulkComposition, double *temperature, 
                 char *phasePtr, int *nCharInName, int *numberPhases, double *phaseProperties, int phaseIndices[]);
/*void getMeltsPhaseParameters(int *failure, char *phaseName, int type, int set, int *index1, int *index2,
  int nIn, double *phaseParam);
void setMeltsPhaseParameters(int *failure, char *phaseName, int type, int set, int *index1, int *index2,
int nIn, double *phaseParam);*/
