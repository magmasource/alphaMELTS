const char *read_write_ver(void) { return "$Id: read_write.c,v 1.6 2008/03/06 17:51:23 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: read_write.c,v $
MELTS Source Code: RCS Revision 1.6  2008/03/06 17:51:23  ghiorso
MELTS Source Code: RCS New fluid fractionation mode and other enhancements.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.5  2007/06/08 17:25:43  ghiorso
MELTS Source Code: RCS Added code to allow regression of Ghiorso EOS parameters
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.4  2007/02/01 02:03:08  ghiorso
MELTS Source Code: RCS Minor updates.
MELTS Source Code: RCS Prior to modifying preclb.c code to read xml data file
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
MELTS Source Code: RCS Revision 1.4  2005/01/25 03:25:03  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.3  2004/12/04 19:10:36  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2004/01/02 19:21:49  cvsaccount
MELTS Source Code: RCS CTserver University of Chicago
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2003/04/28 20:44:46  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:58:12  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
 * Revision 3.10  1997/06/21  22:49:29  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.9  1997/05/03  20:23:08  ghiorso
 * *** empty log message ***
 *
 * Revision 3.8  1997/03/27  17:03:13  ghiorso
 * *** empty log message ***
 *
 * Revision 3.7  1996/09/24  20:33:23  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 3.6  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 3.5  1995/11/23  22:37:42  ghiorso
 * Final implementation of subsolidus fO2 buffering.
 *
 * Revision 3.4  1995/11/01  22:40:27  ghiorso
 * Implementation of subsolidus options after Asimow.
 * Additional implementation of nepheline solid solutions.
 *
 * Revision 3.3  1995/09/04  20:01:28  ghiorso
 * Update to allow display of bulk composition (in grams) in the text entry
 * fields of the main silmin display. Liquid composition is no longer
 * display here, and is available only through the popup selection.
 *
 * Revision 3.2  1995/09/04  17:36:28  ghiorso
 * Corrected read/write of MELTS input file to allow for four decimals
 * in bulk system composition.
 *
 * Revision 3.1  1995/08/18  19:11:12  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Perform input/output functions for crystallization mode of
**      the melts pacakge. (File: READ_WRITE.C)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  August    31, 1991 - Original Version
**                   September  1, 1991 - continue
**      V1.0-2  Mark S. Ghiorso  September 3, 1991
**          (1) Removed includedSolids[].status reference
**          (2) Changed access to tpValues through global macros
**          (3) Changed compositionValues[].textwg -> .name
**          (4) changed statusEntries reference to new style
**          (5) added magma and assimilant initialization
**          (6) added toggle button initialization and setting
**      V1.0-3  Mark S. Ghiorso  September 4, 1991
**          (1) Added reference to silminState global strcture
**          and silminInputData global structure
**          (2) Finished assimilant and magma mixing initialization
**          (3) Added initialization of status_adb and solid_adb entries
**          ->->-> Final version to implement all expected options
**          (4) Since the text widgets set the numeric values of
**          their string displays into the appropriate global
**          structures, this initialization has been removed.
**          (5) Toggle buttons are now initialized using the
**          DecWindows XmToggleButtonSetState.
**      V1.0-4  Mark S. Ghiorso  September 20, 1991
**          (1) Removed initialization of phases[0].* labels
**      V2.0-1  Mark S. Ghiorso  November 14, 1991
**          Conversion to OSF Motif V1.1.1/X11 Release 4
**      V2.0-2  Mark S. Ghiorso  November 23, 1991
**          (1) Converted the update method for the display of text in the
**          statusEntries[STATUS_ADB_INDEX_STATUS].name widget to
**          reflect that it is now a scrolled text widget. The
**          wprintf (widget printf) function is used for this purpose.
**      V2.0-3  Mark S. Ghiorso  November 27, 1991
**          (1) Added vframe, vheader, vlist include files inorder to
**          call VListRemoveAll lines on the new phases vlist widget
**      V2.0-4  Mark S. Ghiorso  November 28, 1991
**          (1) Explicitly initialize compositionValues[].value to zero
**          since the textField valueChanged callback is not invoked
**          for a NULL string
**      V2.0-5  Mark S. Ghiorso  December 3, 1991
**          (1) Same as 2.0-4 for tpValues[].value
**      V2.0-6  Mark S. Ghiorso  December 10, 1991
**          (1) Removed reference to magma_padb
**          (2) reorganized assimilantValues and assimilantUnits structures
**          (3) removed tpValues[].name references and direct setting
**          of text widget strings
**          (4) rewrote code to load assimilantValues structure
**      V2.0-7  Mark S. Ghiorso  January 6, 1992
**          (1) Added functions putInputDataToFile() and
**          putOutputDataToFile()
**          Mark S. Ghiorso  January 7, 1992
**          (1) Continued work on putInputDataToFile()
**          Mark S. Ghiorso  January 8, 1992
**          (1) Completed work on putOutputDatatoFile()
**      V2.0-8  Mark S. Ghiorso  January 15, 1992
**          (1) Minor format corrections to putInputDataToFile()
**          (2) Added assimilant printout to putOutputDataToFile()
**          (3) Added calls to updateAssimilantPADB() to clear the
**          list entries and update the list
**      V2.0-9  Mark S. Ghiorso  January 29, 1992
**          Corrected definition of "output" to be a pointer to FILE
**      V2.0-10 Mark S. Ghiorso  February 19, 1992
**          (1) Corrected minor casting violations for ANSI C compliance
**          (2) Removed global dependence on arg_set
**      V2.0-11 Mark S. Ghiorso  March 14, 1992
**          (1) Corrected error in write file routines to avoid NULL
**          memory reference in the event that "title" is undefined
**          (2) same as 2.0-11.1 for uninitialized tpValues structure
**          (3) same as 2.0-11.1 for included solids structure
**      V2.0-12 Mark S. Ghiorso  March 16, 1992
**           (1) Updated "output file" function to report thermodynamic
**           properties and viscosity data
**      V2.0-13 Mark S. Ghiorso  April 30, 1992
**          (1) Corrected error in reading value of "Final Pressure"
**      V2.0-14 Mark S. Ghiorso  May 4, 1992
**          Added output of oxygen content and extensive thermodynamiic
**          properties of oxygen when reaction path is mu O2 constrained
**      V2.1-1  Mark S. Ghiorso  July 13, 1992
**          Added formula output to solid description
**      V2.1-2  Mark S. Ghiorso  September 29, 1992
**          Converted TextField to Text widgets as a bug workaround
**          for DECWindows Motif V 1.1
**      V2.1-3  Mark S. Ghiorso  June 14, 1993
**          Added code for additional fo2 path constraints
**      V2.1-4  Mark S. Ghiorso  September 21, 1993
**          XtFree -> XmStringFree
**      V2.1-5  Mark S. Ghiorso  September 29, 1993
**          Modified call to realloc to catch zero pointer (SPARC port)
**      V2.2-1  Mark S. Ghiorso  September 30, 1993
**          Added code to export DELTAGRAPH compatible tables for
**          MacIntosh plotting (may be turned off by undef MAKE_TABLES).
**                   October 5, 1993
**          Finished modifications
**      V2.2-2  Mark S. Ghiorso  January 17, 1994
**          (1) Fixed an error on read in Assimilant liquid composition
**          7 -> 12 character length problem
**          (2) corr initialization of
**          assimilantValues[npc+nc+ASSIM_PADB_INDEX_LIQUID_MASS].ns
**          when liquid mass is read in. Note that the code assumes
**          only liquid is assimilated
**      V3.0-1  Mark S. Ghiorso  May 11, 1994
**          (1) Modified calls to *cpmix and *vmix to reflect additional
**          derivatives
**          (2) Added reference to isentropic constraints
**                   June 9, 1994
**          (3) Added initialization of silminState->refEnthalpy,
**          silminState->refEntropy and silminState->refVolume
**                   July 2, 1994
**          (4) Modified initialization of silminState->refEnthalpy,
**          silminState->refEntropy and silminState->refVolume so that
**          when constrained to an fo2 path, these quatities are
**          reset on each pass through the program.
**          (5) Revoked 3.0-1.4
**      V3.1-1  Mark S. Ghiorso  March 4, 1995
**          (1) Modified MELTS_DEMO code for putOutputDataToFile()
**      V3.2-1  Mark S. Ghiorso  March 22, 1995
**          (1) Added code to putInputDataToFile() to output H, S, V
**          increments and related data
**          (2) Added code to getInputDataFromFile() to input H, S, V
**          increments and related data
**      V4.0-1  Paul D. Asimow April 26, 1995
**          (1) Make sensible printouts without liquid
**--
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "silmin.h"           /*SILMIN structures include file      */

#define MAKE_TABLES

#define REALLOC(x, y) (((x) == NULL) ? malloc(y) : realloc((x), (y)))
#define REC   134       /* Maximal record length for input/output files */

#ifndef BATCH_VERSION

#include <Xm/Text.h>
#include <Xm/ToggleBG.h>

#include "vframe.h"
#include "vheader.h"
#include "vlist.h"

#include "interface.h"        /*Specific external declarations      */

#define ABORT(string1, string2) \
    cstring1 = XmStringCreateLtoR(string1, "ISO8859-1"); \
    cstring2 = XmStringCreateLtoR(string2, "ISO8859-1"); \
    cstring3 = XmStringConcat(cstring1, cstring2); \
    i = 0; \
    XtSetArg(arg_set[i], XmNmessageString, cstring3); i++; \
    XtSetValues(message, arg_set, i); \
    XtManageChild(message); \
    XmStringFree(cstring1); \
    XmStringFree(cstring2); \
    XmStringFree(cstring3);

#define READ_ERROR \
    cstring1 = XmStringCreateLtoR( \
    "Error during input file read, offending record:\n", "ISO8859-1"); \
    cstring2 = XmStringCreateLtoR(line, "ISO8859-1"); \
    cstring3 = XmStringConcat(cstring1, cstring2); \
    i = 0; \
    XtSetArg(arg_set[i], XmNmessageString, cstring3); i++; \
    XtSetValues(message, arg_set, i); \
    XtManageChild(message); \
    XmStringFree(cstring1); \
    XmStringFree(cstring2); \
    XmStringFree(cstring3); \
    return GET_INPUT_ERROR_BAD_READ;

#define UPDATE_DISPLAY \
    while (XtAppPending(app_context) == XtIMXEvent) { \
    XtAppNextEvent(app_context, &event); \
    XtDispatchEvent(&event); \
    }

/******************************************************************************
 * Global functions and statically global variables
 ******************************************************************************/

/*
static char *defFileName;

int getInputDataFromFile(char *fileName)
int putInputDataToFile(char *fileName)

*/

#else

/* From here... */
#ifndef LIBXML_STATIC
#include <libxml/encoding.h>
#include <libxml/xmlschemas.h>
#include <libxml/xmlschemastypes.h>
#include <libxml/xmlwriter.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#endif

#include "interface.h"

static SilminState *previousSilminState;

#ifdef RHYOLITE_ADJUSTMENTS
#define RELEASE "rhyolite-MELTS (1.0.2, 1.1.0, 1.2.0) pMELTS (5.6.1) - (" __DATE__ " - " __TIME__ ")"
#else
#define RELEASE "()(p)(x)Melts (MELTS V5.6.0) - (" __DATE__ " - " __TIME__ ")"
#endif

/* ... to here from interface.c */

#endif /* BATCH_VERSION */

/* adapted from GUI version */
/* Storage for display (e.g. tpValues) used as flags instead */
int putInputDataToFile(char *fileName, int tpValues, int includedSolids, int assimilantValues)
{
    int i,j ;
    FILE *output;

/* -> Check validity of file name */

    if (fileName != (char *) NULL) {
    if((output = fopen (fileName, "w")) == NULL) {
        printf("Error in SILMIN file output procedure. Cannot open file: %s\n", fileName);
        return GET_INPUT_ERROR_BAD_FILE;
    }


    /* NEED TO CHECK THIS BIT */

    /* if (silminInputData.name != (char *) NULL) free(silminInputData.name);
    silminInputData.name = (char *) malloc((unsigned) (strlen(fileName)+1)*sizeof(char));
    (void) strcpy(silminInputData.name, fileName); */
    } else if (silminInputData.name != (char *) NULL) {
    if((output = fopen (silminInputData.name, "w")) == NULL) {
        printf("Error in SILMIN file output procedure. Cannot open file: %s\n", silminInputData.name);
        return GET_INPUT_ERROR_BAD_FILE;
    }
    } else {
    printf("Error in SILMIN file output procedure. Cannot open file: no file specified!\n");
    return GET_INPUT_ERROR_BAD_FILE;
    }

/* -> title record */
    if (silminInputData.title != NULL) fprintf(output, "Title: %s\n", silminInputData.title);
    else fprintf(output, "Title: Dummy Title\n");

/* -> initial composition record */
    if (silminState != NULL) {
    for (i=0; i<nc; i++) if ((silminState->bulkComp)[i] != 0.0)
        fprintf(output, "Initial Composition: %s %.4f\n", bulkSystem[i].label, (silminState->bulkComp)[i]*bulkSystem[i].mw);
    }

/* -> T,P records */
    if (tpValues) {
/* -> initial temperature record */
    fprintf(output, "Initial Temperature: %.2f\n", silminState->dspTstart);
/* -> final temperature record */
    fprintf(output, "Final Temperature: %.2f\n",   silminState->dspTstop);
/* -> initial pressure record */
    fprintf(output, "Initial Pressure: %.2f\n",    silminState->dspPstart);
/* -> final pressure record */
    fprintf(output, "Final Pressure: %.2f\n",      silminState->dspPstop);

    if (silminState->isenthalpic) {
/* -> increment pressure record */
        fprintf(output, "Increment Pressure: %.2f\n", silminState->dspPinc);
/* -> increment enthalpy record */
        fprintf(output, "Increment Enthalpy: %.2f\n", silminState->dspHinc);
/* -> dp/dH record */
        fprintf(output, "dp/dH: %.2f\n",          silminState->dspDPDH);

    } else if (silminState->isentropic) {
/* -> increment pressure record */
        fprintf(output, "Increment Pressure: %.2f\n", silminState->dspPinc);
/* -> increment entropy record */
        fprintf(output, "Increment Entropy: %.2f\n",  silminState->dspSinc);
/* -> dp/dS record */
        fprintf(output, "dp/dS: %.2f\n",          silminState->dspDPDS);

    } else if (silminState->isochoric) {
/* -> increment temperature record */
        fprintf(output, "Increment Temperature: %.2f\n", silminState->dspTinc);
/* -> increment volume record */
        fprintf(output, "Increment Volume: %.2f\n",      silminState->dspVinc);
/* -> dt/dV record */
#ifdef PHMELTS_ADJUSTMENTS
        fprintf(output, "dt/dV: %.2f\n",         silminState->dspDTDV);
#else
        fprintf(output, "dt/dV: %.2f\n", (silminState->dspDVDt != 0.0) ? 1.0 / silminState->dspDVDt : 0.0);
#endif
    } else {
/* -> increment temperature record */
        fprintf(output, "Increment Temperature: %.2f\n", silminState->dspTinc);
/* -> increment pressure record */
        fprintf(output, "Increment Pressure: %.2f\n",    silminState->dspPinc);
/* -> dp/dt record */
        fprintf(output, "dp/dt: %.2f\n",         silminState->dspDPDt);
    }

   } else {

/* -> initial temperature record */
    fprintf(output, "Initial Temperature: %.2f\n", silminState->T - 273.15);
/* -> final temperature record */
    fprintf(output, "Final Temperature: %.2f\n",   silminState->T - 273.15);
/* -> initial pressure record */
    fprintf(output, "Initial Pressure: %.2f\n",    silminState->P);
/* -> final pressure record */
    fprintf(output, "Final Pressure: %.2f\n",      silminState->P);

/* -> increment temperature record */
    fprintf(output, "Increment Temperature: 0.00\n");
/* -> increment pressure record */
    fprintf(output, "Increment Pressure: 0.00\n");
/* -> dp/dt record */
    fprintf(output, "dp/dt: 0.00\n");

    }

/* -> log fo2 path record */
    fprintf(output, "log fo2 Path: ");
    if (silminState->fo2Path == FO2_NONE) fprintf(output, "None\n");
    else if (silminState->fo2Path == FO2_QFM)    fprintf(output, "FMQ\n");
    else if (silminState->fo2Path == FO2_COH)    fprintf(output, "COH\n");
    else if (silminState->fo2Path == FO2_NNO)    fprintf(output, "NNO\n");
    else if (silminState->fo2Path == FO2_IW)     fprintf(output, "IW\n");
    else if (silminState->fo2Path == FO2_HM)     fprintf(output, "HM\n");
#ifdef PHMELTS_ADJUSTMENTS
    else if (silminState->fo2Path == FO2_ABS) fprintf(output, "Absolute\n");
#else
    else if (silminState->fo2Path == FO2_QFM_P3) fprintf(output, "+3FMQ\n");
    else if (silminState->fo2Path == FO2_QFM_P2) fprintf(output, "+2FMQ\n");
    else if (silminState->fo2Path == FO2_QFM_P1) fprintf(output, "+1FMQ\n");
    else if (silminState->fo2Path == FO2_QFM_M1) fprintf(output, "-1FMQ\n");
    else if (silminState->fo2Path == FO2_QFM_M2) fprintf(output, "-2FMQ\n");
    else if (silminState->fo2Path == FO2_QFM_M3) fprintf(output, "-3FMQ\n");
    else if (silminState->fo2Path == FO2_QFM_M4) fprintf(output, "-4FMQ\n");
    else if (silminState->fo2Path == FO2_QFM_M5) fprintf(output, "-5FMQ\n");
    else if (silminState->fo2Path == FO2_QFM_M6) fprintf(output, "-6FMQ\n");
    else if (silminState->fo2Path == FO2_QFM_M7) fprintf(output, "-7FMQ\n");
    else if (silminState->fo2Path == FO2_QFM_M8) fprintf(output, "-8FMQ\n");
    else if (silminState->fo2Path == FO2_QFM_M9) fprintf(output, "-9FMQ\n");
    else if (silminState->fo2Path == FO2_QFM_P0_5) fprintf(output, "+0.5FMQ\n");
    else if (silminState->fo2Path == FO2_QFM_P1_5) fprintf(output, "+1.5FMQ\n");
#endif

/* -> suppress a solid phase record */
    if (includedSolids) {
    for (i=0, j=0; i<npc; i++) if (solids[i].type == PHASE) {
        if (silminState->incSolids[j]) fprintf(output, "Suppress: %s\n", solids[i].label);
        j++;
    }
    }

/* -> mode record */
    if (tpValues) {
    if (silminState->fractionateSol)  fprintf(output, "Mode: Fractionate Solids\n");
    if (silminState->fractionateLiq)  fprintf(output, "Mode: Fractionate Liquids\n");
    if (silminState->fractionateFlu)  fprintf(output, "Mode: Fractionate Fluids\n");
#ifdef PHMELTS_ADJUSTMENTS
    if (silminState->incLiquids > 1)    fprintf(output, "Mode: Multiple Liquids\n");
#else
    if (silminState->multipleLiqs)    fprintf(output, "Mode: Multiple Liquids\n");
#endif
    if (silminState->isenthalpic)     fprintf(output, "Mode: Isenthalpic\n");
    if (silminState->isentropic)      fprintf(output, "Mode: Isentropic\n");
    if (silminState->isochoric)       fprintf(output, "Mode: Isochoric\n");
    }

/* -> assimilant phase record */
    if (assimilantValues && (silminState->dspAssimMass > 0.0) &&
        (silminState->dspAssimT > 0.0) && (silminState->dspAssimInc > 0.0)) {

    fprintf(output, "Assimilant: Temperature %.2f\n", silminState->dspAssimT);
    fprintf(output, "Assimilant: Mass %.2f\n",    silminState->dspAssimMass);
    fprintf(output, "Assimilant: Increments %d\n",  silminState->dspAssimInc);
    fprintf(output, "Assimilant: Liquid Mass %.2f\n", silminState->dspAssimLiqM);
    if      (silminState->dspAssimUnits == ASSIM_PADB_UNITS_VOLUME) fprintf(output, "Assimilant: Units Vol %%\n");
    else if (silminState->dspAssimUnits == ASSIM_PADB_UNITS_WEIGHT) fprintf(output, "Assimilant: Units Wt %%\n");

    for (i=0; i<npc; i++) for (j=0; j<silminState->nDspAssimComp[i]; j++)
        fprintf(output, "Assimilant: %s %.2f\n", solids[i].label, (silminState->dspAssimComp[i])[j]);

    for (i=0; i<nc; i++) for (j=0; j<silminState->nDspAssimComp[npc+i]; j++)
        fprintf(output, "Assimilant: %s %.2f\n", bulkSystem[i].label, (silminState->dspAssimComp[npc+i])[j]);
    }

/* -> Close and discard file, display status entries, and return */
    fclose(output);

    printf("Input description written to file %s.\n", silminInputData.name);
    return GET_INPUT_SUCCESS;

}

int putOutputDataToFile(char *fileName)
{
    static FILE *output;
#ifdef MAKE_TABLES
    static FILE *tableLiq;
    static const char *liquidFile = "melts-liquid.tbl";
    static FILE **tableSol;
    static int rowIndex = 0;
#endif
    static double *m, *r, *oxVal;
    double moles,    viscosity, gLiq, hLiq, sLiq, vLiq, cpLiq, dvdtLiq, dvdpLiq, mLiq,
     mass,     totalMass,
     gibbsEnergy,  totalGibbsEnergy,
     enthalpy,     totalEnthalpy,
     entropy,      totalEntropy,
     volume,       totalVolume,
     heatCapacity, totalHeatCapacity,
	      dVolumeDt,    totaldVolumeDt,
	      dVolumeDp,    totaldVolumeDp,
     oxSum,
     fo2Delta;
    int i, j, nl, ns;
    int hasLiquid = (silminState->liquidMass != 0.0);
#ifndef BATCH_VERSION
    XmString  cstring1, cstring2, cstring3;
    Arg arg_set[20];
#endif

#ifndef BATCH_VERSION
    if(output == NULL && (output = fopen (meltsEnviron.OUTPUT_FILE, "w")) == NULL) {
    ABORT("Error in SILMIN file output procedure. Cannot open file:\n", meltsEnviron.OUTPUT_FILE)
    return GET_INPUT_ERROR_BAD_FILE;
    }
#else
    if(output == NULL && (output = fopen ("melts.out", "w")) == NULL) {
    printf("Error in SILMIN file output procedure. Cannot open file: melts.out!\n");
    return FALSE;
    }
#endif

#ifdef PHMELTS_ADJUSTMENTS
    if (previousSilminState == NULL) {
        previousSilminState = allocSilminStatePointer();
        if ((silminState->fractionateSol || silminState->fractionateFlu) && previousSilminState->fracSComp == (double **) NULL) {
        previousSilminState->fracSComp    = (double **) calloc((unsigned) npc, sizeof(double *));
        previousSilminState->nFracCoexist = (int *) calloc((unsigned) npc, sizeof(int));
        }
        if (silminState->fractionateLiq && previousSilminState->fracLComp == (double *) NULL) {
        previousSilminState->fracLComp = (double *) calloc((unsigned) nlc, sizeof(double));
        }
    }
    /* Use fracMass, fracSComp and fracLComp only for most recent fractionation (zeroed in silmin.c)*/
    /* Frac options may have changed at menu (unlike Melts-batch) */
    if (previousSilminState == NULL) {
        previousSilminState = copySilminStateStructure(silminState, previousSilminState);
        previousSilminState->fracMass = 0.0;/* should already be zero but just in case */
    }
    else {
        previousSilminState->fractionateSol = silminState->fractionateSol;
        previousSilminState->fractionateFlu = silminState->fractionateFlu;
        previousSilminState->fractionateLiq = silminState->fractionateLiq;
    }
    if ((previousSilminState->fractionateSol || previousSilminState->fractionateFlu) && previousSilminState->fracSComp == (double **) NULL) {
        previousSilminState->fracSComp    = (double **) calloc((unsigned) npc, sizeof(double *));
        previousSilminState->nFracCoexist = (int *) calloc((unsigned) npc, sizeof(int));
    }
    if (previousSilminState->fractionateLiq && previousSilminState->fracLComp == (double *) NULL) {
        previousSilminState->fracLComp = (double *) calloc((unsigned) nlc, sizeof(double));
    }
#endif


#ifdef MAKE_TABLES
    rowIndex++;
    if (tableLiq == NULL && hasLiquid) {
    if ((tableLiq = fopen (liquidFile, "w")) == NULL) {
#ifndef BATCH_VERSION
        ABORT("Error in SILMIN file output procedure. Cannot open file:\n", (char *) liquidFile)
        return GET_INPUT_ERROR_BAD_FILE;
#else
        printf("Error in SILMIN file output procedure. Cannot open file: %s!\n", (char *) liquidFile);
        return FALSE;
#endif /* BATCH_VERSION */
    }
    fprintf(tableLiq, "Index,T (C),P (kbars),log(10) f O2");
    fprintf(tableLiq, ",liq mass (gm),liq rho (gm/cc)");
    for (i=0; i<nc; i++) fprintf(tableLiq, ",wt%% %s", bulkSystem[i].label);
    fprintf(tableLiq, ",liq G (kJ),liq H (kJ),liq S (J/K),liq V (cc),liq Cp (J/K)");
    for (i=0; i<nlc; i++) fprintf(tableLiq, ",activity %s", liquid[i].label);
    fprintf(tableLiq, ",liq vis (log 10 poise),sol mass (gm),sol rho (gm/cc)");
    fprintf(tableLiq, ",sol G (kJ),sol H (kJ),sol S (J/K),sol V (cc),sol Cp (J/K)");
    fprintf(tableLiq, ",sys G (kJ),sys H (kJ),sys S (J/K),sys V (cc),sys Cp (J/K)");
    fprintf(tableLiq, ",sys dVdT (cc/K),sys dVdP (cc/bar),sys alpha (1/K),sys beta (1/bar)");
    fprintf(tableLiq, ",liq dVdT (cc/K),liq dVdP (cc/bar),liq alpha (1/K),liq beta (1/bar)");
    fprintf(tableLiq, "\n");
    }
    if (tableSol == NULL) {
    tableSol  = (FILE **) calloc((unsigned) npc, sizeof(FILE *));
    }
#endif /* MAKES_TABLES */
    if (r     == NULL) r     = (double *) malloc((unsigned) nlc*sizeof(double));
    if (m     == NULL) m     = (double *) malloc((unsigned) nlc*sizeof(double));
    if (oxVal == NULL) oxVal = (double *) malloc((unsigned) nc*sizeof(double));

    fprintf(output, "\n**********----------**********\nTitle: %s\n\n",
#ifndef BATCH_VERSION
    (silminInputData.title != NULL) ? silminInputData.title : "Dummy Title");
#elif defined(ALPHAMELTS_UPDATE_SYSTEM)
    (silminInputData.title != NULL) ? silminInputData.title : "Dummy Title");
#else
    "Dummy Title");
#endif
    fprintf(output, "T = %.2f (C)  P = %.3f (kbars)  log(10) f O2 = %.2f  ",
    silminState->T-273.15, silminState->P/1000.0, silminState->fo2);

    fo2Delta = silminState->fo2Delta;
    silminState->fo2Delta = 0;

    fprintf(output, "delta HM = %.2f  NNO = %.2f  QFM = %.2f  COH = %.2f  IW = %.2f\n\n",
    silminState->fo2 - getlog10fo2(silminState->T, silminState->P, FO2_HM),
    silminState->fo2 - getlog10fo2(silminState->T, silminState->P, FO2_NNO),
    silminState->fo2 - getlog10fo2(silminState->T, silminState->P, FO2_QFM),
    silminState->fo2 - getlog10fo2(silminState->T, silminState->P, FO2_COH),
    silminState->fo2 - getlog10fo2(silminState->T, silminState->P, FO2_IW));

    silminState->fo2Delta = fo2Delta;

    fprintf(output, "Constraint Flags: ");

    if      (silminState->fo2Path == FO2_HM)     fprintf(output, "fO2 path = HM  ");
    else if (silminState->fo2Path == FO2_QFM)    fprintf(output, "fO2 path = QFM  ");
    else if (silminState->fo2Path == FO2_COH)    fprintf(output, "fO2 path = C-COH  ");
    else if (silminState->fo2Path == FO2_NNO)    fprintf(output, "fO2 path = NNO  ");
    else if (silminState->fo2Path == FO2_IW)     fprintf(output, "fO2 path = IW  ");
    else if (silminState->fo2Path == FO2_QFM_P3) fprintf(output, "fO2 path = QFM+3  ");
    else if (silminState->fo2Path == FO2_QFM_P2) fprintf(output, "fO2 path = QFM+2  ");
    else if (silminState->fo2Path == FO2_QFM_P1) fprintf(output, "fO2 path = QFM+1  ");
    else if (silminState->fo2Path == FO2_QFM_M1) fprintf(output, "fO2 path = QFM-1  ");
    else if (silminState->fo2Path == FO2_QFM_M2) fprintf(output, "fO2 path = QFM-2  ");
    else if (silminState->fo2Path == FO2_QFM_M3) fprintf(output, "fO2 path = QFM-3  ");
    else if (silminState->fo2Path == FO2_QFM_M4) fprintf(output, "fO2 path = QFM-4  ");
    else if (silminState->fo2Path == FO2_QFM_M5) fprintf(output, "fO2 path = QFM-5  ");
    else if (silminState->fo2Path == FO2_QFM_M6) fprintf(output, "fO2 path = QFM-6  ");
    else if (silminState->fo2Path == FO2_QFM_M7) fprintf(output, "fO2 path = QFM-7  ");
    else if (silminState->fo2Path == FO2_QFM_M8) fprintf(output, "fO2 path = QFM-8  ");
    else if (silminState->fo2Path == FO2_QFM_M9) fprintf(output, "fO2 path = QFM-9  ");
    else if (silminState->fo2Path == FO2_QFM_P0_5) fprintf(output, "fO2 path = QFM+0.5  ");
    else if (silminState->fo2Path == FO2_QFM_P1_5) fprintf(output, "fO2 path = QFM+1.5  ");

    if (silminState->fractionateSol) fprintf(output, "Fractionate Solids  ");
    if (silminState->fractionateFlu) fprintf(output, "Fractionate Fluids  ");
    if (silminState->fractionateLiq) fprintf(output, "Fractionate Liquids  ");
#ifdef PHMELTS_ADJUSTMENTS
    if (silminState->incLiquids > 1) fprintf(output, "Multiple Liquids  ");
#else
    if (silminState->multipleLiqs)   fprintf(output, "Multiple Liquids  ");
#endif
    if (silminState->isenthalpic)    fprintf(output, "Isenthalpic Path  ");
    if (silminState->isentropic)     fprintf(output, "Isentropic Path  ");
    if (silminState->isochoric)      fprintf(output, "Isochoric Path  ");
    if (silminState->assimilate)     fprintf(output, "Assimilation  ");

    fprintf(output, "\n\n");

    if (hasLiquid) {
    gLiq = 0.0; hLiq = 0.0; sLiq = 0.0; vLiq = 0.0; cpLiq = 0.0; mLiq = 0.0; dvdtLiq = 0.0; dvdpLiq = 0.0;
    for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
        conLiq(SECOND, THIRD, silminState->T, silminState->P, NULL, silminState->liquidComp[nl], r, NULL, NULL, NULL, NULL);

        gmixLiq (FIRST, silminState->T, silminState->P, r, &gibbsEnergy,  NULL, NULL);
        hmixLiq (FIRST, silminState->T, silminState->P, r, &enthalpy,     NULL);
        smixLiq (FIRST, silminState->T, silminState->P, r, &entropy,      NULL, NULL, NULL);
        vmixLiq (FIRST | FOURTH | FIFTH, silminState->T, silminState->P, r, &volume,
            NULL, NULL, &dVolumeDt, &dVolumeDp, NULL, NULL, NULL, NULL, NULL, NULL);
        cpmixLiq(FIRST, silminState->T, silminState->P, r, &heatCapacity, NULL, NULL);
        visLiq  (FIRST, silminState->T, silminState->P, r, &viscosity);

        for (i=0, moles=0.0; i<nlc; i++) moles +=  (silminState->liquidComp)[nl][i];
        gibbsEnergy *= moles; enthalpy	 *= moles; entropy   *= moles;
        volume	  *= moles; heatCapacity *= moles; dVolumeDt *= moles;
        dVolumeDp   *= moles;

        for (i=0; i<nlc; i++) {
        gibbsEnergy  += (silminState->liquidComp)[nl][i]*(liquid[i].cur).g;
        enthalpy     += (silminState->liquidComp)[nl][i]*(liquid[i].cur).h;
        entropy      += (silminState->liquidComp)[nl][i]*(liquid[i].cur).s;
        volume       += (silminState->liquidComp)[nl][i]*(liquid[i].cur).v;
        heatCapacity += (silminState->liquidComp)[nl][i]*(liquid[i].cur).cp;
        dVolumeDt    += (silminState->liquidComp)[nl][i]*(liquid[i].cur).dvdt;
        dVolumeDp    += (silminState->liquidComp)[nl][i]*(liquid[i].cur).dvdp;
        }

        for (i=0, oxSum=0.0; i<nc; i++) {
    	    for (j=0, oxVal[i]=0.0; j<nlc; j++) oxVal[i] += (liquid[j].liqToOx)[i]*(silminState->liquidComp)[nl][j];
        	oxVal[i] *= bulkSystem[i].mw;
    	    oxSum	 += oxVal[i];
        }
        if (oxSum != 0.0) for (i=0; i<nc; i++) oxVal[i] /= oxSum;

        gLiq += gibbsEnergy; hLiq += enthalpy; sLiq += entropy; vLiq += volume; cpLiq += heatCapacity; mLiq += oxSum;
        dvdtLiq += dVolumeDt; dvdpLiq += dVolumeDp;

        fprintf(output, "%-15.15s  mass = %.2f (gm)  density = %.2f (gm/cc)", "Liquid", oxSum, (volume == 0.0) ? 0.0 : oxSum/(10.0*volume));
        fprintf(output, "  viscosity = %.2f (log 10 poise)", viscosity);
        fprintf(output, "     (analysis in wt %%)\n");

        fprintf(output, " 		");
        fprintf(output, "G = %.2f (J)  ",    gibbsEnergy);
        fprintf(output, "H = %.2f (J)  ",    enthalpy);
        fprintf(output, "S = %.2f (J/K)  ",  entropy);
        fprintf(output, "V = %.2f (cc)  ",   volume*10.0);
        fprintf(output, "Cp = %.2f (J/K)  ", heatCapacity);
        fprintf(output, "\n");

        fprintf(output, "     ");
        for (i=0; i<nc; i++) fprintf(output, "%7.7s", bulkSystem[i].label);
        fprintf(output, "\n");
        fprintf(output, "     ");
        for (i=0; i<nc; i++) fprintf(output, "  %5.2f", oxVal[i]*100.0);
        fprintf(output, "\n");

#ifdef MAKE_TABLES
    fprintf(tableLiq, "%d,%.2f,%.3f,%.3f", rowIndex, silminState->T-273.15, silminState->P/1000.0, silminState->fo2);
    fprintf(tableLiq, ",%.13e,%.4f", oxSum, (volume == 0.0) ? 0.0 : oxSum/(10.0*volume));
    for (i=0; i<nc; i++) fprintf(tableLiq, ",%.4f", oxVal[i]*100.0);
    fprintf(tableLiq, ",%.13e,%.13e,%.13e,%.13e,%.13e", gibbsEnergy/1000.0, enthalpy/1000.0, entropy, volume*10.0, heatCapacity);
    actLiq(FIRST, silminState->T, silminState->P, r, m, NULL, NULL, NULL);
    for (i=0; i<nlc; i++) fprintf(tableLiq, ",%.7e", m[i]);
    fprintf(tableLiq, ",%.4f", viscosity);
#endif
    }
    } else {
    gLiq = 0.0; hLiq = 0.0; sLiq = 0.0; vLiq = 0.0; cpLiq = 0.0; mLiq = 0.0; dvdtLiq = 0.0; dvdpLiq = 0.0;
    }

    for (j=0, totalMass=0.0, totalGibbsEnergy=0.0, totalEnthalpy=0.0,
    totalEntropy=0.0, totalVolume=0.0, totalHeatCapacity=0.0, totaldVolumeDt=0.0, totaldVolumeDp=0.0; j<npc; j++)
    for (ns=0; ns<(silminState->nSolidCoexist)[j]; ns++) {
    if (solids[j].na == 1) {
        mass           = (silminState->solidComp)[j][ns]*solids[j].mw;
        gibbsEnergy    = (silminState->solidComp)[j][ns]*(solids[j].cur).g;
        enthalpy       = (silminState->solidComp)[j][ns]*(solids[j].cur).h;
        entropy        = (silminState->solidComp)[j][ns]*(solids[j].cur).s;
        volume         = (silminState->solidComp)[j][ns]*(solids[j].cur).v;
        heatCapacity       = (silminState->solidComp)[j][ns]*(solids[j].cur).cp;
        dVolumeDt      = (silminState->solidComp)[j][ns]*(solids[j].cur).dvdt;
        dVolumeDp      = (silminState->solidComp)[j][ns]*(solids[j].cur).dvdp;
        totalMass     += (silminState->solidComp)[j][ns]*solids[j].mw;
        totalGibbsEnergy  += (silminState->solidComp)[j][ns]*(solids[j].cur).g;
        totalEnthalpy     += (silminState->solidComp)[j][ns]*(solids[j].cur).h;
        totalEntropy      += (silminState->solidComp)[j][ns]*(solids[j].cur).s;
        totalVolume       += (silminState->solidComp)[j][ns]*(solids[j].cur).v;
        totalHeatCapacity += (silminState->solidComp)[j][ns]*(solids[j].cur).cp;
        totaldVolumeDt    += (silminState->solidComp)[j][ns]*(solids[j].cur).dvdt;
        totaldVolumeDp    += (silminState->solidComp)[j][ns]*(solids[j].cur).dvdp;
        fprintf(output, "\n%-15.15s  mass = %.2f (gm)  density = %.2f (gm/cc)", solids[j].label, mass, (volume == 0.0) ? 0.0 : mass/(10.0*volume));
        fprintf(output, "\n");
        fprintf(output, "         %s\n", solids[j].formula);
        fprintf(output, "         ");
        fprintf(output, "G = %.2f (J)  ", gibbsEnergy);
        fprintf(output, "H = %.2f (J)  ", enthalpy);
        fprintf(output, "S = %.2f (J/K)  ", entropy);
        fprintf(output, "V = %.2f (cc)  ", volume*10.0);
        fprintf(output, "Cp = %.2f (J/K)  ", heatCapacity);
        fprintf(output, "\n");
#ifdef MAKE_TABLES
        if (tableSol[j] == NULL) {
        int len = (int) strlen(solids[j].label);
        char *nameOfFile = (char *) calloc((unsigned) (len+5), sizeof(char));
        strcpy(nameOfFile, solids[j].label);
        for (i=0; i<len; i++) if(nameOfFile[i] == ' ') nameOfFile[i] = '-';
        nameOfFile = strncat(nameOfFile, ".tbl", 4);
        if ((tableSol[j] = fopen (nameOfFile, "w")) == NULL) {
#ifndef BATCH_VERSION
            ABORT("Error in SILMIN file output procedure. Cannot open file:\n", nameOfFile)
            return GET_INPUT_ERROR_BAD_FILE;
#else
            printf("Error in SILMIN file output procedure. Cannot open file: %s!\n", nameOfFile);
            return FALSE;
#endif /* BATCH_VERSION */
        }
        free(nameOfFile);
        fprintf(tableSol[j], "Index,T (C),P (kbars),log(10) f O2");
        fprintf(tableSol[j], ",mass (gm),rho (gm/cc)");
        for (i=0; i<nc; i++) fprintf(tableSol[j], ",wt%% %s", bulkSystem[i].label);
        fprintf(tableSol[j], ",G (kJ),H (kJ),S (J/K),V (cc),Cp (J/K)");
        fprintf(tableSol[j], "\n");
        }
        fprintf(tableSol[j], "%d,%.2f,%.3f,%.3f", rowIndex, silminState->T-273.15, silminState->P/1000.0, silminState->fo2);
        fprintf(tableSol[j], ",%.4f,%.4f", mass, (volume == 0.0) ? 0.0 : mass/(10.0*volume));

        for (i=0, oxSum=0.0; i<nc; i++) {
        oxVal[i]  = (solids[j].solToOx)[i]*bulkSystem[i].mw;
    	oxSum	 += oxVal[i];
        }
        if (oxSum != 0.0) for (i=0; i<nc; i++) oxVal[i] /= oxSum;
        for (i=0; i<nc; i++) fprintf(tableSol[j], ",%.4f", oxVal[i]*100.0);

        fprintf(tableSol[j], ",%.3f,%.3f,%.3f,%.3f,%.3f", gibbsEnergy/1000.0, enthalpy/1000.0, entropy, volume*10.0, heatCapacity);
        fprintf(tableSol[j], "\n");
#endif /* MAKE_TABLES */
    } else {
        char *formula;
        for (i=0, mass=0.0; i<solids[j].na; i++) {
        m[i] = (silminState->solidComp)[j+1+i][ns];
        mass += (silminState->solidComp)[j+1+i][ns]*solids[j+1+i].mw;
        }
        (*solids[j].convert)(SECOND, THIRD, silminState->T, silminState->P, NULL, m, r, NULL, NULL, NULL, NULL, NULL);
        (*solids[j].display)(FIRST, silminState->T, silminState->P, r, &formula);
        (*solids[j].gmix) (FIRST, silminState->T, silminState->P, r, &gibbsEnergy,  NULL, NULL, NULL);
        (*solids[j].hmix) (FIRST, silminState->T, silminState->P, r, &enthalpy);
        (*solids[j].smix) (FIRST, silminState->T, silminState->P, r, &entropy,      NULL, NULL);
        (*solids[j].vmix) (FIRST | FOURTH | FIFTH, silminState->T, silminState->P, r, &volume,
                NULL, NULL, &dVolumeDt, &dVolumeDp, NULL, NULL,  NULL, NULL, NULL);
        (*solids[j].cpmix)(FIRST, silminState->T, silminState->P, r, &heatCapacity, NULL, NULL);
        gibbsEnergy  *= (silminState->solidComp)[j][ns];
        enthalpy     *= (silminState->solidComp)[j][ns];
        entropy      *= (silminState->solidComp)[j][ns];
        volume       *= (silminState->solidComp)[j][ns];
        heatCapacity *= (silminState->solidComp)[j][ns];
        dVolumeDt    *= (silminState->solidComp)[j][ns];
        dVolumeDp    *= (silminState->solidComp)[j][ns];
        for (i=0; i<solids[j].na; i++) {
        gibbsEnergy  += m[i]*(solids[j+1+i].cur).g;
        enthalpy     += m[i]*(solids[j+1+i].cur).h;
        entropy      += m[i]*(solids[j+1+i].cur).s;
        volume       += m[i]*(solids[j+1+i].cur).v;
        heatCapacity += m[i]*(solids[j+1+i].cur).cp;
	        dVolumeDt    += m[i]*(solids[j+1+i].cur).dvdt;
	        dVolumeDp    += m[i]*(solids[j+1+i].cur).dvdp;
        }
        totalMass     += mass;
        totalGibbsEnergy  += gibbsEnergy;
        totalEnthalpy     += enthalpy;
        totalEntropy      += entropy;
        totalVolume       += volume;
        totalHeatCapacity += heatCapacity;
        totaldVolumeDt    += dVolumeDt;
        totaldVolumeDp    += dVolumeDp;
        fprintf(output, "\n%-15.15s  mass = %.2f (gm)  density = %.2f (gm/cc)", solids[j].label, mass, (volume == 0.0) ? 0.0 : mass/(10.0*volume));
        fprintf(output, "     (analysis in mole %%)\n");
        fprintf(output, "         %s\n", formula); free(formula);
        fprintf(output, "         ");
        fprintf(output, "G = %.2f (J)  ", gibbsEnergy);
        fprintf(output, "H = %.2f (J)  ", enthalpy);
        fprintf(output, "S = %.2f (J/K)  ", entropy);
        fprintf(output, "V = %.2f (cc)  ", volume*10.0);
        fprintf(output, "Cp = %.2f (J/K)  ", heatCapacity);
        fprintf(output, "\n");
        fprintf(output, "     ");
        for (i=0; i<solids[j].na; i++) fprintf(output, " %13.13s", solids[j+1+i].label);
        fprintf(output, "\n");
        fprintf(output, "     ");
        for (i=0; i<solids[j].na; i++) fprintf(output, " %13.2f", 100.0*m[i]/(silminState->solidComp)[j][ns]);
        fprintf(output, "\n");
#ifdef MAKE_TABLES
        if (tableSol[j] == NULL) {
        int len = (int) strlen(solids[j].label);
        char *nameOfFile = (char *) calloc((unsigned) (len+5), sizeof(char));
        strcpy(nameOfFile, solids[j].label);
        for (i=0; i<len; i++) if(nameOfFile[i] == ' ') nameOfFile[i] = '-';
        nameOfFile = strncat(nameOfFile, ".tbl", 4);
        if ((tableSol[j] = fopen (nameOfFile, "w")) == NULL) {
#ifndef BATCH_VERSION
            ABORT("Error in SILMIN file output procedure. Cannot open file:\n", nameOfFile)
            return GET_INPUT_ERROR_BAD_FILE;
#else
            printf("Error in SILMIN file output procedure. Cannot open file: %s\n", nameOfFile);
            free(nameOfFile);
            return FALSE;
#endif /* BATCH_VERSION */
        }
        free(nameOfFile);
        fprintf(tableSol[j], "Index,T (C),P (kbars),log(10) f O2");
        fprintf(tableSol[j], ",mass (gm),rho (gm/cc)");
        for (i=0; i<nc; i++) fprintf(tableSol[j], ",wt%% %s", bulkSystem[i].label);
        fprintf(tableSol[j], ",G (kJ),H (kJ),S (J/K),V (cc),Cp (J/K)");
#ifdef PHMELTS_ADJUSTMENTS
        for (i=0; i<solids[j].na; i++) fprintf(tableSol[j], ",%s", solids[j+1+i].label);
#else
        for (i=0; i<solids[j].na; i++) fprintf(tableSol[j], ",%13.13s", solids[j+1+i].label);
#endif
        fprintf(tableSol[j], "\n");
        }
        fprintf(tableSol[j], "%d,%.2f,%.3f,%.3f", rowIndex, silminState->T-273.15, silminState->P/1000.0, silminState->fo2);
        fprintf(tableSol[j], ",%.4f,%.4f", mass, (volume == 0.0) ? 0.0 : mass/(10.0*volume));

        for (i=0, oxSum=0.0; i<nc; i++) {
        int k;
	for (k=0, oxVal[i]=0.0; k<solids[j].na; k++) oxVal[i] += (solids[j+1+k].solToOx)[i]*m[k]*bulkSystem[i].mw;
    	oxSum += oxVal[i];
        }
        if (oxSum != 0.0) for (i=0; i<nc; i++) oxVal[i] /= oxSum;
        for (i=0; i<nc; i++) fprintf(tableSol[j], ",%.4f", oxVal[i]*100.0);

        fprintf(tableSol[j], ",%.3f,%.3f,%.3f,%.3f,%.3f", gibbsEnergy/1000.0, enthalpy/1000.0, entropy, volume*10.0, heatCapacity);
        for (i=0; i<solids[j].na; i++) fprintf(tableSol[j], ",%.6f", m[i]/(silminState->solidComp)[j][ns]);
        fprintf(tableSol[j], "\n");
#endif /* MAKE_TABLES */
    }
    }

    fprintf(output, "\n%-15.15s  mass = %.2f (gm)  density = %.2f (gm/cc)\n", "Total solids", totalMass, (totalVolume == 0.0) ? 0.0 : totalMass/(10.0*totalVolume));
    fprintf(output, "         ");
    fprintf(output, "G = %.2f (J)  ", totalGibbsEnergy);
    fprintf(output, "H = %.2f (J)  ", totalEnthalpy);
    fprintf(output, "S = %.2f (J/K)  ", totalEntropy);
    fprintf(output, "V = %.2f (cc)  ", totalVolume*10.0);
    fprintf(output, "Cp = %.2f (J/K)  ", totalHeatCapacity);
    fprintf(output, "\n");
#ifdef MAKE_TABLES
    if (hasLiquid) {
    fprintf(tableLiq, ",%.13e,%.4f", totalMass, (totalVolume == 0.0) ? 0.0 : totalMass/(10.0*totalVolume));
    fprintf(tableLiq, ",%.13e,%.13e,%.13e,%.13e,%.13e", totalGibbsEnergy/1000.0, totalEnthalpy/1000.0, totalEntropy, totalVolume*10.0, totalHeatCapacity);
    }
#endif

    if (silminState->isenthalpic && (silminState->refEnthalpy == 0.0)) silminState->refEnthalpy = hLiq+totalEnthalpy;
    if (silminState->isentropic  && (silminState->refEntropy  == 0.0)) silminState->refEntropy  = sLiq+totalEntropy;
    if (silminState->isochoric   && (silminState->refVolume   == 0.0)) silminState->refVolume   = vLiq+totalVolume;

#ifdef PHMELTS_ADJUSTMENTS
    if (silminState->fractionateSol || silminState->fractionateFlu || silminState->fractionateLiq)
    fprintf(output, "\nSummary of all fractionated phases: (total mass = %.2f grams)\n", previousSilminState->fracMass);
#else
    if (silminState->fractionateSol || silminState->fractionateFlu || silminState->fractionateLiq)
    fprintf(output, "\nSummary of all fractionated phases: (total mass = %.2f grams)\n", silminState->fracMass);
#endif

    if (silminState->fractionateSol || silminState->fractionateFlu) {
    for (j=0; j<npc; j++) {
        if ( silminState->fractionateSol && !silminState->fractionateFlu && !strcmp((char *) solids[j].label, "fluid")) continue;
        if (!silminState->fractionateSol &&  silminState->fractionateFlu &&  strcmp((char *) solids[j].label, "fluid")) continue;
#ifdef PHMELTS_ADJUSTMENTS
        for (ns=0; ns<(previousSilminState->nFracCoexist)[j]; ns++) {
        if (solids[j].na == 1) {
            mass  	     = (previousSilminState->fracSComp)[j][ns]*solids[j].mw;
            gibbsEnergy	     = (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).g;
            enthalpy	     = (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).h;
            entropy	     = (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).s;
            volume	     = (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).v;
            heatCapacity       = (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).cp;
            dVolumeDt      = (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).dvdt;
            dVolumeDp      = (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).dvdp;
            totalMass	    += (previousSilminState->fracSComp)[j][ns]*solids[j].mw;
            totalGibbsEnergy  += (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).g;
            totalEnthalpy     += (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).h;
            totalEntropy      += (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).s;
            totalVolume	    += (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).v;
            totalHeatCapacity += (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).cp;
            totaldVolumeDt    += (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).dvdt;
            totaldVolumeDp    += (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).dvdp;
#else
        for (ns=0; ns<(silminState->nFracCoexist)[j]; ns++) {
        if (solids[j].na == 1) {
            mass  	     = (silminState->fracSComp)[j][ns]*solids[j].mw;
            gibbsEnergy	     = (silminState->fracSComp)[j][ns]*(solids[j].cur).g;
            enthalpy	     = (silminState->fracSComp)[j][ns]*(solids[j].cur).h;
            entropy	     = (silminState->fracSComp)[j][ns]*(solids[j].cur).s;
            volume	     = (silminState->fracSComp)[j][ns]*(solids[j].cur).v;
            heatCapacity       = (silminState->fracSComp)[j][ns]*(solids[j].cur).cp;
            dVolumeDt      = (silminState->fracSComp)[j][ns]*(solids[j].cur).dvdt;
            dVolumeDp      = (silminState->fracSComp)[j][ns]*(solids[j].cur).dvdp;
            totalMass	    += (silminState->fracSComp)[j][ns]*solids[j].mw;
            totalGibbsEnergy  += (silminState->fracSComp)[j][ns]*(solids[j].cur).g;
            totalEnthalpy     += (silminState->fracSComp)[j][ns]*(solids[j].cur).h;
            totalEntropy      += (silminState->fracSComp)[j][ns]*(solids[j].cur).s;
            totalVolume	    += (silminState->fracSComp)[j][ns]*(solids[j].cur).v;
            totalHeatCapacity += (silminState->fracSComp)[j][ns]*(solids[j].cur).cp;
            totaldVolumeDt    += (silminState->fracSComp)[j][ns]*(solids[j].cur).dvdt;
            totaldVolumeDp    += (silminState->fracSComp)[j][ns]*(solids[j].cur).dvdp;
#endif
            if (mass > 0.0) {
            fprintf(output, "\n%-15.15s  mass = %.2f (gm)  density = %.2f (gm/cc)", solids[j].label, mass, (volume == 0.0) ? 0.0 : mass/(10.0*volume));
            fprintf(output, "\n");
            fprintf(output, "		    %s\n", solids[j].formula);
            fprintf(output, "		    ");
            fprintf(output, "G = %.2f (J)  ", gibbsEnergy);
            fprintf(output, "H = %.2f (J)  ", enthalpy);
            fprintf(output, "S = %.2f (J/K)  ", entropy);
            fprintf(output, "V = %.2f (cc)  ", volume*10.0);
            fprintf(output, "Cp = %.2f (J/K)  ", heatCapacity);
            fprintf(output, "\n");
            }
        } else {
            char *formula;
            for (i=0, mass=0.0; i<solids[j].na; i++) {
#ifdef PHMELTS_ADJUSTMENTS
            m[i] = (previousSilminState->fracSComp)[j+1+i][ns];
            mass += (previousSilminState->fracSComp)[j+1+i][ns]*solids[j+1+i].mw;
#else
            m[i] = (silminState->fracSComp)[j+1+i][ns];
            mass += (silminState->fracSComp)[j+1+i][ns]*solids[j+1+i].mw;
#endif
            }
            if (mass > 0.0) {
            (*solids[j].convert)(SECOND, THIRD, silminState->T, silminState->P, NULL, m, r, NULL, NULL, NULL, NULL, NULL);
            (*solids[j].display)(FIRST, silminState->T, silminState->P, r, &formula);
            (*solids[j].gmix) (FIRST, silminState->T, silminState->P, r, &gibbsEnergy,  NULL, NULL, NULL);
            (*solids[j].hmix) (FIRST, silminState->T, silminState->P, r, &enthalpy);
            (*solids[j].smix) (FIRST, silminState->T, silminState->P, r, &entropy,      NULL, NULL);
            (*solids[j].vmix) (FIRST | FOURTH | FIFTH, silminState->T, silminState->P, r, &volume,
                        NULL, NULL, &dVolumeDt, &dVolumeDp, NULL, NULL,  NULL, NULL, NULL);
            (*solids[j].cpmix)(FIRST, silminState->T, silminState->P, r, &heatCapacity, NULL, NULL);
#ifdef PHMELTS_ADJUSTMENTS
            gibbsEnergy  *= (previousSilminState->fracSComp)[j][ns];
            enthalpy     *= (previousSilminState->fracSComp)[j][ns];
            entropy      *= (previousSilminState->fracSComp)[j][ns];
            volume       *= (previousSilminState->fracSComp)[j][ns];
            heatCapacity *= (previousSilminState->fracSComp)[j][ns];
            dVolumeDt    *= (previousSilminState->fracSComp)[j][ns];
            dVolumeDp    *= (previousSilminState->fracSComp)[j][ns];
#else
            gibbsEnergy  *= (silminState->fracSComp)[j][ns];
            enthalpy     *= (silminState->fracSComp)[j][ns];
            entropy      *= (silminState->fracSComp)[j][ns];
            volume       *= (silminState->fracSComp)[j][ns];
            heatCapacity *= (silminState->fracSComp)[j][ns];
            dVolumeDt    *= (silminState->fracSComp)[j][ns];
            dVolumeDp    *= (silminState->fracSComp)[j][ns];
#endif
            for (i=0; i<solids[j].na; i++) {
                gibbsEnergy  += m[i]*(solids[j+1+i].cur).g;
                enthalpy	 += m[i]*(solids[j+1+i].cur).h;
                entropy	 += m[i]*(solids[j+1+i].cur).s;
                volume	 += m[i]*(solids[j+1+i].cur).v;
                heatCapacity += m[i]*(solids[j+1+i].cur).cp;
                dVolumeDt    += m[i]*(solids[j+1+i].cur).dvdt;
                dVolumeDp    += m[i]*(solids[j+1+i].cur).dvdp;
            }
            totalMass	    += mass;
            totalGibbsEnergy  += gibbsEnergy;
            totalEnthalpy     += enthalpy;
            totalEntropy      += entropy;
            totalVolume	    += volume;
            totalHeatCapacity += heatCapacity;
            totaldVolumeDt    += dVolumeDt;
            totaldVolumeDp    += dVolumeDp;
            fprintf(output, "\n%-15.15s  mass = %.2f (gm)  density = %.2f (gm/cc)", solids[j].label, mass, (volume == 0.0) ? 0.0 : mass/(10.0*volume));
            fprintf(output, "	(analysis in mole %%)\n");
            fprintf(output, "		    %s\n", formula); free(formula);
            fprintf(output, "		    ");
            fprintf(output, "G = %.2f (J)  ", gibbsEnergy);
            fprintf(output, "H = %.2f (J)  ", enthalpy);
            fprintf(output, "S = %.2f (J/K)  ", entropy);
            fprintf(output, "V = %.2f (cc)  ", volume*10.0);
            fprintf(output, "Cp = %.2f (J/K)  ", heatCapacity);
            fprintf(output, "\n");
            fprintf(output, "	");
            for (i=0; i<solids[j].na; i++) fprintf(output, " %13.13s", solids[j+1+i].label);
            fprintf(output, "\n");
            fprintf(output, "	");
#ifdef PHMELTS_ADJUSTMENTS
            for (i=0; i<solids[j].na; i++) fprintf(output, " %13.2f", 100.0*(previousSilminState->fracSComp)[j+1+i][ns]/(previousSilminState->fracSComp)[j][ns]);
#else
            for (i=0; i<solids[j].na; i++) fprintf(output, " %13.2f", 100.0*(silminState->fracSComp)[j+1+i][ns]/(silminState->fracSComp)[j][ns]);
#endif
       fprintf(output, "\n");
            }
        }
        }
    }
    }

    if (silminState->fractionateLiq) {
    char *formula;
    for (i=0, mass=0.0, moles=0.0; i<nlc; i++) {
        double mw;
        for (j=0, mw = 0.0; j<nc; j++) mw += (liquid[i].liqToOx)[j]*bulkSystem[j].mw;
#ifdef PHMELTS_ADJUSTMENTS
        m[i]   = (previousSilminState->fracLComp)[i];
        moles += m[i];
        mass  += (previousSilminState->fracLComp)[i]*mw;
#else
        m[i]   = (silminState->fracLComp)[i];
        moles += m[i];
        mass  += (silminState->fracLComp)[i]*mw;
#endif
    }
    if (mass > 0.0) {
        conLiq  (SECOND, THIRD, silminState->T, silminState->P, NULL, m, r, NULL, NULL, NULL, NULL);
        dispLiq (FIRST, silminState->T, silminState->P, r, &formula);
        gmixLiq (FIRST, silminState->T, silminState->P, r, &gibbsEnergy, NULL, NULL);
        hmixLiq (FIRST, silminState->T, silminState->P, r, &enthalpy, NULL);
        smixLiq (FIRST, silminState->T, silminState->P, r, &entropy, NULL, NULL, NULL);
        vmixLiq (FIRST | FOURTH | FIFTH, silminState->T, silminState->P, r, &volume,
            NULL, NULL, &dVolumeDt, &dVolumeDp, NULL, NULL,  NULL, NULL, NULL, NULL);
        cpmixLiq(FIRST, silminState->T, silminState->P, r, &heatCapacity, NULL, NULL);
        gibbsEnergy  *= moles;
        enthalpy     *= moles;
        entropy	   *= moles;
        volume	   *= moles;
        heatCapacity *= moles;
        dVolumeDt    *= moles;
        dVolumeDp    *= moles;
        for (i=0; i<nlc; i++) {
        gibbsEnergy  += m[i]*(liquid[i].cur).g;
        enthalpy     += m[i]*(liquid[i].cur).h;
        entropy      += m[i]*(liquid[i].cur).s;
        volume       += m[i]*(liquid[i].cur).v;
        heatCapacity += m[i]*(liquid[i].cur).cp;
      	    dVolumeDt    += m[i]*(liquid[i].cur).dvdt;
	        dVolumeDp    += m[i]*(liquid[i].cur).dvdp;
        }
        totalMass 	+= mass;
        totalGibbsEnergy  += gibbsEnergy;
        totalEnthalpy	+= enthalpy;
        totalEntropy	+= entropy;
        totalVolume	+= volume;
        totalHeatCapacity += heatCapacity;
        totaldVolumeDt    += dVolumeDt;
        totaldVolumeDp    += dVolumeDp;
        fprintf(output, "\n%-15.15s  mass = %.2f (gm)  density = %.2f (gm/cc)", "liquid", mass, (volume == 0.0) ? 0.0 : mass/(10.0*volume));
        fprintf(output, "     (analysis in wt %%)\n");
        fprintf(output, " 		%s\n", formula); free(formula);
        fprintf(output, " 		");
        fprintf(output, "G = %.2f (J)  ", gibbsEnergy);
        fprintf(output, "H = %.2f (J)  ", enthalpy);
        fprintf(output, "S = %.2f (J/K)  ", entropy);
        fprintf(output, "V = %.2f (cc)  ", volume*10.0);
        fprintf(output, "Cp = %.2f (J/K)  ", heatCapacity);
        fprintf(output, "\n");
    }
    }

    /* This includes the fractionated material (all counted as "solid") in estimate of system viscosity. */
    if (vLiq > totalVolume) fprintf(output,"\nViscosity of the System: %.2f (log 10 poise)\n", viscosity - 2.0*log10(1.0-2.0*totalVolume/(totalVolume+vLiq)));
    else fprintf(output,"\nViscosity of the System cannot be computed.\n");

    /* Total system including fractionated material */
    /* https://melts.ofm-research.org/Manual/UnixManHtml/output-explanations.html#SYSTEM */
    fprintf(output, "\n%-15.15s  mass = %.2f (gm)  density = %.2f (gm/cc)\n", "System", mLiq+totalMass, (totalVolume+vLiq == 0.0) ? 0.0 :
    (mLiq+totalMass)/(10.0*(vLiq+totalVolume)));
    fprintf(output, "         ");
    fprintf(output, "G = %.2f (J)  ",    gLiq+totalGibbsEnergy);
    fprintf(output, "H = %.2f (J)  ",    hLiq+totalEnthalpy);
    fprintf(output, "S = %.2f (J/K)  ",  sLiq+totalEntropy);
    fprintf(output, "V = %.2f (cc)  ",   (vLiq+totalVolume)*10.0);
    fprintf(output, "Cp = %.2f (J/K)  ", cpLiq+totalHeatCapacity);
    fprintf(output, "\n");
#ifdef MAKE_TABLES
    if (hasLiquid) {
    fprintf(tableLiq, ",%20.13e", (gLiq+totalGibbsEnergy)/1000.0);
    fprintf(tableLiq, ",%20.13e", (hLiq+totalEnthalpy)/1000.0);
    fprintf(tableLiq, ",%20.13e", (sLiq+totalEntropy));
    fprintf(tableLiq, ",%20.13e", (vLiq+totalVolume)*10.0);
    fprintf(tableLiq, ",%20.13e", (cpLiq+totalHeatCapacity));
    fprintf(tableLiq, ",%20.13e", (dvdtLiq+totaldVolumeDt)*10.0);
    fprintf(tableLiq, ",%20.13e", (dvdpLiq+totaldVolumeDp)*10.0);
    fprintf(tableLiq, ",%20.13e", ((vLiq+totalVolume) != 0.0) ?  (dvdtLiq+totaldVolumeDt)/(vLiq+totalVolume) : 0.0);
    fprintf(tableLiq, ",%20.13e", ((vLiq+totalVolume) != 0.0) ? -(dvdpLiq+totaldVolumeDp)/(vLiq+totalVolume) : 0.0);
    fprintf(tableLiq, ",%20.13e", dvdtLiq*10.0);
    fprintf(tableLiq, ",%20.13e", dvdpLiq*10.0);
    fprintf(tableLiq, ",%20.13e", (vLiq != 0.0) ?  dvdtLiq/vLiq : 0.0);
    fprintf(tableLiq, ",%20.13e", (vLiq != 0.0) ? -dvdpLiq/vLiq : 0.0);
    fprintf(tableLiq, "\n");
    }
#endif

    if (silminState->fo2Path != FO2_NONE) {
    double mO2 = -silminState->oxygen;
    for (nl=0; nl<silminState->nLiquidCoexist; nl++) for (i=0; i<nlc; i++) mO2 += (oxygen.liqToOx)[i]*(silminState->liquidComp)[nl][i];
    for (i=0; i<npc; i++) for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
        if (solids[i].na == 1) mO2 += (oxygen.solToOx)[i]*(silminState->solidComp)[i][ns];
        else {
        for (j=0; j<solids[i].na; j++) mO2 += (oxygen.solToOx)[i+1+j]*(silminState->solidComp)[i+1+j][ns];
        }
    }
    fprintf(output, "\n%-15.15s  delta moles = %.6g  delta grams = %.6g\n", "Oxygen", mO2, mO2*31.9988);
    fprintf(output, "         ");
    fprintf(output, "G = %.2f (J)  ",    mO2*(oxygen.cur).g);
    fprintf(output, "H = %.2f (J)  ",    mO2*(oxygen.cur).h);
    fprintf(output, "S = %.2f (J/K)  ",  mO2*(oxygen.cur).s);
    fprintf(output, "V = %.2f (cc)  ",   mO2*10.0*(oxygen.cur).v);
    fprintf(output, "Cp = %.2f (J/K)  ", mO2*(oxygen.cur).cp);
    fprintf(output, "\n");
    }

    if (silminState->assimilate) {
    fprintf(output, "\nSummary of assimilant: ");
    fprintf(output, "(total mass = %.2f grams, temperature = %.2f)\n", silminState->assimMass, silminState->assimT);

    for (j=0; j<npc; j++) if (solids[j].type == PHASE)
    for (ns=0; ns<(silminState->nAssimComp)[j]; ns++) {
        if (solids[j].na == 1) {
        mass = (silminState->assimComp)[j][ns]*solids[j].mw*silminState->assimMass/silminState->dspAssimMass;
        fprintf(output, "\n%-15.15s  mass = %.2f (gm)", solids[j].label, mass);
        fprintf(output, "     (analysis in mole %%)\n");
        } else {
        for (i=0, mass=0.0; i<solids[j].na; i++) mass += (silminState->assimComp)[j+1+i][ns]*solids[j+1+i].mw;
        mass *= silminState->assimMass/silminState->dspAssimMass;
        fprintf(output, "\n%-15.15s  mass = %.2f (gm)", solids[j].label, mass);
        fprintf(output, "     (analysis in mole %%)\n");
        fprintf(output, "     ");
        for (i=0; i<solids[j].na; i++) fprintf(output, " %13.13s", solids[j+1+i].label);
        fprintf(output, "\n");
        fprintf(output, "     ");
        for (i=0; i<solids[j].na; i++) fprintf(output, " %13.2f", 100.0*(silminState->assimComp)[j+1+i][ns]/(silminState->assimComp)[j][ns]);
        fprintf(output, "\n");
        }
    }
    }

#ifndef BATCH_VERSION
    wprintf(statusEntries[STATUS_ADB_INDEX_STATUS].name, "Current state of the system recorded in file %s.\n", meltsEnviron.OUTPUT_FILE);
    return GET_INPUT_SUCCESS;
#elif defined(ALPHAMELTS_UPDATE_SYSTEM)
    if(fflush(NULL)) fprintf(stderr, "Error returned when attempting to flush output files.");
    printf("Current state of the system recorded in melts.out and .tbl files.\n");
    return TRUE;
#else
    if(!fflush(NULL)) fprintf(stderr, "Output files flushed. ");
    else fprintf(stderr, "Error returned when attempting to flush output files. ");
    fprintf(stderr, " Current state of the system recorded in file melts.out.\n");
    return TRUE;
#endif
}

#ifdef BATCH_VERSION
#ifndef LIBXML_STATIC
int putSequenceDataToXmlFile(int active) {
    static xmlTextWriterPtr writer;
    static char *sequenceFile;

    //size_t len = strlen(silminInputData.name) - 4;
    size_t len = strlen(silminInputData.name) - 6; // assume .melts
    char *outputFile = (char *) malloc((size_t) (len+9)*sizeof(char));

    int rc;
    time_t tp;
    char * cOut, *temporary = (char *) malloc((size_t) 40*sizeof(char));
    double gLiq = 0.0, hLiq = 0.0, sLiq = 0.0, vLiq = 0.0, cpLiq = 0.0, mLiq = 0.0, viscosity = 0.0;
    double totalMass=0.0, totalGibbsEnergy=0.0, totalEnthalpy=0.0, totalEntropy=0.0, totalVolume=0.0, totalHeatCapacity=0.0, fracMass=0.0;
    static double *m, *r, *oxVal;
    int i, j;
    double fo2Delta;

    if (m == NULL)     m = (double *) malloc((size_t)      nc*sizeof(double));
    if (r == NULL)     r = (double *) malloc((size_t) (nlc-1)*sizeof(double));
    if (oxVal == NULL) oxVal = (double *) malloc((size_t)      nc*sizeof(double));

    if (active == FALSE) {
    if (writer != NULL) {
        rc = xmlTextWriterEndElement(writer);
        rc = xmlTextWriterEndDocument(writer);
        xmlFreeTextWriter(writer);
        printf("Sequence file name is %s\n", sequenceFile);
    }
    free (outputFile);
    free (temporary);
    return active;
    }

#ifdef PHMELTS_ADJUSTMENTS
    /* Use fracMass, fracSComp and fracLComp only for most recent fractionation (zeroed in silmin.c)*/
    /* Zero previous entries - should only be necessary if have never fractionated */
    /* Frac options may have changed at menu (like MCS) */
    if (previousSilminState == NULL) {
        previousSilminState = allocSilminStatePointer();
        previousSilminState->nLiquidCoexist = 1;
    }
    else {
        previousSilminState->fractionateSol = silminState->fractionateSol; // check
        previousSilminState->fractionateFlu = silminState->fractionateFlu;
        previousSilminState->fractionateLiq = silminState->fractionateLiq;
    }
    if (previousSilminState->fractionateSol || previousSilminState->fractionateFlu) {
        if (previousSilminState->nFracCoexist == NULL) {
            previousSilminState->fracSComp    = (double **) calloc((unsigned) npc, sizeof(double *));
            previousSilminState->nFracCoexist = (int *) calloc((unsigned) npc, sizeof(int));
        }
    }
    if (previousSilminState->fractionateLiq && previousSilminState->fracLComp == (double *) NULL)
        previousSilminState->fracLComp = (double *) calloc((unsigned) nlc, sizeof(double));
#endif

    (void) strncpy(outputFile, silminInputData.name, len);
    (void) strcpy(&outputFile[len], "-out.xml");

    if (writer == NULL) {
    sequenceFile = (char *) malloc((size_t) (len+14)*sizeof(char));
    (void) strncpy(sequenceFile, silminInputData.name, len);
    (void) strcpy(&sequenceFile[len], "-sequence.xml");

    printf("Sequence file name is %s\n", sequenceFile);
    writer = xmlNewTextWriterFilename(sequenceFile, 0);
    rc = xmlTextWriterStartDocument(writer, NULL, "UTF-8", NULL);
    rc = xmlTextWriterStartElement(writer, BAD_CAST "MELTSsequence");
    }

    (void) time(&tp);
    cOut = ctime(&tp);
    len = strlen(cOut);
    cOut[len-1] = '\0';

    rc = xmlTextWriterStartElement(writer, BAD_CAST "MELTSoutput");

    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "inputFile", "%s", silminInputData.name);
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "outputFile","%s", outputFile);
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "title",     "%s", silminInputData.title);
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "time",      "%s", cOut);
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "release",   "%s", RELEASE);
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "buildDate", "%s", __DATE__);
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "buildTime", "%s", __TIME__);

    rc = sprintf(temporary, "%23.16e", silminState->T-273.15); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "temperature", "%s", temporary);
    rc = sprintf(temporary, "%23.16e", silminState->P);    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "pressure",    "%s", temporary);
    rc = sprintf(temporary, "%23.16e", silminState->fo2);      rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "log_fO2",     "%s", temporary);

    fo2Delta = silminState->fo2Delta; silminState->fo2Delta = 0;
    rc = sprintf(temporary, "%23.16e", silminState->fo2 - getlog10fo2(silminState->T, silminState->P, FO2_HM));  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "deltaHM",     "%s", temporary);
    rc = sprintf(temporary, "%23.16e", silminState->fo2 - getlog10fo2(silminState->T, silminState->P, FO2_NNO)); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "deltaNNO",    "%s", temporary);
    rc = sprintf(temporary, "%23.16e", silminState->fo2 - getlog10fo2(silminState->T, silminState->P, FO2_QFM)); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "deltaFMQ",    "%s", temporary);
    rc = sprintf(temporary, "%23.16e", silminState->fo2 - getlog10fo2(silminState->T, silminState->P, FO2_COH)); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "deltaCOH",    "%s", temporary);
    rc = sprintf(temporary, "%23.16e", silminState->fo2 - getlog10fo2(silminState->T, silminState->P, FO2_IW));  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "deltaIW",     "%s", temporary);
    silminState->fo2Delta = fo2Delta;

    if (silminState->liquidMass != 0.0) {
        int nl;

        for (nl=0; nl<silminState->nLiquidCoexist; nl++) {
            double moles, oxSum;
            double gibbsEnergy, enthalpy, entropy, volume, heatCapacity;

            rc = xmlTextWriterStartElement(writer, BAD_CAST "liquid");

            conLiq(SECOND, THIRD, silminState->T, silminState->P, NULL, silminState->liquidComp[nl], r, NULL, NULL, NULL, NULL);

            gmixLiq (FIRST, silminState->T, silminState->P, r, &gibbsEnergy,  NULL, NULL);
            hmixLiq (FIRST, silminState->T, silminState->P, r, &enthalpy,	  NULL);
            smixLiq (FIRST, silminState->T, silminState->P, r, &entropy,	  NULL, NULL, NULL);
            vmixLiq (FIRST, silminState->T, silminState->P, r, &volume,	  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
            cpmixLiq(FIRST, silminState->T, silminState->P, r, &heatCapacity, NULL, NULL);
            visLiq  (FIRST, silminState->T, silminState->P, r, &viscosity);

            for (i=0, moles=0.0; i<nlc; i++) moles +=  (silminState->liquidComp)[nl][i];
            gibbsEnergy  *= moles;
            enthalpy     *= moles;
            entropy      *= moles;
            volume	     *= moles;
            heatCapacity *= moles;

            for (i=0; i<nlc; i++) {
                gibbsEnergy  += (silminState->liquidComp)[nl][i]*(liquid[i].cur).g;
                enthalpy     += (silminState->liquidComp)[nl][i]*(liquid[i].cur).h;
                entropy      += (silminState->liquidComp)[nl][i]*(liquid[i].cur).s;
                volume       += (silminState->liquidComp)[nl][i]*(liquid[i].cur).v;
                heatCapacity += (silminState->liquidComp)[nl][i]*(liquid[i].cur).cp;
            }

            for (i=0, oxSum=0.0; i<nc; i++) {
                for (j=0, oxVal[i]=0.0; j<nlc; j++) oxVal[i] += (liquid[j].liqToOx)[i]*(silminState->liquidComp)[nl][j];
                oxVal[i] *= bulkSystem[i].mw;
                oxSum	   += oxVal[i];
            }
            if (oxSum != 0.0) for (i=0; i<nc; i++) oxVal[i] /= oxSum;

            gLiq += gibbsEnergy; hLiq += enthalpy; sLiq += entropy; vLiq += volume; cpLiq += heatCapacity; mLiq += oxSum;

            rc = sprintf(temporary, "%23.16e", oxSum);    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "mass",    "%s", temporary);
            rc = sprintf(temporary, "%23.16e", (volume == 0.0) ? 0.0 : oxSum/(10.0*volume)); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "density", "%s", temporary);
            rc = sprintf(temporary, "%23.16e", viscosity);    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "viscosity",       "%s", temporary);
            rc = sprintf(temporary, "%23.16e", gibbsEnergy);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "gibbsFreeEnergy", "%s", temporary);
            rc = sprintf(temporary, "%23.16e", enthalpy);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "enthalpy",    "%s", temporary);
            rc = sprintf(temporary, "%23.16e", entropy);      rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "entropy",     "%s", temporary);
            rc = sprintf(temporary, "%23.16e", volume*10.0);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "volume",      "%s", temporary);
            rc = sprintf(temporary, "%23.16e", heatCapacity); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "heatCapacity",    "%s", temporary);

            for (i=0; i<nc; i++) if (oxVal[i] != 0.0) {
                rc = sprintf(temporary, "%23.16e", oxVal[i]*100.0); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST bulkSystem[i].label, "%s", temporary);
            }
            rc = xmlTextWriterEndElement(writer); //liquid
        }
    }

    for (j=0, totalMass=0.0, totalGibbsEnergy=0.0, totalEnthalpy=0.0,
        totalEntropy=0.0, totalVolume=0.0, totalHeatCapacity=0.0; j<npc; j++) {
        int ns;
        for (ns=0; ns<(silminState->nSolidCoexist)[j]; ns++) {
        double oxSum, mass, gibbsEnergy, enthalpy, entropy, volume, heatCapacity;

        rc = xmlTextWriterStartElement(writer, BAD_CAST "solid");

        if (solids[j].na == 1) {
            mass  	     = (silminState->solidComp)[j][ns]*solids[j].mw;
            gibbsEnergy	     = (silminState->solidComp)[j][ns]*(solids[j].cur).g;
            enthalpy	     = (silminState->solidComp)[j][ns]*(solids[j].cur).h;
            entropy	     = (silminState->solidComp)[j][ns]*(solids[j].cur).s;
            volume	     = (silminState->solidComp)[j][ns]*(solids[j].cur).v;
            heatCapacity       = (silminState->solidComp)[j][ns]*(solids[j].cur).cp;
            totalMass	    += (silminState->solidComp)[j][ns]*solids[j].mw;
            totalGibbsEnergy  += (silminState->solidComp)[j][ns]*(solids[j].cur).g;
            totalEnthalpy     += (silminState->solidComp)[j][ns]*(solids[j].cur).h;
            totalEntropy      += (silminState->solidComp)[j][ns]*(solids[j].cur).s;
            totalVolume	    += (silminState->solidComp)[j][ns]*(solids[j].cur).v;
            totalHeatCapacity += (silminState->solidComp)[j][ns]*(solids[j].cur).cp;

            rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "name",    "%s",   solids[j].label);
            rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "formula",     "%s",   solids[j].formula);
            rc = sprintf(temporary, "%23.16e", mass);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "mass", 	   "%s", temporary);
            rc = sprintf(temporary, "%23.16e", (volume == 0.0) ? 0.0 : mass/(10.0*volume)); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "density", "%s", temporary);
            rc = sprintf(temporary, "%23.16e", gibbsEnergy);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "gibbsFreeEnergy", "%s", temporary);
            rc = sprintf(temporary, "%23.16e", enthalpy);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "enthalpy",	     "%s", temporary);
            rc = sprintf(temporary, "%23.16e", entropy);      rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "entropy",	     "%s", temporary);
            rc = sprintf(temporary, "%23.16e", volume*10.0);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "volume",	     "%s", temporary);
            rc = sprintf(temporary, "%23.16e", heatCapacity); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "heatCapacity",    "%s", temporary);

            for (i=0, oxSum=0.0; i<nc; i++) {
                oxVal[i]  = (solids[j].solToOx)[i]*bulkSystem[i].mw;
                oxSum    += oxVal[i];
            }
            if (oxSum != 0.0) for (i=0; i<nc; i++) oxVal[i] /= oxSum;
            for (i=0; i<nc; i++) if (oxVal[i] != 0.0) {
                rc = sprintf(temporary, "%23.16e", oxVal[i]*100.0); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST bulkSystem[i].label, "%s", temporary);
            }
            rc = xmlTextWriterStartElement(writer, BAD_CAST "component");
            rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "name",     "%s",   solids[j].label);
            rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "formula",      "%s",   solids[j].formula);
            rc = sprintf(temporary, "%23.16e", 1.0); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "moleFraction", "%s", temporary);
            rc = xmlTextWriterEndElement(writer);

        } else {
            char *formula;
            for (i=0, mass=0.0; i<solids[j].na; i++) {
                m[i] = (silminState->solidComp)[j+1+i][ns];
                mass += (silminState->solidComp)[j+1+i][ns]*solids[j+1+i].mw;
            }
            (*solids[j].convert)(SECOND, THIRD, silminState->T, silminState->P, NULL, m, r, NULL, NULL, NULL, NULL, NULL);
            (*solids[j].display)(FIRST, silminState->T, silminState->P, r, &formula);
            (*solids[j].gmix) (FIRST, silminState->T, silminState->P, r, &gibbsEnergy,  NULL, NULL, NULL);
            (*solids[j].hmix) (FIRST, silminState->T, silminState->P, r, &enthalpy);
            (*solids[j].smix) (FIRST, silminState->T, silminState->P, r, &entropy,      NULL, NULL);
            (*solids[j].vmix) (FIRST, silminState->T, silminState->P, r, &volume,       NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL);
            (*solids[j].cpmix)(FIRST, silminState->T, silminState->P, r, &heatCapacity, NULL, NULL);
            gibbsEnergy  *= (silminState->solidComp)[j][ns];
            enthalpy     *= (silminState->solidComp)[j][ns];
            entropy      *= (silminState->solidComp)[j][ns];
            volume       *= (silminState->solidComp)[j][ns];
            heatCapacity *= (silminState->solidComp)[j][ns];
            for (i=0; i<solids[j].na; i++) {
                gibbsEnergy  += m[i]*(solids[j+1+i].cur).g;
                enthalpy	 += m[i]*(solids[j+1+i].cur).h;
                entropy	 += m[i]*(solids[j+1+i].cur).s;
                volume	 += m[i]*(solids[j+1+i].cur).v;
                heatCapacity += m[i]*(solids[j+1+i].cur).cp;
            }
            totalMass	    += mass;
            totalGibbsEnergy  += gibbsEnergy;
            totalEnthalpy     += enthalpy;
            totalEntropy      += entropy;
            totalVolume	    += volume;
            totalHeatCapacity += heatCapacity;

            rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "name",    "%s",   solids[j].label);
            rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "formula",     "%s",   formula); free(formula);
            rc = sprintf(temporary, "%23.16e", mass);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "mass", 	   "%s", temporary);
            rc = sprintf(temporary, "%23.16e", (volume == 0.0) ? 0.0 : mass/(10.0*volume)); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "density", "%s", temporary);
            rc = sprintf(temporary, "%23.16e", gibbsEnergy);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "gibbsFreeEnergy", "%s", temporary);
            rc = sprintf(temporary, "%23.16e", enthalpy);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "enthalpy",	     "%s", temporary);
            rc = sprintf(temporary, "%23.16e", entropy);      rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "entropy", 	     "%s", temporary);
            rc = sprintf(temporary, "%23.16e", volume*10.0);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "volume",	     "%s", temporary);
            rc = sprintf(temporary, "%23.16e", heatCapacity); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "heatCapacity",    "%s", temporary);

            for (i=0, oxSum=0.0; i<nc; i++) {
                int k;
                for (k=0, oxVal[i]=0.0; k<solids[j].na; k++) oxVal[i] += (solids[j+1+k].solToOx)[i]*m[k]*bulkSystem[i].mw;
                oxSum += oxVal[i];
            }
            if (oxSum != 0.0) for (i=0; i<nc; i++) oxVal[i] /= oxSum;
            for (i=0; i<nc; i++) if (oxVal[i] != 0.0) {
                rc = sprintf(temporary, "%23.16e", oxVal[i]*100.0); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST bulkSystem[i].label, "%s", temporary);
            }

            for (i=0; i<solids[j].na; i++) {
                rc = xmlTextWriterStartElement(writer, BAD_CAST "component");
                rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "name",     "%s",   solids[j+1+i].label);
                rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "formula",      "%s",   solids[j+1+i].formula);
                rc = sprintf(temporary, "%23.16e", m[i]/(silminState->solidComp)[j][ns]); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "moleFraction", "%s", temporary);
                rc = xmlTextWriterEndElement(writer);
            }
        }

        rc = xmlTextWriterEndElement(writer); //totalSolids
        }
    }

    if (totalMass != 0.0) {
        rc = xmlTextWriterStartElement(writer, BAD_CAST "totalSolids");
        rc = sprintf(temporary, "%23.16e", totalMass);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "mass",		      "%s", temporary);
        rc = sprintf(temporary, "%23.16e", (totalVolume == 0.0) ? 0.0 : totalMass/(10.0*totalVolume)); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "density", "%s", temporary);
        rc = sprintf(temporary, "%23.16e", totalGibbsEnergy);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "gibbsFreeEnergy", "%s", temporary);
        rc = sprintf(temporary, "%23.16e", totalEnthalpy);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "enthalpy",    "%s", temporary);
        rc = sprintf(temporary, "%23.16e", totalEntropy);      rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "entropy",     "%s", temporary);
        rc = sprintf(temporary, "%23.16e", totalVolume*10.0);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "volume",      "%s", temporary);
        rc = sprintf(temporary, "%23.16e", totalHeatCapacity); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "heatCapacity",	  "%s", temporary);
        rc = xmlTextWriterEndElement(writer); //totalSolids
    }

    if (silminState->isenthalpic && (silminState->refEnthalpy == 0.0)) silminState->refEnthalpy = hLiq+totalEnthalpy;
    if (silminState->isentropic  && (silminState->refEntropy  == 0.0)) silminState->refEntropy  = sLiq+totalEntropy;
    if (silminState->isochoric   && (silminState->refVolume   == 0.0)) silminState->refVolume   = vLiq+totalVolume;

#ifdef PHMELTS_ADJUSTMENTS
    /* Web services uses "previous"; Melts-batch uses "current" and subtracts off previous; here "previous" has already been zerod out */
    fracMass = (silminState->fractionateSol || silminState->fractionateFlu || silminState->fractionateLiq) ? (silminState->fracMass) : 0.0;

    if ((silminState->fractionateSol || silminState->fractionateFlu || silminState->fractionateLiq) && (fracMass > 0.0)) {
        rc = xmlTextWriterStartElement(writer, BAD_CAST "fractionate");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "type", BAD_CAST "current");
        rc = sprintf(temporary, "%23.16e", silminState->fracMass); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "mass", "%s", temporary);

        if (silminState->fractionateSol || silminState->fractionateFlu) {
        for (j=0; j<npc; j++) {
            int ns;

            if ( silminState->fractionateSol && !silminState->fractionateFlu && !strcmp((char *) solids[j].label, "fluid")) continue;
            if (!silminState->fractionateSol &&  silminState->fractionateFlu &&  strcmp((char *) solids[j].label, "fluid")) continue;

            for (ns=0; ns<(silminState->nFracCoexist)[j]; ns++) {
                double oxSum, mass, gibbsEnergy, enthalpy, entropy, volume, heatCapacity;
                double tmpMoles = (silminState->fracSComp)[j][ns];
                if (fabs(tmpMoles) < 10.0*DBL_EPSILON) continue;

                rc = xmlTextWriterStartElement(writer, BAD_CAST "solid");

                if (solids[j].na == 1) {
                mass		   = tmpMoles*solids[j].mw;
                gibbsEnergy	   = tmpMoles*(solids[j].cur).g;
                enthalpy	   = tmpMoles*(solids[j].cur).h;
                entropy 	   = tmpMoles*(solids[j].cur).s;
                volume  	   = tmpMoles*(solids[j].cur).v;
                heatCapacity	   = tmpMoles*(solids[j].cur).cp;
                totalMass	  += tmpMoles*solids[j].mw;
                totalGibbsEnergy  += tmpMoles*(solids[j].cur).g;
                totalEnthalpy	  += tmpMoles*(solids[j].cur).h;
                totalEntropy	  += tmpMoles*(solids[j].cur).s;
                totalVolume	  += tmpMoles*(solids[j].cur).v;
                totalHeatCapacity += tmpMoles*(solids[j].cur).cp;

                rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "name",		 "%s",   solids[j].label);
                rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "formula",	 "%s",   solids[j].formula);
                rc = sprintf(temporary, "%23.16e", mass);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "mass",	     "%s", temporary);
                rc = sprintf(temporary, "%23.16e", (volume == 0.0) ? 0.0 : mass/(10.0*volume)); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "density", "%s", temporary);
                rc = sprintf(temporary, "%23.16e", gibbsEnergy);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "gibbsFreeEnergy", "%s", temporary);
                rc = sprintf(temporary, "%23.16e", enthalpy);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "enthalpy",	     "%s", temporary);
                rc = sprintf(temporary, "%23.16e", entropy);      rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "entropy",	     "%s", temporary);
                rc = sprintf(temporary, "%23.16e", volume*10.0);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "volume",      "%s", temporary);
                rc = sprintf(temporary, "%23.16e", heatCapacity); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "heatCapacity",	 "%s", temporary);

                for (i=0, oxSum=0.0; i<nc; i++) {
                    oxVal[i]  = (solids[j].solToOx)[i]*bulkSystem[i].mw;
                    oxSum    += oxVal[i];
                }
                        if (oxSum != 0.0) for (i=0; i<nc; i++) oxVal[i] /= oxSum;
                        for (i=0; i<nc; i++) if (oxVal[i] != 0.0) {
                            rc = sprintf(temporary, "%23.16e", oxVal[i]*100.0); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST bulkSystem[i].label, "%s", temporary);
                        }
                        rc = xmlTextWriterStartElement(writer, BAD_CAST "component");
                        rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "name",     "%s",   solids[j].label);
                        rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "formula",      "%s",   solids[j].formula);
                        rc = sprintf(temporary, "%23.16e", 1.0); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "moleFraction", "%s", temporary);
                        rc = xmlTextWriterEndElement(writer);

                    } else {
                        char *formula;
                        for (i=0, mass=0.0; i<solids[j].na; i++) {
                            double tmpCmpMoles = (silminState->fracSComp)[j+1+i][ns];
                            m[i]  = tmpCmpMoles;
                            mass += tmpCmpMoles*solids[j+1+i].mw;
                        }
                        (*solids[j].convert)(SECOND, THIRD, silminState->T, silminState->P, NULL, m, r, NULL, NULL, NULL, NULL, NULL);
                        (*solids[j].display)(FIRST, silminState->T, silminState->P, r, &formula);
                        (*solids[j].gmix)   (FIRST, silminState->T, silminState->P, r, &gibbsEnergy,  NULL, NULL, NULL);
                        (*solids[j].hmix)   (FIRST, silminState->T, silminState->P, r, &enthalpy);
                        (*solids[j].smix)   (FIRST, silminState->T, silminState->P, r, &entropy,      NULL, NULL);
                        (*solids[j].vmix)   (FIRST, silminState->T, silminState->P, r, &volume,	      NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL);
                        (*solids[j].cpmix)  (FIRST, silminState->T, silminState->P, r, &heatCapacity, NULL, NULL);
                        gibbsEnergy  *= tmpMoles;
                        enthalpy     *= tmpMoles;
                        entropy      *= tmpMoles;
                        volume       *= tmpMoles;
                        heatCapacity *= tmpMoles;
                        for (i=0; i<solids[j].na; i++) {
                            gibbsEnergy  += m[i]*(solids[j+1+i].cur).g;
                            enthalpy     += m[i]*(solids[j+1+i].cur).h;
                            entropy      += m[i]*(solids[j+1+i].cur).s;
                            volume       += m[i]*(solids[j+1+i].cur).v;
                            heatCapacity += m[i]*(solids[j+1+i].cur).cp;
                        }
                        totalMass	  += mass;
                        totalGibbsEnergy  += gibbsEnergy;
                        totalEnthalpy	  += enthalpy;
                        totalEntropy	  += entropy;
                        totalVolume	  += volume;
                        totalHeatCapacity += heatCapacity;

                        rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "name",		 "%s",   solids[j].label);
                        rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "formula",	 "%s",   formula); free(formula);
                        rc = sprintf(temporary, "%23.16e", mass);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "mass",		     "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", (volume == 0.0) ? 0.0 : mass/(10.0*volume)); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "density", "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", gibbsEnergy);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "gibbsFreeEnergy", "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", enthalpy);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "enthalpy",    "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", entropy);      rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "entropy",     "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", volume*10.0);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "volume",      "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", heatCapacity); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "heatCapacity",	 "%s", temporary);

                        for (i=0, oxSum=0.0; i<nc; i++) {
                            int k;
                            for (k=0, oxVal[i]=0.0; k<solids[j].na; k++) oxVal[i] += (solids[j+1+k].solToOx)[i]*m[k]*bulkSystem[i].mw;
                            oxSum += oxVal[i];
                        }
                        if (oxSum != 0.0) for (i=0; i<nc; i++) oxVal[i] /= oxSum;
                        for (i=0; i<nc; i++) if (oxVal[i] != 0.0) {
                            rc = sprintf(temporary, "%23.16e", oxVal[i]*100.0); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST bulkSystem[i].label, "%s", temporary);
                        }
                        for (i=0; i<solids[j].na; i++) {
                            rc = xmlTextWriterStartElement(writer, BAD_CAST "component");
                            rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "name", 	  "%s",   solids[j+1+i].label);
                            rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "formula",  "%s",   solids[j+1+i].formula);
                            rc = sprintf(temporary, "%23.16e", m[i]/tmpMoles);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "moleFraction", "%s", temporary);
                            rc = xmlTextWriterEndElement(writer);
                        }
                    }

                    rc = xmlTextWriterEndElement(writer); //solid
                }
            }
        }

        if (silminState->fractionateLiq) {
            double oxSum, mass, moles, gibbsEnergy, enthalpy, entropy, volume, heatCapacity;
            char *formula;

            for (i=0, mass=0.0, moles=0.0; i<nlc; i++) {
                double mw;
                double tmpMoles = (silminState->fracLComp)[i];
                for (j=0, mw = 0.0; j<nc; j++) mw += (liquid[i].liqToOx)[j]*bulkSystem[j].mw;
                m[i]   = tmpMoles;
                moles += m[i];
                mass  += tmpMoles*mw;
            }

            if (mass > 0.0) {
                rc = xmlTextWriterStartElement(writer, BAD_CAST "liquid");

                conLiq  (SECOND, THIRD, silminState->T, silminState->P, NULL, m, r, NULL, NULL, NULL, NULL);
                dispLiq (FIRST, silminState->T, silminState->P, r, &formula);
                gmixLiq (FIRST, silminState->T, silminState->P, r, &gibbsEnergy,  NULL, NULL);
                hmixLiq (FIRST, silminState->T, silminState->P, r, &enthalpy,     NULL);
                smixLiq (FIRST, silminState->T, silminState->P, r, &entropy,      NULL, NULL, NULL);
                vmixLiq (FIRST, silminState->T, silminState->P, r, &volume,       NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL);
                cpmixLiq(FIRST, silminState->T, silminState->P, r, &heatCapacity, NULL, NULL);
                gibbsEnergy  *= moles;
                enthalpy	 *= moles;
                entropy	 *= moles;
                volume	 *= moles;
                heatCapacity *= moles;
                for (i=0; i<nlc; i++) {
                    gibbsEnergy  += m[i]*(liquid[i].cur).g;
                    enthalpy     += m[i]*(liquid[i].cur).h;
                    entropy	   += m[i]*(liquid[i].cur).s;
                    volume	   += m[i]*(liquid[i].cur).v;
                    heatCapacity += m[i]*(liquid[i].cur).cp;
                }
                totalMass	      += mass;
                totalGibbsEnergy  += gibbsEnergy;
                totalEnthalpy     += enthalpy;
                totalEntropy      += entropy;
                totalVolume       += volume;
                totalHeatCapacity += heatCapacity;

                for (i=0, oxSum=0.0; i<nc; i++) {
                    for (j=0, oxVal[i]=0.0; j<nlc; j++) oxVal[i] += (liquid[j].liqToOx)[i]*m[j];
                    oxVal[i] *= bulkSystem[i].mw;
                    oxSum    += oxVal[i];
                }
                if (oxSum != 0.0) for (i=0; i<nc; i++) oxVal[i] /= oxSum;

                rc = sprintf(temporary, "%23.16e", mass);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "mass",	     "%s", temporary);
                rc = sprintf(temporary, "%23.16e", (volume == 0.0) ? 0.0 : mass/(10.0*volume)); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "density", "%s", temporary);
                rc = sprintf(temporary, "%23.16e", gibbsEnergy);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "gibbsFreeEnergy", "%s", temporary);
                rc = sprintf(temporary, "%23.16e", enthalpy);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "enthalpy",	     "%s", temporary);
                rc = sprintf(temporary, "%23.16e", entropy);      rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "entropy",	     "%s", temporary);
                rc = sprintf(temporary, "%23.16e", volume*10.0);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "volume",	     "%s", temporary);
                rc = sprintf(temporary, "%23.16e", heatCapacity); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "heatCapacity",    "%s", temporary);

                for (i=0; i<nc; i++) if (oxVal[i] != 0.0) {
                    rc = sprintf(temporary, "%23.16e", oxVal[i]*100.0); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST bulkSystem[i].label, "%s", temporary);
                }
                rc = xmlTextWriterEndElement(writer); //liquid
            }
        }

        rc = xmlTextWriterEndElement(writer); //fractionate
    }
#endif
    /* Web services uses "previous"; Melts-batch uses "current" and subtracts off previous */
    if ( (previousSilminState->fractionateSol || previousSilminState->fractionateFlu || previousSilminState->fractionateLiq)
        && (previousSilminState->fracMass > 0.0) ) {
        rc = xmlTextWriterStartElement(writer, BAD_CAST "fractionate");
        rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "type", BAD_CAST "previous");
        rc = sprintf(temporary, "%23.16e", previousSilminState->fracMass); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "mass", "%s", temporary);

        if (previousSilminState->fractionateSol || previousSilminState->fractionateFlu) {
            for (j=0; j<npc; j++) {
                int ns;

                if ( previousSilminState->fractionateSol && !previousSilminState->fractionateFlu && !strcmp((char *) solids[j].label, "fluid")) continue;
                if (!previousSilminState->fractionateSol &&  previousSilminState->fractionateFlu &&  strcmp((char *) solids[j].label, "fluid")) continue;

                for (ns=0; ns<(previousSilminState->nFracCoexist)[j]; ns++) {
                    double oxSum, mass, gibbsEnergy, enthalpy, entropy, volume, heatCapacity;

                    rc = xmlTextWriterStartElement(writer, BAD_CAST "solid");

                    if (solids[j].na == 1) {
                        mass		   = (previousSilminState->fracSComp)[j][ns]*solids[j].mw;
                        gibbsEnergy	   = (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).g;
                        enthalpy	   = (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).h;
                        entropy 	   = (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).s;
                        volume  	   = (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).v;
                        heatCapacity	   = (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).cp;
                        totalMass	  += (previousSilminState->fracSComp)[j][ns]*solids[j].mw;
                        totalGibbsEnergy  += (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).g;
                        totalEnthalpy	  += (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).h;
                        totalEntropy	  += (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).s;
                        totalVolume	  += (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).v;
                        totalHeatCapacity += (previousSilminState->fracSComp)[j][ns]*(solids[j].cur).cp;

                        rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "name",		 "%s",   solids[j].label);
                        rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "formula",	 "%s",   solids[j].formula);
                        rc = sprintf(temporary, "%23.16e", mass);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "mass",	     "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", (volume == 0.0) ? 0.0 : mass/(10.0*volume)); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "density", "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", gibbsEnergy);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "gibbsFreeEnergy", "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", enthalpy);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "enthalpy",    "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", entropy);      rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "entropy",     "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", volume*10.0);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "volume",      "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", heatCapacity); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "heatCapacity",	 "%s", temporary);

                        for (i=0, oxSum=0.0; i<nc; i++) {
                            oxVal[i]  = (solids[j].solToOx)[i]*bulkSystem[i].mw;
                            oxSum    += oxVal[i];
                        }
                        if (oxSum != 0.0) for (i=0; i<nc; i++) oxVal[i] /= oxSum;
                        for (i=0; i<nc; i++) if (oxVal[i] != 0.0) {
                            rc = sprintf(temporary, "%23.16e", oxVal[i]*100.0); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST bulkSystem[i].label, "%s", temporary);
                        }
                        rc = xmlTextWriterStartElement(writer, BAD_CAST "component");
                        rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "name",     "%s",   solids[j].label);
                        rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "formula",      "%s",   solids[j].formula);
                        rc = sprintf(temporary, "%23.16e", 1.0); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "moleFraction", "%s", temporary);
                        rc = xmlTextWriterEndElement(writer);

                    } else {
                        char *formula;
                        for (i=0, mass=0.0; i<solids[j].na; i++) {
                            m[i] = (previousSilminState->fracSComp)[j+1+i][ns];
                            mass += (previousSilminState->fracSComp)[j+1+i][ns]*solids[j+1+i].mw;
                        }
                        (*solids[j].convert)(SECOND, THIRD, silminState->T, silminState->P, NULL, m, r, NULL, NULL, NULL, NULL, NULL);
                        (*solids[j].display)(FIRST, silminState->T, silminState->P, r, &formula);
                        (*solids[j].gmix)   (FIRST, silminState->T, silminState->P, r, &gibbsEnergy,  NULL, NULL, NULL);
                        (*solids[j].hmix)   (FIRST, silminState->T, silminState->P, r, &enthalpy);
                        (*solids[j].smix)   (FIRST, silminState->T, silminState->P, r, &entropy,      NULL, NULL);
                        (*solids[j].vmix)   (FIRST, silminState->T, silminState->P, r, &volume,	      NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL);
                        (*solids[j].cpmix)  (FIRST, silminState->T, silminState->P, r, &heatCapacity, NULL, NULL);
                        gibbsEnergy  *= (previousSilminState->fracSComp)[j][ns];
                        enthalpy     *= (previousSilminState->fracSComp)[j][ns];
                        entropy      *= (previousSilminState->fracSComp)[j][ns];
                        volume       *= (previousSilminState->fracSComp)[j][ns];
                        heatCapacity *= (previousSilminState->fracSComp)[j][ns];
                        for (i=0; i<solids[j].na; i++) {
                            gibbsEnergy  += m[i]*(solids[j+1+i].cur).g;
                            enthalpy     += m[i]*(solids[j+1+i].cur).h;
                            entropy      += m[i]*(solids[j+1+i].cur).s;
                            volume       += m[i]*(solids[j+1+i].cur).v;
                            heatCapacity += m[i]*(solids[j+1+i].cur).cp;
                        }
                        totalMass	  += mass;
                        totalGibbsEnergy  += gibbsEnergy;
                        totalEnthalpy	  += enthalpy;
                        totalEntropy	  += entropy;
                        totalVolume	  += volume;
                        totalHeatCapacity += heatCapacity;

                        rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "name",		 "%s",   solids[j].label);
                        rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "formula",	 "%s",   formula); free(formula);
                        rc = sprintf(temporary, "%23.16e", mass);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "mass",	     "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", (volume == 0.0) ? 0.0 : mass/(10.0*volume)); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "density", "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", gibbsEnergy);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "gibbsFreeEnergy", "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", enthalpy);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "enthalpy",    "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", entropy);      rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "entropy",     "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", volume*10.0);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "volume",      "%s", temporary);
                        rc = sprintf(temporary, "%23.16e", heatCapacity); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "heatCapacity",	 "%s", temporary);

                        for (i=0, oxSum=0.0; i<nc; i++) {
                            int k;
                            for (k=0, oxVal[i]=0.0; k<solids[j].na; k++) oxVal[i] += (solids[j+1+k].solToOx)[i]*m[k]*bulkSystem[i].mw;
                            oxSum += oxVal[i];
                        }
                        if (oxSum != 0.0) for (i=0; i<nc; i++) oxVal[i] /= oxSum;
                        for (i=0; i<nc; i++) if (oxVal[i] != 0.0) {
                            rc =  sprintf(temporary, "%23.16e", oxVal[i]*100.0); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST bulkSystem[i].label, "%s", temporary);
                        }
                        for (i=0; i<solids[j].na; i++) {
                            rc = xmlTextWriterStartElement(writer, BAD_CAST "component");
                            rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "name", 	  "%s",   solids[j+1+i].label);
                            rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "formula",	  "%s",   solids[j+1+i].formula);
                            rc = sprintf(temporary, "%23.16e", m[i]/(previousSilminState->fracSComp)[j][ns]); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "moleFraction", "%s", temporary);
                            rc = xmlTextWriterEndElement(writer);
                        }
                    }

                    rc = xmlTextWriterEndElement(writer); //solid
                }
            }
        }

        if (previousSilminState->fractionateLiq) {
            double oxSum, mass, moles, gibbsEnergy, enthalpy, entropy, volume, heatCapacity;
            char *formula;

            for (i=0, mass=0.0, moles=0.0; i<nlc; i++) {
                double mw;
                for (j=0, mw = 0.0; j<nc; j++) mw += (liquid[i].liqToOx)[j]*bulkSystem[j].mw;
                m[i]   = (previousSilminState->fracLComp)[i];
                moles += m[i];
                mass  += (previousSilminState->fracLComp)[i]*mw;
            }

            if (mass > 0.0) {
                rc = xmlTextWriterStartElement(writer, BAD_CAST "liquid");

                conLiq  (SECOND, THIRD, silminState->T, silminState->P, NULL, m, r, NULL, NULL, NULL, NULL);
                dispLiq (FIRST, silminState->T, silminState->P, r, &formula);
                gmixLiq (FIRST, silminState->T, silminState->P, r, &gibbsEnergy,  NULL, NULL);
                hmixLiq (FIRST, silminState->T, silminState->P, r, &enthalpy,     NULL);
                smixLiq (FIRST, silminState->T, silminState->P, r, &entropy,      NULL, NULL, NULL);
                vmixLiq (FIRST, silminState->T, silminState->P, r, &volume,       NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL);
                cpmixLiq(FIRST, silminState->T, silminState->P, r, &heatCapacity, NULL, NULL);
                gibbsEnergy  *= moles;
                enthalpy	 *= moles;
                entropy	 *= moles;
                volume	 *= moles;
                heatCapacity *= moles;
                for (i=0; i<nlc; i++) {
                    gibbsEnergy  += m[i]*(liquid[i].cur).g;
                    enthalpy     += m[i]*(liquid[i].cur).h;
                    entropy	   += m[i]*(liquid[i].cur).s;
                    volume	   += m[i]*(liquid[i].cur).v;
                    heatCapacity += m[i]*(liquid[i].cur).cp;
                }
                totalMass	      += mass;
                totalGibbsEnergy  += gibbsEnergy;
                totalEnthalpy     += enthalpy;
                totalEntropy      += entropy;
                totalVolume       += volume;
                totalHeatCapacity += heatCapacity;

                for (i=0, oxSum=0.0; i<nc; i++) {
                    for (j=0, oxVal[i]=0.0; j<nlc; j++) oxVal[i] += (liquid[j].liqToOx)[i]*(previousSilminState->fracLComp)[j];
                    oxVal[i] *= bulkSystem[i].mw;
                    oxSum    += oxVal[i];
                }
                if (oxSum != 0.0) for (i=0; i<nc; i++) oxVal[i] /= oxSum;

                rc = sprintf(temporary, "%23.16e", mass);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "mass",       "%s", temporary);
                rc = sprintf(temporary, "%23.16e", (volume == 0.0) ? 0.0 : mass/(10.0*volume)); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "density", "%s", temporary);
                rc = sprintf(temporary, "%23.16e", gibbsEnergy);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "gibbsFreeEnergy", "%s", temporary);
                rc = sprintf(temporary, "%23.16e", enthalpy);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "enthalpy",	     "%s", temporary);
                rc = sprintf(temporary, "%23.16e", entropy);      rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "entropy",	     "%s", temporary);
                rc = sprintf(temporary, "%23.16e", volume*10.0);  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "volume",	     "%s", temporary);
                rc = sprintf(temporary, "%23.16e", heatCapacity); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "heatCapacity",    "%s", temporary);

                for (i=0; i<nc; i++) if (oxVal[i] != 0.0) {
                    rc = sprintf(temporary, "%23.16e", oxVal[i]*100.0); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST bulkSystem[i].label, "%s", temporary);
                }
                rc = xmlTextWriterEndElement(writer); //liquid
            }
        }

        rc = xmlTextWriterEndElement(writer); //fractionate
    }

    rc = xmlTextWriterStartElement(writer, BAD_CAST "system");
    rc = sprintf(temporary, "%23.16e", mLiq+totalMass); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "mass", "%s", temporary);
    rc = sprintf(temporary, "%23.16e", (vLiq+totalVolume == 0.0) ? 0.0 : (mLiq+totalMass)/(10.0*(vLiq+totalVolume)));
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "density", "%s", temporary);

    if (vLiq > totalVolume) {
        rc = sprintf(temporary, "%23.16e", viscosity - 2.0*log10(1.0-2.0*totalVolume/(totalVolume+vLiq)));
        rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "viscosity", "%s", temporary);
    }
    rc = sprintf(temporary, "%23.16e", gLiq+totalGibbsEnergy);   rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "gibbsFreeEnergy", "%s", temporary);
    rc = sprintf(temporary, "%23.16e", hLiq+totalEnthalpy);      rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "enthalpy",    "%s", temporary);
    rc = sprintf(temporary, "%23.16e", sLiq+totalEntropy);       rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "entropy",     "%s", temporary);
    rc = sprintf(temporary, "%23.16e", (vLiq+totalVolume)*10.0); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "volume",	    "%s", temporary);
    rc = sprintf(temporary, "%23.16e", cpLiq+totalHeatCapacity); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "heatCapacity",    "%s", temporary);
    rc = xmlTextWriterEndElement(writer); //system

    if (silminState->fo2Path != FO2_NONE) {
        double mO2 = -silminState->oxygen;
        int nl, ns;
        for (nl=0; nl<silminState->nLiquidCoexist; nl++) for (i=0; i<nlc; i++) mO2 += (oxygen.liqToOx)[i]*(silminState->liquidComp)[nl][i];
        for (i=0; i<npc; i++) for (ns=0; ns<(silminState->nSolidCoexist)[i]; ns++) {
            if (solids[i].na == 1) mO2 += (oxygen.solToOx)[i]*(silminState->solidComp)[i][ns];
            else {
                for (j=0; j<solids[i].na; j++) mO2 += (oxygen.solToOx)[i+1+j]*(silminState->solidComp)[i+1+j][ns];
            }
        }
        rc = xmlTextWriterStartElement(writer, BAD_CAST "oxygen");
        rc = sprintf(temporary, "%23.16e", mO2);             rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "moles",	    "%s", temporary);
        rc = sprintf(temporary, "%23.16e", mO2*31.9988);         rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "mass",	    "%s", temporary);
        rc = sprintf(temporary, "%23.16e", mO2*(oxygen.cur).g);      rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "gibbsFreeEnergy", "%s", temporary);
        rc = sprintf(temporary, "%23.16e", mO2*(oxygen.cur).h);      rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "enthalpy",    "%s", temporary);
        rc = sprintf(temporary, "%23.16e", mO2*(oxygen.cur).s);      rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "entropy",     "%s", temporary);
        rc = sprintf(temporary, "%23.16e", mO2*10.0*(oxygen.cur).v); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "volume",	    "%s", temporary);
        rc = sprintf(temporary, "%23.16e", mO2*(oxygen.cur).cp);     rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "heatCapacity",    "%s", temporary);
        rc = xmlTextWriterEndElement(writer); //oxygen
    }

    if (silminState->assimilate) {
        int ns;
        rc = xmlTextWriterStartElement(writer, BAD_CAST "assimilant");
        rc = sprintf(temporary, "%23.16e", silminState->assimMass); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "mass",	    "%s", temporary);
        rc = sprintf(temporary, "%23.16e", silminState->assimT); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "temperature", "%s", temporary);

        for (j=0; j<npc; j++) if (solids[j].type == PHASE) {
            for (ns=0; ns<(silminState->nAssimComp)[j]; ns++) {
                rc = xmlTextWriterStartElement(writer, BAD_CAST "solid");
                rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "name", "%s", solids[j].label);
                if (solids[j].na == 1) {
                    double mass = (silminState->assimComp)[j][ns]*solids[j].mw*silminState->assimMass/silminState->dspAssimMass;
                    rc = sprintf(temporary, "%23.16e", mass); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "mass", "%s", temporary);
                    rc = xmlTextWriterStartElement(writer, BAD_CAST "component");
                    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "name",     "%s",   solids[j].label);
                    rc = sprintf(temporary, "%23.16e", 1.0); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "moleFraction", "%s", temporary);
                    rc = xmlTextWriterEndElement(writer);

                } else {
                    double mass = 0.0;
                    for (i=0; i<solids[j].na; i++) mass += (silminState->assimComp)[j+1+i][ns]*solids[j+1+i].mw;
                    mass *= silminState->assimMass/silminState->dspAssimMass;
                    rc = sprintf(temporary, "%23.16e", mass); rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "mass", "%s", temporary);

                    for (i=0; i<solids[j].na; i++) {
                        rc = xmlTextWriterStartElement(writer, BAD_CAST "component");
                        rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "name",    "%s", solids[j+1+i].label);
                        rc = sprintf(temporary, "%23.16e", (silminState->assimComp)[j+1+i][ns]/(silminState->assimComp)[j][ns]);
                        rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "molFrac", "%s", temporary);
                        rc = xmlTextWriterEndElement(writer);
                    }

                }
                rc = xmlTextWriterEndElement(writer); //solid
            }
        }
        rc = xmlTextWriterEndElement(writer); //assimilant
    }

    rc = xmlTextWriterEndElement(writer); //MELTSoutput
    /*
        rc = xmlTextWriterEndDocument(writer);
        xmlFreeTextWriter(writer);
    */


    if(xmlTextWriterFlush(writer) >= 0) fprintf(stderr, "XML files flushed.\n");
    else fprintf(stderr, "Error returned when attempting to flush XML files.\n");

    free (outputFile);
    free (temporary);

    if (silminState->fractionateSol || silminState->fractionateFlu || silminState->fractionateLiq)
        previousSilminState = copySilminStateStructure(silminState, previousSilminState);

    return active;

}
#endif
#endif
/* end of file READ_WRITE.C */
