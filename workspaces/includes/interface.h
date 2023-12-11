#ifndef _Interface_h
#define _Interface_h

/*
MELTS Source Code: RCS $Log: interface.h,v $
MELTS Source Code: RCS Revision 1.3  2007/02/21 21:51:18  ghiorso
MELTS Source Code: RCS New regressions options and parameter selection buttons
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2006/10/20 00:59:22  ghiorso
MELTS Source Code: RCS (1) Made initial modifications for thread safe code.
MELTS Source Code: RCS (2) Added support for XML I/O in batch mode
MELTS Source Code: RCS (3) Added support for Melts-batch listener for eventual integration into VIGMCS
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2006/08/15 16:57:36  ghiorso
MELTS Source Code: RCS xMELTS gcc 3.x sources
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2004/01/02 19:21:49  cvsaccount
MELTS Source Code: RCS CTserver University of Chicago
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:48:56  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
 * Revision 3.9  1997/06/21  22:49:51  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.8  1997/05/03  20:23:30  ghiorso
 * *** empty log message ***
 *
 * Revision 3.7  1997/03/27  17:03:33  ghiorso
 * *** empty log message ***
 *
 * Revision 3.6  1996/09/24  20:33:37  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 3.5  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 3.4  1995/09/04  19:59:58  ghiorso
 * Update to allow display of bulk composition (in grams) in the text entry
 * fields of the main silmin display. Liquid composition is no longer
 * display here, and is available only through the popup selection.
 *
 * Revision 3.4  1995/09/04  19:59:58  ghiorso
 * Update to allow display of bulk composition (in grams) in the text entry
 * fields of the main silmin display. Liquid composition is no longer
 * display here, and is available only through the popup selection.
 *
 * Revision 3.3  1995/09/01  23:54:46  ghiorso
 * Modifications made to update interface for V3.x and consolidate
 * Graph Widgets
 *
 * Revision 3.2  1995/08/31  00:31:48  ghiorso
 * Changed help widget to invoke dialogs and WEB browser.
 * Removed graphics menu entry and toggle button entries.
 *
 * Revision 3.1  1995/08/18  17:44:18  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      User Interface include file (file: INTERFACE.H)
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  August 23, 1990   Original Version
**      V1.0-2  Mark S. Ghiorso  August 30, 1991
**              (1) Removed widget creation macro definitions
**              (2) Made menu bar entry and button widgets globally
**                  known
**              (3) Removed hardcopy menu entry and button widgets
**              (4) Removed callback tags and callback declarations that
**                  needn't be global
**      V1.0-3  Mark S. Ghiorso  August 31, 1991
**              Cleaned-up a number of extern/global declarations. All
**              extern variables are now defined globally in the appropriate
**              module
**      V1.0-4  Mark S. Ghiorso  September 1, 1991
**              Added getInputDataFromFile global declaration and
**              made extern declarations of other functions ANSI
**      V1.0-5  Mark S. Ghiorso  September 3, 1991
**              (1) Removed separate storage for status in includedSolids
**                  structure
**              (2) altered magmaValues, assimilantValues and tpValues
**                  structures
**              (3) Defined certain toggle button widgets of the menu bar
**                  as global as in (V1.0-2)
**              (4) altered phases and debugEntries structures and defined
**                  more global macros for debug label widgets
**              (5) as in (4) for statusEntries
**              (6) added global macros for extra assimilant and magma mixing
**                  text widgets
**      V1.0-6  Mark S. Ghiorso  September 4, 1991
**              (1) Added units display information structure for
**                  assimilant_padb
**              (2) Added global macros to describe units display in
**                  solid_adb
**      V1.0-7  Mark S. Ghiorso  September 6, 1991
**              Added workProcData structure, and ANSI function defns for
**              XtWorkProcedures
**      V1.0-8  Mark S. Ghiorso  September 7, 1991
**              (1) Added declaraction of new liquidus work procedure
**              (2) Added cs_strings halt_me_first and empty_string
**      V1.0-9  Mark S. Ghiorso  September 9, 1991
**              (1) Made debug and plot toggles global widgets
**      V1.0-10 Mark S. Ghiorso  September 13, 1991
**              Added client_data pointer to the preclbRank structure
**      V1.0-11 Mark S. Ghiorso  September 20, 1991
**              (1) Added tg_incl_frac to global widget list
**              (2) Removed nyi_callback() declaration
**      V1.0-12 Mark S. Ghiorso  October 18, 1991
**              (1) Added globally known widget for silmin display popup
**                  attached dialog box for graph legend
**      V2.0-1  Mark S. Ghiorso  November 14, 1991
**              Conversion to OSF Motif V1.1.1/X11 Release 4
**              (1) Added declaration of wprintf function which is
**                  temporarily defined in create_minos.c
**      V2.0-2  Mark S. Ghiorso  November 23, 1991
**              (1) Removed debug entries from silmin display and menu bar
**              (2) Removed unreferenced global widgets for the silmin
**                  display and the popup dialog for graph legends
**              (3) Removed postclb_padb1 and postclb_padb2 widget declarations
**              (4) Added more message widgets
**      V2.0-3  Mark S. Ghiorso  November 27, 1991
**              (1) Removed phases[] array and changed graph_legend widget
**                  to phases widget.
**              (2) Added popup menu option to choose mode of display for
**                  log fo2 in the status box
**      V2.0-4  Mark S. Ghiorso  December 10, 1991
**              (1) Removed magma_padb and related buttons and callbacks
**      V2.0-5  Mark S. Ghiorso  December 13, 1991
**              (1) changed caddr_t pointer declaration (R3) to
**                  XtPointer (R4)
**              (2) changed tpValues structures to eliminate widgets
**              (3) altered preclbRank structure to eliminate pseudovariable
**                  in function declaration, i.e...
**                  void (*fillANOVA) (int rank);  => void (*fillANOVA) (int);
**                  (this is necessary for RISC C compliance)
**                               December 18, 1991
**              (4) Removed magma mix status entry
**              (5) item (3) revoked, RISC C is non-ANSI
**      V2.0-6  Mark S. Ghiorso  January 6, 1992
**              (1) Added declaration of two new functions:
**                  putInputDataToFile(), putOutputDataToFile()
**      V2.0-7  Mark S. Ghiorso  February 19, 1992
**              Removed references to all local variables. These became
**              multiply defined in ANSI C implementations
**      V2.0-8  Mark S. Ghiorso  March 27, 1992
**              Corrected prototype definitions for global callback procs
**      V2.0-9  Mark S. Ghiorso  June 14, 1993
**              Added new fo2 path toggle constants
**      V2.0-10 Mark S. Ghiorso March 12, 1994
**              (1) Added external declaration of create_minos().
**      V2.0-11 Mark S. Ghiorso April 23, 1994
**              Corrected declaration of typedef for preclb work_proc structure
**      V3.0-1  Mark S. Ghiorso May 11, 1994
**              (1) Added declaration of global toggle button tg_isentropic
**              (2) Added external declaration of k_mb_tg_options_isentropic
**              (3) Removed Grove and Walker entries from graphics options
**                              June 10, 1994
**              (4) Added new defines for TP_PADB tpValues instances
**              (5) Rearranged TP_PADB entries
**--
*/

/*
 *=============================================================================
 * List of globally known Widgets:
 *
 *     Name:                 Type of widget:
 */

/*
 *=============================================================================
 * List of known colors: white and black are guaranteed to exist.
 *                       on a B/W system red, green and blue are set to black
 */

/*
 *=============================================================================
 * List of globally known compound strings:
 *
 *            Name:
 */

/*
 *=============================================================================
 * Global macros
 *
 */

/*
 *=============================================================================
 * Callback tags that must be declared global:
 */

/*
 *=============================================================================
 * Callback tags associated with launch buttons on the main window:
 */

/*
 *=============================================================================
 * Global functions used as callbacks:
 *
 *   Name:                 Widget: (bt = button gadget)
 */

/*
 *=============================================================================
 * Other external functions and return value macros:
 */
int  getInputDataFromFile(char *fileName);
int  getInputDataFromLine(char *line, int *np, int *ns, double *oxideWt);
int  putInputDataToFile(char *fileName, int tpValues, int includedSolids, int assimilantValues);
int  putOutputDataToFile(char *fileName);

#define GET_INPUT_SUCCESS        0
#define GET_INPUT_ERROR_BAD_FILE 1
#define GET_INPUT_ERROR_BAD_READ 2

void doBatchFractionation(double fracOut);
int doBatchStateChange(void);
void doBatchAssimilation(void);

/*
 *=============================================================================
 *      Work proceedures:
 */

/*
 *=============================================================================
 * Globally accessible variables (silmin structures):
 */

/*
 *=============================================================================
 * Globally accessible variables (preclb structures involving the user
 *                                interface):
 */

/*
 *=============================================================================
 * Globally accessible environment and option variables:
 */

#endif /* _Interface_h */
