#include <stdio.h>
#include <stdlib.h>

#ifdef USEEDITLINE
#include <editline/readline.h>
#else
#include <readline/history.h>
#include <readline/readline.h>
#endif

#include "silmin.h"
#include "alphamelts.h"
#include "interface.h"

/*#define MODE_xMELTS           0
#define MODE__MELTS           1
#define MODE_pMELTS           2
#define MODE__MELTSandCO2     3
#define MODE__MELTSandCO2_H2O 4
#define MODE_DEFAULT          MODE_xMELTS*/

#define QUIT 0
#define MELTS_ERROR printf("Please read in MELTS file (option 1) first.\n"); \
    { char *buf = readline("Press any key to continue.\n"); if (ilog) fprintf(logfile, "\n"); free(buf); }
#define FILE_ERROR(LINE) printf("Could not read file name / path: %s\n", LINE); \
    { char *buf = readline("Press any key to continue.\n"); if (ilog) fprintf(logfile, "\n"); free(buf); }
#define MENU_ERROR(LINE) printf("Unrecognized option: %s\n", LINE); \
    { char *buf = readline("Press any key to continue.\n"); if (ilog) fprintf(logfile, "\n"); free(buf); }
#define PHASE_ERROR(LINE) printf("Unrecognized phase: %s\n", LINE); \
    { char *buf = readline("Press any key to continue.\n"); if (ilog) fprintf(logfile, "\n"); free(buf); }
#define SORRY printf("Sorry, option not yet implemented.\n"); \
    { char *buf = readline("Press any key to continue.\n"); if (ilog) fprintf(logfile, "\n"); free(buf); }

static int iAmInitialized = FALSE;
//extern int guessFlag;

void printmenu(int iprint) {

        printf("\nalphaMELTS main menu (enter 'x' when done):\n");
        if(iprint) {
            printf(" 1. Read MELTS file to set composition and system properties\n");
            printf(" 2. Twiddle starting or continuation parameters\n");
            printf(" 3. Single (batch) calculation\n");
            printf(" 4. Execute (follow path, mineral isograd or melt contour)\n");
            printf(" 5. Set fO2 buffer\n");
            printf(" 6. Set H2O (ppm) or aH2O\n");
            printf(" 7. Impose isenthalpic, isentropic or isochoric conditions\n");
            printf(" 8. Adjust solid phase settings (including fractionation)\n");
            printf(" 9. Adjust liquid or fluid settings (including melt extraction)\n");
            printf("10. Turn phase diagram mode on / off (includes wet liquidus)\n");
            printf("11. Turn geothermal, grid or PT-path file mode on / off\n");
            printf("12. Source mixer / Remixer (includes AFC and flux melting)\n");
            printf("13. Update system properties, or I/O settings\n");
            printf("14. Write out one or more MELTS files\n");
            printf("15. Write thermodynamic output for all phases\n");
            printf("16. Calculate integrated melt and output file(s)\n");
            printf("17. Fit parental melt composition with amoeba\n");
            printf("-1. Turn off menu display for options 1-17\n");
        }
        else {
            printf("-1. Turn on menu display for options 1-17\n");
        }
        printf(" X. QUIT\n");

}

int main(int argc, char** argv) {

    int mainmenu = -2, iprint = 1, ilog = 1, mode = MODE_DEFAULT; /* i.e. FALSE */
    double newestversion = 0.0;
    char* menubuf;
    FILE* logfile = fopen("logfile.txt", "a+");

    if (getenv("ALPHAMELTS_VERSION") != NULL) newestversion = atof(getenv("ALPHAMELTS_VERSION"));

    if (argc > 1) ilog = 0;

#ifdef DEBUG
    printf("Press any key to continue.\n");
    getchar();
    getchar();
#endif

    if (getenv("ALPHAMELTS_CALC_MODE") != NULL) {
        if (!strcmp(getenv("ALPHAMELTS_CALC_MODE"),"MELTS")) mode = MODE__MELTS;
        else if (!strcmp(getenv("ALPHAMELTS_CALC_MODE"),"MELTSandCO2")) mode = MODE__MELTSandCO2;
        else if (!strcmp(getenv("ALPHAMELTS_CALC_MODE"),"MELTSandCO2_H2O")) mode = MODE__MELTSandCO2_H2O;
        else if (!strcmp(getenv("ALPHAMELTS_CALC_MODE"),"pMELTS")) mode = MODE_pMELTS;
    }
    if (!mode) {
        printf("ALPHAMELTS_CALC_MODE not set!\n\n");
        while (!mode) {

            //                       11111111112222222222333333333344444444445555555555666666666677777777778888888888
            if ((menubuf = readline("---> Default calculation mode is rhyolite-MELTS (v. 1.0.2).  Change this? (y or n): ")) != NULL
                    && strlen(menubuf) > 0) add_history(menubuf);

            if (tolower(menubuf[0]) == 'y') {
                if (ilog) fprintf(logfile, "%s\n", menubuf); free(menubuf);
                //                       11111111112222222222333333333344444444445555555555666666666677777777778888888888
                if ((menubuf = readline("     Set calculation mode to rhyolite-MELTS (public release v 1.1.0)? (y or n): ")) != NULL
                        && strlen(menubuf) > 0) add_history(menubuf);
                if (tolower(menubuf[0]) == 'y') { mode = MODE__MELTSandCO2; break; }
                else {
                    if (ilog) fprintf(logfile, "%s\n", menubuf); free(menubuf);
                    if ((menubuf = readline("     Set calculation mode to rhyolite-MELTS (public release v 1.2.0)? (y or n): ")) != NULL
                            &&  strlen(menubuf) > 0) add_history(menubuf);
                    if (tolower(menubuf[0]) == 'y') { mode = MODE__MELTSandCO2_H2O; break; }
                    else {
                        if (ilog) fprintf(logfile, "%s\n", menubuf); free(menubuf);
                        if ((menubuf = readline("     Set calculation mode to pMELTS (public release v 5.6.1)? (y or n): ")) != NULL
                                && strlen(menubuf) > 0) add_history(menubuf);
                        if (tolower(menubuf[0]) == 'y') { mode = MODE_pMELTS; break; }
                        else { mode = MODE_DEFAULT; continue; }
                    }
                }
            } else { mode = MODE__MELTS; break; }
        }
        if (ilog) fprintf(logfile, "%s\n", menubuf); free(menubuf);
    }

    mode = assignGlobalStatics(mode);
    if (!mode) return 1; /* Something went wrong setting the calculation mode */
    else mode = getEnvironmentSettings();
    splashScreen(newestversion);

    printmenu(iprint);

    /* Ideally, turn off file tab completion until are actually looking for a file */
    while ((menubuf = readline("Your choice: ")) != NULL) {

        if (strlen(menubuf) > 0) {
            add_history(menubuf);
            if (tolower(menubuf[0]) == 'x') mainmenu = QUIT;
            else if (!sscanf(menubuf, "%d", &mainmenu) || mainmenu == QUIT) mainmenu = -2;
        }
        if (ilog) fprintf(logfile, "%s\n", menubuf); free(menubuf);

        if (mainmenu == QUIT) {
            /* clean up */


            /* need to think about closing xml here */

            menu_option0();
            break;
        }

        switch (mainmenu) {
        case -1:
            /* Turn off menu display for options 1-17 */
            iprint = !iprint;
            break;
        case 1:
        {
            /*Read input file to set composition of system */
            char *buf, *filename = NULL;
            int success = FALSE, iquit = FALSE;

            if ((buf = readline("MELTS filename: ")) != NULL && strlen(buf) > 0) add_history(buf);
            if (strlen(buf) < 2 && tolower(buf[0]) == 'x') iquit = TRUE;
            else success = ((filename = scanfilename(buf)) != NULL);
            if (ilog) fprintf(logfile, "%s\n", buf); free(buf);

            if (iquit) break;
            else if (!success) { FILE_ERROR(filename) }
            else iAmInitialized = menu_option1(filename);
            free(filename);


            /* need to think about closing xml here */

            if (!iAmInitialized) { MELTS_ERROR }
            // else update the completer

            break;
        }
        case 2:
            /* Twiddle starting or continuation parameters */
            if (!iAmInitialized) { MELTS_ERROR }
            else {
                char *buf;
                int menu;
                double dspTstart = silminState->dspTstart, dspPstart = silminState->dspPstart,
                    dspTstop = silminState->dspTstop, dspPstop = silminState->dspPstop,
                    dspTinc = silminState->dspTinc, dspPinc = silminState->dspPinc;

                printf("Parameters to twiddle (enter 'x' to return):\n");
                printf(" 1. Initial Temperature (%.2f C)\n", dspTstart);
                printf(" 2. Final Temperature (%.2f C)\n", dspTstop);
                printf(" 3. Temperature Increment (%.2f C)\n", dspTinc);
                printf(" 4. Initial Pressure (%.2f bars)\n", dspPstart);
                printf(" 5. Final Pressure (%.2f bars)\n", dspPstop);
                printf(" 6. Pressure Increment (%.2f bars)\n", dspPinc);

                while ((buf = readline("Choose: ")) != NULL) {

                    if (strlen(buf) > 0) {
                        add_history(buf);
                        if (tolower(buf[0]) == 'x') break;
                    }

                    if (!sscanf(buf, "%d", &menu) || menu < 1 || menu > 6) {
                        if (ilog) fprintf(logfile, "%s\n", buf); { MENU_ERROR(buf) }
                        free(buf);
                    }
                    else {
                        if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                        switch(menu) {
                        case 1:
                            if ((buf = readline("Type new Tstart: ")) != NULL && strlen(buf) > 0) add_history(buf);
                            if (tolower(buf[0]) == 'x' || !sscanf(buf, "%lf", &dspTstart)) dspTstart = silminState->dspTstart;
                            else dspTstart = (dspTstart > 0.0) ? dspTstart : silminState->dspTstart;
                            break;
                        case 2:
                            if ((buf = readline("Type new Tstop (or 0 to use current Tstart):")) != NULL
                                && strlen(buf) > 0) add_history(buf);
                            if (tolower(buf[0]) == 'x' || !sscanf(buf, "%lf", &dspTstop)) dspTstop = silminState->dspTstop;
                            else dspTstop = (dspTstop > 0.0) ? dspTstop : dspTstart;
                            break;
                        case 3:
                            if ((buf = readline("Type new Tinc: ")) != NULL && strlen(buf) > 0) add_history(buf);
                            if (tolower(buf[0]) == 'x' || !sscanf(buf, "%lf", &dspTinc)) dspTinc = silminState->dspTinc;
                            break;
                        case 4:
                            if ((buf = readline("Type new Pstart: ")) != NULL && strlen(buf) > 0) add_history(buf);
                            if (tolower(buf[0]) == 'x' || !sscanf(buf, "%lf", &dspPstart)) dspPstart = silminState->dspPstart;
                            else dspPstart = (dspPstart > 0.0) ? dspPstart : silminState->dspPstart;
                            break;
                        case 5:
                            if ((buf = readline("Type new Pstop (or 0 to use current Pstop): ")) != NULL
                                && strlen(buf) > 0) add_history(buf);
                            if (tolower(buf[0]) == 'x' || !sscanf(buf, "%lf", &dspPstop)) dspPstop = silminState->dspPstop;
                            else dspPstop = (dspPstop > 0.0) ? dspPstop : dspPstart;
                            break;
                        case 6:
                            if ((buf = readline("Type new Pinc: ")) != NULL && strlen(buf) > 0) add_history(buf);
                            if (tolower(buf[0]) == 'x' || !sscanf(buf, "%lf", &dspPinc)) dspPinc = silminState->dspPinc;
                            break;
                        default:
                            break;
                        }
                        if (ilog) fprintf(logfile, "%s\n", buf); free(buf);

                    }

                    printf("Parameters to twiddle (enter 'x' to return):\n");
                    printf(" 1. Initial Temperature (%.2f)\n", dspTstart);
                    printf(" 2. Final Temperature (%.2f)\n", dspTstop);
                    printf(" 3. Temperature Increment (%.2f)\n", dspTinc);
                    printf(" 4. Initial Pressure (%.2f bars)\n", dspPstart);
                    printf(" 5. Final Pressure (%.2f bars)\n", dspPstop);
                    printf(" 6. Pressure Increment (%.2f bars)\n", dspPinc);

                }

                if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                menu_option2(dspTstart, dspTstop, dspTinc, dspPstart, dspPstop, dspPinc);

            }
            break;
        case 3:
            /* Single (batch) calculation  without constraints (except fO2) */
            if (!iAmInitialized) { MELTS_ERROR }
            else {
                char *buf;
                int equilibriumGuess = 1, iquit = 0;

                if(!guessFlag) {
                    printf("Starting guess to use (enter 'x' to return):\n0. Norm-calculated subsolidus assemblage\n");
                    printf("1. Superliquidus starting guess at current temperature (%.2f C).\n", silminState->dspTstart);
                    printf("2. Superliquidus starting guess at liquidus temperature.\n");
                    printf("3. Superliquidus starting guess at wet liquidus temperature.\n");
                    if ((buf = readline("Choose: ")) != NULL && strlen(buf) > 0) add_history(buf);
                    if (tolower(buf[0]) == 'x') iquit = 1;
                    else if (!sscanf(buf, "%d", &equilibriumGuess)) equilibriumGuess = silminState->incLiquids;
                    if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                    if (iquit) break;
                }
                menu_option3(equilibriumGuess);
            }
            break;
        case 4:
            /* Execute (follow path, mineral isograd or melt contour) */
            if (!iAmInitialized) { MELTS_ERROR }
            else {
                char *buf;
                int equilibriumGuess = 1, iquit = 0, saveAll = 1, iterMax = 0;

                if (!guessFlag) {
                    printf("Starting guess to use (enter 'x' to return):\n0. Norm-calculated subsolidus assemblage\n");
                    printf("1. Superliquidus starting guess at current temperature (%.2f C).\n", silminState->dspTstart);
                    printf("2. Superliquidus starting guess at liquidus temperature.\n");
                    printf("3. Superliquidus starting guess at wet liquidus temperature.\n");
                    if ((buf = readline("Choose: ")) != NULL && strlen(buf) > 0) add_history(buf);
                    if (tolower(buf[0]) == 'x') iquit = 1;
                    else if (!sscanf(buf, "%d", &equilibriumGuess)) equilibriumGuess = silminState->incLiquids;
                    if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                    if (iquit) break;
                }

                if (guessFlag && (silminState->txtOutput > TEXT_TABLE)) {
                    if ((buf = readline("Record output / fractionations from scratch (1) or append to previous (0)? ")) != NULL
                        && strlen(buf) > 0) add_history(buf);
                    if (tolower(buf[0]) == 'x') iquit = 1;
                    else if (!sscanf(buf, "%d", &saveAll)) saveAll = 1;
                    else saveAll = !saveAll;
                    if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                    if (iquit) break;
                }

                printf("Maximum number of steps to take, or 0 for 'unlimited' (i.e. limited by other criteria);\n");
                if ((buf = readline("enter negative number to limit failed steps instead (e.g. -4 to skip failure 4 times): "))
                    != NULL && strlen(buf) > 0) add_history(buf);
                if (tolower(buf[0]) == 'x') iquit = 1;
                else if (!sscanf(buf, "%d", &iterMax)) iterMax = 0;
                if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                if (iquit) break;

                menu_option4(equilibriumGuess, saveAll, iterMax);

            }
            break;
        case 5:
            /* Set fO2 buffer */
            /* NEED TO FIX UP FO2LIQ IMPLEMENTATION? */

            if (!iAmInitialized) { MELTS_ERROR }
            else {
                char *buf;
                int fo2Path = silminState->fo2Path, menu;
                double fo2Delta = silminState->fo2Delta;
                double fo2Offset = (fo2Path <= 0 && fo2Delta == 0.0) ? silminState->fo2 : fo2Delta;

                printf("Oxygen buffer to impose (enter 'x' to return):\n");
                printf(" 0. Select no buffer%s\n", fo2Path == FO2_NONE ? " (current)" : "");
                printf(" 1. Select Hematite-Magnetite%s\n", fo2Path == FO2_HM ? " (current)" : "");
                printf(" 2. Select Nickel-Nickel Oxide%s\n", fo2Path == FO2_NNO ? " (current)" : "");
                printf(" 3. Select Fayalite-Magnetite-Quartz%s\n", fo2Path == FO2_QFM ? " (current)" : "");
                printf(" 4. Select Graphite-COH Fluid%s\n", fo2Path == FO2_COH ? " (current)" : "");
                printf(" 5. Select Iron-Wustite%s\n", fo2Path == FO2_IW ? " (current)" : "");
                if (fo2Path <= 0) printf(" 6. Use absolute logfo2 offset (%s %f)\n", fo2Path == FO2_ABS ? "current" : "default", fo2Offset);
                else printf(" 6. Use absolute logfo2 offset (default %f from option %d)\n", fo2Offset, fo2Path);
                if (fo2Path <= 0) printf(" 7. Change offset from buffer, if any (current %f).\n", fo2Delta);
                else printf(" 7. Change offset from buffer, if any (current %f from option %d).\n", fo2Delta, fo2Path);

                while ((buf = readline("Choose: ")) != NULL) {

                    if (strlen(buf) > 0) {
                        add_history(buf);
                        if (tolower(buf[0]) == 'x') break;
                    }

                    if (!sscanf(buf, "%d", &menu) || menu < 0 || menu > 7) {
                        if (ilog) fprintf(logfile, "%s\n", buf);  { MENU_ERROR(buf) }
                        free(buf);
                        fo2Path = (silminState->fo2Path < 0 ? 6 : silminState->fo2Path);
                    }
                    else if (menu == 7) {
                        if ((buf = readline("Offset from buffer in log 10 units: ")) != NULL && strlen(buf) > 0) add_history(buf);
                        if (tolower(buf[0]) == 'x' || !sscanf(buf, "%lf", &fo2Delta)) fo2Delta = silminState->fo2Delta;
                        if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                    }
                    else if (menu == 6) {
                        fo2Offset += ((fo2Path != FO2_NONE) ? (getlog10fo2(silminState->T, silminState->P, fo2Path) - silminState->fo2Delta) : 0.0);
                        fo2Path = -1;
                        fo2Delta = fo2Offset;
                        if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                    }
                    else {
                        fo2Offset += ((fo2Path != FO2_NONE) ? (getlog10fo2(silminState->T, silminState->P, fo2Path) - silminState->fo2Delta) : 0.0);
                        fo2Path = menu;
                        fo2Offset -= ((fo2Path != FO2_NONE) ? (getlog10fo2(silminState->T, silminState->P, fo2Path) - silminState->fo2Delta) : 0.0);
                        if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                    }

                    printf("Oxygen buffer to impose (enter 'x' when done):\n");
                    printf(" 0. Select no buffer%s\n", fo2Path == FO2_NONE ? " (current)" : "");
                    printf(" 1. Select Hematite-Magnetite%s\n", fo2Path == FO2_HM ? " (current)" : "");
                    printf(" 2. Select Nickel-Nickel Oxide%s\n", fo2Path == FO2_NNO ? " (current)" : "");
                    printf(" 3. Select Fayalite-Magnetite-Quartz%s\n", fo2Path == FO2_QFM ? " (current)" : "");
                    printf(" 4. Select Graphite-COH Fluid%s\n", fo2Path == FO2_COH ? " (current)" : "");
                    printf(" 5. Select Iron-Wustite%s\n", fo2Path == FO2_IW ? " (current)" : "");
                    if (fo2Path <= 0) printf(" 6. Use absolute logfo2 offset (%s %f)\n", fo2Path == FO2_ABS ? "current" : "default", fo2Offset);
                    else printf(" 6. Use absolute logfo2 offset (default %f from option %d)\n", fo2Offset, fo2Path);
                    if (fo2Path <= 0) printf(" 7. Change offset from buffer, if any (current %f).\n", fo2Delta);
                    else printf(" 7. Change offset from buffer, if any (current %f from option %d).\n", fo2Delta, fo2Path);
                }
                if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                menu_option5(fo2Path, fo2Delta);

            }
            break;

        case 6:
            /* Set H2O (ppm) or aH2O */
            { SORRY }
            break;
        case 7:
            /* Impose entropy, enthalpy, or volume */

            if (!iAmInitialized) { MELTS_ERROR }
            else {
                char *buf;
                int isentropic = silminState->isentropic, isenthalpic = silminState->isenthalpic, isochoric = silminState->isochoric,
                    constraint = (isenthalpic) ? 1 : ((isentropic) ? 2 : ((isochoric) ? 3 : 0)), menu;
                double newRef = (isenthalpic) ? silminState->refEnthalpy :
                    ((isentropic) ? silminState->refEntropy : ((isochoric) ? silminState->refVolume : 0.0)),
                    newStop = (isenthalpic) ? silminState->dspHstop : ((isentropic) ? silminState->dspSstop : ((isochoric) ? silminState->dspVstop : 0.0)),
                    newInc = (isenthalpic) ? silminState->dspHinc : ((isentropic) ? silminState->dspSinc : ((isochoric) ? silminState->dspVinc : 0.0));

                printf("Constraint to impose (enter 'x' when done):\n");
                printf(" 0. None%s\n", (!isenthalpic && !isentropic && !isochoric) ? " (current)" : "");
                printf(" 1. Isenthalpic%s\n", isenthalpic ? " (current)" : "");
                printf(" 2. Isentropic%s\n", isentropic ? " (current)" : "");
                printf(" 3. Isochoric%s\n", isochoric ? " (current)" : "");
                printf(" 4. Change reference value (%f) for current constraint (option %d).\n", newRef, constraint);
                printf(" 5. Change stop value (%f) for current constraint (option %d).\n", newStop, constraint);
                printf(" 6. Change increment (%f) for current constraint (option %d).\n", newInc, constraint);

                while ((buf = readline("Choose: ")) != NULL) {

                    if (strlen(buf) > 0) {
                        add_history(buf);
                        if (tolower(buf[0]) == 'x') break;
                    }

                    if (!sscanf(buf, "%d", &menu) || menu < 0 || menu > 6) {
                        if (ilog) fprintf(logfile, "%s\n", buf); { MENU_ERROR(buf) }
                        free(buf);
                    }
                    else {
                        if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                        if (menu < 4) constraint = menu;
                        switch(menu) {
                        case 0:
                            isenthalpic = FALSE; isentropic = FALSE; isochoric = FALSE;
                            break;
                        case 1:
                            isenthalpic = TRUE; isentropic = FALSE; isochoric = FALSE;
                            if (guessFlag && newRef == 0.0) newRef = silminState->bulkTD.h;
                            break;
                        case 2:
                            isenthalpic = FALSE; isentropic = TRUE; isochoric = FALSE;
                            if (guessFlag && newRef == 0.0) newRef = silminState->bulkTD.s;
                            break;
                        case 3:
                            isenthalpic = FALSE; isentropic = FALSE; isochoric = TRUE;
                            if (guessFlag && newRef == 0.0) newRef = silminState->bulkTD.v;
                            break;
                        case 4:
                   if (isenthalpic && (buf = readline("Type new Hstart (or 0 to calculate later): ")) != NULL && strlen(buf) > 0) add_history(buf);
                            else if (isentropic  && (buf = readline("Type new Sstart (or 0 to calculate later): ")) != NULL && strlen(buf) > 0) add_history(buf);
                            else if (isochoric   && (buf = readline("Type new Vstart (or 0 to calculate later): ")) != NULL && strlen(buf) > 0) add_history(buf);
                            else if       ((buf = readline("Type new reference value: ")) != NULL && strlen(buf) > 0) add_history(buf);
                            if (tolower(buf[0]) == 'x' || !sscanf(buf, "%lf", &newRef)) newRef = (isenthalpic) ? silminState->refEnthalpy :
                                ((isentropic) ? silminState->refEntropy : ((isochoric) ? silminState->refVolume : 0.0));
                            if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                            break;
                        case 5:
                   if (isenthalpic && (buf = readline("Type new Hstop (or 0 to use current Hstart):")) != NULL && strlen(buf) > 0) add_history(buf);
                            else if (isentropic  && (buf = readline("Type new Sstop (or 0 to use current Sstart):")) != NULL && strlen(buf) > 0) add_history(buf);
                            else if (isochoric   && (buf = readline("Type new Vstop (or 0 to use current Vstart):")) != NULL && strlen(buf) > 0) add_history(buf);
                            else if      ((buf = readline("Type new stop value (or 0 to use current start value):")) != NULL && strlen(buf) > 0) add_history(buf);
                            if (tolower(buf[0]) == 'x' || !sscanf(buf, "%lf", &newStop)) newStop = (isenthalpic) ? silminState->dspHstop :
                                ((isentropic) ? silminState->dspSstop : ((isochoric) ? silminState->dspVstop : 0.0));
                            else newStop = (newStop != 0.0) ? newStop : newRef;
                            if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                            break;
                        case 6:
                   if (isenthalpic && (buf = readline("Type new Hinc: ")) != NULL && strlen(buf) > 0) add_history(buf);
                            else if (isentropic  && (buf = readline("Type new Sinc: ")) != NULL && strlen(buf) > 0) add_history(buf);
                            else if (isochoric   && (buf = readline("Type new Vinc: ")) != NULL && strlen(buf) > 0) add_history(buf);
                            else if           ((buf = readline("Type new increment: ")) != NULL && strlen(buf) > 0) add_history(buf);
                            if (tolower(buf[0]) == 'x' || !sscanf(buf, "%lf", &newInc)) newInc = (isenthalpic) ? silminState->dspHinc :
                                ((isentropic) ? silminState->dspSinc : ((isochoric) ? silminState->dspVinc : 0.0));
                            if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                            break;
                        }
                    }

                    printf("Constraint to impose (enter 'x' when done):\n");
                    printf(" 0. None%s\n", (!isenthalpic && !isentropic && !isochoric) ? " (current)" : "");
                    printf(" 1. Isenthalpic%s\n", isenthalpic ? " (current)" : "");
                    printf(" 2. Isentropic%s\n", isentropic ? " (current)" : "");
                    printf(" 3. Isochoric%s\n", isochoric ? " (current)" : "");
                    printf(" 4. Change reference value (%f) for current constraint (option %d).\n", newRef, constraint);
                    printf(" 5. Change stop value (%f) for current constraint (option %d).\n", newStop, constraint);
                    printf(" 6. Change increment (%f) for current constraint (option %d).\n", newInc, constraint);

                }
                if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                menu_option7(constraint, newRef, newStop, newInc);
            }

            break;
        case 8:

            if (!iAmInitialized) { MELTS_ERROR }
            else {
                char *buf;
                int menu;
                while ((buf = readline("Solid phase(s) to change settings for ('all' for all; 'x' when done): ")) != NULL) {

                    int i, j, len = strlen(buf);
                    if (len > 0) {
                        add_history(buf);
                        if (tolower(buf[0]) == 'x') break;
                        else for (i = 0; i<len; i++) buf[i] = tolower(buf[i]);
                    }

                    if (!strncmp(buf, "all", MIN(len, 3))) {
                        i = -1;
                    }
                    // All called fluid now but keep so old batch files work
                    else if (!strncmp(buf, "water", MIN(len, 5)) || !strncmp(buf, "fluid", MIN(len, 5))) {
                        i = npc;
                    }
                    else {
                        for (i=0, j=0; i<npc; i++) {
                            if (solids[i].type == PHASE) {
                                int phaseStrLen = (int) strlen(solids[i].label);
                                if (((len-phaseStrLen)  == 0) && !strncmp(buf, solids[i].label, phaseStrLen)) {
                                    break;
                                }
                                j++;
                            }
                        }
                    }
                    if (i == npc || (i >= 0 && solids[i].nr > 0 && solids[i].convert == NULL)) {
                        if (ilog) fprintf(logfile, "%s\n", buf); { PHASE_ERROR(buf) }
                        free(buf);
                    }
                    else {
                        int incSolids = -2, fracSolids = -2;
                        int massF, ratio, iratio, minType;
                        double minF, newMin;

                        if (ilog) fprintf(logfile, "%s\n", buf); free(buf);

                        if (i < 0) { /* all solids in standard set */

                            massF = silminState->maxF < 1.0;
                            ratio = (!massF && silminState->fracOut < 1.0);
                            iratio = ratio && (silminState->fracSolids != NULL);

                            minType = (massF) ? 1 : ((ratio && !iratio) ? 2 : ((iratio) ? 3 : 0));

                            minF = silminState->maxF;
                            newMin = (massF) ? minF : ((ratio) ? silminState->fracOut : 1.0);

                            if (incSolids < -1) {
                                for (i=0, j=0; i<npc; i++) if (solids[i].type == PHASE) {
                                        int phaseStrLen = strlen(solids[i].label);
                                        if (strncmp(solids[i].label, "fluid", MIN(phaseStrLen, 5)) && (solids[i].inStdSet) && silminState->incSolids[j]) break;
                                        j++;
                                }
                                incSolids = (i != npc) ? -1 : 0;
                                i = -1;
                            }
                            if (fracSolids < -1) fracSolids = (silminState->fractionateSol) ? -1 : 0; // check this!!

                        }
                        else { /* a particular solid */

                            newMin =  (silminState->fracSolids != NULL) ? silminState->fracSolids[i] : silminState->fracOut;

                            massF = silminState->maxF < 1.0;
                            ratio = (!massF && newMin < 1.0);
                            iratio = ratio && (silminState->fracSolids != NULL);

                            minType = (massF) ? 1 : ((ratio && !iratio) ? 2 : ((iratio) ? 3 : 0));

                            minF = silminState->maxF;
                            newMin = (massF) ? minF : ((ratio) ? newMin : silminState->fracOut);

                            if (incSolids < -1) incSolids = silminState->incSolids[j];
                            if (fracSolids < -1)
                                fracSolids = (silminState->fracSolids != NULL) ? (silminState->fracSolids[j-1] != 0.0) : silminState->fractionateSol;;

                        }

                        printf("Setting to adjust (enter 'x' when done):\n");
                        if (i < 0) {
                            printf(" 1. Suppress all solid phases%s\n", incSolids ? "" : " (current)");
                            printf(" 2. Allow all standard solid phases; optionally limit maximum number\n");
                            printf(" 3. Do not fractionate any solid phases%s\n", fracSolids ? "" : " (current)");
                            printf(" 4. Fractionate all standard solid phases%s\n", fracSolids ? " (current)" : "");
                            printf(" 5. Change whether residual solid is measured by solid mass frac (1-F%s), or fractionation coefficient(s)%s, or default%s\n",
                                (massF) ? "; current" : "", (ratio) ? " (current)" : "", (!massF && !ratio) ? " (current)" : "");
                            printf(" 6. Adjust amount of solid %s during fractionation (%f) unless individual values or default bahavior\n", (massF) ? "retained" : "fractionated", newMin);
                        }
                        else {
                            printf(" 1. Suppress phase %s%s\n", solids[i].label, incSolids ? "" : " (current)");
                            if ((incSolids >= 0) && (incSolids < solids[i].na))
                                printf(" 2. Limit number of %s allowed to join assemblage (%d)\n", solids[i].label, incSolids);
                            else
                                printf(" 2. Limit number of %s allowed to join assemblage (unlimited)\n", solids[i].label);
                            printf(" 3. Do not fractionate phase %s%s\n", solids[i].label, fracSolids ? "" : " (current)");
                            printf(" 4. Fractionate phase %s%s\n", solids[i].label, fracSolids ? " (current)" : "");
                            printf(" 5. Change whether residual %s is measured by solid mass frac (1-F%s), or fractionation coefficient(s)%s, or default%s\n",
                                solids[i].label, (massF) ? "; current" : "", (ratio) ? " (current)" : "", (!massF && !ratio) ? " (current)" : "");
                            printf(" 6. Adjust amount of %s %s during fractionation (%f) unless default behavior\n", solids[i].label, (massF) ? "retained" : "fractionated", newMin);
                        }

                        while((buf = readline("Choose: ")) != NULL) {

                            if (strlen(buf) > 0) {
                                add_history(buf);
                                if (tolower(buf[0]) == 'x') break;
                            }

                            if (!sscanf(buf, "%d", &menu) || menu < 1 || menu > 6) {
                                if (ilog) fprintf(logfile, "%s\n", buf); { MENU_ERROR(buf) }
                                free(buf);
                            }
                            else {
                                if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                                switch(menu) {
                                case 1:
                                    incSolids = FALSE;
                                    break;
                                case 2:
                                    if ((buf = readline("Maximum number of single phase allowed (or 0 for unlimited): ")) != NULL
                                        && strlen(buf) > 0) add_history(buf);
                                    if (tolower(buf[0]) == 'x' || !sscanf(buf, "%d", &incSolids)) incSolids = -2; /* unchanged */
                                    if (!incSolids) incSolids = -1; /* unlimited */
                                    if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                                    break;
                                case 3:
                                    fracSolids = FALSE;
                                    break;
                                case 4:
                                    fracSolids = TRUE;
                                    break;
                                case 5:
                                    if (i < 0) {
                                        printf(" 1. Residual solid measured by solid mass frac (1-F%s)\n", (massF) ? "; current" : "");
                                        printf(" 2. Residual solid measured from globally set ratio fractionated%s, relative to solid\n", (ratio && !iratio) ? " (current)" : "");
                                        printf(" 3. Residual solid measured from individual ratios fractionated%s, relative to solid\n", (ratio) ? " (current)" : "");
                                        printf(" 4. Default behavior (1.e-5 mols per phase retained%s)\n", (!massF && !ratio) ? "; current" : "");
                                    }
                                    else {
                                        printf(" 1. Residual %s measured by solid mass frac (1-F%s)\n", solids[i].label, (massF) ? "; current" : "");
                                        printf(" 2. Residual %s measured from globally set ratio extracted%s, relative to %s\n", solids[i].label,
                                            (ratio) ? " (current)" : "", solids[i].label);
                                        printf(" 3. Residual %s measured from individual ratios extracted%s, relative to %s\n", solids[i].label,
                                            (ratio && !iratio) ? " (current)" : "", solids[i].label);
                                        printf(" 4. Default behavior (1.e-5 mols per phase retained%s)\n", (!massF && !ratio) ? "; current" : "");
                                    }
                                    if ((buf = readline("Choose: ")) != NULL && strlen(buf) > 0) add_history(buf);
                                    if (tolower(buf[0]) == 'x' || !sscanf(buf, "%d", &minType)) minType = (massF) ? 1 : ((ratio && !iratio) ? 2 : ((iratio) ? 3 : 0));
                                    switch(minType) {
                                        case 1:
                                            massF = TRUE; ratio = FALSE;
                                            break;
                                        case 2:
                                            massF = FALSE; ratio = TRUE; iratio = FALSE;
                                            break;
                                        case 3:
                                            massF = FALSE; ratio = TRUE; iratio = TRUE;
                                            break;
                                        case 4:
                                            massF = FALSE; ratio = FALSE;
                                            break;
                                        default:
                                            { MENU_ERROR(buf) }
                                            break;
                                    }
                                    if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                                    break;
                                case 6:
                                    {
                                        double oldMin = newMin;
                                        if      ((massF) && ((buf = readline("Mass frac (1-F) of solids to retain during extraction (relative to system): ")) != NULL
                                            && strlen(buf) > 0)) add_history(buf);
                                        else if ((ratio) && ((buf = readline("Ratio of phase to extract (relative to total amount of phase): ")) != NULL
                                            && strlen(buf) > 0)) add_history(buf);
                                        else break;
                                        if (tolower(buf[0]) == 'x' || !sscanf(buf, "%lf", &newMin)) newMin = oldMin;
                                    }
                                    if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                                    break;
                                }
                            }

                            printf("Setting to adjust (enter 'x' when done):\n");
                            if (i < 0) {
                                printf(" 1. Suppress all solid phases%s\n", incSolids ? "" : " (current)");
                                printf(" 2. Allow all standard solid phases; optionally limit maximum number\n");
                                printf(" 3. Do not fractionate any solid phases%s\n", fracSolids ? "" : " (current)");
                                printf(" 4. Fractionate all standard solid phases%s\n", fracSolids ? " (current)" : "");
                                printf(" 5. Change whether residual solid is measured by solid mass frac (1-F%s), or fractionation coefficient(s)%s, or default%s\n",
                                    (massF) ? "; current" : "", (ratio) ? " (current)" : "", (!massF && !ratio) ? " (current)" : "");
                                printf(" 6. Adjust amount of solid %s during fractionation (%f) unless individual values or default behavior\n", (massF) ? "retained" : "fractionated", newMin);
                            }
                            else {
                                printf(" 1. Suppress phase %s%s\n", solids[i].label, incSolids ? "" : " (current)");
                                if ((incSolids >= 0) && (incSolids < solids[i].na))
                                    printf(" 2. Limit number of %s allowed to join assemblage (%d)\n", solids[i].label, incSolids);
                                else
                                    printf(" 2. Limit number of %s allowed to join assemblage (unlimited)\n", solids[i].label);
                                printf(" 3. Do not fractionate phase %s%s\n", solids[i].label, fracSolids ? "" : " (current)");
                                printf(" 4. Fractionate phase %s%s\n", solids[i].label, fracSolids ? " (current)" : "");
                                printf(" 5. Change whether residual %s is measured by solid mass frac (1-F%s), or fractionation coefficient(s)%s, or default%s\n",
                                    solids[i].label, (massF) ? "; current" : "", (ratio) ? " (current)" : "", (!massF && !ratio) ? " (current)" : "");
                                printf(" 6. Adjust amount of %s %s during fractionation (%f) unless default behavior\n", solids[i].label, (massF) ? "retained" : "fractionated", newMin);
                            }
                        }

                        if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                        menu_option8(i, j, incSolids, fracSolids, minType, newMin);

                    }
                }
                if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
            }
            break;
        case 9:
            /* Turn liquid on / off */
            if (!iAmInitialized) { MELTS_ERROR }
            else {
                char *buf;

              	int haveWater = ((calculationMode == MODE__MELTS) || (calculationMode == MODE_pMELTS));
                while ((buf = readline((haveWater) ? "Phase(s) to change settings for ('liquid', 'liquid n', where n>=1, or 'fluid' or 'water'; 'x' when done): " :
                    "Phase(s) to change settings for ('liquid', 'liquid n', where n>=1, or 'fluid'; 'x' when done): ")) != NULL) {

                    int i, j, len = strlen(buf);
                    if (len > 0) {
                        add_history(buf);
                        if (tolower(buf[0]) == 'x') break;
                        else for (i = 0; i<len; i++) buf[i] = tolower(buf[i]);
                    }

                    if (!strncmp(buf, "liquid", MIN(len, 6))) {
                        i = -1; j = 0;
                        if (len > 6) {
                            if (sscanf(&buf[6], "%d", &j) == EOF) j = 0; /* all liquids */
                            else if (j < 1 || j > nlc) j = 0;
                        }
                    }
                    // All called fluid now but keep so old batch files work
                    else if (haveWater && !strncmp(buf, "water", MIN(len, 5))) {
                        for (i=0, j=0; i<npc; i++) {
                            if (solids[i].type == PHASE) {
                                int phaseStrLen = (int) strlen(solids[i].label);
                                if (((len-phaseStrLen)  == 0) && !strncmp("fluid", solids[i].label, phaseStrLen)) {
                                    break;
                                }
                                j++;
                            }
                        }
                    }
                    else if (!strncmp(buf, "fluid", MIN(len, 5))) {
                        for (i=0, j=0; i<npc; i++) {
                            if (solids[i].type == PHASE) {
                                int phaseStrLen = (int) strlen(solids[i].label);
                                if (((len-phaseStrLen)  == 0) && !strncmp("fluid", solids[i].label, phaseStrLen)) {
                                    break;
                                }
                                j++;
                            }
                        }
                    }
                    if (i == npc || (i >= 0 && solids[i].nr > 0 && solids[i].convert == NULL)) {
                        if (ilog) fprintf(logfile, "%s\n", buf); { PHASE_ERROR(buf) }
                        free(buf);
                    }
                    else {
                        int incLiquids = -2, fracLiquids = -2;
                        int massF, volF, ratio, iratio, menu, minType;
                        double minF, minPhi, newMin;

                        if (ilog) fprintf(logfile, "%s\n", buf); free(buf);

                        if (i < 0 && j) { /* a particular liquid (j counted from 1)*/

                            newMin = (silminState->fracLiquids != NULL) ? silminState->fracLiquids[j-1] : silminState->fracOut;

                            volF = (silminState->minLiqPhi != 0.0);
                            massF = (!volF && silminState->minF != 0.0);
                            ratio = (!volF && !massF && newMin < 1.0);
                            iratio = ratio && (silminState->fracLiquids != NULL);

                            minType = (volF) ? 1 : ((massF) ? 2 : ((ratio && !iratio) ? 3 : ((iratio) ? 4 : 0)));

                            minPhi = silminState->minLiqPhi;
                            minF = silminState->minF;
                            newMin = (volF) ? minPhi : ((massF) ? minF : ((ratio) ? newMin : silminState->fracOut));

                            if (incLiquids < -1) incLiquids = (silminState->incLiquids >= j);
                            if (fracLiquids < -1)
                                fracLiquids = (silminState->fracLiquids != NULL) ? (silminState->fracLiquids[j-1] != 0.0) : silminState->fractionateLiq;;

                        }
                        else if (i < 0) { /* all liquids */

                            volF = (silminState->minLiqPhi != 0.0);
                            massF = (!volF && silminState->minF != 0.0);
                            ratio = (!volF && !massF && silminState->fracOut < 1.0);
                            iratio = ratio && (silminState->fracLiquids != NULL);

                            minType = (volF) ? 1 : ((massF) ? 2 : ((ratio && !iratio) ? 3 : ((iratio) ? 4 : 0)));

                            minPhi = silminState->minLiqPhi;
                            minF = silminState->minF;
                            newMin = (volF) ? minPhi : ((massF) ? minF : ((ratio) ? silminState->fracOut : 1.0));

                            if (incLiquids < -1) incLiquids = silminState->incLiquids;
                            if (fracLiquids < -1) fracLiquids = silminState->fractionateLiq;

                        }
                        else { /* fluid phase (j counted from 0) */

                            newMin =  (silminState->fracFluids != NULL) ? silminState->fracFluids[0] : silminState->fracOut;

                            volF = (silminState->minFluPhi != 0.0);
                            massF = (!volF && silminState->maxF != 1.0);
                            ratio = (!volF && !massF && newMin < 1.0);
                            iratio = ratio && (silminState->fracFluids != NULL);

                            minType = (volF) ? 1 : ((massF) ? 2 : ((ratio && !iratio) ? 3 : ((iratio) ? 4 : 0)));

                            minPhi = silminState->minFluPhi;
                            minF = 1.0 - silminState->maxF;
                            newMin = (volF) ? minPhi : ((massF) ? minF : ((ratio) ? newMin : silminState->fracOut));

                            if (incLiquids < -1) incLiquids = silminState->incSolids[j];
                            if (fracLiquids < -1) fracLiquids = silminState->fractionateFlu;

                        }

                        printf("Setting to adjust (enter 'x' when done):\n");
                        if (i < 0 && j) {
                            printf(" 1. Suppress liquid phase %d %s\n", j, incLiquids ? "" : " (current)");
                            if ((incLiquids >= 0) && (incLiquids < nlc))
                                printf(" 2. Allow at least %d liquids to join assemblage (current limit %d)\n", j, incLiquids);
                            else
                                printf(" 2. Allow at least %d liquids to join assemblage (unlimited)\n", j);
                            printf(" 3. Do not extract liquid %d %s\n", j, fracLiquids ? "" : " (current)");
                            printf(" 4. Extract liquid %d %s\n", j, fracLiquids ? " (current)" : "");
                            printf(" 5. Change whether residual liquids are measured by vol frac (phi%s), mass frac (F%s), from ratio extracted%s, or default%s\n",
                                (volF) ? "; current" : "", (massF) ? "; current" : "", (ratio) ? " (current)" : "", (!volF && !massF && !ratio) ? " (current)" : "");
                            printf(" 6. Adjust amount of liquid %d %s during extraction (%f) unless default behavior\n", j, (volF || massF) ? "retained" : "extracted", newMin);
                        }
                        else if (i < 0) {
                            printf(" 1. Suppress liquid phase %s\n", incLiquids ? "" : " (current)");
                            if ((incLiquids >= 0) && (incLiquids < nlc))
                                printf(" 2. Limit number of multiple liquids allowed to join assemblage (%d)\n", incLiquids);
                            else
                                printf(" 2. Limit number of multiple liquids allowed to join assemblage (unlimited)\n");
                            printf(" 3. Do not extract any liquids %s\n", fracLiquids ? "" : " (current)");
                            printf(" 4. Extract all liquids %s\n", fracLiquids ? " (current)" : "");
                            printf(" 5. Change whether residual liquids are measured by vol frac (phi%s), mass frac (F%s), from ratio(s) extracted%s, or default%s\n",
                                (volF) ? "; current" : "", (massF) ? "; current" : "", (ratio) ? " (current)" : "", (!volF && !massF && !ratio) ? " (current)" : "");
                            printf(" 6. Adjust amount of liquid %s during extraction (%f) unless individual values or default behavior\n", (volF || massF) ? "retained" : "extracted", newMin);
                        }
                        else {
                            printf(" 1. Suppress %s phase %s\n", solids[i].label, incLiquids ? "" : " (current)");
                            if ((incLiquids >= 0) && (incLiquids < solids[i].na))
                                printf(" 2. Limit number of %s phase allowed to join assemblage (%d)\n", solids[i].label, incLiquids);
                            else
                                printf(" 2. Limit number of %s phase allowed to join assemblage (unlimited)\n", solids[i].label);
                            printf(" 3. Do not fractionate / degas %s phase %s\n", solids[i].label, fracLiquids ? "" : " (current)");
                            printf(" 4. Fractionate / degas %s phase %s\n", solids[i].label, fracLiquids ? " (current)" : "");
                            printf(" 5. Change whether residual %s is measured by %s vol frac%s, solid mass frac (1-F%s), from ratio extracted%s, or default%s\n",
                                solids[i].label, solids[i].label, (volF) ? " (current)" : "", (massF) ? "; current" : "", (ratio) ? " (current)" : "",
                                (!volF && !massF && !ratio) ? " (current)" : "");
                            printf(" 6. Adjust amount of %s %s during fractionation / degassing (%f) unless default behavior\n", solids[i].label,
                                (volF || massF) ? "retained" : "extracted", newMin);
                        }

                        while((buf = readline("Choose: ")) != NULL) {

                            if (strlen(buf) > 0) {
                                add_history(buf);
                                if (tolower(buf[0]) == 'x') break;
                            }

                            if (!sscanf(buf, "%d", &menu) || menu < 1 || menu > 6) {
                                if (ilog) fprintf(logfile, "%s\n", buf); { MENU_ERROR(buf) }
                                free(buf);
                            }
                            else {
                                if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                                switch(menu) {
                                case 1:
                                    incLiquids = FALSE;
                                    break;
                                case 2:
                                    if ((buf = readline("Maximum number of single phase allowed (or 0 for unlimited): ")) != NULL
                                        && strlen(buf) > 0) add_history(buf);
                                    if (tolower(buf[0]) == 'x' || !sscanf(buf, "%d", &incLiquids)) incLiquids = -2; /* unchanged */
                                    if (!incLiquids) incLiquids = -1; /* unlimited */
                                    else if (incLiquids < j) incLiquids = j; /* allow at least j liquids */
                                    if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                                    break;
                                case 3:
                                    fracLiquids = FALSE;
                                    break;
                                case 4:
                                    fracLiquids = TRUE;
                                    break;
                                case 5:
                                    if (i < 0) {
                                        printf(" 1. Residual liquid measured by vol frac (phi%s), relative to system\n", (volF) ? "; current" : "");
                                        printf(" 2. Residual liquid measured by mass frac (F%s), relative to system\n", (massF) ? "; current" : "");
                                        printf(" 3. Residual liquid measured from globally set ratio extracted%s, relative to liquid\n", (ratio && !iratio) ? " (current)" : "");
                                        printf(" 4. Residual liquid measured from individual ratios extracted%s, relative to liquid\n", (iratio) ? " (current)" : "");
                                        printf(" 5. Default behavior (1.e-5 mols per phase retained%s)\n", (!volF && !massF && !ratio) ? "; current" : "");
                                    }
                                    else {
                                        printf(" 1. Residual %s measured by vol frac%s, relative to system\n", solids[i].label, (volF) ? " (current)" : "");
                                        printf(" 2. Residual %s measured by solid mass frac (1-F%s)\n", solids[i].label, (massF) ? "; current" : "");
                                        printf(" 3. Residual %s measured from globally set ratio extracted%s, relative to %s\n", solids[i].label,
                                            (ratio && !iratio) ? " (current)" : "", solids[i].label);
                                        printf(" 4. Residual %s measured from individual ratios extracted%s, relative to %s\n", solids[i].label,
                                            (iratio) ? " (current)" : "", solids[i].label);
                                        printf(" 5. Default behavior (1.e-5 mols per phase retained%s)\n", (!volF && !massF && !ratio) ? "; current" : "");
                                    }
                                    if ((buf = readline("Choose: ")) != NULL && strlen(buf) > 0) add_history(buf);
                                    if (tolower(buf[0]) == 'x' || !sscanf(buf, "%d", &minType))
                                        minType = (volF) ? 1 : ((massF) ? 2 : ((ratio && !iratio) ? 3 : ((iratio) ? 4 : 0)));
                                    switch(minType) {
                                        case 1:
                                            massF = FALSE; volF = TRUE; ratio = FALSE;
                                            break;
                                        case 2:
                                            massF = TRUE; volF = FALSE; ratio = FALSE;
                                            break;
                                        case 3:
                                            massF = FALSE; volF = FALSE; ratio = TRUE; iratio = FALSE;
                                            break;
                                        case 4:
                                            massF = FALSE; volF = FALSE; ratio = TRUE; iratio = TRUE;
                                            break;
                                        case 5:
                                            massF = FALSE; volF = FALSE; ratio = FALSE;
                                            break;
                                        default:
                                            { MENU_ERROR(buf) }
                                            break;
                                    }
                                    if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                                    break;
                                case 6:
                                    if (i < 0) {
                                        double oldMin = newMin;
                                        if      ((volF) && ((buf = readline("Vol frac (phi) of liquid to retain during extraction (relative to system): ")) != NULL
                                            && strlen(buf) > 0)) add_history(buf);
                                        else if ((massF) && ((buf = readline("Mass frac (F) of liquid to retain during extraction (relative to system): ")) != NULL
                                            && strlen(buf) > 0)) add_history(buf);
                                        else if ((ratio) && ((buf = readline("Ratio of liquid to extract (relative to total amount of liquid): ")) != NULL
                                            && strlen(buf) > 0)) add_history(buf);
                                        else break;
                                        if (tolower(buf[0]) == 'x' || !sscanf(buf, "%lf", &newMin)) newMin = oldMin;
                                    }
                                    else {
                                        double oldMin = newMin;
                                        if      ((volF) && ((buf = readline("Vol frac of phase to retain during fractionation (relative to system): ")) != NULL
                                            && strlen(buf) > 0)) add_history(buf);
                                        else if ((massF) && ((buf = readline("Mass frac (1-F) of solids to retain during extraction (relative to system): ")) != NULL
                                            && strlen(buf) > 0)) add_history(buf);
                                        else if ((ratio) && ((buf = readline("Ratio of phase to extract (relative to total amount of phase): ")) != NULL
                                            && strlen(buf) > 0)) add_history(buf);
                                        else break;
                                        if (tolower(buf[0]) == 'x' || !sscanf(buf, "%lf", &newMin)) newMin = oldMin;
                                    }
                                    if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                                    break;
                                }
                            }

                            printf("Setting to adjust (enter 'x' when done):\n");
                            if (i < 0 && j) {
                                printf(" 1. Suppress liquid phase %d %s\n", j, incLiquids ? "" : " (current)");
                                if ((incLiquids >= 0) && (incLiquids < nlc))
                                    printf(" 2. Allow at least %d liquids to join assemblage (current limit %d)\n", j, incLiquids);
                                else
                                    printf(" 2. Allow at least %d liquids to join assemblage (unlimited)\n", j);
                                printf(" 3. Do not extract liquid %d %s\n", j, fracLiquids ? "" : " (current)");
                                printf(" 4. Extract liquid %d %s\n", j, fracLiquids ? " (current)" : "");
                                printf(" 5. Change whether residual liquids are measured by vol frac (phi%s), mass frac (F%s), from ratio extracted%s, or default%s\n",
                                    (volF) ? "; current" : "", (massF) ? "; current" : "", (ratio) ? " (current)" : "", (!volF && !massF && !ratio) ? " (current)" : "");
                                printf(" 6. Adjust amount of liquid %d %s during extraction (%f) unless default behavior\n", j, (volF || massF) ? "retained" : "extracted", newMin);
                            }
                            else if (i < 0) {
                                printf(" 1. Suppress liquid phase %s\n", incLiquids ? "" : " (current)");
                                if ((incLiquids >= 0) && (incLiquids < nlc))
                                    printf(" 2. Limit number of multiple liquids allowed to join assemblage (%d)\n", incLiquids);
                                else
                                    printf(" 2. Limit number of multiple liquids allowed to join assemblage (unlimited)\n");
                                printf(" 3. Do not extract any liquids %s\n", fracLiquids ? "" : " (current)");
                                printf(" 4. Extract all liquids %s\n", fracLiquids ? " (current)" : "");
                                printf(" 5. Change whether residual liquids are measured by vol frac (phi%s), mass frac (F%s), from ratio(s) extracted%s, or default%s\n",
                                    (volF) ? "; current" : "", (massF) ? "; current" : "", (ratio) ? " (current)" : "", (!volF && !massF && !ratio) ? " (current)" : "");
                                printf(" 6. Adjust amount of liquid %s during extraction (%f) unless individual values or default behavior\n", (volF || massF) ? "retained" : "extracted", newMin);
                            }
                            else {
                                printf(" 1. Suppress %s phase %s\n", solids[i].label, incLiquids ? "" : " (current)");
                                if ((incLiquids >= 0) && (incLiquids < solids[i].na))
                                    printf(" 2. Limit number of %s phase allowed to join assemblage (%d)\n", solids[i].label, incLiquids);
                                else
                                    printf(" 2. Limit number of %s phase allowed to join assemblage (unlimited)\n", solids[i].label);
                                printf(" 3. Do not fractionate / degas %s phase %s\n", solids[i].label, fracLiquids ? "" : " (current)");
                                printf(" 4. Fractionate / degas %s phase %s\n", solids[i].label, fracLiquids ? " (current)" : "");
                                printf(" 5. Change whether residual %s is measured by %s vol frac%s, solid mass frac (1-F%s), from ratio extracted%s, or default%s\n",
                                    solids[i].label, solids[i].label, (volF) ? " (current)" : "", (massF) ? "; current" : "", (ratio) ? " (current)" : "",
                                    (!volF && !massF && !ratio) ? " (current)" : "");
                                printf(" 6. Adjust amount of %s %s during fractionation / degassing (%f) unless default behavior\n", solids[i].label,
                                    (volF || massF) ? "retained" : "extracted", newMin);
                            }
                        }

                        if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                        menu_option9(i, j, incLiquids, fracLiquids, minType, newMin);
                    }
                }
                if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
            }
            break;
        case 10:
            /* Turn phase diagram mode on / off */



            { SORRY }
            break;
        case 11:
            /* Geothermal, grids, or select PT-path file */
            { SORRY }



#ifdef NEVER_DEFINED


            /* Cannot set DPdt = 0.0 !!!! */
            if (modeFlag == GEOTHERMAL) printf("dP/dT = %f bars\n Type new dP/dT (or 0 to use current value): ", silminState->dspDPDt);
            else printf("Pinc = %f bars\nType new pressure increment: ", silminState->dspPinc);
            buffer = readline(buffer, size, lines, count, NULL); sscanf(buffer, "%lf", &newPT);
            if(newPT > 0.0) {
	if (modeFlag == GEOTHERMAL) silminState->dspDPDt = newPT;
	else silminState->dspPinc = newPT;
            }



else if (modeFlag == GEOTHERMAL) {
	  if      (silminState->isenthalpic) printf("dP/dH = %f\nType new dP/dH", silminState->dspDPDH);
	  else if (silminState->isentropic)  printf("dP/dS = %f\nType new dP/dS", silminState->dspDPDS);
	  else if (silminState->isochoric)   printf("dT/dV = %f\nType new dT/dV", 1.0 / silminState->dspDVDt);
	  printf(" (or 0 to use current value): ");
	  buffer = readline(buffer, size, lines, count, NULL);
	  scanf(buffer, "%lf", &newGrad);
	  if (newGrad != 0.0) {
	    if (silminState->isenthalpic)  silminState->dspDPDH = newGrad;
	    else if (silminState->isentropic) silminState->dspDPDS = newGrad;
	    else if (silminState->isochoric) silminState->dspDVDt = 1.0 / newGrad;
	  }
	}



if (modeFlag == PTFILE) {
	  if      (silminState->isenthalpic) printf("New P-H path filename: ");
	  else if (silminState->isentropic)  printf("New P-S path filename: ");
	  else if (silminState->isochoric)   printf("New T-V path filename: ");
	  else printf("New P-T path filename: ");
	  buffer = readline(buffer, size, lines, count, NULL);
	  silminState->PTfile  = scanfilename(silminState->PTfile, buffer);
	}

#endif


            break;
        case 12:
            /* Remixer, source mixer (was Update composition using MELTS file) */
            { SORRY }
            break;
        case 13:
            /* Read in settings via MELTS or other file (was Write out restart file */
            if (!iAmInitialized) { MELTS_ERROR }
            else {
                /*Read input file to set composition of system */
                char *buf, *filename = NULL;
                int success = FALSE, iquit = FALSE, menu;

                // Options to change other I/O settings
                printf("Type of settings file to read in or alter (enter 'x' when done):\n");
                printf(" 1. Settings in MELTS file format\n");
                printf(" 2. Trace data file\n");
                printf(" 3. Other options not yet implemented\n");
                while((buf = readline("Choose: ")) != NULL) {

                    if (strlen(buf) > 0) {
                        add_history(buf);
                        if (tolower(buf[0]) == 'x') break;
                    }

                    if (!sscanf(buf, "%d", &menu) || menu != 1) {
                        { SORRY }
                        break;
                    }
                    if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                    // Will be switch in here

                    if ((buf = readline("MELTS filename: ")) != NULL && strlen(buf) > 0) add_history(buf);
                    if (strlen(buf) < 2 && tolower(buf[0]) == 'x') iquit = TRUE;
                    else success = ((filename = scanfilename(buf)) != NULL);
                    //if (ilog) fprintf(logfile, "%s\n", buf); free(buf);

                    if (iquit) break;
                    else if (!success) { FILE_ERROR(filename) }
                    else iAmInitialized = menu_option13(filename);
                    free(filename);

                    if (!iAmInitialized) { MELTS_ERROR }
                    // else update the completer

                    printf("Type of settings file to read in or alter (enter 'x' when done):\n");
                    printf(" 1. Settings in MELTS file format\n");
                    printf(" 2. Trace data file\n");
                    printf(" 3. Other options not yet implemented\n");

                }
                if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                //if (!iquit) menu_option13(filename);
            }
            break;
        case 14:
            /* Write out MELTS file */

            /*printf("MELTS filename: ");
            scanfilename(filename);
            (void) putInputDataToFile(filename);*/

            if (!iAmInitialized) { MELTS_ERROR }
            else {
                /*Read input file to set composition of system */
                char *buf, *filename = NULL;
                int success = FALSE, iquit = FALSE, menu;


// choose to write out bulk, (one) liquid or residue (but not if no solid / liquid of course)

/*
        if((silminState->liquidMass == 0.0) && (itemp == 0)) {
            printf("WARNING: no liquid present, writing residue composition instead!\n");
            itemp = 1;
        }
        else if((silminState->solidMass == 0.0) && (itemp != 0)) {
            printf("WARNING: no solid present, writing liquid composition instead!\n");
            itemp = 0;
        }




*/





                // Options to change other I/O settings
                printf("Type of settings file to read in or alter (enter 'x' when done):\n");
                printf(" 1. Settings in MELTS file format\n");
                printf(" 2. Trace data file\n");
                printf(" 3. Other options not yet implemented\n");
                while((buf = readline("Choose: ")) != NULL) {

                    if (strlen(buf) > 0) {
                        add_history(buf);
                        if (tolower(buf[0]) == 'x') break;
                    }

                    if (!sscanf(buf, "%d", &menu) || menu != 1) {
                        { SORRY }
                        break;
                    }
                    if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                    // Will be switch in here

                    if ((buf = readline("MELTS filename: ")) != NULL && strlen(buf) > 0) add_history(buf);
                    if (strlen(buf) < 2 && tolower(buf[0]) == 'x') iquit = TRUE;
                    else success = ((filename = scanfilename(buf)) != NULL);
                    //if (ilog) fprintf(logfile, "%s\n", buf); free(buf);

                    if (iquit) break;
                    else if (!success) { FILE_ERROR(filename) }
                    else iAmInitialized = menu_option13(filename);
                    free(filename);

                    if (!iAmInitialized) { MELTS_ERROR }
                    // else update the completer

                    printf("Type of settings file to read in or alter (enter 'x' when done):\n");
                    printf(" 1. Settings in MELTS file format\n");
                    printf(" 2. Trace data file\n");
                    printf(" 3. Other options not yet implemented\n");

                }
                if (ilog) fprintf(logfile, "%s\n", buf); free(buf);
                //if (!iquit) menu_option13(filename);
            }




            break;
        case 15:
            /* Write thermodynamic output for all phases */
            { SORRY }
            break;
        case 16:
            /* Calculate integrated melt and output file(s) */
            { SORRY }
            break;
        case 17:
            /* amoeba */
            { SORRY }
            break;
        default:
            { MENU_ERROR(menubuf) }
            break;
        }

        fflush(stderr);
        fflush(stdout);
        fflush(logfile);
        printmenu(iprint);

    }

    fclose(logfile);
    return 0;
}
