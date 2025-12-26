#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "melts_gsl.h"
#include "silmin.h"
#include "adiabat.h"

#ifndef TRUE
#define TRUE  1
#endif

#define REALLOC(x, y) (((x) == NULL) ? malloc(y) : realloc((x), (y)))
#define REC   134           /* Maximal record length for input/output files */

#define FREE(x) free(x); x = NULL

#define SUCCESS 1
#define FAILURE 0

int *outsw = NULL;

/*
 *=============================================================================
 * Text tools for reading files and command line
 * See interface.c for the getInputDataFromFile
 */

char *fgetstring(char *s, int size, FILE *stream)
{
    /* modified from following url to try to catch Mac line endings
     http://www.mail-archive.com/debian-glibc@lists.debian.org/msg36705.html
     use (threadsafe) getc rather than getc_unlocked for compatibility with Windows.  Similarly for feof */
	char *p;
	int c;

	p=s,c=EOF;
	while(--size>0&&(c=getc(stream))!=EOF) { /* was getc_unlocked */
    *p++=(char)c;
    if (c=='\n'||c=='\r'||c=='\0') { /* Many do not stop at '\0'. */
            break;
    }
	}
	*p='\0'; /* Always mark the end. */
	if(c==EOF) { /* Expected. */
                    if(feof(stream)) { /* was feof_unlocked */
	    if(p==s) { /* Got nothing - not even '\0'. */
	      return (char *)NULL;
	    }
	  } else { /* Other error. */
	    return (char *)NULL;
	  }
	} /*else got something or an empty string which means null by itself*/
	return s;
}

char *scanfilename(char *line) {

    /* for now we will assume filename is big enough */
    int i, len;
    char *filename = NULL;

    len = (int) strlen(line);
    filename = (char *) malloc((size_t) (len+1) * sizeof(char));

    /* strip any leading spaces */
    for (i=0; i<len; i++) if (line[i] != ' ') break;
    ( void ) strncpy(filename, &line[i], len-i);
    filename[len-i] = '\0';
    len = strlen(filename);

    while(TRUE) {

        /* quotes match and not a false positive */
        if (((filename[0] == '"') && (filename[len-1] == '"')) || ((filename[0] == '\'') && (filename[len-1] == '\'') &&
                            strncmp(&filename[len-4], "'\\''", 4))) {
            filename[len-1] = '\0';
            for (i=0; i<len-1; i++) filename[i] = filename[i+1];
            break;
        }
        else if (filename[len-1] == ' ') {
            filename[len-1] = '\0'; /* was '\x20' when might need to call scanf again */
            len = strlen(filename);
        }
        /* not quoted and not ending in space */
        else if ((filename[0] != '"') && (filename[len-1] != '"') && (filename[0] != '\'') && (filename[len-1] != '\'')) {
            break;
        }

    }

        /* need to deal with spaces that are escaped by \? */


#ifndef WINDOWS
    int jj;
    /* we should be good except for any escaped apostrophes etc.... */
    for (i=0; i<len-3; i++) {
        if (!strncmp(&filename[i], "'\\''", 4)) for (jj=i+1; jj<len-2; jj++) filename[jj] = filename[jj+3];
        else  if (filename[i] == '\\') for (jj=i; jj<len; jj++) filename[jj] = filename[jj+1];
    }
    len = strlen(filename);
#endif

    return filename;

}

void setoutput() {

    int i;

    /* will generally be handled by subbasic_read_write */
    if (outsw == ( int * ) NULL) {
        outsw = (int *) calloc(nc, sizeof(int));
    }
    for (i=0; i<nc - 3; i++) {  /* SO3, Cl2O-1, F2O-1 */
        outsw[i] = strlen(bulkSystem[i].label);
        if (calculationMode != MODE_xMELTS) { /* i.e. calculationMode has been set */
            if ((calculationMode == MODE_pMELTS) &&  /* not calibrated */
	            (!strncmp(bulkSystem[i].label, "MnO", 3) || !strncmp(bulkSystem[i].label, "NiO", 3) ||
                !strncmp(bulkSystem[i].label, "CoO", 3) || !strncmp(bulkSystem[i].label, "CO2", 3))) outsw[i] = 0;
            else if ((calculationMode == MODE__MELTS) && !strncmp(bulkSystem[i].label, "CO2", 3)) outsw[i] = 0;
        }
    }
}

void putMultipleDataToFile(char *fileName, char *fileName2, SilminState sBlock[], int z, int writeorappend) {
    /* This is now mostly 'put one data set of multiple sets to file (except some fO2 parts) */
    int i, j, k, l, ns, nl;
    FILE *fp;
    static double *r = NULL, *m;
    static gsl_matrix *solidMasses, *liquidMasses;
    int solidCount = 0, liquidCount = 1;
    static int *solidIDs, rowIndex = 0;
    static double totalMass, mass, dspComp, temp[7];
    char *formula;
    double fo2, fo20, F, rhol, rhos, tprint, chisqr;

    //ThermoData temp;
    double fo2Delta, viscosity;

    if (r==NULL) {
        r = (double *) malloc((size_t) nc*sizeof(double));
        m = (double *) malloc((size_t) nc*sizeof(double));
        solidIDs = (int *) calloc((size_t) npc, sizeof(int));
        writeorappend = TRUE;
    }
    if (writeorappend) rowIndex = 1;
    else rowIndex++;

    //solidMasses = (double **) dmatrix(0,npc,0,z);
    solidMasses = gsl_matrix_alloc((size_t) npc, (size_t) z);
    gsl_matrix_set_zero(solidMasses);
    //liquidMasses = (double **) dmatrix(0,5,0,z);
    liquidMasses = gsl_matrix_alloc((size_t) nlc, (size_t) z);
    gsl_matrix_set_zero(liquidMasses);

    /* Old "Solids Compositions:" */
    for (i=0; i<z; i++) {
        for (j=0, sBlock[i].solidMass = 0.0; j<npc; j++) {
            for (ns=0; ns<(sBlock[i].nSolidCoexist)[j]; ns++) {
                if (solids[j].na == 1) { // ns=0
                    gsl_matrix_set(solidMasses, j, i, (sBlock[i].solidComp)[j][ns]*solids[j].mw);
                    if (solidIDs[j] == 0) {solidIDs[j] = 1; solidCount++;}
                } else {
                    for (l=0, mass = 0.0; l<solids[j].na; l++) mass += (sBlock[i].solidComp)[j+1+l][ns]*solids[j+1+l].mw;
                    gsl_matrix_set(solidMasses, j+ns, i, mass);
                    if (solidIDs[j+ns] == 0) {solidIDs[j+ns] = 1; solidCount++;}
                }
              	sBlock[i].solidMass += gsl_matrix_get(solidMasses, j+ns, i);
            }
        }
    }

    /* Make sure fileName is NULL if do not want main output */
    if(fileName == (char *) NULL) return;

    /*  for ALPHAMELTS_DOUBLE_OUTPUT use format %23.16e (?) */
    /* In older versions 'mode' was true for isobaric, isochoric or DELTAP = 0.0 */
    /* cout was true if ALPHAMELTS_CELSIUS_OUTPUT */
    i = z-1;

    /* Use the same update as silmin here. Note that volume is J/bar because it's mostly for internal calculations */
    (void) getSystemProperties(&(sBlock[i]), TRUE); /* TRUE for derivatives */

    //tprint = (getenv("ALPHAMELTS_CELSIUS_OUTPUT") != NULL) ? sBlock[i].T - 273.15 : sBlock[i].T;
    tprint = sBlock[i].T - 273.15;
    totalMass = sBlock[i].liquidMass + sBlock[i].solidMass;

    /* Old "Liquid Thermodynamic Data:" */
    fo2Delta = silminState->fo2Delta;

    fo2 = sBlock[i].fo2;
    silminState->fo2Delta = sBlock[i].fo2Delta;
    if (((sBlock[i].fo2Path != FO2_NONE) || (sBlock[0].fo2Path != FO2_NONE)) && !((sBlock[i].fo2Path != FO2_NONE) && (sBlock[0].fo2Path != FO2_NONE)))
        fo2 -= getlog10fo2(sBlock[i].T, sBlock[i].P, sBlock[i].fo2Path);
    else if ((sBlock[i].fo2Path == FO2_NONE) && (sBlock[0].fo2Path == FO2_NONE))
        fo2 -= getlog10fo2(sBlock[i].T, sBlock[i].P, FO2_QFM);

    fo20 = sBlock[i].fo2;
    silminState->fo2Delta = sBlock[0].fo2Delta;
    if (sBlock[0].fo2Path != FO2_NONE) fo20 -= getlog10fo2(sBlock[i].T,sBlock[i].P,sBlock[0].fo2Path);

    silminState->fo2Delta = fo2Delta;

    F = sBlock[i].liquidMass/totalMass;
    rhol = (sBlock[i].liquidMass > 0.0) ? 0.1*sBlock[i].liquidMass/sBlock[i].liquidTD.v : 0.0;
    rhos = (sBlock[i].solidMass > 0.0) ? 0.1*sBlock[i].solidMass/sBlock[i].solidTD.v : 0.0;
    /* Amoeba can't be called with liquid fractionation */
    //chisqr = (sBlock[i].liquidout < 0.0) ? -sBlock[i].liquidout : 0.0;
    chisqr = 0.0;

    /* Old "Liquid Compositions:" */
    for (i=0; i<z; i++) {
        if (sBlock[i].liquidMass != 0.0) {
            if (sBlock[i].nLiquidCoexist > liquidCount) liquidCount = sBlock[i].nLiquidCoexist;
            for (nl=0; nl<sBlock[i].nLiquidCoexist; nl++) {
            	for (j=0, mass = 0.0; j<nlc; j++) {
	                //	  liquidMasses[ns][i] += (sBlock[i].liquidComp)[ns][k]*liquid[k].mw;
	                for (k=0; k<nc; k++) mass += (sBlock[i].liquidComp)[nl][j]*(liquid[j].liqToOx)[k]*bulkSystem[k].mw;
            	}
                gsl_matrix_set(liquidMasses, nl, i, mass);
            }
        }
    }

    i = z-1;
    if ((fp = fopen("System_main_tbl.txt", (writeorappend ? "w" : "a"))) == NULL) {
        printf("Output error -- can't open file.\n");
        return;
    }
    else {

        if (writeorappend) {

            /* Put H, S, V (instead of S, H, V) so consistent with .h files etc. */
            if((silminInputData.title != NULL) && (int) strlen(silminInputData.title) > 0)
	            fprintf (fp, "Title:%s\n\n", silminInputData.title);
            fprintf (fp, "System Thermodynamic Data:\n");

            /* 10^6 refers to V in cc */
            fprintf(fp, "index Pressure Temperature mass F phi H S V Cp dVdP*10^6 dVdT*10^6 fO2");

            /* If current and first state are both buffered then print absolute fO2 once */
            if (((sBlock[i].fo2Path != FO2_NONE) || (sBlock[0].fo2Path != FO2_NONE)) && !((sBlock[i].fo2Path != FO2_NONE) && (sBlock[0].fo2Path != FO2_NONE))) {
                if (sBlock[i].fo2Path == FO2_QFM) fprintf(fp, "-(QFM");
                if (sBlock[i].fo2Path == FO2_NNO) fprintf(fp, "-(NNO");
                if (sBlock[i].fo2Path == FO2_IW) fprintf(fp, "-(IW");
                if (sBlock[i].fo2Path == FO2_HM) fprintf(fp, "-(HM");
                if(sBlock[i].fo2Delta > 0.0) fprintf(fp, "+%.1f) fO2", sBlock[i].fo2Delta);
                else if(sBlock[i].fo2Delta < 0.0) fprintf(fp, "%.1f) fO2", sBlock[i].fo2Delta);
                else fprintf(fp, ") fO2");
            }
            /* put relative to QFM if both are unbuffered */
            else if ((sBlock[i].fo2Path == FO2_NONE) && (sBlock[0].fo2Path == FO2_NONE)) {
                fprintf(fp, "-(QFM) fO2");
            }
            else fprintf(fp, "(absolute) fO2");

            /* Have option to change buffer during calculation so print relative to first state too */
            if (sBlock[0].fo2Path != FO2_NONE) {
                if (sBlock[0].fo2Path == FO2_QFM) fprintf(fp, "-(QFM");
                if (sBlock[0].fo2Path == FO2_NNO) fprintf(fp, "-(NNO");
                if (sBlock[0].fo2Path == FO2_IW) fprintf(fp, "-(IW");
                if (sBlock[0].fo2Path == FO2_HM) fprintf(fp, "-(HM");
                if(sBlock[0].fo2Delta > 0.0) fprintf(fp, "+%.1f)", sBlock[0].fo2Delta);
                else if(sBlock[0].fo2Delta < 0.0) fprintf(fp, "%.1f)", sBlock[0].fo2Delta);
                else fprintf(fp, ")");
            } else fprintf(fp, "(absolute)");

            fprintf(fp, " rhol rhos viscosity aH2O chisqr\n");

        }

        // MAKE THIS SYSTEM VISCOSITY?
        viscosity = 0.0; /* just the first liquid for now */
        if (sBlock[i].liquidMass != 0.0) {
            (void) getDspLiquidProperties(&sBlock[i], 0, temp);
            viscosity = temp[6];
        }
        // viscosity - 2.0*log10(1.0-2.0*totalVolume/(totalVolume+vLiq)

        fprintf(fp, "%d %.2f %.2f %.6f %.16f %.6f %9.6f %f %.6f %.6f %.6f %.6f %.3f %.3f %.6f %.6f %.3f n/a n/a\n",
            rowIndex, sBlock[i].P, tprint, totalMass, F, (sBlock[i].liquidTD.v)/(sBlock[i].bulkTD.v),
            sBlock[i].bulkTD.h, sBlock[i].bulkTD.s, sBlock[i].bulkTD.v * 10.0, sBlock[i].bulkTD.cp,
            sBlock[i].bulkTD.dvdp * 1.0e7, sBlock[i].bulkTD.dvdt * 1.0e7,
            fo2, fo20, rhol, rhos, viscosity); //, sBlock[i].aH2O, chisqr);  %.6g %.3f

        fclose(fp);

    }

    if ((fp = fopen("Liquid_comp_tbl.txt", (writeorappend ? "w" : "a"))) == NULL) {
        printf("Output error -- can't open file.\n");
        return;
    }
    else {

        if (writeorappend) {

            if((silminInputData.title != NULL) && (int) strlen(silminInputData.title) > 0) fprintf (fp, "Title:%s\n\n", silminInputData.title);
            fprintf(fp, "Liquid Composition:\n");

            fprintf(fp, "index Pressure Temperature mass ");
            for (k=0; k<nc; k++) if (outsw[k]) fprintf(fp, "%s ", bulkSystem[k].label);

        }

        if (sBlock[i].liquidMass != 0.0) {

            fprintf(fp, "\n%d %.2f %.2f %.6f", rowIndex, sBlock[i].P, tprint, sBlock[i].liquidMass);

            for (k=0; k<nc; k++) if (outsw[k]) {
            	for (j=0, dspComp=0.0; j<nlc; j++)
            	    for (nl=0; nl<sBlock[i].nLiquidCoexist; nl++) dspComp += sBlock[i].liquidComp[nl][j] * liquid[j].liqToOx[k];
                dspComp *= bulkSystem[k].mw;
                fprintf(fp, " %.6g", dspComp*100.0/sBlock[i].liquidMass);
            }

        }
        else fprintf(fp, "\n%d %.2f %.2f %.6f ---", rowIndex, sBlock[i].P, tprint, 0.0);

        fclose(fp);

    }

    if ((fp = fopen("Phase_main_tbl.txt", (writeorappend ? "w" : "a"))) == NULL) {
        printf("Output error -- can't open file.\n");
        return;
    }
    else {

        if (writeorappend) {
            if((silminInputData.title != NULL) && (int) strlen(silminInputData.title) > 0) fprintf (fp, "Title:%s\n\n", silminInputData.title);
        }

        fprintf(fp, "index %d Pressure %.2f Temperature %.2f", rowIndex, sBlock[i].P, tprint);
        for (k=0; k<nc; k++) if (outsw[k]) fprintf(fp, " %s", bulkSystem[k].label);
        fprintf(fp, "\n");

        /* Probably not actually necessary */
        //updateThermoProp(&sBlock[i]);

        if (sBlock[i].liquidMass != 0.0) {
            for (nl=0; nl<sBlock[i].nLiquidCoexist; nl++) {

            	//ThermoData temp;
            	//double viscosity;
                (void) getDspLiquidProperties(&sBlock[i], nl, temp);

                if ((silminState->txtOutput == TEXT_ALPHA_0) || (silminState->txtOutput == TEXT_BOTH_0))
                    fprintf(fp, "liquid_%d %.6f %9.6f %f %.6f %.6f %.3f", nl, gsl_matrix_get(liquidMasses, nl, i), temp[1], temp[2], temp[3], temp[4], temp[6]);
                else
                    fprintf(fp, "liquid%d %.6f %9.6f %f %.6f %.6f %.3f", (nl+1), gsl_matrix_get(liquidMasses, nl, i), temp[1], temp[2], temp[3], temp[4], temp[6]);

                for (k=0; k<nc; k++) if (outsw[k]) {
                    for (j=0, dspComp=0.0; j<nlc; j++) dspComp += sBlock[i].liquidComp[nl][j]*liquid[j].liqToOx[k]*bulkSystem[k].mw;
                    fprintf(fp, " %.6g", dspComp*100.0/gsl_matrix_get(liquidMasses, nl, i));
                }
                fprintf(fp, "\n");
            }
        }

        for (j=0; j<npc; j++) {
            for (ns=0; ns<(sBlock[i].nSolidCoexist)[j]; ns++) {

                //ThermoData temp;
                (void) getDspSolidProperties(&sBlock[i], j, ns, temp);

                if (solids[j].na == 1) {

                    if ((silminState->txtOutput == TEXT_ALPHA_0) || (silminState->txtOutput == TEXT_BOTH_0))
                        fprintf(fp,"%s_%d %f %9.6f %f %.6f %.6f %s", solids[j].label, ns, gsl_matrix_get(solidMasses, j, i), temp[1], temp[2], temp[3], temp[4], solids[j].formula);
                    else
                        fprintf(fp,"%s%d %f %9.6f %f %.6f %.6f %s", solids[j].label, (ns+1), gsl_matrix_get(solidMasses, j, i), temp[1], temp[2], temp[3], temp[4], solids[j].formula);

                    for (k=0; k<nc; k++) if (outsw[k]) {
                        dspComp = (solids[j].solToOx[k])*bulkSystem[k].mw*(silminState->solidComp)[j][0];
                        fprintf(fp, " %g", dspComp*100.0/gsl_matrix_get(solidMasses, j, i));
                    }
                    fprintf(fp, "\n");

                } else {

                    for (l=0; l<solids[j].na; l++) m[l] = (sBlock[i].solidComp)[j+1+l][ns];
                    (*solids[j].convert)(SECOND, THIRD, sBlock[i].T, sBlock[i].P, NULL, m, r, NULL, NULL, NULL, NULL, NULL);
                    (*solids[j].display)(FIRST, sBlock[i].T, sBlock[i].P, r, &formula);
                    if ((silminState->txtOutput == TEXT_ALPHA_0) || (silminState->txtOutput == TEXT_BOTH_0))
                        fprintf(fp,"%s_%d %.6f %9.6f %f %.6f %.6f %s", solids[j].label, ns, gsl_matrix_get(solidMasses, j+ns, i), temp[1], temp[2], temp[3], temp[4], formula);
                    else
                        fprintf(fp,"%s%d %.6f %9.6f %f %.6f %.6f %s", solids[j].label, (ns+1), gsl_matrix_get(solidMasses, j+ns, i), temp[1], temp[2], temp[3], temp[4], formula);
                    free(formula);

                    for (k=0; k<nc; k++) if (outsw[k]) {
                        for (l=0, dspComp=0.0; l<solids[j].na; l++) dspComp += m[l]*(solids[j+1+l].solToOx[k])*bulkSystem[k].mw;
                        fprintf(fp, " %g", dspComp*100.0/gsl_matrix_get(solidMasses, j+ns, i));
                    }
                    fprintf(fp, "\n");
                }

            }
        }

        fclose(fp);

    }

    i = z-1;
    if ((fp = fopen("Solid_comp_tbl.txt", (writeorappend ? "w" : "a"))) == NULL) {
        printf("Output error -- can't open file.\n");
        return;
    }
    else {

        if (writeorappend) {

            if((silminInputData.title != NULL) && (int) strlen(silminInputData.title) > 0) fprintf (fp, "Title:%s\n\n", silminInputData.title);
            fprintf(fp, "Solid Composition:\n");
            fprintf(fp, "index Pressure Temperature mass ");
            for (k=0; k<nc; k++) if (outsw[k]) fprintf(fp, "%s ", bulkSystem[k].label);

        }

        if (sBlock[i].solidMass != 0.0) {

            fprintf(fp, "\n%d %.2f %.2f %.6f", rowIndex, sBlock[i].P, tprint, sBlock[i].solidMass);

            for (k=0; k<nc; k++) if (outsw[k]) {
	            dspComp=sBlock[i].bulkComp[k];
	            if (sBlock[i].liquidMass != 0.0) {
	                for (j=0; j<nlc; j++)
	                    for (nl=0; nl<sBlock[i].nLiquidCoexist; nl++) dspComp -= sBlock[i].liquidComp[nl][j] * liquid[j].liqToOx[k];
                }
                dspComp *= bulkSystem[k].mw;
                dspComp *= 100.0/sBlock[i].solidMass;
                fprintf(fp, " %.6g", (dspComp > 1.0e-8) ? dspComp : 0.0);
            }

        }
        else fprintf(fp, "\n%d %.2f %.2f %.6f ---", rowIndex, sBlock[i].P, tprint, 0.0);

        fclose(fp);

    }

    if ((fp = fopen("Bulk_comp_tbl.txt", (writeorappend ? "w" : "a"))) == NULL) {
        printf("Output error -- can't open file.\n");
        return;
    }
    else {

        if (writeorappend) {

            if((silminInputData.title != NULL) && (int) strlen(silminInputData.title) > 0) fprintf (fp, "Title:%s\n\n", silminInputData.title);
            fprintf(fp, "Bulk Composition:\n");
            fprintf(fp, "index Pressure Temperature mass ");
            for (k=0; k<nc; k++) if (outsw[k]) fprintf(fp, "%s ", bulkSystem[k].label);

        }

        fprintf(fp, "\n%d %.2f %.2f %.6f", rowIndex, sBlock[i].P, tprint, totalMass);

        for (k=0; k<nc; k++) if (outsw[k]) {
            dspComp = sBlock[i].bulkComp[k]*bulkSystem[k].mw;
            fprintf(fp, " %.6g", dspComp*100.0/totalMass);
        }

        fclose(fp);
        printf("Current state of the system recorded in _tbl.txt files.\n");

    }

    gsl_matrix_free(solidMasses);
    gsl_matrix_free(liquidMasses);

}
