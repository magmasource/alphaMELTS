
/* Version of cubicquad for use with release version of alphaMELTS */
/* Uses information stored in memory, not in a file */
/* PDA 06/29/04 */

/*#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "recipes.h"
#include "silmin.h"
#include "adiabat.h"
//#include "trace.h"

#define REC2   4000 /* Maximal record length for input/output files */

void CQint(double x[], double y[], int n, double ixdy[], int sign);
void Dumbint(double x[], double y[], int n, double ixdy[], int sign);
void spline(double x[], double y[], int n, double yp1, double ypn, 
double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, 
double *y);
double qgaus(double a, double b, double xa[], double ya[], double 
y2a[], int n);

int doingTraces = FALSE;

/* n is passed from main and is the number of calculations in the extracted array
   nmajors is the number of major components, ntrace is the number of trace components */
int cubicquad(double **extractions, int n, int nmajors, int ntrace, int *outsw, char **outlb) {
  int i, j, nox, cubit = 1, itemp, Pint = 0, Fint = 0, Phi = (getenv("ALPHAMELTS_INTEGRATE_PHI") != NULL);
  char filename1[300], filename2[300], line[REC2];
  FILE *fp1;
  double *P, *T, *dF, **comp, *F, FB, ZC, PC, *y2;
  double Fmax, PC1D, ZC1D, dFnew;
  double Pbase_2, *compbase_2, *compbase_1;
  double *iFdP, **icompdF, *iPdF;
  double **icompFdP,  *iPFdP;

   
   nox = nmajors;
   if (doingTraces) nox += ntrace;

   if (Phi) {
     cubit = 0;
   }
   else {
     printf("Calculate pressure integral (1) or aggregated melts only (0)? ");
     scanf("%d", &Pint);

     if (Pint) {
       if (!strcmp(getenv("ALPHAMELTS_RUN_MODE"), "isobaric") || !strcmp(getenv("ALPHAMELTS_RUN_MODE"), "isochoric")) {
	 printf("WARNING: pressure integral is not valid for %s calculations.\n", getenv("ALPHAMELTS_RUN_MODE"));
	 Pint = 0;
       }
       else if (Pint) {
	 printf("Use simple rectangle rule (1) or Cubic Gaussian Quadrature (0)? ");
	 scanf("%d", &cubit);
       }
     }
   }
   
   if (!cubit) {
     printf("Use equally spaced melt fractions (1), or do not interpolate (0)? ");
     scanf("%d", &Fint);
     if (Phi) cubit = !Fint;
   }

   if (Fint) {
     printf("Constant %s relative to system %s (1), or starting %s (0)? ",
	    ((Phi) ? "dPhi" : "dF"), ((Phi) ? "volume" : "mass"), ((Phi) ? "volume" : "mass"));
     scanf("%d", &Fint);
     if (!Fint) Fint = -1;
     printf("%s to use for interpolation: ", ((Phi) ? "dPhi" : "dF"));
     scanf("%lf", &dFnew);
   }
   
   printf("Major element filename: ");
   scanfilename(filename1);
   if (doingTraces) {
     printf("Trace element filename: ");
     scanfilename(filename2);
   }
 
   if(n) {

     P = dvector(1,n+1); T = dvector(1,n+1); comp = dmatrix(1, nox+1, 1,n+1);
     dF = dvector(1,n+1); F = dvector(1,n+1);
     
     P[1] = extractions[0][0];
     T[1] = extractions[0][1];
     dF[1] = extractions[0][2];
     if (getenv("ALPHAMELTS_CELSIUS_OUTPUT") != NULL) T[1] -= 273.15;

     if(Pint && dF[1] != 0.0) {
       double Ptemp = P[1];

       P[1] -= (DELTAP < 0.0) ? DELTAP : -1000;
       F[1] = dF[1] = 0.0;

       printf("WARNING: Calculation did not start subsolidus!\n");
       if (DELTAP >= 0.0) printf("WARNING: DELTAP is not negative; using default value of -1000 bars.\n");
       printf("Assuming solidus is between %f and %f for melt extraction.\n", Ptemp, P[1]);
   
     }
     else F[1] = dF[1];
     for (j=1; j<=nox; j++) comp[j][1] = extractions[0][j+2];

     for (i=2;i<=n;i++) {
       P[i] = extractions[i-1][0];
       T[i] = extractions[i-1][1];
       if (getenv("ALPHAMELTS_CELSIUS_OUTPUT") != NULL) T[i] -= 273.15;
       dF[i] = extractions[i-1][2];
       for (j=1; j<=nox; j++) comp[j][i] = extractions[i-1][j+2];
       F[i] = F[i-1]+dF[i];
     }

   }
   else if (getenv("ALPHAMELTS_INTEGRATE_FILE") != NULL) {
      if ((fp1 = fopen(getenv("ALPHAMELTS_INTEGRATE_FILE"), "r")) != NULL) {
 
	if (fgets(line, REC2, fp1) != NULL && strlen(line) > 1) {
	  double Ptemp, Ttemp, Ftemp;
	  sscanf(strtok(line, ",\t "), "%lf", &Ptemp);
	  sscanf(strtok(NULL, ",\t "), "%lf", &Ttemp);
	  sscanf(strtok(NULL, ",\t "), "%lf", &Ftemp);
	  n = (Ftemp != 0.0) ? 1 : 0;
	  if(Pint) n++;
	}

	while (fgets(line, REC2, fp1) != NULL && strlen(line) > 1) {
	  double Ptemp, Ttemp, Ftemp;
	  sscanf(strtok(line, ",\t "), "%lf", &Ptemp);
	  sscanf(strtok(NULL, ",\t "), "%lf", &Ttemp);
	  sscanf(strtok(NULL, ",\t "), "%lf", &Ftemp);
	  if(Pint || (Ftemp != 0.0)) n++;
	}
	fclose(fp1);

	P = dvector(1,n+1); T = dvector(1,n+1); comp = dmatrix(1, nox+1, 1,n+1);
	dF = dvector(1,n+1); F = dvector(1,n+1);
	
	if ((fp1 = fopen(getenv("ALPHAMELTS_INTEGRATE_FILE"), "r")) != NULL) {

	  i = 1;
	  while (i<=n) {

	    fgets(line, REC2, fp1);
	    sscanf(strtok(line, ",\t "), "%lf", &P[i]);
	    sscanf(strtok(NULL, ",\t "), "%lf", &T[i]);
	    sscanf(strtok(NULL, ",\t "), "%lf", &dF[i]);
	    if (getenv("ALPHAMELTS_CELSIUS_OUTPUT") != NULL) T[i] -= 273.15;
	    
	    if(Pint || (dF[i] != 0.0)) {
	      
	      for (j=1; j<=nox; j++) sscanf(strtok(NULL, ",\t "), "%lf", &comp[j][i]);
	      
	      if((i == 1) && Pint && (dF[i] != 0.0)) {
		double Ptemp = P[1];

		i = 3;
		P[2] = P[1]; T[2] = T[1];
		dF[2] = dF[1]; F[2] = dF[2];
		for (j=1; j<=nox; j++) comp[j][2] = comp[j][1];

		P[1] -= (DELTAP < 0.0) ? DELTAP : -1000;
		F[1] = dF[1] = 0.0;

		printf("WARNING: Calculation did not start subsolidus!\n");
		if (DELTAP >= 0.0) printf("WARNING: DELTAP is not negative; assuming a value of -1000 bars!\n");
		printf("Assuming solidus is between %f and %f for melt extraction\n", Ptemp, P[1]);
 
	      }
	      else {
		F[i] = (i > 1) ? F[i-1]+dF[i] : dF[i];
		i++;
	      }
	    }
	  }
	  fclose(fp1);
	
	}

      }
      else {
	printf("Cubicquad could not open integrate file.\n");
	printf("Printing output file headers only.\n");
      }
   }
    
   if(n) {

     PC1D = P[1];

     if (Fint) {

       double *P0 = P, *T0 = T, **comp0 = comp, *F0 = F;
       int n0 = n;

       n = ( int ) floor((F0[n0] - F0[1]) / dFnew);
       
       P = dvector(1,n+1); T = dvector(1,n+1); comp = dmatrix(1, nox+1, 1,n+1);
       F = dvector(1,n+1); y2 = dvector(1,n0);

       F[1] = F0[1];
       for (i=2;i<=n;i++) {
	 if (Fint <= 0) F[i] = F[i-1]+dFnew;
	 else F[i] = F[i-1]+dFnew*(1.0 - F[i-1]);
       }

       P[1] = P0[1];
       spline(F0, P0, n0, 1e32, 1e32, y2);
       for (i=2;i<=n;i++) splint(F0, P0, y2, n0, F[i], &P[i]);

       T[1] = T0[1];
       spline(F0, T0, n0, 1e32, 1e32, y2);
       for (i=2;i<=n;i++) splint(F0, T0, y2, n0, F[i], &T[i]);
       
       for (j=1; j<=nox; j++) {
	 comp[j][1] = comp0[j][1];
	 spline(F0, comp0[j], n0, 1e32, 1e32, y2);
	 for (i=2;i<=n;i++) splint(F0, comp0[j], y2, n0, F[i], &comp[j][i]); 
       }
       
       free_dvector(P0,1,n+1); free_dvector(T0,1,n+1); free_dmatrix(comp0, 1, nox+1, 1,n+1);
       free_dvector(F0,1,n+1); free_dvector(y2,1,n0);

     }

     /* Get Integral CompdF */
     icompdF = dmatrix(1, nox, 1,n+2);
     for (j=1; j<=nox; j++) {
       icompdF[j][1] = F[1]*comp[j][1]; /* if Pint this will be 0.0 anyway, as before */
       if (!cubit) CQint(F,comp[j],n,icompdF[j],1);
       else Dumbint(F,comp[j],n,icompdF[j],1);
     }

     Pint = (Pint && !Phi);
     
     if(Pint) {

       y2 = dvector(1,n);

       /* Get Integral FdP */
       iFdP = dvector(1,n+2); iFdP[1] = 0.0;
       if (!cubit) CQint(P,F,n,iFdP,-1);
       else Dumbint(P,F,n,iFdP,-1);

       /* Get Integral PdF */
       iPdF = dvector(1,n+2); iPdF[1] = 0.0;
       if (!cubit) CQint(F,P,n,iPdF,1);
       else Dumbint(F,P,n,iPdF,1);

       /* Get Integral <Comp>FdP, Integral <P>FdP */
       icompFdP = dmatrix(1, nox, 1,n);
       for (j=1; j<=nox; j++) {
	 icompFdP[j][1] = 0.0;
	 if (!cubit) CQint(P,icompdF[j],n,icompFdP[j],-1);
	 else Dumbint(P,icompdF[j],n,icompFdP[j],-1);
       }
       iPFdP = dvector(1,n); iPFdP[1] = 0.0;
       if (!cubit) CQint(P,iPdF,n,iPFdP,-1);
       else Dumbint(P,iPdF,n,iPFdP,-1);

       /* Get P1D */
       iPdF[1] = P[1];
       for (i=2;i<=n;i++) {
	 if (F[i] != 0.0) iPdF[i] /= F[i];
       }

       /* Get <Comp>2D */
       for (j=1; j<=nox; j++) icompFdP[j][1] = comp[j][1];
       iPFdP[1] = P[1];
       for (i=2;i<=n;i++) {
	 if (iFdP[i] != 0.0) {
	   for (j=1; j<=nox; j++) icompFdP[j][i] /= iFdP[i];
	   iPFdP[i] /= iFdP[i];
	 }
       }


       /* Pick Base of Crust 2-D */
       PC = P[1];
       for (i=2;i<=n;i++) {
	 if (iFdP[i] == P[i]) break;
	 if (iFdP[i] > P[i]) {
	   PC = (iFdP[i-1]*P[i] - iFdP[i]*P[i-1])/((P[i]-P[i-1])-(iFdP[i]-iFdP[i-1]));
	   break;
	 }
       }

       compbase_2 = dvector(1, nox);
       if(PC == P[1]) {
	 printf("WARNING: calculation did not get to base of crust (2-D)\n");
	 printf("Writing final aggregate composition instead\n");
	 PC = FB = ZC = Pbase_2 = 0.0;
	 for (j=1; j<=nox; j++) compbase_2[j] = icompFdP[j][n];
       }
       else {
	 FB = PC/(P[1]-PC);
	 ZC = PC/1000*10.2/2.62;

	 spline(P, iPFdP, n, 1e32, 1e32, y2);
	 splint(P, iPFdP, y2, n, PC, &Pbase_2);

	 for (j=1; j<=nox; j++) {
	   spline(P, icompFdP[j], n, 1e32, 1e32, y2);
	   splint(P, icompFdP[j], y2, n, PC, &compbase_2[j]);
	 }
       }

       /* Pick base of crust 1-D */
       for (i=2; i<=n; i++) {
	 if (P[i] == F[i]*P[1]/(1.0+F[i])) break;
	 if (P[i] <  F[i]*P[1]/(1.0+F[i])) {
	   PC1D = (F[i-1]*P[1]/(1.0+F[i-1])*P[i] - F[i]*P[1]/(1.0+F[i])*P[i-1])/
	     ((P[i]-P[i-1]) - (F[i]*P[1]/(1.0+F[i]) - F[i-1]*P[1]/(1.0+F[i-1])));
	   break;
	 }
       }

     } /* end of Pint */

     /* Get <Comp>1D */
     for (j=1; j<=nox; j++) icompdF[j][1] = comp[j][1];
     for (i=2;i<=n;i++) {
       if (F[i] != 0.0)
	 for (j=1; j<=nox; j++) icompdF[j][i] /= F[i];
     }

     compbase_1 = dvector(1, nox);
     if (PC1D == P[1]) {
       if(Pint && !Phi) {
	 printf("WARNING: calculation did not get to base of crust (1-D)\n");
	 printf("Writing final aggregate composition instead\n");
       }
       PC1D = Fmax = ZC1D = 0.0;       
       for (j=1; j<=nox; j++) compbase_1[j] = icompdF[j][n];
     }
     else {
       ZC1D = PC1D/1000*10.2/2.62;
       
       spline(P, F, n, 1e32, 1e32, y2);
       splint(P, F, y2, n, PC1D, &Fmax);
       
       for (j=1; j<=nox; j++) {
	 spline(P, icompdF[j], n, 1e32, 1e32, y2);
	 splint(P, icompdF[j], y2, n, PC1D, &compbase_1[j]);
       }
     }
     
     if(Pint) free_dvector(y2,1,n);

   }
   else { /* no extractions - just print out the headers */
     PC = FB = ZC = Pbase_2 = 0.0;
     compbase_2 = dvector(1, nox);
     for (j=1; j<=nox; j++) compbase_2[j] = 0.0;
     PC1D = Fmax = ZC1D = 0.0;
     compbase_1 = dvector(1, nox);
     for (j=1; j<=nox; j++) compbase_1[j] = 0.0;
   }
     
   /* output */
   if ((fp1 = fopen(filename1, "w")) == NULL) {
     printf("Cubicquad could not open major element file\n");
     return -1;
   }
   if((silminInputData.title != NULL) && (int) strlen(silminInputData.title) > 0) {
     fprintf (fp1, "Title:%s\n\n", silminInputData.title);
   }
   if(Pint) {
     fprintf(fp1, "Bulk Crust Properties:\nPc FBbase Zcbase <P>base Pc1D Fmax1D Zc1D");
     for (j=0; j<nc; j++) if (outsw[j]) fprintf(fp1, " %s_2D", bulkSystem[j].label);
     for (j=0; j<nc; j++) if (outsw[j]) fprintf(fp1, " %s_1D", bulkSystem[j].label);
     fprintf(fp1, "\n");
     fprintf(fp1, "%g %g %g %g %g %g %g",PC,FB,ZC,Pbase_2, PC1D, Fmax, ZC1D);
     for (j=1; j<=nmajors; j++) fprintf(fp1, " %g", compbase_2[j]);
     for (j=1; j<=nmajors; j++) fprintf(fp1, " %g", compbase_1[j]);
     fprintf(fp1, "\n\n");
     fprintf(fp1, "Aggregated Liquid:\nP T F P2D P1D iFdP");
     for (j=0; j<nc; j++) if (outsw[j]) fprintf(fp1, " %s_2D", bulkSystem[j].label);
     for (j=0; j<nc; j++) if (outsw[j]) fprintf(fp1, " %s_1D", bulkSystem[j].label);
   }
   else {
     fprintf(fp1, "P T F");
     for (j=0; j<nc; j++) if (outsw[j]) fprintf(fp1, " %s", bulkSystem[j].label);
   }
   fprintf(fp1, "\n");
   if(n) {
     for (i=1;i<=n;i++) {
       fprintf(fp1, "%.2f %.2f %g", P[i], T[i], F[i]);
       if (Pint) {
	 fprintf(fp1, " %g %g %g", iPFdP[i], iPdF[i], iFdP[i]);
	 for (j=1; j<=nmajors; j++) fprintf(fp1, " %g", icompFdP[j][i]);
       }
       for (j=1; j<=nmajors; j++) fprintf(fp1, " %g", icompdF[j][i]);
       fprintf(fp1, "\n");
     }
   }
   fclose(fp1);

   if (doingTraces) {
     if ((fp1 = fopen(filename2, "w")) == NULL) {
       printf("Cubicquad could not open trace element file\n");
       return -1;
     }
     if((silminInputData.title != NULL) && (int) strlen(silminInputData.title) > 0) {
       fprintf (fp1, "Title:%s\n\n", silminInputData.title);
     }
     if(Pint) {
       fprintf(fp1, "Bulk Crust Properties:\nPc FBbase Zcbase <P>base Pc1D Fmax1D Zc1D");
       for (j=0; j<ntrace; j++)  fprintf(fp1, " %s_2D", outlb[j]);
       for (j=0; j<ntrace; j++)  fprintf(fp1, " %s_1D", outlb[j]);
       fprintf(fp1, "\n");
       fprintf(fp1, "%g %g %g %g %g %g %g",PC,FB,ZC,Pbase_2, PC1D, Fmax, ZC1D);
       for (j=nmajors+1; j<=nox; j++) fprintf(fp1, " %g", compbase_2[j]);
       for (j=nmajors+1; j<=nox; j++) fprintf(fp1, " %g", compbase_1[j]);
       fprintf(fp1, "\n\n");
       fprintf(fp1, "Aggregated Liquid:\nP T F P2D P1D iFdP");
       for (j=0; j<ntrace; j++)  fprintf(fp1, " %s_2D", outlb[j]);
       for (j=0; j<ntrace; j++)  fprintf(fp1, " %s_1D", outlb[j]);
     }
     else {
       fprintf(fp1, "P T F");
       for (j=0; j<ntrace; j++) fprintf(fp1, " %s", outlb[j]);
     }
     fprintf(fp1, "\n");
     if(n) {
       for (i=1;i<=n;i++) {
	 fprintf(fp1, "%.2f %.2f %g", P[i], T[i], F[i]);
	 if (Pint) {
	   fprintf(fp1, " %g %g %g", iPFdP[i], iPdF[i], iFdP[i]);
	   for (j=nmajors+1; j<=nox; j++) fprintf(fp1, " %g", icompFdP[j][i]);
	 }
	 for (j=nmajors+1; j<=nox; j++) fprintf(fp1, " %g", icompdF[j][i]);
	 fprintf(fp1, "\n");
       }
     }
     fclose(fp1);
   }

   if(n) {
     
     printf("Write MELTS files (1) or return to menu or (0)? ");   
     scanf("%d", &itemp);

     if(itemp) {

       if(Pint) {
	 printf("Filename for 2-D aggregate melt: ");
	 scanfilename(filename1);
	 printf("Filename for 1-D aggregate melt: ");
	 scanfilename(filename2);

	 if ((fp1 = fopen(filename1, "w")) == NULL) {
	   printf("Cubicquad could not open 2-D melts file\n");
	   return -1;
	 }
	 if((silminInputData.title != NULL) && (int) strlen(silminInputData.title) > 0) {
	   fprintf(fp1, "Title:%s 2-D aggregated melt\n", silminInputData.title);
	 }
	 else {
	   fprintf(fp1, "Title:2-D aggregated melt\n");
	 }   

	 for (i=1, j=0; j<nc; j++) {
	   if (outsw[j])
	     fprintf(fp1, "Initial Composition: %s %.6f\n", bulkSystem[j].label, compbase_2[i++]);
	 }
	 if (doingTraces) {
	   for (j=0; j<ntrace; j++) {
	     if(!strcmp(outlb[j], "H2O") && (getenv("ALPHAMELTS_DO_TRACE_H2O") != NULL)) {
	       fprintf(fp1, "Initial Trace: H2O 0.0\n");
	       i++;
	     }
	     else {
	       fprintf(fp1, "Initial Trace: %s %.6f\n", outlb[j], compbase_2[i++]);
	     }
	   }
	 }
	 fclose(fp1);
       }
       else {
	 printf("Filename for aggregate melt: ");
	 scanfilename(filename2);
       }

       if ((fp1 = fopen(filename2, "w")) == NULL) {
	 printf("Cubicquad could not open 1-D melts file\n");
	 return -1;
       }
       if((silminInputData.title != NULL) && (int) strlen(silminInputData.title) > 0) {
	 if(Pint) fprintf(fp1, "Title:%s 1-D aggregated melt\n", silminInputData.title);
	 else fprintf(fp1, "Title:%s Aggregated melt\n", silminInputData.title);
       }
       else {
	 if(Pint) fprintf(fp1, "Title: 1-D aggregated melt\n");
	 else fprintf(fp1, "Title: Aggregated melt\n");
       }

       for (i=1, j=0; j<nc; j++) {
	 if (outsw[j])
	   fprintf(fp1, "Initial Composition: %s %.6f\n", bulkSystem[j].label, compbase_1[i++]);
       }
       if (doingTraces) {
	 for (j=0; j<ntrace; j++) {
	   if(!strcmp(outlb[j], "H2O") && (getenv("ALPHAMELTS_DO_TRACE_H2O") != NULL)) {
	     fprintf(fp1, "Initial Trace: H2O 0.0\n");
	     i++;
	   }
	   else {
	     fprintf(fp1, "Initial Trace: %s %.6f\n", outlb[j], compbase_1[i++]);
	   }
	 }
       }
       fclose(fp1);
     }

     if(Pint) {
       free_dvector(iFdP,1,n+2); free_dvector(iPdF,1,n+2);
       free_dmatrix(icompFdP, 1, nox, 1,n); free_dvector(iPFdP,1,n);
       free_dvector(compbase_2, 1, nox);
     }
     free_dmatrix(icompdF, 1, nox, 1,n+2); free_dvector(compbase_1, 1, nox);
     free_dvector(P,1,n+1); free_dvector(T,1,n+1); free_dmatrix(comp, 1, nox+1, 1,n+1);
     free_dvector(dF,1,n+1); free_dvector(F,1,n+1);
       
   }

   return 0;
}

void CQint(double x[], double y[], int n, double iydx[], int sign) {
   int i;
   double *y2;

   y2 = dvector(1,n);
   spline(x, y, n, 1e32, 1e32, y2);
   for (i=2;i<=n;i++) 
    iydx[i] = iydx[i-1] + sign*qgaus(x[i-1],x[i],x,y,y2,n);
   free_dvector(y2,1,n);
}

void Dumbint(double x[], double y[], int n, double iydx[], int sign) {
   int i;

   /*iydx[1] = 0.0;*/
   for (i=2;i<=n;i++)
     iydx[i] = iydx[i-1] + sign*(x[i]-x[i-1])*y[i];
}

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]) {
   int i, k;
   double p, qn, sig, un, *u;

   u = dvector(1,n-1);
   if (yp1 > 0.99e30)
     y2[1]=u[1]=0.0;
   else {
     y2[1] = -0.5;
     u[1] = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
   }
   for (i=2;i<=n-1;i++) {
     sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
     p=sig*y2[i-1]+2.0;
     y2[i]=(sig-1.0)/p;
     u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
     u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
   }
   if (ypn > 0.99e30)
     qn=un=0.0;
   else {
     qn=0.5;
     un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
   }
   y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
   for (k=n-1;k>=1;k--)
     y2[k]=y2[k]*y2[k+1]+u[k];
   free_dvector(u,1,n-1);
}

void splint(double xa[], double ya[], double y2a[], int n, double x, double *y) {
   static int klo=1, khi=1000;
   int k;
   double h, b, a;

   /* Allow increasing xa or decreasing xa */
   if (khi-klo != 1 || (xa[khi] < x && xa[klo] < x) || (xa[khi] > x && xa[klo] > x)) {
     klo=1;
     khi=n;
     while(khi-klo>1) {
       k=(khi+klo) >> 1;
       if ((xa[k] > x && xa[klo] < x) || (xa[k] < x && xa[klo] > x)) khi=k;
       else klo=k;
     }
   }
   h=xa[khi]-xa[klo];
   if (h==0.0) {
     printf("Splint failure...xa values must be distinct");
     return;
   }
   a = (xa[khi]-x)/h;
   b = (x-xa[klo])/h;
   *y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

double qgaus(double a, double b, double xa[], double ya[], double y2a[], int n)
{
         int j;
         double xr,xm,dx,s,y1,y2;
         static double x[]={0.0,0.1488743389,0.4333953941,
                 0.6794095682,0.8650633666,0.97390652};
         static double w[]={0.0,0.2955242247,0.2692667193,
                 0.2190863625,0.1494513491,0.06667134};

         xm=0.5*(b+a);
         xr=0.5*(b-a);
         s=0;
         for (j=1;j<=5;j++) {
                 dx=xr*x[j];
    splint(xa,ya,y2a,n,xm+dx,&y1);
    splint(xa,ya,y2a,n,xm-dx,&y2);
                 s += w[j]*(y1+y2);
         }
         return s *= xr;
}
