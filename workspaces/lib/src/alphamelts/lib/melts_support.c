const char *melts_support_ver(void) { return "$Id: melts_support.c,v 1.15 2009/04/16 16:35:23 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: melts_support.c,v $
MELTS Source Code: RCS Revision 1.13  2007/11/01 15:44:06  ghiorso
MELTS Source Code: RCS Permit volume regression for solid phases.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.12  2007/10/03 21:33:49  ghiorso
MELTS Source Code: RCS Updated liquid eos thermodynamics.
MELTS Source Code: RCS Added regression of ferric/ferrous parameters from data file.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.11  2007/09/13 16:12:02  ghiorso
MELTS Source Code: RCS (1) Revised standard state liquid properties.
MELTS Source Code: RCS (2) Revised standard state solid properties (removed non-Berman) Cp, and
MELTS Source Code: RCS     removed Saxena EOS treatment.  All EOS parameterizations are Vinet.
MELTS Source Code: RCS     Updated K, K', alpha to conform to Knittle (1995) and Fei (1995)
MELTS Source Code: RCS     except where refitted Berman (1988) makes more sense.
MELTS Source Code: RCS (3) Updated code to allow for fusion entropies of liquid components to
MELTS Source Code: RCS     be adjusted (fusion enthalpies are dependent).
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.10  2007/09/05 17:52:14  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.9  2007/09/03 16:13:20  ghiorso
MELTS Source Code: RCS Corrected error in reading H,S corrections for liquid components in
MELTS Source Code: RCS param_struct_data.h file.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.8  2007/09/02 23:47:35  ghiorso
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.7  2007/08/26 21:38:30  ghiorso
MELTS Source Code: RCS Normalized residuals (for xMELTS calibration) to the number of atoms in the
MELTS Source Code: RCS endmember mineral formula.  Revised residual-statistics.out file accordingly.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.6  2007/06/16 01:01:55  ghiorso
MELTS Source Code: RCS Revised EOS regression to have K', K'', and K''' as parameters.
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.5  2007/06/08 17:25:43  ghiorso
MELTS Source Code: RCS Added code to allow regression of Ghiorso EOS parameters
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.4  2007/05/07 18:23:21  ghiorso
MELTS Source Code: RCS Modifications to LEPR and calibration algorithms following visit by
MELTS Source Code: RCS Marc and Tim.  Mostly work on cpx and opx inclusion and reclassification.
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
MELTS Source Code: RCS Revision 1.4  2004/12/11 22:19:43  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.3  2004/11/21 16:25:32  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.2  2004/10/02 01:33:04  cvsaccount
MELTS Source Code: RCS *** empty log message ***
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2004/01/02 19:21:49  cvsaccount
MELTS Source Code: RCS CTserver University of Chicago
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 5.1  2000/02/15 17:58:12  ghiorso
MELTS Source Code: RCS MELTS 5.0 - xMELTS (associated solutions, multiple liquids)
MELTS Source Code: RCS
 * Revision 3.6  1997/06/21  22:49:44  ghiorso
 * June 1997 MELTS 3.0.x release
 * (prior to new entropy and regression model being introduced)
 *
 * Revision 3.5  1997/05/03  20:23:23  ghiorso
 * *** empty log message ***
 *
 * Revision 3.4  1997/03/27  17:03:27  ghiorso
 * *** empty log message ***
 *
 * Revision 3.3  1996/09/24  20:33:34  ghiorso
 * Version modified for OSF/1 4.0
 *
 * Revision 3.2  1995/12/09  19:26:38  ghiorso
 * Interface revisions for status box and graphics display
 *
 * Revision 3.1  1995/08/18  17:49:28  ghiorso
 * MELTS Version 3 - Initial Entry
 *
 * Revision 3.1  1995/08/18  17:49:28  ghiorso
 * MELTS Version 3 - Initial Entry
 *
*/

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Support routines for converting and initializing data structures
**      and performing onece only conversions  (file: MELTS_SUPPORT.C)
**
**      Functions:
**
**
**  MODIFICATION HISTORY:
**
**      V1.0-1  Mark S. Ghiorso  May 1, 1991
**      V1.0-2  Mark S. Ghiorso  September 6, 1991
**              Removed initialization of defunct variable liquid[].conv
**      V1.0-3  Mark S. Ghiorso  November 2, 1991
**              (1) Altered formulaToMwStoich to allow for decimal points in
**                  numeric entries
**              (2) Made appropriate changes in declaration and matrix
**                  declarations to accomodate this change
**      V1.0-4  Mark S. Ghiorso  April 20, 1992
**              Corrected error in initializing (solids[i].solToLiq)[] array.
**              Changed reference to deallocated invOOt -> invLLt
**      V1.0-5  Mark S. Ghiorso  May 1, 1992
**              Changed liqToOx and solToOx vectors to produce moles of O2
**--
*/

#include "melts_gsl.h"
#include "silmin.h"

typedef struct _elementsData {
   const char *name;
   const char *sym;
   double     atomicWt;
} ElementsData;

const ElementsData elements[] = {
      { "dummy"           , "!" ,   0.0     },
      { "hydrogen"        , "H" ,   1.0079  },
      { "helium"          , "He",   4.00260 },
      { "lithium"         , "Li",   6.94    },
      { "beryllium"       , "Be",   9.01218 },
      { "boron"           , "B" ,  10.81    },
      { "carbon"          , "C" ,  12.011   },
      { "nitrogen"        , "N" ,  14.0067  },
      { "oxygen"          , "O" ,  15.9994  },
      { "fluorine"        , "F" ,  18.998403},
      { "neon"            , "Ne",  20.179   },
      { "sodium"          , "Na",  22.98977 },
      { "magnesium"       , "Mg",  24.305   },
      { "aluminum"        , "Al",  26.98154 },
      { "silicon"         , "Si",  28.0855  },
      { "phosphorous"     , "P" ,  30.97376 },
      { "sulfur"          , "S" ,  32.06    },
      { "chlorine"        , "Cl",  35.453   },
      { "argon"           , "Ar",  39.948   },
      { "potassium"       , "K" ,  39.102   },
      { "calcium"         , "Ca",  40.08    },
      { "scandium"        , "Sc",  44.9559  },
      { "titanium"        , "Ti",  47.90    },
      { "vanadium"        , "V" ,  50.9415  },
      { "chromium"        , "Cr",  51.996   },
      { "manganese"       , "Mn",  54.9380  },
      { "iron"            , "Fe",  55.847   },
      { "cobalt"          , "Co",  58.9332  },
      { "nickel"          , "Ni",  58.71    },
      { "copper"          , "Cu",  63.546   },
      { "zinc"            , "Zn",  65.38    },
      { "gallium"         , "Ga",  69.735   },
      { "germanium"       , "Ge",  72.59    },
      { "arsenic"         , "As",  74.9216  },
      { "selenium"        , "Se",  78.96    },
      { "bromine"         , "Br",  79.904   },
      { "krypton"         , "Kr",  83.80    },
      { "rubidium"        , "Rb",  85.4678  },
      { "strontium"       , "Sr",  87.62    },
      { "yttrium"         , "Y" ,  88.9059  },
      { "zirconium"       , "Zr",  91.22    },
      { "niobium"         , "Nb",  92.9064  },
      { "molybdenum"      , "Mo",  95.94    },
      { "technetium"      , "Tc",  98.9062  },
      { "ruthenium"       , "Ru", 101.07    },
      { "rhodium"         , "Rh", 102.9055  },
      { "palladium"       , "Pd", 106.4     },
      { "silver"          , "Ag", 107.868   },
      { "cadmium"         , "Cd", 112.41    },
      { "indium"          , "In", 114.82    },
      { "tin"             , "Sn", 118.69    },
      { "antimony"        , "Sb", 121.75    },
      { "tellurium"       , "Te", 127.60    },
      { "iodine"          , "I" , 126.9045  },
      { "xenon"           , "Xe", 131.30    },
      { "cesium"          , "Cs", 132.9054  },
      { "barium"          , "Ba", 137.33    },
      { "lantahnum"       , "La", 138.9055  },
      { "cerium"          , "Ce", 140.12    },
      { "praseodymium"    , "Pr", 140.9077  },
      { "neodymium"       , "Nd", 144.24    },
      { "promethium"      , "Pm", 145.      },
      { "samarium"        , "Sm", 150.4     },
      { "europium"        , "Eu", 151.96    },
      { "gadolinium"      , "Gd", 157.25    },
      { "terbium"         , "Tb", 158.9254  },
      { "dysprosium"      , "Dy", 162.50    },
      { "holmium"         , "Ho", 164.9304  },
      { "erbium"          , "Er", 167.26    },
      { "thulium"         , "Tm", 168.9342  },
      { "ytterbium"       , "Yb", 173.04    },
      { "lutetium"        , "Lu", 174.967   },
      { "hafnium"         , "Hf", 178.49    },
      { "tantalum"        , "Ta", 180.9479  },
      { "tungsten"        , "W" , 183.85    },
      { "rhenium"         , "Re", 186.207   },
      { "osmium"          , "Os", 190.2     },
      { "iridium"         , "Ir", 192.22    },
      { "platinum"        , "Pt", 195.09    },
      { "gold"            , "Au", 196.9665  },
      { "mercury"         , "Hg", 200.59    },
      { "thallium"        , "Tl", 204.37    },
      { "lead"            , "Pb", 207.2     },
      { "bismuth"         , "Bi", 208.9804  },
      { "polonium"        , "Po", 209.      },
      { "astatine"        , "At", 210.      },
      { "radon"           , "Rn", 222.      },
      { "francium"        , "Fr", 223.      },
      { "radium"          , "Ra", 226.0254  },
      { "actinium"        , "Ac", 227.      },
      { "thorium"         , "Th", 232.0381  },
      { "protactinium"    , "Pa", 231.0359  },
      { "uranium"         , "U" , 238.029   },
      { "neptunium"       , "Np", 237.0482  },
      { "plutonium"       , "Pu", 244.      },
      { "americium"       , "Am", 243.      },
      { "curium"          , "Cm", 247.      },
      { "berkelium"       , "Bk", 247.      },
      { "californium"     , "Cf", 251.      },
      { "einsteinium"     , "Es", 254.      },
      { "fermium"         , "Fm", 257.      },
      { "mendelevium"     , "Md", 258.      },
      { "nobelium"        , "No", 259.      },
      { "lawrencium"      , "Lw", 260.      },
      { "ruferfordium"    , "Rf", 260.      },
      { "hahnium"         , "Ha", 260.      },
      { "106ium"          , "ZZ", 263.      }
   };
const int ne = (sizeof elements / sizeof (struct _elementsData));

static double **matrix_alloc(int n1, int n2) {
   int i;
   double *m0 = (double *) malloc((size_t) n1*n2*sizeof(double));
   double **m = (double **) malloc((size_t) n1*sizeof(double *));
   for (i=0; i<n1; i++) m[i] = &m0[i*n2];
   return m;
}

static void matrix_free(double **m, int n1, int n2) {
   free(m[0]);
   free(m);
}

void InitComputeDataStruct(void)
{
   int i, j, k, atNumO;
   gsl_vector *S;
   gsl_matrix *invOOt, *OLt, *OSt, *invLLt, *LOt, *LSt, *V;
   gsl_vector_view b, x;
   double **elementsToOxides, **elementsToLiquid, **elementsToSolids;

   /*========================================================================
     Initialize Phase; Compute element contents/molecular weights
     ========================================================================*/

   /* Allocate temporary storage for element abundance matrices              */

   elementsToOxides = matrix_alloc((size_t)  nc, (size_t) ne);  /* matrix: O */
   elementsToLiquid = matrix_alloc((size_t) nls, (size_t) ne);  /* matrix: L */
   elementsToSolids = matrix_alloc((size_t) npc, (size_t) ne);  /* matrix: S */

   /* Form the matrices O, L and S and compute molecular weights             */

   for(i=0; i<nc; i++) {
      for(j=0; j<ne; j++) elementsToOxides[i][j] = 0.0;
      bulkSystem[i].mw =
         formulaToMwStoich((char *) bulkSystem[i].label, elementsToOxides[i]);
      bulkSystem[i].oxToElm = (int *) malloc((size_t) ne*sizeof(int));
      for (j=0; j<ne; j++)
        (bulkSystem[i].oxToElm)[j] = (int) elementsToOxides[i][j];
   }
   if (calculationMode == MODE_xMELTS) {
     for(i=0; i<nls; i++) {
       size_t lenToken   = strcspn(liquid[i].label, " ");
       size_t lenBracket = strcspn(liquid[i].label, "[");
       size_t lenFormula = (lenToken > lenBracket) ? lenBracket : lenToken;
       char *trueFormula = (char *) malloc((size_t) (lenFormula+1)*sizeof(char));
       (void) strncpy(trueFormula, liquid[i].label, lenFormula); trueFormula[lenFormula] = '\0';
       for(j=0; j<ne; j++) elementsToLiquid[i][j] = 0.0;
       (void) formulaToMwStoich(trueFormula, elementsToLiquid[i]);
       free(trueFormula);
     }
   } else { /* MODE__MELTS or MODE_pMELTS */
     for(i=0; i<nls; i++) {
       for(j=0; j<ne; j++) elementsToLiquid[i][j] = 0.0;
       (void) formulaToMwStoich((char *) liquid[i].label, elementsToLiquid[i]);
     }
   }
   for(i=0; i<npc; i++) {
      for(j=0; j<ne; j++) elementsToSolids[i][j] = 0.0;
      solids[i].mw = formulaToMwStoich((char *) solids[i].formula, elementsToSolids[i]);
      for(j=0,solids[i].nAtoms=0.0; j<ne; j++) solids[i].nAtoms += elementsToSolids[i][j];
      if (solids[i].type == PHASE) {
         for (j=i+1; solids[j].type == COMPONENT && j<npc; j++);
         j--;
         solids[i].na = MAX(j-i  , 1);
         solids[i].nr = MAX(j-i-1, 0);
      }
   }

   /*========================================================================
     Compute entries for the and oxygen structure
     ========================================================================*/

   for (atNumO=0; atNumO<ne; atNumO++) {
      if (strcmp(elements[atNumO].name,"oxygen") == 0) break;
   }
   if (atNumO == ne) exit(1);

   oxygen.liqToOx = (double *) malloc((size_t) nlc*sizeof(double));
   for (i=0; i<nlc; i++) {
      (oxygen.liqToOx)[i] = elementsToLiquid[i][atNumO]/2.0;
   }

   oxygen.solToOx = (double *) malloc((size_t) npc*sizeof(double));
   for (i=0; i<npc; i++) {
      (oxygen.solToOx)[i] = elementsToSolids[i][atNumO]/2.0;
   }

   /*========================================================================
     Compute liqToOx and solToOx entries in liquid[] and solids[] structures
     ========================================================================*/

   /* Compute the square (nc-1 by nc-1) matrix  (O O(t))(-1)                 */

   S = gsl_vector_alloc((size_t) nc);
   V = gsl_matrix_alloc((size_t) nc, (size_t) nc);
   invOOt = gsl_matrix_alloc((size_t) nc, (size_t) nc);
   for (i=0; i<nc; i++) {
      for (j=0; j<nc; j++) {
         double sum = 0.0;
         for (k=0; k<ne; k++) sum += elementsToOxides[i][k]*elementsToOxides[j][k];
         gsl_matrix_set(invOOt, i, j, sum);
      }
   }
   gsl_linalg_SV_decomp_jacobi(invOOt, V, S);

   /* Compute the (nc-1 by nls-1) matrix O L(t)                              */

   OLt = gsl_matrix_alloc((size_t) nc, (size_t) nls);
   for (i=0; i<nc; i++) {
      for (j=0; j<nls; j++) {
         double sum = 0.0;
         for (k=0; k<ne; k++) sum += elementsToOxides[i][k]*elementsToLiquid[j][k];
         gsl_matrix_set(OLt, i, j, sum);
      }
   }

   /* Compute (O O(t))(-1) . (O L(t)) to create the liqToOx vectors required
      for initializing the array of Liquid structures                        */

   for (i=0; i<nls; i++) {
      liquid[i].liqToOx = (double *) malloc((size_t) ((calculationMode == MODE_xMELTS) ? nc+1 : nc)*sizeof(double));
      b = gsl_matrix_column(OLt, i);
      x = gsl_vector_view_array(liquid[i].liqToOx, (size_t) nc);

      gsl_linalg_SV_solve(invOOt, V, S, &b.vector, &x.vector);

      for (j=0; j<nc; j++) {
         if (ABS((liquid[i].liqToOx)[j]) < sqrt(TAU)) (liquid[i].liqToOx)[j] = 0.0;
      }
      if (calculationMode == MODE_xMELTS) (liquid[i].liqToOx)[nc] = 0.0;
   }

   /* Free the space held by the matrix O L(t)                               */

   gsl_matrix_free(OLt);

   /* Compute the rectangular (nc-1 by npc-1) matrix O S(t)                  */

   OSt = gsl_matrix_alloc((size_t) nc, (size_t) npc);
   for (i=0; i<nc; i++) {
      for (j=0; j<npc; j++) {
         double sum = 0.0;
         for (k=0; k<ne; k++) sum += elementsToOxides[i][k]*elementsToSolids[j][k];
         gsl_matrix_set(OSt, i, j, sum);

      }
   }

   /* Compute (O O(t))(-1) . (O S(t)) to create the solToOx vectors required
      for initializing the array of Solid structures                         */

   for (i=0; i<npc; i++) {
      solids[i].solToOx = (double *) malloc((size_t) nc*sizeof(double));
      b = gsl_matrix_column(OSt, i);
      x = gsl_vector_view_array(solids[i].solToOx, (size_t) nc);

      gsl_linalg_SV_solve(invOOt, V, S, &b.vector, &x.vector);

      for (j=0; j<nc; j++) {
         if (ABS((solids[i].solToOx)[j]) < sqrt(TAU)) (solids[i].solToOx)[j] = 0.0;
      }
   }

   /* Free the space held by the matrix O S(t)                               */

   gsl_matrix_free(OSt);

   /* Free the space held by the matrix (O O(t))(-1)                         */

   gsl_matrix_free(invOOt);

   gsl_vector_free(S); gsl_matrix_free(V);

   /*========================================================================
     Compute oxToLiq entries in bulkSystem[] structure
     ========================================================================*/

   /* Compute the square (nlc-1 by nlc-1) matrix  (L L(t))(-1)               */
   S = gsl_vector_alloc((size_t) nlc);
   V = gsl_matrix_alloc((size_t) nlc, (size_t) nlc);
   invLLt = gsl_matrix_alloc((size_t) nlc, (size_t) nlc);
   for (i=0; i<nlc; i++) {
      for (j=0; j<nlc; j++) {
         double sum = 0.0;
         for (k=0; k<ne; k++) sum += elementsToLiquid[i][k]*elementsToLiquid[j][k];
         gsl_matrix_set(invLLt, i, j, sum);
      }
   }
   gsl_linalg_SV_decomp_jacobi(invLLt, V, S);

   /* Compute the rectangular (nlc-1 by nc-1) matrix L O(t)                  */

   LOt = gsl_matrix_alloc((size_t) nlc, (size_t) nc);
   for (i=0; i<nlc; i++) {
      for (j=0; j<nc; j++) {
         double sum = 0.0;
         for (k=0; k<ne; k++) sum += elementsToLiquid[i][k]*elementsToOxides[j][k];
         gsl_matrix_set(LOt, i, j, sum);
      }
   }

   /* Compute (L L(t))(-1) . (L O(t)) to create the oxToLiq vectors required
      for initializing the array of BulkSystem structures                    */

   for (i=0; i<nc; i++) {
      bulkSystem[i].oxToLiq = (double *) malloc((size_t) nlc*sizeof(double));
      b = gsl_matrix_column(LOt, i);
      x = gsl_vector_view_array(bulkSystem[i].oxToLiq, (size_t) nlc);

      gsl_linalg_SV_solve(invLLt, V, S, &b.vector, &x.vector);

      for (j=0; j<nlc; j++) {
         if (ABS((bulkSystem[i].oxToLiq)[j]) < sqrt(TAU)) (bulkSystem[i].oxToLiq)[j] = 0.0;
      }
   }

   /* Free the space held by the matrix O L(t)                               */

   gsl_matrix_free(LOt);

   /* Compute the rectangular (nlc-1 by npc-1) matrix L S(t)                 */

   LSt = gsl_matrix_alloc((size_t) nlc, (size_t) npc);
   for (i=0; i<nlc; i++) {
      for (j=0; j<npc; j++) {
         double sum = 0.0;
         for (k=0; k<ne; k++) sum += elementsToLiquid[i][k]*elementsToSolids[j][k];
         gsl_matrix_set(LSt, i, j, sum);
      }
   }

   /* Compute (L L(t))(-1) . (L S(t)) to create the solToLiq vectors required
      for initializing the array of Solid structures                         */

   for (i=0; i<npc; i++) {
      solids[i].solToLiq = (double *) malloc((size_t) nlc*sizeof(double));
      b = gsl_matrix_column(LSt, i);
      x = gsl_vector_view_array(solids[i].solToLiq, (size_t) nlc);

      gsl_linalg_SV_solve(invLLt, V, S, &b.vector, &x.vector);

      for (j=0; j<nlc; j++) {
         if (ABS((solids[i].solToLiq)[j]) < sqrt(TAU)) (solids[i].solToLiq)[j] = 0.0;
      }
   }

   /* Free the space held by the matrix L S(t)                               */

   gsl_matrix_free(LSt);

   /* Free the space held by the matrix (O O(t))(-1)                         */

   gsl_matrix_free(invLLt);

   gsl_matrix_free(V); gsl_vector_free(S);

   /* If mode is xMELTS then pre-process liquid structure entries            */

   if (calculationMode == MODE_xMELTS) {
     int indexSiO2      = -1; /* in liquid components */
     int indexAl2O3     = -1; /* in oxide array       */
     int indexFe2SiO4_6 = -1; /* in liquid species    */
     int indexFeO1_3 = nc; /* convention */
     for (i=0;  i<nls; i++) if (!strcmp(liquid[i].label, "SiO2"))      { indexSiO2      = i; break; }
#ifndef BUILD_MGO_SIO2_VERSION
     for (i=0;  i<nc;  i++) if (!strcmp(bulkSystem[i].label, "Al2O3")) { indexAl2O3     = i; break; }
#ifndef BUILD_SIO2_AL2O3_CAO_NA2O_K2O_VERSION
     for (i=nc; i<nls; i++) if (!strcmp(liquid[i].label, "Fe2SiO4.6")) { indexFe2SiO4_6 = i; break; }
#ifdef TRUE_xMELTS
     if (indexFe2SiO4_6 == -1) printf("ERROR in InitComputeDataStruct. Cannot find species Fe2SiO4.6\n"); exit(0);
#endif
#endif
#ifdef TRUE_xMELTS
     if (indexAl2O3 == -1)     printf("ERROR in InitComputeDataStruct. Cannot find component Al2O3\n"); exit(0);
#endif
#endif

     for (i=0; i<nc; i++) {
       double v0    = bulkSystem[i].gk_v;
       double alpha = (v0 != 0.0) ? bulkSystem[i].gk_dvdt/v0 : 0.0;
       double c     = bulkSystem[i].gk_c;
       double cp    = bulkSystem[i].gk_cp;
       double mw    = bulkSystem[i].mw;
       double v1    = ((mw != 0.0) && (c != 0.0) && (cp != 0.0)) ? -v0*v0*(1000.0/(mw*c*c) + 1673.15*alpha*alpha/(cp)) : 0.0;
       double K     = (v1 != 0.0) ? -v0/v1 : 0.0;
       double Kp    = eosModelParameters[i].Kp;
       double Kpp   = eosModelParameters[i].Kpp;
       double Kppp  = eosModelParameters[i].Kppp;

       eosModelParameters[i].v2 = (K != 0.0) ? (Kp+1.0)*v0/K/K : 0.0;
       eosModelParameters[i].v3 = (K != 0.0) ? (Kpp*K - (2.0*Kp+1.0)*(Kp+1.0))*v0/K/K/K : 0.0;
       eosModelParameters[i].v4 = (K != 0.0) ? (Kppp*K*K - Kpp*K*(4.0+6.0*Kp) + (3.0*Kp+1.0)*(2.0*Kp+1.0)*(Kp+1.0))*v0/K/K/K/K : 0.0;

       bulkSystem[i].gk_d2vdp2 = 0.0; /* entries in liq_struct_data.h are ignored */
       bulkSystem[i].gk_d3vdp3 = 0.0; /* entries in liq_struct_data.h are ignored */
       bulkSystem[i].gk_d4vdp4 = 0.0; /* entries in liq_struct_data.h are ignored */
     }

#ifndef BUILD_MGO_SIO2_VERSION
#ifndef BUILD_SIO2_AL2O3_CAO_NA2O_K2O_VERSION
     { /* Convert the intermediate iron species */
       double v0    = bulkSystem[indexFeO1_3].gk_v;
       double alpha = (v0 != 0.0) ? bulkSystem[indexFeO1_3].gk_dvdt/v0 : 0.0;
       double c     = bulkSystem[indexFeO1_3].gk_c;
       double cp    = bulkSystem[indexFeO1_3].gk_cp;
       double mw    = bulkSystem[indexFeO1_3].mw;
       double v1    = ((mw != 0.0) && (c != 0.0) && (cp != 0.0)) ? -v0*v0*(1000.0/(mw*c*c) + 1673.15*alpha*alpha/(cp)) : 0.0;
       double K     = (v1 != 0.0) ? -v0/v1 : 0.0;
       double Kp    = eosModelParameters[indexFeO1_3].Kp;
       double Kpp   = eosModelParameters[indexFeO1_3].Kpp;
       double Kppp  = eosModelParameters[indexFeO1_3].Kppp;

       eosModelParameters[indexFeO1_3].v2 = (K != 0.0) ? (Kp+1.0)*v0/K/K : 0.0;
       eosModelParameters[indexFeO1_3].v3 = (K != 0.0) ? (Kpp*K - (2.0*Kp+1.0)*(Kp+1.0))*v0/K/K/K : 0.0;
       eosModelParameters[indexFeO1_3].v4 = (K != 0.0) ? (Kppp*K*K - Kpp*K*(4.0+6.0*Kp) + (3.0*Kp+1.0)*(2.0*Kp+1.0)*(Kp+1.0))*v0/K/K/K/K : 0.0;

       bulkSystem[indexFeO1_3].gk_d2vdp2 = 0.0; /* entries in liq_struct_data.h are ignored */
       bulkSystem[indexFeO1_3].gk_d3vdp3 = 0.0; /* entries in liq_struct_data.h are ignored */
       bulkSystem[indexFeO1_3].gk_d4vdp4 = 0.0; /* entries in liq_struct_data.h are ignored */
     }
#endif
#endif

     /* Compute volumetric properties of liquid species                       */
     /* --> Note loop is executed only if the liquid volume is initially zero */

     for (i=0; i<nls; i++) if ((liquid[i].liq.v == 0.0) && (liquid[i].liq.eos_type == EOS_GHIORSO) ) {
       double totalMoles = 0.0;

       for (j=0; j<nc; j++) {
         totalMoles                       += (liquid[i].liqToOx)[j];
         liquid[i].liq.v                  += (liquid[i].liqToOx)[j]*bulkSystem[j].gk_v;
         liquid[i].liq.eos.Ghiorso.alpha  += (liquid[i].liqToOx)[j]*bulkSystem[j].gk_dvdt;
         liquid[i].liq.eos.Ghiorso.c      += (liquid[i].liqToOx)[j]*bulkSystem[j].gk_c;
         liquid[i].liq.eos.Ghiorso.dcdt   += (liquid[i].liqToOx)[j]*bulkSystem[j].gk_dcdt;
         liquid[i].liq.eos.Ghiorso.cp     += (liquid[i].liqToOx)[j]*bulkSystem[j].gk_cp;
         liquid[i].liq.eos.Ghiorso.mw     += (liquid[i].liqToOx)[j]*bulkSystem[j].mw;
         liquid[i].liq.eos.Ghiorso.d2vdp2 += (liquid[i].liqToOx)[j]*(bulkSystem[j].gk_d2vdp2 + eosModelParameters[j].v2);
         liquid[i].liq.eos.Ghiorso.d3vdp3 += (liquid[i].liqToOx)[j]*(bulkSystem[j].gk_d3vdp3 + eosModelParameters[j].v3);
         liquid[i].liq.eos.Ghiorso.d4vdp4 += (liquid[i].liqToOx)[j]*(bulkSystem[j].gk_d4vdp4 + eosModelParameters[j].v4);
       }

       if (totalMoles == 0.0) totalMoles = 1.0;

       liquid[i].liq.eos.Ghiorso.alpha /= (liquid[i].liq.v != 0.0) ? liquid[i].liq.v : 1.0;
       liquid[i].liq.eos.Ghiorso.c     /= totalMoles;
       liquid[i].liq.eos.Ghiorso.dcdt  /= totalMoles;

       if (indexAl2O3 != -1) for (j=0; j<nc; j++) liquid[i].liq.eos.Ghiorso.c +=
	 (liquid[i].liqToOx)[j]*(liquid[i].liqToOx)[indexAl2O3]*bulkSystem[j].gk_cXal2o3/(totalMoles*totalMoles);

       liquid[i].liq.cp = liquid[i].liq.eos.Ghiorso.cp;

     } else if (liquid[i].liq.cp == 0.0) {
       for (j=0; j<nc; j++) liquid[i].liq.cp += (liquid[i].liqToOx)[j]*bulkSystem[j].gk_cp;

       /* special case of species with intermediate oxidation state    */
       /*   One bar volumetric properties are set in liq_struct_data.h */
       /*   Higher order properties are reinitialized here             */
       if (i == indexFe2SiO4_6) {
         liquid[i].liq.eos.Ghiorso.d2vdp2 =     (bulkSystem[indexSiO2  ].gk_d2vdp2 + eosModelParameters[indexSiO2  ].v2)
	                                  + 2.0*(bulkSystem[indexFeO1_3].gk_d2vdp2 + eosModelParameters[indexFeO1_3].v2);
         liquid[i].liq.eos.Ghiorso.d3vdp3 =     (bulkSystem[indexSiO2  ].gk_d3vdp3 + eosModelParameters[indexSiO2  ].v3)
	                                  + 2.0*(bulkSystem[indexFeO1_3].gk_d3vdp3 + eosModelParameters[indexFeO1_3].v3);
         liquid[i].liq.eos.Ghiorso.d4vdp4 =     (bulkSystem[indexSiO2  ].gk_d4vdp4 + eosModelParameters[indexSiO2  ].v4)
	                                  + 2.0*(bulkSystem[indexFeO1_3].gk_d4vdp4 + eosModelParameters[indexFeO1_3].v4);
       }

     } /* end if block and loop on i over all species */

     /* Convert thermodynamic data for liquid species to a common set of      */
     /* liquid properties at a reference T of 1673.15 K and pressure of 1 bar */

     for (i=0; i<nls; i++) if (liquid[i].liq.sfus != 0.0) {
       double t = 1673.15;
       double p = 1.0;
       gibbs(t, p, (char *) liquid[i].label, &(liquid[i].ref), &(liquid[i].liq), &(liquid[i].fus), &(liquid[i].cur));

       /* Reset the ThermoRef structure         */
       liquid[i].ref.h = liquid[i].cur.h;
       liquid[i].ref.s = liquid[i].cur.s;
       if      (liquid[i].ref.cp_type == CP_BERMAN) {
         liquid[i].ref.cp.Berman.k0     = 0.0;
         liquid[i].ref.cp.Berman.k1     = 0.0;
         liquid[i].ref.cp.Berman.k2     = 0.0;
         liquid[i].ref.cp.Berman.k3     = 0.0;
         liquid[i].ref.cp.Berman.Tt     = 0.0;
         liquid[i].ref.cp.Berman.deltah = 0.0;
         liquid[i].ref.cp.Berman.l1     = 0.0;
         liquid[i].ref.cp.Berman.l2     = 0.0;
       } else if (liquid[i].ref.cp_type == CP_SAXENA) {
         liquid[i].ref.cp.Saxena.a = 0.0;
         liquid[i].ref.cp.Saxena.b = 0.0;
         liquid[i].ref.cp.Saxena.c = 0.0;
         liquid[i].ref.cp.Saxena.d = 0.0;
         liquid[i].ref.cp.Saxena.e = 0.0;
         liquid[i].ref.cp.Saxena.g = 0.0;
         liquid[i].ref.cp.Saxena.h = 0.0;
       }

       /* Reset the ThermoLiq structure         */
       liquid[i].tfus     = liquid[i].liq.tfus;
       liquid[i].liq.tfus = t;
       liquid[i].liq.sfus = 0.0;

       /* Rezero the ThermoData *fus* structure */
       liquid[i].fus.g       = 0.0;
       liquid[i].fus.h       = 0.0;
       liquid[i].fus.s       = 0.0;
       liquid[i].fus.v       = 0.0;
       liquid[i].fus.cp      = 0.0;
       liquid[i].fus.dcpdt   = 0.0;
       liquid[i].fus.dvdt    = 0.0;
       liquid[i].fus.dvdp    = 0.0;
       liquid[i].fus.d2vdt2  = 0.0;
       liquid[i].fus.d2vdp2  = 0.0;
       liquid[i].fus.d2vdtdp = 0.0;

     }

     /* Deal with the "dependent" species case in xMELTS, where thermodynamic */
     /* properties are computed from a basis set                              */

     for (i=0; i<nls; i++) {
       char   *next  = NULL, *temp = NULL;
       double  coeff = 0.0;

       (void) strtok(strcpy(temp = (char *) malloc((size_t) (strlen(liquid[i].label)+1)*sizeof(char)), liquid[i].label), " ");

       while ( (next = strtok(NULL, " ")) != NULL) {
         if (!strcmp(next,"=")) continue; /* to next while test */

	 /* Assume coeff,species pairs separated by " ", no space between minus sign and number */
	 if (coeff == 0.0) coeff = atof(next);
	 else {
	   for (j=0; j<nls; j++) if (!strcmp(next,liquid[j].label)) break; /* from loop, j = inddex */
	   if (j == nls) {
	     printf("Fatal error in InitComputeDataStruct[melts_support.c at line %d].\n", __LINE__);
	     printf("  ----->Attempting to balance reaction for species %s\n", liquid[i].label);
	     printf("  ----->with species unrecognizable %s.  Aborting...\n", next);
	     abort();
	   } else { /* found BOTH the coeff and the species identifier, so assign properties */
	     liquid[i].ref.h += coeff*liquid[j].ref.h;
	     liquid[i].ref.s += coeff*liquid[j].ref.s;
	   }
	   coeff = 0.0; /* reset so next coeff,species pair can be found */
	 }

       } /* closes while loop on reaction coefficients/species */

       if(temp != NULL) free(temp);
     } /* end of loop over all species */

#ifdef TRUE_xMELTS
     printf("---> Correct the std state liquid properties and zero the modelParameters? (y or n): ");
     if (tolower(getchar()) == 'y') {
       getchar();
       for (i=0; i<nls; i++) {
         liquid[i].ref.h += modelParameters[nls*(nls-1)/2+i].enthalpy;
         liquid[i].ref.s += modelParameters[nls*(nls-1)/2+i].entropy;
         liquid[i].liq.v += modelParameters[nls*(nls-1)/2+i].volume;

         modelParameters[nls*(nls-1)/2+i].enthalpy = 0.0;
         modelParameters[nls*(nls-1)/2+i].entropy  = 0.0;
         modelParameters[nls*(nls-1)/2+i].volume   = 0.0;
       }
     } else getchar();
#endif

     for (i=0; i<npc; i++) {
       solids[i].ref.h += modelParameters[nls*(nls-1)/2+nls+i].enthalpy;
       solids[i].ref.s += modelParameters[nls*(nls-1)/2+nls+i].entropy;
       solids[i].ref.v += modelParameters[nls*(nls-1)/2+nls+i].volume;

       modelParameters[nls*(nls-1)/2+nls+i].enthalpy = 0.0;
       modelParameters[nls*(nls-1)/2+nls+i].entropy  = 0.0;
       modelParameters[nls*(nls-1)/2+nls+i].volume   = 0.0;
     }

   } /* End of xMELTS specific code */

   /*========================================================================
     Clean up and exit
     ========================================================================*/

   /* Free the space for the temorpary matrices O, nd S.                     */

   matrix_free(elementsToOxides,  nc, ne);  /* matrix: L */
   matrix_free(elementsToLiquid, nls, ne);  /* matrix: L */
   matrix_free(elementsToSolids, npc, ne);  /* matrix: S */

}

double formulaToMwStoich(char *formula, double *stoich)
{
   int i, j, len;
   char c, *sym, *num, *temp;
   double mw = 0.0, mult = 1.0;

   if((len = (int) strlen(formula)) == 0) return 0.0;

   sym  = (char *) malloc((unsigned)       3*sizeof(char));
   temp = (char *) malloc((unsigned) (len+2)*sizeof(char));
   num  = (char *) malloc((unsigned) (len+2)*sizeof(char));
   sym[0] = '\0';
   num[0] = '\0';
   j = 0;

   temp = strcpy(temp, formula);
   temp = strcat(temp," ");

   for (i=0; i<=len; i++) {
      c = temp[i];
      if(isalpha(c) || isspace(c) || c == '(' || c == ')') {
         if( !islower(c) && strlen(sym) > 0) {
            for (j=1; j<ne; j++) {
               if(strcmp(sym, elements[j].sym) == 0) break;
            }
            if (j == ne) {free(temp); free(num); return 0.0;}
            sym[0] = '\0';
         }
         if (j > 0) {
            if (num[0] == '\0') {
               stoich[j] += mult;
               mw += elements[j].atomicWt*mult;
            } else {
               stoich[j] += atof(num)*mult;
               mw += atof(num) * elements[j].atomicWt * mult;
            }
            j = 0;
            num[0] = '\0';
         }
         if (c == '(') {
            int close, finish, k;

            for(close=i+1; temp[close] != ')' && close < len; close++)
                ;
             if(close == len) {free(sym); free(temp); free(num); return 0.0;}
            close++;
            for(finish=close; !isupper(temp[finish]) &&
                (temp[finish] != '(') && (finish < len); finish++);

            if (finish > close) {
               for(k=0; k < (finish-close); k++) {
                  num[k] = temp[close+k];
               }
               num[k] = '\0';
               mult = atof(num);
               for(k=0; k <= (len+1-finish); k++) {
                  temp[close+k] = temp[finish+k];
               }
               len -= strlen(num);
               num[0] = '\0';
            }

         } else if (c == ')') {
            mult = 1.0;
         } else {
            sym = strncat(sym, &temp[i], 1);
         }
      } else if (isdigit(c) || c == '.' || c == '-') { /* decimal amount of element */
         if(num[0] == '\0') {
            for (j=1; j<ne; j++) {
               if(strcmp(sym, elements[j].sym) == 0) break;
            }
            if (j == ne) {free(temp); free(num); return 0.0;}
            sym[0] = '\0';
         }
         num = strncat(num, &temp[i], 1);
      } else {                                  /* Illegal character         */
         free(temp); free(num); return 0.0;
      }
   }

   free(sym);
   free(temp);
   free(num);

   return mw;
}

/* end of file MELTS_SUPPORT.C */
