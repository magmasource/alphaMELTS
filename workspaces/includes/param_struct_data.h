#ifndef _Param_Struct_Data_h
#define _Param_Struct_Data_h

#ifdef USE_NEW_MELTS_WATER_MODEL
#undef USE_NEW_MELTS_WATER_MODEL
#endif

/*
**++
**  FACILITY:  Silicate Melts Regression Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Include file for initializing values of liquid interaction
**         parameters (file: PARAM_STRUCT_DATA.H)
**
**  This file produced on: Tuesday, May 18, 2010
**--
*/
 
/*
 *=============================================================================
 * Initialize global values of interaction parameters
 */

static ModelParameters originalModelParameters[] = {
  { "W(TiO2	     ,SiO2      )",       26266.7,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Al2O3	 ,SiO2      )",      -39120.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Fe2O3	 ,SiO2      )",        8110.3,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(MgCr2O4   ,SiO2      )",       27886.3,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Fe2SiO4   ,SiO2      )",       23660.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(MnSi0.5O2 ,SiO2      )",       18393.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Mg2SiO4   ,SiO2      )",        3421.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(NiSi0.5O2 ,SiO2      )",       25197.4,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CoSi0.5O2 ,SiO2      )",       14802.8,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CaSiO3	 ,SiO2      )",        -863.7,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Na2SiO3   ,SiO2      )",      -99039.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(KAlSiO4   ,SiO2      )",      -33921.7,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Ca3(PO4)2 ,SiO2      )",       61891.6,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CO2	     ,SiO2      )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(S	     ,SiO2      )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Cl	     ,SiO2      )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(F	     ,SiO2      )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#ifdef USE_NEW_MELTS_WATER_MODEL
  { "W(H2O	     ,SiO2      )",       27357.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#else
  { "W(H2O	     ,SiO2      )",       30967.3,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#endif
  { "W(CaCO3     ,SiO2      )",       63281.2,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "W(Al2O3	 ,TiO2      )",      -29449.8,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Fe2O3	 ,TiO2      )",      -84756.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(MgCr2O4   ,TiO2      )",      -72303.4,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Fe2SiO4   ,TiO2      )",        5209.1,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(MnSi0.5O2 ,TiO2      )",      -16123.5,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Mg2SiO4   ,TiO2      )",       -4178.3,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(NiSi0.5O2 ,TiO2      )",        3614.8,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CoSi0.5O2 ,TiO2      )",       -1640.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CaSiO3	 ,TiO2      )",      -35372.5,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Na2SiO3   ,TiO2      )",      -15415.6,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(KAlSiO4   ,TiO2      )",      -48094.6,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Ca3(PO4)2 ,TiO2      )",       25938.8,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CO2	     ,TiO2      )", 	 -19265.7,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(S	     ,TiO2      )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Cl	     ,TiO2      )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(F	     ,TiO2      )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#ifdef USE_NEW_MELTS_WATER_MODEL
  { "W(H2O	     ,TiO2      )",       88199.2,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#else
  { "W(H2O	     ,TiO2      )",       81879.1,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#endif
  { "W(CaCO3     ,TiO2      )",      -79202.7,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "W(Fe2O3	 ,Al2O3     )",      -17089.4,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(MgCr2O4   ,Al2O3     )",      -31770.3,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Fe2SiO4   ,Al2O3     )",      -30509.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(MnSi0.5O2 ,Al2O3     )",      -53874.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Mg2SiO4   ,Al2O3     )",      -32880.3,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(NiSi0.5O2 ,Al2O3     )",        2985.2,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CoSi0.5O2 ,Al2O3     )",       -2677.4,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CaSiO3	 ,Al2O3     )",      -57917.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Na2SiO3   ,Al2O3     )",     -130785.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(KAlSiO4   ,Al2O3     )",      -25859.2,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Ca3(PO4)2 ,Al2O3     )",       52220.8,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CO2	     ,Al2O3     )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(S	     ,Al2O3     )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Cl	     ,Al2O3     )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(F	     ,Al2O3     )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#ifdef USE_NEW_MELTS_WATER_MODEL
  { "W(H2O	     ,Al2O3     )",       11768.4,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#else
  { "W(H2O	     ,Al2O3     )",      -16098.1,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#endif
  { "W(CaCO3     ,Al2O3     )",       46716.1,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  
  { "W(MgCr2O4   ,Fe2O3     )",       21605.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Fe2SiO4   ,Fe2O3     )",     -179064.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(MnSi0.5O2 ,Fe2O3     )",        3907.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Mg2SiO4   ,Fe2O3     )",      -71518.6,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(NiSi0.5O2 ,Fe2O3     )", 	    408.7,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CoSi0.5O2 ,Fe2O3     )",        -223.7,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CaSiO3	 ,Fe2O3     )",       12076.6,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Na2SiO3   ,Fe2O3     )",     -149662.2,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(KAlSiO4   ,Fe2O3     )",       57555.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Ca3(PO4)2 ,Fe2O3     )",       -4213.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CO2	     ,Fe2O3     )", 	  -3186.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(S	     ,Fe2O3     )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Cl	     ,Fe2O3     )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(F	     ,Fe2O3     )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#ifdef USE_NEW_MELTS_WATER_MODEL
  { "W(H2O	     ,Fe2O3     )",       50105.2,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#else
  { "W(H2O	     ,Fe2O3     )",       31405.5,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#endif
  { "W(CaCO3     ,Fe2O3     )",       65508.7,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "W(Fe2SiO4   ,MgCr2O4   )",      -82971.8,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(MnSi0.5O2 ,MgCr2O4   )", 	    182.4,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Mg2SiO4   ,MgCr2O4   )",       46049.2,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(NiSi0.5O2 ,MgCr2O4   )",        -266.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CoSi0.5O2 ,MgCr2O4   )",        -384.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CaSiO3	 ,MgCr2O4   )",       30704.7,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Na2SiO3   ,MgCr2O4   )",      113646.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(KAlSiO4   ,MgCr2O4   )",       75709.1,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Ca3(PO4)2 ,MgCr2O4   )",        5341.8,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CO2	     ,MgCr2O4   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(S	     ,MgCr2O4   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Cl	     ,MgCr2O4   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(F	     ,MgCr2O4   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(H2O	     ,MgCr2O4   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CaCO3     ,MgCr2O4   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "W(MnSi0.5O2 ,Fe2SiO4   )",       -6823.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Mg2SiO4   ,Fe2SiO4   )",      -37256.7,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(NiSi0.5O2 ,Fe2SiO4   )",      -17019.8,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CoSi0.5O2 ,Fe2SiO4   )",      -11746.3,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CaSiO3	 ,Fe2SiO4   )",      -12970.8,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Na2SiO3   ,Fe2SiO4   )",      -90533.8,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(KAlSiO4   ,Fe2SiO4   )",       23649.4,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Ca3(PO4)2 ,Fe2SiO4   )",       87410.3,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CO2	     ,Fe2SiO4   )", 	 -32464.5,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(S	     ,Fe2SiO4   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Cl	     ,Fe2SiO4   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(F	     ,Fe2SiO4   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#ifdef USE_NEW_MELTS_WATER_MODEL
  { "W(H2O	     ,Fe2SiO4   )",       30936.5,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#else
  { "W(H2O	     ,Fe2SiO4   )",       28873.6,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#endif
  { "W(CaCO3     ,Fe2SiO4   )",      -72996.8,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "W(Mg2SiO4   ,MnSi0.5O2 )",      -13040.1,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(NiSi0.5O2 ,MnSi0.5O2 )", 	    785.8,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CoSi0.5O2 ,MnSi0.5O2 )", 	    -50.6,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CaSiO3	 ,MnSi0.5O2 )",        2934.6,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Na2SiO3   ,MnSi0.5O2 )",      -15780.8,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(KAlSiO4   ,MnSi0.5O2 )",       23727.4,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Ca3(PO4)2 ,MnSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CO2	     ,MnSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(S	     ,MnSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Cl	     ,MnSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(F	     ,MnSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(H2O	     ,MnSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CaCO3     ,MnSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "W(NiSi0.5O2 ,Mg2SiO4   )",      -21175.5,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CoSi0.5O2 ,Mg2SiO4   )",      -14994.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CaSiO3	 ,Mg2SiO4   )",      -31731.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Na2SiO3   ,Mg2SiO4   )",      -41876.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(KAlSiO4   ,Mg2SiO4   )",       22323.1,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Ca3(PO4)2 ,Mg2SiO4   )",      -23208.8,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CO2	     ,Mg2SiO4   )", 	 -40853.6,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(S	     ,Mg2SiO4   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Cl	     ,Mg2SiO4   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(F	     ,Mg2SiO4   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#ifdef USE_NEW_MELTS_WATER_MODEL
  { "W(H2O	     ,Mg2SiO4   )",       20909.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#else
  { "W(H2O	     ,Mg2SiO4   )",       35633.7,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#endif
  { "W(CaCO3     ,Mg2SiO4   )",      -24872.6,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "W(CoSi0.5O2 ,NiSi0.5O2 )", 	    258.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CaSiO3	 ,NiSi0.5O2 )",        7027.5,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Na2SiO3   ,NiSi0.5O2 )",       -3647.8,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(KAlSiO4   ,NiSi0.5O2 )",        4261.4,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Ca3(PO4)2 ,NiSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CO2	     ,NiSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(S	     ,NiSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Cl	     ,NiSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(F	     ,NiSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(H2O	     ,NiSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CaCO3     ,NiSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "W(CaSiO3	 ,CoSi0.5O2 )",      -26685.7,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Na2SiO3   ,CoSi0.5O2 )", 	    531.2,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(KAlSiO4   ,CoSi0.5O2 )", 	    265.7,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */	 
  { "W(Ca3(PO4)2 ,CoSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */	
  { "W(CO2	     ,CoSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(S	     ,CoSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Cl	     ,CoSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(F	     ,CoSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(H2O	     ,CoSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CaCO3     ,CoSi0.5O2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
    					     						    		
  { "W(Na2SiO3   ,CaSiO3    )",      -13247.1,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(KAlSiO4   ,CaSiO3    )",       17111.1,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Ca3(PO4)2 ,CaSiO3    )",       37070.3,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CO2	     ,CaSiO3    )", 	  30012.5,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(S	     ,CaSiO3    )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Cl	     ,CaSiO3    )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(F	     ,CaSiO3    )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#ifdef USE_NEW_MELTS_WATER_MODEL
  { "W(H2O	     ,CaSiO3    )",        9715.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#else
  { "W(H2O	     ,CaSiO3    )",       20374.6,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#endif
  { "W(CaCO3     ,CaSiO3    )",       37534.1,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "W(KAlSiO4   ,Na2SiO3   )",        6522.8,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Ca3(PO4)2 ,Na2SiO3   )",       15571.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CO2	     ,Na2SiO3   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(S	     ,Na2SiO3   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Cl	     ,Na2SiO3   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(F	     ,Na2SiO3   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#ifdef USE_NEW_MELTS_WATER_MODEL
  { "W(H2O	     ,Na2SiO3   )",      -82460.1,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#else
  { "W(H2O	     ,Na2SiO3   )",      -96937.6,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#endif
  { "W(CaCO3     ,Na2SiO3   )",     -311011.3,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "W(Ca3(PO4)2 ,KAlSiO4   )",       17100.6,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CO2	     ,KAlSiO4   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(S	     ,KAlSiO4   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Cl	     ,KAlSiO4   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(F	     ,KAlSiO4   )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#ifdef USE_NEW_MELTS_WATER_MODEL
  { "W(H2O	     ,KAlSiO4   )",        1057.2,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#else
  { "W(H2O	     ,KAlSiO4   )",       10374.2,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#endif
  { "W(CaCO3     ,KAlSiO4   )",      -27865.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "W(CO2	     ,Ca3(PO4)2 )", 	  -3472.8,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(S	     ,Ca3(PO4)2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Cl	     ,Ca3(PO4)2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(F	     ,Ca3(PO4)2 )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#ifdef USE_NEW_MELTS_WATER_MODEL
  { "W(H2O	     ,Ca3(PO4)2 )",       44133.3,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#else
  { "W(H2O	     ,Ca3(PO4)2 )",       43451.3,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
#endif
  { "W(CaCO3     ,Ca3(PO4)2 )",        2011.9,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "W(S	     ,CO2	    )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(Cl	     ,CO2	    )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(F	     ,CO2	    )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(H2O	     ,CO2	    )",       23255.4,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CaCO3     ,CO2	    )",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "W(Cl	     ,S	        )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(F	     ,S	        )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(H2O	     ,S	        )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CaCO3     ,S	        )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "W(F	     ,Cl	    )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(H2O	     ,Cl	    )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CaCO3     ,Cl	    )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "W(H2O	     ,F	        )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "W(CaCO3     ,F	        )", 	      0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "W(CaCO3     ,H2O       )", 	   7873.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "SiO2                    ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "TiO2       	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "Al2O3      	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "Fe2O3      	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "MgCr2O4    	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "Fe2SiO4    	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "MnSi0.5O2  	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "Mg2SiO4    	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "NiSi0.5O2  	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "CoSi0.5O2  	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "CaSiO3     	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "Na2SiO3    	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "KAlSiO4    	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "Ca3(PO4)2  	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "CO2        	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "S          	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "Cl         	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "F          	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */
  { "H2O        	         ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  { "CaCO3     	             ",           0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Basis species */

  {"olivine",		    		          0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"tephroite", 	    		          0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"fayalite",  	    		          0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"co-olivine",	    		          0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"ni-olivine",	    		          0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"monticellite",	    		          0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"forsterite",	    		          0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"fayalite",  	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"sphene",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"garnet",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"almandine", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"grossular", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"pyrope",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"melilite",  	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"akermanite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"gehlenite", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"iron-akermanite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"soda-melilite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"orthopyroxene",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"diopside",  	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"clinoenstatite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"hedenbergite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"alumino-buffonite",     		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"buffonite", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"essenite",  	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"jadeite",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"clinopyroxene",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"diopside",  	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"clinoenstatite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"hedenbergite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"alumino-buffonite",     		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"buffonite", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"essenite",  	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"jadeite",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"aegirine",  	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"aenigmatite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"cummingtonite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"cummingtonite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"grunerite", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"amphibole", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"cummingtonite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"grunerite", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"tremolite", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"hornblende",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"pargasite", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"ferropargasite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"magnesiohastingsite",   		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"biotite",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"annite",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"phlogopite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"muscovite", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"feldspar",  	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"albite",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"anorthite", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"sanidine",  	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"quartz",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"tridymite", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"cristobalite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"nepheline", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"na-nepheline",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"k-nepheline",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"vc-nepheline",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"ca-nepheline",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"kalsilite", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"na-nepheline",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"k-nepheline",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"vc-nepheline",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"ca-nepheline",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"leucite",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"leucite",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"analcime",  	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"na-leucite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"corundum",  	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"sillimanite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"rutile",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"perovskite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"spinel",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"chromite",  	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"hercynite", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"magnetite", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"spinel",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"ulvospinel",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"rhm-oxide", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"geikielite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"hematite",  	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"ilmenite",  	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"pyrophanite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"corundum",	    		          0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"ortho-oxide",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"pseudobrookite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"ferropseudobrookite",   		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"karrooite", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"whitlockite",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"apatite",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"water",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"fluid",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"h2oduan",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"co2duan",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"h2oduan",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"co2duan",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"alloy-solid",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"Fe-metal",  	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"Ni-metal",  	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"alloy-liquid",	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"Fe-liquid", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */
  {"Ni-liquid", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"lime",		    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }, /* Solid phase adjustment */

  {"periclase", 	    		  0.0,  0.0,  0.0, FALSE, FALSE, FALSE, FALSE }  /* Solid phase adjustment */
};

ModelParameters *modelParameters = originalModelParameters;

EosModelParameters eosModelParameters[] = {  /* sets nc entries */
                             /*  Kp   Kpp   Kppp */
  { "SiO2		     ", 0.188587,  1.91476e-05,  6.26232e-10, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "TiO2		     ",  -0.8548,  2.37406e-05, -1.69681e-09, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "Al2O3		     ", 0.720748, -1.22963e-05,   1.7757e-10, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "Fe2O3		     ",-0.516207, -2.77847e-05,  9.13701e-10, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 }, /* -0.516207, -2.77847e-05,  9.13701e-10 */ /*  1.0 */
  { "Cr2O3		     ",   2.8153,  1.05776e-05,  3.40393e-10, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "FeO		     ",  10.8578,  0.000326056, -9.59294e-09, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "MnO		     ",  10.5879,  0.000306443, -1.03951e-08, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "MgO		     ", -2.43032,  0.000271175, -3.12788e-08, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "NiO		     ",  28.8151,   0.00300134,  4.01919e-07, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "CoO		     ",  11.2075,  0.000310096, -1.22111e-08, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "CaO		     ",  18.0476,  0.000891667,  -2.0328e-08, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "Na2O		     ",   22.208,   0.00711102,  3.17635e-06, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "K2O		     ",  4.25535,  0.000752301,  8.46257e-08, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "P2O5		     ",   2.6445,  4.97721e-05,  2.73194e-10, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "H2O		     ",-0.653589, -7.46073e-06,  -8.1586e-12, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "CO2		     ", 0.930718,  4.37796e-06, -4.98837e-10, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "SO3		     ",  2.44963,  4.33544e-05,  9.62372e-11, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "Cl2O-1		     ",  19.2377,   0.00129928,    1.225e-07, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "F2O-1		     ",  3.10454,  6.60495e-05,  7.88912e-10, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 },   
  { "FeO1.3                  ",  1.47474, -7.28951e-05, -2.50008e-09, FALSE, FALSE, FALSE, 0.0, 0.0, 0.0 }  /* 1.47474, -7.28951e-05, -2.50008e-09 */ /* 12.0 */
}; 

int    nCN   = 0;    /* Total number of coordination states beyond reference                 */
double fCN[] = { };  /* configurational collapse term for IV->V, IV->IV oxygen CN shift      */
SsCnModelParameters ssCnModelParameters[] = { };  /* sets nc entries, repeated nCN-1 times   */

#endif /* _Param_Struct_Data_h */
