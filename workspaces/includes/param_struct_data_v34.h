#ifndef _Param_Struct_Data_v34_h
#define _Param_Struct_Data_v34_h

/*
MELTS Source Code: RCS $Log: param_struct_data_v34.h,v $
MELTS Source Code: RCS Revision 1.1.1.1  2006/08/15 16:57:36  ghiorso
MELTS Source Code: RCS xMELTS gcc 3.x sources
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2004/01/02 19:21:49  cvsaccount
MELTS Source Code: RCS CTserver University of Chicago
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1.1.1  2001/12/20 03:25:03  ghiorso
MELTS Source Code: RCS Sources for MELTS 5.x (xMELTS)
MELTS Source Code: RCS
*/

/*
**++
**  FACILITY:  Thread Safe Silicate Melts Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Include file for initializing values of liquidinteraction
**         parameters (file: PARAM_STRUCT_DATA_V34.H)
**--
*/

/*
 *=============================================================================
 * Initialize global values of interaction parameters
 */

static ModelParameters meltsModelParameters[] = {
  { "W(TiO2      ,SiO2      )",       26266.7,           0.0,      0.0 },  /*   0 */
  { "W(Al2O3     ,SiO2      )",      -39120.0,           0.0,      0.0 },  /*   1 */
  { "W(Fe2O3     ,SiO2      )",        8110.3,           0.0,      0.0 },  /*   2 */
  { "W(MgCr2O4   ,SiO2      )",       27886.3,           0.0,      0.0 },  /*   3 */
  { "W(Fe2SiO4   ,SiO2      )",       23660.9,           0.0,      0.0 },  /*   4 */
  { "W(MnSi0.5O2 ,SiO2      )",       18393.9,           0.0,      0.0 },  /*   5 */
  { "W(Mg2SiO4   ,SiO2      )",        3421.0,           0.0,      0.0 },  /*   6 */
  { "W(NiSi0.5O2 ,SiO2      )",       25197.4,           0.0,      0.0 },  /*   7 */
  { "W(CoSi0.5O2 ,SiO2      )",       14802.8,           0.0,      0.0 },  /*   8 */
  { "W(CaSiO3    ,SiO2      )",        -863.7,           0.0,      0.0 },  /*   9 */
  { "W(Na2SiO3   ,SiO2      )",      -99039.0,           0.0,      0.0 },  /*  10 */
  { "W(KAlSiO4   ,SiO2      )",      -33921.7,           0.0,      0.0 },  /*  11 */
  { "W(Ca3(PO4)2 ,SiO2      )",       61891.6,           0.0,      0.0 },  /*  12 */
  { "W(CO2       ,SiO2      )",       66600.7,           0.0,      0.0 },  /*  13 */
  { "W(S         ,SiO2      )",           0.0,           0.0,      0.0 },  /*  14 */
  { "W(Cl        ,SiO2      )",           0.0,           0.0,      0.0 },  /*  15 */
  { "W(F         ,SiO2      )",           0.0,           0.0,      0.0 },  /*  16 */
  { "W(H2O       ,SiO2      )",       30967.3,           0.0,      0.0 },  /*  17 */

  { "W(Al2O3     ,TiO2      )",      -29449.8,           0.0,      0.0 },  /*  18 */
  { "W(Fe2O3     ,TiO2      )",      -84756.9,           0.0,      0.0 },  /*  19 */
  { "W(MgCr2O4   ,TiO2      )",      -72303.4,           0.0,      0.0 },  /*  20 */
  { "W(Fe2SiO4   ,TiO2      )",        5209.1,           0.0,      0.0 },  /*  21 */
  { "W(MnSi0.5O2 ,TiO2      )",      -16123.5,           0.0,      0.0 },  /*  22 */
  { "W(Mg2SiO4   ,TiO2      )",       -4178.3,           0.0,      0.0 },  /*  23 */
  { "W(NiSi0.5O2 ,TiO2      )",        3614.8,           0.0,      0.0 },  /*  24 */
  { "W(CoSi0.5O2 ,TiO2      )",       -1640.0,           0.0,      0.0 },  /*  25 */
  { "W(CaSiO3    ,TiO2      )",      -35372.5,           0.0,      0.0 },  /*  26 */
  { "W(Na2SiO3   ,TiO2      )",      -15415.6,           0.0,      0.0 },  /*  27 */
  { "W(KAlSiO4   ,TiO2      )",      -48094.6,           0.0,      0.0 },  /*  28 */
  { "W(Ca3(PO4)2 ,TiO2      )",       25938.8,           0.0,      0.0 },  /*  29 */
  { "W(CO2       ,TiO2      )",       25427.3,           0.0,      0.0 },  /*  30 */
  { "W(S         ,TiO2      )",           0.0,           0.0,      0.0 },  /*  31 */
  { "W(Cl        ,TiO2      )",           0.0,           0.0,      0.0 },  /*  32 */
  { "W(F         ,TiO2      )",           0.0,           0.0,      0.0 },  /*  33 */
  { "W(H2O       ,TiO2      )",       81879.1,           0.0,      0.0 },  /*  34 */

  { "W(Fe2O3     ,Al2O3     )",      -17089.4,           0.0,      0.0 },  /*  35 */
  { "W(MgCr2O4   ,Al2O3     )",      -31770.3,           0.0,      0.0 },  /*  36 */
  { "W(Fe2SiO4   ,Al2O3     )",      -30509.0,           0.0,      0.0 },  /*  37 */
  { "W(MnSi0.5O2 ,Al2O3     )",      -53874.9,           0.0,      0.0 },  /*  38 */
  { "W(Mg2SiO4   ,Al2O3     )",      -32880.3,           0.0,      0.0 },  /*  39 */
  { "W(NiSi0.5O2 ,Al2O3     )",        2985.2,           0.0,      0.0 },  /*  40 */
  { "W(CoSi0.5O2 ,Al2O3     )",       -2677.4,           0.0,      0.0 },  /*  41 */
  { "W(CaSiO3    ,Al2O3     )",      -57917.9,           0.0,      0.0 },  /*  42 */
  { "W(Na2SiO3   ,Al2O3     )",     -130785.0,           0.0,      0.0 },  /*  43 */
  { "W(KAlSiO4   ,Al2O3     )",      -25859.2,           0.0,      0.0 },  /*  44 */
  { "W(Ca3(PO4)2 ,Al2O3     )",       52220.8,           0.0,      0.0 },  /*  45 */
  { "W(CO2       ,Al2O3     )",       27043.7,           0.0,      0.0 },  /*  46 */
  { "W(S         ,Al2O3     )",           0.0,           0.0,      0.0 },  /*  47 */
  { "W(Cl        ,Al2O3     )",           0.0,           0.0,      0.0 },  /*  48 */
  { "W(F         ,Al2O3     )",           0.0,           0.0,      0.0 },  /*  49 */
  { "W(H2O       ,Al2O3     )",      -16098.1,           0.0,      0.0 },  /*  50 */

  { "W(MgCr2O4   ,Fe2O3     )",       21605.9,           0.0,      0.0 },  /*  51 */
  { "W(Fe2SiO4   ,Fe2O3     )",     -179064.9,           0.0,      0.0 },  /*  52 */
  { "W(MnSi0.5O2 ,Fe2O3     )",        3907.9,           0.0,      0.0 },  /*  53 */
  { "W(Mg2SiO4   ,Fe2O3     )",      -71518.6,           0.0,      0.0 },  /*  54 */
  { "W(NiSi0.5O2 ,Fe2O3     )",         408.7,           0.0,      0.0 },  /*  55 */
  { "W(CoSi0.5O2 ,Fe2O3     )",        -223.7,           0.0,      0.0 },  /*  56 */
  { "W(CaSiO3    ,Fe2O3     )",       12076.6,           0.0,      0.0 },  /*  57 */
  { "W(Na2SiO3   ,Fe2O3     )",     -149662.2,           0.0,      0.0 },  /*  58 */
  { "W(KAlSiO4   ,Fe2O3     )",       57555.9,           0.0,      0.0 },  /*  59 */
  { "W(Ca3(PO4)2 ,Fe2O3     )",       -4213.9,           0.0,      0.0 },  /*  60 */
  { "W(CO2       ,Fe2O3     )",        1688.8,           0.0,      0.0 },  /*  61 */
  { "W(S         ,Fe2O3     )",           0.0,           0.0,      0.0 },  /*  62 */
  { "W(Cl        ,Fe2O3     )",           0.0,           0.0,      0.0 },  /*  63 */
  { "W(F         ,Fe2O3     )",           0.0,           0.0,      0.0 },  /*  64 */
  { "W(H2O       ,Fe2O3     )",       31405.5,           0.0,      0.0 },  /*  65 */

  { "W(Fe2SiO4   ,MgCr2O4   )",      -82971.8,           0.0,      0.0 },  /*  66 */
  { "W(MnSi0.5O2 ,MgCr2O4   )",         182.4,           0.0,      0.0 },  /*  67 */
  { "W(Mg2SiO4   ,MgCr2O4   )",       46049.2,           0.0,      0.0 },  /*  68 */
  { "W(NiSi0.5O2 ,MgCr2O4   )",        -266.0,           0.0,      0.0 },  /*  69 */
  { "W(CoSi0.5O2 ,MgCr2O4   )",        -384.0,           0.0,      0.0 },  /*  70 */
  { "W(CaSiO3    ,MgCr2O4   )",       30704.7,           0.0,      0.0 },  /*  71 */
  { "W(Na2SiO3   ,MgCr2O4   )",      113646.0,           0.0,      0.0 },  /*  72 */
  { "W(KAlSiO4   ,MgCr2O4   )",       75709.1,           0.0,      0.0 },  /*  73 */
  { "W(Ca3(PO4)2 ,MgCr2O4   )",        5341.8,           0.0,      0.0 },  /*  74 */
  { "W(CO2       ,MgCr2O4   )",           0.0,           0.0,      0.0 },  /*  75 */
  { "W(S         ,MgCr2O4   )",           0.0,           0.0,      0.0 },  /*  76 */
  { "W(Cl        ,MgCr2O4   )",           0.0,           0.0,      0.0 },  /*  77 */
  { "W(F         ,MgCr2O4   )",           0.0,           0.0,      0.0 },  /*  78 */
  { "W(H2O       ,MgCr2O4   )",           0.0,           0.0,      0.0 },  /*  79 */

  { "W(MnSi0.5O2 ,Fe2SiO4   )",       -6823.9,           0.0,      0.0 },  /*  80 */
  { "W(Mg2SiO4   ,Fe2SiO4   )",      -37256.7,           0.0,      0.0 },  /*  81 */
  { "W(NiSi0.5O2 ,Fe2SiO4   )",      -17019.8,           0.0,      0.0 },  /*  82 */
  { "W(CoSi0.5O2 ,Fe2SiO4   )",      -11746.3,           0.0,      0.0 },  /*  83 */
  { "W(CaSiO3    ,Fe2SiO4   )",      -12970.8,           0.0,      0.0 },  /*  84 */
  { "W(Na2SiO3   ,Fe2SiO4   )",      -90533.8,           0.0,      0.0 },  /*  85 */
  { "W(KAlSiO4   ,Fe2SiO4   )",       23649.4,           0.0,      0.0 },  /*  86 */
  { "W(Ca3(PO4)2 ,Fe2SiO4   )",       87410.3,           0.0,      0.0 },  /*  87 */
  { "W(CO2       ,Fe2SiO4   )",      134680.2,           0.0,      0.0 },  /*  88 */
  { "W(S         ,Fe2SiO4   )",           0.0,           0.0,      0.0 },  /*  89 */
  { "W(Cl        ,Fe2SiO4   )",           0.0,           0.0,      0.0 },  /*  90 */
  { "W(F         ,Fe2SiO4   )",           0.0,           0.0,      0.0 },  /*  91 */
  { "W(H2O       ,Fe2SiO4   )",       28873.6,           0.0,      0.0 },  /*  92 */

  { "W(Mg2SiO4   ,MnSi0.5O2 )",      -13040.1,           0.0,      0.0 },  /*  93 */
  { "W(NiSi0.5O2 ,MnSi0.5O2 )",         785.8,           0.0,      0.0 },  /*  94 */
  { "W(CoSi0.5O2 ,MnSi0.5O2 )",         -50.6,           0.0,      0.0 },  /*  95 */
  { "W(CaSiO3    ,MnSi0.5O2 )",        2934.6,           0.0,      0.0 },  /*  96 */
  { "W(Na2SiO3   ,MnSi0.5O2 )",      -15780.8,           0.0,      0.0 },  /*  97 */
  { "W(KAlSiO4   ,MnSi0.5O2 )",       23727.4,           0.0,      0.0 },  /*  98 */
  { "W(Ca3(PO4)2 ,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /*  99 */
  { "W(CO2       ,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 100 */
  { "W(S         ,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 101 */
  { "W(Cl        ,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 102 */
  { "W(F         ,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 103 */
  { "W(H2O       ,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 104 */

  { "W(NiSi0.5O2 ,Mg2SiO4   )",      -21175.5,           0.0,      0.0 },  /* 105 */
  { "W(CoSi0.5O2 ,Mg2SiO4   )",      -14994.9,           0.0,      0.0 },  /* 106 */
  { "W(CaSiO3    ,Mg2SiO4   )",      -31731.9,           0.0,      0.0 },  /* 107 */
  { "W(Na2SiO3   ,Mg2SiO4   )",      -41876.9,           0.0,      0.0 },  /* 108 */
  { "W(KAlSiO4   ,Mg2SiO4   )",       22323.1,           0.0,      0.0 },  /* 109 */
  { "W(Ca3(PO4)2 ,Mg2SiO4   )",      -23208.8,           0.0,      0.0 },  /* 110 */
  { "W(CO2       ,Mg2SiO4   )",       47406.1,           0.0,      0.0 },  /* 111 */
  { "W(S         ,Mg2SiO4   )",           0.0,           0.0,      0.0 },  /* 112 */
  { "W(Cl        ,Mg2SiO4   )",           0.0,           0.0,      0.0 },  /* 113 */
  { "W(F         ,Mg2SiO4   )",           0.0,           0.0,      0.0 },  /* 114 */
  { "W(H2O       ,Mg2SiO4   )",       35633.7,           0.0,      0.0 },  /* 115 */

  { "W(CoSi0.5O2 ,NiSi0.5O2 )",         258.9,           0.0,      0.0 },  /* 116 */
  { "W(CaSiO3    ,NiSi0.5O2 )",        7027.5,           0.0,      0.0 },  /* 117 */
  { "W(Na2SiO3   ,NiSi0.5O2 )",       -3647.8,           0.0,      0.0 },  /* 118 */
  { "W(KAlSiO4   ,NiSi0.5O2 )",        4261.4,           0.0,      0.0 },  /* 119 */
  { "W(Ca3(PO4)2 ,NiSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 120 */
  { "W(CO2       ,NiSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 121 */
  { "W(S         ,NiSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 122 */
  { "W(Cl        ,NiSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 123 */
  { "W(F         ,NiSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 124 */
  { "W(H2O       ,NiSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 125 */

  { "W(CaSiO3    ,CoSi0.5O2 )",      -26685.7,           0.0,      0.0 },  /* 126 */
  { "W(Na2SiO3   ,CoSi0.5O2 )",         531.2,           0.0,      0.0 },  /* 127 */
  { "W(KAlSiO4   ,CoSi0.5O2 )",         265.7,           0.0,      0.0 },  /* 128 */
  { "W(Ca3(PO4)2 ,CoSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 129 */
  { "W(CO2       ,CoSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 130 */
  { "W(S         ,CoSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 131 */
  { "W(Cl        ,CoSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 132 */
  { "W(F         ,CoSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 133 */
  { "W(H2O       ,CoSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 134 */

  { "W(Na2SiO3   ,CaSiO3    )",      -13247.1,           0.0,      0.0 },  /* 135 */
  { "W(KAlSiO4   ,CaSiO3    )",       17111.1,           0.0,      0.0 },  /* 136 */
  { "W(Ca3(PO4)2 ,CaSiO3    )",       37070.3,           0.0,      0.0 },  /* 137 */
  { "W(CO2       ,CaSiO3    )",       44513.8,           0.0,      0.0 },  /* 138 */
  { "W(S         ,CaSiO3    )",           0.0,           0.0,      0.0 },  /* 139 */
  { "W(Cl        ,CaSiO3    )",           0.0,           0.0,      0.0 },  /* 140 */
  { "W(F         ,CaSiO3    )",           0.0,           0.0,      0.0 },  /* 141 */
  { "W(H2O       ,CaSiO3    )",       20374.6,           0.0,      0.0 },  /* 142 */

  { "W(KAlSiO4   ,Na2SiO3   )",        6522.8,           0.0,      0.0 },  /* 143 */
  { "W(Ca3(PO4)2 ,Na2SiO3   )",       15571.9,           0.0,      0.0 },  /* 144 */
  { "W(CO2       ,Na2SiO3   )",        9846.7,           0.0,      0.0 },  /* 145 */
  { "W(S         ,Na2SiO3   )",           0.0,           0.0,      0.0 },  /* 146 */
  { "W(Cl        ,Na2SiO3   )",           0.0,           0.0,      0.0 },  /* 147 */
  { "W(F         ,Na2SiO3   )",           0.0,           0.0,      0.0 },  /* 148 */
  { "W(H2O       ,Na2SiO3   )",      -96937.6,           0.0,      0.0 },  /* 149 */

  { "W(Ca3(PO4)2 ,KAlSiO4   )",       17100.6,           0.0,      0.0 },  /* 150 */
  { "W(CO2       ,KAlSiO4   )",        1744.4,           0.0,      0.0 },  /* 151 */
  { "W(S         ,KAlSiO4   )",           0.0,           0.0,      0.0 },  /* 152 */
  { "W(Cl        ,KAlSiO4   )",           0.0,           0.0,      0.0 },  /* 153 */
  { "W(F         ,KAlSiO4   )",           0.0,           0.0,      0.0 },  /* 154 */
  { "W(H2O       ,KAlSiO4   )",       10374.2,           0.0,      0.0 },  /* 155 */

  { "W(CO2       ,Ca3(PO4)2 )",           0.0,           0.0,      0.0 },  /* 156 */
  { "W(S         ,Ca3(PO4)2 )",           0.0,           0.0,      0.0 },  /* 157 */
  { "W(Cl        ,Ca3(PO4)2 )",           0.0,           0.0,      0.0 },  /* 158 */
  { "W(F         ,Ca3(PO4)2 )",           0.0,           0.0,      0.0 },  /* 159 */
  { "W(H2O       ,Ca3(PO4)2 )",       43451.3,           0.0,      0.0 },  /* 160 */

  { "W(S         ,CO2       )",           0.0,           0.0,      0.0 },  /* 161 */
  { "W(Cl        ,CO2       )",           0.0,           0.0,      0.0 },  /* 162 */
  { "W(F         ,CO2       )",           0.0,           0.0,      0.0 },  /* 163 */
  { "W(H2O       ,CO2       )",       43128.4,           0.0,      0.0 },  /* 164 */

  { "W(Cl        ,S         )",           0.0,           0.0,      0.0 },  /* 165 */
  { "W(F         ,S         )",           0.0,           0.0,      0.0 },  /* 166 */
  { "W(H2O       ,S         )",           0.0,           0.0,      0.0 },  /* 167 */

  { "W(F         ,Cl        )",           0.0,           0.0,      0.0 },  /* 168 */
  { "W(H2O       ,Cl        )",           0.0,           0.0,      0.0 },  /* 169 */

  { "W(H2O       ,F         )",           0.0,           0.0,      0.0 }   /* 170 */
};

static const ModelParameters pMeltsModelParameters[] = {
  { "W(TiO2      ,Si4O8     )",       15094.7,           0.0,      0.0 },  /* 000 */
  { "W(Al4O6     ,Si4O8     )",     -296975.2,           0.0,      0.0 },  /* 001 */
  { "W(Fe2O3     ,Si4O8     )",     -164027.4,           0.0,      0.0 },  /* 002 */
  { "W(MgCr2O4   ,Si4O8     )",       37459.2,           0.0,      0.0 },  /* 003 */
  { "W(Fe2SiO4   ,Si4O8     )",      -18841.4,           0.0,      0.0 },  /* 004 */
  { "W(MnSi0.5O2 ,Si4O8     )",           0.0,           0.0,      0.0 },  /* 005 */
  { "W(Mg2SiO4   ,Si4O8     )",      -33833.5,           0.0,      0.0 },  /* 006 */
  { "W(NiSi0.5O2 ,Si4O8     )",           0.0,           0.0,      0.0 },  /* 007 */
  { "W(CoSi0.5O2 ,Si4O8     )",           0.0,           0.0,      0.0 },  /* 008 */
  { "W(Ca2Si2O6  ,Si4O8     )",      -34232.9,           0.0,      0.0 },  /* 009 */
  { "W(NaSi0.5O1.,Si4O8     )",      -59822.7,           0.0,      0.0 },  /* 010 */
  { "W(KAlSiO4   ,Si4O8     )",     -102706.5,           0.0,      0.0 },  /* 011 */
  { "W(Ca3(PO4)2 ,Si4O8     )",       37519.9,           0.0,      0.0 },  /* 012 */
  { "W(CO2       ,Si4O8     )",           0.0,           0.0,      0.0 },  /* 013 */
  { "W(SO3       ,Si4O8     )",           0.0,           0.0,      0.0 },  /* 014 */
  { "W(Cl2O-1    ,Si4O8     )",           0.0,           0.0,      0.0 },  /* 015 */
  { "W(F2O-1     ,Si4O8     )",           0.0,           0.0,      0.0 },  /* 016 */
  { "W(H2O       ,Si4O8     )",      -45181.6,           0.0,      0.0 },  /* 017 */

  { "W(Al4O6     ,TiO2      )",     -144804.9,           0.0,      0.0 },  /* 018 */
  { "W(Fe2O3     ,TiO2      )",     -212292.3,           0.0,      0.0 },  /* 019 */
  { "W(MgCr2O4   ,TiO2      )",      -22455.8,           0.0,      0.0 },  /* 020 */
  { "W(Fe2SiO4   ,TiO2      )",        9324.2,           0.0,      0.0 },  /* 021 */
  { "W(MnSi0.5O2 ,TiO2      )",           0.0,           0.0,      0.0 },  /* 022 */
  { "W(Mg2SiO4   ,TiO2      )",       16335.6,           0.0,      0.0 },  /* 023 */
  { "W(NiSi0.5O2 ,TiO2      )",           0.0,           0.0,      0.0 },  /* 024 */
  { "W(CoSi0.5O2 ,TiO2      )",           0.0,           0.0,      0.0 },  /* 025 */
  { "W(Ca2Si2O6  ,TiO2      )",       -9471.5,           0.0,      0.0 },  /* 026 */
  { "W(NaSi0.5O1.,TiO2      )",       22194.2,           0.0,      0.0 },  /* 027 */
  { "W(KAlSiO4   ,TiO2      )",       -3744.0,           0.0,      0.0 },  /* 028 */
  { "W(Ca3(PO4)2 ,TiO2      )",       65544.0,           0.0,      0.0 },  /* 029 */
  { "W(CO2       ,TiO2      )",           0.0,           0.0,      0.0 },  /* 030 */
  { "W(SO3       ,TiO2      )",           0.0,           0.0,      0.0 },  /* 031 */
  { "W(Cl2O-1    ,TiO2      )",           0.0,           0.0,      0.0 },  /* 032 */
  { "W(F2O-1     ,TiO2      )",           0.0,           0.0,      0.0 },  /* 033 */
  { "W(H2O       ,TiO2      )",       70663.0,           0.0,      0.0 },  /* 034 */

  { "W(Fe2O3     ,Al4O6     )",     -393566.0,           0.0,      0.0 },  /* 035 */
  { "W(MgCr2O4   ,Al4O6     )",     -269339.7,           0.0,      0.0 },  /* 036 */
  { "W(Fe2SiO4   ,Al4O6     )",     -200788.1,           0.0,      0.0 },  /* 037 */
  { "W(MnSi0.5O2 ,Al4O6     )",           0.0,           0.0,      0.0 },  /* 038 */
  { "W(Mg2SiO4   ,Al4O6     )",     -192709.0,           0.0,      0.0 },  /* 039 */
  { "W(NiSi0.5O2 ,Al4O6     )",           0.0,           0.0,      0.0 },  /* 040 */
  { "W(CoSi0.5O2 ,Al4O6     )",           0.0,           0.0,      0.0 },  /* 041 */
  { "W(Ca2Si2O6  ,Al4O6     )",     -270700.8,           0.0,      0.0 },  /* 042 */
  { "W(NaSi0.5O1.,Al4O6     )",     -205068.6,           0.0,      0.0 },  /* 043 */
  { "W(KAlSiO4   ,Al4O6     )",     -114506.5,           0.0,      0.0 },  /* 044 */
  { "W(Ca3(PO4)2 ,Al4O6     )",     -176584.1,           0.0,      0.0 },  /* 045 */
  { "W(CO2       ,Al4O6     )",           0.0,           0.0,      0.0 },  /* 046 */
  { "W(SO3       ,Al4O6     )",           0.0,           0.0,      0.0 },  /* 047 */
  { "W(Cl2O-1    ,Al4O6     )",           0.0,           0.0,      0.0 },  /* 048 */
  { "W(F2O-1     ,Al4O6     )",           0.0,           0.0,      0.0 },  /* 049 */
  { "W(H2O       ,Al4O6     )",     -161944.4,           0.0,      0.0 },  /* 050 */

  { "W(MgCr2O4   ,Fe2O3     )",      201536.3,           0.0,      0.0 },  /* 051 */
  { "W(Fe2SiO4   ,Fe2O3     )",     -211493.4,           0.0,      0.0 },  /* 052 */
  { "W(MnSi0.5O2 ,Fe2O3     )",           0.0,           0.0,      0.0 },  /* 053 */
  { "W(Mg2SiO4   ,Fe2O3     )",     -196914.9,           0.0,      0.0 },  /* 054 */
  { "W(NiSi0.5O2 ,Fe2O3     )",           0.0,           0.0,      0.0 },  /* 055 */
  { "W(CoSi0.5O2 ,Fe2O3     )",           0.0,           0.0,      0.0 },  /* 056 */
  { "W(Ca2Si2O6  ,Fe2O3     )",     -146008.1,           0.0,      0.0 },  /* 057 */
  { "W(NaSi0.5O1.,Fe2O3     )",     -123728.7,           0.0,      0.0 },  /* 058 */
  { "W(KAlSiO4   ,Fe2O3     )",     -130847.5,           0.0,      0.0 },  /* 059 */
  { "W(Ca3(PO4)2 ,Fe2O3     )",     -126339.8,           0.0,      0.0 },  /* 060 */
  { "W(CO2       ,Fe2O3     )",           0.0,           0.0,      0.0 },  /* 061 */
  { "W(SO3       ,Fe2O3     )",           0.0,           0.0,      0.0 },  /* 062 */
  { "W(Cl2O-1    ,Fe2O3     )",           0.0,           0.0,      0.0 },  /* 063 */
  { "W(F2O-1     ,Fe2O3     )",           0.0,           0.0,      0.0 },  /* 064 */
  { "W(H2O       ,Fe2O3     )",     -114508.6,           0.0,      0.0 },  /* 065 */

  { "W(Fe2SiO4   ,MgCr2O4   )",      -74759.0,           0.0,      0.0 },  /* 066 */
  { "W(MnSi0.5O2 ,MgCr2O4   )",           0.0,           0.0,      0.0 },  /* 067 */
  { "W(Mg2SiO4   ,MgCr2O4   )",       -3638.5,           0.0,      0.0 },  /* 068 */
  { "W(NiSi0.5O2 ,MgCr2O4   )",           0.0,           0.0,      0.0 },  /* 069 */
  { "W(CoSi0.5O2 ,MgCr2O4   )",           0.0,           0.0,      0.0 },  /* 070 */
  { "W(Ca2Si2O6  ,MgCr2O4   )",       48337.5,           0.0,      0.0 },  /* 071 */
  { "W(NaSi0.5O1.,MgCr2O4   )",      -43302.5,           0.0,      0.0 },  /* 072 */
  { "W(KAlSiO4   ,MgCr2O4   )",      124517.4,           0.0,      0.0 },  /* 073 */
  { "W(Ca3(PO4)2 ,MgCr2O4   )",       13004.3,           0.0,      0.0 },  /* 074 */
  { "W(CO2       ,MgCr2O4   )",           0.0,           0.0,      0.0 },  /* 075 */
  { "W(SO3       ,MgCr2O4   )",           0.0,           0.0,      0.0 },  /* 076 */
  { "W(Cl2O-1    ,MgCr2O4   )",           0.0,           0.0,      0.0 },  /* 077 */
  { "W(F2O-1     ,MgCr2O4   )",           0.0,           0.0,      0.0 },  /* 078 */
  { "W(H2O       ,MgCr2O4   )",         -18.9,           0.0,      0.0 },  /* 079 */

  { "W(MnSi0.5O2 ,Fe2SiO4   )",           0.0,           0.0,      0.0 },  /* 080 */
  { "W(Mg2SiO4   ,Fe2SiO4   )",      -28736.4,           0.0,      0.0 },  /* 081 */
  { "W(NiSi0.5O2 ,Fe2SiO4   )",           0.0,           0.0,      0.0 },  /* 082 */
  { "W(CoSi0.5O2 ,Fe2SiO4   )",           0.0,           0.0,      0.0 },  /* 083 */
  { "W(Ca2Si2O6  ,Fe2SiO4   )",      -28573.8,           0.0,      0.0 },  /* 084 */
  { "W(NaSi0.5O1.,Fe2SiO4   )",       -4723.9,           0.0,      0.0 },  /* 085 */
  { "W(KAlSiO4   ,Fe2SiO4   )",       22245.0,           0.0,      0.0 },  /* 086 */
  { "W(Ca3(PO4)2 ,Fe2SiO4   )",        4909.8,           0.0,      0.0 },  /* 087 */
  { "W(CO2       ,Fe2SiO4   )",           0.0,           0.0,      0.0 },  /* 088 */
  { "W(SO3       ,Fe2SiO4   )",           0.0,           0.0,      0.0 },  /* 089 */
  { "W(Cl2O-1    ,Fe2SiO4   )",           0.0,           0.0,      0.0 },  /* 090 */
  { "W(F2O-1     ,Fe2SiO4   )",           0.0,           0.0,      0.0 },  /* 091 */
  { "W(H2O       ,Fe2SiO4   )",        9769.4,           0.0,      0.0 },  /* 092 */

  { "W(Mg2SiO4   ,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 093 */
  { "W(NiSi0.5O2 ,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 094 */
  { "W(CoSi0.5O2 ,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 095 */
  { "W(Ca2Si2O6  ,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 096 */
  { "W(NaSi0.5O1.,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 097 */
  { "W(KAlSiO4   ,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 098 */
  { "W(Ca3(PO4)2 ,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 099 */
  { "W(CO2       ,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 100 */
  { "W(SO3       ,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 101 */
  { "W(Cl2O-1    ,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 102 */
  { "W(F2O-1     ,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 103 */
  { "W(H2O       ,MnSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 104 */

  { "W(NiSi0.5O2 ,Mg2SiO4   )",           0.0,           0.0,      0.0 },  /* 105 */
  { "W(CoSi0.5O2 ,Mg2SiO4   )",           0.0,           0.0,      0.0 },  /* 106 */
  { "W(Ca2Si2O6  ,Mg2SiO4   )",         574.1,           0.0,      0.0 },  /* 107 */
  { "W(NaSi0.5O1.,Mg2SiO4   )",        9272.3,           0.0,      0.0 },  /* 108 */
  { "W(KAlSiO4   ,Mg2SiO4   )",       36512.7,           0.0,      0.0 },  /* 109 */
  { "W(Ca3(PO4)2 ,Mg2SiO4   )",       -6766.8,           0.0,      0.0 },  /* 110 */
  { "W(CO2       ,Mg2SiO4   )",           0.0,           0.0,      0.0 },  /* 111 */
  { "W(SO3       ,Mg2SiO4   )",           0.0,           0.0,      0.0 },  /* 112 */
  { "W(Cl2O-1    ,Mg2SiO4   )",           0.0,           0.0,      0.0 },  /* 113 */
  { "W(F2O-1     ,Mg2SiO4   )",           0.0,           0.0,      0.0 },  /* 114 */
  { "W(H2O       ,Mg2SiO4   )",       24630.1,           0.0,      0.0 },  /* 115 */

  { "W(CoSi0.5O2 ,NiSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 116 */
  { "W(Ca2Si2O6  ,NiSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 117 */
  { "W(NaSi0.5O1.,NiSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 118 */
  { "W(KAlSiO4   ,NiSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 119 */
  { "W(Ca3(PO4)2 ,NiSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 120 */
  { "W(CO2       ,NiSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 121 */
  { "W(SO3       ,NiSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 122 */
  { "W(Cl2O-1    ,NiSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 123 */
  { "W(F2O-1     ,NiSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 124 */
  { "W(H2O       ,NiSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 125 */

  { "W(Ca2Si2O6  ,CoSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 126 */
  { "W(NaSi0.5O1.,CoSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 127 */
  { "W(KAlSiO4   ,CoSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 128 */
  { "W(Ca3(PO4)2 ,CoSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 129 */
  { "W(CO2       ,CoSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 130 */
  { "W(SO3       ,CoSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 131 */
  { "W(Cl2O-1    ,CoSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 132 */
  { "W(F2O-1     ,CoSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 133 */
  { "W(H2O       ,CoSi0.5O2 )",           0.0,           0.0,      0.0 },  /* 134 */

  { "W(NaSi0.5O1.,Ca2Si2O6  )",        7430.3,           0.0,      0.0 },  /* 135 */
  { "W(KAlSiO4   ,Ca2Si2O6  )",       19927.4,           0.0,      0.0 },  /* 136 */
  { "W(Ca3(PO4)2 ,Ca2Si2O6  )",       88993.1,           0.0,      0.0 },  /* 137 */
  { "W(CO2       ,Ca2Si2O6  )",           0.0,           0.0,      0.0 },  /* 138 */
  { "W(SO3       ,Ca2Si2O6  )",           0.0,           0.0,      0.0 },  /* 139 */
  { "W(Cl2O-1    ,Ca2Si2O6  )",           0.0,           0.0,      0.0 },  /* 140 */
  { "W(F2O-1     ,Ca2Si2O6  )",           0.0,           0.0,      0.0 },  /* 141 */
  { "W(H2O       ,Ca2Si2O6  )",       -1583.7,           0.0,      0.0 },  /* 142 */

  { "W(KAlSiO4   ,NaSi0.5O1.)",       -1102.3,           0.0,      0.0 },  /* 143 */
  { "W(Ca3(PO4)2 ,NaSi0.5O1.)",      -13062.6,           0.0,      0.0 },  /* 144 */
  { "W(CO2       ,NaSi0.5O1.)",           0.0,           0.0,      0.0 },  /* 145 */
  { "W(SO3       ,NaSi0.5O1.)",           0.0,           0.0,      0.0 },  /* 146 */
  { "W(Cl2O-1    ,NaSi0.5O1.)",           0.0,           0.0,      0.0 },  /* 147 */
  { "W(F2O-1     ,NaSi0.5O1.)",           0.0,           0.0,      0.0 },  /* 148 */
  { "W(H2O       ,NaSi0.5O1.)",       13043.1,           0.0,      0.0 },  /* 149 */

  { "W(Ca3(PO4)2 ,KAlSiO4   )",       85064.0,           0.0,      0.0 },  /* 150 */
  { "W(CO2       ,KAlSiO4   )",           0.0,           0.0,      0.0 },  /* 151 */
  { "W(SO3       ,KAlSiO4   )",           0.0,           0.0,      0.0 },  /* 152 */
  { "W(Cl2O-1    ,KAlSiO4   )",           0.0,           0.0,      0.0 },  /* 153 */
  { "W(F2O-1     ,KAlSiO4   )",           0.0,           0.0,      0.0 },  /* 154 */
  { "W(H2O       ,KAlSiO4   )",       35572.8,           0.0,      0.0 },  /* 155 */

  { "W(CO2       ,Ca3(PO4)2 )",           0.0,           0.0,      0.0 },  /* 156 */
  { "W(SO3       ,Ca3(PO4)2 )",           0.0,           0.0,      0.0 },  /* 157 */
  { "W(Cl2O-1    ,Ca3(PO4)2 )",           0.0,           0.0,      0.0 },  /* 158 */
  { "W(F2O-1     ,Ca3(PO4)2 )",           0.0,           0.0,      0.0 },  /* 159 */
  { "W(H2O       ,Ca3(PO4)2 )",       53448.7,           0.0,      0.0 },  /* 160 */

  { "W(SO3       ,CO2       )",           0.0,           0.0,      0.0 },  /* 161 */
  { "W(Cl2O-1    ,CO2       )",           0.0,           0.0,      0.0 },  /* 162 */
  { "W(F2O-1     ,CO2       )",           0.0,           0.0,      0.0 },  /* 163 */
  { "W(H2O       ,CO2       )",           0.0,           0.0,      0.0 },  /* 164 */

  { "W(Cl2O-1    ,SO3       )",           0.0,           0.0,      0.0 },  /* 165 */
  { "W(F2O-1     ,SO3       )",           0.0,           0.0,      0.0 },  /* 166 */
  { "W(H2O       ,SO3       )",           0.0,           0.0,      0.0 },  /* 167 */

  { "W(F2O-1     ,Cl2O-1    )",           0.0,           0.0,      0.0 },  /* 168 */
  { "W(H2O       ,Cl2O-1    )",           0.0,           0.0,      0.0 },  /* 169 */

  { "W(H2O       ,F2O-1     )",           0.0,           0.0,      0.0 }   /* 170 */
};

#endif /* _Param_Struct_Data_v34_h */