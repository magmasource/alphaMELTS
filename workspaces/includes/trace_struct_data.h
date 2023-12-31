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
