# Exercise: isobaric fractional crystallization of a back-arc basalt

from meltsdynamic import MELTSdynamic
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

'''
For this calculation, which involves H2O-CO2 fluid in a mafic system, we will use rhyoliteMELTS 1.2.0
     1. rhyolite-MELTS 1.0.2
     2. pMELTS
     3. rhyolite-MELTS 1.1.0
     4. rhyolite-MELTS 1.2.0
'''

ptpath = MELTSdynamic(4)

'''
Initial Composition: SiO2 51.27
Initial Composition: TiO2 0.96
Initial Composition: Al2O3 14.80
Initial Composition: Fe2O3 1.1
Initial Composition: Cr2O3 0.042
Initial Composition: FeO 8.70
Initial Composition: MnO 0.176
Initial Composition: MgO 8.45
Initial Composition: NiO 0.011
Initial Composition: CaO 12.47
Initial Composition: Na2O 2.25
Initial Composition: K2O 0.047
Initial Composition: P2O5 0.065
Initial Composition: H2O 0.25
Initial Composition: CO2 0.005
'''

bulk = [51.27, 0.96, 14.80, 1.1, 0.042, 8.7, 0.176, 8.45, 0.011,  0.0, 12.47, 2.25, 0.047, 0.065, 0.25, 0.005, 0.0, 0.0, 0.0]

# There are several ways to do set composition
ptpath.engine.setBulkComposition(bulk)


# Here we passed the whole composition vector but it is possible to set individual oxides separately
# ptpath.engine.set('bulkComposition', 'H2O', 0.5);

pressure = 250.0
temperature = 1191.0

# Note that as of Nov 18th, 2020 temperature is in Celcius
ptpath.engine.pressure = pressure
ptpath.engine.temperature = temperature

# There are other way for setting system properties.

# A 1-D array of strings that look just like the .melts file lines (fo2 Offset also supported)
ptpath.engine.setSystemProperties(["Log fO2 Path: -1FMQ", "Mode: Fractionate Solids"])

'''
Execute the path (equivalent to menu option 4)

 Select 1 to get output after equilibration and before fractionation, 2 for output after fractionation
 (either way, bulk composition will be updated after fractionation)
'''
ptpath.engine.calcEquilibriumState(1, 1)

while ptpath.engine.temperature >= 923:

    ptpath = ptpath.addNodeAfter()
    ptpath.engine.temperature = ptpath.engine.temperature - 1

    # Run mode is saved from previous run
    # Select 1 to get output after equilibration and before fractionation
    ptpath.engine.calcEquilibriumState(1, 1)
    print(ptpath.engine.status.message)

    # This means that we do not skip any failures
    if ptpath.engine.status.failed:
        break

    # Can populate X, mu0, mu for all phases or select key phases; can also calculate afterwards (see tutorial.py).
    if "plagioclase1" in ptpath.engine.solidNames:
        ptpath.engine.calcEndMemberProperties("plagioclase1", ptpath.engine.dispComposition["plagioclase1"])

'''
Examine output

By default alphaMELTS for MATLAB automatically its own _tbl.txt files that can be
processed with run-alphamelts.command -x to generate Phase_mass_tbl.txt etc.

Note that ptpath.engine.meltsIndex gives you the line number of the output in the tbl or _tbl.txt
files, wherease ptpath.nodeIndex gives you the location within the MELTSdynamic linked list
(these won't be the same if the library gets reloaded and calculations continue)

Liquidus phase(s)

You can also explore the liquidus phase by looking at the first few iterations / nodes
'''
print(ptpath.First.engine.solidNames)

# These two are almost equivalent *********************
print(ptpath.First.Next.engine.solidNames)
print(ptpath.getNodeProperty(1, 'solidNames'))

temp = ptpath.getListProperty('temperature')

Xan = ptpath.getListProperty('X', "plagioclase1", "CaAl2Si2O8")

plt.scatter(temp, Xan, marker='*', color='black')
plt.show()

cpx = ptpath.getListProperty('mass', "clinopyroxene1");

plt.scatter(temp, Xan, marker='*', color='black')
plt.show()

# Other phases

'''
indices = list(range(ptpath.Last.nodeIndex+1))
sp = ptpath.getListProperty('mass', "spinel1")
indices = np.isfinite(sp[0, :].astype(np.double))
'''

'''
print(ptpath.getNodeProperty(indices(0), 'X', "spinel1"))
print(ptpath.getNodeProperty(indices(), 'X', "spinel1"))
'''

'''
Xfo = ptpath.getListProperty('X', "olivine1", "Mg2SiO4")
plt.scatter(Xan, Xfo, marker='*', color='black')
plt.show()
'''

# Liquid

# We need to do Prev because the last point failed
print(ptpath.endMemberFormulas["bulk"])
print(ptpath.Prev.engine.dispComposition["liquid1"])


"""
F = ptpath.getListProperty('mass', "liquid1") / ptpath.First.engine.mass("bulk");
figure(4)
plot(temp, F)

feot = ptpath.getListProperty('dispComposition', 'liquid1', 'FeO') +...
    0.9*ptpath.getListProperty('dispComposition', 'liquid1', 'Fe2O3');
mgo = ptpath.getListProperty('dispComposition', 'liquid1', 'MgO');
figure(5)
plot(mgo, feot)

"""
