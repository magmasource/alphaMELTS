# %% [markdown]
# The 'tutorial' example, see:
#
# https://magmasource.caltech.edu/forum/index.php?board=12.0
#
# https://magmasource.caltech.edu/forum/index.php/board=31.0
#
# Converted to ipynb using https://jupytext.readthedocs.io/en/latest/install.html

# %%
# Can use to connect an interactive console to notebook
# %connect_info

# %%
from meltsdynamic import MELTSdynamic

# %%
try:
    import numpy as np
    import pandas as pd
except:
    from tinynumpy import tinynumpy as tnp
    np = None
    pd = None

# %%
try:
    import matplotlib.pyplot as plt
    # %matplotlib inline
except:
    plt = None

# %%
liquidus = MELTSdynamic(1)

# %%
pressure = 500.0
temperature = 1200.0

# %%
bulk = [48.68, 1.01, 17.64, 0.89, 0.03, 7.59, 0.0, 9.10, 0.0, 0.0, 12.45, 2.65, 0.03, 0.08, 0.2, 0.0, 0.0, 0.0, 0.0]

# %%
liquidus.engine.setBulkComposition(bulk)
liquidus.engine.pressure = pressure

# %%
# Note that as of Nov 18th, 2020 temperature is in Celcius
liquidus.engine.temperature = temperature

# %%
# %%capture output --no-stdout --no-display
# This should capture all output but a little goes to stdout, which is a bug that will be fixed soon
liquidus.engine.findLiquidus()
temperature = liquidus.engine.temperature

# %%
print(f'temperature = {temperature:.2f}')

# %%
liquidusPhase = liquidus.engine.solidNames
print(liquidusPhase)

# %%
# Note that as of July 4th, 2022 .calcSaturationState() outputs wt% oxides
liquidus.engine.calcEndMemberProperties(liquidusPhase, liquidus.engine.dispComposition[liquidusPhase[0]])
Xan = liquidus.engine.getProperty('X', liquidusPhase, "CaAl2Si2O8")
print(Xan)

# %%
liquidus.engine.setSystemProperties("Mode", "Fractionate Solids")

# %%
# .addNodeAfter() uses .copy(), which only copies "INPUT" or "INPUT AND OUTPUT"
# .copyAndKeepOutput() copies "OUTPUT" also. We can verify that "OUTPUT" is reset
# in ptpath when .calcEquilibriumState() (or .findLiquidus()) is called.
ptpath = liquidus.copyAndKeepOutput()

# %%
# %%capture output --no-stdout --no-display
# Run mode is 1 (isobaric, isothermal) if not defined
# Select output mode 0, to equilibrate but not do any fractionations
ptpath.engine.calcEquilibriumState(1, 0)
print(ptpath.engine.status.message)
print(ptpath.engine.solidNames)

# %%
# %%capture output --no-display
while ptpath.engine.temperature >= 1000:
    ptpath = ptpath.addNodeAfter()
    ptpath.engine.temperature = ptpath.engine.temperature - 3

    # Run mode is saved from previous run
    # Select 1 to get output after equilibration and before fractionation
    ptpath.engine.calcEquilibriumState(1, 1)
    print(ptpath.engine.status.message)

    # Can populate X, mu0 ,mu as go along (should test that the phase exists)
    # This is a Supplemental Calculator-type method so "OUTPUT" is not reset
    if "plagioclase1" in ptpath.engine.solidNames:
        ptpath.engine.calcEndMemberProperties("plagioclase1", ptpath.engine.dispComposition["plagioclase1"])

# %%
# Get the Phase_mass_tbl.txt and Phase_vol_tbl.txt files if standalone alphaMELTS 2 is installed
# !run-alphamelts.pl -x

# %%
temp = ptpath.getListProperty('temperature')
mgo = ptpath.getListProperty('dispComposition', 'liquid1', 'MgO')
mass = ptpath.getListProperty('mass', 'bulk')

# %%
# A silly example to illustrate sorting the list (in this case by temperature ascending)
# There are multiple ways to do this.
if np is not None:
    indices = np.arange(ptpath.Last.nodeIndex+1)
else:
    indices = list(range(ptpath.Last.nodeIndex+1))

# %%
if plt is not None:

    plt.scatter(indices, mgo, marker='*', color='black')

# %%
if np is not None:
    temp = np.array(temp)
    order = indices[temp.argsort()]
else:
    order = [x for _, x in sorted(zip(temp, indices), key=lambda pair: pair[0])]

# %%
ptpath.sortListBy(order)
mgo = ptpath.getListProperty('dispComposition', 'liquid1', 'MgO')

# %%
if plt is not None:
    plt.scatter(indices, mgo, marker='*', color='blue')
    plt.show()

# %%
# flip it back again
ptpath.sortListBy(order)
mgo = ptpath.getListProperty('dispComposition', 'liquid1', 'MgO')

# %%
Xan = ptpath.getListProperty('X', "plagioclase1", "CaAl2Si2O8")

# %%
# If we didn't calculate Xan as we went along and wanted to avoid repeating the equilibration loop.
plag = ptpath.getListProperty('dispComposition', 'plagioclase1')

# %%
# Need to give the feldspars from different runs different names that start with feldspar
phaseList = ['plagioclase_'] * (ptpath.Last.nodeIndex+1)
for i in range(ptpath.Last.nodeIndex+1):
    phaseList[i] += str(i+1)

# %%
# A pandas example
if pd is not None:
    table = pd.DataFrame(data=np.asarray(plag), index=ptpath.endMemberFormulas["bulk"], columns=phaseList)
    print(table)

# %%
# This works because Xan is independent of pressure and temperature but we need to specify them
# start a new list (as ptpath already has a value for Xan)
melts = MELTSdynamic(1)
melts.engine.temperature = ptpath.engine.temperature
melts.engine.pressure = ptpath.engine.pressure

# %% [markdown]
# There is no feldspar in the first ptpath phase assemblage (composition will be all zeroes; check SiO2)
# there are other ways to do this (see meltsplotter.py); .tolist() on phaseList is optional.
# indices = np.isfinite(plag[0, :].astype(np.double))
# phaseList = (np.array(phaseList)[indices])
# melts.engine.calcEndMemberProperties(phaseList, plag[:, indices])

# %%
# There is feldspar in the first ptpath now, as alphaMELTS for Python always uses findWetLiquidus
melts.engine.calcEndMemberProperties(phaseList, plag)

# %%
#To get properties from a single node, use getNodeProperty with None for the index.
#Or equivalently use the engine's getProperty method instead i.e:
#Xan = melts.engine.getProperty('X', phaseList[indices], "CaAl2Si2O8")

# %%
Xan = melts.getNodeProperty(None, 'X', phaseList, "CaAl2Si2O8")

# %%
# Ryan Clark's MELTSplotter, using matplotlib
if plt is not None:
    try:
        from meltsplotter import MELTSplotter
        mplot1 = MELTSplotter()
        mplot1 = mplot1.plotTAS(ptpath)
        mplot2 = MELTSplotter()
        mplot2 = mplot2.harkerPlot(ptpath, 'SiO2', ['CaO', 'MgO', 'K2O', 'TiO2', 'Na2O', 'FeO', 'MnO', 'Al2O3'], False)
        mplot3 = MELTSplotter()
        mplot3 = mplot3.harkerPlot(ptpath, 'MgO', ['SiO2', 'CaO', 'K2O', 'TiO2', 'Na2O', 'FeO', 'MnO', 'Al2O3'], False)
    except:
        print("Could not find MELTSplotter!")

# %%
temp = ptpath.getListProperty('temperature')
sio2 = ptpath.getListProperty('dispComposition', 'liquid1', 'SiO2')
feo = ptpath.getListProperty('dispComposition', 'liquid1', 'FeO')
#mgo = ptpath.getListProperty('dispComposition', 'liquid1', 'MgO')
cao = ptpath.getListProperty('dispComposition', 'liquid1', 'CaO')
na2o = ptpath.getListProperty('dispComposition', 'liquid1', 'Na2O')

# %%
try:
    import pygmt
    #pygmt.show_versions()
    fig = pygmt.Figure()
    fig.basemap(region=[1000, 1220, 0, 100], projection="X-4i/3i", frame=['WeSn+tMelt composition', 'x+lT @.C system', 'y+Lwt%'])
    fig.plot(x=temp, y=sio2, pen="thin,lightblue", style="c0.1")
    fig.plot(x=temp, y=feo, pen="thin,orange", style="c0.1")
    fig.plot(x=temp, y=mgo, pen="thin,magenta", style="c0.1")
    fig.plot(x=temp, y=cao, pen="thin,gray", style="c0.1")
    fig.plot(x=temp, y=na2o, pen="thin,pink", style="c0.1")
    try:
        fig.legend(position='jTL', box=True, spec='legend.txt')
    except:
        print("Could not find legend.txt")
    fig.show()
except:
    print("Could not find pygmt")

# %%
