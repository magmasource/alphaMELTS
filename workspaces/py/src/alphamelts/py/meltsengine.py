from ctypes import c_int, c_double, c_char_p, create_string_buffer, pointer #, byref
try:
    import numpy as np
except:
    import math
    np = None
from tinynumpy import tinynumpy as tnp
from copy import copy, deepcopy
import re

# https://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook/39662359#39662359
try:
    shell = get_ipython().__class__.__name__
    if shell == 'ZMQInteractiveShell':
        from wurlitzer import sys_pipes as pipes # Jupyter notebook or qtconsole
    elif shell == 'TerminalInteractiveShell':
        pipes = None  # Terminal running IPython
    else:
        pipes = None  # Other type (?)
except:
    pipes = None     # Probably standard Python interpreter (or wurlitzer not available)

from meltsstatus import MELTSstatus

class MELTSengine(object):
    ''' Thermodynamic 'engine' where the real work happens.
        Calculation mode is set in parent MELTSdyanmic list: 1. rhyolite-MELTS 1.0.2; 2. pMELTS; 3. rhyolite-MELTS 1.1.0; 4. rhyolite-MELTS 1.2.0.
        For equilibration, run mode need to be one of: 1. isobaric, isothermal; 2. isenthalpic; 3. isentropic; 4. isochoric.
        Output is post-equilibration / pre-fractionation (1) by default, but calculation can also be batch, no fractionation (0), or output can be after fractionation (2).
        phaseNames format is: [phase name in MELTS system][optional integer, like equilibration output][optional underscore, '_'][optional user-defined name, comprising letters and digits] '''

    # MELTSstatus class used to track status of the C library of MELTS functions.
    status = None
    # MELTS calculation mode inherited from MELTSdynamic.
    calculationMode = None

    # User-defined index for this node of the parent MELTSdynamic list.
    nodeId = None
    # User-defined name and/or description for this node of the parent MELTSdynamic list.
    nodeName = None
    # Tracks the 'index' in any text output (may not work for complex calculations where library is reloaded or MELTS dynamic list is adjusted).
    meltsIndex = None
    # Controls whether any constraints are applied in equilibration calculations: 1. none (default); 2. isenthalpic; 3. isentropic; 4. isochoric.
    runMode = None
    # INPUT and OUTPUT: pressure (bars)
    pressure = None
    # INPUT and OUTPUT: temperature (Celsius)
    temperature = None
    # INPUT and OUTPUT: Reference quantity for isenthalpic (J), isentropic (J/k), or isochoric (cc) calculations.
    reference = None
    # OUTPUT: oxygen fugacity (log 10 units) of the system.
    logfO2 = None
    # List of system properties that resemble lines in .melts files.
    systemProperties = None
    #outputFlag # put in systemProperties?

    # INPUT and OUTPUT: composition in grams (or loosely wt%) for regular MELTS calculations; updated after any fractionations regardless of output mode.
    bulkComposition = None

    # OUTPUT: list of liquid phase name(s) from regular MELTS calculations (multiple instanced count from 0); may be INPUT as phaseNames to 'Supplemental calculator'.
    liquidNames = None
    # OUTPUT: list of solid phase names from regular MELTS calculations (multiple instanced count from 0); may be INPUT as phaseNames to 'Supplemental calculator'.
    solidNames = None

    # INPUT: list of phase names for 'Supplemental calculator' calculations (currently also output for calcSaturationState).
    phaseNames = None
    # INPUT: composition in grams (or loosely wt%) for 'Supplemental calculator' calculations.
    phaseComposition = None
    # INPUT: composition in end-member mol frac (or oxide mol frac) for 'Supplemental calculator' calculations.
    molarComposition = None

    # OUTPUT: Dictionary of Gibbs free energy (J) for solid and liquid phases and/or 'bulk'.
    g = None
    # OUTPUT: Dictionary of enthalpy (J) for solid and liquid phases and/or 'bulk'.
    h = None
    # OUTPUT: Dictionary of entropy (J/K) for solid and liquid phases and/or 'bulk'.
    s = None
    # OUTPUT: Dictionary of volume (cc) for solid and liquid phases and/or 'bulk'.
    v = None
    # OUTPUT: Dictionary of heat capacity (J/K) for solid and liquid phases and/or 'bulk'.
    cp = None
    # OUTPUT: Dictionary of partial derivative of heat capacity wrt temperature (J/K/K) for solid and liquid phases and/or 'bulk'.
    dcpdt = None
    # OUTPUT: Dictionary of partial derivative of volume wrt temperature (cc/K) for solid and liquid phases and/or 'bulk'.
    dvdt = None
    # OUTPUT: Dictionary of partial derivative of volume wrt pressure (cc/bar) for solid and liquid phases and/or 'bulk'.
    dvdp = None
    # OUTPUT: Dictionary of 2nd partial derivative of volume wrt temperature (cc/K/K) for solid and liquid phases and/or 'bulk'.
    d2vdt2 = None
    # OUTPUT: Dictionary of 2nd partial derivative of volume wrt temperature, pressure (cc/K/bar) for solid and liquid phases and/or 'bulk'.
    d2vdtdp = None
    # OUTPUT: Dictionary of 2nd partial derivative of volume wrt pressure (cc/bar/bar) for solid and liquid phases and/or 'bulk'.
    d2vdp2 = None
    # OUTPUT: Dictionary of molecular weights of phases from regular MELTS calculations; use to convert to/from moles.
    molwt = None
    # OUTPUT: Dictionary of density (g/cc) for solid and liquid phases and/or 'bulk'.
    rho = None
    # OUTPUT: Dictionary of mass (grams) for solid and liquid phases and/or 'bulk'.
    mass = None
    # OUTPUT: Dictionary of log 10 viscosity (Poise) for one or more liquid phases, and system viscosity (not including previously fractionated material).
    viscosity = None

    # OUTPUT: Dictionary of compositions (wt%) for solid and liquid phases and/or 'bulk'; may be INPUT for 'Supplemental Calculator' calculations.
    dispComposition = None

    # OUTPUT: Dictionary of affinities (J/mol) for solid and liquid phases.
    affinity = None
    # OUTPUT: Dictionary (mol frac) for solid and liquid phases and/or 'bulk'; may be INPUT for molar properties.
    X = None

    # OUTPUT: Dictionary for end-member activities relative to fixed structural and/or ordering state (for 'true' activities, using mu and mu0 or 'activity' is usually preferable).
    activity0 = None
    # OUTPUT: Dictionary for end-member activities relative to pure phase structural and/or ordering state at the given pressure and temperature (calculated from mu and mu0)
    activity = None
    # OUTPUT: Dictionary for end-member (or oxide) pure-phase chemical potentials (J/mol) according to Berman/MELTS thermodynamic database, for solid and liquid phases and/or 'bulk'.
    mu0 = None
    # OUTPUT: Dictionary for end-member (or oxide) chemical potentials (J/mol) for solid and liquid phases and/or 'bulk'.
    mu = None

    def __init__(self, cMode, *args, **kwargs):
        ''' Called by MELTSdynamic to set up the thermodynamic 'engine' where all the real work happens. '''
        name = None
        if 'nodeName' in kwargs:
            name = kwargs['nodeName']
        elif len(args):
            name = args[::]

        if isinstance(cMode, MELTSstatus):
            self.status = cMode
            self.calculationMode = self.status.calculationMode
        elif cMode is not None:
            self.calculationMode = cMode
            self.status = MELTSstatus(cMode)

        if name is not None:
            self.nodeName = name

        self.runMode = 1
        self.reference = 0.0

    def setBulkComposition(self, oxide='', val=None):
        ''' Set the bulk compsoitions for a particular oxide (in grams), or the entire compositional array. '''
        if self.calculationMode != self.status.getCalculationMode():
            if self.status.setCalculationMode(self.calculationMode):
                self.setSystemProperties(self.systemProperties)

        if self.bulkComposition is None:
            self.bulkComposition = tnp.array([0.0] * len(self.status.endmembers["bulk"]))

        try:
            if not oxide.isnumeric() and val is not None:
                iox = self.status.endmembers["bulk"].index(oxide.lower())
                if iox >= 0 and iox < len(self.bulkComposition):
                    self.bulkComposition[iox] = val
                else:
                    print (f'Oxide {oxide:s} not found, Ignoring value\n.')
        except AttributeError:
            assert (len(oxide) == len(self.bulkComposition)), "Array sizes do not match"
            self.bulkComposition[:] = oxide[:]

    def setPhaseComposition(self, oxide='', val=None):
        ''' For a given phase set either a particular oxide (in grams), or the entire compositional array. '''
        if self.calculationMode != self.status.getCalculationMode():
            if self.status.setCalculationMode(self.calculationMode):
                self.setSystemProperties(self.systemProperties)

        if self.phaseComposition is None:
            self.phaseComposition = tnp.array([0.0] * len(self.status.endmembers["bulk"]))
        try:
            if not oxide.isnumeric() and val is not None:
                iox = self.status.endmembers["bulk"].index(oxide.lower())
                if iox >= 0 and iox < len(self.phaseComposition):
                    self.phaseComposition[iox] = val
                else:
                    print (f'Oxide {oxide:s} not found, Ignoring value\n.')
        except AttributeError:
            assert (len(oxide) == len(self.phaseComposition)), "Array sizes do not match"
            self.phaseComposition[:] = oxide[:]

    def setInitialComposition(self, oxide='', val=None):
        ''' Set both bulk composition and phase composition to the passed value(s). '''
        self.setBulkComposition(oxide, val)
        self.setPhaseComposition(oxide, val)

    def setMolarComposition(self, phase='', *args):
        ''' For a given phase set either a particular end-member mole frac (or oxide for 'bulk'), or the entire mole fraction array. '''

        if self.calculationMode != self.status.getCalculationMode():
            if self.status.setCalculationMode(self.calculationMode):
                self.setSystemProperties(self.systemProperties)

        if self.molarComposition is None:
            self.molarComposition = tnp.array([0.0] * len(self.status.endmembers["bulk"]))

        if isinstance(phase, list):
            phase = phase[0]
        try:
            if not phase.isnumeric() and len(args) > 1:
                iend = self.status.endmembers[phase].index(args[0].lower())
                if iend >= 0 and iend < len(self.phaseComposition):
                    self.molarComposition[iend] = args[1]
                else:
                    print (f'Endmember {phase:s} not found, Ignoring value\n.')
        except AttributeError:
            assert (len(phase) <= len(self.molarComposition)), "Array sizes do not match"
            length = len(phase)
            self.molarComposition[::length] = phase[:]
            self.molarComposition[length::] = 0.0

    def getProperty(self, propertyName, *args):
        ''' Get some property (for system, or one or more phases) from engine.
            Returns a scalar, vector or 2-D matrix. '''

        # Assume endmember or oxide present in all
        phaseList = None
        phaseName = None

        if len(args):
            if propertyName == 'bulkComposition' or propertyName == 'phaseComposition' or propertyName == 'molarComposition':
                value = self.__dict__[propertyName]
                endMemberName = args[::]
            else:
                # MELTSmap / dictionary
                phaseList = args[0]
                if np is not None and isinstance(phaseList, np.ndarray):
                    phaseList = phaseList.tolist()
                elif not isinstance(phaseList, list):
                    phaseList = [phaseList]

                if len(args) > 1:
                    endMemberName = args[1]
                    if propertyName == 'dispComposition':
                        phaseName = "bulk"
                    else:
                        phaseName = re.sub(r'(_+.*)|(\d*)', '', phaseList[0])
                    value = tnp.array([float("NaN")] * len(self.status.endmembers[phaseName]))
                else:
                    endMemberName = None
                    if propertyName == 'dispComposition':
                        value = tnp.array([float("NaN")] * len(self.status.endmembers["bulk"]))
                    elif propertyName == 'X' or propertyName == 'activity' or propertyName == 'activity0' or propertyName == 'mu' or propertyName == 'mu0':
                        phaseName = phaseList[0]
                        phaseName = re.sub(r'(_+.*)|(\d*)', '', phaseName)
                        assert all(map(lambda phase: re.match(phaseName, phase), phaseList)), \
                            "The getProperty method cannot be called for multiple phases for that property, as the number of endmembers may vary."
                        value = tnp.array([float("NaN")] * len(self.status.endmembers[phaseName]))
                    else:
                        value = float("NaN")

                try:
                    if self.__dict__[propertyName] is not None:
                        if len(phaseList) > 1:
                            try:
                                if phaseList[0] in self.__dict__[propertyName]:
                                    vlen = len((self.__dict__[propertyName])[phaseList[0]])
                                else:
                                    vlen = len(value)
                                value2 = tnp.empty((len(phaseList), vlen))

                                for i in range(len(phaseList)):
                                    label = phaseList[i]
                                    if label in self.__dict__[propertyName]:
                                        value2[i, :] = (self.__dict__[propertyName])[label][:]
                                    else:
                                        value2[i, :] = value[:]
                                value2 = value2.transpose()
                            except:
                                # value is scalar or nan
                                value2 = tnp.empty((len(phaseList),))

                                for i in range(len(phaseList)):
                                    label = phaseList[i]
                                    if label in self.__dict__[propertyName]:
                                        value2[i] = (self.__dict__[propertyName])[label]
                                    else:
                                        value2[i] = value

                        elif phaseList[0] in self.__dict__[propertyName]:
                            if isinstance(value, tnp.ndarray):
                                value[:] = self.__dict__[propertyName][phaseList[0]]
                            else:
                                value = self.__dict__[propertyName][phaseList[0]]

                        try:
                            endMembers = self.status.endmembers[phaseName]
                            indices = endMembers.index(endMemberName.lower())
                        except (KeyError, AttributeError, ValueError):
                            # phaseName is None because it's not the kind of map you can select end members for
                            # or endMemberName either wasn't passed or is not a string
                            indices = None

                except KeyError:
                    value = None
                    indices = None

            if len(phaseList) > 1 and indices is not None:
                value = value2[indices, :]
            elif len(phaseList) > 1:
                value = value2
            elif indices is not None:
                value = value[indices]

        else:
            value = self.__dict__[propertyName]

        try:
            if (np is None and not math.isnan(self.status.nullValue)) or (np is not None and not np.isnan(self.status.nullValue)):
                if isinstance(value, tnp.ndarray):
                    for i in range(len(value.ravel())):
                        if (np is None and math.isnan(value.ravel()[i])) or np.isnan(value.ravel()[i]):
                            value.ravel()[i] = self.status.nullValue
                else:
                    if (np is None and math.isnan(value)) or np.isnan(value):
                        value = self.status.nullValue
        except TypeError:
            print("Cannot replace NaN with non-numeric value.")

        return value

    def setSystemProperties(self, *args):
        ''' Set system properties using string arrays or pairs of strings that resemble lines in .melts files. '''
        nchars = 132 # REC is 134 and the last two characters will be used for ' \0'

        if self.calculationMode != self.status.getCalculationMode():
            if self.status.setCalculationMode(self.calculationMode):
                self.setSystemProperties(self.systemProperties)

        if isinstance(args[0], list):
            alllines = ''
            lines = args[0]
            nlines = len(lines)
            for i in range(nlines):
                alllines += lines[i].ljust(nchars)
            properties = create_string_buffer(alllines.encode())
        elif len(args)%2:
            print ("Number of values does not match number of properties?")
            return
        else:
            alllines = ''
            lines = []
            nlines = len(args)//2
            for i in range(nlines):
                alllines += f'{args[2*i]:s}: {args[2*i+1]:s}'.ljust(nchars)
                lines.append(f'{args[2*i]:s}: {args[2*i+1]:s}')
            properties = create_string_buffer(alllines.encode())

        try:
            self.status.failed = True
            failure = c_int(self.status.failed)
            nLinesInArray = c_int(nlines)
            nCharInString = c_int(nchars)
            if pipes is not None:
                with pipes():
                    self.status.libalphamelts.setMeltsSystemProperties(pointer(failure), properties, pointer(nCharInString), pointer(nLinesInArray))
            else:
                self.status.libalphamelts.setMeltsSystemProperties(pointer(failure), properties, pointer(nCharInString), pointer(nLinesInArray))
            self.status.failed = failure.value
        except OSError:
            print("OSERROR")

        if self.status.failed:
            if self.status.reload("setSystemProperties"):
                self.setSystemProperties(self.systemProperties)

        self.systemProperties = self.systemProperties or []

        for line in lines:
            if line not in self.systemProperties:
                self.systemProperties.append(line)

    # Eventually have an alphaMELTS-like find boundary ('isograd') too
    def findLiquidus(self, *args):
        ''' Find liquidus for current or passed bulk composition. Gets affinities for solid phases that are not suppressed.
        '''

        self.calcEquilibriumState(0, 0, *args)
        self.phaseNames = self.calcSaturationState()

        affinities = self.getProperty('affinity', self.phaseNames)

        minAffinity = 0.0
        for j in range (len(affinities)):
            if affinities[j] < minAffinity:
                minAffinity =affinities[j]
                liquidusPhase = self.phaseNames[j]

        self.solidNames.append(liquidusPhase+"1")
        self.mass[liquidusPhase+"1"] = 0.0
        self.dispComposition[liquidusPhase+"1"] = self.dispComposition[liquidusPhase].copy()

    def calcEquilibriumState(self, runMode=None, outputFlag=1, name='', *args):
        ''' Single equilibration calculation for current or passed bulk composition with or without constraints
         (equivalent to alphaMELTS menu option 3 or 4). By default output is post-equilibration / pre-fractionation. '''

        if self.calculationMode != self.status.getCalculationMode():
            if self.status.setCalculationMode(self.calculationMode):
                self.setSystemProperties(self.systemProperties)

        if runMode is not None and runMode:
            # reset runMode if this is not a findLiquidus call
            self.runMode = (runMode or self.runMode or 1)
        elif runMode is None:
            runMode = self.runMode

        if len(name):
            self.nodeName = name

        if len(args):
            self.setBulkComposition(args)

        assert runMode in range(0, 5), "Unexpected run mode?"
        assert self.temperature, "System temperature not set?"
        self.pressure = (self.pressure or 1)
        assert self.pressure, "System pressure not set?"
        assert any(self.bulkComposition), "System bulk composition not set?"
        self.reference = self.reference or 0.0

        nCharInName = c_int(20)
        numberPhases = c_int(20)
        indices = [0] * 20
        phaseIndices = (c_int * len(indices)) (*indices)
        phases = create_string_buffer(20*20)
        nCharInString = c_int(100)
        errorString = create_string_buffer(100)
        runMode = c_int(runMode)
        outputFlag = c_int(outputFlag)

        nox = len(self.status.endmembers['bulk'])
        nf = len(self.status.fields) - 3
        nc = nf + nox + 3
        propArr = [0.0] * (nc*20) # numberPhases
        properties = (c_double * len(propArr)) (*propArr)

        bulkComposition = self.bulkComposition
        bulk = (c_double * len(bulkComposition)) (*bulkComposition)
        temperature = c_double(self.temperature+273.15)
        pressure = c_double(self.pressure)
        reference = c_double(self.reference)

        try:
            self.status.failed = True
            failure = c_int(self.status.failed)

            if pipes is not None:
                with pipes():
                    self.status.libalphamelts.driveMeltsProcess(pointer(failure), pointer(runMode), pointer(pressure), bulk, pointer(reference), pointer(temperature), \
                        phases, pointer(nCharInName), pointer(numberPhases), pointer(outputFlag), errorString, pointer(nCharInString), properties, phaseIndices)
            else:
                self.status.libalphamelts.driveMeltsProcess(pointer(failure), pointer(runMode), pointer(pressure), bulk, pointer(reference), pointer(temperature), \
                    phases, pointer(nCharInName), pointer(numberPhases), pointer(outputFlag), errorString, pointer(nCharInString), properties, phaseIndices)
            self.status.failed = failure.value
            self.status.message = errorString.value.decode()
        except OSError:
            print("OSERROR")
            #self.status.message = errorString.value.decode()

        # Engine output is reset
        self.g = {}
        self.h = {}
        self.s = {}
        self.v = {}
        self.cp = {}
        self.dcpdt = {}
        self.dvdt = {}
        self.dvdp = {}
        self.d2vdt2 = {}
        self.d2vdtdp = {}
        self.d2vdp2 = {}
        self.molwt = {}
        self.rho = {}
        self.mass = {}
        self.viscosity = {}
        self.affinity = {}
        self.dispComposition = {}
        self.X = {}
        self.activity = {}
        self.mu0 = {}
        self.mu = {}

        if self.status.failed:
            if self.status.reload("calcEquilibriumState"):
                self.setSystemProperties(self.systemProperties)

        else:

            if runMode:
                # Find Liquidus does not generate file output, but the rest do
                self.status.meltsIndex += 1
                self.meltsIndex = self.status.meltsIndex

            # Get and process the phaseIndices array
            phaseIndex = [0] * numberPhases.value
            for i in range(numberPhases.value):
                phaseIndex[i] = phaseIndices[i]%10 + 1

            # These values are post-fractionation
            self.pressure = pressure.value
            self.reference = reference.value # could actually be entropy, enthalpy or volume
            self.temperature = temperature.value-273.15
            for i in range(nox):
                self.bulkComposition[i] = bulk[i]

            phases = phases.value.decode()
            phases = phases.split()

            self.dispComposition["bulk"] = tnp.array(properties[nf:nf+nox])
            total = self.dispComposition["bulk"].sum()
            for j in range(nox):
                # element-wise operations don't seem to be available in pip-packaged tinynumpy
                self.dispComposition["bulk"][j] *= 100.0/total
            for j in range(nf):
                self.__dict__[self.status.fields[j]]["bulk"] = properties[j]

            phaseProps = tnp.array(properties[nox+nf:nox+nf+3])
            self.__dict__["viscosity"]["bulk"] = phaseProps[2] # borrow the total grams slot
            phaseProps[2] = phaseProps[1] * properties[3]
            for j in range(3):
                self.__dict__[self.status.fields[nf+j]]["bulk"] = phaseProps[j]

            for j in range(nf):
                self.__dict__[self.status.fields[j]]["oxygen"] = properties[nc + j]
            phaseProps = tnp.array(properties[nc + nox+nf:nc + nox+nf+3])
            self.__dict__["logfO2"] = phaseProps[0]  # borrow the mw slot
            phaseProps[0] = 31.998
            for j in range(3):
                self.__dict__[self.status.fields[nf+j]]["oxygen"] = phaseProps[j]

            self.liquidNames = []
            self.solidNames = []
            for i in range(2, numberPhases.value):

                phaseName = phases[i]+str(phaseIndex[i])

                self.dispComposition[phaseName] = tnp.array(properties[i*nc + nf:i*nc + nf+nox])
                total = self.dispComposition[phaseName].sum()
                for j in range(nox):
                    self.dispComposition[phaseName][j] *= 100.0/total
                for j in range(nf):
                    self.__dict__[self.status.fields[j]][phaseName] = properties[i*nc + j]

                phaseProps = tnp.array(properties[i*nc + nox+nf:i*nc + nox+nf+3])
                if phases[i] == "liquid":
                    # At the moment it is assumed there is only one liquid
                    self.liquidNames.append(phaseName)
                    self.viscosity[phaseName] = phaseProps[2] # borrow the total grams slot
                    phaseProps[2] = phaseProps[1] * properties[i*nc + 3]
                else:
                    self.solidNames.append(phaseName)

                for j in range(3):
                    self.__dict__[self.status.fields[nf+j]][phaseName] = phaseProps[j]

    def calcPhaseProperties(self, phaseList=None, *args):
        ''' 'Supplemental Calculator' type calculation to get thermodynamic properties of one or more phases.
            If composition(s) not passed, will use the contents of phaseComposition (grams). '''
        if self.calculationMode != self.status.getCalculationMode():
            if self.status.setCalculationMode(self.calculationMode):
                self.setSystemProperties(self.systemProperties)

        assert self.temperature, "System temperature not set?"
        self.pressure = (self.pressure or 1)
        assert self.pressure, "System pressure not set?"

        if phaseList is None:
            phaseList = ['liquid']
        elif np is not None and isinstance(phaseList, np.ndarray):
            phaseList = phaseList.tolist()
        elif not isinstance(phaseList, list):
            phaseList = [phaseList]

        nox = len(self.status.endmembers['bulk'])
        nf = len(self.status.fields) - 3
        nc = nf + nox + 3

        if len(phaseList) > 1:
            bulk = args[0] # should be (t)np array
            assert (isinstance(bulk, tnp.ndarray) or isinstance(bulk, np.ndarray)) and bulk.shape == (nox, len(phaseList)), \
                "Please supply a row of phase names, and a matrix with compositions as columns"

        self.phaseNames = self.phaseNames or []
        self.phaseNames.extend(phaseList)

        for i in range(len(phaseList)):

            if len(args):
                if len(phaseList) > 1:
                    # if all nan then continue - not yet implemented
                    self.setPhaseComposition(bulk[:,i])
                else:
                    self.setPhaseComposition(*args)

            assert any(self.phaseComposition), "Phase composition not set?"

            propArr = [0.0] * nc
            properties = (c_double * len(propArr)) (*propArr)

            phaseComposition = self.phaseComposition
            bulk = (c_double * len(phaseComposition)) (*phaseComposition)
            temperature = c_double(self.temperature+273.15)
            pressure = c_double(self.pressure)

            phaseName = phaseList[i]
            phase = create_string_buffer(re.sub(r'(_+.*)|(\d*)', '', phaseName).encode())
            visc = c_double(0.0)

            self.status.failed = True
            failure = c_int(self.status.failed)
            try:
                if pipes is not None:
                    with pipes():
                        if phaseName.startswith('bulk') or phaseName.startswith('liquid'):
                            self.status.libalphamelts.getMeltsPhaseProperties(pointer(failure), 'liquid', pointer(temperature), pointer(pressure), bulk, properties)
                            if failure is False:
                                self.status.libalphamelts.getMeltsViscosity(pointer(failure), 'Shaw', pointer(temperature), bulk, pointer(visc))
                        else:
                            self.status.libalphamelts.getMeltsPhaseProperties(pointer(failure), phase, pointer(temperature), pointer(pressure), bulk, properties)
                else:
                    if phaseName.startswith('bulk') or phaseName.startswith('liquid'):
                        self.status.libalphamelts.getMeltsPhaseProperties(pointer(failure), 'liquid', pointer(temperature), pointer(pressure), bulk, properties)
                        if failure is False:
                            self.status.libalphamelts.getMeltsViscosity(pointer(failure), 'Shaw', pointer(temperature), bulk, pointer(visc))
                    else:
                        self.status.libalphamelts.getMeltsPhaseProperties(pointer(failure), phase, pointer(temperature), pointer(pressure), bulk, properties)
                self.status.failed = failure.value
            except OSError:
                print("OSERROR")

            if self.status.failed:
                if self.status.reload("calcPhaseProperties"):
                    self.setSystemProperties(self.systemProperties)
            else:

                self.dispComposition = self.dispComposition or {}

                self.dispComposition[phaseName] = tnp.array(properties[nf:nf+nox])
                total = self.dispComposition[phaseName].sum()
                for j in range(nox):
                    self.dispComposition[phaseName][j] *= 100.0/total
                for j in range(nf):
                    if not self.status.fields[j] in self.__dict__.keys():
                        self.__dict__[self.status.fields[j]] = {}
                    self.__dict__[self.status.fields[j]][phaseName] = properties[j]

                for j in range(3):
                    if not self.status.fields[nf+j] in self.__dict__.keys():
                        self.__dict__[self.status.fields[nf+j]] = {}
                    self.__dict__[self.status.fields[nf+j]][phaseName] = properties[nf+nox+j]
                if phaseName.startswith('bulk') or phaseName.startswith('liquid'):
                    self.viscosity = self.viscosity or {}
                    self.viscosity[phaseName] = visc.value

    def calcMolarProperties(self, phaseList=None, *args):
        ''' 'Supplemental Calculator' type calculation to get thermodynamic properties of one or more phases, or of oxygen.
            If composition(s) not passed, will use the contents of molarComposition (mol frac). '''

        if self.calculationMode != self.status.getCalculationMode():
            if self.status.setCalculationMode(self.calculationMode):
                self.setSystemProperties(self.systemProperties)

        assert self.temperature, "System temperature not set?"
        self.pressure = (self.pressure or 1)
        assert self.pressure, "System pressure not set?"

        if phaseList is None:
            phaseList = ['liquid']
        elif np is not None and isinstance(phaseList, np.ndarray):
            phaseList = phaseList.tolist()
        elif not isinstance(phaseList, list):
            phaseList = [phaseList]

        nox = len(self.status.endmembers['bulk'])
        nf = len(self.status.fields) - 3
        nc = nf + nox + 3

        if len(phaseList) > 1:
            bulk = args[0] # should be (t)np array
            assert (isinstance(bulk, tnp.ndarray) or isinstance(bulk, np.ndarray)) and bulk.shape[0] <= nox and bulk.shape[1] == len(phaseList), \
                "Please supply a row of phase names, and a matrix with compositions as columns"

        self.phaseNames = self.phaseNames or []
        self.phaseNames.extend(phaseList)

        for i in range(len(phaseList)):

            if len(args):
                if len(phaseList) > 1:
                    bulk = args[0]
                    self.setMolarComposition(bulk[:,i])
                else:
                    self.setMolarComposition(phaseList, *args)

            assert any(self.molarComposition), "Phase composition not set?"

            propArr = [0.0] * nc
            properties = (c_double * len(propArr)) (*propArr)

            phaseComposition = self.phaseComposition
            bulk = (c_double * len(phaseComposition)) (*phaseComposition)
            temperature = c_double(self.temperature+273.15)
            pressure = c_double(self.pressure)

            phaseName = phaseList[i]
            phase = create_string_buffer(re.sub(r'(_+.*)|(\d*)', '', phaseName).encode())
            visc = c_double(0.0)

            self.status.failed = True
            failure = c_int(self.status.failed)
            try:
                if pipes is not None:
                    with pipes():
                        if phaseName.startswith('bulk') or phaseName.startswith('liquid'):
                            self.status.libalphamelts.getMeltsMolarProperties(pointer(failure), 'liquid', pointer(temperature), pointer(pressure), bulk, properties)
                            if failure is False:
                                self.status.libalphamelts.getMeltsViscosity(pointer(failure), 'Shaw', pointer(temperature), bulk, pointer(visc))
                        else:
                            self.status.libalphamelts.getMeltsMolarProperties(pointer(failure), phase, pointer(temperature), pointer(pressure), bulk, properties)
                else:
                    if phaseName.startswith('bulk') or phaseName.startswith('liquid'):
                        self.status.libalphamelts.getMeltsMolarProperties(pointer(failure), 'liquid', pointer(temperature), pointer(pressure), bulk, properties)
                        if failure is False:
                            self.status.libalphamelts.getMeltsViscosity(pointer(failure), 'Shaw', pointer(temperature), bulk, pointer(visc))
                    else:
                        self.status.libalphamelts.getMeltsMolarProperties(pointer(failure), phase, pointer(temperature), pointer(pressure), bulk, properties)
                self.status.failed = failure.value
            except OSError:
                print("OSERROR")

            if self.status.failed:
                if self.status.reload("calcMolarProperties"):
                    self.setSystemProperties(self.systemProperties)
            else:

                self.dispComposition = self.dispComposition or {}

                self.dispComposition[phaseName] = tnp.array(properties[nf:nf+nox])
                total = self.dispComposition[phaseName].sum()
                for j in range(nox):
                    self.dispComposition[phaseName][j] *= 100.0/total
                for j in range(nf):
                    if not self.status.fields[j] in self.__dict__.keys():
                        self.__dict__[self.status.fields[j]] = {}
                    self.__dict__[self.status.fields[j]][phaseName] = properties[j]

                for j in range(3):
                    if not self.status.fields[nf+j] in self.__dict__.keys():
                        self.__dict__[self.status.fields[nf+j]] = {}
                    self.__dict__[self.status.fields[nf+j]][phaseName] = properties[nf+nox+j]
                if phaseName.startswith('bulk') or phaseName.startswith('liquid'):
                    self.viscosity = self.viscosity or {}
                    self.viscosity[phaseName] = visc.value

    def calcViscosityFromGRD(self, phaseList=None, *args):
        ''' 'Supplemental Calculator' type calculation to get viscosity of of one or more liquids using:
            Giordano D, Russell JK, Dingwell DB (2008) Viscosity of magmatic liquids: A model. EPSL 271, 123-134.
            If composition(s) not passed, will use the contents of phaseComposition (grams).
            Call this after calcPhaseProperties or calcMolarProperties to avoid the GRD viscosity being overwritten.
            System ('bulk') viscosity is not updated by this method; it always uses Shaw model and crystallinity.'''
        if self.calculationMode != self.status.getCalculationMode():
            if self.status.setCalculationMode(self.calculationMode):
                self.setSystemProperties(self.systemProperties)

        assert self.temperature, "System temperature not set?"
        self.pressure = (self.pressure or 1)
        assert self.pressure, "System pressure not set?"

        if phaseList is None:
            phaseList = ['liquid']
        elif np is not None and isinstance(phaseList, np.ndarray):
            phaseList = phaseList.tolist()
        elif not isinstance(phaseList, list):
            phaseList = [phaseList]

        nox = len(self.status.endmembers['bulk'])
        nf = len(self.status.fields) - 3
        nc = nf + nox + 3

        if len(phaseList) > 1:
            bulk = args[0] # should be (t)np array
            assert (isinstance(bulk, tnp.ndarray) or isinstance(bulk, np.ndarray)) and bulk.shape == (nox, len(phaseList)), \
                "Please supply a row of phase names, and a matrix with compositions as columns"

        self.phaseNames = self.phaseNames or []
        self.phaseNames.extend(phaseList)

        for i in range(len(phaseList)):

            if len(args):
                if len(phaseList) > 1:
                    bulk = args[0] # should be np array
                    self.setPhaseComposition(bulk[:,i])
                else:
                    self.setPhaseComposition(*args)

            assert any(self.phaseComposition), "Phase composition not set?"

            propArr = [0.0] * nc
            properties = (c_double * len(propArr)) (*propArr)

            phaseComposition = self.phaseComposition
            bulk = (c_double * len(phaseComposition)) (*phaseComposition)
            temperature = c_double(self.temperature+273.15)
            pressure = c_double(self.pressure)

            phaseName = phaseList[i]
            phase = create_string_buffer(re.sub(r'(_+.*)|(\d*)', '', phaseName).encode())
            visc = c_double(0.0)

            assert (phaseName.startswith('bulk') or phaseName.startswith('liquid')), "Phase is not liquid in viscosity calc?"

            self.status.failed = True
            failure = c_int(self.status.failed)
            try:
                if pipes is not None:
                    with pipes():
                        self.status.libalphamelts.getMeltsViscosity(pointer(failure), 'GRD', pointer(temperature), bulk, pointer(visc))
                else:
                    self.status.libalphamelts.getMeltsViscosity(pointer(failure), 'GRD', pointer(temperature), bulk, pointer(visc))
                self.status.failed = failure.value
            except OSError:
                print("OSERROR")

            if self.status.failed:
                if self.status.reload("calcViscosityFromGRD"):
                    self.setSystemProperties(self.systemProperties)
            else:
                self.viscosity = self.viscosity or {}
                self.viscosity = tnp.log10(10.0 * visc.value) # PaS -> log10 Poise

    def calcEndMemberProperties(self, phaseList=None, *args):
        ''' 'Supplemental Calculator' type calculation to get end-member thermodynamic properties of one or more phases.
            If composition(s) not passed, will use the contents of phaseComposition (grams).
            If there is an fO2 buffer, liquid composition will be updated according Kress & Carmichael (1991).
        '''

        if self.calculationMode != self.status.getCalculationMode():
            if self.status.setCalculationMode(self.calculationMode):
                self.setSystemProperties(self.systemProperties)

        assert self.temperature, "System temperature not set?"
        self.pressure = (self.pressure or 1)
        assert self.pressure, "System pressure not set?"

        if phaseList is None:
            phaseList = ['liquid']
        elif np is not None and isinstance(phaseList, np.ndarray):
            phaseList = phaseList.tolist()
        elif not isinstance(phaseList, list):
            phaseList = [phaseList]

        nox = len(self.status.endmembers['bulk'])
        fields = ('X', 'act', 'mu0', 'mu')
        nc = len(fields)

        if len(phaseList) > 1:
            bulk = args[0] # should be (t)np array
            assert (isinstance(bulk, tnp.ndarray) or isinstance(bulk, np.ndarray)) and bulk.shape == (nox, len(phaseList)), \
                "Please supply a row of phase names, and a matrix with compositions as columns"

        self.phaseNames = self.phaseNames or []
        self.phaseNames.extend(phaseList)

        for i in range(len(phaseList)):

            if len(args):
                if len(phaseList) > 1:
                    bulk = args[0] # should be np array
                    self.setPhaseComposition(bulk[:,i])
                else:
                    self.setPhaseComposition(*args)

            assert any(self.phaseComposition), "Phase composition not set?"

            nCharInName = c_int(30)
            numberEndMembers = c_int(20)

            endmembers = create_string_buffer(30*20)

            propArr = [0.0] * (nc*20) # numberPhases
            properties = (c_double * len(propArr)) (*propArr)

            phaseComposition = self.phaseComposition
            bulk = (c_double * len(phaseComposition)) (*phaseComposition)
            temperature = c_double(self.temperature+273.15)
            pressure = c_double(self.pressure)

            phaseName = phaseList[i]
            if phaseName.startswith('bulk'):
                phaseName = 'liquid'
            phase = create_string_buffer(re.sub(r'(_+.*)|(\d*)', '', phaseName).encode())
            phaseName = phaseList[i]

            self.status.failed = True
            failure = c_int(self.status.failed)
            try:
                if pipes is not None:
                    with pipes():
                        if phaseName.startswith('bulk'):
                            self.status.libalphamelts.getMeltsOxideProperties(pointer(failure), phase, pointer(temperature), pointer(pressure), \
                                bulk, endmembers, pointer(nCharInName), pointer(numberEndMembers), properties)
                        else:
                            self.status.libalphamelts.getMeltsEndMemberProperties(pointer(failure), phase, pointer(temperature), pointer(pressure), \
                                bulk, endmembers, pointer(nCharInName), pointer(numberEndMembers), properties)
                else:
                    if phaseName.startswith('bulk'):
                        self.status.libalphamelts.getMeltsOxideProperties(pointer(failure), phase, pointer(temperature), pointer(pressure), \
                            bulk, endmembers, pointer(nCharInName), pointer(numberEndMembers), properties)
                    else:
                        self.status.libalphamelts.getMeltsEndMemberProperties(pointer(failure), phase, pointer(temperature), pointer(pressure), \
                            bulk, endmembers, pointer(nCharInName), pointer(numberEndMembers), properties)
                self.status.failed = failure.value
            except OSError:
                print("OSERROR")

            if self.status.failed:
                if self.status.reload("calcEndMemberProperties"):
                    self.setSystemProperties(self.systemProperties)

            else:

                # If phase is liquid and there is an fO2 buffer composition may have changed
                self.dispComposition = self.dispComposition or {}
                self.dispComposition[phaseName] = tnp.array(bulk)
                total = self.dispComposition[phaseName].sum()
                for j in range(nox):
                    self.dispComposition[phaseName][j] *= 100.0/total

                self.X = self.X or {}
                self.activity0 = self.activity0 or {}
                self.activity = self.activity or {}
                self.mu0 = self.mu0 or {}
                self.mu = self.mu or {}

                endMemberProps = tnp.array(properties[:4*numberEndMembers.value])
                endMemberProps = endMemberProps.reshape((numberEndMembers.value,4))

                self.X[phaseName] = endMemberProps[:,0]
                self.activity0[phaseName] = endMemberProps[:,1]
                self.mu0[phaseName] = endMemberProps[:,2]
                self.mu[phaseName] = endMemberProps[:,3]

                self.activity[phaseName] = self.activity0[phaseName].copy()
                for i in range(numberEndMembers.value):
                    # Set undefined mu, mu0 (and activities) to NaN
                    if self.mu0[phaseName][i] != 0.0 and self.mu[phaseName][i] != 0.0:
                        if np is not None:
                            self.activity[phaseName][i] = np.exp((self.mu[phaseName][i] - self.mu0[phaseName][i]) / (8.3143*(self.temperature+273.15)))
                        else:
                            self.activity[phaseName][i] = math.exp((self.mu[phaseName][i] - self.mu0[phaseName][i]) / (8.3143*(self.temperature+273.15)))
                    else:
                        self.activity[phaseName][i] = float("NaN")
                        self.mu[phaseName][i] = float("NaN")
                    if self.activity0[phaseName][i] == 0.0 and self.mu[phaseName][i] != 0.0:
                        self.activity0[phaseName][i] = float("NaN")
                    if self.mu0[phaseName][i] == 0.0 and self.mu[phaseName][i] != 0.0:
                        self.mu0[phaseName][i] = float("NaN")


    def calcSaturationState(self, *args):
        ''' Get the affinity of each phase not in the assemblage. Can be called before or after equilibration. Also called by findLiquidus.
            Currently does not update bulkComposition for fO2 buffer; this will change once subsolidus start is implemented.
         '''

        if self.calculationMode != self.status.getCalculationMode():
            if self.status.setCalculationMode(self.calculationMode):
                self.setSystemProperties(self.systemProperties)

        if len(args):
            self.setBulkComposition(args)

        assert self.temperature, "System temperature not set?"
        self.pressure = (self.pressure or 1)
        assert self.pressure, "System pressure not set?"
        assert any(self.bulkComposition), "System bulk composition not set?"

        nCharInName = c_int(20)
        numberPhases = c_int(100)
        indices = [0] * 100
        phaseIndices = (c_int * len(indices)) (*indices)
        phases = create_string_buffer(100*20)

        nox = len(self.status.endmembers['bulk'])
        nc = nox + 1
        propArr = [0.0] * (nc*100) # numberPhases
        properties = (c_double * len(propArr)) (*propArr)

        bulkComposition = self.bulkComposition
        bulk = (c_double * len(bulkComposition)) (*bulkComposition)
        temperature = c_double(self.temperature+273.15)
        pressure = c_double(self.pressure)

        try:
            self.status.failed = True
            failure = c_int(self.status.failed)
            if pipes is not None:
                with pipes():
                    self.status.libalphamelts.getMeltsSaturationState(pointer(failure), pointer(pressure), \
                        bulk, pointer(temperature), phases, pointer(nCharInName), pointer(numberPhases), properties, phaseIndices)
            else:
                self.status.libalphamelts.getMeltsSaturationState(pointer(failure), pointer(pressure), \
                    bulk, pointer(temperature), phases, pointer(nCharInName), pointer(numberPhases), properties, phaseIndices)
            self.status.failed = failure.value
        except OSError:
            print("OSERROR")

        if self.status.failed:
            if self.status.reload("getSaturationState"):
                self.setSystemProperties(self.systemProperties)
        else:

            self.affinity = self.affinity or {}
            self.dispComposition  = self.dispComposition or {}

            phaseList = []
            self.phaseNames = self.phaseNames or []

            for i in range(2, numberPhases.value): # not bulk system or oxygen

                if properties[i*nc]: # not in assemblage
                    phaseName = self.status.phases[i]
                    phaseList.append(phaseName)
                    self.phaseNames.append(phaseName)

                    if properties[i*nc] > -99999: # not failed
                        self.affinity[phaseName] = properties[i*nc]
                        self.dispComposition[phaseName] = tnp.array(properties[i*nc+1:(i+1)*nc])
                        total = self.dispComposition[phaseName].sum()
                        for j in range(nox):
                            # element-wise operations don't seem to be available in pip-packaged tinynumpy
                            self.dispComposition[phaseName][j] *= 100.0/total
                    else:
                        self.affinity[phaseName] = float("NaN")

        return phaseList

        ''' Will be option to get fO2 by different methods / compare to
             different buffers in here '''

    def copyAndKeepOutput(self, *args):
        ''' Deep copy the node - keeps (duplicates) same engine input and output as previous engine. '''
        if len(args):
            cp = MELTSengine(args[0])
        else:
            cp = MELTSengine(None)
        cp.calculationMode = copy(self.calculationMode)
        # Engine input is duplicated
        cp.nodeId = copy(self.nodeId)
        cp.nodeName = copy(self.nodeName)
        cp.runMode = copy(self.runMode)
        cp.pressure = copy(self.pressure)
        cp.temperature = copy(self.temperature)
        cp.reference = copy(self.reference)
        cp.systemProperties = copy(self.systemProperties)
        cp.bulkComposition = copy(self.bulkComposition)
        cp.phaseComposition = copy(self.phaseComposition)

        # Engine output is duplicated
        cp.g = deepcopy(self.g)
        cp.h = deepcopy(self.h)
        cp.s = deepcopy(self.s)
        cp.v = deepcopy(self.v)
        cp.cp = deepcopy(self.cp)
        cp.dcpdt = deepcopy(self.dcpdt)
        cp.dvdt = deepcopy(self.dvdt)
        cp.dvdp = deepcopy(self.dvdp)
        cp.d2vdt2 = deepcopy(self.d2vdt2)
        cp.d2vdtdp = deepcopy(self.d2vdtdp)
        cp.d2vdp2 = deepcopy(self.d2vdp2)
        cp.molwt = deepcopy(self.molwt)
        cp.rho = deepcopy(self.rho)
        cp.mass = deepcopy(self.mass)
        cp.viscosity = deepcopy(self.viscosity)
        cp.dispComposition = deepcopy(self.dispComposition)
        cp.affinity = deepcopy(self.affinity)
        cp.X = deepcopy(self.X)
        cp.activity = deepcopy(self.activity)
        cp.mu0 = deepcopy(self.mu0)
        cp.mu = deepcopy(self.mu)

        return cp

    def copyElement(self):
        ''' Copy the engine - keeps (duplicates) same engine input as previous engine; engine output is reset. '''
        # Handles point to the same entities, unless reset
        cp = MELTSengine(self.status)
        cp.calculationMode = copy(self.calculationMode)
        # Engine input is duplicated
        cp.nodeId = copy(self.nodeId)
        cp.nodeName = copy(self.nodeName)
        cp.runMode = copy(self.runMode)
        cp.pressure = copy(self.pressure)
        cp.temperature = copy(self.temperature)
        cp.reference = copy(self.reference)
        cp.systemProperties = copy(self.systemProperties)
        cp.bulkComposition = copy(self.bulkComposition)
        cp.phaseComposition = copy(self.phaseComposition)
        # Engine output is reset
        cp.g = {}
        cp.h = {}
        cp.s = {}
        cp.v = {}
        cp.cp = {}
        cp.dcpdt = {}
        cp.dvdt = {}
        cp.dvdp = {}
        cp.d2vdt2 = {}
        cp.d2vdtdp = {}
        cp.d2vdp2 = {}
        cp.molwt = {}
        cp.rho = {}
        cp.mass = {}
        cp.viscosity = {}
        cp.dispComposition = {}
        cp.affinity = {}
        cp.X = {}
        cp.activity = {}
        cp.mu0 = {}
        cp.mu = {}

        return cp
