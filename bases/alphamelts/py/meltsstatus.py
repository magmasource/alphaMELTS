import os, platform
from ctypes import CDLL, cdll, c_int, c_double, c_char_p, c_void_p, create_string_buffer, pointer #, byref
if platform.system() == 'Windows':
    from ctypes import windll

libIsLoaded = False

class MELTSstatus(object):

    libalphamelts = None
    calculationMode = None
    phases = None
    endmembers = None
    molwts = None
    fields = None
    loaded = None
    nodeIndex = 0
    meltsIndex = 0
    phaseIndex = 0
    nullValue = float("NaN")
    failed = None
    console = None
    message = None
    #finished


    def __init__(self, cMode):

        global libIsLoaded

        if platform.system() == 'Darwin':
            dylib = os.path.join(os.path.dirname(__file__), "libalphamelts.dylib")
            dllib = 'libdl.dylib'
        elif platform.system() == 'Linux':
            dylib = os.path.join(os.path.dirname(__file__), "libalphamelts.so")
            dllib = 'libdl.so'
        else:
            dylib = os.path.join(os.path.dirname(__file__), "libalphamelts.dll")
            dllib = 'libdl.dll'
        dynamicLib = CDLL(dylib)

        cMode = cMode or self.calculationMode
        self.calculationMode = cMode or None
        if libIsLoaded:
            windll.kernel32.FreeLibrary.argtypes = (c_void_p,)
            windll.kernel32.FreeLibrary(dynamicLib._handle)
            del dynamicLib
        self.libalphamelts = cdll.LoadLibrary(dylib)

        self.failed = True
        failure = c_int(self.failed)
        nCharInString = c_int(50)
        errorString = create_string_buffer(50)

        self.libalphamelts.getMeltsVersionString(pointer(failure), errorString, pointer(nCharInString))
        self.failed = failure.value
        self.message = errorString.value.decode()

        self.libalphamelts.addConsole()
        self.console = True

        self.calculationMode = 1 if (not self.calculationMode or self.calculationMode < 0) else self.calculationMode
        # Setting the calculation mode initializes the library
        if not self.libalphamelts.setCalculationMode(self.calculationMode):
            # The library was already initialized: this can be normal behavior on Windows / gcc.
            # It is a problem if the library calculation mode is not what Python thinks it is...
            cMode = self.libalphamelts.getCalculationMode()
            if cMode != self.calculationMode:
                print (f'Could not reset MELTS library calculation mode from {cMode:d} to {self.calculationMode:d}!\n Please save your work and restart Python.')
            else:
                print ('MELTS:libraryAlreadyInitialized. Could not re-initialize MELTS library! Please check stderr (console or terminal). If there are any error messages, please save your work and reload Python; otherwise you should be able to continue, with caution...')
        else:
            self.meltsIndex = 0
        self.nodeIndex = 0

        self.getSystemNamesAndWeights()
        self.fields = ('g', 'h', 's', 'v', 'cp', 'dcpdt', 'dvdt', 'dvdp', 'd2vdt2', 'd2vdtdp', 'd2vdp2', 'molwt', 'rho', 'mass')
        self.loaded = True
        libIsLoaded  = True if dllib is None else False

    def getSystemNamesAndWeights(self):

        self.failed = True
        failure = c_int(self.failed)
        nCharInName = c_int(20)
        numberPhases = c_int(100)
        indices = [0] * 100
        phaseIndices = (c_int * len(indices)) (*indices)
        phaseNames = create_string_buffer(100*20)

        self.libalphamelts.getMeltsPhaseNames(pointer(failure), phaseNames, pointer(nCharInName), pointer(numberPhases), phaseIndices)

        self.failed = failure.value
        phaseNames = phaseNames.value.decode()
        phaseNames = phaseNames.lower()
        self.phases = phaseNames.split()
        #self.phases.pop() # 0x0 char

        self.endmembers = {}
        self.molwts = {}
        for phase in self.phases:

            if phase == 'bulk':

                self.failed = True
                failure = c_int(self.failed)
                nCharInName = c_int(20)
                numberOxides = c_int(50)
                oxideNames = create_string_buffer(50*20)

                self.libalphamelts.getMeltsOxideNames(pointer(failure), oxideNames, pointer(nCharInName), pointer(numberOxides))

                numberOxides = c_int(50)
                oxideArr = [0.0] * (50)
                weights = (c_double * len(oxideArr)) (*oxideArr)

                self.libalphamelts.getMeltsOxideWeights(pointer(failure), weights, pointer(numberOxides))

                oxideNames = oxideNames.value.decode()
                oxideNames = oxideNames.lower()
                self.endmembers['bulk'] = oxideNames.split() # system
                self.endmembers['oxygen'] = self.endmembers['bulk'].pop()

                nox = numberOxides.value
                oxideWeights = [0.0] * nox
                for j in range(nox):
                    oxideWeights[j] = weights[j]
                self.molwts['bulk'] = oxideWeights.copy()
                self.molwts['oxygen'] = self.molwts['bulk'].pop()
                nox = nox - 1 # oxygen

            elif phase != 'oxygen':

                self.failed = True
                failure = c_int(self.failed)
                nCharInName = c_int(20)
                numberEndMembers = c_int(50)
                endMemberFormulas = create_string_buffer(50*20)

                endmemberArr = [0.0] * (50)
                weights = (c_double * len(endmemberArr)) (*endmemberArr)

                self.libalphamelts.getMeltsWeightsAndFormulas(pointer(failure), c_char_p(phase.encode()), \
                    weights, endMemberFormulas, pointer(nCharInName), pointer(numberEndMembers))

                endMemberFormulas = endMemberFormulas.value.decode()
                endMemberFormulas = endMemberFormulas.lower()
                self.endmembers[phase] = endMemberFormulas.split()

                na = numberEndMembers.value
                endmemberWeights = [0.0] * na
                for j in range(na):
                    endmemberWeights[j] = weights[j]
                self.molwts[phase] = endmemberWeights.copy()

    def getCalculationMode(self):

        if self.calculationMode is not None:
            cMode = self.libalphamelts.getCalculationMode()
        else:
            cMode = -1
        return cMode

    def setCalculationMode(self, cMode):

        global libIsLoaded

        if cMode == 1:
            print ('Setting calculation mode to rhyolite-MELTS 1.0.2.')
        elif cMode == 2:
            print ('Setting calculation mode to pMELTS')
        elif cMode == 3:
            print ('Setting calculation mode to rhyolite-MELTS 1.1.0.')
        elif cMode == 4:
            print ('Setting calculation mode to rhyolite-MELTS 1.2.0.')
        else:
            cMode = max(1, cMode)
            print ('Unexpected value for calculationMode. Using default (rhyolite-MELTS 1.0.2)')

        if platform.system() == 'Darwin':
            dylib = os.path.join(os.path.dirname(__file__), "libalphamelts.dylib")
            dllib = 'libdl.dylib'
        elif platform.system() == 'Linux':
            dylib = os.path.join(os.path.dirname(__file__), "libalphamelts.so")
            dllib = 'libdl.so'
        else:
            dylib = os.path.join(os.path.dirname(__file__), "libalphamelts.dll")
            dllib = 'libdl.dll'
        dynamicLib = CDLL(dylib)

        if dllib is None:
            windll.kernel32.FreeLibrary.argtypes = (c_void_p,)
            windll.kernel32.FreeLibrary(dynamicLib._handle)
            del dynamicLib
            self.loaded = libIsLoaded = False

            self.libalphamelts = cdll.LoadLibrary(dylib)
            if self.console:
                self.libalphamelts.addConsole()

        self.calculationMode = cMode
        success = self.libalphamelts.setCalculationMode(self.calculationMode)
        cMode = self.libalphamelts.getCalculationMode()

        if not success:
            if cMode != self.calculationMode:
                print (f'Could not reset MELTS library calculation mode from {cMode:d} to {self.calculationMode:d}!\n Please save your work and restart Python.')
            else:
                print ('MELTS:libraryAlreadyInitialized. Could not re-initialize MELTS library! Please check stderr (console or terminal). If there are any error messages, please save your work and restart Python; otherwise you should be able to continue, with caution...')
        else:
            print ('MELTS library has been reloaded and initialized.')
            self.getSystemNamesAndWeights()
            self.meltsIndex = 0

        self.loaded = True
        libIsLoaded  = True if dllib is None else False

        return success

    def reload(self, func):

        global libIsLoaded

        print (f'MELTS call {func:s} failed. Cleaning up...')

        if platform.system() == 'Darwin':
            dylib = os.path.join(os.path.dirname(__file__), "libalphamelts.dylib")
            dllib = 'libdl.dylib'
        elif platform.system() == 'Linux':
            dylib = os.path.join(os.path.dirname(__file__), "libalphamelts.so")
            dllib = 'libdl.so'
        else:
            dylib = os.path.join(os.path.dirname(__file__), "libalphamelts.dll")
            dllib = 'libdl.dll'
        dynamicLib = CDLL(dylib)

        if dllib is None:
            if libIsLoaded:
                windll.kernel32.FreeLibrary.argtypes = (c_void_p,)
                windll.kernel32.FreeLibrary(dynamicLib._handle)
                del dynamicLib
                self.loaded = libIsLoaded = False

            self.libalphamelts = cdll.LoadLibrary(dylib)
            if self.console:
                self.libalphamelts.addConsole()

        success = self.libalphamelts.setCalculationMode(self.calculationMode)

        if not success:
            print ('Could not re-initialize MELTS library after failure! Please check stderr (console or terminal) for failure details, then save your work and restart Python.')
        else:
            print ('MELTS:calculationFailed. MELTS library has been reloaded and initialized. Please check stderr (console or terminal) for failure details.')
            self.meltsIndex = 0

        self.loaded = True
        libIsLoaded  = True if dllib is None else False

        return success
