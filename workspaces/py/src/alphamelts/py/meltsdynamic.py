from tinynumpy import tinynumpy as tnp
from copy import copy, deepcopy

from meltsengine import MELTSengine
from meltsstatus import MELTSstatus

class MELTSdynamic(object):
    ''' Given an integer for the calculation mode, (re-)loads the C library of MELTS functions, and initializes the node/list and engine.
        Calculation mode needs to be set to one of: 1. rhyolite-MELTS 1.0.2; 2. pMELTS; 3. rhyolite-MELTS 1.1.0; 4. rhyolite-MELTS 1.2.0. '''

    # The alphaMELTS for MATLAB/Python version (for reference). Will be updated from the C library of MELTS functions.
    version = '1.2.0beta'
    # User-defined name and/or description for the MELTSdynamic list.
    title = None
    # Integer value used to tell the C library of MELTS functions which MELTS mode to use.
    calculationMode = None
    # The string (equivalent to ALPHAMELTS_CALC_MODE) that corresponds to the MELTS calculation mode.
    modeString = None
    # A list of possible phase names for the current MELTS calculation mode (also includes 'bulk').
    systemNames = None
    # A dictionary of end-member formulae for the possible phases in the current MELTS system (or system oxides for 'bulk')..
    endMemberFormulas = None
    # A dictionary of end-member molecular weights for the possible phases in the current MELTS system (or system oxides for 'bulk')..
    endMemberWeights = None
    #failureHandle

    # The position of the node in the MELTSdynamic list (First has nodeIndex = 0).
    nodeIndex = 0
    # The thermodynamic 'engine' where all the real work is done.
    engine = None
    # The next node in MELTSdynamic list, relative to the current node.
    Next  = None
    # The previous node in MELTSdynamic list, relative to the current node.
    Prev  = None
    # The node at the start of the MELTSdynamic list that the current node is in.
    First = None
    # The node at the end of the MELTSdynamic list that the current node is in.
    Last  = None

    def __init__(self, cMode=None, *args, **kwargs):
        ''' Given an integer for the calculation mode, (re-)loads the C library of MELTS functions,
            and initializes the node/list and engine. Call without arguments to see options. '''

        if 'calculationMode' in kwargs:
            cMode = kwargs['calculationMode']
        elif cMode is None:
            print('Calculation mode needs to be set:')
            print('     1. rhyolite-MELTS 1.0.2     ')
            print('     2. pMELTS                   ')
            print('     3. rhyolite-MELTS 1.1.0     ')
            print('     4. rhyolite-MELTS 1.2.0     ')
            # print('Closing console (if any)...      ')
            # self.engine.status.libalphamelts.closeConsole()
            return None

        self.engine = MELTSengine(cMode, *args, **kwargs)
        self.version = self.engine.status.message

        # Only print if cMode is an integer (i.e. first time)...
        if not isinstance(cMode, MELTSstatus):
            print (f'\nalphaMELTS for Python version {self.version:s}')
            self.printEngineInit(cMode)

        self.systemNames = self.engine.status.phases
        self.endMemberFormulas = self.engine.status.endmembers
        self.endMemberWeights = self.engine.status.molwts

        self.nodeIndex = self.engine.status.nodeIndex
        self.calculationMode = cMode if cMode != 0 else 1

        self.First = self
        self.Last = self

    def printEngineInit(self, cMode):
        ''' Sets the modeString (equivalent to ALPHAMELTS_CALC_MODE) and prints. '''
        if cMode == 1:
            self.modeString = 'MELTS'
            print ('Setting calculation mode to rhyolite-MELTS 1.0.2\n')
        elif cMode == 2:
            self.modeString = 'pMELTS'
            print ('Setting calculation mode to pMELTS\n')
        elif cMode == 3:
            self.modeString = 'MELTSandCO2'
            print ('Setting calculation mode to rhyolite-MELTS 1.1.0\n')
        elif cMode == 4:
            self.modeString = 'MELTSandCO2_H2O'
            print ('Setting calculation mode to rhyolite-MELTS 1.2.0\n')
        else:
            self.modeString = 'MELTS'
            print ('Unexpected value for calculationMode. Using default (rhyolite-MELTS 1.0.2)\n')

    def toggleConsole(self):
        ''' On Windows, use this to open and close the console where output from the
            C library of MELTS functions is directed. Not known to be used by any Python IDE. '''
        if self.engine.status.console:
            self.engine.status.libalphamelts.closeConsole()
        else:
            self.engine.status.libalphamelts.addConsole()

        self.engine.status.console = not(self.engine.status.console)

    def getNodeProperty(self, inode=None, *args):
        ''' Get some property (for system, or one or more phases) from node with given index.
            Returns a scalar, vector or 2-D matrix. '''

        if inode is None:
            inode = self.nodeIndex
        node = self.findIndex(inode)
        assert (node is not None), "Node index not found in list."
        value = node.engine.getProperty(*args)

        return value

    def getListProperty(self, propertyName, *args):
        ''' Get some property from the list. Returns a vector or 2-D array. '''

        if len(args) and propertyName not in ['bulkComposition', 'phaseComposition', 'molarComposition']:
            phaseList = args[0]
            if not isinstance(phaseList, list):
                phaseList = [phaseList]
            assert len(phaseList) == 1, \
                "The getListProperty method cannot be called for multiple phases. Use getNodeProperty or engine.getProperty instead."

        node = self.First
        property = node.engine.getProperty(propertyName, *args)

        if isinstance(property, tnp.ndarray):

            assert len(property.shape) == 1, "Something went wrong. Please report this to the developers."

            values = tnp.empty((node.Last.nodeIndex+1, len(property)))
            values[node.nodeIndex][:] = property

            node = node.Next
            while node is not None:
                property = node.engine.getProperty(propertyName, *args)
                values[node.nodeIndex][:] = property
                node = node.Next
            values = values.transpose()
        else:
            values = tnp.empty((node.Last.nodeIndex+1,))
            values[0] = property

            none = node.Next
            while node is not None:
                values[node.nodeIndex] = node.engine.getProperty(propertyName, *args)
                node = node.Next

        return values

    def addNodeAfter(self):
        ''' Creates a copy of this node (with same engine input) and inserts it after this node.
            Returns the new node and updates MELTSstatus. '''
        newNode = self.copyElement()
        # Number of nodes pointing to MELTSstatus is incremented
        newNode.engine.status.nodeIndex = newNode.engine.status.nodeIndex + 1
        newNode.nodeIndex = newNode.engine.status.nodeIndex
        newNode.insertNodeAfter(self)
        return newNode

    def addNodeBefore(self):
        ''' Creates a copy of this node (with same engine input) and inserts it before this node.
            Returns the new node and updates MELTSstatus. '''
        newNode = self.copyElement()
        # Number of nodes pointing to MELTSstatus is incremented
        newNode.engine.status.nodeIndex = newNode.engine.status.nodeIndex + 1
        newNode.nodeIndex = newNode.engine.status.nodeIndex
        newNode.insertNodeBefore(self)
        return newNode

    def insertNodeAfter(self, nodeBefore):
        ''' Insert this node after nodeBefore. '''
        assert self != nodeBefore, "Inserting node after itself?"

        if self.Next is None or self.Prev is None:
            # Node at start or end of a list
            newNode = self
        elif self.First is not None and (self.First == nodeBefore.First) and \
            self.Last is not None and (self.Last == nodeBefore.First):
            # Node in middle of same list
            self.removeNode()
            newNode = self
        else:
            # If node belongs in the middle of a different linked list, make a deep copy
            # Note that meltsIndex will be for any files generated with nodeBefore's status
            # object and may not be correct for nodeAfter's (flagged as negative).
            newNode = self.copyAndKeepOutput(nodeBefore.engine.status)
            if nodeBefore.engine.meltsIndex is not None:
                newNode.engine.meltsIndex = -abs(nodeBefore.engine.meltsIndex)

        # Insert newNode after nodeBefore
        newNode.Prev = nodeBefore
        if nodeBefore.Next is not None:
            newNode.Next = nodeBefore.Next
            nodeBefore.Next.Prev = newNode
        nodeBefore.Next = newNode
        # Check whether Last needs updating
        if nodeBefore != nodeBefore.Last:
            newNode.Last = nodeBefore.Last
        else:
            if newNode.Last is None:
                newNode.Last = newNode
            nodeBefore.resetLast(newNode.Last)
        # First stays the same
        newNode.First = nodeBefore.First
        if newNode.nodeIndex != nodeBefore.nodeIndex + 1:
            newNode.resetIndex()

    def joinListAfter(self, listBefore):
        ''' Append this list to listBefore. '''
        assert (self.First != listBefore.First and self.Last != listBefore.Last), "Already part of same list?"
        # Insert list after listBefore.
        self.First.insertNodeAfter(listBefore.Last)
        # Fix up Last for rest of list
        self.resetFirst(self.First)

    def insertNodeBefore(self, nodeAfter):
        ''' Insert this node before nodeAfter. '''
        assert self != nodeAfter, "Inserting a node before itself?"

        if self.Next is None or self.Prev is None:
            # Node at start or end of a list
            newNode = self
        elif self.First is not None and (self.First == nodeAfter.First) and \
            self.Last is not None and (self.Last == nodeAfter.Last):
            # Node in middle of same list
            self.removeNode()
            newNode = self
        else:
            # If node belongs in the middle of a different linked list, make a deep copy
            # Note that meltsIndex will be for any files generated with nodeAfter's status
            # object and may not be correct for nodeBefore's (flagged as negative).
            newNode = self.copyAndKeepOutput(nodeAfter.engine.status)
            if nodeAfter.engine.meltsIndex is not None:
                newNode.engine.meltsIndex = -nodeAfter.engine.meltsIndex

        # Insert newNode before nodeAfter
        newNode.Next = nodeAfter
        if nodeAfter.Prev is not None:
            newNode.Prev = nodeAfter.Prev
            nodeAfter.Prev.Next = newNode
        nodeAfter.Prev = newNode
        # Check whether First needs updating
        if nodeAfter != nodeAfter.First:
            newNode.First = nodeAfter.First
        else:
            if newNode.First is None:
                newNode.First = newNode
            nodeAfter.resetFirst(newNode.First)
        # Last stays the same
        newNode.Last = nodeAfter.Last
        # newNode.First is always correct
        if newNode.nodeIndex != nodeAfter.nodeIndex - 1:
            newNode.resetIndex()

    def joinListBefore(self, listAfter):
        ''' Append listAfter to this list. '''
        assert (self.First != listAfter.First and self.Last != listAfter.Last), "Already part of same list?"
        # Insert list before listAfter.
        self.Last.insertNodeBefore(listAfter.First)
        # Fix up Last for rest of list
        self.resetLast(self.Last)

    def resetFirst(self, firstNode):
        ''' Set First to firstNode for all nodes in list. '''
        node = self.Last
        while node is not None:
            node.First = firstNode
            node = node.Prev

    def resetLast(self, lastNode):
        ''' Set Last to lastNode for all nodes in list. '''
        node = self.First
        while node is not None:
            node.Last = lastNode
            node = node.Next

    def resetIndex(self):
        ''' Reset the indices of all nodes in list. '''
        node = self.First
        node.nodeIndex = 0
        while node.Next is not None:
            node.Next.nodeIndex = node.nodeIndex + 1
            node = node.Next

    def findIndex(self, inode):
        ''' Returns the node with matching index. '''
        node = self.First
        while node is not None:
            if node.nodeIndex == inode:
                return node
            node = node.Next
        return None

    def sortListBy(self, orderOut):
        ''' Takes in a list of sorted indices (orderOut) and returns a sorted list. '''
        assert len(orderOut) == self.Last.nodeIndex + 1, "Sort order not the same length as list?"
        orderIn = sorted(orderOut)
        for i in range(len(orderIn)):
            assert orderIn[i] == i, 'Indices (0-n) do not appear exactly once in sort order?'

        if len(orderOut) == 1:
            return self

        for i in range(len(orderIn)):
            orderIn[orderOut[i]] = i
        newList = self.findIndex(orderIn.index(0))

        node = self.First
        for i in range(len(orderOut)):
            node.nodeIndex = orderIn[i]
            node = node.Next

        if newList == newList.First:
            self.resetFirst(newList.Next)
        elif newList == newList.Last:
            self.resetLast(newList.Prev)
        self = self.findIndex(len(orderOut)-1)

        newList.removeNode()
        newList.First = newList
        newList.Last = newList

        nodeBefore = newList
        for i in range(1, len(orderOut)):
            node = self.findIndex(i)
            if node == self.First:
                self.resetFirst(node.Next)
            elif node == self.Last:
                self.resetLast(node.Prev)
            node.removeNode()
            node.insertNodeAfter(nodeBefore)
            nodeBefore = node

        return self

    def clearNode(self):
        ''' Delete node and fix the list so remaining nodes are properly
            connected and global list properties are updated.
            Copy node before if needed to start a new list. '''
        Prev = self.Prev
        Next = self.Next
        self.removeNode()
        if Next is not None:
            Next.resetIndex()
        if Prev is not None:
            Prev.resetIndex()

    def clearList(self):
        ''' Clear list before clearing list variable. '''
        Prev = self.Prev
        Next = self.Next
        self.removeNode()
        while Next is not None:
            node = Next
            Next = node.Next
            node.removeNode()
        while Prev is not None:
            node = Prev
            Prev = node.Prev
            node.removeNode()

    def copyAndKeepOutput(self, *args):
        ''' Deep copy the node - keeps (duplicates) same engine input and output as previous engine. '''
        cp = MELTSdynamic(self.calculationMode)
        if len(args):
            cp.engine = self.engine.copyAndKeepOutput(*args)
        else:
            cp.engine = self.engine.copyAndKeepOutput(cp.engine.status)

        cp.calculationMode = self.calculationMode
        cp.modeString = self.modeString
        cp.systemNames = self.systemNames
        cp.endMemberFormulas = self.endMemberFormulas

        # Reset settings if appropriate
        cp.engine.setSystemProperties(cp.engine.systemProperties)
        # Node links are reset
        cp.Next = None
        cp.Prev = None
        cp.First = cp
        cp.Last = cp
        return cp

    def copyElement(self):
        ''' Copy the node - keeps (duplicates) same engine input as previous engine; engine output is reset. '''
        cp = MELTSdynamic(self.engine.status)
        cp.engine = self.engine.copyElement()
        # Reset settings if appropriate
        if cp.engine.calculationMode != cp.engine.status.getCalculationMode():
            cp.engine.status.setCalculationMode(cp.calculationMode)
            cp.engine.setSystemProperties(cp.engine.systemProperties)
        # Node links are reset
        cp.Next = None
        cp.Prev = None
        cp.First = cp
        cp.Last = cp
        return cp

    def removeNode(self):
        ''' Delete node and fix the list so that remaining nodes are properly connected.
            Does not fix up global properties of list e.g. nodeIndex, to avoid recursion. '''
        prevNode = self.Prev
        nextNode = self.Next
        if prevNode is not None:
            prevNode.Next = nextNode
        if nextNode is not None:
            nextNode.Prev = prevNode
        self.Next = None
        self.Prev = None
        self.First = None
        self.Last = None
