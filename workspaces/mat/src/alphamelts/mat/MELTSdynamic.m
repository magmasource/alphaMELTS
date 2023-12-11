classdef MELTSdynamic < matlab.mixin.Copyable
    % Given an integer for the calculation mode, (re-)loads the C library of MELTS functions, and initializes the node/list and engine.
    % Calculation mode needs to be set to one of: 1. rhyolite-MELTS 1.0.2; 2. pMELTS; 3. rhyolite-MELTS 1.1.0; 4. rhyolite-MELTS 1.2.0.

    properties (SetAccess = protected)

        % The alphaMELTS for MATLAB/Python version (for reference).
        % Will be updated from the C library of MELTS functions.
        version = '1.2.0beta'
        % User-defined name and/or description for the MELTSdynamic list.
        title
        % Integer value used to tell the C library of MELTS functions which MELTS mode to use.
        calculationMode
        % The string (equivalent to ALPHAMELTS_CALC_MODE) that corresponds to calculationMode.
        modeString
        % A string array of possible phase names for the current MELTS calculation mode (also includes 'bulk').
        systemNames
        % A MELTSmap of end-member formulae for the possible phases in the current MELTS system (or system oxides for 'bulk').
        endMemberFormulas
        % A MELTSmap of end-member molecular weights for the possible phases in the current MELTS system (or system oxides for 'bulk').
        endMemberWeights
        % Used to keep track of MELTS calculation success and triggers reloadOnFailure if needed.
        failureHandle
    end

    properties
        % The position of the node in the MELTSdynamic list (First has nodeIndex = 1).
        nodeIndex
        % The thermodynamic 'engine' where all the real work is done.
        engine
    end

    properties(SetAccess = private)
        % The next node in MELTSdynamic list, relative to the current node.
        Next = MELTSdynamic.empty
        % The previous node in MELTSdynamic list, relative to the current node.
        Prev = MELTSdynamic.empty
        % The node at the start of the MELTSdynamic list that the current node is in.
        First = MELTSdynamic.empty
        % The node at the end of the MELTSdynamic list that the current node is in.
        Last = MELTSdynamic.empty
    end

    methods

        function obj = MELTSdynamic(cMode, varargin)
            % Given an integer for the calculation mode, (re-)loads the C library of MELTS functions, and initializes the node/list and engine. Call without arguments to see options.

            if nargin == 0 || isempty(cMode)
                disp (['Calculation mode needs to be set:';...
                    '     1. rhyolite-MELTS 1.0.2     ';...
                    '     2. pMELTS                   ';...
                    '     3. rhyolite-MELTS 1.1.0     ';...
                    '     4. rhyolite-MELTS 1.2.0     ';...
                    'Closing console (if any)...      ']);
                if libisloaded('libalphamelts'); calllib('libalphamelts', 'closeConsole'); end
                obj = MELTSdynamic.empty;
                return
            end

            obj.engine = MELTSengine(cMode, varargin{:});
            obj.version = obj.engine.status.message;

            % Only print if cMode is an integer (i.e. first time)...
            if ~isa(cMode, 'MELTSstatus')
                disp(newline);
                disp(['alphaMELTS for MATLAB version ' obj.version])
                obj.printEngineInit(cMode);
            end

            obj.systemNames = obj.engine.status.phases;
            obj.endMemberFormulas = obj.engine.status.endmembers;
            obj.endMemberWeights = obj.engine.status.molwts;

            obj.nodeIndex = obj.engine.status.nodeIndex;
            obj.calculationMode = obj.engine.calculationMode;
            obj.addListener;

            obj.First = obj;
            obj.Last = obj;

        end

        function printEngineInit(obj, cMode)
            % Sets the modeString (equivalent to ALPHAMELTS_CALC_MODE) and prints.
            switch cMode
                case 1
                    obj.modeString = 'MELTS';
                    disp ('Setting calculation mode to rhyolite-MELTS 1.0.2');
                case 2
                    obj.modeString = 'pMELTS';
                    disp ('Setting calculation mode to pMELTS');
                case 3
                    obj.modeString = 'MELTSandCO2';
                    disp ('Setting calculation mode to rhyolite-MELTS 1.1.0');
                case 4
                    obj.modeString = 'MELTSandCO2_H2O';
                    disp ('Setting calculation mode to rhyolite-MELTS 1.2.0');
                otherwise
                    obj.modeString = 'MELTS';
                    warning ('Unexpected value for calculationMode. Using default (rhyolite-MELTS 1.0.2)');
            end
            disp(newline);

        end

        function toggleConsole(obj)
            % On Windows MATLAB, use this to open and close the console where output from the C library of MELTS functions is directed.
            if obj.engine.status.console
                calllib('libalphamelts', 'closeConsole');
            else
                calllib('libalphamelts', 'addConsole');
            end
            obj.engine.status.console = ~obj.engine.status.console;
        end

        function addListener(obj)
                obj.failureHandle = event.listener(obj.engine.status, 'finished',...
                @(src, evnt) reloadOnFailure(obj, src, evnt));
        end

        function reloadOnFailure(obj, ~, evnt)
            % Tries to reload the C library of MELTS functions if a calculation fails (do not call directly).
            if obj.engine.status.failed
                func = evnt.funcName;
                success = obj.engine.status.reload(func);
                if success
                    % Load previous setSystemProperties
                    obj.engine.setSystemProperties(obj.engine.systemProperties)
                    obj.engine.status.failed = true;
                end
            end
        end

        function value = getNodeProperty(obj, inode, varargin)
            % Get some property (for system, or one or more phases) from node with given index. Returns a scalar, vector or 2-D matrix.

            if isempty(inode); inode = obj.nodeIndex; end
            node = obj.findIndex(inode);
            assert(~isempty(node), "Node index not found in list.");
            value = node.engine.getProperty(varargin{:});

        end

        function value = getListProperty(obj, propertyName, varargin)
            % Get some property from the list. Returns a vector or 2-D matrix.

            if nargin > 2 && ~strcmp(propertyName, 'bulkComposition') && ...
                ~strcmp(propertyName, 'phaseComposition') && ~strcmp(propertyName, 'molarComposition')
                % phaseList is varargin{1}
                phaseList = string(varargin{1});
                assert(length(phaseList) == 1, ...
                    "The getListProperty method cannot be called for multiple phases. Use getNodeProperty or engine.getProperty instead.");
            end

            node = obj.First;
            property = node.engine.getProperty(propertyName, varargin{:});

            assert(iscolumn(property), "Something went wrong. Please report this to the developers.");

            value = NaN(length(property), obj.Last.nodeIndex);
            value(:, node.nodeIndex) = property;
            while ~isempty(node.Next)
                node = node.Next;
                value(:, node.nodeIndex) = node.engine.getProperty(propertyName, varargin{:});
            end

        end

        function newNode = addNodeAfter(nodeBefore)
            % Creates a copy of this node (with same engine input) and inserts it before this node. Returns the new node and udpates MELTSstatus.
            newNode = copy(nodeBefore);
            % Number of nodes pointing to MELTSstatus is incremented;
            newNode.engine.status.nodeIndex = newNode.engine.status.nodeIndex + 1;
            newNode.nodeIndex = newNode.engine.status.nodeIndex;
            insertNodeAfter(newNode, nodeBefore);
        end

        function newNode = addNodeBefore(nodeAfter)
            % Creates a copy of this node (with same engine input) and inserts it before this node. Returns the new node and updates MELTSstatus.
            newNode = copy(nodeAfter);
            % Number of nodes pointing to MELTSstatus is incremented;
            newNode.engine.status.nodeIndex = newNode.engine.status.nodeIndex + 1;
            newNode.nodeIndex = newNode.engine.status.nodeIndex;
            insertNodeBefore(newNode, nodeAfter);
        end

        function insertNodeAfter(node, nodeBefore)
            % Insert this node after nodeBefore.
            assert(node ~= nodeBefore, "Inserting a node after itself?");
            if isempty(node.Next) || isempty(node.Prev)
                % Node at start or end of a list
                newNode = node;
            elseif ~isempty(node.First) && (node.First == nodeBefore.First) &&...
                    ~isempty(node.Last) && (node.Last == nodeBefore.Last)
                % Node in middle of same list
                removeNode(node);
                newNode = node;
            else
                % If node belongs in the middle of a different linked list, make a deep copy
                % Note that meltsIndex will be for any files generated with nodeBefore's status
                % object and may not be correct for nodeAfter's (flagged as negative).
                newNode = node.copyAndKeepOutput(nodeBefore.engine.status);
                if ~isempty(nodeBefore.engine.meltsIndex)
                    newNode.engine.meltsIndex = -abs(nodeBefore.engine.meltsIndex);
                end
            end

            % Insert newNode after nodeBefore.
            newNode.Prev = nodeBefore;
            if ~isempty(nodeBefore.Next)
                newNode.Next = nodeBefore.Next;
                nodeBefore.Next.Prev = newNode;
            end
            nodeBefore.Next = newNode;
            % Check whether Last needs updating
            if nodeBefore ~= nodeBefore.Last
                newNode.Last = nodeBefore.Last;
            else
                if isempty(newNode.Last)
                    newNode.Last = newNode;
                end
                resetLast(nodeBefore, newNode.Last);
            end
            % First stays the same
            newNode.First = nodeBefore.First;
            % newNode.First is always correct
            if newNode.nodeIndex ~= (nodeBefore.nodeIndex + 1)
                resetIndex(newNode);
            end
        end

        function joinListAfter(list, listBefore)
            % Append this list to listBefore.
            assert((list.First ~= listBefore.First) || (list.Last ~= listBefore.Last), "Already part of same list?");
            % Insert list after listBefore.
            insertNodeAfter(list.First, listBefore.Last);
            % Fix up First for rest of list
            resetFirst(list, list.First);
        end

        function insertNodeBefore(node, nodeAfter)
            % Insert this node before nodeAfter.
            assert(node ~= nodeAfter, "Inserting a node before itself?");

            if isempty(node.Next) || isempty(node.Prev)
                % Node at start or end of a list
                newNode = node;
            elseif ~isempty(node.First) && (node.First == nodeAfter.First) &&...
                    ~isempty(node.Last) && (node.Last == nodeAfter.Last)
                % Node in middle of same list
                removeNode(node);
                newNode = node;
            else
                % If node belongs in the middle of a different linked list, make a deep copy
                % Note that meltsIndex will be for any files generated with nodeAfter's status
                % object and may not be correct for nodeBefore's (flagged as negative).
                newNode = node.copyAndKeepOutput(nodeAfter.engine.status);
                if ~isempty(nodeAfter.engine.meltsIndex)
                    newNode.engine.meltsIndex = -abs(nodeAfter.engine.meltsIndex);
                end
            end

            % Insert newNode before nodeAfter.
            newNode.Next = nodeAfter;
            if ~isempty(nodeAfter.Prev)
                newNode.Prev = nodeAfter.Prev;
                nodeAfter.Prev.Next = newNode;
            end
            nodeAfter.Prev = newNode;
            % Check whether First needs updating
            if nodeAfter ~= nodeAfter.First
                newNode.First = nodeAfter.First;
            else
                if isempty(newNode.First)
                    newNode.First = newNode;
                end
                resetFirst(nodeAfter, newNode.First);
            end
            % Last stays the same
            newNode.Last = nodeAfter.Last;
            % newNode.First is always correct
            if newNode.nodeIndex ~= (nodeAfter.nodeIndex - 1)
                resetIndex(newNode);
            end
        end

        function joinListBefore(list, listAfter)
            % Append listAfter to this list.
            assert((list.First ~= listAfter.First) || (list.Last ~= listAfter.Last), "Already part of same list?");
            % Insert list before listAfter.
            insertNodeBefore(list.Last, listAfter.First);
            % Fix up Last for rest of list
            resetLast(list, list.Last);
        end

        function resetFirst(node, firstNode)
            % Set First to firstNode for all nodes in list.
            node = node.Last;
            while ~isempty(node.Prev)
                node.First = firstNode;
                node = node.Prev;
            end
        end

        function resetLast(node, lastNode)
            % Set Last to lastNode for all nodes in list.
            node = node.First;
            while ~isempty(node.Next)
                node.Last = lastNode;
                node = node.Next;
            end
        end

        function resetIndex(node)
            % Reset the indices of all nodes in list.
            node = node.First;
            node.nodeIndex = 1;
            while ~isempty(node.Next)
                node.Next.nodeIndex = node.nodeIndex +1;
                node = node.Next;
            end
        end

        function node = findIndex(list, inode)
            % Returns the mode with matching index.
            node = list.First;
            while ~isempty(node.Next)
                if node.nodeIndex == inode
                    break
                else
                    node = node.Next;
                end
            end
            if node.nodeIndex ~= inode; node = MELTSdynamic.empty; end
        end

        function list = sortListBy(list, orderOut)
            % Takes in a list of sorted indices (orderOut) and returns a sorted list. '''
            assert(length(orderOut) == list.Last.nodeIndex, "Sort order not the same length as list?");
            assert(all(sort(orderOut) == 1:length(orderOut)), "Indices (1-n) do not appear exactly once in sort order?");

            if length(orderOut) > 1

                orderIn(orderOut) = 1:length(orderOut);
                newList = findIndex(list, find(orderIn == 1));

                node = list.First;
                for i = 1:length(orderOut)
                    node.nodeIndex = orderIn(i);
                    node = node.Next;
                end

                if newList == newList.First
                    resetFirst(list, newList.Next);
                elseif newList == newList.Last
                    resetLast(list, newList.Prev);
                end
                list = findIndex(list, length(orderOut));

                removeNode(newList);
                newList.First = newList;
                newList.Last = newList;

                nodeBefore = newList;
                for i = 2:length(orderOut)
                    node = findIndex(list, i);
                    if node == list.First
                        resetFirst(list, node.Next);
                    elseif node == list.Last
                        resetLast(list, node.Prev);
                    end
                    removeNode(node);
                    insertNodeAfter(node, nodeBefore);
                    nodeBefore = node;
                end

            end % Length > 1

        end

        function clearNode(node)
            % Delete node and fix the list so that remaining nodes are properly connected and global list properties are updated. Copy node before if needed to start a new list.
            prev = node.Prev;
            next = node.Next;
            removeNode(node);
            if ~isempty(next)
                resetIndex(next)
            elseif ~isempty(prev)
                resetIndex(prev)
            end
        end

        function clearList(node)
            % Clear the list before clearing list variable
            prev = node.Prev;
            next = node.Next;
            removeNode(node)
            while ~isempty(next)
                node = next;
                next = node.Next;
                removeNode(node);
            end
            while ~isempty(prev)
                node = prev;
                prev = node.Prev;
                removeNode(node)
            end
        end

        function cp = copyAndKeepOutput(obj, varargin)
            % Deep copy of the node - keeps (duplicates) same engine output as previous engine.
            cp = MELTSdynamic(obj.calculationMode);
            if nargin > 1
                cp.engine = obj.engine.copyAndKeepOutput(varargin{:});
            else
                cp.engine = obj.engine.copyAndKeepOutput(cp.engine.status);
            end

            cp.version = obj.version;
            cp.title = obj.title;
            cp.calculationMode = obj.calculationMode;
            cp.modeString = obj.modeString;
            cp.systemNames = obj.systemNames;
            cp.endMemberFormulas = obj.endMemberFormulas;
            cp.failureHandle = obj.failureHandle;

            % Reset settings if appropriate
            cp.engine.setSystemProperties(cp.engine.systemProperties);
            % Node links are reset
            cp.Next = MELTSdynamic.empty;
            cp.Prev = MELTSdynamic.empty;
            cp.First = cp;
            cp.Last = cp;

        end

    end

    methods (Access = protected)

        function cp = copyElement(obj)
            % Copy the node - keeps (duplicates) same engine input as previous engine; engine output is reset.
            % Handles point to the same entities
            cp = copyElement@matlab.mixin.Copyable(obj);
            % Console flag is unaffected
            cp.engine = copy(obj.engine);
            % Reload library if doesn't match the current calculationMode
            if cp.engine.calculationMode ~= cp.engine.status.getCalculationMode
                success = cp.engine.status.setCalculationMode(cp.calculationMode);
                if success; cp.engine.setSystemProperties(cp.engine.systemProperties); end
            end
            % Node links are reset
            cp.Next = MELTSdynamic.empty;
            cp.Prev = MELTSdynamic.empty;
            cp.First = cp;
            cp.Last = cp;
        end

    end

    methods (Access = private)

        function removeNode(node)
            % Delete node and fix the list so that remaining nodes are properly connected. Does not fix up global properties of list e.g. nodeIndex, to avoid recursion.
            if ~isscalar(node)
                error('Node input must be scalar')
            end
            prevNode = node.Prev;
            nextNode = node.Next;
            if ~isempty(prevNode)
                prevNode.Next = nextNode;
            end
            if ~isempty(nextNode)
                nextNode.Prev = prevNode;
            end
            node.Next = MELTSdynamic.empty;
            node.Prev = MELTSdynamic.empty;
            node.First = MELTSdynamic.empty;
            node.Last = MELTSdynamic.empty;
        end

        % Class destructor
        function delete(node)
            clearList(node)
        end

        function saveNode = saveobj(node)
            saveNode = node;
            saveNode.failureHandle = event.listener.empty;
        end

    end

end
