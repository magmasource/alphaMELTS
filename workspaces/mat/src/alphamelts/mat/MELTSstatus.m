classdef MELTSstatus < matlab.mixin.SetGet

    properties (SetAccess = protected)
        calculationMode = []
        phases
        endmembers
        molwts
        fields
    end

    properties
        loaded = false
        nodeIndex = 1
        meltsIndex = []
        phaseIndex = 1
        nullValue = NaN
        failed = false
        console = false
        message
    end

    events
        finished
    end

    methods

        function obj = MELTSstatus(cMode)

            if nargin == 0; cMode = []; end

            obj.calculationMode = cMode;
            if libisloaded('libalphamelts'); unloadlibrary('libalphamelts'); end
            loadlibrary('libalphamelts', 'MELTSdynamic.h', 'mfilename', 'libalphamelts_proto.m');

            nCharInString = 50;
            errorString = blanks(nCharInString);
            [~, errorString, ~] = calllib('libalphamelts', 'getMeltsVersionString',...
                true, errorString, nCharInString);
            obj.set('message', errorString);

            calllib('libalphamelts', 'addConsole');
            obj.console = true;

            if isempty(obj.calculationMode) || obj.calculationMode <= 0; obj.calculationMode = 1; end
            % Setting the calculation mode initializes the library
            success = calllib('libalphamelts', 'setCalculationMode', obj.calculationMode);
            cMode = calllib('libalphamelts', 'getCalculationMode');
            if ~success
                % The library was already initialized: this can be normal behavior on Windows / gcc.
                % It is a problem if the library calculation mode is not what MATLAB thinks it is...
                if cMode ~= obj.calculationMode
                    error(['Could not reset MELTS library calculation mode from %d to %d! ',...
                        'Please save your work and reload MATLAB.'], cMode, obj.calculationMode);
                else
                    warning('MELTS:libraryAlreadyInitialized', ['Could not re-initialize MELTS library! ',...
                        'Please check stderr (console or terminal). ',...
                        'If there are any error messages, please save your work and reload MATLAB; ',...
                        'otherwise you should be able to continue, with caution...']);
                end
            else
                obj.meltsIndex = 0;
            end

            obj.getSystemNamesAndWeights;
            obj.fields = {'g'; 'h'; 's'; 'v'; 'cp'; 'dcpdt'; 'dvdt'; 'dvdp'; 'd2vdt2'; 'd2vdtdp'; 'd2vdp2'; 'molwt'; 'rho'; 'mass'};
            obj.loaded = 1;

        end

        function getSystemNamesAndWeights(obj)

            nCharInName = 20;
            numberPhases = 100; phaseIndices = zeros(1, 100);
            phaseNames = blanks(numberPhases*nCharInName);

            [~, phaseNames, ~, numberPhases]  = calllib('libalphamelts', 'getMeltsPhaseNames',...
                true, phaseNames, nCharInName, numberPhases, phaseIndices);

            obj.phases = split(string(phaseNames));
            obj.phases(numberPhases+1:end) = []; % 0x0 char
            obj.phases = obj.phases';

            nCharInName = 20;
            numberOxides = 50;
            oxideNames = blanks(numberOxides*nCharInName);
            oxideWeights(1:numberOxides) = 0;

            [~, oxideNames, ~, ~]  = calllib('libalphamelts', 'getMeltsOxideNames',...
                true, oxideNames, nCharInName, numberOxides);

            [~, oxideWeights, numberOxides]  = calllib('libalphamelts', 'getMeltsOxideWeights',...
                true, oxideWeights, numberOxides);

            oxideNames = strsplit(string(oxideNames));
            oxideNames(numberOxides+1:end) = []; % ""
            oxideWeights(numberOxides+1:end) = [];
            numberOxides = numberOxides - 1;

            obj.endmembers = MELTSmap();
            obj.endmembers(obj.phases(1)) = oxideNames(1:numberOxides)'; % system
            obj.endmembers(obj.phases(2)) = oxideNames(numberOxides+1); % oxygen

            obj.molwts = MELTSmap();
            obj.molwts(obj.phases(1)) = oxideWeights(1:numberOxides)'; % system
            obj.molwts(obj.phases(2)) = oxideWeights(numberOxides+1); % oxygen

            for i = 3:length(obj.phases)
                numberEndMembers = 50;
                endMemberFormulas = blanks(numberEndMembers*nCharInName);
                endMemberWeights(1:numberEndMembers) = 0;

                [~, ~, endMemberWeights, endMemberFormulas, ~, numberEndMembers] = calllib('libalphamelts', 'getMeltsWeightsAndFormulas',...
                    true, char(obj.phases(i)), endMemberWeights, endMemberFormulas, nCharInName, numberEndMembers);

                endMemberFormulas = strsplit(string(endMemberFormulas));
                endMemberFormulas(numberEndMembers+1:end) = []; % ""
                endMemberWeights(numberEndMembers+1:end) = [];
                obj.endmembers(obj.phases(i)) = endMemberFormulas';
                obj.molwts(obj.phases(i)) = endMemberWeights';
            end

        end

        function cMode = getCalculationMode(obj)

            if ~isempty(obj.calculationMode)
                cMode = calllib('libalphamelts', 'getCalculationMode');
            else
                cMode = -1;
            end

        end

        function success = setCalculationMode(obj, cMode)

            switch cMode
                case 1
                    disp ('Setting calculation mode to rhyolite-MELTS 1.0.2.');
                case 2
                    disp ('Setting calculation mode to pMELTS');
                case 3
                    disp ('Setting calculation mode to rhyolite-MELTS 1.1.0.');
                case 4
                    disp ('Setting calculation mode to rhyolite-MELTS 1.2.0.');
                otherwise
                    cMode = max(1, cMode);
                    warning ('Unexpected value for calculationMode. Using default (rhyolite-MELTS 1.0.2)');
            end

            unloadlibrary('libalphamelts');
            obj.loaded = false;

            loadlibrary('libalphamelts', 'MELTSdynamic.h', 'mfilename', 'libalphamelts_proto.m');
            if obj.console; calllib('libalphamelts', 'addConsole'); end

            obj.calculationMode = cMode;
            success = calllib('libalphamelts', 'setCalculationMode', obj.calculationMode);
            cMode = calllib('libalphamelts', 'getCalculationMode');
            if ~success
                if cMode ~= obj.calculationMode
                    error(['Could not reset MELTS library calculation mode from %d to %d! ',...
                        'Please save your work and reload MATLAB.'], cMode, obj.calculationMode);
                else
                    warning('MELTS:libraryAlreadyInitialized', ['Could not re-initialize MELTS library! ',...
                        'Please check stderr (console or terminal): ',...
                        'if there are any error messages, please save your work and reload MATLAB; ',...
                        'otherwise you should be able to continue, with caution...']);
                end
            else
                disp('MELTS library has been reloaded and initialized.');
                obj.getSystemNamesAndWeights;
                obj.meltsIndex = 0;
            end
            obj.loaded = true;

        end

        function success = reload(obj, func)

            disp(['MELTS call ', func, ' failed. Cleaning up...'])
            if ~obj.console; calllib('libalphamelts', 'closeConsole'); end
            unloadlibrary('libalphamelts');
            obj.loaded = false;

            loadlibrary('libalphamelts', 'MELTSdynamic.h', 'mfilename', 'libalphamelts_proto.m');
            if obj.console; calllib('libalphamelts', 'addConsole'); end
            success = calllib('libalphamelts', 'setCalculationMode', obj.calculationMode);

            if ~success
                error(['Could not re-initialize MELTS library after failure! ',...
                        'Please check stderr (console or terminal) for failure details, ',...
                        'then save your work and reload MATLAB.']);
            else
                warning('MELTS:calculationFailed', ['MELTS library has been reloaded and initialized. ',...
                    'Please check stderr (console or terminal) for failure details.'])
                obj.meltsIndex = 0;
            end
            obj.loaded = true;

        end

    end

end
