classdef MELTSengine < matlab.mixin.SetGet & matlab.mixin.Copyable
    % Thermodynamic 'engine' where the real work happens.
    % Calculation mode is set in parent MELTSdyanmic list: 1. rhyolite-MELTS 1.0.2; 2. pMELTS; 3. rhyolite-MELTS 1.1.0; 4. rhyolite-MELTS 1.2.0.
    % For equilibration, run mode need to be one of: 1. isobaric, isothermal; 2. isenthalpic; 3. isentropic; 4. isochoric.
    % Output is post-equilibration / pre-fractionation (1) by default, but calculation can also be batch, no fractionation (0), or output can be after fractionation (2).
    % phaseNames format is: [phase name in MELTS system][optional integer, like equilibration output][optional underscore, '_'][optional user-defined name, comprising letters and digits]

    properties (SetAccess = protected)
        % MELTSstatus class used to track status of the C library of MELTS functions.
        status
        % MELTS calculation mode inherited from MELTSdynamic list.
        calculationMode = []
    end

    properties
        % User-defined index for this node of the parent MELTSdynamic list.
        nodeId
        % User-defined name and/or description for this node of the parent MELTSdynamic list.
        nodeName
        % Tracks the 'index' in any text output files (may not work for complex calculations where library is reloaded or MELTS dynamic list is adjusted).
        meltsIndex
        % Controls whether any constraints are applied in equilibration calculations: 1. none (default); 2. isenthalpic; 3. isentropic; 4. isochoric.
        runMode
        % INPUT and OUTPUT: pressure (bars)
        pressure
        % INPUT and OUTPUT: temperature (Celsius)
        temperature
        % INPUT and OUTPUT: Reference quantity for isenthalpic (J), isentropic (J/k), or isochoric (cc) calculations.
        reference
        % OUTPUT: oxygen fugacity (log 10 units) of the system.
        logfO2
        % String array of system properties that resemble lines in .melts files.
        systemProperties
        %outputFlag % put in systemProperties?

        % INPUT and OUTPUT: composition in grams (or loosely wt%) for regular MELTS calculations; updated after any fractionations regardless of output mode.
        bulkComposition

        % OUTPUT: string array of liquid phase name(s) from regular MELTS calculations (multiple instanced count from 1); may be INPUT as phaseNames to 'Supplemental calculator'.
        liquidNames
        % OUTPUT: string array of solid phase names from regular MELTS calculations (multiple instanced count from 1); may be INPUT as phaseNames to 'Supplemental calculator'.
        solidNames

        % INPUT: string array of phase names for 'Supplemental calculator' calculations (currently also output for calcSaturationState).
        phaseNames
        % INPUT: composition in grams (or loosely wt%) for 'Supplemental calculator' calculations.
        phaseComposition
        % INPUT: composition in end-member mol frac (or oxide mol frac) for 'Supplemental calculator' calculations.
        molarComposition

        % OUTPUT: MELTSmap of Gibbs free energy (J) for solid and liquid phases and/or 'bulk'.
        g = MELTSmap()
        % OUTPUT: MELTSmap of enthalpy (J) for solid and liquid phases and/or 'bulk'.
        h = MELTSmap()
        % OUTPUT: MELTSmap of entropy (J/K) for solid and liquid phases and/or 'bulk'.
        s = MELTSmap()
        % OUTPUT: MELTSmap of volume (cc) for solid and liquid phases and/or 'bulk'.
        v = MELTSmap()
        % OUTPUT: MELTSmap of heat capacity (J/K) for solid and liquid phases and/or 'bulk'.
        cp = MELTSmap()
        % OUTPUT: MELTSmap of partial derivative of heat capacity wrt temperature (J/K/K) for solid and liquid phases and/or 'bulk'.
        dcpdt = MELTSmap()
        % OUTPUT: MELTSmap of partial derivative of volume wrt temperature (cc/K) for solid and liquid phases and/or 'bulk'.
        dvdt = MELTSmap()
        % OUTPUT: MELTSmap of partial derivative of volume wrt pressure (cc/bar) for solid and liquid phases and/or 'bulk'.
        dvdp = MELTSmap()
        % OUTPUT: MELTSmap of 2nd partial derivative of volume wrt temperature (cc/K/K) for solid and liquid phases and/or 'bulk'.
        d2vdt2 = MELTSmap()
        % OUTPUT: MELTSmap of 2nd partial derivative of volume wrt temperature, pressure (cc/K/bar) for solid and liquid phases and/or 'bulk'.
        d2vdtdp = MELTSmap()
        % OUTPUT: MELTSmap of 2nd partial derivative of volume wrt pressure (cc/bar/bar) for solid and liquid phases and/or 'bulk'.
        d2vdp2 = MELTSmap()
        % OUTPUT: MELTSmap of molecular weights of phases from regular MELTS calculations; use to convert to/from moles.
        molwt = MELTSmap()
        % OUTPUT: MELTSmap of density (g/cc) for solid and liquid phases and/or 'bulk'.
        rho = MELTSmap()
        % OUTPUT: MELTSmap of mass (grams) for solid and liquid phases and/or 'bulk'.
        mass = MELTSmap()
        % OUTPUT: MELTSmap of log 10 viscosity (Poise) for one or more liquid phases, and system viscosity (not including previously fractionated material).
        viscosity = MELTSmap()

        % OUTPUT: MELTSmap of compositions (wt%) for solid and liquid phases and/or 'bulk'; may be INPUT for 'Supplemental Calculator' calculations.
        dispComposition = MELTSmap()

        % OUTPUT: MELTSmap of affinities (J/mol) for solid and liquid phases.
        affinity = MELTSmap()
        % OUTPUT: MELTSmap (mol frac) for solid and liquid phases and/or 'bulk'; may be INPUT for molar properties.
        X = MELTSmap()


        % OUTPUT: Dictionary for end-member activities relative to fixed structural and/or ordering state (for 'true' activities, using mu and mu0 or 'activity' is usually preferable).
        activity0 = MELTSmap()
        % OUTPUT: Dictionary for end-member activities relative to pure phase structural and/or ordering state at the given pressure and temperature (calculated from mu and mu0)
        activity = MELTSmap()
        % OUTPUT: MELTSmap for end-member (or oxide) pure-phase chemical potentials (J/mol) according to Berman/MELTS thermodynamic database, for solid and liquid phases and/or 'bulk'.
        mu0 = MELTSmap()
        % OUTPUT: MELTSmap for end-member (or oxide) chemical potentials (J/mol) for solid and liquid phases and/or 'bulk'.
        mu = MELTSmap()
    end

    methods

        function obj = MELTSengine(cMode, name)
            % Called by MELTSdynamic to set up the thermodynamic 'engine' where all the real work happens.
            if nargin == 0; cMode = []; end
            if nargin > 0 && isa(cMode, 'MELTSstatus')
                obj.status = cMode;
                obj.calculationMode = obj.status.calculationMode;
            elseif ~isempty(cMode)
                obj.calculationMode = cMode;
                obj.status = MELTSstatus(cMode);
            end

            if nargin > 1; obj.nodeName = name; end

        end

        function setBulkComposition(obj, oxide, val)
            % Set the bulk compsoitions for a particular oxide (in grams), or the entire compositional array.
            if obj.calculationMode ~= obj.status.getCalculationMode
                success = obj.status.setCalculationMode(obj.calculationMode);
                if success; obj.setSystemProperties(obj.systemProperties); end
            end

            if isempty(obj.bulkComposition)
                obj.bulkComposition = zeros(size(obj.status.endmembers("bulk")));
            end
            if isnumeric(oxide)
                assert(length(oxide) == length(obj.bulkComposition), "Array sizes do not match");
                obj.bulkComposition(:) = oxide(:);
            elseif nargin > 1 && ~isempty(val)
                iox = find(strcmpi(obj.status.endmembers("bulk"), oxide));
                if iox; obj.bulkComposition(iox) = val;
                else; disp(['Oxide ', oxide, ' not found. Ignoring value.']); end
            end

        end

        function setPhaseComposition(obj, oxide, val)
            % For a given phase set either a particular oxide (in grams), or the entire compositional array.
            if obj.calculationMode ~= obj.status.getCalculationMode
                success = obj.status.setCalculationMode(obj.calculationMode);
                if success; obj.setSystemProperties(obj.systemProperties); end
            end

            if isempty(obj.phaseComposition)
                obj.phaseComposition = zeros(size(obj.status.endmembers("bulk")));
            end
            if isnumeric(oxide)
                assert(length(oxide) == length(obj.phaseComposition), "Array sizes do not match");
                obj.phaseComposition(:) = oxide(:);
            elseif nargin > 1 && ~isempty(val)
                iox = find(strcmpi(obj.status.endmembers("bulk"), oxide));
                if iox; obj.phaseComposition(iox) = val;
                else; disp(['Oxide ', oxide, ' not found. Ignoring value.']); end
            end

        end

        function setInitialComposition(obj, oxide, val)
            % Set both bulk composition and phase composition to the passed value(s).
            if nargin == 2; val = []; end
            obj.setBulkComposition(oxide, val);
            obj.setPhaseComposition(oxide, val);

        end

        function setMolarComposition(obj, varargin)
            % For a given phase set either a particular end-member mole frac (or oxide for 'bulk'), or the entire mole fraction array.
            if obj.calculationMode ~= obj.status.getCalculationMode
                success = obj.status.setCalculationMode(obj.calculationMode);
                if success; obj.setSystemProperties(obj.systemProperties); end
            end

            nox = length(obj.status.endmembers("bulk"));
            if isempty(obj.molarComposition)
                obj.molarComposition = zeros(size(obj.status.endmembers("bulk")));
            end
            if nargin > 1 && isnumeric(varargin{1})
                len = length(varargin{1});
                assert(len <= length(obj.molarComposition), "Array sizes do not match");
                obj.molarComposition(1:len) = varargin{1}(:);
                obj.molarComposition(len+1:nox) = 0.0;
            elseif nargin > 3
                iox = find(strcmpi(obj.status.endmembers(varargin{1}), varargin{2}));
                if iox; obj.molarComposition(iox) = varargin{3};
                else; disp(['Endmember ', varargin{2}, ' not found. Ignoring value.']); end
            end

        end

        function value = getProperty(obj, propertyName, varargin)
            % Get some property (for system, or one or more phases) from engine. Returns a scalar, vector or 2-D matrix.

            % Assume endmember or oxide present in all
            phaseList = [];
            if nargin > 2
                if strcmp(propertyName, 'bulkComposition') || strcmp(propertyName, 'phaseComposition') || ...
                        strcmp(propertyName, 'molarComposition')
                    value = obj.get(propertyName);
                    if isempty(value); value = NaN(size(obj.status.endmembers("bulk"))); end
                    endMemberName = varargin{:};
                else % MELTSmap

                    phaseList = varargin{1};

                    if nargin > 3
                        endMemberName = varargin{2};
                        if strcmp(propertyName, 'dispComposition')
                            phaseName = "bulk";
                        else
                            phaseName = regexprep(char(phaseList(1)), '(_+.*)|(\d*)', '');
                        end
                        value = NaN(size(obj.status.endmembers(phaseName)));
                    else
                        if strcmp(propertyName, 'dispComposition')
                            endMemberName = true(size(obj.status.endmembers("bulk")));
                            value = NaN(size(obj.status.endmembers("bulk")));
                        elseif strcmp(propertyName, 'X') || strcmp(propertyName, 'activity') || strcmp(propertyName, 'activity0') || ...
                            strcmp(propertyName, 'mu') || strcmp(propertyName, 'mu0')
                            phaseName = regexprep(char(phaseList(1)), '(_+.*)|(\d*)', '');
                            assert(all(cell2mat(regexp(phaseList, phaseName))), ...
                                "The getProperty method cannot be called for multiple phases for that property, as the number of endmembers may vary.");
                            endMemberName = true(size(obj.status.endmembers(phaseName)));
                            value = NaN(size(obj.status.endmembers(phaseName)));
                        else
                            endMemberName = 1;
                            value = NaN;
                        end
                    end

                    if ~isstring(phaseList); phaseList = string(phaseList); end
                    if isKey(obj.(propertyName), char(phaseList(1))); value = obj.(propertyName)(phaseList); end

                end

                if isnumeric(endMemberName) || islogical(endMemberName)
                    indices = endMemberName;
                else
                    indices = find(strcmp(obj.status.endmembers(phaseName), endMemberName));
                end

                if length(phaseList) > 1
                    value = value(indices, :);
                else
                    value = value(indices);
                end

            else
                value = obj.get(propertyName);
                if isempty(value); value = NaN; end
            end

            assert(isnumeric(obj.status.nullValue), "Cannot replace NaN with non-numeric value.")
            if ~isnan(obj.status.nullValue)
                value(isnan(value)) = obj.status.nullValue;
            end

        end

        function setSystemProperties(obj, varargin)
            % Set system properties using string arrays or string pairs that resemble lines in .melts files.
            nCharInString = 132; % REC is 134 and the last two characters will be used for ' \0'

            if obj.calculationMode ~= obj.status.getCalculationMode
                success = obj.status.setCalculationMode(obj.calculationMode);
                if success; obj.setSystemProperties(obj.systemProperties); end
            end

            if nargin > 2
                assert(logical(mod(nargin, 2)), "Number of values does not match number of properties?");
                lines = reshape(varargin, [], 2);
            elseif nargin == 2
                lines = varargin{:};
            end

            if ~isstring(lines)
                lines = string(lines);
            end
            if size(lines, 2) > 1
                if size(lines, 2) > 2
                    delimiter = [": "; string(reshape(blanks(size(lines, 2)-2), [], 1))]';
                else
                    delimiter = ": ";
                end
                lines = join(lines, repmat(delimiter, size(lines, 1), 1), 2);
            end
            properties = pad(lines, nCharInString);
            properties = char(join(properties, ''));

            funcData = MELTSevent('setSystemProperties');
            checkSuccess = onCleanup(@() notify(obj.status, 'finished', funcData));

            failure = true;
            obj.status.set('failed', true);
            failure = calllib('libalphamelts', 'setMeltsSystemProperties', failure,...
                properties, nCharInString, length(lines));
            obj.status.set('failed', failure);

            if ~obj.status.failed
                for i = 1:length(lines)
                    if ~any(strcmp(lines(i), obj.systemProperties))
                        obj.systemProperties = [obj.systemProperties; lines(i)];
                    end
                end
            end

        end

        % Eventually have an alphaMELTS-like find boundary ('isograd') too
        function findLiquidus(obj, varargin)
            % Find liquidus for current or passed bulk composition. Gets affinities for solid phases that are not suppressed.

            obj.calcEquilibriumState(0, 0, varargin{:});
            obj.phaseNames = obj.calcSaturationState;

            % There may be 'solid' water or fluid
            affinities = obj.getProperty('affinity', obj.phaseNames);
            liquidusPhase = obj.phaseNames(affinities == min(affinities(affinities ~= 0), [], 'omitnan'));
            obj.solidNames = [obj.solidNames liquidusPhase+"1"];
            obj.mass(liquidusPhase+"1") = 0.0;
            obj.dispComposition(liquidusPhase+"1") = obj.dispComposition(liquidusPhase);

        end

        function calcEquilibriumState(obj, runMode, outputFlag, name, varargin)
            % Single equilibration calculation for current or passed bulk composition with or without constraints (equivalent to alphaMELTS menu option 3 or 4). By default output is post-equilibration / pre-fractionation.
            if obj.calculationMode ~= obj.status.getCalculationMode
                success = obj.status.setCalculationMode(obj.calculationMode);
                if success; obj.setSystemProperties(obj.systemProperties); end
            end

            if (nargin < 2) || isempty(runMode); runMode = obj.runMode; end
            if (nargin < 3) || isempty(outputFlag); outputFlag = 1; end
            if (nargin > 3) && ~isempty(name); obj.nodeName = name; end
            if (nargin > 4); obj.setBulkComposition(varargin{:}); end

            if isempty(runMode); runMode = 1; end
            if isempty(obj.reference); obj.reference = 0.0; end

            % reset runMode if this is not a findLiquidus call
            if runMode; obj.runMode = runMode; end
            assert(any(runMode == 0:4), "Unexpected run mode?");
            assert(any(obj.temperature > 0), "System temperature not set?");
            if obj.pressure == 0; obj.pressure = 1; end
            assert(any(obj.pressure > 0), "System pressure not set?");
            assert(any(obj.bulkComposition), "System bulk composition not set?");

            nCharInName = 20; numberPhases = 20;  phaseIndices = zeros(1, 20);
            phases = blanks(numberPhases*nCharInName);
            nCharInString = 100;
            errorString = blanks(nCharInString);

            nox = length(obj.status.endmembers('bulk'));
            nf = length(obj.status.fields) - 3;
            nc = nf + nox + 3;
            properties(nc*numberPhases) = 0;

            funcData = MELTSevent('calcEquilibriumState');
            checkSuccess = onCleanup(@() notify(obj.status, 'finished', funcData));

            failure = true;
            obj.status.set('failed', failure);

            [failure, ~, press, bulk, enthalpy, temp, phases, ~,...
                numberPhases, ~, errorString, ~, properties, phaseIndices] = ...
                calllib('libalphamelts', 'driveMeltsProcess', failure, runMode,...
                obj.pressure, obj.bulkComposition, obj.reference, obj.temperature+273.15,...
                phases, nCharInName, numberPhases, outputFlag, errorString, nCharInString, properties, phaseIndices);
            obj.status.set('failed', failure);
            obj.status.set('message', errorString);

            % Engine output is reset
            obj.g = MELTSmap();
            obj.h = MELTSmap();
            obj.s = MELTSmap();
            obj.v = MELTSmap();
            obj.cp = MELTSmap();
            obj.dcpdt = MELTSmap();
            obj.dvdt = MELTSmap();
            obj.dvdp = MELTSmap();
            obj.d2vdt2 = MELTSmap();
            obj.d2vdtdp = MELTSmap();
            obj.d2vdp2 = MELTSmap();
            obj.molwt = MELTSmap();
            obj.rho = MELTSmap();
            obj.mass = MELTSmap();
            obj.viscosity = MELTSmap();
            obj.dispComposition = MELTSmap();
            obj.affinity = MELTSmap();
            obj.X = MELTSmap();
            obj.activity = MELTSmap();
            obj.mu0 = MELTSmap();
            obj.mu = MELTSmap();

            if ~obj.status.failed

                % Find Liquidus does not generate file output, but the rest do
                if runMode
                    obj.status.meltsIndex = obj.status.meltsIndex + 1;
                    obj.meltsIndex = obj.status.meltsIndex;
                end

                % Get and process the phaseIndices array
                phaseIndices(numberPhases + 1:end) = [];
                phaseIndices = mod(phaseIndices, 10) + 1;

                % These values are post-fractionation
                obj.pressure = press;
                obj.reference = enthalpy; % could actually be entropy, enthalpy or volume
                obj.temperature = temp-273.15;
                obj.bulkComposition = bulk;

                phases = split(string(phases));
                phases(numberPhases+1:end) = [];
                phases = phases';

                % These values are post-equilibration / pre-fractionation if outputFlag=1;
                % if outputFlag=2 then everything is post-fractionation (may be needed for bookkeeping);
                properties(nc*numberPhases + 1:end) = [];
                properties = reshape(properties, nc, numberPhases);

                obj.dispComposition("bulk") = 100.0*properties(nf+1:nf+nox, 1) / sum(properties(nf+1:nf+nox, 1));
                phaseProps = properties([1:nf nc-2:nc], 1);
                obj.viscosity("bulk") = phaseProps(end); % borrow the total grams slot
                phaseProps(nf+3) = phaseProps(nf+2) * phaseProps(4);
                for j = 1:nf+3; obj.(obj.status.fields{j})("bulk") = phaseProps(j); end

                %obj.dispComposition("oxygen") = 100.0*properties(nf+1:nf+nox, 1) / sum(properties(nf+1:nf+nox, 1));
                phaseProps = properties([1:nf nc-2:nc], 2);
                obj.logfO2 = phaseProps(end-2); % borrow the mw slot
                phaseProps(nf+1) = 31.998;
                for j = 1:nf+3; obj.(obj.status.fields{j})("oxygen") = phaseProps(j); end

                % This can be made to handle multiple liquids and/or alloy-liquids etc.
                hasLiquid = find(strcmpi(phases, 'liquid'));
                hasSolids = find(~strcmpi(phases, 'liquid')); hasSolids(1:2) = []; % system and oxygen

                obj.liquidNames = [];
                obj.solidNames = [];

                if ~isempty(hasLiquid)
                    obj.liquidNames = phases(hasLiquid) + phaseIndices(hasLiquid);
                    obj.dispComposition(obj.liquidNames) = 100.0*properties(nf+1:nf+nox, hasLiquid) ./ sum(properties(nf+1:nf+nox, hasLiquid));
                    phaseProps = properties([1:nf nc-2:nc], hasLiquid);
                    obj.viscosity(obj.liquidNames) = phaseProps(nf+3, :); % borrow the total grams slot
                    phaseProps(nf+3, :) = phaseProps(nf+2, :) .* phaseProps(4, :);
                    for j = 1:nf+3;  obj.(obj.status.fields{j})(obj.liquidNames) = phaseProps(j, :); end
                end

                if ~isempty(hasSolids)
                    obj.solidNames = phases(hasSolids) + phaseIndices(hasSolids);
                    obj.dispComposition(obj.solidNames) =...
                        100.0*properties(nf+1:nf+nox, hasSolids) ./ sum(properties(nf+1:nf+nox, hasSolids));
                    phaseProps = properties([1:nf nc-2:nc], hasSolids);
                    for j = 1:nf+3; obj.(obj.status.fields{j})(obj.solidNames) = phaseProps(j, :); end
                end

            end

        end

        function calcPhaseProperties(obj, phaseList, varargin)
            % 'Supplemental Calculator' type calculation to get thermodynamic properties of one or more phases.
            % If composition(s) not passed, will use the contents of phaseComposition (grams).
            if obj.calculationMode ~= obj.status.getCalculationMode
                success = obj.status.setCalculationMode(obj.calculationMode);
                if success; obj.setSystemProperties(obj.systemProperties); end
            end

            assert(any(obj.temperature > 0), "System temperature not set?");
            if obj.pressure == 0; obj.pressure = 1; end
            assert(any(obj.pressure > 0), "System pressure not set?");

            if (nargin < 2) || isempty(phaseList); phaseList = "liquid";
            elseif ~isstring(phaseList); phaseList = string(phaseList); end

            nox = length(obj.status.endmembers('bulk'));
            nf = length(obj.status.fields) - 3;
            nc = nf + nox + 3;

            if length(phaseList) > 1
                assert(isrow(phaseList) && isnumeric(varargin{1}) && all(size(varargin{1}) == [nox length(phaseList)]), ...
                    "Please supply a row of phase names, and a matrix with compositions as columns")
            end

            obj.phaseNames = [obj.phaseNames phaseList];

            for i = 1:length(phaseList)

                if nargin > 2
                    if length(phaseList) > 1
                        % if all nan then continue - not yet implemented
                        obj.setPhaseComposition(varargin{1}(:,i));
                    else
                        % varargin{1} may be an oxide name or composition
                        obj.setPhaseComposition(varargin{:});
                    end
                end

                assert(any(obj.phaseComposition), "Phase composition not set?");

                phaseName = char(phaseList(i));
                phaseProps(1:nc) = 0;
                visc = 0;

                funcData = MELTSevent('calcPhaseProperties');
                checkSuccess = onCleanup(@() notify(obj.status, 'finished', funcData));

                failure = true;
                obj.status.set('failed', failure);

                if startsWith(phaseName, 'bulk') || startsWith(phaseName, 'liquid')

                    [failure, ~, ~, ~, ~, phaseProps] = calllib('libalphamelts', 'getMeltsPhaseProperties',...
                    true, 'liquid', obj.temperature+273.15, obj.pressure, obj.phaseComposition, phaseProps);

                    if ~failure
                        [failure, ~, ~, ~, visc] = calllib('libalphamelts', 'getMeltsViscosity',...
                        true, 'Shaw', obj.temperature+273.15, obj.phaseComposition, visc);
                    end

                else

                    [failure, ~, ~, ~, ~, phaseProps] = calllib('libalphamelts', 'getMeltsPhaseProperties',...
                    true, regexprep(phaseName, '(_+.*)|(\d*)', ''), obj.temperature+273.15, obj.pressure, obj.phaseComposition, phaseProps);

                end

                obj.status.set('failed', failure);

                if ~obj.status.failed
                    obj.dispComposition(phaseName) = 100.0*phaseProps(nf+1:nf+nox) / sum(phaseProps(nf+1:nf+nox));
                    phaseProps = phaseProps([1:nf nc-2:nc]);
                    for j = 1:nf+3; obj.(obj.status.fields{j})(phaseName) = phaseProps(j); end
                    if visc ~= 0; obj.viscosity(phaseName) = visc; end
                end

            end

        end

        function calcMolarProperties(obj, phaseList, varargin)
            % 'Supplemental Calculator' type calculation to get thermodynamic properties of one or more phases, or of oxygen.
            % If composition(s) not passed, will use the contents of molarComposition (mol frac).
            if obj.calculationMode ~= obj.status.getCalculationMode
                success = obj.status.setCalculationMode(obj.calculationMode);
                if success; obj.setSystemProperties(obj.systemProperties); end
            end

            assert(any(obj.temperature > 0), "System temperature not set?");
            if obj.pressure == 0; obj.pressure = 1; end
            assert(any(obj.pressure > 0), "System pressure not set?");

            if (nargin < 2) || isempty(phaseList); phaseList = "liquid";
            elseif ~isstring(phaseList); phaseList = string(phaseList); end

            nox = length(obj.status.endmembers('bulk'));
            nf = length(obj.status.fields) - 3;
            nc = nf + nox + 3;

            if length(phaseList) > 1
                assert(isrow(phaseList) && isnumeric(varargin{1}) && ...
                    (size(varargin{1}, 1) <= nox) && (size(varargin{1}, 2) == length(phaseList)), ...
                    "Please supply a row of phase names, and a matrix with compositions as columns")
            end

            obj.phaseNames = [obj.phaseNames phaseList];

            for i = 1:length(phaseList)

                if nargin > 2
                    if length(phaseList) > 1
                        obj.setMolarComposition(varargin{1}(:,i));
                    else
                        obj.setMolarComposition(phaseList, varargin{:});
                    end
                end

                assert(any(obj.molarComposition), "Phase composition not set?");

                phaseName = char(phaseList(i));
                phaseProps(1:nc) = 0;
                visc = 0;

                funcData = MELTSevent('calcMolarProperties');
                checkSuccess = onCleanup(@() notify(obj.status, 'finished', funcData));

                failure = true;
                obj.status.set('failed', failure);

                if startsWith(phaseName, 'bulk') || startsWith(phaseName, 'liquid')

                    [failure, ~, ~, ~, ~, phaseProps] = calllib('libalphamelts', 'getMeltsMolarProperties',...
                    true, 'liquid', obj.temperature+273.15, obj.pressure, obj.molarComposition, phaseProps);

                    if ~failure
                        obj.phaseComposition = 100.0*phaseProps(nf+1:nf+nox) / sum(phaseProps(nf+1:nf+nox));
                        [failure, ~, ~, ~, visc] = calllib('libalphamelts', 'getMeltsViscosity',...
                        true, 'Shaw', obj.temperature+273.15, obj.phaseComposition, visc);
                    end

                else

                    [failure, ~, ~, ~, ~, phaseProps] = calllib('libalphamelts', 'getMeltsMolarProperties',...
                    true, regexprep(phaseName, '(_+.*)|(\d*)', ''), obj.temperature+273.15, obj.pressure, obj.molarComposition, phaseProps);

                end

                obj.status.set('failed', failure);

                if ~obj.status.failed
                    obj.dispComposition(phaseName) = 100.0*phaseProps(nf+1:nf+nox) / sum(phaseProps(nf+1:nf+nox));
                    phaseProps = phaseProps([1:nf nc-2:nc]);
                    for j = 1:nf+3; obj.(obj.status.fields{j})(phaseName) = phaseProps(j); end
                    if visc ~= 0; obj.viscosity(phaseName) = visc; end
                end

            end

        end

        function calcViscosityFromGRD(obj, phaseList, varargin)
            % 'Supplemental Calculator' type calculation to get viscosity of one or more liquids. If composition(s) not passed, will use the contents of phaseComposition (grams).
            % System ('bulk') viscosity is not updated by this method; it always uses Shaw model and crystallinity.
            if obj.calculationMode ~= obj.status.getCalculationMode
                success = obj.status.setCalculationMode(obj.calculationMode);
                if success; obj.setSystemProperties(obj.systemProperties); end
            end
            assert(any(obj.temperature > 0), "System temperature not set?");

            if (nargin < 2) || isempty(phaseList); phaseList = "liquid";
            elseif ~isstring(phaseList); phaseList = string(phaseList); end

            nox = length(obj.status.endmembers('bulk'));

            if length(phaseList) > 1
                assert(isrow(phaseList) && isnumeric(varargin{1}) && all(size(varargin{1}) == [nox length(phaseList)]), ...
                    "Please supply a row of phase names, and a matrix with compositions as columns")
            end

            obj.phaseNames = [obj.phaseNames phaseList];

            for i = 1:length(phaseList)

                if nargin > 2
                    if length(phaseList) > 1
                        obj.setPhaseComposition(varargin{1}(:,i));
                    else
                        % varargin{1} may be an oxide name or composition
                        obj.setPhaseComposition(varargin{:});
                    end
                end

                assert(any(obj.phaseComposition), "Phase composition not set?");

                phaseName = char(phaseList(i));
                visc = 0;

                assert((startsWith(phaseName, 'bulk') || startsWith(phaseName, 'liquid')), "Phase is not liquid in viscosity calc?");

                funcData = MELTSevent('calcViscosityFromGRD');
                checkSuccess = onCleanup(@() notify(obj.status, 'finished', funcData));

                failure = true;
                obj.status.set('failed', failure);

                [failure, ~, ~, ~, visc] = calllib('libalphamelts', 'getMeltsViscosity',...
                true, 'GRD', obj.temperature+273.15, obj.phaseComposition, visc);

                obj.status.set('failed', failure);

                if ~obj.status.failed
                    if visc ~= 0; obj.viscosity(phaseName) = log10(10.0 * visc); end % PaS -> log10 Poise
                end

            end

        end

        function calcEndMemberProperties(obj, phaseList, varargin)
            % 'Supplemental Calculator' type calculation to get end-member thermodynamic properties of one or more phases.
            % If composition(s) not passed, will use the contents of phaseComposition (grams).
            % If there is an fO2 buffer, liquid composition will be updated according Kress & Carmichael (1991).
            if obj.calculationMode ~= obj.status.getCalculationMode
                success = obj.status.setCalculationMode(obj.calculationMode);
                if success; obj.setSystemProperties(obj.systemProperties); end
            end

            assert(any(obj.temperature > 0), "System temperature not set?");
            if obj.pressure == 0; obj.pressure = 1; end
            assert(any(obj.pressure > 0), "System pressure not set?");

            if (nargin < 2) || isempty(phaseList); phaseList = "bulk";
            elseif ~isstring(phaseList); phaseList = string(phaseList); end

            nox = length(obj.status.endmembers('bulk'));
            fields = {'X'; 'act'; 'mu0'; 'mu'};
            nc = length(fields);

            funcData = MELTSevent('calcEndMemberProperties');
            checkSuccess = onCleanup(@() notify(obj.status, 'finished', funcData));

            nCharInName = 30; numberEndMembers = 20;
            endMemberFormulas = blanks(numberEndMembers*nCharInName);

            if length(phaseList) > 1
                assert(isrow(phaseList) && isnumeric(varargin{1}) && all(size(varargin{1}) == [nox length(phaseList)]), ...
                    "Please supply a row of phase names, and a matrix with compositions as columns")
            end

            obj.phaseNames = [obj.phaseNames phaseList];

            for i = 1:length(phaseList)

                if nargin > 2
                    if length(phaseList) > 1
                        obj.setPhaseComposition(varargin{1}(:,i));
                    else
                        % varargin{1} may be an oxide name or composition
                        obj.setPhaseComposition(varargin{:});
                    end
                end

                assert(any(obj.phaseComposition), "Phase composition not set?");

                failure = true;
                obj.status.set('failed', failure);

                phaseName = char(phaseList(i)); numberEndMembers = 20;
                properties(1:nc*numberEndMembers) = 0;

                if startsWith(phaseName, 'bulk')

                    [failure, ~, ~, ~, bulk, ~, ~, numberEndMembers, properties] = calllib('libalphamelts',...
                        'getMeltsOxideProperties', failure, 'liquid', ...
                        obj.temperature+273.15, obj.pressure, obj.phaseComposition,...
                        endMemberFormulas, nCharInName, numberEndMembers, properties);

                else

                    [failure, ~, ~, ~, bulk, ~, ~, numberEndMembers, properties] = calllib('libalphamelts',...
                        'getMeltsEndMemberProperties', failure, regexprep(phaseName, '(_+.*)|(\d*)', ''), ...
                        obj.temperature+273.15, obj.pressure, obj.phaseComposition,...
                        endMemberFormulas, nCharInName, numberEndMembers, properties);

                end

                obj.status.set('failed', failure);

                if ~obj.status.failed
                    properties(nc*numberEndMembers + 1:end) = [];
                    endMemberProps = reshape(properties, nc, numberEndMembers);
                    endMemberProps = endMemberProps';

                    % Set undefined mu, mu0 (and activities) to NaN
                    endMemberProps(endMemberProps(:,2) == 0.0 & endMemberProps(:,4) ~= 0.0, 2) = NaN;
                    endMemberProps(endMemberProps(:,3) == 0.0 & endMemberProps(:,4) ~= 0.0, 3) = NaN;
                    endMemberProps(endMemberProps(:,4) == 0.0, 2:4) = NaN;

                    % If phase is liquid and there is an fO2 buffer composition may have changed
                    obj.dispComposition(phaseName) = 100.0*bulk / sum(bulk);

                    obj.X(string(phaseName)) = endMemberProps(:, 1);
                    obj.activity0(string(phaseName)) = endMemberProps(:, 2);
                    obj.activity(string(phaseName)) = exp((endMemberProps(:,4) - endMemberProps(:,3)) / (8.3143*(obj.temperature+273.15)));
                    obj.mu0(string(phaseName)) = endMemberProps(:, 3);
                    obj.mu(string(phaseName)) = endMemberProps(:, 4);
                end

            end

        end

        function phaseList = calcSaturationState(obj, varargin)
            % Get the affinity of each phase not in the assemblage. Can be called before or after equilibration. Also called by findLiquidus.
            % Currently does not update bulkComposition for fO2 buffer; his will change once subsolidus start is implemented.
            if obj.calculationMode ~= obj.status.getCalculationMode
                success = obj.status.setCalculationMode(obj.calculationMode);
                if success; obj.setSystemProperties(obj.systemProperties); end
            end

            if (nargin > 1); obj.setBulkComposition(varargin{:}); end
            assert(any(obj.bulkComposition), "System bulk composition not set?");

            nCharInName = 20; numberPhases = 100; phaseIndices = zeros(1, 100);
            phases = blanks(numberPhases*nCharInName);

            nc = length(obj.status.endmembers('bulk')) + 1;
            properties(1:nc*numberPhases) = 0;

            funcData = MELTSevent('calcSaturationState');
            checkSuccess = onCleanup(@() notify(obj.status, 'finished', funcData));

            failure = true;
            obj.status.set('failed', failure);

            [failure, ~, ~, ~, ~, ~, numberPhases, properties, ~] = ...
                calllib('libalphamelts', 'getMeltsSaturationState', failure, obj.pressure,...
                obj.bulkComposition, obj.temperature+273.15, phases, nCharInName, numberPhases, properties, phaseIndices);
            obj.status.set('failed', failure);

            if ~obj.status.failed
                properties(nc*numberPhases + 1:end) = [];
                properties = reshape(properties, nc, numberPhases);

                phaseList = strings(1, numberPhases);
                j = 0;
                for i = 3:numberPhases % not bulk system or oxygen
                    phaseProps = properties(:, i);
                    if (phaseProps(1)) % not in assemblage
                        j = j+1; phaseName = obj.status.phases(i);
                        phaseList(j) = phaseName;
                        if (phaseProps(1) > -99999) % not failed
                            obj.affinity(phaseName) = phaseProps(1);
                            obj.dispComposition(phaseName) = 100.0*phaseProps(2:end) / sum(phaseProps(2:end));
                        else
                            obj.affinity(phaseName) = NaN;
                        end
                    end
                end
                phaseList(j+1:end) = [];
                obj.phaseNames = [obj.phaseNames phaseList];
            end

        end

        % Will be option to get fO2 by different methods / compare to
        % different buffers in here

        function set(obj, varargin)
            % Extension of set to make sure compositions are stored as columns, and phase lists as rows (so that getListProperty etc. work).
            set@matlab.mixin.SetGet(obj, varargin{:});

            if isrow(obj.bulkComposition); obj.bulkComposition = obj.bulkComposition'; end
            if isrow(obj.phaseComposition); obj.phaseComposition = obj.phaseComposition'; end
            if isrow(obj.molarComposition); obj.molarComposition = obj.molarComposition'; end

            if iscolumn(obj.phaseNames); obj.phaseNames = obj.phaseNames'; end

        end

        function cp = copyAndKeepOutput(obj, varargin)
            % Deep copy of the engine - keeps (duplicates) same engine output as previous engine.
            cp = MELTSengine(varargin{:});
            cp.calculationMode = obj.calculationMode;
            % Engine input is duplicated
            cp.nodeId = obj.nodeId;
            cp.nodeName = obj.nodeName;
            cp.runMode = obj.runMode;
            cp.pressure = obj.pressure;
            cp.temperature = obj.temperature;
            cp.reference = obj.reference;
            cp.systemProperties = obj.systemProperties;
            cp.bulkComposition = obj.bulkComposition;
            cp.phaseComposition = obj.phaseComposition;
            cp.molarComposition = obj.molarComposition;
            % Engine output is duplicated
            cp.g = obj.g;
            cp.h = obj.h;
            cp.s = obj.s;
            cp.v = obj.v;
            cp.cp = obj.cp;
            cp.dcpdt = obj.dcpdt;
            cp.dvdt = obj.dvdt;
            cp.dvdp = obj.dvdp;
            cp.d2vdt2 = obj.d2vdt2;
            cp.d2vdtdp = obj.d2vdtdp;
            cp.d2vdp2 = obj.d2vdp2;
            cp.molwt = obj.molwt;
            cp.rho = obj.rho;
            cp.mass = obj.mass;
            cp.viscosity = obj.viscosity;
            cp.dispComposition = obj.dispComposition;
            cp.affinity = obj.affinity;
            cp.X = obj.X;
            cp.activity = obj.activity;
            cp.mu0 = obj.mu0;
            cp.mu = obj.mu;

        end

    end

    methods (Access = protected)

        function cp = copyElement(obj)
            % Copy the engine - keeps (duplicates) same engine input as previous engine; engine output is reset.
            % Handles point to the same entities, unless reset
            cp = copyElement@matlab.mixin.Copyable(obj);
            % Engine input is duplicated; engine output is reset
            cp.g = MELTSmap();
            cp.h = MELTSmap();
            cp.s = MELTSmap();
            cp.v = MELTSmap();
            cp.cp = MELTSmap();
            cp.dcpdt = MELTSmap();
            cp.dvdt = MELTSmap();
            cp.dvdp = MELTSmap();
            cp.d2vdt2 = MELTSmap();
            cp.d2vdtdp = MELTSmap();
            cp.d2vdp2 = MELTSmap();
            cp.molwt = MELTSmap();
            cp.rho = MELTSmap();
            cp.mass = MELTSmap();
            cp.viscosity = MELTSmap();
            cp.dispComposition = MELTSmap();
            cp.affinity = MELTSmap();
            cp.X = MELTSmap();
            cp.activity = MELTSmap();
            cp.mu0 = MELTSmap();
            cp.mu = MELTSmap();

        end


   end

end
