%% Exercise: isobaric fractional crystallization of a back-arc basalt

% This just suppresses a warning you might get if you do not have the full Xcode installation on Mac
warning('off', 'MATLAB:loadlibrary:cppoutput');

% For this calculation, which involves H2O-CO2 fluid in a mafic system, we will use rhyoliteMELTS 1.2.0
%     1. rhyolite-MELTS 1.0.2
%     2. pMELTS
%     3. rhyolite-MELTS 1.1.0
%     4. rhyolite-MELTS 1.2.0

ptpath = MELTSdynamic(4);

% Setup populates the systemNames (all available phases for rhyoliteMELTS 1.2.0), endMemberFormulas
% and endMemberWeights
display(ptpath.endMemberFormulas("bulk"))

%% The important stuff is all in the MELTSengine

%{
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
%}

bulk = zeros(19, 1);
bulk([1:9 11:16]) = [51.27 0.96 14.80 1.1 0.042 8.7 0.176 8.45 0.011 12.47 2.25 0.047 0.065 0.25 0.005];

% There are several ways to do set composition 
ptpath.engine.set('bulkComposition', bulk);

% Here we passed the whole composition vector but it is possible to set individual oxides separately
%ptpath.engine.set('bulkComposition', 'H2O', 0.5);

pressure = 250.0;
temperature = 1191.0;

% Note that as of Nov 18th, 2020 temperature is in Celcius
ptpath.engine.pressure = pressure;
ptpath.engine.temperature = temperature;

% Either way works for setting system properties.
% A 1-D string array of strings that look just like the .melts file lines:
%ptpath.engine.setSystemProperties(["Log fO2 Path: -1FMQ"; "Mode: Fractionate Solids"]);

% A 2-D string array of property and value pairs (fo2 Offset also supported):
settings = ["Log fO2 Path" "-1FMQ";...
    "Mode" "Fractionate Solids"];
disp(settings)

ptpath.engine.setSystemProperties(settings);

%% Execute the path (equivalent to menu option 4)

% Select 1 to get output after equilibration and before fractionation, 2 for output after fractionation
% (either way, bulk composition will be updated after fractionation)

ptpath.engine.calcEquilibriumState(1, 1);

while ptpath.engine.temperature >= 600
    
    ptpath = ptpath.addNodeAfter;
    % Increment Temperature: -1
    ptpath.engine.temperature = ptpath.engine.temperature - 1;

    ptpath.engine.calcEquilibriumState(1, 1);
    disp(ptpath.engine.status.message);

    % This means that we do not skip any failures; a counter could be set up to allow limited skips
    if ptpath.engine.status.failed; break; end

    % Can populate X, mu0, mu for all phases or select key phases; can also calculate afterwards
    % (see tutorial.m). Passing [] for the nodeIndex gives the current node.
    solidNames = ptpath.engine.solidNames;
    if ~isempty(solidNames)
        solidComps = ptpath.getNodeProperty([], 'dispComposition', solidNames);
        ptpath.engine.calcEndMemberProperties(solidNames, solidComps);
    end
    
end

%% Examine output

% By default alphaMELTS for MATLAB automatically its own _tbl.txt files that can be
% processed with run-alphamelts.command -x to generate Phase_mass_tbl.txt etc.

% Note that ptpath.engine.meltsIndex gives you the line number of the output in the tbl or _tbl.txt 
% files, wherease ptpath.nodeIndex gives you the location within the MELTSdynamic linked list
% (these won't be the same if the library gets reloaded and calculations continue)

%% Liquidus phase(s)

% You can also explore the liquidus phase by looking at the first few iterations / nodes
solidNames1 = ptpath.First.engine.solidNames;
display(solidNames1)

% These two are almost equivalent
solidNames2 = ptpath.First.Next.engine.solidNames;
display(solidNames2)
solidNames2 = ptpath.getNodeProperty(2, 'solidNames');
display(solidNames2)

liquidusPhase = ptpath.getNodeProperty(3, 'solidNames');
display(liquidusPhase)

temp = ptpath.getListProperty('temperature');

Xan = ptpath.getListProperty('X', "plagioclase1", "CaAl2Si2O8");
figure(1)
plot(temp, Xan);

solidNames4 = ptpath.getNodeProperty(4, 'solidNames');
display(solidNames4)

cpx = ptpath.getListProperty('mass', "clinopyroxene1");
figure(2)
plot(temp, cpx)

%% Other phases

indices = 1:ptpath.Last.nodeIndex;
sp = ptpath.getListProperty('mass', "spinel1");
indices = indices(~isnan(sp));

% spinel1 needs to be in "" rather than '' (this is a minor bug)
spFirst = ptpath.getNodeProperty(indices(1), 'X', "spinel1");
spLast = ptpath.getNodeProperty(indices(end), 'X', "spinel1");
display([ptpath.endMemberFormulas("spinel") spFirst spLast])

Xfo = ptpath.getListProperty('X', "olivine1", "Mg2SiO4");
figure(3)
scatter(Xan, Xfo)

%% Liquid

% We need to do Prev because the last point failed
display([ptpath.endMemberFormulas("bulk") ptpath.Prev.engine.dispComposition("liquid1")])

F = ptpath.getListProperty('mass', "liquid1") / ptpath.First.engine.mass("bulk");
figure(4)
plot(temp, F)

feot = ptpath.getListProperty('dispComposition', 'liquid1', 'FeO') +...
    0.9*ptpath.getListProperty('dispComposition', 'liquid1', 'Fe2O3');
mgo = ptpath.getListProperty('dispComposition', 'liquid1', 'MgO');
figure(5)
plot(mgo, feot)