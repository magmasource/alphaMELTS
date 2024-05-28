% The 'tutorial' example, see:
% https://magmasource.caltech.edu/forum/index.php?board=12.0
% https://magmasource.caltech.edu/forum/index.php/board=31.0

warning('off', 'MATLAB:loadlibrary:cppoutput');

liquidus = MELTSdynamic(1);

display(liquidus.endMemberFormulas("bulk"))

pressure = 500.0;
temperature = 1200.0;

%{
bulk( 1) = 48.68;
bulk( 2) =  1.01;
bulk( 3) = 17.64;
bulk( 4) =  0.89;
%bulk( 5) =  0.0425;
bulk( 5) =  0.03;
bulk( 6) =  7.59;
bulk( 7) =  0.0;
bulk( 8) =  9.10;
bulk( 9) =  0.0;
bulk(10) =  0.0;
bulk(11) = 12.45;
bulk(12) =  2.65;
bulk(13) =  0.03;
bulk(14) =  0.08;
bulk(15) =  0.2;
bulk(16) =  0.0;
bulk(17) =  0.0;
bulk(18) =  0.0;
bulk(19) =  0.0;
%}

% Can be a row or a column but  is a column
bulk = zeros(19, 1);
bulk([1:6 8 11:15]) = [48.68 1.01 17.64 0.89 0.03 7.59 9.10 12.45 2.65 0.03 0.08 0.2];

% There are several ways to do this
liquidus.engine.set('bulkComposition', bulk);
liquidus.engine.pressure = pressure;

% Note that as of Nov 18th, 2020 temperature is in Celcius
liquidus.engine.temperature = temperature;

% Find liquidus
liquidus.engine.findLiquidus;
temperature = liquidus.engine.temperature;
disp(liquidus.engine.status.message);
display(temperature);

% Note that as of July 4th, 2022 .calcSaturationState() outputs wt% oxides
liquidusPhase = liquidus.engine.solidNames;
display(liquidusPhase);
liquidus.engine.calcEndMemberProperties(liquidusPhase, liquidus.engine.dispComposition(liquidusPhase))
Xan = liquidus.engine.getProperty('X', liquidusPhase, "CaAl2Si2O8");
display(Xan);

% Tutorial example has no fO2 buffer (this is just illustration of different ways to set properties)
%ptpath.engine.setSystemProperties(["Log fO2 Path: FMQ";...
%    "Mode: Fractionate Solids"]);
%ptpath.engine.setSystemProperties(["Log fO2 Path" "FMQ";...
%    "Mode" "Fractionate Solids"]);

liquidus.engine.setSystemProperties("Mode", "Fractionate Solids");

% Starts a new list (which reloads the library)
ptpath = liquidus.copyAndKeepOutput;

% Run mode is 1 (isobaric, isothermal) if empty
% Select output mode 0, to equilibrate but not do any fractionations
ptpath.engine.calcEquilibriumState(1, 0);
disp(ptpath.engine.status.message);

while ptpath.engine.temperature >= 1000

    ptpath = ptpath.addNodeAfter;
    ptpath.engine.temperature = ptpath.engine.temperature - 3;

    % Select 1 to get output after equilibration and before fractionation, 2 for output after fractionation
    % (either way, bulk composition will be updated after fractionation)
    ptpath.engine.calcEquilibriumState(1, 1);
    disp(ptpath.engine.status.message);

    % Can populate X, mu0 ,mu as go along (should test that the phase exists)
    if any(ptpath.engine.solidNames == "plagioclase1")
        ptpath.engine.calcEndMemberProperties("plagioclase1", ptpath.engine.dispComposition("plagioclase1"))
    end

end

indices = 1:ptpath.Last.nodeIndex;
temp = ptpath.getListProperty('temperature');
mgo = ptpath.getListProperty('dispComposition', 'liquid1', 'MgO');
mass = ptpath.getListProperty('mass', 'bulk');


% A silly example to illustrate sorting the list (in this case by temperature ascending)
% We could have just flipped indices, but instead use sort to pass a sorted indices
figure(1)
plot(indices, mgo)
hold on
[~, order] = sort(temp); % sort by t ascending
ptpath.sortListBy(order);
mgo = ptpath.getListProperty('dispComposition', 'liquid1', 'MgO');
% Note that nodeIndex is reset i.e. ptpath.First.nodeIndex is 1
figure(1)
plot(indices, mgo)

% flip it back again
ptpath.sortListBy(order);
mgo = ptpath.getListProperty('dispComposition', 'liquid1', 'MgO');

Xan = ptpath.getListProperty('X', "plagioclase1", "CaAl2Si2O8");
figure(2)
plot(temp, Xan);

% If we didn't calculate Xan as we went along and wanted to avoid repeating the equilibration loop.

%{
% For single phase composition a row is returned
plag = zeros(19, 75);
plag( 1, :) = ptpath.getListProperty('dispComposition', 'plagioclase1', 'SiO2');
plag( 3, :) = ptpath.getListProperty('dispComposition', 'plagioclase1', 'Al2O3');
plag(11, :) = ptpath.getListProperty('dispComposition', 'plagioclase1', 'CaO');
plag(12, :) = ptpath.getListProperty('dispComposition', 'plagioclase1', 'Na2O');
plag(13, :) = ptpath.getListProperty('dispComposition', 'plagioclase1', 'K2O');
%}

plag = ptpath.getListProperty('dispComposition', 'plagioclase1');

% Need to give the feldspars from different runs different names that start with feldspar
phaseList = 1:ptpath.Last.nodeIndex;
phaseList = "plagioclase_"+phaseList;

% This works because Xan is independent of pressure and temperature but we need to specify them
melts = MELTSdynamic(1); % start a new list (as ptpath already has a value for Xan)
melts.engine.temperature = ptpath.engine.temperature;
melts.engine.pressure = ptpath.engine.pressure;

% There is no feldspar in the first ptpath phase assemblage (composition will be all NaNs; check SiO2)
indices = ~isnan(plag(1, :));
melts.engine.calcEndMemberProperties(phaseList(indices), plag(:, indices));

% To get properties from a single node, use getNodeProperty with [] for the index.
% Or equivalently use the engine's getProperty method instead.
%Xan = NaN(size(temp));
%Xan(indices) = melts.engine.getProperty('X', phaseList(indices), "CaAl2Si2O8");

Xan = NaN(size(temp));
Xan(indices) = melts.getNodeProperty([], 'X', phaseList(indices), "CaAl2Si2O8");

figure(2)
plot(temp, Xan, 'o');

temp = ptpath.getListProperty('temperature')';
sio2 = ptpath.getListProperty('dispComposition', 'liquid1', 'SiO2')';
feo = ptpath.getListProperty('dispComposition', 'liquid1', 'FeO')';
%mgo = ptpath.getListProperty('dispComposition', 'liquid1', 'MgO')';
mgo = mgo';
cao = ptpath.getListProperty('dispComposition', 'liquid1', 'CaO')';
na2o = ptpath.getListProperty('dispComposition', 'liquid1', 'Na2O')';

delete GMT_Plot.pdf

gmt_plot = false;
try
    gmt('begin', 'GMT_Plot', 'pdf')
    gmt_plot = true;
catch
    warning("Could not find gmt")
end

if gmt_plot
    gmt ('basemap', '-R1000/1220/0/100 -JX-4i/3i -BWeSn+t"Melt composition" -Bx+l"T @.C system" -By+L"wt%"')
    gmt ('plot', '-R1000/1220/0/100 -Wthin,lightblue -Sc0.1', [temp sio2])
    gmt ('plot', '-R1000/1220/0/100 -Wthin,orange -Sc0.1', [temp feo])
    gmt ('plot', '-R1000/1220/0/100 -Wthin,magenta -Sc0.1', [temp mgo])
    gmt ('plot', '-R1000/1220/0/100 -Wthin,gray -Sc0.1', [temp cao])
    gmt ('plot', '-R1000/1220/0/100 -Wthin,pink -Sc0.1', [temp na2o])
    %#S [dx1 symbol size fill pen [ dx2 text ]]
    specs = {'S 0.1i s 0.15i lightblue 0.25p 0.2i SiO@-2@-';...
            'S 0.1i s 0.15i orange 0.25p 0.2i FeO';... 
            'S 0.1i s 0.15i magenta 0.25p 0.2i MgO';...
            'S 0.1i s 0.15i gray 0.25p 0.2i CaO';...
            'S 0.1i s 0.15i pink 0.25p 0.2i Na@-2@-O'};
    gmt ('legend', '-DjTL -F', specs);
    gmt ('end')
    open GMT_PLot.pdf
end
