clear; clc;
addpath('libs')
export = 1;      % Auswahl: Exportieren und Speichern der generierten Daten
plot_grid = 0;   % Auswahl: Plotten der Triangulierung mit Koeffizientenfunktion

fprintf("############ Erstelle Trainingsdaten Start ############\n")
fprintf("Startzeit %s\n", datestr(datetime))
if export
    fprintf("Die Daten werden gespeichert\n")
else
    fprintf("Die Daten werden nicht gespeichert\n")
end

%% Funktion rechte Seite
f = @(vert,y) ones(size(vert));   % Rechte Seite der DGL

%% Toleranz zur Auswahl der Eigenwerte
TOL = 100;  

%% Erstelle das Gitter
n = 40;         % 2*n^2 Elemente pro Teilgebiet
N = 4;          % Partition in NxN quadratische Teilgebiete
H = 1/N;        % Schrittweite: Teilgebiete
h = 1/(N*n);    % Schrittweite: Elemente
fprintf("Das Verhaeltnis H/h betraegt %f\n",H/h);
numSD = N^2;    % Anzahl Teilgebiete
xyLim = [0,1];  % Gebiet: Einheitsquadrat

[vert,tri] = genMeshSquare(N,n);            % Erstelle Knoten- und Elementliste
numVert=size(vert,1);   numTri=size(tri,1); % Anzahl Knoten und Dreiecke

% Erstelle Knoten- und Elementlisten pro Teilgebiet und logischen Vektor,
% welche Dreiecke in welchem TG enthalten sind
[vert__sd,tri__sd,l2g__sd,logicalTri__sd] = meshPartSquare(N,vert,tri);

% Markiere Dirichletknoten in logischem Vektor
dirichlet = or(ismember(vert(:,1),xyLim), ismember(vert(:,2),xyLim));

% Structure fuer grid-Variablen
grid_struct = struct('vert__sd',{vert__sd},'tri__sd',{tri__sd},'l2g__sd',{l2g__sd},'dirichlet',{dirichlet});

% Lade vertTris fuer schnellere Berechnung der Koeffizientenfunktion
% Enthaelt fuer jeden Knoten die Nummern der anliegenden Elemente
vertTris = load("./libs/precomputed_vertTris.mat").vertTris;

%% Vorbereitung benoetigte Kanten & TG
% Pruefe fuer jede Teilgebietskante, ob eines der anliegenden TG einen 
% Dirichletknoten enthaelt und damit kein floating TG ist.
% Bei der Erstellung der Trainingsdaten werden nur floating Kanten beruecksichtigt

floatingSD = false(numSD,1); % logischer Vektor
for sd = 1:numSD
    floatingSD(sd) = nnz(dirichlet(l2g__sd{sd})) == 0;
end

% Erstelle Kantenliste mit anliegenden TG-Nummern
edgesSD = [(1:numSD-N)',(N+1:numSD)';...
           setdiff(1:numSD,N:N:numSD)', setdiff(1:numSD,1:N:numSD)'];
edgesSD = sortrows(edgesSD);

% Extrahiere floating Teilgebietskanten
floatingEdges = all(floatingSD(edgesSD),2); % logischer Vektor
validEdges = find(floatingEdges);   % enthaelt Kanteninidzes
nValidEdges = length(validEdges);   % Anzahl floating Kanten
fprintf("In dem gegebenen Gitter sind %i (%4.1f%%) Kanten keine floating Kanten und werden uebersprungen.\n", ...
                                            sum(~floatingEdges),sum(~floatingEdges)/length(floatingEdges)*100)

%% Erstellung der Koeffizientenfunktionverteilung vorbereiten
rng(24); % Setze seed fuer random number generator

% Anzahl an Trainingssamples pro Koeffizentenfunktionen 
nSamplesConstant        = 0;    % Konstante Koeffizienten auf den Teilgebieten (maximal 2)
nSamplesStrips          = 5;    % Streifen
nSamplesBlocks          = 5;    % Bloecke
nSamplesRandBlocks      = 5;    % Random Bloecke
nSamplesStripRandBlocks = 5;    % Strip - Random Bloecke
nSamplesRand            = 5;    % Random Elemente

% Anzahl Trainingssamples gesamt
nCases = nSamplesConstant+nSamplesStrips+nSamplesBlocks+nSamplesRandBlocks+nSamplesStripRandBlocks+nSamplesRand;

% Cell-Arrays, in denen Informationen zwischengespeichert werden
coeffFun_cell = cell(nCases,1); 
parameter_cell = cell(nCases,4); 

% Setze counter
coeffFun_counter = 1;

% Globale Parametergrenzen
rhoBound = 10.^[0,6]; % enthaelt den minimalen und maximalen Koeffizienten
indexShiftBound = 0:2:20; % Array moeglicher Verschiebungen der Elemente 
                          % mit hoeherem Koeffizienten in x- oder y-Richtung

%% Erstellen der verschiedenen Koeffizientenverteilungen
%% Konstante Koeffizientenfunktion auf den Teilgebieten
% Erstelle die Parameterstruktur
param_names = ["affectedSubdomains","rhoMin","rhoMax","indexShiftx","indexShifty"];
fprintf("%s: Insgesamt %i Parameter zur Auswahl.\n","Constant",length(param_names))
affectedSubdomains = [6,10]; % TG-Wahl fuer hoeheren Koeffizienten; decken alle moeglichen Faelle ab
% Samples hier haendisch erstellen
parameter_const = {affectedSubdomains; affectedSubdomains};
parameter_const = [parameter_const, num2cell([rhoBound;circshift(rhoBound,1,2);0,0;0,0])']';
sample_parameters = cell2struct(parameter_const,param_names,1);

% Erstelle in einer Schleife die Koeffizientenfunktionen anhand der
% gewaehlten Parameterkombinationen und speichere diese zwischen.
for sampleID = 1:length(sample_parameters)
    affectedSubdomains  = sample_parameters(sampleID).affectedSubdomains;
    rhoMin              = sample_parameters(sampleID).rhoMin;
    rhoMax              = sample_parameters(sampleID).rhoMax;
    indexShiftx         = sample_parameters(sampleID).indexShiftx;
    indexShifty         = sample_parameters(sampleID).indexShifty;

    coeffFun_cell{coeffFun_counter}    = @(vertices) coeffFun_subdomains(vertices(:,1),vertices(:,2),affectedSubdomains,vert__sd,indexShiftx,indexShifty);
    parameter_cell{coeffFun_counter,1} = "Constant";
    parameter_cell{coeffFun_counter,2} = param_names;
    parameter_cell{coeffFun_counter,3} = [affectedSubdomains,rhoMin,rhoMax,indexShiftx,indexShifty];
    parameter_cell{coeffFun_counter,4} = sample_parameters(sampleID);
    coeffFun_counter = coeffFun_counter + 1;
end

%% Streifen Koeffizientenfunktion
% Lege die Parametergrenzen fest
heightBound = -2:4;  % Breite der Streifen, 0 ist initiale Breite abhaengig von der Anzahl Streifen je TG
nStripsBound = 1:5; % Anzahl Streifen je TG

% Erstelle die Parameterstruktur
param_names = ["rhoMin","rhoMax","height","nStrips","indexShifty"];
fprintf("%s: Insgesamt %i Parameter zur Auswahl.\n","Strip",length(param_names))
% Samples hier mit Hilfsfunktion erstellen
sample_parameters = generateSampleParameters(nSamplesStrips,param_names,rhoBound,heightBound,nStripsBound,indexShiftBound);

% Erstelle in einer Schleife die Koeffizientenfunktionen anhand der
% gewaehlten Parameterkombinationen und speichere diese zwischen.
for sampleID = 1:length(sample_parameters)
    height      = sample_parameters(sampleID).height;
    nStrips     = sample_parameters(sampleID).nStrips;
    rhoMin      = sample_parameters(sampleID).rhoMin;
    rhoMax      = sample_parameters(sampleID).rhoMax;
    indexShifty = sample_parameters(sampleID).indexShifty;

    coeffFun_cell{coeffFun_counter} = @(vertices) coeffFun_strip(vertices(:,2),N,n,height,nStrips,indexShifty);
    parameter_cell{coeffFun_counter,1} = "Strip";
    parameter_cell{coeffFun_counter,2} = param_names;
    parameter_cell{coeffFun_counter,3} = [rhoMin,rhoMax,height,nStrips,indexShifty];
    parameter_cell{coeffFun_counter,4}  = sample_parameters(sampleID);
    coeffFun_counter = coeffFun_counter + 1;
end

%% Bloecke Koeffizientenfunktion
% Lege die Parametergrenzen fest
difBound = -20:2:20;  % Gibt an, wie weit die Bloecke in jedem zweiten TG (spaltenweise ab 2.Spalte) 
                      % voneinander versetzt sind. 0 entspricht keiner Versetzung
prop1Bound = 0:0.2:1; % Gibt den Anteil an Block in jedem zweiten TG (spaltenweise ab 1.Spalte) an
prop2Bound = 0:0.2:1; % Gibt den Anteil an Block in jedem zweiten TG (spaltenweise ab 2.Spalte) an
heightBound = 2:2:38; % Hoehe der Bloecke

% Erstelle die Parameterstruktur
param_names = ["rhoMin","rhoMax","height","dif","prop1","prop2","indexShiftx","indexShifty"];
fprintf("%s: Insgesamt %i Parameter zur Auswahl.\n","Blocks",length(param_names))
% Samples hier mit Hilfsfunktion erstellen
sample_parameters = generateSampleParameters(nSamplesBlocks,param_names,rhoBound,heightBound,difBound,prop1Bound,prop2Bound,indexShiftBound,indexShiftBound);

% Erstelle in einer Schleife die Koeffizientenfunktionen anhand der
% gewaehlten Parameterkombinationen und speichere diese zwischen.
for sampleID = 1:length(sample_parameters)
    height      = sample_parameters(sampleID).height;
    dif         = sample_parameters(sampleID).dif;
    prop1       = sample_parameters(sampleID).prop1;
    prop2       = sample_parameters(sampleID).prop2;
    rhoMin      = sample_parameters(sampleID).rhoMin;
    rhoMax      = sample_parameters(sampleID).rhoMax;
    indexShiftx = sample_parameters(sampleID).indexShiftx;
    indexShifty = sample_parameters(sampleID).indexShifty;

    coeffFun_cell{coeffFun_counter} = @(vertices) coeffFun_block(vertices(:,1), vertices(:,2), N, n, prop1,prop2,dif,height,indexShiftx,indexShifty);
    parameter_cell{coeffFun_counter,1}  = "Blocks";
    parameter_cell{coeffFun_counter,2}  = param_names;
    parameter_cell{coeffFun_counter,3}  = [rhoMin,rhoMax,height,dif,prop1,prop2,indexShiftx,indexShifty];
    parameter_cell{coeffFun_counter,4}  = sample_parameters(sampleID);
    coeffFun_counter = coeffFun_counter + 1;
end

%% Zufalls - Bloecke Koeffizientenfunktion
% Lege die Parametergrenzen fest
widthBound      =  4:4:16;  % Breite der Bloecke mit Faktor der Schrittweite
heightBound     =  4:4:16;  % Hoehe der Bloecke mit Faktor der Schrittweite
varianceBound   =  0:4:16;  % positive Varianz in Breite und Hoehe
nBlocksBound    =  10:5:20; % Anzahl an random erstellten Bloecken

% Erstelle die Parameterstruktur
param_names = ["rhoMin","rhoMax","nBlocks","height","heightVariance","width","widthVariance","indexShiftx","indexShifty"];
fprintf("%s: Insgesamt %i Parameter zur Auswahl.\n","Random Blocks",length(param_names))
% Samples hier mit Hilfsfunktion erstellen
sample_parameters = generateSampleParameters(nSamplesRandBlocks,param_names,rhoBound,nBlocksBound,heightBound,varianceBound,widthBound,varianceBound,indexShiftBound,indexShiftBound);

% Erstelle in einer Schleife die Koeffizientenfunktionen anhand der
% gewaehlten Parameterkombinationen und speichere diese zwischen.
for sampleID = 1:length(sample_parameters)
    nBlocks         = sample_parameters(sampleID).nBlocks;
    height          = sample_parameters(sampleID).height;
    heightVariance  = sample_parameters(sampleID).heightVariance;
    width           = sample_parameters(sampleID).width;
    widthVariance   = sample_parameters(sampleID).widthVariance;
    rhoMin          = sample_parameters(sampleID).rhoMin;
    rhoMax          = sample_parameters(sampleID).rhoMax;
    indexShiftx     = sample_parameters(sampleID).indexShiftx;
    indexShifty     = sample_parameters(sampleID).indexShifty;

    coeffFun_cell{coeffFun_counter} = @(vertices) coeffFun_randomBlocks(vertices(:,1),vertices(:,2),N,n,nBlocks,width:width+widthVariance,height:height+heightVariance,indexShiftx,indexShifty);
    parameter_cell{coeffFun_counter,1}  = "Random Blocks";
    parameter_cell{coeffFun_counter,2}  = param_names;
    parameter_cell{coeffFun_counter,3}  = [rhoMin,rhoMax,nBlocks,height,heightVariance,width,widthVariance,indexShiftx,indexShifty];
    parameter_cell{coeffFun_counter,4}  = sample_parameters(sampleID);
    coeffFun_counter = coeffFun_counter + 1;
end

%% Streifen + Zufalls - Bloecke Koeffizientenfunktion
% Lege die Parametergrenzen fest fuer die Random-Bloecke
widthBound      =  1:1:3; % Breite der Bloecke mit Faktor der Schrittweite
heightBound     =  2:1:4; % Hoehe der Bloecke mit Faktor der Schrittweite
varianceBound   =  0:1:2;  % positive Varianz in Breite und Hoehe
nBlocksBound    = 10:10:50; % Anzahl an random erstellten Bloecken

% Lege die Parametergrenzen fest fuer die Streifen
heightBoundS = -2:2;  % Breite der Streifen, 0 ist dabei eine initiale Breite abhaengig von der Anzahl an Streifen je TG
nStripBound = 1:4; % Gibt die Anzahl Streifen je TG an

% Erstelle die Parameterstruktur
param_names = ["rhoMin","rhoMax","nBlocks","height","heightVariance","width","widthVariance","heightBoundS","nStripBound","indexShifty"];
fprintf("%s: Insgesamt %i Parameter zur Auswahl.\n","StripesNBlocks",length(param_names))
% Samples hier mit Hilfsfunktion erstellen
sample_parameters = generateSampleParameters(nSamplesStripRandBlocks,param_names,rhoBound,nBlocksBound,heightBound,varianceBound,widthBound,varianceBound,heightBoundS,nStripBound,indexShiftBound);

% Erstelle in einer Schleife die Koeffizientenfunktionen anhand der
% gewaehlten Parameterkombinationen und speichere diese zwischen.
for sampleID = 1:length(sample_parameters)
    nBlocks         = sample_parameters(sampleID).nBlocks;
    height          = sample_parameters(sampleID).height;
    heightVariance  = sample_parameters(sampleID).heightVariance;
    width           = sample_parameters(sampleID).width;
    widthVariance   = sample_parameters(sampleID).widthVariance;
    heightBoundS    = sample_parameters(sampleID).heightBoundS;
    nStripBound     = sample_parameters(sampleID).nStripBound;
    rhoMin          = sample_parameters(sampleID).rhoMin;
    rhoMax          = sample_parameters(sampleID).rhoMax;
    indexShifty     = sample_parameters(sampleID).indexShifty;

    coeffFun_cell{coeffFun_counter} = @(vertices) coeffFun_stripRandomBlocks(vertices(:,1),vertices(:,2),N,n,heightBoundS,nStripBound,nBlocks,width:width+widthVariance,height:height+heightVariance,indexShifty);
    parameter_cell{coeffFun_counter,1}  = "StripesNBlocks";
    parameter_cell{coeffFun_counter,2}  = param_names;
    parameter_cell{coeffFun_counter,3}  = [rhoMin,rhoMax,heightBoundS,nStripBound,nBlocks,height,heightVariance,width,widthVariance,indexShifty];
    parameter_cell{coeffFun_counter,4}  = sample_parameters(sampleID);
    coeffFun_counter = coeffFun_counter + 1;
end

%% Zufalls - Elemente Koeffizientenfunktion
% Lege die Parametergrenzen fest
randomPercentageBound   = 0.2:0.05:0.7; % Anteil Elemente mit max Koeffizienten
randomStateBound        = 1:1:10;   % verschiedene random states

% Erstelle die Parameterstruktur
param_names = ["rhoMin","rhoMax","randomPercentage","randomState","indexShiftx","indexShifty"];
fprintf("%s: Insgesamt %i Parameter zur Auswahl.\n","Completely Random",length(param_names))
% Samples hier mit Hilfsfunktion erstellen
sample_parameters = generateSampleParameters(nSamplesRand,param_names,rhoBound,randomPercentageBound,randomStateBound,indexShiftBound,indexShiftBound);

% Erstelle in einer Schleife die Koeffizientenfunktionen anhand der
% gewaehlten Parameterkombinationen und speichere diese zwischen.
for sampleID = 1:length(sample_parameters)
    randomPercentage    = sample_parameters(sampleID).randomPercentage;
    randomState         = sample_parameters(sampleID).randomState;
    rhoMin              = sample_parameters(sampleID).rhoMin;
    rhoMax              = sample_parameters(sampleID).rhoMax;
    indexShifty         = sample_parameters(sampleID).indexShifty;

    coeffFun_cell{coeffFun_counter} = @(vertices) coeffFun_random(vertices(:,1),vertices(:,2),randomPercentage,randomState,indexShiftx,indexShifty);
    parameter_cell{coeffFun_counter,1}  = "Completely Random";
    parameter_cell{coeffFun_counter,2}  = param_names;
    parameter_cell{coeffFun_counter,3}  = [rhoMin,rhoMax,randomPercentage,randomState,indexShiftx,indexShifty];
    parameter_cell{coeffFun_counter,4}  = sample_parameters(sampleID);
    coeffFun_counter = coeffFun_counter + 1;
end

%% Erstelle Trainingsdaten
% Ueberpruefe, dass keine leeren Eintraege in dem Koeffizientenarray
% enthalten sind, welche spaeter zu Fehlern führen koennten.
empty_cells_ind = cellfun('isempty',coeffFun_cell);
coeffFun_cell = coeffFun_cell(~empty_cells_ind);
parameter_cell = parameter_cell(~empty_cells_ind,:);

% Erstelle die Arrays, welche spaeter exportiert und gespeichert werden.
output_cell = cell(nCases * nValidEdges,1);
parameter_output_cell = cell(nCases * nValidEdges,3);

t_casesStart = tic;
for case_id = 1:nCases
    fprintf("#### Starte Fall %5i/%5i: Koeffizientenfunktion %s ####\n",case_id,nCases,parameter_cell{case_id})
    t_coeffFun = tic;
    
    rhoMin = parameter_cell{case_id,4}.rhoMin;
    rhoMax = parameter_cell{case_id,4}.rhoMax;
    markerType = 'verts'; % Koeffizientenfunktion ist knotenweise definiert
    
    % Erstelle die Koeffizientenmatrizen, welche in der Berechnung der
    % FETI-DP Matrizen benoetigt werden
    [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = getCoefficientMatrices(coeffFun_cell{case_id},markerType,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid,vertTris);
    fprintf("Benoetigte Zeit: Aufstellen der Koeffizientenmatrizen: %5fs",toc(t_coeffFun))
    rho_struct = struct('rhoTriSD',{rhoTriSD},'maxRhoVert',{maxRhoVert},'maxRhoVertSD',{maxRhoVertSD});

    t_matrices = tic;
    % Erstelle die FETI-DP Matrizen, welche zur Berechnung der Label benoetigt werden
    [edgesPrimalGlobal,cGamma,~,cLocalPrimal,cB,cBskal,cInner,cK,cDirichlet] = setup_matrices(rho_struct,grid_struct,f);
    fprintf(", Sprungoperators/Steifigkeitsmatrix: %5fs\n", toc(t_matrices))

    fprintf("Kanten:")
    for i = 1:length(validEdges)
        edgeID = validEdges(i);
        % Input fuer neuronales Netz: Koeffizientenverteilung elementweise
        input = generate_input(edgeID,edgesSD,rhoTriSD,vert__sd,tri__sd);
        % Zugehoeriger label: 1 kritische Kante, 0 unktritische Kante
        label = generate_label(edgeID,edgesPrimalGlobal,cGamma,edgesSD,cLocalPrimal,cB,cBskal,cInner,cK,TOL);
        % Koeffizientenverteilung mit zugehoerigem label als Trainingsdaten
        output_cell{(case_id-1)*nValidEdges + i,:} =  [input,label];
        % Zugehoerige Parameter zur Nachvollziehbarkeit
        parameter_output_cell((case_id-1)*nValidEdges + i,:) = parameter_cell(case_id,1:3);
        fprintf("   %2i", edgeID)
    end
    fprintf("\nLabels:")
    fprintf("   %2i",label)
    fprintf("\n")
    elapsedTime = toc(t_casesStart);
    remainingTime = (nCases-case_id)*elapsedTime/case_id;
    
    fprintf("Durchschnittliche Zeit pro Fall bisher %fs. Verbleibende Zeit ca: %.2fm (%s)\n", elapsedTime/case_id, remainingTime/60, datestr(datetime + seconds(remainingTime),'dd.mm.yyyy HH:MM:SS'))
end

%% Daten exportieren
if export
    % Input-Label Kombinationen der einzelnen Faelle abspeichern
    export_time = datestr(datetime,'yyyy-mm-dd-HH-MM-SS');
    file_name = sprintf("./resources/train_data/%s-train_data_dump.csv",export_time);
    fprintf("Speichere Traininsdaten als %s...",file_name)
    output_mat = cell2mat(output_cell);
    writematrix(output_mat,file_name);
    fprintf("Fertig!\n")
    % Parameter der einzelnen Faelle abspeichern
    file_name2 = sprintf("./resources/train_data/%s-parameter_dump.csv",export_time);
    fprintf("Speichere Parameterdaten als %s...",file_name2)
    writecell(parameter_output_cell,file_name2);
    fprintf("Fertig!\n")
end

%% Hilfsfunktion zur Generierung der sample Parameter
function sample_parameters = generateSampleParameters(nRandSamples,param_names,rhoBound,varargin)
% Erstelle eine Auswahl aller moeglichen Kombinationen innerhalb der in
% varargin uebergebenen Parametergrenzen + die Kombinationen aus rhoMax und
% rhoMin.

rhoMin = min(rhoBound); rhoMax = max(rhoBound);

% Erstelle alle moeglichen Kombinationen von Parametern in varargin
parameter_cell = cell(numel(varargin),1);
[parameter_cell{:}] = ndgrid(varargin{:});
parameter_cell = cellfun(@(X) reshape(X,1,[]),parameter_cell,'UniformOutput',false);
parameter_cell(:,2) = parameter_cell(:,1);

% Fuege die rhoMin und rhoMax Parameter hinzu: einmal normal und einmal
% umgedreht
nVari = length(parameter_cell{1});
rhoTemp = {repmat(rhoMin,1,nVari), repmat(rhoMax,1,nVari); repmat(rhoMax,1,nVari), repmat(rhoMin,1,nVari)};
parameter_cell = [rhoTemp; parameter_cell];
fprintf("Wahle zufaellig %i aus %8i moeglichen Parameterkombinationen aus.\n",nRandSamples,2*length(parameter_cell{1}));

% Permutiere die Parameterkombinationen und waehle anschliessend
% nRandSamples aus
random_permutation = randperm(2*size(parameter_cell{1},2));
parameter_mat = cell2mat(parameter_cell);
parameter_mat = parameter_mat(:,random_permutation(1:nRandSamples));

% Wandle die erstellten samples in eine structure um
sample_parameters = cell2struct(num2cell(parameter_mat),param_names,1);
end