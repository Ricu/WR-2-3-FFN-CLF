clear; clc;
addpath('libs')
export = 1;
plot_grid = 0;   % Auswahl: Plotten der Triangulierung mit Kanal-Koeffizientenfunktion

fprintf("############ Erstelle Trainingsdaten Start ############\n")
fprintf("Startzeit %s\n", datestr(datetime))
if export
    fprintf("Die Daten werden gespeichert\n")
else
    fprintf("Die Daten werden nicht gespeichert\n")
end

%% Funktion rechte Seite
f = @(vert,y) ones(size(vert));   % Rechte Seite der DGL

%% Lade vertTris fuer schnellere Berechnung der coeff-funktion
vertTris = load("./libs/precomputed_vertTris.mat").vertTris;
%% Erstelle das Gitter
n = 40;         % 2*n^2 Elemente pro Teilgebiet
N = 4;          % Partition in NxN quadratische Teilgebiete
H = 1/N;
h = 1/(N*n);
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

%% Vorbereitung benoetigte Kanten & TG
% Pruefe ob eines der beteiligten TG einen Dirichletknoten enthaelt.
% Falls ja ist eins der TG kein floating TG und wird daher nicht
% beruecksichtigt
floatingSD = false(numSD,1);
for sd = 1:numSD
    floatingSD(sd) = nnz(dirichlet(l2g__sd{sd})) == 0;
end

edgesSD = [(1:numSD-N)',(N+1:numSD)';...
           setdiff(1:numSD,N:N:numSD)', setdiff(1:numSD,1:N:numSD)'];
edgesSD = sortrows(edgesSD);
floatingEdges = all(floatingSD(edgesSD),2);
validEdges = find(floatingEdges);
fprintf("Fuer das gegebene Gitter werden %i (%4.1f%%) Kanten uebersprungen\n",sum(~floatingEdges),sum(~floatingEdges)/length(floatingEdges)*100)
% Schmeiße alle anderen Kanten raus

%% Koeffizientenfunktion vorbereiten
TOL = 100;  % Toleranz zur Auswahl der Eigenwerte
rng(42);

% Anzahl an Trainingsamples pro Koeffizentenfunktionen 
nSamplesConstant    = 3;
nSamplesStrips      = 50;
nSamplesBlocks      = 1;
nSamplesRandBlocks  = 1;
nSamplesRand        = 1;
nSamples = nSamplesConstant+nSamplesStrips+nSamplesBlocks+nSamplesRandBlocks+nSamplesRand;
coeffFun_cell = cell(nSamples*4,1); % Faktor 10 kommt daher, dass ... ?
coeffFun_counter = 1;
parameter_cell = cell(nSamples*4,4); % Faktor 4 kommt daher, dass ... ?

rhoBound = 10.^[0,6]; % enthaelt den minimalen und maximalen Koeffizienten
indexShiftBound = 0:2:20; % Array moeglicher Verschiebungen der Elemente mit hoeherem Koeffizienten in x- oder y-Richtung

%% Konstante Koeffizientenfunktion
% Test verschiedene Parameter fuer die Kanalfunktion

% Erstelle die Parameterstruktur
param_names = ["affectedSubdomains","rhoMin","rhoMax","indexShiftx","indexShifty"];
fprintf("%s: Insgesamt %i Parameter zur Auswahl.\n","Constant",length(param_names))
affectedSubdomains = [6,10]; % TG-Wahl fuer hoeherem Koeffizienten decken alle moeglichen Faelle ab
% Samples hier haendisch erstellen
parameter_const = {affectedSubdomains; affectedSubdomains};
parameter_const = [parameter_const, num2cell([rhoBound;circshift(rhoBound,1,2);0,0;0,0])']';

sample_parameters = cell2struct(parameter_const,param_names,1);

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
% Test verschiedene Parameter fuer die Kanalfunktion
heightBound = -2:4;  % Breite der Kanaele, 0 ist dabei eine initiale Breite abhaengig von der Anzahl an Kanaelen je TG
nStripsBound = 1:5; % Gibt die Anazhl Kanaele je TG an

% Erstelle die Parameterstruktur
param_names = ["rhoMin","rhoMax","height","nStrips","indexShifty"];
fprintf("%s: Insgesamt %i Parameter zur Auswahl.\n","Strip",length(param_names))
sample_parameters = generateSampleParameters(nSamplesStrips,param_names,rhoBound,heightBound,nStripsBound,indexShiftBound);
% sample_parameters(6) = cell2struct(num2cell([0; 1; 2; 3; 4]),param_names,1);

for sampleID = 1:length(sample_parameters)
    height      = sample_parameters(sampleID).height;
    nStrips     = sample_parameters(sampleID).nStrips;
    rhoMin      = sample_parameters(sampleID).rhoMin;
    rhoMax      = sample_parameters(sampleID).rhoMax;
    indexShifty = sample_parameters(sampleID).indexShifty;

    coeffFun_cell{coeffFun_counter} = @(vertices) coeffFun_canal(vertices(:,2),N,n,height,nStrips,indexShifty);
    parameter_cell{coeffFun_counter,1} = "Strip";
    parameter_cell{coeffFun_counter,2} = param_names;
    parameter_cell{coeffFun_counter,3} = [rhoMin,rhoMax,height,nStrips,indexShifty];
    parameter_cell{coeffFun_counter,4}  = sample_parameters(sampleID);
    coeffFun_counter = coeffFun_counter + 1;
end

%% Bloecke Koeffizientenfunktion
% Teste verschiedene Parameter fuer die Kanalfunktion
difBound = -20:2:20; % Gibt an, wie weit die Bloecke in jedem zweiten TG (spaltenweise ab 2.Spalte) 
                     % voneinander versetzt sind. 0 entspricht keiner Versetzung
prop1Bound = 0:0.2:1; % Gibt den Anteil an Block in jedem zweiten TG (spaltenweise ab 1.Spalte) an
prop2Bound = 0:0.2:1; % Gibt den Anteil an Block in jedem zweiten TG (spaltenweise ab 2.Spalte) an
heightBound = 2:2:38; % Hoehe der Bloecke

% Erstelle die Parameterstruktur
param_names = ["rhoMin","rhoMax","height","dif","prop1","prop2","indexShiftx","indexShifty"];
fprintf("%s: Insgesamt %i Parameter zur Auswahl.\n","Blocks",length(param_names))
sample_parameters = generateSampleParameters(nSamplesBlocks,param_names,rhoBound,heightBound,difBound,prop1Bound,prop2Bound,indexShiftBound,indexShiftBound);

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
% Test verschiedene parameter fuer die Kanalfunktion
widthBound      =  2: 1:10; % Breite der Bloecke mit Faktor der Schrittweite
heightBound     =  2: 1:10; % Hoehe der Bloecke mit Faktor der Schrittweite
varianceBound   =  0: 2:16;  % positive Varianz in Breite und Hoehe
nBlocksBound    = 10:10:90; % Anzahl an random erstellten Bloecken

% Erstelle die Parameterstruktur
param_names = ["rhoMin","rhoMax","nBlocks","height","heightVariance","width","widthVariance","indexShiftx","indexShifty"];
fprintf("%s: Insgesamt %i Parameter zur Auswahl.\n","Random Blocks",length(param_names))
sample_parameters = generateSampleParameters(nSamplesRandBlocks,param_names,rhoBound,nBlocksBound,heightBound,varianceBound,heightBound,varianceBound,indexShiftBound,indexShiftBound);

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

%% Zufalls Koeffizientenfunktion
% Teste verschiedene Parameter fuer die Kanalfunktion
randomPercentageBound = 0.2:0.05:0.7;
randomStateBound = 1:10;

% Erstelle die Parameterstruktur
param_names = ["rhoMin","rhoMax","randomPercentage","randomState","indexShiftx","indexShifty"];
fprintf("%s: Insgesamt %i Parameter zur Auswahl.\n","Completely Random",length(param_names))
sample_parameters = generateSampleParameters(nSamplesRand,param_names,rhoBound,randomPercentageBound,randomStateBound,indexShiftBound,indexShiftBound);

for sampleID = 1:length(sample_parameters)
    randomPercentage    = sample_parameters(sampleID).randomPercentage;
    randomState         = sample_parameters(sampleID).randomState;
    rhoMin              = sample_parameters(sampleID).rhoMin;
    rhoMax              = sample_parameters(sampleID).rhoMax;
    indexShiftx         = sample_parameters(sampleID).indexShiftx;
    indexShifty         = sample_parameters(sampleID).indexShifty;

    coeffFun_cell{coeffFun_counter} = @(vertices) coeffFun_random(vertices(:,1),vertices(:,2),randomPercentage,randomState,indexShiftx,indexShifty);
    parameter_cell{coeffFun_counter,1}  = "Completely Random";
    parameter_cell{coeffFun_counter,2}  = param_names;
    parameter_cell{coeffFun_counter,3}  = [rhoMin,rhoMax,randomPercentage,randomState,indexShiftx,indexShifty];
    parameter_cell{coeffFun_counter,4}  = sample_parameters(sampleID);
    coeffFun_counter = coeffFun_counter + 1;
end

%% Faelle 
empty_cells_ind = cellfun('isempty',coeffFun_cell);
coeffFun_cell = coeffFun_cell(~empty_cells_ind);
parameter_cell = parameter_cell(~empty_cells_ind,:);
n_cases = length(coeffFun_cell);
output_cell = cell(n_cases,1);

t_casesStart = tic;
for case_id = 1:n_cases
    fprintf("#### Starte Fall %5i/%5i: Koeffizientenfunktion %s ####\n",case_id,n_cases,parameter_cell{case_id})
    t_coeffFun = tic;
    % Definiere Koeffizient auf den Elementen (und teilgebietsweise);
    % maximalen Koeffizienten pro Knoten (und teilgebietsweise)
    rhoMin = parameter_cell{case_id,4}.rhoMin;
    rhoMax = parameter_cell{case_id,4}.rhoMax;
    markerType = 'verts';
    [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = getCoefficientMatrices(coeffFun_cell{case_id},markerType,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid,vertTris);
    fprintf("Benoetigte Zeit: Aufstellen der Koeffizientenmatrizen: %5fs",toc(t_coeffFun))
    rho_struct = struct('rhoTriSD',{rhoTriSD},'maxRhoVert',{maxRhoVert},'maxRhoVertSD',{maxRhoVertSD});

    t_matrices = tic;
    [edgesPrimalGlobal,cGamma,edgesSD,cLocalPrimal,cB,cBskal,cInner,cK,cDirichlet] = setup_matrices(rho_struct,grid_struct,f);
    fprintf(", Sprungoperators/Steifigkeitsmatrix: %5fs\n", toc(t_matrices))

    nValidEdges = length(validEdges);
    input = cell(nValidEdges,1);
    label = zeros(nValidEdges,1);

    fprintf("Kanten:")
    for i = 1:length(validEdges)
        edgeID = validEdges(i);
        input{i} = generate_input(edgeID,edgesSD,rhoTriSD,vert__sd,tri__sd);
        label(i) = generate_label(edgeID,edgesPrimalGlobal,cGamma,edgesSD,cLocalPrimal,cB,cBskal,cInner,cK,TOL);
%         fprintf("Kante %2i bzgl. der TG (%2i,%2i) erhaelt das Label %i\n",edgeID,edgesSD(edgeID,1),edgesSD(edgeID,2),label(i))
        fprintf("   %2i", edgeID)
    end
    fprintf("\nLabels:")
    fprintf("   %2i",label)
    fprintf("\n")
    % Fuege neue Daten an den Trainingsdatensatz an
    output_cell{case_id} = [cell2mat(input),label];
    fprintf("Durchschnittliche Zeit pro Fall bisher %fs. Verbleibende Zeit ca: %.2fm\n", toc(t_casesStart)/case_id, (n_cases-case_id)*toc(t_casesStart)/case_id/60)
end

%% Daten exportieren
if export
    % Input-Label Kombinationen der einzelnen Faelle abspeichern
    file_name = sprintf("./resources/train_data/%s-train_data_dump.csv",datestr(datetime,'yyyy-mm-dd-HH-MM-SS'));
    fprintf("Speichere Traininsdaten als %s...",file_name)
    output_mat = cell2mat(output_cell);
    writematrix(output_mat,file_name);
    fprintf("Fertig!\n")
    % Parameter der einzelnen Faelle abspeichern
    file_name2 = sprintf("./resources/train_data/%s-parameter_dump.csv",datestr(datetime,'yyyy-mm-dd-HH-MM-SS'));
    fprintf("Speichere Parameterdaten als %s...",file_name2)
    writecell(parameter_cell(:,1:3),file_name2);
    fprintf("Fertig!\n")
end


function sample_parameters = generateSampleParameters(nRandSamples,param_names,rhoBound,varargin)
rhoMin = min(rhoBound); rhoMax = max(rhoBound);

parameter_cell = cell(numel(varargin),1);
[parameter_cell{:}] = ndgrid(varargin{:});
parameter_cell = cellfun(@(X) reshape(X,1,[]),parameter_cell,'UniformOutput',false);
parameter_cell(:,2) = parameter_cell(:,1);
nVari = length(parameter_cell{1});
rhoTemp = {repmat(rhoMin,1,nVari), repmat(rhoMax,1,nVari); repmat(rhoMax,1,nVari), repmat(rhoMin,1,nVari)};
parameter_cell = [rhoTemp; parameter_cell];
fprintf("Wahle zufaellig %i aus %8i moeglichen Parameterkombinationen aus.\n",nRandSamples,2*length(parameter_cell{1}));

random_permutation = randperm(2*size(parameter_cell{1},2));
parameter_mat = cell2mat(parameter_cell);
parameter_mat = parameter_mat(:,random_permutation(1:nRandSamples));
% parameter_cell = num2cell(parameter_mat);

sample_parameters = cell2struct(num2cell(parameter_mat),param_names,1);
end