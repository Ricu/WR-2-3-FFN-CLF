clear; clc;
addpath('libs')
export = 1;
plot_grid = 1;   % Auswahl: Plotten der Triangulierung mit Kanal-Koeffizientenfunktion

fprintf("############ Erstelle Trainingsdaten Start ############\n")
fprintf("Startzeit %s\n", datestr(datetime))
if export
    fprintf("Die Daten werden gespeichert\n")
else
    fprintf("Die Daten werden nicht gespeichert\n")
end

%% Funktion rechte Seite
f = @(vert,y) ones(size(vert));   % Rechte Seite der DGL

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

%% Koeffizientenfunktion vorbereiten
TOL = 100;  % Toleranz zur Auswahl der Eigenwerte
rng(0);
nRandSamples = 1;
coeffFun_cell = cell(nRandSamples*5,1);
coeffFun_counter = 1;
parameter_cell = cell(nRandSamples*5,4);

rhoBound = 10.^[0,3,6];

%% Konstante Koeffizientenfunktion
% Test verschiedene parameter fuer die Kanalfunktion

% Erstelle die Parameterstruktur
param_names = ["affectedSubdomains","rhoMin","rhoMax"];
affectedSubdomains = [6,10];
% Samples hier haendisch erstellen
parameter_cell = {affectedSubdomains; affectedSubdomains; affectedSubdomains};
parameter_cell = [parameter_cell, num2cell([rhoBound;circshift(rhoBound,1,2)])']';

sample_parameters = cell2struct(parameter_cell,param_names,1);

for sampleID = 1:length(sample_parameters)
    affectedSubdomains = sample_parameters(sampleID).affectedSubdomains;
    rhoMin  = sample_parameters(sampleID).rhoMin;
    rhoMax  = sample_parameters(sampleID).rhoMax;

    coeffFun_cell{coeffFun_counter} = @(vertices) coeffFun_subdomains(vertices(:,1),vertices(:,2),affectedSubdomains,vert__sd);
    parameter_cell{coeffFun_counter,1} = "Constant";
    parameter_cell{coeffFun_counter,2} = param_names;
    parameter_cell{coeffFun_counter,3} = [affectedSubdomains,rhoMin,rhoMax];
    parameter_cell{coeffFun_counter,4}  = sample_parameters(sampleID);
    coeffFun_counter = coeffFun_counter + 1;
end

%% Streifen Koeffizientenfunktion
% Test verschiedene parameter fuer die Kanalfunktion
yOffsetBound = -2:2;
widthBound = -2:2;
nStripsBound = 1:5;

% Erstelle die Parameterstruktur
param_names = ["yOffset","width","nStrips","rhoMin","rhoMax"];
sample_parameters = generateSampleParameters(nRandSamples,param_names,yOffsetBound,widthBound,nStripsBound,rhoBound,rhoBound);
% sample_parameters(6) = cell2struct(num2cell([0; 1; 2; 3; 4]),param_names,1);

for sampleID = 1:length(sample_parameters)
    sample = sample_parameters(sampleID);
    yOffset = sample_parameters(sampleID).yOffset;
    width   = sample_parameters(sampleID).width;
    nStrips = sample_parameters(sampleID).nStrips;
    rhoMin  = sample_parameters(sampleID).rhoMin;
    rhoMax  = sample_parameters(sampleID).rhoMax;

    coeffFun_cell{coeffFun_counter} = @(vertices) coeffFun_canal(vertices(:,2),N,n,yOffset,width,nStrips);
    parameter_cell{coeffFun_counter,1} = "Strip";
    parameter_cell{coeffFun_counter,2} = param_names;
    parameter_cell{coeffFun_counter,3} = [yOffset,width,nStrips,rhoMin,rhoMax];
    parameter_cell{coeffFun_counter,4}  = sample_parameters(sampleID);
    coeffFun_counter = coeffFun_counter + 1;
end

%% Bloecke Koeffizientenfunktion
% Teste verschiedene Parameter fuer die Kanalfunktion
difBound = -4:2:4;
prop1Bound = 0:0.25:1;
prop2Bound = 0:0.25:1;

% Erstelle die Parameterstruktur
param_names = ["yOffset","width","rhoMin","rhoMax","dif","prop1","prop2"];
sample_parameters = generateSampleParameters(nRandSamples,param_names,yOffsetBound,widthBound,rhoBound,rhoBound,difBound,prop1Bound,prop2Bound);

for sampleID = 1:length(sample_parameters)
    yOffset = sample_parameters(sampleID).yOffset;
    width   = sample_parameters(sampleID).width;
    dif     = sample_parameters(sampleID).dif;
    prop1   = sample_parameters(sampleID).prop1;
    prop2   = sample_parameters(sampleID).prop2;
    rhoMin  = sample_parameters(sampleID).rhoMin;
    rhoMax  = sample_parameters(sampleID).rhoMax;

    coeffFun_cell{coeffFun_counter} = @(vertices) coeffFun_block(vertices(:,1), vertices(:,2), N, n, prop1,prop2,dif,yOffset,width);
    parameter_cell{coeffFun_counter,1}  = "Blocks";
    parameter_cell{coeffFun_counter,2}  = param_names;
    parameter_cell{coeffFun_counter,3}  = [yOffset,width,rhoMin,rhoMax,dif,prop1,prop2];
    parameter_cell{coeffFun_counter,4}  = sample_parameters(sampleID);
    coeffFun_counter = coeffFun_counter + 1;
end

%% Zufalls - Bloecke Koeffizientenfunktion
% Test verschiedene parameter fuer die Kanalfunktion
widthBound      =  2: 1: 6; 
heightBound     =  2: 1: 6; 
varianceBound   =  0: 1: 5;
nBlocksBound    = 10:10:70;

% Erstelle die Parameterstruktur
param_names = ["nBlocks","height","heightVariance","width","widthVariance","rhoMin","rhoMax"];
sample_parameters = generateSampleParameters(nRandSamples,param_names,nBlocksBound,heightBound,varianceBound,widthBound,varianceBound,rhoBound,rhoBound);

for sampleID = 1:length(sample_parameters)
    nBlocks         = sample_parameters(sampleID).nBlocks;
    height          = sample_parameters(sampleID).height;
    heightVariance  = sample_parameters(sampleID).heightVariance;
    width           = sample_parameters(sampleID).width;
    widthVariance   = sample_parameters(sampleID).widthVariance;
    rhoMin          = sample_parameters(sampleID).rhoMin;
    rhoMax          = sample_parameters(sampleID).rhoMax;

    coeffFun_cell{coeffFun_counter} = @(vertices) coeffFun_randomBlocks(vertices(:,1),vertices(:,2),N,n,nBlocks,width:width+widthVariance,height:height+heightVariance);
    parameter_cell{coeffFun_counter,1}  = "Random Blocks";
    parameter_cell{coeffFun_counter,2}  = param_names;
    parameter_cell{coeffFun_counter,3}  = [nBlocks,height,heightVariance,width,widthVariance,rhoMin,rhoMax];
    parameter_cell{coeffFun_counter,4}  = sample_parameters(sampleID);
    coeffFun_counter = coeffFun_counter + 1;
end

%% Zufalls Koeffizientenfunktion
% Test verschiedene parameter fuer die Kanalfunktion
randomPercentageBound = 0.2:0.1:0.5;
randomStateBound = 1:5;

% Erstelle die Parameterstruktur
param_names = ["randomPercentage","randomState","rhoMin","rhoMax"];
sample_parameters = generateSampleParameters(nRandSamples,param_names,randomPercentageBound,randomStateBound,rhoBound,rhoBound);

for sampleID = 1:length(sample_parameters)
    randomPercentage    = sample_parameters(sampleID).randomPercentage;
    randomState         = sample_parameters(sampleID).randomState;
    rhoMin  = sample_parameters(sampleID).rhoMin;
    rhoMax  = sample_parameters(sampleID).rhoMax;

    coeffFun_cell{coeffFun_counter} = @(vertices) coeffFun_randomBlocks(vertices(:,1),vertices(:,2),N,n,nBlocks,width:width+widthVariance,height:height+heightVariance);
    parameter_cell{coeffFun_counter,1}  = "Completely Random";
    parameter_cell{coeffFun_counter,2}  = param_names;
    parameter_cell{coeffFun_counter,3}  = [randomPercentage,randomState,rhoMin,rhoMax];
    parameter_cell{coeffFun_counter,4}  = sample_parameters(sampleID);
    coeffFun_counter = coeffFun_counter + 1;
end

%% Faelle 
empty_cells_ind = cellfun('isempty',coeffFun_cell);
coeffFun_cell = coeffFun_cell(~empty_cells_ind);
parameter_cell = parameter_cell(~empty_cells_ind,:);
n_cases = length(coeffFun_cell);
output_cell = cell(n_cases,1);

for case_id = 1:n_cases
    fprintf("#### Starte Fall: Koeffizientenfunktion %s ####\n",parameter_cell{case_id})
    tic
    % Definiere Koeffizient auf den Elementen (und teilgebietsweise);
    % maximalen Koeffizienten pro Knoten (und teilgebietsweise)
    rhoMin = parameter_cell{case_id,4}.rhoMin;
    rhoMax = parameter_cell{case_id,4}.rhoMax;
    markerType = 'verts';
    [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = getCoefficientMatrices(coeffFun_cell{case_id},markerType,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid);
    fprintf("Benoetigte Zeit: Aufstellen der Koeffizientenmatrizen: %5fs ",toc)
    rho_struct = struct('rhoTriSD',{rhoTriSD},'maxRhoVert',{maxRhoVert},'maxRhoVertSD',{maxRhoVertSD});

    tic
    [edgesPrimalGlobal,cGamma,edgesSD,cLocalPrimal,cB,cBskal,cInner,cK,cDirichlet] = setup_matrices(rho_struct,grid_struct,f);
    fprintf(", Aufstellen des Sprungoperators/Steifigkeitsmatrix: %5fs\n", toc)

    numEdges = length(edgesSD);
    input = cell(numEdges,1);
    label = zeros(numEdges,1);

    for edgeID = 1:numEdges
        % Pruefe ob eines der beteiligten TG einen Dirichletknoten enthaelt.
        % Falls ja ist eins der TG kein floating TG und wird daher nicht
        % beruecksichtigt
        if (nnz(cDirichlet{edgesSD(edgeID,1)}) > 0) || (nnz(cDirichlet{edgesSD(edgeID,2)}) > 0)
            label(edgeID) = 2;
            continue
        else
            input{edgeID} = generate_input(edgeID,edgesSD,rhoTriSD,vert__sd,tri__sd);
            label(edgeID) = generate_label(edgeID,edgesPrimalGlobal,cGamma,edgesSD,cLocalPrimal,cB,cBskal,cInner,cK,TOL);
        end
        fprintf("Kante %2i bzgl. der TG (%2i,%2i) erhaelt das Label %i\n",edgeID,edgesSD(edgeID,1),edgesSD(edgeID,2),label(edgeID))
    end
    skipped_edges = nnz(label == 2);
    fprintf("Fuer das gegebene Gitter wurden %i (%4.1f%%) Kanten uebersprungen\n",skipped_edges,skipped_edges/numEdges*100)
    % Fuege neue Daten an den Trainingsdatensatz an
    output_cell{case_id} = [cell2mat(input),label(label ~= 2)];
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


function sample_parameters = generateSampleParameters(nRandSamples,param_names,varargin)
parameter_cell = cell(numel(varargin),1);
[parameter_cell{:}] = ndgrid(varargin{:});
parameter_cell = cellfun(@(X) reshape(X,1,[]),parameter_cell,'UniformOutput',false);
fprintf("Wahle zufaellig %i aus %i moeglichen Parameterkombinationen aus.\n",nRandSamples,length(parameter_cell));

random_permutation = randperm(size(parameter_cell{1},2));
parameter_mat = cell2mat(parameter_cell);
parameter_mat = parameter_mat(:,random_permutation(1:nRandSamples));
% parameter_cell = num2cell(parameter_mat);

sample_parameters = cell2struct(num2cell(parameter_mat),param_names,1);
end