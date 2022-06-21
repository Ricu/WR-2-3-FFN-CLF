clear; clc;
addpath('libs')
export = 0;
plot_grid = 1;   % Auswahl: Plotten der Triangulierung mit Kanal-Koeffizientenfunktion

fprintf("############ Erstelle Testdaten Start ############\n")
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
parameter_cell = cell(nRandSamples*5,3);

rhoBound = 10.^[0,3,6];

%% Streifen Koeffizientenfunktion
% Test verschiedene parameter fuer die Kanalfunktion
yOffsetBound = -2:2;
widthBound = -2:2;
nStripsBound = 1:5;

% Erstelle die Parametervektoren
sample_parameters = generateSampleParameters(nRandSamples,yOffsetBound,widthBound,nStripsBound,rhoBound,rhoBound);
param_names = ["yOffset","width","nStrips","rhoMin","rhoMax"];
position_vec = [-2, sample_parameters(1,:)];
width_vec = [1, sample_parameters(2,:)];
nStrips_vec = [1, sample_parameters(3,:)];
rhoMin_vec = [1, sample_parameters(4,:)];
rhoMax_vec = [10^6, sample_parameters(5,:)];

for sampleID = 1:length(position_vec)
    position = position_vec(sampleID);
    width = width_vec(sampleID);
    nStrips = nStrips_vec(sampleID);
    rhoMin = rhoMin_vec(sampleID);
    rhoMax = rhoMax_vec(sampleID);

    
    coeffFun_cell{coeffFun_counter} = @(vertices) coeffFun_canal(vertices(:,2),N,n,position,width,nStrips);
    parameter_cell{coeffFun_counter,1} = "Strip";
    parameter_cell{coeffFun_counter,2} = param_names;
    parameter_cell{coeffFun_counter,3} = sample_parameters(:,sampleID)';
    coeffFun_counter = coeffFun_counter + 1;
end

%% Bloecke Koeffizientenfunktion
% Teste verschiedene Parameter fuer die Kanalfunktion
difBound = -4:2:4;
prop1Bound = 0:0.25:1;
prop2Bound = 0:0.25:1;

% Erstelle die Parametervektoren
sample_parameters = generateSampleParameters(nRandSamples,yOffsetBound,widthBound,nStripsBound,rhoBound,rhoBound,difBound,prop1Bound,prop2Bound);
position_vec = [-2, sample_parameters(1,:)];
width_vec = [1, sample_parameters(2,:)];
nStrips_vec = [1, sample_parameters(3,:)];
rhoMin_vec = [1, sample_parameters(4,:)];
rhoMax_vec = [10^6, sample_parameters(5,:)];
dif_vec = [-3, sample_parameters(6,:)];
prop1_vec = [0.25, sample_parameters(7,:)];
prop2_vec = [0.25, sample_parameters(8,:)];

for sampleID = 1:length(position_vec)
    position = position_vec(sampleID);
    width = width_vec(sampleID);
    nStrips = nStrips_vec(sampleID);
    rhoMin = rhoMin_vec(sampleID);
    rhoMax = rhoMax_vec(sampleID);
    dif = dif_vec(sampleID);
    prop1 = prop1_vec(sampleID);
    prop2 = prop2_vec(sampleID);

    coeffFun_cell{coeffFun_counter} = @(vertices) coeffFun_block(vertices(:,1), vertices(:,2), N, n, prop1, prop2, dif, position, width, nStrips);
    parameter_cell{coeffFun_counter} = "Blocks";
    coeffFun_counter = coeffFun_counter + 1;
end

%% Zufalls - Bloecke Koeffizientenfunktion
% Test verschiedene parameter fuer die Kanalfunktion
widthBound      =  2: 1: 6; 
heightBound     =  2: 1: 6; 
varianceBound   =  0: 1: 5;
nBlocksBound    = 10:10:70;

% Erstelle die Parametervektoren
sample_parameters = generateSampleParameters(nRandSamples,nBlocksBound,heightBound,varianceBound,widthBound,varianceBound,rhoBound,rhoBound);
nBlocks_vec =           [16,    sample_parameters(1,:)];
height_vec =            [5,     sample_parameters(2,:)];
heightVariance_vec =    [0,     sample_parameters(3,:)];
width_vec =             [5,     sample_parameters(4,:)];
widthVariance_vec =     [1,     sample_parameters(5,:)];
rhoMin_vec =            [1,     sample_parameters(6,:)];
rhoMax_vec =            [10^6,  sample_parameters(7,:)];

for sampleID = 1:length(height_vec)
    nBlocks          = nBlocks_vec(sampleID);
    height          = height_vec(sampleID);
    heightVariance  = heightVariance_vec(sampleID);
    width           = width_vec(sampleID);
    widthVariance   = widthVariance_vec(sampleID);
    rhoMin          = rhoMin_vec(sampleID);
    rhoMax          = rhoMax_vec(sampleID);

    coeffFun_cell{coeffFun_counter} = @(vertices) coeffFun_randomBlocks(vertices(:,1),vertices(:,2),N,n,nBlocks,width:width+widthVariance,height:height+heightVariance);
    parameter_cell{coeffFun_counter} = "Blocks";
    coeffFun_counter = coeffFun_counter + 1;
end


%% Faelle 
empty_cells_ind = cellfun('isempty',coeffFun_cell);
coeffFun_cell = coeffFun_cell(~empty_cells_ind);
parameter_cell = parameter_cell(~empty_cells_ind);
n_cases = length(coeffFun_cell);
output_cell = cell(n_cases,1);

for case_id = 1:n_cases
    fprintf("#### Starte Fall: Koeffizientenfunktion %s ####\n",parameter_cell{case_id})
    tic
    % Definiere Koeffizient auf den Elementen (und teilgebietsweise);
    % maximalen Koeffizienten pro Knoten (und teilgebietsweise)
    [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = getCoefficientMatrices(coeffFun_cell{case_id},rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid);
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
    % Fuege neue Daten an den trainingsdatensatz an
    output_cell{case_id} = [cell2mat(input),label(label ~= 2)];
end

%% Daten exportieren
if export
    file_name = sprintf("./resources/train_data/%s-train_data_dump.csv",datestr(datetime,'yyyy-mm-dd-HH-MM-SS'));
    fprintf("Speichere Traininsdaten als %s...",file_name)
    output_mat = cell2mat(output_cell);
    writematrix(output_mat,file_name);
    fprintf("Fertig!\n")
end


function sample_parameters = generateSampleParameters(nRandSamples,param_names,varargin)
% varNames = cell(1,size(varargin,2)-1);
% for i = 2:length(varargin)+1
%     varNames{i-1} = erase(inputname(i),"Bound");
% end

parameter_grid = cell(numel(varargin),1);
[parameter_grid{:}] = ndgrid(varargin{:});
parameter_grid = cellfun(@(X) reshape(X,1,[]),parameter_grid,'UniformOutput',false);
sample_parameters = cell2struct(parameter_grid,param_names,1);
parameter_grid = cell2mat(parameter_grid);
fprintf("Wahle %i aus %i moeglichen Parameterkombinationen aus\n",nRandSamples,length(parameter_grid));
random_permutation = randperm(size(parameter_grid,2));
sample_parameters = parameter_grid(:,random_permutation(1:nRandSamples));
end