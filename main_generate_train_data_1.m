clear; clc;
addpath('libs')
export = 0;
plot_grid = 0;   % Auswahl: Plotten der Triangulierung mit Kanal-Koeffizientenfunktion

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

%% Koeffizientenfunktion aufstellen
% Definiere rho im Kanal und außerhalb des Kanals
rhoMax = 10^6;
rhoMin = 1;


% Definiere Koeffizient auf den Elementen (und teilgebietsweise);
% maximalen Koeffizienten pro Knoten (und teilgebietsweise)


TOL = 100;  % Toleranz zur Auswahl der Eigenwerte
rng(0);
nRandSamples = 2;
output_cell = cell(nRandSamples*5,1);
output_counter = 1;
rhoBound = 10.^(0:7);

%% Kanal Koeffizientenfunktion
% Test verschiedene parameter fuer die Kanalfunktion
positionBound = -2:2;
widthBound = -2:2;
numberBound = 1:5;

% Erstelle eine zufaellige Auswahl an Parametern
parameter_cell = {positionBound,widthBound,numberBound,rhoBound,rhoBound};
parameter_grid = cell(numel(parameter_cell),1);
[parameter_grid{:}] = ndgrid(parameter_cell{:});
parameter_grid = cellfun(@(X) reshape(X,1,[]),parameter_grid,'UniformOutput',false);
parameter_grid = cell2mat(parameter_grid);
random_permutation = randperm(size(parameter_grid,2));
sample_parameters = parameter_grid(:,random_permutation(1:nRandSamples));

% Erstelle die Parametervektoren
position_vec = [-2, sample_parameters(1,:)];
width_vec = [1, sample_parameters(2,:)];
number_vec = [1, sample_parameters(3,:)];
rhoMin_vec = [1, sample_parameters(4,:)];
rhoMax_vec = [10^6, sample_parameters(5,:)];



for sampleID = 1:length(position_vec)
    position = position_vec(sampleID);
    width = width_vec(sampleID);
    number = number_vec(sampleID);
    rhoMin = rhoMin_vec(sampleID);
    rhoMax = rhoMax_vec(sampleID);

    fprintf("#### Starte Durchlauf: Koeffizientenfunktion Kanal mit position=%2i, width=%2i, number=%2i, rhoMin=%7i, rhoMax=%7i\n",position,width,number,rhoMin,rhoMax)
    tic
    [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_Canal(position,width,number,h,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid);
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
            input{edgeID} = generate_input(edgeID,edgesSD,maxRhoVert,l2g__sd);
            label(edgeID) = generate_label(edgeID,edgesPrimalGlobal,cGamma,edgesSD,cLocalPrimal,cB,cBskal,cInner,cK,TOL);
        end
        fprintf("Kante %2i bzgl. der TG (%2i,%2i) erhaelt das Label %i\n",edgeID,edgesSD(edgeID,1),edgesSD(edgeID,2),label(edgeID))
    end
    skipped_edges = nnz(label == 2);
    fprintf("Fuer das gegebene Gitter wurden %i (%4.1f%%) Kanten uebersprungen\n",skipped_edges,skipped_edges/numEdges*100)
    % Fuege neue Daten an den trainingsdatensatz an
    output_cell{output_counter} = [cell2mat(input),label(label ~= 2)];
    output_counter = output_counter + 1;
end

%% Blockstruktur Koeffizientenfunktion
% Test verschiedene parameter fuer die Kanalfunktion
difBound = -4:2:4;
prop1Bound = 0:0.25:1;
prop2Bound = 0:0.25:1;

% Erstelle eine zufaellige Auswahl an Parametern
parameter_cell = {positionBound,widthBound,numberBound,rhoBound,rhoBound,difBound,prop1Bound,prop2Bound};
parameter_grid = cell(numel(parameter_cell),1);
[parameter_grid{:}] = ndgrid(parameter_cell{:});
parameter_grid = cellfun(@(X) reshape(X,1,[]),parameter_grid,'UniformOutput',false);
parameter_grid = cell2mat(parameter_grid);
random_permutation = randperm(size(parameter_grid,2));
sample_parameters = parameter_grid(:,random_permutation(1:nRandSamples));

% Erstelle die Parametervektoren
position_vec = [-2, sample_parameters(1,:)];
width_vec = [1, sample_parameters(2,:)];
number_vec = [1, sample_parameters(3,:)];
rhoMin_vec = [1, sample_parameters(4,:)];
rhoMax_vec = [10^6, sample_parameters(5,:)];
dif_vec = [-3, sample_parameters(6,:)];
prop1_vec = [0.25, sample_parameters(7,:)];
prop2_vec = [0.25, sample_parameters(8,:)];

for sampleID = 1:length(position_vec)
    position = position_vec(sampleID);
    width = width_vec(sampleID);
    number = number_vec(sampleID);
    rhoMin = rhoMin_vec(sampleID);
    rhoMax = rhoMax_vec(sampleID);
    dif = dif_vec(sampleID);
    prop1 = prop1_vec(sampleID);
    prop2 = prop2_vec(sampleID);

    fprintf("#### Starte Durchlauf: Koeffizientenfunktion Block mit position=%2i, width=%2i, number=%2i, rhoMin=%7i, rhoMax=%7i\n",position,width,number,rhoMin,rhoMax)
    tic
    [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_Block(prop1,prop2,dif,position,width,number,h,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid);
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
            input{edgeID} = generate_input(edgeID,edgesSD,maxRhoVert,l2g__sd);
            label(edgeID) = generate_label(edgeID,edgesPrimalGlobal,cGamma,edgesSD,cLocalPrimal,cB,cBskal,cInner,cK,TOL);
        end
        fprintf("Kante %2i bzgl. der TG (%2i,%2i) erhaelt das Label %i\n",edgeID,edgesSD(edgeID,1),edgesSD(edgeID,2),label(edgeID))
    end
    skipped_edges = nnz(label == 2);
    fprintf("Fuer das gegebene Gitter wurden %i (%4.1f%%) Kanten uebersprungen\n",skipped_edges,skipped_edges/numEdges*100)
    % Fuege neue Daten an den trainingsdatensatz an
    output_cell{output_counter} = [cell2mat(input),label(label ~= 2)];
    output_counter = output_counter + 1;
end


%% Zufalls Koeffizientenfunktion
% Test verschiedene parameter fuer die Kanalfunktion
widthBound = 1:5; % Noch unbenutzt
heightBound = 1:5; % Noch unbenutzt
numberBound = 10:10:70;

% Erstelle eine zufaellige Auswahl an Parametern
parameter_cell = {heightBound,widthBound,numberBound,rhoBound,rhoBound};
parameter_grid = cell(numel(parameter_cell),1);
[parameter_grid{:}] = ndgrid(parameter_cell{:});
parameter_grid = cellfun(@(X) reshape(X,1,[]),parameter_grid,'UniformOutput',false);
parameter_grid = cell2mat(parameter_grid);
random_permutation = randperm(size(parameter_grid,2));
sample_parameters = parameter_grid(:,random_permutation(1:nRandSamples));

% Erstelle die Parametervektoren
height_vec = [-2, sample_parameters(1,:)];
width_vec = [1, sample_parameters(2,:)];
number_vec = [16, sample_parameters(3,:)];
rhoMin_vec = [1, sample_parameters(4,:)];
rhoMax_vec = [10^6, sample_parameters(5,:)];

for sampleID = 1:length(height_vec)
    height = height_vec(sampleID);
    width = width_vec(sampleID);
    number = number_vec(sampleID);
    rhoMin = rhoMin_vec(sampleID);
    rhoMax = rhoMax_vec(sampleID);

    fprintf("#### Starte Durchlauf: Koeffizientenfunktion Zufall mit position=%2i, width=%2i, number=%2i, rhoMin=%7i, rhoMax=%7i\n",position,width,number,rhoMin,rhoMax)
    tic
    [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_Rand(number,h,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid);
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
            input{edgeID} = generate_input(edgeID,edgesSD,maxRhoVert,l2g__sd);
            label(edgeID) = generate_label(edgeID,edgesPrimalGlobal,cGamma,edgesSD,cLocalPrimal,cB,cBskal,cInner,cK,TOL);
        end
        fprintf("Kante %2i bzgl. der TG (%2i,%2i) erhaelt das Label %i\n",edgeID,edgesSD(edgeID,1),edgesSD(edgeID,2),label(edgeID))
    end
    skipped_edges = nnz(label == 2);
    fprintf("Fuer das gegebene Gitter wurden %i (%4.1f%%) Kanten uebersprungen\n",skipped_edges,skipped_edges/numEdges*100)
    % Fuege neue Daten an den trainingsdatensatz an
    output_cell{output_counter} = [cell2mat(input),label(label ~= 2)];
    output_counter = output_counter + 1;
end





%% Daten exportieren
if export
    file_name = sprintf("%s-train_data_dump.csv",datestr(datetime,'yyyy-mm-dd-HH-MM-SS'));
    fprintf("Speichere Traininsdaten als %s...",file_name)
    output_mat = cell2mat(output_cell);
    writematrix(output_mat,file_name);
    fprintf("Fertig!\n")
end