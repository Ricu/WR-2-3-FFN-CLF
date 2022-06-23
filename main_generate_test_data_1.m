clear; clc;
addpath('libs')
export = 1;
plot_grid = true;   % Auswahl: Plotten der Triangulierung mit Kanal-Koeffizientenfunktion

%% Erstelle das Gitter
n = 40;         % 2*n^2 Elemente pro Teilgebiet
N = 4;          % Partition in NxN quadratische Teilgebiete
H = 1/N;
h = 1/(N*n);
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
% Definiere minimales und maximales rho
rhoMin = 1;
rhoMax = 10^6;

%Parameter fuer Hufeisen/Streifen Triangulierung
yStripeLim = [0.1,0.9];
position = 2;
width = 1;
hight = 2;

% Definiere zu testende Koeffizientenverteilung: 1 -  Hufeisen
coeffFun = @(vertices) coeffFun_horseshoe(vert,vert(:,1),vert(:,2),N,n,yStripeLim,position,width,hight);

% Definiere Koeffizient auf den Elementen (und teilgebietsweise);
% maximalen Koeffizienten pro Knoten (und teilgebietsweise)
[rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = getCoefficientMatrices(coeffFun,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid);
rho_struct = struct('rhoTriSD',{rhoTriSD},'maxRhoVert',{maxRhoVert},'maxRhoVertSD',{maxRhoVertSD});

%% Funktion rechte Seite
f = @(vert,y) ones(size(vert));   % Rechte Seite der DGL

%% Kantenliste erstellen
[~,~,edgesSD,~,~,~,~,~,cDirichlet] = setup_matrices(rho_struct,grid_struct,f);
numEdges = length(edgesSD);

%% Inputgenerierung
input_cell = cell(numEdges,1);
for edgeID = 1:numEdges
    % Pruefe ob eines der beteiligten TG einen Dirichletknoten enthaelt.
    % Falls ja ist eins der TG kein floating TG und wird daher nicht
    % beruecksichtigt
    if (nnz(cDirichlet{edgesSD(edgeID,1)}) > 0) || (nnz(cDirichlet{edgesSD(edgeID,2)}) > 0)
        continue
    else
        input_cell{edgeID} = generate_input(edgeID,edgesSD,rhoTriSD,vert__sd,tri__sd);
    end
end
% Fuege neue Daten an den Testdatensatz an
input_mat = cell2mat(input_cell);


%% Daten exportieren
if export
    file_name = sprintf("./resources/test_data/%s-test_data_dump.csv",datestr(datetime,'yyyy-mm-dd-HH-MM-SS'));
    fprintf("Speichere Traininsdaten als %s...",file_name)
    writematrix(input_mat,file_name);
    fprintf("Fertig!\n")
end