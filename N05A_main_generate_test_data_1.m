clear; clc;
addpath('libs')
export = 1;         % Auswahl: Testdaten exporieren und abspeichern
plot_grid = true;   % Auswahl: Plotten der Triangulierung mit Kanal-Koeffizientenfunktion

if export
    fprintf("Die Daten werden gespeichert\n")
else
    fprintf("Die Daten werden nicht gespeichert\n")
end

%% Erstelle das Gitter
n = 40;         % 2*n^2 Elemente pro Teilgebiet
N = 4;          % Partition in NxN quadratische Teilgebiete
H = 1/N;        % Schrittweite: Teilgebiete
h = 1/(N*n);    % Schrittweite: Elemente
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

%% Koeffizientenfunktion aufstellen
% Definiere minimales und maximales rho
rhoMin = 1;
rhoMax = 10^6;

% Parameter fuer die vorgegebene Koeffizientenverteilung
yStripeLim = [0.1,0.9];
xDisplacement = 6;
width = 3;
hight = 5;

% Definiere zu testende Koeffizientenverteilung: 1 -  Hufeisen
coeffFun = @(tri) coeffFun_horseshoe(tri,vert(:,1),vert(:,2),N,n,yStripeLim,xDisplacement,width,hight);
markerType = 'elements';  % Die Koeffizientenverteilung ist elementweise definiert

% Erstelle die Koeffizientenmatrizen, welche in der Berechnung der
% FETI-DP Matrizen benoetigt werden
[~,rhoTriSD,maxRhoVert,maxRhoVertSD] = getCoefficientMatrices(coeffFun,markerType,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid,vertTris);
rho_struct = struct('rhoTriSD',{rhoTriSD},'maxRhoVert',{maxRhoVert},'maxRhoVertSD',{maxRhoVertSD});

%% Funktion rechte Seite
f = @(vert,y) ones(size(vert));   % Rechte Seite der DGL

%% Toleranz zur Auswahl der Eigenwerte
TOL = 100;  

%% Erstelle Testdaten
% Erstelle FETI-DP Matrizen, welche zur Berechnung der Label benoetigt werden
[edgesPrimalGlobal,cGamma,edgesSD,cLocalPrimal,cB,cBskal,cInner,cK,cDirichlet] = setup_matrices(rho_struct,grid_struct,f);

numEdges = length(edgesSD);
input_cell = cell(numEdges,1);
label = zeros(numEdges,1);
for edgeID = 1:numEdges
    % Input fuer neuronales Netz: Koeffizientenverteilung elementweise
    input_cell{edgeID} = generate_input(edgeID,edgesSD,rhoTriSD,vert__sd,tri__sd);
    % Zugehoeriger label: 1 kritische Kante, 0 unktritische Kante
    label(edgeID) = generate_label(edgeID,edgesPrimalGlobal,cGamma,edgesSD,cLocalPrimal,cB,cBskal,cInner,cK,TOL);
end
% Koeffizientenverteilung mit zugehoerigem label als Testdaten
input_mat = [cell2mat(input_cell),label];

%% Daten exportieren
if export
    file_name = sprintf("./resources/test_data/%s-test_data_1_dump.csv",datestr(datetime,'yyyy-mm-dd-HH-MM-SS'));
    fprintf("Speichere Testdaten als %s...",file_name)
    writematrix(input_mat,file_name);
    fprintf("Fertig!\n")
end