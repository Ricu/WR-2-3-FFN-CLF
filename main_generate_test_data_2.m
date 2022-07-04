clear; clc;
addpath('libs')
export = 1;         % Auswahl: Testdaten abspeichern
plot_grid = true;   % Auswahl: Plotten der Triangulierung mit Kanal-Koeffizientenfunktion

%% Lade vertTris fuer schnellere Berechnung der coeff-funktion
vertTris = load("./libs/precomputed_vertTris.mat").vertTris;

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

TOL = 100;  % Toleranz zur Auswahl der Eigenwerte

%% Koeffizientenfunktion aufstellen
%Bilddatei einlesen 
pic = imread('./resources/img/rho_coeff_multiple_stripes.png');
pic_bw = pic(:,:,1)'; %benoetigen nur einen Kanal, da schwarz-weiss Bild
num_pixel = length(pic_bw); %Anzahl Pixel je Dimension

% Definiere minimales und maximales rho
rhoMin = 1;
rhoMax = 10^6;

% Aufstellen der zu testenden Koeffizientenverteilung: 2 -  Streifen
coeffFun = @(vert) coeffFun_image(vert(:,1),vert(:,2),pic_bw,num_pixel);
markerType = 'verts';  % Die Koeffizientenverteilung ist knotenweise definiert

% Definiere Koeffizientenfunktion auf den Elementen (teilgebietsweise);
% maximalen Koeffizienten pro Knoten (und teilgebietsweise)
[~,rhoTriSD,maxRhoVert,maxRhoVertSD] = getCoefficientMatrices(coeffFun,markerType,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid,vertTris);
rho_struct = struct('rhoTriSD',{rhoTriSD},'maxRhoVert',{maxRhoVert},'maxRhoVertSD',{maxRhoVertSD});

%% Funktion rechte Seite
f = @(vert,y) ones(size(vert));   % Rechte Seite der DGL

%% Benoetigte Matrizen aufstellen
[edgesPrimalGlobal,cGamma,edgesSD,cLocalPrimal,cB,cBskal,cInner,cK,cDirichlet] = setup_matrices(rho_struct,grid_struct,f);

%% Inputgenerierung
numEdges = length(edgesSD);
input_cell = cell(numEdges,1);
label = zeros(numEdges,1);
for i = 1:numEdges
    edgeID = edgesSD(i);
    input_cell{i} = generate_input(edgeID,edgesSD,rhoTriSD,vert__sd,tri__sd);
    label(i) = generate_label(edgeID,edgesPrimalGlobal,cGamma,edgesSD,cLocalPrimal,cB,cBskal,cInner,cK,TOL);
end
input_mat = [cell2mat(input_cell),label];

%% Daten exportieren
if export
    file_name = sprintf("./resources/test_data/%s-test_data_2_dump.csv",datestr(datetime,'yyyy-mm-dd-HH-MM-SS'));
    fprintf("Speichere Testdaten als %s...",file_name)
    writematrix(input_mat,file_name);
    fprintf("Fertig!\n")
end
