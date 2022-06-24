clear; clc;
addpath('libs')
export = 0;         % Auswahl: Testdaten abspeichern
plot_grid = true;   % Auswahl: Plotten der Triangulierung mit Kanal-Koeffizientenfunktion

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

%% Koeffizientenfunktion aufstellen
%Bilddatei einlesen 
pic = imread('./resources/img/multiple_circle_bw_512x512px.jpeg');
pic_bw = pic(:,:,1)'; % Nur ein Kanal, da schwarz-weiss Bild
num_pixel = length(pic_bw); % Anzahl Pixel je Dimension

% Definiere minimales und maximales rho
rhoMin = 1;
rhoMax = 10^6;

% Aufstellen der zu testenden Koeffizientenverteilung: 3 -  Kreise
coeffFun = @(tri) coeffFun_pixel(vert(:,1),vert(:,2),tri,rhoMax,rhoMin,pic_bw,num_pixel);
base = 'elements';  % Die Koeffizientenverteilung ist elementweise definiert

% Definiere Koeffizientenfunktion auf den Elementen (teilgebietsweise);
% maximalen Koeffizienten pro Knoten (und teilgebietsweise)
[~,rhoTriSD,maxRhoVert,maxRhoVertSD] = getCoefficientMatrices(coeffFun,base,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid);
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
input_mat = cell2mat(input_cell);

%% Daten exportieren
if export
    file_name = sprintf("./resources/test_data/%s-test_data_3_dump.csv",datestr(datetime,'yyyy-mm-dd-HH-MM-SS'));
    fprintf("Speichere Testdaten als %s...",file_name)
    writematrix(input_mat,file_name);
    fprintf("Fertig!\n")
end