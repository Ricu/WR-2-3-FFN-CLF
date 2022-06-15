clear; clc;
addpath('libs')
export = 0;
plot_grid = true;   % Auswahl: Plotten der Triangulierung mit Kanal-Koeffizientenfunktion

fprintf("############ Erstelle Testdaten Start ############\n")
fprintf("Startzeit %s\n", datestr(datetime))
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

%Parameter fuer die Kanalstruktur
position = -2:2;
width = -2:2;
number = 1:5;

%Parameter fuer die Blockstruktur
dif = -4:2:4;
prop1 = 0:0.25:1;
prop2 = 0:0.25:1;

%Parameter fuer Randomstruktur
numrand = 0:5:100;

%Beispielparameter zum Testen
positionC = position(1);
widthC = width(4);
numberC = number(1);
propB1 = prop1(2);
propB2 = prop2(2);
difB = dif(2);

numRand = numrand(15);

[rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_Rand(numRand,h,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid);


[rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_Block(propB1,propB2,difB,positionC,widthC,numberC,h,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid);

[rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_Canal(positionC,widthC,numberC,h,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid);
% Structure fuer rho-Variablen
rho_struct = struct('rhoTriSD',{rhoTriSD},'maxRhoVert',{maxRhoVert},'maxRhoVertSD',{maxRhoVertSD});
  
%% Definiere constraint-Typ
TOL = 100;  % Toleranz zur Auswahl der Eigenwerte

[edgesPrimalGlobal,cGamma,edgesSD,cLocalPrimal,cB,cBskal,cInner,cK,cDirichlet] = setup_matrices(rho_struct,grid_struct,f);

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

%% Daten exportieren
if export
    file_name = sprintf("%s-train_data_dump.csv",datestr(datetime,'yyyy-mm-dd-HH-MM-SS'));
    fprintf("Speichere Traininsdaten als %s...",file_name)
    output_mat = [cell2mat(input),label(label ~= 2)];
    writematrix(output_mat,file_name);
    fprintf("Fertig!\n")
end