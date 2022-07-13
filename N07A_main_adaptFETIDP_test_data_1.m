clear; clc;
addpath('libs')
plot_grid = false;   % Auswahl: Plotten der Triangulierung mit Kanal-Koeffizientenfunktion

%% Daten importieren
predicted_labels_loc = "./resources/trained_model/predicted_labels_1.csv";
fprintf("Lese predicted lables aus %s...",predicted_labels_loc)
predicted_labels = readmatrix(predicted_labels_loc);
fprintf("Fertig \n")
true_labels_loc = "./resources/test_data/test_data_1_dump.csv";
fprintf("Lese wahre labels aus %s...",true_labels_loc)
true_labels = readmatrix(true_labels_loc);
true_labels = true_labels(:,end);
fprintf("Fertig \n")

%% Nuetzliche Ausgaben
fprintf("Kante    : ")
fprintf("%2i,", 1:length(true_labels)-1); fprintf("%2i\n",length(true_labels));
fprintf("predicted: ")
fprintf("%2i,", predicted_labels(1:end-1)); fprintf("%2i\n",predicted_labels(end));
fprintf("true     : ")
fprintf("%2i,", true_labels(1:end-1)); fprintf("%2i\n",true_labels(end));

%% Definiere zu vergleichende Verfahren
method_type = {'Dirichlet','none';
               'Balancing','non-adaptive'
               'Balancing','adaptive';
               'Balancing','adaptive-improved';
               };       
 numMethods = length(method_type);

%% Initialisiere Parameter fuer PCG
x0 = @(dim) zeros(dim,1);           % Startvektor
tol = 10^(-8);                      % Toleranz fuer die Abbruchbedingung
resid_type = {'vorkonditioniert'};  % Residuum fuer die Abbruchbedingung

% Structure fuer PCG-Parameter
pcg_param = struct('tol', tol, 'x0',x0, 'resid_type',resid_type);

plot_iteration = false; % Auswahl: Plotten der Loesung nach den ersten Iterationen von PCG

%% Funktion rechte Seite
f = @(vert,y) ones(size(vert));   % Rechte Seite der DGL

%% Toleranz zur Auswahl der Eigenwerte bei adaptive 
TOL = 100;  

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

%Parameter fuer die vorgegebene Koeffizientenverteilung
yStripeLim = [0.1,0.9];
position = 6;
width = 3;
hight = 5;

% Definiere zu testende Koeffizientenverteilung: 1 -  Hufeisen
coeffFun = @(tri) coeffFun_horseshoe(tri,vert(:,1),vert(:,2),N,n,yStripeLim,position,width,hight);
markerType = 'elements';  % Die Koeffizientenverteilung ist elementweise definiert

% Erstelle die Koeffizientenmatrizen, welche in der Berechnung der
% FETI-DP Matrizen benoetigt werden
[rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = getCoefficientMatrices(coeffFun,markerType,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid);
rho_struct = struct('rhoTriSD',{rhoTriSD},'maxRhoVert',{maxRhoVert},'maxRhoVertSD',{maxRhoVertSD});

%% Plot der klassifizierten Kanten
plot_classified_edges(coeffFun,markerType,rhoMax,rhoMin,vert,tri,predicted_labels,true_labels,vert__sd)

%% Loesen des Systems mit FETI-DP fuer versch. Verfahren
iters = cell(numMethods,1);
kappa_ests = cell(numMethods,1);
for m = 1:numMethods
    VK = method_type{m,1};
    constraint_type = method_type{m,2};
    pc_param = struct('VK',VK,'constraint_type',constraint_type,'adaptiveTol',TOL);
    
    % Loesen des Systems mit FETI-DP mit entsprechendem VK
    [cu,u_FETIDP_glob,~,iters{m},kappa_ests{m}] = fetidp(grid_struct,f,pc_param,rho_struct,pcg_param,plot_iteration,predicted_labels);
end

%% Ergebnistabelle
rowNames = ["Anzahl Iterationen","Konditionszahl"];
variableNames = (method_type(:,2));
T_results = cell2table([iters';kappa_ests'],"RowNames",rowNames,"VariableNames",variableNames);
disp(T_results)