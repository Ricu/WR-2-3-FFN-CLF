clear; clc;
addpath('libs')
plot_grid = true;   % Auswahl: Plotten der Triangulierung mit Kanal-Koeffizientenfunktion

%% Daten importieren
file_name = sprintf("./resources/trained_model/predicted_labels_2.csv");
fprintf("Lese Testdaten aus %s...",file_name)
predicted_labels = readmatrix(file_name);

%% Definiere zu vergleichende Verfahren
method_type = {'Dirichlet','none';
               'Balancing','non-adaptive'
               'Balancing','adaptive';
               'Balancing','adaptive-improved'
               };       
 numMethods = length(method_type);

%% Initialisiere Parameter fuer PCG
x0 = @(dim) zeros(dim,1);    % Startvektor
tol = 10^(-8);               % Toleranz fuer die Abbruchbedingung

% Residuum fuer die Abbruchbedingung
resid_type = {'vorkonditioniert'};

% Structure fuer PCG-Parameter
pcg_param = struct('tol', tol, 'x0',x0, 'resid_type',resid_type);

%% Funktion rechte Seite
f = @(vert,y) ones(size(vert));   % Rechte Seite der DGL

%% Lade vertTris fuer schnellere Berechnung der coeff-funktion
vertTris = load("./libs/precomputed_vertTris.mat").vertTris;

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

%% Loesen des Systems mit FETI-DP fuer versch. Verfahren
% Referenzloesung auftsllen und vergleichen? 
% diffs = cell(numMethods,1);
iters = cell(numMethods,1);
kappa_ests = cell(numMethods,1);

plot_iteration = false; % Auswahl: Plotten der Loesung nach den ersten Iterationen von PCG

TOL = 100;  % Toleranz zur Auswahl der Eigenwerte bei adaptive 

for m = 1:numMethods
    VK = method_type{m,1};
    constraint_type = method_type{m,2};
    pc_param = struct('VK',VK,'constraint_type',constraint_type,'adaptiveTol',TOL);
    
    % Loesen des Systems mit FETI-DP mit entsprechendem VK
    [cu,u_FETIDP_glob,~,iters{m},kappa_ests{m}] = fetidp(grid_struct,f,pc_param,rho_struct,pcg_param,plot_iteration,predicted_labels);
end

%% Ergebnistabelle
%rowNames = ["Anzahl Iterationen","Konditionszahl","Abweichung von Referenzloesung"];
rowNames = ["Anzahl Iterationen","Konditionszahl"];
variableNames = (method_type(:,2));
T_results = cell2table([iters';kappa_ests'],"RowNames",rowNames,"VariableNames",variableNames);
disp(T_results)