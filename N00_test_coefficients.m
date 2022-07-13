% In dieser Datei koennen die verschiedenen Koeffizientenfunktionen
% ausprobiert und geplottet werden.
addpath('libs')
clc; clear;
plot_grid = 1;

%% Waehle die gewünschte Koeffizientenfunktion
coeff_subdomain = false;
coeff_strip = false;
coeff_block = false;
coeff_randomBlocks = false;
coeff_stripRandomBlocks = true;
coeff_random = false;

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


%% Konstante Koeffizientenfunktion auf den Teilgebieten
if coeff_subdomain == true
    affectedSubdomains = [6,10]; % TG-Wahl fuer hoeheren Koeffizienten; decken alle moeglichen Faelle ab
    indexShiftx = 0;
    indexShifty = 0;

    coeffFun= @(vertices) coeffFun_subdomains(vertices(:,1),vertices(:,2),affectedSubdomains,vert__sd,indexShiftx,indexShifty);
end

%% Streifen Koeffizientenfunktion
% widthBound = -2:4; % Breite der Streifen, 0 ist initiale Breite abhaengig von der Anzahl Streifen je TG
% nStripsBound = 1:5; % Anzahl Streifen je TG
if coeff_strip == true
    width = -2;
    nStrips = 1;
    indexShifty = 0;
    
    coeffFun = @(vertices) coeffFun_canal(vertices(:,2),N,n,width,nStrips,indexShifty);
end

%% Bloecke Koeffizientenfunktion
% difBound = -20:2:20;  % Gibt an, wie weit die Bloecke in jedem zweiten TG (spaltenweise ab 2.Spalte) 
%                       % voneinander versetzt sind. 0 entspricht keiner Versetzung
% prop1Bound = 0:0.2:1; % Gibt den Anteil an Block in jedem zweiten TG (spaltenweise ab 1.Spalte) an
% prop2Bound = 0:0.2:1; % Gibt den Anteil an Block in jedem zweiten TG (spaltenweise ab 2.Spalte) an
% heightBound = 2:2:38; % Hoehe der Bloecke

if coeff_block == true
    height   = 4;
    dif     = 20;
    prop1   = 0.5;
    prop2   = 0.5;
    indexShiftx = 0;
    indexShifty = 0;
    
    coeffFun = @(vertices) coeffFun_block(vertices(:,1), vertices(:,2), N, n, prop1,prop2,dif,height,indexShiftx,indexShifty);
end

%% Zufalls - Bloecke Koeffizientenfunktion
% widthBound      =  4:4:16;  % Breite der Bloecke mit Faktor der Schrittweite
% heightBound     =  4:4:16;  % Hoehe der Bloecke mit Faktor der Schrittweite
% varianceBound   =  0:4:16;  % positive Varianz in Breite und Hoehe
% nBlocksBound    =  10:5:20; % Anzahl an random erstellten Bloecken
if coeff_randomBlocks == true
    nBlocks         = 40;
    height          =  4;
    heightVariance  = 10;
    width           =  4;
    widthVariance   = 10;
    indexShiftx = 0;
    indexShifty = 0;

    coeffFun = @(vertices) coeffFun_randomBlocks(vertices(:,1),vertices(:,2),N,n,nBlocks,width:width+widthVariance,height:height+heightVariance,indexShiftx,indexShifty);
end

%% Streifen + Zufalls - Bloecke Koeffizientenfunktion
% widthBound      =  1:1:3; % Breite der Bloecke mit Faktor der Schrittweite
% heightBound     =  2:1:4; % Hoehe der Bloecke mit Faktor der Schrittweite
% varianceBound   =  0:1:2;  % positive Varianz in Breite und Hoehe
% nBlocksBound    = 10:10:50; % Anzahl an random erstellten Bloecken
% heightBoundS    = -2:2;  % Breite der Streifen, 0 ist dabei eine initiale Breite abhaengig von der Anzahl an Streifen je TG
% nStripBound     = 1:4; % Gibt die Anzahl Streifen je TG an
if coeff_stripRandomBlocks == true
    heightS = 1;
    nStrip = 2;
    nBlocks = 20;
    width = 2;
    widthVariance = 1;
    height = 3;
    heightVariance = 1;
    indexShifty = 1;

    coeffFun = @(vertices) coeffFun_stripRandomBlocks(vertices(:,1),vertices(:,2),N,n,heightS,nStrip,nBlocks,width:width+widthVariance,height:height+heightVariance,indexShifty);
end       

%% Zufalls Koeffizientenfunktion
% randomPercentageBound = 0.2:0.05:0.7;
% randomStateBound = 1:10;
if coeff_random == true
    randomPercentage    = 0.7;
    randomState         = 1;
    indexShiftx         = 0;
    indexShifty         = 0;
    
    coeffFun = @(vertices) coeffFun_random(vertices(:,1),vertices(:,2),randomPercentage,randomState,indexShiftx,indexShifty);
end

%% Plot

rhoMin = 1;
rhoMax = 10^6;
markerType = 'verts';
vertTris = load("./libs/precomputed_vertTris.mat").vertTris;
getCoefficientMatrices(coeffFun,markerType,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid,vertTris);
