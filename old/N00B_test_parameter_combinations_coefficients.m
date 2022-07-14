% In dieser Datei koennen die verschiedenen Koeffizientenfunktionen
% ausprobiert und geplottet werden.
addpath('libs')
clc; clear;
plot_grid = 1;

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


%% Konstante Koeffizientenfunktion
% affectedSubdomains = [6,10];
% indexShiftx = 10;
% indexShifty = 2;

% coeffFun= @(vertices) coeffFun_subdomains(vertices(:,1),vertices(:,2),affectedSubdomains,vert__sd,indexShiftx,indexShifty);

%% Streifen Koeffizientenfunktion
% widthBound = -2:4;
% nStripsBound = 1:5;

% width   = -2;
% nStrips = 1;
% indexShifty = 0;
% 
% coeffFun = @(vertices) coeffFun_canal(vertices(:,2),N,n,width,nStrips,indexShifty);

%% Bloecke Koeffizientenfunktion
% heightBound = 2:2:38;
% difBound = -20:2:20;
% prop1Bound = 0:0.2:1;
% prop2Bound = 0:0.2:1;

% height   = 4;
% dif     = 20;
% prop1   = 0.5;
% prop2   = 0.5;
% indexShiftx = 0;
% indexShifty = 0;
% 
% coeffFun = @(vertices) coeffFun_block(vertices(:,1), vertices(:,2), N, n, prop1,prop2,dif,height,indexShiftx,indexShifty);

%% Zufalls - Bloecke Koeffizientenfunktion
% widthBound      =  2: 1:10; 
% heightBound     =  2: 1:10; 
% varianceBound   =  0: 2:16;
% nBlocksBound    = 10:10:90;

% nBlocks         = 40;
% height          =  4;
% heightVariance  = 10;
% width           =  4;
% widthVariance   = 10;
% indexShiftx = 0;
% indexShifty = 0;
% 
% coeffFun = @(vertices) coeffFun_randomBlocks(vertices(:,1),vertices(:,2),N,n,nBlocks,width:width+widthVariance,height:height+heightVariance,indexShiftx,indexShifty);

%% Zufalls Koeffizientenfunktion
% randomPercentageBound = 0.2:0.05:0.7;
% randomStateBound = 1:10;

% randomPercentage    = 0.7;
% randomState         = 1;
% indexShiftx         = 0;
% indexShifty         = 0;
% 
% coeffFun = @(vertices) coeffFun_random(vertices(:,1),vertices(:,2),randomPercentage,randomState,indexShiftx,indexShifty);

%% Plot
file_path = "./resources/train_data/2022-07-08-09-40-22-parameter_dump.csv";
opts = detectImportOptions(file_path);
params = readmatrix(file_path,opts);

parameters = cell2mat(cellfun(@ str2num,params(8:end,3),'UniformOutput',false));
parameters = unique(parameters,'rows');
parameters = parameters(:,3:end);

fig = figure();
t = tiledlayout('flow');
rhoMin = 1;
rhoMax = 10^6;
markerType = 'verts';
vertTris = load("./libs/precomputed_vertTris.mat").vertTris;
for i = length(parameters):-1:length(parameters) - 9
    nBlocks         = parameters(i,1);
    height          = parameters(i,2);
    heightVariance  = parameters(i,3);
    width           = parameters(i,4);
    widthVariance   = parameters(i,5);
    indexShiftx     = parameters(i,6);
    indexShifty     = parameters(i,7);
    coeffFun = @(vertices) coeffFun_randomBlocks(vertices(:,1),vertices(:,2),N,n,nBlocks,width:width+widthVariance,height:height+heightVariance,indexShiftx,indexShifty);

    
    plotCoefficientMatrices(coeffFun,markerType,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid,vertTris);
end








function [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = plotCoefficientMatrices(f_coeff,markerType,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot,vertTris)
% Input: xCanalLim,yCanalLim: Grenzen des Kanalgebiets in x- und y-Richtung
% Input: rhoMax,rhoMin: rho im Kanal und außerhalb des Kanals
% Input: vert,tri: Knoten- und Elementliste
% Input: logicalTri__sd: Logischer Vektor, welche Dreiecke in welchem TG enthalten sind
% Input: plot: Boolean, ob das Gitter mit Kanal geplottet werden soll

% Output: rhoTri,rhoTriSD: Koeffizient pro Element (und teilgebietsweise)
% Output: indElementsCanal: Logischer Vektor, welche Elemente im Kanal liegen
% Output: maxRhoVert,maxRhoVertSD: maximaler Koeffizient pro Knoten (und teilgebietsweise)

if exist("vertTris","var")
    assert(length(vertTris) == length(vert),'Unpassendes Gitter fuer das geladene vertTris')
    loadedVertTris = 1;
else
    loadedVertTris = 0;
end
numSD = length(logicalTri__sd);
N = sqrt(numSD);
numTri = length(tri);
numVert = length(vert);

 %% Definiere Koeffizientenfunktion auf den Elementen
if strcmp('verts',markerType)
    markedVertices = find(f_coeff(vert)); % Knotenindizes der markierten Knoten
    
    % Idee: direkt die markedElements zurueckgeben lassen: die elementliste in
    % knoten indizes uebersetzen, reshapen. markierung prüfen -> zurueck
    % reshapen -> mittels any die elemente markieren. In der FKT: bei nargout =
    % 1 die markierten Knoten, bei nargout = 2 die markierten elemente
    % zurueckgeben
    
    % Erstelle Vektor welcher die Koeffizientenfunktion pro Element angibt. Ein
    % Element gilt als markiert wenn alle zugehörigen Eckknoten markiert sind.
    % Alle markierten Elemente werden auf rhoMax und der Rest auf rhoMin gesetzt.
    rhoTri = rhoMin*ones(numTri,1);
    markedElements = (sum(ismember(tri,markedVertices),2)==3);
elseif strcmp('elements',markerType)
    rhoTri = rhoMin*ones(numTri,1);
    markedElements = find(f_coeff(tri));
else
    error('Ungueltigen Markertyp angegeben.')
end
rhoTri(markedElements) = rhoMax;

%% Definiere Koeffizientenfunktion auf den Elementen eines TG
rhoTriSD = cell(numSD,1);
for i = 1:numSD
    rhoTriSD{i} = rhoTri(logicalTri__sd{i});  % Koeffizientenfunktion pro Element teilgebietsweise
end

%% Definiere maximalen Koeffizienten pro Knoten
maxRhoVert = zeros(numVert,1);
if ~loadedVertTris
    vertTris = cell(numVert,1);
end
maxRhoVertSD = cell(numVert,1);
% [vertTris2,test] = cellfun(@(x) find(x==tri),num2cell(1:numVert),'UniformOutput',false)';
for i = 1:numVert % Iteriere ueber Knoten
    if ~loadedVertTris
        [vertTris{i},~,~] = find(i == tri);
    end

    maxRhoVert(i) = max(rhoTri(vertTris{i})); % Maximaler Koeffizient pro Knoten
    
    %% Definiere maximalen Koeffizienten pro Knoten eines TG
    for k = 1:numSD % Iteriere ueber TG
        vertTrisSD = logicalTri__sd{k}(vertTris{i}); % Logischer Vektor, welche Dreiecke des Knotens im TG liegen
        maxRhoVertSD{i} = [maxRhoVertSD{i},max(rhoTri(vertTris{i}(vertTrisSD)))]; % Maximaler Koeffizient pro Knoten teilgebietsweise
    end
end

if plot
    %% Plotten des Gitters mit Kanal
%     figure("Name","Triangulierung des Gebiets mit Koeffizientenfunktion");
    gcf
    nexttile
    patch('vertices',vert,'faces',tri,'facecol',[1,1,1],'edgecolor',"#FFFFFF"); 
    hold on; axis equal tight;
    patch('vertices',vert,'faces',tri(markedElements,:),'facecol',"#2b8cbe",'edgecolor',"#5a5a5a");
    for i = 1:N-1
        line([0,1],[i/N,i/N],'LineWidth', 1.5, 'color', 'r')
        line([i/N,i/N],[0,1],'LineWidth', 1.5, 'color', 'r')
    end
end
end