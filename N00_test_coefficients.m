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
%
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
% widthBound      =  2:1:10; 
% heightBound     =  2:1:10; 
% varianceBound   =  0:1: 5;
% nBlocksBound    = 10:10:70;

% nBlocks         = 70;
% height          = 2;
% heightVariance  = 0;
% width           = 2;
% widthVariance   = 0;
% indexShiftx = 0;
% indexShifty = 0;
% 
% coeffFun = @(vertices) coeffFun_randomBlocks(vertices(:,1),vertices(:,2),N,n,nBlocks,width:width+widthVariance,height:height+heightVariance,indexShiftx,indexShifty);

%% Zufalls Koeffizientenfunktion
% randomPercentageBound = 0.2:0.05:0.7;
% randomStateBound = 1:10;

randomPercentage    = 0.7;
randomState         = 1;
indexShiftx         = 0;
indexShifty         = 0;

coeffFun = @(vertices) coeffFun_random(vertices(:,1),vertices(:,2),randomPercentage,randomState,indexShiftx,indexShifty);

%% Plot

rhoMin = 1;
rhoMax = 1;
markerType = 'verts';
getCoefficientMatrices(coeffFun,markerType,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot_grid);