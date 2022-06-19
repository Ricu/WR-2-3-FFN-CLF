clear; clc;
addpath('libs')
plot_grid = true;   % Auswahl: Plotten der Triangulierung mit Bild-Koeffizientenfunktion

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
rhoMin = 1;
rhoMax = 10^6;
%Bilddatei einlesen 
%ii)
pic_ii = imread('./resources/img/rho_coeff_multiple_stripes.png');
pic_ii_bw = pic_ii(:,:,1); %benoetigen nur einen Kanal, da schwarz-weiss Bild
num_pixel_ii = length(pic_ii_bw); %Anzahl Pixel je Dimension

%Erstelle Koeffizientenfunktion fuer Gitter ii)
[rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_pixel(rhoMax,rhoMin,pic_ii_bw,num_pixel_ii,vert,tri,logicalTri__sd,plot_grid);


%iii)
pic_iii = imread('./resources/img/multiple_circle_bw_512x512px.jpeg');
pic_iii_bw = pic_iii(:,:,1); %benoetigen nur einen Kanal, da schwarz-weiss Bild
num_pixel_iii = length(pic_iii_bw); %Anzahl Pixel je Dimension

%Erstelle Koeffizientenfunktion fuer Gitter iii)
[rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_pixel(rhoMax,rhoMin,pic_iii_bw,num_pixel_iii,vert,tri,logicalTri__sd,plot_grid);



