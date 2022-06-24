function [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_Horseshow(yStripeLim,position,width,hight,h,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot)
% Input: yStripeLim: Grenzen der Streifen in y-Richtung
% Input: rhoMax,rhoMin: rho im Kanal und außerhalb des Kanals
% Input: vert,tri: Knoten- und Elementliste
% Input: logicalTri__sd: Logischer Vektor, welche Dreiecke in welchem TG enthalten sind
% Input: plot: Boolean, ob das Gitter mit Kanal geplottet werden soll

% Output: rhoTri,rhoTriSD: Koeffizient pro Element (und teilgebietsweise)
% Output: indElementsCanal: Logischer Vektor, welche Elemente im Kanal liegen
% Output: maxRhoVert,maxRhoVertSD: maximaler Koeffizient pro Knoten (und teilgebietsweise)

numSD = length(logicalTri__sd);
N = sqrt(numSD);
numTri = length(tri);
numVert = length(vert);


%% Definiere Koeffizientenfunktion auf den Elementen
SD_size = 1/N;

%Erstellen der vertikalen Streifen
indySt = (yStripeLim(1) <= vert(:,1) & vert(:,1) <= yStripeLim(2));
indxSt = false(numVert,1);    %initialisiere mit logical false fuer x-Koordinate
for j = 0:N-1
    bool1 = SD_size(j+0.5)+h*position;
    bool2 = SD_size(j+0.5)+h*(1+position+width);
    indxSt = indxSt | (bool1 <= vert(:,1)) & (vert(:,1) <= bool2);
end
indSt = indySt & indxSt;  % Logischer Vektor, welche Knoten in den Streifen liegen


%Erstelle Hufeisen
indHuf = false(numVert,1);    %Logischer Vektor, welche Knoten im Hufeisen liegen
                              %initialisiere mit logical false fuer x- und y-Koordinate
for j = 0:N-1
    for i = 1:N-1
        bool3 = SD_size*j + h*position;
        indx1 = ((bool3 <= vert(:,1)) & (vert(:,1) <= bool3 +h)) | ((bool3 +2*h <= vert(:,1)) & (vert(:,1) <= bool3 +3*h));
        indy1 = (SD_size*i - h*hight <= vert(:,2)) & (vert(:,2) <= SD_size*i + h*hight);
        indHuf = indHuf | (indx1&indy1);
        indx2 = (bool3 +h <= vert(:,1)) & (vert(:,1) <= bool3 +2*h);
        indy2 = (SD_size*i -2*h <= vert(:,2)) & (vert(:,2) <= SD_size*i -h);
        indHuf = indHuf | (indx2&indy2);
    end
end

numVertMax = find(indSt&indHuf); % Knotennummern der Knoten mit maximalem Koeffizeinten

rhoTri = rhoMin*ones(numTri,1); % Koeffizient (zuerst) auf allen Elementen = Koeffizient (gleich) auf Elementen außerhalb des Kanals

rhoTri((sum(ismember(tri,numVertMax),2)==3)) = rhoMax;

indElementsCanal = rhoTri > 1; % Logischer Vektor, welche Elemente hoeheren Koeffizeinten haben

%% Definiere Koeffizientenfunktion auf den Elementen eines TG
rhoTriSD = cell(numSD,1);
for i = 1:numSD
    rhoTriSD{i} = rhoTri(logicalTri__sd{i});  % Koeffizientenfunktion pro Element teilgebietsweise
end

%% Definiere maximalen Koeffizienten pro Knoten
maxRhoVert = zeros(numVert,1);
vertTris = cell(numVert,1); 
maxRhoVertSD = cell(numVert,1);
for i = 1:numVert % Iteriere ueber Knoten
    [vertTris{i},~,~] = find(i == tri);
    maxRhoVert(i) = max(rhoTri(vertTris{i})); % Maximaler Koeffizient pro Knoten
    
    %% Definiere maximalen Koeffizienten pro Knoten eines TG
    for k = 1:numSD % Iteriere ueber TG
        vertTrisSD = logicalTri__sd{k}(vertTris{i}); % Logischer Vektor, welche Dreiecke des Knotens im TG liegen
        maxRhoVertSD{i} = [maxRhoVertSD{i},max(rhoTri(vertTris{i}(vertTrisSD)))]; % Maximaler Koeffizient pro Knoten teilgebietsweise
    end
end

if plot
    %% Plotten des Gitters mit Kanal
    figure("Name","Triangulierung des Gebiets mit Koeffizientenfunktion");
    patch('vertices',vert,'faces',tri,'facecol',[1,1,1],'edgecolor',"#5a5a5a"); 
    hold on; axis equal tight;
    patch('vertices',vert,'faces',tri(indElementsCanal,:),'facecol',"#2b8cbe",'edgecolor',"#5a5a5a");
    for i = 1:N-1
        line([0,1],[i/N,i/N],'LineWidth', 1.5, 'color', 'r')
        line([i/N,i/N],[0,1],'LineWidth', 1.5, 'color', 'r')
    end
    rhoMax = sprintf('\\rho = %.0e',rhoMax);
    rhoMin = sprintf('\\rho = %g',rhoMin);
    legend(rhoMin,rhoMax,'Interface','','','')
    title("Triangulierung mit Koeffizientenfunktion")
end
end
