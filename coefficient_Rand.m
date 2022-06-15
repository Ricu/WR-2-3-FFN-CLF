function [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_Rand(numRand,widthR,heightR,h,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot)
% Input: xCanalLim,yCanalLim: Grenzen des Kanalgebiets in x- und y-Richtung
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
%SD_size = 1/N;
%propStripes = SD_size/(2*numberC+1); %Gibt an in wie viele Teile das TG vom Kanal geteilt wird

%Waehle random Knoten im Gitter und erstelle um diese herum Bloecke
%hoeherer Koeffizientenfunktionen
numCor = sqrt(numVert);
corVec = vert(1:numCor,2);

indy = false(numVert,1);    %initialisiere mit logical false fuer y-Koordinate
indx = false(numVert,1);    %initialisiere mit logical false fuer x-Koordinate
indRand = false(numVert,1); %initialisiere mit logical false fuer Punkt
for r = 1:numRand
    randx = randi([1 numCor]);
    randy = randi([1 numCor]);
    vertx = corVec(randx);
    verty = corVec(randy);

    indx = indx | (vertx - widthR*h <= vert(:,1)) & (vert(:,1) <= vertx + widthR*h);
    indy = indy | (verty - heightR*h <= vert(:,2)) & (vert(:,2) <= verty + heightR*h);

    indRand = indRand | (indx & indy);
end

indVertCanal = indRand;  % Logischer Vektor, welche Knoten im Kanal liegen
numVertCanal = find(indVertCanal); % Knotennummern der Knoten, die im Kanal liegen

rhoTri = rhoMin*ones(numTri,1); % Koeffizient (zuerst) auf allen Elementen = Koeffizient (gleich) auf Elementen außerhalb des Kanals

rhoTri((sum(ismember(tri,numVertCanal),2)==3)) = rhoMax;

% for i=1:numTri % Iteriere ueber die Elemente
%     if (min(abs(tri(i,1)-numVertCanal))==0 && min(abs(tri(i,2)-numVertCanal))==0 && min(abs(tri(i,3)-numVertCanal))==0)   % Alle Knoten des Elements liegen im Kanal und damit das Element selber
%         rhoTri(i)=rhoMax;    % Koeffizient auf Elementen innerhalb des Kanals
%     end
% end
indElementsCanal = rhoTri == rhoMax; % Logischer Vektor, welche Elemente im Kanal liegen

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