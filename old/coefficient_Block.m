function [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_Block(propB1,propB2,difB,positionC,widthC,numberC,h,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot)
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
SD_size = 1/N;

%Bloecke fuer jedes zweite TG spaltenweise ab 1.Spalte
indx = false(numVert,1);   %initialisiere mit logical false fuer x-Koordinate
for j = 0:2:N-2
    bool1 = j*SD_size;
    bool2 = (propB1+j)*SD_size;
    indx = indx | (bool1 <= vert(:,1)) & (vert(:,1) <= bool2);
end
indy = false(numVert,1);  %initialisiere mit logical false fuer y-Koordinate
for i = 0:N-1
    bool3 = (i+0.5)*SD_size + h*(positionC -0.5*widthC);
    bool4 = (i+0.5)*SD_size + h*(positionC +0.5*widthC);
    indy = indy | (bool3 <= vert(:,2)) & (vert(:,2) <= bool4);
end
indBlock1 = (indx&indy);

%Bloecke fuer jedes zweite TG spaltenweise ab 2.Spalte
indx = false(numVert,1);   %initialisiere mit logical false fuer x-Koordinate
for j = 1:2:N-1
    bool1 = (propB2+j)*SD_size;
    bool2 = (1+j)*SD_size;
    indx = indx | (bool1 <= vert(:,1)) & (vert(:,1) <= bool2);
end
indy = false(numVert,1);  %initialisiere mit logical false fuer y-Koordinate
for i = 0:N-1
    bool3 = (i+0.5)*SD_size + h*(positionC -0.5*widthC +difB);
    bool4 = (i+0.5)*SD_size + h*(positionC +0.5*widthC +difB);
    indy = indy | (bool3 <= vert(:,2)) & (vert(:,2) <= bool4);
end
indBlock2 = (indx&indy);

indBlock = indBlock1|indBlock2;
numVertCanal = find(indBlock); % Knotennummern der Knoten, die maximalen Koeffizienten haben

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