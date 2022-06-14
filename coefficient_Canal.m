function [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_Canal(positionC,widthC,numberC,h,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot)
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
propStripes = SD_size/(2*numberC+1); %Gibt an in wie viele Teile das TG vom Kanal geteilt wird
%indx = true(1,numVert);     %Kanal unabh. von x-Koordinate
indy = false(numVert,1);    %initialisiere mit logical false fuer y-Koordinate
for j = 0:N-1
    for i = 1:numberC
        a = (2*i-1)*propStripes + j*SD_size + h*positionC;
        b = 2*i*propStripes + j*SD_size + h*(positionC+widthC);
        indy = indy | (a <= vert(:,2)) & (vert(:,2) <= b); 
    end
end
indVertCanal = indy;  % Logischer Vektor, welche Knoten im Kanal liegen
numVertCanal = find(indVertCanal); % Knotennummern der Knoten, die im Kanal liegen

rhoTri = rhoMin*ones(numTri,1); % Koeffizient (zuerst) auf allen Elementen = Koeffizient (gleich) auf Elementen außerhalb des Kanals

rhoTri((sum(ismember(tri,numVertCanal),2)==3)) = rhoMax;

% for i=1:numTri % Iteriere ueber die Elemente
%     if (min(abs(tri(i,1)-numVertCanal))==0 && min(abs(tri(i,2)-numVertCanal))==0 && min(abs(tri(i,3)-numVertCanal))==0)   % Alle Knoten des Elements liegen im Kanal und damit das Element selber
%         rhoTri(i)=rhoMax;    % Koeffizient auf Elementen innerhalb des Kanals
%     end
% end
indElementsCanal = rhoTri > 1; % Logischer Vektor, welche Elemente im Kanal liegen

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

