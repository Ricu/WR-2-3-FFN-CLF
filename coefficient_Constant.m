function [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_Constant(numSub,vert,tri,logicalTri__sd,plot)
% Input: rhoMax,rhoMin: rho Werte
% Input: affectedSubdomains: markiert in welchen TG rhoMax angenommen wird
% Input: vert,tri: Knoten- und Elementliste
% Input: logicalTri__sd: Logischer Vektor, welche Dreiecke in welchem TG enthalten sind
% Input: plot: Boolean, ob das Gitter mit Kanal geplottet werden soll

% Output: rhoTri,rhoTriSD: Koeffizient pro Element (und teilgebietsweise)
% Output: indElementsrhoMax: Logischer Vektor, welche Elemente in rhoMax liegen
% Output: maxRhoVert,maxRhoVertSD: maximaler Koeffizient pro Knoten (und teilgebietsweise)

numSD = length(logicalTri__sd);
N = sqrt(numSD);
numTri = length(tri);
numVert = length(vert);

%% Definiere Koeffizientenfunktion auf den Elementen
%SD_size = 1/N;
%propStripes = SD_size/(2*numberC+1); %Gibt an in wie viele Teile das TG vom Kanal geteilt wird

rhoTri = ones(numTri,1); % Die Koeffizientenfunktion entspricht rhoMin außerhalb der markierten TG

rho = 10.^(0:8);
for s = 1:numSub
    randSub = randi([1 N^2]);
    randRho = randi([1 length(rho)]);
    RhoSub = rho(randRho);

    rhoTri(logicalTri__sd{randSub}) = RhoSub;
end

%indElementsCanal = rhoTri == rhoMax; % Logischer Vektor, welche Elemente im Kanal liegen

%% Definiere Koeffizientenfunktion auf den Elementen eines TG
rhoTriSD = cell(numSD,1);
for i = 1:numSD
    rhoTriSD{i} = rhoTri(logicalTri__sd{i}); % Koeffizientenfunktion pro Element teilgebietsweise
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
%     hold on; axis equal tight;
%     patch('vertices',vert,'faces',tri(indElementsrhoMax,:),'facecol',"#2b8cbe",'edgecolor',"#5a5a5a");
    for i = 1:N-1
        line([0,1],[i/N,i/N],'LineWidth', 1.5, 'color', 'r')
        line([i/N,i/N],[0,1],'LineWidth', 1.5, 'color', 'r')
    end
    %rhoMax = sprintf('\\rho = %.0e',rhoMax);
    %legend('\rho = 1',rhoMax,'Interface','','','')
    title("Triangulierung mit Koeffizientenfunktion")
end

end

