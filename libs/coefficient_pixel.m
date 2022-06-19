function [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = coefficient_pixel(rhoMax,rhoMin,pic_bw,num_pixel,vert,tri,logicalTri__sd,plot)
% Input: rhoMax,rhoMin: rho Werte
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
rhoTri = rhoMin*ones(numTri,1); % Die Koeffizientenfunktion entspricht rhoMin außerhalb der markierten TG

for p = 1 : length(tri)
    %Berechne Schwerpunkt fuer rechtwinkliges Dreieckselement
    x_centroid = (vert(tri(p,1),1)+vert(tri(p,2),1)+vert(tri(p,3),1))/3;
    y_centroid = (vert(tri(p,1),2)+vert(tri(p,2),2)+vert(tri(p,3),2))/3;
    %Bestimme die Pixelnummer in x- und y-Richtung, indem wir auf den
    %neahsten integer-Wert runden
    err_cor = 0.5*(1/num_pixel); %korrigiere zum Pixelmittelpunkt
    x_pixel = round((x_centroid-err_cor)*num_pixel); 
    y_pixel = round((y_centroid-err_cor)*num_pixel);
    
    if x_pixel < 1
        x_pixel = x_pixel +1;
    else if x_pixel > num_pixel
        x_pixel = x_pixel -1;
    end
    end

    if y_pixel < 1
        y_pixel = y_pixel +1;
    else if y_pixel > num_pixel
        y_pixel = y_pixel -1;
    end
    end
   
    %werte die Pixelnummer im Bild aus
    value_bw = pic_bw(x_pixel,y_pixel);
    if value_bw < 100
        rhoTri(p) = rhoMax;
    end
end

indElementsrhoMax = (rhoTri == rhoMax); % Logischer Vektor, welche Elemente in rhoMax liegen

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
    hold on; axis equal tight;
    patch('vertices',vert,'faces',tri(indElementsrhoMax,:),'facecol',"#2b8cbe",'edgecolor',"#5a5a5a");
    for i = 1:N-1
        line([0,1],[i/N,i/N],'LineWidth', 1.5, 'color', 'r')
        line([i/N,i/N],[0,1],'LineWidth', 1.5, 'color', 'r')
    end
    rhoMax = sprintf('\\rho = %.0e',rhoMax);
    legend('\rho = 1',rhoMax,'Interface','','','')
    title("Triangulierung mit Koeffizientenfunktion")
end

end

