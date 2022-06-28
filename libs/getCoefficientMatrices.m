function [rhoTri,rhoTriSD,maxRhoVert,maxRhoVertSD] = getCoefficientMatrices(f_coeff,markerType,rhoMax,rhoMin,vert,tri,logicalTri__sd,plot,vertTris)
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
    figure("Name","Triangulierung des Gebiets mit Koeffizientenfunktion");
    patch('vertices',vert,'faces',tri,'facecol',[1,1,1],'edgecolor',"#5a5a5a"); 
    hold on; axis equal tight;
    patch('vertices',vert,'faces',tri(markedElements,:),'facecol',"#2b8cbe",'edgecolor',"#5a5a5a");
    for i = 1:N-1
        line([0,1],[i/N,i/N],'LineWidth', 1.5, 'color', 'r')
        line([i/N,i/N],[0,1],'LineWidth', 1.5, 'color', 'r')
    end
    rhoMax = sprintf('\\rho = %.0e',rhoMax);
    rhoMin = sprintf('\\rho = %g',rhoMin);
    legend(rhoMin,rhoMax,'Interface','','','')
    title("Triangulierung mit Koeffizientenfunktion")
    temp = (1:2:2*N)/(2*N);
    yt = reshape(repmat(temp,N,1),N^2,1);
    xt = repmat(temp,1,N);
    str = compose('%g',1:N^2);
    text(xt,yt,str,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
end
end

