function [markedVertices] = coeffFun_subdomains(x,y,affectedSubdomains,vert__sd,indexShiftx,indexShifty)
%Input: x   x-Koordinaten aller Knoten
%Input: y   y-Koordinaten aller Knoten
%Input: affectedSubdomains     Teilgebiete mit maximalem Koeffizienten
%Input: vert__sd    Knotenliste teilgebietsweise
%Input: indexShiftx    Verschiebt die Kaneale in x-Richtung
%Input: indexShifty    Verschiebt die Kaneale in y-Richtung

%Output: markedVertices    Gibt an welchen Knoten der maximale Koeffizient zugewiesen wird

assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')
markedVertices = false(size(x)); %initialisiere mit logical false fuer Punkt
for i = 1:length(affectedSubdomains)
    sd = affectedSubdomains(i);
    indx = x >= min(vert__sd{sd}(:,1)) & x <= max(vert__sd{sd}(:,1));
    indy = y >= min(vert__sd{sd}(:,2)) & y <= max(vert__sd{sd}(:,2));
    markedVertices = markedVertices | (indx & indy);
end

markedMatrix = reshape(markedVertices,sqrt(length(y)),sqrt(length(y)));
markedMatrix= circshift(markedMatrix,indexShiftx,2);
markedMatrix= circshift(markedMatrix,indexShifty,1);
markedVertices = markedMatrix(:);
end

