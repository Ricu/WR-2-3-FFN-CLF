function [markedVertices] = coeffFun_randomBlocks(x,y,N,n,n_blocks,widthBound,heightBound,indexShiftx,indexShifty)
%Input: x   x-Koordinaten aller Knoten
%Input: y   y-Koordinaten aller Knoten
%Input: N   Anzahl Teilgebiete in einer Richtung
%Input: n   Feinheit des Gitters
%Input: n_blocks   Anzahl an random erstellten Bloecken
%Input: prop1    Anteil des Blocks in x-Richtung im zugehoerigen TG (jedes zweite ab 1.Spalte)
%Input: prop2    Anteil des Blocks in x-Richtung im zugehoerigen TG (jedes zweite ab 2.Spalte)
%Input: dif    Versetzung der Bloecke spaltenweise in den TG
%Input: hight    Hoehe der Bloecke
%Input: indexShiftx    Verschiebt die Kaneale in x-Richtung
%Input: indexShifty    Verschiebt die Kaneale in y-Richtung

%Output: markedVertices    Gibt an welchen Knoten der maximale Koeffizient zugewiesen wird

assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')

numCoord = N*n+1;
coordVec = linspace(0,1,numCoord)';
h = 1/(N*n);
rng(0);

markedVertices = false(size(x)); %initialisiere mit logical false fuer Punkt
for current_block = 1:n_blocks
    randx = randi([1 numCoord]);
    randy = randi([1 numCoord]);
    randwidth = randi([1 length(widthBound)]);
    randheight = randi([1 length(heightBound)]);

    vertx = coordVec(randx);
    verty = coordVec(randy);
    block_width = widthBound(randwidth);
    block_height = heightBound(randheight);

    indx = (vertx - block_width  * h <= x) & (x <= vertx + block_width  * h);
    indy = (verty - block_height * h <= y) & (y <= verty + block_height * h);

    markedVertices = markedVertices | (indx & indy);
end

markedMatrix = reshape(markedVertices,sqrt(length(y)),sqrt(length(y)));
markedMatrix= circshift(markedMatrix,indexShiftx,2);
markedMatrix= circshift(markedMatrix,indexShifty,1);
markedVertices = markedMatrix(:);
end

