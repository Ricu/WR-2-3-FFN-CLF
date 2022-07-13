function [markedVertices] = coeffFun_randomBlocks(x,y,N,n,n_blocks,widthBound,heightBound,indexShiftx,indexShifty)
%Input: x   x-Koordinaten aller Knoten
%Input: y   y-Koordinaten aller Knoten
%Input: N   Anzahl Teilgebiete in einer Richtung
%Input: n   Feinheit des Gitters
%Input: n_blocks   Anzahl an random erstellten Bloecken
%Input: widthbound    Array mit moeglichen Breiten der Bloecke
%Input: hightbound    Array mit moeglichen Hoehen der Bloecke
%Input: indexShiftx    Verschiebt die Kaneale in x-Richtung
%Input: indexShifty    Verschiebt die Kaneale in y-Richtung

%Output: markedVertices    Gibt an welchen Knoten der maximale Koeffizient zugewiesen wird

assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')

numCoord = N*n+1; %Anzahl Koordinaten in einer Dimension
coordVec = linspace(0,1,numCoord)';
h = 1/(N*n);

markedVertices = false(size(x));
for current_block = 1:n_blocks
    randx = randi([1 numCoord]);
    randy = randi([1 numCoord]);
    randwidth = randi([1 length(widthBound)]);
    randheight = randi([1 length(heightBound)]);

    vertx = coordVec(randx);
    verty = coordVec(randy);
    block_width = widthBound(randwidth);
    block_height = heightBound(randheight);

    indx = (vertx - h*0.5*block_width <= x) & (x <= vertx + h*0.5*block_width);
    indy = (verty - h*0.5*block_height <= y) & (y <= verty + h*0.5*block_height);

    markedVertices = markedVertices | (indx & indy);
end

markedMatrix = reshape(markedVertices,sqrt(length(y)),sqrt(length(y)));
markedMatrix= circshift(markedMatrix,indexShiftx,2);
markedMatrix= circshift(markedMatrix,indexShifty,1);
markedVertices = markedMatrix(:);
end

