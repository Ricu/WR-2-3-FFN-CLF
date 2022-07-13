function [markedVertices] = coeffFun_stripRandomBlocks(x,y,N,n,height,n_strip,n_blocks,widthBound,heightBound,indexShifty)
%Input: x   x-Koordinaten aller Knoten
%Input: y   y-Koordinaten aller Knoten
%Input: N   Anzahl Teilgebiete in einer Richtung
%Input: n   Feinheit des Gitters
%Input: height  Hoehe der Streifen, Wert von 0 entspricht der initialen Breite und ist abhaengig von der Anzahl an Streifen 
%Input: n_strip    Anzahl an Streifen je Teilgebiet
%Input: n_blocks   Anzahl an random erstellten Bloecken
%Input: widthbound    Array mit moeglichen Breiten der Bloecke
%Input: hightbound    Array mit moeglichen Hoehen der Bloecke
%Input: indexShifty    Verschiebt die Kaneale in y-Richtung

%Output: markedVertices    Gibt an welchen Knoten der maximale Koeffizient zugewiesen wird

SD_size = 1/N;
h = 1/(N*n);

%% Markiere die Knoten innerhalb der Streifen
propStripes = SD_size/(2*n_strip+1); %Gibt an in wie viele Teile das TG vom Kanal geteilt wird
%indx = true(1,numVert);     %Kanal unabh. von x-Koordinate
markedVertices = false(size(y));    %initialisiere mit logical false fuer y-Koordinate
for j = 0:N-1 %iteriere ueber die Teilgebiete in y-Richtung
    for i = 1:n_strip
        a = (2*i-1)*propStripes + j*SD_size - h * 0.5 * height;
        b = 2*i*propStripes + j*SD_size + h* 0.5 * height;
        markedVertices = markedVertices | (a <= y) & (y <= b); 
    end
end

markedVerticesStreifen = markedVertices;

%% Markiere die Knoten innerhalb der randomBlocks

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

    indx = (vertx - h*0.5*block_width <= x) & (x <= vertx + h*0.5*block_width);
    indy = (verty - h*0.5*block_height <= y) & (y <= verty + h*0.5*block_height);

    markedVertices = markedVertices | (indx & indy);
end

markedVerticesBlocks = markedVertices;

markedVertices = markedVerticesStreifen | markedVerticesBlocks;

markedMatrix = reshape(markedVertices,sqrt(length(y)),sqrt(length(y)));
markedMatrix= circshift(markedMatrix,indexShifty,1);
markedVertices = markedMatrix(:);
end
