function [markedVertices] = coeffFun_block(x,y,N,n,prop1,prop2,dif,width,indexShiftx,indexShifty)
%Input: x   x-Koordinaten aller Knoten
%Input: y   y-Koordinaten aller Knoten
%Input: N   Anzahl Teilgebiete in einer Richtung
%Input: n   Feinheit des Gitters

%Input: canal_width    Breite des Kanals, Wert von 0 entspricht der initialen Breite und ist abhaengig von der Anzahl an Kanaelen
%Input: indexShifty    Verschiebt die Kaneale in y-Richtung

%Output: markedVertices    Gibt an welchen Knoten der maximale Koeffizient zugewiesen wird
assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')
SD_size = 1/N;
h = 1/(N*n);
numVert = length(x);

%Bloecke fuer jedes zweite TG spaltenweise ab 1.Spalte
indx = false(numVert,1);   %initialisiere mit logical false fuer x-Koordinate
for j = 0:2:N-2
    bool1 = j*SD_size;
    bool2 = (prop1+j)*SD_size;
    indx = indx | (bool1 <= x) & (x <= bool2);
end
indy = false(numVert,1);  %initialisiere mit logical false fuer y-Koordinate
for i = 0:N-1
    bool3 = (i+0.5)*SD_size - h * 0.5 * width;
    bool4 = (i+0.5)*SD_size + h * 0.5 * width;
    indy = indy | (bool3 <= y) & (y <= bool4);
end
indBlock1 = (indx&indy);

%Bloecke fuer jedes zweite TG spaltenweise ab 2.Spalte
indx = false(numVert,1);   %initialisiere mit logical false fuer x-Koordinate
for j = 1:2:N-1
    bool1 = (prop2+j)*SD_size;
    bool2 = (1+j)*SD_size;
    indx = indx | (bool1 <= x) & (x <= bool2);
end
indy = false(numVert,1);  %initialisiere mit logical false fuer y-Koordinate
for i = 0:N-1
    bool3 = (i+0.5)*SD_size + h*(dif - 0.5 * width);
    bool4 = (i+0.5)*SD_size + h*(dif + 0.5 * width);
    indy = indy | (bool3 <= y) & (y <= bool4);
end
indBlock2 = (indx&indy);

markedVertices = indBlock1|indBlock2;

markedMatrix = reshape(markedVertices,sqrt(length(y)),sqrt(length(y)));
markedMatrix= circshift(markedMatrix,indexShiftx,2);
markedMatrix= circshift(markedMatrix,indexShifty,1);
markedVertices = markedMatrix(:);
end

