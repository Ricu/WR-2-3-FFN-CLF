function [markedVertices] = coeffFun_canal(y,N,n,yOffset,canal_width,n_canals,indexShifty)
%Input: y   y-Koordinaten aller Knoten
%Input: N   Anzahl Teilgebiete in einer Richtung
%Input: n   Feinheit des Gitters
%Input: canal_width    Breite des Kanals, Wert von 0 entspricht der initialen Breite und ist abhaengig von der Anzahl an Kanaelen
%Input: indexShifty    Verschiebt die Kaneale in y-Richtung

%Output: 
SD_size = 1/N;
h = 1/(N*n);
propStripes = SD_size/(2*n_canals+1); %Gibt an in wie viele Teile das TG vom Kanal geteilt wird
%indx = true(1,numVert);     %Kanal unabh. von x-Koordinate
markedVertices = false(size(y));    %initialisiere mit logical false fuer y-Koordinate
for j = 0:N-1 %iteriere ueber die Teilgebiete in y-Richtung
    for i = 1:n_canals
        a = (2*i-1)*propStripes + j*SD_size + h*(yOffset -0.5*canal_width);
        b = 2*i*propStripes + j*SD_size + h*(yOffset+0.5*canal_width);
        markedVertices = markedVertices | (a <= y) & (y <= b); 
    end
end

markedMatrix = reshape(markedVertices,sqrt(length(y)),sqrt(length(y)));
markedMatrix= circshift(markedMatrix,indexShiftx,2);
markedMatrix= circshift(markedMatrix,indexShifty,1);
markedVertices = markedMatrix(:);
end

