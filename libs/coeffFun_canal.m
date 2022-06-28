function [markedVertices] = coeffFun_canal(y,N,n,canal_width,n_canals,indexShiftx,indexShifty)
SD_size = 1/N;
h = 1/(N*n);
propStripes = SD_size/(2*n_canals+1); %Gibt an in wie viele Teile das TG vom Kanal geteilt wird
%indx = true(1,numVert);     %Kanal unabh. von x-Koordinate
markedVertices = false(size(y));    %initialisiere mit logical false fuer y-Koordinate
for j = 0:N-1
    for i = 1:n_canals
        a = (2*i-1)*propStripes + j*SD_size + h;
        b = 2*i*propStripes + j*SD_size + h*canal_width;
        markedVertices = markedVertices | (a <= y) & (y <= b); 
    end
end

markedMatrix = reshape(markedVertices,sqrt(length(y)),sqrt(length(y)));
markedMatrix= circshift(markedMatrix,indexShiftx,2);
markedMatrix= circshift(markedMatrix,indexShifty,1);
markedVertices = markedMatrix(:);
end

