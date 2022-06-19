function [logicalVerticesInCanal] = coeffFun_canal(y,N,n,position,canal_width,n_canals)
SD_size = 1/N;
h = 1/(N*n);
propStripes = SD_size/(2*n_canals+1); %Gibt an in wie viele Teile das TG vom Kanal geteilt wird
%indx = true(1,numVert);     %Kanal unabh. von x-Koordinate
logicalVerticesInCanal = false(size(y));    %initialisiere mit logical false fuer y-Koordinate
for j = 0:N-1
    for i = 1:n_canals
        a = (2*i-1)*propStripes + j*SD_size + h*position;
        b = 2*i*propStripes + j*SD_size + h*(position+canal_width);
        logicalVerticesInCanal = logicalVerticesInCanal | (a <= y) & (y <= b); 
    end
end
end

