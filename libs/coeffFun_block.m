assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')
SD_size = 1/N;
h = 1/(N*n);

%Bloecke fuer jedes zweite TG spaltenweise ab 1.Spalte
indx = false(numVert,1);   %initialisiere mit logical false fuer x-Koordinate
for j = 0:2:N-2
    bool1 = j*SD_size;
    indx = indx | (bool1 <= x) & (x <= bool2);
end
indy = false(numVert,1);  %initialisiere mit logical false fuer y-Koordinate
for i = 0:N-1
    indy = indy | (bool3 <= y) & (y <= bool4);
end
indBlock1 = (indx&indy);

%Bloecke fuer jedes zweite TG spaltenweise ab 2.Spalte
indx = false(numVert,1);   %initialisiere mit logical false fuer x-Koordinate
for j = 1:2:N-1
    bool2 = (1+j)*SD_size;
    indx = indx | (bool1 <= x) & (x <= bool2);
end
indy = false(numVert,1);  %initialisiere mit logical false fuer y-Koordinate
for i = 0:N-1
    indy = indy | (bool3 <= y) & (y <= bool4);
end
indBlock2 = (indx&indy);

markedVertices = indBlock1|indBlock2;
end

