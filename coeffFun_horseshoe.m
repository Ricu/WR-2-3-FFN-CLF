function [markedVertices] = coeffFun_horseshoe(x,y,N,n,yStripeLim,position,width,hight)
assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')
SD_size = 1/N;
h = 1/(N*n);

%Erstellen der vertikalen Streifen
indySt = (yStripeLim(1) <= y) & (y <= yStripeLim(2));
indxSt = false(size(x));    %initialisiere mit logical false fuer x-Koordinate
for j = 0:N-1
    bool1 = SD_size*(j+0.5)+h*position;
    bool2 = SD_size*(j+0.5)+h*(1+position+width);
    indxSt = indxSt | ((bool1 <= x) & (x <= bool2));
end
indSt = indySt & indxSt;  %Logischer Vektor, welche Knoten in den Streifen liegen


%Erstelle Hufeisen
indHuf = false(size(x)); %Logischer Vektor, welche Knoten im Hufeisen liegen
                         %initialisiere mit logical false fuer x- und y-Koordinate
for j = 0:N-1
    for i = 1:N-1
        bool3 = SD_size*j + h*position;
        indx1 = ((bool3 <= x) & (x <= bool3 +h)) | ((bool3 +2*h <= x) & (x <= bool3 +3*h));
        indy1 = (SD_size*i - h*hight <= y) & (y <= SD_size*i + h*hight);
        indHuf = indHuf | (indx1&indy1);
        indx2 = (bool3 +h <= x) & (x <= bool3 +2*h);
        indy2 = (SD_size*i -2*h <= y) & (y <= SD_size*i -h);
        indHuf = indHuf | (indx2&indy2);
    end
end
markedVertices = indSt | indHuf;
end

