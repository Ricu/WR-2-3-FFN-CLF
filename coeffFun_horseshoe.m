function [markedVertices] = coeffFun_horseshoe(vert,x,y,N,n,yStripeLim,position,width,hight)
assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')
SD_size = 1/N;
h = 1/(N*n);

%Erstellen der vertikalen Streifen
indySt = (yStripeLim(1) <= y) & (y <= yStripeLim(2));  %Logischer Vektor, welche Knoten im Streifen liegen
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
for j = 0:N-1   %iteriere ueber TG in x-Richtung
    bool3 = SD_size*j + h*position;
    for i = 1:N-1   %iteriere ueber y-Koordinaten
        bool4 = SD_size*i;
        indx1 = ((bool3 <= x) & (x <= bool3 +h +10^6)) | ((bool3 +2*h <= x) & (x <= bool3 +3*h));
        indy1 = (bool4 - h*hight <= y) & (y <= bool4 + h*hight);
        indHuf = indHuf | (indx1&indy1);
        indx2 = (bool3 <= x) & (x <= bool3 +3*h);
        indy2 = (bool4 -2*h <= y) & (y <= bool4 -h);
        indHuf = indHuf | (indx2&indy2);
    end
end
markedVertices = indSt | indHuf;


%% Plotten des Gitters mit Kanal
figure(2);
scatter(vert(markedVertices,1),vert(markedVertices,2),[],"r") % Dirichletknoten markieren

end

