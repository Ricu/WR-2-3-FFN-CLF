function [markedElements] = coeffFun_horseshoe(tri,x,y,N,n,yStripeLim,xDisplacement,width,height)
%Input: tri    Elementliste mit Knotennummern
%Input: x   x-Koordinaten aller Knoten
%Input: y   y-Koordinaten aller Knoten
%Input: N   Anzahl Teilgebiete in einer Richtung
%Input: n   Feinheit des Gitters
%Input: yStripeLim    Intervall, indem der Kanal in y-Richtung liegt
%Input: xDisplacement     Verschiebung der Figuren in x-Richtung   
%Input: width   Breite der Figuren
%Input: height    Hoehe der Hufeisen-Figuren

%Output: markedElements    Gibt an welchen Elementen der maximale Koeffizient zugewiesen wird

assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')
SD_size = 1/N;
h = 1/(N*n);
numTri = size(tri,1);
markedElements = false(numTri,1); % Logischer Vektor, welche Elemente hoeheren Koeffizienten haben

tol = 10^(-6); %Toleranz gegen Rundungsfehler
 
%Erstellen der vertikalen Streifen
indySt = (yStripeLim(1) <= y) & (y <= yStripeLim(2));  %Logischer Vektor, welche Knoten im Streifen liegen
indxSt = false(size(x));    %initialisiere mit logical false fuer x-Koordinate
for j = 0:N-1
    bool1 = SD_size*(j+0.5)+h*xDisplacement;
    bool2 = SD_size*(j+0.5)+h*(1+xDisplacement+width);
    indxSt = indxSt | ((bool1 <= x) & (x <= bool2));
end
indSt = indySt & indxSt;
vertNum = find(indSt); % Knotennummern der Knoten mit hoehem Koeffizienten
markedElements(sum(ismember(tri,vertNum),2)==3) = true;

%Erstelle Hufeisen
indHuf = false(size(x)); %Logischer Vektor, welche Knoten im Hufeisen liegen
                         %initialisiere mit logical false fuer x- und y-Koordinate
for j = 0:N-1   %iteriere ueber TG in x-Richtung
    for i = 1:N-1   %iteriere ueber y-Koordinaten
        bool3 = SD_size*j + h*xDisplacement;
        bool4 = SD_size*i;
        %linker Streifen des Hufeisens
        elements1 = false(numTri,1);
        indx1 = (bool3 -tol <= x) & (x <= bool3 + width*h +tol) ;
        indy1 = (bool4 - h*height -tol <= y) & (y <= bool4 + h*height +tol);
        ind1 = indx1&indy1;
        vertNum1 = find(ind1); % Knotennummern der Knoten mit hoehem Koeffizienten
        elements1(sum(ismember(tri,vertNum1),2)==3) = true;

        %rechter Streifen des Hufeisens
        elements2 = false(numTri,1);
        indx2 = (bool3 +2*width*h -tol <= x) & (x <= bool3 +3*width*h +tol);
        indy2 = (bool4 - h*height -tol <= y) & (y <= bool4 + h*height +tol);
        ind2 = indx2&indy2;
        vertNum2 = find(ind2); % Knotennummern der Knoten mit hoehem Koeffizienten
        elements2(sum(ismember(tri,vertNum2),2)==3) = true;

        %Verbiundung zwischen Streifen
        elements3 = false(numTri,1);
        indx3 = (bool3 -tol <= x) & (x <= bool3 +3*width*h +tol);
        indy3 = (bool4 -(width+height)*h -tol <= y) & (y <= bool4 -h*height +tol);
        ind3 = indx3&indy3;
        vertNum3 = find(ind3); % Knotennummern der Knoten mit hoehem Koeffizienten
        elements3(sum(ismember(tri,vertNum3),2)==3) = true;


        markedElements = markedElements | (elements1 | elements2 | elements3);
        indHuf = indHuf | ind1 | ind2 | ind3;
    end
end
% markedVertices = indSt | indHuf;

% %% Plotten des Gitters mit Kanal
% figure(2);
% scatter(vert(markedVertices,1),vert(markedVertices,2),[],"r") % Dirichletknoten markieren

end

