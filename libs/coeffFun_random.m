function [markedVertices] = coeffFun_random(x,y,randomPercentage,randomState,indexShiftx,indexShifty)
%Input: x   x-Koordinaten aller Knoten
%Input: y   y-Koordinaten aller Knoten
%Input: randomPercentage    Auswahlkriterium der random values, entspricht 
        % dem ungefaehren Anteil der Knoten, die mit rhoMax markiert werden  
%Input: randomState    random seed des random number generators
%Input: indexShiftx    Verschiebt die Kaneale in x-Richtung
%Input: indexShifty    Verschiebt die Kaneale in y-Richtung

%Output: markedVertices    Gibt an welchen Knoten der maximale Koeffizient zugewiesen wird

assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')

if nargin > 3
    rng(randomState)
end
markedVertices = rand(size(x)) < randomPercentage;

markedMatrix = reshape(markedVertices,sqrt(length(y)),sqrt(length(y)));
markedMatrix= circshift(markedMatrix,indexShiftx,2);
markedMatrix= circshift(markedMatrix,indexShifty,1);
markedVertices = markedMatrix(:);
end

