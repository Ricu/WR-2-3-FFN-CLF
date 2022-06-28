function [markedVertices] = coeffFun_random(x,y,randomPercentage,randomState,indexShiftx,indexShifty)
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

