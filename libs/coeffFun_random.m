function [markedVertices] = coeffFun_random(x,y,randomPercentage,randomState)
assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')

if nargin > 3
    rng(randomState)
end
markedVertices = rand(size(x)) < randomPercentage;

end

