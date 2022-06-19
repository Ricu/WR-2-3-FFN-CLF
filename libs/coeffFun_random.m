function [markedVertices] = coeffFun_random(x,y,random_percentage,random_state)
assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')

if nargin > 3
    rng(random_state)
end
markedVertices = rand(size(x)) < random_percentage;

end

