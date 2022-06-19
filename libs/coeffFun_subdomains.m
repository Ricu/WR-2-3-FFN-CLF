function [markedVertices] = coeffFun_subdomains(x,y,affectedSubdomains,vert__sd)
assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')
markedVertices = false(size(x)); %initialisiere mit logical false fuer Punkt
for i = 1:length(affectedSubdomains)
    sd = affectedSubdomains(i);
    indx = x >= min(vert__sd{sd}(:,1)) & x <= max(vert__sd{sd}(:,1));
    indy = y >= min(vert__sd{sd}(:,2)) & y <= max(vert__sd{sd}(:,2));
    markedVertices = markedVertices | (indx & indy);
end
end

