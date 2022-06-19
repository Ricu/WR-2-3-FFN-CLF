function [input] = generate_input(edgeID,edgesSD,maxRhoVert,l2g__sd)

sd1 = edgesSD(edgeID,1);
sd2 = edgesSD(edgeID,2);
nVert1 = length(l2g__sd{sd1});
nVert2 = length(l2g__sd{sd2});
if abs(sd1-sd2) == 1
%     % Test Vektoren
%     input1 = reshape((1:nVert1)',sqrt(nVert1),sqrt(nVert1))';
%     input2 = reshape((1:nVert2)',sqrt(nVert2),sqrt(nVert2))';
    input1 = reshape(maxRhoVert(l2g__sd{sd1}),sqrt(nVert1),sqrt(nVert1))';
    input2 = reshape(maxRhoVert(l2g__sd{sd2}),sqrt(nVert2),sqrt(nVert2))';
else
%     % Test Vektoren
%     input1 = reshape((1:nVert1)',sqrt(nVert1),sqrt(nVert1));
%     input2 = reshape((1:nVert2)',sqrt(nVert2),sqrt(nVert2));
    input1 = reshape(maxRhoVert(l2g__sd{sd1}),sqrt(nVert1),sqrt(nVert1));
    input2 = reshape(maxRhoVert(l2g__sd{sd2}),sqrt(nVert2),sqrt(nVert2));
%     input1 = maxRhoVert(l2g__sd{sd1});
%     input2 = maxRhoVert(l2g__sd{sd2});
end

% Die Indizes 111 bis 121 sollen die Kantenknoten von TG 1 enthalten
% Die Indizes 122 bis 132 sollen die Kantenknoten von TG 2 enthalten
% Die Sortierung ist wie folgt:
% Die ersten 11 Knoten von TG 1 sind jene welche am weitesten entfernt von
% der Kante sind. Danach nähert man sich reihenweise an. Erreicht man das
% 2. Teilgebiet so entfernt man sich schrittweise von der Kante.

input = [input1(:) ; input2(:)]';
end

