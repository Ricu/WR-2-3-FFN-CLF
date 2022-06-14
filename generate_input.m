function [input] = generate_input(edgeID,edgesSD,maxRhoVert,l2g__sd)
sd1 = edgesSD(edgeID,1);
input1 = maxRhoVert(l2g__sd{sd1});

sd2 = edgesSD(edgeID,2);
input2 = maxRhoVert(l2g__sd{sd2});

input = [input1 ; input2];
end

