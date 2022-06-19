function [markedVertices] = coeffFun_randomBlocks(x,y,n_blocks,h)
assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')
numCor = sqrt(numVert);
corVec = vert(1:numCor,2);

markedVertices = false(size(x),1); %initialisiere mit logical false fuer Punkt
width = 1:5;
height = 1:5;
for current_block = 1:n_blocks
    randx = randi([1 numCor]);
    randy = randi([1 numCor]);
    randwidth = randi([1 length(width)]);
    randheight = randi([1 length(height)]);

    vertx = corVec(randx);
    verty = corVec(randy);
    block_width = width(randwidth);
    block_height = height(randheight);

    indx = (vertx - block_width  * h <= x) & (x <= vertx + block_width  * h);
    indy = (verty - block_height * h <= y) & (y <= verty + block_height * h);

    markedVertices = markedVertices | (indx & indy);
end
end

