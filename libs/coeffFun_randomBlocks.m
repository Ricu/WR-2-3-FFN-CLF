function [markedVertices] = coeffFun_randomBlocks(x,y,N,n,n_blocks)
assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')

numCoord = N*n+1;
coordVec = linspace(0,1,numCoord)';
h = 1/(N*n);

markedVertices = false(size(x)); %initialisiere mit logical false fuer Punkt
width = 1:5;
height = 1:5;
for current_block = 1:n_blocks
    randx = randi([1 numCoord]);
    randy = randi([1 numCoord]);
    randwidth = randi([1 length(width)]);
    randheight = randi([1 length(height)]);

    vertx = coordVec(randx);
    verty = coordVec(randy);
    block_width = width(randwidth);
    block_height = height(randheight);

    indx = (vertx - block_width  * h <= x) & (x <= vertx + block_width  * h);
    indy = (verty - block_height * h <= y) & (y <= verty + block_height * h);

    markedVertices = markedVertices | (indx & indy);
end
end

