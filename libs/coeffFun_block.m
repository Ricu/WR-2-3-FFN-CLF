function [logicalVerticesInBlocks] = coeffFun_block(x,y,N,n,propB1,propB2,difB,position,canal_width,n_canals)
assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')
SD_size = 1/N;
h = 1/(N*n);
propStripes = SD_size/(2*n_canals+1); %Gibt an in wie viele Teile das TG vom Kanal geteilt wird

invindx = false(size(x));  %initialisiere mit logical false fuer inverse x-Koordinate
for j = 0:2:N-2
    left_bound = (propB1+j)*SD_size;
    right_bound = (propB2+j+1)*SD_size;
    invindx = invindx | (left_bound <= x) & (x <= right_bound);
end
indx = not(invindx);

indy = false(size(y));    %initialisiere mit logical false fuer y-Koordinate
for j = 0:N-1
    for i = 1:n_canals
        top_bound = (2*i-1)*propStripes + j*SD_size + h*(position-mod(j,2)*difB);
        bottom_bound = 2*i*propStripes + j*SD_size + h*(position+canal_width-mod(j,2)*difB);
        indy = indy | (top_bound <= y) & (y <= bottom_bound); 
    end
end
logicalVerticesInBlocks = indx & indy;
end

