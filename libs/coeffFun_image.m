function [markedVertices] = coeffFun_image(x,y,pic_bw,num_pixel)
%Input: x   x-Koordinaten aller Knoten
%Input: y   y-Koordinaten aller Knoten
%Input: pic_bw     Bildmatrix mit Pixeleintaegen
%Input: num_pixel    Anzahl Pixel je Dimension

%Output: markedVertices    Gibt an welchen Knoten der maximale Koeffizient zugewiesen wird

assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')

numVert = length(x);
value_bw = zeros(numVert,1);

err_cor = 0.5*(1/num_pixel); %korrigiere zum Pixelmittelpunkt
for v = 1 : numVert
    %Bestimme die Pixelnummer in x- und y-Richtung, indem wir auf den
    %neahsten integer-Wert runden
    x_pixel = round((x(v)-err_cor)*num_pixel); 
    y_pixel = num_pixel -round((y(v)-err_cor)*num_pixel);
    
    if x_pixel < 1
        x_pixel = 1;
    else 
        if x_pixel > num_pixel
            x_pixel = num_pixel;
        end
    end

    if y_pixel < 1
        y_pixel = 1;
    else 
        if y_pixel > num_pixel
            y_pixel = num_pixel;
        end
    end
   
    %werte die Pixelnummer im Bild aus
    value_bw(v) = pic_bw(x_pixel,y_pixel);
end

markedVertices = (value_bw < 100); %enthaelt die Knoten mit rhoMax als Koeffizient
end
