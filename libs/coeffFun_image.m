function [markedVertices] = coeffFun_image(x,y)
assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')

markedElements = false(size(x));
for p = 1 : length(tri)
    %Berechne Schwerpunkt fuer rechtwinkliges Dreieckselement
    x_centroid = (vert(tri(p,1),1)+vert(tri(p,2),1)+vert(tri(p,3),1))/3;
    y_centroid = (vert(tri(p,1),2)+vert(tri(p,2),2)+vert(tri(p,3),2))/3;
    %Bestimme die Pixelnummer in x- und y-Richtung, indem wir auf den
    %neahsten integer-Wert runden
    err_cor = 0.5*(1/num_pixel); %korrigiere zum Pixelmittelpunkt
    x_pixel = round((x_centroid-err_cor)*num_pixel); 
    y_pixel = round((y_centroid-err_cor)*num_pixel);
    
    if x_pixel < 1
        x_pixel = x_pixel +1;
    else if x_pixel > num_pixel
        x_pixel = x_pixel -1;
    end
    end

    if y_pixel < 1
        y_pixel = y_pixel +1;
    else if y_pixel > num_pixel
        y_pixel = y_pixel -1;
    end
    end
   
    %werte die Pixelnummer im Bild aus
    value_bw = pic_bw(x_pixel,y_pixel);
    if value_bw < 100
        markedElements(p) = 1;
    end
end
end

