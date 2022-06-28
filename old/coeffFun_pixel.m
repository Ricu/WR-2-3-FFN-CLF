function [markedElements] = coeffFun_pixel(x,y,tri,rhoMax,rhoMin,pic_bw,num_pixel)
assert(all(size(x) == size(y)),'Die Vektoren x und y haben unterschiedliche Groesse')

numTri = length(tri);
rhoTri = rhoMin*ones(numTri,1); %initialisiere die Koeffizientenfunktion je Element mit rhoMin

for p = 1 : numTri
    %Berechne Schwerpunkt fuer rechtwinkliges Dreieckselement
    x_centroid = (x(tri(p,1))+x(tri(p,2))+x(tri(p,3)))/3;
    y_centroid = (y(tri(p,1))+y(tri(p,2))+y(tri(p,3)))/3;
    %Bestimme die Pixelnummer in x- und y-Richtung, indem wir auf den
    %neahsten integer-Wert runden
    err_cor = 0.5*(1/num_pixel); %korrigiere zum Pixelmittelpunkt
    x_pixel = round((x_centroid-err_cor)*num_pixel); 
    y_pixel = num_pixel -round((y_centroid-err_cor)*num_pixel);
    
    if x_pixel < 1
        x_pixel = x_pixel +1;
    else 
        if x_pixel > num_pixel
            x_pixel = x_pixel -1;
        end
    end

    if y_pixel < 1
        y_pixel = y_pixel +1;
    else 
        if y_pixel > num_pixel
            y_pixel = y_pixel -1;
        end
    end
   
    %werte die Pixelnummer im Bild aus
    value_bw = pic_bw(x_pixel,y_pixel);
    if value_bw < 100
        rhoTri(p) = rhoMax;
    else
        rhoTri(p) = rhoMin;
    end
end
markedElements = (rhoTri == rhoMax); % Logischer Vektor, welche Elemente rhoMax-Wert haben
end

