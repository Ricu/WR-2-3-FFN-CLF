function [input] = generate_input(edgeID,edgesSD,rhoTriSD,vert__sd,tri__sd)
% Input: edgeID: Kantenindex
% Input: edgesSD:  Kantenliste mit angrenzenden Teilgebietsnummern
% Input: rhoTriSD: Koeffizienten pro Element (teilgebietsweise)
% Input: vert__sd: Knotenliste (teilgebietsweise)
% Input: tri__sd: Elementliste (teilgebietsweise)
%
% Output: Input fuer neuronales Netz: Koeffizientenverteilung an Samplepunkten

h = vert__sd{1}(2,2);   % Schrittweite knotenweise
N = sqrt(length(vert__sd)); % Anzahl TG in eine Koordinatenrichtung

% Parameter welche vor Erstellen des Netzwerkes festgelegt werden muessen:
% Bestimmen die Anzahl Samplepunkte und damit die Inputlayergroeße
coverage = 3/4; % Anteil der TG, der von Samplepunkten abgedeckt werden soll
n = 1/(N*h);    % Anzahl Knoten pro TG in eine Koordinatenrichtung

% Cell-Array, in welchem die Koeffizienten an den Samplepunkten der TG gespeichert werden
input_cell = cell(1,2);
for i = 1:2 % Iteriere ueber die angrenzenden TG der Kante
    % Bestimmte die Grenzen des Teilgebiets
    sd = edgesSD(edgeID,i);
    left_sd_bound = min(vert__sd{sd}(:,1));
    right_sd_bound = max(vert__sd{sd}(:,1));
    bottom_sd_bound = min(vert__sd{sd}(:,2));
    top_sd_bound = max(vert__sd{sd}(:,2));
    
    % Verschiebe Grenze des TG, da wir nicht komplettes TG mit
    % Samplepunkten abdecken wollen, sondern nur den coverage Anteil
    if abs(edgesSD(edgeID,1)-edgesSD(edgeID,2)) == 1 
        % Fall 1: TG liegen nebeneinander
        xsize = (right_sd_bound-left_sd_bound); % Breite des aktuellen TG
        if i == 1
            % Fall 1a: TG liegen nebeneinander und betrachten linkes TG 
            left_sd_bound = right_sd_bound - coverage * xsize;
        else
            % Fall 1b: TG liegen nebeneinander und betrachten rechtes TG
            right_sd_bound = left_sd_bound + coverage * xsize;
        end
        % Waehle Samplepunkte an den Mittelpunkten der Elemente
        x1 = linspace(left_sd_bound+ h/3,right_sd_bound - 2/3*h,coverage*n);
        y1 = linspace(bottom_sd_bound + h/3, top_sd_bound - 2/3*h, n);
        x2 = linspace(left_sd_bound+ 2/3*h,right_sd_bound - 1/3*h,coverage*n);
        y2 = linspace(bottom_sd_bound + 2/3*h, top_sd_bound - 1/3*h, n);
    else
        % Fall 2: TG liegen uebereinander
        ysize = (right_sd_bound-left_sd_bound);
        if i == 1
            % Fall 2a: TG liegen uebereinander und betrachten unteres TG
            bottom_sd_bound = top_sd_bound - coverage * ysize;
        else
            % Fall 2b: TG liegen uebereinander und betrachten oberes TG
            top_sd_bound = bottom_sd_bound + coverage * ysize;
        end
        % Waehle Samplepunkte an den Mittelpunkten der Elemente
        x1 = linspace(left_sd_bound+ h/3,right_sd_bound - 2/3*h,n);
        y1 = linspace(bottom_sd_bound + h/3, top_sd_bound - 2/3*h, coverage*n);
        x2 = linspace(left_sd_bound+ 2/3*h,right_sd_bound - 1/3*h,n);
        y2 = linspace(bottom_sd_bound + 2/3*h, top_sd_bound - 1/3*h, coverage*n);
    end

    % Erstelle das Samplegitter
    [XX1,YY1] = meshgrid(x1,y1);
    [XX2,YY2] = meshgrid(x2,y2);

    if ~(abs(edgesSD(edgeID,1)-edgesSD(edgeID,2)) == 1)
        % Falls die TG uebereinander liegen, transponiere die Reihenfolge
        % der Knoten. Dadurch wird eine konsistente Nummerierung bzgl des
        % Abstands zur Kante gewaehrleistet.
        XX1 = XX1'; YY1 = YY1';
        XX2 = XX2'; YY2 = YY2';
    end

    X = [XX1(:); XX2(:)]; Y = [YY1(:); YY2(:)]; % Gitter in Vektor umwandeln
    % In X_element wird fuer jeden Knoten gespeichert, in welchem Element
    % dieser liegt.
    X_element = zeros(length(X),1);
    % Pruefe fuer jeden Knoten in welchem Element er liegt
    for ind= 1:size(tri__sd{sd},1)
        element = tri__sd{sd}(ind,:);
        element_corners_x = vert__sd{sd}(element,1);
        element_corners_y = vert__sd{sd}(element,2);
        in = inpolygon(X,Y,element_corners_x,element_corners_y);
        X_element(in) = ind;
    end
    % Extrahiere den Koeffizienten fuer jeden Knoten, je nachdem in welchem 
    % Element er liegt.
    input_cell{i} = rhoTriSD{sd}(X_element)';

    % Da wir immer nur rho als 1 oder 10^6 waehlen: normalisiere bereits
    % hier auf 0 und 1, um Speicherplatz zu sparen
    input_cell{i}(input_cell{i} == min(rhoTriSD{sd})) = 0;
    input_cell{i}(input_cell{i} == max(rhoTriSD{sd})) = 1;
    
    % Plotbefehle zum Ueberpruefen der Knotennummerierung
%     scatter(X,Y,'filled')
%     str = compose('%g',1:length(X));
%     text(X,Y,str,'HorizontalAlignment','left')
end
input = cell2mat(input_cell);
end

