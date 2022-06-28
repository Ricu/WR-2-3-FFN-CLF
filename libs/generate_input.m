function [input] = generate_input(edgeID,edgesSD,rhoTriSD,vert__sd,tri__sd)
%
%
%

h = vert__sd{1}(2,2);
N = sqrt(length(vert__sd));
n = 1/(N*h);

% Parameter welche vor erstellen des Netzwerkes festgelegt werden muessen!
n_vertx = 40;
n_verty = 40;
coverage = 3/4;

% Cell-Array in welchem die rhoWerte an den Knoten fuer die TG gespeichert
% werden
input_cell = cell(1,2);

% Iteriere ueber die angrenzenden TG
for i = 1:2
    % Bestimmte die Grenzen des Teilgebiets
    sd = edgesSD(edgeID,i);
    left_sd_bound = min(vert__sd{sd}(:,1));
    right_sd_bound = max(vert__sd{sd}(:,1));
    bottom_sd_bound = min(vert__sd{sd}(:,2));
    top_sd_bound = max(vert__sd{sd}(:,2));
    

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
        x1 = linspace(left_sd_bound+ h/3,right_sd_bound - 2/3*h,n);
        y1 = linspace(bottom_sd_bound + h/3, top_sd_bound - 2/3*h, coverage*n);
        x2 = linspace(left_sd_bound+ 2/3*h,right_sd_bound - 1/3*h,n);
        y2 = linspace(bottom_sd_bound + 2/3*h, top_sd_bound - 1/3*h, coverage*n);
    end

    % Gleichmaessig verteilte Punkte erstellen. Erstelle 2 zusaetzliche
    % Punkte welche spaeter rausgeworfen werden, sodass alle Punkte im
    % inneren liegen
%     x = linspace(left_sd_bound,right_sd_bound,n_vertx+2);
%     y = linspace(bottom_sd_bound,top_sd_bound,n_verty+2);

    % Erstelle das Gitter
%     [XX,YY] = meshgrid(x(2:end-1),y(2:end-1));
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
    % In X_element wird fuer jeden Knoten gespeichert in welchem Element
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
    % Extrahiere den rho Wert fuer jeden Knoten entsprechend fuer in
    % welchem Element er liegt.
    input_cell{i} = rhoTriSD{sd}(X_element)';
    
    % Plotbefehle zum ueberpruefen
    scatter(X,Y,'filled')
%     str = compose('%g',1:length(X));
%     text(X,Y,str,'HorizontalAlignment','left')
end

% Die Indizes 111 bis 121 sollen die Kantenknoten von TG 1 enthalten
% Die Indizes 122 bis 132 sollen die Kantenknoten von TG 2 enthalten
% Die Sortierung ist wie folgt:
% Die ersten 11 Knoten von TG 1 sind jene welche am weitesten entfernt von
% der Kante sind. Danach nähert man sich reihenweise an. Erreicht man das
% 2. Teilgebiet so entfernt man sich schrittweise von der Kante.
input = cell2mat(input_cell);
end

