function [input] = generate_input_2(edgeID,edgesSD,rhoTriSD,vert__sd,tri__sd)

sd1 = edgesSD(edgeID,1);
sd2 = edgesSD(edgeID,2);
n_vertx = 40;
n_verty = 40;
input_cell = cell(2,1);
coverage = 3/4;


for i = 1:2
    sd = edgesSD(edgeID,i);
    left_sd_bound = min(vert__sd{sd}(:,1));
    right_sd_bound = max(vert__sd{sd}(:,1));
    bottom_sd_bound = min(vert__sd{sd}(:,2));
    top_sd_bound = max(vert__sd{sd}(:,2));


    if abs(sd1-sd2) == 1
        xsize = (right_sd_bound-left_sd_bound);
        if i == 1
            left_sd_bound = right_sd_bound - coverage * xsize;
        else
            right_sd_bound = left_sd_bound + coverage * xsize;
        end
    else
        ysize = (right_sd_bound-left_sd_bound);
        if i == 1
            bottom_sd_bound = top_sd_bound - coverage * ysize;
        else
            top_sd_bound = bottom_sd_bound + coverage * ysize;
        end
    end

    x = linspace(left_sd_bound,right_sd_bound,n_vertx+2);
    y = linspace(bottom_sd_bound,top_sd_bound,n_verty+2);

    [XX,YY] = meshgrid(x(2:end-1),y(2:end-1));

    if ~(abs(sd1-sd2) == 1)
        XX = XX'; YY = YY';
    end

    X = XX(:); Y = YY(:);
    X_element = zeros(length(X),1);

    for ind= 1:size(tri__sd{sd},1)
        element = tri__sd{sd}(ind,:);
        element_corners_x = vert__sd{sd}(element,1);
        element_corners_y = vert__sd{sd}(element,2);
        in = inpolygon(X,Y,element_corners_x,element_corners_y);
        X_element(in) = ind;
    end
    input_cell{i} = rhoTriSD{sd}(X_element);
    scatter(X,Y)
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

