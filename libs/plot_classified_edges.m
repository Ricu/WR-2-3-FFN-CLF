function plot_classified_edges(predicted_labels,true_labels,vert,tri,N,rhoMax,rhoMin)
%% Ermittlung der (TP, TN, FP, FN) Kanten
numLabels = length(predicted_labels);
confusion_cell = cell(numLabels,1);
for i = 1:numLabels
    if predicted_labels(i) == 1 % als kritisch klassifizierte Kante
        if true_labels(i) == 1  % kritische Kante
            confusion_cell(i) = "TP";
        else % unkritische Kante
           confusion_cell(i) = "FP"; 
        end
    else % als unkritisch klassifizierte Kante
        if true_labels(i) == 1  % kritische Kante
            confusion_cell(i) = "FN";
        else % unkritische Kante
           confusion_cell(i) = "TP"; 
        end
    end
end

%% Plotten des Gitters mit klassifizierten Kanten
figure("Name","Gebiet mit klassifizierten Kanten");
patch('vertices',vert,'faces',tri,'facecol',[1,1,1],'edgecolor',"#5a5a5a");
hold on; axis equal tight;
patch('vertices',vert,'faces',tri(markedElements,:),'facecol',"#2b8cbe",'edgecolor',"#5a5a5a");
for i = 1:N-1
    line([0,1],[i/N,i/N],'LineWidth', 1.5, 'color', 'r')
    line([i/N,i/N],[0,1],'LineWidth', 1.5, 'color', 'r')
end
rhoMax = sprintf('\\rho = %.0e',rhoMax);
rhoMin = sprintf('\\rho = %g',rhoMin);
legend(rhoMin,rhoMax,'Interface','','','')
title("Triangulierung mit Koeffizientenfunktion")
temp = (1:2:2*N)/(2*N);
yt = reshape(repmat(temp,N,1),N^2,1);
xt = repmat(temp,1,N);
str = compose('%g',1:N^2);
text(xt,yt,str,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)

end

