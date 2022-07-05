function plot_classified_edges(f_coeff,markerType,rhoMax,rhoMin,vert,tri,predicted_labels,true_labels,vert__sd)

N = sqrt(length(vert__sd));
numSD = N^2;

edgesSD = [(1:numSD-N)',(N+1:numSD)';...
           setdiff(1:numSD,N:N:numSD)', setdiff(1:numSD,1:N:numSD)'];
edgesSD = sortrows(edgesSD);
nEdges = length(edgesSD);

% edge_info
edge_ends = cell(length(edgesSD),6);
for edgeID = 1:nEdges
    common_vert = intersect(vert__sd{edgesSD(edgeID,1)},vert__sd{edgesSD(edgeID,2)},'rows');
    xmin = min(common_vert(:,1)); xmax = max(common_vert(:,1));
    ymin = min(common_vert(:,2)); ymax = max(common_vert(:,2));
    xLine = [xmin,xmax];
    yLine = [ymin,ymax];
    edge_ends{edgeID,1} = xLine;
    edge_ends{edgeID,2} = yLine;
end

TP = predicted_labels == 1 & true_labels == 1;
TN = predicted_labels == 0 & true_labels == 0;
FP = predicted_labels == 1 & true_labels == 0;
FN = predicted_labels == 0 & true_labels == 1;


%% Definiere Koeffizientenfunktion auf den Elementen
if strcmp('verts',markerType)
    markedVertices = find(f_coeff(vert)); % Knotenindizes der markierten Knoten
    markedElements = (sum(ismember(tri,markedVertices),2)==3);
elseif strcmp('elements',markerType)
    markedElements = find(f_coeff(tri));
else
    error('Ungueltigen Markertyp angegeben.')
end

%% Plotten des Gitters mit Kanal
figure("Name","Triangulierung des Gebiets mit Koeffizientenfunktion");
p1 = patch('vertices',vert,'faces',tri,'facecol',[1,1,1],'edgecolor',"#5a5a5a");
hold on; axis equal tight;
p2 = patch('vertices',vert,'faces',tri(markedElements,:),'facecol',"#2b8cbe",'edgecolor',"#5a5a5a");

title("Triangulierung mit Koeffizientenfunktion")
temp = (1:2:2*N)/(2*N);
yt = reshape(repmat(temp,N,1),N^2,1);
xt = repmat(temp,1,N);
str = compose('%g',1:N^2);
text(xt,yt,str,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)

%% Plotten der klassifizierten Kanten
h1 =  line([0,0],[0,0],'LineWidth', 1.5, 'color', 'k');
h2 =  line([0,0],[0,0],'LineWidth', 1.5, 'color', 'r');
h3 =  line([0,0],[0,0],'LineWidth', 1.5, 'color', 'y');

line(cell2mat(edge_ends(TP | TN,1))',cell2mat(edge_ends(TP | TN,2))','LineWidth', 1.5, 'color', 'k');
line(cell2mat(edge_ends(FN,1))',cell2mat(edge_ends(FN,2))','LineWidth', 1.5, 'color', 'r');
line(cell2mat(edge_ends(FP,1))',cell2mat(edge_ends(FP,2))','LineWidth', 1.5, 'color', 'y');

rhoMax = sprintf('\\rho = %.0e',rhoMax);
rhoMin = sprintf('\\rho = %g',rhoMin);
legend([p1; p2; h1; h2; h3],rhoMin,rhoMax,'TN/TP','FN','FP')
end

