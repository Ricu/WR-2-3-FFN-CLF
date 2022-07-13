function [edgesPrimalGlobal,cGamma,edgesSD,cLocalPrimal,cB,cBskal,cInner,cK,cDirichlet] = setup_matrices(rho_struct,grid_struct,f)
% Input: rho_struct: Structure mit allen Koeffizientenkomponenten
%        Komponenten: rhoTriSD, maxRhoVert, maxRhoVertSD
% Input: grid_struct: Structure mit allen Gitterkomponenten:
%        Komponenten: vert__sd,tri__sd,l2g__sd,dirichlet
% Input: f: rechte Seite der DGL

% Output: edgesPrimalGlobal: Cell-Array: mit primalen Knoten pro TG-Kante (global)
% Output: cGamma: Cell-Array: Interfaceknoten pro TG
% Output: edgesSD:  Kantenliste mit angrenzenden Teilgebietsnummern
% Output: cLocalPrimal: Cell-Array: mit primalen Knoten pro TG-Kante (global)
% Output: cB: Cell-Array: lokale Sprungoperatoren
% Output: cBskal: Cell-Array: skalierte Sprungoperatoren
% Output: cInner: Cell-Array: Innere Knoten pro TG
% Output: cK: Cell-Array: lokale Steifigkeitsmatrizen
% Output: cDirichlet: Cell-Array: Dirichlet Knoten pro TG

%% Structures entpacken
rhoTriSD = rho_struct.rhoTriSD; % Koeffizienten pro Element teilgebietsweise
maxRhoVertSD = rho_struct.maxRhoVertSD; % Maximaler Koeffizient pro Knoten teilgebietsweise

vert__sd = grid_struct.vert__sd;
tri__sd  = grid_struct.tri__sd;
l2g__sd  = grid_struct.l2g__sd;
dirichlet = grid_struct.dirichlet;

numSD = length(vert__sd);       % Anzahl Teilgebiete
numVert = length(dirichlet);    % Anzahl Knoten Global


%% Partition der Knoten in Interfaceknoten und in primale und duale Knoten
% Zaehle Anzahl Teilgebiete, in denen Knoten enthalten ist
multiplicity = zeros(numVert,1);
for i = 1:numSD
    multiplicity(l2g__sd{i}) = multiplicity(l2g__sd{i}) + 1;
end
gamma = (multiplicity > 1) & ~dirichlet;    % Extrahiere Interfaceknoten
primal = (multiplicity == 4) & ~dirichlet;  % Extrahiere primale Knoten
dual = gamma & ~primal;                     % Extrahiere duale Knoten

%% Partition der Knotengruppen teilgebietsweise
cDirichlet = cell(numSD,1);
cInner = cell(numSD,1);
cGamma = cell(numSD,1);
cDual = cell(numSD,1);
cPrimal = cell(numSD,1);
cIDual = cell(numSD,1);
for i = 1:numSD
    cDirichlet{i} = dirichlet(l2g__sd{i});      % Dirichletknoten pro TG
    cGamma{i} = gamma(l2g__sd{i});              % Interfaceknoten pro TG
    cInner{i} = ~(cGamma{i} | cDirichlet{i});   % Innere Knoten pro TG
    cPrimal{i} = primal(l2g__sd{i});            % Primale Knoten pro TG
    cDual{i} = dual(l2g__sd{i});                % Duale Knoten pro TG
    cIDual{i} = cInner{i} | cDual{i};           % Innere+Duale Knoten pro TG
end

%% Mappings
% Mapping: Global -> Interface
mapGamma = zeros(numVert,1);
mapGamma(gamma)=1:nnz(gamma);

% Mapping: Global -> Primal
mapPrimal = zeros(numVert,1);
mapPrimal(primal)=1:nnz(primal);

% Mapping: Global -> Dual
mapDual = zeros(numVert,1);
mapDual(dual)=1:nnz(dual);

% Mappings: Interface lokal -> Interface global
%           Primal lokal -> Primal global
%           Dual lokal -> Dual global
cGammaMap = cell(numSD,1);
cPrimalMap = cell(numSD,1);
cDualMap = cell(numSD,1);
for i = 1:numSD
    cGammaMap{i} = mapGamma((l2g__sd{i}(cGamma{i})));
    cPrimalMap{i} = mapPrimal((l2g__sd{i}(cPrimal{i})));
    cDualMap{i} = mapDual((l2g__sd{i}(cDual{i})));
end

%% Lagrange Multiplikatoren Info
cLM = cell(sum(dual),1);
for i = 1:numSD
    i_ind = find(cDual{i}); % Lokale Knotennummern der dualen Knoten des TG
    for j = 1:length(i_ind)
        s=cDualMap{i}(j);   % Lagrangescher-Multipikator (duale globale Knotennummer)
        cLM{s}=[cLM{s},[i;i_ind(j)]]; % Enthaelt TG-Nummer und lokale Knotennummer des LM
    end
end

%% Lokale Sprungoperatoren: mit und ohne Skalierung
% Initialisierung
n_LM = sum(multiplicity(dual)-1);   % Anzahl Lagrangescher Multiplikatoren
cB = cell(1,numSD);
cBskal=cell(1,numSD);
for i = 1:numSD
    cB{i}=sparse(n_LM,length(vert__sd{i}));
    cBskal{i}=sparse(n_LM,length(vert__sd{i}));
end

% Sprungoperator ohne Skalierung
row_ind_LM = 1; % Jede Zeile von cB repraesentiert einen LM
for i = 1:length(cLM)   % Iteriere ueber LM
    % Erstes zugehoeriges TG des LM
    cB{cLM{i}(1,1)}(row_ind_LM,cLM{i}(2,1)) = 1;   % cLM{i}(1,.) TG-Nummer, cLM{i}(2,.) lokale Knotennummer
    for j = 2:size(cLM{i},2)    % Iteriere ueber restliche zugehoerige TG des LM
        cB{cLM{i}(1,j)}(row_ind_LM,cLM{i}(2,j)) = -1;
        row_ind_LM = row_ind_LM + 1;
    end
end

% Sprungoperator mit Skalierung
row_ind_LM = 1;
for i = 1:length(cLM) % Iteriere ueber LM
    globNum = l2g__sd{cLM{i}(1,1)}(cLM{i}(2,1));    % Globale Knotennummer des dualen Knotens
    sumMaxRhoDual = sum(maxRhoVertSD{globNum}); % Summer ueber maximale Koeffizienten der zugehoerigen TG zu diesem Knoten
    % Verwende zur Skalierung jeweils den Koeffizienten des anderen TG
    cBskal{cLM{i}(1,1)}(row_ind_LM,cLM{i}(2,1)) = (maxRhoVertSD{globNum}(2)/sumMaxRhoDual)*cB{cLM{i}(1,1)}(row_ind_LM,cLM{i}(2,1));
    cBskal{cLM{i}(1,2)}(row_ind_LM,cLM{i}(2,2)) = (maxRhoVertSD{globNum}(1)/sumMaxRhoDual)*cB{cLM{i}(1,2)}(row_ind_LM,cLM{i}(2,2));
    row_ind_LM = row_ind_LM + 1;
end

%% Assembliere die lokalen Steifigkeitsmatrizen und die lokalen Lastvektoren
cK = cell(numSD,1); % Steifigkeitsmatrizen
cb = cell(numSD,1); % Lastvektoren
store = 0;
if store 
    storedMatrices = cell(numSD,1);
    for i = 1:numSD
        [cK{i},~,cb{i},storedMatrices{i}] = assemble(tri__sd{i}, vert__sd{i},1,f,rhoTriSD{i},store);
    end
    save("./libs/localMatrices_N4n40.mat",'storedMatrices')
else
    storedMatrices = load("./libs/localMatrices_N4n40.mat").storedMatrices;
    for i = 1:numSD
        [cK{i},~,cb{i}] = assemble(tri__sd{i}, vert__sd{i},1,f,rhoTriSD{i},store,storedMatrices{i});
    end
end

%% Assembliere globale Steifigkeitsmatrix in primalen Variablen
K_PiPiTilde = sparse(sum(primal),sum(primal));
f_PiTilde = sparse(sum(primal),1);
for i = 1:numSD
    K_PiPiTilde(cPrimalMap{i},cPrimalMap{i}) = K_PiPiTilde(cPrimalMap{i},cPrimalMap{i}) + cK{i}(cPrimal{i},cPrimal{i});
    f_PiTilde(cPrimalMap{i}) = f_PiTilde(cPrimalMap{i}) + cb{i}(cPrimal{i});
end

%% Nebenbedingung Vorarbeit

%% Erstelle Kantenlisten
cEdgesSD = cell(1,1);
for i = 1:length(cLM)
    cEdgesSD{1} = [cEdgesSD{1};cLM{i}(1,:)];
end
edgesSD = unique(cEdgesSD{1},'rows'); % Enthaelt die beiden angrenzenden Teilgebietsnummern pro TG-Kante
numEdges = size(edgesSD,1);           % Anzahl der Kanten

edgesDualGlobalAll = cell(size(edgesSD)); % Enthaelt fuer jedes angrenzende TG die dualen Knoten (GLOBALE Knotennummern)
edgesDual = cell(numEdges,1);             % Enthaelt fuer die TG-Kanten die dualen Knoten (DUALE Knotennummern)
edgesDualGlobal = cell(numEdges,1);       % Enthaelt fuer die TG-Kanten die dualen Knoten (GLOBALE Knotennummern)

edgesPrimalGlobalAll = cell(size(edgesSD)); % Enthaelt fuer jedes angrenzende TG die primalen Knoten (GLOBALE Knotennummern)
edgesPrimalGlobal = cell(numEdges,1);       % Enthaelt fuer die TG-Kanten die primalen Knoten (GLOBALE Knotennummern)

for i = 1:numEdges % Iteriere ueber TG-Kanten
    for j = 1:size(edgesSD,2) % Iteriere ueber angrenzende TG
        SD = edgesSD(i,j);
        edgesDualGlobalAll{i,j} = l2g__sd{SD}(cDual{SD});      % Enthaelt fuer jedes angrenzende TG die dualen Knoten (GLOBALE Knotennummern)
        edgesPrimalGlobalAll{i,j} =  l2g__sd{SD}(cPrimal{SD}); % Enthaelt fuer jedes angrenzende TG die primalen Knoten (GLOBALE Knotennummern)
    end
    edgesDualGlobal{i} = intersect(edgesDualGlobalAll{i,1},edgesDualGlobalAll{i,2}); % Enthaelt fuer die TG-Kanten die dualen Knoten (GLOBALE Knotennummern)
    edgesDual{i} = mapDual(edgesDualGlobal{i});     % Enthaelt fuer die TG-Kanten die dualen Knoten (DUALE Knotennummern)
    edgesPrimalGlobal{i} = intersect(edgesPrimalGlobalAll{i,1},edgesPrimalGlobalAll{i,2}); % Enthaelt fuer die TG-Kanten die primalen Knoten (GLOBALE Knotennummern)
end

% Umkehrabbildung von globaler zu lokaler Nummerierung
g2l__sd = cell(numSD,1);
for sd = 1:numSD
    g2l__sd{sd} = zeros(numVert,1);
    ind = l2g__sd{sd};
    g2l__sd{sd}(ind) = 1:length(l2g__sd{sd});
end

cLocalPrimal = cell(size(edgesSD)); % Enthaelt die auf einer TG-Kante liegenden primalen Knoten des TG in lokaler Nummerierung
for i = 1:numEdges
    for j = 1:2
        sd = edgesSD(i,j);
        cLocalPrimal{i,j} = g2l__sd{sd}(edgesPrimalGlobal{i}); % Enthaelt die auf einer TG-Kante liegenden primalen Knoten des TG in lokaler Nummerierung
    end
end
end

