function label = generate_label(edgeID,edgesPrimalGlobal,cGamma,edgesSD,cLocalPrimal,cB,cBskal,cInner,cK,adaptiveTOL)
% Input: edgeID: Kantenindex
% Input: edgesPrimalGlobal: Cell-Array: mit primalen Knoten pro TG-Kante (global)
% Input: cGamma: Cell-Array: Interfaceknoten pro TG
% Input: edgesSD:  Kantenliste mit angrenzenden Teilgebietsnummern
% Input: cLocalPrimal: Cell-Array: mit primalen Knoten pro TG-Kante (global)
% Input: cB: Cell-Array: lokale Sprungoperatoren
% Input: cBskal: Cell-Array: skalierte Sprungoperatoren
% Input: cInner: Cell-Array: Innere Knoten pro TG
% Input: cK: Cell-Array: lokale Steifigkeitsmatrizen
% Input: adaptiveTOL: Toleranz zur Auswahl der Eigenwerte
%
% Output: label: 0 = unkritische, 1 = kritische Kante 

%% Assemblierungmatrix R_ij
nPrimal = length(edgesPrimalGlobal{edgeID}); % Anzahl primaler Knoten auf der Kante
nGamma = [nnz(cGamma{edgesSD(edgeID,1)}),nnz(cGamma{edgesSD(edgeID,2)})]; % Anzahl Interfaceknoten pro TG
nGammaUnass = sum(nGamma); % Anzahl unassemblierte Interface Knoten beider TG
nRest = nGamma-nPrimal; % Anzahl unassemblierter dualer Knoten beider TG
% Teile die Matrix R in drei Teile auf:
% Knotensortierung in Spalten ist:
% primal, rest(1), rest(2)
% Knotensortierung in Zeilen ist:
% gamma(1), gamma(2) (=primal(1),rest(1),primal(2),rest(2))
P_e = zeros(nGammaUnass,nPrimal); % primal
R_1 = zeros(nGammaUnass,nRest(1)); % rest(1)
R_2 = zeros(nGammaUnass,nRest(2)); % rest(2)

P_e(1:nPrimal,1:nPrimal) = eye(nPrimal); %primal(1),primal
R_1(nPrimal+1 : nGamma(1),:) = eye(nRest(1)); %rest(1),rest(1)

P_e(nGamma(1) + 1:nGamma(1)+nPrimal,1:nPrimal) = eye(nPrimal); %primal(2),primal
R_2(nGamma(1)+nPrimal+1 : nGammaUnass,:) = eye(nRest(2)); % rest(2),rest(2)
R = [P_e, R_1, R_2]; % Zusammensetzen
pi = R*((R'*R)\R');
% pi hat Dimension nGammaUnass x nGammaUnass

%% Sprungoperator und Projektion B und P_D
% P_D wird Dimension dim(B_D_e,2) x dim(B_e,2) haben
% Also muss dim(B_D_e,2) = nGammaUnass entsprechen
% Also muss dim(B_e,2) = nGammaUnass entsprechen
B_D = cell(1,2);
B   = cell(1,2);
n_LM = size(cB{1},1);
for k = 1:2
    sd = edgesSD(edgeID,k);
    B{k} = zeros(n_LM,nGamma(k));
    B_D{k} = zeros(n_LM,nGamma(k));

    % Der folgende Index listet alle lokalen Knotenindizes vom
    % aktuellen Teilgebiet, welche zum Interface gehoeren, aber kein
    % primaler Knoten auf der betrachteten Kante ist. Die primalen
    % Knoten auf der Kante sind lediglich Nullspalten, welche vorne
    % angefuegt werden.
    relevantGamma = setdiff(find(cGamma{sd}),cLocalPrimal{edgeID,k});
    B{k}(:,nPrimal+1:end) = cB{sd}(:,relevantGamma);
    B_D{k}(:,nPrimal+1:end) = cBskal{sd}(:,relevantGamma);
end
B = cell2mat(B);
B_D = cell2mat(B_D);

% Loesche nun alle Zeilen welche nicht zu einem LM auf der
% betrachteten Kante gehoeren.
subset = (full(sum(abs(B),2)) == 2);
B = B(subset,:);
B_D = B_D(subset,:);
P_D = B_D'*B;

%% Schurkomplement S
S_temp = cell(2,1);
% Erstelle nacheinander die lokalen Schurkomplemente
for k = 1:2
    sd = edgesSD(edgeID,k);
    inner_local = cInner{sd};
    % Finde die lokalen Indizes die zu primal und rest gehoeren
    relevantPrimal = cLocalPrimal{edgeID,k};
    relevantGamma = setdiff(find(cGamma{sd}),relevantPrimal);

    % Schreibe die Kombination aus primal und rest in einen 2x2-Cell
    S = cell(2,2);

    % Nur in primalen (primal,primal)
    index_1 = relevantPrimal;
    index_2 = relevantPrimal;

    S{1,1} = cK{sd}(index_1,index_2);
    temp = cK{sd}(inner_local,inner_local) \ cK{sd}(inner_local,index_2);
    S{1,1} = S{1,1} - cK{sd}(index_1,inner_local) * temp;

    % In primalen und restlichen (primal,rest)
    index_1 = relevantPrimal;
    index_2 = relevantGamma;

    S{1,2} = cK{sd}(index_1,index_2);
    temp = cK{sd}(inner_local,inner_local) \ cK{sd}(inner_local,index_2);
    S{1,2} = S{1,2} - cK{sd}(index_1,inner_local) * temp;
    S{2,1} = S{1,2}';

    % Nur in restlichen (rest,rest)
    index_1 = relevantGamma;
    index_2 = relevantGamma;

    S{2,2} = cK{sd}(index_1,index_2);
    temp = cK{sd}(inner_local,inner_local) \ cK{sd}(inner_local,index_2);
    S{2,2} = S{2,2} - cK{sd}(index_1,inner_local) * temp;

    S = cell2mat(S);
    S_temp{k} = S;
end
S = blkdiag(S_temp{:});
% Berechne sigma
sigma = max(diag(S));


%% Stelle den Vector c auf
c = ones(nGammaUnass,1) / norm(ones(nGammaUnass,1));

%% Verallgmeinertes Eigenwertproblem loesen
[eigenvalues] = adaptiveEigenvalues(c, pi, P_D, S,sigma);

% Eigenwerte entsprechend Toleranz auswählen
if nnz(diag(eigenvalues) > adaptiveTOL) > 0
    label = 1;
else
    label = 0;
end

end

%% Hilfsfunktion
function [eigenvalues, eigenvectors] = adaptiveEigenvalues(c, pi, P_D, S,sigma)
% A = Pi*P_D^T*S*P_D*Pi (Korrektur notwendig?)
% Fall 1: c ist NICHT im Kern von S:
% Setze B = B_tilde = Pi * S * Pi + sigma * (I-Pi)
% Fall 2: c ist im Kern von S:
% Setze B = Pibar * B_tilde * Pibar + sigma * (I-Pibar)

% LHS = pi_bar * pi * P_D' * S * P_D * pi * pi_bar;
A =  pi * P_D' * S * P_D * pi;
Btilde = pi * S * pi + sigma*(eye(size(pi))-pi);
if norm(S*c) < 10^(-8)
    pibar = eye(length(c)) - c*c';
    B = pibar * Btilde * pibar + sigma*(eye(size(pibar))-pibar);
else
    B = Btilde;
end
% Korrigiere um evenuelle Asymmetrien
A = 1/2*(A + A');
B = 1/2*(B + B');

% Loese das Eigenwertproblem
[eigenvectors,eigenvalues] = eig(A,B);
end