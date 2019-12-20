function [X_selected, X_selected_labels, R2_adj_max] = stepwise_reg(Y, X_cand, X_cand_labels)
%STEPWISE REGRESSION
% Y = BX + E
% Y = misure effettuate(parametri cinetica)

n_covariate = size(X_cand,2);
n = size(Y,1);

%% Modello 0 regressori

X = ones(n,1);                
B = pinv(X'*X)*X'*Y;          % Stima matrice B dei coefficienti di regressione
Y_hat = X*B;                  % Predizione del modello
E = Y - Y_hat ;               % Residui
SSE = sum((E).^2) ;           % Somma dei quadrati dell'errore (devianza non spiegata)
SST = sum((Y - mean(Y)).^2);  % Somma dei quadrati totale (dev spiegata + n_spiegata)
SSR = SST - SSE;              % Somma dei quadrati spiegata  
R2_std = SSR/SST;             % Quantita' di variabilita' spiegata sul totale

% Man mano che aggiungo regressori R2 si avvicina ad 1
% Il modello migliore e' quello che ha R2 piu' alto e basso n di regressori


%% Forward stepwise

% Regressione multipla semplice (max 7 regressori)
%   Ciclo sui candidati
%   Per ogni candidato calcolo una cifra di riferimento( R2 o R2_adj)
%   Cerco il regressore che ha aumentato l'R2 complessivo del modello (se esiste)
%   Aggiorno R2 max e aggiungo il regressore nel vettore degli X
%   selezionati

X_selected = [];            % Vettore regressori selezionati
X_selected_labels = {};     % Labels regeressori selezionati
% z = 1;                      % Inizializzazione indice labels regressori selezionati


% Il ciclo sarebbe piu' corretto con un while - uso il for per semplicita'
for i = 1 : n_covariate
    
    % Calcolo delle cifre di merito 
    for j = 1 : n_covariate
        X = [X X_cand{j}];          % aggiungo temporaneamente regressore_i 
        B = pinv(X'*X)*X'*Y;        % stimo B con il nuovo regressore
        Y_hat = X*B;                % Predizione del modello
        E = Y - Y_hat ;             % Residui
        SSE = sum((E).^2) ;           
        SST = sum((Y - mean(Y)).^2);  
        SSR = SST - SSE; 
        
        % Se ho una dummy variable k non posso calcolarlo 
        % come il numero di colonne della matrice dei regressori
        % un regressore dummy e' infatti codificato da piu' di una colonna
        if j == 1
            k = size(X,2) - 1 - (size(X_cand{j},2)-1);
        else
            k = size(X,2) - 1;          % numero di regressori utilizzati    
        end
                              
        % Cifre di merito
        R2 = SSR/SST;
        R2_adj(i,j) = R2 - (1-R2)*(k/(n-k-1)); 
        
        % Rimozione regressore temporaneo
        if j == 1
            X(:,end-2:end) = [];
        else
            X(:,end) = [];    
        end       
    end
    
    
    % Check aggiunta regressore al modello
    % Se l'R2 max trovato aggiungendo l'i-esimo regressore
    % e' maggiore di quello ottenuto con un regressore in meno lo aggiungo
    % (Aggiungo il regressore se ha portato un miglioramento)
    if max(R2_adj(i,:)) > R2_std
        
        [R2_std, idx_selected] = max(R2_adj(i,:));     % Aggiornamento R2 max
        
        X = [X X_cand{idx_selected}];                  % Aggiornamento vettore regressori
        
        
        X_selected = X;                                        % Aggiornamento vettore regressori selezionati
        X_selected_labels{i} = X_cand_labels{idx_selected};    % Aggiornamento lista regressori selezionati
        
        R2_adj_max = R2_std;       % R2 adjusted massimo prodotto dall'aggiunta di regressori 
        
    else
        % Se non ci sono stati miglioramenti esco dal ciclo
        break
    end
    
    
end




end

