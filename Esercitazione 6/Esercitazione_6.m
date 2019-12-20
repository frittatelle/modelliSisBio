%% ESERCITAZIONE 6

% Obbiettivo: Costruire a partire dai dati un modello in
% grado di predire i parametri cinetici del C-Peptide usando
% il metodo della regressione lineare.

clear
dati_vc;

%% Data

% misure antropomorfiche (covariate), x
stato_salute = dati(:, 1);
sesso = dati(:, 2);
eta = dati(:, 3);
altezza = dati(:, 4);
peso = dati(:, 5);
BMI = dati(:, 6);
BSA = dati(:, 7);
X = [stato_salute sesso eta altezza peso BMI BSA];
labels_cov = {'stato salute' 'sesso' 'eta' 'altezza' 'peso' 'BMI' 'BSA'};

% parametri, y
volume_distribuzione = dati(:, 8);
t_emivita_corto = dati(:, 9);
t_emivita_lungo = dati(:, 10);
fraction = dati(:, 11);
Y = [volume_distribuzione t_emivita_corto t_emivita_lungo fraction];
labels_par = {'volume distribuzione' 't emivita corto' 't emivita lungo' 'fraction'};

%% Analisi monovariata

% Graficare i dati raccolti in funzione delle singole covariate

r = zeros(length(labels_par),length(labels_cov));

% Per ogni parametro ricavato grafico le relazioni con le diverse covariate
for i = 1:length(labels_par)
    figure(i)
    for j=1:length(labels_cov)
        
        if j<=2
        subplot(3,3,j)
        boxplot( Y(:,i),X(:,j))
        ylabel(labels_par{i})
        xlabel(labels_cov{j})
        grid
        
        else             
        R = corrcoef(Y(:, i), X(:, j));
        r(i, j) = R(1, 2);                
        subplot(3,3,j)
        scatter( X(:, j),Y(:,i),'.')
        hold on 
        h = lsline;
        h.Color = 'r';
        hold off;
        
        ylabel(labels_par{i})
        xlabel(labels_cov{j})
        title(strcat('r = ',num2str(r(i, j))));
        grid                
        end 
    end 
end

%% Regressione multipla 

% Dummy var stato salute (1 vect) --> (3 vect) 
stato_salute_dmy = dummyvar(stato_salute + 1);  

n_soggetti = size(Y,1);
n_parametri = size(Y,2);
X_cand_labels = labels_cov;

% Per ciascun parametro della cinetica trovo il modello di
% regressione opportuno attraverso una stepwise regression.
% Una volta trovati i regressori calcolo gli intervalli di
% confidenza delle stime:
% - Stimo la varianza dell'errore con i residui
% - Calcolo la devizione standard
% - Costruisco gli intervalli di confidenza

for i = 1:n_parametri
      
	X_cand = {stato_salute_dmy sesso eta altezza peso BMI BSA};      % Vettore candidati
    
    % Selezione dei regressori (Forward stepwise regression)
	[X_selected{i}, X_selected_labels{i}, R2_adjusted_max{i}] = stepwise_reg(Y(:,i), X_cand, X_cand_labels);   
    
    
    % Predizione con regressori ottimi
    X = X_selected{i};               % Matrice dei regressori ottimi
    B{i} = pinv(X'*X)*X'*Y(:,i);     % Stima matrice coefficienti di regressione
    Y_hat{i} = X*B{i};               % Predizione del modello
    
    
    % Calcolo cifre di merito
    E{i} = Y(:,i) - Y_hat{i};           % Errore di predizione del modello
    SSE{i} = sum((E{i}).^2);            % Somma dei quadrati dell'errore 
    SST{i} = sum((Y - mean(Y)).^2);     % Somma dei quadrati totale
    SSR{i} = SST{i} - SSE{i};           % Somma dei quadrati di regressione
    
    R2{i} = SSR{i}/SST{i};                     % Proporzione di variabilita' spiegata sul totale
    p{i} = length(X_selected_labels{i});       % numero dei regressori del modello
    s2{i} = SSE{i}/(n_soggetti-p{i}-1);        % varianze dell'errore residuo
    
    sigma_B{i} = pinv(X'*(1/s2{i})*X);                % covarianza delle stime
    CV_B{i} = sqrt(diag(sigma_B{i}))./sigma_B{i};     % CV stime
    std_B{i} = sqrt(diag(sigma_B{i}));                % deviazione standard delle stime
    
    alpha = 0.05;
    
    % CI Beta stimati
    for j = 1:length(B{i})
        CI_B{i} = [B{i} - tcdf(alpha/2,n_soggetti-p{i}-1)*std_B{i}(j)  B{i} + tcdf(alpha/2,n_soggetti-p{i}-1)*std_B{i}(j)];  
    end 
    
    % T-TEST
    % H0: beta k sono significativamente diversi da 0?
    for j = 1:length(B{i})
        t = abs(B{i}(j)/(std_B{i}(j)));
        p_value_B{i}(j) = 1 - tcdf(t,n_soggetti-p{i}-1);
    end
    
    
end

% F-TEST
% Verifico la significativita' statistica dei regressori inclusi nel modello
for i = 1:n_parametri
    F{i}=(SSR{i}/p{i})/(SSE{i}/(n_soggetti - p{i} - 1));
    p_value{i} = 1 - fcdf(F{i},p{i},n_soggetti - p{i} - 1);    
end 


%% My data

stato_salute_L = [1 0 0] ;
sesso_L = 0 ;
eta_L = 22;
altezza_L = 1.80;
peso_L = 72;
BMI_L = peso_L/(altezza_L)^2;
BSA_L = 0.20247*((altezza_L^0.725)*(peso_L^0.425));

dati_L = {stato_salute_L sesso_L eta_L altezza_L peso_L BMI_L BSA_L};

% Volume distribuzione
% Modello: V_dist = B0 + B1*BSA + B2*eta;
b = B{1,1};
volume_distribuzione_L = b(1) + b(2) * BSA_L + b(3) * eta_L;  

% tempo di emivita corto
% Modello: t_emi_short = B0 + B1*BMI + B2*eta;
b = B{1,2};
t_emi_short_L = b(1) + b(2) * BMI_L + b(3) * eta_L; 

% tempo di emivita lungo
% Modello: t_emi_long = B0 + B1*eta + B2*altezza + B3*stato_salute;
b = B{1,3};
t_emi_long_L = b(1) + b(2) * eta_L + b(3) * altezza_L + b(4)*stato_salute_L; 
t_emi_long_L = t_emi_long_L(1);

% fraction
% Modello: fraction = B0 + B1*BMI
b = B{1,4};
fraction_L = b(1) + b(2) * BMI_L;

% Parametri cincetica C-peptide usando i miei dati
Y_L = [volume_distribuzione_L t_emi_short_L t_emi_long_L fraction_L];

% Costruire gli intervalli di confidenza dei parametri della
% cinetica ricavati con i diversi modelli
% Xi * B_hat +- t_n-p-1 * sy_hat
% sy_hat e' la deviazione standard; 
% s_y = X'*(X'*sigma_B^-1*X)^-1 * X

% es. volume di distribuzione
xVolL = [1 BSA_L eta_L];   
s_volDist = xVolL*sigma_B{1}*xVolL';   % stdev volDist
IC_volDist = [volume_distribuzione_L - 1.972*s_volDist volume_distribuzione_L + 1.972*s_volDist];

%% Regressione con interazione regressori

% Rule of thumb: per esserci interazione, i regressori che
% interagiscono devono essere stati selezionati nel modello
% senza interazione

X_cand = {stato_salute_dmy eta altezza BMI BSA eta.*altezza eta.*BSA BMI.*BSA BMI.^2 BSA.^2 eta.^2 altezza.^2 };      % Vettore candidati con interazioni
X_cand_labels = {'stato salute','eta','altezza','BMI','BSA','eta x altezza','eta x BSA','BMI x BSA','BMI^2','BSA^2','eta^2','altezza^2'};

for i = 1:n_parametri
      
	  
    % Selezione dei regressori (Forward stepwise regression)
	[X_selected_interaction{i}, X_selected_interaction_labels{i}, R2_adjusted_max{i}] = stepwise_reg(Y(:,i), X_cand, X_cand_labels);   
    
    
    % Predizione con regressori ottimi
    X = X_selected_interaction{i};               % Matrice dei regressori ottimi
    B_interaction{i} = pinv(X'*X)*X'*Y(:,i);     % Stima matrice coefficienti di regressione
    Y_hat{i} = X*B_interaction{i};               % Predizione del modello
    
    
    % Calcolo cifre di merito
    E{i} = Y(:,i) - Y_hat{i};           % Errore di predizione del modello
    SSE{i} = sum((E{i}).^2);            % Somma dei quadrati dell'errore 
    SST{i} = sum((Y - mean(Y)).^2);     % Somma dei quadrati totale
    SSR{i} = SST{i} - SSE{i};           % Somma dei quadrati di regressione
    
    R2_interaction{i} = SSR{i}/SST{i};                              % Proporzione di variabilita' spiegata sul totale
    p{i} = length(X_selected_interaction_labels{i});    % numero dei regressori del modello
    s2{i} = SSE{i}/(n_soggetti-p{i}-1);                 % varianze dell'errore residuo
    
    sigma_B{i} = pinv(X'*(1/s2{i})*X);                  % covarianza delle stime
    CV_B{i} = sqrt(diag(sigma_B{i}))./sigma_B{i};       % CV stime
    std_B{i} = sqrt(diag(sigma_B{i}));                  % deviazione standard delle stime
    
    alpha = 0.05;
    
    % CI Beta stimati
    for j = 1:length(B_interaction{i})
        CI_B_interaction{i} = [B_interaction{i} - tcdf(alpha/2,n_soggetti-p{i}-1)*std_B{i}(j)  B_interaction{i} + tcdf(alpha/2,n_soggetti-p{i}-1)*std_B{i}(j)];  
    end 
    
    % T-TEST
    % H0: beta k sono significativamente diversi da 0?
    for j = 1:length(B{i})
        t = abs(B_interaction{i}(j)/(std_B{i}(j)));
        p_value_B{i}(j) = 1 - tcdf(t,n_soggetti-p{i}-1);
    end
    
end

% Commento: dopo aver selzionato i regressori 'semplici'
% sono stati inseriti regressori derivati da interazione tra
% i regressori selezionati in precedenza e si e' visto un
% miglioramento delle performance dei modelli.
% Vedi R2 e R2_interaction per dettagli numerici.

%% My data new models

L = CPepModels(stato_salute_L,sesso_L,eta_L,altezza_L,peso_L);
