%% ESERCITAZIONE 7 - DECONVOLUZIONE

% Obbiettivo: stimare la secrezione insulinica con diversi
% metodi di deconvoluzione (raw deconvolution,
% regolarizzazione, approccio bayesiano)

clear;

datiEs7;
t = DATA(:,1);       % tempo[min]
c = DATA(:,2);       % concentrazione[pmol/l]
c = c - yb;          % concentrazione corrected[pmol/l]

stato_salute = [0 0 1];
altezza = 1.64;
peso = 59.6;
eta = 32;
sesso = 1;

% Parametri cinetica C-peptide soggetto
s = CPepModels(stato_salute,sesso,eta,altezza,peso);

% Risposta all'impulso
g = s.A * exp(-s.alfa*t) + s.B * exp(-s.beta*t);

figure(1)
subplot(2,1,1)
simplePlot(t,c,'Concentrazione misurata','tempo[min]','concentrazione[pmol/l]');
subplot(2,1,2)
simplePlot(t,g,'Risposta impulsiva','tempo[min]','concentrazione[pmol/l]');


%% Raw Deconvolution


% Costruzione vettore gk tramite la risoluzione
% dell'integrale di g(t) negli estremi:
% 0 - t1
% ti - ti+1
gk(1) = (s.A/s.alfa) * (1-exp(-s.alfa*t(1))) + (s.B/s.beta)*(1-exp(-s.beta*t(1))); 
for i = 1:length(t)-1    
    gk(i+1) = (s.A/s.alfa) * (exp(-s.alfa*t(i)) - exp(-s.alfa*t(i+1))) + ...
              (s.B/s.beta) * (exp(-s.beta*t(i)) - exp(-s.beta*t(i+1)));  
end
gk = gk';
gk = toeplitz(gk);

% Costruzione di G (matrice di convoluzione)
G = tril(gk);  

% Deconvolution
u_rd = pinv(G) * c;
% Reconvolution
y_rd = G*u_rd;

% Commento: la deconvoluzione e' un problema mal
% condizionato (rumore piccolo sull'uscita = rumore grande
% all'ingresso). Calcolo k per valutare l'entita' del
% malcondizionamento di cui la raw deconvolution non si
% cura.
k = norm(G) * norm(pinv(G));


% Plots
figure(2)

% Secrezione insulinica (stima ingresso)
subplot(2,1,1)
stairs(t,u_rd)
title('Stima secrezione insulinica - Raw deconvolution')
xlabel('tempo[min]')
ylabel('concentrazione[pmol/l]')

% Riconvolution vs Data
subplot(2,1,2)
plot(t,y_rd,t,c,'o')
title('Reconvolution vs Data')
xlabel('tempo[min]')
ylabel('concentrazione[pmol/l]')
legend('Reconvolution','Data')


%% Metodo della regolarizzazione

% Il metodo tiene conto del rumore delle misure. E'
% necessario quindi dare una misura dell'entita' della
% variabilita' dell'errore. Suppongo un CV del 5%.
CV = 0.05;
sigma_v = diag((CV.*(c + yb)).^2); 

% Penalizzo soluzioni poco regolari con la matrice di
% penalita'. Soluzioni poco regolari sono quelle con le
% derivate di punti vicine diverse tra loro
D = [1 -1 zeros(1,length(t) - 2)]';
D = tril(toeplitz(D));
P = D;   % grado = 1

% u_hat e' la soluzione del  problema ai minimi quadrati in
% cui il funzionale di costo e' la norma quadrata dei
% residui. Nel metodo della regolarizzazione appare nel
% funzionale un termine di regolarita' gamma.
gamma = 1000;   % Scelta arbitraria del termine di regolarita'
u_rg = pinv(G' * pinv(sigma_v) * G + gamma * P' * P) * G' *pinv(sigma_v) * c;
y_rg = G * u_rg;
res_rg = c - y_rg;
RSS = norm(res_rg).^2;

%% Criterio di discrepanza di Twomey

tol = 1e-9;
gamma_tw = 10000;
upper_thr = trace(sigma_v) + tol;

while RSS > upper_thr 
   
   gamma_tw = gamma_tw - 1; 
   % Stima di u tw
   u_tw = pinv(G' * pinv(sigma_v) * G + gamma_tw * P' * P) * G' *pinv(sigma_v) * c;
   % Reconvolution tw
   y_tw = G * u_tw;
   % Residuals tw
   res_tw = c - y_tw;
   % RSS tw
   RSS = norm(res_tw).^2;
   
end

figure('Name','Regolarizzazione - Twomey criterion','NumberTitle','off');
subplot(3,1,1)
plot(t,u_tw)
title('Reg - Twomey criterion')
xlabel('tempo[min]')
ylabel('concentrazione[pmol/l]')
subplot(3,1,2)
plot(t,y_tw)
title('Reconvolution')
xlabel('tempo[min]')
ylabel('concentrazione[pmol/l]')
subplot(3,1,3)
plot(t,res_tw)
title('Residuals')
xlabel('tempo[min]')



%% Generalized Cross Validation

% gamma inizializzazione
gamma_gcv = gamma_tw;   
% gamma per minimizzazione funzionale gcv
gamma_gcv = lsqnonlin(@(gamma)gcv(G,c,gamma,P,sigma_v),gamma_gcv);   
% Stima di u con gamma ottimizzato gcv 
u_gcv = pinv(G' * pinv(sigma_v) * G + gamma_gcv * P' * P)*G' * pinv(sigma_v) * c;
% Reconvolution y gcv
y_gcv = G * u_gcv;
% Residuals gcv
res_gcv  = c - y_gcv;

figure('Name','Regolarizzazione - Generalized Cross Validation','NumberTitle','off');
subplot(3,1,1)
plot(t,u_gcv)
title('Reg - Generalized Cross Validation')
xlabel('tempo[min]')
ylabel('concentrazione[pmol/l]')
subplot(3,1,2)
plot(t,y_gcv)
title('Reconvolution')
xlabel('tempo[min]')
ylabel('concentrazione[pmol/l]')
subplot(3,1,3)
plot(t,res_gcv)
title('Residuals')
xlabel('tempo[min]')

%% Maximum Likelihood

% gamma inizializzazione
gamma_ml = gamma_tw;   
% gamma per minimizzazione funzionale ml
gamma_ml = lsqnonlin(@(gamma)ml(G,c,gamma,P,sigma_v),gamma_ml);   
% Stima di u con gamma ottimizzato ml 
u_ml = pinv(G' * pinv(sigma_v) * G + gamma_ml * P' * P)*G' * pinv(sigma_v) * c;
% Reconvolution y gcv
y_ml = G * u_ml;
% Residuals gcv
res_ml  = c - y_ml;

figure('Name','Regolarizzazione - Maximum Likelihood','NumberTitle','off');
subplot(3,1,1)
plot(t,u_ml)
title('Reg - Maximum Likelihood')
xlabel('tempo[min]')
ylabel('concentrazione[pmol/l]')
subplot(3,1,2)
plot(t,y_ml)
title('Reconvolution')
xlabel('tempo[min]')
ylabel('concentrazione[pmol/l]')
subplot(3,1,3)
plot(t,res_ml)
title('Residuals')
xlabel('tempo[min]')


%% Plots

figure()

% Stima secrezione insulinica (u)
subplot(2,1,1)
stairs(t,u_rd)
hold on
plot(t,u_rg)
plot(t,u_tw)
plot(t,u_gcv)
plot(t,u_ml)
hold off
title('Stima secrezione insulinica')
xlabel('tempo[min]')
ylabel('concentrazione[pmol/l]')
legend('Raw Deconvolution','Reg - Twomey criterion','Reg - Generalized Cross Validation','Reg - Maximum Likelihood')

% Reconvolution (y_hat)
subplot(2,1,2)
plot(t,y_rd)
hold on
plot(t,y_rg)
plot(t,y_tw)
plot(t,y_gcv)
plot(t,y_ml)
plot(t,c,'o')
hold off
title('Reconvolution')
xlabel('tempo[min]')
ylabel('concentrazione[pmol/l]')
legend('Raw Deconvolution','Reg - Twomey criterion','Reg - Generalized Cross Validation','Reg - Maximum Likelihood','Data')








