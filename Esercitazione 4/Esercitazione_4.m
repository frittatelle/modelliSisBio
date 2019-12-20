%% ESERCITAZIONE 4 - IDENTIFICAZIONE MODELLI COMPARTIMENTALI LINEARI

% Obbiettivo: Identificare il modello compartimentale del
% C-Peptide, ricavare le relazioni tra modello
% compartimentale ed esponenziale e confrontare le stime
% ottenute

clear

% Soggetto 1
[t,z] = textread('DatiCPsog1.dat','%d %f','headerlines',4);

%% Preprocessing

% Rimozione basale
zb = z(1);
z = z-zb;
% Rimozione misura al min = 1
z(1:2) = [];
t(1:2) = [];

tsim = [0:t(end)];

dose = 49650;

% Data 
data.x = t;
data.z = z;

CV = 0.04;                     % Coefficiente di variazione errore di misura  
sigma_v = (CV* (z+zb)).^2;     % Matrice di covarianza delle misure 

%% Stima LS - Parametri EXP

% Initial guess
p0_exp = [4 0.05 8 0.2];   

% Plot options - iter
data.objFun = @(x)respred_expode(x,t,z,2,0); 
data.title = 'Fitting animation EXP';
data.pauseTime = 1;
options = optimset('plotFcn',@(x,optimValues,state)optimplotfit(x,optimValues,state,data));
% Stima dei parametri (EXP)
% [p_hat_exp,wrss_exp,r_exp,~, ~, ~, jacobian_exp] = lsqnonlin(@(p) respred_expode(p,t,z,2,dose),p0_exp,[],[],options);
[p_hat_exp,wrss_exp,r_exp,~, ~, ~, jacobian_exp] = lsqnonlin(@(p) respred_expode(p,t,z,2,dose),p0_exp,[],[]);
% Fitting con parametri stimati (EXP)
[res_exp,y_exp] = respred_expode(p_hat_exp,t,z,2,dose);

A_hat_exp = p_hat_exp(1);
alfa_hat_exp = p_hat_exp(2);
B_hat_exp = p_hat_exp(3);
beta_hat_exp = p_hat_exp(4);

% Precisione stime EXP
sigma_v = diag(sigma_v);
sigma_p_hat_exp = pinv(jacobian_exp' * pinv(sigma_v) * jacobian_exp);
var_p_hat_exp = diag(sigma_p_hat_exp);
SD_p_hat_exp = sqrt(var_p_hat_exp);
% CI / CV 
for i = 1:length(SD_p_hat_exp)
   CI_exp{i} = [p_hat_exp(i) - SD_p_hat_exp(i) ; p_hat_exp(i) + SD_p_hat_exp(i)];
   CV_exp(i) = 100 * SD_p_hat_exp(i) / p_hat_exp(i);
end


% Plot EXP
res = zeros(size(res_exp));
res_due = 2.*ones(size(res_exp));
res_mindue = -2.*ones(size(res_exp));

figure()

subplot(2,2,1)
plot(t, z, '*', t, y_exp)
grid
title('Stima LS - EXP')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')
subplot(2,2,3)
plot(t,res_exp,'o-',t,res,t,res_due,t,res_mindue)
grid
title('Residui LS - EXP')


%% Stima LS - Parametri ODE

% Initial guess
p0_ode = [0.06 0.05 0.07 3];  

% Plot options - iter
data.objFun = @(x)respred_expode(x,t,z,1,dose); 
data.title = 'Fitting animation ODE';
data.pauseTime = 0.5;
options = optimset('outputfcn',@(x,optimValues,state)optimplotfit(x,optimValues,state,data));

% Stima dei parametri (ODE)
% [p_hat_ode,wrss_ode,r_ode,~, ~, ~, jacobian_ode] = lsqnonlin(@(p) respred_expode(p,t,z,1,dose),p0_ode,[],[],options);
[p_hat_ode,wrss_ode,r_ode,~, ~, ~, jacobian_ode] = lsqnonlin(@(p) respred_expode(p,t,z,1,dose),p0_ode,[],[]);
% Fitting con parametri stimati (ODE)
[res_ode,y_ode] = respred_expode(p_hat_ode,tsim,z,1,dose);

% Parametri ODE
k01 = p_hat_ode(1);
k12 = p_hat_ode(2);
k21 = p_hat_ode(3);
V = p_hat_ode(4);

% Precisione stime ODE
sigma_p_hat_ode = pinv(jacobian_ode' * pinv(sigma_v) * jacobian_ode);
var_p_hat_ode = diag(sigma_p_hat_ode);
SD_p_hat_ode = sqrt(var_p_hat_ode);
% CI / CV 
for i = 1:length(SD_p_hat_ode)
   CI_ode{i} = [p_hat_ode(i) - SD_p_hat_ode(i) ; p_hat_ode(i) + SD_p_hat_ode(i)];
   CV_ode(i) = 100 * SD_p_hat_ode(i) / p_hat_ode(i);
end


%% Simulated annealing

data.pauseTime = 0.1;
data.dose = dose;
objFun = @(x)simple_objective(x,data);
x0 = [4 0.05 8 0.2]; 
% options = optimoptions(@simulannealbnd,'display','iter');
[p_hat_SA,fval,exitFlag,output] = simulannealbnd(objFun,x0,-Inf,Inf);

r = respred_expode(p_hat_SA,t,z,2,0);

% Analisi residui
% POCHI DATI
% h = 1 --> rifiuto
% Runstest randomness
%  h0: random order
[h_random,p_random] = runstest(r);
% Anderson-Darling normalita'
%  h0: campione proveniente da dist normale
[h_normal,p_normal] = adtest(r);

[r,y] = respred_expode(p_hat_SA,t,z,2,0);
% Plot EXP
res = zeros(size(r));
res_due = 2.*ones(size(r));
res_mindue = -2.*ones(size(r));


subplot(2,2,2)
plot(t, z, '*', t, y)
grid
title('Stima Simulated Annealing - EXP')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')
subplot(2,2,4)
plot(t,r,'o-',t,res,t,res_due,t,res_mindue)
grid
title('Residui Simulated Annealing - EXP')




%% Espressione analitica dell'uscita 

% Espressione analitica dell'uscita da ODE per confronto con modello EXP
% 3 modi:
% 1. Carta penna - Laplace Transform - Antitrasformata
% 2. Metodo dei fratti semplici 
% 3. Espressione analitica uscita con dsolve

% Fratti semplici

A = [-(k01+k21),k12;k21,-k12];
B = [dose;0];
C = [1/V 0];
sys = ss(A,B,C,0);

g = tf(sys);
N=cell2mat(g.num);
D=cell2mat(g.den);

[R,P]=residue(N,D);

alfa_hat_fs = - P(2); 
beta_hat_fs = - P(1); 
A_hat_fs = R(2);
B_hat_fs = R(1);

% Dsolve

syms x1(t) x2(t) y(t) k01 k12 k21 V

eqn_x1 = [diff(x1,t) == -(k01+k21).*x1 + k12.*x2];
eqn_x2 = [diff(x2,t) == k21.*x1 - k12.*x2];
eqn_y = [diff(y,t) == x1./V];

eqns = [eqn_x1,eqn_x2,eqn_y];

sol = dsolve(eqns);

k01 = p_hat_ode(1);
k12 = p_hat_ode(2);
k21 = p_hat_ode(3);

beta_hat_dsolve = ((k01 + k12 + k21 + (k01^2 - 2*k01*k12 + 2*k01*k21 + k12^2 + 2*k12*k21 + k21^2)^(1/2)))/2;
alfa_hat_dsolve = ((k01 + k12 + k21 - (k01^2 - 2*k01*k12 + 2*k01*k21 + k12^2 + 2*k12*k21 + k21^2)^(1/2)))/2;

%% Monte Carlo

% Genero un set di dati di numerosita' pari a quella di
% partenza a partire dai parametri stimati

% EXP 
p_hat_exp0 = p_hat_exp;
t = data.x;
P_MC = [];
P_MC = P_MC';
y_MC = y_exp;
y = y_exp;

CV = 0.04;
% maxIter = 10000;
maxIter = 100;
for i = 1:maxIter
   
   noise = normrnd(0, CV .* y);
   y = y_exp + noise;
  
   % Stima dei parametri (EXP,MC)
   p_hat_mc = lsqnonlin(@(p) respred_expode(p,t,y,2,dose),p_hat_exp0,[],[]);
   [~,y_MC] = respred_expode(p_hat_mc,t,y,2,dose);
   
   P_MC = vertcat(P_MC,p_hat_mc);
   
end

% Dalle distribuzioni simulate dei parametri che ho generato
% estraggo media e SD per calcolare la precisione delle stime. 
% Calcolo inoltre i CV per confrontare le stime con il metodo LS
figure()
for i = 1:size(P_MC,2)
   SD_p_hat_MC{i} = std(P_MC(:,i));
   mean_p_hat_MC{i} = mean(P_MC(:,i));
   CI_MC{i} = [mean_p_hat_MC{i} - SD_p_hat_MC{i},mean_p_hat_MC{i} + SD_p_hat_MC{i}];
   CV_MC(i) = 100 * SD_p_hat_MC{i} ./ mean_p_hat_MC{i};
   subplot(2,2,i)
   histfit(P_MC(:,i));
end

% Per qualche motivo i CV mi vengono piu' grandi degli altri
% metodi di identificazione. Distribuzione non normale dei
% parametri?
% MA HA SENSO CALCOLARE I CV IN QUESTO CASO? 
% La simulazione mi serve per costruire gli intervalli di confidenza 
% e per quantificare la precisione delle stime (fatte con altri metodi)

% Test KS per verificare la normalita' dei parametri
% H = 0 --> viene probabilmente da una normale
for i = 1:size(P_MC,2)
   [h,p] = kstest(P_MC(:,i));
   H{i} = h;
   P_VALUE{i} = p;
end


[r_MC,y_MC] = respred_expode(p_hat_mc,t,y,2,dose);

% Plot 
res = zeros(size(r_MC));
res_due = 2.*ones(size(r_MC));
res_mindue = -2.*ones(size(r_MC));

figure()
subplot(2,1,1)
plot(t,y_MC,t,z,'*')
grid
title('Stima Monte Carlo - EXP')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')
subplot(2,1,2)
plot(t,r_MC,'o-',t,res,t,res_due,t,res_mindue)
grid
title('Residui Monte Carlo - EXP')
