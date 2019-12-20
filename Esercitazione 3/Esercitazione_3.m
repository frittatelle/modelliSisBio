%% ESERCITAZIONE 3

% Obbiettivo: Identificare i modelli a una, due e tre
% esponenziali valutando la precisione delle stime e quale
% approccio (LS / WLS) sia migliore basnadosi su diverse cifre di
% merito.

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

CV = 0.04;                     % Coefficiente di variazione errore di misura  
sigma_v = (CV* (z+zb)).^2;     % Matrice di covarianza delle misure 

%% Stima LS

% Initial guess
p0_1 = [10 0.05]; 
p0_2 = [4 0.05 5 0.2];   
p0_3 = [ 3 0.2 2 0.02 8 0.05];    

% Stima dei parametri (LS)
[p_hat1,wrss_1,r_1,~, ~, ~, jacobian1] = lsqnonlin(@(p) respred(p,t,z,1),p0_1);     % 1 exp
[p_hat2,wrss_2,r_2,~, ~, ~, jacobian2] = lsqnonlin(@(p) respred(p,t,z,2),p0_2);     % 2 exp
[p_hat3,wrss_3,r_3,~, ~, ~, jacobian3] = lsqnonlin(@(p) respred(p,t,z,3),p0_3);     % 3 exp

% Fitting con parametri stimati (LS)
[r1,y1] = respred(p_hat1,tsim,z,1);     % 1 exp
[r2,y2] = respred(p_hat2,tsim,z,2);     % 2 exp
[r3,y3] = respred(p_hat3,tsim,z,3);     % 3 exp

sigma_v = diag(sigma_v);

% Matrice di covarianza delle stime (LS)
sigma_p1 = inv(jacobian1'*inv(sigma_v)*jacobian1);   % 1 exp
sigma_p2 = inv(jacobian2'*inv(sigma_v)*jacobian2);   % 2 exp
sigma_p3 = inv(jacobian3'*inv(sigma_v)*jacobian3);   % 3 exp

% Calcolo CV delle stime (LS)
var_p1 = diag(sigma_p1);    % varianze p1
var_p2 = diag(sigma_p2);    % varianze p2
var_p3 = diag(sigma_p3);    % varianze p3

% CV p1
for i = 1:length(var_p1)
    SD_p1(i) = sqrt(var_p1(i));
    CV_p1(i) = 100*SD_p1(i)/p_hat1(i);
    CI_p1{i} = [p_hat1(i) - SD_p1(i),p_hat1(i) + SD_p1(i)];
end

% CV p2
for i = 1:length(var_p2)
    SD_p2(i) = sqrt(var_p2(i));
    CV_p2(i) = 100*SD_p2(i)/p_hat2(i);
    CI_p2{i} = [p_hat2(i) - SD_p2(i),p_hat2(i) + SD_p2(i)];
end

% CV p3
for i = 1:length(var_p3)
    SD_p3(i) = sqrt(var_p3(i));
    CV_p3(i) = 100*SD_p3(i)/p_hat3(i);
    CI_p3{i} = [p_hat3(i) - SD_p3(i),p_hat3(i) + SD_p3(i)];
end


% Indici di merito
N = length(t);          %numero dati
M1 = length(p_hat1);    %numero parametri
M2 = length(p_hat2);
M3 = length(p_hat3);

% 1
AIC1 = (wrss_1/N)+(2*M1/N);          %penalizza modelli con molti parametri che fittano male
SC1 = wrss_1+log(N);
MDL1 = wrss_1/N+log(N)*M1/N;
BIC1 = wrss_1+log(N)*M1;

% 2
AIC2 = (wrss_2/N)+(2*M2/N);
SC2 = wrss_2+log(N);
MDL2 = wrss_2/N+log(N)*M2/N;
BIC2 = wrss_2+log(N)*M2;

% 3
AIC3 = (wrss_3/N)+(2*M3/N);
SC3 = wrss_3+log(N);
MDL3 = wrss_3/N+log(N)*M3/N;
BIC3 = wrss_3+log(N)*M3;


%% Stima WLS 

sigma_v = diag(sigma_v);

% Initial guess
p0_1 = [10 1]; 
p0_2 = [3 0.05 10 0.2];   
p0_3 = [3 0.2 2 0.01 10 0.05];  

% Stima dei parametri (WLS)
[wp_hat1,wrss_1,wr_1,~, ~, ~, w_jacobian1] = lsqnonlin(@(p) respred(p,t,z,1,sigma_v),p0_1);     % 1 exp
[wp_hat2,wrss_2,wr_2,~, ~, ~, w_jacobian2] = lsqnonlin(@(p) respred(p,t,z,2,sigma_v),p0_2);     % 2 exp
[wp_hat3,wrss_3,wr_3,~, ~, ~, w_jacobian3] = lsqnonlin(@(p) respred(p,t,z,3,sigma_v),p0_3);     % 3 exp


% Fitting con parametri stimati (WLS)
[w_r1,w_y1] = respred(wp_hat1,tsim,z,1,sigma_v);     % 1 exp
[w_r2,w_y2] = respred(wp_hat2,tsim,z,2,sigma_v);     % 2 exp
[w_r3,w_y3] = respred(wp_hat3,tsim,z,3,sigma_v);     % 3 exp

% Matrice di covarianza delle stime (WLS)
% non moltiplico per sigma_v inversa perche' gia' pesato nel calcolo della funzione residui
sigma_wp1 = inv(w_jacobian1'*w_jacobian1);   % 1 exp
sigma_wp2 = inv(w_jacobian2'*w_jacobian2);   % 2 exp
sigma_wp3 = inv(w_jacobian3'*w_jacobian3);   % 3 exp

% Calcolo CV delle stime (WLS)
var_wp1 = diag(sigma_wp1);    % varianze p1
var_wp2 = diag(sigma_wp2);    % varianze p2
var_wp3 = diag(sigma_wp3);    % varianze p3

% CV wp1
for i = 1:length(var_wp1)
    SD_wp1(i) = sqrt(var_wp1(i));
    CV_wp1(i) = 100*SD_wp1(i)/wp_hat1(i);
    CI_wp1{i} = [wp_hat1(i) - SD_wp1(i),wp_hat1(i) + SD_wp1(i)];
end

% CV wp2
for i = 1:length(var_wp2)
    SD_wp2(i) = sqrt(var_wp2(i));
    CV_wp2(i) = 100*SD_wp2(i)/wp_hat2(i);
    CI_wp2{i} = [wp_hat2(i) - SD_wp2(i),wp_hat2(i) + SD_wp2(i)];
end

% CV wp3
for i = 1:length(var_wp3)
    SD_wp3(i) = sqrt(var_wp3(i));
    CV_wp3(i) = 100*SD_wp3(i)/wp_hat3(i);
    CI_wp3{i} = [wp_hat3(i) - SD_wp3(i),wp_hat3(i) + SD_wp3(i)];
end


% Indici di merito
N = length(t);
M1w = length(wp_hat1);
M2w = length(wp_hat2);
M3w = length(wp_hat3);

% 1w
AIC1_w = (wrss_1/N)+(2*M1/N);          %penalizza modelli con molti parametri che fittano male
SC1_w = wrss_1+log(N);
MDL1_w = wrss_1/N+log(N)*M1/N;
BIC1_w = wrss_1+log(N)*M1;

% 2w
AIC2_w = (wrss_2/N)+(2*M2/N);
SC2_w = wrss_2+log(N);
MDL2_w = wrss_2/N+log(N)*M2/N;
BIC2_w = wrss_2+log(N)*M2;

% 3w
AIC3_w = (wrss_3/N)+(2*M3/N);
SC3_w = wrss_3+log(N);
MDL3_w = wrss_3/N+log(N)*M3/N;
BIC3_w = wrss_3+log(N)*M3;


%% Plots - Modello 1 exp

res = zeros(size(r_1));
res_due = 2.*ones(size(r_1));
res_mindue = -2.*ones(size(r_1));


figure(1)

% LS
subplot(2,2,1)
plot(t, z, '*', tsim, y1)
grid
title('Stima LS - 1 esponenziale')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')
subplot(2,2,3)
plot(t,r_1,'o-',t,res,t,res_due,t,res_mindue)
grid
title('Residui LS - 1 esponenziale')

% WLS
subplot(2,2,2)
plot(t, z, '*', tsim, w_y1)
grid
title('Stima WLS - 1 esponenziale')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')
subplot(2,2,4)
plot(t,wr_1,'o-',t,res,t,res_due,t,res_mindue)
grid
title('Residui WLS - 1 esponenziale')

%% Plots - Modello 2 exp
figure(2)

% LS
subplot(2,2,1)
plot(t, z, '*', tsim, y2)
grid
title('Stima LS - 2 esponenziali')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')
subplot(2,2,3)
plot(t,r_2,'o-',t,res,t,res_due,t,res_mindue)
grid
title('Residui LS - 2 esponenziali')

% WLS
subplot(2,2,2)
plot(t, z, '*', tsim, w_y2)
grid
title('Stima WLS - 2 esponenziali')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')
subplot(2,2,4)
plot(t,wr_2,'o-',t,res,t,res_due,t,res_mindue)
grid
title('Residui WLS - 2 esponenziali')

%% Plots - Modello 3 exp
figure(3)

% LS
subplot(2,2,1)
plot(t, z, '*', tsim, y3)
grid
title('Stima LS - 3 esponenziali')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')
subplot(2,2,3)
plot(t,r_3,'o-',t,res,t,res_due,t,res_mindue)
grid
title('Residui LS - 3 esponenziali')

% WLS
subplot(2,2,2)
plot(t, z, '*', tsim, w_y3)
grid
title('Stima WLS - 3 esponenziali')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')
subplot(2,2,4)
plot(t,wr_3,'o-',t,res,t,res_due,t,res_mindue)
grid
title('Residui WLS - 3 esponenziali')


%% Plots - Confronto modelli LS/WLS
figure(4)

% LS
subplot(3,2,1)
plot(t, z, '*', tsim, y1)
grid
title('Stima LS - 1 esponenziale')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')

subplot(3,2,3)
plot(t, z, '*', tsim, y2)
grid
title('Stima LS - 2 esponenziali')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')

subplot(3,2,5)
plot(t, z, '*', tsim, y3)
grid
title('Stima LS - 3 esponenziali')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')


% WLS
subplot(3,2,2)
plot(t, z, '*', tsim, w_y1)
grid
title('Stima WLS - 1 esponenziale')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')

subplot(3,2,4)
plot(t, z, '*', tsim, w_y2)
grid
title('Stima WLS - 2 esponenziali')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')

subplot(3,2,6)
plot(t, z, '*', tsim, w_y3)
grid
title('Stima WLS - 3 esponenziali')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')

%% Plots - Confronto modelli LS/WLS semilog
figure(5)

% LS
subplot(3,2,1)
semilogy(t, z, '*', tsim, y1)
grid
title('Stima LS - 1 esponenziale')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')

subplot(3,2,3)
semilogy(t, z, '*', tsim, y2)
grid
title('Stima LS - 2 esponenziali')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')

subplot(3,2,5)
semilogy(t, z, '*', tsim, y3)
grid
title('Stima LS - 3 esponenziali')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')


% WLS
subplot(3,2,2)
semilogy(t, z, '*', tsim, w_y1)
grid
title('Stima WLS - 1 esponenziale')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')

subplot(3,2,4)
semilogy(t, z, '*', tsim, w_y2)
grid
title('Stima WLS - 2 esponenziali')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')

subplot(3,2,6)
semilogy(t, z, '*', tsim, w_y3)
grid
title('Stima WLS - 3 esponenziali')
xlabel('tempo [min]')
ylabel('concentrazione [pmol/ml]')
legend('data', 'model')
