%Modello tricompartimentale

%Parametri modello
V1 = 5;
k21 = 2.2;
k12 = 0.859;
k31 = 0.031;
k13 = 0.008;
k01 = 1.2;
dose = 500;

A = [-(k31+k21+k01) k12 k13 ; k21 -k12 0 ; k31 0 -k13];
B = [dose 0 0]';
C = [1/V1 0 0];
D = 0;

SYS = ss(A,B,C,D);
SYS.TimeUnit = 'hours';


% (a) Valutare il volume di distribuzione, la clearance e il tempo medio di residenza
%--------------------------------------------------------------------

%Risposta impulsiva del sistema

[y,t] = impulse(SYS);

figure(4);
plot(t,y);
title('Risposta Impulsiva sistema tricompartimentale');
xlabel('Tempo [ore]'),ylabel('Concentrazione [mg/l]');

% La percentuale piu' stringente che mi restituisce un risultato e' lo 0.1%
perc = 0.001;
y_max = max(y);
ind_almostzero = find(y<y_max*perc);
tempi_post = t(ind_almostzero);
t_eliminazione = tempi_post(1);


%-------------------------------------------------------------------

AUC = trapz(t,y);
AUMC = trapz(t,t.*y);

% MRT
MRT = AUMC./AUC;

% Clearance
CL = dose./AUC;

tau_par = 1/abs(eig(A));

% Volume di distribuzione
V_d = CL./k01;



