%Modello monocompartimentale - no abs

% dx/dt = Ax(t) + Bu(t)
%  y(t) = Cx(t) + Du(t)

% V = [l]
% k = [ore^-1]
% Dose = [mg] 

%Parametri modello
V1 = 5;
k01 = 1.2;
dose = 500;

A = -k01;
B = dose;
C = 1/V1;
D = 0;

SYS = ss(A,B,C,D);
SYS.TimeUnit = 'hours';

%--------------------------------------------------------------------

%Risposta impulsiva del sistema

%In quanto tempo viene sostanzialmente eliminato il farmaco? 
%A quale parametro del modello `e collegato tale tempo?

[y,t] = impulse(SYS);
% y = valori di concentrazione
% t = tempo trascorso

figure(1);
subplot(3,1,1),plot(t,y);
title('Risposta Impulsiva sistema monocompartimentale');
xlabel('Tempo [ore]'),ylabel('Concentrazione [mg/l]');
% ALTERNATIVA: impulseplot(SYS);

% tempo eliminazione quando è sotto lo 0.1% y_max
% Sono neccessarie 10 emivite per eliminare il 99,9% del farmaco

% La percentuale piu' stringente che mi restituisce un risultato e' lo 0.4%
perc = 0.004;
y_max = max(y);
ind_almostzero = find(y<y_max*perc);
tempi_post = t(ind_almostzero);
t_eliminazione = tempi_post(1);

%--------------------------------------------------------------------

%Tempo di emivita

% Qual e l’emivita (tempo in cui la quantita raggiunge il 50% del valore iniziale)?
%A che parametro del modello `e collegata?

[t_emivita,y_emivita] = half_life(y,t);
hold on
plot(t_emivita,y_emivita,'r.','LineWidth',15);
text(t_emivita,y_emivita,'   emivita');


%--------------------------------------------------------------------

%Costanti di tempo
subplot(3,1,2),semilogy(t,y);
title('Andamento concentrazioni - log, Costanti tempo');
xlabel('Tempo [ore]'),ylabel('Concentrazione [mg/l]');


%--------------------------------------------------------------------

%Concentrazioni

% Qual `e il valore iniziale di concentrazione e la concentrazione massima misurabile? 
% A quali parametri del modello sono legati tali valori?

concentrazione_iniziale = y(1);
concentrazione_max = y(1);

%--------------------------------------------------------------------

%AUC
% Valutare l’area sotto la curva di concentrazione detta comunemente AUC
% (consiglio usare l’istruzione trapz). 
% A quali parametri del modello `e legato tale valore?

AUC = trapz(t,y);
AUMC = trapz(t,t.*y);

tau = 1/abs(eig(A));

AUC_inf = AUC + y(end)*tau;
AUMC_inf = AUMC + t(end)*y(end)*tau+y(end)*tau.^2;

MRT = AUMC_inf./AUC_inf;

%--------------------------------------------------------------------
% Valutare il tempo medio di residenza con un approccio non compartimentale

% Simulazione con dati empirici (pochi punti)
% simulazione su otto ore con 4 campionamenti
T = linspace(0,8,8);
Y = impulse(SYS,T);

% Calcolo AUC su dati empirici
AUC_emp = trapz(T,Y);
AUMC_emp = trapz(T,T.*Y);

tau = 1/abs(eig(A));

AUC_inf_emp = AUC_emp + Y(end)*tau;
AUMC_inf_emp = AUMC_emp + T(end)*Y(end)*tau+Y(end)*tau.^2;

MRT_emp = AUMC_inf_emp./AUC_inf_emp;

% rapporto stima/reale AUC
r_AUC = AUC_inf_emp./AUC;


%--------------------------------------------------------------------

%CL tot
% Calcolare la clearance totale

CL_tot = dose./AUC;


%--------------------------------------------------------------------

%K01 sweep
% Variare il valore dei parametri (uno alla volta) e verificare L’impatto sulle grandezze
% calcolate.

%vettore valori di k01 x sweep
sweep = linspace(k01./5,k01.*5,10);

%plot grafico con k01 diversi
subplot(3,1,3)

%Ciclo sui valori di sweep.
%Definisco un nuovo sistema per ogni k01.
%Plotto la risposta all'impulso.

k01_legend = {};

for i = 1:length(sweep)       
    A = -sweep(i);
    SYS = ss(A,B,C,D);
    SYS.TimeUnit = 'hours';
    
    [y,t] = impulse(SYS);
    plot(t,y);
    
    k01_legend{end+1} = num2str(sweep(i));
    
    hold on   
end

title('k01 Sweep');
xlabel('Tempo [ore]'),ylabel('Concentrazione [mg/l]');
legend(strcat('k01= ',k01_legend));

hold off



