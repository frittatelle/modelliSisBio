%Modello monocompartimentale - abs

%Parametri modello
V2 = 5;
k01 = 1.2;
k02 = 1.2;
k21 = 2.2;
dose = 500;

A = [-(k01+k21) 0; k21 -k02];
B = [dose 0]';
C = [0 1/V2];
D = 0;

SYS = ss(A,B,C,D);
SYS.TimeUnit = 'hours';

%--------------------------------------------------------------------

%Risposta impulsiva del sistema

%In quanto tempo viene sostanzialmente eliminato il farmaco? 
%A quale parametro del modello `e collegato tale tempo?

[y,t] = impulse(SYS);

figure(2);
subplot(3,1,1),plot(t,y);
title('Risposta Impulsiva sistema monocompartimentale (abs)');
xlabel('Tempo [ore]'),ylabel('Concentrazione [mg/l]');

%--------------------------------------------------------------------
%Tempo eliminazione farmaco

% Vorrei tempo eliminazione quando è sotto lo 0.1% y_max
% La percentuale piu' stringente che mi restituisce un risultato e' lo 0.2%

perc = 0.002;

[y_max,ind_y_max] = max(y);

ind_almostzero = find(y<y_max*perc);
ind_almostzero = ind_almostzero(ind_almostzero>ind_y_max);
tempi_post = t(ind_almostzero);
t_eliminazione = tempi_post(1);

%--------------------------------------------------------------------

%Concentrazioni

% Qual `e il valore iniziale di concentrazione e la concentrazione massima misurabile? 
% A quali parametri del modello sono legati tali valori?

hold on

plot(t(ind_y_max),y_max,'r.','LineWidth',15);
text(t(ind_y_max),y_max,'      Concentrazione massima');

concentrazione_iniziale = y(1);
concentrazione_max = y_max;

%--------------------------------------------------------------------

%AUC
% Valutare l’area sotto la curva di concentrazione detta comunemente AUC
% (consiglio usare l’istruzione trapz). 
% A quali parametri del modello `e legato tale valore?

AUC = trapz(t,y);

%--------------------------------------------------------------------

%F
% Quale frazione di farmaco passa dal compartimento 1 al compartimento 2? A quali parametri
% del modello `e legato tale valore?

F = k21./(k01+k21);

%--------------------------------------------------------------------

%k01 Sweep
% Provare a variare il valore di k01 (es. provare con una scala di 10 valori compresi tra k01/100
% e 100k01), lasciando fissi i valori di tutti gli altri parametri, cosa cambia?

%vettore valori di k01 x sweep
sweep = linspace(k01./100,k01.*100,10);

%plot grafico con k01 diversi
subplot(3,1,2)

%Ciclo sui valori di sweep.
%Definisco un nuovo sistema per ogni k01.
%Plotto la risposta all'impulso.

k01_legend = {};

for i = 1:length(sweep)       
    A = [-(sweep(i)+k21) 0; k21 -k02];
    SYS = ss(A,B,C,D);
    SYS.TimeUnit = 'hours';
    
    [y,t] = impulse(SYS);
    plot(t,y);
    
    k01_legend{end+1} = num2str(sweep(i));
    
    hold on   
end

title('k01 Sweep');
xlabel('Tempo [ore]'),ylabel('Concentrazione [mg/ml]');
legend(strcat('k01= ',k01_legend));

hold off

%Torno al modello iniziale
A = [-(k01+k21) 0; k21 -k02];
B = [dose 0]';
C = [0 1/V2];
D = 0;

SYS = ss(A,B,C,D);
SYS.TimeUnit = 'hours';
[y,t] = impulse(SYS);

% Cambiando i valori di k01 si puo` notare un cambiamento dell`andamento della concentrazione.\n');
% In particolare, essendo k01 il parametro che regola l`espulsione da C1,si puo` notare come per valori bassi, 
% la concentrazione max in C2 raggiunga valori piu` alti.;
% La velocita` di espulsione e` bassa e una quantita` maggiore di farmaco raggiunge C2

%--------------------------------------------------------------------

%Costanti di tempo
subplot(3,1,3),semilogy(t,y);
title('Andamento concentrazioni - log, Costanti tempo');
xlabel('Tempo [ore]'),ylabel('Concentrazione [mg/ml]');
    
%--------------------------------------------------------------------
%Valutare la biodisponibilita' attraverso un approccio non compartimentale.

D_iv = 500;
D_oral = dose;
AUC_oral = AUC;
AUC_iv = 83.096;
F_biodisponibilita = (D_iv*AUC_oral)./(D_oral*AUC_iv);

%--------------------------------------------------------------------
%Valutare la costante di assorbimento apparente con un approccio non compartimentale.

% Se l’assorbimento è un first-order process, si assume
% 1/ka=MAT (mean absorbtion time) dove MAT=MRTniMRTiv ossia la differenza tra il tempo medio di
% residenza in somministrazione non istantanea (orale o
% intramuscolare) e quello a seguito di
% somministrazione endovenosa.

AUMC_oral = trapz(t,t.*y);
AUMC_iv = 68.005;

MRT_oral = AUMC_oral./AUC_oral;
MRT_iv = 0.833;

MAT = MRT_oral*MRT_iv;
ka = 1/MAT;



