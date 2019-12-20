%Modello compartimentale a due compartimenti (caso di studio C-peptide)

%Parametri modello
% V = [l]
% k = [minuti^-1]
% Dose = [pmol] 

V1 = 3.29;
k01 = 6.54e-2;
k12 = 5.68e-2;
k21 = 7.16e-2;
dose = 49650;

A = [-(k01+k21) k12; k21 -k12];
B = [dose 0]';
C = [1/V1 0];
D = 0;

SYS = ss(A,B,C,D);
SYS.TimeUnit = 'minutes';

%--------------------------------------------------------------------

%Risposta impulsiva del sistema

%In quanto tempo viene sostanzialmente eliminato il C-peptide? 
%A quali parametri del modello `e collegato tale tempo?

[y,t] = impulse(SYS);

figure(3);
subplot(3,1,1),plot(t,y);
title('Risposta Impulsiva sistema bicompartimentale');
xlabel('Tempo [minuti]'),ylabel('Concentrazione [pmol/l]');

% La percentuale piu' stringente che mi restituisce un risultato e' lo 0.5%
perc = 0.005;
y_max = max(y);
ind_almostzero = find(y<y_max*perc);
tempi_post = t(ind_almostzero);
t_eliminazione = tempi_post(1);


%--------------------------------------------------------------------

%Concentrazioni

% Qual `e il valore iniziale di concentrazione e la concentrazione massima misurabile? 
% A quali parametri del modello sono legati tali valori?

hold on

[y_max,ind_y_max] = max(y);

concentrazione_iniziale = y(1);
concentrazione_max = y_max;

%--------------------------------------------------------------------

%AUC
% Valutare l’area sotto la curva di concentrazione detta comunemente AUC
% (consiglio usare l’istruzione trapz). 
% A quali parametri del modello `e legato tale valore?

AUC = trapz(t,y);

%--------------------------------------------------------------------

%Costanti di tempo
subplot(3,1,3),semilogy(t,y);
title('Andamento concentrazioni - log, Costanti tempo');
xlabel('Tempo [minuti]'),ylabel('Concentrazione [pmol/l]');

%--------------------------------------------------------------------

%CL tot
% Calcolare la clearance totale

CL_tot = dose./AUC;

%--------------------------------------------------------------------

%k21 Sweep
% Provare a variare il valore di k21 (es. provare con una scala di 10 valori compresi tra k21/100
% e 100k12), lasciando fissi i valori di tutti gli altri parametri, cosa cambia?


%vettore valori di k21 x sweep
sweep = linspace(k21./100,k12.*100,10);

%plot grafico con k21 diversi
subplot(3,1,2)

%Ciclo sui valori di sweep.
%Definisco un nuovo sistema per ogni k01.
%Plotto la risposta all'impulso.

k21_legend = {};

for i = 1:length(sweep)       
    A = [-(k01+sweep(i)) k12; sweep(i) -k12];
    SYS = ss(A,B,C,D);
    SYS.TimeUnit = 'hours';
    
    [y,t] = impulse(SYS);
    plot(t,y);
    
    k21_legend{end+1} = num2str(sweep(i));
    
    hold on   
end

title('k21 Sweep');
xlabel('Tempo [ore]'),ylabel('Concentrazione [mg/l]');
legend(strcat('k01= ',k21_legend));

hold off
