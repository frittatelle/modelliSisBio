%Modello monocompartimentale - no abs

% V = [l]
% k = [ore^-1]
% Dose = [mg]

%Parametri modello
V1 = 5;
k01 = 1.2;
dose = 50;

A = -k01;
B = dose;
C = 1/V1;
D = 0;

SYS = ss(A,B,C,D);
SYS.TimeUnit = 'hours';

[y,t] = impulse(SYS);

%--------------------------------------------------------------------

% Si ipotizzi ora che il farmaco venga somministrato
% con una serie di boli (es. 12) distanziati gli uni dagli altri da un tempo tau (es. 6 ore).
% Provare a cambiare l’intervallo di tempo tra due somministrazioni (es. 6, 3, 1, 0.5, 0.2, ore) e il numero
% di somministrazioni mantenendo costante la durata del trattamento (3 giorni).

%   SOL 1 - PULSE TRAIN
%   PulseTrain
%   n_impulsi = numero di impulsi
%   tau = intervallo temporale tra gli impulsi
%   t_eliminazione = tempo di eliminazione del sistema
%   t = vettore risposta impulsiva del sistema
%   durata_sim = durata simulazione in ore

%   tau ha la priorita' sul numero di impulsi
%   ex. vorrei 10 impulsi con tau = 6 in una sim da 24 ore

days = 1;
sim_time = days*24;

U = pulseTrain(50,0.2,5,t,sim_time);
T = linspace(0,sim_time,length(U));
DT = T(2);

[y,t] = lsim(SYS,U,T);
% La risposta deve essere riscalata per il periodo di campionamento
% (lsim)
y = y/DT;

figure(5)
subplot(2,2,1)
plot(t,y);
title('Dosi ripetute - n = 50, tau = 0.2, 24H')
xlabel('Tempo [ore]'),ylabel('Concentrazione [mg/l]');



%   SOL 2 - SUPERPOSITION

%--------------------------------------------------------------------

% Regime

% Ipotizzando un numero infinito di somministrazioni, esiste un livello (medio) di regime?
% Qual `e il tempo necessario affinch`e il sistema arrivi a regime?


% Stima valore regime
perc = 0.63;
% Trovo i valori maggiori del 63% del valore max di concentrazione
ind_y_ss = find(y > max(y)*perc);
y_ss = y(ind_y_ss);
% Valor medio dei dati appena trovati
ss_val = mean(y_ss);
hold on
plot(t,ones(length(t),1)*ss_val);
legend('Concentrazione','SS stimato (63% y_m_a_x)');


% Valore regime da segnale filtrato
subplot(2,2,2)
% Plot segnale filtrato per visualizzare andamento
% e concentrazione a regime
y_lpf = lowpass(y,0.001);
plot(t,y_lpf);
title('Dosi ripetute - LPF')
xlabel('Tempo [ore]'),ylabel('Concentrazione [mg/l]');

% Moda segnale filtrato = concentrazione a regime
y_lpf_mean = mean(y_lpf(74:346));

hold on
plot(t,ones(length(t),1)*y_lpf_mean);
legend('Concentrazione - LPF','SS');

hold off

% Errore di stima
err_stima = ss_val./y_lpf_mean;


%-------------------------------------------------------------------



% Sweep tau,n

days = 3;
sim_time = days*24;

%vettore valori di tau x sweep
sweep_n = [6 3 1 0.5 0.2];

subplot(2,2,3)

tau_legend = {};

for i = 1:length(sweep_n)
    U = pulseTrain(20,sweep_n(i),5,t,sim_time);
    T = linspace(0,sim_time,length(U));
    DT = T(2);
    
    [y,t] = lsim(SYS,U,T);
    y = y/DT;
    
    plot(t,y);
    
    tau_legend{end+1} = num2str(sweep_n(i));
    
    hold on
end

title('tau Sweep - n = 20, days = 3');
xlabel('Tempo [ore]'),ylabel('Concentrazione [mg/l]');
legend(strcat('tau= ',tau_legend));

hold off


%vettore valori di n x sweep
sweep_n = linspace(5,50,10);

subplot(2,2,4)

n_legend = {};

for i = 1:length(sweep_n)
    U = pulseTrain(sweep_n(i),0.5,5,t,sim_time);
    T = linspace(0,sim_time,length(U));
    DT = T(2);
    
    [y,t] = lsim(SYS,U,T);
    y = y/DT;
    
    plot(t,y);
    
    n_legend{end+1} = num2str(sweep_n(i));
    
    hold on
end

title('n Sweep - tau = 0.5, days = 3');
xlabel('Tempo [ore]'),ylabel('Concentrazione [mg/l]');
legend(strcat('n= ',n_legend));

hold off


%-------------------------------------------------------------------

% Volendo raggiungere nel modo piu rapido possibile un livello plasmatico medio di 50 mg/l
% con una variazione del ±10%, volendo mantenere inoltre quei livelli di concentrazione per
% 7 giorni, quale schema di somministrazione proporresti (boli ripetuti con eventuale dose di
% attacco nel primo periodo)?

% AUC/tau = y_ss;

days = 7;
sim_time = 24*days;
% tau, n, dose d'attacco, dose?

y_toxic = 50+50*(0.1);
y_terap = 50-50*(0.1);

% Risposta al bolo per calcolo di AUC
[y_imp,t_imp] = impulse(SYS);

%Concentrazione di steady-state --> tau
y_ss = 50;
AUC = trapz(t_imp,y_imp);
%tau steady state
tau_ss = AUC./y_ss;


% Dose per SS
dose = 50;

B = dose;
SYS = ss(A,B,C,D);
SYS.TimeUnit = 'hours';

[y_imp,t_imp] = impulse(SYS);

tau = 1/abs(eig(A));

% La dose di attacco non mi serve per determinare il valore di regime.
% Puo' sostanzialmente diminuire il tempo con cui raggiungo il regime.
% Il valore di regime e' determinato in pratica dal tempo tra due
% somministrazioni.
% PROBLEMA: probabilmente, a causa del ceil(tau) di pulseTrain non riesco
% mai a raggiungere il valore di regime desiderato. (Sono sempre poco sopra o poco sotto)
% il valore corretto dovrebbe essere compreso tra 0.1 e 0.2: ~0.16.


% Per mostrare lo stesso l'andamento desiderato, suppongo che il tempo di
% eliminazione sia 4.1 tau invece che 5 tau
% Come tau tra una somministrazione e l'altra prendo 0.1

U = pulseTrain(1000,0.1,tau*4.1,t_imp,sim_time,0.1,V1);
T = linspace(0,sim_time,length(U));
DT = T(2);

[y,t] = lsim(SYS,U,T);
y = y/DT;

plot(t,y);
title('Schema somministrazione: n = 1000, tau = 0.1');
xlabel('Tempo [ore]'),ylabel('Concentrazione [mg/l]');

%-------------------------------------------------------------------

% Soluzione alternativa con gensig

%La velocità di eliminazione dipende da k01 costante di eliminazione,
%dalla frequenza del bolo (se sono molto ripetute l'organismo non fa in
%tempo ad eliminarne prima che la dose successiva venga somministrata e
%quindi si assesterà su un livello di concentrazione stazionario (Css))
%e dall'emivita della particolare sostanza.

d = 500;
%Il livello medio di regime (Css) lo calcolo ad esempio con tau=0.2, posso
%agire solo su d e tau poichè il resto dipende da farmaco e paziente
F = 1;         %somministrazione endovena: tutto il farmaco viene assorbito
ke = k01;
Css = (F*d)/(V1*ke*0.2);

%In quanto tempo arriva a regime?
%Fss = 1-e^(-tau_regime*ke) tempo per raggiungere il 99% di Css -->  tau_regime = log(1/(1-Fss))/ke
tau_regime = log(1-(1-0.99*Css))/ke;    %dipende solo da ke

%Dosi ripetute con dose d'attacco per arrivare a 50mg/L con una variazione
%di +-10% (quindi [45-55]) e trattamento di 7 giorni
durata = 7*24;               %durata trattamento
Css2 = 50;
ampiezza_finestra = 10;      %[55-45]
d2=V1*ampiezza_finestra;     %dose terapeutica (quella che viene ripetuta)
tau = d2/(V1*ke*Css2);     %oppure tau=log(55/45)/ke sapendo quanto devono essere cmax[55] e cmin[45] faccio il loro rapporto
n_boli = ceil(durata/tau);   %numero intero di somministrazioni
Tc = tau/10;

%Calcolo dose attacco
dose_attacco = d2/(1-exp(-ke*tau));   %se usassi d2 cambierebbe di molto poco il tempo per arrivare a regime
B = dose_attacco;                     %la dose d'attacco sarà l'ingresso del nuovo sistema
sist = ss(A,B,C,D);
[y_att] = impulse(sist,0:Tc:180);     %prendo un po' piu' di 7*24 cosi vedo la discesa

%Calcolo dosi ripetute
B = d2;
sist = ss(A,B,C,0);
[y,t] = impulse(sist,0:Tc:180);
indice_t = find(t==tau);            %indice della prima somministrazione (dove finisce)
c=zeros(1,length(y));               %vettore di zeri
b=indice_t-1;                       %momento precedente la somministrazione
Y=zeros(length(y),1);               %preparo il vettore y di zeri

for i=2:n_boli               %perche' il primo e' quello della dose di attacco
    for j=1:b
        c(j)=0;                 %niente somministrazione
    end
    for j=(b+1):length(y)
        c(j)=y(j-b);            %metto le concentrazioni y scalate
    end
    Y=Y+c';                 %aggiorno il vettore delle concentrazioni
    b=(indice_t-1)*i;       %ogni quanto devo fare la dose di attacco (11,21,31...)
end

figure(6)

subplot(2,1,1)
plot(t,y_att+Y)
title('Boli ripetuti con dose di attacco (modello monocompartimentale)');
xlabel('Tempo [ore]')
ylabel('Concentrazione [mg/L]')

%Utilizzo il modello monocompartimentale con assorbimento del punto 2
V2 = 5;
k01 = 1.2;
k02 = 1.2;
k21 = 2.2;
ka = k01+k21;
F = k21/(k01+k21);
d3 = 50/F;

A = [-ka 0;k21 -k02];
B = [d3; 0];
C = [0 1/V2];
D = 0;
sist = ss(A,B,C,D);

[y2,t] = impulse(sist,0:Tc:180);   %risposta del sistema al singolo bolo che sommo ripetutamente
indice_t = find(t==tau);           %tau è quello di prima
n_boli = ceil(durata/tau);

Y2=zeros(length(y2),1);
b=indice_t-1;
c=zeros(1,length(y2));

for i=2:n_boli
    for j=1:b
        c(j)=0;
    end
    for j=(b+1):length(y2)
        c(j)=y2(j-b);
    end
    Y2=Y2+c';
    b=(indice_t-1)*i;
end

dose_attacco = (d3*F)/(1-exp(-k02*tau));  %devo tenere conto di F
B = [dose_attacco; 0];
sist = ss(A,B,C,D);
[y_att2] = impulse(sist,0:Tc:180);        %calcolo la risposta impulsiva dovuta alla dose d'attacco, per aggiungerla alla normale somministrazione


subplot(2,1,2)
plot(t,y_att2+Y2)
title('Boli ripetuti con dose di attacco (modello monocompartimentale con assorbimento)')
xlabel('Tempo [ore]')
ylabel('Concentrazione [mg/L]');

