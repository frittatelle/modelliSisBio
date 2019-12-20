%% ESERCITAZIONE 5 - REAZIONI ENZIMATICHE

% Obbiettivo: Simulare l'evoluzione nel tempo di un sistema in cui
% avvengono reazioni enzimatiche imponendo diverse
% condizioni sui parametri di reazione


%% Set 1

k1 = 1;     %[mM^-1 * sec^-1]
k2 = 1;     %[sec^-1]
k_1 = 1;    %[sec^-1]
k_2 = 1;    %[mM^-1 * sec^-1]

s0 = 1;       %[mM]
e0 = 100;     %[mM]
c0 = 0;
p0 = 0;

cond0 = [s0 e0 c0 p0];

t = [0:0.001:0.1];

[t,y] = ode45(@(t,y)reazioneEnzimatica(y,k1,k_1,k2,k_2),t,cond0);
s = y(:,1);
e = y(:,2);
c = y(:,3);
p = y(:,4);

% Plots
simplePlot = @simplePlot;
labels = {' Substrato',' Enzima',' Complesso',' Prodotto'};
tLabel = 'tempo [s]';

% Andamento concentrazioni 
cLabel = 'concentrazione [mMoli]';
cTitle = 'Andamento concentrazione';
figure
for i = 1:length(labels)
   subplot(3,2,i)
   simplePlot(t,y(:,i),strcat(cTitle,labels{i}),tLabel,cLabel);
end

% Velocita' formazione del prodotto 
dp = k2.*c -k_2.*e.*p;
vLabel = 'velocita` [mMoli*s^{-1}]';
vTitle = 'Velocita` formazione Prodotto';
subplot(3,2,5)
simplePlot(t,dp,vTitle,tLabel,vLabel);

% - Crescita/decrescita prodotto/substrato
% Il substrato ha una decrescita di tipo esponenziale,
% mentre il proddotto ha una crescita di tipo sigmoidale(?).
% Questo e' probabilmente dovuto a k-2 diverso da 0.
% L'interazione tra enzima e substrato potrebbe essere
% cooperativa, pertanto il prodotto seguirebbe una cinetica
% di tipo Hill invece che di tipo Michaelis-Menten

% - Regime? Quando? Rapporto?
% Substrato, Composto, Enzima, Prodotto raggiungono tutti una condizione
% di regime, dopo circa 0.06 secondi

% - Enzima a regime
% A regime l'enzima si trova al 99,02% sotto forma di enzima semplice, e
% all' 0.98% sotto forma di composto. Il tutto e' dovuto a
% k-2 diverso da 0. Se esso fosse uguale a 0 mi aspetterei
% l'enzima sia al 100% in forma di enzima


%% Set 2

k1 = 1;     %[mM^-1 * sec^-1]
k2 = 1;     %[sec^-1]
k_1 = 1;    %[sec^-1]
k_2 = 1;    %[mM^-1 * sec^-1]

s0 = 100;   %[mM]
e0 = 1;     %[mM]
c0 = 0;
p0 = 0;

cond0 = [s0 e0 c0 p0];

t = [0:0.001:0.1];

[t,y] = ode45(@(t,y)reazioneEnzimatica(y,k1,k_1,k2,k_2),t,cond0);

s = y(:,1);
e = y(:,2);
c = y(:,3);
p = y(:,4);

% Plots

% Andamento concentrazioni 
figure
for i = 1:length(labels)
   subplot(3,2,i) 
   simplePlot(t,y(:,i),strcat(cTitle,labels{i}),tLabel,cLabel);
end

% Velocita' formazione del prodotto 
dp = k2.*c -k_2.*e.*p;
subplot(3,2,5)
simplePlot(t,dp,vTitle,tLabel,vLabel);

% - Crescita/decrescita prodotto/substrato
% Dopo una prima fase di crescita esponenziale, il prodotto
% cresce linearmente. Il substrato dopo una fase di
% decrescita esponenziale, decresce linearmente.

% - Regime? Quando? Rapporto?
% Solo enzima e composto, 0, 1.

% - Enzima a regime
% L'enzima a regime si trova sotto forma di composto .

% - Differenze con 5.1
% Enzima e composto mantengono lo stesso andamento.
% Substrato e prodotto invece, seguono all'inizio una
% cinetica di tipo espnonenziale per poi crescere o
% decrescere linearmente.


%% Set 3

k1 = 75e3;
k_1 = 75;
k2 = 600e3;
k_2 = 0;

s0= 100;
e0 = 1; 
c0 = 0;
p0 = 0;

cond0 = [s0 e0 c0 p0];

t = [0:0.00000001:0.0005];

[t,y] = ode45(@(t,y)reazioneEnzimatica(y,k1,k_1,k2,k_2),t,cond0);
s = y(:,1);
e = y(:,2);
c = y(:,3);
p = y(:,4);

% Plots

% Andamento concentrazioni 
figure
for i = 1:length(labels)
   subplot(3,2,i)  
   simplePlot(t,y(:,i),strcat(cTitle,labels{i}),tLabel,cLabel);
end

% Velocita' formazione del prodotto 
dp = k2.*c -k_2.*e.*p;
subplot(3,2,5)
simplePlot(t,dp,'Velocita` formazione Prodotto',tLabel,vLabel);

% - Che tipo di crescita del prodotto e decrescita del substrato si osserva?
% Il prodotto segue un andamento di crescita che può essere visto come
% derivata di una composizione di un andamento esponenziale nel primo
% tratto e lineare nel secondo.
% Il substrato mantiene un andamento simile nella forma, ma nella prima
% parte, quella che possiamo chiamare di transizione vi è una
% decrescita lineare

% - Si raggiunge un valore di regime? In quanto tempo? In che rapporto si trovano le diverse concentrazioni?
% Raggiungo un valore di regime per tutti e quattro le entità,ì: prodotto = 100, composto, substrato = 0, enzima = 1

% - A regime in che forma si trova l’enzima?
% A regima l'enzima di trova libero

% - Su che scala temporale posso valutare l’andamento a “regime”? 
% - Quante fasi transitorie posso individuare?

%  - Verificare che tau0 sia il tempo per cui dc/dt sia
%  uguale a 0 (variazione nulla --> c costante)
tau0 = 1 / (k1*s0 + k_1 + k2);
dc =  k1 * e .* s - (k_1+k2) * c + k_2 * p .* e;
idx = find(dc < e-5);
tau0_hat = t(idx(1));
% tau0 e tau0_hat differiscono di un ordine di grandezza
% (molto piccolo). Accettabile(?)

% - Verificare la legge di Michaelis-Menten al variare di s0
Vmax = k2 * e0;
km = (k_1 + k2) / k1;
s0 = linspace(0,200);
V = (Vmax * s0)./(km + s0);
figure()
simplePlot(s0,V,'Michaelis - Menten','substrato[mM]','velocita`[mM/t]');
hold on
plot(km,Vmax/2,'o')
text(km + 5,Vmax/2 + 0.1,'Km,V_{max/2}','FontSize',11)
hold off


%% Set 4

k1  = 75e+3;   % [mM^-1*sec^-1]
k_1 = 75;        % [sec^-1]
k2  = 600e+3;  % [sec^-1]
k_2 = 0;         % [mM^-1*sec^-1]
k3  = 75e+2;   % [mM^-1*sec^-1]
k_3 = 75;        % [sec^-1]
k4  = 600e+3;  % [sec^-1]
k_4 = 0;         % [mM^-1*sec^-1]

s0= 100;
e0 = 1; 
c0 = 0;
p0 = 0;

i0   = 0:100:1000;
c2_0 = 0;
p2_0 = 0;


for i = 1:length(i0)

      cond0 = [s0, e0, c0, p0, i0(i), c2_0, p2_0];
      tsim = 0:0.00001:0.001;
      
      [t,y] = ode45(@(t,y)reazioneEnzimaticaInibitore(y,k1,k_1,k2,k_2,k3,k_3,k4,k_4), tsim ,cond0);

      S{i} = y(:,1);
      C{i} = y(:,2);
      E{i} = y(:,3);
      P{i} = y(:,4);
      I{i} = y(:,5);
      C2{i} = y(:,6);
      P2{i} = y(:,7);
      dP{i} = k2 * C{i} - k_2 * P{i}.*E{i};
      
end

figure()

subplot(4,2,1)
for i = 1:10
   plot(tsim,S{i})
   i0_label{i} = strcat('i0 = ',num2str(i * 100));
   hold on
end
hold off 
legend(i0_label{:})
title('Andamento substrato')
xlabel('[sec]')
ylabel('[mM-1]')

subplot(4,2,2)
for i = 1:10
   plot(tsim,C{i})
   hold on
end
hold off 
legend(i0_label{:})
title('Andamento complesso')
xlabel('[sec]')
ylabel('[mM-1]')

subplot(4,2,3)
for i = 1:10
   plot(tsim,E{i})
   hold on
end
hold off 
legend(i0_label{:})
title('Andamento enzima')
xlabel('[sec]')
ylabel('[mM-1]')

subplot(4,2,4)
for i = 1:10
   plot(tsim,P{i})
   hold on
end
hold off 
legend(i0_label{:})
title('Andamento prodotto')
xlabel('[sec]')
ylabel('[mM-1]')

subplot(4,2,5)
for i = 1:10
   plot(tsim,I{i})
   hold on
end
hold off 
legend(i0_label{:})
title('Andamento enzima 2')
xlabel('[sec]')
ylabel('[mM-1]')

subplot(4,2,5)
for i = 1:10
   plot(tsim,C2{i})
   hold on
end
hold off 
legend(i0_label{:})
title('Andamento complesso 2')
xlabel('[sec]')
ylabel('[mM-1]')

subplot(4,2,6)
for i = 1:10
   plot(tsim,P2{i})
   hold on
end
hold off 
legend(i0_label{:})
title('Andamento prodotto 2')
xlabel('[sec]')
ylabel('[mM-1]')

subplot(4,2,7)
for i = 1:10
   plot(tsim,dP{i})
   hold on
end
hold off 
legend(i0_label{:})
title('Velocita` formazione prodotto')
xlabel('[sec]')
ylabel('[mM-1]')



i0 = 25;
k1 = 100*10^3 : 300*10^3 : 130*10^4;

figure()
for j = 1:length(k1)
    
     [t,x]=ode45(@(t,x)reazioneEnzimaticaInibitore(x,k1(j),k_1,k2,k_2,k3,k_3,k4,k_4), [0 0.001] ,[s0 c0 e0 p0 i0 c2_0 p2_0]); %risolvo il sistema differenziale   
    
    s = x(:,1);
    c = x(:,2);
    e = x(:,3);
    p = x(:,4);
    i = x(:,5);
    c2 = x(:,6);
    p2= x(:,7);    
    
    dp = k2*c - k_2*p.*e;
    
    legenda{j}=(['k1 = ', num2str(k1(j))]);
   
    
    subplot(4,2,1)
    plot(t,s)
    legend(legenda)
    title('Substrato')
    xlabel('[sec]')
    ylabel('[mM-1]')
    grid
hold on
    subplot(4,2,2)
    plot(t,c)
    legend(legenda)
    title('Complesso')
    xlabel('[sec]')
    ylabel('[mM-1]')
    grid
hold on
    subplot(4,2,3)
    plot(t,e)
    legend(legenda)
    title('Enzima')
    xlabel('[sec]')
    ylabel('[mM-1]')
    grid
hold on
    subplot(4,2,4)
    plot(t,p)
    legend(legenda)
    title('Prodotto')
    xlabel('[sec]')
    ylabel('[mM-1]')
    grid
hold on
    subplot(4,2,5)
    plot(t,c2)
    legend(legenda)
    title('Complesso2')
    xlabel('[sec]')
    ylabel('[mM-1]')
    grid
hold on
    subplot(4,2,6)
    plot(t,i)
    legend(legenda)
    title('Inibitore')
    xlabel('[sec]')
    ylabel('[mM-1]')
    grid
hold on
    subplot(4,2,7)
    plot(t,p2)
    legend(legenda)
    title('Prodotto2')
    xlabel('[sec]')
    ylabel('[mM-1]')
    grid
hold on
    subplot(4,2,8)
    plot(t,dp)
    legend(legenda)
    title('Velocità formazione prodotto')
    xlabel('[sec]')
    ylabel('[mM-1]')
    grid
hold on

end













