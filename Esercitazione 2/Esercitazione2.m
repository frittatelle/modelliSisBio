%% ESERCITAZIONE 2, compartimento non lineare con assorbimento

% Obbiettivo: Simulare un modello a compartimenti con
% cinetica non lineare. Valutare poi le dosi ripetute

 clear all
 close all
 clc
 
V2=5;                             %[l]
k01 = 1.2;                        %[ore^-1]
k02 = 1.2;                        %[ore^-1]
Vmax = 110;                       %[mg/ore^-1]
km = 50;                          %[mg]
bolo = 500;                       %[mg]

% Bolo iniziale, somministrato a t=0 è uguale a imporre q1(0)=bolo, come
% condizione iniziale

[t,q]=ode45(@(t,q)modelloMM(t,q,Vmax,km),[0:0.02:8],[bolo 0]); %risolvo il sistema differenziale

for i=1:length(q) %Calcolo i valori di k21 e di F21 al variare della quantità nel compartimento 1
    k21(i) =Vmax/(km + q(i,1));
    F21(i)=k21(i)*q(i,1);
end

c=q(:,2)/V2;

figure

subplot(1,2,1)
plot(t,q)
title('Andamento delle quantità con Michaelis-Menten')
xlabel('Tempo [ore]')
ylabel('Quantità [mg]')
legend('comp 1','comp 2')
grid

subplot(1,2,2)
plot(t,c)
title('Andamento C=q2/V2')
xlabel('Tempo [ore]')
ylabel('Concentrazione [mg/l]')
grid

%Mi aspetto che per t>t* tale che q1(t*)=0, il valore di k21 'saturerà' a
%Vmax/km=2.2

figure

subplot(3,1,1)
plot(t,k21)
title('Andamento della velocità di assorbimento')
xlabel('Tempo [ore]')
ylabel('[ore^-1]')
grid

%Mi aspetto F21(q1=0)=0 e per q1(t) che tende all'infinito, F21 dovrebbe
%tenedere a Vmax, verifico dal grafico

subplot(3,1,2)
plot(q(:,1),F21)
title('Andamento del flusso d ingresso nel compartimento 2 : F(q1)')
xlabel('Quantità comp1 [mg]')
ylabel('Flusso [mg/ore]')
grid

subplot(3,1,3)
plot(t,F21)
title('Andamento del flusso d ingresso nel compartimento 2 : F(t)')
xlabel('t[ore]')
ylabel('Flusso [mg/ore]')
grid

%% a) Confronto con modello lineare, k21_lin=2.2

[c_lin,t_lin,q_lin]=modello_lineare(k01,k02,2.2,bolo,V2,t); 

% Lineare, equazione = (bolo/V2)*e^(-3.4*t)
% Non lineare, equazione = (bolo/V2)*e^-(1.2+k21), non si vede poichè
% max(k21) = 2.2
% k21 Michaelis Menten, allora per valori di Q1(t) (NL) > 110/k21_lin -50
% Il sistema non lineare avrà un assorbimento più lento

figure

subplot(3,1,1)
plot(t_lin,q_lin(:,1),t,q(:,1))
title('Andamento delle quantità del compartimento 1')
legend('Modello lineare','Modello non lineare')
xlabel('Tempo [ore]')
ylabel('Quantità [mg]')
grid

subplot(3,1,2)
plot(t_lin,q_lin(:,2),t,q(:,2))
title('Andamento delle quantità del compartimento 2')
legend('Modello lineare','Modello non lineare')
xlabel('Tempo [ore]')
ylabel('Quantità [mg]')
grid

subplot(3,1,3)
plot(t_lin,c_lin, t, c)
title('Andamento della concentrazione del compartimento 2')
legend('Modello lineare','Modello non lineare')
xlabel('Tempo [ore]')
ylabel('Concentrazione [mg]')
grid

%% b) Cosa varia al variare di Vmax e km

Vmax2=Vmax*2;
km2=km/2;

[t_vmax,q_vmax]=ode45(@(t,q)modelloMM(t,q,Vmax2,km),[0 8],[bolo 0]); %modello con bolo iniziale è uguale a definire una condizione iniziale della quantità di q1=bolo
[t_km,q_km]=ode45(@(t,q)modelloMM(t,q,Vmax,km2),[0 8],[bolo 0]); %modello con bolo iniziale è uguale a definire una condizione iniziale della quantità di q1=bolo

for i=1:length(q_vmax)
k21_vmax(i) =Vmax2/(km + q_vmax(i,1));
F21_vmax(i)=k21_vmax(i)*q_vmax(i,1);
end

for i=1:length(q_km)
k21_km(i) =Vmax/(km2 + q_km(i,1));
F21_km(i)=k21_km(i)*q_km(i,1);
end

figure

subplot(2,2,1)
plot(t,q(:,1), t_vmax,q_vmax(:,1), t_km, q_km(:,1))
title('Andamento delle quantità del compartimento 1')
legend('Vmax=110, km=50','Vmax=2*Vmax, km=50', 'Vmax=110, km=km/2')
xlabel('Tempo [ore]')
ylabel('Quantità [mg]')
grid

subplot(2,2,2)
plot( t,q(:,2),t_vmax,q_vmax(:,2), t_km, q_km(:,2))
title('Andamento delle quantità del compartimento 2')
legend('Vmax=110, km=50','Vmax=2*Vmax, km=50', 'Vmax=110, km=km/2')
xlabel('Tempo [ore]')
ylabel('Quantità [mg]')
grid

subplot(2,2,3)
plot(t,k21,t_vmax,k21_vmax, t_km, k21_km)
title('Andamento della velocità di assorbimento')
legend('Vmax=110, km=50','Vmax=2*Vmax, km=50', 'Vmax=110, km=km/2')
xlabel('Tempo [ore]')
ylabel('[ore^-1]')
grid

subplot(2,2,4)
plot(q(:,1),F21,q_vmax(:,1),F21_vmax, q_km(:,1), F21_km)
title('Andamento del flusso d ingresso nel compartimento 2')
legend('Vmax=110, km=50','Vmax=2*Vmax, km=50', 'Vmax=110, km=km/2')
xlabel('Quantità [mg]')
ylabel('Flusso [mg/ore]')
grid

clear q k21 F21

%% c) Cosa cambia riducendo a 50mg la dose

bolo=50;

[t,q]=ode45(@(t,q)modelloMM(t,q,Vmax,km),t,[bolo 0]);
[c_lin,t_lin,q_lin]=modello_lineare(k01,k02,2.2,bolo,V2,t); %2.2 valore di k21_lin

for i=1:length(q)
k21(i) =Vmax/(km + q(i,1));
F21(i)=k21(i)*q(i,1);
end

c=q(:,2)/V2;

figure

subplot(2,2,1)
plot(t,q)
title('Andamento delle quantità con Michaelis-Menten (dose=50mg)')
xlabel('Tempo [ore]')
ylabel('Quantità [mg]')
grid

subplot(2,2,2)
plot(t,c)
title('Andamento C=q2/V2 (dose=50mg)')
xlabel('Tempo [ore]')
ylabel('Concentrazione [mg/l]')
grid

subplot(2,2,3)
plot(t,k21)
title('Andamento della velocità di assorbimento')
xlabel('Tempo [ore]')
ylabel('[ore^-1]')
grid

subplot(2,2,4)
plot(q(:,1),F21)
title('Andamento del flusso d ingresso nel compartimento 2')
xlabel('Quantità [mg]')
ylabel('Flusso [mg/ore]')
grid

figure

subplot(3,1,1)
plot(t_lin,q_lin(:,1),t,q(:,1))
title('Andamento delle quantità del compartimento 1')
legend('Modello lineare','Modello non lineare')
xlabel('Tempo [ore]')
ylabel('Quantità [mg]')
grid

subplot(3,1,2)
plot(t_lin,q_lin(:,2),t,q(:,2))
title('Andamento delle quantità del compartimento 2')
legend('Modello lineare','Modello non lineare')
xlabel('Tempo [ore]')
ylabel('Quantità [mg]')
grid

subplot(3,1,3)
plot(t_lin,c_lin, t, c)
title('Andamento della concentrazione del compartimento 2')
legend('Modello lineare','Modello non lineare')
xlabel('Tempo [ore]')
ylabel('Concentrazione [mg]')
grid

%% d) Considerare ora una non linearità di Hill, con q che può assumere diversi valori

clear t q k21 F21

bolo=500;

p_max=5;


for p=1:p_max

    clear k21 F21
    
    [t,q]=ode45(@(t,q)modelloHill(t,q,Vmax,km,p),[0 8],[bolo 0]);
    
    figure(7)
    
    legenda{p}=(['q = ', num2str(p)]);
    
    subplot(3,1,1)
    plot(t,q(:,1))
    title('Andamento delle quantità, comp 1, con non linearità di Hill')
    legend(legenda)
    xlabel('Tempo [ore]')
    ylabel('Quantità [mg]')
    grid
    
 hold on
 
    subplot(3,1,2)
    plot(t,q(:,2))
    title('Andamento delle quantità, comp 2, con non linearità di Hill')
    legend(legenda)
    xlabel('Tempo [ore]')
    ylabel('Quantità [mg]')
    grid
    
 hold on
 
    c=q(:,2)/V2;
    
    subplot(3,1,3)
    plot(t,c)
    title('Andamento C=q2/V2')
    legend(legenda)
    xlabel('Tempo [ore]')
    ylabel('Concentrazione [mg/l]')
    grid

    hold on
    
    for i=1:length(q)
        k21(i) =Vmax*(q(i,1).^(p-1))/(km.^(p) + q(i,1).^(p));
        biodisp_hill(p,i) = (k21(i)./(k21(i) + k01));
        F21(i)=k21(i)*q(i,1); % Per NL di Hill, tende a 0, per NL di M-M satura
    end

    
    figure(8)
    
    subplot(2,1,1)
    plot(t,k21)
    title('Andamento della velocità di assorbimento')
    legend(legenda)
    xlabel('Tempo [ore]')
    ylabel('[ore^-1]')
    grid
    
    %Velocità di assorbimento cambia drasticamente il suo grafico passando
    %da q=1 a q=2, siccome, per Q1=0, nel caso di q=1 avrò 

hold on
    subplot(2,1,2)
    plot(q(:,1),F21) %Si intersecano tutti nel punto in cui Q=km
    title('Andamento del flusso d ingresso nel compartimento 2')
    legend(legenda)
    xlabel('Quantità [mg]')
    ylabel('Flusso [mg/ore]')
    grid
    
    hold on

end

hold off
hold off

%Pensare se grafici esplicativi, nel caso organizzarli in modo diverso


