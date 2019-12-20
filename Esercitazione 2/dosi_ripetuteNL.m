%% Dosi ripetute

clear all
close all

%Creo impulsi ripetute, ogni 4 ore

%% Lineare

k01=1.2;
k02=1.2;
k21=2.2;
bolo=500;
V2=5;

A=[-(k01+k21), 0; k21 , -k02];
B=[bolo;0];
C=[0, 1/V2];
D=0;

sys=ss(A,B,C,D);

tau=4; %Ogni 4 ore do dose ripetuta
Ts=tau/50;
Tf=3*24;

[u,t]=gensig('pulse', 4, Tf+10, Ts);
u((Tf/Ts):end)=0;

[c,t,q]=lsim(sys,u,t);

massimo_valore(1)=max(c);

figure(1)

plot(t,c/Ts)
hold on
%% Dose ripetuta, Michaelis Menten/Hill 

Vmax = 110; %[mg/ore^-1]
km = 50; %[mg]

n_boli=ceil(Tf/tau);

% Mi ricordo che Hill, con q=1 è una M-M 

p_max=3;

for p=1:p_max
    
    ci1_new=0;
    ci2_new=0;
    
    q_rip=[];
    t_rip=[];
    
    for i=1:n_boli 
            
        [t,q]=ode23s(@(t,q)modelloHill(t,q,Vmax,km,p),[tau*(i-1):Ts:(tau*i - Ts)],[bolo+ci1_new ci2_new]);

        q1=q(:,1);
        q2=q(:,2);

        ci1_new = q1(end); %Nuova condizione iniziale 1
        ci2_new = q2(end); %Nuova condizione iniziale 2

        q_rip = vertcat(q_rip, q(:,2));
        t_rip = vertcat(t_rip, t);
    
    end
    
    plot(t_rip,q_rip/V2) %Grafico dosi ripetute
    hold on
    
end

title('Risposta del sistema a dosi ripetute (dose=500mg, tau=4 ore)')
xlabel('tempo[ore]')
ylabel('concentrazione[mg/l]')
legend('lineare', 'q=1', 'q=2', 'q=3')
hold off

% All'aumentare di q la cinetica non lineare di Hill tende
% alla linearita'. In particolare Vmax = Vlin e di
% conseguenza Cmax = Clin