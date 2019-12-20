function dq = modelloHill(t,q,Vmax,km,p)
%creo modello

V2=5; %[l]
k01 = 1.2; %[ore?1]
k02 = 1.2; %[ore?1]
% Vmax = 110; %[mg/ore]
% km = 50; %[mg]
k21 =(Vmax*q(1).^(p-1))/(km.^(p) + q(1).^(p));

dq=zeros(2,1); %Perchè ode vuole vettore colonna

dq(1) = -(k01+k21)*q(1);
dq(2) = k21*q(1) - k02*q(2);
end
