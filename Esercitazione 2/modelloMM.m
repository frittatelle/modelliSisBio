function dq = modelloMM(t,q,Vmax,km)
%creo modello

V2=5; %[l]
k01 = 1.2; %[ore?1]
k02 = 1.2; %[ore?1]
k21 =Vmax/(km + q(1));

dq=zeros(2,1); %Perchè ode vuole vettore colonna

dq(1) = -(k01+k21)*q(1);
dq(2) = k21*q(1) - k02*q(2);
end
