function[c_lin,t_lin,q_lin] = modello_lineare(k01,k02,k21,bolo,V2,t)

A=[-(k01+k21), 0; k21 , -k02];
B=[bolo;0];
C=[0, 1/V2];
D=0;

sys=ss(A,B,C,D);

[c_lin,t_lin,q_lin]=impulse(sys,t);