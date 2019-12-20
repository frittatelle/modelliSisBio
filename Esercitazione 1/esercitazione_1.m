%% ESERCITAZIONE 1 - Simulazione di modelli compartimentali

% Obbiettivo: Simulare diversi modelli compartimentali
% lineari plottando l'andamento delle concentrazioni dei
% compartimenti accessibili


% Funzioni utili del Control System Toolbox:
% sist=ss(A,B,C,D) Sistema LTI definito dalle matrici A,B,C,D.
% [A,B,C,D] = ssdata(sist) Matrici A,B,C,D del il sistema LTI sist.
% size(sist) Dimensioni dei vettori di stato, ingresso e uscita.
% p = eig(sist) Autovalori del sistema.
% mu = dcgain(sist) Guadagno statico del sistema.
% [y,t,x] = initial(sist,x0) Movimenti liberi del sistema generati da x(0)=x0.
% [y,t,x] = lsim(sist,u,tu,x0) Movimenti generati dall�ingresso u, definito negli
% istanti tu, a partire dallo stato iniziale x0.
% [y,t,x] = step(sist) Movimenti forzati generati da uno scalino
% di ampiezza unitaria.
% [y,t,x] = impulse(sist) Movimenti forzati generati da un impulso unitario.


clear;


%%
%Modello monocompartimentale - no abs

%provare step
%AUC pochi punti
%provare initial
%calcolare eig
%contributi exp con i fratti semplici

monocomp;

%%
%Modello monocompartimentale - abs
monocomp_abs;

%%
%Modello compartimentale a due compartimenti (caso di studio C-peptide)

%trasformare pmol

bicomp_cpeptide;

%%
%Modello tricompartimentale

%da fare

tricomp;

%%
%Modello monocompartimentale - no abs - dosi ripetute

monocomp_dosi_ripetute;
