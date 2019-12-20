function [vect] = pulseTrain(n_impulsi,tau,t_eliminazione,t,durata_sim,tau_atck,scaling_f)
%PULSETRAIN Crea un treno di impulsi 
%   n_impulsi = numero di impulsi 
%   tau = intervallo temporale tra gli impulsi
%   t_eliminazione = tempo di eliminazione del sistema
%   t = vettore risposta impulsiva del sistema
%   durata_sim = durata simulazione in ore
%   tau_atck = intervallo temporale tra il primo impulso e il secondo
%   scaling_f = fattore di moltiplicativo per dose d'attacco

tau_eliminazione = t_eliminazione/5;

%Da quanti campioni e' formato un tau?
%n elementi di t < tau
n_samples_tau = length(find(t<tau_eliminazione));

if nargin > 5
    vect = [ones(1,1).*scaling_f zeros(1,ceil(tau_atck*n_samples_tau))];
else
    vect = [ones(1,1) zeros(1,round(tau*n_samples_tau))];
end

n_impulsi = n_impulsi-1;

if n_impulsi>1
    for i = 1:n_impulsi
        if i<n_impulsi
            vect = horzcat(vect,[ones(1,1) zeros(1,round(tau*n_samples_tau))]);
        elseif i == n_impulsi
            %se e' l'ultimo impulso non aggiungo gli zeri alla fine
            vect = horzcat(vect,[ones(1,1)]);
        end
    end
end

durata_sim_samples = durata_sim*n_samples_tau;

%aggiusto la lunghezza di vect in funIone della durata della simulazione
if (durata_sim_samples<length(vect))
    vect = vect(1:durata_sim_samples);
elseif (durata_sim_samples>length(vect))    
    rem_samples = durata_sim_samples - length(vect);    
    vect = horzcat(vect,[zeros(1,rem_samples)]);
end


end

