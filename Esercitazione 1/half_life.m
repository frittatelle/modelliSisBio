function [t_emivita,y_emivita,ind_emivita] = half_life(y,t)

    %Concentrazione iniziale
    y_start = y(1);
    %Concentrazione al t di emivita
    y_emivita = y_start.*0.5;
    
    %Indice t emivita
    
    ind_emivita = find(abs(y-y_emivita) < 0.2);
    
    t_emivita = t(ind_emivita);
    
end

