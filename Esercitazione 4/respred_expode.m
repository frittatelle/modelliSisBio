function [r,y] = respred_expode(p,t,z,mod,dose,sigma_v)
%RESIDUI
% Residui, predizioni con modelli esponenziali

switch mod
  
    % Stima parametri modello ODE
    case 1 
        k01 = p(1);
        k12 = p(2);
        k21 = p(3);
        V = p(4);

        A = [-(k01+k21),k12; k21,-k12];
        B = [dose; 0];
        C = [1/V,0];
        D = 0;
        
        sys = ss(A,B,C,D);
        [y,t_imp] = impulse(sys, 0:t(end));
        y = interp1(t_imp, y, t);

        % Residui
        r = z-y;
        
        % Residui pesati
        if nargin>5
            r = (z-y)./sqrt(sigma_v);
        end
          
        
    % Stima parametri modello esponenziale    
    case 2   
        
        A = p(1);
        alfa = p(2);
        B = p(3);
        beta = p(4);
               
        y = A*exp(-alfa*t) + B*exp(-beta*t);
        
        % Residui
        r = z-y;
        
        % Residui pesati
        if nargin>5
            r = (z-y)./sqrt(sigma_v);
        end      
            
end
    

end

