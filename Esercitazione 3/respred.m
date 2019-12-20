function [r,y] = respred(p,t,z,mod_exp,sigma_v)
%RESIDUI
% Residui, predizioni con modelli esponenziali

switch mod_exp
    case 1
        y = p(1)*exp(-p(2)*t);
    case 2
        y = p(1)*exp(-p(2)*t) + p(3)*exp(-p(4)*t);
    case 3
        y = p(1)*exp(-p(2)*t) + p(3)*exp(-p(4)*t) + p(5)*exp(-p(6)*t);
end

% Residui
r = z-y;

% Residui pesati 
if nargin>4
    r = (z-y)./sqrt(sigma_v);
end


end

