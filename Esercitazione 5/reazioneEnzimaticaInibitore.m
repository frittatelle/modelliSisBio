function dx = reazioneEnzimaticaInibitore(x,k1,k_1,k2,k_2,k3,k_3,k4,k_4)

%creo modello

dx=zeros(7,1);

dx(1) = -k1*x(3)*x(1) + k_1*x(2) ;
dx(2) = k1*x(3)*x(1) - (k_1+k2)*x(2) + k_2*x(4)*x(3);
dx(3) = -k1*x(3)*x(1) + (k_1+k2)*x(2) - k_2*x(3)*x(4) - k3*x(3)*x(5) + (k_3+k4)*x(6) - k_4*x(7)*x(3);
dx(4) = k2*x(2) - k_2*x(4)*x(3);
dx(5) = -k3*x(3)*x(5) + k_3*dx(6);
dx(6) = k3*x(3)*x(5) - (k_3+k4)*x(6) + k_2*x(7)*x(3);
dx(7) = k4*x(6) - k_4*x(7)*x(3);

% dx(1)=ds;
% dx(2)=dc;
% dx(3)=de;
% dx(4)=dp;
% dx(5)=di;
% dx(6)=dc2;
% dx(7)=dp2;

end
