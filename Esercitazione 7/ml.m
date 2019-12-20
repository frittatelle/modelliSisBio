function costFun = ml(G,c,gamma,P,sigma_v)

u_ml = (inv(G'*pinv(sigma_v)*G + gamma*P'*P))*G'*inv(sigma_v)*c;
y_ml = G * u_ml;
res_ml = c - y_ml;

n = length(c);
H = G*(G'*(sigma_v)^-1*G + gamma*P'*P)^-1*G'*(sigma_v)^-1;
q = trace(H);

costFun = abs(gamma - (q*res_ml'*(sigma_v)^-1*res_ml)/((n-q)*u_ml'*(P'*P)*u_ml));

end
