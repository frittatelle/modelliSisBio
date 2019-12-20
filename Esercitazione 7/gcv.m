function costFun = gcv(G,c,gamma,P,sigma_v)

u_gcv = (inv(G'*inv(sigma_v)*G + gamma*P'*P))*G'*inv(sigma_v)*c;
y_gcv = G*u_gcv;
res_cgv = c-y_gcv;

n = length(c);
H = G*(G'*(sigma_v)^-1*G + gamma*P'*P)^-1*G'*(sigma_v)^-1;
q = trace(H);

costFun = (n*res_cgv'*(sigma_v)^-1*res_cgv)/((n-q)^2);

end

