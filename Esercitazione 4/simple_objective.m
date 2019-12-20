function y = simple_objective(x,data)

t = data.x;
z = data.z;
dose = data.dose;

[r,~] = respred_expode(x,t,z,2,dose);

y = sum(r.^2);


end

