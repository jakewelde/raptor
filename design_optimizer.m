

syms time 

C = sym('C',[1 48]);
z_0 = sym('z_0',[3 1]);
z_d0 = sym('z_d0',[3 1]);
syms J0 J1

r = reshape(C,[8,6]).'*basis.';
% - (z_0 + z_d0*time + 1/2*time^2*g_*e3);
rd = diff(r,time);
rdd = diff(rd,time);
rddd = diff(rdd,time);
rdddd = diff(rddd,time);


% integrand = J0*(r.'*r) + J1*(rd.'*rd);
Q0 = diag([10000 10000 10000 .1 .1 .1]);
Q1 = diag([100 100 100 1 1 1]);
J0 = 1;
J1 = .001;

integrand = J0*(r.'*Q0*r);
% + J1*(rddd.'*Q1*rddd);

syms t_apex
syms delta_t
% cost = simplify(int(integrand,time,t_apex - delta_t,t_apex + delta_t));
cost = simplify(int(integrand,time,0,t_apex));
% compute_cost = matlabFunction(cost);
cost_function = matlabFunction(cost,'Vars',{C,t_apex,delta_t,z_0});

disp('derived cost function.')