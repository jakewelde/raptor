global cost_function

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

integrand = J0*(rdddd.'*Q0*rdddd) + J1*(rddd.'*Q1*rddd);

syms t_apex
syms delta_t
% cost = simplify(int(integrand,time,t_apex - delta_t,t_apex + delta_t));
cost = simplify(int(integrand,time,0,t_apex));
compute_cost = matlabFunction(cost);

% cost_function = @(C,J0,J1,t_apex,delta_t,z0) compute_cost(C(1),C(2),C(3),C(4),C(5),C(6),C(7),C(8),C(9),C(10),C(11),C(12),C(13),C(14),C(15),C(16),C(17),C(18),C(19),C(20),C(21),C(22),C(23),C(24),J0,J1,delta_t,t_apex,z0(1),z0(2),z0(3),z0(4),z0(5),z0(6));
% cost_function = @(C,J0,J1,t_apex,delta_t,z0) compute_cost(C(4),C(5),C(6),C(7),C(8),C(12),C(13),C(14),C(15),C(16),C(20),C(21),C(22),C(23),C(24),J0,J1,t_apex);
% cost_function = @(C,t_apex,z0) compute_cost(C(2),C(3),C(4),C(5),C(6),C(7),C(8),C(10),C(11),C(12),C(13),C(14),C(15),C(16),C(18),C(19),C(20),C(21),C(22),C(23),C(24),C(26),C(27),C(28),C(29),C(30),C(31),C(32),C(34),C(35),C(36),C(37),C(38),C(39),C(40),C(42),C(43),C(44),C(45),C(46),C(47),C(48),t_apex);
cost_function = @(C,t_apex,z0) compute_cost(...
    C(4),C(5),C(6),C(7),C(8),...
    C(12),C(13),C(14),C(15),C(16),...
    C(20),C(21),C(22),C(23),C(24),...
    C(28),C(29),C(30),C(31),C(32),...
    C(36),C(37),C(38),C(39),C(40),...
    C(44),C(45),C(46),C(47),C(48),...
 t_apex);

disp('derived cost function.')