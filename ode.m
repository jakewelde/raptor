function xdot = ode(x,u)
% 
% velocity_cascade = zeros(3,22);
% velocity_cascade(1:3,3+9+2+(1:3)) = eye(3);

global mg_ mq_ g_
global e3
global hat
global compute_Rq_state
global compute_M_state compute_B_state compute_a_state

[~, Rg, th1, th2, xs_d, w, th1d, th2d] = state_from_vector(x);

Rq = compute_Rq_state(Rg,th1,th2);

xdot = [
  xs_d; % first three rows
  reshape(Rg*hat(w),[9 1]); % d/dt Rq = Rq Om_hat
  th1d;
  th2d;
  u(1)/(mg_+mq_)*Rq*e3 - g_*e3; % acceleration of system center of mass is due to gravity and thrust
  compute_M_state(Rg,th1,th2) \ ( compute_a_state(Rg,w,th1,th2,th1d,th2d) + compute_B_state(Rg,th1,th2)*u)
];

% M(x) x' = a(x) + B(x) u

end