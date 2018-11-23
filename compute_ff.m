% M x'' = B u + a

%% Compute numerical dynamics at this point in state space
M = [
 (mg_+mq_)   zeros(1,6);
  zeros(6,1) compute_M_state(Rg,Rq);
];
B = [
    e3.'*(1/(mg_+mq_)*Rq*e3) zeros(1,5);
    compute_B_state(Rg,Rq);
];
a = [
    -g_;
    compute_a_state(Rg,Rq,Om,w);
];

%% Determine desired accelerations from differential flatness
x_dd_des3 = 0;
Om_d_des = [0;1.5; 0];
w_d_des = [0;1.5; 0];
acc_des = [x_dd_des3; Om_d_des; w_d_des];

%% Use feedback linearization to compute necessary inputs

REF = rref([B (M*acc_des - a)]);

if(REF(end,:) == zeros(1,7))
    u_ff = double(REF(1:6,7));
else
    u_ff = zeros(6,1);
    disp('Error: encountered inconsistent system while computing feed forward control');
end
