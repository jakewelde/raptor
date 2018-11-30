% compute the nominal control inputs based on a trajectory

% TODO :
% - polynomial trajectory planner (matrix inversion for starters)
% - euler angles to rotation matrix trajectory
% - TVLQR

%% Numerical Dynamics

% Compute numerical function that represents the evolution of the system's
% state at the current point in state space, but is still a function of the
% inputs of the system

% M x'' = B u + a(x)

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

%% Differential Flatness

% using a minimal representation of the trajectory in a flat-ish space, we
% can plan a kinematic trajectory and then use that trajectory to determine
% the necessary accelerations in state space to follow it 

t = segment_dt*j;
flat_state = compute_derivatives(C,t);

% determine these based on trajectory 
x_dd_des = [flat_state(3);0;0]; % COM acceleration (2nd derivative)
x_ddd_des = [flat_state(4);0;0]; % COM jerk (3rd derivative)
x_dddd_des = [flat_state(5);0;0]; % COM snap (4th derivative)

% x_dd_des = [0;0;0]; % COM acceleration (2nd derivative)
% x_ddd_des = [0;0;0]; % COM jerk (3rd derivative)
% x_dddd_des = [0;0;0]; % COM snap (4th derivative)

% x_dd_des = [0;flat_state(3);0]; % COM acceleration (2nd derivative)
% x_ddd_des = [0;flat_state(4);0]; % COM jerk (3rd derivative)
% x_dddd_des = [0;flat_state(5);0]; % COM snap (4th derivative)
% x_dd_des = [0;0;flat_state(3)]; % COM acceleration (2nd derivative)
% x_ddd_des = [0;0;flat_state(4)]; % COM jerk (3rd derivative)
% x_dddd_des = [0;0;flat_state(5)]; % COM snap (4th derivative)

% 15*cos(t*(4*pi))

if(t < .5)
    w_d_des = [4;7*pi/2; 0]; % Gripper angular acceleration 
else
    w_d_des = [-4;-7*pi/2; 0]; % Gripper angular acceleration 
end
% sol = solve(compute_thrust_constraint_state(eye(3),Om_d_des,zeros(3),zeros(3),zeros(3),zeros(3)) == 0,Om_d_des)
F = compute_F_state(Rg, Rq,x_dd_des);
d = compute_d_state(Rg,Rq,Om,w,w_d_des,x_dd_des,x_ddd_des,x_dddd_des);

REF = rref([F d]);


H = [F d];
[U,S,V] = svd(H);

if(S(end) ~= 0)
%     disp('Warning: encountered inconsistent system in feedback linearization. Using previous control input.');
%     disp('Singular Value Ratio:')
%     disp(singular_vals(1) / singular_vals(end));
    S(end) = 0;
    H = U*S*V.';
end
REF = rref(H);
Om_d_des = double(REF(1:3,4));

% if(REF(end,:) == zeros(1,4))
%     Om_d_des = double(REF(1:3,4));
% else
%     Om_d_des = zeros(3,1);
%     s = svd(G)
%     disp('Error: encountered inconsistent system in differential flatness');
%     disp(s(1) / s(end));
% 
% end

% fictitious input
acc_des = [x_dd_des(3); Om_d_des; w_d_des];

% TODO : Later will need to compute desired velocities and positions based
% on differential flatness in order to do TVLQR

%% Feedback Linearization

% using the accelerations in state space computed by the differential
% flatness planner (creating a fictitious input), we can compute the
% necessary control inputs to achieve those accelerations using the
% dynamics at the current point in state
% space

G = [B (M*acc_des - a)];
[U,S,V] = svd(G);

% if(REF(end,:) ~= zeros(1,7))
if(S(end) ~= 0)
%     disp('Warning: encountered inconsistent system in feedback linearization. Using previous control input.');
%     disp('Singular Value Ratio:')
%     disp(singular_vals(1) / singular_vals(end));
    S(end) = 0;
    G = U*S*V.';
end
REF = rref(G);
u_ff = double(REF(1:6,7));

