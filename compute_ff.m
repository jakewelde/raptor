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

if(nargin(angular_derivatives) > 0)
    angular = angular_derivatives(t);
else
    angular = angular_derivatives();
end
if(nargin(desired_orientation) > 0)
    Rg_des = desired_orientation(t); 
else
    Rg_des = desired_orientation();
end


w_des     = angular(1,:).';
w_d_des   = angular(2,:).';
w_dd_des  = angular(3,:).';
w_ddd_des = angular(4,:).';

flat_state = compute_derivatives(C,t);

% determine these based on trajectory 
% xe_des = [flat_state(1);0;0];
% xe_d_des = [flat_state(2);0;0];
% xe_dd_des = [flat_state(3);0;0]; % COM acceleration (2nd derivative)
% xe_ddd_des = [flat_state(4);0;0]; % COM jerk (3rd derivative)
% xe_dddd_des = [flat_state(5);0;0]; % COM snap (4th derivative)

xe_des      = [0;flat_state(1);0];
xe_d_des    = [0;flat_state(2);0];
xe_dd_des   = [0;flat_state(3);0]; % COM acceleration (2nd derivative)
xe_ddd_des  = [0;flat_state(4);0]; % COM jerk (3rd derivative)
xe_dddd_des = [0;flat_state(5);0]; % COM snap (4th derivative)


% xs_dd_des = xe_dd_des-Ls_*Rg_des*(hat(w_d_des)+hat(w_des)^2)*e1;
% % xs_dd_rec(:,j) = xs_dd_des;
% 
% xs_ddd_des = xe_ddd_des-Ls_*Rg_des*(hat(w_dd_des) + hat(w_d_des)*hat(w_des)+2*hat(w_des)*hat(w_d_des)+hat(w_des)^3)*e1;
% % 
% xs_dddd_des = xe_dd_des-Ls_*Rg_des*(hat(w_ddd_des) + hat(w_dd_des)*hat(w_des) + 3*hat(w_d_des)^2 + 2*hat(w_des)*hat(w_dd_des) + 3*hat(w_des)*hat(w_d_des) + hat(w_des)*hat(w_dd_des) + hat(w_des)*hat(w_d_des)*hat(w_des) + 2 * hat(w_des)^2*hat(w_d_des) + hat(w_des)^4)*e1;

% xs_dd_des



% xs_dd_des = xe_dd_des;
% xs_ddd_des = xe_ddd_des;
% xs_dddd_des = xe_dddd_des;

wh = hat(w_des);
wdh = hat(w_d_des);
wddh = hat(w_dd_des);
wdddh = hat(w_ddd_des);

xs_des = xe_des-Ls_*Rg_des*e1;
xs_d_des = xe_d_des-Ls_*Rg_des*wh*e1;
xs_dd_des = xe_dd_des-Ls_*Rg_des*(wdh+wh^2)*e1;
xs_ddd_des = xe_ddd_des-Ls_*Rg_des*(wddh + 3*wh*wdh + wh^3)*e1;
xs_dddd_des = xe_dddd_des-Ls_*Rg_des*(wdddh + 4*wh*wddh + 6*wh^2*wdh + 3*wdh^2 + wh^4)*e1;


xq_des = xe_des-Le_*Rg_des*e1;
xq_rec(:,j) = xq_des;
xe_rec(:,j) = xe_des;


xs_rec(1:3,j) = xs_des;
xs_rec(4:6,j) = xs_d_des;
xs_rec(7:9,j) = xs_dd_des;
xs_rec(10:12,j) = xs_ddd_des;
xs_rec(13:15,j) = xs_dddd_des;

%% 

F = compute_F_state(Rg, Rq,xs_dd_des);
d = compute_d_state(Rg_des,Rq,Om,w_des,w_d_des,xs_dd_des,xs_ddd_des,xs_dddd_des);

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

% fictitious input
acc_des = [xs_dd_des(3); Om_d_des; w_d_des];

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
    S(end) = 0;
    G = U*S*V.';
end
REF = rref(G);
u_ff = double(REF(1:6,7));

