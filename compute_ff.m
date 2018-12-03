% compute the nominal control inputs based on a trajectory

%% Evaluate Planned Trajectory in Flat Space

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

flat_state = zeros(5,3);

flat_state(:,1) = compute_derivatives(trajectory.x,t,total_dt);
flat_state(:,2) = compute_derivatives(trajectory.y,t,total_dt);
flat_state(:,3) = compute_derivatives(trajectory.z,t,total_dt);

     xe_des = flat_state(1,:).'; % end effector position     (0th derivative)
   xe_d_des = flat_state(2,:).'; % end effector velocity     (1st derivative)
  xe_dd_des = flat_state(3,:).'; % end effector acceleration (2nd derivative)
 xe_ddd_des = flat_state(4,:).'; % end effector jerk         (3rd derivative)
xe_dddd_des = flat_state(5,:).'; % end effector snap         (4th derivative)

wh = hat(w_des);
wdh = hat(w_d_des);
wddh = hat(w_dd_des);
wdddh = hat(w_ddd_des);

%% Differential Flatness

% using the trajectory in flat space, apply the diffeomorphism to compute
% the trajectory in state space

% find system center of mass trajectory
xs_des = xe_des-Ls_*Rg_des*e1;
xs_d_des = xe_d_des-Ls_*Rg_des*wh*e1;
xs_dd_des = xe_dd_des-Ls_*Rg_des*(wdh+wh^2)*e1;
xs_ddd_des = xe_ddd_des-Ls_*Rg_des*(wddh + 3*wh*wdh + wh^3)*e1;
xs_dddd_des = xe_dddd_des-Ls_*Rg_des*(wdddh + 4*wh*wddh + 6*wh^2*wdh + 3*wdh^2 + wh^4)*e1;

% record spatial trajectories for visualization
xe_rec(1:3,j) = xe_des;
xe_rec(4:6,j) = xe_d_des;
xs_rec(1:3,j) = xs_des;
xs_rec(4:6,j) = xs_d_des;

% record angular trajectories for visualization
w_rec(:,j) = w_des;

% find angular acceleration of quadrotor body
F = compute_F_state(Rg, Rq,xs_dd_des);
d = compute_d_state(Rg_des,Rq,Om,w_des,w_d_des,xs_dd_des,xs_ddd_des,xs_dddd_des);

% necessary due to numerical issues
H = [F d];
[U,S,V] = svd(H);
if(S(end) ~= 0)
    S(end) = 0;
    H = U*S*V.';
end
REF = rref(H);
Om_d_des = double(REF(1:3,4));


%% Feedback Linearization

acc_des = [xs_dd_des(3); Om_d_des; w_d_des];

% Compute numerical function that represents the evolution of the system's
% state at the current point in state space, but is still a function of the
% inputs of the system
M = [
  (mg_+mq_)   zeros(1,6);
  zeros(6,1) compute_M_state(Rg,Rq);
];
B = [
    e3.'*(Rq*e3) zeros(1,5);
    compute_B_state(Rg,Rq);
];
a = [
    -(mg_+mq_)*g_;
    compute_a_state(Rg,Rq,Om,w);
];

% M x'' = B u + a(x)

% TODO : replace all this with the nominal computation from desired state,
% not the current state

G = [B (M*acc_des - a)];
[U,S,V] = svd(G);

% necessary due to numerical issues
if(S(end) ~= 0)
    S(end) = 0;
    G = U*S*V.';
end
REF = rref(G);

u_ff = double(REF(1:6,7));