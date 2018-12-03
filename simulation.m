%% Set Initial Conditions

Rq0 = eye(3); % hover 
Rg0 = axisangle(e2,pi/2); % arm downwards

x0 = vector_from_state(...
    [0;0;Ls_],Rq0,Rg0,...
    [0;0;0],[0;0;0],Rg0.'*[0;0;0]...
);

%% Configure Simulation Parameters
segment_dt = .0001;
total_dt = .5;
n = floor(total_dt/segment_dt);
state = zeros(n,size(x0,1));
state(1,:) = x0;
current_state = x0;

us = zeros(n,6);
[xs, Rq, Rg, xs_d, Om, w] = state_from_vector(x0);

%% Plan Trajectory

trajectory_planning;

trajectory.x = find_coefficients([0;0;0;0],[0;0;0;0],total_dt);
trajectory.y = find_coefficients([0;0;0;0],[0;0;0;0],total_dt);
trajectory.z = find_coefficients([0;0;0;0],[0;0;0;0],total_dt);
% trajectory.x = find_coefficients([0;0;0;0],[.1;0;0;0],total_dt);
% trajectory.y = find_coefficients([0;0;0;0],[.05;0;0;0],total_dt);
% trajectory.z = find_coefficients([0;0;0;0],[-.02;0;0;0],total_dt);

C_R = find_coefficients([0;0;0;0],[pi/3;0;0;0],total_dt);
C_S = find_coefficients([pi/2;0;0;0],[.9*pi/2;0;0;0],total_dt);
C_W = find_coefficients([0;0;0;0],[pi/10;0;0;0],total_dt);

roll_derivatives = compute_derivatives(C_R,time,total_dt);
swing_derivatives = compute_derivatives(C_S,time,total_dt);
wrist_derivatives = compute_derivatives(C_W,time,total_dt);
Rg_des = axisangle(e2,roll_derivatives(1))*axisangle(e2,swing_derivatives(1))*axisangle(e1,wrist_derivatives(1));

w_des = unhat(Rg_des.' * diff(Rg_des,time));
w_d_des = diff(w_des,time);
w_dd_des = diff(w_d_des,time);
w_ddd_des = diff(w_dd_des,time);

desired_orientation = matlabFunction(Rg_des);
angular_derivatives = matlabFunction([
       w_des.';
     w_d_des.';
    w_dd_des.';
   w_ddd_des.';
]);

%% Dynamic Simulation

xe_rec = zeros(6,n); % records the planned trajectory of end effector
xs_rec = zeros(6,n); % records the trajectory computed with diff. flatness

w_rec = zeros(3,n);  % records the planned gripper ang. vel
Om_rec = zeros(3,n); % records the ang. vel. computed with diff. flatness

percent_done = -1;

for j=1:n
    
    % progress indicator
    percent = floor((j / n)*100);
    if(percent > percent_done)
      fprintf('simulating dynamics: %d%% done.\n',percent);
      percent_done = percent;
    end

    % compute feedforward control
    compute_ff;

    % integrate dynamics 
    tspan=segment_dt*(j-1)+[0 segment_dt];
    [~,qs] = ode45(@(t,x) ode(x,u_ff),tspan,current_state);   
    current_state = qs(end,:)';
    
    % reorthonormalize rotation matrices (project back onto manifold)
    [xs, Rq, Rg, xs_d, Om, w] = state_from_vector(current_state);
    [U, ~, V] = svd(Rq);
    Rq = U * V';
    [U, ~, V] = svd(Rg);
    Rg = U * V';

    % record state and control inputs for plotting
    state(j,:) = vector_from_state(xs, Rq, Rg, xs_d, Om, w);
    us(j,:) = u_ff.';
 
end

%% Visualization

plot_sim_results