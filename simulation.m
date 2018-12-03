%% Set Initial Conditions

Rg0 = axisangle(e2,pi/2);

x0 = vector_from_state(...
    [0;0;Ls_],Rq0,Rg0,...
    [0;0;0],[0;0;0],Rg0.'*[0;0;0]...
);

%% Configure Simulation Parameters
segment_dt = .0001;
total_dt = 1;
n = floor(total_dt/segment_dt);
state = zeros(n,size(x0,1));
state(1,:) = x0;
current_state = x0;

us = zeros(n,6);
[xs, Rq, Rg, xs_d, Om, w] = state_from_vector(x0);

%% Plan Trajectory

planning_matrix;

init_state = [0;0;0;0];
final_state = [.5;0;0;0];

C = D*[init_state;final_state];

init_angle = [pi/2;0;0;0];
final_angle = [7*pi/8;0;0;0];
C_R = D*[init_angle;final_angle];

init_wrist = [0;0;0;0];
final_wrist = [pi/6;0;0;0];
C_W = D*[init_wrist;final_wrist];

Rg_des = axisangle(e2,basis*C_R)*axisangle(e1,basis*C_W);

w_des = unhat(Rg_des.' * diff(Rg_des,time));
w_d_des = diff(w_des,time);
w_dd_des = diff(w_d_des,time);
w_ddd_des = diff(w_dd_des,time);

desired_orientation = matlabFunction(Rg_des);
angular_derivatives = matlabFunction([
   w_des w_d_des w_dd_des w_ddd_des
].');

%% Dynamic Simulation

xs_rec = zeros(15,n);
xe_rec = zeros(3,n);
xq_rec = zeros(3,n);

for j=1:n
    
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

%     compute_ff;
 
end

%% Plotting

plot_sim_results