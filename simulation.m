%% Initial state defined here
set_initial_conditions;

%% Configure simulation parameters
segment_dt = .001;
total_dt = 1;
n = floor(total_dt/segment_dt);
state = zeros(n,size(x0,1));
state(1,:) = x0;
current_state = x0;

us = zeros(n,6);
[xs, Rq, Rg, xs_d, Om, w] = state_from_vector(x0);

%% Plan trajectory

planning_matrix;

init_state = [0;3;0;0];
final_state = [4;5;0;0];

C = D*[init_state;final_state];

% figure(3);

% derivatives = zeros(5,n);
% for j=1:n
%    derivatives(:,j) = compute_derivatives(C,j*segment_dt);
% end

% plot((1:n)*segment_dt,derivatives)

%% Dynamic Simulation

u_ff = zeros(6,1); % TODO : eliminate

for j=1:n
    
    % compute feedforward control
    
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

    compute_ff;
%     u_ff = zeros(6,1);
 
end

%% Plotting

plot_sim_results