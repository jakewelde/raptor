%% Set Initial Conditions

Rq0 = eye(3); % hover 
Rg0 = axisangle(e2,pi/2); % arm downwards

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


%% Ball trajectory

z0 = [-.5; 0; -9; .5; -.7; 13];

ball_position = zeros(3,n);
ball_velocity = zeros(3,n);

time = 0:segment_dt:total_dt;

t_apex = z0(6)/g_;
z_apex = z0(1:3) + z0(4:6)*t_apex - 1/2*g_*t_apex^2*e3;
z_d_apex = z0(4:6) - g_*t_apex*e3;

for i=1:size(time,2)
   t = time(i);
   ball_position(:,i) = z0(1:3) + z0(4:6)*t - 1/2*g_*t^2*e3;
end


%% Plan Trajectory

% trajectory.x = find_coefficients_intermediate([0;0;0;0],[z_apex(1);z_d_apex(1);0;0],t_apex,total_dt);
% trajectory.y = find_coefficients_intermediate([0;0;0;0],[z_apex(2);z_d_apex(2);0;0],t_apex,total_dt);
% trajectory.z = find_coefficients_intermediate([0;0;0;0],[z_apex(3);z_d_apex(3);0;0],t_apex,total_dt);
trajectory.x = find_coefficients([0;0;0;0],[0;0;0;0],total_dt);
trajectory.y = find_coefficients([0;0;0;0],[0;0;0;0],total_dt);
trajectory.z = find_coefficients([0;0;0;0],[0;0;0;0],total_dt);
% trajectory.x = find_coefficients([0;0;0;0],[z_apex(1);z_d_apex(1);0;0],total_dt);
% trajectory.y = find_coefficients([0;0;0;0],[z_apex(2);z_d_apex(2);0;0],total_dt);
% trajectory.z = find_coefficients([0;0;0;0],[z_apex(3);z_d_apex(3);0;0],total_dt);

trajectory.a = find_coefficients([0;0;0;0],[pi/6;-pi/8;0;0],total_dt); 
trajectory.b = find_coefficients([pi/2;0;0;0],[pi/2;3;0;0],total_dt);
trajectory.g = find_coefficients([0;0;0;0],[pi/8;0;0;0],total_dt);
% trajectory.a = find_coefficients([0;0;0;0],[pi/3;0;0;0],total_dt); 
% trajectory.b = find_coefficients([pi/2;0;0;0],[pi/7;0;0;0],total_dt);
% trajectory.g = find_coefficients([0;0;0;0],[pi/2;0;0;0],total_dt);

stacked = [
    trajectory.x; trajectory.y; trajectory.z;
    trajectory.a; trajectory.b; trajectory.g;
];

figure(3)
clf;
names = {'x','y','z','\alpha','\beta','\gamma'};
for coord=1:6
    derivatives = zeros(5,n);
    for i=1:n
        derivatives(:,i) = compute_derivatives(stacked((coord-1)*8+(1:8)),i*segment_dt,total_dt);
    end
    subplot(2,3,coord);
    plot(segment_dt*(1:n),derivatives.')
    title(names{coord});
end
drawnow;

%% Dynamic Simulation

xe_rec = zeros(6,n); % records the planned trajectory of end effector
xs_rec = zeros(6,n); % records the trajectory computed with diff. flatness

w_rec = zeros(3,n);  % records the planned gripper ang. vel
Om_rec = zeros(3,n); % records the ang. vel. computed with diff. flatness

percent_done = -1;

options = odeset('AbsTol',1e-8,'RelTol',1e-4, 'MaxStep',0.00001);

for j=1:n
    
    % progress indicator
    percent = floor((j / n)*100);
    if(percent > percent_done)
      fprintf('simulating dynamics: %d%% done.\n',percent);
      percent_done = percent;
    end

    t = segment_dt * j;
    
    % compute feedforward control
    [u_ff, xe_des, xs_des, w_des, Om_des] = compute_control(stacked, t, total_dt);

    record_nominal_trajectory;
    
    % integrate dynamics 
    tspan=segment_dt*(j-1)+[0 segment_dt];
    [~,qs] = ode45(@(t,x) ode(x,u_ff),tspan,current_state,options);   
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