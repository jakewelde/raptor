%% Configure Simulation Parameters
segment_dt = .001;
total_dt = 1;
n = floor(total_dt/segment_dt);

%% Ball trajectory

% z0 = [-2.1; 1.8; -5.3; 2; -2; g_];
% z0 = [-1.5; 3.9; -4.7; 2; -2; g_];

% z0 = [-.5; .9; -4.7; 1; -1; g_];

% z0 = [-.5; .9; -4.7; 3; -3; g_];

v0 = 1;
z0 = [.75-v0; 0; -g_/2; v0; 0; g_]; 


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

trajectory.x = find_coefficients([0;0;0;0],[z_apex(1);z_d_apex(1)/2;0;0],total_dt);
trajectory.y = find_coefficients([0;0;0;0],[z_apex(2);z_d_apex(2)/2;0;0],total_dt);
trajectory.z = find_coefficients([0;0;0;0],[z_apex(3);z_d_apex(3)/2;0;0],total_dt);
trajectory.a = find_coefficients([0;0;0;0],[0;0;0;0],total_dt); % wrist pronation a: + pi, - pi
trajectory.b = find_coefficients([pi/2;0;0;0],[pi/2;0;0;0],total_dt); % swing      b: 0, pi
trajectory.g = find_coefficients([0;0;0;0],[0;0;0;0],total_dt); % yaw               -2*pi,2*pi

% aggressive example
% trajectory.x = find_coefficients([0;0;0;0],[.3;0;0;0],total_dt);
% trajectory.y = find_coefficients([0;0;0;0],[.2;0;0;0],total_dt);
% trajectory.z = find_coefficients([0;0;0;0],[.8;0;0;0],total_dt);
% trajectory.a = find_coefficients([0;0;0;0],[pi/4;0;0;0],total_dt); 
% trajectory.b = find_coefficients([pi/2;0;0;0],[pi/2+pi/4;0;0;0],total_dt);
% trajectory.g = find_coefficients([0;0;0;0],[pi/4;0;0;0],total_dt);

% conservative example
% trajectory.x = find_coefficients([0;0;0;0],[.03;0;0;0],total_dt);
% trajectory.y = find_coefficients([0;0;0;0],[.05;0;0;0],total_dt);
% trajectory.z = find_coefficients([0;0;0;0],[.4;0;0;0],total_dt);
% trajectory.a = find_coefficients([0;0;0;0],[pi/10;0;0;0],total_dt); 
% trajectory.b = find_coefficients([pi/2;0;0;0],[pi/2+pi/4;0;0;0],total_dt);
% trajectory.g = find_coefficients([0;0;0;0],[pi/10;0;0;0],total_dt);



stacked = [
    trajectory.x; trajectory.y; trajectory.z;
    trajectory.a; trajectory.b; trajectory.g;
];

[stacked, minval, retcode] = trajectory_optimization(z0,stacked);
minval
retcode

show_trajectory;
pause


% 
% figure(3)
% clf;
% names = {'x','y','z','\alpha','\beta','\gamma'};
% for coord=1:6
%     derivatives = zeros(5,n);
%     for i=1:n
%         derivatives(:,i) = compute_derivatives(stacked((coord-1)*8+(1:8)),i*segment_dt,total_dt);
%     end
%     subplot(2,3,coord);
%     plot(segment_dt*(1:n),derivatives.')
%     title(names{coord});
% end
% drawnow;

%% Tracking

[~,current_state] = compute_control(stacked,0,total_dt);

state = zeros(n,size(current_state,1));
state_des = zeros(size(state));

state(1,:) = current_state;

us = zeros(n,6);
[xs, Rg, th1, th2, xs_d, w, th1d, th2d] = state_from_vector(current_state);


%% Dynamic Simulation

xe_rec = zeros(6,n); % records the planned trajectory of end effector
xs_rec = zeros(6,n); % records the trajectory computed with diff. flatness

w_rec = zeros(3,n);  % records the planned gripper ang. vel
% Om_rec = zeros(3,n); % records the ang. vel. computed with diff. flatness

percent_done = -1;

% options = odeset('AbsTol',1e-8,'RelTol',1e-4, 'MaxStep',0.00001);

for j=1:n
    
    % progress indicator
    percent = floor((j / n)*100);
    if(percent > percent_done)
      fprintf('simulating dynamics: %d%% done.\n',percent);
      percent_done = percent;
    end

    t = segment_dt * j;
    
    % compute feedforward control
    [u_ff, current_state_des] = compute_control(stacked, t, total_dt,current_state);
    
    % integrate dynamics 
    tspan=segment_dt*(j-1)+[0 segment_dt];
    [~,qs] = ode45(@(t,x) ode(x,u_ff),tspan,current_state);   
    current_state = qs(end,:)';
    
    % reorthonormalize rotation matrices (project back onto manifold)
    [xs, Rg, th1, th2, xs_d, w, th1d, th2d]  = state_from_vector(current_state);
    [U, ~, V] = svd(Rg);
    Rg = U * V';

    % record state and control inputs for plotting
    state(j,:) = vector_from_state(xs, Rg, th1, th2, xs_d, w, th1d, th2d);
    state_des(j,:) = current_state_des;
    us(j,:) = u_ff.';
 
end

%% Visualization

plot_sim_results