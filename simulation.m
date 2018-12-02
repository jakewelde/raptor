%% Initial state defined here
set_initial_conditions;

%% Configure simulation parameters
segment_dt = .0001;
total_dt = 1;
n = floor(total_dt/segment_dt);
state = zeros(n,size(x0,1));
state(1,:) = x0;
current_state = x0;

us = zeros(n,6);
[xs, Rq, Rg, xs_d, Om, w] = state_from_vector(x0);

%% Plan trajectory

planning_matrix;

% init_state = [pi/8;0;0;0];
% final_state = [pi-pi/8;0;0;0];
% init_state = [0;3;0;0];
% final_state = [4;5;0;0];

init_state = [0;0;0;0];
final_state = [.5;0;0;0];


C = D*[init_state;final_state];

% figure(3);

derivatives = zeros(5,n);
for j=1:n
   derivatives(:,j) = compute_derivatives(C,j*segment_dt);
end


plot((1:n)*segment_dt,derivatives(1,:));

shg;

init_angle = [pi/2;0;0;0];
% final_angle = [pi/2;0;0;0];
final_angle = [7*pi/8;0;0;0];
C_R = D*[init_angle;final_angle];

% init_wrist = [0;0;0;0];
% final_wrist = [pi/2;0;0;0];
% C_W = D*[init_wrist;final_wrist];

Rg_des = axisangle(e2,basis*C_R);
%  *axisangle(e1,basis*C_W);

w_des = unhat(Rg_des.' * diff(Rg_des,time));
w_d_des = diff(w_des,time);
w_dd_des = diff(w_d_des,time);
w_ddd_des = diff(w_dd_des,time);

desired_orientation = matlabFunction(Rg_des);
angular_derivatives = matlabFunction([
   w_des w_d_des w_dd_des w_ddd_des
].');


% ang_vel = zeros(3,n);
% for j=1:10:n
%    deriv = angular_derivatives(j*segment_dt);
%    ang_vel(:,j) = deriv(1,:);
% end

% clf;
% plot((1:n)*segment_dt,ang_vel);
% legend('Om1','Om2','Om3');

% decim = 10;
% 
% Rg_x_des = zeros(3,n/decim);
% Rg_y_des = zeros(3,n/decim);
% Rg_z_des = zeros(3,n/decim);
% for j=1:n/decim
%     Rg_des_t = desired_orientation(j*segment_dt*decim);
%     Rg_x_des(:,j) = Rg_des_t(:,1);
%     Rg_y_des(:,j) = Rg_des_t(:,2);
%     Rg_z_des(:,j) = Rg_des_t(:,3);
% end

% zo = zeros(1,n/decim);
% figure(4);
% clf;
% hold on;
% quiver3(zo,zo,zo,Rg_x_des(1,:),Rg_x_des(2,:),Rg_x_des(3,:),'color',[1 0 0]);
% quiver3(zo,zo,zo,Rg_y_des(1,:),Rg_y_des(2,:),Rg_y_des(3,:),'color',[0 1 0]);
% quiver3(zo,zo,zo,Rg_z_des(1,:),Rg_z_des(2,:),Rg_z_des(3,:),'color',[0 0 1]);
% xlabel('x'); ylabel('y'); zlabel('z');
% hold off;
% pause

%% Dynamic Simulation

u_ff = zeros(6,1); % TODO : eliminate

xs_rec = zeros(15,n);
xe_rec = zeros(3,n);
xq_rec = zeros(3,n);

for j=1:n
    
    % compute feedforward control
    compute_ff;

    % integrate dynamics 
    tspan=segment_dt*(j-1)+[0 segment_dt];
    [~,qs] = ode113(@(t,x) ode(x,u_ff),tspan,current_state);   
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