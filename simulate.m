Rq0 = eye(3);

% T1 direction
% Rg0 = axisangle(e1,pi/6)*axisangle(e2,pi/2);
% T2 direction
Rg0 = axisangle(e2,pi/4);
% Mix
Rg0 = axisangle(e1,pi/4)*axisangle(e2,pi/4);
% Down
% Rg0 = axisangle(e2,pi/2);


% normal
x0 = vector_from_state(...
    [0;0;0],Rq0,Rg0,...
    zeros(3,1),[0;0;0],Rg0.'*[0;0;0]...
);

% twist
% x0 = vector_from_state(...
%     [0;0;0],Rq0,Rg0,...
%     zeros(3,1),[0;0;.1],Rg0.'*[0;0;.1]...
% );

% tumble
% x0 = vector_from_state(...
%     [0;0;0],Rq0,Rg0,...
%     [0;0;5],[0;0;.1],Rg0.'*[0;0;.1]...
% );


u = zeros(6,1);
u(1) = (mg_+mq_)*g_;


segment_dt = .001;
total_dt = 1;

n = floor(total_dt/segment_dt);
state = zeros(n,size(x0,1));

% syms x_dd_des3
% Om_d_des = sym('Om_d_des',[3 1]);
% w_d_des = sym('w_d_des',[3 1]);

u_ff = zeros(6,1);
u_ff = u;

us = zeros(n,6);

for j=1:n
    tspan=segment_dt*(j-1)+[0 segment_dt];
    [~,qs] = ode45(@(t,x) ode(x,u_ff),tspan,x0);   
    us(j,:) = u_ff.';
    x0 = qs(end,:)';
    [xs, Rq, Rg, xs_d, Om, w] = state_from_vector(x0);
    [U, ~, V] = svd(Rq);
    Rq = U * V';
    [U, ~, V] = svd(Rg);
    Rg = U * V';
    state(j,:) = vector_from_state(xs, Rq, Rg, xs_d, Om, w);
%     U(1)
%     xdot = ode(state(j,:)',u)
    

%     M = compute_M_state(Rg,Rq);
%     a = compute_a_state(Rg,Rq,Om,w);
%     B = compute_B_state(Rg,Rq);
%     M2 = [
%      (mg_+mq_)  zeros(1,6);
%       zeros(6,1)       M
%     ];
%     B2 = [
%         e3.'*((mg_+mq_)*Rq*e3) zeros(1,5);
%         B
%     ];
%     a2 = [
%         -g_;
%         a
%     ];
%     


    % M x'' = B u + a

%     x_dd_des3 = 0;
%     Om_d_des = [0;0; 0];
%     w_d_des = zeros([3 1]);
    
%     kp = 10;
%     Om_d_des(3) = -kp*(Om(3)-.5);
%     acc_des = [x_dd_des3; Om_d_des; w_d_des];
%     acc_des(5) = Om_d_des.'*Rq.'*Rg*e1 + Om.'*(-hat(Om)*Rq.'*Rg+Rq.'*Rg*hat(w))*e1;
% 
%     REF = rref([B2 (M2*acc_des - a2)]);
%     
%     if(REF(end,:) == zeros(1,7))
%         u_ff = double(REF(1:6,7));
%     else
%         u_ff = zeros(6,1);
%         disp('Error: Inconsistent system');
%     end


%     u_ff = u;
%     
%     u(5) = 0;
    
%     inv(M)*[B(:,1)]
%     ode(state(j,:)',u_ff)

% dampen in g2 direction only controller
%     kd = 1;
%     u_ff(6) =  -kd * w(2);
%     u_ff(3) =  -kd * w(2);

    
end

%% Plotting


ts = segment_dt*(1:n);
shg;

xg = zeros(n,3);
xq = zeros(n,3);
Oms = zeros(n,3);
ws = zeros(n,3);
constraint_check = zeros(n,1);
energy = zeros(n,1);
input_energy = zeros(n,1);
energy_in = 0;

Jq_ = diag([Jqx_ Jqy_ Jqz_]); % quad inertia
Jg_ = diag([Jgx_ Jgy_ Jgz_]); % gripper inertia

for i=1:n
    [xs,Rq,Rg,xs_d,Om,w] = state_from_vector(state(i,:).');
    xg(i,:) = xs+Lg_*mq_/(mg_+mq_)*Rg*e1;
    xq(i,:) = xg(i,:)' - Lg_*Rg*e1;
    Oms(i,:) = .2*Rq*Om;
    ws(i,:) = .2*Rg*w;
    constraint_check(i,:) = (Rg*e2).'*(Rq*e1);
    energy_in = energy_in + (xs_d * segment_dt).'*(us(i,1) * Rq*e3);
    translational_energy = (mq_+mg_)*xs(3)*g_ + 1/2*(mq_+mg_)*xs_d.'*xs_d;
    rotational_energy = 1/2*(Om.'*Jq_*Om + w.'*Jg_*w);
    energy(i) =  rotational_energy + translational_energy - energy_in;
end

figure(1);
subplot(5,2,1);
plot(ts,state(:,1),ts,state(:,2),ts,state(:,3));
title('system center of mass position');
legend('x','y','z');
ylabel('position [m]')
subplot(5,2,3);
plot(ts,state(:,22),ts,state(:,23),ts,state(:,24));
title('system center of mass velocity');
legend('x','y','z');
ylabel('velocity [m]')
subplot(5,2,5);
plot(ts,constraint_check);
title('f(x) = 0');
legend('rotation');
ylabel('constraint')
subplot(5,2,7);
plot(ts,state(:,25),ts,state(:,26),ts,state(:,27));
title('quad angular velocity');
legend('Om1','Om2','Om3');
ylabel('velocity [m]')
subplot(5,2,9);
plot(ts,state(:,28),ts,state(:,29),ts,state(:,30));
title('gripper angular velocity');
legend('w1','w2','w3');
ylabel('velocity [m]')

subplot(5,2,10);
plot(ts,energy);
title('energy');
ylabel('energy [J]')

subplot(5,2,[2 4 6 8]);

range = 1:20:n;

% for range=1:5:n
    cla;
    hold on;
    view([160 40])

    quiver3(xq(range,1),xq(range,2),xq(range,3),state(range,4),state(range,5),state(range,6),'AutoScale','off','color',[1 0 0]);
    quiver3(xq(range,1),xq(range,2),xq(range,3),state(range,7),state(range,8),state(range,9),'AutoScale','off','color',[0 1 0]);
    quiver3(xq(range,1),xq(range,2),xq(range,3),state(range,10),state(range,11),state(range,12),'AutoScale','off','color',[0 0 1]);

    quiver3(xg(range,1),xg(range,2),xg(range,3),state(range,13),state(range,14),state(range,15),'AutoScale','off','color',[1 0 0]);
    quiver3(xg(range,1),xg(range,2),xg(range,3),state(range,16),state(range,17),state(range,18),'AutoScale','off','color',[0 1 0]);
    quiver3(xg(range,1),xg(range,2),xg(range,3),state(range,19),state(range,20),state(range,21),'AutoScale','off','color',[0 0 1]);
   
%     quiver3(xq(range,1),xq(range,2),xq(range,3),Oms(range,1),Oms(range,2),Oms(range,3),'AutoScale','off','color',[1 0 1]);
%     quiver3(xg(range,1),xg(range,2),xg(range,3),ws(range,1),ws(range,2),ws(range,3),'AutoScale','off','color',[0 1 1]);
    
    plot3(state(range,1),state(range,2),state(range,3),'x');
    axis equal
    axis([-2 4 -2 2 -2 2])
    xlabel('x');
    ylabel('y');
    zlabel('z');
    hold off;
    drawnow;
% end

% figure(2);
% clf;

