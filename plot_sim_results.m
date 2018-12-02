%% Initilalize plotting variables
n = size(state,1);
ts = segment_dt*(1:n);

show_profile_plots = true;


figure(5);
clf;
plot(...
    ts,xs_rec(1,:),...
    ts,xs_rec(4,:),...
    ts,xs_rec(7,:),...
    ts,xs_rec(10,:),...
    ts,xs_rec(13,:)...
)
legend('position','velocity','acceleration','jerk','snap');



xg = zeros(n,3);
ve = zeros(n,3);
xq = zeros(n,3);
Oms = zeros(n,3);
ws = zeros(n,3);
constraint_check = zeros(n,1);
energy = zeros(n,1);
input_energy = zeros(n,1);
power_in = 0;

Jq_ = diag([Jqx_ Jqy_ Jqz_]); % quad inertia
Jg_ = diag([Jgx_ Jgy_ Jgz_]); % gripper inertia

%% Compute constraint checks, energy, and reformat variables
global Le_ Lg_

for i=1:n
    [xs,Rq,Rg,xs_d,Om,w] = state_from_vector(state(i,:).');
    xg(i,:) = xs+Lg_*mq_/(mg_+mq_)*Rg*e1;
    xq(i,:) = xg(i,:)' - Lg_*Rg*e1;
    ve(i,:) = xs_d + Rg*cross(w,Rg.'*(xg(i,:).'-xs)*Le_/Lg_);
    xe(i,:) = xs + Ls_*Rg*e1;
    Oms(i,:) = .2*Rq*Om;
    ws(i,:) = .2*Rg*w;
    constraint_check(i,:) = (Rg*e2).'*(Rq*e1);
    power_in = (...
        us(i,1) * xs_d.'*Rq*e3 + ...     % f * v
        us(i,2) * Om.'*e1+...            % M1 * Om1
        us(i,3) * Om.'*e2+...            % M2 * Om2 
        us(i,4) * Om.'*e3+...            % M3 * Om3
        us(i,5) * (Rg*w).'*(Rq*e1)+...   % T1 * component of w in b1 direction
        us(i,6) * w.'*e2 + ...           % T2 * w2
        -us(i,5) * Om.'*e1+...           % T1 * Om1
        -us(i,6) * (Rq*Om).'*(Rg*e2)...  % T2 * component of Om in g2 direction
    );
    % TODO : NEED TO *FIX* TORQUE CONTRIBUTIONS TO ENERGY 
    input_energy(i+1) = input_energy(i) + segment_dt*power_in;
    translational_energy = (mq_+mg_)*xs(3)*g_ + 1/2*(mq_+mg_)*xs_d.'*xs_d;
    rotational_energy = 1/2*(Om.'*Jq_*Om + w.'*Jg_*w);
    energy(i) =  rotational_energy + translational_energy - input_energy(i);
end

%% Plot trajectories

figure(2);
shg;

if(show_profile_plots)
    subplot(4,2,1);
    plot(ts,state(:,1),ts,state(:,2),ts,state(:,3));
    title('system center of mass position');
    legend('x','y','z');
    ylabel('position [m]')
    subplot(4,2,3);
    plot(ts,state(:,22),ts,state(:,23),ts,state(:,24));
    title('system center of mass velocity');
    legend('x','y','z');
    ylabel('velocity [m]')
    subplot(4,2,5);
    plot(ts,state(:,25),ts,state(:,26),ts,state(:,27));
    title('quad angular velocity');
    legend('Om1','Om2','Om3');
    ylabel('velocity [1/s]')
    subplot(4,2,7);
    plot(ts,state(:,28),ts,state(:,29),ts,state(:,30));
    title('gripper angular velocity');
    legend('w1','w2','w3');
    ylabel('velocity [1/s]')
    legend('x','y','z');

    subplot(4,2,2);
    plot(ts,us);
    title('control efforts');
    legend('f','M1','M2','M3','T1','T2');

    subplot(4,2,4);
    plot(ts,ve);
    title('end effector velocity in world');
    legend('x','y','z');
    ylabel('velocity[m/s]')
    
    subplot(4,2,6);
    plot(ts,energy-energy(1));
    title('change in energy');
    ylabel('energy [J]')

    subplot(4,2,8);
    plot(ts,constraint_check);
    title('f(x) = 0');
    legend('rotation');
    ylabel('constraint')    
    
end



%% Visualize Robot
% state(:,4:21) = state(:,4:21);
figure(1)
% range = 1:20:n;
if(exist('az') && exist('el'))
    view(az,el)
end

step = 50;
for range=[1:step:n n]
    [az,el]=view;
    cla;
    hold on;
    view(az,el);
    d_trail = 50;
    pts = 1:d_trail:range;
    plot3(xe(pts,1),xe(pts,2),xe(pts,3),'.g','MarkerSize',15)
    plot3(xq(pts,1),xq(pts,2),xq(pts,3),'.b','MarkerSize',15)
    plot3(state(pts,1),state(pts,2),state(pts,3),'.m','MarkerSize',15)
    
    plot3(xs_rec(1,:),xs_rec(2,:),xs_rec(3,:),'.')
    plot3(xe_rec(1,:),xe_rec(2,:),xe_rec(3,:),'.')
    plot3(xq_rec(1,:),xq_rec(2,:),xq_rec(3,:),'.')

%     quiver3(xq(range,1),xq(range,2),xq(range,3),state(range,4),state(range,5),state(range,6),'AutoScale','off','color',[1 0 0]);
%     quiver3(xq(range,1),xq(range,2),xq(range,3),state(range,7),state(range,8),state(range,9),'AutoScale','off','color',[0 1 0]);
%     quiver3(xq(range,1),xq(range,2),xq(range,3),state(range,10),state(range,11),state(range,12),'AutoScale','off','color',[0 0 1]);
% 
%     quiver3(xg(range,1),xg(range,2),xg(range,3),state(range,13),state(range,14),state(range,15),'AutoScale','off','color',[1 0 0]);
%     quiver3(xg(range,1),xg(range,2),xg(range,3),state(range,16),state(range,17),state(range,18),'AutoScale','off','color',[0 1 0]);
%     quiver3(xg(range,1),xg(range,2),xg(range,3),state(range,19),state(range,20),state(range,21),'AutoScale','off','color',[0 0 1]);

%     quiver3(xq(range,1),xq(range,2),xq(range,3),Oms(range,1),Oms(range,2),Oms(range,3),'AutoScale','off','color',[1 0 1]);
%     quiver3(xg(range,1),xg(range,2),xg(range,3),ws(range,1),ws(range,2),ws(range,3),'AutoScale','off','color',[0 1 1]);

%     plot3(state(range,1),state(range,2),state(range,3),'x');

    draw_robot(state(range,:));
    axis equal    
    bounds = [min(state(:,1:3))-3*Lg_; max(state(:,1:3))+3*Lg_];
    axis(bounds(:));
    xlabel('x');
    ylabel('y');
    zlabel('z');
    hold off;
    drawnow;
    
end
