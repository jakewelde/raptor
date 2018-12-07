%% Initilalize plotting variables
n = size(state,1);
ts = segment_dt*(1:n);

show_profile_plots = true;

xg = zeros(n,3);
xe = zeros(n,3);
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
global Le_ Lg_ Ls_

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

    input_energy(i+1) = input_energy(i) + segment_dt*power_in;
    translational_energy = (mq_+mg_)*xs(3)*g_ + 1/2*(mq_+mg_)*xs_d.'*xs_d;
    rotational_energy = 1/2*(Om.'*Jq_*Om + w.'*Jg_*w);
    energy(i) =  rotational_energy + translational_energy - input_energy(i);
end

%% Plot trajectories

if(show_profile_plots)
    
    figure(2);
    clf;
    subplot(4,2,1);
    cla;
    hold on;
%     plot(ts,state(:,1),'r',ts,state(:,2),'g',ts,state(:,3),'b');
%     plot(ts,xs_rec(1,:),'r--',ts,xs_rec(2,:),'g--',ts,xs_rec(3,:),'b--');
    plot(ts,state(:,1).'-xs_rec(1,:),'r',ts,state(:,2).'-xs_rec(2,:),'g',ts,state(:,3).'-xs_rec(3,:),'b');
    hold off;
    title('system center of mass error');
    ylabel('position [m]')
    
    subplot(4,2,3);
    cla;
    hold on;
%     plot(ts,state(:,22),'r',ts,state(:,23),'g',ts,state(:,24),'b');
%     plot(ts,xs_rec(4,:),'r--',ts,xs_rec(5,:),'g--',ts,xs_rec(6,:),'b--');
    plot(ts,state(:,22).'-xs_rec(4,:),'r',ts,state(:,23).'-xs_rec(5,:),'g',ts,state(:,24).'-xs_rec(6,:),'b');
    hold off;
    title('system center of mass velocity');
    ylabel('velocity [m/s]')
    
    subplot(4,2,5);
    cla;
    hold on;
%     plot(ts,state(:,25),'r',ts,state(:,26),'g',ts,state(:,27),'b');
%     plot(ts,Om_rec(1,:),'r--',ts,Om_rec(2,:),'g--',ts,Om_rec(3,:),'b--');
    plot(ts,state(:,25).'-Om_rec(1,:),'r',ts,state(:,26).'-Om_rec(2,:),'g',ts,state(:,27).'-Om_rec(3,:),'b');
    hold off;
    
    title('quad angular velocity error');
    ylabel('velocity [1/s]')
    
    subplot(4,2,7);
    cla;
    hold on;
%     plot(ts,state(:,28),'r',ts,state(:,29),'g',ts,state(:,30),'b');
%     plot(ts,w_rec(1,:),'r--',ts,w_rec(2,:),'g--',ts,w_rec(3,:),'b--');
    plot(ts,state(:,28).'-w_rec(1,:),'r',ts,state(:,29).'-w_rec(2,:),'g',ts,state(:,30).'-w_rec(3,:),'b');
    hold off;
    title('gripper angular velocity error');
    ylabel('velocity [1/s]')

    subplot(4,2,2);
    cla;
    hold on;
%     plot(ts,xe(:,1),'r',ts,xe(:,2),'g',ts,xe(:,3),'b');
%     plot(ts,xe_rec(1,:),'r--',ts,xe_rec(2,:),'g--',ts,xe_rec(3,:),'b--');
    plot(ts,xe(:,1).'-xe_rec(1,:),'r',ts,xe(:,2).'-xe_rec(2,:),'g',ts,xe(:,3).'-xe_rec(3,:),'b');
    hold off;
    title('end effector error in world');
    ylabel('position [m]')

    subplot(4,2,4);
    cla;
    hold on;
%     plot(ts,ve(:,1),'r',ts,ve(:,2),'g',ts,ve(:,3),'b');
%     plot(ts,xe_rec(4,:),'r--',ts,xe_rec(5,:),'g--',ts,xe_rec(6,:),'b--');
    plot(ts,ve(:,1).'-xe_rec(4,:),'r',ts,ve(:,2).'-xe_rec(5,:),'g',ts,ve(:,3).'-xe_rec(6,:),'b');
    hold off;
    title('end effector velocity in world');
    ylabel('velocity [m/s]')

   
    subplot(4,2,6);
    plot(ts,us);
    title('control efforts');
    legend('f','M1','M2','M3','T1','T2');
    ylabel('[N] or [N \cdot m]')

    subplot(4,2,8);
    plot(ts,constraint_check);
    title('g_2 \cdot b_1 = 0');
    ylabel('[radians]')    
    
end

shg;
drawnow;


%% Visualize Robot



figure(1)
if(exist('az') && exist('el'))
    view(az,el)
end

step = 900;
for range=[1:step:n n]
    [az,el]=view;
    cla;
    hold on;
    view(az,el);
    d_trail = 70;
    pts = [1:d_trail:range range];
    plot3(xe(pts,1),xe(pts,2),xe(pts,3),'.','MarkerSize',15,'color',[.9 .2 .6])
    plot3(xe_rec(1,:),xe_rec(2,:),xe_rec(3,:),'color',[.9 .2 .6])
    plot3(xs_rec(1,:),xs_rec(2,:),xs_rec(3,:),'color',[.9 .2 .6])
    plot3(ball_position(1,pts),ball_position(2,pts),ball_position(3,pts),'k');
    plot3(ball_position(1,pts(end)),ball_position(2,pts(end)),ball_position(3,pts(end)),...
        'o','MarkerEdgeColor','k',...
        'MarkerFaceColor','k',...
        'MarkerSize',6);
    plot3(z_apex(1),z_apex(2),z_apex(3),'rx');

    draw_robot(state(range,:));
    axis equal    
%     bounds = [min([(min(state(:,1:3))-3*Lg_).'  min(ball_position.').']); max([(max(state(:,1:3))+3*Lg_).'   max(ball_position.').'])];
    bounds = [min(state(:,1:3)-3*Lg_); max(state(:,1:3)+3*Lg_)];
%     bounds = [min([state(:,1:3)-3*Lg_;ball_position.']); max([state(:,1:3)+3*Lg_;ball_position.'])];
    axis(bounds(:));
    xlabel('x');
    ylabel('y');
    zlabel('z');
    hold off;
    drawnow;
    
end
