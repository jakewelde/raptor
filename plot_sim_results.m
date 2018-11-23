%% Initilalize plotting variables
n = size(state,1);
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

%% Compute constraint checks, energy, and reformat variables

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

%% Plot trajectories
% 
% figure(1);
% subplot(5,2,1);
% plot(ts,state(:,1),ts,state(:,2),ts,state(:,3));
% title('system center of mass position');
% legend('x','y','z');
% ylabel('position [m]')
% subplot(5,2,3);
% plot(ts,state(:,22),ts,state(:,23),ts,state(:,24));
% title('system center of mass velocity');
% legend('x','y','z');
% ylabel('velocity [m]')
% subplot(5,2,5);
% plot(ts,constraint_check);
% title('f(x) = 0');
% legend('rotation');
% ylabel('constraint')
% subplot(5,2,7);
% plot(ts,state(:,25),ts,state(:,26),ts,state(:,27));
% title('quad angular velocity');
% legend('Om1','Om2','Om3');
% ylabel('velocity [1/s]')
% subplot(5,2,9);
% plot(ts,state(:,28),ts,state(:,29),ts,state(:,30));
% title('gripper angular velocity');
% legend('w1','w2','w3');
% ylabel('velocity [1/s]')
% 
% subplot(5,2,10);
% plot(ts,energy);
% title('energy');
% ylabel('energy [J]')
% 
% subplot(5,2,8);
% plot(ts,us);
% title('control efforts');
% legend('f','M1','M2','M3','T1','T2');
% 


%% Visualize Robot

% subplot(5,2,[2 4 6]);
range = 1:20:n;
% for range=1:5:n
    cla;
    hold on;
    view([160 40])

    quiver3(xq(range,1),xq(range,2),xq(range,3),state(range,4),state(range,5),state(range,6),'AutoScale','off','color',[1 0 0]);
    quiver3(xq(range,1),xq(range,2),xq(range,3),state(range,7),state(range,8),state(range,9),'AutoScale','off','color',[0 1 0]);
    quiver3(xq(range,1),xq(range,2),xq(range,3),state(range,10),state(range,11),state(range,12),'AutoScale','off','color',[0 0 1]);

    quiver3(xg(range,1),xg(range,2),xg(range,3),state(range,13),state(range,14),state(range,15),'AutoScale','off','color',[1 1 0]);
    quiver3(xg(range,1),xg(range,2),xg(range,3),state(range,16),state(range,17),state(range,18),'AutoScale','off','color',[0 1 1]);
    quiver3(xg(range,1),xg(range,2),xg(range,3),state(range,19),state(range,20),state(range,21),'AutoScale','off','color',[1 0 1]);

%     quiver3(xq(range,1),xq(range,2),xq(range,3),Oms(range,1),Oms(range,2),Oms(range,3),'AutoScale','off','color',[1 0 1]);
%     quiver3(xg(range,1),xg(range,2),xg(range,3),ws(range,1),ws(range,2),ws(range,3),'AutoScale','off','color',[0 1 1]);

    plot3(state(range,1),state(range,2),state(range,3),'x');
    axis equal
    
    bounds = [min(state(:,1:3))-2*Lg_; max(state(:,1:3))+2*Lg_];
    axis(bounds(:));
    title(total_dt);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    hold off;
    drawnow;
% end
