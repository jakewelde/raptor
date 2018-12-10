cx = stacked(1:8);
cy = stacked(9:16);
cz = stacked(17:24);
ca = stacked(25:32);
cb = stacked(33:40);
cg = stacked(41:48);

delta_t = .1;
time_range = 0:.01:1+delta_t;
positions = zeros(6,size(time_range,2));

states = zeros(size(time_range,2),22);
us = zeros(size(time_range,2),6);

for i=1:size(time_range,2)
    positions(:,i) = [
        [1 0 0 0 0]*compute_derivatives(cx,time_range(i),1);
        [1 0 0 0 0]*compute_derivatives(cy,time_range(i),1);
        [1 0 0 0 0]*compute_derivatives(cz,time_range(i),1)
        [1 0 0 0 0]*compute_derivatives(ca,time_range(i),1)
        [1 0 0 0 0]*compute_derivatives(cb,time_range(i),1)
        [1 0 0 0 0]*compute_derivatives(cg,time_range(i),1)
    ];

    [u, state_des] = compute_control(stacked,time_range(i),total_dt);
    states(i,:) = state_des.';
    us(i,:) = u.';
end
figure(3);
clf;
subplot(3,2,[1 3 5]);
hold on;
plot3(positions(1,:),positions(2,:),positions(3,:))

plot3(states(:,1),states(:,2),states(:,3))

plot3(z_apex(1),z_apex(2),z_apex(3),'x');
plot3(positions(1,floor((1-delta_t)/.01)),...
    positions(2,floor((1-delta_t)/.01)),...
    positions(3,floor((1-delta_t)/.01)),'o');

plot3(positions(1,end),...
    positions(2,end),...
    positions(3,end),'o');

plot3(positions(1,1),...
    positions(2,1),...
    positions(3,1),'o');

vec = z_d_apex/norm(z_d_apex);
quiver3(z_apex(1),z_apex(2),z_apex(3),vec(1),vec(2),vec(3));
hold off;
axis equal;

xlabel('x'); ylabel('y'); zlabel('z');
title(retcode)

subplot(3,2,2);
plot(time_range,positions(4:6,:));
title('Euler Angles');
subplot(3,2,4);
plot(time_range,states(:,13).',time_range,states(:,14).');
title('Joint Angles');
subplot(3,2,6);
plot(time_range,us.');
title('Control Inputs');


