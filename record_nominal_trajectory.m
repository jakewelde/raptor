% record spatial trajectories for visualization
xe_rec(1:3,j) = xe_des(1:3);
xe_rec(4:6,j) = xe_des(4:6);
xs_rec(1:3,j) = xs_des(1:3);
xs_rec(4:6,j) = xs_des(4:6);

% record angular trajectories for visualization
w_rec(:,j) = w_des;
Om_rec(:,j) = Om_des;