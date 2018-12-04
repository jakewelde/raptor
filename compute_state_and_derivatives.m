function [xs_des, xe_des, Rq_des, Rg_des, xs_d_des, xe_d_des, Om_des, w_des, xs_dd_des, Om_d_des, w_d_des] = compute_state_and_derivatives(C_stacked,t,total_dt)

    global Ls_ e1 e3 g_
    global compute_derivatives hat
    global compute_F_state compute_d_state compute_L_state compute_o_state

    trajectory.x = C_stacked(1:8);
    trajectory.y = C_stacked(9:16);
    trajectory.z = C_stacked(17:24);
    trajectory.a = C_stacked(25:32);
    trajectory.b = C_stacked(33:40);
    trajectory.g = C_stacked(41:48);
    
    angular = compute_angular_des(trajectory.a,trajectory.b,trajectory.g,t,total_dt);
    Rg_des = compute_Rg_des(trajectory.a,trajectory.b,trajectory.g,t,total_dt);

    w_des     = angular(1,:).';
    w_d_des   = angular(2,:).';
    w_dd_des  = angular(3,:).';
    w_ddd_des = angular(4,:).';

    flat_state = zeros(5,3);

    flat_state(:,1) = compute_derivatives(trajectory.x,t,total_dt);
    flat_state(:,2) = compute_derivatives(trajectory.y,t,total_dt);
    flat_state(:,3) = compute_derivatives(trajectory.z,t,total_dt);

         xe_des = flat_state(1,:).'; % end effector position     (0th derivative)
       xe_d_des = flat_state(2,:).'; % end effector velocity     (1st derivative)
      xe_dd_des = flat_state(3,:).'; % end effector acceleration (2nd derivative)
     xe_ddd_des = flat_state(4,:).'; % end effector jerk         (3rd derivative)
    xe_dddd_des = flat_state(5,:).'; % end effector snap         (4th derivative)

    wh = hat(w_des);
    wdh = hat(w_d_des);
    wddh = hat(w_dd_des);
    wdddh = hat(w_ddd_des);

    %% Differential Flatness

    % using the trajectory in flat space, apply the diffeomorphism to compute
    % the trajectory in state space

    % find system center of mass trajectory
    xs_des = xe_des-Ls_*Rg_des*e1;
    xs_d_des = xe_d_des-Ls_*Rg_des*wh*e1;
    xs_dd_des = xe_dd_des-Ls_*Rg_des*(wdh+wh^2)*e1;
    xs_ddd_des = xe_ddd_des-Ls_*Rg_des*(wddh + 3*wh*wdh + wh^3)*e1;
    xs_dddd_des = xe_dddd_des-Ls_*Rg_des*(wdddh + 4*wh*wddh + 6*wh^2*wdh + 3*wdh^2 + wh^4)*e1;




    b3 = (xs_dd_des + g_ * e3)/norm(xs_dd_des + g_ * e3);

    g2 = Rg_des(:,2);
    g3 = Rg_des(:,3);

    % A is full rank unless b3 g2 g3 lie in the same plane, which can only
    % occur when g1 lies in the plane spanned by b1 and b2
    A = [
        transpose(b3);
        transpose(g2);
        transpose(g3);
    ];

    % TODO : is this valid in all orientations?

    b1 = A\[0;0;1];
    b1 = b1 / norm(b1);
    b2 = cross(b3, b1);
    Rq_des = [b1 b2 b3];

    L = compute_L_state(Rg_des,Rq_des);
    o = compute_o_state(Rg_des, Rq_des, w_des, xs_dd_des, xs_ddd_des);

    Om_des = L \ o;

    % find angular acceleration of quadrotor body
    F = compute_F_state(Rg_des, Rq_des,xs_dd_des);
    d = compute_d_state(Rg_des,Rq_des,Om_des,w_des,w_d_des,xs_dd_des,xs_ddd_des,xs_dddd_des);

    % necessary due to numerical issues
    H = [F d];
    [U,S,V] = svd(H);
    if(S(end) ~= 0)
        S(end) = 0;
        H = U*S*V.';
    end
    REF = rref(H);
    Om_d_des = double(REF(1:3,4));

end