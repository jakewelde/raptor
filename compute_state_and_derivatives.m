function [xe_des, xe_d_des, xs_dd_des, w_d_des, th1_dd_des, th2_dd_des, state_des] = compute_state_and_derivatives(C_stacked,t,total_dt)

    %% Trajectory Evaluation 
    
    % Compute the higher derivatives of flat state using coefficients and
    % also transform the Euler Angle trajectory and derivatives to the
    % space of SO(3), angular velocity, angular acceleration, ...
    
    global Ls_ e1 e2 e3 g_ mg_ mq_
    global compute_derivatives hat
    global compute_Rq_state compute_Om_state;
    global compute_thd_from_Om12_state compute_thdd_from_Om_d12_state
    
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
    
    components = Rg_des.'*b3; 
    % == [ 
    %  cos(th1)*sin(th2); 
    % -sin(th1); 
    %  cos(th1)*cos(th2);
    % ];

    th2_des = atan2(components(1),components(3));
    sinth1 = components(2);
    costh1 = components(3) / cos(th2_des);
    th1_des = atan2(sinth1,costh1);

    Rq_des = compute_Rq_state(Rg_des,th1_des,th2_des);
    
    Om_des12 = 1 / (norm(xs_dd_des+g_*e3)) * [-(Rq_des*e2).'; (Rq_des*e1).']*xs_ddd_des;
    
    thds = compute_thd_from_Om12_state(Rg_des,w_des,th1_des,th2_des,Om_des12(1),Om_des12(2));
    th1_d_des = thds(1);
    th2_d_des = thds(2);

    
    Om_des = compute_Om_state(Rg_des,th1_des,th2_des,th1_d_des,th2_d_des,w_des);
    

    thrust    = (Rq_des*e3).' * (mg_+mq_) * (xs_dd_des + g_*e3);
    thrust_d  = (Rq_des*e3).' * (mg_+mq_) * xs_ddd_des;
    thrust_dd = (Rq_des*e3).' * (mg_+mq_) * xs_dddd_des - thrust*e3.'*hat(Om_des)^2*e3;

    Om_d_des12 = ...
    [0 -1 0; 1 0 0]*Rq_des.'*(...
    -Rq_des * ( hat(Om_des)^2 ) * e3 + ...
    1/thrust * (...
    (mg_+mq_) * xs_dddd_des + ...
    -thrust_d * 2 * Rq_des * hat(Om_des) * e3  + ...
    -thrust_dd * Rq_des * e3...
    ));

    thdds = compute_thdd_from_Om_d12_state(Rg_des,w_des,w_d_des,th1_des,th2_des,th1_d_des,th2_d_des,Om_d_des12(1),Om_d_des12(2));
    th1_dd_des = thdds(1);
    th2_dd_des = thdds(2);
    
    
    state_des = vector_from_state(xs_des, Rg_des, th1_des, th2_des, xs_d_des, w_des, th1_d_des, th2_d_des);

end