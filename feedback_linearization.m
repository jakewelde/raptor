function u = feedback_linearization(xs_des, xe_des, Rg_des, th1_des, th2_des, xs_d_des, xe_d_des, w_des, th1_d_des, th2_d_des, xs_dd_des, w_d_des, th1_dd_des, th2_dd_des, state)

    global e3
    global mg_ mq_ g_
    global compute_M_state compute_B_state compute_a_state
    global compute_Rq_state;
    global unhat
    
    if(exist('state','var'))
        [xs, Rg, th1, th2, xs_d, w, th1_d, th2_d] = state_from_vector(state);
        Rq_des = compute_Rq_state(Rg,th1,th2);
        xe = [
          (Rq_des*e3).'*(xs - xs_des);
          (Rq_des*e3).'*(xs_d - xs_d_des);
          1/2*unhat(Rg_des.'*Rg - Rg.'*Rg_des); % gripper rotation error   
          w - Rg.'*Rg_des*w_des;
          th1-th1_des;
          th1_d-th1_d_des;
          th2-th2_des;
          th2_d-th2_d_des;
        ];
    else
        xe = zeros(12,1);
        Rg = Rg_des;
        th1 = th1_des;
        th2 = th2_des;
        th1_d = th1_d_des;
        th2_d = th2_d_des;
        w = w_des;
        Rq_des = compute_Rq_state(Rg,th1,th2);
    end
        
    acc_nominal = [(Rq_des*e3).'*xs_dd_des; w_d_des; th1_dd_des; th2_dd_des];

    % Compute numerical function that represents the evolution of the system's
    % state at the current point in state space, but is still a function of the
    % inputs of the system
    M = [
      (mg_+mq_)   zeros(1,5);
      zeros(5,1) compute_M_state(Rg,th1,th2)
    ];
    B = [
        1 zeros(1,5);
        compute_B_state(Rg,th1,th2);
    ];
    a = [
        (Rq_des*e3).'*(-(mg_+mq_)*g_*e3);
        compute_a_state(Rg,w,th1,th2,th1_d,th2_d);
    ];

    % M x'' = B u + a(x)

    % x'' = inv(M)*(B u + a(x))
    % x_e = x - x_nom
    % xe' = x' - x_nom'
    % xe'' = x'' - x_nom''
    % sub in dynamics for x''
    % xe'' = inv(M)*(B u + a(x)) - x_nom''
    % xe'' = inv(M)*a(x) - x_nom'' + inv(M)*B*u
    % Pick u so that xe'' = A x and A is stable
    % Propose u = inv(B)*( -a(x) + M*x_nom'' + M*x_e_des'' ) 
    %  sub back in, get  
    % xe'' = inv(M)*a(x) - x_nom'' + inv(M)*B*inv(B)*( -a(x) + M*x_nom'' + M*x_e_des'' ) 
    % xe'' = inv(M)*a(x) - x_nom'' - inv(M)*(a(x)) + x_nom'' + x_e_des'' 
    % xe'' = inv(M)*a(x) - inv(M)*(a(x)) + x_nom'' - x_nom'' + x_e_des'' 
    % xe'' = x_e_des'' 
    % 
    % So to stabilize xe, just choose  
    % x_e_des'' = A x_e
    % where A has negative eigenvalues. Or choose a similar nonlinear
    % system.
    % 
    % Plugging that back in,
    % u = inv(B)*( -a(x) + M*x_nom'' + M*A x_e ) 
    % u = inv(B)*( -a(x) + M*x_nom'' + M*A (x - x_nom )) 
    
    
    % Nominal case:
    % 
    % M x_nom'' = B u_ff + a(x)
    % u_ff = inv(B) * (-a(x) + M x_nom'')

    
    Kx = 1e5;
    Kv = 2*sqrt(Kx);

    Krg = 1e5;
    Kw = 2*sqrt(Krg);

    Kp = 1e5;
    Kd = 2*sqrt(Kp);
    K = -[
        Kx  Kv  0   0   0   0   0   0   0   0   0   0;
        0   0   Krg 0   0   Kw  0   0   0   0   0   0;
        0   0   0   Krg 0   0   Kw  0   0   0   0   0;
        0   0   0   0   Krg 0   0   Kw  0   0   0   0;
        0   0   0   0   0   0   0   0   Kp  Kd  0   0;
        0   0   0   0   0   0   0   0   0   0   Kp  Kd;
    ];

    u_ff = B \ (- a + M*(acc_nominal));
    u_fb = B \ (M*K*xe);
%     u_fb
    
    u = u_ff + u_fb;
    
%     if(exist('state','var')) % real feedback scenario, not nominal
%         saturation = [g_*(mq_+mg_)*2 10 10 1 1 1].';
%         u = max(min(u,saturation),-saturation);
%     end
    
end
