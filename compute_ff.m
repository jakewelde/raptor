function u_ff = compute_ff(Rq_des, Rg_des, Om_des, w_des, xs_dd_des,Om_d_des, w_d_des)

    global e3

    global mg_ mq_ g_
    global compute_M_state compute_B_state compute_a_state

    % Feedback Linearization

    acc_des = [xs_dd_des(3); Om_d_des; w_d_des];

    % Compute numerical function that represents the evolution of the system's
    % state at the current point in state space, but is still a function of the
    % inputs of the system
    M = [
      (mg_+mq_)   zeros(1,6);
      zeros(6,1) compute_M_state(Rg_des,Rq_des);
    ];
    B = [
        e3.'*(Rq_des*e3) zeros(1,5);
        compute_B_state(Rg_des,Rq_des);
    ];
    a = [
        -(mg_+mq_)*g_;
        compute_a_state(Rg_des,Rq_des,Om_des,w_des);
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
    
    G = [B (M*acc_des - a)];
    [U,S,V] = svd(G);

    % necessary due to numerical issues
    if(S(end) ~= 0)
        S(end) = 0;
        G = U*S*V.';
    end
    REF = rref(G);

    u_ff = double(REF(1:6,7));
    
end
