function u = feedback_linearization(xs_des, Rq_des, Rg_des, xs_d_des, Om_des, w_des, xs_dd_des,Om_d_des, w_d_des, state)

    global e3

    global mg_ mq_ g_
    global compute_M_state compute_B_state compute_a_state

    % Feedback Linearization

    [xs, Rq, Rg, xs_d, Om, w] = state_from_vector(state);

    acc_des = [(Rq*e3).'*xs_dd_des; Om_d_des; w_d_des];

    % Compute numerical function that represents the evolution of the system's
    % state at the current point in state space, but is still a function of the
    % inputs of the system
    M = [
      (mg_+mq_)   zeros(1,6);
      zeros(6,1) compute_M_state(Rg,Rq);
    ];
    B = [
        1 zeros(1,5);
        compute_B_state(Rg,Rq);
    ];
    a = [
        (Rq*e3).'*(-(mg_+mq_)*g_*e3);
        compute_a_state(Rg,Rq,Om,w);
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
    
%     if(nargin < 7)
%         xe = zeros(12,1); 
%     else
        [~,~,~,~,Om,w] = state_from_vector(state);
        xe = [
          (Rq*e3).'*(xs - xs_des);
          (Rq*e3).'*(xs_d - xs_d_des);
          zeros(3,1); % body rotation error   
          zeros(3,1); % gripper rotation error   
          Om - Om_des;
          w - w_des
        ];

%     end
%     -4*sqrt(1e17)
    Kp = 90;
    Kd = 2*sqrt(Kp);
    K = -[
%         [-1e17 -2*sqrt(1e17) zeros(1,12)];
        [Kp Kd zeros(1,12)];
        zeros(6,14)
    ];
%     disp((M*K*xe).')

%         -eye(3) -eye(3) zeros(3,6);
%         zeros(3,6) -eye(3) -eye(3)

%     G = [B (- a + M*(acc_des))];
%     B \ M*K*xe
    G = [B (- a + M*(acc_des + K*xe))];
    [U,S,V] = svd(G);
    
    % necessary due to numerical issues
    if(S(end) ~= 0)
        S(end) = 0;
        G = U*S*V.';
    end
    REF = rref(G);
    
    u = double(REF(1:6,7));  
%     disp([1 zeros(1,6)]*(M*K*xe));
%     disp(u(1)-g_*(mg_+mq_));
%     disp([rank(M) rank(B) size(B)])
    
end
