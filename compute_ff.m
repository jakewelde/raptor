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
