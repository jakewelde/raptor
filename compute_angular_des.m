function angular_derivs = compute_angular_des(C_A, C_B, C_G, t, t_f)
    global compute_derivatives compute_angular_derivatives
    m = [
        compute_derivatives(C_A,t,t_f).';
        compute_derivatives(C_B,t,t_f).';
        compute_derivatives(C_G,t,t_f).'
    ];

    angular_derivs = compute_angular_derivatives(...
        m(1,1),m(1,2),m(1,3),m(1,4),m(1,5),...
        m(2,1),m(2,2),m(2,3),m(2,4),m(2,5),...
               m(3,2),m(3,3),m(3,4),m(3,5)...
    );
end