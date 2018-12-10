function [value] = inequality_constraints(C_stacked,t,t_apex,i)
    
    % all values must be less than zero to be admissible  

    global e3 saturation_ compute_Rq_state

    
    [u, state_des] = compute_control(C_stacked.',t,t_apex);
    
    [xs, Rg, th1, th2, xs_d, w, th1_d, th2_d] = state_from_vector(state_des);
    
    
    angle_from_vertical = acos(e3.'*compute_Rq_state(Rg,th1,th2)*e3);
    
    values = [
        (u - saturation_);
        (-u - saturation_);
        (th1 - pi/2);
        (- th1 -pi/2);
        (th2 - pi);
        (- th2);
        (angle_from_vertical - pi/3);
    ];


    value = values(i);

end