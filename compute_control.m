function [u, xe_rec, xs_rec, w_rec] = compute_control(C_stacked,t,total_dt,state)
    
    [xs_des, xe_des, Rg_des, th1_des, th2_des, xs_d_des, xe_d_des, w_des, th1_d_des, th2_d_des, xs_dd_des, w_d_des, th1_dd_des, th2_dd_des]...
        = compute_state_and_derivatives(C_stacked,t,total_dt);
    u = feedback_linearization(xs_des, xe_des, Rg_des, th1_des, th2_des, xs_d_des, xe_d_des, w_des, th1_d_des, th2_d_des, xs_dd_des, w_d_des, th1_dd_des, th2_dd_des, state);
    
    if (nargout > 1)
        xe_rec = [xe_des; xe_d_des];
        xs_rec = [xs_des; xs_d_des];

        w_rec = w_des;
%         Om_rec = Om_des;
    end

end