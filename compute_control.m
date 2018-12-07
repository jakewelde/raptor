function [u, xe_rec, xs_rec, w_rec, Om_rec] = compute_control(C_stacked,t,total_dt,state)
    
    [xs_des, xe_des, Rq_des, Rg_des, xs_d_des, xe_d_des, Om_des, w_des, xs_dd_des, Om_d_des, w_d_des]...
        = compute_state_and_derivatives(C_stacked,t,total_dt);
%     if(nargin > 3)
%         u = feedback_linearization(Rq_des, Rg_des, Om_des, w_des, xs_dd_des,Om_d_des, w_d_des,state);
%     else
    u = feedback_linearization(xs_des, Rq_des, Rg_des, xs_d_des, Om_des, w_des, xs_dd_des,Om_d_des, w_d_des,state);
%     end
    
    if (nargout > 1)
        xe_rec = [xe_des; xe_d_des];
        xs_rec = [xs_des; xs_d_des];

        w_rec = w_des;
        Om_rec = Om_des;
    end

end