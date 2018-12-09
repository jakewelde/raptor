function [C_stacked, min, retcode] = trajectory_optimization(z0,C0,delta_t)
    
    % TODO : nondimensionalize! 

    global g_
    global e3
    global cost_function
    
    t_apex = z0(6)/g_;

    addpath('/usr/local/lib/matlab/')
   
    zf = (z0(1:3) + z0(4:6) - 1/2*g_*e3);
    zdf = (z0(4:6) - 1*g_*e3);
    
    opt.algorithm = NLOPT_LN_COBYLA;
    opt.lower_bounds = -inf*ones(1,48);
    opt.min_objective = @(C) cost_function(C,t_apex,delta_t,z0);
    % opt.fc = { (@(x) myconstraint(x,2,0)), (@(x) myconstraint(x,-1,1)) };
    % opt.fc_tol = [1e-8, 1e-8];
    opt.h = {
        @(C) C(1), @(C) C(9), @(C) C(17), @(C) C(25), @(C) C(33)-pi/2, @(C) C(41),...
        @(C) C(2), @(C) C(10), @(C) C(18), @(C) C(26), @(C) C(34), @(C) C(42),...
        @(C) 2*C(3), @(C) 2*C(11), @(C) 2*C(19), @(C) 2*C(27), @(C) 2*C(35), @(C) 2*C(43),... 
        @(C) 6*C(4), @(C) 6*C(12), @(C) 6*C(20), @(C) 6*C(28), @(C) 6*C(36), @(C) 6*C(44),...
        @(C) sum(C(1:8)) - zf(1), @(C) sum(C(9:16)) - zf(2), @(C) sum(C(17:24)) - zf(3),...
        @(C) C(1:8)*(0:7).' - zdf(1), @(C) C(9:16)*(0:7).' - zdf(2), @(C) C(17:24)*(0:7).' - zdf(3),...
    };
    opt.xtol_rel = 1e-5;
    opt.verbose = 0;

    [C_stacked, min, retcode] = nlopt_optimize(opt, C0);

end