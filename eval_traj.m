function val = eval_traj(C,t)
    val = t.^(0:7)*C.';
end