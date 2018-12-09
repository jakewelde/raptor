function val = compute_objective(C_stacked,t_apex,delta_t,z0)
   trajectory.x = C_stacked(1:8);
   trajectory.y = C_stacked(9:16);
   trajectory.z = C_stacked(17:24);
   trajectory.a = C_stacked(25:32);
   trajectory.b = C_stacked(33:40);
   trajectory.g = C_stacked(41:48);

end