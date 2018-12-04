function Rg_des = compute_Rg_des(C_A, C_B, C_G,t,t_f) 
    global compute_Rg_angles;
    basis = (t/t_f).^(0:7);
    a = basis*C_A;
	b = basis*C_B;
	g = basis*C_G;
    Rg_des = compute_Rg_angles(a,b,g);
end