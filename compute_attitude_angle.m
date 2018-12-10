function projection = compute_attitude_angle(C,t,tf)
    
    C = C.';

    global compute_com_acceleration

    cx = C(1:8);
    cy = C(9:16);
    cz = C(17:24);
    
%     ca = C(25:32);
    cb = C(33:40);
    cc = C(41:48);

    global Le_ Lg_ mg_ mq_ e3 g_
    
    dt = [ 0, 1, 2*t, 3*t^2, 4*t^3, 5*t^4, 6*t^5, 7*t^6];
    ddt = [ 0, 0, 2, 6*t, 12*t^2, 20*t^3, 30*t^4, 42*t^5];
    t = t.^(0:7);

    xe_dd_des = [
        ddt*cx;
        ddt*cy;
        ddt*cz;
    ];
    
    b = t*cb;   bd = dt*cb;     bdd = ddt*cb;
    c = t*cc;   cd = dt*cc;     cdd = ddt*cc;
 
    xs_dd_des = compute_com_acceleration(Le_,Lg_,b,bd,bdd,c,cd,cdd,mg_,mq_,xe_dd_des(1),xe_dd_des(2),xe_dd_des(3));
    
    projection = acos(e3.'*(xs_dd_des + g_*e3)  / norm(xs_dd_des + g_*e3)); 
    
end