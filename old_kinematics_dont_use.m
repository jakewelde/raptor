%% OLD Angular kinematics for gripper

syms angle1(t) angle2(t)

Ra = diff(axisangle(e1,angle1(t))*axisangle(e2,angle2(t)),t);

   Ra_d = simplify(diff(Ra,t));
  Ra_dd = simplify(diff(Ra_d,t));
 Ra_ddd = simplify(diff(Ra_dd,t));
Ra_dddd = simplify(diff(Ra_ddd,t));

syms th1 th2 
syms th1d th2d 
syms th1dd th2dd 
syms th1ddd th2ddd 
syms th1dddd th2dddd

differential = [
    angle1 angle2;
    diff(angle1,t) diff(angle2,t);
    diff(angle1,t,2) diff(angle2,t,2);
    diff(angle1,t,3) diff(angle2,t,3);
    diff(angle1,t,4) diff(angle2,t,4);
];

variant = [
   th1 th2; 
   th1d th2d;
   th1dd th2dd;
   th1ddd th2ddd;
   th1dddd th2dddd;
];

     Ra = subs(Ra,differential,variant)
   Ra_d = subs(Ra_d,differential,variant)
  Ra_dd = subs(Ra_dd,differential,variant)
 Ra_ddd = subs(Ra_ddd,differential,variant)
Ra_dddd = subs(Ra_dddd,differential,variant)

Rq_d = simplify(Rq*hat(Om));
Rq_dd = simplify(Rq_d*hat(Om) + Rq*hat(Om_d));
Rq_ddd = simplify(Rq_dd*hat(Om) + Rq_d*hat(Om_d) + Rq_d*hat(Om_d) + Rq*hat(Om_dd));
Rq_dddd = simplify(Rq_ddd*hat(Om) + Rq_dd*hat(Om_d) + ...
2*(Rq_dd*hat(Om_d) + Rq_d*hat(Om_dd)) + ...
Rq_d*hat(Om_dd) + Rq*hat(Om_ddd));

Rg = Rq*Ra;

Rg_d = simplify(Rq_d * Ra + Rq * Ra_d);

Rg_dd = simplify(Rq_dd * Ra + 2*(Rq_d * Ra_d) + Rq * Ra_dd);

Rg_ddd = ...
simplify(Rq_ddd * Ra + Rq_dd * Ra_d + ...
2*(Rq_dd * Ra_d + Rq_d * Ra_dd) + ...
Rq_d * Ra_dd + Rq * Ra_ddd);

Rg_dddd = ...
simplify(Rq_dddd * Ra + Rq_ddd * Ra_d + ...
Rq_ddd * Ra_d + Rq_dd * Ra_dd + ...
2*(Rq_ddd * Ra_d + Rq_dd * Ra_dd + Rq_dd * Ra_dd + Rq_d * Ra_ddd) + ...
Rq_dd * Ra_dd + Rq_d * Ra_ddd + ...
Rq_d * Ra_ddd + Rq * Ra_dddd);

w_hat   = simplify(Rg.' * Rg_d);
w_d_hat = simplify(Rg_d.' * Rg_d + Rg.' * Rg_dd);
w_dd_hat = simplify(Rg_dd.' * Rg_d + 2*(Rg_d.' * Rg_dd) + Rg.' * Rg_ddd);
w_ddd_hat = simplify(...
Rg_ddd.' * Rg_d + Rg_dd.' * Rg_dd + ...
2*(Rg_dd.' * Rg_dd + Rg_d.' * Rg_ddd) + ...
Rg_d.' * Rg_ddd + Rg.' * Rg_dddd);

    w = unhat(w_hat)
  w_d = unhat(w_d_hat)
 w_dd = unhat(w_dd_hat)
w_ddd = unhat(w_ddd_hat)
return;
