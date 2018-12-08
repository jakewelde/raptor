%% 

w = sym('w',[3 1]);
w_d = sym('w_d',[3 1]);
w_dd = sym('w_dd',[3 1]);
w_ddd = sym('w_ddd',[3 1]);

Rg = sym('Rg',[3 3]); % rotation from quad to world

Rg_d = simplify(Rg*hat(w));
Rg_dd = simplify(Rg_d*hat(w) + Rg*hat(w_d));
Rg_ddd = simplify(Rg_dd*hat(w) + Rg_d*hat(w_d) + Rg_d*hat(w_d) + Rg*hat(w_dd));
Rg_dddd = simplify(Rg_ddd*hat(w) + Rg_dd*hat(w_d) + ...
2*(Rg_dd*hat(w_d) + Rg_d*hat(w_dd)) + ...
Rg_d*hat(w_dd) + Rg*hat(w_ddd));

syms angle1(t) angle2(t)

Ra = axisangle(e2,angle2(t))*axisangle(e1,angle1(t));

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

Rq = Rg*Ra;

Rq_d = Rg_d*Ra + Rg*Ra_d;

Rq_dd = Rg_dd*Ra + 2*Rg_d*Ra_d + Rg*Ra_dd;

Om_hat   = simplify(Rq.' * Rq_d);
Om_d_hat = simplify(Rq_d.' * Rq_d + Rq.' * Rq_dd);

Om = unhat(Om_hat)
Om_d = unhat(Om_d_hat)

% Are some of these useful? Unclear!

compute_Rq_all = matlabFunction(Rq);
compute_Om_all = matlabFunction(Om);
compute_Om_d_all = matlabFunction(Om_d);

compute_Rq_state = @(Rg,th1,th2) compute_Rq_all(Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),th1,th2);
compute_Om_state = @(Rg,th1,th2,th1d,th2d,w) compute_Om_all(Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),th1,th2,th1d,th2d,w(1),w(2),w(3));
compute_Om_d_state = @(Rg,th1,th2,th1d,th2d,th1dd,th2dd,w,w_d) compute_Om_d_all(Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),th1,th2,th1d,th2d,th1dd,th2dd,w(1),w(2),w(3),w_d(1),w_d(2),w_d(3));

