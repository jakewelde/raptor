%% Methodology of Derivation

% Create an equation of the form M(x) [Om_d; w_d] = a(x) + B(x)*u so that
% we can find [Om_d; w_d].
%
% The first three rows of this vector equation will come from the Euler
% equations for the two bodies. They can each be solved for the moments
% between them and set equal, giving three equations. 
%
% Then, two equations will be found by projecting the moments between the
% two bodies onto the arm joint axes. This will be equal to the control
% inputs that are the torques on those axes. To get this equation in terms
% of Om_d, use the angular dynamics of one of the bodies solved for the
% moments, as computed above, in the projection.
% 
% The final constraint comes from the kinematics of the robot. Becuase
% there are only two joints between the two bodies, there is a constraint
% on their relative angular velocity since it cannot be any arbitrary
% direction. In particular, consider the difference between their angular
% velocities, expressed in a common frame. The projection of this vector
% onto the direction that is instantaneously perpendicular to both motor
% axes must be zero, since there is no freedom of movement between the two
% bodies along that direction. In order to get a constraint on
% acceleration, differentiate this velocity constraint.

clear;

%% Model Parameters

disp('Generating preliminaries...');

syms Jqx Jqy Jqz Jgx Jgy Jgz;
Jq = diag([Jqx Jqy Jqz]); % quad inertia
Jg = diag([Jgx Jgy Jgz]); % gripper inertia

syms mq; % quad mass
syms mg; % gripper mass
syms Lg; % distance from center of actuation to gripper COM
syms Le; % distance from center of actuation to end effector

Ls = mq/(mg+mq)*Lg - Le; % distance from end effector to system center of mass

syms g;

%% State variables

% POSITION
xs = sym('xs',[3 1]); % position of center of mass
Rg = sym('Rg',[3 3]); % rotation from gripper to world
syms th1 th2          % joint angles

% VELOCITY
xs_d = sym('xs_d',[3 1]); % velocity of center of mass
w = sym('w',[3 1]);       % angular velocity of quadrotor body relative to world, expressed in quadrotor frame
syms thd1 thd2            % joint velocities

% ACCELERATION (not truly part of state but used below)
xs_dd = sym('xs_dd',[3 1]); % acceleration of center of mass
w_d = sym('w_d',[3 1]);     % angular acceleration of quadrotor body relative to world, expressed in quadrotor frame
syms thdd1 thdd2            % joint accelerations

% Higher derivatives

w_dd = sym('w_dd',[3 1]);
w_ddd = sym('w_ddd',[3 1]);

%% Control Inputs

syms f % net thrust
syms M1 M2 M3 
Mq = [M1;M2;M3]; % rotor thrust moments around quadrotor body axes
syms T1 % torque at joint 1 - closer to quadrotor body, rotates around b1 
syms T2 % torque at joint 2 - closer to gripper body, rotates around g2 (which is rotated by joint 1 so is not constant in quadrotor body frame)

u = [f; M1; M2; M3; T1; T2];

%% Operators
global hat unhat
hat = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
unhat = @(v_hat) [v_hat(3,2); v_hat(1,3); v_hat(2,1)];

%% Unit Vectors

global e1 e2 e3
e1 = [1; 0; 0];
e2 = [0; 1; 0];
e3 = [0; 0; 1];


%% Derive angular kinematics for quad based on minimal state representation

disp('Deriving angular quantities...');

Rg_d = simplify(Rg*hat(w));
Rg_dd = simplify(Rg_d*hat(w) + Rg*hat(w_d));
Rg_ddd = simplify(Rg_dd*hat(w) + Rg_d*hat(w_d) + Rg_d*hat(w_d) + Rg*hat(w_dd));
Rg_dddd = simplify(Rg_ddd*hat(w) + Rg_dd*hat(w_d) + ...
2*(Rg_dd*hat(w_d) + Rg_d*hat(w_dd)) + ...
Rg_d*hat(w_dd) + Rg*hat(w_ddd));

syms angle1(t) angle2(t)

Ra = axisangle(e2,angle2(t))*axisangle(-e1,angle1(t));

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

     Ra = subs(Ra,differential,variant);
   Ra_d = subs(Ra_d,differential,variant);
  Ra_dd = subs(Ra_dd,differential,variant);

Rq = Rg*Ra;

Rq_d = Rg_d*Ra + Rg*Ra_d;

Rq_dd = Rg_dd*Ra + 2*Rg_d*Ra_d + Rg*Ra_dd;

Om_hat   = simplify(Rq.' * Rq_d);
Om_d_hat = simplify(Rq_d.' * Rq_d + Rq.' * Rq_dd);

Om = unhat(Om_hat);
Om_d = unhat(Om_d_hat);

% Are some of these useful? Unclear!

compute_Rq_all = matlabFunction(Rq);
compute_Om_all = matlabFunction(Om);
compute_Om_d_all = matlabFunction(Om_d);

global compute_Rq_state compute_Om_state compute_Om_d_state

compute_Rq_state = @(Rg,th1,th2) compute_Rq_all(Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),th1,th2);
compute_Om_state = @(Rg,th1,th2,th1d,th2d,w) compute_Om_all(Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),th1,th2,th1d,th2d,w(1),w(2),w(3));
compute_Om_d_state = @(Rg,th1,th2,th1d,th2d,th1dd,th2dd,w,w_d) compute_Om_d_all(Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),th1,th2,th1d,th2d,th1dd,th2dd,w(1),w(2),w(3),w_d(1),w_d(2),w_d(3));


%% Kinematics

disp('Computing translational dynamics...');

r = -Lg*Rg*e1; % vector from gripper center of mass to center of actuation

xg = xs + (Lg*mq)/(mg+mq)*Rg*e1; % gripper position using system enter of mass and gripper orientation
xs_dd = 1/(mq+mg)*(f*Rq*e3) - g*e3; % sum of external forces on entire system 
xg_dd = xs_dd + Lg*mq/(mq+mg)*Rg*(hat(w)^2+hat(w_d))*e1; % acceleration kinematics from differentiating (A) above

%% Dynamics

% Mr is reaction moments on the gripper in the gripper frame, at the center
% of actuation. Includes the effect of input torques on the two joints as
% well as torque due to the kinematic constraints.

Fr = mg*xg_dd + mg*g*e3; % solving sum of forces on gripper subsystem for the reaction force on the gripper in the world frame

% Solve angular dynamics of the two bodies for the moments between them
Mr_gripper = Rg*(Jg*w_d + cross(w,Jg*w))-cross(r,Fr);
Mr_quad = Rq*(Mq - cross(Om,Jq*Om) - Jq*Om_d);

% set these moments equal
equations = [
    Mr_quad - Mr_gripper; % == 0
    (Rq*e1).'*(Mr_gripper) - T1 % == 0;
    (Rg*e2).'*(Mr_gripper) - T2 % == 0;
];

disp('Arranging in form M(q) q'''' = a(q,q'') + B(q,q'') u');

acc = [w_d; th1dd; th2dd];
[M,affine] = equationsToMatrix(equations,acc);
[B,a] = equationsToMatrix(affine,u);
a = -a;

% M(x) x' = a(x) + B(x) u

%% Dynamics Solution 

% check that the dynamics are actually of the form necessary to use
% jacobian() etc to massage the equations like we did above

disp('Functionalizing...');

if ~eval(jacobian(M(:),acc) == 0)
   disp('Error : Dynamics are nonlinear in acceleration, cannot solve dynamics');
   return;
end

if ~eval(jacobian(B(:),u) == 0)
   disp('Error : B depends on u, cannot solve dynamics');
   return;
end

if eval(jacobian(a,u) == 0)
   disp('The system is control affine');
else
   disp('Error: a depends on u');
   return;
end

% turn these symbolic expressions into MATLAB functions that we can
% use with numeric inputs and outputs instead of symbolics

compute_M_all = matlabFunction(M);
compute_B_all = matlabFunction(B);
compute_a_all = matlabFunction(a);

%% Differential Flatness

disp('Computing differential flatness relationships')

syms Om_des1 Om_des2
syms Om_d_des1 Om_d_des2

[J1, j1] = equationsToMatrix([Om(1); Om(2)] == [Om_des1; Om_des2],[th1d; th2d]);
compute_thd_from_Om12_all = matlabFunction(inv(J1)*j1);

[J2, j2] = equationsToMatrix([Om_d(1); Om_d(2)] == [Om_d_des1;Om_d_des2],[th1dd;th2dd]);
compute_thdd_from_Om_d12_all = matlabFunction(inv(J2)*j2);



% 
% % Om_d_des = F \ d
% 
% xs_dd_des = sym('x_dd_des',[3 1]);
% xs_ddd_des = sym('x_ddd_des',[3 1]);
% xs_dddd_des = sym('x_dddd_des',[3 1]);
% 
% 
% thrust    = (Rq*e3).' * (mg+mq) * (xs_dd_des + g*e3);
% thrust_d  = (Rq*e3).' * (mg+mq) * xs_ddd_des;
% thrust_dd = (Rq*e3).' * (mg+mq) * xs_dddd_des - thrust*e3.'*hat(Om_des)^2*e3;
% 
% dynamics_thrust_constraint = ...
%      -(mg+mq) * xs_dddd_des + ...
%        thrust * Rq * ( hat(Om_d) + hat(Om)^2 ) * e3 + ...
%      thrust_d * 2 * Rq * hat(Om) * e3  + ...
%     thrust_dd * Rq * e3; % == 0
% 
% 
% 
% 
% 
% thrust_vectoring_constraints = [
%     dynamics_thrust_constraint(1);
%     dynamics_thrust_constraint(2);
% ];
% [F,d] = equationsToMatrix(thrust_vectoring_constraints == 0,[Om_d(1); Om_d(2)])
% 



% 
% F = simplify(jacobian(thrust_vectoring_constraints,Om_d));
% d = -simplify(expand(thrust_vectoring_constraints - F * Om_d));
% 
% compute_F_all = matlabFunction(F);
% compute_d_all = matlabFunction(d);
% 
% 
% quad_velocity_constraints = [
%     [e1.'; e2.'] * Om - 1 / (norm(x_dd_des+g*e3)) * [-(Rq*e2).'; (Rq*e1).']*x_ddd_des; % == 0 
% ];
% % 
% % % L Om = o
% % 
% L = simplify(jacobian(quad_velocity_constraints,Om));
% o = -simplify(expand(quad_velocity_constraints - L * Om));
% 
% compute_L_all = matlabFunction(L);
% compute_o_all = matlabFunction(o);

%% Euler Angle planning

% utilities for planning trajectories on SO(3) using Euler angles.

% trajectory_planning;

% order = 8;
% syms t_f
% C_R = sym('C_R',[order 1]);
% C_S = sym('C_S',[order 1]);
% C_W = sym('C_W',[order 1]);
% syms t 
% 
% roll_derivatives = compute_derivatives(C_R,t,t_f);
% swing_derivatives = compute_derivatives(C_S,t,t_f);
% wrist_derivatives = compute_derivatives(C_W,t,t_f);
% 
% Rg_des = axisangle(e2,roll_derivatives(1))*axisangle(e2,swing_derivatives(1))*axisangle(e1,wrist_derivatives(1));
% 
% w_d_des = diff(w_des,time);
% w_dd_des = diff(w_d_des,time);
% w_ddd_des = diff(w_dd_des,time);
% 
% desired_orientation = matlabFunction(Rg_des);
% angular_derivatives = matlabFunction([
%        w_des.';
%      w_d_des.';
%     w_dd_des.';
%    w_ddd_des.';
% ]);
% 
% compute_Rg_des = @(C_R, C_S, C_W, t, t_f) desired_orientation(C_R(1),C_R(2),C_R(3),C_R(4),C_R(5),C_R(6),C_R(7),C_R(8),C_S(1),C_S(2),C_S(3),C_S(4),C_S(5),C_S(6),C_S(7),C_S(8),C_W(1),C_W(2),C_W(3),C_W(4),C_W(5),C_W(6),C_W(7),C_W(8),t_f,t);
% compute_angulars_des = @(C_R, C_S, C_W, t, t_f) angular_derivatives(C_R(1),C_R(2),C_R(3),C_R(4),C_R(5),C_R(6),C_R(7),C_R(8),C_S(1),C_S(2),C_S(3),C_S(4),C_S(5),C_S(6),C_S(7),C_S(8),C_W(1),C_W(2),C_W(3),C_W(4),C_W(5),C_W(6),C_W(7),C_W(8),t_f,t);


trajectory_planning;

syms alpha(t) beta(t) gamma(t)

Rg_des = axisangle(e1,alpha(t))*axisangle(e2,beta(t))*axisangle(e1,gamma(t));
w_des = unhat(Rg_des.' * diff(Rg_des,t));
w_d_des = diff(w_des,t);
w_dd_des = diff(w_d_des,t);
w_ddd_des = diff(w_dd_des,t);

angular_block = [
       w_des.';
     w_d_des.';
    w_dd_des.';
   w_ddd_des.';
];

syms a ad add addd adddd
syms b bd bdd bddd bdddd
syms c cd cdd cddd cdddd

sym_derivatives = [
    alpha  diff(alpha,t)  diff(alpha,t,2)  diff(alpha,t,3)  diff(alpha,t,4) 
    beta  diff(beta,t)  diff(beta,t,2)  diff(beta,t,3)  diff(beta,t,4)
    gamma  diff(gamma,t)  diff(gamma,t,2)  diff(gamma,t,3)  diff(gamma,t,4)
];

symbol_derivatives = [
    a ad add addd adddd
	b bd bdd bddd bdddd
	c cd cdd cddd cdddd
];

substituted_block = simplify(subs(angular_block,sym_derivatives,symbol_derivatives));

global compute_angular_derivatives compute_Rg_angles
compute_angular_derivatives = matlabFunction(substituted_block);
compute_Rg_angles = matlabFunction(simplify(subs(Rg_des,sym_derivatives,symbol_derivatives)));


%% Numerical Physical Parameters

global Jqx_ Jqy_ Jqz_ Jgx_ Jgy_ Jgz_  
global Jq_ Jg_
global mq_ mg_ g_
global Lg_ Le_ Ls_

% Quad Inertia
Jqx_ = .001; % [kg*m^2]
Jqy_ = .001; % [kg*m^2]
Jqz_ = .001; % [kg*m^2]

% Gripper Inertia
Jgx_ = .001; % [kg*m^2]
Jgy_ = .001; % [kg*m^2]
Jgz_ = .001; % [kg*m^2]

mq_ = .5;  % [kg] quad mass
mg_ = .15; % [kg] gripper mass
Lg_ = .5;  % [m] distance from center of actuation to gripper COM
Le_ = .6;  % [m] distance from center of actuation to end effector

g_ = 9.81; % [m/s^2] acceleration due to gravity

Ls_ = (-(mg_*Lg_)/(mg_+mq_)+Le_);

Jq_ = diag([Jqx_ Jqy_ Jqz_]); % quad inertia
Jg_ = diag([Jgx_ Jgy_ Jgz_]); % gripper inertia

%% Function Parameter Substitution

% Substitute physical parameters of the robot into the functions we are
% creating to use with numerical inputs and outputs so that we don't have
% to supply these parameters every time we call the functions.

global compute_M_state compute_B_state compute_a_state
global compute_thd_from_Om12_state compute_thdd_from_Om_d12_state
% global compute_F_state compute_d_state compute_L_state compute_o_state

% dynamics
compute_M_state = @(Rg,th1,th2) compute_M_all(Jgx_,Jgy_,Jgz_,Jqx_,Jqy_,Jqz_,Lg_,Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),mg_,mq_,th1,th2);
compute_B_state = @(Rg,th1,th2) compute_B_all(Lg_,Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),mg_,mq_,th1,th2);
compute_a_state = @(Rg,w,th1,th2,th1d,th2d) compute_a_all(Jgx_,Jgy_,Jgz_,Jqx_,Jqy_,Jqz_,Lg_,Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),g_,mg_,mq_,th1,th2,th1d,th2d,w(1),w(2),w(3));

% differential flatness
compute_thd_from_Om12_state = @(Rg,w,th1,th2,Om_des1,Om_des2) compute_thd_from_Om12_all(Om_des1,Om_des2,Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),th1,th2,w(1),w(2),w(3));
compute_thdd_from_Om_d12_state = @(Rg,w,w_d,th1,th2,th1d,th2d,Om_d_des1,Om_d_des2) compute_thdd_from_Om_d12_all(Om_d_des1,Om_d_des2,Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),th1,th2,th1d,th2d,w(1),w(2),w(3),w_d(1),w_d(2),w_d(3));

