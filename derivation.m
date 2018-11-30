clear;

%% Model Parameters

syms Jqx Jqy Jqz Jgx Jgy Jgz;
Jq = diag([Jqx Jqy Jqz]); % quad inertia
Jg = diag([Jgx Jgy Jgz]); % gripper inertia

syms mq; % quad mass
syms mg; % gripper mass
syms Lg; % distance from center of actuation to gripper COM

syms g;

%% State variables

% POSITION
xs = sym('xs',[3 1]); % position of center of mass
Rq = sym('Rq',[3 3]); % rotation from quad to world
Rg = sym('Rg',[3 3]); % rotation from gripper to world

% VELOCITY
xs_d = sym('xs_d',[3 1]); % velocity of center of mass
Om = sym('Om',[3 1]); % angular velocity of quadrotor body relative to world, expressed in quadrotor frame
w = sym('w',[3 1]); % angular velocity of gripper relative to world, expressed in gripper frame

% ACCELERATION (not truly part of state)
xs_dd = sym('xs_dd',[3 1]); % acceleration of center of mass
Om_d = sym('Om_d',[3 1]); % angular acceleration of quadrotor body relative to world, expressed in quadrotor frame
w_d = sym('w_d',[3 1]); % angular acceleration of gripper relative to world, expressed in gripper frame

%% Control Inputs

syms f % net thrust
syms M1 M2 M3 
Mq = [M1;M2;M3]; % rotor thrust moments around quadrotor body axes
syms T1 % torque at joint 1 - closer to quadrotor body, rotates around b1 
syms T2 % torque at joint 2 - closer to gripper body, rotates around g2 (which is rotated by joint 1 so is not constant in quadrotor body frame)

u = [f; M1; M2; M3; T1; T2];

%% Operators
global hat
hat = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];

%% Unit Vectors

global e1 e2 e3
e1 = [1; 0; 0];
e2 = [0; 1; 0];
e3 = [0; 0; 1];

%% Kinematics

r = -Lg*Rg*e1; % vector from gripper center of mass to center of actuation

xq = sym('xq',[3 1]); % position of quadrotor body center of mass
xg = sym('xg',[3 1]); % position of gripper center of mass

xs_exp = subs((mg*xg + mq*xq)/(mg+mq),xq,xg+r);
sol = solve(xs == xs_exp,xg);
xg = [sol.xg1;sol.xg2;sol.xg3]; % (A)

xs_dd = 1/(mq+mg)*(f*Rq*e3) - g*e3; % sum of external forces on entire system 
xg_dd = xs_dd + Lg*mq/(mq+mg)*Rg*(hat(w)^2+hat(w_d))*e1; % acceleration kinematics from differentiating (A) above

%% Dynamics

% Mr is reaction moments on the gripper in the gripper frame, at the center of actuation. Includes the effect of input torques on the two joints.

Fr = mg*xg_dd + mg*g*e3; % solving sum of forces on gripper subsystem for the reaction force on the gripper in the world frame

Mr_quad = Rg.'*Rq*(Mq - cross(Om,Jq*Om) - Jq*Om_d);
Mr_gripper = Jg*w_d + cross(w,Jg*w)-Rg.'*cross(r,Fr);

combination = Mr_quad - Mr_gripper; % == 0
acc = [Om_d; w_d];

top_3_rows = jacobian(combination, acc);
top_affine_part = -simplify(combination - top_3_rows*acc);

%% Kinematic Constraints

% OLD (INCORRECT) CONSTRAINTS
% kinematic_constraint_exp = Om_d.'*Rq.'*Rg*e1 + Om.'*(-hat(Om)*Rq.'*Rg+Rq.'*Rg*hat(w))*e1 - e1.'*w_d % == 0
% kinematic_constraint_exp = e2.'*(-hat(w)*Rg.'*Rq*hat(Om) + Rg.'*Rq.'*(hat(Om_d)+hat(Om)*hat(Om))-hat(w_d)*Rg.'*Rq+hat(w)*hat(w)*Rg.'*Rq-hat(w)*Rg.'*Rq*hat(Om))*e1 % == 0

% NEW CORRECT VELOCITY BASED KINEMATIC CONSTRAINT
kinematic_constraint_exp = (cross(Rq*e1,Rg*hat(w)*e2)+cross(Rq*hat(Om)*e1,Rg*e2)).'*(Rg*w-Rq*Om)+cross(Rq*e1,Rg*e2).'*( Rg*hat(w)*w - Rq*hat(Om)*Om + Rg*w_d - Rq*Om_d ); % == 0;

torque_constraint_1_exp = (Rg*Mr_gripper).'*Rq*e1 - T1; % == 0
torque_constraint_2_exp = (Mr_gripper).'*e2 - T2; % == 0

K = [
 kinematic_constraint_exp;
 torque_constraint_1_exp;
 torque_constraint_2_exp;
];

bottom_3_rows = jacobian(K, acc);
bottom_affine_part = -simplify(expand(K) - expand(bottom_3_rows*acc));

%% Dynamics Solution 
M = [
  top_3_rows;
  bottom_3_rows
];
affine = [
  top_affine_part;
  bottom_affine_part;
];
B = jacobian(affine,u);
a = simplify(affine - B*u);

% M(x) x' = a(x) + B(x) u

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

compute_M_all = matlabFunction(M);
compute_B_all = matlabFunction(B);
compute_a_all = matlabFunction(a);

%% Differential Flatness

x_dd_des = sym('x_dd_des',[3 1]);
x_ddd_des = sym('x_ddd_des',[3 1]);
x_dddd_des = sym('x_dddd_des',[3 1]);

thrust = (mg+mq)*e3.'*x_dd_des + g/(e3.'*(1/(mg+mq)*Rq*e3));
thrust_d = (Rq*e3).'*(mg+mq)*x_ddd_des;
thrust_dd = (Rq*e3).'*(mg+mq)*x_dddd_des - thrust*e3.'*hat(Om)^2*e3;

% Om_d = sym('Om_d_des',[3 1]);
dynamics_thrust_constraint = ...
    -(mg+mq)*x_dddd_des + ...
    thrust*(Rq*hat(Om_d)*e3+Rq*hat(Om)*hat(Om)*e3) + ...
    thrust_d*Rq*hat(Om)*e3 + ...
    thrust_dd*Rq*e3; % == 0

% F * Om_d_des - g = 0

thrust_vectoring_constraints = [
    kinematic_constraint_exp;
    dynamics_thrust_constraint
];

F = simplify(jacobian(thrust_vectoring_constraints,Om_d));
d = -simplify(expand(thrust_vectoring_constraints - F * Om_d));

compute_F_all = matlabFunction(F);
compute_d_all = matlabFunction(d);

%% Parameter Substitution

global Jqx_ Jqy_ Jqz_ Jgx_ Jgy_ Jgz_  mq_ mg_ Lg_ Le_ g_

% Quad Inertia
Jqx_ = .005; % [kg*m^2]
Jqy_ = .005; % [kg*m^2]
Jqz_ = .010; % [kg*m^2]

% Gripper Inertia
Jgx_ = .001; % [kg*m^2]
Jgy_ = .007; % [kg*m^2]
Jgz_ = .007; % [kg*m^2]

mq_ = .5;  % [kg] quad mass
mg_ = .35; % [kg] gripper mass
Lg_ = .5;  % [m] distance from center of actuation to gripper COM
Le_ = .6;  % [m] distance from center of actuation to end effector

g_ = 9.81; % [m/s^2] acceleration due to gravity

%% ode parameter substitution
compute_M_state = @(Rg,Rq) compute_M_all(Jgx_,Jgy_,Jgz_,Jqx_,Jqy_,Jqz_,Lg_,Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),Rq(1,1),Rq(1,2),Rq(1,3),Rq(2,1),Rq(2,2),Rq(2,3),Rq(3,1),Rq(3,2),Rq(3,3),mg_,mq_);
compute_B_state = @(Rg,Rq) compute_B_all(Lg_,Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),Rq(1,1),Rq(1,2),Rq(1,3),Rq(2,1),Rq(2,2),Rq(2,3),Rq(3,1),Rq(3,2),Rq(3,3),mg_,mq_);
compute_a_state = @(Rg,Rq,Om,w) compute_a_all(Jgx_,Jgy_,Jgz_,Jqx_,Jqy_,Jqz_,Lg_,Om(1),Om(2),Om(3),Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),Rq(1,1),Rq(1,2),Rq(1,3),Rq(2,1),Rq(2,2),Rq(2,3),Rq(3,1),Rq(3,2),Rq(3,3),mg_,mq_,w(1),w(2),w(3));
compute_F_state = @(Rg,Rq,x_dd_des) compute_F_all(Rg(1,2),Rg(2,2),Rg(3,2),Rq(1,1),Rq(1,2),Rq(1,3),Rq(2,1),Rq(2,2),Rq(2,3),Rq(3,1),Rq(3,2),Rq(3,3),g_,mg_,mq_,x_dd_des(3));
compute_d_state = @(Rg,Rq,Om,w,w_d,x_dd_des,x_ddd_des,x_dddd_des) compute_d_all(Om(1),Om(2),Om(3),Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),Rq(1,1),Rq(1,2),Rq(1,3),Rq(2,1),Rq(2,2),Rq(2,3),Rq(3,1),Rq(3,2),Rq(3,3),g_,mg_,mq_,w(1),w(2),w(3),w_d(1),w_d(3),x_dd_des(3),x_ddd_des(1),x_ddd_des(2),x_ddd_des(3),x_dddd_des(1),x_dddd_des(2),x_dddd_des(3));


% @(Om1,Om2,Om3,Rq1_1,Rq1_2,Rq1_3,Rq2_1,Rq2_2,Rq2_3,Rq3_1,Rq3_2,Rq3_3,g,mg,mq,x_dd_des3,x_ddd_des1,x_ddd_des2,x_ddd_des3,x_dddd_des1,x_dddd_des2,x_dddd_des3)

% compute_thrust_constraint_state = @(Rq,Om,Om_d_des,x_dd_des,x_ddd_des,x_dddd_des) compute_thrust_constraint_all(Om(1),Om(2),Om(3),Om_d_des(1),Om_d_des(2),Rq(1,1),Rq(1,2),Rq(1,3),Rq(2,1),Rq(2,2),Rq(2,3),Rq(3,1),Rq(3,2),Rq(3,3),g_,mg_,mq_,x_dd_des(3),x_ddd_des(1),x_ddd_des(2),x_ddd_des(3),x_dddd_des(1),x_dddd_des(2),x_dddd_des(3))

%% ode45 interface

velocity_cascade = zeros(3,(3+9+9+3+3+3));
velocity_cascade(1:3,3+9+9+(1:3)) = eye(3);

ode = @(x,u) [
  velocity_cascade*x; % first three rows
  reshape(sp(x(4:(4+8)),hat(x(25:27))),[9 1]);
  reshape(sp(x(13:(13+8)),hat(x(28:30))),[9 1]);
  u(1)/(mg_+mq_)*sp(x(4:(4+8)),e3) - g_*e3;
  compute_M_state(sp(x(13:13+8),eye(3)),sp(x(4:4+8),eye(3))) \ (compute_B_state(sp(x(13:13+8),eye(3)),sp(x(4:4+8),eye(3)))* u + compute_a_state(sp(x(13:13+8),eye(3)),sp(x(4:4+8),eye(3)),x(25:27),x(28:30)))
];

% M(x) [Om_d w_d] = (B(x) u + a(x))

% x = [
%   x_s
%   Rq
%   Rg
%   x_s_d
%   Om
%   w
% ]
% u = [
%   f
%   M1
%   M2
%   M3
%   T_1
%   T_2
% ]

%% Differential Flatness
