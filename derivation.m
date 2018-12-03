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
Rq = sym('Rq',[3 3]); % rotation from quad to world
Rg = sym('Rg',[3 3]); % rotation from gripper to world

% VELOCITY
xs_d = sym('xs_d',[3 1]); % velocity of center of mass
Om = sym('Om',[3 1]); % angular velocity of quadrotor body relative to world, expressed in quadrotor frame
w = sym('w',[3 1]); % angular velocity of gripper relative to world, expressed in gripper frame

% ACCELERATION (not truly part of state but used below)
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
global hat unhat
hat = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
unhat = @(v_hat) [v_hat(3,2); v_hat(1,3); v_hat(2,1)];

%% Unit Vectors

global e1 e2 e3
e1 = [1; 0; 0];
e2 = [0; 1; 0];
e3 = [0; 0; 1];

%% Kinematics

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
Mr_quad = Rg.'*Rq*(Mq - cross(Om,Jq*Om) - Jq*Om_d);
Mr_gripper = Jg*w_d + cross(w,Jg*w)-Rg.'*cross(r,Fr);

% set these moments equal
combination = Mr_quad - Mr_gripper; % == 0

% massage this equation into the form Ax=b
acc = [Om_d; w_d];
top_3_rows = jacobian(combination, acc);
top_affine_part = -simplify(combination - top_3_rows*acc);

%% Kinematic Constraints

% constraint on accelerations found by differentiating the constraint on
% angular velocities due to missing degree of freedom between the bodies
kinematic_constraint_exp = (cross(Rq*e1,Rg*hat(w)*e2)+cross(Rq*hat(Om)*e1,Rg*e2)).'*(Rg*w-Rq*Om)+cross(Rq*e1,Rg*e2).'*( Rg*hat(w)*w - Rq*hat(Om)*Om + Rg*w_d - Rq*Om_d ); % == 0;

% torque between bodies projected onto actuator axes equals the torque
% applied by those actuators
torque_constraint_1_exp = (Rg*Mr_quad).'*Rq*e1 - T1; % == 0
torque_constraint_2_exp = (Mr_quad).'*e2 - T2; % == 0

K = [
 kinematic_constraint_exp;
 torque_constraint_1_exp;
 torque_constraint_2_exp;
];

% massage this equation into the form Ax=b
bottom_3_rows = jacobian(K, acc);
bottom_affine_part = -simplify(expand(K) - expand(bottom_3_rows*acc));

%% Dynamics Solution 

% Combine the results from above into one vector equation that gives six
% different constraints on the angular accelerations, which should allow us
% to solve for them.

% M(x) x' = a(x) + B(x) u

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

% check that the dynamics are actually of the form necessary to use
% jacobian() etc to massage the equations like we did above

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

x_dd_des = sym('x_dd_des',[3 1]);
x_ddd_des = sym('x_ddd_des',[3 1]);
x_dddd_des = sym('x_dddd_des',[3 1]);

thrust    = (Rq*e3).' * (mg+mq) * (x_dd_des + g*e3);
thrust_d  = (Rq*e3).' * (mg+mq) * x_ddd_des;
thrust_dd = (Rq*e3).' * (mg+mq) * x_dddd_des - thrust*e3.'*hat(Om)^2*e3;

dynamics_thrust_constraint = ...
     -(mg+mq) * x_dddd_des + ...
       thrust * Rq * ( hat(Om_d) + hat(Om)^2 ) * e3 + ...
     thrust_d * 2 * Rq * hat(Om) * e3  + ...
    thrust_dd * Rq * e3; % == 0

% F * Om_d_des - g = 0

thrust_vectoring_constraints = [
    kinematic_constraint_exp;
    dynamics_thrust_constraint
];

F = simplify(jacobian(thrust_vectoring_constraints,Om_d));
d = -simplify(expand(thrust_vectoring_constraints - F * Om_d));

compute_F_all = matlabFunction(F);
compute_d_all = matlabFunction(d);

%% Numerical Physical Parameters

global Jqx_ Jqy_ Jqz_ Jgx_ Jgy_ Jgz_  mq_ mg_ Lg_ Le_ g_ Ls_

% Quad Inertia
Jqx_ = .001; % [kg*m^2]
Jqy_ = .001; % [kg*m^2]
Jqz_ = .001; % [kg*m^2]

% Gripper Inertia
Jgx_ = .001; % [kg*m^2]
Jgy_ = .001; % [kg*m^2]
Jgz_ = .001; % [kg*m^2]

mq_ = .5;  % [kg] quad mass
mg_ = .5; % [kg] gripper mass
Lg_ = .5;  % [m] distance from center of actuation to gripper COM
Le_ = .6;  % [m] distance from center of actuation to end effector

g_ = 10; % [m/s^2] acceleration due to gravity

Ls_ = (-(mg_*Lg_)/(mg_+mq_)+Le_);

%% Function Parameter Substitution

% Substitute physical parameters of the robot into the functions we are
% creating to use with numerical inputs and outputs so that we don't have
% to supply these parameters every time we call the functions.

% dynamics
compute_M_state = @(Rg,Rq) compute_M_all(Jgx_,Jgy_,Jgz_,Jqx_,Jqy_,Jqz_,Lg_,Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),Rq(1,1),Rq(1,2),Rq(1,3),Rq(2,1),Rq(2,2),Rq(2,3),Rq(3,1),Rq(3,2),Rq(3,3),mg_,mq_);
compute_B_state = @(Rg,Rq) compute_B_all(Lg_,Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),Rq(1,1),Rq(1,2),Rq(1,3),Rq(2,1),Rq(2,2),Rq(2,3),Rq(3,1),Rq(3,2),Rq(3,3),mg_,mq_);
compute_a_state = @(Rg,Rq,Om,w) compute_a_all(Jgx_,Jgy_,Jgz_,Jqx_,Jqy_,Jqz_,Lg_,Om(1),Om(2),Om(3),Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),Rq(1,1),Rq(1,2),Rq(1,3),Rq(2,1),Rq(2,2),Rq(2,3),Rq(3,1),Rq(3,2),Rq(3,3),mg_,mq_,w(1),w(2),w(3));

% differential flatness
compute_F_state = @(Rg,Rq,x_dd_des) compute_F_all(Rg(1,2),Rg(2,2),Rg(3,2),Rq(1,1),Rq(1,2),Rq(1,3),Rq(2,1),Rq(2,2),Rq(2,3),Rq(3,1),Rq(3,2),Rq(3,3),g_,mg_,mq_,x_dd_des(1),x_dd_des(2),x_dd_des(3));
compute_d_state = @(Rg,Rq,Om,w,w_d,x_dd_des,x_ddd_des,x_dddd_des) compute_d_all(Om(1),Om(2),Om(3),Rg(1,1),Rg(1,2),Rg(1,3),Rg(2,1),Rg(2,2),Rg(2,3),Rg(3,1),Rg(3,2),Rg(3,3),Rq(1,1),Rq(1,2),Rq(1,3),Rq(2,1),Rq(2,2),Rq(2,3),Rq(3,1),Rq(3,2),Rq(3,3),g_,mg_,mq_,w(1),w(2),w(3),w_d(1),w_d(3),x_dd_des(1),x_dd_des(2),x_dd_des(3),x_ddd_des(1),x_ddd_des(2),x_ddd_des(3),x_dddd_des(1),x_dddd_des(2),x_dddd_des(3));

%% ode45 Interface

% Take the functions we've generated above, which deal only with the
% angular accelerations, and combine them with the translational dynamics
% and the position-velocity integration relationships to create the full
% dynamics of the robot as a system of first order ODE's in order to
% provide to MATLAB ode solver.

% derivative of center of mass position is center of mass velocity
velocity_cascade = zeros(3,(3+9+9+3+3+3));
velocity_cascade(1:3,3+9+9+(1:3)) = eye(3);

ode = @(x,u) [
  velocity_cascade*x; % first three rows
  reshape(sp(x(4:(4+8)),hat(x(25:27))),[9 1]); % d/dt Rq = Rq Om_hat
  reshape(sp(x(13:(13+8)),hat(x(28:30))),[9 1]); % d/dt Rg = Rg w_hat
  u(1)/(mg_+mq_)*sp(x(4:(4+8)),e3) - g_*e3; % acceleration of system center of mass is due to gravity and thrust
  compute_M_state(sp(x(13:13+8),eye(3)),sp(x(4:4+8),eye(3))) \ (compute_B_state(sp(x(13:13+8),eye(3)),sp(x(4:4+8),eye(3)))* u + compute_a_state(sp(x(13:13+8),eye(3)),sp(x(4:4+8),eye(3)),x(25:27),x(28:30)))
];

% x = [x_s; Rq; Rg; xs_d; Om; w]