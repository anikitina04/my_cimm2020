close all
clc
clear

% Four-bar linkage kinematic analysis
end_time = 10;
time_step = 0.001;
l_ground = [0.4; 0.2];
l_crank = 0.2;
l_coupler = 0.5;
l_follower = 0.4;
phi_2 = 0;

%% Coordinates
% ground
q0 = [0; 0; 0];

% crank
q1 = zeros(3,1);
q1(1) = -0.5*l_crank*sin(phi_2);
q1(2) = 0.5*l_crank*cos(phi_2);
q1(3) = phi_2;

% coupler
alpha = phi_2 + pi/2 - atan(l_ground(1)/l_ground(2));
c = sqrt((l_ground(1)^2) + (l_ground(2)^2));
e = sqrt(l_crank^2 + c^2 - (2*l_crank*c*cos(alpha)));
beta1 = acos((l_crank^2+e^2-c^2) / (2*l_crank*e));
beta2 = acos((l_coupler^2+e^2-l_follower^2) / (2*l_coupler*e));
beta = beta1 + beta2;
q2 = zeros(3, 1);
q2(3) = phi_2 + beta - (pi/2);
q2(1) = -l_crank*sin(phi_2) + 0.5*l_coupler*cos(q2(3));
q2(2) = l_crank*cos(phi_2) + 0.5*l_coupler*sin(q2(3));

% follower
delta = acos((l_follower^2+l_coupler^2-e^2) / (2*l_follower*l_coupler));
q3 = zeros(3,1);
q3(3) = q2(3) + delta - pi/2;
q3(1) = q2(1) + 0.5*l_coupler*cos(q2(3)) + 0.5*l_follower*sin(q3(3));
q3(2) = q2(2) + 0.5*l_coupler*sin(q2(3)) + 0.5*l_follower*cos(q3(3));

% ground 2
q4 = [l_ground; 0];

q_0 = [q0; q1; q2; q3; q4]; % initial coordinates
qp_0 = zeros(length(q_0), 1);

%% We need two constraint types (geometric ones)
% - revolute
% - simple constraints

%% Revolute joints
% 1 connects ground (frame) and crank
revolute(1).i = 1;
revolute(1).j = 2;
revolute(1).s_i = [0; 0];
revolute(1).s_j = [0; -0.1];

% 2 connects crank and coupler
revolute(2).i = 2;
revolute(2).j = 3;
revolute(2).s_i = [0; 0.1];
revolute(2).s_j = [-0.25; 0];

% 3 connects coupler and follower
revolute(3).i = 3;
revolute(3).j = 4;
revolute(3).s_i = [0.25; 0];
revolute(3).s_j = [0; 0.2];

% 4 connects follower and ground (frame)
revolute(4).i = 4;
revolute(4).j = 5;
revolute(4).s_i = [0; -0.2];
revolute(4).s_j = [0; 0];

%% Simple constraints

% Three simple joints to fix left side
simple(1).i = 1;
simple(1).k = 1;
simple(1).c_k = 0;

simple(2).i = 1;
simple(2).k = 2;
simple(2).c_k = 0;

simple(3).i = 1;
simple(3).k = 3;
simple(3).c_k = 0;

% Three simple joints to fix right side
simple(4).i = 5;
simple(4).k = 1;
simple(4).c_k = l_ground(1);

simple(5).i = 5;
simple(5).k = 2;
simple(5).c_k = l_ground(2);

simple(6).i = 5;
simple(6).k = 3;
simple(6).c_k = 0;

%% Add driving constraints
driving.i = 2;
driving.k = 3;
driving.d_k = @(t,q) pi/9 - 1.5 * t;
driving.d_k_t = @(t,q) -1.5;
driving.d_k_tt = @(t,q) 0;

%% Add gravity to each body
body(1).m = 1;
body(2).m = 2;
body(3).m = 3;
body(1).Ic = body(1).m * 0.2^2 / 12; % mass moment of inertia along center of mass in kgm2
body(2).Ic = body(2).m * 0.5^2 / 12;
body(3).Ic = body(3).m * 0.2^2 / 12;
M = mass_matrix(body);
grav = [0; -9.81]; % gravitational acceleration
% q0 = system_coordinates(body);
F = @(t, q, qd) force_vector(grav, sforce, body, q_0);

%% Solve constraint equation using NR for position and velocity
C_fun = @(t, q) constraint(revolute, simple, driving, t, q);
Cq_fun = @(t, q) constraint_dq(revolute, simple, driving, t, q);
Ct_fun = @(t, q) constraint_dt(revolute, simple, driving, t, q);
Ctt_fun = @(t, q, qd) constraint_dtt(revolute, simple, driving, t, q, qd);
[T, Q, QP, QA] = pos_vel_acc_NR(C_fun, Cq_fun, Ct_fun, Ctt_fun, end_time, q_0, time_step);

%% Dynamic solution
% alpha = 10;
% beta = 10;
% 
% LHS = @(t, q) [M, Cq_fun(t, q)'; Cq_fun(t, q), zeros(size(Cq_fun(t,q)))];
% RHS = @(t, q, qd) [F(t, q, qd); Ctt_fun(t, q, qd) - 2*alpha*Ct_fun(t, q) - beta^2 * C_fun(t, q)];
% acc_fun = @(t, q, qd) LHS(t, q)\RHS(t, q, qd);
% acc_f_ode = @(t,y) [y(length(M)+1:end) ;...
%                     acc_f(t, y(1:length(M)), y(length(M)+1:end))];
% opts = odeset('AbsTol', 1e-9, 'RelTol', 1e-6);
% y0 = [q_0,qp_0];
% disp('start ode45 solver')
% [t,solution_ode] = ode45(acc_f_ode,[0,end_time],y0,opts);


%% Some verification plots
figure
plot(Q(:, 4), Q(:, 5), ...
    Q(:, 7), Q(:, 8), ...
    Q(:, 10), Q(:, 11), ...
    0, 0, '*', 'LineWidth', 2);
axis equal
legend("Frame", "Coupler", "Follower")
title("Position")

%% Some verification plots
figure
plot(QP(:, 4), QP(:, 5), ...
    QP(:, 7), QP(:, 8), ...
    QP(:, 10), QP(:, 11), ...
    0, 0, '*', 'LineWidth', 2);
axis equal
legend("Frame", "Coupler", "Follower")
title("Velocity")

%% Some verification plots
figure
plot(QA(:, 4), QA(:, 5), ...
    QA(:, 7), QA(:, 8), ...
    QA(:, 10), QA(:, 11), ...
    0, 0, '*', 'LineWidth', 2);
axis equal
legend("Frame", "Coupler", "Follower")
title("Acceleration")
