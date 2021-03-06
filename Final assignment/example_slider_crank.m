close all

% Slider crank kinematic analysis
%% Coordinates
% ground
q1 = [0; 0; 0];
% crank
q2 = [-0.1 * cosd(30)
    0.1 * sind(30)
    -deg2rad(30)];
% link
h_B = 0.2 * sind(30); % y coordinate of point B
phi_l = asin(h_B / 0.5); % link's angle
q3 = [-0.2 * cosd(30) - 0.3 * cos(phi_l)
    h_B - 0.3 * sin(phi_l)
    phi_l];
% slider
q4 = [-0.2 * cosd(30) - 0.5 * cos(phi_l)
    0
    0];

q_0 = [q1; q2; q3; q4]; % initial coordinates
qp_0 = zeros(length(q_0), 1);

%% We need two constraint types (geometric ones)
% - revolute
% - simple constraints

%% Revolute joints
% 1 connects ground and crank
revolute(1).i = 1;
revolute(1).j = 2;
revolute(1).s_i = [0; 0];
revolute(1).s_j = [0.1; 0];

% 2 connects crank and link
revolute(2).i = 2;
revolute(2).j = 3;
revolute(2).s_i = [-0.1; 0];
revolute(2).s_j = [0.3; 0];

% 3 connects link and slider
revolute(3).i = 3;
revolute(3).j = 4;
revolute(3).s_i = [-0.2; 0];
revolute(3).s_j = [0; 0];

% % Check revolute joint constraints
% r = revolute(3);
% C_r_i = revolute_joint(r.i, r.j, r.s_i, r.s_j, q_0)

%% Simple constraints

% Three simple joints to fix the ground origin
simple(1).i = 1;
simple(1).k = 1;
simple(1).c_k = 0;

simple(2).i = 1;
simple(2).k = 2;
simple(2).c_k = 0;

simple(3).i = 1;
simple(3).k = 3;
simple(3).c_k = 0;

% slider - use simple joints instead of translational
simple(4).i = 4;
simple(4).k = 2;
simple(4).c_k = 0;

simple(5).i = 4;
simple(5).k = 3;
simple(5).c_k = 0;

% % check simple constraints
% for s = simple
%     C_s_i = simple_joint(s.i, s.k, s.c_k, q_0)
% end

%% Add some driving constraints
driving.i = 2;
driving.k = 3;
driving.d_k = @(t) -pi/6 - 1.2 * t;
driving.d_k_t = @(t) -1.2;
driving.d_k_tt = @(t) 0;

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
end_time = 10;
time_step = 0.001;
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

%% Some verification plots
figure
plot(QP(:, 4), QP(:, 5), ...
    QP(:, 7), QP(:, 8), ...
    QP(:, 10), QP(:, 11), ...
    0, 0, '*', 'LineWidth', 2);
axis equal

%% Some verification plots
figure
plot(QA(:, 4), QA(:, 5), ...
    QA(:, 7), QA(:, 8), ...
    QA(:, 10), QA(:, 11), ...
    0, 0, '*', 'LineWidth', 2);
axis equal