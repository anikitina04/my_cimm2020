% Assignement D
clear
close all

a = 0.1; % m
b = 0.2; % m
omega = 1; % rad/s, also derivative(phi)
phi = @(t) pi/6 + omega.*t;
t = 0:0.001:4*pi;
x = zeros(2, length(t));
x_vel = zeros(2, length(t));
eps = 1e-8;

for n = 1:length(t)
fun = @(x) constr(x, a, b, phi(t(n)));
jacobian = @(x) jacob(x, b);

x(:,n) = newton_raphson (fun, jacobian, x(:,n), eps);

x0 = x(:,n);

fun_vel = @(x) constr_vel(x, a, b, omega, x0, phi(t(n)));

jacobian_vel = @(x) jacob(x0, b);

x_vel(:,n) = newton_raphson (fun_vel, jacobian_vel, x0, eps);

end

figure
plot(t, x(1, :), t, x(2, :))
legend("Angle theta", "Distance D")
figure
plot(t, x_vel(1, :), t, x_vel(2, :))
legend("Rotational velocity", "Velocity of point B")

function result = constr(x0, a, b, phi)
    theta = x0(1);
    d = x0(2);
    result = [a.*cos(phi) + b.*cos(theta) - d;
              a.*sin(phi) - b*sin(theta)];
end

function result = constr_vel(x0_vel, a, b, omega, x0, phi)
    theta_d = x0_vel(1);
    d_d = x0_vel(2);
    theta = x0(1);
    result = [-a.*omega.*sin(phi) - b.*sin(theta).*theta_d - d_d;
              a.*omega.*cos(phi) - b*cos(theta).*theta_d];
end

function result = jacob(u, b)
    theta = u(1);
    result = [-b.*sin(theta), -1;
              -b.*cos(theta), 0];
end