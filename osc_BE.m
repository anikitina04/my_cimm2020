%% Exercise 4.14
% Solve an oscillating system using BE scheme

omega = 2;
period = 2 * pi / omega;
dt = period / 20;
% dt = period / 2000; % gives better result
T = 3 * period;
t = 0:dt:T;

% Preallocation
u = zeros(1, length(t));
v = zeros(1, length(t));

% Initial conditions
X_0 = 2;
u(1) = X_0;
v(1) = 0;

for ii = 2:length(t)
    u(ii) = (dt*v(ii-1) + u(ii-1)) / (1 + (dt*omega)^2);
    v(ii) = (-dt*omega^2*u(ii-1) + v(ii-1)) / (1 + (dt*omega)^2);
end

plot(t, u, t, X_0*cos(omega*t))
legend("Numerical", "Exact")
