%% Exercise 4.10
% A script for caclulating sum of potential and kinetic energies of the
% system using Forward Euler method.

omega = 2;
period = 2 * pi / omega;
dt = period / 20;
T = 3 * period;
t = 0:dt:T;

% Preallocation
u = zeros(1, length(t));
v = zeros(1, length(t));

% Initial conditions
X_0 = 2;
u(1) = X_0;
v(1) = 0;

for ii = 1:length(t)-1
   u(ii+1) = u(ii) + dt.*v(ii);
   v(ii+1) = v(ii) - dt.*omega.^2.*u(ii);
end

[U, K] = osc_energy(u, v, omega);
plot(t, U+K)
xlabel("Time")
ylabel("Sum of potential and kinetic energies")
title("Sum of potential and kinetic energis using FE method")
