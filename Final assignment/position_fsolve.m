function [T, Q] = position_fsolve(C_fun, tend, q_0, dt)

N_t = floor(round(tend/dt));
T = linspace(0, N_t*dt, N_t+1)';
Q = zeros(N_t+1, length(q_0));

Q(1, :) = fsolve(@(q) C_fun(0, q), q_0)';
% Step equations forward in time
for n = 1 : N_t
    Q(n + 1, :) = fsolve(@(q) C_fun(T(n + 1), q), Q(n, :)')';
end