function [T, Q] = position_NR(C_fun, Cq_fun, tend, q_0, dt)

N_t = floor(round(tend/dt));
T = linspace(0, N_t*dt, N_t+1)';
Q = zeros(N_t+1, length(q_0));

% [x, iteration_counter] = NR_method(F, J, x, eps)

[qi, icnt] = NR_method(@(q) C_fun(0, q), @(q) Cq_fun(0, q), q_0, 1e-8);
Q(1, :) = qi';
% Step equations forward in time
for n = 1 : N_t
    [qi, icnt] = NR_method(@(q) C_fun(T(n + 1), q), ...
    @(q) Cq_fun(T(n + 1), q), ...
    qi, 1e-8);
    disp(icnt)
    Q(n + 1, :) = qi';
end