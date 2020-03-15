function [t, u] = ode_FE(f, tspan, U_0)
    dt = tspan(1);
    T = tspan(2);
    N_t = floor(T/dt);
    u = zeros(N_t+1, length(U_0));
    t = linspace(0, N_t*dt, length(u))';
    u(1, :) = U_0';
    for n = 1:N_t
        u(n+1, :) = u(n, :) + dt.*f(u(n, :), t(n))';
    end
end