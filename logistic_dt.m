%% Exercise 4.4
% Find an appropriate timestep for population model using logistic model
% optimum dt = 0.01953125

clc
clear
close all

f = @(u, t) 0.1*(1 - u/500)*u;
U_0 = 100;

dt = 20; T = 100;
tspan = [dt, T];
[t_prev, u_prev] = ode_FE(f, tspan, U_0);

k = 1;

while 1
    dt_k = 2^(-k) * dt;
    tspan_cur = [dt_k, T];
    [t_cur, u_cur] = ode_FE(f, tspan_cur, U_0);

    plot_label = strcat("Smaller timestep: ", num2str(dt_k));
    plot(t_prev, u_prev, "b-")
    hold on
    plot(t_cur, u_cur, "r--")
    title(plot_label)
    xlabel("t"); ylabel("N(t)")
    legend("dt_k", "dt_k-1")
    hold off
    
    k = k+1;
    
    continue_loop = input("Do you want to try smaller timestep? (y/n)\n", "s");
    
    if continue_loop == "n"
        break
    else
        t_prev = t_cur;
        u_prev = u_cur;
    end
end
