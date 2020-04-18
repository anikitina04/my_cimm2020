function Ct = driving_joint_dt(d_k_t, t)
% i - body id
% k = 1, 2, 3 for x, y and phi
% d_k - function of time that define coordinate value
% q - coordinate vector

Ct = -d_k_t(t);