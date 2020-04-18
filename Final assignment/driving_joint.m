function C = driving_joint(i, k, d_k, t, q)
% i - body id
% k = 1, 2, 3 for x, y and phi
% d_k - function of time that define coordinate value
% q - coordinate vector

idx_i = body_idx(i);

C = q(idx_i(k)) - d_k(t);